"""
Data mixin for GenBankBrowser.

Handles data operations including:
- Loading GenBank files
- Project save/load (.bcproject format)
- Record sorting and filtering
- State reset
"""
import tkinter as tk
from tkinter import filedialog
import os
import json
import zipfile
import shutil
import tempfile
import logging
from config.constants import *
from core.data_manager import load_genbank_files

logger = logging.getLogger(__name__)

class DataMixin:
    """
    Handles data loading, state management, and sorting for GenBankBrowser.
    """

    def open_files(self):
        """Opens file dialog to load GenBank files (Appends to current session)."""
        files = filedialog.askopenfilenames(title="Select GenBank Files", 
                                            filetypes=[("GenBank", "*.gb *.gbk"), ("All Files", "*.*")])
        if not files: return

        self.status_label.config(text="Loading...")
        self.update_idletasks()
        
        # Note: We do NOT call _reset_app_state() here anymore, to allow appending.
        # If the user wants to clear, they should restart or we can add a "Clear" button later.

        # Load NEW files
        new_records, new_cache, new_pmap, new_span = load_genbank_files(files, self.next_import_id)
        self.next_import_id += 1  # Increment for next import
        
        if not new_records:
            self.status_label.config(text=f"No valid clusters found in selection.")
            return

        # Merge with existing data
        self.gbk_records.extend(new_records)
        self.cached_cds_data.update(new_cache)
        self.protein_map.update(new_pmap)
        
        if new_span > self.max_cluster_span:
            self.max_cluster_span = new_span
        if self.max_cluster_span < 1: self.max_cluster_span = 1

        self.status_label.config(text=f"Loaded {len(self.gbk_records)} clusters (Added {len(new_records)}).")
        
        self._update_sorted_cache()
        self._update_layout_map()
        self.update_button_states()
        self._full_redraw_sequence()

    def _reset_app_state(self):
        """Clears all data and state variables for a fresh load."""
        # Cleanup previous temp dir if it exists
        if hasattr(self, 'project') and self.project.temp_dir:
            try:
                shutil.rmtree(self.project.temp_dir)
            except Exception as e:
                logger.warning(f"Failed to cleanup temp dir: {e}")
            self.project.temp_dir = None

        # Clear data storage
        self.gbk_records.clear()
        self.sorted_records_cache.clear()
        self.protein_map.clear()
        self.cached_cds_data.clear()
        self.visible_canvas_tag_map.clear()
        
        # Clear grouped state using dataclass methods
        self.selection.clear()
        self.alignment.clear()
        self.view.clear()
        self.clipboard.clear()
        self.metadata.clear()
        
        # Sync clipboard UI with state (hide frame, reset button)
        self.clipboard_content.pack_forget()
        self.btn_clipboard_toggle.config(text="+")
        
        # Reset drag state
        self.drag_data = None
        
        # Clear project path
        self.project.clear()
        
        # Redraw canvases to show empty state
        self._draw_clipboard()
        self._draw_special_row()
        self.org_canvas.delete("all")
        self.gene_canvas.delete("all")
        self.status_label.config(text="No files loaded.")
        
        # Reset cluster count label
        self.cluster_count_label.config(text="")
        
        # Reset entry fields to defaults
        self.core_cutoff_entry.delete(0, tk.END)
        self.core_cutoff_entry.insert(0, "40")
        self.cutoff_entry.delete(0, tk.END)
        self.cutoff_entry.insert(0, "40")
        self.bench_cutoff_entry.delete(0, tk.END)
        self.bench_cutoff_entry.insert(0, "40")
        self.limit_entry.delete(0, tk.END)
        self.limit_entry.insert(0, "8")
        
        # Reset checkboxes to defaults
        self.var_hide_genes.set(True)
        self.var_highlight_genes.set(True)
        self.var_limit_range.set(True)
        self.var_inspect_monochrome.set(True)
        self.var_show_identity.set(True)
        self.var_color_by_gene.set(True)
        
        # Reset zoom to default
        self.zoom_scale.set(DEFAULT_ZOOM)
        
        # Reset filter text
        if hasattr(self, 'filter_entry'):
            self.filter_entry.delete(0, tk.END)

    def _update_sorted_cache(self):
        """Sorts records: Favorites -> Then based on Sort Mode."""
        # Initial list (exclude pinned)
        records = [r for r in self.gbk_records if r is not self.pinned_record]

        # Filter by text input (Organism Name)
        query = self.filter_text.get().lower().strip()
        if query:
            # Keep only records where organism name contains query
            records = [r for r in records if query in r.annotations.get('organism', r.id).lower()]
        
        favs = [r for r in records if r.uid in self.favorite_uids]
        others = [r for r in records if r.uid not in self.favorite_uids]
        
        if self.sort_mode == "alpha":
             others.sort(key=lambda r: r.annotations.get('organism', r.id))
        elif self.sort_mode == "similarity":
             others.sort(key=lambda r: self.cluster_similarity_scores.get(r.uid, 0), reverse=True)

        self.sorted_records_cache = favs + others
        
        # Update cluster count in header
        self._update_cluster_count()

    def on_filter_change(self, *args):
        """Callback for filter text change."""
        self._update_sorted_cache()
        self._update_layout_map()
        self._full_redraw_sequence()

    def _update_cluster_count(self):
        """Updates the cluster count label in the Gene clusters header."""
        count = len(self.sorted_records_cache)
        self.cluster_count_label.config(text=f" ({count})")

    def _update_layout_map(self):
        """Calculates Y-position and height for every row (handling expansions)."""
        self.row_layout_map = []
        current_y = 0.0
        scale = self.zoom_scale.get() / 1000.0

        limit_count = float('inf')
        if self.var_limit_range.get():
            try: limit_count = int(self.limit_entry.get())
            except ValueError: pass

        for i, rec in enumerate(self.sorted_records_cache):
            h = DEFAULT_ROW_HEIGHT
            
            # If Expanded, calculate extra height required for labels
            if rec.uid in self.expanded_rows:
                h = self._calculate_expanded_height(rec, scale, limit_count)

            self.row_layout_map.append({'y': current_y, 'h': h})
            current_y += h

    def _calculate_expanded_height(self, rec, scale, limit_count):
        """Determines height of an expanded row based on stacked labels."""
        cds_list = self.cached_cds_data.get(rec.uid, [])
        levels = [] 
        is_flipped = rec.uid in self.flipped_rows
        
        min_start = min((c['start'] for c in cds_list), default=0)
        max_end = max((c['end'] for c in cds_list), default=0)
        span_bp = max_end - min_start
        
        # Filter based on Core range if applicable
        indices_to_render = self._get_render_indices(cds_list, limit_count)

        # 1. Collect all labels first
        labels_to_stack = []

        for idx in indices_to_render:
            cds = cds_list[idx]
            prod = cds['product']

            # If 'Hide other' is on and gene is in list, don't allocate space
            if self.var_hide_genes.get() and any(k in prod.lower() for k in self.hide_genes_list):
                continue

            # Aligned genes (hits) are drawn as WHITE blocks without text labels in expanded view.
            hit = self.current_alignment_results.get(cds['unique_tag'])
            if hit: continue 
            
            # Calculate visual position to detect collisions
            start, end = cds['start'], cds['end']
            
            if is_flipped:
                rel_s = span_bp - (end - min_start)
            else:
                rel_s = start - min_start
            
            # Font calculation needs to match drawing logic
            is_bold = self.var_highlight_genes.get() and any(k in prod.lower() for k in self.highlight_genes_list)
            font_obj = self.label_font_bold if is_bold else self.label_font
            txt_w = font_obj.measure(prod) + 6 
            
            # Match Renderer Logic: Labels are centered on the gene
            gene_width = (end - start) * scale
            x_start_gene = rel_s * scale
            x_center = x_start_gene + (gene_width / 2)
            
            label_x1 = x_center - (txt_w / 2)
            label_x2 = x_center + (txt_w / 2)
            
            labels_to_stack.append((label_x1, label_x2))

        # 2. Sort by visual x_start for Left-to-Right stacking
        labels_to_stack.sort(key=lambda x: x[0])

        max_level = 0
        for x_start, x_end in labels_to_stack:
            found_lvl = -1
            for lvl_idx, lvl_end in enumerate(levels):
                if x_start > lvl_end + 5: 
                    found_lvl = lvl_idx
                    levels[lvl_idx] = x_end
                    break
            
            if found_lvl == -1:
                if len(levels) < MAX_LABEL_STACK:
                    found_lvl = len(levels)
                    levels.append(x_end)
                else:
                    continue 
            
            if found_lvl > max_level: max_level = found_lvl
        
        extra_h = (max_level + 1) * LABEL_STEP_HEIGHT if labels_to_stack else 0
        # Dynamic height: y_mid is fixed at bottom (h-45). Labels stack up from y_mid-15.
        # Top = h - 45 - 15 - extra_h - 12 (font). We want Top >= 10.
        # h >= extra_h + 82. Using 85 for safety.
        return max(DEFAULT_ROW_HEIGHT, extra_h + 85)

    def _get_render_indices(self, cds_list, limit_count):
        return self.renderer._get_render_indices(cds_list, limit_count)

    def set_sort_mode(self, mode):
        self.sort_mode = mode
        self._update_sorted_cache()
        self._update_layout_map()
        self._full_redraw_sequence()

    # --- Project Save/Load ---

    def save_project_as(self):
        """Opens dialog to save project to a new file."""
        filepath = filedialog.asksaveasfilename(
            title="Save Project As",
            defaultextension=".bcproject",
            filetypes=[("BisCEET Project", "*.bcproject"), ("All Files", "*.*")]
        )
        if not filepath: return
        self.save_project(filepath)

    def save_project(self, filepath=None):
        """Saves current state to a bundled ZIP project file.
        
        If filepath is provided, saves to that location.
        If no filepath but current_project_path exists, does a quick save.
        Otherwise prompts for a new location (first-time save).
        """
        if not filepath:
            if self.current_project_path:
                # Quick save to existing project
                filepath = self.current_project_path
            else:
                # First-time save, prompt for location
                self.save_project_as()
                return

        # (modules already imported at top of file)

        # Create a temporary directory to build the project structure
        with tempfile.TemporaryDirectory() as temp_dir:
            sequences_dir = os.path.join(temp_dir, "sequences")
            os.makedirs(sequences_dir, exist_ok=True)

            # 1. Collect File Paths & Copy Files
            # We copy all referenced files into the 'sequences' folder
            # For duplicate imports, we save each with a unique name
            
            # Track which (file_path, import_instance) combinations we've saved
            saved_file_instances = {}
            saved_file_paths = []

            for rec in self.gbk_records:
                if hasattr(rec, 'file_path') and hasattr(rec, 'import_instance'):
                    if os.path.exists(rec.file_path):
                        import_inst = rec.import_instance
                        base_filename = os.path.basename(rec.file_path)
                        
                        # Create unique filename for this instance
                        if import_inst == 0:
                            filename = base_filename
                        else:
                            # Add instance number for duplicates: cluster.gb -> cluster_inst1.gb
                            name, ext = os.path.splitext(base_filename)
                            filename = f"{name}_inst{import_inst}{ext}"
                        
                        # Copy file with unique name
                        dest_path = os.path.join(sequences_dir, filename)
                        if not os.path.exists(dest_path):  # Only copy if not already there
                            shutil.copy2(rec.file_path, dest_path)
                        
                        rel_path = f"sequences/{filename}"
                        
                        # Track this instance
                        key = (rec.file_path, import_inst)
                        if key not in saved_file_instances:
                            saved_file_instances[key] = rel_path
                            saved_file_paths.append(rel_path)
                    else:
                        logger.warning(f"Source file not found: {rec.file_path}")

            # 2. Collect Record Metadata (UID mapping)
            # Map UIDs to (Relative File, Index, Import Instance)
            record_meta = {}
            for rec in self.gbk_records:
                if hasattr(rec, 'file_path') and hasattr(rec, 'index_in_file'):
                    import_inst = getattr(rec, 'import_instance', 0)
                    key = (rec.file_path, import_inst)
                    
                    # Get the relative path for this specific instance
                    rel_path = saved_file_instances.get(key, rec.file_path)
                    
                    record_meta[rec.uid] = {
                        'file': rel_path,
                        'idx': rec.index_in_file,
                        'import_instance': import_inst
                    }

            # 3. Build State Dictionary
            state = {
                'version': '2.0', # Bump version for bundled format
                'file_paths': saved_file_paths, # These are now relative paths like "sequences/foo.gb"
                'record_meta': record_meta,
                
                # UI State
                'pinned_uid': self.pinned_record.uid if self.pinned_record else None,
                'selected_rows': [r.uid for r in self.selected_rows],
                'favorite_rows': [r.uid for r in self.favorite_rows],
                'flipped_rows': list(self.flipped_rows),
                'expanded_rows': list(self.expanded_rows),
                'row_x_offsets': self.row_x_offsets,
                'global_x_offset': self.global_x_offset,
                'sort_mode': self.sort_mode,
                'zoom_val': self.zoom_scale.get(),
                
                # User Data
                'cluster_notes': self.cluster_notes,
                'cluster_products': self.cluster_products,
                'gene_alt_names': self.gene_alt_names,
                
                # Alignment & Colors
                'query_gene_colors': self.query_gene_colors,
                'alignment_results': self.current_alignment_results,
                'similarity_scores': self.cluster_similarity_scores,
                
                # Selections & Bench
                'selected_gene_tags': list(self.selected_gene_canvas_tags),
                'pending_alignment_tags': list(self.pending_alignment_tags),
                'gene_bench_data': self.clipboard_data,
                
                # Custom Protein Sequences (for gene bench custom proteins)
                # Extract sequences for custom proteins (those starting with "custom_")
                'custom_protein_sequences': {
                    tag: str(seq) for tag, seq in self.protein_map.items() if tag.startswith('custom_')
                },
                
                # Toggles
                'hide_genes': self.var_hide_genes.get(),
                'highlight_genes': self.var_highlight_genes.get(),
                'limit_range': self.var_limit_range.get(),
                'limit_val': self.limit_entry.get()
            }

            # Write project.json
            json_path = os.path.join(temp_dir, "project.json")
            with open(json_path, 'w') as f:
                json.dump(state, f, indent=2)

            # Create ZIP Archive
            try:
                with zipfile.ZipFile(filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    # Add project.json
                    zipf.write(json_path, arcname="project.json")
                    # Add sequences
                    for root, dirs, files in os.walk(sequences_dir):
                        for file in files:
                            abs_file = os.path.join(root, file)
                            rel_archive_path = os.path.relpath(abs_file, temp_dir)
                            zipf.write(abs_file, arcname=rel_archive_path)
                            
                self.current_project_path = filepath  # Remember this path for quick save
                self.status_label.config(text=f"Project saved (Bundled): {os.path.basename(filepath)}")
            except Exception as e:
                tk.messagebox.showerror("Save Error", f"Failed to create project archive:\n{e}")

    def load_project(self):
        """Loads a project from a .bcproject file (ZIP or JSON)."""
        filepath = filedialog.askopenfilename(
            title="Load Project",
            filetypes=[("BisCEET Project", "*.bcproject"), ("All Files", "*.*")]
        )
        if not filepath: return
        
        # (modules already imported at top of file)
        
        state = None
        temp_extract_dir = None 

        try:
            # Check if ZIP
            if zipfile.is_zipfile(filepath):
                # It's a bundled project
                temp_extract_dir = tempfile.mkdtemp()
                with zipfile.ZipFile(filepath, 'r') as zip_ref:
                    zip_ref.extractall(temp_extract_dir)
                
                json_path = os.path.join(temp_extract_dir, "project.json")
                if not os.path.exists(json_path):
                    raise Exception("Invalid project archive: missing project.json")
                
                with open(json_path, 'r') as f:
                    state = json.load(f)
                
                # Fix up file paths in state to point to extracted files
                abs_file_paths = []
                for rel_path in state.get('file_paths', []):
                    # Sanitize path just in case
                    abs_path = os.path.abspath(os.path.join(temp_extract_dir, rel_path))
                    abs_file_paths.append(abs_path)
                
                # Update state to use these absolute paths for loading
                state['file_paths'] = abs_file_paths
                
                # Also update record_meta to match these absolute paths for UID reconstruction
                for uid, meta in state.get('record_meta', {}).items():
                    rel = meta['file']
                    meta['file'] = os.path.abspath(os.path.join(temp_extract_dir, rel))

            else:
                # Legacy JSON format
                with open(filepath, 'r') as f:
                    state = json.load(f)

        except Exception as e:
            tk.messagebox.showerror("Load Error", f"Failed to parse project file:\n{e}")
            if temp_extract_dir: shutil.rmtree(temp_extract_dir)
            return

        # 1. Load Files
        file_paths = state.get('file_paths', [])
        # Verify files exist
        missing = [p for p in file_paths if not os.path.exists(p)]
        if missing:
            tk.messagebox.showerror("Missing Files", f"Cannot find the following files:\n" + "\n".join(missing[:5]))
            if temp_extract_dir: shutil.rmtree(temp_extract_dir)
            return

        self.status_label.config(text="Loading project files...")
        self.update_idletasks()
        
        # Reset state (this cleans up old temp dir)
        self._reset_app_state()
        
        # Store new temp dir so it persists
        if temp_extract_dir:
            self.project_temp_dir = temp_extract_dir

        # Load fresh
        self.gbk_records, self.cached_cds_data, self.protein_map, self.max_cluster_span = load_genbank_files(file_paths)
        if self.max_cluster_span < 1: self.max_cluster_span = 1

        # 2. Restore UIDs (Critical Step)
        # Map (File, Index, Import Instance) -> Saved UID
        saved_meta = state.get('record_meta', {})
        # Create lookup: (file_path, index, import_instance) -> uid
        lookup = {}
        max_import_id = 0  # Track the highest import_instance for next_import_id
        
        for uid, meta in saved_meta.items():
            # Ensure key is normalized and includes import_instance
            import_inst = meta.get('import_instance', 0)
            key = (os.path.abspath(meta['file']), meta['idx'], import_inst)
            lookup[key] = uid
            max_import_id = max(max_import_id, import_inst)
            
        # Iterate loaded records and update UIDs + import_instance
        new_cache = {}
        
        for rec in self.gbk_records:
            if hasattr(rec, 'file_path') and hasattr(rec, 'index_in_file'):
                # Try to find a match by checking all possible import_instances
                matched = False
                for uid, meta in saved_meta.items():
                    if (os.path.abspath(rec.file_path) == os.path.abspath(meta['file']) and
                        rec.index_in_file == meta['idx']):
                        # Found a match - restore the import_instance
                        original_uid = uid
                        rec.import_instance = meta.get('import_instance', 0)
                        
                        # Move cache data to old UID
                        if rec.uid in self.cached_cds_data:
                            data = self.cached_cds_data.pop(rec.uid)
                            
                            # Update Record UID
                            rec.uid = original_uid
                            
                            # Update CDS List & Protein Map
                            for cds in data:
                                old_tag = cds['unique_tag']
                                new_tag = f"gene_{original_uid}_{cds['start']}_{cds['end']}"
                                cds['unique_tag'] = new_tag
                                
                                if old_tag in self.protein_map:
                                    self.protein_map[new_tag] = self.protein_map.pop(old_tag)
                            
                            new_cache[original_uid] = data
                        
                        # Remove from saved_meta so we don't match it again
                        saved_meta = {k: v for k, v in saved_meta.items() if k != uid}
                        matched = True
                        break
                
                if not matched:
                    # No match found
                    new_cache[rec.uid] = self.cached_cds_data.get(rec.uid, [])
        
        self.cached_cds_data = new_cache
        
        # Update next_import_id to continue from where we left off
        self.next_import_id = max_import_id + 1

        # 3. Restore State
        
        # Pinned Record
        pinned_uid = state.get('pinned_uid')
        if pinned_uid:
            self.pinned_record = next((r for r in self.gbk_records if r.uid == pinned_uid), None)
            
        # Lists - restore from saved UIDs and populate UID sets
        sel_uids = set(state.get('selected_rows', []))
        self.selected_rows = [r for r in self.gbk_records if r.uid in sel_uids]
        self.selected_uids = sel_uids  # Populate UID set
        
        fav_uids = set(state.get('favorite_rows', []))
        self.favorite_rows = [r for r in self.gbk_records if r.uid in fav_uids]
        self.favorite_uids = fav_uids  # Populate UID set
        
        self.flipped_rows = set(state.get('flipped_rows', []))
        self.expanded_rows = set(state.get('expanded_rows', []))
        
        self.row_x_offsets = state.get('row_x_offsets', {})
        self.global_x_offset = state.get('global_x_offset', 0.0)
        self.sort_mode = state.get('sort_mode', 'default')
        
        # User Data
        self.cluster_notes = state.get('cluster_notes', {})
        self.cluster_products = state.get('cluster_products', {})
        self.gene_alt_names = state.get('gene_alt_names', {})
        
        # Alignment & Colors
        self.query_gene_colors = state.get('query_gene_colors', {})
        self.current_alignment_results = state.get('alignment_results', {})
        self.cluster_similarity_scores = state.get('similarity_scores', {})
        
        # Selections
        self.selected_gene_canvas_tags = set(state.get('selected_gene_tags', []))
        self.pending_alignment_tags = set(state.get('pending_alignment_tags', []))
        
        # Gene Bench
        self.clipboard_data = state.get('gene_bench_data', [])
        
        # Restore Custom Protein Sequences
        # Custom proteins are not in GenBank files, so their sequences must be restored from saved state
        from Bio.Seq import Seq
        custom_seqs = state.get('custom_protein_sequences', {})
        for tag, seq_str in custom_seqs.items():
            self.protein_map[tag] = Seq(seq_str)
        
        # Toggles & Controls
        self.var_hide_genes.set(state.get('hide_genes', True))
        self.var_highlight_genes.set(state.get('highlight_genes', True))
        self.var_limit_range.set(state.get('limit_range', True))
        
        limit_val = state.get('limit_val', "8")
        self.limit_entry.delete(0, tk.END)
        self.limit_entry.insert(0, limit_val)
        
        zoom = state.get('zoom_val', DEFAULT_ZOOM)
        self.zoom_scale.set(zoom)
        self.current_scale_val = zoom # Important for zoom logic

        # 4. Finalize
        self.current_project_path = filepath  # Remember this path for quick save
        self.status_label.config(text=f"Project loaded: {os.path.basename(filepath)}")
        self._update_sorted_cache()
        self._update_layout_map()
        self.update_button_states()
        self._full_redraw_sequence()
        if self.is_clipboard_open: self._draw_clipboard()

    # ==========================================
    #           EXPORT FUNCTIONS
    # ==========================================
    
    def export_marked_clusters(self):
        """Export marked clusters to a user-selected folder."""
        self._export_clusters(include_reference=False)
    
    def export_marked_and_reference(self):
        """Export marked clusters and reference cluster to a user-selected folder."""
        self._export_clusters(include_reference=True)
    
    def _export_clusters(self, include_reference: bool):
        """
        Export GenBank files to a folder.
        
        Args:
            include_reference: If True, also export the reference cluster
        """
        from tkinter import messagebox
        
        # Check if there are marked clusters
        if not self.selected_rows:
            messagebox.showinfo("Info", "No clusters are marked. Right-click clusters and select 'Mark Cluster' first.")
            return
        
        # Check reference if requested
        if include_reference and not self.pinned_record:
            messagebox.showinfo("Info", "No reference cluster is set.")
            return
        
        # Ask for destination folder
        dest_folder = filedialog.askdirectory(title="Select Destination Folder for Export")
        if not dest_folder:
            return  # User cancelled
        
        exported_count = 0
        errors = []
        
        # Collect records to export
        records_to_export = []
        
        # Add marked clusters
        for rec in self.gbk_records:
            if rec.uid in self.selected_uids:
                records_to_export.append(rec)
        
        # Add reference if requested
        if include_reference and self.pinned_record:
            # Avoid duplicates
            if self.pinned_record.uid not in self.selected_uids:
                records_to_export.append(self.pinned_record)
        
        # Export each record
        for rec in records_to_export:
            source_path = getattr(rec, 'file_path', None)
            if source_path and os.path.isfile(source_path):
                try:
                    filename = os.path.basename(source_path)
                    dest_path = os.path.join(dest_folder, filename)
                    
                    # Handle filename conflicts
                    counter = 1
                    base, ext = os.path.splitext(filename)
                    while os.path.exists(dest_path):
                        dest_path = os.path.join(dest_folder, f"{base}_{counter}{ext}")
                        counter += 1
                    
                    shutil.copy2(source_path, dest_path)
                    exported_count += 1
                except Exception as e:
                    errors.append(f"{rec.id}: {str(e)}")
                    logger.error(f"Failed to export {rec.id}: {e}")
            else:
                errors.append(f"{rec.id}: Source file not found")
        
        # Show result
        if errors:
            error_msg = f"Exported {exported_count} files.\n\nErrors:\n" + "\n".join(errors[:5])
            if len(errors) > 5:
                error_msg += f"\n... and {len(errors) - 5} more errors"
            messagebox.showwarning("Export Complete with Errors", error_msg)
        else:
            messagebox.showinfo("Export Complete", f"Successfully exported {exported_count} GenBank files to:\n{dest_folder}")
        
        self.status_label.config(text=f"Exported {exported_count} files to {os.path.basename(dest_folder)}")
