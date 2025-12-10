"""
Action mixin for GenBankBrowser.

Handles user actions including:
- Context menus for clusters and genes
- Sequence alignment orchestration  
- Clipboard (Gene Bench) operations
- External links (BLAST, GenBank)
"""
import tkinter as tk
from tkinter import messagebox
import webbrowser
import urllib.parse
import multiprocessing
import threading
import concurrent.futures
import logging
from typing import Dict, List, Tuple

from config.constants import *
from core.alignment import process_alignment_chunk
from ui.dialogs import open_note_editor, open_product_editor, prompt_custom_protein, open_gene_name_editor

logger = logging.getLogger(__name__)

class ActionMixin:
    """
    Handles user actions, context menus, alignment logic, and clipboard operations.
    """

    def on_row_right_click(self, event):
        """Shows context menu for a cluster row or specific gene."""
        if self.drag_data: return
        
        # Check if click is on the gene canvas (vs name canvas)
        clicked_gene_data = None
        rec = None
        row_idx = -1
        
        abs_y = self.org_canvas.canvasy(event.y)
        for i, layout in enumerate(self.row_layout_map):
            if layout['y'] <= abs_y < layout['y'] + layout['h']:
                row_idx = i
                break
        
        if row_idx != -1 and (0 <= row_idx < len(self.sorted_records_cache)):
            rec = self.sorted_records_cache[row_idx]
            
            # If clicked on Gene Canvas, check for gene hit
            if event.widget == self.gene_canvas:
                cx = self.gene_canvas.canvasx(event.x)
                cy = self.gene_canvas.canvasy(event.y)
                
                # Find overlapping items (more precise than find_closest for dense maps)
                items = self.gene_canvas.find_overlapping(cx-1, cy-1, cx+1, cy+1)
                if items:
                    # Check items for gene tags
                    prefix = f"cv_{row_idx}_"
                    for item_id in items:
                        tags = self.gene_canvas.gettags(item_id)
                        for t in tags:
                            if t.startswith(prefix):
                                try:
                                    start_val = int(t[len(prefix):])
                                    # Find matching CDS
                                    cds_list = self.cached_cds_data.get(rec.uid, [])
                                    for cds in cds_list:
                                        if cds['start'] == start_val:
                                            clicked_gene_data = cds
                                            break
                                except ValueError: pass
                        if clicked_gene_data: break
        
        if row_idx == -1: return

        self.row_menu.delete(0, "end")
        
        # --- 1. GENE ACTIONS (Top) ---
        if clicked_gene_data:
            pid = clicked_gene_data.get('protein_id', 'unknown')
            self.row_menu.add_command(label="Add to Gene bench", 
                                      command=lambda: self.add_to_clipboard(rec.uid, clicked_gene_data))
            self.row_menu.add_command(label="Set alternative gene name", 
                                      command=lambda: self.open_gene_name_editor(clicked_gene_data['unique_tag']))
            self.row_menu.add_command(label="Copy Sequence", command=lambda: self._copy_seq(clicked_gene_data['unique_tag'], pid))
            self.row_menu.add_command(label="Blast Protein", command=lambda: self._blast_seq(clicked_gene_data['unique_tag'], pid))
            self.row_menu.add_separator()
        
        # --- 2. CLUSTER ACTIONS ---
        def toggle_selected(r):
            self.selection.toggle_selected(r)
            self._update_sorted_cache()
            self._update_layout_map()
            self.update_button_states()
            self._full_redraw_sequence()
            
        def toggle_favorite(r):
            self.selection.toggle_favorite(r)
            self._update_sorted_cache()
            self._update_layout_map()
            self.update_button_states()
            self._full_redraw_sequence()

        is_sel = self.selection.is_selected(rec)
        
        # --- 2. VIEW / LAYOUT ACTIONS (Moved to top) ---
        is_expanded = rec.uid in self.expanded_rows
        self.row_menu.add_command(label="Close Inspect View" if is_expanded else "Inspect View", command=lambda: self.toggle_row_expansion(rec.uid))
        self.row_menu.add_command(label="Inspect View - All", command=lambda: self.bulk_expand(True))
        # Only show Close Inspect View - All if at least one cluster is expanded
        if self.expanded_rows:
            self.row_menu.add_command(label="Close Inspect View - All", command=lambda: self.bulk_expand(False))
        
        self.row_menu.add_separator()
        
        # --- 3. CLUSTER ACTIONS ---
        self.row_menu.add_command(label="Unmark Cluster" if is_sel else "Mark Cluster", 
                                  command=lambda: toggle_selected(rec))
        self.row_menu.add_command(label="Unfavorite" if self.selection.is_favorite(rec) else "Favorite", 
                                  command=lambda: toggle_favorite(rec))
        
        # RENAMED: Pin Cluster -> Set Cluster as Reference
        self.row_menu.add_command(label="Set Cluster as Reference", command=lambda: self._confirm_and_pin_cluster(rec))
        
        self.row_menu.add_separator()
        
        # --- 4. EDIT ACTIONS ---
        self.row_menu.add_command(label="Edit Note", command=lambda: self.open_note_editor(rec.uid))
        self.row_menu.add_command(label="Edit Product", command=lambda: self.open_product_editor(rec.uid))
        
        self.row_menu.add_separator()

        # --- 5. FLIP ACTION ---
        self.row_menu.add_command(label="Flip", command=lambda: self.on_cluster_flip(row_idx))
        
        # GenBank link (if accession is available)
        if hasattr(rec, 'genbank_accession') and rec.genbank_accession:
            self.row_menu.add_command(label="Open cluster in GenBank", 
                                      command=lambda: self.open_genbank_entry(rec.genbank_accession))
        
        self.row_menu.add_separator()

        # --- 5. SORT ACTIONS ---
        self.row_menu.add_command(label="Sort A-Z", command=lambda: self.set_sort_mode("alpha"))
        
        state_sim = "normal" if self.cluster_similarity_scores else "disabled"
        self.row_menu.add_command(label="Sort by Similarity", 
                                  command=lambda: self.set_sort_mode("similarity"), 
                                  state=state_sim)
        
        self.row_menu.add_separator()
        
        # --- 6. REMOVE CLUSTER ---
        self.row_menu.add_command(label="Remove Cluster from Project", 
                                  command=lambda: self.remove_cluster(rec))

        self.row_menu.tk_popup(event.x_root, event.y_root)

    def on_pinned_right_click(self, event):
        """Shows context menu for the pinned (query) cluster."""
        if not self.pinned_record: return
        menu = tk.Menu(self, tearoff=0)
        
        cx, cy = self.special_gene_canvas.canvasx(event.x), self.special_gene_canvas.canvasy(event.y)
        clicked_gene_data = None
        items = self.special_gene_canvas.find_overlapping(cx-1, cy-1, cx+1, cy+1)
        
        gene_tag = None
        if items:
            for item_id in items:
                tags = self.special_gene_canvas.gettags(item_id)
                for t in tags:
                    if t.startswith(TAG_PREFIX_REFERENCE):
                        gene_tag = t
                        try:
                            start_val = int(t[11:])
                            cds_list = self.cached_cds_data.get(self.pinned_record.uid, [])
                            for cds in cds_list:
                                if cds['start'] == start_val:
                                    clicked_gene_data = cds
                                    break
                        except: pass
                        break
                if gene_tag: break
        
        if gene_tag and clicked_gene_data:
            is_sel = gene_tag in self.selected_gene_canvas_tags
            menu.add_command(label="Deselect Gene" if is_sel else "Select Gene", 
                             command=lambda: self.on_gene_select(gene_tag))
            menu.add_command(label="Set alternative gene name", 
                             command=lambda: self.open_gene_name_editor(clicked_gene_data['unique_tag']))
            menu.add_command(label="Visually Align Clusters", 
                             command=lambda: self.align_clusters_visual(gene_tag))
            menu.add_separator()
            
        menu.add_command(label="Select All Genes", command=self.select_all_pinned_genes)
        menu.add_command(label="Deselect All Genes", command=self.deselect_all_pinned_genes)
        menu.add_separator()
        menu.add_command(label="Flip Cluster", command=lambda: self.on_cluster_flip("special"))

        if gene_tag and clicked_gene_data:
            menu.add_separator()
            pid = clicked_gene_data.get('protein_id', 'unknown')
            menu.add_command(label="Copy Sequence", command=lambda: self._copy_seq(clicked_gene_data['unique_tag'], pid))
            menu.add_command(label="Blast Protein", command=lambda: self._blast_seq(clicked_gene_data['unique_tag'], pid))

        menu.tk_popup(event.x_root, event.y_root)

    def on_pinned_name_right_click(self, event):
        """Shows context menu for the pinned (query) name panel."""
        if not self.pinned_record: return
        menu = tk.Menu(self, tearoff=0)
        
        rec = self.pinned_record
        
        # Add standard row actions
        def toggle_fav():
            self.selection.toggle_favorite(rec)
            self._full_redraw_sequence()

        menu.add_command(label="Unfavorite" if self.selection.is_favorite(rec) else "Favorite", 
                         command=toggle_fav)
        
        menu.add_command(label="Edit Note", command=lambda: self.open_note_editor(rec.uid))
        menu.add_command(label="Edit Product", command=lambda: self.open_product_editor(rec.uid))
        menu.add_separator()
        menu.add_command(label="Flip Cluster", command=lambda: self.on_cluster_flip("special"))
        
        menu.tk_popup(event.x_root, event.y_root)

    def _confirm_and_pin_cluster(self, rec):
        """Confirms with user before replacing existing reference cluster."""
        if self.pinned_record is not None:
            result = messagebox.askyesno(
                "Replace Reference Cluster",
                "Do you want to replace the current reference cluster?"
            )
            if not result:
                return
        self.pin_cluster(rec)

    def pin_cluster(self, rec):
        """Pins a cluster to the top 'Query' row."""
        self.pinned_record = rec
        # Remove from selected using dataclass helper
        self.selection.remove_selected(rec)
        self.selected_gene_canvas_tags.clear()
        # Clear alignment state
        self.alignment.clear()
        self.query_gene_colors.clear()
        self.pending_alignment_tags.clear()
        self.flipped_rows.clear()
        self.row_x_offsets.clear()
        self.global_x_offset = 0.0
        self.expanded_rows.clear()
        # Reset sort mode to default on new reference
        self.sort_mode = "default"
        self._update_sorted_cache()
        self._update_layout_map()
        self.update_button_states()
        self._full_redraw_sequence()

    def remove_cluster(self, rec):
        """Removes a cluster from the project."""
        # Remove from main list
        self.gbk_records = [r for r in self.gbk_records if r is not rec]
        
        # Remove from selection state using dataclass helpers
        self.selection.remove_selected(rec)
        self.selection.remove_favorite(rec)
        
        # Unpin if this was the pinned cluster
        if self.pinned_record is rec:
            self.pinned_record = None
            self.selected_gene_canvas_tags.clear()
            self.alignment.clear()
        
        # Remove from view state
        self.view.flipped_rows.discard(rec.uid)
        self.view.expanded_rows.discard(rec.uid)
        self.view.row_x_offsets.pop(rec.uid, None)
        
        # Remove CDS data
        self.cached_cds_data.pop(rec.uid, None)
        
        # Remove user metadata
        self.cluster_notes.pop(rec.uid, None)
        self.cluster_products.pop(rec.uid, None)
        
        # Remove protein sequences for this cluster's genes
        cds_list = self.cached_cds_data.get(rec.uid, [])
        for cds in cds_list:
            self.protein_map.pop(cds['unique_tag'], None)
        
        # Update display
        self._update_sorted_cache()
        self._update_layout_map()
        self.update_button_states()
        self._full_redraw_sequence()
        
        self.status_label.config(text=f"Removed cluster. {len(self.gbk_records)} clusters remaining.")


    def on_cluster_flip(self, row_idx_or_special):
        """Toggles the flip state of a cluster."""
        key = "special"
        if row_idx_or_special != "special":
             key = self.sorted_records_cache[row_idx_or_special].uid
        
        if key in self.flipped_rows: self.flipped_rows.remove(key)
        else: self.flipped_rows.add(key)
        
        if key in self.expanded_rows: self._full_redraw_sequence()
        else: 
            if key == "special": self._draw_special_row()
            else: self._draw_visible_rows_only()

    def on_gene_select(self, canvas_tag):
        """Toggles selection of a gene for alignment."""
        if canvas_tag in self.selected_gene_canvas_tags:
            self.selected_gene_canvas_tags.remove(canvas_tag)
            if canvas_tag in self.pending_alignment_tags:
                self.pending_alignment_tags.remove(canvas_tag)
            self.current_alignment_results = {k:v for k,v in self.current_alignment_results.items() if v[0] != canvas_tag}
        else:
            self.selected_gene_canvas_tags.add(canvas_tag)
            self.pending_alignment_tags.add(canvas_tag)
            if canvas_tag not in self.query_gene_colors:
                c = AVAILABLE_COLORS[len(self.query_gene_colors) % len(AVAILABLE_COLORS)]
                self.query_gene_colors[canvas_tag] = c
        self.update_button_states()
        self._full_redraw_sequence()

    def select_all_pinned_genes(self):
        if not self.pinned_record: return
        pinned_cds = self.cached_cds_data.get(self.pinned_record.uid, [])
        for cds in pinned_cds:
            canvas_tag = f"cv_special_{cds['start']}"
            if canvas_tag not in self.selected_gene_canvas_tags:
                self.selected_gene_canvas_tags.add(canvas_tag)
                self.pending_alignment_tags.add(canvas_tag)
                if canvas_tag not in self.query_gene_colors:
                    c = AVAILABLE_COLORS[len(self.query_gene_colors) % len(AVAILABLE_COLORS)]
                    self.query_gene_colors[canvas_tag] = c
        self.update_button_states()
        self._full_redraw_sequence()

    def deselect_all_pinned_genes(self):
        to_remove = [t for t in self.selected_gene_canvas_tags if t.startswith(TAG_PREFIX_REFERENCE)]
        for t in to_remove:
            self.selected_gene_canvas_tags.remove(t)
            if t in self.pending_alignment_tags:
                self.pending_alignment_tags.remove(t)
            # Also remove from results if it was a query
            self.current_alignment_results = {k:v for k,v in self.current_alignment_results.items() if v[0] != t}
            
        self.update_button_states()
        self._full_redraw_sequence()

    def align_clusters_visual(self, gene_tag):
        """Aligns visible clusters to center the selected gene or its homolog, auto-flipping if needed."""
        # 1. Identify target gene tag
        target_unique_tag = self.visible_canvas_tag_map.get(gene_tag)
        if not target_unique_tag: return

        # Get Pinned Cluster Data
        pinned_uid = self.pinned_record.uid
        pinned_cds_list = self.cached_cds_data.get(pinned_uid, [])
        if not pinned_cds_list: return

        # Find specific gene in pinned cluster
        pinned_cds = next((c for c in pinned_cds_list if c['unique_tag'] == target_unique_tag), None)
        if not pinned_cds: return
        
        # --- Calculate Pinned Gene's Relative Center ---
        p_min = min(c['start'] for c in pinned_cds_list)
        p_max = max(c['end'] for c in pinned_cds_list)
        p_span = p_max - p_min
        p_center_bp = (pinned_cds['start'] + pinned_cds['end']) / 2
        
        is_pinned_flipped = "special" in self.flipped_rows
        
        # Calculate center relative to visual start (accounts for flip)
        if is_pinned_flipped:
            p_rel_center = p_span + p_min - p_center_bp
        else:
            p_rel_center = p_center_bp - p_min
            
        p_offset = self.row_x_offsets.get("special", 0)

        # 2. Iterate through visible rows to align
        for rec in self.sorted_records_cache:
            cds_list = self.cached_cds_data.get(rec.uid, [])
            
            # Find homolog mapping to the pinned gene (query)
            homolog_cds = None
            for cds in cds_list:
                hit = self.current_alignment_results.get(cds['unique_tag'])
                if hit and hit[0] == gene_tag: 
                    homolog_cds = cds
                    break
            
            if homolog_cds:
                # --- Calculate Homolog Center & Auto-Flip ---
                h_min = min(c['start'] for c in cds_list)
                h_max = max(c['end'] for c in cds_list)
                h_span = h_max - h_min
                h_center_bp = (homolog_cds['start'] + homolog_cds['end']) / 2
                
                # Auto-Flip: Match visual orientation to pinned gene
                p_eff_strand = pinned_cds['strand'] * (-1 if is_pinned_flipped else 1)
                should_flip = (p_eff_strand != homolog_cds['strand'])
                
                if should_flip:
                    self.flipped_rows.add(rec.uid)
                elif rec.uid in self.flipped_rows:
                    self.flipped_rows.remove(rec.uid)
                
                # Calculate homolog relative center with new flip state
                if should_flip:
                    h_rel_center = h_span + h_min - h_center_bp
                else:
                    h_rel_center = h_center_bp - h_min
                
                # Set offset to align visual centers
                self.row_x_offsets[rec.uid] = p_offset + p_rel_center - h_rel_center

        self._full_redraw_sequence()

    def update_button_states(self) -> None:
        """Updates the enabled/disabled state of mapping buttons based on selections."""
        # Map Reference Genes button enabled if there are selected reference tags
        reference_tags = [t for t in self.selected_gene_canvas_tags if t.startswith(TAG_PREFIX_REFERENCE)]
        self.align_button.config(state="normal" if reference_tags else "disabled")
        self.core_align_button.config(state="normal" if self.pinned_record and self.selected_rows else "disabled")
        # Gene Bench button enabled if there are selected clip_ tags
        selected_clip_tags = [t for t in self.selected_gene_canvas_tags if t.startswith(TAG_PREFIX_CLIPBOARD)]
        self.bench_align_button.config(state="normal" if selected_clip_tags else "disabled")

    def run_alignment(self) -> None:
        """Starts the pairwise alignment process for REFERENCE cluster genes only."""
        self._run_gene_mapping(
            tag_prefix=TAG_PREFIX_REFERENCE,
            cutoff_entry=self.cutoff_entry,
            button=self.align_button,
            status_prefix="Aligning",
            empty_message="No reference genes selected. Click on genes in the Reference cluster to select them.",
            use_canvas_map=True  # Reference genes use visible_canvas_tag_map
        )

    def run_gene_bench_alignment(self) -> None:
        """Starts the pairwise alignment process for SELECTED Gene Bench genes."""
        self._run_gene_mapping(
            tag_prefix=TAG_PREFIX_CLIPBOARD,
            cutoff_entry=self.bench_cutoff_entry,
            button=self.bench_align_button,
            status_prefix="Mapping Gene Bench",
            empty_message="No Gene Bench genes selected. Use right-click > Select All Genes or click on genes to select.",
            use_canvas_map=False  # Gene Bench genes extract unique_tag from tag directly
        )

    def _run_gene_mapping(
        self,
        tag_prefix: str,
        cutoff_entry: tk.Entry,
        button: tk.Button,
        status_prefix: str,
        empty_message: str,
        use_canvas_map: bool
    ) -> None:
        """
        Consolidated mapping function for both reference and gene bench genes.
        
        Args:
            tag_prefix: Canvas tag prefix to filter selected genes (e.g., 'cv_special_' or 'clip_')
            cutoff_entry: Entry widget containing the cutoff percentage
            button: Button to disable during mapping
            status_prefix: Status message prefix
            empty_message: Message to show if no matching genes selected
            use_canvas_map: If True, use visible_canvas_tag_map to get unique_tag; else extract from tag
        """
        # Get matching tags
        matching_tags = [t for t in self.selected_gene_canvas_tags if t.startswith(tag_prefix)]
        
        if not matching_tags:
            messagebox.showinfo("Info", empty_message)
            return
            
        # Parse cutoff
        try: 
            cutoff = float(cutoff_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid Cutoff value. Please enter a number.")
            return
        
        self.status_label.config(text=f"{status_prefix}... (0%)")
        self.update_idletasks()
        button.config(state="disabled")
        
        # Build query data
        query_data: Dict[str, str] = {}
        prefix_len = len(tag_prefix)
        
        for tag in matching_tags:
            # Get unique_tag based on source type
            if use_canvas_map:
                ut = self.visible_canvas_tag_map.get(tag)
            else:
                ut = tag[prefix_len:]  # Extract unique_tag from tag
                
            if ut and ut in self.protein_map:
                query_data[tag] = self.protein_map[ut]
                # Assign color if not already assigned
                if tag not in self.query_gene_colors:
                    self.query_gene_colors[tag] = AVAILABLE_COLORS[len(self.query_gene_colors) % len(AVAILABLE_COLORS)]
                # Clear pending tag since we're mapping now
                self.pending_alignment_tags.discard(tag)
        
        if not query_data:
            self.status_label.config(text=f"No valid sequences in selected genes.")
            button.config(state="normal")
            return
        
        # Build target list from all clusters (except pinned)
        target_list: List[Tuple[str, str]] = []
        records = [r for r in self.gbk_records if r is not self.pinned_record]
        for rec in records:
            cds_list = self.cached_cds_data.get(rec.uid, [])
            for cds in cds_list:
                ut = cds['unique_tag']
                seq = self.protein_map.get(ut)
                if seq:
                    target_list.append((ut, seq))

        if not target_list:
            self.status_label.config(text="No clusters to align against.")
            button.config(state="normal")
            return

        self._launch_multiprocess(target_list, query_data, cutoff, self._finish_align, status_prefix)

    def run_core_gene_alignment(self):
        """Starts the 'Find Homologs' process to identify core genes."""
        if not self.pinned_record or not self.selected_rows: return
            
        try: 
            cutoff = float(self.core_cutoff_entry.get())
        except ValueError: 
            messagebox.showerror("Error", "Invalid Cutoff value. Please enter a number.")
            return

        self.status_label.config(text="Finding Core... (0%)")
        self.update_idletasks()
        self.core_align_button.config(state="disabled")

        query_data = {}
        pinned_cds = self.cached_cds_data.get(self.pinned_record.uid, [])
        
        # Always clear selection and search all genes in pinned cluster
        self.selected_gene_canvas_tags.clear()
        self.query_gene_colors.clear()
        self.pending_alignment_tags.clear()

        for i, cds in enumerate(pinned_cds):
            canvas_tag = f"cv_special_{cds['start']}"
            seq = self.protein_map.get(cds['unique_tag'])
            if seq: query_data[canvas_tag] = seq

        target_list = []
        for rec in self.selected_rows:
            cds_list = self.cached_cds_data.get(rec.uid, [])
            for cds in cds_list:
                ut = cds['unique_tag']
                seq = self.protein_map.get(ut)
                if seq: target_list.append((ut, seq))

        self._launch_multiprocess(target_list, query_data, cutoff, self._finish_core, "Finding Core")

    def _launch_multiprocess(self, target_list, query_data, cutoff, callback, status_prefix="Processing"):
        if not target_list: return
        try: num_workers = max(1, multiprocessing.cpu_count())
        except NotImplementedError: num_workers = 2
        
        # Create more chunks than workers to ensure smoother progress updates
        # Target at least 4 chunks per worker
        chunk_size = max(1, len(target_list) // (num_workers * 4))
        chunks = [target_list[i:i + chunk_size] for i in range(0, len(target_list), chunk_size)]
        
        threading.Thread(target=self._multiprocess_handler, args=(chunks, query_data, cutoff, callback, status_prefix), daemon=True).start()

    def _multiprocess_handler(self, chunks, query_data, cutoff, callback, status_prefix):
        combined_results = {}
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Use imported process_alignment_chunk
            futures = [executor.submit(process_alignment_chunk, chunk, query_data, cutoff) for chunk in chunks]
            for future in concurrent.futures.as_completed(futures):
                try:
                    res = future.result()
                    combined_results.update(res)
                except Exception as e:
                    logger.error(f"Worker failed: {e}")
                
                # Update Progress
                completed_count = sum(1 for f in futures if f.done())
                pct = int((completed_count / len(chunks)) * 100)
                self.master.after(0, lambda p=pct: self.status_label.config(text=f"{status_prefix}... ({p}%)"))
        self.master.after(0, callback, combined_results)

    def _finish_align(self, res):
        """Finishes alignment and merges results using higher-similarity-wins logic."""
        # Filter: Keep only the best hit per query gene for each cluster
        filtered_res = self._filter_best_hits_per_cluster(res)
        
        # Merge with existing results using higher-similarity-wins logic
        self._merge_alignment_results(filtered_res)
        
        # Update cluster similarity scores based on current alignment results
        self.cluster_similarity_scores.clear()
        for rec in self.gbk_records:
            cds_list = self.cached_cds_data.get(rec.uid, [])
            scores = [self.current_alignment_results[c['unique_tag']][1] for c in cds_list if c['unique_tag'] in self.current_alignment_results]
            self.cluster_similarity_scores[rec.uid] = sum(scores) / len(scores) if scores else 0.0

        self.status_label.config(text=f"Done. Found {len(filtered_res)} best hits ({len(self.current_alignment_results)} total mapped).")
        self.update_button_states()
        self._full_redraw_sequence()

    def _merge_alignment_results(self, new_results: Dict[str, Tuple[str, float]]):
        """Merge new alignment results, keeping only higher similarity mappings.
        
        First removes all existing mappings from the query genes being remapped,
        then adds new results using higher-similarity-wins logic.
        """
        # Collect all query tags in this batch
        query_tags_in_batch = set(query_tag for query_tag, _ in new_results.values())
        
        # Remove all existing mappings from these query genes
        # This ensures that if cutoff is raised, old low-similarity hits are cleared
        targets_to_remove = [
            target_tag for target_tag, (query_tag, _) in self.current_alignment_results.items()
            if query_tag in query_tags_in_batch
        ]
        for target_tag in targets_to_remove:
            del self.current_alignment_results[target_tag]
        
        # Now add new results using higher-similarity-wins logic
        for target_tag, (query_tag, similarity) in new_results.items():
            if target_tag in self.current_alignment_results:
                existing_query, existing_sim = self.current_alignment_results[target_tag]
                if similarity > existing_sim:
                    # New mapping wins - replace
                    self.current_alignment_results[target_tag] = (query_tag, similarity)
            else:
                # No existing mapping - add new
                self.current_alignment_results[target_tag] = (query_tag, similarity)

    def _filter_best_hits_per_cluster(self, raw_results: Dict[str, List[Tuple[str, float]]]) -> Dict[str, Tuple[str, float]]:
        """Ensures that for a given cluster, each Query Gene maps to at most ONE Target Gene (the best one)."""
        final_results = {}

        # Iterate over each cluster
        for rec in self.gbk_records:
            cds_list = self.cached_cds_data.get(rec.uid, [])
            
            # Group hits in this cluster by Query Tag
            # Key: Query Tag, Value: List of (Target Tag, Score)
            hits_by_query = {} 
            
            for cds in cds_list:
                t_tag = cds['unique_tag']
                if t_tag in raw_results:
                    # Iterate through ALL hits for this target
                    for q_tag, score in raw_results[t_tag]:
                        if q_tag not in hits_by_query:
                            hits_by_query[q_tag] = []
                        hits_by_query[q_tag].append((t_tag, score))
            
            # For each query, pick the best target in this cluster
            for q_tag, candidates in hits_by_query.items():
                # Sort by score descending
                candidates.sort(key=lambda x: x[1], reverse=True)
                best_t_tag, best_score = candidates[0]
                
                # Add to final results
                final_results[best_t_tag] = (q_tag, best_score)
                
        return final_results

    def _finish_core(self, results):
        tag_to_rec_id = {}
        for rec in self.selected_rows:
            for cds in self.cached_cds_data.get(rec.uid, []):
                tag_to_rec_id[cds['unique_tag']] = rec.uid
        
        query_hits_per_record = {} 
        for t_ut, hits in results.items():
            rec_id = tag_to_rec_id.get(t_ut)
            if rec_id:
                for q_ct, _ in hits:
                    if q_ct not in query_hits_per_record: query_hits_per_record[q_ct] = set()
                    query_hits_per_record[q_ct].add(rec_id)
        
        core_genes = set()
        num_required = len(self.selected_rows)
        for q_ct, hit_recs in query_hits_per_record.items():
            if len(hit_recs) == num_required: core_genes.add(q_ct)
                
        self.selected_gene_canvas_tags = core_genes
        # CRITICAL FIX: ensure these are marked as pending so '*' appears
        self.pending_alignment_tags = set(core_genes)
        
        self.query_gene_colors.clear()
        for i, q_ct in enumerate(sorted(list(core_genes))):
             self.query_gene_colors[q_ct] = AVAILABLE_COLORS[i % len(AVAILABLE_COLORS)]
             
        self.current_alignment_results = {} 
        self.selected_rows.clear()
        self.selected_uids.clear()
        self.status_label.config(text=f"Core analysis done. Found {len(core_genes)} core genes.")
        self.update_button_states()
        self._full_redraw_sequence()

    def _copy_seq(self, unique_tag, protein_id):
        if unique_tag in self.protein_map:
            seq_str = str(self.protein_map[unique_tag])
            # Format: >identifier [line break] [sequence]
            text_to_copy = f">{protein_id}\n{seq_str}"
            self.clipboard_clear()
            self.clipboard_append(text_to_copy)
            messagebox.showinfo("Copied", f"Sequence for {protein_id} copied to clipboard.")

    def _blast_seq(self, unique_tag, protein_id):
        if unique_tag in self.protein_map:
            seq_str = str(self.protein_map[unique_tag])
            # Format query: >protein_id\nSEQUENCE
            query = f">{protein_id}\n{seq_str}"
            encoded_query = urllib.parse.quote(query)
            url = f"{NCBI_BLAST_URL}?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY={encoded_query}"
            webbrowser.open(url)

    def open_genbank_entry(self, accession):
        """Opens the NCBI GenBank page for the given accession number."""
        if accession:
            url = f"{NCBI_GENBANK_URL}/{accession}"
            webbrowser.open(url)

    # --- DIALOG WRAPPERS ---
    def open_note_editor(self, uid):
        open_note_editor(self, uid)

    def open_product_editor(self, uid):
        open_product_editor(self, uid)

    def open_gene_name_editor(self, unique_tag):
        open_gene_name_editor(self, unique_tag)

    # --- CLIPBOARD METHODS ---
    
    def toggle_clipboard(self):
        self.is_clipboard_open = not self.is_clipboard_open
        if self.is_clipboard_open:
            self.clipboard_content.pack(fill=tk.X)
            self.btn_clipboard_toggle.config(text="-")
            self._draw_clipboard()
        else:
            self.clipboard_content.pack_forget()
            self.btn_clipboard_toggle.config(text="+")

    def add_to_clipboard(self, rec_uid, cds_data):
        # Check duplicates based on unique_tag
        for item in self.clipboard_data:
            if item['cds']['unique_tag'] == cds_data['unique_tag']:
                return # Already present
                
        self.clipboard_data.append({'rec_uid': rec_uid, 'cds': cds_data})
        if self.is_clipboard_open:
            self._draw_clipboard()
            
    def on_clipboard_click(self, event):
        pass # Selection disabled on left click

    def on_clipboard_right_click(self, event):
        cx = self.clipboard_canvas.canvasx(event.x)
        cy = self.clipboard_canvas.canvasy(event.y)
        
        # Check if clicked on background (to add custom protein)
        menu = tk.Menu(self, tearoff=0)
        
        # Check if an item was clicked
        item = self.clipboard_canvas.find_closest(cx, cy, halo=2)
        gene_tag = None
        unique_tag = None
        
        if item:
            tags = self.clipboard_canvas.gettags(item[0])
            for t in tags:
                if t.startswith(TAG_PREFIX_CLIPBOARD):
                    gene_tag = t
                    unique_tag = t[5:] # Remove clip_
                    break
        
        if unique_tag:
            # Find CDS data
            cds_data = None
            for cd in self.clipboard_data:
                if cd['cds']['unique_tag'] == unique_tag:
                    cds_data = cd['cds']
                    break
            
            is_sel = gene_tag in self.selected_gene_canvas_tags
            menu.add_command(label="Deselect Gene" if is_sel else "Select Gene", 
                             command=lambda: self.on_gene_select(gene_tag))
            menu.add_separator()

            menu.add_command(label="Remove from Gene bench", command=lambda: self._remove_from_clipboard(unique_tag))
            menu.add_separator()
            
            # Copy and BLAST work for both regular and custom genes (sequence is in protein_map)
            if unique_tag in self.protein_map:
                pid = cds_data.get('protein_id', 'unknown') if cds_data else 'unknown'
                product_name = cds_data.get('product', 'Custom Protein') if cds_data else 'Custom Protein'
                menu.add_command(label="Copy Sequence", command=lambda: self._copy_seq(unique_tag, pid))
                menu.add_command(label="Blast Protein", command=lambda: self._blast_seq(unique_tag, pid))
            
            menu.add_separator()
        
        # Always show option to add custom protein
        menu.add_command(label="Add Custom Protein", command=self._prompt_custom_protein)
        menu.add_separator()
        menu.add_command(label="Select All Genes", command=self.select_all_clipboard_genes)
        menu.add_command(label="Deselect All Genes", command=self.deselect_all_clipboard_genes)
        
        menu.tk_popup(event.x_root, event.y_root)

    def select_all_clipboard_genes(self):
        """Selects all genes currently in the clipboard (Gene Bench)."""
        for item in self.clipboard_data:
            cds = item['cds']
            canvas_tag = f"clip_{cds['unique_tag']}"
            if canvas_tag not in self.selected_gene_canvas_tags:
                self.selected_gene_canvas_tags.add(canvas_tag)
                self.pending_alignment_tags.add(canvas_tag)
                if canvas_tag not in self.query_gene_colors:
                    c = AVAILABLE_COLORS[len(self.query_gene_colors) % len(AVAILABLE_COLORS)]
                    self.query_gene_colors[canvas_tag] = c
        self.update_button_states()
        self._draw_clipboard()

    def deselect_all_clipboard_genes(self):
        """Deselects all genes currently in the clipboard (Gene Bench)."""
        to_remove = [t for t in self.selected_gene_canvas_tags if t.startswith(TAG_PREFIX_CLIPBOARD)]
        for t in to_remove:
            self.selected_gene_canvas_tags.remove(t)
            if t in self.pending_alignment_tags:
                self.pending_alignment_tags.remove(t)
            # Also remove from results if it was a query
            self.current_alignment_results = {k:v for k,v in self.current_alignment_results.items() if v[0] != t}
            
        self.update_button_states()
        self._draw_clipboard()
        self._full_redraw_sequence()  # Redraw to update cluster colors

    def _remove_from_clipboard(self, unique_tag):
        self.clipboard_data = [x for x in self.clipboard_data if x['cds']['unique_tag'] != unique_tag]
        # Also deselect if selected
        clip_tag = f"clip_{unique_tag}"
        if clip_tag in self.selected_gene_canvas_tags:
            self.selected_gene_canvas_tags.remove(clip_tag)
        if clip_tag in self.pending_alignment_tags:
            self.pending_alignment_tags.remove(clip_tag)
        
        # Remove alignment results that were created by this gene
        targets_to_remove = [
            target_tag for target_tag, (query_tag, _) in self.current_alignment_results.items()
            if query_tag == clip_tag
        ]
        for target_tag in targets_to_remove:
            del self.current_alignment_results[target_tag]
        
        # Update cluster similarity scores
        if targets_to_remove:
            self.cluster_similarity_scores.clear()
            for rec in self.gbk_records:
                cds_list = self.cached_cds_data.get(rec.uid, [])
                scores = [self.current_alignment_results[c['unique_tag']][1] for c in cds_list if c['unique_tag'] in self.current_alignment_results]
                self.cluster_similarity_scores[rec.uid] = sum(scores) / len(scores) if scores else 0.0
        
        self.update_button_states()
        self._draw_clipboard()
        self._full_redraw_sequence()  # Redraw to update cluster colors

    # --- Custom Protein Logic ---
    def _prompt_custom_protein(self):
        prompt_custom_protein(self)
