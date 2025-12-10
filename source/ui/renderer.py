"""
Cluster renderer for BisCEET.

Handles all canvas drawing operations including:
- Gene arrows with directional indicators
- Organism name panels
- Stacked gene labels
- Core gene highlighting
- Clipboard canvas rendering
"""
import tkinter as tk
from tkinter import ttk
from config.constants import *
from utils.utils import draw_arrow_poly, draw_stacked_labels, draw_core_bracket

class ClusterRenderer:
    """
    Handles all canvas drawing operations for the main window.
    Separates rendering logic from the main application class.
    """
    def __init__(self, app):
        self.app = app
        self._group_color_cache = {}  # Maps group_id -> color
    
    def _get_group_color(self, query_tag: str) -> str:
        """
        Get color based on the query tag's source group.
        
        Groups:
        - 'reference' = all reference cluster genes (cv_special_*)
        - 'organism_{rec_uid}' = genes from that organism in Gene Bench  
        - 'custom' = custom protein genes
        
        Returns a consistent color for the group.
        """
        # Determine group from query tag
        if query_tag.startswith(TAG_PREFIX_REFERENCE):
            group_id = "reference"
        elif query_tag.startswith(TAG_PREFIX_CLIPBOARD):
            # Extract unique_tag from clip_tag and find the organism
            ut = query_tag[len(TAG_PREFIX_CLIPBOARD):]
            # Look up in clipboard_data to find the organism
            group_id = "custom"  # Default if not found
            for item in self.app.clipboard_data:
                if item['cds']['unique_tag'] == ut:
                    group_id = f"organism_{item['rec_uid']}"
                    break
        else:
            group_id = "unknown"
        
        # Assign color if not cached
        if group_id not in self._group_color_cache:
            color_idx = len(self._group_color_cache) % len(AVAILABLE_COLORS)
            self._group_color_cache[group_id] = AVAILABLE_COLORS[color_idx]
        
        return self._group_color_cache[group_id]
    
    def clear_group_colors(self):
        """Clear the group color cache (call when reference changes)."""
        self._group_color_cache.clear()

    def draw_visible_rows_only(self):
        """
        Draws only the rows currently visible in the viewport to optimize performance.
        Clears canvases and redraws name panels and gene tracks.
        """
        self.app.org_canvas.delete("all")
        self.app.gene_canvas.delete("all")
        for b in self.app.visible_buttons: b.destroy()
        self.app.visible_buttons.clear()

        if not self.app.sorted_records_cache: return

        # Viewport Culling
        top_y = self.app.gene_canvas.canvasy(0)
        bot_y = self.app.gene_canvas.canvasy(self.app.gene_canvas.winfo_height())

        start_idx = 0
        for i, layout in enumerate(self.app.row_layout_map):
            if layout['y'] + layout['h'] > top_y:
                start_idx = i
                break
        
        scale = self.app.zoom_scale.get() / 1000.0
        
        sr_w = self.app.gene_canvas.cget("scrollregion").split()[2] if self.app.gene_canvas.cget("scrollregion") else "50000"
        bg_width = float(sr_w) + 5000
        
        for i in range(start_idx, len(self.app.sorted_records_cache)):
            layout = self.app.row_layout_map[i]
            if layout['y'] > bot_y: break 

            rec = self.app.sorted_records_cache[i]
            y_top = layout['y']
            y_bot = y_top + layout['h']
            y_mid = y_bot - (DEFAULT_ROW_HEIGHT / 2)

            # Draw Name Panel
            bg = COLOR_BG_SELECTED if rec.uid in self.app.selected_uids else COLOR_BG_WHITE
            self.draw_name_panel(self.app.org_canvas, rec, y_top, y_bot, y_mid, bg)

            # Draw Genes
            self.app.gene_canvas.create_rectangle(0, y_top, bg_width, y_bot, fill=bg, outline="#e0e0e0")
            self.draw_cluster_features(self.app.gene_canvas, rec, i, y_mid, scale)

    def _get_render_indices(self, cds_list, limit_count):
        """
        Calculates which genes to render when 'Truncate' is active.
        Shows genes around the aligned core, hiding distant neighbors.
        """
        if limit_count == float('inf'):
            return range(len(cds_list))
            
        aligned_indices = [idx for idx, c in enumerate(cds_list) if self.app.current_alignment_results.get(c['unique_tag'])]
        
        if not aligned_indices: 
            return range(len(cds_list))

        min_core = min(aligned_indices)
        max_core = max(aligned_indices)
        
        valid = []
        for idx in range(len(cds_list)):
            if min_core - limit_count <= idx <= max_core + limit_count:
                valid.append(idx)
        return valid

    def draw_name_panel(self, canvas, rec, y_top, y_bot, y_mid, bg, is_reference=False):
        """Draws the left-side panel with organism name, favorites star, and product info."""
        # Draw background
        canvas.create_rectangle(0, y_top, ROW_NAME_WIDTH, y_bot, fill=bg, outline="#e0e0e0")
        
        current_y = y_mid
        
        # Reference cluster text is now handled by the title bar, so we don't draw it here.

        txt = rec.annotations.get('organism', rec.id)
        
        # Create a unique tag for the organism name text
        name_tag = f"org_name_{rec.uid}"
        
        if rec.uid in self.app.favorite_uids:
            canvas.create_text(10, current_y, text="★", fill=COLOR_TEXT_FAVORITE, font=("tahoma", 10, "bold"))
            canvas.create_text(25, current_y, text=txt, anchor="w", width=200, tags=name_tag)
        else:
            canvas.create_text(5, current_y, text=txt, anchor="w", width=230, tags=name_tag)
        
        # Bind tooltip to show DEFINITION if available
        if hasattr(rec, 'definition') and rec.definition:
            def show_def_tooltip(event, definition=rec.definition):
                self.app.button_tooltip.showtip(definition, event.x_root, event.y_root)
            
            def hide_def_tooltip(event):
                self.app.button_tooltip.hidetip()
            
            canvas.tag_bind(name_tag, '<Enter>', show_def_tooltip)
            canvas.tag_bind(name_tag, '<Leave>', hide_def_tooltip)
            
        # Product Name (Blue)
        prod_name = self.app.cluster_products.get(rec.uid, "")
        if prod_name:
             canvas.create_text(25 if rec.uid in self.app.favorite_uids else 5, 
                                         current_y + 14, text=prod_name, fill=COLOR_PRODUCT_NAME, 
                                         font=("tahoma", 7), anchor="w")

        # --- Round Note Button Logic ---
        has_note = bool(self.app.cluster_notes.get(rec.uid, "").strip())
        if has_note:
            # Draw Circle
            r = 8  # Radius
            x_center = ROW_NAME_WIDTH - 15
            y_center = current_y # Vertically centered
            
            # Tags for binding
            btn_tag = f"note_btn_{rec.uid}"
            
            # Circle background
            canvas.create_oval(x_center - r, y_center - r, x_center + r, y_center + r, 
                                        fill="#ffeb3b", outline="#d4ac0d", tags=btn_tag)
            # Icon (simple 'N' or could use unicode)
            canvas.create_text(x_center, y_center, text="N", font=("tahoma", 7, "bold"), tags=btn_tag)
            
            # Bindings for the canvas items
            note_txt = self.app.cluster_notes[rec.uid]
            
            def on_click(e, uid=rec.uid):
                self.app.open_note_editor(uid)
            
            def on_enter(e, t=note_txt):
                self.app.button_tooltip.showtip(t, e.x_root, e.y_root)
                
            def on_leave(e):
                self.app.button_tooltip.hidetip()
                
            canvas.tag_bind(btn_tag, "<Button-1>", on_click)
            canvas.tag_bind(btn_tag, "<Enter>", on_enter)
            canvas.tag_bind(btn_tag, "<Leave>", on_leave)

    def draw_cluster_features(self, canvas, record, row_idx_or_special, y_mid, scale):
        """
        Draws the genes (arrows) for a single cluster.
        Handles coordinate scaling, flipping, coloring, and labeling.
        """
        is_special = (row_idx_or_special == "special")
        # Use UUID
        cds_list = self.app.cached_cds_data.get(record.uid, [])
        if not cds_list: return

        min_start = min(c['start'] for c in cds_list)
        span_bp = max(c['end'] for c in cds_list) - min_start
        
        # Use UUID
        state_key = "special" if is_special else record.uid
        is_flipped = state_key in self.app.flipped_rows
        is_expanded = (state_key in self.app.expanded_rows) and not is_special
        
        drag_off_bp = self.app.row_x_offsets.get(state_key, 0)
        x_base = 10 + ((drag_off_bp + self.app.global_x_offset) * scale)
        
        canvas.create_line(x_base, y_mid, x_base + (span_bp * scale), y_mid, dash=(2, 2), fill="grey")
        
        limit_count = float('inf')
        if self.app.var_limit_range.get():
            try: limit_count = int(self.app.limit_entry.get())
            except ValueError: pass
            
        indices = self._get_render_indices(cds_list, limit_count) if is_expanded else range(len(cds_list))

        draw_items = []
        for idx in indices:
            cds = cds_list[idx]
            start, end = cds['start'], cds['end']
            
            if is_flipped:
                rel_s = span_bp - (end - min_start)
                strand = -cds['strand']
            else:
                rel_s = start - min_start
                strand = cds['strand']
            
            x1 = x_base + (rel_s * scale)
            x2 = x_base + ((rel_s + (end - start)) * scale)
            draw_items.append({'x1': x1, 'x2': x2, 'strand': strand, 'data': cds})

        draw_items.sort(key=lambda item: item['x1'])

        label_tasks = []
        core_x_coords = []
        text_above = True 
        alt_name_count = 0

        for item in draw_items:
            x1, x2, strand, cds = item['x1'], item['x2'], item['strand'], item['data']
            canvas_tag = f"cv_{row_idx_or_special}_{cds['start']}"
            unique_tag = cds['unique_tag']
            
            color = COLOR_DEFAULT_GENE
            hit = self.app.current_alignment_results.get(unique_tag)
            
            if is_special:
                 self.app.visible_canvas_tag_map[canvas_tag] = unique_tag
                 if canvas_tag in self.app.selected_gene_canvas_tags:
                    # Use group color or individual gene color based on toggle
                    if self.app.var_color_by_gene.get():
                        color = self.app.query_gene_colors.get(canvas_tag, COLOR_QUERY_DEFAULT)
                    else:
                        color = self._get_group_color(canvas_tag)
            else:
                if hit:
                    core_x_coords.extend([x1, x2])
                    # In expanded (Inspect) mode, use white if Monochrome is on, else use color
                    if is_expanded:
                        if self.app.var_inspect_monochrome.get():
                            color = COLOR_CORE_HIGHLIGHT  # White
                        else:
                            # Group or gene color
                            if self.app.var_color_by_gene.get():
                                color = self.app.query_gene_colors.get(hit[0], COLOR_QUERY_DEFAULT)
                            else:
                                color = self._get_group_color(hit[0])
                    else:
                        # Group or gene color
                        if self.app.var_color_by_gene.get():
                            color = self.app.query_gene_colors.get(hit[0], COLOR_QUERY_DEFAULT)
                        else:
                            color = self._get_group_color(hit[0])
                elif is_expanded and self.app.var_highlight_genes.get() and any(k in cds['product'].lower() for k in self.app.highlight_genes_list):
                    color = COLOR_HIGHLIGHT_MATCH

            # Use imported utility function
            poly = draw_arrow_poly(canvas, x1, x2, y_mid, strand, color, canvas_tag)
            
            tip = f"Product: {cds['product']}\nSize: {cds['end']-cds['start']} bp"
            alt_name = self.app.gene_alt_names.get(unique_tag)
            if alt_name:
                tip += f"\nAlt. name: {alt_name}"
                
                # Draw Alt Name over gene if Reference Row
                if is_special:
                     # Stagger logic: 3 steps (Base, +8, +16)
                     # Base is 20px up.
                     step = alt_name_count % 3
                     y_offset = 20 + (step * 8)
                     
                     text_y = y_mid - y_offset
                     center_x = (x1+x2)/2
                     
                     # Draw connecting line (from text bottom to 3px above gene top)
                     # Assuming text height ~12px, bottom is text_y + 6. 
                     # Gene top is approx y_mid - 6. Gap of 3px -> y_mid - 9.
                     canvas.create_line(center_x, text_y + 6, center_x, y_mid - 9, fill="black", width=1)
                     
                     canvas.create_text(center_x, text_y, text=alt_name, font=("tahoma", 8, "bold"), fill="black")
                     alt_name_count += 1
            
            # Bind hover events for tooltip AND border highlighting
            canvas.tag_bind(poly, "<Enter>", lambda e, t=tip, c=canvas, i=poly: self.app._on_gene_enter(e, t, c, i))
            canvas.tag_bind(poly, "<Leave>", lambda e, c=canvas, i=poly: self.app._on_gene_leave(e, c, i))

            # Draw '*' for pending selections
            if canvas_tag in self.app.pending_alignment_tags:
                 canvas.create_text((x1+x2)/2, y_mid + 19, text="*", font=("tahoma", 12, "bold"))

            # In Inspect mode, show labels for non-mapped genes only (mapped genes have no labels)
            if is_expanded and not hit:
                should_show = True
                if self.app.var_hide_genes.get() and any(k in cds['product'].lower() for k in self.app.hide_genes_list):
                    should_show = False

                if should_show:
                    label_tasks.append((x1, x2, y_mid, cds['product'], color == COLOR_HIGHLIGHT_MATCH))

            # Show identity % if toggle is on
            if not is_special and not is_expanded and hit and hit[1] >= 20 and self.app.var_show_identity.get():
                lbl_y = y_mid - 18 if text_above else y_mid + 22
                canvas.create_text((x1+x2)/2, lbl_y, text=f"{int(hit[1])}%", font=("tahoma", 7))
                text_above = not text_above

        if label_tasks:
            # Use imported utility function
            draw_stacked_labels(canvas, label_tasks, self.app.label_font, self.app.label_font_bold)

        if is_expanded and core_x_coords and not is_special:
             # Use imported utility function
             draw_core_bracket(canvas, min(core_x_coords), max(core_x_coords), y_mid)

    def draw_special_row(self):
        """Redraws the pinned 'Query' cluster at the top of the window."""
        self.app.special_gene_canvas.delete("all")
        self.app.special_org_frame.delete("all") # Now a Canvas
        self.app.visible_canvas_tag_map.clear()

        if not self.app.pinned_record: 
            # Title bar handles the label now
            self.app.special_gene_canvas.configure(scrollregion=(0,0,0,0))
            return

        rec = self.app.pinned_record
        
        # Draw Name Panel using generalized function
        # y_top=0, y_bot=DEFAULT_ROW_HEIGHT, y_mid=DEFAULT_ROW_HEIGHT/2
        self.draw_name_panel(self.app.special_org_frame, rec, 0, DEFAULT_ROW_HEIGHT, DEFAULT_ROW_HEIGHT/2, COLOR_BG_WHITE, is_reference=True)
        
        scale = self.app.zoom_scale.get() / 1000.0
        self.draw_cluster_features(self.app.special_gene_canvas, rec, "special", DEFAULT_ROW_HEIGHT/2, scale)

    def draw_clipboard(self):
        """Renders the clipboard contents (saved genes/proteins) in the bottom panel."""
        self.app.clipboard_canvas.delete("all")
        if not self.app.clipboard_data: return
        
        # Group by record
        grouped = {}
        for item in self.app.clipboard_data:
            rid = item['rec_uid']
            if rid not in grouped: grouped[rid] = []
            grouped[rid].append(item['cds'])
            
        current_x = 20
        # Move genes lower to make room for product name: y_mid was 40, now 50
        y_mid = (CLIPBOARD_HEIGHT - 20) / 2 + 10
        scale = self.app.zoom_scale.get() / 1000.0
        
        for rid, cds_list in grouped.items():
            # Find record name and product name
            rec_name = "Unknown Cluster"
            prod_name = ""
            if rid == "CUSTOM_PROTEINS":
                rec_name = "Custom Proteins"
            else:
                for r in self.app.gbk_records:
                    if r.uid == rid:
                        rec_name = r.annotations.get('organism', r.id)
                        prod_name = self.app.cluster_products.get(rid, "")
                        break
            
            # Draw Group Label (organism name) - moved higher (y=8)
            self.app.clipboard_canvas.create_text(current_x, 8, text=rec_name, anchor="nw", font=("tahoma", 8, "bold"))
            
            # Draw Product Name in blue below organism name (y=20)
            if prod_name:
                self.app.clipboard_canvas.create_text(current_x, 20, text=prod_name, anchor="nw", font=("tahoma", 7), fill=COLOR_PRODUCT_NAME)
            
            # Draw genes
            # Simple linear layout: 10px padding between genes
            gene_x = current_x
            
            for cds in cds_list:
                length_bp = cds['end'] - cds['start']
                width_px = length_bp * scale
                if width_px < 5: width_px = 5
                
                x1 = gene_x
                x2 = gene_x + width_px
                
                # Determine Color
                canvas_tag = f"clip_{cds['unique_tag']}"
                
                # Map for alignment queries
                self.app.visible_canvas_tag_map[canvas_tag] = cds['unique_tag']
                
                color = COLOR_DEFAULT_GENE
                if canvas_tag in self.app.selected_gene_canvas_tags:
                    # Use group color or individual gene color based on toggle
                    if self.app.var_color_by_gene.get():
                        color = self.app.query_gene_colors.get(canvas_tag, COLOR_QUERY_DEFAULT)
                    else:
                        color = self._get_group_color(canvas_tag)
                
                # Use imported utility function
                poly = draw_arrow_poly(self.app.clipboard_canvas, x1, x2, y_mid, 1, color, canvas_tag) # Always strand 1 for simplicity
                
                tip = f"Product: {cds['product']}\nSize: {cds['end']-cds['start']} bp"
                alt_name = self.app.gene_alt_names.get(cds['unique_tag'])
                if alt_name:
                    tip += f"\nAlt. name: {alt_name}"
                self.app.clipboard_canvas.tag_bind(poly, "<Enter>", lambda e, t=tip, c=self.app.clipboard_canvas, i=poly: self.app._on_gene_enter(e, t, c, i))
                self.app.clipboard_canvas.tag_bind(poly, "<Leave>", lambda e, c=self.app.clipboard_canvas, i=poly: self.app._on_gene_leave(e, c, i))

                # Draw '*' for pending selections
                if canvas_tag in self.app.pending_alignment_tags:
                     self.app.clipboard_canvas.create_text((x1+x2)/2, y_mid + 19, text="*", font=("tahoma", 12, "bold"))

                gene_x += width_px + 5
            
            group_width = gene_x - current_x
            # Account for product name width when calculating step
            text_width = max(len(rec_name), len(prod_name)) * 7
            step = max(group_width, text_width) + 30
            
            # Divider line
            self.app.clipboard_canvas.create_line(current_x + step - 15, 5, current_x + step - 15, CLIPBOARD_HEIGHT-30, fill="#ddd")
            
            current_x += step
            
        self.app.clipboard_canvas.configure(scrollregion=(0, 0, current_x, CLIPBOARD_HEIGHT-20))
