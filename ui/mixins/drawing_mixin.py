"""
Drawing mixin for GenBankBrowser.

Handles rendering orchestration including:
- Full redraw sequences
- Visible row rendering
- Row expansion toggling
"""
import tkinter as tk
from config.constants import *


class DrawingMixin:
    """
    Handles rendering calls and layout updates for GenBankBrowser.
    """

    def _full_redraw_sequence(self):
        """Triggers a full update of the layout and redraws visible elements."""
        self._draw_special_row()
        self._draw_clipboard() # Update clipboard mapping
        self._update_layout_map()

        total_h = self.row_layout_map[-1]['y'] + self.row_layout_map[-1]['h'] if self.row_layout_map else 0
        scale = self.zoom_scale.get() / 1000.0
        
        # Determine scroll width
        max_drag = max(self.row_x_offsets.values(), default=0) if self.row_x_offsets else 0
        total_drag_px = (max_drag + max(0, self.global_x_offset)) * scale
        total_w = (self.max_cluster_span * scale) + total_drag_px + CANVAS_SCROLL_PADDING 
        
        self.gene_canvas.configure(scrollregion=(0, 0, total_w, total_h))
        self.special_gene_canvas.configure(scrollregion=(0, 0, total_w, DEFAULT_ROW_HEIGHT))
        self.org_canvas.configure(scrollregion=(0, 0, ROW_NAME_WIDTH, total_h))
        
        self._draw_visible_rows_only()

    def _draw_visible_rows_only(self):
        self.renderer.draw_visible_rows_only()

    def _draw_name_panel(self, rec, y_top, y_bot, y_mid, bg):
        self.renderer.draw_name_panel(rec, y_top, y_bot, y_mid, bg)

    def _draw_cluster_features(self, canvas, record, row_idx_or_special, y_mid, scale):
        self.renderer.draw_cluster_features(canvas, record, row_idx_or_special, y_mid, scale)

    def _draw_special_row(self):
        self.renderer.draw_special_row()

    def _draw_clipboard(self):
        self.renderer.draw_clipboard()

    def toggle_row_expansion(self, rec_uid):
        if rec_uid in self.expanded_rows: self.expanded_rows.remove(rec_uid)
        else: self.expanded_rows.add(rec_uid)
        self._full_redraw_sequence()

    def bulk_expand(self, expand: bool):
        if expand:
            # Expand all records (except pinned), including hidden ones
            records = [r for r in self.gbk_records if r is not self.pinned_record]
            self.expanded_rows = {r.uid for r in records}
        else: self.expanded_rows.clear()
        self._full_redraw_sequence()
