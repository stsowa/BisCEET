"""
Event mixin for GenBankBrowser.

Handles UI events including:
- Mouse wheel scrolling and zooming
- Canvas resize events
- Drag and pan interactions
- Cluster flipping
"""
import sys
from config.constants import *

class EventMixin:
    """
    Handles mouse and keyboard events (Scroll, Zoom, Drag, Pan, Hover).
    """

    def on_vscroll(self, *args):
        """Handles vertical scrolling for both canvases."""
        self.org_canvas.yview(*args)
        self.gene_canvas.yview(*args)
        self._draw_visible_rows_only()

    def on_hscroll(self, *args):
        """Handles horizontal scrolling for gene canvases."""
        self.gene_canvas.xview(*args)
        self.special_gene_canvas.xview(*args)

    def on_mousewheel(self, event):
        delta = -event.delta if sys.platform == "darwin" else -int(event.delta / 120)
        self.org_canvas.yview_scroll(delta, "units")
        self.gene_canvas.yview_scroll(delta, "units")
        self._draw_visible_rows_only()
        return "break"

    def _on_canvas_resize(self, event):
        if hasattr(self, 'v_scroll') and hasattr(self, 'query_spacer_right'):
            w = self.v_scroll.winfo_width()
            if w > 1: 
                self.query_spacer_right.config(width=w)
                if hasattr(self, 'clipboard_spacer'):
                    self.clipboard_spacer.config(width=w)
        self._draw_visible_rows_only()

    def on_zoom(self, val):
        """Handles zoom slider changes, centering on the current view."""
        new_scale_val = float(val)
        view_w = self.gene_canvas.winfo_width()
        center_x_px = self.gene_canvas.canvasx(view_w / 2)
        
        if not hasattr(self, 'current_scale_val'):
            self.current_scale_val = DEFAULT_ZOOM
            
        old_scale = self.current_scale_val / 1000.0
        new_scale = new_scale_val / 1000.0
        self.current_scale_val = new_scale_val

        center_bp = (center_x_px - 10) / old_scale if old_scale > 0 else 0
        
        x, y = self.zoom_scale.winfo_rootx(), self.zoom_scale.winfo_rooty()
        # Tooltip removed for performance
        if self.zoom_tooltip.tip_window: self.zoom_tooltip.hidetip()
        
        self._full_redraw_sequence()
        
        new_center_x_px = 10 + (center_bp * new_scale)
        sr = self.gene_canvas.cget("scrollregion")
        if sr:
            try:
                coords = [float(x) for x in sr.split()]
                total_width = coords[2] - coords[0]
                if total_width > 0:
                    target_left_px = new_center_x_px - (view_w / 2)
                    fraction = (target_left_px - coords[0]) / total_width
                    self.gene_canvas.xview_moveto(fraction)
                    self.special_gene_canvas.xview_moveto(fraction)
            except: pass

    # --- Drag / Pan Handlers ---
    def on_drag_start(self, event, c_type):
        """Initiates dragging of a cluster row."""
        self.is_drag_active = False
        row = "special"
        state_key = "special"

        if c_type == 'main':
            abs_y = self.gene_canvas.canvasy(event.y)
            found_row = -1
            for i, layout in enumerate(self.row_layout_map):
                if layout['y'] <= abs_y < layout['y'] + layout['h']:
                    found_row = i
                    break
            if found_row != -1: 
                row = found_row
                # Use UUID
                state_key = self.sorted_records_cache[found_row].uid
            else: return
            
        self.drag_data = {
            "row": row, 
            "start_x": self.gene_canvas.canvasx(event.x),
            "state_key": state_key,
            "orig": self.row_x_offsets.get(state_key, 0)
        }

    def on_drag_move(self, event):
        if not self.drag_data: return
        self.is_drag_active = True
        cur_x = self.gene_canvas.canvasx(event.x)
        scale = self.zoom_scale.get() / 1000.0
        delta_px = cur_x - self.drag_data["start_x"]
        delta_bp = delta_px / scale if scale > 0 else 0
        
        key = self.drag_data["state_key"]
        self.row_x_offsets[key] = self.drag_data["orig"] + delta_bp
        
        if self.drag_data["row"] == "special": self._draw_special_row()
        else: self._draw_visible_rows_only()

    def on_drag_release(self, e):
        self.drag_data = None
        self.is_drag_active = False

    def on_pan_start(self, event):
        """Initiates global panning of the view."""
        self.pan_start_x = event.x
        self.pan_start_y = event.y
        self.is_panning = True

    def on_pan_move(self, event):
        if not self.is_panning: return
        dx_px = event.x - self.pan_start_x
        scale = self.zoom_scale.get() / 1000.0
        if scale > 0:
            dx_bp = dx_px / scale
            self.global_x_offset += dx_bp
            self.pan_start_x = event.x
            self._full_redraw_sequence()

    def on_pan_release(self, event):
        self.is_panning = False

    def on_double_click_main(self, event):
        abs_y = self.gene_canvas.canvasy(event.y)
        for i, layout in enumerate(self.row_layout_map):
            if layout['y'] <= abs_y < layout['y'] + layout['h']:
                self.on_cluster_flip(i)
                break

    # --- Hover Handlers ---
    def _on_gene_enter(self, event, text, canvas, item_id):
        self.gene_tooltip.showtip(text, event.x_root, event.y_root)
        canvas.itemconfig(item_id, width=3)

    def _on_gene_leave(self, event, canvas, item_id):
        self.gene_tooltip.hidetip()
        canvas.itemconfig(item_id, width=1)

    # Clipboard Panning
    def on_clip_pan_start(self, event):
        self.clip_pan_start_x = event.x

    def on_clip_pan_move(self, event):
        if not hasattr(self, 'clip_pan_start_x'): return
        dx = self.clip_pan_start_x - event.x
        self.clipboard_canvas.xview_scroll(int(dx/2), "units") # Slower scroll
        self.clip_pan_start_x = event.x

    def on_clip_mousewheel(self, event):
        delta = -event.delta if sys.platform == "darwin" else -int(event.delta / 120)
        self.clipboard_canvas.xview_scroll(delta, "units")
        return "break"
