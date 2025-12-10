# utils.py
"""
Utility functions and classes for UI rendering and file handling.
"""
import tkinter as tk
import os
from config.constants import *

class ToolTip:
    """Creates a pop-up tooltip for a given tkinter widget."""
    def __init__(self, widget: tk.Widget):
        self.widget = widget
        self.tip_window = None
        self.tip_label = None
        self.hide_timer = None
        self.text = ""

    def showtip(self, text: str, event_x_root: int, event_y_root: int):
        self.text = text
        x = event_x_root + 20
        y = event_y_root + 10

        if self.tip_window:
            if self.tip_label: 
                self.tip_label.config(text=self.text)
            self.tip_window.wm_geometry(f"+{int(x)}+{int(y)}")
            return

        try:
            self.tip_window = tk.Toplevel(self.widget)
            self.tip_window.wm_overrideredirect(True)
            self.tip_window.wm_geometry(f"+{int(x)}+{int(y)}")
            self.tip_label = tk.Label(
                self.tip_window, text=self.text, justify=tk.LEFT,
                background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                font=("tahoma", "8", "normal")
            )
            self.tip_label.pack(ipadx=1)
        except tk.TclError:
            self.tip_window = None

    def hidetip(self):
        if self.hide_timer:
            self.widget.after_cancel(self.hide_timer)
            self.hide_timer = None
        if self.tip_window:
            try: 
                self.tip_window.destroy()
            except tk.TclError: 
                pass
            self.tip_window = None

def load_keywords(filename):
    """
    Loads keywords from a file, ignoring empty lines.
    Checks both relative path and current working directory.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(base_dir, filename)
    
    keywords = []
    if os.path.exists(full_path):
        try:
            with open(full_path, 'r') as f:
                for line in f:
                    k = line.strip().lower()
                    if k: keywords.append(k)
        except Exception as e:
            print(f"Warning: Could not read {full_path}: {e}")
    elif os.path.exists(filename): # Fallback to cwd
            try:
                with open(filename, 'r') as f:
                    for line in f:
                        k = line.strip().lower()
                        if k: keywords.append(k)
            except Exception: pass
            
    return keywords

def draw_arrow_poly(canvas, x1, x2, y_mid, strand, color, tags):
    """
    Draws a gene arrow on the canvas.
    
    Args:
        strand: 1 (forward), -1 (reverse), or 0 (unknown).
    """
    h = GENE_ARROW_HEIGHT
    w = x2 - x1
    head = min(w * GENE_ARROW_HEAD_RATIO, GENE_ARROW_HEAD_MAX)
    y1, y2 = y_mid - h/2, y_mid + h/2
    
    if w < 2: 
        return canvas.create_line(x1, y1, x1, y2, fill=color, width=2, tags=tags)
    
    if strand == 1:
        pts = [x1, y1, x2-head, y1, x2, y_mid, x2-head, y2, x1, y2]
    elif strand == -1:
        pts = [x1+head, y1, x2, y1, x2, y2, x1+head, y2, x1, y_mid]
    else:
        pts = [x1, y1, x2, y1, x2, y2, x1, y2]
        
    return canvas.create_polygon(pts, fill=color, outline=COLOR_ARROW_OUTLINE, width=1, tags=tags)

def draw_stacked_labels(canvas, tasks, label_font, label_font_bold):
    """
    Draws stacked text labels with connector lines to avoid overlap.
    Uses a simple greedy algorithm to find the first available vertical level.
    """
    levels = []
    for x1, x2, y_mid, text, is_bold in tasks:
        font_obj = label_font_bold if is_bold else label_font
        text_w = font_obj.measure(text) + 6
        cx = (x1 + x2) / 2
        lx1, lx2 = cx - (text_w / 2), cx + (text_w / 2)
        
        found_lvl = -1
        for i, lvl_end in enumerate(levels):
            if lx1 > lvl_end + 5: 
                found_lvl = i
                levels[i] = lx2
                break
        if found_lvl == -1:
            if len(levels) < MAX_LABEL_STACK:
                found_lvl = len(levels)
                levels.append(lx2)
            else:
                continue 

        lbl_y = y_mid - 15 - (found_lvl * LABEL_STEP_HEIGHT)
        canvas.create_line(cx, y_mid - 10, cx, lbl_y + 5, fill=COLOR_CONNECTOR_LINE)
        canvas.create_text(cx, lbl_y, text=text, font=font_obj, fill="black", anchor="s")

def draw_core_bracket(canvas, min_x, max_x, y_mid):
    """Draws a bracket indicating the core gene range."""
    by = y_mid + 21
    canvas.create_line(min_x, by, max_x, by, width=2, fill="black")
    canvas.create_line(min_x, by, min_x, by - 5, width=2, fill="black")
    canvas.create_line(max_x, by, max_x, by - 5, width=2, fill="black")
