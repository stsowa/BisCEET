# main_window.py
"""
Main application window for BisCEET.

This module contains the GenBankBrowser class which is the core UI class,
combining multiple mixins for actions, data management, drawing, and events.
"""
import tkinter as tk
from tkinter import filedialog, ttk, messagebox, font
import multiprocessing
from typing import List, Dict, Set, Tuple, Optional

# Third-party imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# Local imports
from config.constants import *
from utils.utils import ToolTip, load_keywords
from core.data_manager import load_genbank_files
from ui.renderer import ClusterRenderer
from ui.state import (
    AlignmentState, SelectionState, ViewState, 
    ClipboardState, UserMetadata, ProjectState
)

# Mixins
from ui.mixins.data_mixin import DataMixin
from ui.mixins.drawing_mixin import DrawingMixin
from ui.mixins.event_mixin import EventMixin
from ui.mixins.action_mixin import ActionMixin

class GenBankBrowser(DataMixin, DrawingMixin, EventMixin, ActionMixin, ttk.Frame):
    """
    Main application class for BisCEET.
    Manages UI layout, data state, user interactions, and orchestrates rendering/alignment.
    """

    # ==========================================
    #           SECTION 1: INITIALIZATION
    # ==========================================

    def __init__(self, master: tk.Tk):
        super().__init__(master)
        master.title("BisCEET")
        
        # Fix for resizing: Set minsize and allow root window to expand
        master.geometry(f"{WINDOW_WIDTH}x{WINDOW_HEIGHT}")
        master.minsize(800, 600)
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        master.resizable(True, True) # Explicitly enable resizing

        # --- Data Storage ---
        self.gbk_records: List[SeqRecord] = []
        self.sorted_records_cache: List[SeqRecord] = [] 
        self.protein_map: Dict[str, Seq] = {} 
        self.cached_cds_data: Dict[str, List[Dict]] = {}  # Keys are Record UUIDs
        self.max_cluster_span = 1
        self.visible_canvas_tag_map: Dict[str, str] = {} 

        # --- Grouped State (using dataclasses) ---
        self.selection = SelectionState()      # Pinned, selected, favorites
        self.alignment = AlignmentState()       # Alignment results, colors
        self.view = ViewState()                 # Flipped, expanded, offsets
        self.clipboard = ClipboardState()       # Gene bench data
        self.metadata = UserMetadata()          # Notes, products, alt names
        self.project = ProjectState()           # Project file tracking

        # --- Drag & Interaction State ---
        self.drag_data: Optional[Dict] = None
        self.is_drag_active = False 
        self.is_panning = False
        self.pan_start_x = 0
        
        # --- Editor State ---
        self.editor_window = None  # Tracks active editor popup

        # --- Filters ---
        self.hide_genes_list = load_keywords("config/hide_genes.txt")
        self.highlight_genes_list = load_keywords("config/highlight_genes.txt")

        # --- Fonts & Helpers ---
        self.label_font = font.Font(family="Tahoma", size=8)
        self.label_font_bold = font.Font(family="Tahoma", size=8, weight="bold")
        
        self.renderer = ClusterRenderer(self)

        # --- Filter State ---
        self.filter_text = tk.StringVar()
        self.filter_text.trace_add("write", self.on_filter_change)

        # --- Build UI ---
        self.pack(fill=tk.BOTH, expand=True)
        self._create_widgets()
        self.zoom_tooltip = ToolTip(master)
        # Helper tooltip for buttons
        self.button_tooltip = ToolTip(master) 
        self.visible_buttons: List[tk.Widget] = []
        self.row_menu = tk.Menu(self, tearoff=0)
        
        # --- Keyboard Shortcuts ---
        master.bind('<Control-s>', lambda e: self.save_project())
        master.bind('<Control-e>', lambda e: self.bulk_expand(True))
        master.bind('<Control-r>', lambda e: self.bulk_expand(False))
        
        # Zoom level shortcuts
        master.bind('<Control-Key-1>', lambda e: self.zoom_scale.set(5))
        master.bind('<Control-Key-2>', lambda e: self.zoom_scale.set(6))
        master.bind('<Control-Key-3>', lambda e: self.zoom_scale.set(7))
        master.bind('<Control-Key-4>', lambda e: self.zoom_scale.set(8))
        master.bind('<Control-Key-5>', lambda e: self.zoom_scale.set(9))
        master.bind('<Control-Key-6>', lambda e: self.zoom_scale.set(10))
        master.bind('<Control-Key-7>', lambda e: self.zoom_scale.set(13))
        master.bind('<Control-Key-8>', lambda e: self.zoom_scale.set(17))
        master.bind('<Control-Key-9>', lambda e: self.zoom_scale.set(25))

    # ==========================================
    #    BACKWARD-COMPATIBLE PROPERTY ALIASES
    # ==========================================
    # These properties delegate to the new dataclass instances,
    # allowing existing code to work without modification.
    
    # --- Selection State Aliases ---
    @property
    def pinned_record(self):
        return self.selection.pinned_record
    
    @pinned_record.setter
    def pinned_record(self, value):
        self.selection.pinned_record = value
    
    @property
    def selected_rows(self):
        return self.selection.selected_rows
    
    @selected_rows.setter
    def selected_rows(self, value):
        self.selection.selected_rows = value
    
    @property
    def selected_uids(self):
        return self.selection.selected_uids
    
    @selected_uids.setter
    def selected_uids(self, value):
        self.selection.selected_uids = value
    
    @property
    def favorite_rows(self):
        return self.selection.favorite_rows
    
    @favorite_rows.setter
    def favorite_rows(self, value):
        self.selection.favorite_rows = value
    
    @property
    def favorite_uids(self):
        return self.selection.favorite_uids
    
    @favorite_uids.setter
    def favorite_uids(self, value):
        self.selection.favorite_uids = value
    
    @property
    def selected_gene_canvas_tags(self):
        return self.selection.selected_gene_tags
    
    @selected_gene_canvas_tags.setter
    def selected_gene_canvas_tags(self, value):
        self.selection.selected_gene_tags = value
    
    # --- Alignment State Aliases ---
    @property
    def current_alignment_results(self):
        return self.alignment.results
    
    @current_alignment_results.setter
    def current_alignment_results(self, value):
        self.alignment.results = value
    
    @property
    def cluster_similarity_scores(self):
        return self.alignment.similarity_scores
    
    @cluster_similarity_scores.setter
    def cluster_similarity_scores(self, value):
        self.alignment.similarity_scores = value
    
    @property
    def pending_alignment_tags(self):
        return self.alignment.pending_tags
    
    @pending_alignment_tags.setter
    def pending_alignment_tags(self, value):
        self.alignment.pending_tags = value
    
    @property
    def query_gene_colors(self):
        return self.alignment.query_colors
    
    @query_gene_colors.setter
    def query_gene_colors(self, value):
        self.alignment.query_colors = value
    
    # --- View State Aliases ---
    @property
    def flipped_rows(self):
        return self.view.flipped_rows
    
    @flipped_rows.setter
    def flipped_rows(self, value):
        self.view.flipped_rows = value
    
    @property
    def expanded_rows(self):
        return self.view.expanded_rows
    
    @expanded_rows.setter
    def expanded_rows(self, value):
        self.view.expanded_rows = value
    
    @property
    def row_x_offsets(self):
        return self.view.row_x_offsets
    
    @row_x_offsets.setter
    def row_x_offsets(self, value):
        self.view.row_x_offsets = value
    
    @property
    def global_x_offset(self):
        return self.view.global_x_offset
    
    @global_x_offset.setter
    def global_x_offset(self, value):
        self.view.global_x_offset = value
    
    @property
    def sort_mode(self):
        return self.view.sort_mode
    
    @sort_mode.setter
    def sort_mode(self, value):
        self.view.sort_mode = value
    
    @property
    def row_layout_map(self):
        return self.view.row_layout_map
    
    @row_layout_map.setter
    def row_layout_map(self, value):
        self.view.row_layout_map = value
    
    # --- Clipboard State Aliases ---
    @property
    def clipboard_data(self):
        return self.clipboard.data
    
    @clipboard_data.setter
    def clipboard_data(self, value):
        self.clipboard.data = value
    
    @property
    def is_clipboard_open(self):
        return self.clipboard.is_open
    
    @is_clipboard_open.setter
    def is_clipboard_open(self, value):
        self.clipboard.is_open = value
    
    # --- Metadata Aliases ---
    @property
    def cluster_notes(self):
        return self.metadata.cluster_notes
    
    @cluster_notes.setter
    def cluster_notes(self, value):
        self.metadata.cluster_notes = value
    
    @property
    def cluster_products(self):
        return self.metadata.cluster_products
    
    @cluster_products.setter
    def cluster_products(self, value):
        self.metadata.cluster_products = value
    
    @property
    def gene_alt_names(self):
        return self.metadata.gene_alt_names
    
    @gene_alt_names.setter
    def gene_alt_names(self, value):
        self.metadata.gene_alt_names = value
    
    # --- Project State Aliases ---
    @property
    def current_project_path(self):
        return self.project.current_path
    
    @current_project_path.setter
    def current_project_path(self, value):
        self.project.current_path = value
    
    @property
    def next_import_id(self):
        return self.project.next_import_id
    
    @next_import_id.setter
    def next_import_id(self, value):
        self.project.next_import_id = value

    # ==========================================
    #           SECTION 2: WIDGET CREATION
    # ==========================================

    def _create_widgets(self):
        """Initializes all UI components."""
        self._create_control_bar()
        self._create_query_row()
        self._create_clipboard_row() 
        self._create_main_area()
        self.gene_tooltip = ToolTip(self.gene_canvas)

    def _create_control_bar(self):
        """Creates the top toolbar with buttons and filters."""
        control_frame = ttk.Frame(self, padding=(10, 10, 10, 5))
        control_frame.pack(fill=tk.X)

        # 1. Main Menu (File operations and app info)
        menu_btn = ttk.Menubutton(control_frame, text="Menu")
        main_menu = tk.Menu(menu_btn, tearoff=0)
        menu_btn.config(menu=main_menu)
        
        main_menu.add_command(label="New Project", command=self._reset_app_state)
        main_menu.add_separator()
        main_menu.add_command(label="Import gene clusters (.gb/.gbk)", command=self.open_files)
        main_menu.add_separator()
        main_menu.add_command(label="Load Project...", command=self.load_project)
        main_menu.add_command(label="Save Project", command=lambda: self.save_project(None))
        main_menu.add_command(label="Save Project As...", command=self.save_project_as)
        main_menu.add_separator()
        main_menu.add_command(label="Export Marked Clusters...", command=self.export_marked_clusters)
        main_menu.add_command(label="Export Marked Clusters and Reference...", command=self.export_marked_and_reference)
        main_menu.add_separator()
        main_menu.add_command(label="About", command=self._show_about_dialog)
        main_menu.add_command(label="Exit", command=self.master.quit)
        
        menu_btn.pack(side=tk.LEFT, padx=(0, 10))

        self.status_label = ttk.Label(control_frame, text="No files loaded.", anchor=tk.W)
        self.status_label.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)

        # 2. Toggles & Filters (moved up since mapping controls go to title boxes)
        self.var_hide_genes = tk.BooleanVar(value=True)
        self.var_highlight_genes = tk.BooleanVar(value=True)
        self.var_limit_range = tk.BooleanVar(value=True)
        self.var_inspect_monochrome = tk.BooleanVar(value=True)  # White genes in Inspect view
        self.var_show_identity = tk.BooleanVar(value=True)       # Show identity % on mapped genes
        self.var_color_by_gene = tk.BooleanVar(value=True)       # Color by gene (vs by group)
        
        toggle_frame = ttk.Frame(control_frame)
        ttk.Checkbutton(toggle_frame, text="Show %", variable=self.var_show_identity, command=self._full_redraw_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(toggle_frame, text="Hide other", variable=self.var_hide_genes, command=self._full_redraw_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(toggle_frame, text="Highlight", variable=self.var_highlight_genes, command=self._full_redraw_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(toggle_frame, text="Monotone Inspect", variable=self.var_inspect_monochrome, command=self._full_redraw_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(toggle_frame, text="Color by Gene", variable=self.var_color_by_gene, command=self._full_redraw_sequence).pack(side=tk.LEFT, padx=5)
        
        limit_subframe = ttk.Frame(toggle_frame)
        ttk.Checkbutton(limit_subframe, text="Truncate:", variable=self.var_limit_range, command=self._full_redraw_sequence).pack(side=tk.LEFT)
        self.limit_entry = ttk.Entry(limit_subframe, width=3)
        self.limit_entry.insert(0, "8")
        self.limit_entry.bind("<Return>", lambda e: self._full_redraw_sequence())
        self.limit_entry.pack(side=tk.LEFT, padx=(2,0))
        limit_subframe.pack(side=tk.LEFT, padx=10)
        toggle_frame.pack(side=tk.LEFT, padx=10)

        # 5. Zoom Control
        zoom_frame = ttk.Frame(control_frame)
        ttk.Label(zoom_frame, text="Zoom:").pack(side=tk.LEFT, padx=5)
        self.zoom_scale = ttk.Scale(zoom_frame, from_=5.0, to=25.0, value=DEFAULT_ZOOM, orient=tk.HORIZONTAL, command=self.on_zoom)
        self.zoom_scale.pack(side=tk.LEFT)
        zoom_frame.pack(side=tk.LEFT, padx=10)

    def _create_query_row(self):
        """Fixed top row for pinned cluster."""
        # Container to align with Gene Bench (padding 10, 0, 10, 0)
        query_container = ttk.Frame(self, padding=(10, 0, 10, 0))
        query_container.pack(fill=tk.X)

        # Spacer on right to match scrollbar width
        self.query_spacer_right = ttk.Frame(query_container, width=16)
        self.query_spacer_right.pack(side=tk.RIGHT, fill=tk.Y)

        # Inner frame for content (Left side, expands) - Excludes scrollbar
        query_inner_frame = ttk.Frame(query_container)
        query_inner_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Header Bar for Reference Cluster with mapping controls
        header_frame = tk.Frame(query_inner_frame, borderwidth=1, relief="solid", bg=COLOR_HEADER_BG)
        header_frame.pack(fill=tk.X)
        tk.Label(header_frame, text=" Reference cluster", font=("tahoma", 9, "bold"), bg=COLOR_HEADER_BG).pack(side=tk.LEFT)
        
        # Spacer to match Gene Bench toggle button width (aligns controls)
        tk.Label(header_frame, text="", width=3, bg=COLOR_HEADER_BG).pack(side=tk.RIGHT)
        
        # Map to Gene Clusters (right side, aligned with Gene Bench)
        self.cutoff_entry = ttk.Entry(header_frame, width=4)
        self.cutoff_entry.insert(0, "40")
        self.cutoff_entry.pack(side=tk.RIGHT, padx=(2, 5))
        tk.Label(header_frame, text="Cutoff %:", font=("tahoma", 8), bg=COLOR_HEADER_BG).pack(side=tk.RIGHT)
        self.align_button = ttk.Button(header_frame, text="Map to Gene Clusters", command=self.run_alignment, state="disabled")
        self.align_button.pack(side=tk.RIGHT, padx=(10, 5))
        
        # Select Homologs button (also right side, left of Map button)
        self.core_cutoff_entry = ttk.Entry(header_frame, width=4)
        self.core_cutoff_entry.insert(0, "40")
        self.core_cutoff_entry.pack(side=tk.RIGHT, padx=(2, 20))  # Extra padding on right to separate from Map button
        tk.Label(header_frame, text="Cutoff %:", font=("tahoma", 8), bg=COLOR_HEADER_BG).pack(side=tk.RIGHT)
        self.core_align_button = ttk.Button(header_frame, text="Select Homologs from Marked Clusters", command=self.run_core_gene_alignment, state="disabled")
        self.core_align_button.pack(side=tk.RIGHT, padx=(10, 5))

        # Query Frame (Row Content) - No padding to attach to header
        query_frame = ttk.Frame(query_inner_frame, padding=(0, 0, 0, 5))
        query_frame.pack(fill=tk.X)
        query_frame.grid_columnconfigure(1, weight=1)
        
        # Name Panel
        # WRAPPED in Frame to match gene_frame border/height exactly
        self.special_org_frame_container = ttk.Frame(query_frame, borderwidth=1, relief="solid", height=DEFAULT_ROW_HEIGHT, width=ROW_NAME_WIDTH)
        self.special_org_frame_container.grid(row=0, column=0, sticky="ns")
        self.special_org_frame_container.pack_propagate(False) # Enforce size
        
        self.special_org_frame = tk.Canvas(self.special_org_frame_container, bg=COLOR_BG_WHITE, highlightthickness=0, borderwidth=0)
        self.special_org_frame.pack(fill=tk.BOTH, expand=True)

        # Canvas Panel
        self.special_gene_frame = ttk.Frame(query_frame, borderwidth=1, relief="solid", height=DEFAULT_ROW_HEIGHT)
        self.special_gene_frame.grid(row=0, column=1, sticky="ew", padx=(1,0))
        
        self.special_gene_canvas = tk.Canvas(self.special_gene_frame, bg=COLOR_BG_WHITE, highlightthickness=0, borderwidth=0, height=DEFAULT_ROW_HEIGHT-2)
        self.special_gene_canvas.pack(fill=tk.BOTH, expand=True)

        # Bindings
        self._bind_canvas_events(self.special_gene_canvas, is_main=False)
        self.special_org_frame.bind("<Button-3>", self.on_pinned_name_right_click)

    def _create_clipboard_row(self):
        """Collapsible clipboard area below pinned row."""
        self.clipboard_container = ttk.Frame(self, padding=(10, 0, 10, 0))
        self.clipboard_container.pack(fill=tk.X)
        
        # Spacer to match main scrollbar width (Right side)
        self.clipboard_spacer = ttk.Frame(self.clipboard_container, width=16)
        self.clipboard_spacer.pack(side=tk.RIGHT, fill=tk.Y)

        # Inner frame for content (Left side, expands)
        self.clipboard_frame = ttk.Frame(self.clipboard_container)
        self.clipboard_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Header Bar with Gene Bench mapping controls
        header_frame = tk.Frame(self.clipboard_frame, borderwidth=1, relief="solid", bg=COLOR_HEADER_BG)
        header_frame.pack(fill=tk.X)
        
        tk.Label(header_frame, text=" Gene bench", font=("tahoma", 9, "bold"), bg=COLOR_HEADER_BG).pack(side=tk.LEFT)
        
        self.btn_clipboard_toggle = tk.Button(header_frame, text="+", width=3, command=self.toggle_clipboard, bd=1, relief="raised", font=("tahoma", 8))
        self.btn_clipboard_toggle.pack(side=tk.RIGHT)
        
        # Gene Bench mapping controls
        self.bench_cutoff_entry = ttk.Entry(header_frame, width=4)
        self.bench_cutoff_entry.insert(0, "40")
        self.bench_cutoff_entry.pack(side=tk.RIGHT, padx=(2, 5))
        tk.Label(header_frame, text="Cutoff %:", font=("tahoma", 8), bg=COLOR_HEADER_BG).pack(side=tk.RIGHT)
        self.bench_align_button = ttk.Button(header_frame, text="Map to Gene Clusters", command=self.run_gene_bench_alignment, state="disabled")
        self.bench_align_button.pack(side=tk.RIGHT, padx=(10, 5))
        
        # Content Frame (Hidden by default)
        self.clipboard_content = ttk.Frame(self.clipboard_frame, borderwidth=1, relief="solid", height=CLIPBOARD_HEIGHT)
        # We don't pack it yet
        
        # Use Grid layout to position canvas and scrollbar
        self.clipboard_content.grid_columnconfigure(0, weight=1)
        self.clipboard_content.grid_rowconfigure(0, weight=1)

        self.clipboard_canvas = tk.Canvas(self.clipboard_content, bg=COLOR_CLIPBOARD_BG, height=CLIPBOARD_HEIGHT-20, highlightthickness=0)
        self.clipboard_canvas.grid(row=0, column=0, sticky="nsew")
        
        # Horizontal Scrollbar for Clipboard
        self.clipboard_h_scroll = ttk.Scrollbar(self.clipboard_content, orient=tk.HORIZONTAL, command=self.clipboard_canvas.xview)
        self.clipboard_h_scroll.grid(row=1, column=0, sticky="ew")
        self.clipboard_canvas.configure(xscrollcommand=self.clipboard_h_scroll.set)
        
        # Clipboard Bindings
        self.clipboard_canvas.bind("<Button-3>", self.on_clipboard_right_click)

    def _create_main_area(self):
        """Creates the main scrolling canvas area for clusters."""
        self.main_frame = ttk.Frame(self, padding=(10, 5, 10, 10))
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Header Bar
        header_frame = tk.Frame(self.main_frame, borderwidth=1, relief="solid", bg=COLOR_HEADER_BG)
        header_frame.grid(row=0, column=0, columnspan=2, sticky="ew")
        
        # Header title with count
        header_title_frame = tk.Frame(header_frame, bg=COLOR_HEADER_BG)
        header_title_frame.pack(side=tk.LEFT)
        tk.Label(header_title_frame, text=" Gene clusters", font=("tahoma", 9, "bold"), bg=COLOR_HEADER_BG).pack(side=tk.LEFT)
        self.cluster_count_label = tk.Label(header_title_frame, text="", font=("tahoma", 9), bg=COLOR_HEADER_BG)
        self.cluster_count_label.pack(side=tk.LEFT)
        
        # Filter Box
        ttk.Entry(header_frame, textvariable=self.filter_text, width=20).pack(side=tk.RIGHT, padx=2)
        tk.Label(header_frame, text="Filter: ", font=("tahoma", 8), bg=COLOR_HEADER_BG).pack(side=tk.RIGHT)

        # Configure grid weights
        self.main_frame.grid_rowconfigure(1, weight=1) # Canvas row gets weight
        self.main_frame.grid_columnconfigure(0, minsize=ROW_NAME_WIDTH, weight=0)
        self.main_frame.grid_columnconfigure(1, weight=1)

        self.org_canvas = tk.Canvas(self.main_frame, bg=COLOR_BG_WHITE, highlightthickness=0, borderwidth=0, width=ROW_NAME_WIDTH)
        self.org_canvas.grid(row=1, column=0, sticky="ns", pady=(0,10))

        self.gene_canvas = tk.Canvas(self.main_frame, bg=COLOR_BG_WHITE, highlightthickness=0, borderwidth=0)
        self.gene_canvas.grid(row=1, column=1, sticky="nsew", padx=(1,0), pady=(0,10))

        # Scrollbars
        self.v_scroll = ttk.Scrollbar(self.main_frame, orient=tk.VERTICAL, command=self.on_vscroll)
        self.v_scroll.grid(row=1, column=2, sticky="ns", pady=(0,10))
        self.h_scroll = ttk.Scrollbar(self.main_frame, orient=tk.HORIZONTAL, command=self.on_hscroll)
        self.h_scroll.grid(row=2, column=1, sticky="ew", padx=(0, 10), pady=(0,10))
        
        # FIX: Only bind scrollbar set to gene_canvas to prevent conflicts with org_canvas
        self.org_canvas.configure(yscrollcommand=None) 
        self.gene_canvas.configure(yscrollcommand=self.v_scroll.set, xscrollcommand=self.h_scroll.set)
        self.special_gene_canvas.configure(xscrollcommand=self.h_scroll.set)
        
        # Bindings
        self._bind_canvas_events(self.gene_canvas, is_main=True)
        self.org_canvas.bind("<MouseWheel>", self.on_mousewheel)
        self.org_canvas.bind("<Button-3>", self.on_row_right_click)
        self.gene_canvas.bind("<Configure>", self._on_canvas_resize)

    def _bind_canvas_events(self, canvas, is_main=True):
        """Binds mouse events for scrolling, dragging, and clicking."""
        canvas.bind("<MouseWheel>", self.on_mousewheel)
        
        # Dragging (Cluster Shift)
        row_type = 'main' if is_main else 'special'
        canvas.bind("<Button-1>", lambda e: self.on_drag_start(e, row_type))
        canvas.bind("<B1-Motion>", self.on_drag_move)
        canvas.bind("<ButtonRelease-1>", self.on_drag_release)
        
        # Panning (Middle Mouse)
        canvas.bind("<Button-2>", self.on_pan_start)
        canvas.bind("<B2-Motion>", self.on_pan_move)
        canvas.bind("<ButtonRelease-2>", self.on_pan_release)

        if is_main:
            canvas.bind("<Button-3>", self.on_row_right_click)
            canvas.bind("<Double-Button-1>", self.on_double_click_main)
        else:
            canvas.bind("<Button-3>", self.on_pinned_right_click)
            canvas.bind("<Double-Button-1>", lambda e: self.on_cluster_flip("special"))

    def _show_about_dialog(self):
        """Shows the About dialog with version and copyright information."""
        about_win = tk.Toplevel(self.master)
        about_win.title("About BisCEET")
        about_win.resizable(False, False)
        about_win.transient(self.master)
        about_win.grab_set()
        
        # Center on parent
        about_win.geometry("+%d+%d" % (self.master.winfo_x() + 100, self.master.winfo_y() + 100))
        
        about_text = """BisCEET v1.0.0

Biosynthetic Cluster Environment Examination Tool

A Python/Tkinter application for visualizing, comparing, 
and analyzing gene clusters from GenBank files.

© 2024 BisCEET Contributors
Licensed under the MIT License

https://github.com/yourusername/BisCEET"""
        
        frame = ttk.Frame(about_win, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(frame, text=about_text, justify=tk.CENTER).pack(pady=10)
        ttk.Button(frame, text="OK", command=about_win.destroy).pack(pady=10)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    root = tk.Tk()
    app = GenBankBrowser(root)
    root.mainloop()
