# constants.py
"""
Global configuration and constants for the BisCEET application.
Includes window settings, colors, and biological parameters.
"""

# ==========================================
#        CONFIGURATION & CONSTANTS
# ==========================================

# --- Window & Layout ---
WINDOW_WIDTH = 1200
WINDOW_HEIGHT = 800
DEFAULT_ROW_HEIGHT = 90  # Height of each cluster row (increased for product names)
ROW_NAME_WIDTH = 250     # Width of the left-side organism name panel
DEFAULT_ZOOM = 10.0      # Pixels per kbp (1000 bp)
LABEL_STEP_HEIGHT = 15   # Vertical spacing for stacked gene labels
MAX_LABEL_STACK = 20     # Max tiers for labels before they are hidden
CLIPBOARD_HEIGHT = 100   # Height of the clipboard panel

# --- Gene Arrow Dimensions ---
GENE_ARROW_HEIGHT = 20           # Height of gene arrow polygons
GENE_ARROW_HEAD_MAX = 20         # Maximum arrow head size in pixels
GENE_ARROW_HEAD_RATIO = 0.4      # Arrow head as ratio of gene width

# --- Canvas Scroll Settings ---
CANVAS_SCROLL_PADDING = 500      # Extra padding for scroll region width

# --- Colors ---
COLOR_BG_WHITE = "white"
COLOR_BG_SELECTED = "#d1e7dd"  # Light green background for selected rows
COLOR_BG_HOVER = "#f8f9fa"
COLOR_TEXT_FAVORITE = "#f0ad4e" # Star color for favorites
COLOR_DEFAULT_GENE = "#a0a0a0"  # Default grey for unaligned genes
COLOR_QUERY_DEFAULT = "#808080" # Fallback color for query genes
COLOR_CORE_HIGHLIGHT = "white"  # Core genes (homologs found in all) are white
COLOR_HIGHLIGHT_MATCH = "#27272a"  # Highlighted genes (by keyword) - dark grey
COLOR_ARROW_OUTLINE = "black"
COLOR_CONNECTOR_LINE = "#555555" # Line connecting label to gene
COLOR_CLIPBOARD_BG = "#FFFFFF"
COLOR_NOTE_ACTIVE = "#ffeb3b"    # Yellow indicator for active notes
COLOR_PRODUCT_NAME = "blue"
COLOR_HEADER_BG = "#E9EBF2"  # Light blue background for section headers

# --- Palette for Query Genes ---
# Distinct colors used to highlight different query genes and their homologs
AVAILABLE_COLORS = [
    "#f2b218", "#9d60f2", "#08b21e", "#cc3d55", "#67a1e5", "#789916", "#e818f2", 
    "#60f2c1", "#b24808", "#483dcc", "#81e567", "#991658", "#18c5f2", "#f2e560", 
    "#7208b2", "#3dcc6c", "#e56c67", "#163799", "#8ef218", "#f260da", "#08b29c", 
    "#cc903d", "#8b67e5", "#179916", "#f21858", "#60b6f2", "#9db208", "#b33dcc", 
    "#67e5ab", "#993616", "#1822f2", "#91f260", "#b20873", "#3dc0cc", "#e5ca67", 
    "#571699", "#18f244", "#f2606d", "#0849b2", "#9dcc3d", "#e567e0", "#169977", 
    "#f27a18", "#7860f2", "#1fb208", "#cc3d79", "#67c1e5", "#999816", "#b118f2", 
    "#60f29c"
]

# --- External URLs ---
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
NCBI_GENBANK_URL = "https://www.ncbi.nlm.nih.gov/nuccore"

# Valid characters for protein sequences (IUPAC extended)
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

# --- Canvas Tag Prefixes ---
# Used to identify gene sources in selected_gene_canvas_tags
TAG_PREFIX_REFERENCE = "cv_special_"  # Reference cluster genes
TAG_PREFIX_CLIPBOARD = "clip_"         # Gene Bench genes
TAG_PREFIX_CLUSTER = "cv_"             # General cluster genes (not reference)

