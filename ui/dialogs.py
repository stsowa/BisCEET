"""
Dialog functions for BisCEET.

Provides popup dialog windows for:
- Note editing
- Product editing
- Gene name editing
- Custom protein sequence input
"""
import tkinter as tk
from tkinter import ttk, messagebox
import re
import uuid
from Bio.Seq import Seq
from config.constants import VALID_AMINO_ACIDS

def open_note_editor(parent, uid):
    """
    Opens a popup window to edit user notes for a specific cluster.
    Notes are saved to the parent application state.
    """
    if parent.editor_window and tk.Toplevel.winfo_exists(parent.editor_window):
         parent.editor_window.lift()
         return
         
    top = tk.Toplevel(parent)
    parent.editor_window = top
    top.title("Edit Note")
    top.geometry("400x300")
    
    # Save button packed first (at bottom)
    def save():
        content = txt.get("1.0", "end-1c")
        parent.cluster_notes[uid] = content
        parent._full_redraw_sequence() # Refresh all, including pinned row
        top.destroy()
        
    btn_frame = ttk.Frame(top)
    btn_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=5)
    ttk.Button(btn_frame, text="Save (Ctrl+Enter)", command=save).pack()

    # Text widget fills remainder
    txt = tk.Text(top, wrap="word", font=("tahoma", 10))
    txt.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=5)
    
    existing = parent.cluster_notes.get(uid, "")
    txt.insert("1.0", existing)
    txt.focus_set() # Auto-focus
    top.bind('<Control-Return>', lambda e: save()) # Ctrl+Enter to save note

def open_product_editor(parent, uid):
    """
    Opens a popup to rename the 'Product' label for a cluster.
    Useful for annotating clusters with specific compound names.
    """
    if parent.editor_window and tk.Toplevel.winfo_exists(parent.editor_window):
         parent.editor_window.lift()
         return

    top = tk.Toplevel(parent)
    parent.editor_window = top
    top.title("Edit Product Name")
    top.geometry("300x100")
    
    ttk.Label(top, text="Product Name:").pack(pady=(10,5))
    entry = ttk.Entry(top, width=40)
    entry.pack(padx=10)
    
    existing = parent.cluster_products.get(uid, "")
    entry.insert(0, existing)
    
    def save():
        parent.cluster_products[uid] = entry.get().strip()
        parent._full_redraw_sequence()
        top.destroy()
        
    ttk.Button(top, text="Save", command=save).pack(pady=10)
    entry.focus_set() # Auto-focus
    top.bind('<Return>', lambda e: save()) # Enter to save

def open_gene_name_editor(parent, unique_tag):
    """
    Opens a popup to set an alternative name for a specific gene.
    """
    if parent.editor_window and tk.Toplevel.winfo_exists(parent.editor_window):
         parent.editor_window.lift()
         return

    top = tk.Toplevel(parent)
    parent.editor_window = top
    top.title("Set Alternative Gene Name")
    top.geometry("300x100")
    
    ttk.Label(top, text="Alternative Name:").pack(pady=(10,5))
    entry = ttk.Entry(top, width=40)
    entry.pack(padx=10)
    
    existing = parent.gene_alt_names.get(unique_tag, "")
    entry.insert(0, existing)
    
    def save():
        val = entry.get().strip()
        if val:
            parent.gene_alt_names[unique_tag] = val
        elif unique_tag in parent.gene_alt_names:
            del parent.gene_alt_names[unique_tag]
            
        parent._full_redraw_sequence()
        top.destroy()
        
    ttk.Button(top, text="Save", command=save).pack(pady=10)
    entry.focus_set() # Auto-focus
    top.bind('<Return>', lambda e: save()) # Enter to save

def prompt_custom_protein(parent):
    """
    Opens a dialog to manually add a protein sequence to the clipboard.
    Validates the sequence against standard amino acid codes.
    """
    if parent.editor_window and tk.Toplevel.winfo_exists(parent.editor_window):
         parent.editor_window.lift()
         return

    top = tk.Toplevel(parent)
    parent.editor_window = top
    top.title("Add Custom Protein")
    top.geometry("400x300")
    
    ttk.Label(top, text="Protein Name:").pack(pady=(10,0))
    name_entry = ttk.Entry(top, width=40)
    name_entry.pack(pady=5)
    
    ttk.Label(top, text="Protein Sequence:").pack(pady=(10,0))
    seq_text = tk.Text(top, height=8, width=50, wrap="char")
    seq_text.pack(pady=5, padx=10)
    
    def save_custom():
        name = name_entry.get().strip()
        raw_seq = seq_text.get("1.0", "end-1c")
        
        # 1. Clean Sequence (remove whitespace, newlines)
        clean_seq = re.sub(r'\s+', '', raw_seq).upper()
        
        if not name:
            messagebox.showwarning("Input Error", "Please enter a protein name.")
            return
        if not clean_seq:
            messagebox.showwarning("Input Error", "Please enter a protein sequence.")
            return
            
        # 2. Validate Sequence
        invalid_chars = set(clean_seq) - VALID_AMINO_ACIDS
        if invalid_chars:
            messagebox.showwarning("Invalid Sequence", 
                                   f"Sequence contains invalid characters: {', '.join(invalid_chars)}\n"
                                   "Only standard amino acids (and * for stop) are allowed.")
            return
        
        # 3. Create Dummy CDS Data
        custom_uid = str(uuid.uuid4())
        unique_tag = f"custom_{custom_uid}"
        
        cds_data = {
            'start': 0,
            'end': len(clean_seq) * 3, # Fake BP length
            'strand': 1,
            'product': name,
            'protein_id': name,  # Use the user's name as the ID
            'unique_tag': unique_tag
        }
        
        # 4. Store Sequence
        parent.protein_map[unique_tag] = Seq(clean_seq)
        
        # 5. Add to Clipboard
        parent.clipboard_data.append({'rec_uid': 'CUSTOM_PROTEINS', 'cds': cds_data})
        parent._draw_clipboard()
        top.destroy()
        
    ttk.Button(top, text="Save to Gene Bench", command=save_custom).pack(pady=10)
    name_entry.focus_set() # Auto-focus
    top.bind('<Return>', lambda e: save_custom()) # Enter to save
