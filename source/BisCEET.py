#!/usr/bin/env python3
"""
BisCEET - Biosynthetic Cluster Environment Examination Tool

A Python/Tkinter application for visualizing, comparing, and analyzing 
gene clusters from GenBank files.

Entry point for the application. Initializes the main window and starts 
the Tkinter event loop.

Usage:
    python BisCEET.py
"""
import tkinter as tk
import multiprocessing
import os
from ui.main_window import GenBankBrowser

if __name__ == "__main__":
    # Ensure multiprocessing works correctly on Windows (frozen executable support)
    multiprocessing.freeze_support()
    
    root = tk.Tk()
    
    # Set custom window icon
    icon_path = os.path.join(os.path.dirname(__file__), "icon", "icon.png")
    if os.path.exists(icon_path):
        icon = tk.PhotoImage(file=icon_path)
        root.iconphoto(True, icon)
    
    app = GenBankBrowser(root)
    root.mainloop()

