# state.py
"""
Dataclass definitions for grouping related application state.
Reduces complexity of the main GenBankBrowser class by organizing
40+ instance variables into focused, cohesive state containers.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple, Optional, Any


@dataclass
class AlignmentState:
    """Holds all state related to sequence alignment operations."""
    results: Dict[str, Tuple[str, float]] = field(default_factory=dict)
    similarity_scores: Dict[str, float] = field(default_factory=dict)
    pending_tags: Set[str] = field(default_factory=set)
    query_colors: Dict[str, str] = field(default_factory=dict)
    
    def clear(self):
        """Reset all alignment state."""
        self.results.clear()
        self.similarity_scores.clear()
        self.pending_tags.clear()
        self.query_colors.clear()


@dataclass  
class SelectionState:
    """Holds state for row and gene selection."""
    pinned_record: Optional[Any] = None  # SeqRecord, but avoid circular import
    selected_rows: List[Any] = field(default_factory=list)
    selected_uids: Set[str] = field(default_factory=set)
    favorite_rows: List[Any] = field(default_factory=list)
    favorite_uids: Set[str] = field(default_factory=set)
    selected_gene_tags: Set[str] = field(default_factory=set)
    
    def clear(self):
        """Reset all selection state."""
        self.pinned_record = None
        self.selected_rows.clear()
        self.selected_uids.clear()
        self.favorite_rows.clear()
        self.favorite_uids.clear()
        self.selected_gene_tags.clear()
    
    def add_selected(self, record) -> bool:
        """Add a record to selected rows. Returns True if added."""
        if record.uid not in self.selected_uids:
            self.selected_uids.add(record.uid)
            self.selected_rows.append(record)
            return True
        return False
    
    def remove_selected(self, record):
        """Remove a record from selected rows."""
        self.selected_uids.discard(record.uid)
        self.selected_rows = [r for r in self.selected_rows if r.uid != record.uid]
    
    def toggle_selected(self, record) -> bool:
        """Toggle selection state. Returns True if now selected."""
        if record.uid in self.selected_uids:
            self.remove_selected(record)
            return False
        else:
            self.add_selected(record)
            return True
    
    def add_favorite(self, record) -> bool:
        """Add a record to favorites. Returns True if added."""
        if record.uid not in self.favorite_uids:
            self.favorite_uids.add(record.uid)
            self.favorite_rows.append(record)
            return True
        return False
    
    def remove_favorite(self, record):
        """Remove a record from favorites."""
        self.favorite_uids.discard(record.uid)
        self.favorite_rows = [r for r in self.favorite_rows if r.uid != record.uid]
    
    def toggle_favorite(self, record) -> bool:
        """Toggle favorite state. Returns True if now favorited."""
        if record.uid in self.favorite_uids:
            self.remove_favorite(record)
            return False
        else:
            self.add_favorite(record)
            return True
    
    def is_selected(self, record) -> bool:
        """Check if record is selected."""
        return record.uid in self.selected_uids
    
    def is_favorite(self, record) -> bool:
        """Check if record is a favorite."""
        return record.uid in self.favorite_uids


@dataclass
class ViewState:
    """Holds state for visual display options."""
    flipped_rows: Set[str] = field(default_factory=set)
    expanded_rows: Set[str] = field(default_factory=set)
    row_x_offsets: Dict[str, float] = field(default_factory=dict)
    global_x_offset: float = 0.0
    sort_mode: str = "default"
    row_layout_map: List[Dict[str, float]] = field(default_factory=list)
    
    def clear(self):
        """Reset all view state."""
        self.flipped_rows.clear()
        self.expanded_rows.clear()
        self.row_x_offsets.clear()
        self.global_x_offset = 0.0
        self.sort_mode = "default"
        self.row_layout_map.clear()
    
    def is_flipped(self, key: str) -> bool:
        """Check if a row is flipped."""
        return key in self.flipped_rows
    
    def toggle_flip(self, key: str) -> bool:
        """Toggle flip state. Returns True if now flipped."""
        if key in self.flipped_rows:
            self.flipped_rows.remove(key)
            return False
        else:
            self.flipped_rows.add(key)
            return True
    
    def is_expanded(self, key: str) -> bool:
        """Check if a row is expanded."""
        return key in self.expanded_rows
    
    def toggle_expand(self, key: str) -> bool:
        """Toggle expand state. Returns True if now expanded."""
        if key in self.expanded_rows:
            self.expanded_rows.remove(key)
            return False
        else:
            self.expanded_rows.add(key)
            return True


@dataclass
class ClipboardState:
    """Holds state for the gene bench/clipboard."""
    data: List[Dict] = field(default_factory=list)
    is_open: bool = False
    
    def clear(self):
        """Reset clipboard state."""
        self.data.clear()
        self.is_open = False
    
    def add_gene(self, rec_uid: str, cds_data: Dict) -> bool:
        """Add a gene to clipboard. Returns True if added (not duplicate)."""
        for item in self.data:
            if item['cds']['unique_tag'] == cds_data['unique_tag']:
                return False  # Already present
        self.data.append({'rec_uid': rec_uid, 'cds': cds_data})
        return True
    
    def remove_gene(self, unique_tag: str):
        """Remove a gene from clipboard by unique_tag."""
        self.data = [x for x in self.data if x['cds']['unique_tag'] != unique_tag]


@dataclass
class UserMetadata:
    """Holds user-created annotations."""
    cluster_notes: Dict[str, str] = field(default_factory=dict)
    cluster_products: Dict[str, str] = field(default_factory=dict)
    gene_alt_names: Dict[str, str] = field(default_factory=dict)
    
    def clear(self):
        """Reset all user metadata."""
        self.cluster_notes.clear()
        self.cluster_products.clear()
        self.gene_alt_names.clear()


@dataclass
class ProjectState:
    """Holds project-level state."""
    current_path: Optional[str] = None
    next_import_id: int = 0
    temp_dir: Optional[str] = None
    
    def clear(self):
        """Reset project state (except temp_dir which needs special handling)."""
        self.current_path = None
        self.next_import_id = 0
