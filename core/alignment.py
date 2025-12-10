# alignment.py
"""
Handles pairwise sequence alignment logic using Biopython.
Designed to run in worker processes for performance.
"""
from Bio import Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from typing import List, Dict, Tuple

def process_alignment_chunk(targets: List[Tuple[str, Seq]], 
                            query_data: Dict[str, Seq], 
                            cutoff: float) -> Dict[str, Tuple[str, float]]:
    """
    Worker function to align a chunk of target sequences against query sequences.
    
    Args:
        targets: List of (unique_tag, sequence) tuples.
        query_data: Dict mapping query_tag -> sequence.
        cutoff: Minimum percent identity to report a hit.
        
    Returns:
        dict: { target_tag: (best_query_tag, percent_identity) }
    """
    # Re-instantiate aligner per process to ensure thread safety
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -20
    aligner.extend_gap_score = -5
    
    local_results = {}
    
    for t_ut, t_seq in targets:
        hits = []
        len_t = len(t_seq)
        
        for q_ct, q_seq in query_data.items():
            len_q = len(q_seq)
            
            # Heuristic: Skip if length difference is > 50%
            if abs(len_q - len_t) > max(len_q, len_t) * 0.5:
                continue
            
            try:
                alignment = next(aligner.align(q_seq, t_seq))
                
                # Biopython version compatibility check
                try:
                    matches = alignment.counts().identities
                except AttributeError:
                    s1, s2 = alignment[0], alignment[1]
                    matches = sum(1 for a, b in zip(s1, s2) if a == b)
                
                pct = (matches / len_q) * 100
                
                if pct >= cutoff:
                    hits.append((q_ct, pct))
                        
            except Exception:
                continue
        
        if hits:
            local_results[t_ut] = hits
            
    return local_results
