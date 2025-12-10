# data_manager.py
"""
Handles loading and parsing of GenBank files.
Extracts CDS features, translations, and manages unique IDs for genes.
"""
import os
import uuid
import logging
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def parse_features(record):
    """
    Extracts CDS features and translation from a GenBank record.
    
    Args:
        record: Bio.SeqRecord object.
        
    Returns:
        tuple: (cds_list, l_start, l_end, protein_map)
    """
    cds_list = []
    l_start = []
    l_end = []
    protein_map = {}
    
    for feat in record.features:
        if feat.type.upper() == 'CDS':
            translation = feat.qualifiers.get('translation', [None])[0]
            if translation:
                start, end = int(feat.location.start), int(feat.location.end)
                l_start.append(start)
                l_end.append(end)
                # Use UUID in gene tag to differentiate duplicates
                unique_tag = f"gene_{record.uid}_{start}_{end}"
                protein_map[unique_tag] = Seq(translation)
                cds_list.append({
                    'start': start, 'end': end,
                    'strand': feat.location.strand,
                    'product': feat.qualifiers.get('product', ['unknown'])[0],
                    'protein_id': feat.qualifiers.get('protein_id', ['unknown'])[0],
                    'unique_tag': unique_tag
                })
    return cds_list, l_start, l_end, protein_map

def load_genbank_files(files, import_id=0):
    """
    Loads and parses a list of GenBank files.
    
    Args:
        files (list): List of file paths.
        import_id (int): Unique identifier for this import session to distinguish duplicates.
        
    Returns:
        tuple: (gbk_records, cached_cds_data, protein_map, global_max_span)
    """
    gbk_records = []
    cached_cds_data = {}
    protein_map = {}
    global_max_span = 0
    
    for f in files:
        try:
            for idx, record in enumerate(SeqIO.parse(f, "genbank")):
                # Ensure ID exists (visual only)
                if not record.id or record.id == "<unknown id>":
                    record.id = os.path.basename(f)
                
                # Assign internal UUID for state tracking
                record.uid = str(uuid.uuid4())
                
                # Store source metadata for project saving/loading
                record.file_path = os.path.abspath(f)
                record.index_in_file = idx
                record.import_instance = import_id  # NEW: Track which import session this is from
                
                # Extract GenBank accession (from annotations or id)
                # Try to get from annotations first (most reliable)
                accessions = record.annotations.get('accessions', [])
                if accessions:
                    record.genbank_accession = accessions[0]
                elif hasattr(record, 'id') and record.id and record.id != "<unknown id>":
                    # Fall back to record.id which often contains the accession
                    record.genbank_accession = record.id
                else:
                    record.genbank_accession = None
                
                # Extract DEFINITION field (description)
                record.definition = record.description if hasattr(record, 'description') else None
                
                gbk_records.append(record)
                
                # Extract CDS and Translation
                cds_list, l_start, l_end, p_map = parse_features(record)
                
                # Store using UUID
                cached_cds_data[record.uid] = cds_list
                protein_map.update(p_map)
                
                # Update max span for scrolling
                if l_start:
                    span = max(l_end) - min(l_start)
                    if span > global_max_span: global_max_span = span
        except Exception as e:
            logger.warning(f"Skipping file {f} due to error: {e}")
            
    return gbk_records, cached_cds_data, protein_map, global_max_span

