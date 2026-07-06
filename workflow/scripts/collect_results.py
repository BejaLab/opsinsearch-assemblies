import csv
import os
import re
from pathlib import Path
from Bio import SeqIO
from opsintools.classes import Hmmer

# Extract variables from the Snakemake object
fasta_in = snakemake.input['fasta']
blast_in = snakemake.input['blast']
opsinmap_dir = snakemake.input['opsinmap']
output_file = snakemake.output[0]

positions_param = snakemake.params['pos']
hmm_score_thresh = float(snakemake.params['hmm_score'])
blast_evalue_thresh = float(snakemake.params['blast_evalue'])
mapping_score_thresh = float(snakemake.params['mapping_score'])

def parse_blast(blast_files):
    blast_results = { }
    family_re = re.compile(pattern = r'/([A-Za-z0-9]{1,2}:[^/]+)/')
    
    for i, b_file in enumerate(blast_files):
        blast_results[i] = { }
        seen = set()
        with open(file = b_file, mode = 'r') as f:
            reader = csv.reader(f, delimiter = '\t')
            for row in reader:
                if len(row) < 13:
                    raise ValueError(f"Expected 13 columns in BLAST file {b_file}, but found {len(row)} at row: {row}")
                    
                qseqid_raw = row[0]
                qseqid = qseqid_raw.split(sep = '/')[0]
                
                if qseqid in seen:
                    continue
                seen.add(qseqid)
                
                sseqid = row[1]
                pident = row[2]
                length = row[3]
                evalue = row[10]
                bitscore = row[11]
                stitle = row[12]
                
                family = ""
                if i == 0:
                    fam_match = family_re.search(string = stitle)
                    if fam_match:
                        family = fam_match.group(1)
                
                blast_results[i][qseqid] = { 'sseqid': sseqid, 'stitle': stitle, 'family': family, 'length': length, 'evalue': evalue, 'bitscore': bitscore, 'pident': pident }
                
    return blast_results

def parse_residues(mapping_file, positions, score_thresh):
    residue_data = { }
    with open(file = mapping_file, mode = 'r') as f:
        reader = csv.DictReader(f, delimiter = '\t')
        for row in reader:
            query = row['query']
            ref = row['ref']
            ref_pos = int(row['ref_pos'])
            query_res = row['query_res']
            query_score = float(row['query_score'])
            
            if query_score >= score_thresh:
                if ref in positions and ref_pos in positions[ref]:
                    idx = positions[ref].index(ref_pos)
                    if query not in residue_data:
                        residue_data[query] = { }
                    residue_data[query][idx] = query_res
                
    return residue_data

def parse_hmm(hmm_files):
    hmm_best = { }
    
    for hmm_file in hmm_files:
        hmmer = Hmmer.Hmmer(search = hmm_file, profile_file = None, profile_cons = None)
        
        for match in hmmer.matches:
            seq_name = match['seq_name']
            score = match['full']['score']
            evalue = match['full']['evalue']
            hmm_name = match['hmm_name']
            
            if seq_name not in hmm_best or score > hmm_best[seq_name]['score']:
                hmm_best[seq_name] = { 'score': score, 'hmm_name': hmm_name, 'evalue': evalue }
                
    return hmm_best

primary_ref = list(positions_param.keys())[0]

mapping_in = os.path.join(opsinmap_dir, 'mapping.tsv')
if not os.path.exists(path = mapping_in):
    raise FileNotFoundError(f"Mapping file not found at expected path: {mapping_in}")

hmm_files_paths = list(Path(opsinmap_dir).rglob(pattern = 'hmmsearch.txt'))
if not hmm_files_paths:
    raise FileNotFoundError(f"No hmmsearch.txt files found in subdirectories of {opsinmap_dir}")
hmm_in = [str(p) for p in hmm_files_paths]

blast_data = parse_blast(blast_files = blast_in)
residue_data = parse_residues(mapping_file = mapping_in, positions = positions_param, score_thresh = mapping_score_thresh)
hmm_data = parse_hmm(hmm_files = hmm_in)

target_orfs = set(hmm_data.keys()) | set(residue_data.keys())
for i in blast_data:
    target_orfs.update(blast_data[i].keys())

fasta_index = SeqIO.index(filename = fasta_in, format = "fasta")

filtered_orfs = set()
for orf in target_orfs:
    if orf not in fasta_index:
        raise KeyError(f"Sequence '{orf}' referenced in results but missing from FASTA index.")
        
    pass_hmm = False
    pass_blast = False
    
    if orf in hmm_data and float(hmm_data[orf]['score']) >= hmm_score_thresh:
        pass_hmm = True
        
    if 0 in blast_data and orf in blast_data[0] and float(blast_data[0][orf]['evalue']) <= blast_evalue_thresh:
        pass_blast = True
        
    if pass_hmm or pass_blast:
        filtered_orfs.add(orf)

res_headers = [f"Res_{pos}" for pos in positions_param[primary_ref]]

blast_headers = []
for i, b_file in enumerate(blast_in):
    stem = Path(b_file).stem
    blast_headers.append(f"{stem}:")
    
    if i == 0:
        blast_headers.extend(["Name", "Desc", "Family", "Length", "Evalue", "Bitscore", "Pident"])
    else:
        blast_headers.extend(["Name", "Desc", "Length", "Evalue", "Bitscore", "Pident"])

header = ["ORF_Name", "Contig", "Start", "End", "Strand", "ORF_Length", "Sequence", "HMM_Profile", "HMM_Score"] + res_headers + blast_headers
header_re = re.compile(pattern = r'\[(\d+)\s+-\s+(\d+)\]')

all_rows = []

for orf in filtered_orfs:
    record = fasta_index[orf]
    match = header_re.search(string = record.description)
    if not match:
        raise ValueError(f"Could not parse coordinates from FASTA description: {record.description}")
        
    start = int(match.group(1))
    end = int(match.group(2))
    strand = '+' if start < end else '-'
    contig = record.id.rsplit(sep = '_', maxsplit = 1)[0]
    seq_len = len(record.seq)
    
    row = [
        orf,
        contig,
        start,
        end,
        strand,
        seq_len,
        str(record.seq)
    ]
    
    hmm_match = hmm_data.get(orf, { 'hmm_name': '', 'score': '', 'evalue': '' })
    row.extend([hmm_match['hmm_name'], hmm_match['score']])
    
    res_dict = residue_data.get(orf, { })
    for idx in range(len(positions_param[primary_ref])):
        row.append(res_dict.get(idx, ""))
        
    for i in range(len(blast_in)):
        b_match = blast_data[i].get(orf, { })
        
        row.append("") 
        
        if i == 0:
            row.extend([
                b_match.get('sseqid', ''),
                b_match.get('stitle', ''),
                b_match.get('family', ''),
                b_match.get('length', ''),
                b_match.get('evalue', ''),
                b_match.get('bitscore', ''),
                b_match.get('pident', '')
            ])
        else:
            row.extend([
                b_match.get('sseqid', ''),
                b_match.get('stitle', ''),
                b_match.get('length', ''),
                b_match.get('evalue', ''),
                b_match.get('bitscore', ''),
                b_match.get('pident', '')
            ])
            
    all_rows.append(row)

all_rows.sort(key = lambda r: float(r[8]) if r[8] != "" else float('-inf'), reverse = True)

with open(file = output_file, mode = 'w', newline = '') as out_f:
    writer = csv.writer(out_f, delimiter = '\t')
    writer.writerow(header)
    writer.writerows(all_rows)
