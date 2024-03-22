from Bio import SeqIO
import os
import subprocess
from collections import defaultdict
from Bio.Blast import NCBIXML

def blast(sequence, genome, output_dir):
    command1 = f"makeblastdb -in {genome} -dbtype nucl -out {os.path.splitext(os.path.basename(genome))[0]}"
    subprocess.check_call(command1, shell=True)
    command2 = f"tblastn -db {os.path.splitext(os.path.basename(genome))[0]} -query {sequence} -out {output_dir}/{os.path.splitext(os.path.basename(genome))[0]}.out -outfmt 5"
    subprocess.check_call(command2, shell=True)
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.nhr")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.nin")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.nsq")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.ndb")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.njs")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.not")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.ntf")
    os.remove(f"{os.path.splitext(os.path.basename(genome))[0]}.nto")


def blast_result_analysis(blast_file):
    genome_protein_scaffold = defaultdict(list)
    sequence = os.path.splitext(os.path.basename(blast_file))[0]
    result_handle = open(f"{blast_file}")
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for description in blast_record.descriptions:
            scaffold_id = description.title.split(' ', 2)[1]
            query = blast_record.query
            break
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                start = hsp.sbjct_start
                end = hsp.sbjct_end
                break
            break
        genome_protein_scaffold[sequence].append([query, scaffold_id, start, end])
    return genome_protein_scaffold


def extract_gene_cluster(genome_protein_scaffold, genome, output, length : int):
    for key in genome_protein_scaffold.keys():
        scaffold_query_start_end = defaultdict(list)
        for data in genome_protein_scaffold[key]:
            scaffold_query_start_end[data[1]].append([data[0], data[2], data[3]])
        scaffold_ls = list(scaffold_query_start_end.keys())
        genome_id = os.path.splitext(os.path.basename(genome))[0]
        id_scaffold = {}
        if key == genome_id:
            for seq_record in SeqIO.parse(f"{genome}", "fasta"):
                if seq_record.id in scaffold_ls:
                    for query in scaffold_query_start_end[seq_record.id]:
                        if query[1] - length <= 0 and query[2] + length >= len(seq_record.seq):
                            id_scaffold[f'{query[0]}'] = seq_record.seq
                        elif query[1] - length <= 0 and query[2] + length < len(seq_record.seq):
                            id_scaffold[f'{query[0]}'] = seq_record.seq[:query[2] + length]
                        elif query[1] - length > 0 and query[2] + length >= len(seq_record.seq):
                            id_scaffold[f'{query[0]}'] = seq_record.seq[query[1] - length:]
                        else:
                            id_scaffold[f'{query[0]}'] = seq_record.seq[query[1] - length:query[2] + length]
        with open(f'{output}/{key}.fasta', 'w') as f:
            for id in id_scaffold.keys():
                f.write(f'>{id}\n')
                f.write(f'{id_scaffold[id]}\n')


