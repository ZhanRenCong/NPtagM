import os
import subprocess
import shutil
import enzyme_extraction
from io import StringIO
import gene_clusters_extraction
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from Bio.Blast import NCBIXML


def find_core_gene_by_antismash(sequence_file, taxon, annotation_software, core, output_dir, fungi_species):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    if taxon == 'bacteria':
        command = f'antismash --taxon bacteria -c {core} --output-dir {output_dir}/{name} {sequence_file} --genefinding-tool prodigal'
        subprocess.check_call(command, shell=True)
    elif taxon == 'fungi':
        if annotation_software == 'glimmerhmm':
            command = f'antismash --taxon fungi --genefinding-tool glimmerhmm -c {core} --output-dir {output_dir}/{name} {sequence_file}'
            subprocess.check_call(command, shell=True)
        elif annotation_software == 'augustus':
            command_1 = f'augustus --gff3=on --species={fungi_species} {sequence_file} > {name}.gff'
            subprocess.check_call(command_1, shell=True)
            command_2 = f'antismash --taxon fungi --genefinding-gff3 {name}.gff -c {core} --output-dir {output_dir}/{name} {sequence_file}'
            subprocess.check_call(command_2, shell=True)
            os.remove(f'{name}.gff')


def extract_protein_sequence(genbank_file):
    record = SeqIO.read(genbank_file, "genbank")
    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    protein_sequences = []
    for cds_feature in cds_features:
        if "translation" in cds_feature.qualifiers:
            protein_sequence = cds_feature.qualifiers["translation"][0]
            protein_sequences.append(protein_sequence)
    id = record.description
    return protein_sequences, id


def extract_gene_cluster_from_antismash_result(sequence_file, antismash_dir, gbk_dir):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    a = 0
    for file in os.listdir(f'{antismash_dir}/{name}'):
        if os.path.splitext(file)[1] == '.gbk':
            if os.path.splitext(file)[0] == name:
                continue
            else:
                a += 1
                proteins_sequence, id_ = extract_protein_sequence(f'{antismash_dir}/{name}/{file}')
                shutil.copy(f'{antismash_dir}/{name}/{file}', f'{gbk_dir}/{id_}_{a}.region001.gbk')


def glimmerhmm_annotation(sequence_file, output_dir):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    annotation = {}
    for seq_record in SeqIO.parse(f"{sequence_file}", "fasta"):
        with open(f'{seq_record.id}.fasta', 'w') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')
        command = f"glimmerhmm {seq_record.id}.fasta train_crypto -g -o {seq_record.id}.gff"
        subprocess.check_call(command, shell=True)
        f = open(f'{seq_record.id}.gff', 'r')
        output = f.read()
        handle = StringIO(output)
        features = enzyme_extraction.get_features_from_file(handle)[f"{seq_record.id}"]
        genome = SeqIO.read(f'{seq_record.id}.fasta', 'fasta')
        a = 1
        for feature in features:
            seq = ''
            if feature.location.strand == 1:
                for part in feature.location.parts:
                    seq += genome.seq[part.start:part.end]
            elif feature.location.strand == -1:
                for part in feature.location.parts:
                    seq += genome.seq[part.start:part.end].reverse_complement()
            annotation[f'{seq_record.id}_{a}'] = seq.translate()
            a += 1
        os.remove(f'{seq_record.id}.gff')
        os.remove(f'{seq_record.id}.fasta')
    with open(f'{output_dir}/{name}.faa', 'w') as f:
        for key in annotation.keys():
            f.write(f'>{key}\n')
            f.write(f'{annotation[key]}\n')


def prodigal_annotation(sequence_file, output_dir):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    annotation = {}
    for seq_record in SeqIO.parse(f"{sequence_file}", "fasta"):
        with open(f'{seq_record.id}.fasta', 'w') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')
        command = f'prodigal -i {seq_record.id}.fasta -f gff -o {seq_record.id}.gff'
        subprocess.check_call(command, shell=True)
        f = open(f'{seq_record.id}.gff', 'r')
        output = f.read()
        handle = StringIO(output)
        genome = SeqIO.read(f'{seq_record.id}.fasta', 'fasta')
        features = enzyme_extraction.get_features_from_file(handle)[f"{seq_record.id}"]
        a = 1
        for feature in features:
            seq = ''
            if feature.location.strand == 1:
                for part in feature.location.parts:
                    seq += genome.seq[part.start:part.end]
            elif feature.location.strand == -1:
                for part in feature.location.parts:
                    seq += genome.seq[part.start:part.end].reverse_complement()
            annotation[f'{seq_record.id}_{a}'] = seq.translate()
            a += 1
        os.remove(f'{seq_record.id}.fasta')
        os.remove(f'{seq_record.id}.gff')
    with open(f'{output_dir}/{name}.faa', 'w') as f:
        for key in annotation.keys():
            f.write(f'>{key}\n')
            f.write(f'{annotation[key]}\n')

def extract_protein_sequences(gff_file, output_file):
    name = os.path.splitext(os.path.basename(gff_file))[0]
    with open(gff_file) as gff_handle:
        in_protein_sequence = False
        protein_sequence = ""
        gene_name = 0
        for line in gff_handle:
            if line.startswith("# protein sequence = "):
                gene_name += 1
                in_protein_sequence = True
                line_ = line.split('[')[1].strip()
                line__ = line_.replace(']', '')
                protein_sequence += line__
            elif in_protein_sequence and "Evidence" in line:
                with open(output_file, "a") as output_handle:
                    output_handle.write(f">{name}_{gene_name}\n{protein_sequence}\n")
                in_protein_sequence = False
                protein_sequence = ""
            elif in_protein_sequence and line.strip():  # 非空行
                line_ = line.replace('# ', '')
                line__ = line_.replace(']', '')
                protein_sequence += line__.strip()


def augustus_annotation(sequence_file, output_dir, species):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    annotation = {}
    for seq_record in SeqIO.parse(f"{sequence_file}", "fasta"):
        with open(f'{seq_record.id}.fasta', 'w') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')
        command = f"augustus --gff3=on --species={species} {seq_record.id}.fasta > {seq_record.id}.gff"
        subprocess.check_call(command, shell=True)
        extract_protein_sequences(f'{seq_record.id}.gff', f'{seq_record.id}.faa')
        if os.path.exists(f"{seq_record.id}.faa") is True:
            for seq_record_ in SeqIO.parse(f"{seq_record.id}.faa", "fasta"):
                annotation[seq_record_.id] = seq_record_.seq
            os.remove(f'{seq_record.id}.faa')
        os.remove(f'{seq_record.id}.fasta')
        os.remove(f'{seq_record.id}.gff')
    with open(f'{output_dir}/{name}.faa', 'w') as f:
        for key in annotation.keys():
            f.write(f'>{key}\n')
            f.write(f'{annotation[key]}\n')


def gene_cluster_annotation(sequence_file, taxon, output_dir, annotation_software, fungi_species):
    if taxon == 'bacteria':
        prodigal_annotation(sequence_file, output_dir)
    elif taxon == 'fungi':
        if annotation_software == 'glimmerhmm':
            glimmerhmm_annotation(sequence_file, output_dir)
        elif annotation_software == 'augustus':
            augustus_annotation(sequence_file, output_dir, fungi_species)


def create_gbk_from_gff(gff_file, sequence_file, output_file):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    sequence_record = SeqIO.read(sequence_file, "fasta")
    gbk_record = SeqRecord(sequence_record.seq, id=sequence_record.id, name=sequence_record.name,
                           description=sequence_record.description, annotations={'molecule_type': 'DNA'})
    f = open(gff_file, 'r')
    output = f.read()
    handle = StringIO(output)
    genome = SeqIO.read(sequence_file, 'fasta')
    features = enzyme_extraction.get_features_from_file(handle)[name]
    for feature in features:
        seq = ''
        if feature.location.strand == 1:
            for part in feature.location.parts:
                seq += genome.seq[part.start:part.end]
        elif feature.location.strand == -1:
            for part in feature.location.parts:
                seq += genome.seq[part.start:part.end].reverse_complement()
        feature.qualifiers['translation'] = seq.translate()
        gbk_record.features.append(feature)
    # 写入GenBank文件
    with open(output_file, "w") as output_handle:
        SeqIO.write(gbk_record, output_handle, "genbank")
    os.remove(gff_file)
    os.remove(sequence_file)


def find_core_gene_by_users_file(sequence_file, taxon, output_dir, annotation_software, species, core_gene_folder,
                                 blast_evalue_CE, hmm_evalue_CE, BGC_length):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    gene_cluster_annotation(sequence_file, taxon, output_dir, annotation_software, species)
    if os.path.exists(f'{output_dir}/{name}.faa') is True:
        for file in os.listdir(core_gene_folder):
            if os.path.splitext(file)[1] == '.hmm':
                enzyme_extraction.hmm(f'{core_gene_folder}/{file}', f'{output_dir}/{name}.faa',
                                  f'{output_dir}/blast_hmm_result_2', hmm_evalue_CE)
            elif os.path.splitext(file)[1] == '.fasta':
                enzyme_extraction.blast(f'{core_gene_folder}/{file}', f'{output_dir}/{name}.faa',
                                    f'{output_dir}/blast_hmm_result_2', blast_evalue_CE)
            else:
                (f'{file} is not provided as right format')
        enzyme_extraction.extract_sequence_from_blast_hmm_result(f'./{output_dir}/blast_hmm_result_2',
                                                             f'./{output_dir}/core_enzyme_sequence',
                                                             f'{output_dir}/{name}.faa')
    os.remove(f'{output_dir}/{name}.faa')

    # extract gene clusters
    if os.path.getsize(f'./{output_dir}/core_enzyme_sequence/{name}.fasta') != 0:
        gene_clusters_extraction.blast(f'./{output_dir}/core_enzyme_sequence/{name}.fasta', sequence_file,
                                   f'./{output_dir}/gene_cluster_blast_2')
        genome_protein_scaffold = gene_clusters_extraction.blast_result_analysis(f'./{output_dir}/gene_cluster_blast_2/{name}.out')
        gene_clusters_extraction.extract_gene_cluster(genome_protein_scaffold, sequence_file, f'./{output_dir}/gene_cluster_2', BGC_length)
        # annotate gene_cluster
        a = 1
        for seq_record in SeqIO.parse(f"./{output_dir}/gene_cluster_2/{name}.fasta", "fasta"):
            with open(f'{seq_record.id}.fasta', 'w') as f:
                f.write(f'>{seq_record.id}\n')
                f.write(f'{seq_record.seq}\n')
            if taxon == 'bacteria':
                command = f'prodigal -i {seq_record.id}.fasta -f gff -o {seq_record.id}.gff'
                subprocess.check_call(command, shell=True)
            elif taxon == 'fungi':
                if annotation_software == 'glimmerhmm':
                    command = f"glimmerhmm {seq_record.id}.fasta train_crypto -g -o {seq_record.id}.gff"
                    subprocess.check_call(command, shell=True)
                elif annotation_software == 'augustus':
                    command = f"augustus --gff3=on --species={species} {seq_record.id}.fasta > {seq_record.id}.gff"
                    subprocess.check_call(command, shell=True)
            new_id = ''
            data = seq_record.id.split('_')
            for i in data[0:-1]:
                new_id += i + '_'
            new_id = new_id[0:-1]
            create_gbk_from_gff(f'{seq_record.id}.gff', f'{seq_record.id}.fasta',
                                f'./{output_dir}/BGC_gbk/{new_id}_{a}.region001.gbk')
            a += 1



def get_result_from_blast_hmm_result(result_file):
    gene_list = []
    if '_hmm.out' in result_file:
        with open(result_file, 'r') as f:
            for line in f.readlines():
                if line[0] == '>':
                    data = line.split(" ")
                    gene_list.append(data[1])
    elif '_blast.out' in result_file:
        result_handle = open(result_file)
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for description in blast_record.descriptions:
                gene_list.append(description.title.split(' ', 2)[1])
    gene_list = list(set(gene_list))
    return gene_list


def check_tailoring_enzyme_in_new_BGC(query, outputdir, hmm_evalue_TE, blast_evalue_TE):
    gene_cluster_data = defaultdict(dict)
    new_gene_cluster_data = defaultdict(dict)
    for file in os.listdir(f'{outputdir}/BGC_gbk'):
        proteins_sequence, useless = extract_protein_sequence(f'{outputdir}/BGC_gbk/{file}')
        id_ = os.path.splitext(os.path.basename(file))[0]
        a = 1
        with open(f'{id_}.fasta', 'w') as f:
            for i in proteins_sequence:
                f.write(f'>{id_}_{a}\n')
                f.write(f"{i}\n")
                a += 1
        gene_cluster_data[id_]['gene_list'] = []
        for file_ in os.listdir(query):
            if os.path.splitext(file_)[1] == '.hmm':
                enzyme_extraction.hmm(f'{query}/{file_}', f'{id_}.fasta', outputdir, hmm_evalue_TE)
                gene_list = get_result_from_blast_hmm_result(
                    f'{outputdir}/{id_}_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
                for gene in gene_list:
                    gene_cluster_data[id_]['gene_list'].append(gene)
                os.remove(f'{outputdir}/{id_}_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
            elif os.path.splitext(file_)[1] == '.fasta':
                enzyme_extraction.blast(f'{query}/{file_}', f'{id_}.fasta', outputdir, blast_evalue_TE)
                gene_list = get_result_from_blast_hmm_result(
                    f'{outputdir}/{id_}_{os.path.splitext(os.path.basename(file_))[0]}_blast.out')
                for gene in gene_list:
                    gene_cluster_data[id_]['gene_list'].append(gene)
                os.remove(f'{outputdir}/{id_}_{os.path.splitext(os.path.basename(file_))[0]}_blast.out')
        os.remove(f'{id_}.fasta')
    tailoring_enzyme_list = []
    for key in gene_cluster_data.keys():
        if len(gene_cluster_data[key]['gene_list']) != 0:
            new_gene_cluster_data[key]['gene_list'] = gene_cluster_data[key]['gene_list']
            shutil.copy(f'{outputdir}/BGC_gbk/{key}.gbk', f'{outputdir}/BGC_annotation/{key}.gbk')
            new_id = ''
            data = key.split('_')
            for i in data[0:-1]:
                new_id += i + '_'
            new_id = new_id[0:-1]
            tailoring_enzyme_list.append(new_id)
    tailoring_enzyme_list = list(set(tailoring_enzyme_list))
    with open(f'{outputdir}/tailoring_enzyme.fasta', 'w') as f:
        for genome in os.listdir(f'{outputdir}/tailoring_enzyme_sequence'):
            genome_id = os.path.splitext(os.path.basename(genome))[0]
            for enzyme in tailoring_enzyme_list:
                enzyme_genome = ''
                data_ = enzyme.split('_')
                for i_ in data_[0:-1]:
                    enzyme_genome += i_ + '_'
                enzyme_genome = enzyme_genome[0:-1]
                if enzyme_genome == genome_id:
                    for seq_record in SeqIO.parse(f"{outputdir}/tailoring_enzyme_sequence/{genome}", "fasta"):
                        if seq_record.id == enzyme:
                            f.write(f'>{seq_record.id}\n')
                            f.write(f'{seq_record.seq}\n')
    return new_gene_cluster_data


def statistics(new_gene_cluster_data, outputdir):
    for gbk in os.listdir(f'{outputdir}/BGC_annotation'):
        name = os.path.splitext(gbk)[0]
        record = SeqIO.read(f'{outputdir}/BGC_annotation/{gbk}', "genbank")
        features = [feature for feature in record.features if feature.type == "protocluster"]
        if features != []:
            new_gene_cluster_data[name]['type'] = features[0].qualifiers['category']
        else:
            new_gene_cluster_data[name]['type'] = "user's BGC definition"
    with open(f'{outputdir}/BGC_statistics.txt', 'w') as f:
        f.write('BGC_id\ttarget tailroing enzyme number\tBGC type\n')
        for key in new_gene_cluster_data.keys():
            f.write(f'{key}\t{len(new_gene_cluster_data[key]["gene_list"])}\t{new_gene_cluster_data[key]["type"]}\n')


def extract_sequence_from_blast_hmm_result(blast_hmm_result):
    gene_list = []
    sequence = {}
    if '_hmm.out' in blast_hmm_result:
        with open(blast_hmm_result, 'r') as f:
            for line in f.readlines():
                if line[0] == '>':
                    data = line.split(" ")
                    gene_list.append(data[1])
    elif '_blast.out' in blast_hmm_result:
        result_handle = open(blast_hmm_result)
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for description in blast_record.descriptions:
                gene_list.append(description.title.split(' ', 2)[1])
        gene_list = list(set(gene_list))
    for seq_record in SeqIO.parse("MIBIG/MIBIG.fasta", "fasta"):
        if seq_record.id in gene_list:
            sequence[seq_record.id] = seq_record.seq
    return sequence


def check_core_enzyme_in_new_BGC(core_gene, outputdir, gene_cluster_list, hmm_evalue_CE, blast_evalue_CE):
    BGC_list = []
    for BGC in gene_cluster_list:
        proteins_sequence, useless = extract_protein_sequence(f'MIBIG/MIBIG/{BGC}.gbk')
        a = 1
        with open(f'{BGC}.fasta', 'w') as f:
            for i in proteins_sequence:
                f.write(f'>{BGC}_{a}\n')
                f.write(f"{i}\n")
                a += 1
        for file_ in os.listdir(core_gene):
            if os.path.splitext(file_)[1] == '.hmm':
                enzyme_extraction.hmm(f'{core_gene}/{file_}', f'{BGC}.fasta', outputdir, hmm_evalue_CE)
                gene_list = get_result_from_blast_hmm_result(
                    f'{outputdir}/{BGC}_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
                if len(gene_list) > 0:
                    for i in gene_list:
                        enzyme_BGC = ''
                        data_ = i.split('_')
                        for i_ in data_[0:-1]:
                            enzyme_BGC += i_ + '_'
                        enzyme_BGC = enzyme_BGC[0:-1]
                        BGC_list.append(enzyme_BGC)
                os.remove(f'{outputdir}/{BGC}_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
            elif os.path.splitext(file_)[1] == '.fasta':
                enzyme_extraction.blast(f'{core_gene}/{file_}', f'{BGC}.fasta', outputdir, blast_evalue_CE)
                gene_list = get_result_from_blast_hmm_result(
                    f'{outputdir}/{BGC}_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
                if len(gene_list) > 0:
                    for i in gene_list:
                        enzyme_BGC = ''
                        data_ = i.split('_')
                        for i_ in data_[0:-1]:
                            enzyme_BGC += i_ + '_'
                        enzyme_BGC = enzyme_BGC[0:-1]
                        BGC_list.append(enzyme_BGC)
                os.remove(f'{outputdir}/{BGC}_{os.path.splitext(os.path.basename(file_))[0]}_blast.out')
        os.remove(f'{BGC}.fasta')
    BGC_list = list(set(BGC_list))
    return BGC_list


def find_tailoring_enzyme_in_MIBIG(query, outputdir, hmm_evalue_TE, blast_evalue_TE, core_gene, hmm_evalue_CE,
                                   blast_evalue_CE):
    total_sequence = {}
    for file_ in os.listdir(query):
        if os.path.splitext(file_)[1] == '.hmm':
            enzyme_extraction.hmm(f'{query}/{file_}', f'MIBIG/MIBIG.fasta', outputdir, hmm_evalue_TE)
            sequence = extract_sequence_from_blast_hmm_result(
                f'{outputdir}/MIBIG_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
            for key in sequence.keys():
                total_sequence[key] = sequence[key]
            os.remove(f'{outputdir}/MIBIG_{os.path.splitext(os.path.basename(file_))[0]}_hmm.out')
        elif os.path.splitext(file_)[1] == '.fasta':
            enzyme_extraction.blast(f'{query}/{file_}', f'MIBIG/MIBIG.fasta', outputdir, blast_evalue_TE)
            sequence = extract_sequence_from_blast_hmm_result(
                f'{outputdir}/MIBIG_{os.path.splitext(os.path.basename(file_))[0]}_blast.out')
            for key in sequence.keys():
                total_sequence[key] = sequence[key]
            os.remove(f'{outputdir}/MIBIG_{os.path.splitext(os.path.basename(file_))[0]}_blast.out')
    tailoring_enzyme_list = list(total_sequence.keys())
    TE_gene_cluster_list = []
    for key in tailoring_enzyme_list:
        enzyme_BGC = ''
        data_ = key.split('_')
        for i_ in data_[0:-1]:
            enzyme_BGC += i_ + '_'
        enzyme_BGC = enzyme_BGC[0:-1]
        TE_gene_cluster_list.append(enzyme_BGC)
    gene_cluster_list = list(set(TE_gene_cluster_list))
    if core_gene == None:
        for i in gene_cluster_list:
            shutil.copy(f'MIBIG/MIBIG/{i}.gbk', f'{outputdir}/BGC_annotation/{i}.region001.gbk')
        with open(f'{outputdir}/tailoring_enzyme_MIBIG.fasta', 'w') as f:
            for seq_record in SeqIO.parse("MIBIG/MIBIG.fasta", "fasta"):
                if seq_record.id in tailoring_enzyme_list:
                    f.write(f'>{seq_record.id}\n')
                    f.write(f'>{seq_record.seq}\n')

    else:
        BGC_list = check_core_enzyme_in_new_BGC(core_gene, outputdir, gene_cluster_list, hmm_evalue_CE, blast_evalue_CE)
        new_tailoring_enzyme_list = []
        with open(f'{outputdir}/tailoring_enzyme_MIBIG.fasta', 'w') as f:
            for i in gene_cluster_list:
                if i in BGC_list:
                    shutil.copy(f'MIBIG/MIBIG/{i}.gbk', f'{outputdir}/BGC_annotation/{i}.region001.gbk')
                    for j in tailoring_enzyme_list:
                        if i in j:
                            new_tailoring_enzyme_list.append(j)
        sequence_list = []
        with open(f'{outputdir}/tailoring_enzyme_MIBIG.fasta', 'w') as f:
            for seq_record in SeqIO.parse("MIBIG/MIBIG.fasta", "fasta"):
                if seq_record.id in new_tailoring_enzyme_list:
                    if seq_record.seq not in sequence_list:
                        f.write(f'>{seq_record.id}\n')
                        f.write(f'{seq_record.seq}\n')
                        sequence_list.append(seq_record.seq)



