import os
import subprocess
from collections import defaultdict
from Bio.Blast import NCBIXML
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from io import StringIO

MODIFY_LOCATIONS_BY_PHASE = False

def generate_details_from_subfeature(sub_feature,
                                     existing_qualifiers,
                                     locations,
                                     trans_locations):
    """ Finds the locations of a subfeature and any mismatching qualifiers

        Arguments:
            sub_feature: the GFF subfeature to work on
            existing_qualifiers: a dict of any existing qualifiers from other
                                 subfeatures
            locations: a list of any existing FeatureLocations from other
                       subfeatures
            trans_locations: a list of any existing FeatureLocations for
                             translations

        Returns:
            a set of qualifiers from the subfeature for which an existing
            qualifier existed but had a different value
    """
    mismatching_qualifiers = set()
    start = sub_feature.location.start.real
    end = sub_feature.location.end.real
    if MODIFY_LOCATIONS_BY_PHASE:
        phase = int(sub_feature.qualifiers.get('phase', [0])[0])
        if sub_feature.strand == 1:
            start += phase
        else:
            end -= phase
    locations.append(FeatureLocation(start, end, strand=sub_feature.strand))
    # Make sure CDSs lengths are multiple of three. Otherwise extend to next full codon.
    # This only applies for translation.
    modulus = (end - start) % 3
    if modulus and sub_feature.strand == 1:
        end += 3 - modulus
    elif modulus and sub_feature.strand == -1:
        start -= 3 - modulus
    trans_locations.append(FeatureLocation(start, end, strand=sub_feature.strand))
    # For split features (CDSs), the final feature will have the same qualifiers as the children ONLY if
    # they're the same, i.e.: all children have the same "protein_ID" (key and value).
    for qual in sub_feature.qualifiers:
        if qual not in existing_qualifiers:
            existing_qualifiers[qual] = sub_feature.qualifiers[qual]
        elif existing_qualifiers[qual] != sub_feature.qualifiers[qual]:
            mismatching_qualifiers.add(qual)
    return mismatching_qualifiers

def check_sub(feature):
    """ Recursively checks a GFF feature for any subfeatures and generates any
        appropriate SeqFeature instances from them.
    """
    new_features = []
    locations = []  # type: List[FeatureLocation]
    trans_locations = []  # type: List[FeatureLocation]
    qualifiers = {}  # type: Dict[str, List[str]]
    mismatching_qualifiers = set()  # type: Set[str]
    for sub in feature.sub_features:
        if sub.sub_features:  # If there are sub_features, go deeper
            new_features.extend(check_sub(sub))
        elif sub.type == 'CDS':
            sub_mismatch = generate_details_from_subfeature(sub, qualifiers,
                                                            locations, trans_locations)
            mismatching_qualifiers.update(sub_mismatch)

    for qualifier in mismatching_qualifiers:
        del qualifiers[qualifier]
    if 'Parent' in qualifiers:
        del qualifiers['Parent']

    # if nothing to work on
    if not new_features and not locations:
        return []

    # Only works in tip of the tree, when there's no new_feature built yet. If there is,
    # it means the script just came out of a check_sub and it's ready to return.
    if not new_features:
        new_loc = locations[0]
        # construct a compound location if required
        if len(locations) > 1:
            locations = sorted(locations, key=lambda x: x.start.real)
            trans_locations = sorted(trans_locations, key=lambda x: x.start.real)
            if locations[0].strand == 1:
                new_loc = CompoundLocation(locations)
            else:
                new_loc = CompoundLocation(list(reversed(locations)))
                trans_locations = list(reversed(trans_locations))
        new_feature = SeqFeature(new_loc)
        new_feature.qualifiers = qualifiers
        new_feature.type = 'CDS'
        new_features.append(new_feature)

    return new_features

def get_features_from_file(handle):
    """ Generates new SeqFeatures from a GFF file.

        Arguments:
            handle: a file handle/stream with the GFF contents

        Returns:
            a dictionary mapping record ID to a list of SeqFeatures for that record
    """
    gff_records = list(GFF.parse(handle))

    results = {}
    for gff_record in gff_records:
        features = []
        for feature in gff_record.features:
            if feature.type == 'CDS':
                new_features = [feature]
            else:
                new_features = check_sub(feature)
                if not new_features:
                    continue

            name = feature.id
            locus_tag = feature.qualifiers.get("locus_tag")

            for qtype in ["gene", "name", "Name"]:
                if qtype in feature.qualifiers:
                    name_tmp = feature.qualifiers[qtype][0]
                    # Assume name/Name to be sane if they don't contain a space
                    if " " in name_tmp:
                        continue
                    name = name_tmp
                    break

            for i, new_feature in enumerate(new_features):
                variant = name
                if len(new_features) > 1:
                    variant = "{0}_{1}".format(name, i)
                new_feature.qualifiers['gene'] = [variant]
                if locus_tag is not None:
                    new_feature.qualifiers["locus_tag"] = locus_tag
                features.append(new_feature)
        results[gff_record.id] = features
    return results

def glimmerhmm_annotation(sequence_file, output_dir):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    annotation = {}
    for seq_record in SeqIO.parse(f"{sequence_file}", "fasta"):
        with open(f'{seq_record.id}.fasta', 'w') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')
        command = f"glimmerhmm {seq_record.id}.fasta train_crypto -g"
        process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, unused_err = process.communicate()
        output = output.decode("utf-8")
        handle = StringIO(output)
        features_ = get_features_from_file(handle)
        genome = SeqIO.read(f'{seq_record.id}.fasta', 'fasta')
        if list(features_.keys()) != []:
            features = features_ [f"{seq_record.id}"]
            for feature in features:
                seq = ''
                if feature.location.strand == 1:
                    for part in feature.location.parts:
                        seq += genome.seq[part.start:part.end]
                elif feature.location.strand == -1:
                    for part in feature.location.parts:
                        seq += genome.seq[part.start:part.end].reverse_complement()
                annotation[feature.qualifiers['gene'][0]] = seq.translate()
        os.remove(f'{seq_record.id}.fasta')
    with open(f'{output_dir}/{name}.faa', 'w') as f:
        for key in annotation.keys():
            f.write(f'>{key}\n')
            f.write(f'{annotation[key]}\n')


def extract_protein_sequences(gff_file, output_file):
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
                    output_handle.write(f">g{gene_name}\n{protein_sequence}\n")
                in_protein_sequence = False
                protein_sequence = ""
            elif in_protein_sequence and line.strip():  # 非空行
                line_ = line.replace('# ', '')
                line__ = line_.replace(']', '')
                protein_sequence += line__.strip()

def augustus_annotation(sequence_file, output_dir, species):
    name = os.path.splitext(os.path.basename(sequence_file))[0]
    command = f"augustus --gff3=on --species={species} {sequence_file} > {name}.gff"
    subprocess.check_call(command, shell=True)
    extract_protein_sequences(f'{name}.gff', f'{output_dir}/{name}.faa')
    os.remove(f'{name}.gff')

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
        features_ = get_features_from_file(handle)
        genome = SeqIO.read(f'{seq_record.id}.fasta', 'fasta')
        if list(features_.keys()) != []:
            features = features_ [f"{seq_record.id}"]
            for feature in features:
                seq = ''
                if feature.location.strand == 1:
                    for part in feature.location.parts:
                        seq += genome.seq[part.start:part.end]
                elif feature.location.strand == -1:
                    for part in feature.location.parts:
                        seq += genome.seq[part.start:part.end].reverse_complement()
                annotation[feature.qualifiers['gene'][0]] = seq.translate()
        os.remove(f'{seq_record.id}.fasta')
        os.remove(f'{seq_record.id}.gff')
    with open(f'{output_dir}/{name}.faa', 'w') as f:
        for key in annotation.keys():
            f.write(f'>{key}\n')
            f.write(f'{annotation[key]}\n')

def blast(sequence, protein_file, output_dir, evalue):
    command1 = f"makeblastdb -in {protein_file} -dbtype prot -out {os.path.splitext(protein_file)[0]}"
    subprocess.check_call(command1, shell=True)
    command2 = f"blastp -db {os.path.splitext(protein_file)[0]} -query {sequence} -out {output_dir}/{os.path.splitext(os.path.basename(protein_file))[0]}_{os.path.splitext(os.path.basename(sequence))[0]}_blast.out -outfmt 5 -evalue {evalue}"
    subprocess.check_call(command2, shell=True)
    os.remove(f"{os.path.splitext(protein_file)[0]}.phr")
    os.remove(f"{os.path.splitext(protein_file)[0]}.pin")
    os.remove(f"{os.path.splitext(protein_file)[0]}.psq")
    os.remove(f"{os.path.splitext(protein_file)[0]}.pdb")
    os.remove(f"{os.path.splitext(protein_file)[0]}.pjs")
    os.remove(f"{os.path.splitext(protein_file)[0]}.pot")
    os.remove(f"{os.path.splitext(protein_file)[0]}.ptf")
    os.remove(f"{os.path.splitext(protein_file)[0]}.pto")

def hmm(hmm_file, protein_file , output_dir, evalue):
    command = f'hmmsearch -E {evalue} {hmm_file} {protein_file} > {output_dir}/{os.path.splitext(os.path.basename(protein_file))[0]}_{os.path.splitext(os.path.basename(hmm_file))[0]}_hmm.out'
    subprocess.check_call(command, shell=True)

def extract_sequence_from_blast_hmm_result(blast_hmm_result_dir, output_dir, protein_file):
    gene_list = []
    name = os.path.splitext(os.path.basename(protein_file))[0]
    for blast in os.listdir(blast_hmm_result_dir):
        name_ = os.path.splitext(os.path.basename(blast))[0]
        list_ = name_.split('_')
        name_1 = ''
        for i in list_[0:-2]:
            name_1 += f'{i}_'
        name_1 = name_1[0:-1]
        if name_1 == name:
            if '_hmm.out' in blast:
                with open(f'{blast_hmm_result_dir}/{blast}') as f:
                    for line in f.readlines():
                        if line[0] == '>':
                            data = line.split(" ")
                            gene_list.append(data[1])
            elif '_blast.out' in blast:
                result_handle = open(f"{blast_hmm_result_dir}/{blast}")
                blast_records = NCBIXML.parse(result_handle)
                for blast_record in blast_records:
                    for description in blast_record.descriptions:
                        gene_list.append(description.title.split(' ', 2)[1])
        gene_list = list(set(gene_list))
        f = open(f'{output_dir}/{name}.fasta', 'w')
        for seq_record in SeqIO.parse(f"{protein_file}", "fasta"):
            if seq_record.id in gene_list:
                f.write(f">{seq_record.id}\n")
                f.write(f"{seq_record.seq}\n")
        f.close()

def CD_Hit(protein_file_dir, output_dir, cdhit_cutoff):
    id_sequence = {}
    for protein in os.listdir(protein_file_dir):
        name = os.path.splitext(os.path.basename(protein))[0]
        a = 1
        for seq_record in SeqIO.parse(f'{protein_file_dir}/{protein}', "fasta"):
            id_sequence[f'{name}_{a}'] = seq_record.seq
            a += 1
    with open(f'{output_dir}/CD_Hit.fasta', 'w') as f:
        for id in id_sequence.keys():
            f.write(f">{id}\n")
            f.write(f"{id_sequence[id]}\n")
    command = f'cd-hit -i {output_dir}/CD_Hit.fasta -o {output_dir}/CD_Hit_{cdhit_cutoff}.fasta -c {cdhit_cutoff}'
    subprocess.check_call(command, shell=True)
    genome_sequence = defaultdict(list)
    for seq_record in SeqIO.parse(f'{output_dir}/CD_Hit_{cdhit_cutoff}.fasta', "fasta"):
        list_ = seq_record.id.split('_')
        name_ = ''
        for i in list_[0:-1]:
            name_ += f'{i}_'
        name_ = name_[0:-1]
        genome_sequence[name_].append(seq_record.seq)
    for genome in genome_sequence.keys():
        with open(f'{output_dir}/tailoring_enzyme_sequence/{genome}.fasta', 'w') as f:
            a = 1
            for sequence in genome_sequence[genome]:
                f.write(f">{genome}_{a}\n")
                f.write(f"{sequence}\n")
                a += 1
    os.remove(f'{output_dir}/CD_Hit_{cdhit_cutoff}.fasta')
    os.remove(f'{output_dir}/CD_Hit.fasta')
    os.remove(f'{output_dir}/CD_Hit_{cdhit_cutoff}.fasta.clstr')

