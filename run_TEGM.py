import argparse
import os
import enzyme_extraction
import shutil
import gene_clusters_extraction
import find_core_gene
from multiprocessing import Pool
import itertools

parser = argparse.ArgumentParser(description='python run_TEGM.py')
parser.add_argument('-query', help='A folder which containing all query files. This algorithm accepted two types of files: '
                                   'the protein family of target tailoring enzyme (hmm file, XXX.hmm) or '
                                   'all target tailoring enzyme protein sequences in one sequence file (fasta file, XXX.fasta)', required=True)
parser.add_argument('-genome', help='A folder which containing all genome files used as database.'
                                    'This algorithm accepted fasta file as input', required=True)
parser.add_argument('-o', help='result and intermediate file will store here', required=True)
parser.add_argument('-taxon', help='The genome is from bacteria or Fungi', required=True)
parser.add_argument('-annotation_software', default='glimmerhmm', help='augustus/glimmerhmm. Default: glimmerhmm.'
                                                                       'If the genome is from fungi, two different software could be used for genome annotation.'
                                                                       'glimmerhmm is faster with lower accuracy.'
                                                                       'augustus is slower with higher accuracy.', required=False)
parser.add_argument('-fungi_species', default='aspergillus_oryzae', help='Default: aspergillus_oryzae.'
                                                                       'If the genome is from fungi and the user choose augustus to annotate genomes. '
                                                                        'They could further choose a reference species. Better choose a reference species which is taxonomically related'
                                                                        ' to the genomes. All available reference species could be found in'
                                                                        ' https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md', required=False)
parser.add_argument('-protein',default=None, help='The annotation of genomes will take a long time. If you have the protein annotation file of genomes,'
                                                  'you could also put them into one folder to save the annotation time. This algorithm will find tailoring enzymes in this folder directly.'
                                                  'The file in this folder should have a same name as it in genome folder', required=False)
parser.add_argument('-cdhit_cutoff', default = 0.9, help='Default: 0.9. A lot of potential tailoring enzymes will be extracted from genomes.'
                                                         'To save time for further biosynthetic gene clusters detection. CD-hit is used to de-replicate by sequence identity.'
                                                         'The user could choose the identity cutoff here', required=False)
parser.add_argument('-blast_evalue_TE', default='1E-15', help='Default:1E-5. The evalue used in extraction of homologous tailoring enzymes by blastp', required=False)
parser.add_argument('-hmm_evalue_TE', default='0.01', help='Default:0.01. The evalue used in extraction of homologous tailoring enzymes by hmm', required=False)
parser.add_argument('-blast_evalue_CE', default='1E-15', help='Default:1E-5. The evalue used in extraction of homologous core enzymes by blastp', required=False)
parser.add_argument('-hmm_evalue_CE', default='0.01', help='Default:0.01. The evalue used in extraction of homologous core enzymes by hmm', required=False)
parser.add_argument('-core_gene', default=None, help='Default:all types of natural products biosynthesis core gene. The user could also use their own core gene. '
                                                     'If the user want to use their own core gene, a folder which containing all core gene query files need to be provided.'
                                                     'This algorithm accepted two types of files: '
                                                     'the protein family of target core enzyme (hmm file, XXX.hmm) or '
                                                     'all target core enzyme protein sequences in one sequence file (fasta file, XXX.fasta)', required=False)
parser.add_argument('-core_gene_BGC_length', default=20, help='Default:20 kb. If the user define their own core genes. They could also define the length of BGC. Default length is 20 kb, which'
                                                              'means genes within 20 kb from the core gene will be determined as potential tailoring enzymes. This parameter should be less then tailoring_gene_BGC_length parameter.', required=False)
parser.add_argument('-tailoring_gene_BGC_length', default=50, help='Default:50 kb. If the user want to extract a longer gene cluster, they could change this parameter to a bigger number', required=False)
parser.add_argument('-core', default=4, help='Default: 2. How many CPUs to run in pararllel.', required=False)

args = parser.parse_args()


if not os.path.exists(args.o):
    os.mkdir(args.o)

if not os.path.exists(f'./{args.o}/genome_annotation'):
    os.mkdir(f'./{args.o}/genome_annotation')

if not os.path.exists(f'./{args.o}/blast_hmm_result'):
    os.mkdir(f'./{args.o}/blast_hmm_result')

if not os.path.exists(f'./{args.o}/tailoring_enzyme_sequence'):
    os.mkdir(f'./{args.o}/tailoring_enzyme_sequence')

if not os.path.exists(f'./{args.o}/gene_cluster_blast'):
    os.mkdir(f'./{args.o}/gene_cluster_blast')

if not os.path.exists(f'./{args.o}/gene_cluster'):
    os.mkdir(f'./{args.o}/gene_cluster')

if not os.path.exists(f'./{args.o}/antismash'):
    os.mkdir(f'./{args.o}/antismash')

if not os.path.exists(f'./{args.o}/core_enzyme_sequence'):
    os.mkdir(f'./{args.o}/core_enzyme_sequence')

if not os.path.exists(f'./{args.o}/blast_hmm_result_2'):
    os.mkdir(f'./{args.o}/blast_hmm_result_2')

if not os.path.exists(f'./{args.o}/gene_cluster_2'):
    os.mkdir(f'./{args.o}/gene_cluster_2')

if not os.path.exists(f'./{args.o}/gene_cluster_blast_2'):
    os.mkdir(f'./{args.o}/gene_cluster_blast_2')

if not os.path.exists(f'./{args.o}/BGC_gbk'):
    os.mkdir(f'./{args.o}/BGC_gbk')

if not os.path.exists(f'./{args.o}/BGC_annotation'):
    os.mkdir(f'./{args.o}/BGC_annotation')
    

def final_tailoring_enzyme_extraction(genome, query_folder, output_dir, taxon, annotation_software, fungi_species, protein_dir, blast_evalue, hmm_evalue):
    # annotate genome
    print('start')
    name = os.path.splitext(os.path.basename(genome))[0]
    print(name)
    if os.path.exists(f'{output_dir}/genome_annotation/{name}.faa'):
        print(f'{name} already have annotation file in {output_dir}/genome_annotation directory, skip the genome annotation step')
    else:
        if protein_dir:
            for protein in os.listdir(protein_dir):
                name_1 = os.path.splitext(os.path.basename(protein))[0]
                if name == name_1:
                    shutil.copy(f'{protein_dir}/{protein}', f'{output_dir}/genome_annotation/{name}.faa')
        else:
            assert taxon != 'bacteria' or 'fungi', 'taxon need to be bacteria or fungi'
            if taxon == 'bacteria':
                enzyme_extraction.prodigal_annotation(genome, f'{output_dir}/genome_annotation')
            elif taxon == 'fungi':
                assert annotation_software != 'glimmerhmm' or 'augustus', 'software need to be glimmerhmm or augustus'
                if annotation_software == 'glimmerhmm':
                    enzyme_extraction.glimmerhmm_annotation(genome, f'{output_dir}/genome_annotation')
                elif annotation_software == 'augustus':
                    enzyme_extraction.augustus_annotation(genome, f'{output_dir}/genome_annotation', fungi_species)

    # extract tailoring enzymes from annotation file
    for file in os.listdir(query_folder):
        if os.path.splitext(file)[1] == '.hmm':
            enzyme_extraction.hmm(f'{query_folder}/{file}', f'{output_dir}/genome_annotation/{name}.faa',
                                          f'{output_dir}/blast_hmm_result', hmm_evalue)
        elif os.path.splitext(file)[1] == '.fasta':
            enzyme_extraction.blast(f'{query_folder}/{file}', f'{output_dir}/genome_annotation/{name}.faa',
                                            f'{output_dir}/blast_hmm_result', blast_evalue)
        else:
            print(f'{file} is not provided as right format')

    enzyme_extraction.extract_sequence_from_blast_hmm_result(f'./{output_dir}/blast_hmm_result',
                                                             f'./{output_dir}/tailoring_enzyme_sequence',
                                                             f'{output_dir}/genome_annotation/{name}.faa')
    print('finish') 
                                                             
def final_BGC_extraction(genome, output_dir, tailoring_gene_BGC_length):
    name = os.path.splitext(os.path.basename(genome))[0]
    # extract gene clusters
    gene_clusters_extraction.blast(f'./{output_dir}/tailoring_enzyme_sequence/{name}.fasta', genome, f'./{output_dir}/gene_cluster_blast')
    genome_protein_scaffold = gene_clusters_extraction.blast_result_analysis(f'./{output_dir}/gene_cluster_blast/{name}.out')
    gene_clusters_extraction.extract_gene_cluster(genome_protein_scaffold, genome, f'./{output_dir}/gene_cluster', 1000 * int(tailoring_gene_BGC_length))


def annotate_BGC_extraction(genome, taxon, annotation_software, output_dir, core, fungi_species, core_gene, blast_evalue_CE, hmm_evalue_CE, core_gene_BGC_length):
    name = os.path.splitext(os.path.basename(genome))[0]
    # use antismash/users' file to find core genes in gene clusters
    if core_gene == None:
        if os.path.exists(f'{output_dir}/antismash/{name}'):
            print(f'{name} already have antismash result in {output_dir}/antismash directory, skip the antismash step.')
            find_core_gene.extract_gene_cluster_from_antismash_result(f'{args.o}/gene_cluster/{name}.fasta', f'{args.o}/antismash',  f'{args.o}/BGC_gbk')
        else:
            find_core_gene.find_core_gene_by_antismash(f'{args.o}/gene_cluster/{name}.fasta', taxon, annotation_software, core, f'./{output_dir}/antismash', fungi_species)
            find_core_gene.extract_gene_cluster_from_antismash_result(f'{args.o}/gene_cluster/{name}.fasta', f'{args.o}/antismash', f'{args.o}/BGC_gbk')
    else:
        find_core_gene.find_core_gene_by_users_file(f'{args.o}/gene_cluster/{name}.fasta', taxon, output_dir, annotation_software, fungi_species, core_gene, blast_evalue_CE, hmm_evalue_CE, 1000 * int(core_gene_BGC_length))

arg_list = []

def wrapper_1(args):
    func, genome, query_folder, output_dir, taxon, annotation_software, fungi_species, protein_dir, blast_evalue, hmm_evalue = args
    return func(genome, query_folder, output_dir, taxon, annotation_software, fungi_species, protein_dir, blast_evalue, hmm_evalue)

for genome_ in os.listdir(args.genome):
    arg = (final_tailoring_enzyme_extraction, f'{args.genome}/{genome_}', args.query, args.o,
                                      args.taxon, args.annotation_software, args.fungi_species, args.protein,
                                      args.blast_evalue_TE, args.hmm_evalue_TE)
    arg_list.append(arg)

pool = Pool(processes=int(args.core))
pool.map(wrapper_1, arg_list)
pool.close()
pool.join()

enzyme_extraction.CD_Hit(f'./{args.o}/tailoring_enzyme_sequence', f'./{args.o}', args.cdhit_cutoff)

arg_list_1 = []
arg_list_2 = []

def wrapper(args):
    func, genome, output_dir, tailoring_gene_BGC_length = args
    return func(genome, output_dir, tailoring_gene_BGC_length)

def wrapper_2(args):
    func, genome, taxon, annotation_software, output_dir, core, fungi_species, core_gene, blast_evalue_CE, hmm_evalue_CE, BGC_length = args
    return func(genome, taxon, annotation_software, output_dir, core, fungi_species, core_gene, blast_evalue_CE, hmm_evalue_CE, BGC_length)

for genome in os.listdir(args.genome):
    name = os.path.splitext(os.path.basename(genome))[0]
    if os.path.getsize(f'./{args.o}/tailoring_enzyme_sequence/{name}.fasta') != 0:
        arg = (final_BGC_extraction, f'{args.genome}/{genome}', args.o, args.tailoring_gene_BGC_length)
        arg_list_1.append(arg)
        arg_ = (annotate_BGC_extraction, f'{args.genome}/{genome}', args.taxon, args.annotation_software, args.o,
                args.core, args.fungi_species, args.core_gene, args.blast_evalue_CE, args.hmm_evalue_CE,
                args.core_gene_BGC_length)
        arg_list_2.append(arg_)


pool = Pool(processes=int(args.core))
pool.map(wrapper, arg_list_1)
pool.close()
pool.join()


pool = Pool(processes=int(args.core))
pool.map(wrapper_2, arg_list_2)
pool.close()
pool.join()


# statistics
new_gene_cluster_data = find_core_gene.check_tailoring_enzyme_in_new_BGC(args.query, args.o, args.hmm_evalue_TE, args.blast_evalue_TE)
find_core_gene.statistics(new_gene_cluster_data, args.o)


shutil.rmtree(f'{args.o}/blast_hmm_result')
shutil.rmtree(f'{args.o}/blast_hmm_result_2')
shutil.rmtree(f'{args.o}/core_enzyme_sequence')
shutil.rmtree(f'{args.o}/gene_cluster')
shutil.rmtree(f'{args.o}/gene_cluster_2')
shutil.rmtree(f'{args.o}/gene_cluster_blast')
shutil.rmtree(f'{args.o}/gene_cluster_blast_2')
shutil.rmtree(f'{args.o}/tailoring_enzyme_sequence')
shutil.rmtree(f'{args.o}/BGC_gbk')

find_core_gene.find_tailoring_enzyme_in_MIBIG(args.query, args.o, args.hmm_evalue_TE, args.blast_evalue_TE, args.core_gene, args.hmm_evalue_CE, args.blast_evalue_CE)

