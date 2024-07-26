This python program is designed to extract natural products tailoring enzymes from microorganism genomes.

Installation：
This program operates depending on the Python environment (Python Version 3.0 or above) and can be run on Linux Ubuntu systems.
1. Install conda:
   The installation guide of conda could be found in https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent
   If you have already installed conda, this step can be skipped.
2. Install antiSMASH through bioconda
   The users need to first install the antiSMASH using Bioconda. Because our algorithm shared same environment with antiSMASH. So this step is also created the right environment for NPtagM. The installation of antiSMASH through Docker or other ways are not OK.
   The installation guide could be found in https://docs.antismash.secondarymetabolites.org/install/.
   2.1 install Bioconda:
       After installing conda, perform a one-time set up of Bioconda with the following commands. This will modify your ~/.condarc file:
       conda config --add channels defaults
       conda config --add channels bioconda
       conda config --add channels conda-forge
       conda config --set channel_priority strict
   2.2 install antiSMASH with the following:
       conda create -n antismash antismash
       conda activate antismash
       download-antismash-databases
       conda deactivate
   2.3 activate environment
       conda activate antismash
3. install AUGUSTUS with the following command:
   sudo apt install augustus augustus-data augustus-doc
4. Install CD-Hit with the following command:
   sudo apt install cd-hit
5. Get the NPtagM package with the following command:
   wget https://github.com/ZhanRenCong/NPtagM
6. Enter the NPtagM directory:
   cd NPtagM

Usage:
You could find all parameters by command line:
python run_TEGM.py -h

We provided some test files. You could try command line:
python run_TEGM.py -genome test/fungi_genome -query test/P450 -taxon fungi -o test_out
If you only want to extract terpenoid biosynthetic tailoring enzymes. Try command line:
python run_TEGM.py -genome test/fungi_genome -query test/P450 -taxon fungi -core_gene test/core_gene -o test_out
