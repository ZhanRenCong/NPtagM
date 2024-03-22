This python program is designed to extract natural products tailoring enzymes from microorganism genomes.

Installation：
This program operates depending on the Python environment (Python Version 3.0 or above) and can be run on Linux systems.
The users need to first install the antiSMASH using Bioconda. The installation guide could be found in https://docs.antismash.secondarymetabolites.org/install/.
After the installation of antiSMASH, first activate the environment by command line:
conda activate antismash
Then install AUGUSTUS by:
sudo apt install augustus augustus-data augustus-doc
Finally install CD-Hit by:
sudo apt install cd-hit

Usage:
You could find all parameters by command line:
python run_TEGM.py -h

We provided some test files. You could try command line:
python run_TEGM.py -genome test/fungi_genome -query test/P450 -taxon fungi -out test_out
If you only want to extract terpenoid biosynthetic tailoring enzymes. Try command line:
python run_TEGM.py -genome test/fungi_genome -query test/P450 -taxon fungi -core_gene test/core_gene -out test_out
