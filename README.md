# Probabilistic Reconstruction of Ancestral Gene Orders Using Intermediate Genomes - PragIG

PragIG is a software for ancestral reconstruction of gene orders based on the concept of Intermediate Genomes (Feijao, 2015; Feijao and Araújo, 2016).
It extends this concept with a probabilistic approach and calculates several Intermediate Genomes and their probability of being the ancestor.
PragIG was implemented in the scope of a master thesis at the Bielefeld University.

# Installation
PragIG is implemented in Python 2.7. Aside from dependencies stated in the requirements.txt, the BioPython package is used by PragIG as well.

```
git clone https://github.com/klamkiew/PragIG.git
cd PragIG
pip install -r requirements.txt
```

BioPython can be found here:
```
http://biopython.org/wiki/Download
```

# Running PragIG

In order to run PragIG the script `pragig.py` is called. All other scripts in the folder are then imported and used by this main script.
```
pragig.py -a ALPHA -e EPSILON -r Repetition <GENOME_FILE> <INPUT_TREE> <OUTPUT_FOLDER>
```
The genome file has to be formatted in the GRIMM format. Other formats are not supported, yet.
Furthermore, only the NEWICK format is supported for the input tree at the moment.
PragIG will create a subfolder in the OUTPUT_FOLDER, that is based on the parameters.

The alpha parameter sets the strictness of a filter used by PragIG that discards sampled genomes that are not likely to be the ancestor.
The epsilon parameter is used to weight adjacencies that are part of an intermediate genome but not conserved in the extant genomes.
The number of sampled genomes for each internal node is determined with the repetition parameter.

# Testing PragIG

There are five datasets in this repository that can be used to test PragIG.
All this datasets have 12 extant genomes with 1000 genes for each genome.
The number of DCJ operations is given in the tree. Each tree has a different diameter, starting by 500 and ending by 2.500.

```
folder = Test_Trees/tree12_n1000_sc.0.5_chr1
src/pragig.py $folder/extant_genomes.txt $folder/sim_tree.nwk $folder
```
This will create a subfolder in the provided test folder with the results of PragIG.

## References

* M. D. Braga and J. Stoye. “The solution space of sorting by DCJ”. In: J. Comput. Biol. 17.9 (Sept. 2010), pp. 1145–1165.
* Pedro Feijão. “Reconstruction of ancestral gene orders using intermediate genomes”. In: BMC Bioinformatics 16.14 (2015).
* Pedro Feijão and Eloi Araujo. “Fast ancestral gene order reconstruction of genomes with unequal gene content”. In: unpublished (2016)
* P. Biller, L. Gueguen, and E. Tannier. “Moments of genome evolution by Double Cut-and-Join”. In: BMC Bioinformatics 16 Suppl 14 (2015),S7.