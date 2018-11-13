# PhyCA

----

**PhyCA** is *Julia* package that implements the inference of the *phyletic couplings* presented in the paper *"A multi-scale coevolutionary approach to predict protein-protein
interactions"* by Giancarlo Croce, Thomas Gueudré, Maria Virginia Ruiz Cuevas, Victoria Keidel, Matteo Figliuzzi, Hendrik Szurmant, Martin Weigt.

The **phylogenetic profiling** is a classic bioinformatics technique in which the joint presence or joint absence of two traits across large numbers of species is used to infer a meaningful biological connections

see [Phylogenetic profiling](https://en.wikipedia.org/wiki/Phylogenetic_profiling).

We revisit the classical ideas of **phylogenetic profiling** by introducing the novel concept of **phyletic couplings**, which can be estimated via a global statistical modelling approach ( taking inspiration by *direct coupling analysis* [DCA](https://en.wikipedia.org/wiki/Direct_coupling_analysis) )

The following Figure shows a  schematic representation of the inference of phylogenetic couplings

![figure_method](https://github.com/GiancarloCroce/PhyloDCA/blob/master/figure_1.png)


## Installation
To install the package run *'julia'* in the terminal and type the command

```
    julia>Pkg.clone("https://github.com/GiancarloCroce/PhyCA")
```

---

## Input file 

The first step is to prepare the input file in the right format. Two formats are supported at the moment:

### 1) Genomes in terms of protein families 

In this case you need to construct a file with the composition of genomes in terms of protein families: the first column of the file must be the name of the species, while the columns contains the proteins families included in the genome.

For example the file *"test/phylo_data_ecoli.txt"* has the following structure:

    ACAM1       PF00011 PF00011 PF00023 PF00027 PF00034 PF00011 PF00042 PF00043
    ACCPU       PF00011 PF00015 PF00027 PF00034 PF00037
    ....
    ZINIC       PF00037 PF00109 PF00111 PF00115 PF00116 PF00146


### 2) Phylogenetic profile matrix
It is also possible to give as input file directly the **phylogenetic
profile matrix** : a binary matrix *P* whose entries capture the presence
*(Pij = 1)* or absence *(Pij = 0)* of a domain across genomes, with *i = 1, . . . , M*
(the number of genomes) and *j = 1, . . . , N* (the number of domains). 

Consequently, each domain (the columns of the PPM) is represented by a long binary number with a digit for each genome.

See, as an example, the file *"test/phylo_matrix_ecoli.txt"*.


## Usage

A real documentation is not available yet, but we report here some usage examples to get started.

To run the program type 'julia' in the terminal and load the module:

```
    julia> using PhyCA
```

The software provides two main functions ```phyca(filename_data::String, PhyloDCA.PhylogenticDistance)``` if the input file *"filename_data"* is in the first format  and ```phyca_matrix(filename_matrix::String, PhyloDCA.PhylogenticDistance)``` if the input file *"filename_matrix"* is a Phylogenetic Profile Matrix.


Next is to decide which Phylogenetic Distances we want to use for the analysis (a list of all supported Phylogenetic Distances is in the next section). 

For example if we want to use the "phyletic couplings inferred with mean field DCA", then run

```
    julia> ecoli_results = phyca("phylo_data_ecoli.txt",mfDCA()) 
```



### Output
The output "ecoli_results" is a type PhyCA.PhyloOut  with 6 fields:

* **list_domains**: a list of all proteins families
* **list_species**: a list of all species 
* **PhyloProfile**: the phylogenetic profile matrix
* **PhyloDistance**: the distance matrix between protein domains
* **result_sorted**: a (String, String, Float)  vector containing the candidate candidate domain-domain connections in descending order
* **result_unsort**: a (String, String, Float)  vector containing the candidate candidate domain-domain connections not sorted


### Supported Phylogenetic Distances
For the sake of comparison several Phylogenetic distances have been included in the code:

* Hamming distance [ Hamming() ] 
* Pearson Correlation [ Correlation() ] 
* pValue of the exact Fisher Test  [ pValue() ] 	
* Phylogenetic couplings from the Mean Field DCA [ mfDCA() ]
* Phylogenetic couplings from the pseudo-likelihood DCA [ plmDCA() ]

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
