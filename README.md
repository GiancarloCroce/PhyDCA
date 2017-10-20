# PhyloDCA

PhyloDCA is a Julia-based tool for predicting functional and/or physical
interactions among proteins domains from genome sequence.
See article ...

## Getting Started

Several methods and datasets have been integrated already and inclusion of
others is under way.

### Supported Phylogenetic Distances
* Hamming distance
* Pearson Correlation
* Mutual Information
* pValue of the exact Fisher Test
* Phylogenetic couplings of the Mean Field DCA
* Phylogenetic couplings of the pseudo-likelihood DCA

```
Give examples
```

### Installing

It requires the installation of:

* [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl). 
```
julia> Pkg.add("NLopt")
```
* [GaussDCA](https://github.com/carlobaldassi/GaussDCA.jl)
julia> Pkg.clone("https://github.com/carlobaldassi/GaussDCA.jl")
```
* [PlmDCA](https://github.com/pagnani/PlmDCA.jl)
julia> Pkg.clone(GaussDCA"https://github.com/pagnani/PlmDCA.jl")
```

## Usage
bla bla bla


## Output
## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```
## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
