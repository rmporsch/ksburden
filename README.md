# ksburden

Most commonly used rare variant tests do not take the position of individual variants into account when testing a gene or a genomic region.
The `ksburden` test aims to make use of hypothesized differences in the distribution of rare variants between cases and controls.


### Prerequisites

Currently `ksburden` only runs on Linux machine machines.
No particular requirements are needed prior installation with the exception of `cmake`.


### Installing

Installation is simple. Just clone this repository then run:

```
cmake .
make
```

The necessary libraries are then downloaded and the program is compiled.

## Usage

### Gene Based Testing

After compilation the program `ksburden` can be found in the `bin` directory.
It is fairly simple and requires only three inputs as shown below

| Flag | Function |
| ---- | -------- |
| --genotypes | variant file in vcf format |
| --variant | list of variants in tab format |
| --ped | pedigree file in plink format |

The variant file is a tab separated text file with 3 columns:

1. Chromosome
2. Position
3. Gene name

In addition a number of other options can be specified. 

| Flag | Function|
| ---- | ------- |
| --out | output file name |
| --iter | number of iterations |
| --ks | perform the KS test in addition to ksburden |
| --burden | perform the burden test in addition to ksburden |
| --CMC | perform the CMC test in addition to ksburden |

### Simulations

In addition to the `ksburden` test I also included a simulation framework.
The framework simulates different mutation scenarios and estimates statistical power for various tests.

For example the code below takes a vcf file `skat.matrix.vcf` and simulates 100 times a case-control sample with a life-time disease risk of 10%. 
It then estimates statistical power for `ksburden`, `burden` and `cmc`.
A full description of all flags can be found below.

```bash
simmat="skat.matrix.vcf"

sim --simmat=$simmat \
	--lifetimerisk=0.10 \
	--powerIter=10 \
	--iter=100 \
	--threads=5 \
	--tests="ksburden,burden,cmc" \
	--minEffect=0.010 \
	--maxEffect=0.010 \
	--effectSteps=0.001 \
	--percentageSteps=0.1 \
	--maxPercentage=1.0 \
	--verbose=9 \
	--path="path/to/somulation/folder/" \
	--out="sometest.txt" 
```

Required Flags:

| Flag | Function |
| ---- | -------- |
| --genotypes | variant file in vcf or plink format (see --simmat) |
| --variant | list of variants in tab format |
| --gene | name of gene to do simulation on |
| --simmat | 'vcf' or 'plink' for files in the corresponding format |


Optional flags:

| Flag | Function |
| ------ | ------- |
| --threads | number of threads to use |
| --subjects | total number of subjects |
| --powerIter | number of iteration to calculate power |
| --verbose | verbose level |
| --iter | number of iterations |
| --lifetimerisk | life time risk |
| --tests | models to perform, should be a string separated by commas |
| --minEffect | minimal effect size |
| --maxEffect | maximal effect size |
| --effectSteps | effect size steps |
| --percentageSteps | percentage of causal variants steps in percentage |
| --maxPercentage | maximal percent of gene covered by causal variants |
| --causalVar | number of causal mutations |
| --outPath | output path of the simulations |

## Built With

* [Armadillo](http://arma.sourceforge.net/) - C++ linear algebra library
* [OpenMP](http://openmp.org/) - Shared Memory multiprocessing programming library
* [Easylogging++](https://github.com/muflihun/easyloggingpp) - C++ logging library
* [LibStatGen](https://github.com/statgen/libStatGen) - C++ library for statistical genetics


## Authors

* **Robert Milan Porsch**  

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

