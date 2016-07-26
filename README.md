# Installation

1. git clone the repository and go into the `ksburden` folder
2. use `cmake .`
3. use `make`

The resulting executable programs should be found in `./bin`

# Usage

## Genome Wide Testing

The program `ksburden` can be found in the `bin` directory.
It is fairly simple and requires only three inputs as shown below

Required Flags:

| Flag | Function |
| ---- | -------- |
| --genotypes | variant file in vcf format |
| --variant | list of variants in tab format |
| --ped | pedigree file in plink format |

The variant file is a tab separated text file with 3 columns:

1. Chromosome
2. Position
3. Gene name

## Simulation

The simulation 

Required Flags:

| Flag | Function |
| ---- | -------- |
| --genotypes | variant file in vcf format |
| --variant | list of variants in tab format |
| --gene | name of gene to do simualtion on |

Further on can also import a genotype matrix in csv format

| Flag | Function |
| ---- | -------- |
| --simmat | genotypematrix in csv format (no headers) |

Optional flags:

| Flag | Function |
| ------ | ------- |
| --threads | number of threds to use |
| --subjects | total number of subjects |
| --powerIter | number of iteration to calculate power |
| --verbose | verbose level |
| --iter | number of iterations |
| --lifetimerisk | life time risk |
| --tests | models to perform, should a string seperated by commas |
| --minEffect | minimal effect size |
| --maxEffect | maximal effect size |
| --effectSteps | effect size steps |
| --percentageSteps | percentage of causal variants steps in percentage |
| --maxPercentage | maximal percent of gene covered by causal variants |
| --causalVar | number of causal mutations |
| --outPath | output path of the simulations |
 
# Examples 

one could execute the following command:

```bash
simmat="skat.matrix"

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

