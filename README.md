# Installation

1. git clone the repository and go into the `ksburden` folder
2. use `cmake .`
3. use `make`

The resulting executable programs should be found in `./bin`

# Usage

## Simulation

Required Flags:

| Flag | Function |
| ---- | -------- |
| --genotypes | variant file in vcf format |
| --variant | list of variants in tab format |
| --gene | name of gene to do simualtion on |

Optional flags:

| Flag | Function |
| ------ | ------- |
| --threads | number of threds to use |
| --subjects | total number of subjects |
| --powerIter | number of iteration to calculate power |
| --verbose | verbose level |
| --iter | number of iterations |
| --lifetimerisk | life time risk |
| --causalVar | number of causal mutations |
| --outPath | output path of the simulations |



# Troubleshooting

If you get the following error message:
```
/usr/bin/ld: lib/libStatGen/src/libStatGen/libStatGen.a(VcfFileReader.o): relocation R_X86_64_32 against `.rodata.str1.1' can not be used when making a shared object; recompile with -fPIC
lib/libStatGen/src/libStatGen/libStatGen.a: error adding symbols: Bad value
```

you will need to recompile libstatgen with `-fPIC`.
This can be done by:
```
cd lib/libStatGen/src/libStatGen/
vi Makefiles/Makefile.include
```

here you will need to add the line

```
GFLAGS += -fPIC
```
