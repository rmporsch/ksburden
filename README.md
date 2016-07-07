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
