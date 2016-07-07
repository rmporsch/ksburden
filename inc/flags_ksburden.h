#ifndef ksburdenflags_H
#define ksburdenflags_H

#include <gflags/gflags.h>

DECLARE_string(out);
DECLARE_string(vcf);
DECLARE_string(variant);
DECLARE_string(ped);
DECLARE_string(out);

DECLARE_int32(threads);
DECLARE_int32(iter);
DECLARE_int32(verbose);

DECLARE_bool(ks);
DECLARE_bool(burden);
DECLARE_bool(CMC);

void sanity_check(int argc, char* argv[]);

#endif
