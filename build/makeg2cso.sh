#!/bin/sh
out=$ALICE_ROOT/lib/tgt_$ALICE_TARGET/libg2c_sh
cd /tmp
rm -rf g2c ; mkdir g2c ; cd g2c
ar -x /sw/lib/libg2c.a
export MACOSX_DEPLOYMENT_TARGET=10.3
unset LD_PREBIND
rm -f ${out}.so ${out}.dylib
g++ -bundle -undefined dynamic_lookup -o ${out}.so `ls *.o | grep -v buggy`
g++ -dynamiclib -undefined dynamic_lookup -single_module -o ${out}.dylib `ls *.o | grep -v buggy`

