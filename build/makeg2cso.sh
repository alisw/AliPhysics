#!/bin/sh
out=$ALICE_ROOT/lib/tgt_$ALICE_TARGET/libg2c_sh
cd /tmp
rm -rf g2c ; mkdir g2c ; cd g2c
ar -x /usr/lib/gcc/darwin/default/libg2c.a
export MACOSX_DEPLOYMENT_TARGET=10.2
rm -f ${out}.so ${out}.dylib
g++ -bundle -flat_namespace -undefined suppress -o ${out}.so `ls *.o | grep -v buggy`
g++ -dynamiclib -flat_namespace -undefined suppress -single_module -o ${out}.dylib `ls *.o | grep -v buggy`

