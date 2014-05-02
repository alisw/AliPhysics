#!/bin/sh
fwd=${ALICE_ROOT}/PWGLF/FORWARD/analysis2
fwd=$ANA_SRC

# 1: file
# 2: data type 
# 3: year
# 4: period
# 5: pass 
# 6: run
aliroot -l -b -q ${fwd}/qa/RunQA.C\(\"$1\",\"$2\",$3,\"$4\",\"$5\",$6\)
cp ${fwd}/qa/style.css .
cp ${fwd}/qa/script.css .
cp ${fwd}/qa/fmd_favicon.png . 
cp ${fwd}/qa/fmd_logo.png .

#
# EOF
# 

