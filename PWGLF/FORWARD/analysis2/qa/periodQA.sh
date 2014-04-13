#!/bin/sh
fwd=${ALICE_ROOT}/PWGLF/FORWARD/analysis2
fwd=$ANA_SRC

# 1: file
# 2: data type 
# 3: year
# 4: period
# 5: pass 
aliroot -l -b -q ${fwd}/qa/PeriodQA.C\(\"$1\",\"\",0,\"\",\"\"\)
cp ${fwd}/qa/style.css .
cp ${fwd}/qa/script.js .
cp ${fwd}/qa/fmd_favicon.png . 
cp ${fwd}/qa/fmd_logo.png .


#
# EOF
# 

