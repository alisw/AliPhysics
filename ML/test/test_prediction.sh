#!/bin/bash

DIRPATH="test_extBDT"
mkdir ${DIRPATH}

curl http://personalpages.to.infn.it/~fecchio/test_extBDT/test_xgboost_pt8_12.model -o ${DIRPATH}/test_xgboost_pt8_12.model
curl http://personalpages.to.infn.it/~fecchio/test_extBDT/test_tree_pt8_12.root -o ${DIRPATH}/test_tree_pt8_12.root
curl http://personalpages.to.infn.it/~fecchio/test_extBDT/xgboost_pred.txt -o ${DIRPATH}/xgboost_pred.txt

root -q -b -l ../macros/test_AliEsternalBDT.cc\(\"${DIRPATH}\"\)
