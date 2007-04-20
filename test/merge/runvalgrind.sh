#!/bin/sh
rm -f */*.root */*.dat */*.log */fort* */hough */hlt */*~
cd ./backgr
valgrind --error-limit=no --leak-check=full aliroot -b -q sim.C\(1\) 2>&1 | tee sim.log
valgrind --error-limit=no --leak-check=full aliroot -b -q rec.C      2>&1 | tee rec.log
#chmod a-w *.root
#cd ../signal
#valgrind --error-limit=no --leak-check=full aliroot -b -q sim.C\(3\) 2>&1 | tee sim.log
#valgrind --error-limit=no --leak-check=full aliroot -b -q rec.C      2>&1 | tee rec.log


