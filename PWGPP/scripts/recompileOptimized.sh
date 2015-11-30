#
# $ALICE_PHYSICS/../src/PWGPP/scripts/recompileOptimized.sh
#
# shell script recompile specified files with another compiling option
#
# Examples:
 
#  1.) recompile TMatrix and TLinear Fitter  in hghly optimized mode
#         ( touch $ROOTSYS/../src/math/matrix/inc/TMatrix.h ; source $ALICE_PHYSICS/../src/PWGPP/scripts/recompileOptimized.sh;  recompileOptimized  "c++ -g -O5" $ROOTSYS/../src/math/minuit/src/TLinearFitter.cxx  20)
#  2.) recompile AliNDLocalRegression and TStatToolkit in highly optimized mode
#          ( touch $ALICE_ROOT/../src/STAT/AliNDLocalRegression.cxx;  cd $ALICE_ROOT/../build; source $ALICE_PHYSICS/../src/PWGPP/scripts/recompileOptimized.sh;  recompileOptimized  "c++ -g -O5" $ALICE_ROOT/../src/STAT/TStatToolkit.cxx  20 )
#


recompileOptimized(){
    flag=$1 
    inputFile=$2    
    nP=$3
    echo recompile $inputFile with option $flag instead of default 
    touch $inputFile
    make -j $nP  CXX="$flag" VERBOSE=1

}
