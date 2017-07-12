#
# test $Ali_Physics/PWGPP/scripts/utilities.sh
#
enable_testing()

# AliTreePlayer test
add_test( "testAliTreePlayer"
           bash -c "( export AliRoot_SRC=${PROJECT_SOURCE_DIR}; source $AliRoot_SRC/STEER/Utilities/alilog4bash.sh; source $AliRoot_SRC/STAT/test/statTest.sh; testAliTreePlayer)")
# command to be executed can be tested from command line
#  ctest -R testAliTreePlayer --verbose
# AliTMinuitToolkitTest
add_test( "testAliTMinuitToolkitTest"
           bash -c "( export AliRoot_SRC=${PROJECT_SOURCE_DIR}; source $AliRoot_SRC/STEER/Utilities/alilog4bash.sh; source $AliRoot_SRC/STAT/test/statTest.sh; testAliTMinutiToolkitTestLinear )")
# command to be executed can be tested from command line
#  ctest -R  testAliTMinuitToolkitTest --verbose
