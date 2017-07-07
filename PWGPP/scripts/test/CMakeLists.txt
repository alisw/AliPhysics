#
# test $Ali_Physics/PWGPP/scripts/utilities.sh
#
enable_testing()
add_test( "utilitiesTest"
        bash -c "(source ${CMAKE_CURRENT_SOURCE_DIR}/utilities.test.sh; test_all|tee test_all.log; exit   `cat test_all.log | grep -c error`)")
# command to be executed can be tested from command line 
# bash "-c" "(source /data/alicesw/sw/SOURCES/AliPhysics/alicesw/0/PWGPP/scripts/test/utilities.test.sh; test_all|tee test_all.log; exit   `cat test_all.log | grep -c error`)"
# or  ctest --verbose