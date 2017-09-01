#
#
enable_testing()
add_test( "func_STEER_Utilities"
        bash -c "(source ${CMAKE_CURRENT_SOURCE_DIR}/utilities.test.sh; test_all|tee test_all.log; exit   `cat test_all.log | grep -c error`)")
