enable_testing()
add_test( "utilitiesTest"
        bash -c "(source ${CMAKE_CURRENT_SOURCE_DIR}/utilities.test.sh; test_all|tee test_all.log; exit   `cat test_all.log | grep -c error`)")

# /bin/bash "-c" "(source /data/alicesw/sw/SOURCES/AliPhysics/alicesw/0/PWGPP/scripts/test/utilities.test.sh; test_all|tee test_all.log; exit   `cat test_all.log | grep -c error`)"
#  /bin/bash "-c" "(source /data/alicesw/sw/SOURCES/AliPhysics/alicesw/0/PWGPP/scripts/test/utilities.test.sh; value=`test_all|  grep -c error`; exit  $value)"
#bash -c "(source ${CMAKE_CURRENT_SOURCE_DIR}/utilities.test.sh; test_all;)")
# (pwd; ls; source ${CMAKE_CURRENT_SOURCE_DIR}/utilities.test.sh; test_all))
