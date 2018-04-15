enable_testing()
add_test("func_TPCbase_AliTPCSpaceCharge3DDriftLine"
         bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/TPC/TPCbase/test/spacechargeTest.sh testAliTPCSpaceCharge3DDriftLine")
