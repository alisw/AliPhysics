enable_testing()
add_test("func_STAT_AliTreePlayer"
         bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/STAT/test/statTest.sh testAliTreePlayer")
add_test("func_STAT_AliTMinuitToolkitTest"
         bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/STAT/test/statTest.sh testAliTMinutiToolkitTestLinear")
add_test("func_STAT_AliDrawStyleTest"
         bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/STAT/test/statTest.sh testAliDrawStyleTest")
add_test("func_STAT_AliPainterTest"
    bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/STAT/test/statTest.sh testAliPainterTest")
add_test("func_STAT_AliParserTest"
        bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/STAT/test/statTest.sh testAliParserTest")
add_test("func_STAT_TPCVolumeDemoTest"
    bash -c "ALIROOT_SOURCE=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/STAT/test/statTest.sh testTPCVolumeDemoTest")
