#if defined(__CLING__)
#include "AliAnalysisTaskMKTest.h"
#endif

AliAnalysisTaskMKTest* AddTaskMKTest(const char* name = "TaskMKTest", const char* outfile = 0) {
  return AliAnalysisTaskMKTest::AddTaskMKTest(name, outfile);
}
