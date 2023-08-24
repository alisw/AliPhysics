#if defined(__CLING__)
#include "AliMultDepSpecOnTheFlyAnalysisTask.h"
#endif

AliMultDepSpecOnTheFlyAnalysisTask* AddTask_mkrueger_MultDepSpecOnTheFlyMC(const string& dataSet,
                                                                           TString options = "")
{
  return AliMultDepSpecOnTheFlyAnalysisTask::AddTaskMultDepSpecOnTheFlyMC(dataSet, options);
}
