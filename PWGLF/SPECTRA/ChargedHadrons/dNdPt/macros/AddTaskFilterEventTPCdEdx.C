#if defined(__CLING__)
#include "AliAnalysisTaskFilterEventTPCdEdx.h"
#endif

AliAnalysisTaskFilterEventTPCdEdx* AddTaskFilterEventTPCdEdx(const char* name = "TaskFilterEventTPCdEdx", const char* treefile = "AliESDs_filtered.root", const char* outfile = 0) {
  return AliAnalysisTaskFilterEventTPCdEdx::AddTaskFilterEventTPCdEdx(name, treefile, outfile);
}
