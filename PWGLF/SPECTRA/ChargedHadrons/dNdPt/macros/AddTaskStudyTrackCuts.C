#if defined(__CLING__)
#include "AliAnalysisTaskTrackCuts.h"
#endif

 AliAnalysisTaskTrackCuts* AddTaskStudyTrackCuts(const char* name = "TrackCuts", const char* outfile = 0, int _cutMode=100) {
  return AliAnalysisTaskTrackCuts::AddTaskTrackCuts(name, outfile, _cutMode);
}
