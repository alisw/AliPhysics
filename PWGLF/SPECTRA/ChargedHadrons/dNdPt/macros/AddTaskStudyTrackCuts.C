#if defined(__CLING__)
#include "AliAnalysisTaskSpectraTrackCuts.h"
#endif

 AliAnalysisTaskSpectraTrackCuts* AddTaskStudyTrackCuts(const char* name = "TrackCuts", const char* outfile = 0, int _cutMode=100) {
  return AliAnalysisTaskSpectraTrackCuts::AddTaskSpectraTrackCuts(name, outfile, _cutMode);
}
