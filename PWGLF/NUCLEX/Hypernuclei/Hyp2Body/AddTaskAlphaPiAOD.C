#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>

#include "AliAnalysisTaskAlphaPiAOD.h"
#endif

AliAnalysisTaskAlphaPiAOD *AddTaskAlphaPiAOD(bool isMC, TString tskname = "Hyperhydrogen", TString suffix = "") {
    return AliAnalysisTaskAlphaPiAOD::AddTask(isMC, tskname, suffix);
}
