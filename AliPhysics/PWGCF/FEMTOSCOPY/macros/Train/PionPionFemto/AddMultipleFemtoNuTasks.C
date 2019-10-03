///
/// \file PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/AddMultipleTasks.C
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#include <TROOT.h>
#include <TClass.h>
#include <iostream>

#include <TString.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>

#define ROOT6_MACROLOAD_WORKAROUND

#include "AliAnalysisTaskFemtoNu.h"

AliAnalysisTask* AddMultipleFemtoNuTasks(TString container,
                                         TString configuration,
                                         TString params,
                                         TString subwagon_suffix)
{
  return AliAnalysisTaskFemtoNu::BuildForMacro(container, configuration, params, subwagon_suffix);
}
