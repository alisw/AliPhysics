/**************************************************************************
* Copyright(c) 2019, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// =================================================================================================
// AliUniFlowCorrTask - ALICE Unified Flow framework : CorrTask
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2019
// =================================================================================================

#ifndef UNIFLOWCORRTASK_CXX
#define UNIFLOWCORRTASK_CXX

#include <vector>
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "AliUniFlowCorrTask.h"

ClassImp(AliUniFlowCorrTask);

// ============================================================================
AliUniFlowCorrTask::AliUniFlowCorrTask() :
  TObject{},
  fbDoRefs{0},
  fbDoPOIs{0},
  fiNumHarm{0},
  fiNumGaps{0},
  fMaxWeightPower{0},
  fMaxHarm{0},
  fsName{},
  fsLabel{},
  fiHarm{},
  fdGaps{}
{}
// ============================================================================
AliUniFlowCorrTask::AliUniFlowCorrTask(Bool_t doRFPs, Bool_t doPOIs, std::vector<Int_t> harms, std::vector<Double_t> gaps) :
  TObject{},
  fbDoRefs{doRFPs},
  fbDoPOIs{doPOIs},
  fiNumHarm{0},
  fiNumGaps{0},
  fMaxWeightPower{0},
  fMaxHarm{0},
  fsName{},
  fsLabel{},
  fiHarm{harms},
  fdGaps{gaps}
{
  // constructor of CorrTask

  fiNumHarm = harms.size();
  fiNumGaps = gaps.size();

  if(fiNumHarm < 2) { return; }

  fMaxWeightPower = harms.size();
  Int_t maxHarm = 0;
  Int_t maxHarmNeg = 0;
  for(Int_t i = 0; i < fMaxWeightPower; i++)
  {
    Int_t harm = fiHarm[i];
    if(harm > 0) maxHarm += harm;
    else maxHarmNeg += harm;
    Int_t abs = TMath::Abs(maxHarmNeg);
    if(abs > maxHarm) maxHarm = abs;
  }
  fMaxHarm = maxHarm;

  // generating name
  TString sName = Form("<<%d>>(%d",fiNumHarm,fiHarm[0]);
  for(Int_t i(1); i < fiNumHarm; ++i) { sName += Form(",%d",fiHarm[i]); }
  sName += ")";

  if(fiNumGaps > 0) {
    sName += Form("_%dsub(%.2g",fiNumGaps+1,fdGaps[0]);
    for(Int_t i(1); i < fiNumGaps; ++i) { sName += Form(",%.2g",fdGaps[i]); }
    sName += ")";
  }

  // generating label
  TString sLabel = Form("<<%d>>_{%d",fiNumHarm,fiHarm[0]);
  for(Int_t i(1); i < fiNumHarm; ++i) { sLabel += Form(",%d",fiHarm[i]); }
  sLabel += "}";

  if(fiNumGaps > 0) {
    sLabel += Form(" %dsub(|#Delta#eta| > %.2g",fiNumGaps+1,fdGaps[0]);
    for(Int_t i(1); i < fiNumGaps; ++i) { sLabel += Form(", |#Delta#eta| > %.2g",fdGaps[i]); }
    sLabel += ")";
  }

  fsName = sName;
  fsLabel = sLabel;
}
// ============================================================================
void AliUniFlowCorrTask::PrintTask() const
{
  printf("AliUniFlowCorrTask::Print():\n");
  printf("# fsName:\t %s\n", fsName.Data());
  printf("# fsLabel:\t %s\n", fsLabel.Data());
  printf("# fbDoRefs:\t %d\n", fbDoRefs);
  printf("# fbDoPOIs:\t %d\n", fbDoPOIs);
  printf("# fiNumHarm:\t %d\n", fiNumHarm);
  printf("# fiHarm:\t { "); for(Int_t i(0); i < fiNumHarm; ++i) { printf("%d ",fiHarm[i]); }  printf("}\n");
  printf("# fiNumGaps:\t %d\n", fiNumGaps);
  printf("# fdGaps:\t { "); for(Int_t i(0); i < fiNumGaps; ++i) { printf("%f ",fdGaps[i]); } printf("}\n");
  printf("############################################\n");
}

#endif
