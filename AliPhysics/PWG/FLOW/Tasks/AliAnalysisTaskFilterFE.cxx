/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////
// AliAnalysisTaskFilterFE:
//
// (re)tag RFP and POI in flowEvent in order to 
// reuse it in train
////////////////////////////////////////////////////

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimpleCuts.h"
#include "AliAnalysisTaskFilterFE.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskFilterFE)

//________________________________________________________________________
AliAnalysisTaskFilterFE::AliAnalysisTaskFilterFE() :
  AliAnalysisTaskSE(),
  fCutsRFP(NULL),
  fCutsPOI(NULL),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fFlowEvent(NULL)
{
  // Constructor
  AliDebug(2,"AliAnalysisTaskFilterFE::AliAnalysisTaskFilterFE()");
}

//________________________________________________________________________
AliAnalysisTaskFilterFE::AliAnalysisTaskFilterFE(const char *name, AliFlowTrackSimpleCuts *cutsRFP, AliFlowTrackSimpleCuts *cutsPOI) :
  AliAnalysisTaskSE(name),
  fCutsRFP(cutsRFP),
  fCutsPOI(cutsPOI),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fFlowEvent(NULL)
{
  // Constructor
  AliDebug(2,"AliAnalysisTaskFilterFE::AliAnalysisTaskFilterFE(const char *name, ...)");
  DefineInput( 0, AliFlowEventSimple::Class());
  DefineOutput(1, AliFlowEventSimple::Class());
}

//________________________________________________________________________
AliAnalysisTaskFilterFE::~AliAnalysisTaskFilterFE()
{
  // Destructor
  delete fCutsRFP;
  delete fCutsPOI;
}
//________________________________________________________________________
void AliAnalysisTaskFilterFE::UserCreateOutputObjects()
{
  // Called at every worker node to initialize
  AliDebug(2,"AliAnalysisTaskFilterFE::CreateOutputObjects()");
  PostData(1,fFlowEvent);
}

//________________________________________________________________________
void AliAnalysisTaskFilterFE::UserExec(Option_t *)
{
  // Main loop
  fFlowEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0)); // from TaskSE
  if (!fFlowEvent) return;
  if(fCutsRFP) fFlowEvent->TagRP(fCutsRFP);
  if(fCutsPOI) fFlowEvent->TagPOI(fCutsPOI);
  fFlowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);
  PostData(1,fFlowEvent);
}
