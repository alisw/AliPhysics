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

// AliAnalysisTaskTemplate
// 
// author: Redmer A. Bertens, Utrecht University, 2013
// rbertens@cern.ch, rbertens@nikhef.nl, rbertens@uu.nl
//
// template task for flow analysis in the flow package
// can be used as a basis for a general flow analysis task
// see corresponding task in $ALICE_PHYSICS/PWG/FLOW/Base/

// aliroot include
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskTemplate.h"
#include "AliFlowAnalysisTemplate.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliLog.h"

using std::endl;
using std::cout;

ClassImp(AliAnalysisTaskTemplate)

//_____________________________________________________________________________
AliAnalysisTaskTemplate::AliAnalysisTaskTemplate() : AliAnalysisTaskSE(), 
    fFlowTask(NULL),
    fOutputList(NULL),
    fUsePhiWeights(kFALSE),
    fListWeights(NULL),
    fApplyCorrectionForNUA(kFALSE),
    fHarmonic(2) { /* constructor for ROOT IO */ }
//_____________________________________________________________________________
AliAnalysisTaskTemplate::AliAnalysisTaskTemplate(const char *name, Bool_t usePhiWeights) : AliAnalysisTaskSE(name), 
    fFlowTask(NULL),
    fOutputList(NULL),
    fUsePhiWeights(usePhiWeights),
    fListWeights(NULL),
    fApplyCorrectionForNUA(kFALSE),
    fHarmonic(2)
{
    // constructor
    DefineInput(0, AliFlowEventSimple::Class());
    if(usePhiWeights) {
        DefineInput(1, TList::Class()); 
    }
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskTemplate::~AliAnalysisTaskTemplate()
{
    // destructor
}
//_____________________________________________________________________________
void AliAnalysisTaskTemplate::UserCreateOutputObjects() 
{
    // create output objects
    fFlowTask = new AliFlowAnalysisTemplate();
    if (fApplyCorrectionForNUA) fFlowTask->SetApplyCorrectionForNUA(fApplyCorrectionForNUA);
    fFlowTask->SetHarmonic(fHarmonic);
    fFlowTask-> Init();
    // connect the output list of the flow task to the output slot of this task
    if( !(fOutputList = fFlowTask->GetHistList()) ) {
        return;
    } else fOutputList->SetOwner(kTRUE);
    PostData(1,fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskTemplate::UserExec(Option_t *) 
{
    if(dynamic_cast<AliFlowEventSimple*>(GetInputData(0))) {
        fFlowTask->Make(static_cast<AliFlowEventSimple*>(GetInputData(0)));
        PostData(1,fOutputList);
    }
} 
//_____________________________________________________________________________
void AliAnalysisTaskTemplate::Terminate(Option_t *) 
{
    // terminate
}
//_____________________________________________________________________________
