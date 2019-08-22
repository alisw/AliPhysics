/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnalysisTaskFlowExample
 *
 * example flow class (empty)
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskFlowExample.h"

class AliAnalysisTaskFlowExample;

using namespace std;

ClassImp(AliAnalysisTaskFlowExample)

AliAnalysisTaskFlowExample::AliAnalysisTaskFlowExample() : AliAnalysisTaskSE(),
    fAOD(0), fOutputList(0), fHistPhiEta(0)
{}
//_____________________________________________________________________________
AliAnalysisTaskFlowExample::AliAnalysisTaskFlowExample(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fHistPhiEta(0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskFlowExample::~AliAnalysisTaskFlowExample()
{
    if(fOutputList) {
        delete fOutputList;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowExample::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    fHistPhiEta = new TH2F("fHistPhiEta", "fHistPhiEta; phi; eta", 100, -0.5, 7, 100, -1.5, 1.5);
    fOutputList->Add(fHistPhiEta);
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowExample::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    Int_t iTracks(fAOD->GetNumberOfTracks());
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !track->TestFilterBit(1)) continue;
        fHistPhiEta->Fill(track->Phi(), track->Eta());
    }
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowExample::Terminate(Option_t *)
{}
//_____________________________________________________________________________
