/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//    Analysis task for creating a reduced data tree                     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <TChain.h>
#include <TH1D.h>
#include <TH2I.h>
#include <TFile.h>
#include <TBits.h>
#include <TRandom.h>
#include <TTimeStamp.h>

#include <AliAnalysisTaskSE.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDHeader.h>
#include <AliAODHeader.h>
#include <AliTriggerAnalysis.h>
#include <AliESDtrack.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliMultEstimator.h>
#include <AliCentrality.h>
#include "TGeoGlobalMagField.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronVarCuts.h"
#include "AliAnalysisTaskComputeLumi.h"

#include <iostream>
#include <vector>
#include <algorithm>
using std::cout;
using std::endl;
using std::flush;

ClassImp(AliAnalysisTaskComputeLumi)

//_________________________________________________________________________________
AliAnalysisTaskComputeLumi::AliAnalysisTaskComputeLumi() :
  AliAnalysisTaskSE(),
  fEventFilter(0x0),
  fHistogramList(0x0),
  fTriggerAnalysis(0x0),
  fTrigAliasVsCentV0M_before(0x0),
  fL0TriggerInputsVsCent_before(0x0),
  fL0TriggerInputsVsTrigAlias_before(0x0),
  fVtxVsCentV0_before(0x0),
  fTrigAliasVsCentV0M_acc(0x0),
  fL0TriggerInputsVsCent_acc(0x0),
  fL0TriggerInputsVsTrigAlias_acc(0x0),
  fVtxVsCentV0_acc(0x0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskComputeLumi::AliAnalysisTaskComputeLumi(const char *name) :
  AliAnalysisTaskSE(name),
  fEventFilter(0x0),
  fHistogramList(0x0),
  fTriggerAnalysis(0x0),
  fTrigAliasVsCentV0M_before(0x0),
  fL0TriggerInputsVsCent_before(0x0),
  fL0TriggerInputsVsTrigAlias_before(0x0),
  fVtxVsCentV0_before(0x0),
  fTrigAliasVsCentV0M_acc(0x0),
  fL0TriggerInputsVsCent_acc(0x0),
  fL0TriggerInputsVsTrigAlias_acc(0x0),
  fVtxVsCentV0_acc(0x0)
{
  //
  // Constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());   // list of output histograms
}


//_________________________________________________________________________________
void AliAnalysisTaskComputeLumi::UserCreateOutputObjects()
{
  //
  // create the output TList
  //
  
  // TList for histograms
  fHistogramList = new TList();
  fHistogramList->SetOwner();

  const Char_t* aliasNames[32] = {
          "MB/INT1", "INT7", "MUON", "HighMult/HighMultSPD", 
          "EMC1", "CINT5/INT5", "CMUS5/MUSPB/INT7inMUON", "MuonSingleHighPt7/MUSH7/MUSHPB", 
          "MuonLikeLowPt7/MUL7/MuonLikePB", "MuonUnlikeLowPt7/MUU7/MuonUnlikePB", "EMC7/EMC8", "MUS7/MuonSingleLowPt7", 
          "PHI1", "PHI7/PHI8/PHOSPb", "EMCEJE", "EMCEGA", 
          "Central/HighMultV0", "SemiCentral", "DG/DG5", "ZED", 
          "SPI7/SPI", "INT8", "MuonSingleLowPt8", "MuonSingleHighPt8",
          "MuonLikeLowPt8", "MuonUnlikeLowPt8", "MuonUnlikeLowPt0/INT6", "UserDefined", 
          "TRD", "MuonCalo/CaloOnly", "FastOnly", "N/A"
  };
  
  fTrigAliasVsCentV0M_before = new TH2I("TrigAliasVsCent_before","Trig alias vs cent V0M, before;centV0M;", 50, 0.0, 100., 33, -1.0, 32.0);
  fTrigAliasVsCentV0M_acc = new TH2I("TrigAliasVsCent_acc","Trig alias vs cent V0M, accepted;centV0M;", 50, 0.0, 100., 33, -1.0, 32.0);
  for (Int_t i=0; i<32; i++) {
    fTrigAliasVsCentV0M_before->GetYaxis()->SetBinLabel(i+2,aliasNames[i]);
    fTrigAliasVsCentV0M_acc->GetYaxis()->SetBinLabel(i+2,aliasNames[i]);
  }
  fTrigAliasVsCentV0M_before->GetYaxis()->SetBinLabel(1,"Total");
  fTrigAliasVsCentV0M_acc->GetYaxis()->SetBinLabel(1,"Total");
  fHistogramList->Add(fTrigAliasVsCentV0M_before);
  fHistogramList->Add(fTrigAliasVsCentV0M_acc);
  
  fL0TriggerInputsVsCent_before = new TH2I("L0InputVsCent_before","L0 input vs cent V0M, before;centV0M;", 50, 0.0, 100., 33, -1.0, 32.0);
  fL0TriggerInputsVsCent_acc = new TH2I("L0InputVsCent_acc","L0 input vs cent V0M, accepted;centV0M;", 50, 0.0, 100., 33, -1.0, 32.0);
  for (Int_t i=0; i<32; i++) {
    fL0TriggerInputsVsCent_before->GetYaxis()->SetBinLabel(i+2,Form("L0input %d", i+1));
    fL0TriggerInputsVsCent_acc->GetYaxis()->SetBinLabel(i+2,Form("L0input %d", i+1));
  }
  fL0TriggerInputsVsCent_before->GetYaxis()->SetBinLabel(1,"Total");
  fL0TriggerInputsVsCent_acc->GetYaxis()->SetBinLabel(1,"Total");
  fHistogramList->Add(fL0TriggerInputsVsCent_before);
  fHistogramList->Add(fL0TriggerInputsVsCent_acc); 
  
  fL0TriggerInputsVsTrigAlias_before = new TH2I("L0InputVsTrigAlias_before","L0 input vs trigger alias, before", 33, -1.0, 32., 33, -1.0, 32.0);
  fL0TriggerInputsVsTrigAlias_acc = new TH2I("L0InputVsTrigAlias_accepted","L0 input vs trigger alias, accepted", 33, -1.0, 32., 33, -1.0, 32.0);
  for (Int_t i=0; i<32; i++) {
    fL0TriggerInputsVsTrigAlias_before->GetXaxis()->SetBinLabel(i+2,aliasNames[i]);
    fL0TriggerInputsVsTrigAlias_before->GetYaxis()->SetBinLabel(i+2,Form("L0input %d", i+1));
    fL0TriggerInputsVsTrigAlias_acc->GetXaxis()->SetBinLabel(i+2,aliasNames[i]);
    fL0TriggerInputsVsTrigAlias_acc->GetYaxis()->SetBinLabel(i+2,Form("L0input %d", i+1));
  }
  fL0TriggerInputsVsTrigAlias_before->GetXaxis()->SetBinLabel(1,"Total");
  fL0TriggerInputsVsTrigAlias_acc->GetXaxis()->SetBinLabel(1,"Total");
  fL0TriggerInputsVsTrigAlias_before->GetYaxis()->SetBinLabel(1,"Total");
  fL0TriggerInputsVsTrigAlias_acc->GetYaxis()->SetBinLabel(1,"Total");
  fHistogramList->Add(fL0TriggerInputsVsTrigAlias_before);
  fHistogramList->Add(fL0TriggerInputsVsTrigAlias_acc);
  
  fVtxVsCentV0_before = new TH2I("VtxVsCentV0_before", "Vtx Z vs CentV0; centV0; vtx. Z (cm)", 50, 0.0, 100.0, 60, -15.0, 15.0);
  fVtxVsCentV0_acc = new TH2I("VtxVsCentV0_acc", "Vtx Z vs CentV0; centV0; vtx. Z (cm)", 50, 0.0, 100.0, 60, -15.0, 15.0);
  fHistogramList->Add(fVtxVsCentV0_before);
  fHistogramList->Add(fVtxVsCentV0_acc);
  
  fTriggerAnalysis = new AliTriggerAnalysis();
  
  PostData(1, fHistogramList);
}

//_________________________________________________________________________________
void AliAnalysisTaskComputeLumi::UserExec(Option_t *option)
{
  //
  // Main loop. Called for every event
  //
  //cout << "***************************** AliAnalysisTaskReducedTreeMaker::UserExec()  IN" << endl;  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) return;

  // Was event selected ?
  UInt_t isPhysSel = AliVEvent::kAny;
    
  if((isESD && inputHandler->GetEventSelection()) || isAOD){
    isPhysSel = inputHandler->IsEventSelected();
  }
    
  // Get centrality object
  const Int_t nCentEstimators = 3;
  const Char_t* estimatorNames[nCentEstimators] = {"V0M", "ZNA", "CL1"};
  Double_t percentileEstimators[nCentEstimators] = {999., 999., 999.};
  AliVEvent* event = InputEvent();
  AliMultSelection* multSelection = (AliMultSelection*)event->FindListObject("MultSelection"); // new centrality framework
  if(!multSelection) 
    AliInfo("No centrality object found");
  for(Int_t i=0; i<nCentEstimators; ++i) 
    percentileEstimators[i] = multSelection->GetMultiplicityPercentile(Form("%s",estimatorNames[i]));
  
  const AliVVertex *primVtx = event->GetPrimaryVertex();
  Float_t vtxZ = -999.0;
  if (primVtx) { 
    vtxZ = primVtx->GetZ();
  }
  UInt_t inputsL0 = ((AliESDEvent*)event)->GetHeader()->GetL0TriggerInputs();
  
  fVtxVsCentV0_before->Fill(percentileEstimators[0], vtxZ);
  
  fTrigAliasVsCentV0M_before->Fill(percentileEstimators[0], -0.5);
  fL0TriggerInputsVsTrigAlias_before->Fill(-0.5, -0.5);
  fL0TriggerInputsVsCent_before->Fill(percentileEstimators[0], -0.5);
  for (Int_t i=0; i<32; i++) {
    if (isPhysSel & (UInt_t(1)<<i)) {
      fTrigAliasVsCentV0M_before->Fill(percentileEstimators[0], 0.5+i);        
      for (Int_t j=0; j<32; j++) {
        if (inputsL0 & (UInt_t(1)<<j)) {
          fL0TriggerInputsVsTrigAlias_before->Fill(0.5+i, 0.5+j);
        }      
      }
    }
    if (inputsL0 & (UInt_t(1)<<i)) {
      fL0TriggerInputsVsCent_before->Fill(percentileEstimators[0], 0.5+i);        
    }
  }
  for (Int_t i=0; i<32; i++) {
    if (isPhysSel & (UInt_t(1)<<i)) {
      fL0TriggerInputsVsTrigAlias_before->Fill(0.5+i, -0.5);
    }
    if (inputsL0 & (UInt_t(1)<<i)) {
      fL0TriggerInputsVsTrigAlias_before->Fill(-0.5, 0.5+i);      
    }
  }
  
  if(fEventFilter && !fEventFilter->IsSelected(InputEvent())) {
    return;        
  }
  
  fVtxVsCentV0_acc->Fill(percentileEstimators[0], vtxZ);
  
  fTrigAliasVsCentV0M_acc->Fill(percentileEstimators[0], -0.5);
  fL0TriggerInputsVsTrigAlias_acc->Fill(-0.5, -0.5);
  fL0TriggerInputsVsCent_acc->Fill(percentileEstimators[0], -0.5);
  for (Int_t i=0; i<32; i++) {
    if (isPhysSel & (UInt_t(1)<<i)) {
      fTrigAliasVsCentV0M_acc->Fill(percentileEstimators[0], 0.5+i);        
      for (Int_t j=0; j<32; j++) {
        if (inputsL0 & (UInt_t(1)<<j)) {
          fL0TriggerInputsVsTrigAlias_acc->Fill(0.5+i, 0.5+j);
        }      
      }
    }
    if (inputsL0 & (UInt_t(1)<<i)) {
      fL0TriggerInputsVsCent_acc->Fill(percentileEstimators[0], 0.5+i);        
    }
  }
  for (Int_t i=0; i<32; i++) {
    if (isPhysSel & (UInt_t(1)<<i)) {
      fL0TriggerInputsVsTrigAlias_acc->Fill(0.5+i, -0.5);
    }
    if (inputsL0 & (UInt_t(1)<<i)) {
      fL0TriggerInputsVsTrigAlias_acc->Fill(-0.5, 0.5+i);      
    }
  }
  
  PostData(1, fHistogramList);
}


//_________________________________________________________________________________
void AliAnalysisTaskComputeLumi::FinishTaskOutput()
{
  //
  // Finish Task 
  //
  PostData(1, fHistogramList);
}
