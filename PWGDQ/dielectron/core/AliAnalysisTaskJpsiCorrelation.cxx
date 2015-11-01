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
//                        Analysis task for J/psi - hadron correlations  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
using namespace std;

#include <TChain.h>
#include <TH1D.h>

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliTriggerAnalysis.h>
#include <AliESDtrackCuts.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliAnalysisTaskMultiDielectron.h"

#include "AliAnalysisTaskJpsiCorrelation.h"
ClassImp(AliAnalysisTaskJpsiCorrelation)

//_________________________________________________________________________________
AliAnalysisTaskJpsiCorrelation::AliAnalysisTaskJpsiCorrelation() :
  AliAnalysisTaskMultiDielectron(),
  fTreesList(),
  fESD(0x0),
  fESDTrackCuts(0x0),
  fIdxDielectron(0),
  fNjpsiPerEvent(0),
  fSign(0),
  fJpsiM(0.0),
  fJpsiPt(0.0),
  fJpsiPhi(0.0),
  fJpsiTheta(0.0),
  fJpsiY(0.0),
  fTrackPt(0.0),
  fTrackPhi(0.0),
  fTrackTheta(0.0),
  fTrackEta(0.0),
  fMultiDieleOutputs(0)
{
  //
  // Constructor
  //
  cout << "AliAnalysisTaskJpsiCorrelation::AliAnalysisTaskJpsiCorrelation()" << endl;
}

//_________________________________________________________________________________
AliAnalysisTaskJpsiCorrelation::AliAnalysisTaskJpsiCorrelation(const char *name) :
  AliAnalysisTaskMultiDielectron(name),
  fTreesList(),
  fESD(0x0),
  fESDTrackCuts(0x0),
  fIdxDielectron(0),
  fNjpsiPerEvent(0),
  fSign(0),
  fJpsiM(0.0),
  fJpsiPt(0.0),
  fJpsiPhi(0.0),
  fJpsiTheta(0.0),
  fJpsiY(0.0),
  fTrackPt(0.0),
  fTrackPhi(0.0),
  fTrackTheta(0.0),
  fTrackEta(0.0),
  fMultiDieleOutputs(0)
{
  //
  // Constructor
  //
  /*
  // Done in the AliAnalysisTaskMultiDielectron(name) constructor
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1D::Class());
  */
  // Add other outputs here
  fMultiDieleOutputs = GetNoutputs();  // number of output slots in the mother analysis task
  DefineOutput(fMultiDieleOutputs, TList::Class());
}


//_________________________________________________________________________________
void AliAnalysisTaskJpsiCorrelation::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //
  cout << "AliAnalysisTaskJpsiCorrelation::UserCreateOutputObjects()" << endl;
  AliAnalysisTaskMultiDielectron::UserCreateOutputObjects();

  if (!fTreesList.IsEmpty()) return; //already initialised
  fTreesList.SetOwner();

  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    TTree* tree=new TTree(die->GetName(), die->GetTitle());
    tree->Branch("idxDielectron",&fIdxDielectron,"fIdxDielectron/I");
    tree->Branch("nJpsiPerEvent",&fNjpsiPerEvent,"fNjpsiPerEvent/I");
    tree->Branch("sign",&fSign,"fSign/I");
    tree->Branch("jpsiM",&fJpsiM,"fJpsiM/D");
    tree->Branch("jpsiPt",&fJpsiPt,"fJpsiPt/D");
    tree->Branch("jpsiPhi",&fJpsiPhi,"fJpsiPhi/D");
    tree->Branch("jpsiTheta",&fJpsiTheta,"fJpsiTheta/D");
    tree->Branch("jpsiY",&fJpsiY,"fJpsiY/D");
    tree->Branch("trackPt",&fTrackPt,"fTrackPt/D");
    tree->Branch("trackPhi",&fTrackPhi,"fTrackPhi/D");
    tree->Branch("trackTheta",&fTrackTheta,"fTrackTheta/D");
    tree->Branch("trackEta",&fTrackEta,"fTrackEta/D");
    fTreesList.Add(tree);
  }
  
  /*
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3, fEventStat);
  */
  PostData(fMultiDieleOutputs, &fTreesList);
}

//_________________________________________________________________________________
void AliAnalysisTaskJpsiCorrelation::UserExec(Option_t *option)
{
  //
  // Main loop. Called for every event
  //
  //cout << "AliAnalysisTaskJpsiCorrelation::UserExec(Option_t *option)" << endl;
  if(fTreesList.IsEmpty()) return;
  AliAnalysisTaskMultiDielectron::UserExec(option);

  //cout << "ndielectrons = " << fListDielectron.GetEntries() << endl;    

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliError("Could not get ESDInputHandler");
    return;
  } 

  fESD = esdH->GetEvent();
  
  if (!fESD) {
    AliError("fESD not available");
    return;
  }
    
  //Process event in all AliDielectron instances
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  Int_t idie=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    //cout << "idie = " << idie << endl;
    die->Process(fESD);
    for(Int_t isign=0; isign<=2; ++isign) {
      //cout << "isign = " << isign << endl;
      Int_t nCandidates=0;
      if(die->GetPairArray(isign))
	nCandidates = die->GetPairArray(isign)->GetEntriesFast();
      if(nCandidates<1) continue;
      cout << "candidates = " << nCandidates << endl;
      
      fIdxDielectron = idie;
      fSign = isign;
      fNjpsiPerEvent = nCandidates;
      
      for(Int_t iCandidate=0; iCandidate<nCandidates; ++iCandidate) {
        TIter nextPair(die->GetPairArray(isign));
	AliDielectronPair* pair = NULL;
	while ((pair = static_cast<AliDielectronPair*>(nextPair()))) {

	  fJpsiM     = pair->M();
	  fJpsiPt    = pair->Pt();
	  fJpsiPhi   = pair->Phi();
	  fJpsiTheta = pair->Theta();
	  fJpsiY     = pair->Y();
	  cout << "jpsiM = " << fJpsiM << endl;
	  cout << "jpsiPt = " << fJpsiPt << endl;
	  cout << "jpsiPhi = " << fJpsiPhi << endl;
	  cout << "jpsiTheta = " << fJpsiTheta << endl;
	  cout << "jpsiY = " << fJpsiY << endl;
	  
	  AliESDtrack* d1 = static_cast<AliESDtrack*>(pair->GetFirstDaughter());
	  AliESDtrack* d2 = static_cast<AliESDtrack*>(pair->GetSecondDaughter());
	  
	  for (Int_t idx = 0; idx < fESD->GetNumberOfTracks(); idx++) {
	    AliESDtrack* esdTrack = fESD->GetTrack(idx);
	    if(!esdTrack) continue;
            if(esdTrack==d1) continue;
            if(esdTrack==d2) continue;
            if(!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
	    
	    fTrackPt = esdTrack->Pt();
	    fTrackPhi = esdTrack->Phi();
	    fTrackTheta = esdTrack->Theta();
	    fTrackEta = esdTrack->Eta();

	    cout << "trackPt = " << fTrackPt << endl;
	    cout << "trackPhi = " << fTrackPhi << endl;
	    cout << "trackTheta = " << fTrackTheta << endl;
	    cout << "trackEta = " << fTrackEta << endl;

	    ((TTree*)fTreesList.At(idie))->Fill();
	  }
	}
      }
    }
    
    ++idie;
  }
  
  PostData(fMultiDieleOutputs, &fTreesList);
}

//_________________________________________________________________________________
void AliAnalysisTaskJpsiCorrelation::FinishTaskOutput()
{
  //
  // Finish function, called after all events 
  //
  cout << "AliAnalysisTaskJpsiCorrelation::FinishTaskOutput()" << endl;
  AliAnalysisTaskMultiDielectron::FinishTaskOutput();
}
