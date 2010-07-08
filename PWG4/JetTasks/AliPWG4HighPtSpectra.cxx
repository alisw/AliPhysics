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

//-----------------------------------------------------------------------
// Example of task (running locally, on AliEn and CAF),
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------


#ifndef ALIPWG4HIGHPTSPECTRA_CXX
#define ALIPWG4HIGHPTSPECTRA_CXX

#include "AliPWG4HighPtSpectra.h"

#include "TVector3.h"
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TChain.h"

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"

#include "AliLog.h"

#include "AliStack.h"
#include "TParticle.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliCFContainer.h"

//#include "$ALICE_ROOT/PWG4/JetTasks/AliAnalysisHelperJetTasks.h"

//#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtSpectra)

//__________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra() : AliAnalysisTask("AliPWG4HighPtSpectra", ""), 
  fReadAODData(0),
  fCFManagerPos(0x0),
  fCFManagerNeg(0x0),
  fESD(0),
  fTrackCuts(0),
  fTrackCutsTPConly(0),
  fHistList(0),
  fNEventAll(0),
  fNEventSel(0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra(const Char_t* name) :
  AliAnalysisTask(name,""),
  fReadAODData(0),
  fCFManagerPos(0x0),
  fCFManagerNeg(0x0),
  fESD(0),
  fTrackCuts(),
  fTrackCutsTPConly(0),
  fHistList(0),
  fNEventAll(0),
  fNEventSel(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  AliDebug(2,Form("AliPWG4HighPtSpectra","Calling Constructor"));
  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0,TList::Class());
  // Output slot #1, #2 writes into a AliCFContainer
  DefineOutput(1,AliCFContainer::Class());
  DefineOutput(2,AliCFContainer::Class());
  // Output slot #3 writes into a AliESDtrackCuts
  DefineOutput(3, AliESDtrackCuts::Class());
  DefineOutput(4, AliESDtrackCuts::Class());
}

//________________________________________________________________________
void AliPWG4HighPtSpectra::LocalInit()
{
  //
  // Only called once at beginning
  //
  PostData(3,fTrackCuts);
  PostData(4,fTrackCutsTPConly);
}

//________________________________________________________________________
void AliPWG4HighPtSpectra::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  // Called once
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::ConnectInputData \n"));
  //  cout << "cout >> AliPWG4HighPtSpectra::ConnectInputData" << endl;
  printf(">> AliPWG4HighPtSpectra::ConnectInputData \n");

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    AliDebug(2,Form("ERROR: Could not read chain from input slot 0"));
  } else {
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
    } else {
      fESD = esdH->GetEvent();
    }
  }
  
}
//_________________________________________________
void AliPWG4HighPtSpectra::Exec(Option_t *)
{
  //
  // Main loop function
  //
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::Exec \n"));  

  // All events without selection
  fNEventAll->Fill(0.);

  if (!fESD) {
    AliDebug(2,Form("ERROR: fESD not available"));
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }

  Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!isSelected) { //Select collison candidates
    AliDebug(2,Form(" Trigger Selection: event REJECTED ... "));
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }

  // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
  // This handler can return the current MC event
  
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  //  AliMCEventHandler* eventHandler = (AliMCEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  AliStack* stack = 0x0;
  AliMCEvent* mcEvent = 0x0;
  
  if(eventHandler) {
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      AliDebug(2,Form("ERROR: Could not retrieve MC event"));
      PostData(0,fHistList);
      PostData(1,fCFManagerPos->GetParticleContainer());
      PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
    }
    
    AliDebug(2,Form("MC particles: %d", mcEvent->GetNumberOfTracks()));
    
    stack = mcEvent->Stack();                //Particles Stack
    
    AliDebug(2,Form("MC particles stack: %d", stack->GetNtrack()));
  }
  
  const AliESDVertex *vtx = fESD->GetPrimaryVertex();
  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors()));
  // Need vertex cut
  TString vtxName(vtx->GetName());
  if(vtx->GetNContributors() < 2 || (vtxName.Contains("TPCVertex")) ) {
    // SPD vertex
    vtx = fESD->GetPrimaryVertexSPD();
    if(vtx->GetNContributors()<2) {
      vtx = 0x0;
      // Post output data
      PostData(0,fHistList);
      PostData(1,fCFManagerPos->GetParticleContainer());
      PostData(2,fCFManagerNeg->GetParticleContainer());
      return;
    }
  }
  
  double primVtx[3];
  vtx->GetXYZ(primVtx);
  if(TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2]>10.)){
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }
  
  if(!fESD->GetNumberOfTracks() || fESD->GetNumberOfTracks()<2){ 
    // Post output data
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }
  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d", nTracks));

  if(!fTrackCuts) { 
    // Post output data
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }

  // Selected events for analysis
  fNEventSel->Fill(0.);
  
  
  Double_t containerInputRec[5] ;
  Double_t containerInputTPConly[5];
  Double_t containerInputMC[5];
  //Now go to rec level
  for (Int_t iTrack = 0; iTrack<nTracks; iTrack++) 
    {   
      if(!fESD->GetTrack(iTrack) ) continue;
      AliESDtrack* track = fESD->GetTrack(iTrack);
      if(!(AliExternalTrackParam *)track->GetTPCInnerParam()) continue;
      AliExternalTrackParam *trackTPC = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!track || !trackTPC) continue;

      Float_t dca2D, dcaZ;
      track->GetImpactParameters(dca2D,dcaZ);
      Float_t dca2DTPC, dcaZTPC;
      track->GetImpactParametersTPC(dca2DTPC,dcaZTPC); 
      Float_t chi2PerClusterTPC = -1.;
      Float_t nClustersTPC = track->GetTPCNcls();//track->GetTPCclusters(0);
      if(nClustersTPC>0.) chi2PerClusterTPC = track->GetTPCchi2()/(2.*nClustersTPC-5.);
      Float_t chi2PerClusterTPCIter1 = -1.;
      Float_t nClustersTPCIter1 = track->GetTPCNclsIter1();   
      if(nClustersTPCIter1>0.) chi2PerClusterTPCIter1 = track->GetTPCchi2Iter1()/(2.*nClustersTPCIter1-5.);

      //fill the container
      containerInputRec[0] = track->Pt();
      containerInputRec[1] = track->Phi();
      containerInputRec[2] = track->Eta();
      containerInputRec[3] = dca2D;
      containerInputRec[4] = chi2PerClusterTPC;

      //Store TPC Inner Params for TPConly tracks
      containerInputTPConly[0] = trackTPC->Pt();
      containerInputTPConly[1] = trackTPC->Phi();
      containerInputTPConly[2] = trackTPC->Eta();
      containerInputTPConly[3] = dca2DTPC/10.; //Divide by 10 in order to store in same container. Should be corrected back when looking at output.
      containerInputTPConly[4] = chi2PerClusterTPCIter1;//TPC;

      AliESDtrack* trackTPCESD = fTrackCutsTPConly->GetTPCOnlyTrack(fESD, iTrack);
      if(trackTPCESD) {
	if (fTrackCutsTPConly->AcceptTrack(trackTPCESD)) {
	  if(trackTPC->GetSign()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputTPConly,kStepReconstructedTPCOnly);
	  if(trackTPC->GetSign()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputTPConly,kStepReconstructedTPCOnly);
	}
      }

      if (fTrackCuts->AcceptTrack(track)) {
	if(track->GetSign()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);
	if(track->GetSign()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);

  	
	//Only fill the MC containers if MC information is available
	if(eventHandler) {
	  Int_t label = TMath::Abs(track->GetLabel());
	  TParticle *particle = stack->Particle(label) ;
	  if(!particle) continue;

	  containerInputMC[0] = particle->Pt();      
	  containerInputMC[1] = particle->Phi();      
	  containerInputMC[2] = particle->Eta();  
	  containerInputMC[3] = 0.0;      
	  containerInputMC[4] = 0.0;  

	  //Container with primaries
	  if(stack->IsPhysicalPrimary(label)) {
	    if(particle->GetPDG()->Charge()>0.) {
	      fCFManagerPos->GetParticleContainer()->Fill(containerInputMC,kStepReconstructedMC);
	    }
	    if(particle->GetPDG()->Charge()<0.) {
	      fCFManagerNeg->GetParticleContainer()->Fill(containerInputMC,kStepReconstructedMC);
	    }
	  }

	  //Container with secondaries
	  if (!stack->IsPhysicalPrimary(label) ) {
	    if(particle->GetPDG()->Charge()>0.) {
	      fCFManagerPos->GetParticleContainer()->Fill(containerInputRec,kStepSecondaries);
	    }
	    if(particle->GetPDG()->Charge()<0.) {
	      fCFManagerNeg->GetParticleContainer()->Fill(containerInputRec,kStepSecondaries);
	    }
	  }
	}
	
      }//trackCuts

      delete trackTPCESD;
    }//track loop
  

  //Fill MC containters if particles are findable
  if(eventHandler) {
    for(int iPart = 1; iPart<(mcEvent->GetNumberOfPrimaries()); iPart++)//stack->GetNprimary();
      {
	AliMCParticle *mcPart  = (AliMCParticle*)mcEvent->GetTrack(iPart);
	if(!mcPart) continue;
	
	//fill the container
	containerInputMC[0] = mcPart->Pt();
	containerInputMC[1] = mcPart->Phi();      
	containerInputMC[2] = mcPart->Eta();  
	containerInputMC[3] = 0.0;
	containerInputMC[4] = 0.0;

	if(stack->IsPhysicalPrimary(iPart)) {
	  if(mcPart->Charge()>0. && fCFManagerPos->CheckParticleCuts(kStepMCAcceptance,mcPart)) fCFManagerPos->GetParticleContainer()->Fill(containerInputMC,kStepMCAcceptance);
	  if(mcPart->Charge()<0. && fCFManagerNeg->CheckParticleCuts(kStepMCAcceptance,mcPart)) fCFManagerNeg->GetParticleContainer()->Fill(containerInputMC,kStepMCAcceptance);
	}
      }
  }
  
  PostData(0,fHistList);
  PostData(1,fCFManagerPos->GetParticleContainer());
  PostData(2,fCFManagerNeg->GetParticleContainer());
  
}


//___________________________________________________________________________
void AliPWG4HighPtSpectra::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

//___________________________________________________________________________
void AliPWG4HighPtSpectra::CreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  AliDebug(2,Form("CreateOutputObjects","CreateOutputObjects of task %s", GetName()));

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 

  //slot #1
  OpenFile(0);
  fHistList = new TList();
  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistList->Add(fNEventSel);

  TH1::AddDirectory(oldStatus);   

}

#endif
