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
// Efficiency between different steps of the procedure.
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
#include "TProfile.h"
#include "TList.h"
#include "TChain.h"
#include "TKey.h"
#include "TSystem.h"
#include "TFile.h"

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"

#include "AliLog.h"

#include "AliStack.h"
#include "TParticle.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliCFContainer.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliESDVertex.h"
//#include "$ALICE_ROOT/PWG4/JetTasks/AliAnalysisHelperJetTasks.h"

//#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtSpectra)

//__________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra() : AliAnalysisTask("AliPWG4HighPtSpectra", ""), 
  fReadAODData(0),
  fNoPythiaInfo(0),
  fCFManagerPos(0x0),
  fCFManagerNeg(0x0),
  fESD(0x0),
  fAOD(0x0),
  fMC(0x0),
  fStack(0x0),
  fArrayMCAOD(0x0),
  fVtx(0x0),
  fTriggerMask(AliVEvent::kMB),
  fIsPbPb(0),
  fCentClass(10),
  fTrackType(0),
  fTrackCuts(0x0),
  fTrackCutsReject(0x0),
  fFilterMask(0),
  fbSelectHIJING(kFALSE),
  fSigmaConstrainedMax(100.),
  fAvgTrials(1),
  fHistList(0),
  fNEventAll(0),
  fNEventSel(0),
  fNEventReject(0),
  fh1Centrality(0x0),
  fh1Xsec(0),
  fh1Trials(0),
  fh1PtHard(0),
  fh1PtHardTrials(0),
  fPtRelUncertainty1PtPrim(0x0),
  fPtRelUncertainty1PtSec(0x0)
{
  //
  //Default ctor
  //
}

//___________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra(const Char_t* name) :
  AliAnalysisTask(name,""),
  fReadAODData(0),
  fNoPythiaInfo(0),
  fCFManagerPos(0x0),
  fCFManagerNeg(0x0),
  fESD(0x0),
  fAOD(0x0),
  fMC(0x0),
  fStack(0x0),
  fArrayMCAOD(0x0),
  fVtx(0x0),
  fTriggerMask(AliVEvent::kMB),
  fIsPbPb(0),
  fCentClass(10),
  fTrackType(0),
  fTrackCuts(0x0),
  fTrackCutsReject(0x0),
  fFilterMask(0),
  fbSelectHIJING(kFALSE),
  fSigmaConstrainedMax(100.),
  fAvgTrials(1),
  fHistList(0),
  fNEventAll(0),
  fNEventSel(0),
  fNEventReject(0),
  fh1Centrality(0x0),
  fh1Xsec(0),
  fh1Trials(0),
  fh1PtHard(0),
  fh1PtHardTrials(0),
  fPtRelUncertainty1PtPrim(0x0),
  fPtRelUncertainty1PtSec(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  AliDebug(2,Form("AliPWG4HighPtSpectra Calling Constructor"));
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
}

//________________________________________________________________________
void AliPWG4HighPtSpectra::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  // Called once
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::ConnectInputData \n"));

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    AliDebug(2,Form( "ERROR: Could not read chain from input slot 0 \n"));
    return;
  }

  if(fReadAODData) {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
    if (!aodH) {
      AliDebug(2,Form("ERROR: Could not get AODInputHandler"));
      return;
    } else
      fAOD = aodH->GetEvent();
  } else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
    if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
      return;
    } else
      fESD = (AliESDEvent*)esdH->GetEvent();
  }

  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    AliDebug(2,Form( "ERROR: Could not retrieve MC event handler \n"));
  }
  else
    fMC = eventHandler->MCEvent();

}

//________________________________________________________________________
Bool_t AliPWG4HighPtSpectra::SelectEvent()
{
  //
  // Decide if event should be selected for analysis
  //

  // Checks following requirements:
  // - fESD available
  // - trigger info from AliPhysicsSelection
  // - number of reconstructed tracks > 1
  // - primary vertex reconstructed
  // - z-vertex < 10 cm

  Bool_t selectEvent = kTRUE;

  //fESD or fAOD object available?
  if (!fESD && !fAOD) {
    AliDebug(2,Form("ERROR: fInputEvent not available\n"));
    fNEventReject->Fill("noESD/AOD",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Trigger
  UInt_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fTriggerMask != AliVEvent::kAny && !(isSelected&fTriggerMask)) { //Select collison candidates
    AliDebug(2,Form(" Trigger Selection: event REJECTED ... "));
    fNEventReject->Fill("Trigger",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Check if number of reconstructed tracks is larger than 1
  if(!fReadAODData) {
    if(!fESD->GetNumberOfTracks() || fESD->GetNumberOfTracks()<2)  {
      fNEventReject->Fill("NTracks<2",1);
      selectEvent = kFALSE;
      return selectEvent;
    }
  }

  //Check if vertex is reconstructed
  if(fESD)
    fVtx = fESD->GetPrimaryVertexSPD();
  if(fAOD)
    fVtx = fAOD->GetPrimaryVertex();

  if(!fVtx) {
    fNEventReject->Fill("noVTX",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  if(!fReadAODData){
    const AliESDVertex* esdVtx = fESD->GetPrimaryVertexSPD();
    if(!esdVtx->GetStatus()) {
      fNEventReject->Fill("VtxStatus",1);
      selectEvent = kFALSE;
      return selectEvent;
    }
  }

  // Need vertex cut
  //  TString vtxName(fVtx->GetName());
  if(fVtx->GetNContributors()<2) {
    fNEventReject->Fill("NCont<2",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Check if z-vertex < 10 cm
  double primVtx[3];
  fVtx->GetXYZ(primVtx);
  if(TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2]>10.)){
    fNEventReject->Fill("ZVTX>10",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Centrality selection should only be done in case of PbPb
  if(IsPbPb()) {
    Float_t cent = 0.;
    Int_t calccent = 0;
    if(fReadAODData)
      calccent = CalculateCentrality(fAOD);
    else
      calccent = CalculateCentrality(fESD);
    if(fCentClass!=calccent && fCentClass!=10) {
      fNEventReject->Fill("cent",1);
      selectEvent = kFALSE;
      return selectEvent;
    }
    else {
      if(fReadAODData) {
        if(dynamic_cast<AliAODEvent*>(fAOD)->GetCentrality()) {
          cent = dynamic_cast<AliAODEvent*>(fAOD)->GetCentrality()->GetCentralityPercentile("V0M");
        }
      }
      else {
        if(dynamic_cast<AliESDEvent*>(fESD)->GetCentrality()) {
          cent = dynamic_cast<AliESDEvent*>(fESD)->GetCentrality()->GetCentralityPercentile("V0M");
        }
      }

      if(cent>90.) {
	fNEventReject->Fill("cent>90",1);
	selectEvent = kFALSE;
	return selectEvent;
      }
      fh1Centrality->Fill(cent);
    }
  }

  return selectEvent;

}

//________________________________________________________________________
Int_t AliPWG4HighPtSpectra::CalculateCentrality(AliVEvent *event)
{
  Float_t cent = 999;

  if(event){
    if(event->GetCentrality()){
      cent = event->GetCentrality()->GetCentralityPercentile("V0M");
    }
  }

  if(cent<0) return 5;
  if(cent>80)return 4;
  if(cent>50)return 3;
  if(cent>30)return 2;
  if(cent>10)return 1;
  return 0;
}

//_________________________________________________
void AliPWG4HighPtSpectra::Exec(Option_t *)
{
  //
  // Main loop function
  //
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::Exec \n"));  

  AliVEvent* event;
  if(fReadAODData)
    event = fAOD;
  else {
    event = fESD;
    if(!fMC) {
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
        AliDebug(2,Form( "ERROR: Could not retrieve MC event handler \n"));
      }
      else
        fMC = eventHandler->MCEvent();
    }
  }

  // All events without selection
  fNEventAll->Fill(0.);

  if(!SelectEvent()) {
    // Post output data
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }

  if(fReadAODData){
    //Open the MC particles in the case of AODanalysis
    fArrayMCAOD = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  }
  else {
    //In case of ESD event: MCEvent available? if yes: get MC particles stack
    if(fMC) {
      AliDebug(2,Form("MC particles: %d", fMC->GetNumberOfTracks()));
      fStack = fMC->Stack();                //Particles Stack
      AliDebug(2,Form("MC particles stack: %d", fStack->GetNtrack()));
    }
  }

  Int_t nTracks = event->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d", nTracks));

  if((!fTrackCuts && !fReadAODData) || (fFilterMask==0 && fReadAODData)) { //if the TrackCuts are missing in ESD analysis or the FilterBit is missing in AOD analysis
    fNEventReject->Fill("noTrackCuts",1);
    // Post output data
    PostData(0,fHistList);
    PostData(1,fCFManagerPos->GetParticleContainer());
    PostData(2,fCFManagerNeg->GetParticleContainer());
    return;
  }

  // Selected events for analysis
  fNEventSel->Fill(0.);
  
  const Int_t nvar = 4;
  
  Double_t containerInputRec[nvar]       = {0.,0.,0.,0.};
  Double_t containerInputMC[nvar]        = {0.,0.,0.,0.};
  Double_t containerInputRecMC[nvar]     = {0.,0.,0.,0.}; //reconstructed yield as function of MC variable

  //Now go to rec level
  if(fReadAODData) {//AOD analysis
    for (Int_t iTrack = 0; iTrack<nTracks; iTrack++) 
    {
      //Get track for analysis
      AliVTrack *track = 0x0;

      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
      if(!aodtrack) AliFatal("Not a standard AOD");
      if(!aodtrack)
        continue;
      if((!aodtrack->TestFilterBit(fFilterMask)) || 
         ((!fCFManagerPos->CheckParticleCuts(kStepReconstructed,aodtrack)) && (aodtrack->Charge()>0.)) ||
         ((!fCFManagerNeg->CheckParticleCuts(kStepReconstructed,aodtrack)) && (aodtrack->Charge()<0))    )
        continue;
      else {
        track = static_cast<AliVTrack*>(aodtrack);
      }
      //fill the container
      containerInputRec[0] = track->Pt();
      containerInputRec[1] = track->Phi();
      containerInputRec[2] = track->Eta();
      containerInputRec[3] = track->GetTPCNcls();

      if(track->Charge()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);
      if(track->Charge()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);

      if(fArrayMCAOD) {
        Int_t label = TMath::Abs(track->GetLabel());
	if(label>fArrayMCAOD->GetEntries()) {
	  continue;
	}
	AliAODMCParticle *particle = (AliAODMCParticle*) fArrayMCAOD->At(label);
	if(!particle) {
	  continue;
	}
	//Only select particles generated by HIJING if requested
	if(fbSelectHIJING) {
	  if(!IsHIJINGParticle(label)) {
	    continue;
	  }
	}
	containerInputRecMC[0] = particle->Pt();      
	containerInputRecMC[1] = particle->Phi();      
	containerInputRecMC[2] = particle->Eta();  
	containerInputRecMC[3] = track->GetTPCNcls();

	//Container with primaries
	if(particle->IsPhysicalPrimary()) {
	  if(particle->Charge()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputRecMC,kStepReconstructedMC);
	  if(particle->Charge()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputRecMC,kStepReconstructedMC);
	  //Fill pT resolution plots for primaries
	  //fPtRelUncertainty1PtPrim->Fill(containerInputRec[0],containerInputRec[0]*TMath::Sqrt(track->GetSigma1Pt2())); //This has not been implemented in AOD analysis, since they are also produced by the AddTaskPWG4HighPtTrackQA.C macro
	}

	//Container with secondaries
	if (!particle->IsPhysicalPrimary() ) {
	  if(particle->Charge()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputRecMC,kStepSecondaries);
	  if(particle->Charge()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputRecMC,kStepSecondaries);
	  //Fill pT resolution plots for primaries
	  //fPtRelUncertainty1PtSec->Fill(containerInputRec[0],containerInputRec[0]*TMath::Sqrt(track->GetSigma1Pt2())); //This has not been implemented in AOD analysis, since they are also produced by the AddTaskPWG4HighPtTrackQA.C macro
	}
      }
    }//track loop

    //Fill MC containers if particles are findable
    if(fArrayMCAOD) {
      int noPart = fArrayMCAOD->GetEntriesFast();

      for(int iPart = 1; iPart<noPart; iPart++) {
	AliAODMCParticle *mcPart = (AliAODMCParticle*) fArrayMCAOD->At(iPart);
        if(!mcPart) continue;

        //Only select particles generated by HIJING if requested
        if(fbSelectHIJING) {
	  if(!IsHIJINGParticle(iPart))
	    continue;
        }

        Int_t pdg = mcPart->PdgCode();

        // select charged pions, protons, kaons , electrons, muons
        if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) ==
	   321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13){
  
	  //fill the container
	  containerInputMC[0] = mcPart->Pt();
	  containerInputMC[1] = mcPart->Phi();      
	  containerInputMC[2] = mcPart->Eta();  
	  //      AliESDtrack *esdtrack = fESD->GetTrack(mcPart->GetLabel());
	  containerInputMC[3] = 159.;
	  
	  if(mcPart->IsPhysicalPrimary()) {
	    if(mcPart->Charge()>0. && fCFManagerPos->CheckParticleCuts(kStepMCAcceptance,mcPart)) fCFManagerPos->GetParticleContainer()->Fill(containerInputMC,kStepMCAcceptance);
	    if(mcPart->Charge()<0. && fCFManagerNeg->CheckParticleCuts(kStepMCAcceptance,mcPart)) fCFManagerNeg->GetParticleContainer()->Fill(containerInputMC,kStepMCAcceptance);
	  }
        }
      }
    }

  }//end of AOD analysis
  else {//ESD analysis
    for (Int_t iTrack = 0; iTrack<nTracks; iTrack++) 
    {
      //Get track for analysis
      AliESDtrack *track = 0x0;
      AliESDtrack *esdtrack = fESD->GetTrack(iTrack);
      if(!esdtrack)
	continue;
      

      if(fTrackType==4) {
	if (!(fTrackCuts->AcceptTrack(esdtrack)))
	  continue;
      }

      if(fTrackType==1)
	track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
      else if(fTrackType==2 || fTrackType==4) {
	track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(fESD),esdtrack->GetID());
	if(!track)
	  continue;
	
	AliExternalTrackParam exParam;
	Bool_t relate = track->RelateToVertexTPC(fESD->GetPrimaryVertexSPD(),fESD->GetMagneticField(),kVeryBig,&exParam);
	if( !relate ) {
	  if(track) delete track;
	  continue;
	}
	track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
      }
      else if(fTrackType==5 || fTrackType==6) {
	if(fTrackCuts->AcceptTrack(esdtrack)) {
	  continue;
	}
	else {
	  if( !(fTrackCutsReject->AcceptTrack(esdtrack)) && fTrackCuts->AcceptTrack(esdtrack) ) {

	    if(fTrackType==5) {
	      //use TPConly constrained track
	      track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
	      if(!track)
		continue;
	      
	      AliExternalTrackParam exParam;
	      Bool_t relate = track->RelateToVertexTPC(fESD->GetPrimaryVertexSPD(),fESD->GetMagneticField(),kVeryBig,&exParam);
	      if( !relate ) {
		if(track) delete track;
		continue;
	      }
	      track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
	    }
	    else if(fTrackType==6) {
	      //use global constrained track
	      track = new AliESDtrack(*esdtrack);
	      track->Set(esdtrack->GetConstrainedParam()->GetX(),esdtrack->GetConstrainedParam()->GetAlpha(),esdtrack->GetConstrainedParam()->GetParameter(),esdtrack->GetConstrainedParam()->GetCovariance());
	    }
	  }
	}
      }
      else if(fTrackType==7) {
	//use global constrained track
	track = new AliESDtrack(*esdtrack);
      }
      else
	track = esdtrack;
    
      if(!track)
	continue;
      

      if(fTrackType==2 || fTrackType==4 || fTrackType==5) {
	//Cut on chi2 of constrained fit
	if(track->GetConstrainedChi2TPC() > fSigmaConstrainedMax*fSigmaConstrainedMax && fSigmaConstrainedMax>0.) {
	  if(track) delete track;
	  continue;
	}
      }
      if (!(fTrackCuts->AcceptTrack(track)) && fTrackType!=4 && fTrackType!=5 && fTrackType!=6) {
	if(fTrackType==1 || fTrackType==2 || fTrackType==7) {
	  if(track) delete track;
	}
	continue;
      }

      if(fTrackType==7) {
	if(fTrackCutsReject ) {
	  if(fTrackCutsReject->AcceptTrack(track) ) {
	    if(track) delete track;
	    continue;
	  }
	}
	if(esdtrack->GetConstrainedParam()) 
	  track->Set(esdtrack->GetConstrainedParam()->GetX(),esdtrack->GetConstrainedParam()->GetAlpha(),esdtrack->GetConstrainedParam()->GetParameter(),esdtrack->GetConstrainedParam()->GetCovariance());
      }

      if(!track) {
	if(fTrackType==1 || fTrackType==2 || fTrackType==4 || fTrackType==5 || fTrackType==6 || fTrackType==7) {
	  if(track) delete track;
	}
	continue;
      }

      //fill the container
      containerInputRec[0] = track->Pt();
      containerInputRec[1] = track->Phi();
      containerInputRec[2] = track->Eta();
      containerInputRec[3] = track->GetTPCNclsIter1();

      if(track->GetSign()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);
      if(track->GetSign()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);

  	
      //Only fill the MC containers if MC information is available
      if(fMC) {
	Int_t label = TMath::Abs(track->GetLabel());
	if(label>fStack->GetNtrack()) {
	  if(fTrackType==1 || fTrackType==2 || fTrackType==7)
	    delete track;
	  continue;
	}
	//Only select particles generated by HIJING if requested
	if(fbSelectHIJING) {
	  if(!IsHIJINGParticle(label)) {
	    if(fTrackType==1 || fTrackType==2 || fTrackType==7)
	      delete track;
	    continue;
	  }
	}
	TParticle *particle = fStack->Particle(label) ;
	if(!particle) {
	  if(fTrackType==1 || fTrackType==2 || fTrackType==7)
	    delete track;
	  continue;
	}
	containerInputRecMC[0] = particle->Pt();      
	containerInputRecMC[1] = particle->Phi();      
	containerInputRecMC[2] = particle->Eta();  
	containerInputRecMC[3] = track->GetTPCNclsIter1();

	//Container with primaries
	if(fStack->IsPhysicalPrimary(label)) {
	  if(particle->GetPDG()->Charge()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputRecMC,kStepReconstructedMC);
	  if(particle->GetPDG()->Charge()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputRecMC,kStepReconstructedMC);
	  //Fill pT resolution plots for primaries
	  fPtRelUncertainty1PtPrim->Fill(containerInputRec[0],containerInputRec[0]*TMath::Sqrt(track->GetSigma1Pt2()));
	}

	//Container with secondaries
	if (!fStack->IsPhysicalPrimary(label) ) {
	  if(particle->GetPDG()->Charge()>0.) fCFManagerPos->GetParticleContainer()->Fill(containerInputMC,kStepSecondaries);
	  if(particle->GetPDG()->Charge()<0.) fCFManagerNeg->GetParticleContainer()->Fill(containerInputMC,kStepSecondaries);
	  //Fill pT resolution plots for primaries
	  fPtRelUncertainty1PtSec->Fill(containerInputRec[0],containerInputRec[0]*TMath::Sqrt(track->GetSigma1Pt2()));
	}
      }

      if(fTrackType==1  || fTrackType==2 || fTrackType==4 || fTrackType==5 || fTrackType==6 || fTrackType==7) {
	if(track) delete track;
      }
    }//track loop

    //Fill MC containers if particles are findable
    if(fMC) {
      for(int iPart = 1; iPart<(fMC->GetNumberOfPrimaries()); iPart++) {
        AliMCParticle *mcPart  = (AliMCParticle*)fMC->GetTrack(iPart);
        if(!mcPart) continue;
  
        //Only select particles generated by HIJING if requested
        if(fbSelectHIJING) {
	  if(!IsHIJINGParticle(iPart))
	    continue;
        }
  
        Int_t pdg = mcPart->PdgCode();
        
        // select charged pions, protons, kaons , electrons, muons
        if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) ==
	   321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13){
  
	  //fill the container
	  containerInputMC[0] = mcPart->Pt();
	  containerInputMC[1] = mcPart->Phi();      
	  containerInputMC[2] = mcPart->Eta();  
	  //      AliESDtrack *esdtrack = fESD->GetTrack(mcPart->GetLabel());
	  containerInputMC[3] = 159.;
	  
	  if(fStack->IsPhysicalPrimary(iPart)) {
	    if(mcPart->Charge()>0. && fCFManagerPos->CheckParticleCuts(kStepMCAcceptance,mcPart)) fCFManagerPos->GetParticleContainer()->Fill(containerInputMC,kStepMCAcceptance);
	    if(mcPart->Charge()<0. && fCFManagerNeg->CheckParticleCuts(kStepMCAcceptance,mcPart)) fCFManagerNeg->GetParticleContainer()->Fill(containerInputMC,kStepMCAcceptance);
	  }
        }
      }
    }
  }//end of ESD analysis

  PostData(0,fHistList);
  PostData(1,fCFManagerPos->GetParticleContainer());
  PostData(2,fCFManagerNeg->GetParticleContainer());
  
}

//________________________________________________________________________
Bool_t AliPWG4HighPtSpectra::PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials)
{
  //
  // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // Copied from AliAnalysisTaskJetSpectrum2
  //

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if(file.Contains("root_archive.zip#")){
    Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    file.Replace(pos+1,20,"");
  }
  else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  //  Printf("%s",file.Data());
   
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec){
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec){
      // not a severe condition but indicates that we have no information
      return kFALSE;
    }
    else{
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if(!key){
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list){
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree){
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliPWG4HighPtSpectra::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // Copied from AliAnalysisTaskJetSpectrum2
  // 

  if(fNoPythiaInfo)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      //      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    fh1Xsec->Fill("<#sigma>",xsection);
    fh1Trials->Fill("#sum{ntrials}",ftrials);
  } 
  return kTRUE;
}

//________________________________________________________________________
AliGenPythiaEventHeader*  AliPWG4HighPtSpectra::GetPythiaEventHeader(AliMCEvent *mcEvent)
{
  
  if(!mcEvent)return 0;
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  if(!pythiaGenHeader){
    // cocktail ??
    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    
    if (!genCocktailHeader) {
      //      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Unknown header type (not Pythia or Cocktail)");
      //      AliWarning(Form("%s %d: Unknown header type (not Pythia or Cocktail)",(char*)__FILE__,__LINE__));
      return 0;
    }
    TList* headerList = genCocktailHeader->GetHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }
    if(!pythiaGenHeader){
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Pythia event header not found");
      return 0;
    }
  }
  return pythiaGenHeader;

}

//________________________________________________________________________
AliGenHijingEventHeader*  AliPWG4HighPtSpectra::GetHijingEventHeader(AliMCEvent *mcEvent)
{
  
  if(!mcEvent)return 0;
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
  if(!hijingGenHeader){
    // cocktail ??
    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    
    if (!genCocktailHeader) {
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Unknown header type (not Pythia or Cocktail)");
      //      AliWarning(Form("%s %d: Unknown header type (not Pythia or Cocktail)",(char*)__FILE__,__LINE__));
      return 0;
    }
    TList* headerList = genCocktailHeader->GetHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(i));
      if (hijingGenHeader)
        break;
    }
    if(!hijingGenHeader){
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Hijing event header not found");
      return 0;
    }
  }
  return hijingGenHeader;

}

//________________________________________________________________________
AliGenHijingEventHeader*  AliPWG4HighPtSpectra::GetHijingEventHeaderAOD(AliAODEvent *aodEvent)
{
  AliGenHijingEventHeader* hijingGenHeader = 0x0;
  if(aodEvent) {
    AliAODMCHeader* mcHeader = (AliAODMCHeader*) aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    TList* headerList = mcHeader->GetCocktailHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(i));
      if (hijingGenHeader)
        break;
    }
  }
  return hijingGenHeader;
}

//___________________________________________________________________________
void AliPWG4HighPtSpectra::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

//___________________________________________________________________________
void AliPWG4HighPtSpectra::CreateOutputObjects() 
{
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  AliDebug(2,Form("CreateOutputObjects CreateOutputObjects of task %s", GetName()));

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 

  //slot #1
  OpenFile(0);
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistList->Add(fNEventSel);

  fNEventReject = new TH1F("fNEventReject","Reason events are rejectected for analysis",20,0,20);
  //Set labels
  fNEventReject->Fill("noESD/AOD",0);
  fNEventReject->Fill("Trigger",0);
  fNEventReject->Fill("NTracks<2",0);
  fNEventReject->Fill("noVTX",0);
  fNEventReject->Fill("VtxStatus",0);
  fNEventReject->Fill("NCont<2",0);
  fNEventReject->Fill("ZVTX>10",0);
  fNEventReject->Fill("cent",0);
  fNEventReject->Fill("cent>90",0);
  fHistList->Add(fNEventReject);

  fh1Centrality = new TH1F("fh1Centrality","fh1Centrality; Centrality %",100,0,100);
  fHistList->Add(fh1Centrality);

  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fHistList->Add(fh1Xsec);

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fHistList->Add(fh1Trials);

  fh1PtHard       = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fHistList->Add(fh1PtHard);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);
  fHistList->Add(fh1PtHardTrials);

  Int_t fgkNPtBins = 100;
  Float_t kMinPt   = 0.;
  Float_t kMaxPt   = 100.;
  Double_t *binsPt = new Double_t[fgkNPtBins+1];
  for(Int_t i=0; i<=fgkNPtBins; i++) binsPt[i]=(Double_t)kMinPt + (kMaxPt-kMinPt)/fgkNPtBins*(Double_t)i ;

  Int_t fgkNRel1PtUncertaintyBins=50;
  Float_t fgkRel1PtUncertaintyMin = 0.;
  Float_t fgkRel1PtUncertaintyMax = 1.;

  Double_t *binsRel1PtUncertainty=new Double_t[fgkNRel1PtUncertaintyBins+1];
  for(Int_t i=0; i<=fgkNRel1PtUncertaintyBins; i++) binsRel1PtUncertainty[i]=(Double_t)fgkRel1PtUncertaintyMin + (fgkRel1PtUncertaintyMax-fgkRel1PtUncertaintyMin)/fgkNRel1PtUncertaintyBins*(Double_t)i ;

  fPtRelUncertainty1PtPrim = new TH2F("fPtRelUncertainty1PtPrim","fPtRelUncertainty1PtPrim",fgkNPtBins,binsPt,fgkNRel1PtUncertaintyBins,binsRel1PtUncertainty);
  fHistList->Add(fPtRelUncertainty1PtPrim);

  fPtRelUncertainty1PtSec = new TH2F("fPtRelUncertainty1PtSec","fPtRelUncertainty1PtSec",fgkNPtBins,binsPt,fgkNRel1PtUncertaintyBins,binsRel1PtUncertainty);
  fHistList->Add(fPtRelUncertainty1PtSec);

  TH1::AddDirectory(oldStatus);   

  OpenFile(1);
  OpenFile(2);

  PostData(0,fHistList);
  PostData(1,fCFManagerPos->GetParticleContainer());
  PostData(2,fCFManagerNeg->GetParticleContainer());

  if(binsPt)                delete [] binsPt;
  if(binsRel1PtUncertainty) delete [] binsRel1PtUncertainty;

}

//________________________________________________________________________
Bool_t AliPWG4HighPtSpectra::IsHIJINGParticle(Int_t label)
{
  //
  // Return kTRUE in case particle is from HIJING event
  //

  AliGenHijingEventHeader* hijingHeader;
  if(fReadAODData)
    hijingHeader = GetHijingEventHeaderAOD(fAOD);
  else
    hijingHeader = GetHijingEventHeader(fMC);

  Int_t nproduced = hijingHeader->NProduced();

  if(fReadAODData) {
    AliAODMCParticle* mom = (AliAODMCParticle*) fArrayMCAOD->At(label);
    Int_t iMom = label;
    Int_t iParent = mom->GetMother();
    while(iParent!=-1){
      iMom = iParent;
      mom = (AliAODMCParticle*) fArrayMCAOD->At(iMom);
      iParent = mom->GetMother();
    }

    if(iMom<nproduced) {
      return kTRUE;
    } else {
      return kFALSE;
    }
  }
  else { //ESD analysis
    TParticle * mom = fStack->Particle(label);
    Int_t iMom = label;
    Int_t iParent = mom->GetFirstMother();
    while(iParent!=-1){
      iMom = iParent;
      mom = fStack->Particle(iMom);
      iParent = mom->GetFirstMother();
    }
  
    if(iMom<nproduced) {
      return kTRUE;
    } else {
      return kFALSE;
    }
  }
}

#endif
