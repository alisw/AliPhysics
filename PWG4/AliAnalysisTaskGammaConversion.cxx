/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// root
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <Riostream.h>

// analysis
//#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGammaConversion.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliStack.h"
#include "AliLog.h"
//#include "TLorentzVector.h"
#include "AliKFVertex.h"

ClassImp(AliAnalysisTaskGammaConversion)


AliAnalysisTaskGammaConversion::AliAnalysisTaskGammaConversion():
  AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fStack(NULL),
  fOutputContainer(NULL),
  fHistograms(NULL),
  fDoMCTruth(kFALSE),
  fMCAllGammas(),
  fMCPi0s(),
  fMCEtas(),
  fMCGammaChi_c(),
  fKFReconstructedGammas(),
  fElectronMass(-1),
  fGammaMass(-1),
  fPi0Mass(-1),
  fEtaMass(-1),
  fGammaWidth(-1),
  fPi0Width(-1),
  fEtaWidth(-1),
  fCalculateBackground(kFALSE)
{
  // Default constructor
  // Common I/O in slot 0
  DefineInput (0, TChain::Class());
  DefineOutput(0, TTree::Class());

  // Your private output
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaConversion::AliAnalysisTaskGammaConversion(const char* name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fStack(NULL),
  fOutputContainer(0x0),
  fHistograms(NULL),
  fDoMCTruth(kFALSE),
  fMCAllGammas(),
  fMCPi0s(),
  fMCEtas(),
  fMCGammaChi_c(),
  fKFReconstructedGammas(),
  fElectronMass(-1),
  fGammaMass(-1),
  fPi0Mass(-1),
  fEtaMass(-1),
  fGammaWidth(-1),
  fPi0Width(-1),
  fEtaWidth(-1),
  fCalculateBackground(kFALSE)
{
  // Common I/O in slot 0
  DefineInput (0, TChain::Class());
  DefineOutput(0, TTree::Class());
  
  // Your private output
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaConversion::~AliAnalysisTaskGammaConversion() 
{
  // Remove all pointers
 
  if(fOutputContainer){
    fOutputContainer->Clear() ; 
    delete fOutputContainer ;
  }
  if(fHistograms){
    delete fHistograms;
  }
  if(fV0Reader){
    delete fV0Reader;
  }
}


void AliAnalysisTaskGammaConversion::Init()
{
  // Initialization
}


void AliAnalysisTaskGammaConversion::Exec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  
  ConnectInputData("");
  
  //clear vectors
  fMCAllGammas.clear();
  fMCPi0s.clear();
  fMCEtas.clear();
  fMCGammaChi_c.clear();
  /*
  for(UInt_t ikf=0;ikf<fKFReconstructedGammas.size();ikf++){
    delete fKFReconstructedGammas[ikf];
    fKFReconstructedGammas[ikf]=NULL;
  }
  */
  fKFReconstructedGammas.clear();

  //Clear the data in the v0Reader
  fV0Reader->UpdateEventByEventData();
  // Process the v0 information
  if(fDoMCTruth){
    ProcessMCData();
  }

  // Process the v0 information
  ProcessV0s();

  if(fCalculateBackground){//calculate background if flag is set
    CalculateBackground();
  }

  // Process reconstructed gammas
  ProcessGammasForNeutralMesonAnalysis();

  PostData(1, fOutputContainer);
  
}

void AliAnalysisTaskGammaConversion::ConnectInputData(Option_t */*option*/){

  if(fV0Reader == NULL){
    // Write warning here cuts and so on are default if this ever happens
  }
  fV0Reader->Initialize();
}

void AliAnalysisTaskGammaConversion::ProcessMCData(){
  
  fStack = fV0Reader->GetMCStack();

  for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
    TParticle* particle = (TParticle *)fStack->Particle(iTracks);

    if (!particle) {
      //print warning here
      continue;
    }
    
    if(particle->Pt()<fV0Reader->GetPtCut()){
      continue;
    }
    
    if(TMath::Abs(particle->Eta())> fV0Reader->GetEtaCut()){
      continue;
    }

    if(particle->R()>fV0Reader->GetMaxRCut()){ // cuts on distance from collision point
      continue;
    }

    Double_t tmpPhi=particle->Phi();
    if(particle->Phi()> TMath::Pi()){
      tmpPhi = particle->Phi()-(2*TMath::Pi());
    }

    
    //process the gammas
    if (particle->GetPdgCode()== 22){
      fMCAllGammas.push_back(particle);
      if(particle->GetMother(0)>-1){ //Means we have a mother
	if( fStack->Particle(particle->GetMother(0))->GetPdgCode() != 22 ){//Checks for a non gamma mother.
	  if(fHistograms->fMC_Gamma_Energy){fHistograms->fMC_Gamma_Energy->Fill(particle->Energy());}
	  if(fHistograms->fMC_Gamma_Pt){fHistograms->fMC_Gamma_Pt->Fill(particle->Pt());}
	  if(fHistograms->fMC_Gamma_Eta){fHistograms->fMC_Gamma_Eta->Fill(particle->Eta());}
	  
	  /*	  Double_t tmpPhi=particle->Phi();
	  if(particle->Phi()> TMath::Pi()){
	    tmpPhi = particle->Phi()-(2*TMath::Pi());
	    }*/
	  if(fHistograms->fMC_Gamma_Phi){fHistograms->fMC_Gamma_Phi->Fill(tmpPhi);}

	  //adding the conversion points from all gammas with e+e- daughters
	  if(particle->GetNDaughters() == 2){
	    TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
	    TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
	    
	    if(daughter0->R()>fV0Reader->GetMaxRCut() || daughter1->R()>fV0Reader->GetMaxRCut()){
	      continue;
	    }
	    
	    if(daughter0->GetPdgCode() == -11 && daughter1->GetPdgCode() == 11 ||
	       daughter0->GetPdgCode() == 11 && daughter1->GetPdgCode() == -11){

	      // begin Mapping 
	      Int_t rBin    = fHistograms->GetRBin(daughter0->R());
	      Int_t phiBin  = fHistograms->GetPhiBin(daughter0->Phi());
	      
	      if(fHistograms->fMC_Mapping[phiBin][rBin] != NULL){fHistograms->fMC_Mapping[phiBin][rBin]->Fill(daughter0->Vz(), particle->Eta());}
	      if(fHistograms->fMC_Mapping_Phi[phiBin] != NULL){fHistograms->fMC_Mapping_Phi[phiBin]->Fill(daughter0->Vz(), particle->Eta());}
	      if(fHistograms->fMC_Mapping_R[rBin] != NULL){fHistograms->fMC_Mapping_R[rBin]->Fill(daughter0->Vz(), particle->Eta());}
	      //end mapping


	      if(fHistograms->fMC_EP_R){fHistograms->fMC_EP_R->Fill(daughter0->R());}
	      if(fHistograms->fMC_EP_Z_R){fHistograms->fMC_EP_Z_R->Fill(daughter0->Vz(),daughter0->R());}
	      if(fHistograms->fMC_EP_X_Y){fHistograms->fMC_EP_X_Y->Fill(daughter0->Vx(),daughter0->Vy());}
	      if(fHistograms->fMC_EP_OpeningAngle){fHistograms->fMC_EP_OpeningAngle->Fill(GetMCOpeningAngle(daughter0, daughter1));}

	    }
	  }
	}
	if( fStack->Particle(particle->GetMother(0))->GetPdgCode()==10441 ||//chi_c0 
	    fStack->Particle(particle->GetMother(0))->GetPdgCode()==20443 ||//psi_2S
	    fStack->Particle(particle->GetMother(0))->GetPdgCode()==445  //chi_c2
	){ 
	  fMCGammaChi_c.push_back(particle);
         }
      }
      else{//means we have a primary particle
	if(fHistograms->fMC_DirectGamma_Energy){fHistograms->fMC_DirectGamma_Energy->Fill(particle->Energy());}
	if(fHistograms->fMC_DirectGamma_Pt){fHistograms->fMC_DirectGamma_Pt->Fill(particle->Pt());}
	if(fHistograms->fMC_DirectGamma_Eta){fHistograms->fMC_DirectGamma_Eta->Fill(particle->Eta());}
	/*
	Double_t tmpPhi=particle->Phi();
	if(particle->Phi()> TMath::Pi()){
	  tmpPhi = particle->Phi()-(2*TMath::Pi());
	}
	*/
	if(fHistograms->fMC_DirectGamma_Phi){fHistograms->fMC_DirectGamma_Phi->Fill(tmpPhi);}

	//adding the conversion points from all gammas with e+e- daughters
	if(particle->GetNDaughters() == 2){
	  TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
	  TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
	  if(daughter0->GetPdgCode() == -11 && daughter1->GetPdgCode() == 11 ||
	     daughter0->GetPdgCode() == 11 && daughter1->GetPdgCode() == -11){
	    
	    if(fHistograms->fMC_EP_R){fHistograms->fMC_EP_R->Fill(daughter0->R());}
	    if(fHistograms->fMC_EP_Z_R){fHistograms->fMC_EP_Z_R->Fill(daughter0->Vz(),daughter0->R());}
	    if(fHistograms->fMC_EP_X_Y){fHistograms->fMC_EP_X_Y->Fill(daughter0->Vx(),daughter0->Vy());}
	    if(fHistograms->fMC_EP_OpeningAngle){fHistograms->fMC_EP_OpeningAngle->Fill(GetMCOpeningAngle(daughter0,daughter1));}

	  }
	}

      }
    }
    else if (TMath::Abs(particle->GetPdgCode())== 11){ // Means we have an electron or a positron
      if(particle->GetMother(0)>-1){ // means we have a mother
	if( fStack->Particle(particle->GetMother(0))->GetPdgCode()==22 ){ // Means we have a gamma mother
	  if(particle->GetPdgCode() == 11){//electron 
	    if(fHistograms->fMC_E_Energy){fHistograms->fMC_E_Energy->Fill(particle->Energy());}
	    if(fHistograms->fMC_E_Pt){fHistograms->fMC_E_Pt->Fill(particle->Pt());}
	    if(fHistograms->fMC_E_Eta){fHistograms->fMC_E_Eta->Fill(particle->Eta());}
	    /*	    Double_t tmpPhi=particle->Phi();
	    if(particle->Phi()> TMath::Pi()){
	      tmpPhi = particle->Phi()-(2*TMath::Pi());
	      }*/
	    if(fHistograms->fMC_E_Phi){fHistograms->fMC_E_Phi->Fill(tmpPhi);}
	  }
	  if(particle->GetPdgCode() == -11){//positron 
	    if(fHistograms->fMC_P_Energy){fHistograms->fMC_P_Energy->Fill(particle->Energy());}
	    if(fHistograms->fMC_P_Pt){fHistograms->fMC_P_Pt->Fill(particle->Pt());}
	    if(fHistograms->fMC_P_Eta){fHistograms->fMC_P_Eta->Fill(particle->Eta());}
	    /*
	    Double_t tmpPhi=particle->Phi();
	    if(particle->Phi()> TMath::Pi()){
	      tmpPhi = particle->Phi()-(2*TMath::Pi());
	    }
	    */
	    if(fHistograms->fMC_P_Phi){fHistograms->fMC_P_Phi->Fill(tmpPhi);}
	  }
	}
      }
    }
    else if(particle->GetNDaughters() == 2){

      TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
      TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
      if(daughter0->GetPdgCode() == 22 && daughter1->GetPdgCode() == 22){//check for gamma gamma daughters
	
	if(particle->GetPdgCode()==111){//Pi0
	  if( iTracks >= fStack->GetNprimary()){
	    
	    if(fHistograms->fMC_Pi0Secondaries_Eta){fHistograms->fMC_Pi0Secondaries_Eta->Fill(particle->Eta());}
	    /*
	      Double_t tmpPhi=particle->Phi();
	      if(particle->Phi()> TMath::Pi()){
	      tmpPhi = particle->Phi()-(2*TMath::Pi());
	      }
	    */
	    if(fHistograms->fMC_Pi0Secondaries_Phi){fHistograms->fMC_Pi0Secondaries_Phi->Fill(tmpPhi);}
	    if(fHistograms->fMC_Pi0Secondaries_Pt){fHistograms->fMC_Pi0Secondaries_Pt->Fill(particle->Pt());}
	    if(fHistograms->fMC_Pi0Secondaries_Energy){fHistograms->fMC_Pi0Secondaries_Energy->Fill(particle->Energy());}
	    if(fHistograms->fMC_Pi0Secondaries_R){fHistograms->fMC_Pi0Secondaries_R->Fill(particle->R());}
	    if(fHistograms->fMC_Pi0Secondaries_Z_R){fHistograms->fMC_Pi0Secondaries_Z_R->Fill(particle->Vz(),particle->R());}
	    if(fHistograms->fMC_Pi0Secondaries_OpeningAngleGamma){fHistograms->fMC_Pi0Secondaries_OpeningAngleGamma->Fill(GetMCOpeningAngle(daughter0,daughter1));}
	    if(fHistograms->fMC_Pi0Secondaries_X_Y){fHistograms->fMC_Pi0Secondaries_X_Y->Fill(particle->Vx(),particle->Vy());}//only fill from one daughter to avoid multiple filling
	  }
	  else{
	    if(fHistograms->fMC_Pi0_Eta){fHistograms->fMC_Pi0_Eta->Fill(particle->Eta());}
	    /*
	    Double_t tmpPhi=particle->Phi();
	    if(particle->Phi()> TMath::Pi()){
	      tmpPhi = particle->Phi()-(2*TMath::Pi());
	    }
	    */
	    if(fHistograms->fMC_Pi0_Phi){fHistograms->fMC_Pi0_Phi->Fill(tmpPhi);}
	    if(fHistograms->fMC_Pi0_Pt){fHistograms->fMC_Pi0_Pt->Fill(particle->Pt());}
	    if(fHistograms->fMC_Pi0_Energy){fHistograms->fMC_Pi0_Energy->Fill(particle->Energy());}
	    if(fHistograms->fMC_Pi0_R){fHistograms->fMC_Pi0_R->Fill(particle->R());}
	    if(fHistograms->fMC_Pi0_Z_R){fHistograms->fMC_Pi0_Z_R->Fill(particle->Vz(),particle->R());}
	    if(fHistograms->fMC_Pi0_OpeningAngleGamma){fHistograms->fMC_Pi0_OpeningAngleGamma->Fill(GetMCOpeningAngle(daughter0,daughter1));}
	    if(fHistograms->fMC_Pi0_X_Y){fHistograms->fMC_Pi0_X_Y->Fill(particle->Vx(),particle->Vy());}//only fill from one daughter to avoid multiple filling
	  }
	}
	else if(particle->GetPdgCode()==221){//Eta
	  if(fHistograms->fMC_Eta_Eta){fHistograms->fMC_Eta_Eta->Fill(particle->Eta());}
	  /*
	    Double_t tmpPhi=particle->Phi();
	  if(particle->Phi()> TMath::Pi()){
	    tmpPhi = particle->Phi()-(2*TMath::Pi());
	  }
	  */
	  if(fHistograms->fMC_Eta_Phi){fHistograms->fMC_Eta_Phi->Fill(tmpPhi);}
	  if(fHistograms->fMC_Eta_Pt){fHistograms->fMC_Eta_Pt->Fill(particle->Pt());}
	  if(fHistograms->fMC_Eta_Energy){fHistograms->fMC_Eta_Energy->Fill(particle->Energy());}
	  if(fHistograms->fMC_Eta_R){fHistograms->fMC_Eta_R->Fill(particle->R());}
	  if(fHistograms->fMC_Eta_Z_R){fHistograms->fMC_Eta_Z_R->Fill(particle->Vz(),particle->R());}
	  if(fHistograms->fMC_Eta_OpeningAngleGamma){fHistograms->fMC_Eta_OpeningAngleGamma->Fill(GetMCOpeningAngle(daughter0,daughter1));}
	  if(fHistograms->fMC_Eta_X_Y){fHistograms->fMC_Eta_X_Y->Fill(particle->Vx(),particle->Vy());}//only fill from one daughter to avoid multiple filling
	}
	
	//the match data should be filled no matter which mother the gamma-gamma comes from
	if(fHistograms->fMC_Match_Gamma_R){fHistograms->fMC_Match_Gamma_R->Fill(particle->R());}
	if(fHistograms->fMC_Match_Gamma_Z_R){fHistograms->fMC_Match_Gamma_Z_R->Fill(particle->Vz(),particle->R());}
	if(fHistograms->fMC_Match_Gamma_X_Y){fHistograms->fMC_Match_Gamma_X_Y->Fill(particle->Vx(),particle->Vy());}
	if(fHistograms->fMC_Match_Gamma_Mass){fHistograms->fMC_Match_Gamma_Mass->Fill(particle->GetCalcMass());}
	if(fHistograms->fMC_Match_Gamma_OpeningAngle){fHistograms->fMC_Match_Gamma_OpeningAngle->Fill(GetMCOpeningAngle(daughter0,daughter1));}
	if(fHistograms->fMC_Match_Gamma_Energy){fHistograms->fMC_Match_Gamma_Energy->Fill(particle->Energy());}
	if(fHistograms->fMC_Match_Gamma_Pt){fHistograms->fMC_Match_Gamma_Pt->Fill(particle->Pt());}
	if(fHistograms->fMC_Match_Gamma_Eta){fHistograms->fMC_Match_Gamma_Eta->Fill(particle->Eta());}
	/*
	Double_t tmpPhi=particle->Phi();
	if(particle->Phi()> TMath::Pi()){
	  tmpPhi = particle->Phi()-(2*TMath::Pi());
	}
	*/
	if(fHistograms->fMC_Match_Gamma_Phi){fHistograms->fMC_Match_Gamma_Phi->Fill(tmpPhi);}
      }
    }
  }
}

void AliAnalysisTaskGammaConversion::ProcessV0s(){
  Int_t nSurvivingV0s=0;
  while(fV0Reader->NextV0()){
    nSurvivingV0s++;
    //-------------------------- filling v0 information -------------------------------------
    if(fHistograms->fESD_EP_OpeningAngle){fHistograms->fESD_EP_OpeningAngle->Fill(fV0Reader->GetOpeningAngle());}    
    if(fHistograms->fESD_EP_R){fHistograms->fESD_EP_R->Fill(fV0Reader->GetXYRadius());}
    if(fHistograms->fESD_EP_Z_R){fHistograms->fESD_EP_Z_R->Fill(fV0Reader->GetZ(),fV0Reader->GetXYRadius());}
    if(fHistograms->fESD_EP_X_Y){fHistograms->fESD_EP_X_Y->Fill(fV0Reader->GetX(),fV0Reader->GetY());}
    
    
    if(fHistograms->fESD_E_Energy){fHistograms->fESD_E_Energy->Fill(fV0Reader->GetNegativeTrackEnergy());}
    if(fHistograms->fESD_E_Pt){fHistograms->fESD_E_Pt->Fill(fV0Reader->GetNegativeTrackPt());}
    if(fHistograms->fESD_E_Eta){fHistograms->fESD_E_Eta->Fill(fV0Reader->GetNegativeTrackEta());}
    if(fHistograms->fESD_E_Phi){fHistograms->fESD_E_Phi->Fill(fV0Reader->GetNegativeTrackPhi());}
    
    if(fHistograms->fESD_P_Energy){fHistograms->fESD_P_Energy->Fill(fV0Reader->GetPositiveTrackEnergy());}
    if(fHistograms->fESD_P_Pt){fHistograms->fESD_P_Pt->Fill(fV0Reader->GetPositiveTrackPt());}
    if(fHistograms->fESD_P_Eta){fHistograms->fESD_P_Eta->Fill(fV0Reader->GetPositiveTrackEta());}
    if(fHistograms->fESD_P_Phi){fHistograms->fESD_P_Phi->Fill(fV0Reader->GetPositiveTrackPhi());}
    
    if(fHistograms->fESD_Gamma_Energy){fHistograms->fESD_Gamma_Energy->Fill(fV0Reader->GetMotherCandidateEnergy());}
    if(fHistograms->fESD_Gamma_Pt){fHistograms->fESD_Gamma_Pt->Fill(fV0Reader->GetMotherCandidatePt());}
    if(fHistograms->fESD_Gamma_Eta){fHistograms->fESD_Gamma_Eta->Fill(fV0Reader->GetMotherCandidateEta());}
    if(fHistograms->fESD_Gamma_Phi){fHistograms->fESD_Gamma_Phi->Fill(fV0Reader->GetMotherCandidatePhi());}


    // begin mapping
    Int_t rBin    = fHistograms->GetRBin(fV0Reader->GetXYRadius());
    Int_t phiBin  = fHistograms->GetPhiBin(fV0Reader->GetNegativeTrackPhi());
    Double_t motherCandidateEta= fV0Reader->GetMotherCandidateEta();
    if(fHistograms->fESD_Mapping[phiBin][rBin] != NULL){fHistograms->fESD_Mapping[phiBin][rBin]->Fill(fV0Reader->GetZ(),motherCandidateEta);}
    if(fHistograms->fESD_Mapping_Phi[phiBin] != NULL){fHistograms->fESD_Mapping_Phi[phiBin]->Fill(fV0Reader->GetZ(),motherCandidateEta);}
    if(fHistograms->fESD_Mapping_R[rBin] != NULL){fHistograms->fESD_Mapping_R[rBin]->Fill(fV0Reader->GetZ(),motherCandidateEta);}
    // end mapping
    

    fKFReconstructedGammas.push_back(*fV0Reader->GetMotherCandidateKFCombination());

    //----------------------------------- checking for "real" conversions (MC match) --------------------------------------
    if(fDoMCTruth){
      if(fV0Reader->HasSameMCMother() == kFALSE){
	continue;
      }

      TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
      TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();
      if(negativeMC->GetPdgCode()!=11 || positiveMC->GetPdgCode()!=-11){
	continue;
      }

      if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
	if(fHistograms->fESD_Match_Gamma_X_Y){fHistograms->fESD_Match_Gamma_X_Y->Fill(fV0Reader->GetX(),fV0Reader->GetY());}
	if(fHistograms->fESD_Match_Gamma_OpeningAngle){fHistograms->fESD_Match_Gamma_OpeningAngle->Fill(fV0Reader->GetOpeningAngle());}
	if(fHistograms->fESD_Match_Gamma_Pt){fHistograms->fESD_Match_Gamma_Pt->Fill(fV0Reader->GetMotherCandidatePt());}
	if(fHistograms->fESD_Match_Gamma_Energy){fHistograms->fESD_Match_Gamma_Energy->Fill(fV0Reader->GetMotherCandidateEnergy());}
	if(fHistograms->fESD_Match_Gamma_Eta){fHistograms->fESD_Match_Gamma_Eta->Fill(fV0Reader->GetMotherCandidateEta());}

	if(fHistograms->fESD_Match_Gamma_Phi){fHistograms->fESD_Match_Gamma_Phi->Fill(fV0Reader->GetMotherCandidatePhi());}
	if(fHistograms->fESD_Match_Gamma_Mass){fHistograms->fESD_Match_Gamma_Mass->Fill(fV0Reader->GetMotherCandidateMass());}
	if(fHistograms->fESD_Match_Gamma_Width){fHistograms->fESD_Match_Gamma_Width->Fill(fV0Reader->GetMotherCandidateWidth());}
	if(fHistograms->fESD_Match_Gamma_Chi2){fHistograms->fESD_Match_Gamma_Chi2->Fill(fV0Reader->GetMotherCandidateChi2());}
	if(fHistograms->fESD_Match_Gamma_NDF){fHistograms->fESD_Match_Gamma_NDF->Fill(fV0Reader->GetMotherCandidateNDF());}
	if(fHistograms->fESD_Match_Gamma_R){fHistograms->fESD_Match_Gamma_R->Fill(fV0Reader->GetXYRadius());}
	if(fHistograms->fESD_Match_Gamma_Z_R){fHistograms->fESD_Match_Gamma_Z_R->Fill(fV0Reader->GetZ(),fV0Reader->GetXYRadius());}

	//resolution
	Double_t mc_pt   = fV0Reader->GetMotherMCParticle()->Pt();
	Double_t esd_pt  = fV0Reader->GetMotherCandidatePt();
	Double_t res_dPt = ((esd_pt - mc_pt)/mc_pt)*100;
	if(fHistograms->fResolution_dPt != NULL){fHistograms->fResolution_dPt->Fill(mc_pt,res_dPt);}
	if(fHistograms->fResolution_MC_Pt != NULL){fHistograms->fResolution_MC_Pt->Fill(mc_pt);}
	if(fHistograms->fResolution_ESD_Pt != NULL){fHistograms->fResolution_ESD_Pt->Fill(esd_pt);}
	
	Double_t res_dZ = ((fV0Reader->GetZ() -fV0Reader->GetNegativeMCParticle()->Vz())/fV0Reader->GetNegativeMCParticle()->Vz())*100;
	if(fHistograms->fResolution_dZ != NULL){fHistograms->fResolution_dZ->Fill(fV0Reader->GetNegativeMCParticle()->Vz(),res_dZ);}
	if(fHistograms->fResolution_MC_Z != NULL){fHistograms->fResolution_MC_Z->Fill(fV0Reader->GetNegativeMCParticle()->Vz());}
	if(fHistograms->fResolution_ESD_Z != NULL){fHistograms->fResolution_ESD_Z->Fill(fV0Reader->GetZ());}
	
	Double_t res_dR = ((fV0Reader->GetXYRadius() - fV0Reader->GetNegativeMCParticle()->R())/fV0Reader->GetNegativeMCParticle()->R())*100;
	if(fHistograms->fResolution_dR != NULL){fHistograms->fResolution_dR->Fill(fV0Reader->GetNegativeMCParticle()->R(),res_dR);}
	if(fHistograms->fResolution_MC_R != NULL){fHistograms->fResolution_MC_R->Fill(fV0Reader->GetNegativeMCParticle()->R());}
	if(fHistograms->fResolution_ESD_R != NULL){fHistograms->fResolution_ESD_R->Fill(fV0Reader->GetXYRadius());}
	if(fHistograms->fResolution_dR_dPt != NULL){fHistograms->fResolution_dR_dPt->Fill(res_dR,res_dPt);}
      }
    }
  }
  if(fHistograms->fNumberOfSurvivingV0s != NULL){fHistograms->fNumberOfSurvivingV0s->Fill(nSurvivingV0s);}
  if(fHistograms->fNumberOfV0s != NULL){fHistograms->fNumberOfV0s->Fill(fV0Reader->GetNumberOfV0s());}
}

void AliAnalysisTaskGammaConversion::ProcessGammasForNeutralMesonAnalysis(){

  for(UInt_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammas.size();firstGammaIndex++){
    for(UInt_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fKFReconstructedGammas.size();secondGammaIndex++){
      AliKFParticle * twoGammaDecayCandidateDaughter0 = &fKFReconstructedGammas[firstGammaIndex];
      AliKFParticle * twoGammaDecayCandidateDaughter1 = &fKFReconstructedGammas[secondGammaIndex];

      // Pi0's
      AliKFParticle *pi0Candidate = new AliKFParticle(*twoGammaDecayCandidateDaughter0,*twoGammaDecayCandidateDaughter1);

      //      pi0Candidate->SetMassConstraint(fPi0Mass,fPi0Width);

      Double_t massPi0 =0.;
      Double_t widthPi0 = 0.;
      Double_t chi2Pi0 =10000.;	
      pi0Candidate->GetMass(massPi0,widthPi0);
      if(pi0Candidate->GetNDF()>0){
	chi2Pi0 = pi0Candidate->GetChi2()/pi0Candidate->GetNDF();
	//	if(chi2Pi0>0 && chi2Pi0<fChi2Cut){//TODO  find this out
	if(chi2Pi0>0 && chi2Pi0<10000){//TODO  find this out see line above

	  /*
	  TVector3 vertexDaughter0(twoGammaDecayCandidateDaughter0->Px(),twoGammaDecayCandidateDaughter0->Py(),twoGammaDecayCandidateDaughter0->Pz());
	  TVector3 vertexDaughter1(twoGammaDecayCandidateDaughter1->Px(),twoGammaDecayCandidateDaughter1->Py(),twoGammaDecayCandidateDaughter1->Pz());
	  TVector3 vertexPi0Candidate(pi0Candidate->Px(),pi0Candidate->Py(),pi0Candidate->Pz());
	  
	  Double_t openingAnglePi0 = vertexDaughter0.Angle(vertexDaughter1);
	  
	  Double_t radiusPi0 = sqrt( vertexPi0Candidate.x()*vertexPi0Candidate.x() + vertexPi0Candidate.y()*vertexPi0Candidate.y() );
	  */

	  TVector3 vectorPi0Candidate(pi0Candidate->Px(),pi0Candidate->Py(),pi0Candidate->Pz());

	  Double_t openingAnglePi0 = twoGammaDecayCandidateDaughter0->GetAngle(*twoGammaDecayCandidateDaughter1);

	  //Double_t vtx000[3] = {0,0,0};
	  //   NOT SURE IF THIS DOES WHAT WE WANT: remember to ask Sergey  Double_t radiusPi0 = pi0Candidate->GetDistanceFromVertex(vtx000);
	  //Calculating by hand the radius
	  Double_t tmpX= pi0Candidate->GetX();
	  Double_t tmpY= pi0Candidate->GetY();
	  
	  Double_t radiusPi0 = TMath::Sqrt(tmpX*tmpX + tmpY*tmpY);

	  if(fHistograms->fESD_Pi0_OpeningAngleGamma){fHistograms->fESD_Pi0_OpeningAngleGamma->Fill(openingAnglePi0);}
	  if(fHistograms->fESD_Pi0_Energy){fHistograms->fESD_Pi0_Energy->Fill(pi0Candidate->GetE());}
	  if(fHistograms->fESD_Pi0_Pt){fHistograms->fESD_Pi0_Pt->Fill(sqrt(pi0Candidate->GetPx()*pi0Candidate->GetPx()+pi0Candidate->GetPy()*pi0Candidate->GetPy()));}
	  if(fHistograms->fESD_Pi0_Eta){fHistograms->fESD_Pi0_Eta->Fill(vectorPi0Candidate.Eta());}
	  if(fHistograms->fESD_Pi0_Phi){fHistograms->fESD_Pi0_Phi->Fill(vectorPi0Candidate.Phi());}
	  if(fHistograms->fESD_Pi0_Mass){fHistograms->fESD_Pi0_Mass->Fill(massPi0);}
	  if(fHistograms->fESD_Pi0_R){fHistograms->fESD_Pi0_R->Fill(radiusPi0);}
	  if(fHistograms->fESD_Pi0_Z_R){fHistograms->fESD_Pi0_Z_R->Fill(tmpY,radiusPi0);}
	  if(fHistograms->fESD_Pi0_X_Y){fHistograms->fESD_Pi0_X_Y->Fill(tmpX,tmpY);}	  
	}
      }
      delete pi0Candidate;

      //Eta's
      AliKFParticle *etaCandidate = new AliKFParticle(*twoGammaDecayCandidateDaughter0,*twoGammaDecayCandidateDaughter1);

      //      etaCandidate->SetMassConstraint(fEtaMass,fEtaWidth);
      Double_t massEta =0.;
      Double_t widthEta = 0.;
      Double_t chi2Eta =10000.;	
      etaCandidate->GetMass(massEta,widthEta);
      if(etaCandidate->GetNDF()>0){
	chi2Eta = etaCandidate->GetChi2()/etaCandidate->GetNDF();
	//	if(chi2Eta>0 && chi2Eta<fChi2Cut){
	if(chi2Eta>0 && chi2Eta<100000){// TODO find this out se line above
	  /*
 	  TVector3 vertexDaughter0(twoGammaDecayCandidateDaughter0->Px(),twoGammaDecayCandidateDaughter0->Py(),twoGammaDecayCandidateDaughter0->Pz());
	  TVector3 vertexDaughter1(twoGammaDecayCandidateDaughter1->Px(),twoGammaDecayCandidateDaughter1->Py(),twoGammaDecayCandidateDaughter1->Pz());
	  TVector3 vertexEtaCandidate(etaCandidate->Px(),etaCandidate->Py(),etaCandidate->Pz());
	  Double_t openingAngleEta = vertexDaughter0.Angle(vertexDaughter1);
	  
	  Double_t radiusEta = sqrt( etaCandidate->GetX()*etaCandidate->GetX() + vertexEtaCandidate.y()*vertexEtaCandidate.y() );
	  */

	  //Calculating by hand the radius
	  Double_t tmpX= etaCandidate->GetX();
	  Double_t tmpY= etaCandidate->GetY();
	  
	  Double_t radiusEta = TMath::Sqrt(tmpX*tmpX + tmpY*tmpY);
	  
	  Double_t openingAngleEta = twoGammaDecayCandidateDaughter0->GetAngle(*twoGammaDecayCandidateDaughter1);

	  TVector3 vectorEtaCandidate(etaCandidate->Px(),etaCandidate->Py(),etaCandidate->Pz());

	  if(fHistograms->fESD_Eta_OpeningAngleGamma){fHistograms->fESD_Eta_OpeningAngleGamma->Fill(openingAngleEta);}
	  if(fHistograms->fESD_Eta_Energy){fHistograms->fESD_Eta_Energy->Fill(etaCandidate->GetE());}
	  if(fHistograms->fESD_Eta_Pt){fHistograms->fESD_Eta_Pt->Fill(sqrt(etaCandidate->GetPx()*etaCandidate->GetPx()+etaCandidate->GetPy()*etaCandidate->GetPy()));}
	  if(fHistograms->fESD_Eta_Eta){fHistograms->fESD_Eta_Eta->Fill(vectorEtaCandidate.Eta());}
	  if(fHistograms->fESD_Eta_Phi){fHistograms->fESD_Eta_Phi->Fill(vectorEtaCandidate.Phi());}
	  if(fHistograms->fESD_Eta_Mass){fHistograms->fESD_Eta_Mass->Fill(massEta);}
	  if(fHistograms->fESD_Eta_R){fHistograms->fESD_Eta_R->Fill(radiusEta);}
	  if(fHistograms->fESD_Eta_Z_R){fHistograms->fESD_Eta_Z_R->Fill(etaCandidate->GetZ(),radiusEta);}
	  if(fHistograms->fESD_Eta_X_Y){fHistograms->fESD_Eta_X_Y->Fill(tmpX,tmpY);}	  
	}
      }
      delete etaCandidate;
    }
  }

}

void AliAnalysisTaskGammaConversion::CalculateBackground(){

  for(UInt_t iCurrent=0;iCurrent<fV0Reader->fCurrentEventGoodV0s.size();iCurrent++){
    AliKFParticle * currentEventGoodV0 = &fV0Reader->fCurrentEventGoodV0s.at(iCurrent);
    for(UInt_t iPrevious=0;iPrevious<fV0Reader->fCurrentEventGoodV0s.size();iPrevious++){
      AliKFParticle * previousEventGoodV0 = &fV0Reader->fPreviousEventGoodV0s.at(iPrevious);

      AliKFParticle *backgroundCandidate = new AliKFParticle(*currentEventGoodV0,*previousEventGoodV0);

      Double_t massBG =0.;
      Double_t widthBG = 0.;
      Double_t chi2BG =10000.;	
      backgroundCandidate->GetMass(massBG,widthBG);
      if(backgroundCandidate->GetNDF()>0){
	chi2BG = backgroundCandidate->GetChi2()/backgroundCandidate->GetNDF();
	//	if(chi2Pi0>0 && chi2Pi0<fChi2Cut){//TODO  find this out
	if(chi2BG>0 && chi2BG<fV0Reader->GetChi2Cut()){//TODO  find this out see line above

	  TVector3 vectorBGCandidate(backgroundCandidate->Px(),backgroundCandidate->Py(),backgroundCandidate->Pz());

	  Double_t openingAngleBG = currentEventGoodV0->GetAngle(*previousEventGoodV0);

	  //Calculating by hand the radius (find a better way)
	  Double_t tmpX= backgroundCandidate->GetX();
	  Double_t tmpY= backgroundCandidate->GetY();
	  
	  Double_t radiusBG = TMath::Sqrt(tmpX*tmpX + tmpY*tmpY);

	  if(fHistograms->fESD_Background_OpeningAngleGamma){fHistograms->fESD_Background_OpeningAngleGamma->Fill(openingAngleBG);}
	  if(fHistograms->fESD_Background_Energy){fHistograms->fESD_Background_Energy->Fill(backgroundCandidate->GetE());}
	  if(fHistograms->fESD_Background_Pt){fHistograms->fESD_Background_Pt->Fill(sqrt(backgroundCandidate->GetPx()*backgroundCandidate->GetPx()+backgroundCandidate->GetPy()*backgroundCandidate->GetPy()));}
	  if(fHistograms->fESD_Background_Eta){fHistograms->fESD_Background_Eta->Fill(vectorBGCandidate.Eta());}
	  if(fHistograms->fESD_Background_Phi){fHistograms->fESD_Background_Phi->Fill(vectorBGCandidate.Phi());}
	  if(fHistograms->fESD_Background_Mass){fHistograms->fESD_Background_Mass->Fill(massBG);}
	  if(fHistograms->fESD_Background_R){fHistograms->fESD_Background_R->Fill(radiusBG);}
	  if(fHistograms->fESD_Background_Z_R){fHistograms->fESD_Background_Z_R->Fill(tmpY,radiusBG);}
	  if(fHistograms->fESD_Background_X_Y){fHistograms->fESD_Background_X_Y->Fill(tmpX,tmpY);}
	}
      }
      delete backgroundCandidate;   
    }
  }
}

void AliAnalysisTaskGammaConversion::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(1,"Do nothing in Terminate");
}

void AliAnalysisTaskGammaConversion::UserCreateOutputObjects()
{
  // Create the output container

  fOutputContainer = fHistograms->GetOutputContainer();
  fOutputContainer->SetName(GetName()) ;  
}

Double_t AliAnalysisTaskGammaConversion::GetMCOpeningAngle(TParticle* daughter0,TParticle* daughter1){
  //helper function
  TVector3 V3D0(daughter0->Px(),daughter0->Py(),daughter0->Pz());
  TVector3 V3D1(daughter1->Px(),daughter1->Py(),daughter1->Pz());
  return V3D0.Angle(V3D1);
}

