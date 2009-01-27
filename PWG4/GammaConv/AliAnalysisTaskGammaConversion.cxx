/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.1                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

// root
#include <TChain.h>

// analysis
#include "AliAnalysisTaskGammaConversion.h"
#include "AliStack.h"
#include "AliLog.h"
#include <vector>

class AliKFVertex;
class AliAODHandler;
class AliAODEvent;
class ALiESDEvent;
class AliMCEvent;
class AliMCEventHandler;
class AliESDInputHandler;
class AliAnalysisManager;
class Riostream;
class TFile;
class TInterpreter;
class TSystem;
class TROOT;

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
  fMCGammaChic(),
  fKFReconstructedGammas(),
  fElectronMass(-1),
  fGammaMass(-1),
  fPi0Mass(-1),
  fEtaMass(-1),
  fGammaWidth(-1),
  fPi0Width(-1),
  fEtaWidth(-1),
  fCalculateBackground(kFALSE),
  fWriteNtuple(kFALSE),
  fGammaNtuple(NULL),
  fNeutralMesonNtuple(NULL),
  fTotalNumberOfAddedNtupleEntries(0)
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
  fMCGammaChic(),
  fKFReconstructedGammas(),
  fElectronMass(-1),
  fGammaMass(-1),
  fPi0Mass(-1),
  fEtaMass(-1),
  fGammaWidth(-1),
  fPi0Width(-1),
  fEtaWidth(-1),
  fCalculateBackground(kFALSE),
  fWriteNtuple(kFALSE),
  fGammaNtuple(NULL),
  fNeutralMesonNtuple(NULL),
  fTotalNumberOfAddedNtupleEntries(0)
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
  AliLog::SetGlobalLogLevel(AliLog::kError);
}


void AliAnalysisTaskGammaConversion::Exec(Option_t */*option*/)
{
  // Execute analysis for current event
  
  ConnectInputData("");
  
  //clear vectors
  fMCAllGammas.clear();
  fMCPi0s.clear();
  fMCEtas.clear();
  fMCGammaChic.clear();

  fKFReconstructedGammas.clear();

  //Clear the data in the v0Reader
  fV0Reader->UpdateEventByEventData();

  // Process the MC information
  if(fDoMCTruth){
    ProcessMCData();
  }

  // Process the v0 information
  ProcessV0s();

  //calculate background if flag is set
  if(fCalculateBackground){
    CalculateBackground();
  }

  // Process reconstructed gammas
  ProcessGammasForNeutralMesonAnalysis();

  PostData(1, fOutputContainer);
  
}

void AliAnalysisTaskGammaConversion::ConnectInputData(Option_t */*option*/){
  // see header file for documentation

  if(fV0Reader == NULL){
    // Write warning here cuts and so on are default if this ever happens
  }
  fV0Reader->Initialize();
}

void AliAnalysisTaskGammaConversion::ProcessMCData(){
  // see header file for documentation
  
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
	  fHistograms->FillHistogram("MC_Gamma_Energy", particle->Energy());
	  fHistograms->FillHistogram("MC_Gamma_Pt", particle->Pt());
	  fHistograms->FillHistogram("MC_Gamma_Eta", particle->Eta());
	  
	  fHistograms->FillHistogram("MC_Gamma_Phi", tmpPhi);

	  //adding the conversion points from all gammas with e+e- daughters
	  if(particle->GetNDaughters() >= 2){
	    TParticle* daughter0 = NULL;
	    TParticle* daughter1 = NULL;
	    
	    for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	      TParticle *tmpDaughter = fStack->Particle(daughterIndex);
	      if(tmpDaughter->GetUniqueID() == 5){
		if(tmpDaughter->GetPdgCode() == 11){
                  daughter0 = tmpDaughter;
		}
		else if(tmpDaughter->GetPdgCode() == -11){
		  daughter1 = tmpDaughter;
		}
	      }
	    }

	    if(daughter0 == NULL || daughter1 == NULL){ // means we do not have two daughters from pair production
	      continue;
	    }

	    if(daughter0->R()>fV0Reader->GetMaxRCut() || daughter1->R()>fV0Reader->GetMaxRCut()){
	      continue;
	    }
	    
	    if((daughter0->GetPdgCode() == -11 && daughter1->GetPdgCode()) == 11 ||
	       (daughter0->GetPdgCode() == 11 && daughter1->GetPdgCode() == -11)){

	      // begin Mapping 
	      Int_t rBin    = fHistograms->GetRBin(daughter0->R());
	      Int_t phiBin  = fHistograms->GetPhiBin(daughter0->Phi());
	      
	      TString nameMCMappingPhiR="";
	      nameMCMappingPhiR.Form("MC_EP_Mapping-Phi%02d-R%02d",phiBin,rBin);
	      fHistograms->FillHistogram(nameMCMappingPhiR, daughter0->Vz(), particle->Eta());
	      
	      TString nameMCMappingPhi="";
	      nameMCMappingPhi.Form("MC_EP_Mapping-Phi%02d",phiBin);
	      fHistograms->FillHistogram(nameMCMappingPhi, particle->Eta());
	      
	      TString nameMCMappingR="";
	      nameMCMappingR.Form("MC_EP_Mapping-R%02d",rBin);
	      fHistograms->FillHistogram(nameMCMappingR, particle->Eta());
	      
	      TString nameMCMappingPhiInR="";
	      nameMCMappingPhiInR.Form("MC_EP_Mapping_Phi_R-%02d",rBin);
	      fHistograms->FillHistogram(nameMCMappingPhiInR, tmpPhi);
	      //end mapping

	      fHistograms->FillHistogram("MC_EP_R",daughter0->R());
	      fHistograms->FillHistogram("MC_EP_ZR",daughter0->Vz(),daughter0->R());
	      fHistograms->FillHistogram("MC_EP_XY",daughter0->Vx(),daughter0->Vy());
	      fHistograms->FillHistogram("MC_EP_OpeningAngle",GetMCOpeningAngle(daughter0, daughter1));
	    }// end if((daughter0->GetPdgCode() == -11 && daughter1->GetPdgCode()) == 11 ||....... approx 20 lines above
	  }// end if(particle->GetNDaughters() >= 2){
	} // end if( fStack->Particle(particle->GetMother(0))->GetPdgCode() != 22 )
	if( fStack->Particle(particle->GetMother(0))->GetPdgCode()==10441 ||//chic0 
	    fStack->Particle(particle->GetMother(0))->GetPdgCode()==20443 ||//psi2S
	    fStack->Particle(particle->GetMother(0))->GetPdgCode()==445  //chic2
	    ){ 
	  fMCGammaChic.push_back(particle);
	}
      }// end if(particle->GetMother(0)>-1)
      else{//means we have a primary particle
	fHistograms->FillHistogram("MC_DirectGamma_Energy",particle->Energy());
	fHistograms->FillHistogram("MC_DirectGamma_Pt", particle->Pt());
	fHistograms->FillHistogram("MC_DirectGamma_Eta", particle->Eta());
	fHistograms->FillHistogram("MCDirectGammaPhi", tmpPhi);

	//adding the conversion points from all gammas with e+e- daughters
	if(particle->GetNDaughters() == 2){
	  TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
	  TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
	  if((daughter0->GetPdgCode() == -11 && daughter1->GetPdgCode() == 11) ||
	     (daughter0->GetPdgCode() == 11 && daughter1->GetPdgCode() == -11)){
	    
	    fHistograms->FillHistogram("MC_EP_R",daughter0->R());
	    fHistograms->FillHistogram("MC_EP_ZR",daughter0->Vz(),daughter0->R());
	    fHistograms->FillHistogram("MC_EP_XY",daughter0->Vx(),daughter0->Vy());
	    fHistograms->FillHistogram("MC_EP_OpeningAngle",GetMCOpeningAngle(daughter0, daughter1));

	  }
	}
      }// end else
    }// end if (particle->GetPdgCode()== 22){
    else if (TMath::Abs(particle->GetPdgCode())== 11){ // Means we have an electron or a positron
      if(particle->GetMother(0)>-1){ // means we have a mother
	if( fStack->Particle(particle->GetMother(0))->GetPdgCode()==22 ){ // Means we have a gamma mother
	  if(particle->GetPdgCode() == 11){//electron 
	    fHistograms->FillHistogram("MC_E_Energy", particle->Energy());
	    fHistograms->FillHistogram("MC_E_Pt", particle->Pt());
	    fHistograms->FillHistogram("MC_E_Eta", particle->Eta());
	    fHistograms->FillHistogram("MC_E_Phi", tmpPhi);
	  }
	  if(particle->GetPdgCode() == -11){//positron 
	    fHistograms->FillHistogram("MC_P_Energy", particle->Energy());
	    fHistograms->FillHistogram("MC_P_Pt", particle->Pt());
	    fHistograms->FillHistogram("MC_P_Eta", particle->Eta());
	    fHistograms->FillHistogram("MC_P_Phi", tmpPhi);
	  }
	}
      }
    } // end else if (TMath::Abs(particle->GetPdgCode())== 11)
    else if(particle->GetNDaughters() == 2){

      TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
      TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
      if(daughter0->GetPdgCode() == 22 && daughter1->GetPdgCode() == 22){//check for gamma gamma daughters
	
	if(particle->GetPdgCode()==111){//Pi0
	  
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_Phi", tmpPhi);

	  if( iTracks >= fStack->GetNprimary()){
	    
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_Eta", particle->Eta());

	    fHistograms->FillHistogram("MC_Pi0_Secondaries_Phi", tmpPhi);
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt", particle->Pt());
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_Energy", particle->Energy());
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_R", particle->R());
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_ZR", particle->Vz(),particle->R());
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_XY", particle->Vx(),particle->Vy());//only fill from one daughter to avoid multiple filling
	  }
	  else{
	    fHistograms->FillHistogram("MC_Pi0_Eta", particle->Eta());

	    fHistograms->FillHistogram("MC_Pi0_Phi", tmpPhi);
	    fHistograms->FillHistogram("MC_Pi0_Pt", particle->Pt());
	    fHistograms->FillHistogram("MC_Pi0_Energy", particle->Energy());
	    fHistograms->FillHistogram("MC_Pi0_R", particle->R());
	    fHistograms->FillHistogram("MC_Pi0_ZR", particle->Vz(),particle->R());
	    fHistograms->FillHistogram("MC_Pi0_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	    fHistograms->FillHistogram("MC_Pi0_XY", particle->Vx(), particle->Vy());//only fill from one daughter to avoid multiple filling
	  }
	}
	else if(particle->GetPdgCode()==221){//Eta
	  fHistograms->FillHistogram("MC_Eta_Eta", particle->Eta());

	  fHistograms->FillHistogram("MC_Eta_Phi",tmpPhi);
	  fHistograms->FillHistogram("MC_Eta_Pt", particle->Pt());
	  fHistograms->FillHistogram("MC_Eta_Energy", particle->Energy());
	  fHistograms->FillHistogram("MC_Eta_R", particle->R());
	  fHistograms->FillHistogram("MC_Eta_ZR", particle->Vz(),particle->R());
	  fHistograms->FillHistogram("MC_Eta_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	  fHistograms->FillHistogram("MC_Eta_XY", particle->Vx(), particle->Vy());//only fill from one daughter to avoid multiple filling
	}
	
	//the match data should be filled no matter which mother the gamma-gamma comes from
	fHistograms->FillHistogram("MC_Match_Gamma_R", particle->R());
	fHistograms->FillHistogram("MC_Match_Gamma_ZR", particle->Vz(),particle->R());
	fHistograms->FillHistogram("MC_Match_Gamma_XY", particle->Vx(),particle->Vy());
	fHistograms->FillHistogram("MC_Match_Gamma_Mass", particle->GetCalcMass());
	fHistograms->FillHistogram("MC_Match_Gamma_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	fHistograms->FillHistogram("MC_Match_Gamma_Energy", particle->Energy());
	fHistograms->FillHistogram("MC_Match_Gamma_Pt", particle->Pt());
	fHistograms->FillHistogram("MC_Match_Gamma_Eta", particle->Eta());
	fHistograms->FillHistogram("MC_Match_Gamma_Phi",tmpPhi);
      }
    }// end else if(particle->GetNDaughters() == 2)
  }// end for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++)
} // end ProcessMCData

void AliAnalysisTaskGammaConversion::FillNtuple(){

  if(fGammaNtuple == NULL){
    return;
  }
  Int_t numberOfV0s = fV0Reader->GetNumberOfV0s();
  for(Int_t i=0;i<numberOfV0s;i++){
    Float_t values[27] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    AliESDv0 * cV0 = fV0Reader->GetV0(i);
    Double_t negPID=0;
    Double_t posPID=0;
    fV0Reader->GetPIDProbability(negPID,posPID);
    values[0]=cV0->GetOnFlyStatus();
    values[1]=fV0Reader->CheckForPrimaryVertex();
    values[2]=negPID;
    values[3]=posPID;
    values[4]=fV0Reader->GetX();
    values[5]=fV0Reader->GetY();
    values[6]=fV0Reader->GetZ();
    values[7]=fV0Reader->GetXYRadius();
    values[8]=fV0Reader->GetMotherCandidateNDF();
    values[9]=fV0Reader->GetMotherCandidateChi2();
    values[10]=fV0Reader->GetMotherCandidateEnergy();
    values[11]=fV0Reader->GetMotherCandidateEta();
    values[12]=fV0Reader->GetMotherCandidatePt();
    values[13]=fV0Reader->GetMotherCandidateMass();
    values[14]=fV0Reader->GetMotherCandidateWidth();
    //    values[15]=fV0Reader->GetMotherMCParticle()->Pt();   MOVED TO THE END, HAS TO BE CALLED AFTER HasSameMother NB: still has the same entry in the array
    values[16]=fV0Reader->GetOpeningAngle();
    values[17]=fV0Reader->GetNegativeTrackEnergy();
    values[18]=fV0Reader->GetNegativeTrackPt();
    values[19]=fV0Reader->GetNegativeTrackEta();
    values[20]=fV0Reader->GetNegativeTrackPhi();
    values[21]=fV0Reader->GetPositiveTrackEnergy();
    values[22]=fV0Reader->GetPositiveTrackPt();
    values[23]=fV0Reader->GetPositiveTrackEta();
    values[24]=fV0Reader->GetPositiveTrackPhi();
    values[25]=fV0Reader->HasSameMCMother();
    if(values[25] != 0){
      values[26]=fV0Reader->GetMotherMCParticlePDGCode();
      values[15]=fV0Reader->GetMotherMCParticle()->Pt();
    }
    fTotalNumberOfAddedNtupleEntries++;
    fGammaNtuple->Fill(values);
  }
  fV0Reader->ResetV0IndexNumber();
  
}

void AliAnalysisTaskGammaConversion::ProcessV0s(){
  // see header file for documentation

  if(fWriteNtuple == kTRUE){
    FillNtuple();
  }

  Int_t nSurvivingV0s=0;
  while(fV0Reader->NextV0()){
    nSurvivingV0s++;
    //-------------------------- filling v0 information -------------------------------------
    fHistograms->FillHistogram("ESD_EP_OpeningAngle", fV0Reader->GetOpeningAngle());    
    fHistograms->FillHistogram("ESD_EP_R", fV0Reader->GetXYRadius());
    fHistograms->FillHistogram("ESD_EP_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
    fHistograms->FillHistogram("ESD_EP_XY", fV0Reader->GetX(),fV0Reader->GetY());
    
    
    fHistograms->FillHistogram("ESD_E_Energy", fV0Reader->GetNegativeTrackEnergy());
    fHistograms->FillHistogram("ESD_E_Pt", fV0Reader->GetNegativeTrackPt());
    fHistograms->FillHistogram("ESD_E_Eta", fV0Reader->GetNegativeTrackEta());
    fHistograms->FillHistogram("ESD_E_Phi", fV0Reader->GetNegativeTrackPhi());
    
    fHistograms->FillHistogram("ESD_P_Energy", fV0Reader->GetPositiveTrackEnergy());
    fHistograms->FillHistogram("ESD_P_Pt", fV0Reader->GetPositiveTrackPt());
    fHistograms->FillHistogram("ESD_P_Eta", fV0Reader->GetPositiveTrackEta());
    fHistograms->FillHistogram("ESD_P_Phi", fV0Reader->GetPositiveTrackPhi());
    
    fHistograms->FillHistogram("ESD_Gamma_Energy", fV0Reader->GetMotherCandidateEnergy());
    fHistograms->FillHistogram("ESD_Gamma_Pt", fV0Reader->GetMotherCandidatePt());
    fHistograms->FillHistogram("ESD_Gamma_Eta", fV0Reader->GetMotherCandidateEta());
    fHistograms->FillHistogram("ESD_Gamma_Phi", fV0Reader->GetMotherCandidatePhi());


    // begin mapping
    Int_t rBin    = fHistograms->GetRBin(fV0Reader->GetXYRadius());
    Int_t phiBin  = fHistograms->GetPhiBin(fV0Reader->GetNegativeTrackPhi());
    Double_t motherCandidateEta= fV0Reader->GetMotherCandidateEta();
    
    TString nameESDMappingPhiR="";
    nameESDMappingPhiR.Form("ESD_EP_Mapping-Phi%02d-R%02d",phiBin,rBin);
    fHistograms->FillHistogram(nameESDMappingPhiR, fV0Reader->GetZ(), motherCandidateEta);

    TString nameESDMappingPhi="";
    nameESDMappingPhi.Form("ESD_EP_Mapping-Phi%02d",phiBin);
    fHistograms->FillHistogram(nameESDMappingPhi, fV0Reader->GetZ(), motherCandidateEta);

    TString nameESDMappingR="";
    nameESDMappingR.Form("ESD_EP_Mapping-R%02d",rBin);
    fHistograms->FillHistogram(nameESDMappingR, fV0Reader->GetZ(), motherCandidateEta);  

    TString nameESDMappingPhiInR="";
    nameESDMappingPhiInR.Form("ESD_EP_Mapping_Phi_R-%02d",rBin);
    fHistograms->FillHistogram(nameESDMappingPhiInR, fV0Reader->GetMotherCandidatePhi());
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
	fHistograms->FillHistogram("ESD_Match_Gamma_XY", fV0Reader->GetX(),fV0Reader->GetY());
	fHistograms->FillHistogram("ESD_Match_Gamma_OpeningAngle", fV0Reader->GetOpeningAngle());
	fHistograms->FillHistogram("ESD_Match_Gamma_Pt", fV0Reader->GetMotherCandidatePt());
	fHistograms->FillHistogram("ESD_Match_Gamma_Energy", fV0Reader->GetMotherCandidateEnergy());
	fHistograms->FillHistogram("ESD_Match_Gamma_Eta", fV0Reader->GetMotherCandidateEta());

	fHistograms->FillHistogram("ESD_Match_Gamma_Phi", fV0Reader->GetMotherCandidatePhi());
	fHistograms->FillHistogram("ESD_Match_Gamma_Mass", fV0Reader->GetMotherCandidateMass());
	fHistograms->FillHistogram("ESD_Match_Gamma_Width", fV0Reader->GetMotherCandidateWidth());
	fHistograms->FillHistogram("ESD_Match_Gamma_Chi2", fV0Reader->GetMotherCandidateChi2());
	fHistograms->FillHistogram("ESD_Match_Gamma_NDF", fV0Reader->GetMotherCandidateNDF());
	fHistograms->FillHistogram("ESD_Match_Gamma_R", fV0Reader->GetXYRadius());
	fHistograms->FillHistogram("ESD_Match_Gamma_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());

	//resolution
	Double_t mcpt   = fV0Reader->GetMotherMCParticle()->Pt();
	Double_t esdpt  = fV0Reader->GetMotherCandidatePt();
	Double_t resdPt = 0;
	if(mcpt != 0){
	  resdPt = ((esdpt - mcpt)/mcpt)*100;
	}

	fHistograms->FillHistogram("Resolution_dPt", mcpt, resdPt);
	fHistograms->FillHistogram("Resolution_MC_Pt", mcpt);
	fHistograms->FillHistogram("Resolution_ESD_Pt", esdpt);
	
	Double_t resdZ = 0;
	if(fV0Reader->GetNegativeMCParticle()->Vz() != 0){
	  resdZ = ((fV0Reader->GetZ() -fV0Reader->GetNegativeMCParticle()->Vz())/fV0Reader->GetNegativeMCParticle()->Vz())*100;
	}
	
	fHistograms->FillHistogram("Resolution_dZ", fV0Reader->GetNegativeMCParticle()->Vz(), resdZ);
	fHistograms->FillHistogram("Resolution_MC_Z", fV0Reader->GetNegativeMCParticle()->Vz());
	fHistograms->FillHistogram("Resolution_ESD_Z", fV0Reader->GetZ());
	
	Double_t resdR = 0;
	if(fV0Reader->GetNegativeMCParticle()->R() != 0){
	  resdR = ((fV0Reader->GetXYRadius() - fV0Reader->GetNegativeMCParticle()->R())/fV0Reader->GetNegativeMCParticle()->R())*100;
	}
	fHistograms->FillHistogram("Resolution_dR", fV0Reader->GetNegativeMCParticle()->R(), resdR);
	fHistograms->FillHistogram("Resolution_MC_R", fV0Reader->GetNegativeMCParticle()->R());
	fHistograms->FillHistogram("Resolution_ESD_R", fV0Reader->GetXYRadius());
	fHistograms->FillHistogram("Resolution_dR_dPt", resdR, resdPt);
      }
    }
  }
  fHistograms->FillHistogram("NumberOfSurvivingV0s", nSurvivingV0s);
  fHistograms->FillHistogram("NumberOfV0s", fV0Reader->GetNumberOfV0s());
}

void AliAnalysisTaskGammaConversion::ProcessGammasForNeutralMesonAnalysis(){
  // see header file for documentation

  for(UInt_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammas.size();firstGammaIndex++){
    for(UInt_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fKFReconstructedGammas.size();secondGammaIndex++){
      AliKFParticle * twoGammaDecayCandidateDaughter0 = &fKFReconstructedGammas[firstGammaIndex];
      AliKFParticle * twoGammaDecayCandidateDaughter1 = &fKFReconstructedGammas[secondGammaIndex];

      
      AliKFParticle *twoGammaCandidate = new AliKFParticle(*twoGammaDecayCandidateDaughter0,*twoGammaDecayCandidateDaughter1);

      Double_t massTwoGammaCandidate =0.;
      Double_t widthTwoGammaCandidate = 0.;
      Double_t chi2TwoGammaCandidate =10000.;	
      twoGammaCandidate->GetMass(massTwoGammaCandidate,widthTwoGammaCandidate);
      if(twoGammaCandidate->GetNDF()>0){
	chi2TwoGammaCandidate = twoGammaCandidate->GetChi2()/twoGammaCandidate->GetNDF();
	if(chi2TwoGammaCandidate>0 && chi2TwoGammaCandidate<fV0Reader->GetChi2CutMeson()){

	  TVector3 vectorTwoGammaCandidate(twoGammaCandidate->Px(),twoGammaCandidate->Py(),twoGammaCandidate->Pz());

	  Double_t openingAngleTwoGammaCandidate = twoGammaDecayCandidateDaughter0->GetAngle(*twoGammaDecayCandidateDaughter1);

	  //Calculating by hand the radius
	  Double_t tmpX= twoGammaCandidate->GetX();
	  Double_t tmpY= twoGammaCandidate->GetY();
	  
	  Double_t radiusTwoGammaCandidate = TMath::Sqrt(tmpX*tmpX + tmpY*tmpY);

	  fHistograms->FillHistogram("ESD_TwoGammaCombination_GammaDaughter_OpeningAngle", openingAngleTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_Energy", twoGammaCandidate->GetE());
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_Pt", sqrt(twoGammaCandidate->GetPx()*twoGammaCandidate->GetPx()+twoGammaCandidate->GetPy()*twoGammaCandidate->GetPy()));
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_Eta", vectorTwoGammaCandidate.Eta());
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_Phi", vectorTwoGammaCandidate.Phi());
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_Mass", massTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_R", radiusTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_ZR", tmpY, radiusTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_TwoGammaCombination_XY", tmpX, tmpY);
	  fHistograms->FillHistogram("InvMass_vs_Pt_Spectra",massTwoGammaCandidate ,sqrt(twoGammaCandidate->GetPx()*twoGammaCandidate->GetPx()+twoGammaCandidate->GetPy()*twoGammaCandidate->GetPy()));
	}
      }
      delete twoGammaCandidate;
    }
  }
}

void AliAnalysisTaskGammaConversion::CalculateBackground(){
  // see header file for documentation

  vector<AliKFParticle> vectorCurrentEventGoodV0s = fV0Reader->GetCurrentEventGoodV0s();
  vector<AliKFParticle> vectorPreviousEventGoodV0s = fV0Reader->GetPreviousEventGoodV0s();
  for(UInt_t iCurrent=0;iCurrent<vectorCurrentEventGoodV0s.size();iCurrent++){
    AliKFParticle * currentEventGoodV0 = &vectorCurrentEventGoodV0s.at(iCurrent);
    for(UInt_t iPrevious=0;iPrevious<vectorPreviousEventGoodV0s.size();iPrevious++){
      AliKFParticle * previousGoodV0 = &vectorPreviousEventGoodV0s.at(iPrevious);

      AliKFParticle *backgroundCandidate = new AliKFParticle(*currentEventGoodV0,*previousGoodV0);

      Double_t massBG =0.;
      Double_t widthBG = 0.;
      Double_t chi2BG =10000.;	
      backgroundCandidate->GetMass(massBG,widthBG);
      if(backgroundCandidate->GetNDF()>0){
	chi2BG = backgroundCandidate->GetChi2()/backgroundCandidate->GetNDF();
	if(chi2BG>0 && chi2BG<fV0Reader->GetChi2CutMeson()){

	  TVector3 vectorBGCandidate(backgroundCandidate->Px(),backgroundCandidate->Py(),backgroundCandidate->Pz());

	  Double_t openingAngleBG = currentEventGoodV0->GetAngle(*previousGoodV0);

	  //Calculating by hand the radius (find a better way)
	  Double_t tmpX= backgroundCandidate->GetX();
	  Double_t tmpY= backgroundCandidate->GetY();
	  
	  Double_t radiusBG = TMath::Sqrt(tmpX*tmpX + tmpY*tmpY);

	  fHistograms->FillHistogram("ESD_Background_GammaDaughter_OpeningAngle", openingAngleBG);
	  fHistograms->FillHistogram("ESD_Background_Energy", backgroundCandidate->GetE());
	  fHistograms->FillHistogram("ESD_Background_Pt", sqrt(backgroundCandidate->GetPx()*backgroundCandidate->GetPx()+backgroundCandidate->GetPy()*backgroundCandidate->GetPy()));
	  fHistograms->FillHistogram("ESD_Background_Eta", vectorBGCandidate.Eta());
	  fHistograms->FillHistogram("ESD_Background_Phi", vectorBGCandidate.Phi());
	  fHistograms->FillHistogram("ESD_Background_Mass", massBG);
	  fHistograms->FillHistogram("ESD_Background_R", radiusBG);
	  fHistograms->FillHistogram("ESD_Background_ZR", tmpY, radiusBG);
	  fHistograms->FillHistogram("ESD_Background_XY", tmpX, tmpY);
	  fHistograms->FillHistogram("Background_InvMass_vs_Pt_Spectra",massBG,sqrt(backgroundCandidate->GetPx()*backgroundCandidate->GetPx()+backgroundCandidate->GetPy()*backgroundCandidate->GetPy()));
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
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer = new TList();
  }
  
  //Adding the histograms to the output container
  fHistograms->GetOutputContainer(fOutputContainer);

  
  if(fWriteNtuple){
    if(fGammaNtuple == NULL){
      fGammaNtuple = new TNtuple("V0ntuple","V0ntuple","OnTheFly:HasVertex:NegPIDProb:PosPIDProb:X:Y:Z:R:MotherCandidateNDF:MotherCandidateChi2:MotherCandidateEnergy:MotherCandidateEta:MotherCandidatePt:MotherCandidateMass:MotherCandidateWidth:MCMotherCandidatePT:EPOpeningAngle:ElectronEnergy:ElectronPt:ElectronEta:ElectronPhi:PositronEnergy:PositronPt:PositronEta:PositronPhi:HasSameMCMother:MotherMCParticlePIDCode",50000);
    }
    if(fNeutralMesonNtuple == NULL){
      fNeutralMesonNtuple = new TNtuple("NeutralMesonNtuple","NeutralMesonNtuple","test");
    }
    TList * ntupleTList = new TList();
    ntupleTList->SetName("Ntuple");
    ntupleTList->Add((TNtuple*)fGammaNtuple);
    fOutputContainer->Add(ntupleTList);
  }

  fOutputContainer->SetName(GetName());
}

Double_t AliAnalysisTaskGammaConversion::GetMCOpeningAngle(TParticle* daughter0, TParticle* daughter1) const{
  //helper function
  TVector3 v3D0(daughter0->Px(),daughter0->Py(),daughter0->Pz());
  TVector3 v3D1(daughter1->Px(),daughter1->Py(),daughter1->Pz());
  return v3D0.Angle(v3D1);
}
