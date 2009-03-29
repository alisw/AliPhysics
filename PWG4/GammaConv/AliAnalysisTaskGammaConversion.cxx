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
  fIsTrueReconstructedGammas(),
  electronv1(),
  electronv2(),
  fElectronMass(-1),
  fGammaMass(-1),
  fPi0Mass(-1),
  fEtaMass(-1),
  fGammaWidth(-1),
  fPi0Width(-1),
  fEtaWidth(-1),
  fMinOpeningAngleGhostCut(0.),
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
  fIsTrueReconstructedGammas(),
  electronv1(),
  electronv2(),
  fElectronMass(-1),
  fGammaMass(-1),
  fPi0Mass(-1),
  fEtaMass(-1),
  fGammaWidth(-1),
  fPi0Width(-1),
  fEtaWidth(-1),
  fMinOpeningAngleGhostCut(0.),
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
  fIsTrueReconstructedGammas.clear();
  electronv1.clear();
  electronv2.clear();
	
  //Clear the data in the v0Reader
  fV0Reader->UpdateEventByEventData();

  
  // Process the MC information
  if(fDoMCTruth){
    ProcessMCData();
  }
  
  //Process the v0 information with no cuts
  ProcessV0sNoCut();
  
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

  if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
    return; // aborts if the primary vertex does not have contributors.
  }

  for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
    TParticle* particle = (TParticle *)fStack->Particle(iTracks);
		
    if (!particle) {
      //print warning here
      continue;
    }

    if(TMath::Abs(particle->Eta())> fV0Reader->GetEtaCut() )	continue;
					
    if(particle->R()>fV0Reader->GetMaxRCut())	continue; // cuts on distance from collision point
		
    Double_t tmpPhi=particle->Phi();
		
    if(particle->Phi()> TMath::Pi()){
      tmpPhi = particle->Phi()-(2*TMath::Pi());
    }
		
    Double_t rapidity;
    if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
      rapidity=0;
    }
    else{
      rapidity = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())));
    }	
		
    //process the gammas
    if (particle->GetPdgCode() == 22){
			
      if(particle->GetMother(0) >-1 && fStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
	continue; // no photon as mothers!
      }

      if(particle->GetMother(0) >= fStack->GetNprimary()){
	continue; // the gamma has a mother, and it is not a primary particle
      }

      fMCAllGammas.push_back(particle);
			
      fHistograms->FillHistogram("MC_allGamma_Energy", particle->Energy());
      fHistograms->FillHistogram("MC_allGamma_Pt", particle->Pt());
      fHistograms->FillHistogram("MC_allGamma_Eta", particle->Eta());
      fHistograms->FillHistogram("MC_allGamma_Phi", tmpPhi);
      fHistograms->FillHistogram("MC_allGamma_Rapid", rapidity);
			
			
      if(particle->GetMother(0) < 0){   // direct gamma
	fHistograms->FillHistogram("MC_allDirectGamma_Energy",particle->Energy());
	fHistograms->FillHistogram("MC_allDirectGamma_Pt", particle->Pt());
	fHistograms->FillHistogram("MC_allDirectGamma_Eta", particle->Eta());
	fHistograms->FillHistogram("MC_allDirectGamma_Phi", tmpPhi);
	fHistograms->FillHistogram("MC_allDirectGamma_Rapid", rapidity);				
      }
			
			
      // looking for conversion (electron + positron from pairbuilding (= 5) )
      TParticle* ePos = NULL;
      TParticle* eNeg = NULL;
			
      if(particle->GetNDaughters() >= 2){
	for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	  TParticle *tmpDaughter = fStack->Particle(daughterIndex);
	  if(tmpDaughter->GetUniqueID() == 5){
	    if(tmpDaughter->GetPdgCode() == 11){
	      eNeg = tmpDaughter;
	    }
	    else if(tmpDaughter->GetPdgCode() == -11){
	      ePos = tmpDaughter;
	    }
	  }
	}
      }
			
			
      if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
	continue;
      }
			
			
      Double_t ePosPhi = ePos->Phi();
      if(ePos->Phi()> TMath::Pi()) ePosPhi = ePos->Phi()-(2*TMath::Pi());
			
      Double_t eNegPhi = eNeg->Phi();
      if(eNeg->Phi()> TMath::Pi()) eNegPhi = eNeg->Phi()-(2*TMath::Pi());
			
			
      if(ePos->Pt()<fV0Reader->GetPtCut() || eNeg->Pt()<fV0Reader->GetPtCut()){
	continue; // no reconstruction below the Pt cut
      }
					
      if(TMath::Abs(ePos->Eta())> fV0Reader->GetEtaCut() || TMath::Abs(eNeg->Eta())> fV0Reader->GetEtaCut()){
	continue;
      }	
				
      if(ePos->R()>fV0Reader->GetMaxRCut()){
	continue; // cuts on distance from collision point
      }
      
      
      if((TMath::Abs(ePos->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue()  > ePos->R()){
	continue;               // line cut to exclude regions where we do not reconstruct
      }		
      		
      fHistograms->FillHistogram("MC_ConvGamma_Energy", particle->Energy());
      fHistograms->FillHistogram("MC_ConvGamma_Pt", particle->Pt());
      fHistograms->FillHistogram("MC_ConvGamma_Eta", particle->Eta());
      fHistograms->FillHistogram("MC_ConvGamma_Phi", tmpPhi);
      fHistograms->FillHistogram("MC_ConvGamma_Rapid", rapidity);
      fHistograms->FillHistogram("MC_ConvGamma_Pt_Eta", particle->Pt(),particle->Eta());
			
      fHistograms->FillHistogram("MC_E_Energy", eNeg->Energy());
      fHistograms->FillHistogram("MC_E_Pt", eNeg->Pt());
      fHistograms->FillHistogram("MC_E_Eta", eNeg->Eta());
      fHistograms->FillHistogram("MC_E_Phi", eNegPhi);
			
      fHistograms->FillHistogram("MC_P_Energy", ePos->Energy());
      fHistograms->FillHistogram("MC_P_Pt", ePos->Pt());
      fHistograms->FillHistogram("MC_P_Eta", ePos->Eta());
      fHistograms->FillHistogram("MC_P_Phi", ePosPhi);
			
			
			
      //cout << "filled histos for converted gamma, ePos, eNeg" << endl;
			
      // begin Mapping 
      Int_t rBin    = fHistograms->GetRBin(ePos->R());
      Int_t phiBin  = fHistograms->GetPhiBin(particle->Phi());
			
      TString nameMCMappingPhiR="";
      nameMCMappingPhiR.Form("MC_Conversion_Mapping-Phi%02d-R%02d",phiBin,rBin);
      fHistograms->FillHistogram(nameMCMappingPhiR, ePos->Vz(), particle->Eta());
			
      TString nameMCMappingPhi="";
      nameMCMappingPhi.Form("MC_Conversion_Mapping-Phi%02d",phiBin);
      fHistograms->FillHistogram(nameMCMappingPhi, particle->Eta());
			
      TString nameMCMappingR="";
      nameMCMappingR.Form("MC_Conversion_Mapping-R%02d",rBin);
      fHistograms->FillHistogram(nameMCMappingR, particle->Eta());
			
      TString nameMCMappingPhiInR="";
      nameMCMappingPhiInR.Form("MC_Conversion_Mapping_Phi_R-%02d",rBin);
      fHistograms->FillHistogram(nameMCMappingPhiInR, tmpPhi);
      //end mapping
			
      fHistograms->FillHistogram("MC_Conversion_R",ePos->R());
      fHistograms->FillHistogram("MC_Conversion_ZR",ePos->Vz(),ePos->R());
      fHistograms->FillHistogram("MC_Conversion_XY",ePos->Vx(),ePos->Vy());
      fHistograms->FillHistogram("MC_Conversion_OpeningAngle",GetMCOpeningAngle(ePos, eNeg));
			
      //cout << "mapping is done" << endl;
			
			
      if(particle->GetMother(0) < 0){ // no mother = direct gamma, still inside converted
	fHistograms->FillHistogram("MC_ConvDirectGamma_Energy",particle->Energy());
	fHistograms->FillHistogram("MC_ConvDirectGamma_Pt", particle->Pt());
	fHistograms->FillHistogram("MC_ConvDirectGamma_Eta", particle->Eta());
	fHistograms->FillHistogram("MC_ConvDirectGamma_Phi", tmpPhi);
	fHistograms->FillHistogram("MC_ConvDirectGamma_Rapid", rapidity);
				
      } // end direct gamma
      else{   // mother exits 
	if( fStack->Particle(particle->GetMother(0))->GetPdgCode()==10441 ||//chic0 
	    fStack->Particle(particle->GetMother(0))->GetPdgCode()==20443 ||//psi2S
	    fStack->Particle(particle->GetMother(0))->GetPdgCode()==445  //chic2
	    ){ 
	  fMCGammaChic.push_back(particle);
	}
      }  // end if mother exits
    } // end if particle is a photon
		
    if(particle->GetNDaughters() == 2){
			
      TParticle* daughter0 = (TParticle*)fStack->Particle(particle->GetFirstDaughter());
      TParticle* daughter1 = (TParticle*)fStack->Particle(particle->GetLastDaughter());
			
      if(daughter0->GetPdgCode() != 22 || daughter1->GetPdgCode() != 22) continue; //check for gamma gamma daughters
			
			
			
      // check for conversions now -> have to pass eta and line cut!
      Bool_t daughter0Electron = kFALSE;
      Bool_t daughter0Positron = kFALSE;
      Bool_t daughter1Electron = kFALSE;
      Bool_t daughter1Positron = kFALSE;
			
			
			
      if(daughter0->GetNDaughters() >= 2){
	for(Int_t TrackIndex=daughter0->GetFirstDaughter();TrackIndex<=daughter0->GetLastDaughter();TrackIndex++){
	  TParticle *tmpDaughter = fStack->Particle(TrackIndex);
	  if(tmpDaughter->GetUniqueID() == 5){
	    if(tmpDaughter->GetPdgCode() == 11){
	      if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() ){
		if( ( TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue()  < tmpDaughter->R() ){
		  daughter0Electron = kTRUE;
		}
								
	      }
	    }
	    else if(tmpDaughter->GetPdgCode() == -11){
	      if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() ){
		if( ( TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue()  < tmpDaughter->R() ){
		  daughter0Positron = kTRUE;
		}
								
	      }
							
	    }
	  }
	}
      }
			
			
	
      if(daughter1->GetNDaughters() >= 2){
	for(Int_t TrackIndex=daughter1->GetFirstDaughter();TrackIndex<=daughter1->GetLastDaughter();TrackIndex++){
	  TParticle *tmpDaughter = fStack->Particle(TrackIndex);
	  if(tmpDaughter->GetUniqueID() == 5){
	    if(tmpDaughter->GetPdgCode() == 11){
	      if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() ){
		if( ( TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue()  < tmpDaughter->R() ){
		  daughter1Electron = kTRUE;
		}
								
	      }
	    }
	    else if(tmpDaughter->GetPdgCode() == -11){
	      if( TMath::Abs(tmpDaughter->Eta()) <= fV0Reader->GetEtaCut() ){
		if( ( TMath::Abs(tmpDaughter->Vz()) * fV0Reader->GetLineCutZRSlope()) - fV0Reader->GetLineCutZValue()  < tmpDaughter->R() ){
		  daughter1Positron = kTRUE;
		}
								
	      }
							
	    }
	  }
	}
      }
			
												
			
			
      if(particle->GetPdgCode()==111){     //Pi0
	if( iTracks >= fStack->GetNprimary()){
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_Eta", particle->Eta());
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_Rapid", rapidity);
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_Phi", tmpPhi);
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt", particle->Pt());
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_Energy", particle->Energy());
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_R", particle->R());
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_ZR", particle->Vz(),particle->R());
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	  fHistograms->FillHistogram("MC_Pi0_Secondaries_XY", particle->Vx(),particle->Vy());//only fill from one daughter to avoid multiple filling
					
	  if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
	    fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
	    if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
	      fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
	      fHistograms->FillHistogram("MC_Pi0_Secondaries_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
	    }
	  }
	}
	else{
	  fHistograms->FillHistogram("MC_Pi0_Eta", particle->Eta());	
	  fHistograms->FillHistogram("MC_Pi0_Rapid", rapidity);
	  fHistograms->FillHistogram("MC_Pi0_Phi", tmpPhi);
	  fHistograms->FillHistogram("MC_Pi0_Pt", particle->Pt());
	  fHistograms->FillHistogram("MC_Pi0_Energy", particle->Energy());
	  fHistograms->FillHistogram("MC_Pi0_R", particle->R());
	  fHistograms->FillHistogram("MC_Pi0_ZR", particle->Vz(),particle->R());
	  fHistograms->FillHistogram("MC_Pi0_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	  fHistograms->FillHistogram("MC_Pi0_XY", particle->Vx(), particle->Vy());//only fill from one daughter to avoid multiple filling
					
	  if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
	    fHistograms->FillHistogram("MC_Pi0_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
	    fHistograms->FillHistogram("MC_Pi0_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
	    if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
	      fHistograms->FillHistogram("MC_Pi0_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
	      fHistograms->FillHistogram("MC_Pi0_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
	    }
	  }
	}
      }
			
      if(particle->GetPdgCode()==221){   //Eta
	fHistograms->FillHistogram("MC_Eta_Eta", particle->Eta());
	fHistograms->FillHistogram("MC_Eta_Rapid", rapidity);
	fHistograms->FillHistogram("MC_Eta_Phi",tmpPhi);
	fHistograms->FillHistogram("MC_Eta_Pt", particle->Pt());
	fHistograms->FillHistogram("MC_Eta_Energy", particle->Energy());
	fHistograms->FillHistogram("MC_Eta_R", particle->R());
	fHistograms->FillHistogram("MC_Eta_ZR", particle->Vz(),particle->R());
	fHistograms->FillHistogram("MC_Eta_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
	fHistograms->FillHistogram("MC_Eta_XY", particle->Vx(), particle->Vy());//only fill from one daughter to avoid multiple filling
				
	if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
	  fHistograms->FillHistogram("MC_Eta_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
	  fHistograms->FillHistogram("MC_Eta_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
	  if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
	    fHistograms->FillHistogram("MC_Eta_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
	    fHistograms->FillHistogram("MC_Eta_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
	  }
					
	}
				
      }
			
      // all motherparticles with 2 gammas as daughters
      fHistograms->FillHistogram("MC_Mother_R", particle->R());
      fHistograms->FillHistogram("MC_Mother_ZR", particle->Vz(),particle->R());
      fHistograms->FillHistogram("MC_Mother_XY", particle->Vx(),particle->Vy());
      fHistograms->FillHistogram("MC_Mother_Mass", particle->GetCalcMass());
      fHistograms->FillHistogram("MC_Mother_GammaDaughter_OpeningAngle", GetMCOpeningAngle(daughter0,daughter1));
      fHistograms->FillHistogram("MC_Mother_Energy", particle->Energy());
      fHistograms->FillHistogram("MC_Mother_Pt", particle->Pt());
      fHistograms->FillHistogram("MC_Mother_Eta", particle->Eta());
      fHistograms->FillHistogram("MC_Mother_Rapid", rapidity);
      fHistograms->FillHistogram("MC_Mother_Phi",tmpPhi);
      fHistograms->FillHistogram("MC_Mother_InvMass_vs_Pt",particle->GetMass(),particle->Pt());			
      if(TMath::Abs(daughter0->Eta()) <= fV0Reader->GetEtaCut() && TMath::Abs(daughter1->Eta()) <= fV0Reader->GetEtaCut() ){
	fHistograms->FillHistogram("MC_Mother_Pt_Eta_withinAcceptance", particle->Pt(),particle->Eta());
	fHistograms->FillHistogram("MC_Mother_Pt_Rapid_withinAcceptance", particle->Pt(),rapidity);
	fHistograms->FillHistogram("MC_Mother_InvMass_vs_Pt_withinAcceptance",particle->GetMass(),particle->Pt());			
	if(daughter0Electron && daughter0Positron && daughter1Electron && daughter1Positron){
	  fHistograms->FillHistogram("MC_Mother_Pt_Eta_ConvGamma_withinAcceptance", particle->Pt(),particle->Eta());
	  fHistograms->FillHistogram("MC_Mother_Pt_Rapid_ConvGamma_withinAcceptance", particle->Pt(),rapidity);
	  fHistograms->FillHistogram("MC_Mother_InvMass_vs_Pt_ConvGamma_withinAcceptance",particle->GetMass(),particle->Pt());			

	}
				
				
      }
			
      //cout << "mother histos are filled" << endl;
			
    } // end if(particle->GetNDaughters() == 2)
		
  }// end for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++)
	
  //cout << "right before the end of processMCdata" << endl;
	
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

void AliAnalysisTaskGammaConversion::ProcessV0sNoCut(){

  Int_t numberOfV0s = fV0Reader->GetNumberOfV0s();
  for(Int_t i=0;i<numberOfV0s;i++){
    /*AliESDv0 * cV0 = */fV0Reader->GetV0(i);

    if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
      return;
    }
    
    if(fDoMCTruth){
      
      if(fV0Reader->HasSameMCMother() == kFALSE){
	continue;
      }
		
      TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
      TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();

      if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
	continue;
      }
      if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
	continue;
      }
	
      if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
      
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Pt", fV0Reader->GetMotherCandidatePt());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Energy", fV0Reader->GetMotherCandidateEnergy());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Eta", fV0Reader->GetMotherCandidateEta());				
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Phi", fV0Reader->GetMotherCandidatePhi());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Mass", fV0Reader->GetMotherCandidateMass());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Width", fV0Reader->GetMotherCandidateWidth());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Chi2", fV0Reader->GetMotherCandidateChi2());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_NDF", fV0Reader->GetMotherCandidateNDF());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Rapid", fV0Reader->GetMotherCandidateRapidity());
	fHistograms->FillHistogram("ESD_NoCutConvGamma_Pt_Eta", fV0Reader->GetMotherCandidatePt(),fV0Reader->GetMotherCandidateEta());
	
	fHistograms->FillHistogram("ESD_NoCutConversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
	fHistograms->FillHistogram("ESD_NoCutConversion_R", fV0Reader->GetXYRadius());
	fHistograms->FillHistogram("ESD_NoCutConversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
	fHistograms->FillHistogram("ESD_NoCutConversion_OpeningAngle", fV0Reader->GetOpeningAngle());
	
	/*
	  ESD_NoCutConvGamma_Pt
	  ESD_NoCutConvGamma_Energy
	  ESD_NoCutConvGamma_Eta
	  ESD_NoCutConvGamma_Phi
	  ESD_NoCutConvGamma_Mass
	  ESD_NoCutConvGamma_Width
	  ESD_NoCutConvGamma_Chi2
	  ESD_NoCutConvGamma_NDF
	  ESD_NoCutConvGamma_PtvsEta
	  ESD_NoCutConversion_XY
	  ESD_NoCutConversion_R
	  ESD_NoCutConversion_ZR
	  ESD_NoCutConversion_OpeningAngle
	*/
      }
    }
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
    fHistograms->FillHistogram("ESD_Conversion_R", fV0Reader->GetXYRadius());
    fHistograms->FillHistogram("ESD_Conversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
    fHistograms->FillHistogram("ESD_Conversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
    fHistograms->FillHistogram("ESD_Conversion_OpeningAngle", fV0Reader->GetOpeningAngle());    
		
    fHistograms->FillHistogram("ESD_E_Energy", fV0Reader->GetNegativeTrackEnergy());
    fHistograms->FillHistogram("ESD_E_Pt", fV0Reader->GetNegativeTrackPt());
    fHistograms->FillHistogram("ESD_E_Eta", fV0Reader->GetNegativeTrackEta());
    fHistograms->FillHistogram("ESD_E_Phi", fV0Reader->GetNegativeTrackPhi());
		
    fHistograms->FillHistogram("ESD_P_Energy", fV0Reader->GetPositiveTrackEnergy());
    fHistograms->FillHistogram("ESD_P_Pt", fV0Reader->GetPositiveTrackPt());
    fHistograms->FillHistogram("ESD_P_Eta", fV0Reader->GetPositiveTrackEta());
    fHistograms->FillHistogram("ESD_P_Phi", fV0Reader->GetPositiveTrackPhi());
		
    fHistograms->FillHistogram("ESD_ConvGamma_Energy", fV0Reader->GetMotherCandidateEnergy());
    fHistograms->FillHistogram("ESD_ConvGamma_Pt", fV0Reader->GetMotherCandidatePt());
    fHistograms->FillHistogram("ESD_ConvGamma_Eta", fV0Reader->GetMotherCandidateEta());
    fHistograms->FillHistogram("ESD_ConvGamma_Phi", fV0Reader->GetMotherCandidatePhi());
    fHistograms->FillHistogram("ESD_ConvGamma_Mass", fV0Reader->GetMotherCandidateMass());
    fHistograms->FillHistogram("ESD_ConvGamma_Width", fV0Reader->GetMotherCandidateWidth());
    fHistograms->FillHistogram("ESD_ConvGamma_Chi2", fV0Reader->GetMotherCandidateChi2());
    fHistograms->FillHistogram("ESD_ConvGamma_NDF", fV0Reader->GetMotherCandidateNDF());
    fHistograms->FillHistogram("ESD_ConvGamma_Rapid", fV0Reader->GetMotherCandidateRapidity());
    fHistograms->FillHistogram("ESD_ConvGamma_Pt_Eta", fV0Reader->GetMotherCandidatePt(),fV0Reader->GetMotherCandidateEta());
		
		
    // begin mapping
    Int_t rBin    = fHistograms->GetRBin(fV0Reader->GetXYRadius());
    Int_t phiBin  = fHistograms->GetPhiBin(fV0Reader->GetNegativeTrackPhi());
    Double_t motherCandidateEta= fV0Reader->GetMotherCandidateEta();
		
    TString nameESDMappingPhiR="";
    nameESDMappingPhiR.Form("ESD_Conversion_Mapping-Phi%02d-R%02d",phiBin,rBin);
    fHistograms->FillHistogram(nameESDMappingPhiR, fV0Reader->GetZ(), motherCandidateEta);
		
    TString nameESDMappingPhi="";
    nameESDMappingPhi.Form("ESD_Conversion_Mapping-Phi%02d",phiBin);
    fHistograms->FillHistogram(nameESDMappingPhi, fV0Reader->GetZ(), motherCandidateEta);
		
    TString nameESDMappingR="";
    nameESDMappingR.Form("ESD_Conversion_Mapping-R%02d",rBin);
    fHistograms->FillHistogram(nameESDMappingR, fV0Reader->GetZ(), motherCandidateEta);  
		
    TString nameESDMappingPhiInR="";
    nameESDMappingPhiInR.Form("ESD_Conversion_Mapping_Phi_R-%02d",rBin);
    fHistograms->FillHistogram(nameESDMappingPhiInR, fV0Reader->GetMotherCandidatePhi());
    // end mapping
		
    fKFReconstructedGammas.push_back(*fV0Reader->GetMotherCandidateKFCombination());
    electronv1.push_back(fV0Reader->GetCurrentV0()->GetPindex());
    electronv2.push_back(fV0Reader->GetCurrentV0()->GetNindex());

		
    //----------------------------------- checking for "real" conversions (MC match) --------------------------------------
    if(fDoMCTruth){
			
      if(fV0Reader->HasSameMCMother() == kFALSE){
	fIsTrueReconstructedGammas.push_back(kFALSE);
	continue;
      }
      
		
      TParticle * negativeMC = (TParticle*)fV0Reader->GetNegativeMCParticle();
      TParticle * positiveMC = (TParticle*)fV0Reader->GetPositiveMCParticle();

      if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
	fIsTrueReconstructedGammas.push_back(kFALSE);
	continue;
      }
      if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
	fIsTrueReconstructedGammas.push_back(kFALSE);
	continue;
      }
	
      if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22){
	fIsTrueReconstructedGammas.push_back(kTRUE);
				
	fHistograms->FillHistogram("ESD_TrueConvGamma_Pt", fV0Reader->GetMotherCandidatePt());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Energy", fV0Reader->GetMotherCandidateEnergy());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Eta", fV0Reader->GetMotherCandidateEta());				
	fHistograms->FillHistogram("ESD_TrueConvGamma_Phi", fV0Reader->GetMotherCandidatePhi());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Mass", fV0Reader->GetMotherCandidateMass());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Width", fV0Reader->GetMotherCandidateWidth());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Chi2", fV0Reader->GetMotherCandidateChi2());
	fHistograms->FillHistogram("ESD_TrueConvGamma_NDF", fV0Reader->GetMotherCandidateNDF());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Pt_Eta", fV0Reader->GetMotherCandidatePt(),fV0Reader->GetMotherCandidateEta());
	fHistograms->FillHistogram("ESD_TrueConvGamma_Rapid", fV0Reader->GetMotherCandidateRapidity());
	fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLength", /*fV0Reader->GetNegativeTrackLength()*/fV0Reader->GetNegativeNTPCClusters());
	fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLength", /*fV0Reader->GetPositiveTrackLength()*/fV0Reader->GetPositiveNTPCClusters());
	fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLengthVSInvMass",/*fV0Reader->GetNegativeTrackLength()*/fV0Reader->GetNegativeNTPCClusters(),fV0Reader->GetMotherCandidateMass());
	fHistograms->FillHistogram("ESD_TrueConvGamma_TrackLengthVSInvMass",/*fV0Reader->GetPositiveTrackLength()*/fV0Reader->GetPositiveNTPCClusters(),fV0Reader->GetMotherCandidateMass());
	
	fHistograms->FillHistogram("ESD_TrueConversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
	fHistograms->FillHistogram("ESD_TrueConversion_R", fV0Reader->GetXYRadius());
	fHistograms->FillHistogram("ESD_TrueConversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
	fHistograms->FillHistogram("ESD_TrueConversion_OpeningAngle", fV0Reader->GetOpeningAngle());

	
	/*
	  fHistograms->FillHistogram("ESD_TrueConversion_XY", fV0Reader->GetX(),fV0Reader->GetY());
	  fHistograms->FillHistogram("ESD_TrueConversion_R", fV0Reader->GetXYRadius());
	  fHistograms->FillHistogram("ESD_TrueConversion_ZR", fV0Reader->GetZ(),fV0Reader->GetXYRadius());
	  fHistograms->FillHistogram("ESD_TrueConversion_OpeningAngle", fV0Reader->GetOpeningAngle());
	*/



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
      }//if(fV0Reader->GetMotherMCParticle()->GetPdgCode() == 22)
      else{
	fIsTrueReconstructedGammas.push_back(kFALSE);
      }
    }//if(fDoMCTruth)
  }//while(fV0Reader->NextV0)
  fHistograms->FillHistogram("ESD_NumberOfSurvivingV0s", nSurvivingV0s);
  fHistograms->FillHistogram("ESD_NumberOfV0s", fV0Reader->GetNumberOfV0s());
	
  //cout << "nearly at the end of doMCTruth" << endl;
	
}

void AliAnalysisTaskGammaConversion::ProcessGammasForNeutralMesonAnalysis(){
  // see header file for documentation
	
  for(UInt_t firstGammaIndex=0;firstGammaIndex<fKFReconstructedGammas.size();firstGammaIndex++){
    for(UInt_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fKFReconstructedGammas.size();secondGammaIndex++){
			
      AliKFParticle * twoGammaDecayCandidateDaughter0 = &fKFReconstructedGammas[firstGammaIndex];
      AliKFParticle * twoGammaDecayCandidateDaughter1 = &fKFReconstructedGammas[secondGammaIndex];
      
      if(electronv1[firstGammaIndex]==electronv1[secondGammaIndex] || electronv1[firstGammaIndex]==electronv2[secondGammaIndex]){
	continue;
      }
      if(electronv2[firstGammaIndex]==electronv1[secondGammaIndex] || electronv2[firstGammaIndex]==electronv2[secondGammaIndex]){
	continue;
      }

      /*
	if(fIsTrueReconstructedGammas[firstGammaIndex] == kFALSE || fIsTrueReconstructedGammas[secondGammaIndex] == kFALSE){
	continue;
	}
      */

      AliKFParticle *twoGammaCandidate = new AliKFParticle(*twoGammaDecayCandidateDaughter0,*twoGammaDecayCandidateDaughter1);
			
      Double_t massTwoGammaCandidate = 0.;
      Double_t widthTwoGammaCandidate = 0.;
      Double_t chi2TwoGammaCandidate =10000.;	
      twoGammaCandidate->GetMass(massTwoGammaCandidate,widthTwoGammaCandidate);
      if(twoGammaCandidate->GetNDF()>0){
	chi2TwoGammaCandidate = twoGammaCandidate->GetChi2()/twoGammaCandidate->GetNDF();
				
	if(chi2TwoGammaCandidate>0 && chi2TwoGammaCandidate<fV0Reader->GetChi2CutMeson()){					
					
	  TVector3 MomentumVectorTwoGammaCandidate(twoGammaCandidate->GetPx(),twoGammaCandidate->GetPy(),twoGammaCandidate->GetPz());
	  TVector3 SpaceVectorTwoGammaCandidate(twoGammaCandidate->GetX(),twoGammaCandidate->GetY(),twoGammaCandidate->GetZ());
					
	  Double_t openingAngleTwoGammaCandidate = twoGammaDecayCandidateDaughter0->GetAngle(*twoGammaDecayCandidateDaughter1);					
	  Double_t rapidity;
	  if(twoGammaCandidate->GetE() - twoGammaCandidate->GetPz() == 0 || twoGammaCandidate->GetE() + twoGammaCandidate->GetPz() == 0){
	    rapidity=0;
	  }
	  else{
	    rapidity = 0.5*(TMath::Log((twoGammaCandidate->GetE() +twoGammaCandidate->GetPz()) / (twoGammaCandidate->GetE()-twoGammaCandidate->GetPz())));
	  }
					
	  if(openingAngleTwoGammaCandidate < fMinOpeningAngleGhostCut) continue;   // minimum opening angle to avoid using ghosttracks

	  fHistograms->FillHistogram("ESD_Mother_GammaDaughter_OpeningAngle", openingAngleTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_Mother_Energy", twoGammaCandidate->GetE());
	  fHistograms->FillHistogram("ESD_Mother_Pt", MomentumVectorTwoGammaCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Mother_Eta", MomentumVectorTwoGammaCandidate.Eta());
	  fHistograms->FillHistogram("ESD_Mother_Rapid", rapidity);					
	  fHistograms->FillHistogram("ESD_Mother_Phi", SpaceVectorTwoGammaCandidate.Phi());
	  fHistograms->FillHistogram("ESD_Mother_Mass", massTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_Mother_R", SpaceVectorTwoGammaCandidate.Pt());    // Pt in Space == R!!!
	  fHistograms->FillHistogram("ESD_Mother_ZR", twoGammaCandidate->GetZ(), SpaceVectorTwoGammaCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Mother_XY", twoGammaCandidate->GetX(), twoGammaCandidate->GetY());
	  fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt",massTwoGammaCandidate ,MomentumVectorTwoGammaCandidate.Pt());
	}
      }
      delete twoGammaCandidate;
			
      //cout << "nearly at the end of processgamma for neutral meson ..." << endl;
			
			
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
					
	  TVector3 MomentumVectorbackgroundCandidate(backgroundCandidate->GetPx(),backgroundCandidate->GetPy(),backgroundCandidate->GetPz());
	  TVector3 SpaceVectorbackgroundCandidate(backgroundCandidate->GetX(),backgroundCandidate->GetY(),backgroundCandidate->GetZ());
					
	  Double_t openingAngleBG = currentEventGoodV0->GetAngle(*previousGoodV0);

	  Double_t rapidity;
	  if(backgroundCandidate->GetE() - backgroundCandidate->GetPz() == 0 || backgroundCandidate->GetE() + backgroundCandidate->GetPz() == 0) rapidity=0;
	  else rapidity = 0.5*(TMath::Log((backgroundCandidate->GetE() +backgroundCandidate->GetPz()) / (backgroundCandidate->GetE()-backgroundCandidate->GetPz())));

					
					
					
	  if(openingAngleBG < fMinOpeningAngleGhostCut ) continue;   // minimum opening angle to avoid using ghosttracks
			
					
	  fHistograms->FillHistogram("ESD_Background_GammaDaughter_OpeningAngle", openingAngleBG);
	  fHistograms->FillHistogram("ESD_Background_Energy", backgroundCandidate->GetE());
	  fHistograms->FillHistogram("ESD_Background_Pt",  MomentumVectorbackgroundCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Background_Eta", MomentumVectorbackgroundCandidate.Eta());
	  fHistograms->FillHistogram("ESD_Background_Rapidity", rapidity);
	  fHistograms->FillHistogram("ESD_Background_Phi", SpaceVectorbackgroundCandidate.Phi());
	  fHistograms->FillHistogram("ESD_Background_Mass", massBG);
	  fHistograms->FillHistogram("ESD_Background_R", SpaceVectorbackgroundCandidate.Pt());  // Pt in Space == R!!!!
	  fHistograms->FillHistogram("ESD_Background_ZR", backgroundCandidate->GetZ(), SpaceVectorbackgroundCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Background_XY", backgroundCandidate->GetX(), backgroundCandidate->GetY());
	  fHistograms->FillHistogram("ESD_Background_InvMass_vs_Pt",massBG,MomentumVectorbackgroundCandidate.Pt());
	}
      }
      delete backgroundCandidate;   
      //cout << "nearly at the end of background" << endl;
			
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
