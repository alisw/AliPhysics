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
#include "TNtuple.h"

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
  fESDEvent(NULL),	
  fOutputContainer(NULL),
  fHistograms(NULL),
  fDoMCTruth(kFALSE),
  fMCAllGammas(),
  fMCPi0s(),
  fMCEtas(),
  fMCGammaChic(),
  fKFReconstructedGammas(),
  fIsTrueReconstructedGammas(),
  fElectronv1(),
  fElectronv2(),
  fCurrentEventPosElectron(),
  fPreviousEventPosElectron(),
  fCurrentEventNegElectron(),
  fPreviousEventNegElectron(),
  fKFReconstructedGammasCut(),		 	
  fPreviousEventTLVNegElectron(),
  fPreviousEventTLVPosElectron(),	
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
  fESDEvent(NULL),	
  fOutputContainer(0x0),
  fHistograms(NULL),
  fDoMCTruth(kFALSE),
  fMCAllGammas(),
  fMCPi0s(),
  fMCEtas(),
  fMCGammaChic(),
  fKFReconstructedGammas(),
  fIsTrueReconstructedGammas(),
  fElectronv1(),
  fElectronv2(),
  fCurrentEventPosElectron(),
  fPreviousEventPosElectron(),
  fCurrentEventNegElectron(),
  fPreviousEventNegElectron(),
  fKFReconstructedGammasCut(),	
  fPreviousEventTLVNegElectron(),
  fPreviousEventTLVPosElectron(),
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
  fElectronv1.clear();
  fElectronv2.clear();
  fCurrentEventPosElectron.clear();
  fCurrentEventNegElectron.clear();	
  fKFReconstructedGammasCut.clear(); 
	
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

  CheckV0Efficiency();
  
  //Process reconstructed gammas electrons for Chi_c Analysis
  ProcessGammaElectronsForChicAnalysis();

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

    ///////////////////////Begin Chic Analysis/////////////////////////////


    if(particle->GetPdgCode() == 443){//Is JPsi

      if(particle->GetNDaughters()==2){
	if(TMath::Abs(fStack->Particle(particle->GetFirstDaughter())->GetPdgCode()) == 11 &&
	   TMath::Abs(fStack->Particle(particle->GetLastDaughter())->GetPdgCode()) == 11){
	  TParticle* daug0 = fStack->Particle(particle->GetFirstDaughter());
	  TParticle* daug1 = fStack->Particle(particle->GetLastDaughter());
	  if(TMath::Abs(daug0->Eta()) < 0.9 && TMath::Abs(daug1->Eta()) < 0.9)
	    fHistograms->FillTable("Table_Electrons",3);//e+ e-  from J/Psi inside acceptance

	  if( TMath::Abs(daug0->Eta()) < 0.9){
	    if(daug0->GetPdgCode() == -11)
	      fHistograms->FillTable("Table_Electrons",1);//e+  from J/Psi inside acceptance
	    else
	      fHistograms->FillTable("Table_Electrons",2);//e-   from J/Psi inside acceptance

	  }
	  if(TMath::Abs(daug1->Eta()) < 0.9){
	    if(daug1->GetPdgCode() == -11)
	      fHistograms->FillTable("Table_Electrons",1);//e+  from J/Psi inside acceptance
	    else
	      fHistograms->FillTable("Table_Electrons",2);//e-   from J/Psi inside acceptance
	  }
	}
      }
    }
    //              const int CHI_C0   = 10441;
    //              const int CHI_C1   = 20443;
    //              const int CHI_C2   = 445
    if(particle->GetPdgCode() == 22){//gamma from JPsi
      if(particle->GetMother(0) > -1){
	if(fStack->Particle(particle->GetMother(0))->GetPdgCode() == 10441 ||
	   fStack->Particle(particle->GetMother(0))->GetPdgCode() == 20443 ||
	   fStack->Particle(particle->GetMother(0))->GetPdgCode() == 445){
	  if(TMath::Abs(particle->Eta()) < 1.2)
	    fHistograms->FillTable("Table_Electrons",17);// gamma from chic inside accptance
	}
      }
    }
    if(particle->GetPdgCode() == 10441 || particle->GetPdgCode() == 20443 || particle->GetPdgCode() == 445){
      if( particle->GetNDaughters() == 2){
	TParticle* daug0 = fStack->Particle(particle->GetFirstDaughter());
	TParticle* daug1 = fStack->Particle(particle->GetLastDaughter());

	if( (daug0->GetPdgCode() == 443 || daug0->GetPdgCode() == 22) && (daug1->GetPdgCode() == 443 || daug1->GetPdgCode() == 22) ){
	  if( daug0->GetPdgCode() == 443){
	    TParticle* daugE0 = fStack->Particle(daug0->GetFirstDaughter());
	    TParticle* daugE1 = fStack->Particle(daug0->GetLastDaughter());
	    if( TMath::Abs(daug1->Eta()) < 1.2 && TMath::Abs(daugE0->Eta()) < 0.9 && TMath::Abs(daugE1->Eta()) < 0.9 )
	      fHistograms->FillTable("Table_Electrons",18);

	  }//if
	  else if (daug1->GetPdgCode() == 443){
	    TParticle* daugE0 = fStack->Particle(daug1->GetFirstDaughter());
	    TParticle* daugE1 = fStack->Particle(daug1->GetLastDaughter());
	    if( TMath::Abs(daug0->Eta()) < 1.2 && TMath::Abs(daugE0->Eta()) < 0.9 && TMath::Abs(daugE1->Eta()) < 0.9 )
	      fHistograms->FillTable("Table_Electrons",18);
	  }//else if
	}//gamma o Jpsi
      }//GetNDaughters
    }


    /////////////////////End Chic Analysis////////////////////////////


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
  //Fills the ntuple with the different values

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
  // Process all the V0's without applying any cuts to it

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
    fElectronv1.push_back(fV0Reader->GetCurrentV0()->GetPindex());
    fElectronv2.push_back(fV0Reader->GetCurrentV0()->GetNindex());

		
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
      
      if(fElectronv1[firstGammaIndex]==fElectronv1[secondGammaIndex] || fElectronv1[firstGammaIndex]==fElectronv2[secondGammaIndex]){
	continue;
      }
      if(fElectronv2[firstGammaIndex]==fElectronv1[secondGammaIndex] || fElectronv2[firstGammaIndex]==fElectronv2[secondGammaIndex]){
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
					
	  TVector3 momentumVectorTwoGammaCandidate(twoGammaCandidate->GetPx(),twoGammaCandidate->GetPy(),twoGammaCandidate->GetPz());
	  TVector3 spaceVectorTwoGammaCandidate(twoGammaCandidate->GetX(),twoGammaCandidate->GetY(),twoGammaCandidate->GetZ());
					
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
	  fHistograms->FillHistogram("ESD_Mother_Pt", momentumVectorTwoGammaCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Mother_Eta", momentumVectorTwoGammaCandidate.Eta());
	  fHistograms->FillHistogram("ESD_Mother_Rapid", rapidity);					
	  fHistograms->FillHistogram("ESD_Mother_Phi", spaceVectorTwoGammaCandidate.Phi());
	  fHistograms->FillHistogram("ESD_Mother_Mass", massTwoGammaCandidate);
	  fHistograms->FillHistogram("ESD_Mother_R", spaceVectorTwoGammaCandidate.Pt());    // Pt in Space == R!!!
	  fHistograms->FillHistogram("ESD_Mother_ZR", twoGammaCandidate->GetZ(), spaceVectorTwoGammaCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Mother_XY", twoGammaCandidate->GetX(), twoGammaCandidate->GetY());
	  fHistograms->FillHistogram("ESD_Mother_InvMass_vs_Pt",massTwoGammaCandidate ,momentumVectorTwoGammaCandidate.Pt());
	  fHistograms->FillHistogram("ESD_Mother_InvMass",massTwoGammaCandidate);
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
	  fHistograms->FillHistogram("ESD_Background_InvMass",massBG);
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

Double_t AliAnalysisTaskGammaConversion::GetMCOpeningAngle(TParticle* const daughter0, TParticle* const daughter1) const{
  //helper function
  TVector3 v3D0(daughter0->Px(),daughter0->Py(),daughter0->Pz());
  TVector3 v3D1(daughter1->Px(),daughter1->Py(),daughter1->Pz());
  return v3D0.Angle(v3D1);
}

void AliAnalysisTaskGammaConversion::CheckV0Efficiency(){
  // see header file for documentation

  vector<Int_t> indexOfGammaParticle;

  fStack = fV0Reader->GetMCStack();

  if(fV0Reader->CheckForPrimaryVertex() == kFALSE){
    return; // aborts if the primary vertex does not have contributors.
  }

  for (Int_t iTracks = 0; iTracks < fStack->GetNprimary(); iTracks++) {
    TParticle* particle = (TParticle *)fStack->Particle(iTracks);
    if(particle->GetPdgCode()==22){     //Gamma
      if(particle->GetNDaughters() >= 2){
	TParticle* electron=NULL;
	TParticle* positron=NULL; 
	for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	  TParticle *tmpDaughter = fStack->Particle(daughterIndex);
	  if(tmpDaughter->GetUniqueID() == 5){
	    if(tmpDaughter->GetPdgCode() == 11){
	      electron = tmpDaughter;
	    }
	    else if(tmpDaughter->GetPdgCode() == -11){
	      positron = tmpDaughter;
	    }
	  }
	}
	if(electron!=NULL && positron!=0){
	  if(electron->R()<160){
	    indexOfGammaParticle.push_back(iTracks);
	  }
	}
      }
    }
  }

  Int_t nFoundGammas=0;
  Int_t nNotFoundGammas=0;

  Int_t numberOfV0s = fV0Reader->GetNumberOfV0s();
  for(Int_t i=0;i<numberOfV0s;i++){
    fV0Reader->GetV0(i);
    
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
      //TParticle * v0Gamma = fV0Reader->GetMotherMCParticle();
      for(UInt_t mcIndex=0;mcIndex<indexOfGammaParticle.size();mcIndex++){
	if(negativeMC->GetFirstMother()==indexOfGammaParticle[mcIndex]){
	  nFoundGammas++;
	}
	else{
	  nNotFoundGammas++;
	}
      }
    }
  }
  //  cout<<"Found: "<<nFoundGammas<<"  of: "<<indexOfGammaParticle.size()<<endl;
}


void AliAnalysisTaskGammaConversion::ProcessGammaElectronsForChicAnalysis(){
  // see header file for documantation

  fESDEvent = fV0Reader->GetESDEvent();


  vector <AliESDtrack*> vESDeNegTemp(0);
  vector <AliESDtrack*> vESDePosTemp(0);
  vector <AliESDtrack*> vESDxNegTemp(0);
  vector <AliESDtrack*> vESDxPosTemp(0);
  vector <AliESDtrack*> vESDeNegNoJPsi(0);
  vector <AliESDtrack*> vESDePosNoJPsi(0); 



  fHistograms->FillTable("Table_Electrons",0);//Count number of Events

  for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
    AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);

    if(!curTrack){
      //print warning here
      continue;
    }

    double p[3];if(!curTrack->GetConstrainedPxPyPz(p))continue;
    double r[3];curTrack->GetConstrainedXYZ(r);

    TVector3 rXYZ(r);

    fHistograms->FillTable("Table_Electrons",4);//Count number of ESD tracks

    Bool_t flagKink       =  kTRUE;
    Bool_t flagTPCrefit   =  kTRUE;
    Bool_t flagTRDrefit   =  kTRUE;
    Bool_t flagITSrefit   =  kTRUE;
    Bool_t flagTRDout     =  kTRUE;
    Bool_t flagVertex     =  kTRUE;


    //Cuts ---------------------------------------------------------------

    if(curTrack->GetKinkIndex(0) > 0){
      fHistograms->FillHistogram("Table_Electrons",5);//Count kink
      flagKink = kFALSE;
    }

    ULong_t trkStatus = curTrack->GetStatus();

    ULong_t tpcRefit = (trkStatus & AliESDtrack::kTPCrefit);

    if(!tpcRefit){
      fHistograms->FillHistogram("Table_Electrons",9);//Count not TPCrefit
      flagTPCrefit = kFALSE;
    }

    ULong_t itsRefit = (trkStatus & AliESDtrack::kITSrefit);
    if(!itsRefit){
      fHistograms->FillHistogram("Table_Electrons",10);//Count not ITSrefit
      flagITSrefit = kFALSE;
    }

    ULong_t trdRefit = (trkStatus & AliESDtrack::kTRDrefit);

    if(!trdRefit){
      fHistograms->FillHistogram("Table_Electrons",8); //Count not TRDrefit
      flagTRDrefit = kFALSE;
    }

    ULong_t trdOut = (trkStatus & AliESDtrack::kTRDout);

    if(!trdOut) {
      fHistograms->FillHistogram("Table_Electrons",7); //Count not TRDout
      flagTRDout = kFALSE;
    }

    double nsigmaToVxt = GetSigmaToVertex(curTrack);

    if(nsigmaToVxt > 3){
      fHistograms->FillHistogram("Table_Electrons",6); //Count Tracks with number of sigmas > 3
      flagVertex = kFALSE;
    }

    if(! (flagKink && flagTPCrefit && flagITSrefit && flagTRDrefit && flagTRDout && flagVertex ) ) continue;
    fHistograms->FillHistogram("Table_Electrons",11);//Count Tracks passed Cuts


    Stat_t pid, weight;
    GetPID(curTrack, pid, weight);

    if(pid!=0){
      fHistograms->FillHistogram("Table_Electrons",12); //Count Tracks with pid != 0
    }

    if(pid == 0){
      fHistograms->FillHistogram("Table_Electrons",13); //Count Tracks with pid != 0
    }




    Int_t labelMC = TMath::Abs(curTrack->GetLabel());
    TParticle* curParticle = fStack->Particle(labelMC);




    TLorentzVector curElec;
    curElec.SetXYZM(p[0],p[1],p[2],fElectronMass);




    if(curTrack->GetSign() > 0){

      vESDxPosTemp.push_back(curTrack);

      if( pid == 0){

	fHistograms->FillHistogram("ESD_ElectronPosNegPt",curElec.Pt());
	fHistograms->FillHistogram("ESD_ElectronPosPt",curElec.Pt());
	fHistograms->FillHistogram("MC_ElectronPosNegPt",curParticle->Pt());
	fHistograms->FillHistogram("ESD_ElectronPosNegEta",curElec.Eta());
	fHistograms->FillHistogram("MC_ElectronPosNegEta",curParticle->Eta());
	vESDePosTemp.push_back(curTrack);



      }

    }
    else {
      vESDxNegTemp.push_back(curTrack);

      if( pid == 0){

	fHistograms->FillHistogram("ESD_ElectronPosNegPt",curElec.Pt());
	fHistograms->FillHistogram("ESD_ElectronNegPt",curElec.Pt());
	fHistograms->FillHistogram("MC_ElectronPosNegPt",curParticle->Pt());
	fHistograms->FillHistogram("ESD_ElectronPosNegEta",curElec.Eta());
	fHistograms->FillHistogram("MC_ElectronPosNegEta",curParticle->Eta());
	vESDeNegTemp.push_back(curTrack);




      }
    }

  }


  Bool_t ePosJPsi = kFALSE;
  Bool_t eNegJPsi = kFALSE;		
  Bool_t ePosPi0  = kFALSE;
  Bool_t eNegPi0  = kFALSE;
	
  UInt_t iePosJPsi=0,ieNegJPsi=0,iePosPi0=0,ieNegPi0=0;
 
  for(UInt_t iNeg=0; iNeg < vESDeNegTemp.size(); iNeg++){
    if(fStack->Particle(TMath::Abs(vESDeNegTemp[iNeg]->GetLabel()))->GetPdgCode() == 11)
      if(fStack->Particle(TMath::Abs(vESDeNegTemp[iNeg]->GetLabel()))->GetMother(0) > -1){
	Int_t labelMother = fStack->Particle(TMath::Abs(vESDeNegTemp[iNeg]->GetLabel()))->GetMother(0);
	TParticle* partMother = fStack ->Particle(labelMother);
	if (partMother->GetPdgCode() == 111){
	  ieNegPi0 = iNeg;
	  eNegPi0 = kTRUE;
	}
	if(partMother->GetPdgCode() == 443){ //Mother JPsi
	  fHistograms->FillTable("Table_Electrons",14);
	  ieNegJPsi = iNeg;
	  eNegJPsi = kTRUE;
	}
	else{	
	  vESDeNegNoJPsi.push_back(vESDeNegTemp[iNeg]);
	  //		cout<<"ESD No Positivo JPsi "<<endl;
	}

      }
  }	

  for(UInt_t iPos=0; iPos < vESDePosTemp.size(); iPos++){
    if(fStack->Particle(TMath::Abs(vESDePosTemp[iPos]->GetLabel()))->GetPdgCode() == -11)
      if(fStack->Particle(TMath::Abs(vESDePosTemp[iPos]->GetLabel()))->GetMother(0) > -1){
	Int_t labelMother = fStack->Particle(TMath::Abs(vESDePosTemp[iPos]->GetLabel()))->GetMother(0);
	TParticle* partMother = fStack ->Particle(labelMother);
	if (partMother->GetPdgCode() == 111){
	  iePosPi0 = iPos;
	  ePosPi0 = kTRUE;
	}
	if(partMother->GetPdgCode() == 443){ //Mother JPsi
	  fHistograms->FillTable("Table_Electrons",15);
	  iePosJPsi = iPos;
	  ePosJPsi = kTRUE;
	}
	else{
	  vESDePosNoJPsi.push_back(vESDePosTemp[iPos]);
	  //		cout<<"ESD No Negativo JPsi "<<endl;
	}

      }
  }
	
  if( eNegJPsi && ePosJPsi ){
    TVector3 tempeNegV,tempePosV;
    tempeNegV.SetXYZ(vESDeNegTemp[ieNegJPsi]->Px(),vESDeNegTemp[ieNegJPsi]->Py(),vESDeNegTemp[ieNegJPsi]->Pz());			
    tempePosV.SetXYZ(vESDePosTemp[iePosJPsi]->Px(),vESDePosTemp[iePosJPsi]->Py(),vESDePosTemp[iePosJPsi]->Pz());
    fHistograms->FillTable("Table_Electrons",16);
    fHistograms->FillHistogram("ESD_ElectronPosNegJPsiAngle",tempeNegV.Angle(tempePosV));	
    fHistograms->FillHistogram("MC_ElectronPosNegJPsiAngle",GetMCOpeningAngle(fStack->Particle(TMath::Abs(vESDeNegTemp[ieNegJPsi]->GetLabel())),
									      fStack->Particle(TMath::Abs(vESDePosTemp[iePosJPsi]->GetLabel()))));	
  }
	
  if( eNegPi0 && ePosPi0 ){
    TVector3 tempeNegV,tempePosV;
    tempeNegV.SetXYZ(vESDeNegTemp[ieNegPi0]->Px(),vESDeNegTemp[ieNegPi0]->Py(),vESDeNegTemp[ieNegPi0]->Pz());
    tempePosV.SetXYZ(vESDePosTemp[iePosPi0]->Px(),vESDePosTemp[iePosPi0]->Py(),vESDePosTemp[iePosPi0]->Pz());
    fHistograms->FillHistogram("ESD_ElectronPosNegPi0Angle",tempeNegV.Angle(tempePosV));
    fHistograms->FillHistogram("MC_ElectronPosNegPi0Angle",GetMCOpeningAngle(fStack->Particle(TMath::Abs(vESDeNegTemp[ieNegPi0]->GetLabel())),
									     fStack->Particle(TMath::Abs(vESDePosTemp[iePosPi0]->GetLabel()))));   
  }
  	 

  FillAngle("ESD_eNegePosAngleBeforeCut",GetTLorentzVector(vESDeNegTemp),GetTLorentzVector(vESDePosTemp));

  CleanWithAngleCuts(vESDeNegTemp,vESDePosTemp,fKFReconstructedGammas);
	
  vector <TLorentzVector> vCurrentTLVeNeg = GetTLorentzVector(fCurrentEventNegElectron);
  vector <TLorentzVector> vCurrentTLVePos = GetTLorentzVector(fCurrentEventPosElectron);


  FillAngle("ESD_eNegePosAngleAfterCut",vCurrentTLVeNeg,vCurrentTLVePos);

 


  //FillAngle("ESD_eNegePosAngleAfterCut",CurrentTLVeNeg,CurrentTLVePos);


  FillElectronInvMass("ESD_InvMass_ePluseMinus",vCurrentTLVeNeg,vCurrentTLVePos);
  FillElectronInvMass("ESD_InvMass_xPlusxMinus",GetTLorentzVector(vESDxNegTemp),GetTLorentzVector(vESDxPosTemp));

       

  FillGammaElectronInvMass("ESD_InvMass_GammaePluseMinusChiC","ESD_InvMass_GammaePluseMinusChiCDiff",
			   fKFReconstructedGammasCut,vCurrentTLVeNeg,vCurrentTLVePos);

  FillGammaElectronInvMass("ESD_InvMass_GammaePluseMinusPi0","ESD_InvMass_GammaePluseMinusPi0Diff",
			   fKFReconstructedGammasCut,vCurrentTLVeNeg,vCurrentTLVePos);

  //BackGround

  //Like Sign e+e-
  ElectronBackground("ESD_ENegBackground",vCurrentTLVeNeg);
  ElectronBackground("ESD_EPosBackground",vCurrentTLVePos);
  ElectronBackground("ESD_EPosENegBackground",vCurrentTLVeNeg);
  ElectronBackground("ESD_EPosENegBackground",vCurrentTLVePos);

  //        Like Sign e+e- no JPsi
  ElectronBackground("ESD_EPosENegNoJPsiBG",GetTLorentzVector(vESDeNegNoJPsi));
  ElectronBackground("ESD_EPosENegNoJPsiBG",GetTLorentzVector(vESDePosNoJPsi));

  //Mixed Event

  if( fCurrentEventPosElectron.size() > 0 && fCurrentEventNegElectron.size() > 0 && fKFReconstructedGammasCut.size() > 0 ){
    FillGammaElectronInvMass("ESD_EPosENegGammaBackgroundMX","ESD_EPosENegGammaBackgroundMXDiff",
			     fKFReconstructedGammasCut,fPreviousEventTLVNegElectron,fPreviousEventTLVPosElectron);
    fPreviousEventTLVNegElectron = vCurrentTLVeNeg;
    fPreviousEventTLVPosElectron = vCurrentTLVePos;

  }

  /*
  //Photons P
  Double_t vtx[3];
  vtx[0]=0;vtx[1]=0;vtx[2]=0;
  for(UInt_t i=0;i<fKFReconstructedGammasChic.size();i++){

  //      if(fMCGammaChicTempCut[i]->GetMother(0) < 0) continue;



  Int_t tempLabel = fStack->Particle(fMCGammaChicTempCut[i]->GetMother(0))->GetPdgCode();
  //      cout<<" Label Pedro Gonzalez " <<tempLabel <<endl;

  //      cout<<" Label Distance"<<fKFReconstructedGammasChic[i].GetDistanceFromVertex(vtx)<<endl;

  if( tempLabel == 10441 || tempLabel == 20443 || tempLabel == 445 )

  fHistograms->FillHistogram("ESD_PhotonsMomentum",fKFReconstructedGammasChic[i].GetMomentum());


  }


  */


}

void AliAnalysisTaskGammaConversion::FillAngle(TString histoName,vector <TLorentzVector> tlVeNeg, vector <TLorentzVector> tlVePos){
  //see header file for documentation
  for( UInt_t iNeg=0; iNeg < tlVeNeg.size(); iNeg++){
    for (UInt_t iPos=0; iPos < tlVePos.size(); iPos++){
      fHistograms->FillHistogram(histoName.Data(),tlVeNeg[iNeg].Vect().Angle(tlVePos[iPos].Vect()));
    }
  }
}
void AliAnalysisTaskGammaConversion::FillElectronInvMass(TString histoName, vector <TLorentzVector> eNeg, vector <TLorentzVector> ePos){
  //see header file for documentation
  for( UInt_t n=0; n < eNeg.size(); n++){

    TLorentzVector en = eNeg.at(n);
    for (UInt_t p=0; p < ePos.size(); p++){
      TLorentzVector ep = ePos.at(p);
      TLorentzVector np = ep + en;
      fHistograms->FillHistogram(histoName.Data(),np.M());
    }
  }

}

void AliAnalysisTaskGammaConversion::FillGammaElectronInvMass(TString histoMass,TString histoDiff,vector <AliKFParticle> fKFGammas,
							      vector <TLorentzVector> tlVeNeg,vector<TLorentzVector> tlVePos)
{
  //see header file for documentation

  for( UInt_t iNeg=0; iNeg < tlVeNeg.size(); iNeg++ ){

    for (UInt_t iPos=0; iPos < tlVePos.size(); iPos++){

      TLorentzVector xy = tlVePos[iPos] + tlVeNeg[iNeg];

      for (UInt_t iGam=0; iGam < fKFGammas.size(); iGam++){

	AliKFParticle * gammaCandidate = &fKFGammas[iGam];
	TLorentzVector g;

	g.SetXYZM(gammaCandidate->GetPx(),gammaCandidate->GetPy(),gammaCandidate->GetPz(),fGammaMass);
	TLorentzVector xyg = xy + g;
	fHistograms->FillHistogram(histoMass.Data(),xyg.M());
	fHistograms->FillHistogram(histoDiff.Data(),(xyg.M()-xy.M()));
      }
    }
  }

}
void AliAnalysisTaskGammaConversion::ElectronBackground(TString hBg, vector <TLorentzVector> e)
{
  // see header file for documentation
  for(UInt_t i=0; i < e.size(); i++)
    {
      for (UInt_t j=i+1; j < e.size(); j++)
	{
	  TLorentzVector ee = e[i] + e[j];

	  fHistograms->FillHistogram(hBg.Data(),ee.M());
	}
    }
}


void AliAnalysisTaskGammaConversion::CleanWithAngleCuts(vector <AliESDtrack*> negativeElectrons,
							vector <AliESDtrack*> positiveElectrons, vector <AliKFParticle> gammas){
  // see header file for documentation

  UInt_t  sizeN = negativeElectrons.size();
  UInt_t  sizeP = positiveElectrons.size();
  UInt_t  sizeG = gammas.size();



  vector <Bool_t> xNegBand(sizeN);
  vector <Bool_t> xPosBand(sizeP);
  vector <Bool_t> gammaBand(sizeG);


  for(UInt_t iNeg=0; iNeg < sizeN; iNeg++) xNegBand[iNeg]=kTRUE;
  for(UInt_t iPos=0; iPos < sizeP; iPos++) xPosBand[iPos]=kTRUE;
  for(UInt_t iGam=0; iGam < sizeG; iGam++) gammaBand[iGam]=kTRUE;


  for(UInt_t iPos=0; iPos < sizeP; iPos++){
	
    Double_t aP[3]; positiveElectrons[iPos]->GetConstrainedPxPyPz(aP); 

    TVector3 ePosV(aP[0],aP[1],aP[2]);

    for(UInt_t iNeg=0; iNeg < sizeN; iNeg++){
	
      Double_t aN[3]; negativeElectrons[iNeg]->GetConstrainedPxPyPz(aN); 
      TVector3 eNegV(aN[0],aN[1],aN[2]);

      if(ePosV.Angle(eNegV) < 0.05){ //e+e- from gamma
	xPosBand[iPos]=kFALSE;
	xNegBand[iNeg]=kFALSE;
      }

      for(UInt_t iGam=0; iGam < sizeG; iGam++){
	AliKFParticle* gammaCandidate = &gammas[iGam];
	TVector3 gammaCandidateVector(gammaCandidate->Px(),gammaCandidate->Py(),gammaCandidate->Pz());
	if(ePosV.Angle(gammaCandidateVector) < 0.05 || eNegV.Angle(gammaCandidateVector) < 0.05)
	  gammaBand[iGam]=kFALSE;
      }
    }
  }




  for(UInt_t iPos=0; iPos < sizeP; iPos++){
    if(xPosBand[iPos]){
      fCurrentEventPosElectron.push_back(positiveElectrons[iPos]);
    }
  }
  for(UInt_t iNeg=0;iNeg < sizeN; iNeg++){
    if(xNegBand[iNeg]){
      fCurrentEventNegElectron.push_back(negativeElectrons[iNeg]);
    }
  }
  for(UInt_t iGam=0; iGam < sizeG; iGam++){
    if(gammaBand[iGam]){
      fKFReconstructedGammasCut.push_back(gammas[iGam]);
    }
  }
}


void  AliAnalysisTaskGammaConversion::GetPID(AliESDtrack *track, Stat_t &pid, Stat_t &weight)
{
  // see header file for documentation
  pid = -1;
  weight = -1;

  double wpart[5];
  double wpartbayes[5];

  //get probability of the diffenrent particle types
  track->GetESDpid(wpart);

  // Tentative particle type "concentrations"
  double c[5]={0.01, 0.01, 0.85, 0.10, 0.05};

  //Bayes' formula
  double rcc = 0.;
  for (int i = 0; i < 5; i++)
    {
      rcc+=(c[i] * wpart[i]);
    }



  for (int i=0; i<5; i++) {
    if( rcc!=0){
      wpartbayes[i] = c[i] * wpart[i] / rcc;
    }
  }



  Float_t max=0.;
  int ipid=-1;
  //find most probable particle in ESD pid
  //0:Electron - 1:Muon - 2:Pion - 3:Kaon - 4:Proton
  for (int i = 0; i < 5; i++)
    {
      if (wpartbayes[i] > max)
        {
          ipid = i;
          max = wpartbayes[i];
        }
    }

  pid = ipid;
  weight = max;
}
double AliAnalysisTaskGammaConversion::GetSigmaToVertex(AliESDtrack* t)
{
  // Calculates the number of sigma to the vertex.

  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  t->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-x**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0)
    return -1;

  double d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // stupid rounding problem screws up everything:
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  if (TMath::Exp(-d * d / 2) < 1e-10)
    return 1000;


  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return d;
}
vector <TLorentzVector> AliAnalysisTaskGammaConversion::GetTLorentzVector(vector <AliESDtrack*> esdTrack){

  vector <TLorentzVector> tlVtrack(0);

  for(UInt_t itrack=0; itrack < esdTrack.size(); itrack++){
    double P[3]; esdTrack[itrack]->GetConstrainedPxPyPz(P);
    TLorentzVector currentTrack;
    currentTrack.SetXYZM(P[0],P[1],P[2],fElectronMass);
    tlVtrack.push_back(currentTrack);
  }

  return tlVtrack;
}

