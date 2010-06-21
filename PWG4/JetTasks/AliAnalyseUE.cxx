/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: A.Abrahantes, E.Lopez, S.Vallero                               *
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
#include <TROOT.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAnalyseUE.h"
#include "AliAnalysisTaskUE.h"
#include "AliAnalysisTask.h"
#include "AliHistogramsUE.h"

#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliKFVertex.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliStack.h"

////////////////////////////////////////////////
//--------------------------------------------- 
// Class for transverse regions analysis
//---------------------------------------------
////////////////////////////////////////////////


using namespace std;

ClassImp(AliAnalyseUE)

//-------------------------------------------------------------------
AliAnalyseUE::AliAnalyseUE() :
  TObject(),
  //fTaskUE(0),
  fkAOD(0x0),            
  fDebug(0),
  fSimulateChJetPt(kFALSE),
  fAnaType(1),         
  fAreaReg(1.5393), // Pi*0.7*0.7
  fConeRadius(0.7),
  fFilterBit(0xFF),
  fRegionType(1),
  fUseChargeHadrons(kFALSE),
  fUseChPartJet(kFALSE),
  fUsePositiveCharge(kTRUE),
  fUseSingleCharge(kFALSE),
  fOrdering(1),
  fJet1EtaCut(0.2),
  fJet2DeltaPhiCut(2.616),    // 150 degrees
  fJet2RatioPtCut(0.8),
  fJet3PtCut(15.),
  fTrackEtaCut(0.9),
  fTrackPtCut(0.),
  fHistos(0x0),
  fSumPtRegionPosit(0.),
  fSumPtRegionNegat(0.),
  fSumPtRegionForward(0.),
  fSumPtRegionBackward(0.),
  fMaxPartPtRegion(0.),
  fNTrackRegionPosit(0),
  fNTrackRegionNegat(0),
  fNTrackRegionForward(0),
  fNTrackRegionBackward(0),
  fSettingsTree(0x0)

{
  // constructor
}


//-------------------------------------------------------------------
AliAnalyseUE::AliAnalyseUE(const AliAnalyseUE & original) :
  TObject(original),
  //fTaskUE(original.fTaskUE),
  fkAOD(original.fkAOD),            
  fDebug(original.fDebug),
  fSimulateChJetPt(original.fSimulateChJetPt),
  fAnaType(original.fAnaType),
  fAreaReg(original.fAreaReg),
  fConeRadius(original.fConeRadius),
  fFilterBit(original.fFilterBit),
  fRegionType(original.fRegionType),
  fUseChargeHadrons(original.fUseChargeHadrons),
  fUseChPartJet(original.fUseChPartJet),
  fUsePositiveCharge(original.fUsePositiveCharge),
  fUseSingleCharge(original.fUseSingleCharge),
  fOrdering(original.fOrdering),
  fJet1EtaCut(original.fJet1EtaCut),
  fJet2DeltaPhiCut(original.fJet2DeltaPhiCut),
  fJet2RatioPtCut(original.fJet2RatioPtCut),
  fJet3PtCut(original.fJet3PtCut),
  fTrackEtaCut(original.fTrackEtaCut),
  fTrackPtCut(original.fTrackPtCut),
  fHistos(original.fHistos),
  fSumPtRegionPosit(original.fSumPtRegionPosit),
  fSumPtRegionNegat(original.fSumPtRegionNegat),
  fSumPtRegionForward(original.fSumPtRegionForward),
  fSumPtRegionBackward(original.fSumPtRegionBackward),
  fMaxPartPtRegion(original.fMaxPartPtRegion),
  fNTrackRegionPosit(original.fNTrackRegionPosit),
  fNTrackRegionNegat(original.fNTrackRegionNegat),
  fNTrackRegionForward(original.fNTrackRegionForward),
  fNTrackRegionBackward(original.fNTrackRegionBackward),
  fSettingsTree(original.fSettingsTree)
{
  //copy constructor	
}

//-------------------------------------------------------------------
AliAnalyseUE & AliAnalyseUE::operator = (const AliAnalyseUE & /*source*/)
{
  // assignment operator
  return *this;
}


//-------------------------------------------------------------------
AliAnalyseUE::~AliAnalyseUE(){

  //clear memory
  delete[] fkAOD;
  fkAOD = NULL;

  
  
}


//-------------------------------------------------------------------
void AliAnalyseUE::AnalyseMC(TVector3 *jetVect,AliMCEvent *mcEvent, AliGenPythiaEventHeader *pythiaGenHeader,Int_t conePosition, Bool_t useAliStack, Bool_t constrainDistance, Double_t minDistance){

  // Execute the analysis in case of MC input
  fSumPtRegionPosit = 0.;
  fSumPtRegionNegat = 0.;
  fSumPtRegionForward = 0.;
  fSumPtRegionBackward = 0.;
  fMaxPartPtRegion = 0.;
  fNTrackRegionPosit = 0;
  fNTrackRegionNegat = 0;
  fNTrackRegionForward = 0;
  fNTrackRegionBackward = 0;

  static Double_t const  kPI     = TMath::Pi();
  static Double_t const  k270rad = 270.*kPI/180.;

  //Get Jets from MC header
  Int_t nPythiaGenJets = pythiaGenHeader->NTriggerJets();
  AliAODJet pythiaGenJets[4];
  TVector3 jetVectnew[4];
  Int_t iCount = 0;
  for(int ip = 0;ip < nPythiaGenJets;++ip){
  	if (iCount>3) break;
  	Float_t p[4];
  	pythiaGenHeader->TriggerJet(ip,p);
  	TVector3 tempVect(p[0],p[1],p[2]);
  	if ( TMath::Abs(tempVect.Eta())>fJet1EtaCut ) continue;
	pythiaGenJets[iCount].SetPxPyPzE(p[0],p[1],p[2],p[3]);
	jetVectnew[iCount].SetXYZ(pythiaGenJets[iCount].Px(), pythiaGenJets[iCount].Py(), pythiaGenJets[iCount].Pz());
	iCount++;
	}

  if (!iCount) return;// no jet in eta acceptance
    
  //Search the index of the nearest MC jet to the leading jet reconstructed from the input data
  Int_t index = 0;
  if (constrainDistance){
  	Float_t deltaR = 0.;
	Float_t dRTemp = 0.;
	for (Int_t i=0; i<iCount; i++){
		if (!i) {
			dRTemp = jetVectnew[i].DeltaR(jetVect[0]);
			index = i;
			}

		deltaR = jetVectnew[i].DeltaR(jetVect[0]);
		if (deltaR < dRTemp){
			index = i;
			dRTemp = deltaR;
			}
		}
   
  	if (jetVectnew[index].DeltaR(jetVect[0]) > minDistance) return;
	}

  //Let's add some taste to jet and simulate pt of charged alone 
  //eta and phi are kept as original
  //Play a Normal Distribution
  Float_t random = 1.;  
  if (fSimulateChJetPt){
  	while(1){
  		random = gRandom->Gaus(0.6,0.25);
  		if (random > 0. && random < 1. && 
  		(random * jetVectnew[index].Pt()>6.)) break;
  		}
  	}
    
  //Set new Pt & Fill histogram accordingly
  Double_t maxPtJet1 = random * jetVectnew[index].Pt();  
    
  fHistos->FillHistogram("hEleadingPt", maxPtJet1 );    
    
  if (useAliStack){//Try Stack Information to perform UE analysis
    
  	AliStack* mcStack = mcEvent->Stack();//Load Stack
  	Int_t nTracksMC = mcStack->GetNtrack();
  	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
  		//Cuts
  		if(!(mcStack->IsPhysicalPrimary(iTracks))) continue;
        
  		TParticle* mctrk = mcStack->Particle(iTracks);
        
  		Double_t charge = mctrk->GetPDG()->Charge();
		Double_t pT = mctrk->Pt();
		Double_t eta = mctrk->Eta();
		Int_t pdgCode = mctrk->GetPdgCode();

		if (!TrackMCSelected(charge, pT, eta, pdgCode))continue;

  		TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());
  		Double_t deltaPhi = jetVectnew[index].DeltaPhi(partVect)+k270rad;
  		if( deltaPhi > 2.*TMath::Pi() )  deltaPhi-= 2.*TMath::Pi();
  		fHistos->FillHistogram("hdNdEtaPhiDist",deltaPhi, maxPtJet1 ); 
        
  		fHistos->FillHistogram("hFullRegPartPtDistVsEt", mctrk->Pt(), maxPtJet1 ); 
        
		//We are not interested on stack organization but don't loose track of info
        
		TVector3 tempVector =  jetVectnew[0];
        	jetVectnew[0] = jetVectnew[index];
        	jetVectnew[index] = tempVector;
        
        	Int_t region = IsTrackInsideRegion( jetVectnew, &partVect, conePosition );  
        
        	if (region == 1) {
        		if( fMaxPartPtRegion < mctrk->Pt() ) fMaxPartPtRegion = mctrk->Pt();
          		fSumPtRegionPosit += mctrk->Pt();
          		fNTrackRegionPosit++;
          		fHistos->FillHistogram("hTransRegPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        		}
        	if (region == -1) {
          		if( fMaxPartPtRegion < mctrk->Pt() ) fMaxPartPtRegion = mctrk->Pt();
          		fSumPtRegionNegat += mctrk->Pt();
          		fNTrackRegionNegat++;
          		fHistos->FillHistogram("hTransRegPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        		}
        	if (region == 2){ //forward
          		fSumPtRegionForward += mctrk->Pt();
          		fNTrackRegionForward++;
          		fHistos->FillHistogram("hRegForwardPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        		}
        	if (region == -2){ //backward
          		fSumPtRegionBackward += mctrk->Pt();
          		fNTrackRegionBackward++;
          		fHistos->FillHistogram("hRegBackwardPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        		}
    		} //end loop on stack particles     
    }else{//Try mc Particle

      TClonesArray* farray = (TClonesArray*)fkAOD->FindListObject("mcparticles");
       
      Int_t ntrks = farray->GetEntries();
      if (fDebug>1) AliInfo(Form("In UE MC analysis tracks %d \n",ntrks));
      for (Int_t i =0 ; i < ntrks; i++){   
      	AliAODMCParticle* mctrk = (AliAODMCParticle*)farray->At(i);
        //Cuts
        if (!(mctrk->IsPhysicalPrimary())) continue;
        //if (!(mctrk->IsPrimary())) continue;
        
  	Double_t charge = mctrk->Charge(); 
	Double_t pT = mctrk->Pt();
	Double_t eta = mctrk->Eta();
	Int_t pdgCode = mctrk->GetPdgCode();

	if (!TrackMCSelected(charge, pT, eta, pdgCode))continue;
        
        TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());

        Double_t deltaPhi = jetVectnew[index].DeltaPhi(partVect)+k270rad;
        if( deltaPhi > 2.*TMath::Pi() )  deltaPhi-= 2.*TMath::Pi();
        fHistos->FillHistogram("hdNdEtaPhiDist", deltaPhi, maxPtJet1 );

	fHistos->FillHistogram("hFullRegPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        
        //We are not interested on stack organization but don't loose track of info
        TVector3 tempVector =  jetVectnew[0];
        jetVectnew[0] = jetVectnew[index];
        jetVectnew[index] = tempVector;
        
        Int_t region = IsTrackInsideRegion( jetVectnew, &partVect, conePosition );  
        
        if (region == 1) { //right
          if( fMaxPartPtRegion < mctrk->Pt() ) fMaxPartPtRegion = mctrk->Pt();
          fSumPtRegionPosit += mctrk->Pt();
          fNTrackRegionPosit++;
	  fHistos->FillHistogram("hTransRegPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        }
        if (region == -1) { //left
          if( fMaxPartPtRegion < mctrk->Pt() ) fMaxPartPtRegion = mctrk->Pt();
          fSumPtRegionNegat += mctrk->Pt();
          fNTrackRegionNegat++;
          fHistos->FillHistogram("hTransRegPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        }
        if (region == 2){ //forward
          fSumPtRegionForward += mctrk->Pt();
          fNTrackRegionForward++;
          fHistos->FillHistogram("hRegForwardPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        }
        if (region == -2){ //backward
          fSumPtRegionBackward += mctrk->Pt();
          fNTrackRegionBackward++;
	  fHistos->FillHistogram("hRegBackwardPartPtDistVsEt", mctrk->Pt(), maxPtJet1 );
        }
        
      }//end loop AliAODMCParticle tracks
   }  
}



//-------------------------------------------------------------------
Bool_t AliAnalyseUE::AnaTypeSelection(TVector3 *jetVect ){

  // Cut events by jets topology
  // anaType:
  //     1 = inclusive,
  //         - Jet1 |eta| < jet1EtaCut
  //     2 = back to back inclusive
  //         - fulfill case 1
  //         - |Jet1.Phi - Jet2.Phi| > jet2DeltaPhiCut
  //         - Jet2.Pt/Jet1Pt > jet2RatioPtCut
  //     3 = back to back exclusive
  //         - fulfill case 2
  //         - Jet3.Pt < jet3PtCut

  Double_t eta=jetVect[0].Eta();
  if( TMath::Abs(eta) > fJet1EtaCut) {
  	if( fDebug > 1 ) AliInfo("\n   Skipping Event...Jet1 |eta| > fJet1EtaCut");
  	return kFALSE;
  	}
  // back to back inclusive
  if( fAnaType > 1 && fAnaType < 4 && jetVect[1].Pt() < 0. ) {
  	if( fDebug > 1 ) AliInfo("\n   Skipping Event... no second Jet found");
  	return kFALSE;
  	}
  if( fAnaType > 1 && fAnaType < 4 && jetVect[1].Pt() > 0. ) {
  	if( TMath::Abs(jetVect[0].DeltaPhi(jetVect[1])) < fJet2DeltaPhiCut ||
  	jetVect[1].Pt()/jetVect[0].Pt() < fJet2RatioPtCut ) {
  		if( fDebug > 1 ) AliInfo("\n   Skipping Event... |Jet1.Phi - Jet2.Phi| < fJet2DeltaPhiCut");
  		return kFALSE;
  		}
  	}
  // back to back exclusive
  if( fAnaType > 2 && fAnaType < 4 && jetVect[2].Pt() > 0. ) {
  	if( jetVect[2].Pt() > fJet3PtCut ) {
      		if( fDebug > 1 ) AliInfo("\n   Skipping Event... Jet3.Pt > fJet3PtCut ");
      		return kFALSE;
    		}
  	}
  return kTRUE;  
}


//-------------------------------------------------------------------
TList* AliAnalyseUE::CreateHistograms(Int_t bins, Double_t min, Double_t max){
  
  //Initialize histograms from class AliHistogramsUE
  fHistos = new AliHistogramsUE();
  TList* list = new TList();
  list = fHistos->CreateHistos(bins, min, max, fTrackEtaCut);


  return list;
}



//-------------------------------------------------------------------
void AliAnalyseUE::FillLeadingJet( Double_t  w){

 fHistos->FillHistogram("hEleadingPt",w);

}


//-------------------------------------------------------------------
void AliAnalyseUE::FillRegions(Bool_t isNorm2Area,  TVector3 *jetVect){
  
  // Fill the different topological regions
  Double_t maxPtJet1 = jetVect[0].Pt();
  static Double_t const  kPI     = TMath::Pi();
  static Double_t const  k120rad = 120.*kPI/180.;
  Double_t const kMyTolerance = 0.0000001;

  //Area for Normalization 
  // Forward and backward
  Double_t normArea = 1.;
  // Transverse
  if (isNorm2Area) {
  	SetRegionArea(jetVect);
    	normArea =  2.*fTrackEtaCut*k120rad ;
  	} else fAreaReg = 1.;
  
  Double_t avePosRegion = (fNTrackRegionPosit) ? fSumPtRegionPosit/fNTrackRegionPosit : 0.;
  Double_t aveNegRegion = (fNTrackRegionNegat) ? fSumPtRegionNegat/fNTrackRegionNegat : 0.;
  if( avePosRegion > aveNegRegion ) {
     FillAvePartPtRegion( maxPtJet1, avePosRegion/fAreaReg, aveNegRegion/fAreaReg );
  } else {
     FillAvePartPtRegion( maxPtJet1, aveNegRegion/fAreaReg, avePosRegion/fAreaReg );
  }
  
  //How quantities will be sorted before Fill Min and Max Histogram
  //  1=Plots will be CDF-like
  //  2=Plots will be Marchesini-like
  //  3=Minimum zone is selected as the one having lowest pt per track 
  if( fOrdering == 1 ) {
    if( fSumPtRegionPosit > fSumPtRegionNegat ) {
      FillSumPtRegion( maxPtJet1, fSumPtRegionPosit/fAreaReg, fSumPtRegionNegat/fAreaReg );
    } else {
      FillSumPtRegion( maxPtJet1, fSumPtRegionNegat/fAreaReg, fSumPtRegionPosit/fAreaReg );
    }
    if (fNTrackRegionPosit > fNTrackRegionNegat ) {
      FillMultRegion( maxPtJet1, fNTrackRegionPosit/fAreaReg, fNTrackRegionNegat/fAreaReg, fSumPtRegionNegat/fAreaReg );
    } else {
      FillMultRegion( maxPtJet1, fNTrackRegionNegat/fAreaReg, fNTrackRegionPosit/fAreaReg, fSumPtRegionPosit/fAreaReg );
    }
  } else if( fOrdering == 2 ) {
    if (fSumPtRegionPosit > fSumPtRegionNegat) {
      FillSumPtRegion( maxPtJet1, fSumPtRegionPosit/fAreaReg, fSumPtRegionNegat/fAreaReg );
      FillMultRegion( maxPtJet1, fNTrackRegionPosit/fAreaReg, fNTrackRegionNegat/fAreaReg, fSumPtRegionNegat/fAreaReg );
    } else {
      FillSumPtRegion( maxPtJet1, fSumPtRegionNegat/fAreaReg, fSumPtRegionPosit/fAreaReg );
      FillMultRegion( maxPtJet1, fNTrackRegionNegat/fAreaReg, fNTrackRegionPosit/fAreaReg, fSumPtRegionPosit/fAreaReg );
    }
  } else if( fOrdering == 3 ){
     if (avePosRegion > aveNegRegion) {
        FillSumPtRegion( maxPtJet1, fSumPtRegionPosit/fAreaReg, fSumPtRegionNegat/fAreaReg );
        FillMultRegion( maxPtJet1, fNTrackRegionPosit/fAreaReg, fNTrackRegionNegat/fAreaReg, fSumPtRegionNegat/fAreaReg );
     }else{
        FillSumPtRegion( maxPtJet1, fSumPtRegionNegat/fAreaReg, fSumPtRegionPosit/fAreaReg );
        FillMultRegion( maxPtJet1, fNTrackRegionNegat/fAreaReg, fNTrackRegionPosit/fAreaReg, fSumPtRegionPosit/fAreaReg );
     }
  }
  fHistos->FillHistogram("hRegionMaxPartPtMaxVsEt",maxPtJet1, fMaxPartPtRegion);  
  
  // Compute pedestal like magnitudes
  fHistos->FillHistogram("hRegionDiffSumPtVsEt",maxPtJet1, (TMath::Abs(fSumPtRegionPosit-fSumPtRegionNegat)/(2.0*fAreaReg))+kMyTolerance);
  fHistos->FillHistogram("hRegionAveSumPtVsEt", maxPtJet1, (fSumPtRegionPosit+fSumPtRegionNegat)/(2.0*fAreaReg));

  // Transverse as a whole
  fHistos->FillHistogram("hRegTransMult", maxPtJet1, fNTrackRegionPosit + fNTrackRegionNegat, (fNTrackRegionPosit + fNTrackRegionNegat)/(2.0*fAreaReg));
 fHistos->FillHistogram("hRegTransSumPtVsMult",maxPtJet1, fNTrackRegionPosit + fNTrackRegionNegat , (fSumPtRegionNegat + fSumPtRegionPosit)/(2.0 *fAreaReg));

  // Fill Histograms for Forward and away side w.r.t. leading jet direction
  // Pt dependence
  //fHistos->FillHistogram("hRegForwardSumPtVsEt",maxPtJet1, fSumPtRegionForward/normArea );
  //fHistos->FillHistogram("hRegForwardMultVsEt",maxPtJet1, fNTrackRegionForward/normArea );
  //fHistos->FillHistogram("hRegBackwardSumPtVsEt",maxPtJet1, fSumPtRegionBackward/normArea );
  //fHistos->FillHistogram("hRegBackwardMultVsEt",maxPtJet1, fNTrackRegionBackward/normArea);
  
  // Multiplicity dependence
  fHistos->FillHistogram("hRegForwardMult", maxPtJet1, fNTrackRegionForward, fNTrackRegionForward/normArea);
  fHistos->FillHistogram("hRegForwardSumPtvsMult", maxPtJet1, fNTrackRegionForward,fSumPtRegionForward/normArea);
  fHistos->FillHistogram("hRegBackwardMult", maxPtJet1, fNTrackRegionBackward, fNTrackRegionBackward/normArea );
  fHistos->FillHistogram("hRegBackwardSumPtvsMult", maxPtJet1, fNTrackRegionBackward,fSumPtRegionBackward/normArea);
}

//-------------------------------------------------------------------
void AliAnalyseUE::FillTrials(const char *namex, Double_t  w){

 fHistos->GetTrials()->Fill(namex,w);

}

//-------------------------------------------------------------------
void AliAnalyseUE::FillVertex(Double_t  w){

 fHistos->FillHistogram("hVertexMult",w);

}

//-------------------------------------------------------------------
void AliAnalyseUE::FillXsec(const char *namex, Double_t  w){

 fHistos->GetXsec()->Fill(namex,w);

}


//-------------------------------------------------------------------
void AliAnalyseUE::FindMaxMinRegions(TVector3 *jetVect, Int_t conePosition){
  
  // Identify the different topological zones
  fSumPtRegionPosit = 0.;
  fSumPtRegionNegat = 0.;
  fSumPtRegionForward = 0.;
  fSumPtRegionBackward = 0.;
  fMaxPartPtRegion = 0.;
  fNTrackRegionPosit = 0;
  fNTrackRegionNegat = 0;
  fNTrackRegionForward = 0;
  fNTrackRegionBackward = 0;
  static Double_t const  kPI     = TMath::Pi();
  static Double_t const  kTWOPI  = 2.*kPI;
  static Double_t const  k270rad = 270.*kPI/180.;
  Double_t const kMyTolerance = 0.0000001;

    Int_t nTracks = fkAOD->GetNTracks();
    if (fDebug > 1) AliInfo(Form(" ==== AOD tracks = %d \n ",nTracks));
    
    for (Int_t ipart=0; ipart<nTracks; ++ipart) {
      
    AliAODTrack* part = fkAOD->GetTrack( ipart );
    if (fDebug > 1) AliInfo(Form(" ==== AOD track = %d pt = %f charge = %d \n ",ipart,part->Pt(),part->Charge()));
    // track selection
    if (! TrackSelected(part)) continue;
      
      
    TVector3 partVect(part->Px(), part->Py(), part->Pz());
    Bool_t isFlagPart = kTRUE;
    Double_t deltaPhi = jetVect[0].DeltaPhi(partVect)+k270rad;
    if( deltaPhi > kTWOPI )  deltaPhi-= kTWOPI;
    if (fAnaType != 4 ) fHistos->FillHistogram("hdNdEtaPhiDist",deltaPhi, jetVect[0].Pt());
    else if (TMath::Abs(deltaPhi-k270rad) >= kMyTolerance && TMath::Abs(jetVect[0].Eta()-partVect.Eta()) >= kMyTolerance){
    	fHistos->FillHistogram("hdNdEtaPhiDist",deltaPhi, jetVect[0].Pt());
    	isFlagPart = kFALSE;
    	}
      
    fHistos->FillHistogram("hFullRegPartPtDistVsEt", part->Pt(), jetVect[0].Pt());  
    
    Int_t region = IsTrackInsideRegion( jetVect, &partVect, conePosition );  
    if (region == 1) {
    	if( fMaxPartPtRegion < part->Pt() ) fMaxPartPtRegion = part->Pt();
        fSumPtRegionPosit += part->Pt();
        fNTrackRegionPosit++;
      	fHistos->FillHistogram("hTransRegPartPtDistVsEt",part->Pt(), jetVect[0].Pt());
	}
    if (region == -1) {
    	if( fMaxPartPtRegion < part->Pt() ) fMaxPartPtRegion = part->Pt();
        fSumPtRegionNegat += part->Pt();
        fNTrackRegionNegat++;
      	fHistos->FillHistogram("hTransRegPartPtDistVsEt",part->Pt(), jetVect[0].Pt());
      	}
    if (region == 2){ //forward
    	fSumPtRegionForward += part->Pt();
        fNTrackRegionForward++;
      	fHistos->FillHistogram("hRegForwardPartPtDistVsEt",part->Pt(), jetVect[0].Pt());
      	}
    if (region == -2){ //backward
    	fSumPtRegionBackward += part->Pt();
        fNTrackRegionBackward++;
      	fHistos->FillHistogram("hRegBackwardPartPtDistVsEt",part->Pt(), jetVect[0].Pt());
      	}
    }//end loop AOD tracks

}


//-------------------------------------------------------------------
TList* AliAnalyseUE::GetHistograms(){

  TList *list = fHistos->GetListOfHistos();

  return list;

}


//-------------------------------------------------------------------
TVector3 AliAnalyseUE::GetOrderedClusters(TString aodBranch, Bool_t chargedJets, Double_t chJetPtMin){ 

  // jets from AOD, on-the-fly or leading particle
  Double_t maxPtJet1 = 0.; 
  Int_t    index1 = -1;
  Double_t maxPtJet2 = 0.;   // jet 2 need for back to back inclusive
  Int_t    index2 = -1;
  Double_t maxPtJet3 = 0.;   // jet 3 need for back to back exclusive
  Int_t    index3 = -1;
  TVector3 jetVect[3];
  
  jetVect[0].SetPtEtaPhi(-1.,-1.,-1.);
  jetVect[1].SetPtEtaPhi(-1.,-1.,-1.);
  jetVect[2].SetPtEtaPhi(-1.,-1.,-1.);
  
  Int_t nJets = 0;
  //TClonesArray* fArrayJets;
  TObjArray* arrayJets;
  // 1) JETS FROM AOD BRANCH (standard, non-standard or delta)
  if (!chargedJets && fAnaType != 4 ) { 
  	AliInfo(" ==== Read AODs  !");
  	AliInfo(Form(" ====  Reading Branch: %s  ", aodBranch.Data()));
  	arrayJets = (TObjArray*)fkAOD->GetList()->FindObject(aodBranch.Data());
	if (!arrayJets){
	       AliFatal(" No jet-array! ");
	       return *jetVect;
	       }

	// Find Leading Jets 1,2,3 
  	// (could be skipped if Jets are sort by Pt...)
    	nJets=arrayJets->GetEntries();
    	for( Int_t i=0; i<nJets; ++i ) {
      		AliAODJet* jet = (AliAODJet*)arrayJets->At(i);
      		Double_t jetPt = jet->Pt();//*1.666; // FIXME Jet Pt Correction ?????!!!
 
      		if( jetPt > maxPtJet1 ) {
	     		maxPtJet3 = maxPtJet2; index3 = index2;
	     		maxPtJet2 = maxPtJet1; index2 = index1;
	     		maxPtJet1 = jetPt; index1 = i;
      			} else if( jetPt > maxPtJet2 ) {
	     		maxPtJet3 = maxPtJet2; index3 = index2;
	     		maxPtJet2 = jetPt; index2 = i;
      			} else if( jetPt > maxPtJet3 ) {
	     		maxPtJet3 = jetPt; index3 = i;
      			}
    		}

    	if( index1 != -1 ) {
      		AliAODJet *jet =(AliAODJet*) arrayJets->At(index1);
      		if(jet)jetVect[0].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
    		}
    	if( index2 != -1 ) {
      		AliAODJet* jet= (AliAODJet*) arrayJets->At(index2);
      		if(jet)jetVect[1].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
    		}
    	if( index3 != -1 ) {
       		AliAODJet* jet = (AliAODJet*) arrayJets->At(index3);
      		if(jet)jetVect[2].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
    		}
    
  }


  // 2) ON-THE-FLY CDF ALGORITHM
  if (chargedJets){ 
    // Printf(" ==== Run CDF algorithm on the fly  !");
  	arrayJets = FindChargedParticleJets(chJetPtMin);
        if( arrayJets ) {
        	nJets = arrayJets->GetEntriesFast();
        	if( nJets > 0 ) {
	     		index1 = 0;
	     		AliAODJet* jet = (AliAODJet*)arrayJets->At(0);
	     		maxPtJet1 = jet->Pt();
	     		jetVect[0].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      			}
      		if( nJets > 1 ) {
	     		index2 = 1;
	     		AliAODJet* jet = (AliAODJet*)arrayJets->At(1);
        		maxPtJet2 = jet->Pt();
	     		jetVect[1].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      			}
      		if( nJets > 2 ) {
	     		index3 = 2;
	     		AliAODJet* jet = (AliAODJet*)arrayJets->At(2);
        		maxPtJet3 = jet->Pt();
	     		jetVect[2].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      			}
      
      		arrayJets->Delete();
      		delete arrayJets;
    		}
  	}
  

  // 3) LEADING PARTICLE
  if( fAnaType == 4 ){
  	TObjArray* tracks = SortChargedParticles();
    	if( tracks ) {
      		nJets = tracks->GetEntriesFast();
      		if( nJets > 0 ) {
        		index1 = 0;
        		AliAODTrack* jet = (AliAODTrack*)tracks->At(0);
        		maxPtJet1 = jet->Pt();
        		jetVect[0].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      			}
      		tracks->Clear();
      		delete tracks; 
    		}

  	}
  //fHistos->FillHistoNJets(nJets);
  if (fHistos ) fHistos->FillHistogram("hNJets",nJets);
   
  return *jetVect;

}


//-------------------------------------------------------------------
void AliAnalyseUE::Initialize(AliAnalysisTaskUE& taskUE){
   
  //Get principal settings from current instance of UE analysis-task
  fAnaType = taskUE.GetAnaTopology();         
  fkAOD = taskUE.GetAOD();           
  fConeRadius = taskUE.GetConeRadius();
  fDebug = taskUE.GetDebugLevel();
  fFilterBit = taskUE.GetFilterBit();
  fJet1EtaCut = taskUE.GetJet1EtaCut();
  fJet2DeltaPhiCut = taskUE.GetJet2DeltaPhiCut();
  fJet2RatioPtCut = taskUE.GetJet2RatioPtCut();
  fJet3PtCut = taskUE.GetJet3PtCut();
  fOrdering = taskUE.GetPtSumOrdering() ;
  fRegionType = taskUE.GetRegionType();
  fSimulateChJetPt = taskUE.GetSimulateChJetPt();
  fTrackEtaCut = taskUE.GetTrackEtaCut(); 
  fTrackPtCut = taskUE.GetTrackPtCut();
  fUseChargeHadrons = taskUE.GetUseChargeHadrons();
  fUseChPartJet = taskUE.GetUseChPartJet();
  fUsePositiveCharge = taskUE.GetUseNegativeChargeType();
  fUseSingleCharge = taskUE.GetUseSingleCharge();
  
  //Write settings to output list
  fSettingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  fSettingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  fSettingsTree->Branch("fConeRadius", &fConeRadius,"Rad/D");
  fSettingsTree->Branch("fJet1EtaCut", &fJet1EtaCut, "LeadJetEtaCut/D");
  fSettingsTree->Branch("fJet2DeltaPhiCut", &fJet2DeltaPhiCut, "DeltaPhi/D");
  fSettingsTree->Branch("fJet2RatioPtCut", &fJet2RatioPtCut, "Jet2Ratio/D");
  fSettingsTree->Branch("fJet3PtCut", &fJet3PtCut, "Jet3PtCut/D");
  fSettingsTree->Branch("fTrackPtCut", &fTrackPtCut, "TrackPtCut/D");
  fSettingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  fSettingsTree->Branch("fAnaType", &fAnaType, "Ana/I");        
  fSettingsTree->Branch("fRegionType", &fRegionType,"Reg/I");
  fSettingsTree->Branch("fOrdering", &fOrdering,"OrderMeth/I");
  fSettingsTree->Branch("fUseChPartJet", &fUseChPartJet,"UseChPart/O");
  fSettingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  fSettingsTree->Branch("fUseSingleCharge", &fUseSingleCharge,"UseSingleCh/O");
  fSettingsTree->Branch("fUsePositiveCharge", &fUsePositiveCharge,"UsePositiveCh/O");
  fSettingsTree->Fill();
  (fHistos->GetListOfHistos())->Add(fSettingsTree);
 
}

//-------------------------------------------------------------------
void AliAnalyseUE::Initialize(Int_t anaType,AliAODEvent* aod,Double_t coneRadius, Int_t debug, Int_t filterBit, Double_t jet1EtaCut, Double_t jet2DeltaPhiCut, Double_t jet2RatioPtCut, Double_t jet3PtCut, Int_t ordering, Int_t regionType,Bool_t simulateChJetPt, Double_t trackEtaCut, Double_t trackPtCut, Bool_t useChargeHadrons, Bool_t useChPartJet, Bool_t useNegativeChargeType, Bool_t useSingleCharge ){
   
  //Get principal settings from generic analysis-task
  fAnaType = anaType;         
  fkAOD = aod;           
  fConeRadius = coneRadius;
  fDebug = debug;
  fFilterBit = filterBit;
  fJet1EtaCut = jet1EtaCut;
  fJet2DeltaPhiCut = jet2DeltaPhiCut;
  fJet2RatioPtCut = jet2RatioPtCut;
  fJet3PtCut = jet3PtCut;
  fOrdering = ordering ;
  fRegionType = regionType;
  fSimulateChJetPt = simulateChJetPt;
  fTrackEtaCut = trackEtaCut; 
  fTrackPtCut = trackPtCut;
  fUseChargeHadrons = useChargeHadrons;
  fUseChPartJet = useChPartJet;
  fUsePositiveCharge = useNegativeChargeType;
  fUseSingleCharge = useSingleCharge;
  
  //Write settings to output list
  fSettingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  fSettingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  fSettingsTree->Branch("fConeRadius", &fConeRadius,"Rad/D");
  fSettingsTree->Branch("fJet1EtaCut", &fJet1EtaCut, "LeadJetEtaCut/D");
  fSettingsTree->Branch("fJet2DeltaPhiCut", &fJet2DeltaPhiCut, "DeltaPhi/D");
  fSettingsTree->Branch("fJet2RatioPtCut", &fJet2RatioPtCut, "Jet2Ratio/D");
  fSettingsTree->Branch("fJet3PtCut", &fJet3PtCut, "Jet3PtCut/D");
  fSettingsTree->Branch("fTrackPtCut", &fTrackPtCut, "TrackPtCut/D");
  fSettingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  fSettingsTree->Branch("fAnaType", &fAnaType, "Ana/I");        
  fSettingsTree->Branch("fRegionType", &fRegionType,"Reg/I");
  fSettingsTree->Branch("fOrdering", &fOrdering,"OrderMeth/I");
  fSettingsTree->Branch("fUseChPartJet", &fUseChPartJet,"UseChPart/O");
  fSettingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  fSettingsTree->Branch("fUseSingleCharge", &fUseSingleCharge,"UseSingleCh/O");
  fSettingsTree->Branch("fUsePositiveCharge", &fUsePositiveCharge,"UsePositiveCh/O");
  fSettingsTree->Fill();
  (fHistos->GetListOfHistos())->Add(fSettingsTree);
 
}

Bool_t  AliAnalyseUE::VertexSelection(AliAODEvent *aod, Int_t tracks, Double_t zed ){

  //Require 1 vertex (no TPC stand-alone) with a minimum number of tracks and z-coordinate in a limited range
  Int_t nVertex = aod->GetNumberOfVertices();
  if( nVertex > 0 ) { // Only one vertex (reject pileup)
  	AliAODVertex* vertex = (AliAODVertex*)aod->GetPrimaryVertex();
  	Int_t nTracksPrim = vertex->GetNContributors();
  	Double_t zVertex = vertex->GetZ();
  	if (fDebug > 1) AliInfo(Form(" Vertex in = %f with %d particles by  %s data ...",zVertex,nTracksPrim,vertex->GetName()));
  	// Select a quality vertex by number of tracks?
  	if( nTracksPrim < tracks || TMath::Abs(zVertex) > zed ) {
  		if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
  		return kFALSE;
  		}
  	if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED...");
  	} else {
  		if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
  		return kFALSE;
  		}

  return kTRUE;
}

// PRIVATE METHODS **************************************************

TObjArray*  AliAnalyseUE::FindChargedParticleJets( Double_t chJetPtMin )
{
  // Return a TObjArray of "charged particle jets"
  
  // Charged particle jet definition from reference:
  // "Charged jet evolution and the underlying event
  //  in proton-antiproton collisions at 1.8 TeV"
  //  PHYSICAL REVIEW D 65 092002, CDF Collaboration
  
  // We defined "jets" as circular regions in eta-phi space with
  // radius defined by R = sqrt( (eta-eta0)^2 +(phi-phi0)^2 ).
  // Our jet algorithm is as follows:
  //   1- Order all charged particles according to their pT .
  //   2- Start with the highest pT particle and include in the jet all
  //      particles within the radius R=0.7 considering each particle
  //      in the order of decreasing pT and recalculating the centroid
  //      of the jet after each new particle is added to the jet .
  //   3- Go to the next highest pT particle not already included in
  //      a jet and add to the jet all particles not already included in
  //      a jet within R=0.7.
  //   4- Continue until all particles are in a jet.
  // We defined the transverse momentum of the jet to be
  // the scalar pT sum of all the particles within the jet, where pT
  // is measured with respect to the beam axis
  
  //  1 - Order all charged particles according to their pT .
  Int_t nTracks = fkAOD->GetNTracks();
  if( !nTracks ) return 0;
  TObjArray tracks(nTracks);
  
  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
  	AliAODTrack* part = fkAOD->GetTrack( ipart );
  	if( !part->TestFilterBit(fFilterBit) ) continue; // track cut selection
  	if( !part->Charge() ) continue;
  	if( part->Pt() < fTrackPtCut ) continue;
  	tracks.AddLast(part);
  	}
  QSortTracks( tracks, 0, tracks.GetEntriesFast() );
  
  nTracks = tracks.GetEntriesFast();
  if( !nTracks ) return 0;

  TObjArray *jets = new TObjArray(nTracks);
  TIter itrack(&tracks);
  while( nTracks ) {
  	// 2- Start with the highest pT particle ...
  	Float_t px,py,pz,pt; 
  	AliAODTrack* track = (AliAODTrack*)itrack.Next();
  	if( !track ) continue;
  	px = track->Px();
  	py = track->Py();
  	pz = track->Pz();
  	pt = track->Pt(); // Use the energy member to store Pt
  	jets->AddLast( new TLorentzVector(px, py, pz, pt) );
  	tracks.Remove( track );
  	TLorentzVector* jet = (TLorentzVector*)jets->Last();
  	jet->SetPtEtaPhiE( 1., jet->Eta(), jet->Phi(), pt );
  	// 3- Go to the next highest pT particle not already included...
  	AliAODTrack* track1;
  	Double_t fPt = jet->E();
  	while ( (track1  = (AliAODTrack*)(itrack.Next())) ) {
  		Double_t tphi = track1->Phi(); // here Phi is from 0 <-> 2Pi
  		if (tphi > TMath::Pi()) tphi -= 2. * TMath::Pi(); // convert to  -Pi <-> Pi
  		Double_t dphi = TVector2::Phi_mpi_pi(jet->Phi()-tphi);
  		Double_t r = TMath::Sqrt( (jet->Eta()-track1->Eta())*(jet->Eta()-track1->Eta()) +dphi*dphi );
  		if( r < fConeRadius ) {
  			fPt   = jet->E()+track1->Pt();  // Scalar sum of Pt
  			// recalculating the centroid
  			Double_t eta = jet->Eta()*jet->E()/fPt + track1->Eta()*track1->Pt()/fPt;
  			Double_t phi = jet->Phi()*jet->E()/fPt + tphi*track1->Pt()/fPt;
  			jet->SetPtEtaPhiE( 1., eta, phi, fPt );
  			tracks.Remove( track1 );
  			}
  		}
    
  	tracks.Compress();
  	nTracks = tracks.GetEntries();
  	//   4- Continue until all particles are in a jet.
  	itrack.Reset();
  	} // end while nTracks
  
  // Convert to AODjets....
  Int_t njets = jets->GetEntriesFast();
  TObjArray* aodjets = new TObjArray(njets);
  aodjets->SetOwner(kTRUE);
  for(Int_t ijet=0; ijet<njets; ++ijet) {
  	TLorentzVector* jet = (TLorentzVector*)jets->At(ijet);
  	if (jet->E() < chJetPtMin) continue;
  	Float_t px, py,pz,en; // convert to 4-vector
  	px = jet->E() * TMath::Cos(jet->Phi());  // Pt * cos(phi)
  	py = jet->E() * TMath::Sin(jet->Phi());  // Pt * sin(phi)
  	pz = jet->E() / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-jet->Eta())));
  	en = TMath::Sqrt(px * px + py * py + pz * pz);

  	aodjets->AddLast( new AliAODJet(px, py, pz, en) );
  	}
  jets->Delete();
  delete jets;
  
  // Order jets according to their pT .
  QSortTracks( *aodjets, 0, aodjets->GetEntriesFast() );
  
  // debug
  if (fDebug>3) AliInfo(Form(" %d Charged jets found\n",njets));
  
  return aodjets;
}


//____________________________________________________________________
void AliAnalyseUE::FillAvePartPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin  )
{

  // Fill average particle Pt of control regions
  
  // Max cone
  fHistos->FillHistogram("hRegionAvePartPtMaxVsEt", leadingE, ptMax);
  // Min cone
  fHistos->FillHistogram("hRegionAvePartPtMinVsEt", leadingE, ptMin);
  // MAke distributions for UE comparison with MB data
  fHistos->FillHistogram("hMinRegAvePt", ptMin);

}

//____________________________________________________________________
void AliAnalyseUE::FillMultRegion( Double_t leadingE, Double_t nTrackPtmax, Double_t nTrackPtmin, Double_t ptMin  )
{

  // Fill Nch multiplicity of control regions
  
  // Max cone
  fHistos->FillHistogram("hRegionMultMaxVsEt", leadingE, nTrackPtmax);
  fHistos->FillHistogram("hRegionMultMax",  nTrackPtmax);
  
  // Min cone
  fHistos->FillHistogram("hRegionMultMinVsEt", leadingE, nTrackPtmin );
  fHistos->FillHistogram("hRegionMultMin", nTrackPtmin);
  
  // MAke distributions for UE comparison with MB data
  fHistos->FillHistogram("hMinRegSumPtvsMult", nTrackPtmin,ptMin);

}

//____________________________________________________________________
void AliAnalyseUE::FillSumPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin  )
{
  // Fill sumPt of control regions
  
  // Max cone
  fHistos->FillHistogram("hRegionSumPtMaxVsEt", leadingE, ptMax);
  
  // Min cone
  fHistos->FillHistogram("hRegionSumPtMinVsEt", leadingE, ptMin);
  
  // MAke distributions for UE comparison with MB data
  fHistos->FillHistogram("hMinRegSumPt", ptMin);
  fHistos->FillHistogram("hMinRegSumPtJetPtBin", leadingE, ptMin);
  fHistos->FillHistogram("hMaxRegSumPtJetPtBin", leadingE, ptMax);

}

//-------------------------------------------------------------------
Int_t AliAnalyseUE::IsTrackInsideRegion(TVector3 *jetVect, TVector3 *partVect, Int_t conePosition) 
{  
  // return de region in delta phi
  // -1 negative delta phi 
  //  1 positive delta phi
  //  0 outside region
  static const Double_t k60rad  = 60.*TMath::Pi()/180.;
  static const Double_t k120rad = 120.*TMath::Pi()/180.;
 
  Int_t region = 0;
  if( fRegionType == 1 ) {
  	if( TMath::Abs(partVect->Eta()) > fTrackEtaCut ) return 0;
  	// transverse regions
  	if (jetVect[0].DeltaPhi(*partVect) < -k60rad && jetVect[0].DeltaPhi(*partVect) > -k120rad ) region = -1; //left
  	if (jetVect[0].DeltaPhi(*partVect) > k60rad && jetVect[0].DeltaPhi(*partVect) < k120rad ) region = 1;    //right
  	if (TMath::Abs(jetVect[0].DeltaPhi(*partVect)) < k60rad ) region = 2;    //forward
  	if (TMath::Abs(jetVect[0].DeltaPhi(*partVect)) > k120rad ) region = -2; //backward
    
  	} else if( fRegionType == 2 ) {
  	// Cone regions
  	Double_t deltaR = 0.;
    
  	TVector3 positVect,negatVect;
  	if (conePosition==1){
  		positVect.SetMagThetaPhi(1, 2.*atan(exp(-jetVect[0].Eta())), jetVect[0].Phi()+TMath::PiOver2());
  		negatVect.SetMagThetaPhi(1, 2.*atan(exp(-jetVect[0].Eta())), jetVect[0].Phi()-TMath::PiOver2());
  	}else if (conePosition==2){
  		if(fAnaType<2) AliFatal("Prevent error in Analysis type there might be only 1 jet. To avoid overflow better Correct UE config");
  		positVect.SetMagThetaPhi(1, 2.*atan(exp(-(jetVect[0].Eta()+jetVect[1].Eta())/2.)), jetVect[0].Phi()+TMath::PiOver2());
  		negatVect.SetMagThetaPhi(1, 2.*atan(exp(-(jetVect[0].Eta()+jetVect[1].Eta())/2.)), jetVect[0].Phi()-TMath::PiOver2());
  	}else if (conePosition==3){
  		if(fAnaType<2) AliFatal("Prevent error in Analysis type there might be only 1 jet. To avoid overflow better Correct UE config");
  		Double_t weightEta = jetVect[0].Eta() * jetVect[0].Pt()/(jetVect[0].Pt() + jetVect[1].Pt()) + 
  		jetVect[1].Eta() * jetVect[1].Pt()/(jetVect[0].Pt() + jetVect[1].Pt());
  		//Double_t weightEta = jetVect[0].Eta() * jetVect[0].Mag()/(jetVect[0].Mag() + jetVect[1].Mag()) + 
  		// jetVect[1].Eta() * jetVect[1].Mag()/(jetVect[0].Mag() + jetVect[1].Mag());
  		positVect.SetMagThetaPhi(1, 2.*atan(exp(-weightEta)), jetVect[0].Phi()+TMath::PiOver2());
  		negatVect.SetMagThetaPhi(1, 2.*atan(exp(-weightEta)), jetVect[0].Phi()-TMath::PiOver2());
  		}
  if (TMath::Abs(positVect.DeltaPhi(*partVect)) < fConeRadius ) { 
  	region = 1;  
  	deltaR = positVect.DrEtaPhi(*partVect);
  } else if (TMath::Abs(negatVect.DeltaPhi(*partVect)) < fConeRadius) { 
  	region = -1;  
  	deltaR = negatVect.DrEtaPhi(*partVect);
  	}
    
  if (deltaR > fConeRadius) region = 0;
    
  } else
  	AliError("Unknow region type");
        
  	return region;
  }



//-------------------------------------------------------------------
void  AliAnalyseUE::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
  // Sort array of TObjArray of tracks by Pt using a quicksort algorithm.
  
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;   // QSortTracks(j + 1, last);
    } else {
      QSortTracks(a, j + 1, last);
      last = j;        // QSortTracks(first, j);
    }
  }
}




//-------------------------------------------------------------------
void AliAnalyseUE::SetRegionArea(TVector3 *jetVect)
{
  // Set region area
  Double_t areaCorrFactor=0.;
  Double_t deltaEta = 0.;
  if (fRegionType==1) fAreaReg = 2.*fTrackEtaCut*60.*TMath::Pi()/180.;
  else if (fRegionType==2){ 
  	deltaEta = 0.9-TMath::Abs(jetVect[0].Eta());
    	if (deltaEta>fConeRadius) fAreaReg = TMath::Pi()*fConeRadius*fConeRadius;
    	else{
      		areaCorrFactor = fConeRadius*fConeRadius*TMath::ACos(deltaEta/fConeRadius) -
      		(deltaEta*fConeRadius)*TMath::Sqrt( 1. - deltaEta*deltaEta/(fConeRadius*fConeRadius));
      		fAreaReg=TMath::Pi()*fConeRadius*fConeRadius-areaCorrFactor;
    		}
  	}else AliWarning("Unknown Region Type");
  if (fDebug>10) AliInfo(Form("\n dEta=%5.3f Angle =%5.3f Region Area = %5.3f Corr Factor=%5.4f \n",deltaEta,TMath::ACos(deltaEta/fConeRadius),fAreaReg,areaCorrFactor));
}


//____________________________________________________________________
TObjArray*  AliAnalyseUE::SortChargedParticles()
{
  //  return an array with all charged particles ordered according to their pT .
  Int_t nTracks = fkAOD->GetNTracks();
  if( !nTracks ) return 0;
  TObjArray* tracks = new TObjArray(nTracks);

  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
  	AliAODTrack* part = fkAOD->GetTrack( ipart );
  	if( !part->TestFilterBit( fFilterBit ) ) continue; // track cut selection
  	if( !part->Charge() ) continue;
  	if( part->Pt() < fTrackPtCut ) continue;
  	tracks->AddLast( part );
  	}
  QSortTracks( *tracks, 0, tracks->GetEntriesFast() );

  nTracks = tracks->GetEntriesFast();
  if( !nTracks ) return 0;

  return tracks;
  }


//____________________________________________________________________
const Bool_t AliAnalyseUE::TrackMCSelected(Double_t charge, Double_t pT, Double_t eta, Int_t pdgCode){
  
  // MC track selection 
  Double_t const kMyTolerance = 0.0000001;
  //if (charge == 0. || charge == -99.) return kFALSE;
  if (charge < kMyTolerance || charge + 99 < kMyTolerance) return kFALSE;
        
  if ( fUseSingleCharge ) { // Charge selection
  	if ( fUsePositiveCharge && charge < 0.) return kFALSE; // keep Positives
  	if ( !fUsePositiveCharge && charge > 0.) return kFALSE; // keep Negatives
  	}
        
  //Kinematics cuts on particle
  if ((pT < fTrackPtCut) || (TMath::Abs(eta) > fTrackEtaCut )) return kFALSE;
        
  Bool_t isHadron = TMath::Abs(pdgCode)==211 ||
  	TMath::Abs(pdgCode)==2212 ||
  	TMath::Abs(pdgCode)==321;
        
  if ( fUseChargeHadrons && !isHadron ) return kFALSE;
        
  return kTRUE;

}

//____________________________________________________________________
const Bool_t AliAnalyseUE::TrackSelected(AliAODTrack* part){

  // Real track selection
  if ( !part->TestFilterBit(fFilterBit) ) return kFALSE; // track cut selection
  if ( !part->IsPrimaryCandidate()) return kFALSE; // reject whatever is not linked to collision point
  // PID Selection: Reject everything but hadrons
  Bool_t isHadron = part->GetMostProbablePID()==AliAODTrack::kPion || 
  	part->GetMostProbablePID()==AliAODTrack::kKaon || 
  	part->GetMostProbablePID()==AliAODTrack::kProton;
  if ( fUseChargeHadrons && !isHadron ) return kFALSE;
      
  if ( !part->Charge() ) return kFALSE; //Only charged
  if ( fUseSingleCharge ) { // Charge selection
  	if ( fUsePositiveCharge && part->Charge() < 0.) return kFALSE; // keep Positives
  	if ( !fUsePositiveCharge && part->Charge() > 0.) return kFALSE; // keep Negatives
 	} 
    
  if ( part->Pt() < fTrackPtCut ) return kFALSE;
  if( TMath::Abs(part->Eta()) > fTrackEtaCut ) return kFALSE;

  return kTRUE;
}


//____________________________________________________________________
void AliAnalyseUE::WriteSettings(){

  // Print analysis settings
  if (fDebug>5){
	AliInfo("All Analysis Settings in Saved Tree");
	fSettingsTree->Scan();
  	}

}
