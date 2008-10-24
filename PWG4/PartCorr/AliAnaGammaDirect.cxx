/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes iGetEntriesFast(s hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: AliAnaGammaDirect.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
// Class for the prompt gamma analysis, isolation cut
//
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TParticle.h>
#include <TH2.h>
#include <TList.h>
#include "Riostream.h"

// --- Analysis system --- 
#include "AliAnaGammaDirect.h" 
#include "AliLog.h"
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"

ClassImp(AliAnaGammaDirect)
  
//____________________________________________________________________________
  AliAnaGammaDirect::AliAnaGammaDirect() : 
    AliAnaPartCorrBaseClass(), fDetector(""), fMakeIC(0),  fReMakeIC(0), 
    fMakeSeveralIC(0), fMakeInvMass(0),
    fhPtGamma(0),fhPhiGamma(0),fhEtaGamma(0), fhConeSumPt(0),
    //Several IC
    fNCones(0),fNPtThresFrac(0), fConeSizes(),  fPtThresholds(),  fPtFractions(), 
    fhPtThresIsolated(), fhPtFracIsolated(),  fhPtSumIsolated(),
    //MC
    fhPtPrompt(0),fhPhiPrompt(0),fhEtaPrompt(0), 
    fhPtThresIsolatedPrompt(), fhPtFracIsolatedPrompt(),  fhPtSumIsolatedPrompt(),
    fhPtFragmentation(0),fhPhiFragmentation(0),fhEtaFragmentation(0), 
    fhPtThresIsolatedFragmentation(), fhPtFracIsolatedFragmentation(),  fhPtSumIsolatedFragmentation(),
    fhPtPi0Decay(0),fhPhiPi0Decay(0),fhEtaPi0Decay(0), 
    fhPtThresIsolatedPi0Decay(), fhPtFracIsolatedPi0Decay(),  fhPtSumIsolatedPi0Decay(),
    fhPtOtherDecay(0),fhPhiOtherDecay(0),fhEtaOtherDecay(0), 
    fhPtThresIsolatedOtherDecay(), fhPtFracIsolatedOtherDecay(),  fhPtSumIsolatedOtherDecay(),
    fhPtConversion(0),fhPhiConversion(0),fhEtaConversion(0), 
    fhPtThresIsolatedConversion(), fhPtFracIsolatedConversion(),  fhPtSumIsolatedConversion(),
    fhPtUnknown(0),fhPhiUnknown(0),fhEtaUnknown(0), 
    fhPtThresIsolatedUnknown(), fhPtFracIsolatedUnknown(),  fhPtSumIsolatedUnknown()
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

  for(Int_t i = 0; i < 5 ; i++){ 
    fConeSizes[i] = 0 ; 
    fhPtSumIsolated[i] = 0 ;  
    
    fhPtSumIsolatedPrompt[i] = 0 ;  
    fhPtSumIsolatedFragmentation[i] = 0 ;  
    fhPtSumIsolatedPi0Decay[i] = 0 ;  
    fhPtSumIsolatedOtherDecay[i] = 0 ;  
    fhPtSumIsolatedConversion[i] = 0 ;  
    fhPtSumIsolatedUnknown[i] = 0 ;  
    
    for(Int_t j = 0; j < 5 ; j++){ 
      fhPtThresIsolated[i][j] = 0 ;  
      fhPtFracIsolated[i][j] = 0 ; 
      
      fhPtThresIsolatedPrompt[i][j] = 0 ;  
      fhPtThresIsolatedFragmentation[i][j] = 0 ; 
      fhPtThresIsolatedPi0Decay[i][j] = 0 ;  
      fhPtThresIsolatedOtherDecay[i][j] = 0 ;  
      fhPtThresIsolatedConversion[i][j] = 0 ;  
      fhPtThresIsolatedUnknown[i][j] = 0 ;  
  
      fhPtFracIsolatedPrompt[i][j] = 0 ;  
      fhPtFracIsolatedFragmentation[i][j] = 0 ;  
      fhPtFracIsolatedPi0Decay[i][j] = 0 ;  
      fhPtFracIsolatedOtherDecay[i][j] = 0 ;  
      fhPtFracIsolatedConversion[i][j] = 0 ;
      fhPtFracIsolatedUnknown[i][j] = 0 ;  
 
    }  
  } 
  
  for(Int_t i = 0; i < 5 ; i++){ 
    fPtFractions[i]=  0 ; 
    fPtThresholds[i]= 0 ; 
  } 


}

//____________________________________________________________________________
AliAnaGammaDirect::AliAnaGammaDirect(const AliAnaGammaDirect & g) : 
  AliAnaPartCorrBaseClass(g), fDetector(g.fDetector),
  fMakeIC(g.fMakeIC),   fReMakeIC(g.fReMakeIC), 
  fMakeSeveralIC(g.fMakeSeveralIC),  fMakeInvMass(g.fMakeInvMass),
  fhPtGamma(g.fhPtGamma),fhPhiGamma(g.fhPhiGamma),
  fhEtaGamma(g.fhEtaGamma), fhConeSumPt(g.fhConeSumPt),  
  //Several IC
  fNCones(g.fNCones),fNPtThresFrac(g.fNPtThresFrac), fConeSizes(), fPtThresholds(),  fPtFractions(), 
  fhPtThresIsolated(),  fhPtFracIsolated(), fhPtSumIsolated(),
  //MC
  fhPtPrompt(g.fhPtPrompt),fhPhiPrompt(g.fhPhiPrompt),fhEtaPrompt(g.fhEtaPrompt), 
  fhPtThresIsolatedPrompt(), fhPtFracIsolatedPrompt(),  fhPtSumIsolatedPrompt(),
  fhPtFragmentation(g.fhPtFragmentation),fhPhiFragmentation(g.fhPhiFragmentation),fhEtaFragmentation(g.fhEtaFragmentation), 
  fhPtThresIsolatedFragmentation(), fhPtFracIsolatedFragmentation(),  fhPtSumIsolatedFragmentation(),
  fhPtPi0Decay(g.fhPtPi0Decay),fhPhiPi0Decay(g.fhPhiPi0Decay),fhEtaPi0Decay(g.fhEtaPi0Decay), 
  fhPtThresIsolatedPi0Decay(), fhPtFracIsolatedPi0Decay(),  fhPtSumIsolatedPi0Decay(),
  fhPtOtherDecay(g.fhPtOtherDecay),fhPhiOtherDecay(g.fhPhiOtherDecay),fhEtaOtherDecay(g.fhEtaOtherDecay), 
  fhPtThresIsolatedOtherDecay(), fhPtFracIsolatedOtherDecay(),  fhPtSumIsolatedOtherDecay(),
  fhPtConversion(g. fhPtConversion),fhPhiConversion(g.fhPhiConversion),fhEtaConversion(g.fhEtaConversion), 
  fhPtThresIsolatedConversion(), fhPtFracIsolatedConversion(),  fhPtSumIsolatedConversion(),
  fhPtUnknown(g.fhPtUnknown),fhPhiUnknown(g.fhPhiUnknown),fhEtaUnknown(g.fhEtaUnknown), 
  fhPtThresIsolatedUnknown(), fhPtFracIsolatedUnknown(),  fhPtSumIsolatedUnknown()
{
  // cpy ctor
  
  //Several IC
  for(Int_t i = 0; i < fNCones ; i++){
    fConeSizes[i] =  g.fConeSizes[i];
    fhPtSumIsolated[i] = g.fhPtSumIsolated[i]; 

   fhPtSumIsolatedPrompt[i] = g.fhPtSumIsolatedPrompt[i]; 
   fhPtSumIsolatedFragmentation[i] = g.fhPtSumIsolatedFragmentation[i]; 
   fhPtSumIsolatedPi0Decay[i] = g.fhPtSumIsolatedPi0Decay[i]; 
   fhPtSumIsolatedOtherDecay[i] = g.fhPtSumIsolatedOtherDecay[i]; 
   fhPtSumIsolatedConversion[i] = g.fhPtSumIsolatedConversion[i]; 
   fhPtSumIsolatedUnknown[i] = g.fhPtSumIsolatedUnknown[i]; 

    for(Int_t j = 0; j < fNPtThresFrac ; j++){
      fhPtThresIsolated[i][j] = g.fhPtThresIsolated[i][j]; 
      fhPtFracIsolated[i][j] = g.fhPtFracIsolated[i][j];

      fhPtThresIsolatedPrompt[i][j] = g.fhPtThresIsolatedPrompt[i][j]; 
      fhPtThresIsolatedFragmentation[i][j] = g.fhPtThresIsolatedFragmentation[i][j]; 
      fhPtThresIsolatedPi0Decay[i][j] = g.fhPtThresIsolatedPi0Decay[i][j]; 
      fhPtThresIsolatedOtherDecay[i][j] = g.fhPtThresIsolatedOtherDecay[i][j]; 
      fhPtThresIsolatedConversion[i][j] = g.fhPtThresIsolatedConversion[i][j]; 
      fhPtThresIsolatedUnknown[i][j] = g.fhPtThresIsolatedUnknown[i][j]; 
 
      fhPtFracIsolatedPrompt[i][j] = g.fhPtFracIsolatedPrompt[i][j]; 
      fhPtFracIsolatedFragmentation[i][j] = g.fhPtFracIsolatedFragmentation[i][j]; 
      fhPtFracIsolatedPi0Decay[i][j] = g.fhPtFracIsolatedPi0Decay[i][j]; 
      fhPtFracIsolatedOtherDecay[i][j] = g.fhPtFracIsolatedOtherDecay[i][j]; 
      fhPtFracIsolatedConversion[i][j] = g.fhPtFracIsolatedConversion[i][j]; 
      fhPtFracIsolatedUnknown[i][j] = g.fhPtFracIsolatedUnknown[i][j]; 

    } 
  }
  
  for(Int_t i = 0; i < fNPtThresFrac ; i++){
    fPtFractions[i]=   g.fPtFractions[i];
    fPtThresholds[i]=   g.fPtThresholds[i];
  }

}

//_________________________________________________________________________
AliAnaGammaDirect & AliAnaGammaDirect::operator = (const AliAnaGammaDirect & g)
{
  // assignment operator
  
  if(&g == this) return *this;

  fMakeIC = g.fMakeIC ;
  fReMakeIC = g.fReMakeIC ;
  fMakeSeveralIC = g.fMakeSeveralIC ;
  fMakeInvMass = g.fMakeInvMass ;
  fDetector = g.fDetector ;

  fhPtGamma = g.fhPtGamma ; 
  fhPhiGamma = g.fhPhiGamma ;
  fhEtaGamma = g.fhEtaGamma ;
  fhConeSumPt = g.fhConeSumPt ;

  fhPtPrompt = g.fhPtPrompt;
  fhPhiPrompt = g.fhPhiPrompt;
  fhEtaPrompt = g.fhEtaPrompt; 
  fhPtFragmentation = g.fhPtFragmentation;
  fhPhiFragmentation = g.fhPhiFragmentation;
  fhEtaFragmentation = g.fhEtaFragmentation; 
  fhPtPi0Decay = g.fhPtPi0Decay;
  fhPhiPi0Decay = g.fhPhiPi0Decay;
  fhEtaPi0Decay = g.fhEtaPi0Decay; 
  fhPtOtherDecay = g.fhPtOtherDecay;
  fhPhiOtherDecay = g.fhPhiOtherDecay;
  fhEtaOtherDecay = g.fhEtaOtherDecay; 
  fhPtConversion = g. fhPtConversion;
  fhPhiConversion = g.fhPhiConversion;
  fhEtaConversion = g.fhEtaConversion; 
  fhPtUnknown = g.fhPtUnknown;
  fhPhiUnknown = g.fhPhiUnknown;
  fhEtaUnknown = g.fhEtaUnknown; 

  //Several IC
  fNCones = g.fNCones ;
  fNPtThresFrac = g.fNPtThresFrac ; 
   
  for(Int_t i = 0; i < fNCones ; i++){
    fConeSizes[i] =  g.fConeSizes[i];
    fhPtSumIsolated[i] = g.fhPtSumIsolated[i] ;

   fhPtSumIsolatedPrompt[i] = g.fhPtSumIsolatedPrompt[i]; 
   fhPtSumIsolatedFragmentation[i] = g.fhPtSumIsolatedFragmentation[i]; 
   fhPtSumIsolatedPi0Decay[i] = g.fhPtSumIsolatedPi0Decay[i]; 
   fhPtSumIsolatedOtherDecay[i] = g.fhPtSumIsolatedOtherDecay[i]; 
   fhPtSumIsolatedConversion[i] = g.fhPtSumIsolatedConversion[i]; 
   fhPtSumIsolatedUnknown[i] = g.fhPtSumIsolatedUnknown[i]; 

    for(Int_t j = 0; j < fNPtThresFrac ; j++){
      fhPtThresIsolated[i][j] = g.fhPtThresIsolated[i][j] ;
      fhPtFracIsolated[i][j] = g.fhPtFracIsolated[i][j] ;

      fhPtThresIsolatedPrompt[i][j] = g.fhPtThresIsolatedPrompt[i][j]; 
      fhPtThresIsolatedFragmentation[i][j] = g.fhPtThresIsolatedFragmentation[i][j]; 
      fhPtThresIsolatedPi0Decay[i][j] = g.fhPtThresIsolatedPi0Decay[i][j]; 
      fhPtThresIsolatedOtherDecay[i][j] = g.fhPtThresIsolatedOtherDecay[i][j]; 
      fhPtThresIsolatedConversion[i][j] = g.fhPtThresIsolatedConversion[i][j]; 
      fhPtThresIsolatedUnknown[i][j] = g.fhPtThresIsolatedUnknown[i][j]; 
 
      fhPtFracIsolatedPrompt[i][j] = g.fhPtFracIsolatedPrompt[i][j]; 
      fhPtFracIsolatedFragmentation[i][j] = g.fhPtFracIsolatedFragmentation[i][j]; 
      fhPtFracIsolatedPi0Decay[i][j] = g.fhPtFracIsolatedPi0Decay[i][j]; 
      fhPtFracIsolatedOtherDecay[i][j] = g.fhPtFracIsolatedOtherDecay[i][j]; 
      fhPtFracIsolatedConversion[i][j] = g.fhPtFracIsolatedConversion[i][j]; 
      fhPtFracIsolatedUnknown[i][j] = g.fhPtFracIsolatedUnknown[i][j]; 

    }
  }
  
  for(Int_t i = 0; i < fNPtThresFrac ; i++){
    fPtThresholds[i]=   g.fPtThresholds[i];
    fPtFractions[i]=   g.fPtFractions[i];
  }

  return *this;
  
}

//____________________________________________________________________________
AliAnaGammaDirect::~AliAnaGammaDirect() 
{
  //dtor
  //do not delete histograms

 delete [] fConeSizes ; 
 delete [] fPtThresholds ; 
 delete [] fPtFractions ; 

}
//_________________________________________________________________________
Bool_t AliAnaGammaDirect::CheckInvMass(const Int_t icalo,const TLorentzVector mom, Double_t *vertex, TClonesArray * pl){
  //Search if there is a companion decay photon to the candidate 
  // and discard it in such case
  TLorentzVector mom2 ;
  for(Int_t jcalo = 0; jcalo < pl->GetEntriesFast(); jcalo++){
    if(icalo==jcalo) continue ;
    AliAODCaloCluster * calo =  dynamic_cast<AliAODCaloCluster*> (pl->At(jcalo));
  
    //Cluster selection, not charged, with photon id and in fidutial cut, fill TLorentz
    if(!SelectCluster(calo, vertex, mom2)) continue ;

    //Select good pair (good phit, pt cuts, aperture and invariant mass)
    if(GetNeutralMesonSelection()->SelectPair(mom, mom2)){
      if(GetDebug()>1)printf("Selected gamma pair: pt %f, phi %f, eta%f",(mom+mom2).Pt(), (mom+mom2).Phi(), (mom+mom2).Eta());
      return kTRUE ;
    }
  }//loop

  return kFALSE;

}
//_________________________________________________________________________
Int_t AliAnaGammaDirect::CheckOrigin(const Int_t label){
  //Play with the MC stack if available
  //Check origin of the candidates, good for PYTHIA
  
  AliStack * stack =  GetMCStack() ;
  if(!stack) AliFatal("Stack is not available, check analysis settings in configuration file, STOP!!");
  
  if(label >= 0 && label <  stack->GetNtrack()){
    //Mother
    TParticle * mom = stack->Particle(label);
    Int_t mPdg = TMath::Abs(mom->GetPdgCode());
    Int_t mStatus =  mom->GetStatusCode() ;
    Int_t iParent =  mom->GetFirstMother() ;
    if(label < 8 ) printf("Mother is parton %d\n",iParent);
    
    //GrandParent
    TParticle * parent = new TParticle ;
    Int_t pPdg = -1;
    Int_t pStatus =-1;
    if(iParent > 0){
      parent = stack->Particle(iParent);
      pPdg = TMath::Abs(parent->GetPdgCode());
      pStatus = parent->GetStatusCode();  
    }
    else
      printf("Parent with label %d\n",iParent);
    
    //return tag
    if(mPdg == 22){
      if(mStatus == 1){
	if(iParent < 8) {
	  if(pPdg == 22) return kPrompt;
	  else  return kFragmentation;
	}
	else if(pStatus == 11){
	  if(pPdg == 111) return kPi0Decay ;
	  else if (pPdg == 321)  return kEtaDecay ;
	  else  return kOtherDecay ;
	}
      }//Status 1 : Pythia generated
      else if(mStatus == 0){
	if(pPdg ==22 || pPdg ==11) return kConversion ;
	if(pPdg == 111) return kPi0Decay ;
	else if (pPdg == 221)  return kEtaDecay ;
	else  return kOtherDecay ;
      }//status 0 : geant generated
    }//Mother gamma
    else if(mPdg == 111)  return kPi0 ;
    else if(mPdg == 221)  return kEta ;
    else if(mPdg ==11){
      if(mStatus == 0) return kConversion ;
      else return kElectron ;
    }
    else return kUnknown;
  }//Good label value
  else{
    if(label < 0 ) printf("*** bad label or no stack ***:  label %d \n", label);
    if(label >=  stack->GetNtrack()) printf("*** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
    return kUnknown;
  }//Bad label
  
  return kUnknown;
  
}

//________________________________________________________________________
TList *  AliAnaGammaDirect::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("DirectGammaHistos") ; 

  //Histograms of highest gamma identified in Event
  fhPtGamma  = new TH1F("hPtGamma","Number of #gamma over calorimeter",240,0,120); 
  fhPtGamma->SetYTitle("N");
  fhPtGamma->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhPtGamma) ; 
  
  fhPhiGamma  = new TH2F
    ("hPhiGamma","#phi_{#gamma}",200,0,120,200,0,7); 
  fhPhiGamma->SetYTitle("#phi");
  fhPhiGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhiGamma) ; 
  
  fhEtaGamma  = new TH2F
    ("hEtaGamma","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
  fhEtaGamma->SetYTitle("#eta");
  fhEtaGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhEtaGamma) ;

  if(!fMakeSeveralIC){
    fhConeSumPt  = new TH2F
      ("hConePtSum","#Sigma p_{T}  in cone ",200,0,120,100,0,100);
    fhConeSumPt->SetYTitle("#Sigma p_{T}");
    fhConeSumPt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhConeSumPt) ;
  }
  
  if(IsDataMC()){
    
    fhPtPrompt  = new TH1F("hPtPrompt","Number of #gamma over calorimeter",240,0,120); 
    fhPtPrompt->SetYTitle("N");
    fhPtPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPrompt) ; 
    
    fhPhiPrompt  = new TH2F
      ("hPhiPrompt","#phi_{#gamma}",200,0,120,200,0,7); 
    fhPhiPrompt->SetYTitle("#phi");
    fhPhiPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiPrompt) ; 
    
    fhEtaPrompt  = new TH2F
      ("hEtaPrompt","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
    fhEtaPrompt->SetYTitle("#eta");
    fhEtaPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaPrompt) ;
    
    fhPtFragmentation  = new TH1F("hPtFragmentation","Number of #gamma over calorimeter",240,0,120); 
    fhPtFragmentation->SetYTitle("N");
    fhPtFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtFragmentation) ; 
    
    fhPhiFragmentation  = new TH2F
      ("hPhiFragmentation","#phi_{#gamma}",200,0,120,200,0,7); 
    fhPhiFragmentation->SetYTitle("#phi");
    fhPhiFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiFragmentation) ; 
    
    fhEtaFragmentation  = new TH2F
      ("hEtaFragmentation","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
    fhEtaFragmentation->SetYTitle("#eta");
    fhEtaFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaFragmentation) ;
    
    fhPtPi0Decay  = new TH1F("hPtPi0Decay","Number of #gamma over calorimeter",240,0,120); 
    fhPtPi0Decay->SetYTitle("N");
    fhPtPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPi0Decay) ; 
    
    fhPhiPi0Decay  = new TH2F
      ("hPhiPi0Decay","#phi_{#gamma}",200,0,120,200,0,7); 
    fhPhiPi0Decay->SetYTitle("#phi");
    fhPhiPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiPi0Decay) ; 
    
    fhEtaPi0Decay  = new TH2F
      ("hEtaPi0Decay","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
    fhEtaPi0Decay->SetYTitle("#eta");
    fhEtaPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaPi0Decay) ;
    
    fhPtOtherDecay  = new TH1F("hPtOtherDecay","Number of #gamma over calorimeter",240,0,120); 
    fhPtOtherDecay->SetYTitle("N");
    fhPtOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtOtherDecay) ; 
    
    fhPhiOtherDecay  = new TH2F
      ("hPhiOtherDecay","#phi_{#gamma}",200,0,120,200,0,7); 
    fhPhiOtherDecay->SetYTitle("#phi");
    fhPhiOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiOtherDecay) ; 
    
    fhEtaOtherDecay  = new TH2F
      ("hEtaOtherDecay","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
    fhEtaOtherDecay->SetYTitle("#eta");
    fhEtaOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaOtherDecay) ;
    
    fhPtConversion  = new TH1F("hPtConversion","Number of #gamma over calorimeter",240,0,120); 
    fhPtConversion->SetYTitle("N");
    fhPtConversion->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtConversion) ; 
    
    fhPhiConversion  = new TH2F
      ("hPhiConversion","#phi_{#gamma}",200,0,120,200,0,7); 
    fhPhiConversion->SetYTitle("#phi");
    fhPhiConversion->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiConversion) ; 
    
    fhEtaConversion  = new TH2F
      ("hEtaConversion","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
    fhEtaConversion->SetYTitle("#eta");
    fhEtaConversion->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaConversion) ;
    
    fhPtUnknown  = new TH1F("hPtUnknown","Number of #gamma over calorimeter",240,0,120); 
    fhPtUnknown->SetYTitle("N");
    fhPtUnknown->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtUnknown) ; 
    
    fhPhiUnknown  = new TH2F
      ("hPhiUnknown","#phi_{#gamma}",200,0,120,200,0,7); 
    fhPhiUnknown->SetYTitle("#phi");
    fhPhiUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiUnknown) ; 
    
    fhEtaUnknown  = new TH2F
      ("hEtaUnknown","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
    fhEtaUnknown->SetYTitle("#eta");
    fhEtaUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaUnknown) ;
  }//Histos with MC
  
  if(fMakeSeveralIC){
    char name[128];
    char title[128];
    for(Int_t icone = 0; icone<fNCones; icone++){
      sprintf(name,"hPtSumIsolated_Cone_%d",icone);
      sprintf(title,"Candidate cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
      fhPtSumIsolated[icone]  = new TH2F(name, title,240,0,120,120,0,10);
      fhPtSumIsolated[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
      fhPtSumIsolated[icone]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtSumIsolated[icone]) ; 
    
      if(IsDataMC()){
	sprintf(name,"hPtSumIsolatedPrompt_Cone_%d",icone);
	sprintf(title,"Candidate Prompt cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
	fhPtSumIsolatedPrompt[icone]  = new TH2F(name, title,240,0,120,120,0,10);
	fhPtSumIsolatedPrompt[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
	fhPtSumIsolatedPrompt[icone]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtSumIsolatedPrompt[icone]) ; 

	sprintf(name,"hPtSumIsolatedFragmentation_Cone_%d",icone);
	sprintf(title,"Candidate Fragmentation cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
	fhPtSumIsolatedFragmentation[icone]  = new TH2F(name, title,240,0,120,120,0,10);
	fhPtSumIsolatedFragmentation[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
	fhPtSumIsolatedFragmentation[icone]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtSumIsolatedFragmentation[icone]) ; 

	sprintf(name,"hPtSumIsolatedPi0Decay_Cone_%d",icone);
	sprintf(title,"Candidate Pi0Decay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
	fhPtSumIsolatedPi0Decay[icone]  = new TH2F(name, title,240,0,120,120,0,10);
	fhPtSumIsolatedPi0Decay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
	fhPtSumIsolatedPi0Decay[icone]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtSumIsolatedPi0Decay[icone]) ; 

	sprintf(name,"hPtSumIsolatedOtherDecay_Cone_%d",icone);
	sprintf(title,"Candidate OtherDecay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
	fhPtSumIsolatedOtherDecay[icone]  = new TH2F(name, title,240,0,120,120,0,10);
	fhPtSumIsolatedOtherDecay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
	fhPtSumIsolatedOtherDecay[icone]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtSumIsolatedOtherDecay[icone]) ; 

	sprintf(name,"hPtSumIsolatedConversion_Cone_%d",icone);
	sprintf(title,"Candidate Conversion cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
	fhPtSumIsolatedConversion[icone]  = new TH2F(name, title,240,0,120,120,0,10);
	fhPtSumIsolatedConversion[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
	fhPtSumIsolatedConversion[icone]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtSumIsolatedConversion[icone]) ; 

	sprintf(name,"hPtSumIsolatedUnknown_Cone_%d",icone);
	sprintf(title,"Candidate Unknown cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
	fhPtSumIsolatedUnknown[icone]  = new TH2F(name, title,240,0,120,120,0,10);
	fhPtSumIsolatedUnknown[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
	fhPtSumIsolatedUnknown[icone]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtSumIsolatedUnknown[icone]) ; 

      }//Histos with MC

      for(Int_t ipt = 0; ipt<fNPtThresFrac;ipt++){ 
	sprintf(name,"hPtThresIsol_Cone_%d_Pt%d",icone,ipt);
	sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,240,0,120);
	fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 

	sprintf(name,"hPtFracIsol_Cone_%d_Pt%d",icone,ipt);
	sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	fhPtFracIsolated[icone][ipt]  = new TH1F(name, title,240,0,120);
	fhPtFracIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtFracIsolated[icone][ipt]) ; 

	if(IsDataMC()){
	  sprintf(name,"hPtThresIsolPrompt_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtThresIsolatedPrompt[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtThresIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtThresIsolatedPrompt[icone][ipt]) ; 
	  
	  sprintf(name,"hPtFracIsolPrompt_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtFracIsolatedPrompt[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtFracIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtFracIsolatedPrompt[icone][ipt]) ; 

	  sprintf(name,"hPtThresIsolFragmentation_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtThresIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtThresIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtThresIsolatedFragmentation[icone][ipt]) ; 
	  
	  sprintf(name,"hPtFracIsolFragmentation_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtFracIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtFracIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtFracIsolatedFragmentation[icone][ipt]) ; 

	  sprintf(name,"hPtThresIsolPi0Decay_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtThresIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtThresIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtThresIsolatedPi0Decay[icone][ipt]) ; 
	  
	  sprintf(name,"hPtFracIsolPi0Decay_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtFracIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtFracIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtFracIsolatedPi0Decay[icone][ipt]) ; 

	  sprintf(name,"hPtThresIsolOtherDecay_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtThresIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtThresIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtThresIsolatedOtherDecay[icone][ipt]) ; 
	  
	  sprintf(name,"hPtFracIsolOtherDecay_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtFracIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtFracIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtFracIsolatedOtherDecay[icone][ipt]) ;

	  sprintf(name,"hPtThresIsolConversion_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtThresIsolatedConversion[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtThresIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtThresIsolatedConversion[icone][ipt]) ; 
	  
	  sprintf(name,"hPtFracIsolConversion_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtFracIsolatedConversion[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtFracIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtFracIsolatedConversion[icone][ipt]) ;

	  sprintf(name,"hPtThresIsolUnknown_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Unknown p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtThresIsolatedUnknown[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtThresIsolatedUnknown[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtThresIsolatedUnknown[icone][ipt]) ; 
	  
	  sprintf(name,"hPtFracIsolUnknown_Cone_%d_Pt%d",icone,ipt);
	  sprintf(title,"Isolated candidate Unknown p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
	  fhPtFracIsolatedUnknown[icone][ipt]  = new TH1F(name, title,240,0,120);
	  fhPtFracIsolatedUnknown[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
	  outputContainer->Add(fhPtFracIsolatedUnknown[icone][ipt]) ;  
 
	}//Histos with MC

      }//icone loop
    }//ipt loop
  }

  //Keep neutral meson selection histograms if requiered
  //Setting done in AliNeutralMesonSelection
  if(fMakeInvMass && GetNeutralMesonSelection()){
    TList * nmsHistos = GetNeutralMesonSelection()->GetCreateOutputObjects() ;
    cout<<"NMSHistos "<< nmsHistos<<endl;
    if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
      for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
  }

  return outputContainer ;

}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  //Search for the isolated photon in fDetector with pt > GetMinPt()

  //Fill CaloClusters if working with ESDs
  //if(GetReader()->GetDataType() == AliCaloTrackReader::kESD) ConnectAODCaloClusters(); 

  Int_t n = 0, nfrac = 0;
  Bool_t isolated = kFALSE ; 
  Float_t coneptsum = 0 ;
  TClonesArray * pl = new TClonesArray; 

  //Get vertex for photon momentum calculation
  Double_t vertex[]={0,0,0} ; //vertex ;
  if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);

  //Select the detector of the photon
  if(fDetector == "PHOS")
    pl = GetAODPHOS();
  else if (fDetector == "EMCAL")
    pl = GetAODEMCAL();
  //cout<<"Number of entries "<<pl->GetEntriesFast()<<endl;
  
  //Fill AODCaloClusters and AODParticleCorrelation with PHOS aods
  TLorentzVector mom ;
  for(Int_t icalo = 0; icalo < pl->GetEntriesFast(); icalo++){
    AliAODCaloCluster * calo =  dynamic_cast<AliAODCaloCluster*> (pl->At(icalo));
  
    //Cluster selection, not charged, with photon id and in fidutial cut
    if(!SelectCluster(calo,vertex,mom)) continue ;
    
    //If too small pt, skip
    if(mom.Pt() < GetMinPt() || mom.Pt() > GetMaxPt() ) continue ; 

    //Play with the MC stack if available
    Int_t tag =-1;
    if(IsDataMC()){
      //Check origin of the candidates
      tag = CheckOrigin(calo->GetLabel(0));
      if(GetDebug() > 0) printf("Origin of candidate %d\n",tag);
    }//Work with stack also   

    //Check invariant mass
    Bool_t decay = kFALSE ;
    if(fMakeInvMass) decay = CheckInvMass(icalo,mom,vertex,pl);
    if(decay) continue ;

    //Create AOD for analysis
    AliAODParticleCorrelation ph = AliAODParticleCorrelation(mom);
    ph.SetLabel(calo->GetLabel(0));
    ph.SetPdg(AliCaloPID::kPhoton);
    ph.SetDetector(fDetector);
    ph.SetTag(tag);  
    if(fMakeIC){
      n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
      GetIsolationCut()->MakeIsolationCut((TSeqCollection*)GetAODCTS(), (TSeqCollection*)pl, 
					  vertex, kTRUE, &ph,icalo,-1,
					  n,nfrac,coneptsum, isolated);
      if(isolated){
	//cout<<"Isolated : E "<<mom.E()<<" ; Phi"<<mom.Phi()<< " ; Eta "<<mom.Eta()<<endl;
	AddAODParticleCorrelation(ph);
      }
    }
    else //Fill all if not isolation requested
      AddAODParticleCorrelation(ph);

  }//loop
  
  if(GetDebug() > 1) printf("End fill AODs ");  
  
}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  Int_t n = 0, nfrac = 0;
  Bool_t isolated = kFALSE ; 
  Float_t coneptsum = 0 ;
  //cout<<"MakeAnalysisFillHistograms"<<endl;

  //Get vertex for photon momentum calculation
  Double_t v[]={0,0,0} ; //vertex ;
  if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(v);

  //Loop on stored AOD photons
  Int_t naod = GetAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("histo aod branch entries %d\n", naod);
  for(Int_t iaod = 0; iaod < naod ; iaod++){
   AliAODParticleCorrelation* ph =  dynamic_cast<AliAODParticleCorrelation*> (GetAODBranch()->At(iaod));
   Int_t pdg = ph->GetPdg();

   //Only isolate photons in detector fDetector
   if(pdg != AliCaloPID::kPhoton || ph->GetDetector() != fDetector) continue;
   
   if(fMakeSeveralIC) {
     //Analysis of multiple IC at same time
     MakeSeveralICAnalysis(ph,v);
     continue;
   }
   else if(fReMakeIC){
     //In case a more strict IC is needed in the produced AOD
     n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
     GetIsolationCut()->MakeIsolationCut((TSeqCollection*)ph->GetRefIsolationConeTracks(), 
					 (TSeqCollection*)ph->GetRefIsolationConeClusters(), 
					 v, kFALSE, ph,-1,-1,
					 n,nfrac,coneptsum, isolated);
   }

   //Fill histograms (normal case)
   if(!fReMakeIC || (fReMakeIC && isolated)){
     
     //printf("Isolated Gamma: pt %f, phi %f, eta %f", ph->Pt(),ph->Phi(),ph->Eta()) ;
     
     //Fill prompt gamma histograms 
     Float_t ptcluster = ph->Pt();
     Float_t phicluster = ph->Phi();
     Float_t etacluster = ph->Eta();
     
     fhPtGamma   ->Fill(ptcluster);
     fhPhiGamma ->Fill(ptcluster,phicluster);
     fhEtaGamma ->Fill(ptcluster,etacluster);
     fhConeSumPt->Fill(ptcluster,coneptsum);

     if(IsDataMC()){
       Int_t tag =ph->GetTag();
           
       if(tag == kPrompt){
	 fhPtPrompt   ->Fill(ptcluster);
	 fhPhiPrompt ->Fill(ptcluster,phicluster);
	 fhEtaPrompt ->Fill(ptcluster,etacluster);
       }
       else if(tag==kFragmentation)
	 {
	   fhPtFragmentation   ->Fill(ptcluster);
	   fhPhiFragmentation ->Fill(ptcluster,phicluster);
	   fhEtaFragmentation ->Fill(ptcluster,etacluster);
	 }
       else if(tag==kPi0Decay)
	 {
	   fhPtPi0Decay   ->Fill(ptcluster);
	   fhPhiPi0Decay ->Fill(ptcluster,phicluster);
	   fhEtaPi0Decay ->Fill(ptcluster,etacluster);
	 }
       else if(tag==kEtaDecay || tag==kOtherDecay)
	 {
	   fhPtOtherDecay   ->Fill(ptcluster);
	   fhPhiOtherDecay ->Fill(ptcluster,phicluster);
	   fhEtaOtherDecay ->Fill(ptcluster,etacluster);
	 }
       else if(tag==kConversion)
	 {
	   fhPtConversion   ->Fill(ptcluster);
	   fhPhiConversion ->Fill(ptcluster,phicluster);
	   fhEtaConversion ->Fill(ptcluster,etacluster);
	 }
       else{
	 
	 fhPtUnknown   ->Fill(ptcluster);
	 fhPhiUnknown ->Fill(ptcluster,phicluster);
	 fhEtaUnknown ->Fill(ptcluster,etacluster);
       }
     }//Histograms with MC
     
   }
  }// aod loop
  
}

//____________________________________________________________________________
void AliAnaGammaDirect::InitParameters()
{
  
  //Initialize the parameters of the analysis.

  fDetector = "PHOS" ;
  fMakeIC = kTRUE;
  fReMakeIC = kFALSE ;
  fMakeSeveralIC = kFALSE ;
  fMakeInvMass = kFALSE ;

 //----------- Several IC-----------------
  fNCones           = 5 ; 
  fNPtThresFrac         = 6 ; 
  fConeSizes[0] = 0.1; fConeSizes[1] = 0.2; fConeSizes[2] = 0.3; fConeSizes[3] = 0.4; fConeSizes[4] = 0.5;
  fPtThresholds[0]=0.; fPtThresholds[1]=1.; fPtThresholds[2]=2.; fPtThresholds[3]=3.; fPtThresholds[4]=4.;fPtThresholds[5]=5.;
  fPtFractions[0]=0.05; fPtFractions[1]=0.075; fPtFractions[2]=0.1; fPtFractions[3]=1.25; fPtFractions[4]=1.5;fPtFractions[5]=2.;

}

//__________________________________________________________________
void  AliAnaGammaDirect::MakeSeveralICAnalysis(AliAODParticleCorrelation* ph, Double_t v[3]) 
{
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  Float_t ptC = ph->Pt();
  Float_t phiC = ph->Phi();
  Float_t etaC = ph->Eta();

  fhPtGamma->Fill(ptC);
  fhPhiGamma->Fill(ptC,ph->Phi());
  fhEtaGamma->Fill(ptC,ph->Eta());
  Int_t tag =ph->GetTag();

  if(IsDataMC()){
    if(tag == kPrompt){
      fhPtPrompt   ->Fill(ptC);
      fhPhiPrompt ->Fill(ptC,phiC);
      fhEtaPrompt ->Fill(ptC,etaC);
    }
    else if(tag==kFragmentation)
      {
	fhPtFragmentation   ->Fill(ptC);
	fhPhiFragmentation ->Fill(ptC,phiC);
	fhEtaFragmentation ->Fill(ptC,etaC);
      }
    else if(tag==kPi0Decay)
      {
	fhPtPi0Decay   ->Fill(ptC);
	fhPhiPi0Decay ->Fill(ptC,phiC);
	fhEtaPi0Decay ->Fill(ptC,etaC);
      }
    else if(tag==kEtaDecay || tag==kOtherDecay)
      {
	fhPtOtherDecay   ->Fill(ptC);
	fhPhiOtherDecay ->Fill(ptC,phiC);
	fhEtaOtherDecay ->Fill(ptC,etaC);
      }
    else if(tag==kConversion)
      {
	fhPtConversion   ->Fill(ptC);
	fhPhiConversion ->Fill(ptC,phiC);
	fhEtaConversion ->Fill(ptC,etaC);
      }
    else{
      
      fhPtUnknown   ->Fill(ptC);
      fhPhiUnknown ->Fill(ptC,phiC);
      fhEtaUnknown ->Fill(ptC,etaC);
    }
  }//Histograms with MC
  //Keep original setting used when filling AODs, reset at end of analysis  
  Float_t ptthresorg = GetIsolationCut()->GetPtThreshold();
  Float_t ptfracorg = GetIsolationCut()->GetPtFraction();
  Float_t rorg = GetIsolationCut()->GetConeSize();

  Float_t coneptsum = 0 ;  
  Int_t n[10][10];//[fNCones][fNPtThresFrac];
  Int_t nfrac[10][10];//[fNCones][fNPtThresFrac];
  Bool_t  isolated   = kFALSE;

  for(Int_t icone = 0; icone<fNCones; icone++){
    GetIsolationCut()->SetConeSize(fConeSizes[icone]);
    coneptsum = 0 ;
    for(Int_t ipt = 0; ipt<fNPtThresFrac ;ipt++){
      n[icone][ipt]=0;
      nfrac[icone][ipt]=0;
      GetIsolationCut()->SetPtThreshold(fPtThresholds[ipt]);
      GetIsolationCut()->MakeIsolationCut((TSeqCollection*)ph->GetRefIsolationConeTracks(), 
					  (TSeqCollection*)ph->GetRefIsolationConeClusters(), 
					  v, kFALSE, ph,-1,-1,
					  n[icone][ipt],nfrac[icone][ipt],coneptsum, isolated);
      if(n[icone][ipt] == 0) {
	fhPtThresIsolated[icone][ipt]->Fill(ptC);
	if(IsDataMC()){
	  if(tag==kPrompt) fhPtThresIsolatedPrompt[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kConversion) fhPtThresIsolatedConversion[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kFragmentation) fhPtThresIsolatedFragmentation[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kPi0Decay) fhPtThresIsolatedPi0Decay[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kOtherDecay || tag==kEtaDecay) fhPtThresIsolatedOtherDecay[icone][ipt]->Fill(ptC,coneptsum) ;
	  else  fhPtThresIsolatedUnknown[icone][ipt]->Fill(ptC,coneptsum) ;
	}
      }
      if(nfrac[icone][ipt] == 0) {
	fhPtFracIsolated[icone][ipt]->Fill(ptC);
	if(IsDataMC()){
	  if(tag==kPrompt) fhPtFracIsolatedPrompt[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kConversion) fhPtFracIsolatedConversion[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kFragmentation) fhPtFracIsolatedFragmentation[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kPi0Decay) fhPtFracIsolatedPi0Decay[icone][ipt]->Fill(ptC,coneptsum) ;
	  else if(tag==kOtherDecay || tag==kEtaDecay) fhPtFracIsolatedOtherDecay[icone][ipt]->Fill(ptC,coneptsum) ;
	  else  fhPtFracIsolatedUnknown[icone][ipt]->Fill(ptC,coneptsum) ;
	}
      }
    }//pt thresh loop
    fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
    if(IsDataMC()){
      if(tag==kPrompt) fhPtSumIsolatedPrompt[icone]->Fill(ptC,coneptsum) ;
      else if(tag==kConversion) fhPtSumIsolatedConversion[icone]->Fill(ptC,coneptsum) ;
      else if(tag==kFragmentation) fhPtSumIsolatedFragmentation[icone]->Fill(ptC,coneptsum) ;
      else if(tag==kPi0Decay) fhPtSumIsolatedPi0Decay[icone]->Fill(ptC,coneptsum) ;
      else if(tag==kOtherDecay || tag==kEtaDecay) fhPtSumIsolatedOtherDecay[icone]->Fill(ptC,coneptsum) ;
      else  fhPtSumIsolatedUnknown[icone]->Fill(ptC,coneptsum) ;
    }

  }//cone size loop

  //Reset original parameters for AOD analysis
  GetIsolationCut()->SetPtThreshold(ptthresorg);
  GetIsolationCut()->SetPtFraction(ptfracorg);
  GetIsolationCut()->SetConeSize(rorg);

}

//__________________________________________________________________
void AliAnaGammaDirect::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  
  printf("Min Gamma pT       =     %2.2f\n",  GetMinPt()) ;
  printf("Max Gamma pT       =     %3.2f\n",  GetMaxPt()) ;

//   if(IsCaloPIDOn())printf("Check PID \n") ;
//   if(IsCaloPIDRecalculationOn())printf("Recalculate PID \n") ;
//   if(IsFidutialCutOn())printf("Check Fidutial cut \n") ;
  printf("Make Isolation     =     %d \n",  fMakeIC) ;
  printf("ReMake Isolation   =     %d \n",  fReMakeIC) ;
  printf("Make Several Isolation = %d \n",  fMakeSeveralIC) ;
  
  if(fMakeSeveralIC){
    printf("N Cone Sizes       =     %d\n", fNCones) ; 
    printf("Cone Sizes          =    \n") ;
    for(Int_t i = 0; i < fNCones; i++)
      printf("  %1.2f;",  fConeSizes[i]) ;
    printf("    \n") ;
    
    printf("N pT thresholds/fractions = %d\n", fNPtThresFrac) ;
    printf(" pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtThresholds[i]) ;
    
    printf("    \n") ;
    
    printf(" pT fractions         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtFractions[i]) ;
    
  }  
  
  printf("    \n") ;
  
} 

//____________________________________________________________________________
Bool_t  AliAnaGammaDirect::SelectCluster(AliAODCaloCluster * calo, Double_t vertex[3], TLorentzVector & mom){
   //Select cluster depending on its pid and acceptance selections
   
   //Skip matched clusters with tracks
   if(calo->GetNTracksMatched() > 0) return kFALSE ;
   
   //Check PID
   calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
   Int_t pdg = AliCaloPID::kPhoton;   
   if(IsCaloPIDOn()){
     //Get most probable PID, 2 options check PID weights (in MC this option is mandatory)
     //or redo PID, recommended option for EMCal.
     if(!IsCaloPIDRecalculationOn() || GetReader()->GetDataType() == AliCaloTrackReader::kMC )
       pdg = GetCaloPID()->GetPdg(fDetector,calo->PID(),mom.E());//PID with weights
     else
       pdg = GetCaloPID()->GetPdg(fDetector,mom,calo->GetM02(),0,0,0,0);//PID recalculated
     
     if(GetDebug() > 1) printf("PDG of identified particle %d\n",pdg);
     
     //If it does not pass pid, skip
     if(pdg != AliCaloPID::kPhoton) return kFALSE ;
   }
   
   //Check acceptance selection
   if(IsFidutialCutOn()){
     Bool_t in = GetFidutialCut()->IsInFidutialCut(mom,fDetector) ;
     if(! in ) return kFALSE ;
   }
   
   if(GetDebug() > 1) printf("Correlation photon selection cuts passed: pT %3.2f, pdg %d\n",mom.Pt(), pdg);
   
   return kTRUE;
   
 }
