/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: AliAnaParticleIsolation.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
// Class for analysis of particle isolation
// Input is selected particles put in AOD branch (AliAODPWG4ParticleCorrelation)
//
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TClonesArray.h>
#include <TList.h>
//#include <TObjString.h>
#include <TH2F.h>
//#include <Riostream.h>
#include <TClass.h>

// --- Analysis system --- 
#include "AliAnaParticleIsolation.h" 
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliMCAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliAODCaloCluster.h"

ClassImp(AliAnaParticleIsolation)
  
//____________________________________________________________________________
  AliAnaParticleIsolation::AliAnaParticleIsolation() : 
    AliAnaPartCorrBaseClass(), fCalorimeter(""), 
    fReMakeIC(0), fMakeSeveralIC(0), fMakeInvMass(0),
    fhPtIso(0),fhPhiIso(0),fhEtaIso(0), fhConeSumPt(0), fhPtInCone(0),
    //Several IC
    fNCones(0),fNPtThresFrac(0), fConeSizes(),  fPtThresholds(),  fPtFractions(), 
    //MC
    fhPtIsoPrompt(0),fhPhiIsoPrompt(0),fhEtaIsoPrompt(0), 
    fhPtThresIsolatedPrompt(), fhPtFracIsolatedPrompt(),  fhPtSumIsolatedPrompt(),
    fhPtIsoFragmentation(0),fhPhiIsoFragmentation(0),fhEtaIsoFragmentation(0), 
    fhPtThresIsolatedFragmentation(), fhPtFracIsolatedFragmentation(),  fhPtSumIsolatedFragmentation(),
    fhPtIsoPi0Decay(0),fhPhiIsoPi0Decay(0),fhEtaIsoPi0Decay(0),
    fhPtThresIsolatedPi0Decay(), fhPtFracIsolatedPi0Decay(),  fhPtSumIsolatedPi0Decay(),
    fhPtIsoOtherDecay(0),fhPhiIsoOtherDecay(0),fhEtaIsoOtherDecay(0), 
    fhPtThresIsolatedOtherDecay(), fhPtFracIsolatedOtherDecay(),  fhPtSumIsolatedOtherDecay(),
    fhPtIsoConversion(0),fhPhiIsoConversion(0),fhEtaIsoConversion(0), 
    fhPtThresIsolatedConversion(), fhPtFracIsolatedConversion(),  fhPtSumIsolatedConversion(),
    fhPtIsoUnknown(0),fhPhiIsoUnknown(0),fhEtaIsoUnknown(0), 
    fhPtThresIsolatedUnknown(), fhPtFracIsolatedUnknown(),  fhPtSumIsolatedUnknown(),
    //Histograms settings
    fHistoNPtSumBins(0),    fHistoPtSumMax(0.),    fHistoPtSumMin(0.),
    fHistoNPtInConeBins(0), fHistoPtInConeMax(0.), fHistoPtInConeMin(0.)
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
/*
//____________________________________________________________________________
AliAnaParticleIsolation::AliAnaParticleIsolation(const AliAnaParticleIsolation & g) : 
  AliAnaPartCorrBaseClass(g), fCalorimeter(g.fCalorimeter),
  fReMakeIC(g.fReMakeIC), fMakeSeveralIC(g.fMakeSeveralIC),  fMakeInvMass(g.fMakeInvMass),
  fhPtIso(g.fhPtIso),fhPhiIso(g.fhPhiIso),fhEtaIso(g.fhEtaIso), 
  fhConeSumPt(g.fhConeSumPt), fhPtInCone(g.fhPtInCone),
  //Several IC
  fNCones(g.fNCones),fNPtThresFrac(g.fNPtThresFrac), fConeSizes(), fPtThresholds(),  fPtFractions(), 
  fhPtThresIsolated(),  fhPtFracIsolated(), fhPtSumIsolated(),
  //MC
  fhPtIsoPrompt(g.fhPtIsoPrompt),fhPhiIsoPrompt(g.fhPhiIsoPrompt),fhEtaIsoPrompt(g.fhEtaIsoPrompt), 
  fhPtThresIsolatedPrompt(), fhPtFracIsolatedPrompt(),  fhPtSumIsolatedPrompt(),
  fhPtIsoFragmentation(g.fhPtIsoFragmentation),fhPhiIsoFragmentation(g.fhPhiIsoFragmentation),fhEtaIsoFragmentation(g.fhEtaIsoFragmentation),
  fhPtThresIsolatedFragmentation(), fhPtFracIsolatedFragmentation(),  fhPtSumIsolatedFragmentation(),
  fhPtIsoPi0Decay(g.fhPtIsoPi0Decay),fhPhiIsoPi0Decay(g.fhPhiIsoPi0Decay),fhEtaIsoPi0Decay(g.fhEtaIsoPi0Decay), 
  fhPtThresIsolatedPi0Decay(), fhPtFracIsolatedPi0Decay(),  fhPtSumIsolatedPi0Decay(),
  fhPtIsoOtherDecay(g.fhPtIsoOtherDecay),fhPhiIsoOtherDecay(g.fhPhiIsoOtherDecay),fhEtaIsoOtherDecay(g.fhEtaIsoOtherDecay), 
  fhPtThresIsolatedOtherDecay(), fhPtFracIsolatedOtherDecay(),  fhPtSumIsolatedOtherDecay(),
  fhPtIsoConversion(g. fhPtIsoConversion),fhPhiIsoConversion(g.fhPhiIsoConversion),fhEtaIsoConversion(g.fhEtaIsoConversion), 
  fhPtThresIsolatedConversion(), fhPtFracIsolatedConversion(),  fhPtSumIsolatedConversion(),
  fhPtIsoUnknown(g.fhPtIsoUnknown),fhPhiIsoUnknown(g.fhPhiIsoUnknown),fhEtaIsoUnknown(g.fhEtaIsoUnknown), 
  fhPtThresIsolatedUnknown(), fhPtFracIsolatedUnknown(),  fhPtSumIsolatedUnknown(),
  //Histograms
  fHistoNPtSumBins(g.fHistoNPtSumBins), fHistoPtSumMax(g.fHistoPtSumMax), fHistoPtSumMin(g.fHistoPtSumMax),
  fHistoNPtInConeBins(g.fHistoNPtInConeBins), fHistoPtInConeMax(g.fHistoPtInConeMax), fHistoPtInConeMin(g.fHistoPtInConeMin)
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
AliAnaParticleIsolation & AliAnaParticleIsolation::operator = (const AliAnaParticleIsolation & g)
{
  // assignment operator
  
  if(&g == this) return *this;
  
  fReMakeIC      = g.fReMakeIC ;
  fMakeSeveralIC = g.fMakeSeveralIC ;
  fMakeInvMass   = g.fMakeInvMass ;
  fCalorimeter   = g.fCalorimeter ;
	
  fhConeSumPt   = g.fhConeSumPt ;
  fhPtInCone    = g.fhPtInCone;
  
  fhPtIso  = g.fhPtIso ; 
  fhPhiIso = g.fhPhiIso ;
  fhEtaIso = g.fhEtaIso ;
  
  fhPtIsoPrompt  = g.fhPtIsoPrompt;
  fhPhiIsoPrompt = g.fhPhiIsoPrompt;
  fhEtaIsoPrompt = g.fhEtaIsoPrompt; 
  fhPtIsoFragmentation  = g.fhPtIsoFragmentation;
  fhPhiIsoFragmentation = g.fhPhiIsoFragmentation;
  fhEtaIsoFragmentation = g.fhEtaIsoFragmentation; 
  fhPtIsoPi0Decay  = g.fhPtIsoPi0Decay;
  fhPhiIsoPi0Decay = g.fhPhiIsoPi0Decay;
  fhEtaIsoPi0Decay = g.fhEtaIsoPi0Decay; 
  fhPtIsoOtherDecay  = g.fhPtIsoOtherDecay;
  fhPhiIsoOtherDecay = g.fhPhiIsoOtherDecay;
  fhEtaIsoOtherDecay = g.fhEtaIsoOtherDecay; 
  fhPtIsoConversion  = g. fhPtIsoConversion;
  fhPhiIsoConversion = g.fhPhiIsoConversion;
  fhEtaIsoConversion = g.fhEtaIsoConversion; 
  fhPtIsoUnknown  = g.fhPtIsoUnknown;
  fhPhiIsoUnknown = g.fhPhiIsoUnknown;
  fhEtaIsoUnknown = g.fhEtaIsoUnknown; 
  
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
  
  
  fHistoNPtSumBins    = g.fHistoNPtSumBins;
  fHistoPtSumMax      = g.fHistoPtSumMax; 
  fHistoPtSumMin      = g.fHistoPtSumMax;
  fHistoNPtInConeBins = g.fHistoNPtInConeBins;
  fHistoPtInConeMax   = g.fHistoPtInConeMax;
  fHistoPtInConeMin   = g.fHistoPtInConeMin;	
  
  return *this;
  
}
*/

//____________________________________________________________________________
AliAnaParticleIsolation::~AliAnaParticleIsolation() 
{
  //dtor
  //do not delete histograms
  
  //delete [] fConeSizes ; 
  //delete [] fPtThresholds ; 
  //delete [] fPtFractions ; 

}

//_________________________________________________________________________
Bool_t AliAnaParticleIsolation::CheckInvMass(const Int_t iaod, const AliAODPWG4Particle * part1) const
{
  // Search if there is a companion decay photon to the candidate 
  // and discard it in such case
  // Use it only if isolation candidates are photons
  // Make sure that no selection on photon pt is done in the input aod photon list.
  TLorentzVector mom1 = *(part1->Momentum());
  TLorentzVector mom2 ;
  for(Int_t jaod = 0; jaod < GetInputAODBranch()->GetEntriesFast(); jaod++){
    if(iaod == jaod) continue ;
    AliAODPWG4ParticleCorrelation * part2 =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(jaod));
    mom2 = *(part2->Momentum());
    //Select good pair (good phi, pt cuts, aperture and invariant mass)
    if(GetNeutralMesonSelection()->SelectPair(mom1, mom2)){
      if(GetDebug() > 1)printf("AliAnaParticleIsolation::CheckInvMass() - Selected gamma pair: pt %f, phi %f, eta%f\n",(mom1+mom2).Pt(), (mom1+mom2).Phi()*180./3.1416, (mom1+mom2).Eta());
      return kTRUE ;
    }
  }//loop
  
  return kFALSE;
}

//________________________________________________________________________
TList *  AliAnaParticleIsolation::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("IsolatedParticleHistos") ; 
  
  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	
  
  Int_t nptsumbins    = fHistoNPtSumBins;
  Float_t ptsummax    = fHistoPtSumMax;
  Float_t ptsummin    = fHistoPtSumMin;	
  Int_t nptinconebins = fHistoNPtInConeBins;
  Float_t ptinconemax = fHistoPtInConeMax;
  Float_t ptinconemin = fHistoPtInConeMin;
  
  if(!fMakeSeveralIC){
    
    fhConeSumPt  = new TH2F
      ("hConePtSum","#Sigma p_{T} in isolation cone ",nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPt->SetYTitle("#Sigma p_{T}");
    fhConeSumPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhConeSumPt) ;
    
    fhPtInCone  = new TH2F
      ("hPtInCone","p_{T} in isolation cone ",nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInCone->SetYTitle("p_{T in cone} (GeV/c)");
    fhPtInCone->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtInCone) ;
    
    fhPtIso  = new TH1F("hPt","Isolated Number of particles",nptbins,ptmin,ptmax); 
    fhPtIso->SetYTitle("N");
    fhPtIso->SetXTitle("p_{T}(GeV/c)");
    outputContainer->Add(fhPtIso) ; 
    
    fhPhiIso  = new TH2F
      ("hPhi","Isolated Number of particles",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiIso->SetYTitle("#phi");
    fhPhiIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPhiIso) ; 
    
    fhEtaIso  = new TH2F
      ("hEta","Isolated Number of particles",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaIso->SetYTitle("#eta");
    fhEtaIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhEtaIso) ;
    
    if(IsDataMC()){
      
      fhPtIsoPrompt  = new TH1F("hPtMCPrompt","Isolated Number of #gamma prompt",nptbins,ptmin,ptmax); 
      fhPtIsoPrompt->SetYTitle("N");
      fhPtIsoPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPrompt) ; 
      
      fhPhiIsoPrompt  = new TH2F
	("hPhiMCPrompt","Isolated Number of #gamma prompt",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPrompt->SetYTitle("#phi");
      fhPhiIsoPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPrompt) ; 
      
      fhEtaIsoPrompt  = new TH2F
	("hEtaMCPrompt","Isolated Number of #gamma prompt ",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPrompt->SetYTitle("#eta");
      fhEtaIsoPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPrompt) ;
      
      fhPtIsoFragmentation  = new TH1F("hPtMCFragmentation","Isolated Number of #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoFragmentation->SetYTitle("N");
      fhPtIsoFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoFragmentation) ; 
      
      fhPhiIsoFragmentation  = new TH2F
	("hPhiMCFragmentation","Isolated Number of #gamma fragmentation",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoFragmentation->SetYTitle("#phi");
      fhPhiIsoFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoFragmentation) ; 
      
      fhEtaIsoFragmentation  = new TH2F
	("hEtaMCFragmentation","Isolated Number of #gamma fragmentation",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoFragmentation->SetYTitle("#eta");
      fhEtaIsoFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoFragmentation) ;
      
      fhPtIsoPi0Decay  = new TH1F("hPtMCPi0Decay","Isolated Number of #gamma from #pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtIsoPi0Decay->SetYTitle("N");
      fhPtIsoPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPi0Decay) ; 
      
      fhPhiIsoPi0Decay  = new TH2F
	("hPhiMCPi0Decay","Isolated Number of #gamma from #pi^{0} decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPi0Decay->SetYTitle("#phi");
      fhPhiIsoPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPi0Decay) ; 
      
      fhEtaIsoPi0Decay  = new TH2F
	("hEtaMCPi0Decay","Isolated Number of #gamma from #pi^{0} decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPi0Decay->SetYTitle("#eta");
      fhEtaIsoPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPi0Decay) ;
      
      fhPtIsoOtherDecay  = new TH1F("hPtMCOtherDecay","Isolated Number of #gamma from non #pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtIsoOtherDecay->SetYTitle("N");
      fhPtIsoOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoOtherDecay) ; 
      
      fhPhiIsoOtherDecay  = new TH2F
	("hPhiMCOtherDecay","Isolated Number of #gamma from non #pi^{0} decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoOtherDecay->SetYTitle("#phi");
      fhPhiIsoOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoOtherDecay) ; 
      
      fhEtaIsoOtherDecay  = new TH2F
	("hEtaMCOtherDecay","Isolated Number of #gamma non #pi^{0} decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoOtherDecay->SetYTitle("#eta");
      fhEtaIsoOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoOtherDecay) ;
      
      fhPtIsoConversion  = new TH1F("hPtMCConversion","Isolated Number of #gamma converted",nptbins,ptmin,ptmax); 
      fhPtIsoConversion->SetYTitle("N");
      fhPtIsoConversion->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoConversion) ; 
      
      fhPhiIsoConversion  = new TH2F
	("hPhiMCConversion","Isolated Number of #gamma converted",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoConversion->SetYTitle("#phi");
      fhPhiIsoConversion->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoConversion) ; 
      
      fhEtaIsoConversion  = new TH2F
	("hEtaMCConversion","Isolated Number of #gamma converted",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoConversion->SetYTitle("#eta");
      fhEtaIsoConversion->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoConversion) ;
      
      fhPtIsoUnknown  = new TH1F("hPtMCUnknown","Isolated Number of non #gamma particles",nptbins,ptmin,ptmax); 
      fhPtIsoUnknown->SetYTitle("N");
      fhPtIsoUnknown->SetXTitle("p_{T}(GeV/c)");
      outputContainer->Add(fhPtIsoUnknown) ; 
      
      fhPhiIsoUnknown  = new TH2F
	("hPhiMCUnknown","Isolated Number of non #gamma particles",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoUnknown->SetYTitle("#phi");
      fhPhiIsoUnknown->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPhiIsoUnknown) ; 
      
      fhEtaIsoUnknown  = new TH2F
	("hEtaMCUnknown","Isolated Number of non #gamma particles",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoUnknown->SetYTitle("#eta");
      fhEtaIsoUnknown->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhEtaIsoUnknown) ;
    }//Histos with MC
    
  }
  
  if(fMakeSeveralIC){
		char name[128];
		char title[128];
		for(Int_t icone = 0; icone<fNCones; icone++){
		  sprintf(name,"hPtSum_Cone_%d",icone);
		  sprintf(title,"Candidate cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		  fhPtSumIsolated[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		  fhPtSumIsolated[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		  fhPtSumIsolated[icone]->SetXTitle("p_{T} (GeV/c)");
		  outputContainer->Add(fhPtSumIsolated[icone]) ; 
		  
		  if(IsDataMC()){
		    sprintf(name,"hPtSumPrompt_Cone_%d",icone);
		    sprintf(title,"Candidate Prompt cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedPrompt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedPrompt[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedPrompt[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedPrompt[icone]) ; 
		    
		    sprintf(name,"hPtSumFragmentation_Cone_%d",icone);
		    sprintf(title,"Candidate Fragmentation cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedFragmentation[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedFragmentation[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedFragmentation[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedFragmentation[icone]) ; 
		    
		    sprintf(name,"hPtSumPi0Decay_Cone_%d",icone);
		    sprintf(title,"Candidate Pi0Decay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedPi0Decay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedPi0Decay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedPi0Decay[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedPi0Decay[icone]) ; 
		    
		    sprintf(name,"hPtSumOtherDecay_Cone_%d",icone);
		    sprintf(title,"Candidate OtherDecay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedOtherDecay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedOtherDecay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedOtherDecay[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedOtherDecay[icone]) ; 
		    
		    sprintf(name,"hPtSumConversion_Cone_%d",icone);
		    sprintf(title,"Candidate Conversion cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedConversion[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedConversion[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedConversion[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedConversion[icone]) ; 
		    
		    sprintf(name,"hPtSumUnknown_Cone_%d",icone);
		    sprintf(title,"Candidate Unknown cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedUnknown[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedUnknown[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedUnknown[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedUnknown[icone]) ; 
		    
		  }//Histos with MC
		  
		  for(Int_t ipt = 0; ipt<fNPtThresFrac;ipt++){ 
		    sprintf(name,"hPtThres_Cone_%d_Pt%d",icone,ipt);
		    sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		    fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		    fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 
		    
		    sprintf(name,"hPtFrac_Cone_%d_Pt%d",icone,ipt);
		    sprintf(title,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		    fhPtFracIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		    fhPtFracIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtFracIsolated[icone][ipt]) ; 
		    
		    if(IsDataMC()){
		      sprintf(name,"hPtThresMCPrompt_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedPrompt[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedPrompt[icone][ipt]) ; 
		      
		      sprintf(name,"hPtFracMCPrompt_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedPrompt[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedPrompt[icone][ipt]) ; 
		      
		      sprintf(name,"hPtThresMCFragmentation_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedFragmentation[icone][ipt]) ; 
		      
		      sprintf(name,"hPtFracMCFragmentation_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedFragmentation[icone][ipt]) ; 
		      
		      sprintf(name,"hPtThresMCPi0Decay_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedPi0Decay[icone][ipt]) ; 
		      
		      sprintf(name,"hPtFracMCPi0Decay_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedPi0Decay[icone][ipt]) ; 
		      
		      sprintf(name,"hPtThresMCOtherDecay_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedOtherDecay[icone][ipt]) ; 
		      
		      sprintf(name,"hPtFracMCOtherDecay_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedOtherDecay[icone][ipt]) ;
		      
		      sprintf(name,"hPtThresMCConversion_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedConversion[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedConversion[icone][ipt]) ; 
		      
		      sprintf(name,"hPtFracMCConversion_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedConversion[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedConversion[icone][ipt]) ;
		      
		      sprintf(name,"hPtThresMCUnknown_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Unknown p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedUnknown[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedUnknown[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedUnknown[icone][ipt]) ; 
		      
		      sprintf(name,"hPtFracMCUnknown_Cone_%d_Pt%d",icone,ipt);
		      sprintf(title,"Isolated candidate Unknown p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedUnknown[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
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
    if(GetNeutralMesonSelection()->AreNeutralMesonSelectionHistosKept())
      for(Int_t i = 0; i < nmsHistos->GetEntries(); i++) outputContainer->Add(nmsHistos->At(i)) ;
	delete nmsHistos;
  }
  
  //Save parameters used for analysis
//  TString parList ; //this will be list of parameters used for this analysis.
//  char onePar[255] ;
//  
//  sprintf(onePar,"--- AliAnaParticleIsolation ---\n") ;
//  parList+=onePar ;	
//  sprintf(onePar,"Calorimeter: %s\n",fCalorimeter.Data()) ;
//  parList+=onePar ;
//  sprintf(onePar,"fReMakeIC =%d (Flag for reisolation during histogram filling) \n",fReMakeIC) ;
//  parList+=onePar ;
//  sprintf(onePar,"fMakeSeveralIC=%d (Flag for isolation with several cuts at the same time ) \n",fMakeSeveralIC) ;
//  parList+=onePar ;
//  sprintf(onePar,"fMakeInvMass=%d (Flag for rejection of candidates with a pi0 inv mass pair) \n",fMakeInvMass) ;
//  parList+=onePar ;
//  
//  if(fMakeSeveralIC){
//    sprintf(onePar,"fNCones =%d (Number of cone sizes) \n",fNCones) ;
//    parList+=onePar ;
//    sprintf(onePar,"fNPtThresFrac=%d (Flag for isolation with several cuts at the same time ) \n",fNPtThresFrac) ;
//    parList+=onePar ;
//    
//    for(Int_t icone = 0; icone < fNCones ; icone++){
//      sprintf(onePar,"fConeSizes[%d]=%1.2f (isolation cone size) \n",icone, fConeSizes[icone]) ;
//      parList+=onePar ;	
//    }
//    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++){
//      sprintf(onePar,"fPtThresholds[%d]=%1.2f (isolation pt threshold) \n",ipt, fPtThresholds[ipt]) ;
//      parList+=onePar ;	
//    }
//    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++){
//      sprintf(onePar,"fPtFractions[%d]=%1.2f (isolation pt fraction threshold) \n",ipt, fPtFractions[ipt]) ;
//      parList+=onePar ;	
//    }		
//  }
//  
//  //Get parameters set in base class.
//  parList += GetBaseParametersList() ;
//  
//  //Get parameters set in IC class.
//  if(!fMakeSeveralIC)parList += GetIsolationCut()->GetICParametersList() ;
//  
//  TObjString *oString= new TObjString(parList) ;
//  outputContainer->Add(oString);
  
  
  return outputContainer ;
  
}

//__________________________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  //Search for the isolated photon in fCalorimeter with pt > GetMinPt()
  
  if(!GetInputAODBranch()){
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - No input particles in AOD with name branch < %s >, STOP \n",GetInputAODName().Data());
    abort();
  }
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation")){
	printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s> \n",GetInputAODBranch()->GetClass()->GetName());
	abort();
  }
	
  Int_t n = 0, nfrac = 0;
  Bool_t isolated = kFALSE ; 
  Float_t coneptsum = 0 ;
  TObjArray * pl = 0x0; ; 
  
  //Select the calorimeter for candidate isolation with neutral particles
  if(fCalorimeter == "PHOS")
    pl = GetAODPHOS();
  else if (fCalorimeter == "EMCAL")
    pl = GetAODEMCAL();
  
  //Loop on AOD branch, filled previously in AliAnaPhoton
  TLorentzVector mom ;
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Input aod branch entries %d\n", naod);
  for(Int_t iaod = 0; iaod < naod; iaod++){
    AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    //If too small or too large pt, skip
    if(aodinput->Pt() < GetMinPt() || aodinput->Pt() > GetMaxPt() ) continue ; 
    
    //Check invariant mass, if pi0, skip.
    Bool_t decay = kFALSE ;
    if(fMakeInvMass) decay = CheckInvMass(iaod,aodinput);
    if(decay) continue ;

    //After cuts, study isolation
    n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
    GetIsolationCut()->MakeIsolationCut(GetAODCTS(),pl,GetReader(), kTRUE, aodinput, GetAODObjArrayName(), n,nfrac,coneptsum, isolated);
    aodinput->SetIsolated(isolated);
    if(GetDebug() > 1 && isolated) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() : Particle %d IS ISOLATED \n",iaod);  
	  
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - End fill AODs \n");  
  
}

//__________________________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  Int_t n = 0, nfrac = 0;
  Bool_t isolated = kFALSE ; 
  Float_t coneptsum = 0 ;
  //Loop on stored AOD 
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Histo aod branch entries %d\n", naod);

  //Get vertex for photon momentum calculation
  Double_t vertex[]={0,0,0} ; //vertex ;
  Double_t vertex2[]={0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
	  GetReader()->GetVertex(vertex);
	  if(GetReader()->GetSecondInputAODTree()) GetReader()->GetSecondInputAODVertex(vertex2);
  }	
	
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* aod =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));

    Bool_t isolation   = aod->IsIsolated();              
    Float_t ptcluster  = aod->Pt();
    Float_t phicluster = aod->Phi();
    Float_t etacluster = aod->Eta();
    //Recover reference arrays with clusters and tracks
	TObjArray * refclusters = aod->GetObjArray(GetAODObjArrayName()+"Clusters");
    TObjArray * reftracks = aod->GetObjArray(GetAODObjArrayName()+"Tracks");
  
    if(fMakeSeveralIC) {
      //Analysis of multiple IC at same time
      MakeSeveralICAnalysis(aod);
      continue;
    }
    else if(fReMakeIC){
      //In case a more strict IC is needed in the produced AOD
      n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
      GetIsolationCut()->MakeIsolationCut(reftracks, refclusters, GetReader(), kFALSE, aod, "", n,nfrac,coneptsum, isolated);
      fhConeSumPt->Fill(ptcluster,coneptsum);    
	  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Energy Sum in Isolation Cone %2.2f\n", coneptsum);    
	}
    
    //Fill pt distribution of particles in cone
    //Tracks
    coneptsum=0;
	if(reftracks){  
        for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++){
            AliAODTrack* track = (AliAODTrack *) reftracks->At(itrack);
            fhPtInCone->Fill(ptcluster,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
            coneptsum+=track->Pt();
        }
	}
	  
    //CaloClusters
	if(refclusters){    
		TLorentzVector mom ;
        for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++){
            AliAODCaloCluster* calo = (AliAODCaloCluster *) refclusters->At(icalo);
			Int_t input = 0;
			if     (fCalorimeter == "EMCAL" && GetReader()->GetAODEMCALNormalInputEntries() <= icalo) input = 1 ;
			else if(fCalorimeter == "PHOS"  && GetReader()->GetAODPHOSNormalInputEntries()  <= icalo) input = 1;
			
			//Get Momentum vector, 
			if     (input == 0) calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
			else if(input == 1) calo->GetMomentum(mom,vertex2);//Assume that come from vertex in straight line  
			
            fhPtInCone->Fill(ptcluster, mom.Pt());
            coneptsum+=mom.Pt();
       }
    }
	  
	if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d Energy Sum in Isolation Cone %2.2f\n", iaod, coneptsum);
    
	if(!fReMakeIC) fhConeSumPt->Fill(ptcluster,coneptsum);
    
    if(isolation){    
		
	  if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d ISOLATED, fill histograms\n", iaod);
		
      fhPtIso  ->Fill(ptcluster);
      fhPhiIso ->Fill(ptcluster,phicluster);
      fhEtaIso ->Fill(ptcluster,etacluster);
      
      if(IsDataMC()){
	Int_t tag =aod->GetTag();
	
	if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)){
	  fhPtIsoPrompt  ->Fill(ptcluster);
	  fhPhiIsoPrompt ->Fill(ptcluster,phicluster);
	  fhEtaIsoPrompt ->Fill(ptcluster,etacluster);
	}
	else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
	  {
	    fhPtIsoFragmentation  ->Fill(ptcluster);
	    fhPhiIsoFragmentation ->Fill(ptcluster,phicluster);
	    fhEtaIsoFragmentation ->Fill(ptcluster,etacluster);
	  }
	else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
	  {
	    fhPtIsoPi0Decay  ->Fill(ptcluster);
	    fhPhiIsoPi0Decay ->Fill(ptcluster,phicluster);
	    fhEtaIsoPi0Decay ->Fill(ptcluster,etacluster);
	  }
	else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
	  {
	    fhPtIsoOtherDecay  ->Fill(ptcluster);
	    fhPhiIsoOtherDecay ->Fill(ptcluster,phicluster);
	    fhEtaIsoOtherDecay ->Fill(ptcluster,etacluster);
	  }
	else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
	  {
	    fhPtIsoConversion  ->Fill(ptcluster);
	    fhPhiIsoConversion ->Fill(ptcluster,phicluster);
	    fhEtaIsoConversion ->Fill(ptcluster,etacluster);
	  }
	else
	  {
	    fhPtIsoUnknown  ->Fill(ptcluster);
	    fhPhiIsoUnknown ->Fill(ptcluster,phicluster);
	    fhEtaIsoUnknown ->Fill(ptcluster,etacluster);
	  }
      }//Histograms with MC
      
    }//Isolated histograms
    
  }// aod loop
  
}

//____________________________________________________________________________
void AliAnaParticleIsolation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("IsolationCone");  
  AddToHistogramsName("AnaIsolation_");

  fCalorimeter = "PHOS" ;
  fReMakeIC = kFALSE ;
  fMakeSeveralIC = kFALSE ;
  fMakeInvMass = kFALSE ;
  
 //----------- Several IC-----------------
  fNCones           = 5 ; 
  fNPtThresFrac     = 5 ; 
  fConeSizes[0]    = 0.1;  fConeSizes[1]     = 0.2;  fConeSizes[2]    = 0.3; fConeSizes[3]    = 0.4;  fConeSizes[4]    = 0.5;
  fPtThresholds[0] = 1.;   fPtThresholds[1] = 2.;    fPtThresholds[2] = 3.;  fPtThresholds[3] = 4.;   fPtThresholds[4] = 5.; 
  fPtFractions[0]  = 0.05; fPtFractions[1]  = 0.075; fPtFractions[2]  = 0.1; fPtFractions[3]  = 1.25; fPtFractions[4]  = 1.5; 

//------------- Histograms settings -------
  fHistoNPtSumBins = 100 ;
  fHistoPtSumMax   = 50 ;
  fHistoPtSumMin   = 0.  ;

  fHistoNPtInConeBins = 100 ;
  fHistoPtInConeMax   = 50 ;
  fHistoPtInConeMin   = 0.  ;

}

//__________________________________________________________________
void  AliAnaParticleIsolation::MakeSeveralICAnalysis(AliAODPWG4ParticleCorrelation* ph) 
{
  //Isolation Cut Analysis for both methods and different pt cuts and cones
  Float_t ptC = ph->Pt();	
  Int_t tag   = ph->GetTag();
  
  //Keep original setting used when filling AODs, reset at end of analysis  
  Float_t ptthresorg = GetIsolationCut()->GetPtThreshold();
  Float_t ptfracorg  = GetIsolationCut()->GetPtFraction();
  Float_t rorg       = GetIsolationCut()->GetConeSize();
  
  Float_t coneptsum = 0 ;  
  Int_t n[10][10];//[fNCones][fNPtThresFrac];
  Int_t nfrac[10][10];//[fNCones][fNPtThresFrac];
  Bool_t  isolated   = kFALSE;
  
  //Loop on cone sizes
  for(Int_t icone = 0; icone<fNCones; icone++){
    GetIsolationCut()->SetConeSize(fConeSizes[icone]);
    coneptsum = 0 ;

    //Loop on ptthresholds
    for(Int_t ipt = 0; ipt<fNPtThresFrac ;ipt++){
      n[icone][ipt]=0;
      nfrac[icone][ipt]=0;
      GetIsolationCut()->SetPtThreshold(fPtThresholds[ipt]);
      GetIsolationCut()->MakeIsolationCut(ph->GetObjArray(GetAODObjArrayName()+"Tracks"), 
				  ph->GetObjArray(GetAODObjArrayName()+"Clusters"),
				  GetReader(), kFALSE, ph, "",n[icone][ipt],nfrac[icone][ipt],coneptsum, isolated);

      //Normal ptThreshold cut
      if(n[icone][ipt] == 0) {
	fhPtThresIsolated[icone][ipt]->Fill(ptC);
	if(IsDataMC()){
	  if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)) fhPtThresIsolatedPrompt[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) fhPtThresIsolatedConversion[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtThresIsolatedFragmentation[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)) fhPtThresIsolatedPi0Decay[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay)) fhPtThresIsolatedOtherDecay[icone][ipt]->Fill(ptC) ;
	  else  fhPtThresIsolatedUnknown[icone][ipt]->Fill(ptC) ;
	}
      }
      
      //Pt threshold on pt cand/ pt in cone fraction
      if(nfrac[icone][ipt] == 0) {
	fhPtFracIsolated[icone][ipt]->Fill(ptC);
	if(IsDataMC()){
	  if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)) fhPtFracIsolatedPrompt[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) fhPtFracIsolatedConversion[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtFracIsolatedFragmentation[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)) fhPtFracIsolatedPi0Decay[icone][ipt]->Fill(ptC) ;
	  else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay)) fhPtFracIsolatedOtherDecay[icone][ipt]->Fill(ptC) ;
	  else  fhPtFracIsolatedUnknown[icone][ipt]->Fill(ptC) ;
	}
      }
    }//pt thresh loop
    
    //Sum in cone histograms
    fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
    if(IsDataMC()){
      if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)) fhPtSumIsolatedPrompt[icone]->Fill(ptC,coneptsum) ;
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) fhPtSumIsolatedConversion[icone]->Fill(ptC,coneptsum) ;
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtSumIsolatedFragmentation[icone]->Fill(ptC,coneptsum) ;
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)) fhPtSumIsolatedPi0Decay[icone]->Fill(ptC,coneptsum) ;
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay)) fhPtSumIsolatedOtherDecay[icone]->Fill(ptC,coneptsum) ;
      else  fhPtSumIsolatedUnknown[icone]->Fill(ptC,coneptsum) ;
    }
    
  }//cone size loop
  
  //Reset original parameters for AOD analysis
  GetIsolationCut()->SetPtThreshold(ptthresorg);
  GetIsolationCut()->SetPtFraction(ptfracorg);
  GetIsolationCut()->SetConeSize(rorg);
  
}

//__________________________________________________________________
void AliAnaParticleIsolation::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  
  printf("ReMake Isolation          = %d \n",  fReMakeIC) ;
  printf("Make Several Isolation    = %d \n",  fMakeSeveralIC) ;
  printf("Calorimeter for isolation = %s \n",  fCalorimeter.Data()) ;
  
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
  
  printf("Histograms: %3.1f < pT sum < %3.1f,  Nbin = %d\n", fHistoPtSumMin,  fHistoPtSumMax,  fHistoNPtSumBins);
  printf("Histograms: %3.1f < pT in cone < %3.1f, Nbin = %d\n", fHistoPtInConeMin, fHistoPtInConeMax, fHistoNPtInConeBins);
  
  printf("    \n") ;
  
} 

