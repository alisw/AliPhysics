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

//-Yaxian Mao (add the possibility for different IC method with different pt range, 01/10/2010)
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TClonesArray.h>
#include <TList.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TClass.h>

// --- Analysis system --- 
#include "AliAnaParticleIsolation.h" 
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliMCAnalysisUtils.h"
#include "AliVTrack.h"
#include "AliVCluster.h"

ClassImp(AliAnaParticleIsolation)
  
//______________________________________________________________________________
  AliAnaParticleIsolation::AliAnaParticleIsolation() : 
    AliAnaPartCorrBaseClass(),        fCalorimeter(""), 
    fReMakeIC(0),                     fMakeSeveralIC(0),
    // Several IC
    fNCones(0),                       fNPtThresFrac(0), 
    fConeSizes(),                     fPtThresholds(),                 fPtFractions(), 
    // Histograms
    fhPtIso(0),                       fhPhiIso(0),                     fhEtaIso(0), 
    fhPtNoIso(0),                     fhPtDecayIso(0),                 fhPtDecayNoIso(0), 
    fhConeSumPt(0),                   fhPtInCone(0),
    fhFRConeSumPt(0),                 fhPtInFRCone(0),
    // MC histograms
    fhPtIsoPrompt(0),                 fhPhiIsoPrompt(0),               fhEtaIsoPrompt(0), 
    fhPtThresIsolatedPrompt(),        fhPtFracIsolatedPrompt(),        fhPtSumIsolatedPrompt(),
    fhPtIsoFragmentation(0),          fhPhiIsoFragmentation(0),        fhEtaIsoFragmentation(0), 
    fhPtThresIsolatedFragmentation(), fhPtFracIsolatedFragmentation(), fhPtSumIsolatedFragmentation(),
    fhPtIsoPi0Decay(0),               fhPhiIsoPi0Decay(0),             fhEtaIsoPi0Decay(0),
    fhPtThresIsolatedPi0Decay(),      fhPtFracIsolatedPi0Decay(),      fhPtSumIsolatedPi0Decay(),
    fhPtIsoEtaDecay(0),               fhPhiIsoEtaDecay(0),             fhEtaIsoEtaDecay(0),
    fhPtThresIsolatedEtaDecay(),      fhPtFracIsolatedEtaDecay(),      fhPtSumIsolatedEtaDecay(),
    fhPtIsoOtherDecay(0),             fhPhiIsoOtherDecay(0),           fhEtaIsoOtherDecay(0), 
    fhPtThresIsolatedOtherDecay(),    fhPtFracIsolatedOtherDecay(),    fhPtSumIsolatedOtherDecay(),
    fhPtIsoConversion(0),             fhPhiIsoConversion(0),           fhEtaIsoConversion(0), 
    fhPtThresIsolatedConversion(),    fhPtFracIsolatedConversion(),    fhPtSumIsolatedConversion(),
    fhPtIsoUnknown(0),                fhPhiIsoUnknown(0),              fhEtaIsoUnknown(0), 
    fhPtThresIsolatedUnknown(),       fhPtFracIsolatedUnknown(),       fhPtSumIsolatedUnknown(),
    fhPtNoIsoPi0Decay(0),             fhPtNoIsoEtaDecay(0),            fhPtNoIsoOtherDecay(0),
    fhPtNoIsoPrompt(0),               fhPtIsoMCPhoton(0),              fhPtNoIsoMCPhoton(0),
    //Histograms settings
    fHistoNPtSumBins(0),              fHistoPtSumMax(0.),              fHistoPtSumMin(0.),
    fHistoNPtInConeBins(0),           fHistoPtInConeMax(0.),           fHistoPtInConeMin(0.)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  	
  for(Int_t i = 0; i < 5 ; i++){ 
    fConeSizes[i]      = 0 ; 
    fhPtSumIsolated[i] = 0 ;  
    
    fhPtSumIsolatedPrompt[i]        = 0 ;  
    fhPtSumIsolatedFragmentation[i] = 0 ;  
    fhPtSumIsolatedPi0Decay[i]      = 0 ;  
    fhPtSumIsolatedEtaDecay[i]      = 0 ;  
    fhPtSumIsolatedOtherDecay[i]    = 0 ;  
    fhPtSumIsolatedConversion[i]    = 0 ;  
    fhPtSumIsolatedUnknown[i]       = 0 ;  
    
    for(Int_t j = 0; j < 5 ; j++){ 
      fhPtThresIsolated[i][j] = 0 ;  
      fhPtFracIsolated[i][j]  = 0 ; 
      
      fhPtThresIsolatedPrompt[i][j]       = 0 ;  
      fhPtThresIsolatedFragmentation[i][j]= 0 ; 
      fhPtThresIsolatedPi0Decay[i][j]     = 0 ; 
      fhPtThresIsolatedEtaDecay[i][j]     = 0 ;  
      fhPtThresIsolatedOtherDecay[i][j]   = 0 ;  
      fhPtThresIsolatedConversion[i][j]   = 0 ;  
      fhPtThresIsolatedUnknown[i][j]      = 0 ;  
  
      fhPtFracIsolatedPrompt[i][j]        = 0 ;  
      fhPtFracIsolatedFragmentation[i][j] = 0 ;  
      fhPtFracIsolatedPi0Decay[i][j]      = 0 ;  
      fhPtFracIsolatedEtaDecay[i][j]      = 0 ;  
      fhPtFracIsolatedOtherDecay[i][j]    = 0 ;  
      fhPtFracIsolatedConversion[i][j]    = 0 ;
      fhPtFracIsolatedUnknown[i][j]       = 0 ;  
 
    }  
  } 
  
  for(Int_t i = 0; i < 5 ; i++){ 
    fPtFractions[i] = 0 ; 
    fPtThresholds[i]= 0 ; 
  } 

}

//______________________________________________________
TObjString *  AliAnaParticleIsolation::GetAnalysisCuts()
{ 
	//Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar, buffersize,"--- AliAnaParticleIsolation ---\n") ;
  parList+=onePar ;	
  snprintf(onePar, buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fReMakeIC =%d (Flag for reisolation during histogram filling) \n",fReMakeIC) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fMakeSeveralIC=%d (Flag for isolation with several cuts at the same time ) \n",fMakeSeveralIC) ;
  parList+=onePar ;
  
  if(fMakeSeveralIC){
    snprintf(onePar, buffersize,"fNCones =%d (Number of cone sizes) \n",fNCones) ;
    parList+=onePar ;
    snprintf(onePar, buffersize,"fNPtThresFrac=%d (Flag for isolation with several cuts at the same time ) \n",fNPtThresFrac) ;
    parList+=onePar ;
    
    for(Int_t icone = 0; icone < fNCones ; icone++){
      snprintf(onePar, buffersize,"fConeSizes[%d]=%1.2f (isolation cone size) \n",icone, fConeSizes[icone]) ;
      parList+=onePar ;	
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++){
      snprintf(onePar, buffersize,"fPtThresholds[%d]=%1.2f (isolation pt threshold) \n",ipt, fPtThresholds[ipt]) ;
      parList+=onePar ;	
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++){
      snprintf(onePar, buffersize,"fPtFractions[%d]=%1.2f (isolation pt fraction threshold) \n",ipt, fPtFractions[ipt]) ;
      parList+=onePar ;	
    }		
  }
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in IC class.
  if(!fMakeSeveralIC)parList += GetIsolationCut()->GetICParametersList() ;
  
  return new TObjString(parList) ;
	
}

//________________________________________________________
TList *  AliAnaParticleIsolation::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("IsolatedParticleHistos") ; 
  
  Int_t   nptbins  = GetHistoPtBins();
  Int_t   nphibins = GetHistoPhiBins();
  Int_t   netabins = GetHistoEtaBins();
  Float_t ptmax    = GetHistoPtMax();
  Float_t phimax   = GetHistoPhiMax();
  Float_t etamax   = GetHistoEtaMax();
  Float_t ptmin    = GetHistoPtMin();
  Float_t phimin   = GetHistoPhiMin();
  Float_t etamin   = GetHistoEtaMin();	
  
  Int_t   nptsumbins    = fHistoNPtSumBins;
  Float_t ptsummax      = fHistoPtSumMax;
  Float_t ptsummin      = fHistoPtSumMin;	
  Int_t   nptinconebins = fHistoNPtInConeBins;
  Float_t ptinconemax   = fHistoPtInConeMax;
  Float_t ptinconemin   = fHistoPtInConeMin;
  
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
    
    fhFRConeSumPt  = new TH2F
    ("hFRConePtSum","#Sigma p_{T} in the froward region isolation cone ",nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhFRConeSumPt->SetYTitle("#Sigma p_{T}");
    fhFRConeSumPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhFRConeSumPt) ;
    
    fhPtInFRCone  = new TH2F
    ("hPtInFRCone","p_{T} in froward region isolation cone ",nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInFRCone->SetYTitle("p_{T in cone} (GeV/c)");
    fhPtInFRCone->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtInFRCone) ;    
    
    fhPtIso  = new TH1F("hPt","Number of isolated particles",nptbins,ptmin,ptmax); 
    fhPtIso->SetYTitle("N");
    fhPtIso->SetXTitle("p_{T}(GeV/c)");
    outputContainer->Add(fhPtIso) ; 
    
    fhPhiIso  = new TH2F
    ("hPhi","Number of isolated particles",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiIso->SetYTitle("#phi");
    fhPhiIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPhiIso) ; 
    
    fhEtaIso  = new TH2F
    ("hEta","Number of isolated particles",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaIso->SetYTitle("#eta");
    fhEtaIso->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhEtaIso) ;
    
    fhPtNoIso  = new TH1F("hPtNoIso","Number of not isolated leading particles",nptbins,ptmin,ptmax); 
    fhPtNoIso->SetYTitle("N");
    fhPtNoIso->SetXTitle("p_{T}(GeV/c)");
    outputContainer->Add(fhPtNoIso) ;
    
    fhPtDecayIso  = new TH1F("hPtDecayIso","Number of isolated #pi^{0} decay particles",nptbins,ptmin,ptmax); 
    fhPtNoIso->SetYTitle("N");
    fhPtNoIso->SetXTitle("p_{T}(GeV/c)");
    outputContainer->Add(fhPtDecayIso) ;
    
    fhPtDecayNoIso  = new TH1F("hPtDecayNoIso","Number of not isolated leading pi0 decay particles",nptbins,ptmin,ptmax); 
    fhPtNoIso->SetYTitle("N");
    fhPtNoIso->SetXTitle("p_{T}(GeV/c)");
    outputContainer->Add(fhPtDecayNoIso) ;
    
    if(IsDataMC()){
      
      fhPtIsoPrompt  = new TH1F("hPtMCPrompt","Number of isolated prompt #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoPrompt->SetYTitle("N");
      fhPtIsoPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPrompt) ; 
      
      fhPhiIsoPrompt  = new TH2F
      ("hPhiMCPrompt","Number of isolated prompt #gamma",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPrompt->SetYTitle("#phi");
      fhPhiIsoPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPrompt) ; 
      
      fhEtaIsoPrompt  = new TH2F
      ("hEtaMCPrompt","Number of isolated prompt #gamma",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPrompt->SetYTitle("#eta");
      fhEtaIsoPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPrompt) ;
      
      fhPtIsoFragmentation  = new TH1F("hPtMCFragmentation","Number of isolated #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoFragmentation->SetYTitle("N");
      fhPtIsoFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoFragmentation) ; 
      
      fhPhiIsoFragmentation  = new TH2F
      ("hPhiMCFragmentation","Number of isolated fragmentation #gamma",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoFragmentation->SetYTitle("#phi");
      fhPhiIsoFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoFragmentation) ; 
      
      fhEtaIsoFragmentation  = new TH2F
      ("hEtaMCFragmentation","Number of isolated fragmentation #gamma",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoFragmentation->SetYTitle("#eta");
      fhEtaIsoFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoFragmentation) ;
      
      fhPtIsoPi0Decay  = new TH1F("hPtMCPi0Decay","Number of isolated #gamma from #pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtIsoPi0Decay->SetYTitle("N");
      fhPtIsoPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoPi0Decay) ; 
      
      fhPhiIsoPi0Decay  = new TH2F
      ("hPhiMCPi0Decay","Number of isolated #gamma from #pi^{0} decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoPi0Decay->SetYTitle("#phi");
      fhPhiIsoPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoPi0Decay) ; 
      
      fhEtaIsoPi0Decay  = new TH2F
      ("hEtaMCPi0Decay","Number of isolated #gamma from #pi^{0} decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoPi0Decay->SetYTitle("#eta");
      fhEtaIsoPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoPi0Decay) ;
      
      fhPtIsoEtaDecay  = new TH1F("hPtMCEtaDecay","Number of isolated #gamma from #eta decay",nptbins,ptmin,ptmax); 
      fhPtIsoEtaDecay->SetYTitle("N");
      fhPtIsoEtaDecay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoEtaDecay) ; 
      
      fhPhiIsoEtaDecay  = new TH2F
      ("hPhiMCEtaDecay","Number of isolated #gamma from #eta decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoEtaDecay->SetYTitle("#phi");
      fhPhiIsoEtaDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoEtaDecay) ; 
      
      fhEtaIsoEtaDecay  = new TH2F
      ("hEtaMCEtaDecay","Number of isolated #gamma from #eta decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoEtaDecay->SetYTitle("#eta");
      fhEtaIsoEtaDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoEtaDecay) ;
      
      fhPtIsoOtherDecay  = new TH1F("hPtMCOtherDecay","Number of isolated #gamma from non-#pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtIsoOtherDecay->SetYTitle("N");
      fhPtIsoOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoOtherDecay) ; 
      
      fhPhiIsoOtherDecay  = new TH2F
      ("hPhiMCOtherDecay","Number of isolated #gamma from non-#pi^{0} decay",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoOtherDecay->SetYTitle("#phi");
      fhPhiIsoOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoOtherDecay) ; 
      
      fhEtaIsoOtherDecay  = new TH2F
      ("hEtaMCOtherDecay","Number of isolated #gamma non-#pi^{0} decay",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoOtherDecay->SetYTitle("#eta");
      fhEtaIsoOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoOtherDecay) ;
      
      fhPtIsoConversion  = new TH1F("hPtMCConversion","Number of isolated converted #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoConversion->SetYTitle("N");
      fhPtIsoConversion->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtIsoConversion) ; 
      
      fhPhiIsoConversion  = new TH2F
      ("hPhiMCConversion","Number of isolated converted #gamma",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoConversion->SetYTitle("#phi");
      fhPhiIsoConversion->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhPhiIsoConversion) ; 
      
      fhEtaIsoConversion  = new TH2F
      ("hEtaMCConversion","Number of isolated converted #gamma",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoConversion->SetYTitle("#eta");
      fhEtaIsoConversion->SetXTitle("p_{T #gamma} (GeV/c)");
      outputContainer->Add(fhEtaIsoConversion) ;
      
      fhPtIsoUnknown  = new TH1F("hPtMCUnknown","Number of isolated non-#gamma particles",nptbins,ptmin,ptmax); 
      fhPtIsoUnknown->SetYTitle("N");
      fhPtIsoUnknown->SetXTitle("p_{T}(GeV/c)");
      outputContainer->Add(fhPtIsoUnknown) ; 
      
      fhPhiIsoUnknown  = new TH2F
      ("hPhiMCUnknown","Number of isolated non-#gamma particles",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiIsoUnknown->SetYTitle("#phi");
      fhPhiIsoUnknown->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPhiIsoUnknown) ; 
      
      fhEtaIsoUnknown  = new TH2F
      ("hEtaMCUnknown","Number of isolated non-#gamma particles",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhEtaIsoUnknown->SetYTitle("#eta");
      fhEtaIsoUnknown->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhEtaIsoUnknown) ;
      
      fhPtNoIsoPi0Decay  = new TH1F
      ("hPtNoIsoPi0Decay","Number of not isolated leading #gamma from #pi^{0} decay",nptbins,ptmin,ptmax); 
      fhPtNoIsoPi0Decay->SetYTitle("N");
      fhPtNoIsoPi0Decay->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtNoIsoPi0Decay) ;
      
      fhPtNoIsoEtaDecay  = new TH1F
      ("hPtNoIsoEtaDecay","Number of not isolated leading #gamma from eta decay",nptbins,ptmin,ptmax); 
      fhPtNoIsoEtaDecay->SetYTitle("N");
      fhPtNoIsoEtaDecay->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtNoIsoEtaDecay) ;
      
      fhPtNoIsoOtherDecay  = new TH1F
      ("hPtNoIsoOtherDecay","Number of not isolated leading #gamma from other decay",nptbins,ptmin,ptmax); 
      fhPtNoIsoOtherDecay->SetYTitle("N");
      fhPtNoIsoOtherDecay->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtNoIsoOtherDecay) ;
      
      fhPtNoIsoPrompt  = new TH1F
      ("hPtNoIsoPrompt","Number of not isolated leading prompt #gamma",nptbins,ptmin,ptmax); 
      fhPtNoIsoPrompt->SetYTitle("N");
      fhPtNoIsoPrompt->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtNoIsoPrompt) ;
      
      fhPtIsoMCPhoton  = new TH1F
      ("hPtIsoMCPhoton","Number of isolated leading  #gamma",nptbins,ptmin,ptmax); 
      fhPtIsoMCPhoton->SetYTitle("N");
      fhPtIsoMCPhoton->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtIsoMCPhoton) ;
      
      fhPtNoIsoMCPhoton  = new TH1F
      ("hPtNoIsoMCPhoton","Number of not isolated leading #gamma",nptbins,ptmin,ptmax); 
      fhPtNoIsoMCPhoton->SetYTitle("N");
      fhPtNoIsoMCPhoton->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtNoIsoMCPhoton) ;
    }//Histos with MC
    
  }
  
  if(fMakeSeveralIC){
    const Int_t buffersize = 255;
		char name[buffersize];
		char title[buffersize];
		for(Int_t icone = 0; icone<fNCones; icone++){
		  snprintf(name, buffersize,"hPtSum_Cone_%d",icone);
		  snprintf(title, buffersize,"Candidate cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		  fhPtSumIsolated[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		  fhPtSumIsolated[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		  fhPtSumIsolated[icone]->SetXTitle("p_{T} (GeV/c)");
		  outputContainer->Add(fhPtSumIsolated[icone]) ; 
		  
		  if(IsDataMC()){
		    snprintf(name, buffersize,"hPtSumPrompt_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate Prompt cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedPrompt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedPrompt[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedPrompt[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedPrompt[icone]) ; 
		    
		    snprintf(name, buffersize,"hPtSumFragmentation_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate Fragmentation cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedFragmentation[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedFragmentation[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedFragmentation[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedFragmentation[icone]) ; 
		    
		    snprintf(name, buffersize,"hPtSumPi0Decay_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate Pi0Decay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedPi0Decay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedPi0Decay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedPi0Decay[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedPi0Decay[icone]) ; 
		    
        snprintf(name, buffersize,"hPtSumEtaDecay_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate EtaDecay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedEtaDecay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedEtaDecay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedEtaDecay[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedEtaDecay[icone]) ;         
        
		    snprintf(name, buffersize,"hPtSumOtherDecay_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate OtherDecay cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedOtherDecay[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedOtherDecay[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedOtherDecay[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedOtherDecay[icone]) ; 
		    
		    snprintf(name, buffersize,"hPtSumConversion_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate Conversion cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedConversion[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedConversion[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedConversion[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedConversion[icone]) ; 
		    
		    snprintf(name, buffersize,"hPtSumUnknown_Cone_%d",icone);
		    snprintf(title, buffersize,"Candidate Unknown cone sum p_{T} for cone size %d vs candidate p_{T}",icone);
		    fhPtSumIsolatedUnknown[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
		    fhPtSumIsolatedUnknown[icone]->SetYTitle("#Sigma p_{T} (GeV/c)");
		    fhPtSumIsolatedUnknown[icone]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtSumIsolatedUnknown[icone]) ; 
		    
		  }//Histos with MC
		  
		  for(Int_t ipt = 0; ipt<fNPtThresFrac;ipt++){ 
		    snprintf(name, buffersize,"hPtThres_Cone_%d_Pt%d",icone,ipt);
		    snprintf(title, buffersize,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		    fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		    fhPtThresIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtThresIsolated[icone][ipt]) ; 
		    
		    snprintf(name, buffersize,"hPtFrac_Cone_%d_Pt%d",icone,ipt);
		    snprintf(title, buffersize,"Isolated candidate p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		    fhPtFracIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		    fhPtFracIsolated[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		    outputContainer->Add(fhPtFracIsolated[icone][ipt]) ; 
		    
		    if(IsDataMC()){
		      snprintf(name, buffersize,"hPtThresMCPrompt_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedPrompt[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedPrompt[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCPrompt_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Prompt p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedPrompt[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedPrompt[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedPrompt[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtThresMCFragmentation_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedFragmentation[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCFragmentation_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Fragmentation p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedFragmentation[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedFragmentation[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedFragmentation[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtThresMCPi0Decay_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedPi0Decay[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCPi0Decay_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Pi0Decay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedPi0Decay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedPi0Decay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedPi0Decay[icone][ipt]) ; 
		      
          snprintf(name, buffersize,"hPtThresMCEtaDecay_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate EtaDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedEtaDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedEtaDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedEtaDecay[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCEtaDecay_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate EtaDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedEtaDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedEtaDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedEtaDecay[icone][ipt]) ; 
          
          
		      snprintf(name, buffersize,"hPtThresMCOtherDecay_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedOtherDecay[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCOtherDecay_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate OtherDecay p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedOtherDecay[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedOtherDecay[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedOtherDecay[icone][ipt]) ;
		      
		      snprintf(name, buffersize,"hPtThresMCConversion_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedConversion[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedConversion[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCConversion_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Conversion p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedConversion[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedConversion[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedConversion[icone][ipt]) ;
		      
		      snprintf(name, buffersize,"hPtThresMCUnknown_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Unknown p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtThresIsolatedUnknown[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtThresIsolatedUnknown[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtThresIsolatedUnknown[icone][ipt]) ; 
		      
		      snprintf(name, buffersize,"hPtFracMCUnknown_Cone_%d_Pt%d",icone,ipt);
		      snprintf(title, buffersize,"Isolated candidate Unknown p_{T} distribution for cone size %d and p_{T}^{th} %d",icone,ipt);
		      fhPtFracIsolatedUnknown[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
		      fhPtFracIsolatedUnknown[icone][ipt]->SetXTitle("p_{T} (GeV/c)");
		      outputContainer->Add(fhPtFracIsolatedUnknown[icone][ipt]) ;  
		      
		    }//Histos with MC
		    
		  }//icone loop
		}//ipt loop
  }
  
  
  return outputContainer ;
  
}

//____________________________________________
void AliAnaParticleIsolation::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("IsolationCone");  
  AddToHistogramsName("AnaIsolation_");
  
  fCalorimeter = "PHOS" ;
  fReMakeIC = kFALSE ;
  fMakeSeveralIC = kFALSE ;
  
  //----------- Several IC-----------------
  fNCones          = 5 ; 
  fNPtThresFrac    = 5 ; 
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

//__________________________________________________
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
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0 ;
  TObjArray * pl    = 0x0; ; 
  
  //Select the calorimeter for candidate isolation with neutral particles
  if(fCalorimeter == "PHOS")
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  //Loop on AOD branch, filled previously in AliAnaPhoton, find leading particle to do isolation only with it
  Double_t ptLeading = 0. ;
  Int_t    idLeading = -1 ;
  TLorentzVector mom ;
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - Input aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod; iaod++){
    AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    //If too small or too large pt, skip
    if(aodinput->Pt() < GetMinPt() || aodinput->Pt() > GetMaxPt() ) continue ; 
    
    //check if it is low pt trigger particle, then adjust the isolation method
    if(aodinput->Pt() < GetIsolationCut()->GetPtThreshold() || aodinput->Pt() < GetIsolationCut()->GetSumPtThreshold()) 
      continue ; //trigger should not come from underlying event
    
    //vertex cut in case of mixing
    if (GetMixedEvent()) {
      Int_t evt=-1;
      Int_t id =-1;
      if     (aodinput->GetCaloLabel(0)  >= 0 ){
        id=aodinput->GetCaloLabel(0); 
        if(id >= 0 )evt= GetMixedEvent()-> EventIndexForCaloCluster(id) ;
      }
      else if(aodinput->GetTrackLabel(0) >= 0 ){
        id=aodinput->GetTrackLabel(0);
        if(id >= 0 )evt= GetMixedEvent()->EventIndex(id) ;
      }
      else continue;
      
      if (TMath::Abs(GetVertex(evt)[2]) > GetZvertexCut()) 
        return ;
    }
    
    //find the leading particles with highest momentum
    if ((aodinput->Pt())>ptLeading) {
      ptLeading = aodinput->Pt() ;
      idLeading = iaod ;
    }
    aodinput->SetLeadingParticle(kFALSE);
  }//finish searching for leading trigger particle
  
  // Check isolation of leading particle
  if(idLeading < 0) return;
  
  AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(idLeading));
  aodinput->SetLeadingParticle(kTRUE);
    
  //After cuts, study isolation
  n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
  GetIsolationCut()->MakeIsolationCut(GetCTSTracks(),pl,GetReader(), kTRUE, aodinput, GetAODObjArrayName(), n,nfrac,coneptsum, isolated);
  aodinput->SetIsolated(isolated);
  
  if(GetDebug() > 1) {
    if(isolated)printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() : Particle %d IS ISOLATED \n",idLeading);
    printf("AliAnaParticleIsolation::MakeAnalysisFillAOD() - End fill AODs \n");  
  }
  
}

//_________________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  Int_t   n = 0, nfrac = 0;
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0 ;
  
  //Loop on stored AOD 
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Histo aod branch entries %d\n", naod);
  
  //Get vertex for photon momentum calculation
  Double_t vertex[]={0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
	  GetReader()->GetVertex(vertex);
  }	
	
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4ParticleCorrelation* aod =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(!aod->IsLeadingParticle()) continue; // Try to isolate only leading cluster or track
    
    Bool_t  isolation  = aod->IsIsolated(); 
    Bool_t  decay      = aod->IsTagged();
    Float_t pt         = aod->Pt();
    Float_t phi        = aod->Phi();
    Float_t eta        = aod->Eta();
    Float_t conesize   = GetIsolationCut()->GetConeSize();
    
    //Recover reference arrays with clusters and tracks
    TObjArray * refclusters = aod->GetObjArray(GetAODObjArrayName()+"Clusters");
    TObjArray * reftracks   = aod->GetObjArray(GetAODObjArrayName()+"Tracks");
    //If too small or too large pt, skip
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ; 
    
    // --- In case of redoing isolation from delta AOD ----
    if(fMakeSeveralIC) {
      //Analysis of multiple IC at same time
      MakeSeveralICAnalysis(aod);
    }
    else if(fReMakeIC){
      //In case a more strict IC is needed in the produced AOD
      n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0;
      GetIsolationCut()->MakeIsolationCut(reftracks, refclusters, GetReader(), kFALSE, aod, "", n,nfrac,coneptsum, isolated);
      fhConeSumPt->Fill(pt,coneptsum);    
      if(GetDebug() > 0) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Energy Sum in Isolation Cone %2.2f\n", coneptsum);    
    }
    // --- -------------------------------------------- ----
    
    //Fill pt distribution of particles in cone
    //Tracks
    coneptsum = 0;
    Double_t sumptFR = 0. ;
    TObjArray * trackList   = GetCTSTracks() ;
    for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++){
      AliVTrack* track = (AliVTrack *) trackList->At(itrack);
      //fill the histograms at forward range
      if(!track){
        printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Track not available?");
        continue;
      }
      Double_t dPhi = phi - track->Phi() + TMath::PiOver2();
      Double_t dEta = eta - track->Eta();
      Double_t arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < conesize){
        fhPtInFRCone->Fill(pt,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptFR+=track->Pt();
      }    
      dPhi = phi - track->Phi() - TMath::PiOver2();
      arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < conesize){
        fhPtInFRCone->Fill(pt,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        sumptFR+=track->Pt();
      }      
    }
    fhFRConeSumPt->Fill(pt,sumptFR);
    if(reftracks){  
      for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++){
        AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
        fhPtInCone->Fill(pt,TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py()));
        coneptsum+=track->Pt();
      }
    }
    
    //CaloClusters
    if(refclusters){    
      TLorentzVector mom ;
      for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++){
        AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
        calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
        
        fhPtInCone->Fill(pt, mom.Pt());
        coneptsum+=mom.Pt();
      }
    }
	  
    if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d Energy Sum in Isolation Cone %2.2f\n", iaod, coneptsum);
    
    if(!fReMakeIC) fhConeSumPt->Fill(pt,coneptsum);
    
    if(isolation){    
      
      if(GetDebug() > 1) printf("AliAnaParticleIsolation::MakeAnalysisFillHistograms() - Particle %d ISOLATED, fill histograms\n", iaod);
      
      fhPtIso  ->Fill(pt);
      fhPhiIso ->Fill(pt,phi);
      fhEtaIso ->Fill(pt,eta);
      if (decay) fhPtDecayIso->Fill(pt);
      
      if(IsDataMC()){
        Int_t tag =aod->GetTag();
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
        {
          fhPtIsoMCPhoton  ->Fill(pt);
        }        
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)){
          fhPtIsoPrompt  ->Fill(pt);
          fhPhiIsoPrompt ->Fill(pt,phi);
          fhEtaIsoPrompt ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
        {
          fhPtIsoFragmentation  ->Fill(pt);
          fhPhiIsoFragmentation ->Fill(pt,phi);
          fhEtaIsoFragmentation ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhPtIsoPi0Decay  ->Fill(pt);
          fhPhiIsoPi0Decay ->Fill(pt,phi);
          fhEtaIsoPi0Decay ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))
        {
          fhPtIsoEtaDecay  ->Fill(pt);
          fhPhiIsoEtaDecay ->Fill(pt,phi);
          fhEtaIsoEtaDecay ->Fill(pt,eta);
        }        
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
        {
          fhPtIsoOtherDecay  ->Fill(pt);
          fhPhiIsoOtherDecay ->Fill(pt,phi);
          fhEtaIsoOtherDecay ->Fill(pt,eta);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
        {
          fhPtIsoConversion  ->Fill(pt);
          fhPhiIsoConversion ->Fill(pt,phi);
          fhEtaIsoConversion ->Fill(pt,eta);
        }
        else // anything else
        {
          fhPtIsoUnknown  ->Fill(pt);
          fhPhiIsoUnknown ->Fill(pt,phi);
          fhEtaIsoUnknown ->Fill(pt,eta);
        }
      }//Histograms with MC
      
    }//Isolated histograms
    
    if(!isolation)
    {
      fhPtNoIso  ->Fill(pt);
      if (decay) fhPtDecayNoIso->Fill(pt);
      
      if(IsDataMC()){
        Int_t tag =aod->GetTag();
        
        if     (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
        {
          fhPtNoIsoMCPhoton->Fill(pt);
        }
        
        if     (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhPtNoIsoPi0Decay->Fill(pt);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))
        {
          fhPtNoIsoEtaDecay->Fill(pt);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
        {
          fhPtNoIsoOtherDecay->Fill(pt);
        }        
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))
        {
          fhPtNoIsoPrompt->Fill(pt);
        }
        
      }
    }
    
  }// aod loop
  
}


//_____________________________________________________________________________________
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
          if     ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtThresIsolatedPrompt[icone][ipt]       ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtThresIsolatedConversion[icone][ipt]   ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtThresIsolatedFragmentation[icone][ipt]->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtThresIsolatedPi0Decay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtThresIsolatedEtaDecay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtThresIsolatedOtherDecay[icone][ipt]   ->Fill(ptC) ;
          else  fhPtThresIsolatedUnknown[icone][ipt]->Fill(ptC) ;
        }
      }
      
      //Pt threshold on pt cand/ pt in cone fraction
      if(nfrac[icone][ipt] == 0) {
        fhPtFracIsolated[icone][ipt]->Fill(ptC);
        if(IsDataMC()){
          if     ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtFracIsolatedPrompt[icone][ipt]       ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtFracIsolatedConversion[icone][ipt]   ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtFracIsolatedFragmentation[icone][ipt]->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtFracIsolatedPi0Decay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtFracIsolatedEtaDecay[icone][ipt]     ->Fill(ptC) ;
          else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtFracIsolatedOtherDecay[icone][ipt]   ->Fill(ptC) ;
          else  fhPtFracIsolatedUnknown[icone][ipt]->Fill(ptC) ;
        }
      }
    }//pt thresh loop
    
    //Sum in cone histograms
    fhPtSumIsolated[icone]->Fill(ptC,coneptsum) ;
    if(IsDataMC()){
      if     ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))        fhPtSumIsolatedPrompt[icone]       ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))    fhPtSumIsolatedConversion[icone]   ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)) fhPtSumIsolatedFragmentation[icone]->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))      fhPtSumIsolatedPi0Decay[icone]     ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))      fhPtSumIsolatedEtaDecay[icone]     ->Fill(ptC,coneptsum) ;
      else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))    fhPtSumIsolatedOtherDecay[icone]   ->Fill(ptC,coneptsum) ;
      else  fhPtSumIsolatedUnknown[icone]->Fill(ptC,coneptsum) ;
    }
    
  }//cone size loop
  
  //Reset original parameters for AOD analysis
  GetIsolationCut()->SetPtThreshold(ptthresorg);
  GetIsolationCut()->SetPtFraction(ptfracorg);
  GetIsolationCut()->SetConeSize(rorg);
  
}

//_____________________________________________________________
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

