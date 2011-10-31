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

//_________________________________________________________________________
//
// Conversions pairs analysis
// Check if cluster comes from a conversion in the material in front of the calorimeter
// Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
// Input are selected clusters with AliAnaPhoton
//
//
//-- Author: Gustavo Conesa (LPSC-IN2P3-CNRS)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system --- 
#include <TH2F.h>
#include <TH3D.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include "TParticle.h"
#include "TDatabasePDG.h"

// --- Analysis system --- 
#include "AliAnaPhotonConvInCalo.h" 
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnaPhotonConvInCalo)

//________________________________________
AliAnaPhotonConvInCalo::AliAnaPhotonConvInCalo() : 
AliAnaPartCorrBaseClass(),   
fRemoveConvertedPair(kFALSE), 
fAddConvertedPairsToAOD(kFALSE), 
fMassCut(0),                  
fConvAsymCut(1.),                  fConvDEtaCut(2.),
fConvDPhiMinCut(-1.),              fConvDPhiMaxCut(7.), 

// Histograms
fhPtPhotonConv(0),                 fhEtaPhiPhotonConv(0),          fhEtaPhi05PhotonConv(0),
fhConvDeltaEta(0),                 fhConvDeltaPhi(0),              fhConvDeltaEtaPhi(0), 
fhConvAsym(0),                     fhConvPt(0),
fhConvDistEta(0),                  fhConvDistEn(0),                fhConvDistMass(0),     
fhConvDistEtaCutEta(0),            fhConvDistEnCutEta(0),          fhConvDistMassCutEta(0),
fhConvDistEtaCutMass(0),           fhConvDistEnCutMass(0), 
fhConvDistEtaCutAsy(0),            fhConvDistEnCutAsy(0),

// MC histograms
fhPtConversionTagged(0),           fhPtAntiNeutronTagged(0),       
fhPtAntiProtonTagged(0),           fhPtUnknownTagged(0),

fhConvDeltaEtaMCConversion(0),     fhConvDeltaPhiMCConversion(0),  fhConvDeltaEtaPhiMCConversion(0),
fhConvAsymMCConversion(0),         fhConvPtMCConversion(0),           
fhConvDispersionMCConversion(0),   fhConvM02MCConversion(0),

fhConvDeltaEtaMCAntiNeutron(0),    fhConvDeltaPhiMCAntiNeutron(0), fhConvDeltaEtaPhiMCAntiNeutron(0), 
fhConvAsymMCAntiNeutron(0),        fhConvPtMCAntiNeutron(0), 
fhConvDispersionMCAntiNeutron(0),  fhConvM02MCAntiNeutron(0),
fhConvDeltaEtaMCAntiProton(0),     fhConvDeltaPhiMCAntiProton(0),  fhConvDeltaEtaPhiMCAntiProton(0),  
fhConvAsymMCAntiProton(0),         fhConvPtMCAntiProton(0),  
fhConvDispersionMCAntiProton(0),   fhConvM02MCAntiProton(0),
fhConvDeltaEtaMCString(0),         fhConvDeltaPhiMCString(0),      fhConvDeltaEtaPhiMCString(0),      
fhConvAsymMCString(0),             fhConvPtMCString(0),      
fhConvDispersionMCString(0),       fhConvM02MCString(0),
fhConvDistMCConversion(0),         fhConvDistMCConversionCuts(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
}

//_________________________________________________
TObjString *  AliAnaPhotonConvInCalo::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPhotonConvInCalo---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Conversion Selection: fConvAsymCut %1.2f, fConvDEtaCut %1.2f fConvDPhiCut (%1.2f,%1.2f)\n",
           fConvAsymCut, fConvDEtaCut, fConvDPhiMinCut, fConvDPhiMaxCut) ;
  parList+=onePar ; 
  
  return new TObjString(parList) ;
}

//___________________________________________________
TList *  AliAnaPhotonConvInCalo::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("PhotonConvInCaloHistos") ; 
	
  Int_t nptbins  = GetHistoPtBins();  Float_t ptmax  = GetHistoPtMax();  Float_t ptmin  = GetHistoPtMin(); 
  Int_t nphibins = GetHistoPhiBins(); Float_t phimax = GetHistoPhiMax(); Float_t phimin = GetHistoPhiMin(); 
  Int_t netabins = GetHistoEtaBins(); Float_t etamax = GetHistoEtaMax(); Float_t etamin = GetHistoEtaMin();	
  
  fhPtPhotonConv  = new TH1F("hPtPhotonConv","Number of #gamma over calorimeter, conversion",nptbins,ptmin,ptmax); 
  fhPtPhotonConv->SetYTitle("N");
  fhPtPhotonConv->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhPtPhotonConv) ; 
  
  fhEtaPhiPhotonConv  = new TH2F
  ("hEtaPhiPhotonConv","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
  fhEtaPhiPhotonConv->SetYTitle("#phi (rad)");
  fhEtaPhiPhotonConv->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiPhotonConv) ;
  if(GetMinPt() < 0.5){
    fhEtaPhi05PhotonConv  = new TH2F
    ("hEtaPhi05PhotonConv","#eta vs #phi, E > 0.5",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi05PhotonConv->SetYTitle("#phi (rad)");
    fhEtaPhi05PhotonConv->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi05PhotonConv) ;
  }
  
  fhConvDeltaEta  = new TH2F
  ("hConvDeltaEta","#Delta #eta of selected conversion pairs",100,0,fMassCut,netabins*2,-0.5,0.5); 
  fhConvDeltaEta->SetYTitle("#Delta #eta");
  fhConvDeltaEta->SetXTitle("Pair Mass (GeV/c^2)");
  outputContainer->Add(fhConvDeltaEta) ;
  
  fhConvDeltaPhi  = new TH2F
  ("hConvDeltaPhi","#Delta #phi of selected conversion pairs",100,0,fMassCut,nphibins*2,-0.5,0.5); 
  fhConvDeltaPhi->SetYTitle("#Delta #phi");
  fhConvDeltaPhi->SetXTitle("Pair Mass (GeV/c^2)");
  outputContainer->Add(fhConvDeltaPhi) ;
  
  fhConvDeltaEtaPhi  = new TH2F
  ("hConvDeltaEtaPhi","#Delta #eta vs #Delta #phi of selected conversion pairs",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
  fhConvDeltaEtaPhi->SetYTitle("#Delta #phi");
  fhConvDeltaEtaPhi->SetXTitle("#Delta #eta");
  outputContainer->Add(fhConvDeltaEtaPhi) ;
  
  fhConvAsym  = new TH2F
  ("hConvAsym","Asymmetry of selected conversion pairs",100,0,fMassCut,100,0,1); 
  fhConvAsym->SetYTitle("Asymmetry");
  fhConvAsym->SetXTitle("Pair Mass (GeV/c^2)");
  outputContainer->Add(fhConvAsym) ;  
  
  fhConvPt  = new TH2F
  ("hConvPt","p_{T} of selected conversion pairs",100,0,fMassCut,100,0.,10.); 
  fhConvPt->SetYTitle("Pair p_{T} (GeV/c)");
  fhConvPt->SetXTitle("Pair Mass (GeV/c^2)");
  outputContainer->Add(fhConvPt) ;
  
  fhConvDistEta  = new TH2F
  ("hConvDistEta","distance to conversion vertex",100,-0.7,0.7,100,0.,5.); 
  fhConvDistEta->SetXTitle("#eta");
  fhConvDistEta->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEta) ;
  
  fhConvDistEn  = new TH2F
  ("hConvDistEn","distance to conversion vertex",nptbins,ptmin,ptmax,100,0.,5.); 
  fhConvDistEn->SetXTitle("E (GeV)");
  fhConvDistEn->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEn) ;
  
  fhConvDistMass  = new TH2F
  ("hConvDistMass","distance to conversion vertex",100,0,fMassCut,100,0.,5.); 
  fhConvDistMass->SetXTitle("m (GeV/c^2)");
  fhConvDistMass->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistMass) ;
  
  fhConvDistEtaCutEta  = new TH2F
  ("hConvDistEtaCutEta","distance to conversion vertex, dEta < 0.05",100,-0.7,0.7,100,0.,5.); 
  fhConvDistEtaCutEta->SetXTitle("#eta");
  fhConvDistEtaCutEta->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEtaCutEta) ;
  
  fhConvDistEnCutEta  = new TH2F
  ("hConvDistEnCutEta","distance to conversion vertex, dEta < 0.05",nptbins,ptmin,ptmax,100,0.,5.); 
  fhConvDistEnCutEta->SetXTitle("E (GeV)");
  fhConvDistEnCutEta->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEnCutEta) ;
  
  fhConvDistMassCutEta  = new TH2F
  ("hConvDistMassCutEta","distance to conversion vertex, dEta < 0.05",100,0,fMassCut,100,0.,5.); 
  fhConvDistMassCutEta->SetXTitle("m (GeV/c^2)");
  fhConvDistMassCutEta->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistMassCutEta) ;
  
  fhConvDistEtaCutMass  = new TH2F
  ("hConvDistEtaCutMass","distance to conversion vertex, dEta < 0.05, m < 10 MeV",100,-0.7,0.7,100,0.,5.); 
  fhConvDistEtaCutMass->SetXTitle("#eta");
  fhConvDistEtaCutMass->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEtaCutMass) ;
  
  fhConvDistEnCutMass  = new TH2F
  ("hConvDistEnCutMass","distance to conversion vertex, dEta < 0.05, m < 10 MeV",nptbins,ptmin,ptmax,100,0.,5.); 
  fhConvDistEnCutMass->SetXTitle("E (GeV)");
  fhConvDistEnCutMass->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEnCutMass) ;
  
  fhConvDistEtaCutAsy  = new TH2F
  ("hConvDistEtaCutAsy","distance to conversion vertex, dEta < 0.05, m < 10 MeV, A < 0.1",100,-0.7,0.7,100,0.,5.); 
  fhConvDistEtaCutAsy->SetXTitle("#eta");
  fhConvDistEtaCutAsy->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEtaCutAsy) ;
  
  fhConvDistEnCutAsy  = new TH2F
  ("hConvDistEnCutAsy","distance to conversion vertex, dEta < 0.05, m < 10 MeV, A < 0.1",nptbins,ptmin,ptmax,100,0.,5.); 
  fhConvDistEnCutAsy->SetXTitle("E (GeV)");
  fhConvDistEnCutAsy->SetYTitle(" distance (m)");
  outputContainer->Add(fhConvDistEnCutAsy) ;
  
  if(IsDataMC()){
    
    fhPtConversionTagged  = new TH1F("hPtMCConversionTagged","Number of converted #gamma over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
    fhPtConversionTagged->SetYTitle("N");
    fhPtConversionTagged->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtConversionTagged) ; 
    
    
    fhPtAntiNeutronTagged  = new TH1F("hPtMCAntiNeutronTagged","Number of AntiNeutron id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
    fhPtAntiNeutronTagged->SetYTitle("N");
    fhPtAntiNeutronTagged->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtAntiNeutronTagged) ; 
    
    fhPtAntiProtonTagged  = new TH1F("hPtMCAntiProtonTagged","Number of AntiProton id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
    fhPtAntiProtonTagged->SetYTitle("N");
    fhPtAntiProtonTagged->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtAntiProtonTagged) ; 
    
    fhPtUnknownTagged  = new TH1F("hPtMCUnknownTagged","Number of Unknown id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
    fhPtUnknownTagged->SetYTitle("N");
    fhPtUnknownTagged->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtUnknownTagged) ;     
    
    fhConvDeltaEtaMCConversion  = new TH2F
    ("hConvDeltaEtaMCConversion","#Delta #eta of selected conversion pairs from real conversions",100,0,fMassCut,netabins,-0.5,0.5); 
    fhConvDeltaEtaMCConversion->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaEtaMCConversion) ;
    
    fhConvDeltaPhiMCConversion  = new TH2F
    ("hConvDeltaPhiMCConversion","#Delta #phi of selected conversion pairs from real conversions",100,0,fMassCut,nphibins,-0.5,0.5); 
    fhConvDeltaPhiMCConversion->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaPhiMCConversion) ;
    
    fhConvDeltaEtaPhiMCConversion  = new TH2F
    ("hConvDeltaEtaPhiMCConversion","#Delta #eta vs #Delta #phi of selected conversion pairs, from real conversions",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
    fhConvDeltaEtaPhiMCConversion->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCConversion->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCConversion) ;
    
    fhConvAsymMCConversion  = new TH2F
    ("hConvAsymMCConversion","Asymmetry of selected conversion pairs from real conversions",100,0,fMassCut,100,0,1); 
    fhConvAsymMCConversion->SetYTitle("Asymmetry");
    fhConvAsymMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvAsymMCConversion) ;
    
    fhConvPtMCConversion  = new TH2F
    ("hConvPtMCConversion","p_{T} of selected conversion pairs from real conversions",100,0,fMassCut,100,0.,10.); 
    fhConvPtMCConversion->SetYTitle("Pair p_{T} (GeV/c)");
    fhConvPtMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvPtMCConversion) ;    
    
    fhConvDispersionMCConversion  = new TH2F
    ("hConvDispersionMCConversion","p_{T} of selected conversion pairs from real conversions",100,0.,1.,100,0.,1.); 
    fhConvDispersionMCConversion->SetYTitle("Dispersion cluster 1");
    fhConvDispersionMCConversion->SetXTitle("Dispersion cluster 2");
    outputContainer->Add(fhConvDispersionMCConversion) ;   
    
    fhConvM02MCConversion  = new TH2F
    ("hConvM02MCConversion","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
    fhConvM02MCConversion->SetYTitle("M02 cluster 1");
    fhConvM02MCConversion->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCConversion) ;           
    
    fhConvDeltaEtaMCAntiNeutron  = new TH2F
    ("hConvDeltaEtaMCAntiNeutron","#Delta #eta of selected conversion pairs from anti-neutrons",100,0,fMassCut,netabins,-0.5,0.5); 
    fhConvDeltaEtaMCAntiNeutron->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaEtaMCAntiNeutron) ;
    
    fhConvDeltaPhiMCAntiNeutron  = new TH2F
    ("hConvDeltaPhiMCAntiNeutron","#Delta #phi of selected conversion pairs from anti-neutrons",100,0,fMassCut,nphibins,-0.5,0.5); 
    fhConvDeltaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaPhiMCAntiNeutron) ;
    
    fhConvDeltaEtaPhiMCAntiNeutron  = new TH2F
    ("hConvDeltaEtaPhiMCAntiNeutron","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-neutrons",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
    fhConvDeltaEtaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCAntiNeutron->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCAntiNeutron) ;    
    
    fhConvAsymMCAntiNeutron  = new TH2F
    ("hConvAsymMCAntiNeutron","Asymmetry of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0,1); 
    fhConvAsymMCAntiNeutron->SetYTitle("Asymmetry");
    fhConvAsymMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvAsymMCAntiNeutron) ;
    
    fhConvPtMCAntiNeutron  = new TH2F
    ("hConvPtMCAntiNeutron","p_{T} of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0.,10.); 
    fhConvPtMCAntiNeutron->SetYTitle("Pair p_{T} (GeV/c)");
    fhConvPtMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvPtMCAntiNeutron) ;    
    
    fhConvDispersionMCAntiNeutron  = new TH2F
    ("hConvDispersionMCAntiNeutron","p_{T} of selected conversion pairs from anti-neutrons",100,0.,1.,100,0.,1.); 
    fhConvDispersionMCAntiNeutron->SetYTitle("Dispersion cluster 1");
    fhConvDispersionMCAntiNeutron->SetXTitle("Dispersion cluster 2");
    outputContainer->Add(fhConvDispersionMCAntiNeutron) ;       
    
    fhConvM02MCAntiNeutron  = new TH2F
    ("hConvM02MCAntiNeutron","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
    fhConvM02MCAntiNeutron->SetYTitle("M02 cluster 1");
    fhConvM02MCAntiNeutron->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCAntiNeutron) ;  
    
    fhConvDeltaEtaMCAntiProton  = new TH2F
    ("hConvDeltaEtaMCAntiProton","#Delta #eta of selected conversion pairs from anti-protons",100,0,fMassCut,netabins,-0.5,0.5); 
    fhConvDeltaEtaMCAntiProton->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaEtaMCAntiProton) ;
    
    fhConvDeltaPhiMCAntiProton  = new TH2F
    ("hConvDeltaPhiMCAntiProton","#Delta #phi of selected conversion pairs from anti-protons",100,0,fMassCut,nphibins,-0.5,0.5); 
    fhConvDeltaPhiMCAntiProton->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaPhiMCAntiProton) ;
    
    fhConvDeltaEtaPhiMCAntiProton  = new TH2F
    ("hConvDeltaEtaPhiMCAntiProton","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-protons",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
    fhConvDeltaEtaPhiMCAntiProton->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCAntiProton->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCAntiProton) ;    
    
    fhConvAsymMCAntiProton  = new TH2F
    ("hConvAsymMCAntiProton","Asymmetry of selected conversion pairs from anti-protons",100,0,fMassCut,100,0,1); 
    fhConvAsymMCAntiProton->SetYTitle("Asymmetry");
    fhConvAsymMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvAsymMCAntiProton) ;
    
    fhConvPtMCAntiProton  = new TH2F
    ("hConvPtMCAntiProton","p_{T} of selected conversion pairs from anti-protons",100,0,fMassCut,100,0.,10.); 
    fhConvPtMCAntiProton->SetYTitle("Pair p_{T} (GeV/c)");
    fhConvPtMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvPtMCAntiProton) ;
    
    fhConvDispersionMCAntiProton  = new TH2F
    ("hConvDispersionMCAntiProton","p_{T} of selected conversion pairs from anti-protons",100,0.,1.,100,0.,1.); 
    fhConvDispersionMCAntiProton->SetYTitle("Dispersion cluster 1");
    fhConvDispersionMCAntiProton->SetXTitle("Dispersion cluster 2");
    outputContainer->Add(fhConvDispersionMCAntiProton) ;       
    
    fhConvM02MCAntiProton  = new TH2F
    ("hConvM02MCAntiProton","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
    fhConvM02MCAntiProton->SetYTitle("M02 cluster 1");
    fhConvM02MCAntiProton->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCAntiProton) ;       
    
    fhConvDeltaEtaMCString  = new TH2F
    ("hConvDeltaEtaMCString","#Delta #eta of selected conversion pairs from string",100,0,fMassCut,netabins,-0.5,0.5); 
    fhConvDeltaEtaMCString->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCString->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaEtaMCString) ;
    
    fhConvDeltaPhiMCString  = new TH2F
    ("hConvDeltaPhiMCString","#Delta #phi of selected conversion pairs from string",100,0,fMassCut,nphibins,-0.5,0.5); 
    fhConvDeltaPhiMCString->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCString->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaPhiMCString) ;
    
    fhConvDeltaEtaPhiMCString  = new TH2F
    ("hConvDeltaEtaPhiMCString","#Delta #eta vs #Delta #phi of selected conversion pairs from string",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
    fhConvDeltaEtaPhiMCString->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCString->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCString) ;    
    
    fhConvAsymMCString  = new TH2F
    ("hConvAsymMCString","Asymmetry of selected conversion pairs from string",100,0,fMassCut,100,0,1); 
    fhConvAsymMCString->SetYTitle("Asymmetry");
    fhConvAsymMCString->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvAsymMCString) ;
    
    fhConvPtMCString  = new TH2F
    ("hConvPtMCString","p_{T} of selected conversion pairs from string",100,0,fMassCut,100,0.,10.); 
    fhConvPtMCString->SetYTitle("Pair p_{T} (GeV/c)");
    fhConvPtMCString->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvPtMCString) ;
    
    fhConvDispersionMCString  = new TH2F
    ("hConvDispersionMCString","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
    fhConvDispersionMCString->SetYTitle("Dispersion cluster 1");
    fhConvDispersionMCString->SetXTitle("Dispersion cluster 2");
    outputContainer->Add(fhConvDispersionMCString) ;       
    
    fhConvM02MCString  = new TH2F
    ("hConvM02MCString","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
    fhConvM02MCString->SetYTitle("M02 cluster 1");
    fhConvM02MCString->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCString) ; 
    
    fhConvDistMCConversion  = new TH2F
    ("hConvDistMCConversion","calculated conversion distance vs real vertes for MC conversion",100,0.,5.,100,0.,5.); 
    fhConvDistMCConversion->SetYTitle("distance");
    fhConvDistMCConversion->SetXTitle("vertex R");
    outputContainer->Add(fhConvDistMCConversion) ; 
    
    fhConvDistMCConversionCuts  = new TH2F
    ("hConvDistMCConversionCuts","calculated conversion distance vs real vertes for MC conversion, deta < 0.05, m < 10 MeV, asym < 0.1",100,0.,5.,100,0.,5.); 
    fhConvDistMCConversionCuts->SetYTitle("distance");
    fhConvDistMCConversionCuts->SetXTitle("vertex R");
    outputContainer->Add(fhConvDistMCConversionCuts) ; 
        
  }
  
  return outputContainer ;

}

//_______________________________________
void AliAnaPhotonConvInCalo::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaPhotonConvInCalo_");
  
  fMassCut                = 0.03; //30 MeV
  fRemoveConvertedPair    = kFALSE;
  fAddConvertedPairsToAOD = kFALSE;
	
}

//_____________________________________________
void  AliAnaPhotonConvInCalo::MakeAnalysisFillAOD() 
{
  //Do conversion photon analysis and fill aods
  
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPhotonConvInCalo::MakeAnalysisFillAOD() - aod branch entries %d\n", naod);
  
  //List to be used in conversion analysis, to tag the cluster as candidate for conversion
  Bool_t * indexConverted = new Bool_t[naod];
  for (Int_t i = 0; i < naod; i++) indexConverted[i] = kFALSE;
	
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* calo =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    
    Bool_t bConverted = kFALSE;
    Int_t id2 = -1;
    
    //Check if set previously as converted couple, if so skip its use.
    if (indexConverted[iaod]) continue;
    
    // Second cluster loop
    AliAODPWG4Particle* calo2 = 0;
    for(Int_t jaod = iaod + 1 ; jaod < naod ; jaod++) {
      //Check if set previously as converted couple, if so skip its use.
      if (indexConverted[jaod]) continue;
      //printf("Check Conversion indeces %d and %d\n",iaod,jaod);
      calo2 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(jaod));
      
      //................................................
      //Get mass of pair, if small, take this pair.
      Float_t pairM     = calo->GetPairMass(calo2);
      //printf("\t both in calo, mass %f, cut %f\n",pairM,fMassCut);
      if(pairM < fMassCut){  
        calo->SetTagged(kFALSE);
        id2 = calo2->GetCaloLabel(0);
        Float_t asymmetry = TMath::Abs(calo->E()-calo2->E())/(calo->E()+calo2->E());
        Float_t dPhi      = (calo->Momentum())->Phi()-(calo2->Momentum())->Phi();
        Float_t dEta      = (calo->Momentum())->Eta()-(calo2->Momentum())->Eta();  
        
        //...............................................
        //Fill few histograms with kinematics of the pair
        //FIXME, move all this to MakeAnalysisFillHistograms ...
        
        fhConvDeltaEta   ->Fill( pairM, dPhi      );
        fhConvDeltaPhi   ->Fill( pairM, dEta      );
        fhConvAsym       ->Fill( pairM, asymmetry );
        fhConvDeltaEtaPhi->Fill( dEta , dPhi      );
        fhConvPt         ->Fill( pairM, (calo->Momentum())->Pt()+(calo2->Momentum())->Pt());          
        
        //Estimate conversion distance, T. Awes, M. Ivanov
        //Under the assumption that the pair has zero mass, and that each electron 
        //of the pair has the same momentum, they will each have the same bend radius 
        //given by R=p/(qB) = p / (300 B) with p in [MeV/c], B in [Tesla] and R in [m]. 
        //With nominal ALICE magnet current of 30kA B=0.5T, and so with E_cluster=p,  
        //R = E/1.5 [cm].  Under these assumptions, the distance from the conversion 
        //point to the MCEal can be related to the separation distance, L=2y, on the MCEal 
        //as d = sqrt(R^2 -(R-y)^2) = sqrt(2Ry - y^2). And since R>>y we can write as 
        //d = sqrt(E*L/1.5) where E is the cluster energy and L is the distance in cm between 
        //the clusters.
        
        TObjArray * clusters    = 0; 
        if(calo->GetDetector() == "EMCAL"){
          clusters = GetEMCALClusters();
        }
        else{
          clusters = GetPHOSClusters();
        }
        
        Int_t iclus = -1;
        AliVCluster *cluster1 = FindCluster(clusters,calo ->GetCaloLabel(0),iclus); 
        AliVCluster *cluster2 = FindCluster(clusters,calo2->GetCaloLabel(0),iclus); 

        Float_t pos1[3];
        cluster1->GetPosition(pos1); 
        Float_t pos2[3];
        cluster2->GetPosition(pos2); 
        Float_t clustDist = TMath::Sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+
                                        (pos1[1]-pos2[1])*(pos1[1]-pos2[1])+
                                        (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
        
        Float_t convDist  = TMath::Sqrt(calo->E() *clustDist*0.01/0.15);
        Float_t convDist2 = TMath::Sqrt(calo2->E()*clustDist*0.01/0.15);
        //printf("l = %f, e1 = %f, d1=%f, e2 = %f, d2=%f\n",clustDist,calo->E(),convDist,calo2->E(),convDist2);
        if(GetDebug() > 2)
          printf("AliAnaPhotonConvInCalo::MakeAnalysisFillAOD(): Pair with mass %2.3f < %2.3f, %1.2f < dPhi %2.2f < %2.2f, dEta %f < %2.2f, asymmetry %2.2f< %2.2f; \n    cluster1 id %d, e %2.3f  SM %d, eta %2.3f, phi %2.3f ; \n    cluster2 id %d, e %2.3f, SM %d,eta %2.3f, phi %2.3f\n",
                 pairM,fMassCut,fConvDPhiMinCut, dPhi, fConvDPhiMaxCut, dEta, fConvDEtaCut, asymmetry, fConvAsymCut,
                 calo->GetCaloLabel(0),calo->E(),GetCaloUtils()->GetModuleNumber(calo,GetReader()->GetInputEvent()), calo->Eta(), calo->Phi(),
                 id2, calo2->E(), GetCaloUtils()->GetModuleNumber(calo2,GetReader()->GetInputEvent()),calo2->Eta(), calo2->Phi());
        
        fhConvDistEta ->Fill(calo ->Eta(),convDist );
        fhConvDistEta ->Fill(calo2->Eta(),convDist2);
        fhConvDistEn  ->Fill(calo ->E(),  convDist );
        fhConvDistEn  ->Fill(calo2->E(),  convDist2);        
        fhConvDistMass->Fill(pairM, convDist );
        //dEta cut
        if(dEta<0.05){
          fhConvDistEtaCutEta ->Fill(calo->Eta(), convDist );
          fhConvDistEtaCutEta ->Fill(calo2->Eta(),convDist2);
          fhConvDistEnCutEta  ->Fill(calo->E(),   convDist );
          fhConvDistEnCutEta  ->Fill(calo2->E(),  convDist2);        
          fhConvDistMassCutEta->Fill(pairM, convDist );
          //mass cut
          if(pairM<0.01){//10 MeV
            fhConvDistEtaCutMass ->Fill(calo ->Eta(), convDist );
            fhConvDistEtaCutMass ->Fill(calo2->Eta(), convDist2);
            fhConvDistEnCutMass  ->Fill(calo ->E(),   convDist );
            fhConvDistEnCutMass  ->Fill(calo2->E(),   convDist2);        
            // asymmetry cut
            if(asymmetry<0.1){
              fhConvDistEtaCutAsy ->Fill(calo ->Eta(), convDist );
              fhConvDistEtaCutAsy ->Fill(calo2->Eta(), convDist2);
              fhConvDistEnCutAsy  ->Fill(calo ->E(),   convDist );
              fhConvDistEnCutAsy  ->Fill(calo2->E(),   convDist2); 
            }//asymmetry cut
          }//mass cut            
        }//dEta cut
        
        //...............................................
        //Select pairs in a eta-phi window
        if(TMath::Abs(dEta) < fConvDEtaCut    && 
           TMath::Abs(dPhi) < fConvDPhiMaxCut &&
           TMath::Abs(dPhi) > fConvDPhiMinCut && 
           asymmetry        < fConvAsymCut       ){
          indexConverted[iaod] = kTRUE;
          indexConverted[jaod] = kTRUE; 
          bConverted           = kTRUE;          
        }
        //printf("Accepted? %d\n",bConverted);
        //...........................................
        //Fill more histograms, simulated data
        //FIXME, move all this to MakeAnalysisFillHistograms ...
        if(IsDataMC()){
          
          //Check the origin of the pair, look for conversion, antinucleons or jet correlations (strings)
          Int_t ancPDG    = 0;
          Int_t ancStatus = 0;
          TLorentzVector momentum;
          TVector3 prodVertex;
          Int_t ancLabel  = GetMCAnalysisUtils()->CheckCommonAncestor(cluster1->GetLabel(), cluster2->GetLabel(), 
                                                                      GetReader(), ancPDG, ancStatus, momentum, prodVertex);
          
          // printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() - Common ancestor label %d, pdg %d, name %s, status %d; \n",
          //                          ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
          
          Int_t tag1 = calo ->GetTag();
          Int_t tag2 = calo2->GetTag();
          if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCConversion)){
            if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) && (ancPDG==22 || TMath::Abs(ancPDG)==11) && ancLabel > -1){
              fhConvDeltaEtaMCConversion   ->Fill( pairM, dEta      );
              fhConvDeltaPhiMCConversion   ->Fill( pairM, dPhi      );
              fhConvAsymMCConversion       ->Fill( pairM, asymmetry );
              fhConvDeltaEtaPhiMCConversion->Fill( dEta , dPhi      );
              fhConvPtMCConversion         ->Fill( pairM, calo->Pt()+calo2->Pt());
              fhConvDispersionMCConversion ->Fill( cluster1->GetDispersion(), cluster2->GetDispersion());
              fhConvM02MCConversion        ->Fill( cluster1->GetM02(), cluster2->GetM02());
              fhConvDistMCConversion       ->Fill( convDist , prodVertex.Mag() );
              fhConvDistMCConversion       ->Fill( convDist2, prodVertex.Mag() );
              
              if(dEta<0.05 && pairM<0.01 && asymmetry<0.1){
                fhConvDistMCConversionCuts->Fill( convDist , prodVertex.Mag() );
                fhConvDistMCConversionCuts->Fill( convDist2, prodVertex.Mag() );
              }
              
            }              
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCAntiNeutron)){
            if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiNeutron) && ancPDG==-2112 && ancLabel > -1){
              fhConvDeltaEtaMCAntiNeutron    ->Fill( pairM, dEta      );
              fhConvDeltaPhiMCAntiNeutron    ->Fill( pairM, dPhi      );
              fhConvAsymMCAntiNeutron        ->Fill( pairM, asymmetry );
              fhConvDeltaEtaPhiMCAntiNeutron ->Fill( dEta , dPhi      );
              fhConvPtMCAntiNeutron          ->Fill( pairM, calo->Pt()+calo2->Pt());
              fhConvDispersionMCAntiNeutron  ->Fill( cluster1->GetDispersion(), cluster2->GetDispersion());
              fhConvM02MCAntiNeutron         ->Fill( cluster1->GetM02(), cluster2->GetM02());
            }
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCAntiProton)){
            if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiProton) && ancPDG==-2212 && ancLabel > -1){
              fhConvDeltaEtaMCAntiProton    ->Fill( pairM, dEta      );
              fhConvDeltaPhiMCAntiProton    ->Fill( pairM, dPhi      );
              fhConvAsymMCAntiProton        ->Fill( pairM, asymmetry );
              fhConvDeltaEtaPhiMCAntiProton ->Fill( dEta , dPhi      );
              fhConvPtMCAntiProton          ->Fill( pairM, calo->Pt()+calo2->Pt());
              fhConvDispersionMCAntiProton  ->Fill( cluster1->GetDispersion(), cluster2->GetDispersion());
              fhConvM02MCAntiProton         ->Fill( cluster1->GetM02(), cluster2->GetM02());
            }
          }
          
          //Pairs coming from fragmenting pairs.
          if(ancPDG < 22 && ancLabel > 7 && (ancStatus == 11 || ancStatus == 12) ){
            fhConvDeltaEtaMCString    ->Fill( pairM, dPhi);
            fhConvDeltaPhiMCString    ->Fill( pairM, dPhi);
            fhConvAsymMCString        ->Fill( pairM, TMath::Abs(calo->E()-calo2->E())/(calo->E()+calo2->E()) );
            fhConvDeltaEtaPhiMCString ->Fill( dPhi,  dPhi);
            fhConvPtMCString          ->Fill( pairM, calo->Pt()+calo2->Pt());
            fhConvDispersionMCString  ->Fill( cluster1->GetDispersion(), cluster2->GetDispersion());
            fhConvM02MCString         ->Fill( cluster1->GetM02(), cluster2->GetM02());
          }
          
        }// Data MC
        
        break;
      }
      
    }//Mass loop
    
    //..........................................................................................................
    //Pair selected as converted, remove both clusters or recombine them into a photon and put them in the AOD
    if(bConverted){ 
      //Add to AOD
      if(fAddConvertedPairsToAOD){
        //Create AOD of pair analysis
        TLorentzVector mpair = *(calo->Momentum())+*(calo2->Momentum());
        AliAODPWG4Particle aodpair = AliAODPWG4Particle(mpair);
        aodpair.SetLabel(calo->GetLabel());
        
        //printf("Index %d, Id %d\n",iaod, calo->GetID());
        //Set the indeces of the original caloclusters  
        aodpair.SetCaloLabel(calo->GetCaloLabel(0),id2);
        aodpair.SetDetector(calo->GetDetector());
        aodpair.SetIdentifiedParticleType(calo->GetIdentifiedParticleType());
        aodpair.SetTag(calo ->GetTag());
        aodpair.SetTagged(kTRUE);
        //Add AOD with pair object to aod branch
        AddAODParticle(aodpair);
        //printf("\t \t both added pair\n");
      }
      
      //Do not add the current calocluster
      if(!fRemoveConvertedPair) 
      {
        //printf("TAGGED\n");
        //Tag this cluster as likely conversion
        calo->SetTagged(kTRUE);
      }
    }//converted pair
    
  }// main loop
  
  // Remove entries identified as conversion electrons
  // Revise if this is OK
  if(fRemoveConvertedPair || fAddConvertedPairsToAOD){
    for(Int_t iaod = 0; iaod < naod ; iaod++)
      if(indexConverted[iaod])GetOutputAODBranch()->RemoveAt(iaod);
    GetOutputAODBranch()->Compress();
  }
  
  delete [] indexConverted;
	
  if(GetDebug() > 1) printf("AliAnaPhotonConvInCalo::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//____________________________________________________
void  AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() 
{
  //Fill histograms
  
  //-------------------------------------------------------------------
  // Access MC information in stack if requested, check that it exists.	
  AliStack         * stack       = 0x0;
  TParticle        * primary     = 0x0;   
  TClonesArray     * mcparticles = 0x0;
  AliAODMCParticle * aodprimary  = 0x0; 
  
  if(IsDataMC()){
    
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;
      if(!stack) {
        printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called? STOP\n");
        abort();
      }
      
    }
    else if(GetReader()->ReadAODMCParticles()){
      
      //Get the list of MC particles
      mcparticles = GetReader()->GetAODMCParticles(0);
      if(!mcparticles && GetDebug() > 0) 	{
        printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
      }	
    }
  }// is data and MC
 
  //----------------------------------
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    //Int_t pdg = ph->GetIdentifiedParticleType();
    
    if(ph->IsTagged()){
      
      if(GetDebug() > 2) 
        printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() - ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
      //................................
      //Fill photon histograms 
      Float_t ptcluster  = ph->Pt();
      Float_t phicluster = ph->Phi();
      Float_t etacluster = ph->Eta();
      Float_t ecluster   = ph->E();
      
      fhPtPhotonConv->Fill(ptcluster);
      if(ecluster > 0.5)        fhEtaPhiPhotonConv  ->Fill(etacluster, phicluster);
      else if(GetMinPt() < 0.5) fhEtaPhi05PhotonConv->Fill(etacluster, phicluster);
      
      
      //.......................................
      //Play with the MC data if available
      if(IsDataMC()){
        
        
        //....................................................................
        // Access MC information in stack if requested, check that it exists.
        Int_t label =ph->GetLabel();
        if(label < 0) {
          if(GetDebug() > 1) printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
          continue;
        }
        
        Float_t eprim   = 0;
        Float_t ptprim  = 0;
        if(GetReader()->ReadStack()){
          
          if(label >=  stack->GetNtrack()) {
            if(GetDebug() > 2)  printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
            continue ;
          }
          
          primary = stack->Particle(label);
          if(!primary){
            printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
            continue;
          }
          eprim   = primary->Energy();
          ptprim  = primary->Pt();		
          
        }
        else if(GetReader()->ReadAODMCParticles()){
          //Check which is the input
          if(ph->GetInputFileIndex() == 0){
            if(!mcparticles) continue;
            if(label >=  mcparticles->GetEntriesFast()) {
              if(GetDebug() > 2)  printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
                                         label, mcparticles->GetEntriesFast());
              continue ;
            }
            //Get the particle
            aodprimary = (AliAODMCParticle*) mcparticles->At(label);
            
          }
          
          if(!aodprimary){
            printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
            continue;
          }
          
          eprim   = aodprimary->E();
          ptprim  = aodprimary->Pt();
          
        }
        
        Int_t tag =ph->GetTag();
        
        if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
        {
          
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
          {
            
            fhPtConversionTagged ->Fill(ptcluster);
            
          }			
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron))
        {
          
          fhPtAntiNeutronTagged ->Fill(ptcluster);
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton))
        {
          fhPtAntiProtonTagged ->Fill(ptcluster);
          
        } 
        
        else {
          fhPtUnknownTagged ->Fill(ptcluster);
          
        }
        
      }//Histograms with MC
    }// tagged by conversion
  }// aod loop
  
}


//________________________________________________________
void AliAnaPhotonConvInCalo::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  
  printf("Add conversion pair to AOD           = %d\n",fAddConvertedPairsToAOD);
  printf("Conversion pair mass cut             = %f\n",fMassCut);
  printf("Conversion selection cut : A < %1.2f; %1.3f < Dphi < %1.3f; Deta < %1.3f\n",
         fConvAsymCut,fConvDPhiMinCut, fConvDPhiMaxCut, fConvDEtaCut);
  
  printf("    \n") ;
	
} 
