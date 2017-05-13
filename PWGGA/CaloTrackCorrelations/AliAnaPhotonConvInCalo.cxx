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

// --- ROOT system --- 
#include <TH2F.h>
#include <TH3D.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include "TDatabasePDG.h"

// --- Analysis system --- 
#include "AliAnaPhotonConvInCalo.h" 
#include "AliCaloTrackReader.h"
#include "AliMCEvent.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVParticle.h"

/// \cond CLASSIMP
ClassImp(AliAnaPhotonConvInCalo) ;
/// \endcond

//________________________________________________
/// Default constructor. Initialize parameters.
//________________________________________________
AliAnaPhotonConvInCalo::AliAnaPhotonConvInCalo() :
AliAnaCaloTrackCorrBaseClass(),   
fRemoveConvertedPair(kFALSE), 
fAddConvertedPairsToAOD(kFALSE), 
fFillClusterConvDistHisto(kFALSE),
fMassCut(0),                       fMassCutTight(0),                  
fConvAsymCut(1.),                  fConvDEtaCut(2.),
fConvDPhiMinCut(-1.),              fConvDPhiMaxCut(7.), 
fMomentum(),                       fProdVertex(),
// Histograms
fhPtPhotonConv(0),                 
fhConvDeltaEta(0),                 fhConvDeltaPhi(0),              fhConvDeltaEtaPhi(0), 
fhConvAsym(0),                     fhConvPt(0),
fhConvDistEta(0),                  fhConvDistPhi(0),               fhConvDistEn(0),                fhConvDistMass(0),     
fhConvDistEtaCutEta(0),            fhConvDistPhiCutEta(0),         fhConvDistEnCutEta(0),          fhConvDistMassCutEta(0),
fhConvDistEtaCutMass(0),           fhConvDistPhiCutMass(0),        fhConvDistEnCutMass(0), 
fhConvDistEtaCutAsy(0),            fhConvDistPhiCutAsy(0),         fhConvDistEnCutAsy(0),          fhConvDistMassCutAsy(0),
fhConvDistEtaCutAll(0),            fhConvDistPhiCutAll(0),         fhConvDistEnCutAll(0),

// MC histograms
fhPtConversionTagged(0),           fhPtAntiNeutronTagged(0),       
fhPtAntiProtonTagged(0),           fhPtUnknownTagged(0),

fhConvDeltaEtaMCConversion(0),     fhConvDeltaPhiMCConversion(0),  fhConvDeltaEtaPhiMCConversion(0),
fhConvAsymMCConversion(0),         fhConvPtMCConversion(0),           
//fhConvDispersionMCConversion(0),   
fhConvM02MCConversion(0),

fhConvDeltaEtaMCAntiNeutron(0),    fhConvDeltaPhiMCAntiNeutron(0), fhConvDeltaEtaPhiMCAntiNeutron(0), 
fhConvAsymMCAntiNeutron(0),        fhConvPtMCAntiNeutron(0), 
//fhConvDispersionMCAntiNeutron(0),  
fhConvM02MCAntiNeutron(0),
fhConvDeltaEtaMCAntiProton(0),     fhConvDeltaPhiMCAntiProton(0),  fhConvDeltaEtaPhiMCAntiProton(0),  
fhConvAsymMCAntiProton(0),         fhConvPtMCAntiProton(0),  
//fhConvDispersionMCAntiProton(0),   
fhConvM02MCAntiProton(0),
fhConvDeltaEtaMCString(0),         fhConvDeltaPhiMCString(0),      fhConvDeltaEtaPhiMCString(0),      
fhConvAsymMCString(0),             fhConvPtMCString(0),      
//fhConvDispersionMCString(0),       
fhConvM02MCString(0),
fhConvDistMCConversion(0),         fhConvDistMCConversionCuts(0)
{
  InitParameters();
  
  for(Int_t ibin = 0; ibin < 6; ibin++)
  {
    fhEtaPhiPhotonConv      [ibin]  = 0;
    fhEtaPhiPhotonConvPaired[ibin]  = 0;
    fhConvPtMCConversionRcut[ibin]  = 0;
  //fhConvPtRcut            [ibin]  = 0;
  }
  
}

//_____________________________________________________
/// Save parameters used for analysis.
//_____________________________________________________
TObjString *  AliAnaPhotonConvInCalo::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPhotonConvInCalo---:") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Conversion pair mass cut = %1.2f and tight mass cut = %1.2f; ",fMassCut,fMassCutTight);
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Conversion Selection: fConvAsymCut %1.2f, fConvDEtaCut %1.2f fConvDPhiCut (%1.2f,%1.2f); ",
           fConvAsymCut, fConvDEtaCut, fConvDPhiMinCut, fConvDPhiMaxCut) ;
  parList+=onePar ; 
  snprintf(onePar,buffersize,"Add conversion pair to AOD = %d, remove partners? %d. ",fAddConvertedPairsToAOD, fRemoveConvertedPair);
  parList+=onePar ;	

  return new TObjString(parList) ;
}

//_______________________________________________________
/// Create histograms to be saved in output file and
/// store them in outputContainer.
//_______________________________________________________
TList *  AliAnaPhotonConvInCalo::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ;
  outputContainer->SetName("PhotonConvInCaloHistos") ;
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins(); Float_t phimax = GetHistogramRanges()->GetHistoPhiMax(); Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins(); Float_t etamax = GetHistogramRanges()->GetHistoEtaMax(); Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  
  Int_t ndist = 250;
  Int_t mindist = 0;
  Int_t maxdist = 500;
  
  TString region[] = {"ITS","TPC","TRD","TOF","Top EMCal","In EMCal"};

  fhPtPhotonConv  = new TH1F("hPtPhotonConv","Number of #gamma over calorimeter, conversion",nptbins,ptmin,ptmax);
  fhPtPhotonConv->SetYTitle("N");
  fhPtPhotonConv->SetXTitle("#it{p}_{T #gamma}(GeV/#it{c})");
  outputContainer->Add(fhPtPhotonConv) ;
  
  Float_t ebin[] = {0.3,0.4,0.5,0.6,0.75,1,2};
  
  for(Int_t iebin = 0; iebin < 6; iebin++)
  {
    fhEtaPhiPhotonConv[iebin]  = new TH2F
    (Form("hEtaPhiPhotonConv_ebin%d",iebin),Form("pair #eta vs #phi, %2.2f<#it{E}<%2.2f GeV/#it{c}",ebin[iebin],ebin[iebin+1]),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiPhotonConv[iebin]->SetYTitle("#phi (rad)");
    fhEtaPhiPhotonConv[iebin]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiPhotonConv[iebin]) ;

    fhEtaPhiPhotonConvPaired[iebin]  = new TH2F
    (Form("hEtaPhiPhotonConvPaired_ebin%d",iebin),Form("cluster #eta vs #phi, %2.2f<#it{E}<%2.2f GeV/#it{c}",ebin[iebin],ebin[iebin+1]),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiPhotonConvPaired[iebin]->SetYTitle("#phi (rad)");
    fhEtaPhiPhotonConvPaired[iebin]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiPhotonConvPaired[iebin]) ;
  }
  
  fhConvDeltaEta  = new TH2F
  ("hConvDeltaEta","#Delta #eta of selected conversion pairs",100,0,fMassCut,netabins*2,-0.5,0.5);
  fhConvDeltaEta->SetYTitle("#Delta #eta");
  fhConvDeltaEta->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
  outputContainer->Add(fhConvDeltaEta) ;
  
  fhConvDeltaPhi  = new TH2F
  ("hConvDeltaPhi","#Delta #phi of selected conversion pairs",100,0,fMassCut,nphibins*2,-0.5,0.5);
  fhConvDeltaPhi->SetYTitle("#Delta #phi");
  fhConvDeltaPhi->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
  outputContainer->Add(fhConvDeltaPhi) ;
  
  fhConvDeltaEtaPhi  = new TH2F
  ("hConvDeltaEtaPhi","#Delta #eta vs #Delta #phi of selected conversion pairs",netabins,-0.5,0.5,nphibins,-0.5,0.5);
  fhConvDeltaEtaPhi->SetYTitle("#Delta #phi");
  fhConvDeltaEtaPhi->SetXTitle("#Delta #eta");
  outputContainer->Add(fhConvDeltaEtaPhi) ;
  
  fhConvAsym  = new TH2F
  ("hConvAsym","Asymmetry of selected conversion pairs",100,0,fMassCut,100,0,1);
  fhConvAsym->SetYTitle("Asymmetry");
  fhConvAsym->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
  outputContainer->Add(fhConvAsym) ;
  
  fhConvPt  = new TH2F
  ("hConvPt","#it{p}_{T} of selected conversion pairs",100,0,fMassCut,100,0.,10.);
  fhConvPt->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
  fhConvPt->SetXTitle("Pair #it{M} (GeV/#it{c}^{2})");
  outputContainer->Add(fhConvPt) ;

//  for(Int_t iR = 0; iR < 6; iR++)
//  {
//    fhConvPtRcut[iR]  = new TH2F
//    (Form("hConvPt_R%d",iR),
//     Form("#it{p}_{T} of selected conversion pairs in %s?",region[iR].Data()),
//     100,0,fMassCut,100,0.,10.);
//    fhConvPtRcut[iR]->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
//    fhConvPtRcut[iR]->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
//    outputContainer->Add(fhConvPtRcut[iR]) ;
//  }    

  if(fFillClusterConvDistHisto)
  {
    fhConvDistEta  = new TH2F
    ("hConvDistEta","distance to conversion vertex",netabins,etamin,etamax,ndist,mindist,maxdist);
    fhConvDistEta->SetXTitle("#eta");
    fhConvDistEta->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEta) ;
    
    fhConvDistPhi  = new TH2F
    ("hConvDistPhi","distance to conversion vertex",nphibins,phimin,phimax,ndist,mindist,maxdist);
    fhConvDistPhi->SetXTitle("#phi (rad)");
    fhConvDistPhi->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistPhi) ;
    
    fhConvDistEn  = new TH2F
    ("hConvDistEn","distance to conversion vertex",nptbins,ptmin,ptmax,ndist,mindist,maxdist);
    fhConvDistEn->SetXTitle("#it{E} (GeV)");
    fhConvDistEn->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEn) ;
    
    fhConvDistMass  = new TH2F
    ("hConvDistMass","distance to conversion vertex",100,0,fMassCut,ndist,mindist,maxdist);
    fhConvDistMass->SetXTitle("#it{M} (GeV/#it{c}^{2})");
    fhConvDistMass->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistMass) ;
    
    fhConvDistEtaCutEta  = new TH2F
    ("hConvDistEtaCutEta",Form("distance to conversion vertex, #Delta #eta < %2.2f",fConvDEtaCut),netabins,etamin,etamax,ndist,mindist,maxdist);
    fhConvDistEtaCutEta->SetXTitle("#eta");
    fhConvDistEtaCutEta->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEtaCutEta) ;
    
    fhConvDistPhiCutEta  = new TH2F
    ("hConvDistPhiCutEta",Form("distance to conversion vertex, #Delta #eta < %2.2f",fConvDEtaCut),nphibins,phimin,phimax,ndist,mindist,maxdist);
    fhConvDistPhiCutEta->SetXTitle("#phi (rad)");
    fhConvDistPhiCutEta->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistPhiCutEta) ;
    
    fhConvDistEnCutEta  = new TH2F
    ("hConvDistEnCutEta",Form("distance to conversion vertex, #Delta #eta < %2.2f",fConvDEtaCut),nptbins,ptmin,ptmax,ndist,mindist,maxdist);
    fhConvDistEnCutEta->SetXTitle("#it{E} (GeV)");
    fhConvDistEnCutEta->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEnCutEta) ;
    
    fhConvDistMassCutEta  = new TH2F
    ("hConvDistMassCutEta",Form("distance to conversion vertex, #Delta #eta < %2.2f",fConvDEtaCut),100,0,fMassCut,ndist,mindist,maxdist);
    fhConvDistMassCutEta->SetXTitle("#it{M} (GeV/#it{c}^{2})");
    fhConvDistMassCutEta->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistMassCutEta) ;
    
    fhConvDistEtaCutMass  = new TH2F
    ("hConvDistEtaCutMass",Form("distance to conversion vertex, #it{M} < %0.2f MeV/#it{c}^{2}",fMassCutTight),netabins,etamin,etamax,ndist,mindist,maxdist);
    fhConvDistEtaCutMass->SetXTitle("#eta");
    fhConvDistEtaCutMass->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEtaCutMass) ;
    
    fhConvDistPhiCutMass  = new TH2F
    ("hConvDistPhiCutMass",Form("distance to conversion vertex,  #it{M} < %0.2f MeV/#it{c}^{2}",fMassCutTight),nphibins,phimin,phimax,ndist,mindist,maxdist);
    fhConvDistPhiCutMass->SetXTitle("#phi (rad)");
    fhConvDistPhiCutMass->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistPhiCutMass) ;
    
    fhConvDistEnCutMass  = new TH2F
    ("hConvDistEnCutMass",Form("distance to conversion vertex, #it{M} < %0.2f MeV/#it{c}^{2}",fMassCutTight),nptbins,ptmin,ptmax,ndist,mindist,maxdist);
    fhConvDistEnCutMass->SetXTitle("#it{E} (GeV)");
    fhConvDistEnCutMass->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEnCutMass) ;
    
    fhConvDistEtaCutAsy  = new TH2F
    ("hConvDistEtaCutAsy",Form("distance to conversion vertex, #it{A} < %2.2f",fConvAsymCut),netabins,etamin,etamax,ndist,mindist,maxdist);
    fhConvDistEtaCutAsy->SetXTitle("#eta");
    fhConvDistEtaCutAsy->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEtaCutAsy) ;
    
    fhConvDistPhiCutAsy  = new TH2F
    ("hConvDistPhiCutAsy",Form("distance to conversion vertex,  #it{A} < %2.2f",fConvAsymCut),nphibins,phimin,phimax,ndist,mindist,maxdist);
    fhConvDistPhiCutAsy->SetXTitle("#phi (rad)");
    fhConvDistPhiCutAsy->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistPhiCutAsy) ;
    
    fhConvDistEnCutAsy  = new TH2F
    ("hConvDistEnCutAsy",Form("distance to conversion vertex, #it{A} < %2.2f",fConvAsymCut),nptbins,ptmin,ptmax,ndist,mindist,maxdist);
    fhConvDistEnCutAsy->SetXTitle("#it{E} (GeV)");
    fhConvDistEnCutAsy->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEnCutAsy) ;
    
    fhConvDistMassCutAsy  = new TH2F
    ("hConvDistMassCutAsy",Form("distance to conversion vertex, #it{A} < %2.2f",fConvDEtaCut),100,0,fMassCut,ndist,mindist,maxdist);
    fhConvDistMassCutAsy->SetXTitle("#it{M} (GeV/#it{c}^{2})");
    fhConvDistMassCutEta->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistMassCutAsy) ;
    
    fhConvDistEtaCutAll  = new TH2F
    ("hConvDistEtaCutAll",Form("distance to conversion vertex, #it{M} < %0.2f MeV/#it{c}^{2}, #Delta #eta < %2.2f, #it{A} < %2.2f",
                               fConvDEtaCut,fMassCutTight,fConvAsymCut),netabins,etamin,etamax,ndist,mindist,maxdist);
    fhConvDistEtaCutAll->SetXTitle("#eta");
    fhConvDistEtaCutAll->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEtaCutAll) ;
    
    fhConvDistPhiCutAll  = new TH2F
    ("hConvDistPhiCutAll",Form("distance to conversion vertex,  #it{M} < %0.2f MeV/#it{c}^{2}, #Delta #eta < %2.2f, #it{A} < %2.2f",
                               fConvDEtaCut,fMassCutTight,fConvAsymCut),nphibins,phimin,phimax,ndist,mindist,maxdist);
    fhConvDistPhiCutAll->SetXTitle("#phi (rad)");
    fhConvDistPhiCutAll->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistPhiCutAll) ;
    
    fhConvDistEnCutAll  = new TH2F
    ("hConvDistEnCutAll",Form("distance to conversion vertex, #it{M} < %0.2f MeV/#it{c}^{2}, #Delta #eta < %2.2f, #it{A} < %2.2f",
                              fConvDEtaCut,fMassCutTight,fConvAsymCut),nptbins,ptmin,ptmax,ndist,mindist,maxdist);
    fhConvDistEnCutAll->SetXTitle("#it{E} (GeV)");
    fhConvDistEnCutAll->SetYTitle(" distance (cm)");
    outputContainer->Add(fhConvDistEnCutAll) ;
  }
  
  if(IsDataMC())
  {
    fhPtConversionTagged  = new TH1F("hPtMCConversionTagged","Number of converted #gamma over calorimeter, tagged as converted",nptbins,ptmin,ptmax);
    fhPtConversionTagged->SetYTitle("N");
    fhPtConversionTagged->SetXTitle("#it{p}_{T #gamma}(GeV/#it{c})");
    outputContainer->Add(fhPtConversionTagged) ;
    
    fhPtAntiNeutronTagged  = new TH1F("hPtMCAntiNeutronTagged","Number of AntiNeutron id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax);
    fhPtAntiNeutronTagged->SetYTitle("N");
    fhPtAntiNeutronTagged->SetXTitle("#it{p}_{T #gamma}(GeV/#it{c})");
    outputContainer->Add(fhPtAntiNeutronTagged) ;
    
    fhPtAntiProtonTagged  = new TH1F("hPtMCAntiProtonTagged","Number of AntiProton id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax);
    fhPtAntiProtonTagged->SetYTitle("N");
    fhPtAntiProtonTagged->SetXTitle("#it{p}_{T #gamma}(GeV/#it{c})");
    outputContainer->Add(fhPtAntiProtonTagged) ;
    
    fhPtUnknownTagged  = new TH1F("hPtMCUnknownTagged","Number of Unknown id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax);
    fhPtUnknownTagged->SetYTitle("N");
    fhPtUnknownTagged->SetXTitle("#it{p}_{T #gamma}(GeV/#it{c})");
    outputContainer->Add(fhPtUnknownTagged) ;
    
    if(fFillClusterConvDistHisto)
    {
      fhConvDistMCConversion  = new TH2F
      ("hConvDistMCConversion","calculated conversion distance vs real vertes for MC conversion",ndist,mindist,maxdist,ndist,mindist,maxdist);
      fhConvDistMCConversion->SetXTitle("distance");
      fhConvDistMCConversion->SetYTitle("vertex R");
      outputContainer->Add(fhConvDistMCConversion) ;
      
      fhConvDistMCConversionCuts  = new TH2F
      ("hConvDistMCConversionCuts",
       Form("calculated conversion distance vs real vertes for MC conversion, #Delta #eta < %2.2f, #it{M} < 10 MeV/#it{c}^{2}, #it{A} < %2.2f",
            fConvDEtaCut,fConvAsymCut),
       ndist,mindist,maxdist,ndist,mindist,maxdist);
      fhConvDistMCConversionCuts->SetXTitle("distance");
      fhConvDistMCConversionCuts->SetYTitle("vertex R");
      outputContainer->Add(fhConvDistMCConversionCuts) ;
    }
    
    fhConvDeltaEtaMCConversion  = new TH2F
    ("hConvDeltaEtaMCConversion","#Delta #eta of selected conversion pairs from real conversions",100,0,fMassCut,netabins,-0.5,0.5);
    fhConvDeltaEtaMCConversion->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCConversion->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaEtaMCConversion) ;
    
    fhConvDeltaPhiMCConversion  = new TH2F
    ("hConvDeltaPhiMCConversion","#Delta #phi of selected conversion pairs from real conversions",100,0,fMassCut,nphibins,-0.5,0.5);
    fhConvDeltaPhiMCConversion->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCConversion->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaPhiMCConversion) ;
    
    fhConvDeltaEtaPhiMCConversion  = new TH2F
    ("hConvDeltaEtaPhiMCConversion","#Delta #eta vs #Delta #phi of selected conversion pairs, from real conversions",netabins,-0.5,0.5,nphibins,-0.5,0.5);
    fhConvDeltaEtaPhiMCConversion->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCConversion->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCConversion) ;
    
    fhConvAsymMCConversion  = new TH2F
    ("hConvAsymMCConversion","Asymmetry of selected conversion pairs from real conversions",100,0,fMassCut,100,0,1);
    fhConvAsymMCConversion->SetYTitle("Asymmetry");
    fhConvAsymMCConversion->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvAsymMCConversion) ;
    
    fhConvPtMCConversion  = new TH2F
    ("hConvPtMCConversion","#it{p}_{T} of selected conversion pairs from real conversions",100,0,fMassCut,nptbins,ptmin,ptmax);
    fhConvPtMCConversion->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
    fhConvPtMCConversion->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvPtMCConversion) ;
    
    for(Int_t iR = 0; iR < 6; iR++)
    {
      fhConvPtMCConversionRcut[iR]  = new TH2F
      (Form("hConvPtMCConversion_R%d",iR),
       Form("#it{p}_{T} of selected conversion pairs from real conversions in %s",region[iR].Data()),
       100,0,fMassCut,nptbins,ptmin,ptmax);
      fhConvPtMCConversionRcut[iR]->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
      fhConvPtMCConversionRcut[iR]->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
      outputContainer->Add(fhConvPtMCConversionRcut[iR]) ;
    }    
    
//    fhConvDispersionMCConversion  = new TH2F
//    ("hConvDispersionMCConversion","#it{p}_{T} of selected conversion pairs from real conversions",100,0.,1.,100,0.,1.);
//    fhConvDispersionMCConversion->SetYTitle("Dispersion cluster 1");
//    fhConvDispersionMCConversion->SetXTitle("Dispersion cluster 2");
//    outputContainer->Add(fhConvDispersionMCConversion) ;
    
    fhConvM02MCConversion  = new TH2F
    ("hConvM02MCConversion","#it{p}_{T} of selected conversion pairs from real conversion",100,0.,1.,100,0.,1.);
    fhConvM02MCConversion->SetYTitle("M02 cluster 1");
    fhConvM02MCConversion->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCConversion) ;
    
    fhConvDeltaEtaMCAntiNeutron  = new TH2F
    ("hConvDeltaEtaMCAntiNeutron","#Delta #eta of selected conversion pairs from anti-neutrons",100,0,fMassCut,netabins,-0.5,0.5);
    fhConvDeltaEtaMCAntiNeutron->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCAntiNeutron->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaEtaMCAntiNeutron) ;
    
    fhConvDeltaPhiMCAntiNeutron  = new TH2F
    ("hConvDeltaPhiMCAntiNeutron","#Delta #phi of selected conversion pairs from anti-neutrons",100,0,fMassCut,nphibins,-0.5,0.5);
    fhConvDeltaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCAntiNeutron->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaPhiMCAntiNeutron) ;
    
    fhConvDeltaEtaPhiMCAntiNeutron  = new TH2F
    ("hConvDeltaEtaPhiMCAntiNeutron","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-neutrons",netabins,-0.5,0.5,nphibins,-0.5,0.5);
    fhConvDeltaEtaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCAntiNeutron->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCAntiNeutron) ;
    
    fhConvAsymMCAntiNeutron  = new TH2F
    ("hConvAsymMCAntiNeutron","Asymmetry of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0,1);
    fhConvAsymMCAntiNeutron->SetYTitle("Asymmetry");
    fhConvAsymMCAntiNeutron->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvAsymMCAntiNeutron) ;
    
    fhConvPtMCAntiNeutron  = new TH2F
    ("hConvPtMCAntiNeutron","#it{p}_{T} of selected conversion pairs from anti-neutrons",100,0,fMassCut,nptbins,ptmin,ptmax);
    fhConvPtMCAntiNeutron->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
    fhConvPtMCAntiNeutron->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvPtMCAntiNeutron) ;
    
//    fhConvDispersionMCAntiNeutron  = new TH2F
//    ("hConvDispersionMCAntiNeutron","#it{p}_{T} of selected conversion pairs from anti-neutrons",100,0.,1.,100,0.,1.);
//    fhConvDispersionMCAntiNeutron->SetYTitle("Dispersion cluster 1");
//    fhConvDispersionMCAntiNeutron->SetXTitle("Dispersion cluster 2");
//    outputContainer->Add(fhConvDispersionMCAntiNeutron) ;
    
    fhConvM02MCAntiNeutron  = new TH2F
    ("hConvM02MCAntiNeutron","#it{p}_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.);
    fhConvM02MCAntiNeutron->SetYTitle("M02 cluster 1");
    fhConvM02MCAntiNeutron->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCAntiNeutron) ;
    
    fhConvDeltaEtaMCAntiProton  = new TH2F
    ("hConvDeltaEtaMCAntiProton","#Delta #eta of selected conversion pairs from anti-protons",100,0,fMassCut,netabins,-0.5,0.5);
    fhConvDeltaEtaMCAntiProton->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCAntiProton->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaEtaMCAntiProton) ;
    
    fhConvDeltaPhiMCAntiProton  = new TH2F
    ("hConvDeltaPhiMCAntiProton","#Delta #phi of selected conversion pairs from anti-protons",100,0,fMassCut,nphibins,-0.5,0.5);
    fhConvDeltaPhiMCAntiProton->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCAntiProton->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaPhiMCAntiProton) ;
    
    fhConvDeltaEtaPhiMCAntiProton  = new TH2F
    ("hConvDeltaEtaPhiMCAntiProton","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-protons",netabins,-0.5,0.5,nphibins,-0.5,0.5);
    fhConvDeltaEtaPhiMCAntiProton->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCAntiProton->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCAntiProton) ;
    
    fhConvAsymMCAntiProton  = new TH2F
    ("hConvAsymMCAntiProton","Asymmetry of selected conversion pairs from anti-protons",100,0,fMassCut,100,0,1);
    fhConvAsymMCAntiProton->SetYTitle("Asymmetry");
    fhConvAsymMCAntiProton->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvAsymMCAntiProton) ;
    
    fhConvPtMCAntiProton  = new TH2F
    ("hConvPtMCAntiProton","#it{p}_{T} of selected conversion pairs from anti-protons",100,0,fMassCut,nptbins,ptmin,ptmax);
    fhConvPtMCAntiProton->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
    fhConvPtMCAntiProton->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvPtMCAntiProton) ;
    
//    fhConvDispersionMCAntiProton  = new TH2F
//    ("hConvDispersionMCAntiProton","#it{p}_{T} of selected conversion pairs from anti-protons",100,0.,1.,100,0.,1.);
//    fhConvDispersionMCAntiProton->SetYTitle("Dispersion cluster 1");
//    fhConvDispersionMCAntiProton->SetXTitle("Dispersion cluster 2");
//    outputContainer->Add(fhConvDispersionMCAntiProton) ;
    
    fhConvM02MCAntiProton  = new TH2F
    ("hConvM02MCAntiProton","#it{p}_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.);
    fhConvM02MCAntiProton->SetYTitle("M02 cluster 1");
    fhConvM02MCAntiProton->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCAntiProton) ;
    
    fhConvDeltaEtaMCString  = new TH2F
    ("hConvDeltaEtaMCString","#Delta #eta of selected conversion pairs from string",100,0,fMassCut,netabins,-0.5,0.5);
    fhConvDeltaEtaMCString->SetYTitle("#Delta #eta");
    fhConvDeltaEtaMCString->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaEtaMCString) ;
    
    fhConvDeltaPhiMCString  = new TH2F
    ("hConvDeltaPhiMCString","#Delta #phi of selected conversion pairs from string",100,0,fMassCut,nphibins,-0.5,0.5);
    fhConvDeltaPhiMCString->SetYTitle("#Delta #phi");
    fhConvDeltaPhiMCString->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvDeltaPhiMCString) ;
    
    fhConvDeltaEtaPhiMCString  = new TH2F
    ("hConvDeltaEtaPhiMCString","#Delta #eta vs #Delta #phi of selected conversion pairs from string",netabins,-0.5,0.5,nphibins,-0.5,0.5);
    fhConvDeltaEtaPhiMCString->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhiMCString->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhiMCString) ;
    
    fhConvAsymMCString  = new TH2F
    ("hConvAsymMCString","Asymmetry of selected conversion pairs from string",100,0,fMassCut,100,0,1);
    fhConvAsymMCString->SetYTitle("Asymmetry");
    fhConvAsymMCString->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvAsymMCString) ;
    
    fhConvPtMCString  = new TH2F
    ("hConvPtMCString","#it{p}_{T} of selected conversion pairs from string",100,0,fMassCut,nptbins,ptmin,ptmax);
    fhConvPtMCString->SetYTitle("Pair #it{p}_{T} (GeV/#it{c})");
    fhConvPtMCString->SetXTitle("Pair Mass (GeV/#it{c}^{2})");
    outputContainer->Add(fhConvPtMCString) ;
    
//    fhConvDispersionMCString  = new TH2F
//    ("hConvDispersionMCString","#it{p}_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.);
//    fhConvDispersionMCString->SetYTitle("Dispersion cluster 1");
//    fhConvDispersionMCString->SetXTitle("Dispersion cluster 2");
//    outputContainer->Add(fhConvDispersionMCString) ;
    
    fhConvM02MCString  = new TH2F
    ("hConvM02MCString","#it{p}_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.);
    fhConvM02MCString->SetYTitle("M02 cluster 1");
    fhConvM02MCString->SetXTitle("M02 cluster 2");
    outputContainer->Add(fhConvM02MCString) ;    
  }
  
  return outputContainer ;
}

//___________________________________________
/// Initialize the parameters of the analysis.
//___________________________________________
void AliAnaPhotonConvInCalo::InitParameters()
{
  AddToHistogramsName("AnaPhotonConvInCalo_");
  
  fConvAsymCut            = 0.2   ;
  fMassCut                = 0.04  ; // 40 MeV/c^2
  fMassCutTight           = 0.03  ; // 30 MeV/c^2
  fConvDEtaCut            = 0.04  ;
  fConvDPhiMinCut         = 0     ;            
  fConvDPhiMaxCut         = 0.05  ;
  fRemoveConvertedPair    = kFALSE;
  fAddConvertedPairsToAOD = kFALSE;
}

//________________________________________________________
/// Fill histograms for selected pairs.
/// Fill aod with selected pair, remove other pairs
//________________________________________________________
void  AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms()
{
  // Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  AliDebug(1,Form("AOD branch entries %d", naod));
  
  // List to be used in conversion analysis, to tag the cluster as candidate for conversion
  Bool_t * indexConverted = new Bool_t[naod];
  for (Int_t i = 0; i < naod; i++) indexConverted[i] = kFALSE;
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* calo =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    
    //Bool_t bConverted = kFALSE;
    
    // Check if set previously as converted couple, if so skip its use.
    //if (indexConverted[iaod]) continue;
        
    // Second cluster loop
    AliAODPWG4Particle* calo2 = 0;
    for(Int_t jaod = iaod + 1 ; jaod < naod ; jaod++)
    {
      // Check if set previously as converted couple, if so skip its use.
      //if (indexConverted[jaod]) continue;
      
      //printf("Check Conversion indeces %d and %d\n",iaod,jaod);
      calo2 =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(jaod));
      
      //
      // Add both clusters
      //
      fMomentum = *(calo->Momentum())+*(calo2->Momentum());
      Float_t ptConvPair  = fMomentum.Pt();
      Float_t phiConvPair = fMomentum.Phi();
      Float_t etaConvPair = fMomentum.Eta();
    //Float_t eConvPair   = fMomentum.E();
      Float_t pairM       = calo->GetPairMass(calo2);
      //printf("\t both in calo, mass %f, cut %f\n",pairM,fMassCut);
      
      //................................................
      // Get mass of pair, if small, take this pair.
      if(pairM < fMassCut)
      {
        Float_t phi  = calo->Phi();
        if(phi < 0) phi +=TMath::TwoPi();
        Float_t phi2 = calo2->Phi();
        if(phi2 < 0) phi2+=TMath::TwoPi();
        
        calo->SetTagged(kFALSE);
        Float_t asymmetry = TMath::Abs(calo->E()-calo2->E())/(calo->E()+calo2->E());
        Float_t dPhi      = phi-phi2;
        Float_t dEta      = calo->Eta() - calo2->Eta();
        
        //...............................................
        // Fill few histograms with kinematics of the pair
        // FIXME, move all this to MakeAnalysisFillHistograms ...
        
        fhConvDeltaEta   ->Fill( pairM, dEta      , GetEventWeight());
        fhConvDeltaPhi   ->Fill( pairM, dPhi      , GetEventWeight());
        fhConvAsym       ->Fill( pairM, asymmetry , GetEventWeight());
        fhConvDeltaEtaPhi->Fill( dEta , dPhi      , GetEventWeight());
        fhConvPt         ->Fill( pairM, ptConvPair, GetEventWeight());
        
        Float_t convDist  = -1;
        Float_t convDist2 = -1;
        if(fFillClusterConvDistHisto)
        {
          // Estimate conversion distance, T. Awes, M. Ivanov
          // Under the assumption that the pair has zero mass, and that each electron
          // of the pair has the same momentum, they will each have the same bend radius
          // given by R=p/(qB) = p / (300 B) with p in [MeV/c], B in [Tesla] and R in [m].
          // With nominal ALICE magnet current of 30kA B=0.5T, and so with E_cluster=p,
          // R = E/1.5 [cm].  Under these assumptions, the distance from the conversion
          // point to the EMCal can be related to the separation distance, L=2y, on the EMCal
          // as d = sqrt(R^2 -(R-y)^2) = sqrt(2Ry - y^2). And since R>>y we can write as
          // d = sqrt(E*L/1.5) where E is the cluster energy and L is the distance in cm between
          // the clusters.
          
          TObjArray * clusters    = 0;
          if(calo->GetDetectorTag() == kEMCAL) clusters = GetEMCALClusters();
          else                                 clusters = GetPHOSClusters ();
          
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
          
          Float_t convDist  = TMath::Sqrt(calo ->E()*clustDist*0.01/0.15) * 100.; // cm
          Float_t convDist2 = TMath::Sqrt(calo2->E()*clustDist*0.01/0.15) * 100.; // cm
                                                                                  //printf("l = %f, e1 = %f, d1=%f, e2 = %f, d2=%f\n",clustDist,calo->E(),convDist,calo2->E(),convDist2);
          AliDebug(2,Form("Pair with mass %2.3f < %2.3f, %1.2f < dPhi %2.2f < %2.2f, dEta %f < %2.2f, asymmetry %2.2f< %2.2f; \n"
                          " cluster1 id %d, e %2.3f  SM %d, eta %2.3f, phi %2.3f ; \n"    
                          " cluster2 id %d, e %2.3f, SM %d, eta %2.3f, phi %2.3f   \n",
                          pairM,fMassCut,fConvDPhiMinCut, dPhi, fConvDPhiMaxCut, dEta, fConvDEtaCut, asymmetry, fConvAsymCut,
                          calo ->GetCaloLabel(0), calo ->E(), GetCaloUtils()->GetModuleNumber(calo ,GetReader()->GetInputEvent()), calo ->Eta(), phi,
                          calo2->GetCaloLabel(0), calo2->E(), GetCaloUtils()->GetModuleNumber(calo2,GetReader()->GetInputEvent()), calo2->Eta(), phi2));
          
          
//          Int_t convDistR1 = -1;
//          if      ( convDist  > 430-75.  ) convDistR1 = 0;
//          else if ( convDist  < 430-275. ) convDistR1 = 1;
//          else if ( convDist  < 430-375. ) convDistR1 = 2;
//          else if ( convDist  < 430-400. ) convDistR1 = 3;
//          else if ( convDist  < 0        ) convDistR1 = 4;
//          else                             convDistR1 = 5;
//          
//          Int_t convDistR2 = -1;
//          if      ( convDist2 > 430-75.  ) convDistR2 = 0;
//          else if ( convDist2 < 430-275. ) convDistR2 = 1;
//          else if ( convDist2 < 430-375. ) convDistR2 = 2;
//          else if ( convDist2 < 430-400. ) convDistR2 = 3;
//          else if ( convDist2 < 0        ) convDistR2 = 4;
//          else                             convDistR2 = 5;
//          
//          if(convDistR1 >= 0)
//            fhConvPtRcut[convDistR1]->Fill( pairM, ptConvPair, GetEventWeight());
//          if(convDistR2 >= 0)
//            fhConvPtRcut[convDistR2]->Fill( pairM, ptConvPair, GetEventWeight());
          
          fhConvDistEta ->Fill(calo ->Eta(),convDist , GetEventWeight());
          fhConvDistEta ->Fill(calo2->Eta(),convDist2, GetEventWeight());        
          fhConvDistPhi ->Fill(phi ,        convDist , GetEventWeight());
          fhConvDistPhi ->Fill(phi2,        convDist2, GetEventWeight());
          fhConvDistEn  ->Fill(calo ->E(),  convDist , GetEventWeight());
          fhConvDistEn  ->Fill(calo2->E(),  convDist2, GetEventWeight());
          fhConvDistMass->Fill(pairM,       convDist , GetEventWeight());
          fhConvDistMass->Fill(pairM,       convDist2, GetEventWeight());
          
          //dEta cut
          if(TMath::Abs(dEta) < fConvDEtaCut)
          {
            fhConvDistEtaCutEta ->Fill(calo ->Eta(), convDist , GetEventWeight());
            fhConvDistEtaCutEta ->Fill(calo2->Eta(), convDist2, GetEventWeight());
            fhConvDistPhiCutEta ->Fill(phi ,         convDist , GetEventWeight());
            fhConvDistPhiCutEta ->Fill(phi2,         convDist2, GetEventWeight());
            fhConvDistEnCutEta  ->Fill(calo ->E(),   convDist , GetEventWeight());
            fhConvDistEnCutEta  ->Fill(calo2->E(),   convDist2, GetEventWeight());
            fhConvDistMassCutEta->Fill(pairM,        convDist , GetEventWeight());
            fhConvDistMassCutEta->Fill(pairM,        convDist2, GetEventWeight());
          }
          
          //mass cut
          if(pairM < fMassCutTight)
          {
            fhConvDistEtaCutMass ->Fill(calo ->Eta(), convDist , GetEventWeight());
            fhConvDistEtaCutMass ->Fill(calo2->Eta(), convDist2, GetEventWeight());
            fhConvDistPhiCutMass ->Fill(phi ,         convDist , GetEventWeight());
            fhConvDistPhiCutMass ->Fill(phi2,         convDist2, GetEventWeight());
            fhConvDistEnCutMass  ->Fill(calo ->E(),   convDist , GetEventWeight());
            fhConvDistEnCutMass  ->Fill(calo2->E(),   convDist2, GetEventWeight());
          }
          
          // asymmetry cut
          if(asymmetry < fConvAsymCut)
          {
            fhConvDistEtaCutAsy ->Fill(calo ->Eta(), convDist , GetEventWeight());
            fhConvDistEtaCutAsy ->Fill(calo2->Eta(), convDist2, GetEventWeight());
            fhConvDistPhiCutAsy ->Fill(phi ,         convDist , GetEventWeight());
            fhConvDistPhiCutAsy ->Fill(phi2,         convDist2, GetEventWeight());
            fhConvDistEnCutAsy  ->Fill(calo ->E(),   convDist , GetEventWeight());
            fhConvDistEnCutAsy  ->Fill(calo2->E(),   convDist2, GetEventWeight());
            fhConvDistMassCutAsy->Fill(pairM,        convDist , GetEventWeight());
            fhConvDistMassCutAsy->Fill(pairM,        convDist2, GetEventWeight());
          }// asymmetry cut
          
          // All cuts
          if(TMath::Abs(dEta) < fConvDEtaCut && pairM < fMassCutTight && asymmetry < fConvAsymCut)
          {
            fhConvDistEtaCutAll ->Fill(calo ->Eta(), convDist , GetEventWeight());
            fhConvDistEtaCutAll ->Fill(calo2->Eta(), convDist2, GetEventWeight());
            fhConvDistPhiCutAll ->Fill(phi ,         convDist , GetEventWeight());
            fhConvDistPhiCutAll ->Fill(phi2,         convDist2, GetEventWeight());
            fhConvDistEnCutAll  ->Fill(calo ->E(),   convDist , GetEventWeight());
            fhConvDistEnCutAll  ->Fill(calo2->E(),   convDist2, GetEventWeight());
          }
        }
        
        //...........................................
        // Fill more histograms, simulated data
        Int_t ancPDG    = 0;
        Int_t ancStatus = 0;
        Int_t ancLabel  =-1;
        if(IsDataMC())
        {
          // Check the origin of the pair, look for conversion, antinucleons or jet correlations (strings)
          
          ancLabel  = GetMCAnalysisUtils()->CheckCommonAncestor(calo->GetLabel(), calo2->GetLabel(),
                                                                GetMC(), ancPDG, ancStatus, fMomentum, fProdVertex);
          
          // printf("AliAnaPhotonConvInCalo::MakeAnalysisFillHistograms() - Common ancestor label %d, pdg %d, name %s, status %d; \n",
          //                          ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
          
          Int_t   tag1  = calo ->GetTag();
          Int_t   tag2  = calo2->GetTag();
//          Float_t disp1 = cluster1->GetDispersion();
//          Float_t disp2 = cluster2->GetDispersion();
          Float_t l0cl1 = calo ->GetM02();
          Float_t l0cl2 = calo2->GetM02();
         
          if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCConversion))
          {
            if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) && (ancPDG==22 || TMath::Abs(ancPDG)==11) && ancLabel > -1)
            {
              fhConvDeltaEtaMCConversion   ->Fill( pairM, dEta      , GetEventWeight());
              fhConvDeltaPhiMCConversion   ->Fill( pairM, dPhi      , GetEventWeight());
              fhConvAsymMCConversion       ->Fill( pairM, asymmetry , GetEventWeight());
              fhConvDeltaEtaPhiMCConversion->Fill( dEta , dPhi      , GetEventWeight());
              fhConvPtMCConversion         ->Fill( pairM, ptConvPair, GetEventWeight());
            //fhConvDispersionMCConversion ->Fill( disp1, disp2     , GetEventWeight());
              fhConvM02MCConversion        ->Fill( l0cl1, l0cl2     , GetEventWeight());
              
              Float_t prodR = TMath::Sqrt(fProdVertex.X()*fProdVertex.X()+fProdVertex.Y()*fProdVertex.Y());

              if(fFillClusterConvDistHisto)
              {
                fhConvDistMCConversion       ->Fill( convDist , prodR, GetEventWeight());
                fhConvDistMCConversion       ->Fill( convDist2, prodR, GetEventWeight());
              }
              
              Int_t convR = -1;
              if      ( prodR < 75.  ) convR = 0;
              else if ( prodR < 275. ) convR = 1;
              else if ( prodR < 375. ) convR = 2;
              else if ( prodR < 400. ) convR = 3;
              else if ( prodR < 430. ) convR = 4;
              else                     convR = 5;

              if ( convR >= 0 )
                fhConvPtMCConversionRcut[convR]->Fill( pairM, ptConvPair, GetEventWeight());

              if ( fFillClusterConvDistHisto && dEta < fConvDEtaCut && pairM < fMassCutTight && asymmetry < fConvAsymCut )
              {
                fhConvDistMCConversionCuts->Fill( convDist , prodR, GetEventWeight());
                fhConvDistMCConversionCuts->Fill( convDist2, prodR, GetEventWeight());
              }
              
            }
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCAntiNeutron))
          {
            if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiNeutron) && ancPDG==-2112 && ancLabel > -1)
            {
              fhConvDeltaEtaMCAntiNeutron    ->Fill( pairM, dEta      , GetEventWeight());
              fhConvDeltaPhiMCAntiNeutron    ->Fill( pairM, dPhi      , GetEventWeight());
              fhConvAsymMCAntiNeutron        ->Fill( pairM, asymmetry , GetEventWeight());
              fhConvDeltaEtaPhiMCAntiNeutron ->Fill( dEta , dPhi      , GetEventWeight());
              fhConvPtMCAntiNeutron          ->Fill( pairM, ptConvPair, GetEventWeight());
            //fhConvDispersionMCAntiNeutron  ->Fill( disp1, disp2     , GetEventWeight());
              fhConvM02MCAntiNeutron         ->Fill( l0cl1, l0cl2     , GetEventWeight());
            }
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag1,AliMCAnalysisUtils::kMCAntiProton))
          {
            if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiProton) && ancPDG==-2212 && ancLabel > -1)
            {
              fhConvDeltaEtaMCAntiProton    ->Fill( pairM, dEta      , GetEventWeight());
              fhConvDeltaPhiMCAntiProton    ->Fill( pairM, dPhi      , GetEventWeight());
              fhConvAsymMCAntiProton        ->Fill( pairM, asymmetry , GetEventWeight());
              fhConvDeltaEtaPhiMCAntiProton ->Fill( dEta , dPhi      , GetEventWeight());
              fhConvPtMCAntiProton          ->Fill( pairM, ptConvPair, GetEventWeight());
            //fhConvDispersionMCAntiProton  ->Fill( disp1, disp2     , GetEventWeight());
              fhConvM02MCAntiProton         ->Fill( l0cl1, l0cl2     , GetEventWeight());
            }
          }
          
          // Pairs coming from fragmenting pairs.
          if( ancPDG < 22 && ancLabel > 7 && (ancStatus == 11 || ancStatus == 12) )
          {
            fhConvDeltaEtaMCString    ->Fill( pairM, dEta      , GetEventWeight());
            fhConvDeltaPhiMCString    ->Fill( pairM, dPhi      , GetEventWeight());
            fhConvAsymMCString        ->Fill( pairM, asymmetry , GetEventWeight());
            fhConvDeltaEtaPhiMCString ->Fill( dEta,  dPhi      , GetEventWeight());
            fhConvPtMCString          ->Fill( pairM, ptConvPair, GetEventWeight());
          //fhConvDispersionMCString  ->Fill( disp1, disp2     , GetEventWeight());
            fhConvM02MCString         ->Fill( l0cl1, l0cl2     , GetEventWeight());
          }
          
        }// Data MC
        
        //...............................................
        //...............................................
        // Select pairs in a eta/phi/asymmetry windows
        //
        if(TMath::Abs(dEta) < fConvDEtaCut    &&
           TMath::Abs(dPhi) < fConvDPhiMaxCut &&
           TMath::Abs(dPhi) > fConvDPhiMinCut &&
           asymmetry        < fConvAsymCut       )
        {
          indexConverted[iaod] = kTRUE;
          indexConverted[jaod] = kTRUE;
          //bConverted           = kTRUE;
          
          Float_t ebin[] = {0.3,0.4,0.5,0.6,0.75,1,2};
          Int_t bin1  = -1;
          Int_t bin2  = -1;
          Int_t bin12 = -1;
          for(Int_t iebin = 0; iebin < 6; iebin++)
          {
            if( calo ->Pt() > ebin[iebin] && calo ->Pt() <= ebin[iebin+1] ) bin1  = iebin ;
            if( calo2->Pt() > ebin[iebin] && calo2->Pt() <= ebin[iebin+1] ) bin2  = iebin ;
            if( ptConvPair  > ebin[iebin] && ptConvPair  <= ebin[iebin+1] ) bin12 = iebin ;
          }
          
          if(bin1  > -1) fhEtaPhiPhotonConvPaired[bin1 ]->Fill(calo ->Eta(), phi        , GetEventWeight());
          if(bin2  > -1) fhEtaPhiPhotonConvPaired[bin2 ]->Fill(calo2->Eta(), phi2       , GetEventWeight());
          if(bin12 > -1) fhEtaPhiPhotonConv      [bin12]->Fill(etaConvPair , phiConvPair, GetEventWeight());

          fhPtPhotonConv->Fill(ptConvPair, GetEventWeight());
          
          if(IsDataMC())
          {
            //....................................................................
            // Access MC information in stack if requested, check that it exists.
            
            Int_t label =calo->GetLabel();
            if ( label < 0 )
            {
              AliDebug(1,Form("*** bad label ***:  label %d", label));
              continue;
            }
            
            if ( label >=  GetMC()->GetNumberOfTracks() )
            {
              AliDebug(1,Form("*** large label ***:  label %d, n tracks %d", label, GetMC()->GetNumberOfTracks()));
              continue ;
            }
            
            //Float_t eprim   = 0;
            //Float_t ptprim  = 0;
                     
            // Get the particle
            AliVParticle* primary = GetMC()->GetTrack(label);
            
            if(!primary)
            {
              AliDebug(2,Form("*** no primary ***:  label %d", label));
              continue;
            }
            //eprim   = aodprimary->E();
            //ptprim  = aodprimary->Pt();
            
            Int_t tag = calo ->GetTag();
            if(ancLabel >=0 )
            {
              if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
              {
                if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
                {
                  fhPtConversionTagged ->Fill(ptConvPair, GetEventWeight());
                }
              }
              else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron))
              {
                fhPtAntiNeutronTagged ->Fill(ptConvPair, GetEventWeight());
              }
              else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton))
              {
                fhPtAntiProtonTagged ->Fill(ptConvPair, GetEventWeight());
              }
            }
            else
            {
              fhPtUnknownTagged ->Fill(ptConvPair, GetEventWeight());
            }
          } // MC
          
          //..........................................................................................................
          // Pair selected as converted, remove both clusters or recombine them into a photon and put them in the AOD
          // Add to AOD
          if(fAddConvertedPairsToAOD)
          {            
            // Create AOD of pair analysis
            AliAODPWG4Particle aodpair = AliAODPWG4Particle(fMomentum);
            aodpair.SetLabel(calo->GetLabel());
            
            // Set the indeces of the original caloclusters
            aodpair.SetCaloLabel(calo->GetCaloLabel(0),calo2->GetCaloLabel(0));
            aodpair.SetDetectorTag(calo->GetDetectorTag());
            aodpair.SetIdentifiedParticleType(calo->GetIdentifiedParticleType());
            aodpair.SetTag(calo->GetTag());
            aodpair.SetTagged(kTRUE);
            
            // Add AOD with pair object to aod branch
            AddAODParticle(aodpair);
            //printf("added pair: naod %d new %d\n",naod,GetOutputAODBranch()->GetEntriesFast());
          }
          
          // Do not add the current calocluster
          if(!fRemoveConvertedPair)
          {
            //printf("TAGGED\n");
            
            // Tag this cluster as likely conversion
            calo ->SetTagged(kTRUE);
            calo2->SetTagged(kTRUE);
          }
        } // selected converted pair
      } // mass cut
      
     } // Mass loop
    
  }// main loop
  
  // Remove entries identified as conversion electrons
  // Revise if this is OK
  if(fRemoveConvertedPair || fAddConvertedPairsToAOD)
  {
    for(Int_t iaod = 0; iaod < naod ; iaod++)
    {
      if(indexConverted[iaod]) GetOutputAODBranch()->RemoveAt(iaod);
    }
    GetOutputAODBranch()->Compress();
  }
  
  delete [] indexConverted;
  
  AliDebug(1,Form("End fill AODs, with %d entries",GetOutputAODBranch()->GetEntriesFast()));
}


//____________________________________________________________
/// Print some relevant parameters set for the analysis.
//____________________________________________________________
void AliAnaPhotonConvInCalo::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Add conversion pair to AOD = %d, removed? %d\n",fAddConvertedPairsToAOD, fRemoveConvertedPair);
  printf("Conversion pair mass cut   = %1.2f and %1.2f\n",fMassCut,fMassCutTight);
  printf("Conversion selection cut : A < %1.2f; %1.3f < Dphi < %1.3f; Deta < %1.3f\n",
         fConvAsymCut,fConvDPhiMinCut, fConvDPhiMaxCut, fConvDEtaCut);
  
  printf("    \n") ;
}

