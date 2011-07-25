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
/* $Id: AliAnaPhoton.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
// Produces input for other analysis classes like AliAnaPi0, 
// AliAnaParticleHadronCorrelation ... 
//
// -- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TH3D.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>
#include "TParticle.h"
#include "TDatabasePDG.h"

// --- Analysis system --- 
#include "AliAnaPhoton.h" 
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODMCParticle.h"
#include "AliMixedEvent.h"


ClassImp(AliAnaPhoton)
  
//____________________________________________________________________________
  AliAnaPhoton::AliAnaPhoton() : 
    AliAnaPartCorrBaseClass(), fCalorimeter(""), 
    fMinDist(0.),fMinDist2(0.),fMinDist3(0.),fRejectTrackMatch(0),
    fTimeCutMin(-1), fTimeCutMax(9999999), fNCellsCut(0), fFillSSHistograms(kFALSE),
    fCheckConversion(kFALSE), fRemoveConvertedPair(kFALSE), fAddConvertedPairsToAOD(kFALSE), fMassCut(0),
    fConvAsymCut(1.), fConvDEtaCut(2.),fConvDPhiMinCut(-1.), fConvDPhiMaxCut(7.), 
    // Histograms
    fhNtraNclu(0), fhNCellsE(0),  
    //Shower shape histograms
    fhNCellsLam0LowE(0),  fhNCellsLam1LowE(0),  fhNCellsDispLowE(0),  
    fhNCellsLam0HighE(0), fhNCellsLam1HighE(0), fhNCellsDispHighE(0),
    fhEtaLam0(0),         fhPhiLam0(0), 
    fhLam1Lam0LowE(0),    fhLam1Lam0HighE(0), 
    fhLam0E(0),           fhLam1E(0),
    // Spectra histograms
    fhEPhoton(0),      fhPtPhoton(0),  fhPhiPhoton(0),  fhEtaPhoton(0),  fhEtaPhiPhoton(0), fhEtaPhi05Photon(0),
    // Conversion histograms
    fhPtPhotonConv(0), fhEtaPhiPhotonConv(0),fhEtaPhi05PhotonConv(0),
    fhConvDeltaEta(0), fhConvDeltaPhi(0),    fhConvDeltaEtaPhi(0), fhConvAsym(0),     fhConvPt(0),
    fhConvDistEta(0),  fhConvDistEn(0),      fhConvDistMass(0),     
    fhConvDistEtaCutEta(0),  fhConvDistEnCutEta(0),       fhConvDistMassCutEta(0),
    fhConvDistEtaCutMass(0), fhConvDistEnCutMass(0), fhConvDistEtaCutAsy(0), fhConvDistEnCutAsy(0), 
    //MC
    fhDeltaE(0), fhDeltaPt(0),fhRatioE(0), fhRatioPt(0),fh2E(0),fh2Pt(0),
    fhPtMCPhoton(0),fhPhiMCPhoton(0),fhEtaMCPhoton(0), 
    fhPtPrompt(0),fhPhiPrompt(0),fhEtaPrompt(0), 
    fhPtFragmentation(0),fhPhiFragmentation(0),fhEtaFragmentation(0), 
    fhPtISR(0),fhPhiISR(0),fhEtaISR(0), 
    fhPtPi0Decay(0),fhPhiPi0Decay(0),fhEtaPi0Decay(0), 
    fhPtOtherDecay(0),  fhPhiOtherDecay(0),  fhEtaOtherDecay(0), 
    fhPtConversion(0),  fhPhiConversion(0),  fhEtaConversion(0),fhEtaPhiConversion(0),fhEtaPhi05Conversion(0),
    fhPtAntiNeutron(0), fhPhiAntiNeutron(0), fhEtaAntiNeutron(0),
    fhPtAntiProton(0),  fhPhiAntiProton(0),  fhEtaAntiProton(0), 
    fhPtUnknown(0),     fhPhiUnknown(0),     fhEtaUnknown(0),
    fhPtConversionTagged(0),        fhPtAntiNeutronTagged(0),       fhPtAntiProtonTagged(0),           fhPtUnknownTagged(0),
    fhConvDeltaEtaMCConversion(0),  fhConvDeltaPhiMCConversion(0),  fhConvDeltaEtaPhiMCConversion(0),  fhConvAsymMCConversion(0),  fhConvPtMCConversion(0),  fhConvDispersionMCConversion(0), fhConvM02MCConversion(0),
    fhConvDeltaEtaMCAntiNeutron(0), fhConvDeltaPhiMCAntiNeutron(0), fhConvDeltaEtaPhiMCAntiNeutron(0), fhConvAsymMCAntiNeutron(0), fhConvPtMCAntiNeutron(0), fhConvDispersionMCAntiNeutron(0),fhConvM02MCAntiNeutron(0),
    fhConvDeltaEtaMCAntiProton(0),  fhConvDeltaPhiMCAntiProton(0),  fhConvDeltaEtaPhiMCAntiProton(0),  fhConvAsymMCAntiProton(0),  fhConvPtMCAntiProton(0),  fhConvDispersionMCAntiProton(0), fhConvM02MCAntiProton(0),
    fhConvDeltaEtaMCString(0),      fhConvDeltaPhiMCString(0),      fhConvDeltaEtaPhiMCString(0),      fhConvAsymMCString(0),      fhConvPtMCString(0),      fhConvDispersionMCString(0),     fhConvM02MCString(0),
    fhConvDistMCConversion(0),      fhConvDistMCConversionCuts(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}//____________________________________________________________________________
AliAnaPhoton::~AliAnaPhoton() 
{
  //dtor

}

//__________________________________________________________________
Bool_t  AliAnaPhoton::ClusterSelected(AliVCluster* calo, TLorentzVector mom) 
{
  //Select clusters if they pass different cuts
  if(GetDebug() > 2) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; Before selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(),
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //.......................................
  //If too small or big energy, skip it
  if(mom.E() < GetMinEnergy() || mom.E() > GetMaxEnergy() ) return kFALSE ; 
  if(GetDebug() > 2) printf("\t Cluster %d Pass E Cut \n",calo->GetID());
  
  //.......................................
  // TOF cut, BE CAREFUL WITH THIS CUT
  Double_t tof = calo->GetTOF()*1e9;
  if(tof < fTimeCutMin || tof > fTimeCutMax) return kFALSE;
  if(GetDebug() > 2)  printf("\t Cluster %d Pass Time Cut \n",calo->GetID());
  
  //.......................................
  if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) return kFALSE;
  if(GetDebug() > 2) printf("\t Cluster %d Pass NCell Cut \n",calo->GetID());
  
  //.......................................
  //Check acceptance selection
  if(IsFiducialCutOn()){
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
    if(! in ) return kFALSE ;
  }
  if(GetDebug() > 2) printf("Fiducial cut passed \n");
  
  //.......................................
  //Skip matched clusters with tracks
  if(fRejectTrackMatch){
    if(IsTrackMatched(calo)) {
      if(GetDebug() > 2) printf("\t Reject track-matched clusters\n");
      return kFALSE ;
    }
    else  
      if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  }// reject matched clusters
  
  //.......................................
  //Check Distance to Bad channel, set bit.
  Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  if(distBad < fMinDist) {//In bad channel (PHOS cristal size 2.2x2.2 cm), EMCAL ( cell units )
    return kFALSE ;
  }
  else if(GetDebug() > 2) printf("\t Bad channel cut passed %4.2f > %2.2f \n",distBad, fMinDist);
  //printf("Cluster %d Pass Bad Dist Cut \n",icalo);

  if(GetDebug() > 0) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; After  selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(), 
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //All checks passed, cluster selected
  return kTRUE;
    
}

//________________________________________________________________________
TObjString *  AliAnaPhoton::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPhoton ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRejectTrackMatch: %d\n",fRejectTrackMatch) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"Conversion Selection: fConvAsymCut %1.2f, fConvDEtaCut %1.2f fConvDPhiCut (%1.2f,%1.2f)\n",
           fConvAsymCut, fConvDEtaCut, fConvDPhiMinCut, fConvDPhiMaxCut) ;
  parList+=onePar ; 
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  parList += GetCaloPID()->GetPIDParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
  
  return new TObjString(parList) ;
}


//________________________________________________________________________
TList *  AliAnaPhoton::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("PhotonHistos") ; 
	
  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	
  
  fhNtraNclu  = new TH2F ("hNtracksNcluster","# of tracks vs # of clusters", 500,0,500, 500,0,500); 
  fhNtraNclu->SetXTitle("# of tracks");
  fhNtraNclu->SetYTitle("# of clusters");
  outputContainer->Add(fhNtraNclu);
  
  
  fhNCellsE  = new TH2F ("hNCellsE","# of cells in cluster vs E of clusters", 100,0, 20, 20,0,20); 
  fhNCellsE->SetXTitle("E (GeV)");
  fhNCellsE->SetYTitle("# of cells in cluster");
  outputContainer->Add(fhNCellsE);  

  if(fFillSSHistograms){
    
    fhNCellsLam0LowE  = new TH2F ("hNCellsLam0LowE","# of cells in cluster vs #lambda_{0}, E < 2 GeV", 20,0, 20, 200,0,1); 
    fhNCellsLam0LowE->SetXTitle("N Cells");
    fhNCellsLam0LowE->SetYTitle("#lambda_{0}");
    outputContainer->Add(fhNCellsLam0LowE);  
    
    
    fhNCellsLam0HighE  = new TH2F ("hNCellsLam0HighE","# of cells in cluster vs #lambda_{0}, E > 2 GeV", 20,0, 20, 200,0,1); 
    fhNCellsLam0HighE->SetXTitle("N Cells");
    fhNCellsLam0HighE->SetYTitle("#lambda_{0}");
    outputContainer->Add(fhNCellsLam0HighE);  
    
    fhNCellsLam1LowE  = new TH2F ("hNCellsLam1LowE","# of cells in cluster vs #lambda_{1}, E < 2 GeV", 20,0, 20, 200,0,1); 
    fhNCellsLam1LowE->SetXTitle("N Cells");
    fhNCellsLam1LowE->SetYTitle("#lambda_{0}");
    outputContainer->Add(fhNCellsLam1LowE);  
    
    fhNCellsLam1HighE  = new TH2F ("hNCellsLam1HighE","# of cells in cluster vs #lambda_{1}, E > 2 GeV", 20,0, 20, 200,0,1); 
    fhNCellsLam1HighE->SetXTitle("N Cells");
    fhNCellsLam1HighE->SetYTitle("#lambda_{0}");
    outputContainer->Add(fhNCellsLam1HighE);  
    
    
    fhNCellsDispLowE  = new TH2F ("hNCellsDispLowE","# of cells in cluster vs dispersion, E < 2 GeV", 20,0, 20, 200,0,2); 
    fhNCellsDispLowE->SetXTitle("N Cells");
    fhNCellsDispLowE->SetYTitle("dispersion");
    outputContainer->Add(fhNCellsDispLowE);  
    
    fhNCellsDispHighE  = new TH2F ("hNCellsDispHighE","# of cells in cluster vs dispersion, E > 2 GeV", 20,0, 20, 200,0,2); 
    fhNCellsDispHighE->SetXTitle("N Cells");
    fhNCellsDispHighE->SetYTitle("dispersion");
    outputContainer->Add(fhNCellsDispHighE);  
    
    fhEtaLam0  = new TH2F ("hEtaLam0","#eta vs #lambda_{0}",200,-0.8,0.8,  200,0,1); 
    fhEtaLam0->SetYTitle("#lambda_{0}");
    fhEtaLam0->SetXTitle("#eta");
    outputContainer->Add(fhEtaLam0);  
    
    fhPhiLam0  = new TH2F ("hPhiLam0","#phi vs #lambda_{0}", 200,80*TMath::DegToRad(),120*TMath::DegToRad(),  200,0,1); 
    fhPhiLam0->SetYTitle("#lambda_{0}");
    fhPhiLam0->SetXTitle("#phi");
    outputContainer->Add(fhPhiLam0);  
    
    
    fhLam1Lam0LowE  = new TH2F ("hLam1Lam0LowE","#lambda_{0} vs #lambda_{1} in cluster of E < 2 GeV",  200,0,1,  200,0,1); 
    fhLam1Lam0LowE->SetYTitle("#lambda_{0}");
    fhLam1Lam0LowE->SetXTitle("#lambda_{1}");
    outputContainer->Add(fhLam1Lam0LowE);  
    
    
    fhLam1Lam0HighE  = new TH2F ("hLam1Lam0HighE","#lambda_{0} vs #lambda_{1} in cluster of E > 2 GeV",  200,0,1,  200,0,1); 
    fhLam1Lam0HighE->SetYTitle("#lambda_{0}");
    fhLam1Lam0HighE->SetXTitle("#lambda_{1}");
    outputContainer->Add(fhLam1Lam0HighE);  
    
    
    fhLam0E  = new TH2F ("hLam0E","#lambda_{0} vs E",  100,0, 20,  200,0,1); 
    fhLam0E->SetYTitle("#lambda_{0}");
    fhLam0E->SetXTitle("E (GeV)");
    outputContainer->Add(fhLam0E);  
    
    fhLam1E  = new TH2F ("hLam1E","#lambda_{1} vs E",  100,0, 20,  200,0,1); 
    fhLam1E->SetYTitle("#lambda_{1}");
    fhLam1E->SetXTitle("E (GeV)");
    outputContainer->Add(fhLam1E);  
  }
  
  fhEPhoton  = new TH1F("hEPhoton","Number of #gamma over calorimeter vs energy",nptbins,ptmin,ptmax); 
  fhEPhoton->SetYTitle("N");
  fhEPhoton->SetXTitle("E_{#gamma}(GeV)");
  outputContainer->Add(fhEPhoton) ;   
  
  fhPtPhoton  = new TH1F("hPtPhoton","Number of #gamma over calorimeter vs p_{T}",nptbins,ptmin,ptmax); 
  fhPtPhoton->SetYTitle("N");
  fhPtPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhPtPhoton) ; 
  
  fhPhiPhoton  = new TH2F
    ("hPhiPhoton","#phi_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
  fhPhiPhoton->SetYTitle("#phi (rad)");
  fhPhiPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhiPhoton) ; 
  
  fhEtaPhoton  = new TH2F
    ("hEtaPhoton","#eta_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
  fhEtaPhoton->SetYTitle("#eta");
  fhEtaPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhEtaPhoton) ;
  
  fhEtaPhiPhoton  = new TH2F
  ("hEtaPhiPhoton","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
  fhEtaPhiPhoton->SetYTitle("#phi (rad)");
  fhEtaPhiPhoton->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiPhoton) ;
  if(GetMinPt() < 0.5){
    fhEtaPhi05Photon  = new TH2F
    ("hEtaPhi05Photon","#eta vs #phi, E > 0.5",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi05Photon->SetYTitle("#phi (rad)");
    fhEtaPhi05Photon->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi05Photon) ;
  }
  
  //Conversion
  if(fCheckConversion){
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
    
  }
  
  if(IsDataMC()){
    fhDeltaE  = new TH1F ("hDeltaE","MC - Reco E ", 200,-50,50); 
    fhDeltaE->SetXTitle("#Delta E (GeV)");
    outputContainer->Add(fhDeltaE);
                
    fhDeltaPt  = new TH1F ("hDeltaPt","MC - Reco p_{T} ", 200,-50,50); 
    fhDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
    outputContainer->Add(fhDeltaPt);

    fhRatioE  = new TH1F ("hRatioE","Reco/MC E ", 200,0,2); 
    fhRatioE->SetXTitle("E_{reco}/E_{gen}");
    outputContainer->Add(fhRatioE);
    
    fhRatioPt  = new TH1F ("hRatioPt","Reco/MC p_{T} ", 200,0,2); 
    fhRatioPt->SetXTitle("p_{T, reco}/p_{T, gen}");
    outputContainer->Add(fhRatioPt);    

    fh2E  = new TH2F ("h2E","E distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fh2E->SetXTitle("E_{rec} (GeV)");
    fh2E->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fh2E);          
    
    fh2Pt  = new TH2F ("h2Pt","p_T distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fh2Pt->SetXTitle("p_{T,rec} (GeV/c)");
    fh2Pt->SetYTitle("p_{T,gen} (GeV/c)");
    outputContainer->Add(fh2Pt);
   
    fhPtMCPhoton  = new TH1F("hPtMCPhoton","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtMCPhoton->SetYTitle("N");
    fhPtMCPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtMCPhoton) ; 
    
    fhPhiMCPhoton  = new TH2F
      ("hPhiMCPhoton","#phi_{#gamma}, #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiMCPhoton->SetYTitle("#phi");
    fhPhiMCPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiMCPhoton) ; 
    
    fhEtaMCPhoton  = new TH2F
      ("hEtaMCPhoton","#eta_{#gamma}, #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaMCPhoton->SetYTitle("#eta");
    fhEtaMCPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaMCPhoton) ;
    
    fhPtPrompt  = new TH1F("hPtMCPrompt","Number of prompt #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtPrompt->SetYTitle("N");
    fhPtPrompt->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPrompt) ; 
    
    fhPhiPrompt  = new TH2F
      ("hPhiMCPrompt","#phi_{#gamma}, prompt #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiPrompt->SetYTitle("#phi");
    fhPhiPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiPrompt) ; 
    
    fhEtaPrompt  = new TH2F
      ("hEtaMCPrompt","#eta_{#gamma}, prompt #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaPrompt->SetYTitle("#eta");
    fhEtaPrompt->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaPrompt) ;
    
    fhPtFragmentation  = new TH1F("hPtMCFragmentation","Number of fragmentation #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtFragmentation->SetYTitle("N");
    fhPtFragmentation->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtFragmentation) ; 
    
    fhPhiFragmentation  = new TH2F
      ("hPhiMCFragmentation","#phi_{#gamma}, fragmentation #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiFragmentation->SetYTitle("#phi");
    fhPhiFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiFragmentation) ; 
    
    fhEtaFragmentation  = new TH2F
      ("hEtaMCFragmentation","#eta_{#gamma}, fragmentation #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaFragmentation->SetYTitle("#eta");
    fhEtaFragmentation->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaFragmentation) ;
    
    fhPtISR  = new TH1F("hPtMCISR","Number of initial state radiation #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtISR->SetYTitle("N");
    fhPtISR->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtISR) ; 
    
    fhPhiISR  = new TH2F
      ("hPhiMCISR","#phi_{#gamma} initial state radiation",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiISR->SetYTitle("#phi");
    fhPhiISR->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiISR) ; 
    
    fhEtaISR  = new TH2F
      ("hEtaMCISR","#eta_{#gamma} initial state radiation",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaISR->SetYTitle("#eta");
    fhEtaISR->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaISR) ;
    
    fhPtPi0Decay  = new TH1F("hPtMCPi0Decay","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtPi0Decay->SetYTitle("N");
    fhPtPi0Decay->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPi0Decay) ; 
    
    fhPhiPi0Decay  = new TH2F
      ("hPhiMCPi0Decay","#phi_{#gamma}, #pi^{0} decay #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiPi0Decay->SetYTitle("#phi");
    fhPhiPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiPi0Decay) ; 
    
    fhEtaPi0Decay  = new TH2F
      ("hEtaMCPi0Decay","#eta_{#gamma}, #pi^{0} #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaPi0Decay->SetYTitle("#eta");
    fhEtaPi0Decay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaPi0Decay) ;
    
    fhPtOtherDecay  = new TH1F("hPtMCOtherDecay","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtOtherDecay->SetYTitle("N");
    fhPtOtherDecay->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtOtherDecay) ; 
    
    fhPhiOtherDecay  = new TH2F
      ("hPhiMCOtherDecay","#phi_{#gamma}, other decay #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiOtherDecay->SetYTitle("#phi");
    fhPhiOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiOtherDecay) ; 
    
    fhEtaOtherDecay  = new TH2F
      ("hEtaMCOtherDecay","#eta_{#gamma}, other decay #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaOtherDecay->SetYTitle("#eta");
    fhEtaOtherDecay->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaOtherDecay) ;
    
    fhPtConversion  = new TH1F("hPtMCConversion","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtConversion->SetYTitle("N");
    fhPtConversion->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtConversion) ; 
    
    fhPhiConversion  = new TH2F
      ("hPhiMCConversion","#phi_{#gamma}, conversion #gamma in MC",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiConversion->SetYTitle("#phi");
    fhPhiConversion->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiConversion) ; 
    
    fhEtaConversion  = new TH2F
      ("hEtaMCConversion","#eta_{#gamma}, conversion #gamma in MC",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaConversion->SetYTitle("#eta");
    fhEtaConversion->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaConversion) ;
    
    fhEtaPhiConversion  = new TH2F
    ("hEtaPhiConversion","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiConversion->SetYTitle("#phi (rad)");
    fhEtaPhiConversion->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiConversion) ;
    
    fhEtaPhi05Conversion  = new TH2F
    ("hEtaPhi05Conversion","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi05Conversion->SetYTitle("#phi (rad)");
    fhEtaPhi05Conversion->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi05Conversion) ;
    
    fhPtAntiNeutron  = new TH1F("hPtMCAntiNeutron","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtAntiNeutron->SetYTitle("N");
    fhPtAntiNeutron->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtAntiNeutron) ; 
    
    fhPhiAntiNeutron  = new TH2F
    ("hPhiMCAntiNeutron","#phi_{#gamma}, unknown origin",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiAntiNeutron->SetYTitle("#phi");
    fhPhiAntiNeutron->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiAntiNeutron) ; 
    
    fhEtaAntiNeutron  = new TH2F
    ("hEtaMCAntiNeutron","#eta_{#gamma}, unknown origin",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaAntiNeutron->SetYTitle("#eta");
    fhEtaAntiNeutron->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaAntiNeutron) ;
        
    fhPtAntiProton  = new TH1F("hPtMCAntiProton","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtAntiProton->SetYTitle("N");
    fhPtAntiProton->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtAntiProton) ; 
    
    fhPhiAntiProton  = new TH2F
    ("hPhiMCAntiProton","#phi_{#gamma}, unknown origin",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiAntiProton->SetYTitle("#phi");
    fhPhiAntiProton->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiAntiProton) ; 
    
    fhEtaAntiProton  = new TH2F
    ("hEtaMCAntiProton","#eta_{#gamma}, unknown origin",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaAntiProton->SetYTitle("#eta");
    fhEtaAntiProton->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaAntiProton) ;
    
    fhPtUnknown  = new TH1F("hPtMCUnknown","Number of #gamma over calorimeter",nptbins,ptmin,ptmax); 
    fhPtUnknown->SetYTitle("N");
    fhPtUnknown->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtUnknown) ; 
    
    fhPhiUnknown  = new TH2F
      ("hPhiMCUnknown","#phi_{#gamma}, unknown origin",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
    fhPhiUnknown->SetYTitle("#phi");
    fhPhiUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhPhiUnknown) ; 
    
    fhEtaUnknown  = new TH2F
      ("hEtaMCUnknown","#eta_{#gamma}, unknown origin",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
    fhEtaUnknown->SetYTitle("#eta");
    fhEtaUnknown->SetXTitle("p_{T #gamma} (GeV/c)");
    outputContainer->Add(fhEtaUnknown) ;
	
    if(fCheckConversion){  
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
    
  }//Histos with MC
    
  //Store calo PID histograms
  if(fRejectTrackMatch){
    TList * caloPIDHistos = GetCaloPID()->GetCreateOutputObjects() ;
    for(Int_t i = 0; i < caloPIDHistos->GetEntries(); i++) {
      outputContainer->Add(caloPIDHistos->At(i)) ;
    }
    delete caloPIDHistos;
  }
  
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaPhoton::Init()
{
  
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPhoton::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPhoton::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}


//____________________________________________________________________________
void AliAnaPhoton::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaPhoton_");

  fCalorimeter = "EMCAL" ;
  fMinDist     = 2.;
  fMinDist2    = 4.;
  fMinDist3    = 5.;
  fMassCut     = 0.03; //30 MeV
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  fNCellsCut   = 0;
	
  fRejectTrackMatch       = kTRUE ;
  fCheckConversion        = kFALSE;
  fRemoveConvertedPair    = kFALSE;
  fAddConvertedPairsToAOD = kFALSE;
	
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillAOD() 
{
  //Do photon analysis and fill aods
  
  //Get the vertex 
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  //Select the Calorimeter of the photon
  TObjArray * pl = 0x0; 
  if(fCalorimeter == "PHOS")
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  if(!pl) {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }

  //Init arrays, variables, get number of clusters
  TLorentzVector mom, mom2 ;
  Int_t nCaloClusters = pl->GetEntriesFast();
  //List to be used in conversion analysis, to tag the cluster as candidate for conversion
  Bool_t * indexConverted = 0x0;
  if(fCheckConversion){
    indexConverted = new Bool_t[nCaloClusters];
    for (Int_t i = 0; i < nCaloClusters; i++) 
      indexConverted[i] = kFALSE;
	}
  
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);

  //----------------------------------------------------
  // Fill AOD with PHOS/EMCAL AliAODPWG4Particle objects
  //----------------------------------------------------
  // Loop on clusters
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++){    
	  
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));	
    //printf("calo %d, %f\n",icalo,calo->E());

    //Get the index where the cluster comes, to retrieve the corresponding vertex
    Int_t evtIndex = 0 ; 
    if (GetMixedEvent()) {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
      //Get the vertex and check it is not too large in z
      if(TMath::Abs(GetVertex(evtIndex)[2])> GetZvertexCut()) continue;
    }

    //Cluster selection, not charged, with photon id and in fiducial cut
	  
    //Input from second AOD?
    //Int_t input = 0;
    //    if (fCalorimeter == "EMCAL" && GetReader()->GetEMCALClustersNormalInputEntries() <= icalo) 
    //      input = 1 ;
    //    else if(fCalorimeter == "PHOS"  && GetReader()->GetPHOSClustersNormalInputEntries()  <= icalo) 
    //      input = 1;
	  
    //Get Momentum vector, 
    //if (input == 0) 
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;}//Assume that come from vertex in straight line
    else{
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }

    //    else if(input == 1) 
    //      calo->GetMomentum(mom,vertex2);//Assume that come from vertex in straight line  
    
    //--------------------------------------
    // Cluster selection
    //--------------------------------------
    
    if(!ClusterSelected(calo,mom)) continue;
    
    //----------------------------
    //Create AOD for analysis
    //----------------------------
    AliAODPWG4Particle aodph = AliAODPWG4Particle(mom);
    
    //...............................................
    //Set the indeces of the original caloclusters (MC, ID), and calorimeter  
    Int_t label = calo->GetLabel();
    aodph.SetLabel(label);
    //aodph.SetInputFileIndex(input);   
    aodph.SetCaloLabel(calo->GetID(),-1);
    aodph.SetDetector(fCalorimeter);
    //printf("Index %d, Id %d, iaod %d\n",icalo, calo->GetID(),GetOutputAODBranch()->GetEntriesFast());

    //...............................................
    //Set bad channel distance bit
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if     (distBad > fMinDist3) aodph.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodph.SetDistToBad(1) ; 
    else                         aodph.SetDistToBad(0) ;
    //printf("DistBad %f Bit %d\n",distBad, aodph.DistToBad());
    
    //...............................................
    //Set number of cells in this cluster
    //Temporary patch FIXME
    aodph.SetBtag(calo->GetNCells());
    // MEFIX
    
    //-------------------------------------
    //PID selection or bit setting
    //-------------------------------------
    // MC
    if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
      //Get most probable PID, check PID weights (in MC this option is mandatory)
      aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());	 
      //If primary is not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;
    }	
    //...............................................
    // Data, PID check on
    else if(IsCaloPIDOn()){
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for EMCal.		
      if(!IsCaloPIDRecalculationOn())
        aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      else
        aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());
      
      //If cluster does not pass pid, not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;			
      
    }
    //...............................................
    // Data, PID check off
    else{
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph, GetCaloUtils());
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PID Bits set \n");		
    }
    
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",aodph.Pt(), aodph.GetIdentifiedParticleType());
    
    //--------------------------------------------------------------------------------------
    //Play with the MC stack if available
    //--------------------------------------------------------------------------------------

    //Check origin of the candidates
    if(IsDataMC()){
      aodph.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(), aodph.GetInputFileIndex()));
      if(GetDebug() > 0)
        printf("AliAnaPhoton::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also   
    
    //--------------------------------------------------------------------------------------
    // Conversions pairs analysis
    // Check if cluster comes from a conversion in the material in front of the calorimeter
    // Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
    //--------------------------------------------------------------------------------------

    // Do analysis only if there are more than one cluster
    if( nCaloClusters > 1 && fCheckConversion){
      Bool_t bConverted = kFALSE;
      Int_t id2 = -1;
		  
      //Check if set previously as converted couple, if so skip its use.
      if (indexConverted[icalo]) continue;
		  
      // Second cluster loop
      for(Int_t jcalo = icalo + 1 ; jcalo < nCaloClusters ; jcalo++) {
        //Check if set previously as converted couple, if so skip its use.
        if (indexConverted[jcalo]) continue;
        //printf("Check Conversion indeces %d and %d\n",icalo,jcalo);
        AliVCluster * calo2 =  (AliVCluster*) (pl->At(jcalo));              //Get cluster kinematics
        
        
        //Mixed event, get index of event
        Int_t evtIndex2 = 0 ; 
        if (GetMixedEvent()) {
          evtIndex2=GetMixedEvent()->EventIndexForCaloCluster(calo2->GetID()) ; 
          
        }      
        
        //Get kinematics of second cluster
        if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
          calo2->GetMomentum(mom2,GetVertex(evtIndex2)) ;}//Assume that come from vertex in straight line
        else{
          Double_t vertex[]={0,0,0};
          calo2->GetMomentum(mom2,vertex) ;
        }
        
        //--------------------------------------
        // Cluster selection
        //--------------------------------------
        
        if(!ClusterSelected(calo2,mom2)) continue;  
        
        //................................................
        // Get TOF of each cluster in pair, calculate difference if small, 
        // take this pair. Way to reject clusters from hadrons (or pileup?)
        Double_t t12diff = calo2->GetTOF()-calo->GetTOF()*1e9;
        if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
        
        //................................................
        //Get mass of pair, if small, take this pair.
        Float_t pairM     = (mom+mom2).M();
        //printf("\t both in calo, mass %f, cut %f\n",pairM,fMassCut);
        if(pairM < fMassCut){  
          aodph.SetTagged(kFALSE);
          id2 = calo2->GetID();
          indexConverted[icalo]=kTRUE;
          indexConverted[jcalo]=kTRUE; 
          Float_t asymmetry = TMath::Abs(mom.E()-mom2.E())/(mom.E()+mom2.E());
          Float_t dPhi      = mom.Phi()-mom2.Phi();
          Float_t dEta      = mom.Eta()-mom2.Eta();  
          
          //...............................................
          //Fill few histograms with kinematics of the pair
          //FIXME, move all this to MakeAnalysisFillHistograms ...
          
          fhConvDeltaEta   ->Fill( pairM, dPhi      );
          fhConvDeltaPhi   ->Fill( pairM, dEta      );
          fhConvAsym       ->Fill( pairM, asymmetry );
          fhConvDeltaEtaPhi->Fill( dEta , dPhi      );
          fhConvPt         ->Fill( pairM, (mom+mom2).Pt());          
          
          //Estimate conversion distance, T. Awes, M. Ivanov
          //Under the assumption that the pair has zero mass, and that each electron 
          //of the pair has the same momentum, they will each have the same bend radius 
          //given by R=p/(qB) = p / (300 B) with p in [MeV/c], B in [Tesla] and R in [m]. 
          //With nominal ALICE magnet current of 30kA B=0.5T, and so with E_cluster=p,  
          //R = E/1.5 [cm].  Under these assumptions, the distance from the conversion 
          //point to the EMCal can be related to the separation distance, L=2y, on the EMCal 
          //as d = sqrt(R^2 -(R-y)^2) = sqrt(2Ry - y^2). And since R>>y we can write as 
          //d = sqrt(E*L/1.5) where E is the cluster energy and L is the distance in cm between 
          //the clusters.
          Float_t pos1[3];
          calo->GetPosition(pos1); 
          Float_t pos2[3];
          calo2->GetPosition(pos2); 
          Float_t clustDist = TMath::Sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+
                                  (pos1[1]-pos2[1])*(pos1[1]-pos2[1])+
                                  (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
          
          Float_t convDist  = TMath::Sqrt(mom.E() *clustDist*0.01/0.15);
          Float_t convDist2 = TMath::Sqrt(mom2.E()*clustDist*0.01/0.15);
          //printf("l = %f, e1 = %f, d1=%f, e2 = %f, d2=%f\n",clustDist,mom.E(),convDist,mom2.E(),convDist2);
          if(GetDebug() > 2)
            printf("AliAnaPhoton::MakeAnalysisFillAOD(): Pair with mass %2.3f < %2.3f, %1.2f < dPhi %2.2f < %2.2f, dEta %f < %2.2f, asymmetry %2.2f< %2.2f; \n    cluster1 id %d, e %2.3f  SM %d, eta %2.3f, phi %2.3f ; \n    cluster2 id %d, e %2.3f, SM %d,eta %2.3f, phi %2.3f\n",
                   pairM,fMassCut,fConvDPhiMinCut, dPhi, fConvDPhiMaxCut, dEta, fConvDEtaCut, asymmetry, fConvAsymCut,
                   calo->GetID(),calo->E(),GetCaloUtils()->GetModuleNumber(calo), mom.Eta(), mom.Phi(),
                   id2, calo2->E(), GetCaloUtils()->GetModuleNumber(calo2),mom2.Eta(), mom2.Phi());
          
          fhConvDistEta ->Fill(mom .Eta(), convDist );
          fhConvDistEta ->Fill(mom2.Eta(), convDist2);
          fhConvDistEn  ->Fill(mom .E(), convDist );
          fhConvDistEn  ->Fill(mom2.E(), convDist2);        
          fhConvDistMass->Fill((mom+mom2).M(), convDist );
          //dEta cut
          if(dEta<0.05){
            fhConvDistEtaCutEta ->Fill(mom .Eta(), convDist );
            fhConvDistEtaCutEta ->Fill(mom2.Eta(), convDist2);
            fhConvDistEnCutEta  ->Fill(mom .E(),   convDist );
            fhConvDistEnCutEta  ->Fill(mom2.E(),   convDist2);        
            fhConvDistMassCutEta->Fill((mom+mom2).M(), convDist );
            //mass cut
            if(pairM<0.01){//10 MeV
              fhConvDistEtaCutMass ->Fill(mom .Eta(), convDist );
              fhConvDistEtaCutMass ->Fill(mom2.Eta(), convDist2);
              fhConvDistEnCutMass  ->Fill(mom .E(),   convDist );
              fhConvDistEnCutMass  ->Fill(mom2.E(),   convDist2);        
              // asymmetry cut
              if(asymmetry<0.1){
                fhConvDistEtaCutAsy ->Fill(mom .Eta(), convDist );
                fhConvDistEtaCutAsy ->Fill(mom2.Eta(), convDist2);
                fhConvDistEnCutAsy  ->Fill(mom .E(),   convDist );
                fhConvDistEnCutAsy  ->Fill(mom2.E(),   convDist2); 
              }//asymmetry cut
            }//mass cut            
          }//dEta cut
          
          //...............................................
          //Select pairs in a eta-phi window
          if(TMath::Abs(dEta) < fConvDEtaCut    && 
             TMath::Abs(dPhi) < fConvDPhiMaxCut &&
             TMath::Abs(dPhi) > fConvDPhiMinCut && 
             asymmetry        < fConvAsymCut       ){
              bConverted = kTRUE;          
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
            Int_t ancLabel  = GetMCAnalysisUtils()->CheckCommonAncestor(calo->GetLabel(), calo2->GetLabel(), 
                                                                        GetReader(), ancPDG, ancStatus, momentum, prodVertex);
            
            // printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Common ancestor label %d, pdg %d, name %s, status %d; \n",
            //                          ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
            
            Int_t tag2 = GetMCAnalysisUtils()->CheckOrigin(calo2->GetLabels(),calo2->GetNLabels(),GetReader(), 0);
            if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCConversion)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) && (ancPDG==22 || TMath::Abs(ancPDG)==11) && ancLabel > -1){
                fhConvDeltaEtaMCConversion   ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCConversion   ->Fill( pairM, dPhi      );
                fhConvAsymMCConversion       ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCConversion->Fill( dEta , dPhi      );
                fhConvPtMCConversion         ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCConversion ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCConversion        ->Fill( calo->GetM02(), calo2->GetM02());
                fhConvDistMCConversion       ->Fill( convDist , prodVertex.Mag() );
                fhConvDistMCConversion       ->Fill( convDist2, prodVertex.Mag() );
                
                if(dEta<0.05 && pairM<0.01 && asymmetry<0.1){
                  fhConvDistMCConversionCuts->Fill( convDist , prodVertex.Mag() );
                  fhConvDistMCConversionCuts->Fill( convDist2, prodVertex.Mag() );
                }

              }              
            }
            else if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCAntiNeutron)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiNeutron) && ancPDG==-2112 && ancLabel > -1){
                fhConvDeltaEtaMCAntiNeutron    ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCAntiNeutron    ->Fill( pairM, dPhi      );
                fhConvAsymMCAntiNeutron        ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCAntiNeutron ->Fill( dEta , dPhi      );
                fhConvPtMCAntiNeutron          ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCAntiNeutron  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCAntiNeutron         ->Fill( calo->GetM02(), calo2->GetM02());
              }
            }
            else if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCAntiProton)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiProton) && ancPDG==-2212 && ancLabel > -1){
                fhConvDeltaEtaMCAntiProton    ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCAntiProton    ->Fill( pairM, dPhi      );
                fhConvAsymMCAntiProton        ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCAntiProton ->Fill( dEta , dPhi      );
                fhConvPtMCAntiProton          ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCAntiProton  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCAntiProton         ->Fill( calo->GetM02(), calo2->GetM02());
              }
            }
            
            //Pairs coming from fragmenting pairs.
            if(ancPDG < 22 && ancLabel > 7 && (ancStatus == 11 || ancStatus == 12) ){
              fhConvDeltaEtaMCString    ->Fill( pairM, dPhi);
              fhConvDeltaPhiMCString    ->Fill( pairM, dPhi);
              fhConvAsymMCString        ->Fill( pairM, TMath::Abs(mom.E()-mom2.E())/(mom.E()+mom2.E()) );
              fhConvDeltaEtaPhiMCString ->Fill( dPhi, dPhi );
              fhConvPtMCString          ->Fill( pairM, (mom+mom2).Pt());
              fhConvDispersionMCString  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
              fhConvM02MCString         ->Fill( calo->GetM02(), calo2->GetM02());
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
          TLorentzVector mpair = mom+mom2;
          AliAODPWG4Particle aodpair = AliAODPWG4Particle(mpair);
          aodpair.SetLabel(aodph.GetLabel());
          //aodpair.SetInputFileIndex(input);
          
          //printf("Index %d, Id %d\n",icalo, calo->GetID());
          //Set the indeces of the original caloclusters  
          aodpair.SetCaloLabel(calo->GetID(),id2);
          aodpair.SetDetector(fCalorimeter);
          aodpair.SetIdentifiedParticleType(aodph.GetIdentifiedParticleType());
          aodpair.SetTag(aodph.GetTag());
          aodpair.SetTagged(kTRUE);
          //Add AOD with pair object to aod branch
          AddAODParticle(aodpair);
          //printf("\t \t both added pair\n");
        }
        
        //Do not add the current calocluster
        if(fRemoveConvertedPair) continue;
        else {
          //printf("TAGGED\n");
          //Tag this cluster as likely conversion
          aodph.SetTagged(kTRUE);
        }
      }//converted pair
    }//check conversion
    //printf("\t \t added single cluster %d\n",icalo);
	  
    //FIXME, this to MakeAnalysisFillHistograms ...
    fhNCellsE->Fill(aodph.E(),calo->GetNCells());
    if(fFillSSHistograms){
      if(calo->E()<2){
        fhNCellsLam0LowE->Fill(calo->GetNCells(),calo->GetM02());
        fhNCellsLam1LowE->Fill(calo->GetNCells(),calo->GetM20());
        fhNCellsDispLowE->Fill(calo->GetNCells(),calo->GetDispersion());
        fhLam1Lam0LowE  ->Fill(calo->GetM20(),   calo->GetM02());
      }
      else {
        fhNCellsLam0HighE->Fill(calo->GetNCells(),calo->GetM02());
        fhNCellsLam1HighE->Fill(calo->GetNCells(),calo->GetM20());
        fhNCellsDispHighE->Fill(calo->GetNCells(),calo->GetDispersion());
        fhLam1Lam0HighE  ->Fill(calo->GetM20(),   calo->GetM02());
      }
      
      fhEtaLam0->Fill(aodph.Eta(), calo->GetM02());
      fhPhiLam0->Fill(aodph.Phi(), calo->GetM02());
      fhLam0E  ->Fill(aodph.E(),   calo->GetM02());
      fhLam1E  ->Fill(aodph.E(),   calo->GetM20());
    }
    
    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  
  delete [] indexConverted;
	
  if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillHistograms() 
{
  //Fill histograms
  
  //-------------------------------------------------------------------
	// Access MC information in stack if requested, check that it exists.	
	AliStack * stack = 0x0;
	TParticle * primary = 0x0;   
	TClonesArray * mcparticles0 = 0x0;
	//TClonesArray * mcparticles1 = 0x0;
	AliAODMCParticle * aodprimary = 0x0; 
	if(IsDataMC()){
		
		if(GetReader()->ReadStack()){
			stack =  GetMCStack() ;
			if(!stack) {
				printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called? STOP\n");
				abort();
			}
      
		}
		else if(GetReader()->ReadAODMCParticles()){
      
			//Get the list of MC particles
			mcparticles0 = GetReader()->GetAODMCParticles(0);
			if(!mcparticles0 && GetDebug() > 0) 	{
				printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
			}	
      //			if(GetReader()->GetSecondInputAODTree()){
      //				mcparticles1 = GetReader()->GetAODMCParticles(1);
      //				if(!mcparticles1 && GetDebug() > 0) 	{
      //					printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
      //				}
      //			}		
			
		}
	}// is data and MC
	
  
  // Get vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  //fhVertex->Fill(v[0],v[1],v[2]);  
  if(TMath::Abs(v[2]) > GetZvertexCut()) return ; // done elsewhere for Single Event analysis, but there for mixed event
  
  //----------------------------------
	//Loop on stored AOD photons
	Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  fhNtraNclu->Fill(GetReader()->GetTrackMultiplicity(), naod);
	if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
	
	for(Int_t iaod = 0; iaod < naod ; iaod++){
	  AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
	  Int_t pdg = ph->GetIdentifiedParticleType();
	  
	  if(GetDebug() > 3) 
	    printf("AliAnaPhoton::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n", ph->GetIdentifiedParticleType(),ph->GetTag(), (ph->GetDetector()).Data()) ;
	  
	  //If PID used, fill histos with photons in Calorimeter fCalorimeter
	  if(IsCaloPIDOn() && pdg != AliCaloPID::kPhoton) continue; 
	  if(ph->GetDetector() != fCalorimeter) continue;
	  
	  if(GetDebug() > 2) 
	    printf("AliAnaPhoton::MakeAnalysisFillHistograms() - ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
	  
    //................................
	  //Fill photon histograms 
	  Float_t ptcluster  = ph->Pt();
	  Float_t phicluster = ph->Phi();
	  Float_t etacluster = ph->Eta();
	  Float_t ecluster   = ph->E();
	  
    fhEPhoton   ->Fill(ecluster);
	  fhPtPhoton  ->Fill(ptcluster);
	  fhPhiPhoton ->Fill(ptcluster,phicluster);
	  fhEtaPhoton ->Fill(ptcluster,etacluster);    
    if(ecluster > 0.5)         fhEtaPhiPhoton ->Fill(etacluster, phicluster);
    else if(GetMinPt() < 0.5) fhEtaPhi05Photon->Fill(etacluster, phicluster);
    
    if(fCheckConversion &&ph->IsTagged()){
      fhPtPhotonConv->Fill(ptcluster);
      if(ecluster > 0.5)        fhEtaPhiPhotonConv  ->Fill(etacluster, phicluster);
      else if(GetMinPt() < 0.5) fhEtaPhi05PhotonConv->Fill(etacluster, phicluster);
    }
    
    //.......................................
	  //Play with the MC data if available
	  if(IsDataMC()){
	    
	    Int_t tag =ph->GetTag();

	    if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
      {
        fhPtMCPhoton  ->Fill(ptcluster);
        fhPhiMCPhoton ->Fill(ptcluster,phicluster);
        fhEtaMCPhoton ->Fill(ptcluster,etacluster);
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion))
        {
          fhPtConversion  ->Fill(ptcluster);
          fhPhiConversion ->Fill(ptcluster,phicluster);
          fhEtaConversion ->Fill(ptcluster,etacluster);
          if(ph->IsTagged()) fhPtConversionTagged ->Fill(ptcluster);
          if(ptcluster > 0.5)fhEtaPhiConversion   ->Fill(etacluster,phicluster);
          else               fhEtaPhi05Conversion ->Fill(etacluster,phicluster);
        }			
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt)){
          fhPtPrompt  ->Fill(ptcluster);
          fhPhiPrompt ->Fill(ptcluster,phicluster);
          fhEtaPrompt ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
        {
          fhPtFragmentation  ->Fill(ptcluster);
          fhPhiFragmentation ->Fill(ptcluster,phicluster);
          fhEtaFragmentation ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
        {
          fhPtISR  ->Fill(ptcluster);
          fhPhiISR ->Fill(ptcluster,phicluster);
          fhEtaISR ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
        {
          fhPtPi0Decay  ->Fill(ptcluster);
          fhPhiPi0Decay ->Fill(ptcluster,phicluster);
          fhEtaPi0Decay ->Fill(ptcluster,etacluster);
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
        {
          fhPtOtherDecay  ->Fill(ptcluster);
          fhPhiOtherDecay ->Fill(ptcluster,phicluster);
          fhEtaOtherDecay ->Fill(ptcluster,etacluster);
        }
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron))
      {
        fhPtAntiNeutron  ->Fill(ptcluster);
        fhPhiAntiNeutron ->Fill(ptcluster,phicluster);
        fhEtaAntiNeutron ->Fill(ptcluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtAntiNeutronTagged ->Fill(ptcluster);

      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton))
      {
        fhPtAntiProton  ->Fill(ptcluster);
        fhPhiAntiProton ->Fill(ptcluster,phicluster);
        fhEtaAntiProton ->Fill(ptcluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtAntiProtonTagged ->Fill(ptcluster);

      }      
	    else{
	      fhPtUnknown  ->Fill(ptcluster);
	      fhPhiUnknown ->Fill(ptcluster,phicluster);
	      fhEtaUnknown ->Fill(ptcluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtUnknownTagged ->Fill(ptcluster);

	      
        //		 printf(" AliAnaPhoton::MakeAnalysisFillHistograms() - Label %d, pT %2.3f Unknown, bits set: ",
        //					ph->GetLabel(),ph->Pt());
        //		  for(Int_t i = 0; i < 20; i++) {
        //			  if(GetMCAnalysisUtils()->CheckTagBit(tag,i)) printf(" %d, ",i);
        //		  }
        //		  printf("\n");
        
	    }
	    
	    //....................................................................
	    // Access MC information in stack if requested, check that it exists.
	    Int_t label =ph->GetLabel();
	    if(label < 0) {
	      printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
	      continue;
	    }
	    
	    Float_t eprim   = 0;
	    Float_t ptprim  = 0;
	    if(GetReader()->ReadStack()){
	      
	      if(label >=  stack->GetNtrack()) {
          if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
          continue ;
	      }
	      
	      primary = stack->Particle(label);
	      if(!primary){
          printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
	      }
	      eprim   = primary->Energy();
	      ptprim  = primary->Pt();		
	      
	    }
	    else if(GetReader()->ReadAODMCParticles()){
	      //Check which is the input
	      if(ph->GetInputFileIndex() == 0){
          if(!mcparticles0) continue;
          if(label >=  mcparticles0->GetEntriesFast()) {
            if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
                                       label, mcparticles0->GetEntriesFast());
            continue ;
          }
          //Get the particle
          aodprimary = (AliAODMCParticle*) mcparticles0->At(label);
          
	      }
//	      else {//Second input
//          if(!mcparticles1) continue;
//          if(label >=  mcparticles1->GetEntriesFast()) {
//            if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
//                                       label, mcparticles1->GetEntriesFast());
//            continue ;
//          }
//          //Get the particle
//          aodprimary = (AliAODMCParticle*) mcparticles1->At(label);
//          
//	      }//second input
	      
	      if(!aodprimary){
          printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
	      }
	      
	      eprim   = aodprimary->E();
	      ptprim  = aodprimary->Pt();
	      
	    }
	    
	    fh2E     ->Fill(ecluster, eprim);
	    fh2Pt    ->Fill(ptcluster, ptprim);     
	    fhDeltaE ->Fill(eprim-ecluster);
	    fhDeltaPt->Fill(ptprim-ptcluster);     
	    if(eprim > 0)  fhRatioE  ->Fill(ecluster/eprim);
	    if(ptprim > 0) fhRatioPt ->Fill(ptcluster/ptprim); 		
	    
	  }//Histograms with MC
	  
	}// aod loop
	
}


//__________________________________________________________________
void AliAnaPhoton::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
  printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
  printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  printf("Reject clusters with a track matched = %d\n",fRejectTrackMatch);
  printf("Check Pair Conversion                = %d\n",fCheckConversion);
  printf("Add conversion pair to AOD           = %d\n",fAddConvertedPairsToAOD);
  printf("Conversion pair mass cut             = %f\n",fMassCut);
  printf("Conversion selection cut : A < %1.2f; %1.3f < Dphi < %1.3f; Deta < %1.3f\n",
         fConvAsymCut,fConvDPhiMinCut, fConvDPhiMaxCut, fConvDEtaCut);
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %d \n", fNCellsCut);
  printf("    \n") ;
	
} 
