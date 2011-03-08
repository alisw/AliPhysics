/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: $ */

//_________________________________________________________________________
// Class to collect two-photon invariant mass distributions for
// extracting raw pi0 yield.
// Input is produced by AliAnaPhoton (or any other analysis producing output AliAODPWG4Particles), 
// it will do nothing if executed alone
//
//-- Author: Dmitri Peressounko (RRC "KI") 
//-- Adapted to PartCorr frame by Lamia Benhabib (SUBATECH)
//-- and Gustavo Conesa (INFN-Frascati)
//_________________________________________________________________________


// --- ROOT system ---
#include "TH3.h"
#include "TH2D.h"
//#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TDatabasePDG.h"

//---- AliRoot system ----
#include "AliAnaPi0.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliNeutralMesonSelection.h"
#include "AliMixedEvent.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnaPi0)

//________________________________________________________________________________________________________________________________________________  
AliAnaPi0::AliAnaPi0() : AliAnaPartCorrBaseClass(),
fDoOwnMix(kFALSE),fNCentrBin(0),//fNZvertBin(0),fNrpBin(0),
fNmaxMixEv(0), fCalorimeter(""),
fNModules(12), fUseAngleCut(kFALSE), fUseAngleEDepCut(kFALSE),fAngleCut(0), fAngleMaxCut(7.),fEventsList(0x0), fMultiCutAna(kFALSE), fMultiCutAnaSim(kFALSE),
fNPtCuts(0),fNAsymCuts(0), fNCellNCuts(0),fNPIDBits(0),  fMakeInvPtPlots(kFALSE), fSameSM(kFALSE),
fUseTrackMultBins(kFALSE),fUsePhotonMultBins(kFALSE),fUseAverCellEBins(kFALSE), fUseAverClusterEBins(kFALSE),
fUseAverClusterEDenBins(0), //fUseAverClusterPairRBins(0), fUseAverClusterPairRWeightBins(0), fUseEMaxBins(0),
fFillBadDistHisto(kFALSE),
fhAverTotECluster(0), fhAverTotECell(0), fhAverTotECellvsCluster(0),
fhEDensityCluster(0), fhEDensityCell(0), fhEDensityCellvsCluster(0),
//fhClusterPairDist(0), fhClusterPairDistWeight(0),  fhAverClusterPairDist(0), fhAverClusterPairDistWeight(0),
//fhAverClusterPairDistvsAverE(0), fhAverClusterPairDistWeightvsAverE(0),fhAverClusterPairDistvsN(0), fhAverClusterPairDistWeightvsN(0),
//fhMaxEvsClustMult(0), fhMaxEvsClustEDen(0),
fhReMod(0x0),    fhReSameSideEMCALMod(0x0), fhReSameSectorEMCALMod(0x0), fhReDiffPHOSMod(0x0), 
fhMiMod(0x0),    fhMiSameSideEMCALMod(0x0), fhMiSameSectorEMCALMod(0x0), fhMiDiffPHOSMod(0x0),
fhReConv(0x0),   fhMiConv(0x0),    fhReConv2(0x0),  fhMiConv2(0x0),
fhRe1(0x0),      fhMi1(0x0),       fhRe2(0x0),      fhMi2(0x0),      fhRe3(0x0),      fhMi3(0x0),
fhReInvPt1(0x0), fhMiInvPt1(0x0),  fhReInvPt2(0x0), fhMiInvPt2(0x0), fhReInvPt3(0x0), fhMiInvPt3(0x0),
fhRePtNCellAsymCuts(0x0), fhRePtNCellAsymCutsSM0(0x0), fhRePtNCellAsymCutsSM1(0x0), fhRePtNCellAsymCutsSM2(0x0), fhRePtNCellAsymCutsSM3(0x0), fhMiPtNCellAsymCuts(0x0),
fhRePIDBits(0x0),fhRePtMult(0x0), fhRePtAsym(0x0), fhRePtAsymPi0(0x0),fhRePtAsymEta(0x0),  
fhEvents(0x0),   fhCentrality(0x0),
fhRealOpeningAngle(0x0),fhRealCosOpeningAngle(0x0), fhMixedOpeningAngle(0x0),fhMixedCosOpeningAngle(0x0),
fhPrimPi0Pt(0x0), fhPrimPi0AccPt(0x0), fhPrimPi0Y(0x0), fhPrimPi0AccY(0x0), fhPrimPi0Phi(0x0), fhPrimPi0AccPhi(0x0),
fhPrimPi0OpeningAngle(0x0), fhPrimPi0CosOpeningAngle(0x0),
fhPrimEtaPt(0x0), fhPrimEtaAccPt(0x0), fhPrimEtaY(0x0), fhPrimEtaAccY(0x0), fhPrimEtaPhi(0x0), fhPrimEtaAccPhi(0x0),
fhPrimPi0PtOrigin(0x0),  fhPrimEtaPtOrigin(0x0), 
fhMCOrgMass(),fhMCOrgAsym(),  fhMCOrgDeltaEta(),fhMCOrgDeltaPhi(),
fhMCPi0MassPtRec(), fhMCPi0MassPtTrue(), fhMCPi0PtTruePtRec(), fhMCEtaMassPtRec(), fhMCEtaMassPtTrue(), fhMCEtaPtTruePtRec(),
fhMCPi0PtOrigin(0x0),  fhMCEtaPtOrigin(0x0)
{
//Default Ctor
 InitParameters();
 
}

//________________________________________________________________________________________________________________________________________________
AliAnaPi0::~AliAnaPi0() {
  // Remove event containers
  
  if(fDoOwnMix && fEventsList){
    for(Int_t ic=0; ic<fNCentrBin; ic++){
      for(Int_t iz=0; iz<GetNZvertBin(); iz++){
        for(Int_t irp=0; irp<GetNRPBin(); irp++){
          fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp]->Delete() ;
          delete fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp] ;
        }
      }
    }
    delete[] fEventsList; 
    fEventsList=0 ;
  }
	
}

//________________________________________________________________________________________________________________________________________________
void AliAnaPi0::InitParameters()
{
//Init parameters when first called the analysis
//Set default parameters
  SetInputAODName("PWG4Particle");
  
  AddToHistogramsName("AnaPi0_");
  fNModules = 12; // set maximum to maximum number of EMCAL modules
  fNCentrBin = 1;
//  fNZvertBin = 1;
//  fNrpBin    = 1;
  fNmaxMixEv = 10;
 
  fCalorimeter  = "PHOS";
  fUseAngleCut = kFALSE;
  fUseAngleEDepCut = kFALSE;
  fAngleCut    = 0.; 
  fAngleMaxCut = TMath::Pi(); 

  fMultiCutAna = kFALSE;
  
  fNPtCuts = 3;
  fPtCuts[0] = 0.; fPtCuts[1] = 0.3;   fPtCuts[2] = 0.5;
  for(Int_t i = fNPtCuts; i < 10; i++)fPtCuts[i] = 0.;
  
  fNAsymCuts = 4;
  fAsymCuts[0] = 1.;  fAsymCuts[1] = 0.8; fAsymCuts[2] = 0.6;   fAsymCuts[3] = 0.1;    
  for(Int_t i = fNAsymCuts; i < 10; i++)fAsymCuts[i] = 0.;

  fNCellNCuts = 3;
  fCellNCuts[0] = 0; fCellNCuts[1] = 1;   fCellNCuts[2] = 2;   
  for(Int_t i = fNCellNCuts; i < 10; i++)fCellNCuts[i]  = 0;

  fNPIDBits = 2;
  fPIDBits[0] = 0;   fPIDBits[1] = 2; //  fPIDBits[2] = 4; fPIDBits[3] = 6;// check, no cut,  dispersion, neutral, dispersion&&neutral
  for(Int_t i = fNPIDBits; i < 10; i++)fPIDBits[i] = 0;

}


//________________________________________________________________________________________________________________________________________________
TObjString * AliAnaPi0::GetAnalysisCuts()
{  
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliAnaPi0 ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Number of bins in Centrality:  %d \n",fNCentrBin) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Number of bins in Z vert. pos: %d \n",GetNZvertBin()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Number of bins in Reac. Plain: %d \n",GetNRPBin()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Depth of event buffer: %d \n",fNmaxMixEv) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Pair in same Module: %d ; TrackMult as centrality: %d; PhotonMult as centrality: %d; cluster E as centrality: %d; cell as centrality: %d; Fill InvPt histos %d\n",
           fSameSM, fUseTrackMultBins, fUsePhotonMultBins, fUseAverClusterEBins, fUseAverCellEBins, fMakeInvPtPlots) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Select pairs with their angle: %d, edep %d, min angle %2.3f, max angle %2.3f,\n",fUseAngleCut, fUseAngleEDepCut,fAngleCut,fAngleMaxCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," Asymmetry cuts: n = %d, asymmetry < ",fNAsymCuts) ;
  for(Int_t i = 0; i < fNAsymCuts; i++) snprintf(onePar,buffersize,"%s %2.2f;",onePar,fAsymCuts[i]);
  parList+=onePar ;
  snprintf(onePar,buffersize," PID selection bits: n = %d, PID bit =\n",fNPIDBits) ;
  for(Int_t i = 0; i < fNPIDBits; i++) snprintf(onePar,buffersize,"%s %d;",onePar,fPIDBits[i]);
  parList+=onePar ;
  snprintf(onePar,buffersize,"Cuts: \n") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Z vertex position: -%f < z < %f \n",GetZvertexCut(),GetZvertexCut()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Calorimeter: %s \n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Number of modules: %d \n",fNModules) ;
  parList+=onePar ;
  if(fMultiCutAna){
    snprintf(onePar, buffersize," pT cuts: n = %d, pt > ",fNPtCuts) ;
    for(Int_t i = 0; i < fNPtCuts; i++) snprintf(onePar,buffersize,"%s %2.2f;",onePar,fPtCuts[i]);
    parList+=onePar ;
    snprintf(onePar,buffersize, " N cell in cluster cuts: n = %d, nCell > ",fNCellNCuts) ;
    for(Int_t i = 0; i < fNCellNCuts; i++) snprintf(onePar,buffersize,"%s %d;",onePar,fCellNCuts[i]);
    parList+=onePar ;
  }
  
  return new TObjString(parList) ;	
}

//________________________________________________________________________________________________________________________________________________
TList * AliAnaPi0::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  //create event containers
  fEventsList = new TList*[fNCentrBin*GetNZvertBin()*GetNRPBin()] ;
	
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t iz=0; iz<GetNZvertBin(); iz++){
      for(Int_t irp=0; irp<GetNRPBin(); irp++){
        fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp] = new TList() ;
        fEventsList[ic*GetNZvertBin()*GetNRPBin()+iz*GetNRPBin()+irp]->SetOwner(kFALSE);
      }
    }
  }
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName(GetName()); 
	
  fhReMod                = new TH2D*[fNModules]   ;
  fhMiMod                = new TH2D*[fNModules]   ;

  if(fCalorimeter == "PHOS"){
    fhReDiffPHOSMod        = new TH2D*[fNModules]   ;  
    fhMiDiffPHOSMod        = new TH2D*[fNModules]   ;
  }
  else{
    fhReSameSectorEMCALMod = new TH2D*[fNModules/2] ;
    fhReSameSideEMCALMod   = new TH2D*[fNModules-2] ;  
    fhMiSameSectorEMCALMod = new TH2D*[fNModules/2] ;
    fhMiSameSideEMCALMod   = new TH2D*[fNModules-2] ;
  }
  
  
  fhRe1 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  fhMi1 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(fFillBadDistHisto){
    fhRe2 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
    fhRe3 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
    fhMi2 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
    fhMi3 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  }
  if(fMakeInvPtPlots) {
    fhReInvPt1 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
    fhMiInvPt1 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
    if(fFillBadDistHisto){
      fhReInvPt2 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
      fhReInvPt3 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
      fhMiInvPt2 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
      fhMiInvPt3 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
    }
  } 
  
  const Int_t buffersize = 255;
  char key[buffersize] ;
  char title[buffersize] ;
  
  Int_t nptbins   = GetHistoPtBins();
  Int_t nphibins  = GetHistoPhiBins();
  Int_t netabins  = GetHistoEtaBins();
  Float_t ptmax   = GetHistoPtMax();
  Float_t phimax  = GetHistoPhiMax();
  Float_t etamax  = GetHistoEtaMax();
  Float_t ptmin   = GetHistoPtMin();
  Float_t phimin  = GetHistoPhiMin();
  Float_t etamin  = GetHistoEtaMin();	
	
  Int_t nmassbins = GetHistoMassBins();
  Int_t nasymbins = GetHistoAsymmetryBins();
  Float_t massmax = GetHistoMassMax();
  Float_t asymmax = GetHistoAsymmetryMax();
  Float_t massmin = GetHistoMassMin();
  Float_t asymmin = GetHistoAsymmetryMin();
  Int_t ntrmbins  = GetHistoTrackMultiplicityBins();
  Int_t ntrmmax   = GetHistoTrackMultiplicityMax();
  Int_t ntrmmin   = GetHistoTrackMultiplicityMin(); 
  
  fhAverTotECluster = new TH1F("hAverTotECluster","hAverTotECluster",200,0,50) ;
  fhAverTotECluster->SetXTitle("E_{cluster, aver. SM} (GeV)");
  outputContainer->Add(fhAverTotECluster) ;
  
  fhAverTotECell    = new TH1F("hAverTotECell","hAverTotECell",200,0,50) ;
  fhAverTotECell->SetXTitle("E_{cell, aver. SM} (GeV)");
  outputContainer->Add(fhAverTotECell) ;
  
  fhAverTotECellvsCluster    = new TH2F("hAverTotECellvsCluster","hAverTotECellvsCluster",200,0,50,200,0,50) ;
  fhAverTotECellvsCluster->SetYTitle("E_{cell, aver. SM} (GeV)");
  fhAverTotECellvsCluster->SetXTitle("E_{cluster, aver. SM} (GeV)");
  outputContainer->Add(fhAverTotECellvsCluster) ;
  
  fhEDensityCluster = new TH1F("hEDensityCluster","hEDensityCluster",200,0,50) ;
  fhEDensityCluster->SetXTitle("#Sigma E_{cluster} / N_{cluster} (GeV)");
  outputContainer->Add(fhEDensityCluster) ;
  
  fhEDensityCell    = new TH1F("hEDensityCell","hEDensityCell",200,0,50) ;
  fhEDensityCell->SetXTitle("#Sigma E_{cell} / N_{cell} (GeV)");
  outputContainer->Add(fhEDensityCell) ;
  
  fhEDensityCellvsCluster    = new TH2F("hEDensityCellvsCluster","hEDensityCellvsCluster",200,0,50,200,0,50) ;
  fhEDensityCellvsCluster->SetYTitle("#Sigma E_{cell} / N_{cell} (GeV)");
  fhEDensityCellvsCluster->SetXTitle("#Sigma E_{cluster} / N_{cluster} (GeV)");
  outputContainer->Add(fhEDensityCellvsCluster) ;
  
//  fhClusterPairDist = new TH1F("hClusterPairDist","Distance between clusters",250,0,750) ;
//  fhClusterPairDist->SetXTitle("#sqrt{(x_{1}-x_{2})^2+(z_{1}-z_{2})^2} (cm)");
//  outputContainer->Add(fhClusterPairDist) ;
//  
//  fhClusterPairDistWeight = new TH1F("hClusterPairDistWeighted","Distance between clusters, weighted by pair energy",200,0,400) ;
//  fhClusterPairDistWeight->SetXTitle("#sqrt{(x_{1}E_{1}-x_{2}E_{2})^{2}+(z_{1}E_{1}-z_{2}E_{2})^{2}}/ (E_{1}+E_{2}) (cm)");
//  outputContainer->Add(fhClusterPairDistWeight) ;
//   
//  fhAverClusterPairDist = new TH1F("hAverClusterPairDist","Average distance between clusters",250,0,750) ;
//  fhAverClusterPairDist->SetXTitle("#Sigma (#sqrt{(x_{1}-x_{2})^{2}+(z_{1}-z_{2})^{2}}) / N_{pairs} (cm)");
//  outputContainer->Add(fhAverClusterPairDist) ;
//  
//  fhAverClusterPairDistWeight = new TH1F("hAverClusterPairDistWeighted","Average distance between clusters, weighted by pair energy",200,0,400) ;
//  fhAverClusterPairDistWeight->SetXTitle("#Sigma (#sqrt{(x_{1}E_{1}-x_{2}E_{2})^{2}+(z_{1}E_{1}-z_{2}E_{2})^{2}}/ (E_{1}+E_{2})) / N_{pairs} (cm)");
//  outputContainer->Add(fhAverClusterPairDistWeight) ;
//  
//  fhAverClusterPairDistvsAverE = new TH2F("hAverClusterPairDistvsAverE","Average distance between clusters",250,0,750,200,0,50) ;
//  fhAverClusterPairDistvsAverE->SetXTitle("#Sigma (#sqrt{(x_{1}-x_{2})^{2}+(z_{1}-z_{2})^{2}}) / N_{pairs} (cm)");
//  fhAverClusterPairDistvsAverE->SetYTitle("#Sigma E_{cluster} / N_{cluster} (GeV)");
//  outputContainer->Add(fhAverClusterPairDistvsAverE) ;
//  
//  fhAverClusterPairDistWeightvsAverE = new TH2F("hAverClusterPairDistWeightedvsAverE","Average distance between clusters, weighted by pair energy",200,0,400,200,0,50) ;
//  fhAverClusterPairDistWeightvsAverE->SetXTitle("#Sigma (#sqrt{(x_{1}E_{1}-x_{2}E_{2})^2+(z_{1}E_{1}-z_{2}E_{2})^2}/ (E_{1}+E_{2})) / N_{pairs} (cm/GeV)");
//  fhAverClusterPairDistWeightvsAverE->SetYTitle("#Sigma E_{cluster} / N_{cluster} (GeV)");
//  outputContainer->Add(fhAverClusterPairDistWeightvsAverE) ;
  
//  fhAverClusterPairDistvsN = new TH2F("hAverClusterPairDistvsN","Average distance between clusters",250,0,750,200,0,50) ;
//  fhAverClusterPairDistvsN->SetXTitle("#Sigma (#sqrt{(x_{1}-x_{2})^{2}+(z_{1}-z_{2})^{2}}) / N_{pairs} (cm)");
//  fhAverClusterPairDistvsN->SetYTitle("N_{cluster}");
//  outputContainer->Add(fhAverClusterPairDistvsN) ;
//  
//  fhAverClusterPairDistWeightvsN = new TH2F("hAverClusterPairDistWeightedvsN","Average distance between clusters, weighted by pair energy",200,0,400,200,0,50) ;
//  fhAverClusterPairDistWeightvsN->SetXTitle("#Sigma (#sqrt{(x_{1}E_{1}-x_{2}E_{2})^{2}+(z_{1}E_{1}-z_{2}E_{2})^{2}}/ (E_{1}+E_{2})) / N_{pairs} (cm)");
//  fhAverClusterPairDistWeightvsN->SetYTitle("N_{cluster}");
//  outputContainer->Add(fhAverClusterPairDistWeightvsN) ;
  
//  fhMaxEvsClustMult = new TH2F("hMaxEvsClustMult","",nptbins,ptmin,ptmax,50,0,50) ;
//  fhMaxEvsClustMult->SetXTitle("E_{max}");
//  fhMaxEvsClustMult->SetYTitle("N_{cluster}");
//  outputContainer->Add(fhMaxEvsClustMult) ;
//  
//  fhMaxEvsClustEDen = new TH2F("hMaxEvsClustEDen","",nptbins,ptmin,ptmax,200,0,50) ;
//  fhMaxEvsClustEDen->SetXTitle("E_{max}");
//  fhMaxEvsClustEDen->SetYTitle("#Sigma E_{cluster} / N_{cluster} (GeV)");
//  outputContainer->Add(fhMaxEvsClustEDen) ;
  
  fhReConv = new TH2D("hReConv","Real Pair with one recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
  fhReConv->SetXTitle("p_{T} (GeV/c)");
  fhReConv->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
  outputContainer->Add(fhReConv) ;
  
  fhReConv2 = new TH2D("hReConv2","Real Pair with 2 recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
  fhReConv2->SetXTitle("p_{T} (GeV/c)");
  fhReConv2->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
  outputContainer->Add(fhReConv2) ;
  
  if(fDoOwnMix){
    fhMiConv = new TH2D("hMiConv","Mixed Pair with one recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhMiConv->SetXTitle("p_{T} (GeV/c)");
    fhMiConv->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
    outputContainer->Add(fhMiConv) ;
  
    fhMiConv2 = new TH2D("hMiConv2","Mixed Pair with 2 recombined conversion ",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhMiConv2->SetXTitle("p_{T} (GeV/c)");
    fhMiConv2->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
    outputContainer->Add(fhMiConv2) ;
  }
  
  for(Int_t ic=0; ic<fNCentrBin; ic++){
      for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
        for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
          Int_t index = ((ic*fNPIDBits)+ipid)*fNAsymCuts + iasym;
          //printf("cen %d, pid %d, asy %d, Index %d\n",ic,ipid,iasym,index);
          //Distance to bad module 1
          snprintf(key, buffersize,"hRe_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                   ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
          fhRe1[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRe1[index]->SetXTitle("p_{T} (GeV/c)");
          fhRe1[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          //printf("name: %s\n ",fhRe1[index]->GetName());
          outputContainer->Add(fhRe1[index]) ;
          
          if(fFillBadDistHisto){
            //Distance to bad module 2
            snprintf(key, buffersize,"hRe_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhRe2[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhRe2[index]->SetXTitle("p_{T} (GeV/c)");
            fhRe2[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhRe2[index]) ;
            
            //Distance to bad module 3
            snprintf(key, buffersize,"hRe_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 3",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhRe3[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhRe3[index]->SetXTitle("p_{T} (GeV/c)");
            fhRe3[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhRe3[index]) ;
          }
          
          //Inverse pT 
          if(fMakeInvPtPlots){
            //Distance to bad module 1
            snprintf(key, buffersize,"hReInvPt_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhReInvPt1[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhReInvPt1[index]->SetXTitle("p_{T} (GeV/c)");
            fhReInvPt1[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhReInvPt1[index]) ;
            
            if(fFillBadDistHisto){
              //Distance to bad module 2
              snprintf(key, buffersize,"hReInvPt_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhReInvPt2[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhReInvPt2[index]->SetXTitle("p_{T} (GeV/c)");
              fhReInvPt2[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
              outputContainer->Add(fhReInvPt2[index]) ;
              
              //Distance to bad module 3
              snprintf(key, buffersize,"hReInvPt_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 3",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhReInvPt3[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhReInvPt3[index]->SetXTitle("p_{T} (GeV/c)");
              fhReInvPt3[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
              outputContainer->Add(fhReInvPt3[index]) ;
            }
          }
          if(fDoOwnMix){
            //Distance to bad module 1
            snprintf(key, buffersize,"hMi_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
            snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                     ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
            fhMi1[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMi1[index]->SetXTitle("p_{T} (GeV/c)");
            fhMi1[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhMi1[index]) ;
            if(fFillBadDistHisto){
              //Distance to bad module 2
              snprintf(key, buffersize,"hMi_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhMi2[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMi2[index]->SetXTitle("p_{T} (GeV/c)");
              fhMi2[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
              outputContainer->Add(fhMi2[index]) ;
              
              //Distance to bad module 3
              snprintf(key, buffersize,"hMi_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 3",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhMi3[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMi3[index]->SetXTitle("p_{T} (GeV/c)");
              fhMi3[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
              outputContainer->Add(fhMi3[index]) ;
            }
            //Inverse pT
            if(fMakeInvPtPlots){
              //Distance to bad module 1
              snprintf(key, buffersize,"hMiInvPt_cen%d_pidbit%d_asy%d_dist1",ic,ipid,iasym) ;
              snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 1",
                       ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
              fhMiInvPt1[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
              fhMiInvPt1[index]->SetXTitle("p_{T} (GeV/c)");
              fhMiInvPt1[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
              outputContainer->Add(fhMiInvPt1[index]) ;
              if(fFillBadDistHisto){
                //Distance to bad module 2
                snprintf(key, buffersize,"hMiInvPt_cen%d_pidbit%d_asy%d_dist2",ic,ipid,iasym) ;
                snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f, dist bad 2",
                         ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
                fhMiInvPt2[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
                fhMiInvPt2[index]->SetXTitle("p_{T} (GeV/c)");
                fhMiInvPt2[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
                outputContainer->Add(fhMiInvPt2[index]) ;
                
                //Distance to bad module 3
                snprintf(key, buffersize,"hMiInvPt_cen%d_pidbit%d_asy%d_dist3",ic,ipid,iasym) ;
                snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for centrality=%d, PID bit=%d and asymmetry %1.2f,dist bad 3",
                         ic,fPIDBits[ipid], fAsymCuts[iasym]) ;
                fhMiInvPt3[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
                fhMiInvPt3[index]->SetXTitle("p_{T} (GeV/c)");
                fhMiInvPt3[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
                outputContainer->Add(fhMiInvPt3[index]) ;
              }
            }
          } 
        }
      }
    }
  
  fhRePtAsym = new TH2D("hRePtAsym","Asymmetry vs pt, for pairs",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
  fhRePtAsym->SetXTitle("p_{T} (GeV/c)");
  fhRePtAsym->SetYTitle("Asymmetry");
  outputContainer->Add(fhRePtAsym);
  
  fhRePtAsymPi0 = new TH2D("hRePtAsymPi0","Asymmetry vs pt, for pairs close to #pi^{0} mass",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
  fhRePtAsymPi0->SetXTitle("p_{T} (GeV/c)");
  fhRePtAsymPi0->SetYTitle("Asymmetry");
  outputContainer->Add(fhRePtAsymPi0);

  fhRePtAsymEta = new TH2D("hRePtAsymEta","Asymmetry vs pt, for pairs close to #eta mass",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
  fhRePtAsymEta->SetXTitle("p_{T} (GeV/c)");
  fhRePtAsymEta->SetYTitle("Asymmetry");
  outputContainer->Add(fhRePtAsymEta);
  
  if(fMultiCutAna){
    
    fhRePIDBits         = new TH2D*[fNPIDBits];
    for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
      snprintf(key,   buffersize,"hRe_pidbit%d",ipid) ;
      snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for PIDBit=%d",fPIDBits[ipid]) ;
      fhRePIDBits[ipid] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhRePIDBits[ipid]->SetXTitle("p_{T} (GeV/c)");
      fhRePIDBits[ipid]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhRePIDBits[ipid]) ;
    }// pid bit loop
    
    fhRePtNCellAsymCuts    = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    fhRePtNCellAsymCutsSM0 = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    fhRePtNCellAsymCutsSM1 = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    fhRePtNCellAsymCutsSM2 = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    fhRePtNCellAsymCutsSM3 = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    fhMiPtNCellAsymCuts    = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    for(Int_t ipt=0; ipt<fNPtCuts; ipt++){
      for(Int_t icell=0; icell<fNCellNCuts; icell++){
        for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%1.2f ",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
          //printf("ipt %d, icell %d, iassym %d, index %d\n",ipt, icell, iasym, index);
          fhRePtNCellAsymCuts[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRePtNCellAsymCuts[index]->SetXTitle("p_{T} (GeV/c)");
          fhRePtNCellAsymCuts[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhRePtNCellAsymCuts[index]) ;
                    
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d_SM0",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%1.2f, SM 0 ",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          fhRePtNCellAsymCutsSM0[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRePtNCellAsymCutsSM0[index]->SetXTitle("p_{T} (GeV/c)");
          fhRePtNCellAsymCutsSM0[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhRePtNCellAsymCutsSM0[index]) ;
          
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d_SM1",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%1.2f, SM 1 ",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          fhRePtNCellAsymCutsSM1[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRePtNCellAsymCutsSM1[index]->SetXTitle("p_{T} (GeV/c)");
          fhRePtNCellAsymCutsSM1[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhRePtNCellAsymCutsSM1[index]) ;
          
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d_SM2",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%1.2f, SM 2 ",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          fhRePtNCellAsymCutsSM2[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRePtNCellAsymCutsSM2[index]->SetXTitle("p_{T} (GeV/c)");
          fhRePtNCellAsymCutsSM2[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhRePtNCellAsymCutsSM2[index]) ;
          
          snprintf(key,   buffersize,"hRe_pt%d_cell%d_asym%d_SM3",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%1.2f, SM 3 ",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          fhRePtNCellAsymCutsSM3[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhRePtNCellAsymCutsSM3[index]->SetXTitle("p_{T} (GeV/c)");
          fhRePtNCellAsymCutsSM3[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhRePtNCellAsymCutsSM3[index]) ;
          
          snprintf(key,   buffersize,"hMi_pt%d_cell%d_asym%d",ipt,icell,iasym) ;
          snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]) ;
          fhMiPtNCellAsymCuts[index] = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMiPtNCellAsymCuts[index]->SetXTitle("p_{T} (GeV/c)");
          fhMiPtNCellAsymCuts[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhMiPtNCellAsymCuts[index]) ;
          
        }
      }
    }
    
    fhRePtMult = new TH3D*[fNAsymCuts] ;
    for(Int_t iasym = 0; iasym<fNAsymCuts; iasym++){
      fhRePtMult[iasym] = new TH3D(Form("hRePtMult_asym%d",iasym),Form("(p_{T},C,M)_{#gamma#gamma}, A<%1.2f",fAsymCuts[iasym]),
                                   nptbins,ptmin,ptmax,ntrmbins,ntrmmin,ntrmmax,nmassbins,massmin,massmax);
      fhRePtMult[iasym]->SetXTitle("p_{T} (GeV/c)");
      fhRePtMult[iasym]->SetYTitle("Track multiplicity");
      fhRePtMult[iasym]->SetZTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhRePtMult[iasym]) ;
    }
    
  }// multi cuts analysis
  
  fhEvents=new TH3D("hEvents","Number of events",fNCentrBin,0.,1.*fNCentrBin,
                    GetNZvertBin(),0.,1.*GetNZvertBin(),GetNRPBin(),0.,1.*GetNRPBin()) ;
  
  fhEvents->SetXTitle("Centrality bin");
  fhEvents->SetYTitle("Z vertex bin bin");
  fhEvents->SetZTitle("RP bin");
  outputContainer->Add(fhEvents) ;
	
  fhCentrality=new TH1D("hCentralityBin","Number of events in centrality bin",fNCentrBin*10,0.,1.*fNCentrBin) ;
  fhCentrality->SetXTitle("Centrality bin");
  outputContainer->Add(fhCentrality) ;

  fhRealOpeningAngle  = new TH2D
  ("hRealOpeningAngle","Angle between all #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,300,0,TMath::Pi()); 
  fhRealOpeningAngle->SetYTitle("#theta(rad)");
  fhRealOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
  outputContainer->Add(fhRealOpeningAngle) ;
  
  fhRealCosOpeningAngle  = new TH2D
  ("hRealCosOpeningAngle","Cosinus of angle between all #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,100,0,1); 
  fhRealCosOpeningAngle->SetYTitle("cos (#theta) ");
  fhRealCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
  outputContainer->Add(fhRealCosOpeningAngle) ;
	
  if(fDoOwnMix){
    
    fhMixedOpeningAngle  = new TH2D
    ("hMixedOpeningAngle","Angle between all #gamma pair vs E_{#pi^{0}}, Mixed pairs",nptbins,ptmin,ptmax,300,0,TMath::Pi()); 
    fhMixedOpeningAngle->SetYTitle("#theta(rad)");
    fhMixedOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhMixedOpeningAngle) ;
    
    fhMixedCosOpeningAngle  = new TH2D
    ("hMixedCosOpeningAngle","Cosinus of angle between all #gamma pair vs E_{#pi^{0}}, Mixed pairs",nptbins,ptmin,ptmax,100,0,1); 
    fhMixedCosOpeningAngle->SetYTitle("cos (#theta) ");
    fhMixedCosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhMixedCosOpeningAngle) ;
    
  }
  
  //Histograms filled only if MC data is requested 	
  if(IsDataMC()){
    //Pi0
    fhPrimPi0Pt     = new TH1D("hPrimPi0Pt","Primary pi0 pt",nptbins,ptmin,ptmax) ;
    fhPrimPi0AccPt  = new TH1D("hPrimPi0AccPt","Primary pi0 pt with both photons in acceptance",nptbins,ptmin,ptmax) ;
    fhPrimPi0Pt   ->SetXTitle("p_{T} (GeV/c)");
    fhPrimPi0AccPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimPi0Pt) ;
    outputContainer->Add(fhPrimPi0AccPt) ;
    
    fhPrimPi0Y      = new TH2D("hPrimPi0Rapidity","Rapidity of primary pi0",nptbins,ptmin,ptmax, netabins,etamin,etamax) ;
    fhPrimPi0Y   ->SetYTitle("Rapidity");
    fhPrimPi0Y   ->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimPi0Y) ;
    
    fhPrimPi0AccY   = new TH2D("hPrimPi0AccRapidity","Rapidity of primary pi0",nptbins,ptmin,ptmax, netabins,etamin,etamax) ; 
    fhPrimPi0AccY->SetYTitle("Rapidity");
    fhPrimPi0AccY->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimPi0AccY) ;
    
    fhPrimPi0Phi    = new TH2D("hPrimPi0Phi","Azimuthal of primary pi0",nptbins,ptmin,ptmax, nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ; 
    fhPrimPi0Phi->SetYTitle("#phi (deg)");
    fhPrimPi0Phi->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimPi0Phi) ;
    
    fhPrimPi0AccPhi = new TH2D("hPrimPi0AccPhi","Azimuthal of primary pi0 with accepted daughters",nptbins,ptmin,ptmax, nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ; 
    fhPrimPi0AccPhi->SetYTitle("#phi (deg)");
    fhPrimPi0AccPhi->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimPi0AccPhi) ;
    
    //Eta
    fhPrimEtaPt     = new TH1D("hPrimEtaPt","Primary eta pt",nptbins,ptmin,ptmax) ;
    fhPrimEtaAccPt  = new TH1D("hPrimEtaAccPt","Primary eta pt with both photons in acceptance",nptbins,ptmin,ptmax) ;
    fhPrimEtaPt   ->SetXTitle("p_{T} (GeV/c)");
    fhPrimEtaAccPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimEtaPt) ;
    outputContainer->Add(fhPrimEtaAccPt) ;
    
    fhPrimEtaY      = new TH2D("hPrimEtaRapidity","Rapidity of primary eta",nptbins,ptmin,ptmax,netabins,etamin,etamax) ; 
    fhPrimEtaY->SetYTitle("Rapidity");
    fhPrimEtaY->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimEtaY) ;
    
    fhPrimEtaAccY   = new TH2D("hPrimEtaAccRapidity","Rapidity of primary eta",nptbins,ptmin,ptmax, netabins,etamin,etamax) ; 
    fhPrimEtaAccY->SetYTitle("Rapidity");
    fhPrimEtaAccY->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimEtaAccY) ;
    
    fhPrimEtaPhi    = new TH2D("hPrimEtaPhi","Azimuthal of primary eta",nptbins,ptmin,ptmax, nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ; 
    fhPrimEtaPhi->SetYTitle("#phi (deg)");
    fhPrimEtaPhi->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimEtaPhi) ;
    
    fhPrimEtaAccPhi = new TH2D("hPrimEtaAccPhi","Azimuthal of primary eta with accepted daughters",nptbins,ptmin,ptmax, nphibins,phimin*TMath::RadToDeg(),phimax*TMath::RadToDeg()) ; 
    fhPrimEtaAccPhi->SetYTitle("#phi (deg)");
    fhPrimEtaAccPhi->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPrimEtaAccPhi) ;
        
    
    //Prim origin
    //Pi0
    fhPrimPi0PtOrigin     = new TH2D("hPrimPi0PtOrigin","Primary pi0 pt vs origin",nptbins,ptmin,ptmax,11,0,11) ;
    fhPrimPi0PtOrigin->SetXTitle("p_{T} (GeV/c)");
    fhPrimPi0PtOrigin->SetYTitle("Origin");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances ");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
    fhPrimPi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
    outputContainer->Add(fhPrimPi0PtOrigin) ;
    
    fhMCPi0PtOrigin     = new TH2D("hMCPi0PtOrigin","Reconstructed pair from generated pi0 pt vs origin",nptbins,ptmin,ptmax,11,0,11) ;
    fhMCPi0PtOrigin->SetXTitle("p_{T} (GeV/c)");
    fhMCPi0PtOrigin->SetYTitle("Origin");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(5 ,"#rho");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(6 ,"#omega");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(7 ,"K");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(8 ,"Other");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(9 ,"#eta");
    fhMCPi0PtOrigin->GetYaxis()->SetBinLabel(10 ,"#eta prime");
    outputContainer->Add(fhMCPi0PtOrigin) ;    
    
    //Eta
    fhPrimEtaPtOrigin     = new TH2D("hPrimEtaPtOrigin","Primary pi0 pt vs origin",nptbins,ptmin,ptmax,7,0,7) ;
    fhPrimEtaPtOrigin->SetXTitle("p_{T} (GeV/c)");
    fhPrimEtaPtOrigin->SetYTitle("Origin");
    fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
    fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(5 ,"Other");
    fhPrimEtaPtOrigin->GetYaxis()->SetBinLabel(6 ,"#eta prime ");

    outputContainer->Add(fhPrimEtaPtOrigin) ;
    
    fhMCEtaPtOrigin     = new TH2D("hMCEtaPtOrigin","Reconstructed pair from generated pi0 pt vs origin",nptbins,ptmin,ptmax,7,0,7) ;
    fhMCEtaPtOrigin->SetXTitle("p_{T} (GeV/c)");
    fhMCEtaPtOrigin->SetYTitle("Origin");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(1 ,"Status 21");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(2 ,"Quark");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(3 ,"qq Resonances");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(4 ,"Resonances");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(5 ,"Other");
    fhMCEtaPtOrigin->GetYaxis()->SetBinLabel(6 ,"#eta prime");

    outputContainer->Add(fhMCEtaPtOrigin) ;
     
    
    fhPrimPi0OpeningAngle  = new TH2D
    ("hPrimPi0OpeningAngle","Angle between all primary #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,100,0,0.5); 
    fhPrimPi0OpeningAngle->SetYTitle("#theta(rad)");
    fhPrimPi0OpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhPrimPi0OpeningAngle) ;
    
    fhPrimPi0CosOpeningAngle  = new TH2D
    ("hPrimPi0CosOpeningAngle","Cosinus of angle between all primary #gamma pair vs E_{#pi^{0}}",nptbins,ptmin,ptmax,100,-1,1); 
    fhPrimPi0CosOpeningAngle->SetYTitle("cos (#theta) ");
    fhPrimPi0CosOpeningAngle->SetXTitle("E_{ #pi^{0}} (GeV)");
    outputContainer->Add(fhPrimPi0CosOpeningAngle) ;
    
    for(Int_t i = 0; i<13; i++){
      fhMCOrgMass[i] = new TH2D(Form("hMCOrgMass_%d",i),Form("mass vs pt, origin %d",i),nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMCOrgMass[i]->SetXTitle("p_{T} (GeV/c)");
      fhMCOrgMass[i]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhMCOrgMass[i]) ;
      
      fhMCOrgAsym[i]= new TH2D(Form("hMCOrgAsym_%d",i),Form("asymmetry vs pt, origin %d",i),nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax) ;
      fhMCOrgAsym[i]->SetXTitle("p_{T} (GeV/c)");
      fhMCOrgAsym[i]->SetYTitle("A");
      outputContainer->Add(fhMCOrgAsym[i]) ;
      
      fhMCOrgDeltaEta[i] = new TH2D(Form("hMCOrgDeltaEta_%d",i),Form("#Delta #eta of pair vs pt, origin %d",i),nptbins,ptmin,ptmax,netabins,-1.4,1.4) ;
      fhMCOrgDeltaEta[i]->SetXTitle("p_{T} (GeV/c)");
      fhMCOrgDeltaEta[i]->SetYTitle("#Delta #eta");
      outputContainer->Add(fhMCOrgDeltaEta[i]) ;
      
      fhMCOrgDeltaPhi[i]= new TH2D(Form("hMCOrgDeltaPhi_%d",i),Form("#Delta #phi of pair vs p_{T}, origin %d",i),nptbins,ptmin,ptmax,nphibins,-0.7,0.7) ;
      fhMCOrgDeltaPhi[i]->SetXTitle("p_{T} (GeV/c)");
      fhMCOrgDeltaPhi[i]->SetYTitle("#Delta #phi (rad)");
      outputContainer->Add(fhMCOrgDeltaPhi[i]) ;
      
    }
    
    if(fMultiCutAnaSim){
      fhMCPi0MassPtTrue  = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCPi0MassPtRec   = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCPi0PtTruePtRec = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCEtaMassPtRec   = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCEtaMassPtTrue  = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCEtaPtTruePtRec = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      for(Int_t ipt=0; ipt<fNPtCuts; ipt++){
        for(Int_t icell=0; icell<fNCellNCuts; icell++){
          for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
            Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
            
            fhMCPi0MassPtRec[index] = new TH2D(Form("hMCPi0MassPtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                               Form("Reconstructed Mass vs reconstructed p_T of true #pi^{0} cluster pairs for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                               nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMCPi0MassPtRec[index]->SetXTitle("p_{T, reconstructed} (GeV/c)");
            fhMCPi0MassPtRec[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhMCPi0MassPtRec[index]) ;    
            
            fhMCPi0MassPtTrue[index] = new TH2D(Form("hMCPi0MassPtTrue_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                Form("Reconstructed Mass vs generated p_T of true #pi^{0} cluster pairs for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMCPi0MassPtTrue[index]->SetXTitle("p_{T, generated} (GeV/c)");
            fhMCPi0MassPtTrue[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhMCPi0MassPtTrue[index]) ;
            
            fhMCPi0PtTruePtRec[index] = new TH2D(Form("hMCPi0PtTruePtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                 Form("Generated vs reconstructed p_T of true #pi^{0} cluster pairs, 0.01 < rec. mass < 0.17 MeV/c^{2} for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                 nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
            fhMCPi0PtTruePtRec[index]->SetXTitle("p_{T, generated} (GeV/c)");
            fhMCPi0PtTruePtRec[index]->SetYTitle("p_{T, reconstructed} (GeV/c)");
            outputContainer->Add(fhMCPi0PtTruePtRec[index]) ;
            
            fhMCEtaMassPtRec[index] = new TH2D(Form("hMCEtaMassPtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                               Form("Reconstructed Mass vs reconstructed p_T of true #eta cluster pairs for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                               nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMCEtaMassPtRec[index]->SetXTitle("p_{T, generated} (GeV/c)");
            fhMCEtaMassPtRec[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhMCEtaMassPtRec[index]) ;
            
            fhMCEtaMassPtTrue[index] = new TH2D(Form("hMCEtaMassPtTrue_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                Form("Reconstructed Mass vs generated p_T of true #eta cluster pairs for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
            fhMCEtaMassPtTrue[index]->SetXTitle("p_{T, generated} (GeV/c)");
            fhMCEtaMassPtTrue[index]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
            outputContainer->Add(fhMCEtaMassPtTrue[index]) ;
            
            fhMCEtaPtTruePtRec[index] = new TH2D(Form("hMCEtaPtTruePtRec_pt%d_cell%d_asym%d",ipt,icell,iasym),
                                                 Form("Generated vs reconstructed p_T of true #eta cluster pairs, 0.01 < rec. mass < 0.17 MeV/c^{2} for pt >%2.2f, ncell>%d and asym >%1.2f",fPtCuts[ipt],fCellNCuts[icell], fAsymCuts[iasym]),
                                                 nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
            fhMCEtaPtTruePtRec[index]->SetXTitle("p_{T, generated} (GeV/c)");
            fhMCEtaPtTruePtRec[index]->SetYTitle("p_{T, reconstructed} (GeV/c)");
            outputContainer->Add(fhMCEtaPtTruePtRec[index]) ;
          }
        }
      }  
    }//multi cut ana
    else {
      fhMCPi0MassPtTrue  = new TH2D*[1];
      fhMCPi0PtTruePtRec = new TH2D*[1];
      fhMCEtaMassPtTrue  = new TH2D*[1];
      fhMCEtaPtTruePtRec = new TH2D*[1];
      
      fhMCPi0MassPtTrue[0] = new TH2D("hMCPi0MassPtTrue","Reconstructed Mass vs generated p_T of true #pi^{0} cluster pairs",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMCPi0MassPtTrue[0]->SetXTitle("p_{T, generated} (GeV/c)");
      fhMCPi0MassPtTrue[0]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhMCPi0MassPtTrue[0]) ;
      
      fhMCPi0PtTruePtRec[0]= new TH2D("hMCPi0PtTruePtRec","Generated vs reconstructed p_T of true #pi^{0} cluster pairs, 0.01 < rec. mass < 0.17 MeV/c^{2}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
      fhMCPi0PtTruePtRec[0]->SetXTitle("p_{T, generated} (GeV/c)");
      fhMCPi0PtTruePtRec[0]->SetYTitle("p_{T, reconstructed} (GeV/c)");
      outputContainer->Add(fhMCPi0PtTruePtRec[0]) ;
      
      fhMCEtaMassPtTrue[0] = new TH2D("hMCEtaMassPtTrue","Reconstructed Mass vs generated p_T of true #eta cluster pairs",nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMCEtaMassPtTrue[0]->SetXTitle("p_{T, generated} (GeV/c)");
      fhMCEtaMassPtTrue[0]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhMCEtaMassPtTrue[0]) ;
      
      fhMCEtaPtTruePtRec[0]= new TH2D("hMCEtaPtTruePtRec","Generated vs reconstructed p_T of true #eta cluster pairs, 0.01 < rec. mass < 0.17 MeV/c^{2}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax) ;
      fhMCEtaPtTruePtRec[0]->SetXTitle("p_{T, generated} (GeV/c)");
      fhMCEtaPtTruePtRec[0]->SetYTitle("p_{T, reconstructed} (GeV/c)");
      outputContainer->Add(fhMCEtaPtTruePtRec[0]) ;
    }
  }
  
  TString pairnamePHOS[] = {"(0-1)","(0-2)","(1-2)","(0-3)","(0-4)","(1-3)","(1-4)","(2-3)","(2-4)","(3-4)"};
  for(Int_t imod=0; imod<fNModules; imod++){
    //Module dependent invariant mass
    snprintf(key, buffersize,"hReMod_%d",imod) ;
    snprintf(title, buffersize,"Real m_{#gamma#gamma} distr. for Module %d",imod) ;
    fhReMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
    fhReMod[imod]->SetXTitle("p_{T} (GeV/c)");
    fhReMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
    outputContainer->Add(fhReMod[imod]) ;
    if(fCalorimeter=="PHOS"){
      snprintf(key, buffersize,"hReDiffPHOSMod_%d",imod) ;
      snprintf(title, buffersize,"Real pairs PHOS, clusters in different Modules: %s",(pairnamePHOS[imod]).Data()) ;
      fhReDiffPHOSMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhReDiffPHOSMod[imod]->SetXTitle("p_{T} (GeV/c)");
      fhReDiffPHOSMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhReDiffPHOSMod[imod]) ;
    }
    else{//EMCAL
      if(imod<fNModules/2){
        snprintf(key, buffersize,"hReSameSectorEMCAL_%d",imod) ;
        snprintf(title, buffersize,"Real pairs EMCAL, clusters in same sector, SM(%d,%d)",imod*2,imod*2+1) ;
        fhReSameSectorEMCALMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhReSameSectorEMCALMod[imod]->SetXTitle("p_{T} (GeV/c)");
        fhReSameSectorEMCALMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
        outputContainer->Add(fhReSameSectorEMCALMod[imod]) ;
      }
      if(imod<fNModules-2){
        snprintf(key, buffersize,"hReSameSideEMCAL_%d",imod) ;
        snprintf(title, buffersize,"Real pairs EMCAL, clusters in same side SM(%d,%d)",imod, imod+2) ;
        fhReSameSideEMCALMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhReSameSideEMCALMod[imod]->SetXTitle("p_{T} (GeV/c)");
        fhReSameSideEMCALMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
        outputContainer->Add(fhReSameSideEMCALMod[imod]) ;
      }
    }//EMCAL
    
    if(fDoOwnMix){ 
      snprintf(key, buffersize,"hMiMod_%d",imod) ;
      snprintf(title, buffersize,"Mixed m_{#gamma#gamma} distr. for Module %d",imod) ;
      fhMiMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
      fhMiMod[imod]->SetXTitle("p_{T} (GeV/c)");
      fhMiMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
      outputContainer->Add(fhMiMod[imod]) ;
      
      if(fCalorimeter=="PHOS"){
        snprintf(key, buffersize,"hMiDiffPHOSMod_%d",imod) ;
        snprintf(title, buffersize,"Mixed pairs PHOS, clusters in different Modules: %s",(pairnamePHOS[imod]).Data()) ;
        fhMiDiffPHOSMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
        fhMiDiffPHOSMod[imod]->SetXTitle("p_{T} (GeV/c)");
        fhMiDiffPHOSMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
        outputContainer->Add(fhMiDiffPHOSMod[imod]) ;
      }//PHOS
      else{//EMCAL
        if(imod<fNModules/2){
          snprintf(key, buffersize,"hMiSameSectorEMCALMod_%d",imod) ;
          snprintf(title, buffersize,"Mixed pairs EMCAL, clusters in same sector, SM(%d,%d)",imod*2,imod*2+1) ;
          fhMiSameSectorEMCALMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMiSameSectorEMCALMod[imod]->SetXTitle("p_{T} (GeV/c)");
          fhMiSameSectorEMCALMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhMiSameSectorEMCALMod[imod]) ;
        }
        if(imod<fNModules-2){
          snprintf(key, buffersize,"hMiSameSideEMCALMod_%d",imod) ;
          snprintf(title, buffersize,"Mixed pairs EMCAL, clusters in same side SM(%d,%d)",imod, imod+2) ;
          fhMiSameSideEMCALMod[imod]  = new TH2D(key,title,nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
          fhMiSameSideEMCALMod[imod]->SetXTitle("p_{T} (GeV/c)");
          fhMiSameSideEMCALMod[imod]->SetYTitle("m_{#gamma,#gamma} (GeV/c^{2})");
          outputContainer->Add(fhMiSameSideEMCALMod[imod]) ;
        }
      }//EMCAL      
     }
  }
      
//  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++){
//  
//    printf("Histogram %d, name: %s\n ",i, outputContainer->At(i)->GetName());
//  
//  }
  
  return outputContainer;
}

//_________________________________________________________________________________________________________________________________________________
void AliAnaPi0::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set for the analysis
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Number of bins in Centrality:  %d \n",fNCentrBin) ;
  printf("Number of bins in Z vert. pos: %d \n",GetNZvertBin()) ;
  printf("Number of bins in Reac. Plain: %d \n",GetNRPBin()) ;
  printf("Depth of event buffer: %d \n",fNmaxMixEv) ;
  printf("Pair in same Module: %d \n",fSameSM) ;
  printf("Cuts: \n") ;
// printf("Z vertex position: -%2.3f < z < %2.3f \n",GetZvertexCut(),GetZvertexCut()) ; //It crashes here, why?
  printf("Number of modules:             %d \n",fNModules) ;
  printf("Select pairs with their angle: %d, edep %d, min angle %2.3f, max angle %2.3f \n",fUseAngleCut, fUseAngleEDepCut, fAngleCut, fAngleMaxCut) ;
  printf("Asymmetry cuts: n = %d, \n",fNAsymCuts) ;
  printf("\tasymmetry < ");
  for(Int_t i = 0; i < fNAsymCuts; i++) printf("%2.2f ",fAsymCuts[i]);
  printf("\n");
  
  printf("PID selection bits: n = %d, \n",fNPIDBits) ;
  printf("\tPID bit = ");
  for(Int_t i = 0; i < fNPIDBits; i++) printf("%d ",fPIDBits[i]);
  printf("\n");
  
  if(fMultiCutAna){
    printf("pT cuts: n = %d, \n",fNPtCuts) ;
    printf("\tpT > ");
    for(Int_t i = 0; i < fNPtCuts; i++) printf("%2.2f ",fPtCuts[i]);
    printf("GeV/c\n");
    
    printf("N cell in cluster cuts: n = %d, \n",fNCellNCuts) ;
    printf("\tnCell > ");
    for(Int_t i = 0; i < fNCellNCuts; i++) printf("%d ",fCellNCuts[i]);
    printf("\n");

  }
  printf("------------------------------------------------------\n") ;
} 

//_____________________________________________________________
void AliAnaPi0::FillAcceptanceHistograms(){
  //Fill acceptance histograms if MC data is available
  
  if(GetReader()->ReadStack()){	
    AliStack * stack = GetMCStack();
    if(stack){
      for(Int_t i=0 ; i<stack->GetNprimary(); i++){
        TParticle * prim = stack->Particle(i) ;
        Int_t pdg = prim->GetPdgCode();
        if( pdg == 111 || pdg == 221){
          Double_t pi0Pt = prim->Pt() ;
          //printf("pi0, pt %2.2f\n",pi0Pt);
          if(prim->Energy() == TMath::Abs(prim->Pz()))  continue ; //Protection against floating point exception	  
          Double_t pi0Y  = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
          Double_t phi   = TMath::RadToDeg()*prim->Phi() ;
          if(pdg == 111){
            if(TMath::Abs(pi0Y) < 1.0){
              fhPrimPi0Pt->Fill(pi0Pt) ;
            }
            fhPrimPi0Y  ->Fill(pi0Pt, pi0Y) ;
            fhPrimPi0Phi->Fill(pi0Pt, phi) ;
          }
          else if(pdg == 221){
            if(TMath::Abs(pi0Y) < 1.0){
              fhPrimEtaPt->Fill(pi0Pt) ;
            }
            fhPrimEtaY  ->Fill(pi0Pt, pi0Y) ;
            fhPrimEtaPhi->Fill(pi0Pt, phi) ;
          }
          
          //Origin of meson
          Int_t momindex  = prim->GetFirstMother();
          if(momindex < 0) continue;
          TParticle* mother = stack->Particle(momindex);
          Int_t mompdg    = TMath::Abs(mother->GetPdgCode());
          Int_t momstatus = mother->GetStatusCode();
          if(pdg == 111){
            if     (momstatus  == 21)fhPrimPi0PtOrigin->Fill(pi0Pt,0.5);//parton
            else if(mompdg     < 22 ) fhPrimPi0PtOrigin->Fill(pi0Pt,1.5);//quark
            else if(mompdg     > 2100  && mompdg   < 2210) fhPrimPi0PtOrigin->Fill(pi0Pt,2.5);// resonances
            else if(mompdg    == 221) fhPrimPi0PtOrigin->Fill(pi0Pt,8.5);//eta
            else if(mompdg    == 331) fhPrimPi0PtOrigin->Fill(pi0Pt,9.5);//eta prime
            else if(mompdg    == 213) fhPrimPi0PtOrigin->Fill(pi0Pt,4.5);//rho
            else if(mompdg    == 223) fhPrimPi0PtOrigin->Fill(pi0Pt,5.5);//omega
            else if(mompdg    >= 310   && mompdg    <= 323) fhPrimPi0PtOrigin->Fill(pi0Pt,6.5);//k0S, k+-,k*
            else if(mompdg    == 130) fhPrimPi0PtOrigin->Fill(pi0Pt,6.5);//k0L
            else if(momstatus == 11 || momstatus  == 12 ) fhPrimPi0PtOrigin->Fill(pi0Pt,3.5);//resonances   
            else                      fhPrimPi0PtOrigin->Fill(pi0Pt,7.5);//other?
          }//pi0
          else {
            if     (momstatus == 21 ) fhPrimEtaPtOrigin->Fill(pi0Pt,0.5);//parton
            else if(mompdg    < 22  ) fhPrimEtaPtOrigin->Fill(pi0Pt,1.5);//quark
            else if(mompdg    > 2100  && mompdg   < 2210) fhPrimEtaPtOrigin->Fill(pi0Pt,2.5);//qq resonances
            else if(mompdg    == 331) fhPrimEtaPtOrigin->Fill(pi0Pt,5.5);//eta prime
            else if(momstatus == 11 || momstatus  == 12 ) fhPrimEtaPtOrigin->Fill(pi0Pt,3.5);//resonances
            else fhPrimEtaPtOrigin->Fill(pi0Pt,4.5);//stable, conversions?
            //printf("Other Meson pdg %d, Mother %s, pdg %d, status %d\n",pdg, TDatabasePDG::Instance()->GetParticle(mompdg)->GetName(),mompdg, momstatus );          
          }
          
          
          //Check if both photons hit Calorimeter
          if(prim->GetNDaughters()!=2) return; //Only interested in 2 gamma decay
          Int_t iphot1=prim->GetFirstDaughter() ;
          Int_t iphot2=prim->GetLastDaughter() ;
          if(iphot1>-1 && iphot1<stack->GetNtrack() && iphot2>-1 && iphot2<stack->GetNtrack()){
            TParticle * phot1 = stack->Particle(iphot1) ;
            TParticle * phot2 = stack->Particle(iphot2) ;
            if(phot1 && phot2 && phot1->GetPdgCode()==22 && phot2->GetPdgCode()==22){
              //printf("2 photons: photon 1: pt %2.2f, phi %3.2f, eta %1.2f; photon 2: pt %2.2f, phi %3.2f, eta %1.2f\n",
              //	phot1->Pt(), phot1->Phi()*180./3.1415, phot1->Eta(), phot2->Pt(), phot2->Phi()*180./3.1415, phot2->Eta());
              
              TLorentzVector lv1, lv2;
              phot1->Momentum(lv1);
              phot2->Momentum(lv2);
              
              Bool_t inacceptance = kFALSE;
              if(fCalorimeter == "PHOS"){
                if(GetPHOSGeometry() && GetCaloUtils()->IsPHOSGeoMatrixSet()){
                  Int_t mod ;
                  Double_t x,z ;
                  if(GetPHOSGeometry()->ImpactOnEmc(phot1,mod,z,x) && GetPHOSGeometry()->ImpactOnEmc(phot2,mod,z,x)) 
                    inacceptance = kTRUE;
                  if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
                else{
                  
                  if(GetFiducialCut()->IsInFiducialCut(lv1,fCalorimeter) && GetFiducialCut()->IsInFiducialCut(lv2,fCalorimeter)) 
                    inacceptance = kTRUE ;
                  if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
                
              }	   
              else if(fCalorimeter == "EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
                if(GetEMCALGeometry()){
                  
                  Int_t absID1=0;
                  Int_t absID2=0;
                  
                  GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(phot1->Eta(),phot1->Phi(),absID1);
                  GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(phot2->Eta(),phot2->Phi(),absID2);

                  if( absID1 >= 0 && absID2 >= 0) 
                    inacceptance = kTRUE;
                  
//                  if(GetEMCALGeometry()->Impact(phot1) && GetEMCALGeometry()->Impact(phot2)) 
//                    inacceptance = kTRUE;
                  if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
                else{
                  if(GetFiducialCut()->IsInFiducialCut(lv1,fCalorimeter) && GetFiducialCut()->IsInFiducialCut(lv2,fCalorimeter)) 
                    inacceptance = kTRUE ;
                  if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
              }	  
              
              if(inacceptance){
                if(pdg==111){
                  fhPrimPi0AccPt ->Fill(pi0Pt) ;
                  fhPrimPi0AccPhi->Fill(pi0Pt, phi) ;
                  fhPrimPi0AccY  ->Fill(pi0Pt, pi0Y) ;
                  Double_t angle  = lv1.Angle(lv2.Vect());
                  fhPrimPi0OpeningAngle   ->Fill(pi0Pt,angle);
                  fhPrimPi0CosOpeningAngle->Fill(pi0Pt,TMath::Cos(angle));
                }
                else if(pdg==221){
                  fhPrimEtaAccPt ->Fill(pi0Pt) ;
                  fhPrimEtaAccPhi->Fill(pi0Pt, phi) ;
                  fhPrimEtaAccY  ->Fill(pi0Pt, pi0Y) ;
                }
              }//Accepted
            }// 2 photons      
          }//Check daughters exist
        }// Primary pi0 or eta
      }//loop on primaries	
    }//stack exists and data is MC
  }//read stack
  else if(GetReader()->ReadAODMCParticles()){
    
    TClonesArray * mcparticles = GetReader()->GetAODMCParticles(0);
    if(mcparticles){
      Int_t nprim = mcparticles->GetEntriesFast();
      for(Int_t i=0; i < nprim; i++)
      {
        AliAODMCParticle * prim = (AliAODMCParticle *) mcparticles->At(i);   
        // Only generator particles
        if( prim->GetStatus() == 0) break;

        Int_t pdg = prim->GetPdgCode();
        if( pdg == 111 || pdg == 221){
          Double_t pi0Pt = prim->Pt() ;
          //printf("pi0, pt %2.2f\n",pi0Pt);
          if(prim->E() == TMath::Abs(prim->Pz()))  continue ; //Protection against floating point exception	  
          Double_t pi0Y  = 0.5*TMath::Log((prim->E()-prim->Pz())/(prim->E()+prim->Pz())) ;
          Double_t phi   = TMath::RadToDeg()*prim->Phi() ;
          if(pdg == 111){
            if(TMath::Abs(pi0Y) < 0.5){
              fhPrimPi0Pt->Fill(pi0Pt) ;
            }
            fhPrimPi0Y  ->Fill(pi0Pt, pi0Y) ;
            fhPrimPi0Phi->Fill(pi0Pt, phi) ;
          }
          else if(pdg == 221){
            if(TMath::Abs(pi0Y) < 0.5){
              fhPrimEtaPt->Fill(pi0Pt) ;
            }
            fhPrimEtaY  ->Fill(pi0Pt, pi0Y) ;
            fhPrimEtaPhi->Fill(pi0Pt, phi) ;
          }
          
          //Origin of meson
          Int_t momindex  = prim->GetMother();
          if(momindex < 0) continue;
          AliAODMCParticle* mother = (AliAODMCParticle *) mcparticles->At(momindex);
          Int_t mompdg    = TMath::Abs(mother->GetPdgCode());
          Int_t momstatus = mother->GetStatus();
            if(pdg == 111){
              if     (momstatus  == 21) fhPrimPi0PtOrigin->Fill(pi0Pt,0.5);//parton
              else if(mompdg     < 22 ) fhPrimPi0PtOrigin->Fill(pi0Pt,1.5);//quark
              else if(mompdg     > 2100  && mompdg   < 2210) fhPrimPi0PtOrigin->Fill(pi0Pt,2.5);// resonances
              else if(mompdg    == 221) fhPrimPi0PtOrigin->Fill(pi0Pt,8.5);//eta
              else if(mompdg    == 331) fhPrimPi0PtOrigin->Fill(pi0Pt,9.5);//eta prime
              else if(mompdg    == 213) fhPrimPi0PtOrigin->Fill(pi0Pt,4.5);//rho
              else if(mompdg    == 223) fhPrimPi0PtOrigin->Fill(pi0Pt,5.5);//omega
              else if(mompdg    >= 310   && mompdg    <= 323) fhPrimPi0PtOrigin->Fill(pi0Pt,6.5);//k0S, k+-,k*
              else if(mompdg    == 130) fhPrimPi0PtOrigin->Fill(pi0Pt,6.5);//k0L
              else if(momstatus == 11 || momstatus  == 12 ) fhPrimPi0PtOrigin->Fill(pi0Pt,3.5);//resonances   
              else                      fhPrimPi0PtOrigin->Fill(pi0Pt,7.5);//other?
            }//pi0
            else {
              if     (momstatus == 21 ) fhPrimEtaPtOrigin->Fill(pi0Pt,0.5);//parton
              else if(mompdg    < 22  ) fhPrimEtaPtOrigin->Fill(pi0Pt,1.5);//quark
              else if(mompdg    > 2100  && mompdg   < 2210) fhPrimEtaPtOrigin->Fill(pi0Pt,2.5);//qq resonances
              else if(mompdg    == 331) fhPrimEtaPtOrigin->Fill(pi0Pt,5.5);//eta prime
              else if(momstatus == 11 || momstatus  == 12 ) fhPrimEtaPtOrigin->Fill(pi0Pt,3.5);//resonances
              else fhPrimEtaPtOrigin->Fill(pi0Pt,4.5);//stable, conversions?
              //printf("Other Meson pdg %d, Mother %s, pdg %d, status %d\n",pdg, TDatabasePDG::Instance()->GetParticle(mompdg)->GetName(),mompdg, momstatus );          
            }
          
          //Check if both photons hit Calorimeter
          if(prim->GetNDaughters()!=2) return; //Only interested in 2 gamma decay
          Int_t iphot1=prim->GetDaughter(0) ;
          Int_t iphot2=prim->GetDaughter(1) ;
          if(iphot1>-1 && iphot1<nprim && iphot2>-1 && iphot2<nprim){
            AliAODMCParticle * phot1 = (AliAODMCParticle *) mcparticles->At(iphot1);   
            AliAODMCParticle * phot2 = (AliAODMCParticle *) mcparticles->At(iphot2);   
            if(phot1 && phot2 && phot1->GetPdgCode()==22 && phot2->GetPdgCode()==22){
              TLorentzVector lv1, lv2;
              lv1.SetPxPyPzE(phot1->Px(),phot1->Py(),phot1->Pz(),phot1->E());
              lv2.SetPxPyPzE(phot2->Px(),phot2->Py(),phot2->Pz(),phot2->E());

              Bool_t inacceptance = kFALSE;
              if(fCalorimeter == "PHOS"){
                if(GetPHOSGeometry() && GetCaloUtils()->IsPHOSGeoMatrixSet()){
                  Int_t mod ;
                  Double_t x,z ;
                  Double_t vtx []={phot1->Xv(),phot1->Yv(),phot1->Zv()};
                  Double_t vtx2[]={phot2->Xv(),phot2->Yv(),phot2->Zv()};
                  if(GetPHOSGeometry()->ImpactOnEmc(vtx, phot1->Theta(),phot1->Phi(),mod,z,x) && 
                     GetPHOSGeometry()->ImpactOnEmc(vtx2,phot2->Theta(),phot2->Phi(),mod,z,x)) 
                    inacceptance = kTRUE;
                  if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
                else{
                  
                  if(GetFiducialCut()->IsInFiducialCut(lv1,fCalorimeter) && GetFiducialCut()->IsInFiducialCut(lv2,fCalorimeter)) 
                    inacceptance = kTRUE ;
                  if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
                
              }	   
              else if(fCalorimeter == "EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
                if(GetEMCALGeometry()){

                  Int_t absID1=0;
                  Int_t absID2=0;

                  //TVector3 vtx(phot1->Xv(),phot1->Yv(),phot1->Zv());
                  //TVector3 vimpact(0,0,0);
                  
                  //GetEMCALGeometry()->ImpactOnEmcal(vtx,phot1->Theta(),phot1->Phi(),absID1,vimpact);
                  //TVector3 vtx2(phot2->Xv(),phot2->Yv(),phot2->Zv());
                  //TVector3 vimpact2(0,0,0);
                  //GetEMCALGeometry()->ImpactOnEmcal(vtx2,phot2->Theta(),phot2->Phi(),absID2,vimpact2);
                 
                  GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(phot1->Eta(),phot1->Phi(),absID1);
                  GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(phot2->Eta(),phot2->Phi(),absID2);

//                  if(TMath::Abs(phot1->Eta()) < 0.7 && phot1->Phi() > 80*TMath::DegToRad() && phot1->Phi() < 120*TMath::DegToRad() ) 
//                    printf("Phot1 ccepted? %d\n",absID1);
//                  if(TMath::Abs(phot2->Eta()) < 0.7 && phot2->Phi() > 80*TMath::DegToRad() && phot2->Phi() < 120*TMath::DegToRad() ) 
//                    printf("Phot2 accepted? %d\n",absID2);

                  if( absID1 >= 0 && absID2 >= 0) 
                    inacceptance = kTRUE;
                  
//                  if(pdg==111 && inacceptance) printf("2 photons: photon 1: absId %d, pt %2.2f, phi %3.2f, eta %1.2f; photon 2: absId %d, pt %2.2f, phi %3.2f, eta %1.2f\n",
//                                                       absID1,phot1->Pt(), phot1->Phi()*TMath::RadToDeg(), phot1->Eta(),  
//                                                       absID2,phot2->Pt(), phot2->Phi()*TMath::RadToDeg(), phot2->Eta());
                  
                  
                  
                  if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
                else{
                  if(GetFiducialCut()->IsInFiducialCut(lv1,fCalorimeter) && GetFiducialCut()->IsInFiducialCut(lv2,fCalorimeter)) 
                    inacceptance = kTRUE ;
                  if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                }
              }	  
              
              if(inacceptance){
                if(pdg==111){
  //                printf("ACCEPTED pi0: pt %2.2f, phi %3.2f, eta %1.2f\n",pi0Pt,phi,pi0Y);
                  fhPrimPi0AccPt ->Fill(pi0Pt) ;
                  fhPrimPi0AccPhi->Fill(pi0Pt, phi) ;
                  fhPrimPi0AccY  ->Fill(pi0Pt, pi0Y) ;
                  Double_t angle  = lv1.Angle(lv2.Vect());
                  fhPrimPi0OpeningAngle   ->Fill(pi0Pt,angle);
                  fhPrimPi0CosOpeningAngle->Fill(pi0Pt,TMath::Cos(angle));
                }
                else if(pdg==221){
                  fhPrimEtaAccPt ->Fill(pi0Pt) ;
                  fhPrimEtaAccPhi->Fill(pi0Pt, phi) ;
                  fhPrimEtaAccY  ->Fill(pi0Pt, pi0Y) ;
                }
              }//Accepted
            }// 2 photons      
          }//Check daughters exist
        }// Primary pi0 or eta
      }//loop on primaries	
    }//stack exists and data is MC
    
    
  }	// read AOD MC
}

//_____________________________________________________________
void AliAnaPi0::FillMCVersusRecDataHistograms(const Int_t index1,  const Int_t index2, 
                                              const Float_t pt1,   const Float_t pt2, 
                                              const Int_t ncell1,  const Int_t ncell2,
                                              const Double_t mass, const Double_t pt, const Double_t asym, 
                                              const Double_t deta, const Double_t dphi){
  //Do some MC checks on the origin of the pair, is there any common ancestor and if there is one, who?
  //Adjusted for Pythia, need to see what to do for other generators.
  //Array of histograms ordered as follows: 0-Photon, 1-electron, 2-pi0, 3-eta, 4-a-proton, 5-a-neutron, 6-stable particles, 
  // 7-other decays, 8-string, 9-final parton, 10-initial parton, intermediate, 11-colliding proton, 12-unrelated
  
  Int_t ancPDG    = 0;
  Int_t ancStatus = 0;
  TLorentzVector ancMomentum;
  Int_t ancLabel  = GetMCAnalysisUtils()->CheckCommonAncestor(index1, index2, 
                                                              GetReader(), ancPDG, ancStatus,ancMomentum);
  
  Int_t momindex  = -1;
  Int_t mompdg    = -1;
  Int_t momstatus = -1;
  if(GetDebug() > 1) printf("AliAnaPi0::FillMCVersusRecDataHistograms() - Common ancestor label %d, pdg %d, name %s, status %d; \n",
                            ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
    
    if(ancLabel > -1){
      if(ancPDG==22){//gamma
        fhMCOrgMass[0]->Fill(pt,mass);
        fhMCOrgAsym[0]->Fill(pt,asym);
        fhMCOrgDeltaEta[0]->Fill(pt,deta);
        fhMCOrgDeltaPhi[0]->Fill(pt,dphi);
      }              
      else if(TMath::Abs(ancPDG)==11){//e
        fhMCOrgMass[1]->Fill(pt,mass);
        fhMCOrgAsym[1]->Fill(pt,asym);
        fhMCOrgDeltaEta[1]->Fill(pt,deta);
        fhMCOrgDeltaPhi[1]->Fill(pt,dphi);
      }          
      else if(ancPDG==111){//Pi0
        fhMCOrgMass[2]->Fill(pt,mass);
        fhMCOrgAsym[2]->Fill(pt,asym);
        fhMCOrgDeltaEta[2]->Fill(pt,deta);
        fhMCOrgDeltaPhi[2]->Fill(pt,dphi);
        if(fMultiCutAnaSim){
          for(Int_t ipt=0; ipt<fNPtCuts; ipt++){          
            for(Int_t icell=0; icell<fNCellNCuts; icell++){
              for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
                Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                if(pt1    >  fPtCuts[ipt]      && pt2    >  fPtCuts[ipt]        && 
                   asym   <  fAsymCuts[iasym]                                   && 
                   ncell1 >= fCellNCuts[icell] && ncell2 >= fCellNCuts[icell]){ 
                  fhMCPi0MassPtRec [index]->Fill(pt,mass);
                  fhMCPi0MassPtTrue[index]->Fill(ancMomentum.Pt(),mass);
                  if(mass < 0.17 && mass > 0.1) fhMCPi0PtTruePtRec[index]->Fill(ancMomentum.Pt(),pt);
                }//pass the different cuts
              }// pid bit cut loop
            }// icell loop
          }// pt cut loop
        }//Multi cut ana sim
        else {
          fhMCPi0MassPtTrue[0]->Fill(ancMomentum.Pt(),mass);
          if(mass < 0.17 && mass > 0.1) {
            fhMCPi0PtTruePtRec[0]->Fill(ancMomentum.Pt(),pt); 
   
            if(GetReader()->ReadStack()){	
              TParticle* ancestor = GetMCStack()->Particle(ancLabel);
              momindex  = ancestor->GetFirstMother();
              if(momindex < 0) return;
              TParticle* mother = GetMCStack()->Particle(momindex);
              mompdg    = TMath::Abs(mother->GetPdgCode());
              momstatus = mother->GetStatusCode();         
            }
            else {
              TClonesArray * mcparticles = GetReader()->GetAODMCParticles(0);
              AliAODMCParticle* ancestor = (AliAODMCParticle *) mcparticles->At(ancLabel);
              momindex  = ancestor->GetMother();
              if(momindex < 0) return;
              AliAODMCParticle* mother = (AliAODMCParticle *) mcparticles->At(momindex);
              mompdg    = TMath::Abs(mother->GetPdgCode());
              momstatus = mother->GetStatus();  
            }            
            
            if     (momstatus  == 21) fhMCPi0PtOrigin->Fill(pt,0.5);//parton
            else if(mompdg     < 22 ) fhMCPi0PtOrigin->Fill(pt,1.5);//quark
            else if(mompdg     > 2100  && mompdg   < 2210) fhMCPi0PtOrigin->Fill(pt,2.5);// resonances
            else if(mompdg    == 221) fhMCPi0PtOrigin->Fill(pt,8.5);//eta
            else if(mompdg    == 331) fhMCPi0PtOrigin->Fill(pt,9.5);//eta prime
            else if(mompdg    == 213) fhMCPi0PtOrigin->Fill(pt,4.5);//rho
            else if(mompdg    == 223) fhMCPi0PtOrigin->Fill(pt,5.5);//omega
            else if(mompdg    >= 310   && mompdg    <= 323) fhMCPi0PtOrigin->Fill(pt,6.5);//k0S, k+-,k*
            else if(mompdg    == 130) fhMCPi0PtOrigin->Fill(pt,6.5);//k0L
            else if(momstatus == 11 || momstatus  == 12 ) fhMCPi0PtOrigin->Fill(pt,3.5);//resonances   
            else                      fhMCPi0PtOrigin->Fill(pt,7.5);//other?
            
          }//pi0 mass region
          
        }
      }
      else if(ancPDG==221){//Eta
        fhMCOrgMass[3]->Fill(pt,mass);
        fhMCOrgAsym[3]->Fill(pt,asym);
        fhMCOrgDeltaEta[3]->Fill(pt,deta);
        fhMCOrgDeltaPhi[3]->Fill(pt,dphi);
        if(fMultiCutAnaSim){
          for(Int_t ipt=0; ipt<fNPtCuts; ipt++){          
            for(Int_t icell=0; icell<fNCellNCuts; icell++){
              for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
                Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                if(pt1    >  fPtCuts[ipt]      && pt2    >  fPtCuts[ipt]        && 
                   asym   <  fAsymCuts[iasym]                                   && 
                   ncell1 >= fCellNCuts[icell] && ncell2 >= fCellNCuts[icell]){ 
                  fhMCEtaMassPtRec [index]->Fill(pt,mass);
                  fhMCEtaMassPtTrue[index]->Fill(ancMomentum.Pt(),mass);
                  if(mass < 0.65 && mass > 0.45) fhMCEtaPtTruePtRec[index]->Fill(ancMomentum.Pt(),pt);
                }//pass the different cuts
              }// pid bit cut loop
            }// icell loop
          }// pt cut loop
        } //Multi cut ana sim
        else {
          fhMCEtaMassPtTrue[0]->Fill(ancMomentum.Pt(),mass);
          if(mass < 0.65 && mass > 0.45) fhMCEtaPtTruePtRec[0]->Fill(ancMomentum.Pt(),pt); 
          
          if(GetReader()->ReadStack()){	
            TParticle* ancestor = GetMCStack()->Particle(ancLabel);
            momindex  = ancestor->GetFirstMother();
            if(momindex < 0) return;
            TParticle* mother = GetMCStack()->Particle(momindex);
            mompdg    = TMath::Abs(mother->GetPdgCode());
            momstatus = mother->GetStatusCode();         
          }
          else {
            TClonesArray * mcparticles = GetReader()->GetAODMCParticles(0);
            AliAODMCParticle* ancestor = (AliAODMCParticle *) mcparticles->At(ancLabel);
            momindex  = ancestor->GetMother();
            if(momindex < 0) return;
            AliAODMCParticle* mother = (AliAODMCParticle *) mcparticles->At(momindex);
            mompdg    = TMath::Abs(mother->GetPdgCode());
            momstatus = mother->GetStatus();  
          }          
          
          if     (momstatus == 21 ) fhMCEtaPtOrigin->Fill(pt,0.5);//parton
          else if(mompdg    < 22  ) fhMCEtaPtOrigin->Fill(pt,1.5);//quark
          else if(mompdg    > 2100  && mompdg   < 2210) fhMCEtaPtOrigin->Fill(pt,2.5);//qq resonances
          else if(mompdg    == 331) fhMCEtaPtOrigin->Fill(pt,5.5);//eta prime
          else if(momstatus == 11 || momstatus  == 12 ) fhMCEtaPtOrigin->Fill(pt,3.5);//resonances
          else fhMCEtaPtOrigin->Fill(pt,4.5);//stable, conversions?
          //printf("Other Meson pdg %d, Mother %s, pdg %d, status %d\n",pdg, TDatabasePDG::Instance()->GetParticle(mompdg)->GetName(),mompdg, momstatus );  
        }// eta mass region
      }
      else if(ancPDG==-2212){//AProton
        fhMCOrgMass[4]->Fill(pt,mass);
        fhMCOrgAsym[4]->Fill(pt,asym);
        fhMCOrgDeltaEta[4]->Fill(pt,deta);
        fhMCOrgDeltaPhi[4]->Fill(pt,dphi);
      }   
      else if(ancPDG==-2112){//ANeutron
        fhMCOrgMass[5]->Fill(pt,mass);
        fhMCOrgAsym[5]->Fill(pt,asym);
        fhMCOrgDeltaEta[5]->Fill(pt,deta);
        fhMCOrgDeltaPhi[5]->Fill(pt,dphi);
      }       
      else if(TMath::Abs(ancPDG)==13){//muons
        fhMCOrgMass[6]->Fill(pt,mass);
        fhMCOrgAsym[6]->Fill(pt,asym);
        fhMCOrgDeltaEta[6]->Fill(pt,deta);
        fhMCOrgDeltaPhi[6]->Fill(pt,dphi);
      }                   
      else if (TMath::Abs(ancPDG) > 100 && ancLabel > 7) {
        if(ancStatus==1){//Stable particles, converted? not decayed resonances
          fhMCOrgMass[6]->Fill(pt,mass);
          fhMCOrgAsym[6]->Fill(pt,asym);
          fhMCOrgDeltaEta[6]->Fill(pt,deta);
          fhMCOrgDeltaPhi[6]->Fill(pt,dphi);  
        }
        else{//resonances and other decays, more hadron conversions?
          fhMCOrgMass[7]->Fill(pt,mass);
          fhMCOrgAsym[7]->Fill(pt,asym);
          fhMCOrgDeltaEta[7]->Fill(pt,deta);
          fhMCOrgDeltaPhi[7]->Fill(pt,dphi);
        }
      }
      else {//Partons, colliding protons, strings, intermediate corrections
        if(ancStatus==11 || ancStatus==12){//String fragmentation
          fhMCOrgMass[8]->Fill(pt,mass);
          fhMCOrgAsym[8]->Fill(pt,asym);
          fhMCOrgDeltaEta[8]->Fill(pt,deta);
          fhMCOrgDeltaPhi[8]->Fill(pt,dphi);
        }
        else if (ancStatus==21){
          if(ancLabel < 2) {//Colliding protons
            fhMCOrgMass[11]->Fill(pt,mass);
            fhMCOrgAsym[11]->Fill(pt,asym);
            fhMCOrgDeltaEta[11]->Fill(pt,deta);
            fhMCOrgDeltaPhi[11]->Fill(pt,dphi);
          }//colliding protons  
          else if(ancLabel < 6){//partonic initial states interactions
            fhMCOrgMass[9]->Fill(pt,mass);
            fhMCOrgAsym[9]->Fill(pt,asym);
            fhMCOrgDeltaEta[9]->Fill(pt,deta);
            fhMCOrgDeltaPhi[9]->Fill(pt,dphi);
          }
          else if(ancLabel < 8){//Final state partons radiations?
            fhMCOrgMass[10]->Fill(pt,mass);
            fhMCOrgAsym[10]->Fill(pt,asym);
            fhMCOrgDeltaEta[10]->Fill(pt,deta);
            fhMCOrgDeltaPhi[10]->Fill(pt,dphi);
          }
          else {
            printf("AliAnaPi0::FillMCVersusRecDataHistograms() - Check ** Common ancestor label %d, pdg %d, name %s, status %d; \n",
                   ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
          }
        }//status 21
        else {
          printf("AliAnaPi0::FillMCVersusRecDataHistograms() - Check *** Common ancestor label %d, pdg %d, name %s, status %d; \n",
                 ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
        }
      }////Partons, colliding protons, strings, intermediate corrections
    }//ancLabel > -1 
    else { //ancLabel <= -1
      //printf("Not related at all label = %d\n",ancLabel);
      fhMCOrgMass[12]->Fill(pt,mass);
      fhMCOrgAsym[12]->Fill(pt,asym);
      fhMCOrgDeltaEta[12]->Fill(pt,deta);
      fhMCOrgDeltaPhi[12]->Fill(pt,dphi);
    }
}  

//____________________________________________________________________________________________________________________________________________________
void AliAnaPi0::MakeAnalysisFillHistograms() 
{
  //Process one event and extract photons from AOD branch 
  // filled with AliAnaPhoton and fill histos with invariant mass
  
  //In case of simulated data, fill acceptance histograms
  if(IsDataMC())FillAcceptanceHistograms();
  
  //if (GetReader()->GetEventNumber()%10000 == 0) 
  // printf("--- Event %d ---\n",GetReader()->GetEventNumber());
  
  //Init some variables
//Int_t   iRun     = (GetReader()->GetInputEvent())->GetRunNumber() ;
  Int_t   nPhot    = GetInputAODBranch()->GetEntriesFast() ;
  Int_t   nClus    = 0;
  Int_t   nCell    = 0;
  Float_t eClusTot = 0;
  Float_t eCellTot = 0;
  Float_t eDenClus = 0;
  Float_t eDenCell = 0;
//  Int_t   ncomb    = 0;
//  Float_t rtmp     = 0;
//  Float_t rtmpw    = 0;
//  Float_t rxz      = 0;
//  Float_t rxzw     = 0;
//  Float_t pos1[3];
//  Float_t pos2[3];
//  Float_t emax     = 0;
  
  if(GetDebug() > 1) 
    printf("AliAnaPi0::MakeAnalysisFillHistograms() - Photon entries %d\n", nPhot);

  //If less than photon 2 entries in the list, skip this event
  if(nPhot < 2 ) return ; 

  // Count the number of clusters and cells, in case multiplicity bins dependent on such numbers
  // are requested
  if(fCalorimeter=="EMCAL"){ 
    nClus = GetEMCALClusters()  ->GetEntriesFast();
    nCell = GetEMCALCells()->GetNumberOfCells();
    for(Int_t icl=0; icl < nClus; icl++) {
      Float_t e1 = ((AliVCluster*)GetEMCALClusters()->At(icl))->E();
      eClusTot +=  e1;
//      if(e1 > emax) emax = e1;
//      ((AliVCluster*)GetEMCALClusters()->At(icl))->GetPosition(pos1);
//      for(Int_t icl2=icl+1; icl2 < nClus; icl2++) {
//        Float_t e2 = ((AliVCluster*)GetEMCALClusters()->At(icl2))->E();
//        ((AliVCluster*)GetEMCALClusters()->At(icl2))->GetPosition(pos2);
//        rtmp  =  TMath::Sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0]) + (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
//        rtmpw =  TMath::Sqrt((pos1[0]*e1-pos2[0]*e2)*(pos1[0]*e1-pos2[0]*e2) + (pos1[2]*e1-pos2[2]*e2)*(pos1[2]*e1-pos2[2]*e2))/(e1+e2);
//        rxz  += rtmp;  
//        rxzw += rtmpw;
//        ncomb++;
//        fhClusterPairDist      ->Fill(rtmp);
//        fhClusterPairDistWeight->Fill(rtmpw);
//        //printf("Distance: %f; weighted  %f\n ",rtmp,rtmp/(e1+((AliVCluster*)GetEMCALClusters()->At(icl2))->E()));
//
//      }// second cluster loop
    }// first cluster
    
    for(Int_t jce=0; jce < nCell; jce++) eCellTot +=  GetEMCALCells()->GetAmplitude(jce);
  }
  else {                     
    nClus = GetPHOSClusters()->GetEntriesFast();
    nCell = GetPHOSCells()   ->GetNumberOfCells();
    for(Int_t icl=0; icl < nClus; icl++) {
      Float_t e1 = ((AliVCluster*)GetPHOSClusters()->At(icl))->E();
      eClusTot +=  e1;
//      ((AliVCluster*)GetPHOSClusters()->At(icl))->GetPosition(pos1);
//      for(Int_t icl2=icl+1; icl2 < nClus; icl2++) {
//        Float_t e2 = ((AliVCluster*)GetPHOSClusters()->At(icl2))->E();
//        ((AliVCluster*)GetPHOSClusters()->At(icl2))->GetPosition(pos2);
//        rtmp  = TMath::Sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0]) + (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
//        rtmpw =  TMath::Sqrt((pos1[0]*e1-pos2[0]*e2)*(pos1[0]*e1-pos2[0]*e2) + (pos1[2]*e1-pos2[2]*e2)*(pos1[2]*e1-pos2[2]*e2))/(e1+e2);
//        rxz  += rtmp;  
//        rxzw += rtmpw;
//        ncomb++;
//        fhClusterPairDist      ->Fill(rtmp);
//        fhClusterPairDistWeight->Fill(rtmpw);
//      }// second cluster loop
    }// first cluster
    for(Int_t jce=0; jce < nCell; jce++) eCellTot +=  GetPHOSCells()->GetAmplitude(jce);
  }
  if(GetDebug() > 1) 
    printf("AliAnaPi0::MakeAnalysisFillHistograms() - # Clusters %d, sum cluster E per SM %f,# Cells %d, sum cell E per SM %f\n", nClus,eClusTot,nCell,eCellTot);

  //Fill histograms with "energy density", ncell and nclust will be > 0 since there are at least 2 "photons"
  eDenClus = eClusTot/nClus;
  eDenCell = eCellTot/nCell;
  fhEDensityCluster      ->Fill(eDenClus);
  fhEDensityCell         ->Fill(eDenCell);
  fhEDensityCellvsCluster->Fill(eDenClus, eDenCell);
  //Fill the average number of cells or clusters per SM
  eClusTot /=fNModules;
  eCellTot /=fNModules;
  fhAverTotECluster      ->Fill(eClusTot);
  fhAverTotECell         ->Fill(eCellTot);
  fhAverTotECellvsCluster->Fill(eClusTot, eCellTot);
  //printf("Average Cluster: E %f, density %f;  Average Cell E %f, density  %f\n ",eClusTot,eDenClus,eCellTot,eDenCell);

//  //Average weighted pair distance
//  rxz  /= ncomb;                             
//  rxzw /= ncomb;                             
//
//  fhAverClusterPairDist             ->Fill(rxz );
//  fhAverClusterPairDistWeight       ->Fill(rxzw);
//  fhAverClusterPairDistvsAverE      ->Fill(rxz ,eDenClus);
//  fhAverClusterPairDistWeightvsAverE->Fill(rxzw,eDenClus);
//  fhAverClusterPairDistvsN          ->Fill(rxz ,nClus);
//  fhAverClusterPairDistWeightvsN    ->Fill(rxzw,nClus);
//  
//  //emax
//  fhMaxEvsClustEDen->Fill(emax,eDenClus);
//  fhMaxEvsClustMult->Fill(emax,nPhot);
  
  //printf("Average Distance: %f; weighted  %f\n ",rxz,rxzw);

  
  //Init variables
  Int_t module1         = -1;
  Int_t module2         = -1;
  Double_t vert[]       = {0.0, 0.0, 0.0} ; //vertex 
  Int_t evtIndex1       = 0 ; 
  Int_t currentEvtIndex = -1; 
  Int_t curCentrBin     = 0 ; 
  Int_t curRPBin        = 0 ; 
  Int_t curZvertBin     = 0 ;
  
  //---------------------------------
  //First loop on photons/clusters
  //---------------------------------
  for(Int_t i1=0; i1<nPhot-1; i1++){
    AliAODPWG4Particle * p1 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i1)) ;
    //printf("AliAnaPi0::MakeAnalysisFillHistograms() : cluster1 id %d\n",p1->GetCaloLabel(0));

    // get the event index in the mixed buffer where the photon comes from 
    // in case of mixing with analysis frame, not own mixing
    evtIndex1 = GetEventIndex(p1, vert) ; 
    //printf("charge = %d\n", track->Charge());
    if ( evtIndex1 == -1 )
      return ; 
    if ( evtIndex1 == -2 )
      continue ; 
    
    //printf("z vertex %f < %f\n",vert[2],GetZvertexCut());
    if(TMath::Abs(vert[2]) > GetZvertexCut()) continue ;   //vertex cut
  
    
    //----------------------------------------------------------------------------
    // Get the multiplicity bin. Different cases: centrality (PbPb), 
    // average cluster multiplicity, average cell multiplicity, track multiplicity 
    // default is centrality bins
    //----------------------------------------------------------------------------
    if (evtIndex1 != currentEvtIndex) {
      if(fUseTrackMultBins){ // Track multiplicity bins
        //printf("track  mult %d\n",GetTrackMultiplicity());
        curCentrBin = (GetTrackMultiplicity()-1)/5; 
        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
        //printf("track mult bin %d\n",curCentrBin);
      }
      else if(fUsePhotonMultBins){ // Photon multiplicity bins
        //printf("photon  mult %d cluster mult %d\n",nPhot, nClus);
        curRPBin = nPhot-2; 
        if(curRPBin > GetNRPBin() -1) curRPBin=GetNRPBin()-1;
        //printf("photon mult bin %d\n",curRPBin);        
      }
      else if(fUseAverClusterEBins){ // Cluster average energy bins
        //Bins for pp, if needed can be done in a more general way
        curCentrBin = (Int_t) eClusTot/10 * fNCentrBin; 
        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
        //printf("cluster E average %f, bin %d \n",eClusTot,curCentrBin);
      }
      else if(fUseAverCellEBins){ // Cell average energy bins
        //Bins for pp, if needed can be done in a more general way
        curCentrBin = (Int_t) eCellTot/10*fNCentrBin; 
        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
        //printf("cell E average %f, bin %d \n",eCellTot,curCentrBin);
      }
      else if(fUseAverClusterEDenBins){ // Energy density bins
        //Bins for pp, if needed can be done in a more general way
        curCentrBin = (Int_t) eDenClus/10*fNCentrBin; 
        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
        //printf("cluster Eden average %f, bin %d \n",eDenClus,curCentrBin);
      } 
//      else if(fUseAverClusterPairRBins){ // Cluster average distance bins
//        //Bins for pp, if needed can be done in a more general way
//        curCentrBin = rxz/650*fNCentrBin; 
//        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
//        //printf("cluster pair R average %f, bin %d \n",rxz,curCentrBin);
//      }
//      else if(fUseAverClusterPairRWeightBins){ // Cluster average distance bins
//        //Bins for pp, if needed can be done in a more general way
//        curCentrBin = rxzw/350*fNCentrBin; 
//        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
//        //printf("cluster pair rW average %f, bin %d \n",rxzw,curCentrBin);
//      }    
//      else if(fUseEMaxBins){ // Cluster average distance bins
//        //Bins for pp, if needed can be done in a more general way
//        curCentrBin = emax/20*fNCentrBin; 
//        if(curCentrBin > fNCentrBin-1) curCentrBin=fNCentrBin-1;
//        //printf("cluster pair rW average %f, bin %d \n",rxzw,curCentrBin);
//      }     
      else { //Event centrality
          curCentrBin = GetEventCentrality();
        //printf("curCentrBin %d\n",curCentrBin);
      }
      
      if (curCentrBin < 0 || curCentrBin >= fNCentrBin){ 
        printf("AliAnaPi0::MakeAnalysisFillHistograms() - Centrality bin <%d> not expected, n bins <%d> , return\n",curCentrBin,fNCentrBin);
        return;
      }
        
      //Reaction plane bin
      curRPBin    = 0 ;
       
      //Get vertex z bin
      curZvertBin = (Int_t)(0.5*GetNZvertBin()*(vert[2]+GetZvertexCut())/GetZvertexCut()) ;
      
      //Fill event bin info
      fhEvents    ->Fill(curCentrBin+0.5,curZvertBin+0.5,curRPBin+0.5) ;
      fhCentrality->Fill(curCentrBin);
      currentEvtIndex = evtIndex1 ; 
      if(GetDebug() > 1) 
        printf("AliAnaPi0::MakeAnalysisFillHistograms() - Centrality %d, Vertex Bin %d, RP bin %d \n",curCentrBin,curRPBin,curZvertBin);
    }
    
    //printf("AliAnaPi0::MakeAnalysisFillHistograms(): Photon 1 Evt %d  Vertex : %f,%f,%f\n",evtIndex1, GetVertex(evtIndex1)[0] ,GetVertex(evtIndex1)[1],GetVertex(evtIndex1)[2]);
    
    //Get the momentum of this cluster
    TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
    
    //Get (Super)Module number of this cluster
    module1 = GetModuleNumber(p1);
    
    //---------------------------------
    //Second loop on photons/clusters
    //---------------------------------
    for(Int_t i2=i1+1; i2<nPhot; i2++){
      AliAODPWG4Particle * p2 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i2)) ;
      
      //In case of mixing frame, check we are not in the same event as the first cluster
      Int_t evtIndex2 = GetEventIndex(p2, vert) ; 
      if ( evtIndex2 == -1 )
        return ; 
      if ( evtIndex2 == -2 )
        continue ;    
      if (GetMixedEvent() && (evtIndex1 == evtIndex2))
        continue ;
      
      //printf("AliAnaPi0::MakeAnalysisFillHistograms(): Photon 2 Evt %d  Vertex : %f,%f,%f\n",evtIndex2, GetVertex(evtIndex2)[0] ,GetVertex(evtIndex2)[1],GetVertex(evtIndex2)[2]);
 
      //Get the momentum of this cluster
      TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
      //Get module number
      module2       = GetModuleNumber(p2);
      
      //---------------------------------
      // Get pair kinematics
      //---------------------------------
      Double_t m    = (photon1 + photon2).M() ;
      Double_t pt   = (photon1 + photon2).Pt();
      Double_t deta = photon1.Eta() - photon2.Eta();
      Double_t dphi = photon1.Phi() - photon2.Phi();
      Double_t a    = TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
      
      if(GetDebug() > 2)
        printf(" E: photon1 %f, photon2 %f; Pair: pT %f, mass %f, a %f\n", p1->E(), p2->E(), (photon1 + photon2).E(),m,a);
      
      //--------------------------------
      // Opening angle selection
      //--------------------------------
      //Check if opening angle is too large or too small compared to what is expected	
      Double_t angle   = photon1.Angle(photon2.Vect());
      if(fUseAngleEDepCut && !GetNeutralMesonSelection()->IsAngleInWindow((photon1+photon2).E(),angle+0.05)) {
        if(GetDebug() > 2)
          printf("AliAnaPi0::MakeAnalysisFillHistograms() -Real pair angle %f not in E %f window\n",angle, (photon1+photon2).E());
        continue;
      }

      if(fUseAngleCut && (angle < fAngleCut || angle > fAngleMaxCut)) {
        if(GetDebug() > 2)
          printf("AliAnaPi0::MakeAnalysisFillHistograms() - Real pair cut %f < angle %f < cut %f\n",fAngleCut, angle, fAngleMaxCut);
        continue;
      }

      //-------------------------------------------------------------------------------------------------
      //Fill module dependent histograms, put a cut on assymmetry on the first available cut in the array
      //-------------------------------------------------------------------------------------------------
      if(a < fAsymCuts[0]){
        if(module1==module2 && module1 >=0 && module1<fNModules)
          fhReMod[module1]->Fill(pt,m) ;

        if(fCalorimeter=="EMCAL"){
          
          // Same sector
          Int_t j=0;
          for(Int_t i = 0; i < fNModules/2; i++){
            j=2*i;
            if((module1==j && module2==j+1) || (module1==j+1 && module2==j)) fhReSameSectorEMCALMod[i]->Fill(pt,m) ;
          }
          
          // Same side
          for(Int_t i = 0; i < fNModules-2; i++){
            if((module1==i && module2==i+2) || (module1==i+2 && module2==i)) fhReSameSideEMCALMod[i]->Fill(pt,m); 
          }
        }//EMCAL
        else {//PHOS
          if((module1==0 && module2==1) || (module1==1 && module2==0)) fhReDiffPHOSMod[0]->Fill(pt,m) ; 
          if((module1==0 && module2==2) || (module1==2 && module2==0)) fhReDiffPHOSMod[1]->Fill(pt,m) ; 
          if((module1==1 && module2==2) || (module1==2 && module2==1)) fhReDiffPHOSMod[2]->Fill(pt,m) ;
        }//PHOS
      }
      
      //In case we want only pairs in same (super) module, check their origin.
      Bool_t ok = kTRUE;
      if(fSameSM && module1!=module2) ok=kFALSE;
      if(ok){
        
        //Check if one of the clusters comes from a conversion 
        if     (p1->IsTagged() && p2->IsTagged()) fhReConv2->Fill(pt,m);
        else if(p1->IsTagged() || p2->IsTagged()) fhReConv ->Fill(pt,m);
        
        //Fill histograms for different bad channel distance, centrality, assymmetry cut and pid bit
        for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
          if((p1->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton)) && (p2->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton))){ 
            for(Int_t iasym=0; iasym < fNAsymCuts; iasym++){
              if(a < fAsymCuts[iasym]){
                Int_t index = ((curCentrBin*fNPIDBits)+ipid)*fNAsymCuts + iasym;
                //printf("index %d :(cen %d * nPID %d + ipid %d)*nasym %d + iasym %d\n",index,curCentrBin,fNPIDBits,ipid,fNAsymCuts,iasym);
                fhRe1     [index]->Fill(pt,m);
                if(fMakeInvPtPlots)fhReInvPt1[index]->Fill(pt,m,1./pt) ;
                if(fFillBadDistHisto){
                  if(p1->DistToBad()>0 && p2->DistToBad()>0){
                    fhRe2     [index]->Fill(pt,m) ;
                    if(fMakeInvPtPlots)fhReInvPt2[index]->Fill(pt,m,1./pt) ;
                    if(p1->DistToBad()>1 && p2->DistToBad()>1){
                      fhRe3     [index]->Fill(pt,m) ;
                      if(fMakeInvPtPlots)fhReInvPt3[index]->Fill(pt,m,1./pt) ;
                    }// bad 3
                  }// bad2
                }// Fill bad dist histos
              }//assymetry cut
            }// asymmetry cut loop
          }// bad 1
        }// pid bit loop
        
        //Fill histograms with opening angle
        fhRealOpeningAngle   ->Fill(pt,angle);
        fhRealCosOpeningAngle->Fill(pt,TMath::Cos(angle));
        
        //Fill histograms with pair assymmetry
        fhRePtAsym->Fill(pt,a);
        if(m > 0.10 && m < 0.17) fhRePtAsymPi0->Fill(pt,a);
        if(m > 0.45 && m < 0.65) fhRePtAsymEta->Fill(pt,a);
        
        //-------------------------------------------------------
        //Get the number of cells needed for multi cut analysis.
        //-------------------------------------------------------        
        Int_t ncell1 = 0;
        Int_t ncell2 = 0;
        if(fMultiCutAna || (IsDataMC() && fMultiCutAnaSim)){

          AliVEvent * event = GetReader()->GetInputEvent();
          if(event){
            for(Int_t iclus = 0; iclus < event->GetNumberOfCaloClusters(); iclus++){
              AliVCluster *cluster = event->GetCaloCluster(iclus);
              
              Bool_t is = kFALSE;
              if     (fCalorimeter == "EMCAL" && GetReader()->IsEMCALCluster(cluster)) is = kTRUE;
              else if(fCalorimeter == "PHOS"  && GetReader()->IsPHOSCluster (cluster)) is = kTRUE;
              
              if(is){
                if      (p1->GetCaloLabel(0) == cluster->GetID()) ncell1 = cluster->GetNCells();
                else if (p2->GetCaloLabel(0) == cluster->GetID()) ncell2 = cluster->GetNCells();
              } // PHOS or EMCAL cluster as requested in analysis
              
              if(ncell2 > 0 &&  ncell1 > 0) break; // No need to continue the iteration
              
            }
            //printf("e 1: %2.2f, e 2: %2.2f, ncells: n1 %d, n2 %d\n", p1->E(), p2->E(),ncell1,ncell2);
          }
        }
        
        //---------
        // MC data
        //---------
        //Do some MC checks on the origin of the pair, is there any common ancestor and if there is one, who?
        if(IsDataMC()) FillMCVersusRecDataHistograms(p1->GetLabel(), p2->GetLabel(),p1->Pt(), p2->Pt(),ncell1, ncell2, m, pt, a,deta, dphi);    
        
        //-----------------------
        //Multi cuts analysis 
        //-----------------------
        if(fMultiCutAna){
          //Histograms for different PID bits selection
          for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
            
            if(p1->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton)    && 
               p2->IsPIDOK(fPIDBits[ipid],AliCaloPID::kPhoton))   fhRePIDBits[ipid]->Fill(pt,m) ;
            
            //printf("ipt %d, ipid%d, name %s\n",ipt, ipid, fhRePtPIDCuts[ipt*fNPIDBitsBits+ipid]->GetName());
          } // pid bit cut loop
          
          //Several pt,ncell and asymmetry cuts
          for(Int_t ipt=0; ipt<fNPtCuts; ipt++){          
            for(Int_t icell=0; icell<fNCellNCuts; icell++){
              for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
                Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                if(p1->Pt() >   fPtCuts[ipt]      && p2->Pt() > fPtCuts[ipt]        && 
                   a        <   fAsymCuts[iasym]                                    && 
                   ncell1   >=  fCellNCuts[icell] && ncell2   >= fCellNCuts[icell]){
                    fhRePtNCellAsymCuts[index]->Fill(pt,m) ;
                  //printf("ipt %d, icell%d, iasym %d, name %s\n",ipt, icell, iasym,  fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym]->GetName());
                  if(module1==module2){
                    if     (module1==0)  fhRePtNCellAsymCutsSM0[index]->Fill(pt,m) ;
                    else if(module1==1)  fhRePtNCellAsymCutsSM1[index]->Fill(pt,m) ;
                    else if(module1==2)  fhRePtNCellAsymCutsSM2[index]->Fill(pt,m) ;
                    else if(module1==3)  fhRePtNCellAsymCutsSM3[index]->Fill(pt,m) ;
                    else printf("AliAnaPi0::FillHistograms() - WRONG SM NUMBER\n");
                  }
                }
              }// pid bit cut loop
            }// icell loop
          }// pt cut loop
          for(Int_t iasym = 0; iasym < fNAsymCuts; iasym++){
            if(a < fAsymCuts[iasym])fhRePtMult[iasym]->Fill(pt,GetTrackMultiplicity(),m) ;
          }
        }// multiple cuts analysis
      }// ok if same sm
    }// second same event particle
  }// first cluster
        
  //-------------------------------------------------------------
  // Mixing
  //-------------------------------------------------------------
  if(fDoOwnMix){
    //printf("Cen bin %d, RP bin %d, e aver %f, mult %d\n",curCentrBin,curRPBin, eClusTot, nPhot);
    //Recover events in with same characteristics as the current event
    TList * evMixList=fEventsList[curCentrBin*GetNZvertBin()*GetNRPBin()+curZvertBin*GetNRPBin()+curRPBin] ;
    Int_t nMixed = evMixList->GetSize() ;
    for(Int_t ii=0; ii<nMixed; ii++){  
      TClonesArray* ev2= (TClonesArray*) (evMixList->At(ii));
      Int_t nPhot2=ev2->GetEntriesFast() ;
      Double_t m = -999;
      if(GetDebug() > 1) 
        printf("AliAnaPi0::MakeAnalysisFillHistograms() - Mixed event %d photon entries %d, centrality bin %d\n", ii, nPhot2, curCentrBin);
      
      //---------------------------------
      //First loop on photons/clusters
      //---------------------------------      
      for(Int_t i1=0; i1<nPhot; i1++){
        AliAODPWG4Particle * p1 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i1)) ;
        if(fSameSM && GetModuleNumber(p1)!=module1) continue;
        
        //Get kinematics of cluster and (super) module of this cluster
        TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
        module1 = GetModuleNumber(p1);
        
        //---------------------------------
        //First loop on photons/clusters
        //---------------------------------        
        for(Int_t i2=0; i2<nPhot2; i2++){
          AliAODPWG4Particle * p2 = (AliAODPWG4Particle*) (ev2->At(i2)) ;
          
          //Get kinematics of second cluster and calculate those of the pair
          TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
          m           = (photon1+photon2).M() ; 
          Double_t pt = (photon1 + photon2).Pt();
          Double_t a  = TMath::Abs(p1->E()-p2->E())/(p1->E()+p2->E()) ;
          
          //Check if opening angle is too large or too small compared to what is expected
          Double_t angle   = photon1.Angle(photon2.Vect());
          if(fUseAngleEDepCut && !GetNeutralMesonSelection()->IsAngleInWindow((photon1+photon2).E(),angle+0.05)){ 
            if(GetDebug() > 2)
              printf("AliAnaPi0::MakeAnalysisFillHistograms() -Mix pair angle %f not in E %f window\n",angle, (photon1+photon2).E());
            continue;
          }
          if(fUseAngleCut && (angle < fAngleCut || angle > fAngleMaxCut)) {
            if(GetDebug() > 2)
              printf("AliAnaPi0::MakeAnalysisFillHistograms() -Mix pair angle %f < cut %f\n",angle,fAngleCut);
            continue; 
          
          } 
          
          if(GetDebug() > 2)
            printf("AliAnaPi0::MakeAnalysisFillHistograms() - Mixed Event: pT: photon1 %2.2f, photon2 %2.2f; Pair: pT %2.2f, mass %2.3f, a %f2.3\n",
                   p1->Pt(), p2->Pt(), pt,m,a);	
          
          //In case we want only pairs in same (super) module, check their origin.
          module2 = GetModuleNumber(p2);
          
          //-------------------------------------------------------------------------------------------------
          //Fill module dependent histograms, put a cut on assymmetry on the first available cut in the array
          //-------------------------------------------------------------------------------------------------          
          if(a < fAsymCuts[0]){
            if(module1==module2 && module1 >=0 && module1<fNModules)
              fhMiMod[module1]->Fill(pt,m) ;

            if(fCalorimeter=="EMCAL"){
              
              // Same sector
              Int_t j=0;
              for(Int_t i = 0; i < fNModules/2; i++){
                j=2*i;
                if((module1==j && module2==j+1) || (module1==j+1 && module2==j)) fhMiSameSectorEMCALMod[i]->Fill(pt,m) ;
              }
              
              // Same side
              for(Int_t i = 0; i < fNModules-2; i++){
                if((module1==i && module2==i+2) || (module1==i+2 && module2==i)) fhMiSameSideEMCALMod[i]->Fill(pt,m); 
              }
            }//EMCAL
            else {//PHOS
              if((module1==0 && module2==1) || (module1==1 && module2==0)) fhMiDiffPHOSMod[0]->Fill(pt,m) ; 
              if((module1==0 && module2==2) || (module1==2 && module2==0)) fhMiDiffPHOSMod[1]->Fill(pt,m) ; 
              if((module1==1 && module2==2) || (module1==2 && module2==1)) fhMiDiffPHOSMod[2]->Fill(pt,m) ;
            }//PHOS
            
            
          }
          
          Bool_t ok = kTRUE;
          if(fSameSM && module1!=module2) ok=kFALSE;
          if(ok){
            
            //Check if one of the clusters comes from a conversion 
            if     (p1->IsTagged() && p2->IsTagged()) fhMiConv2->Fill(pt,m);
            else if(p1->IsTagged() || p2->IsTagged()) fhMiConv ->Fill(pt,m);

            //Fill histograms for different bad channel distance, centrality, assymmetry cut and pid bit
            for(Int_t ipid=0; ipid<fNPIDBits; ipid++){ 
              if((p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton))){ 
                for(Int_t iasym=0; iasym < fNAsymCuts; iasym++){
                  if(a < fAsymCuts[iasym]){
                    Int_t index = ((curCentrBin*fNPIDBits)+ipid)*fNAsymCuts + iasym;
                    fhMi1     [index]->Fill(pt,m) ;
                    if(fMakeInvPtPlots)fhMiInvPt1[index]->Fill(pt,m,1./pt) ;
                    if(fFillBadDistHisto){
                      if(p1->DistToBad()>0 && p2->DistToBad()>0){
                        fhMi2     [index]->Fill(pt,m) ;
                        if(fMakeInvPtPlots)fhMiInvPt2[index]->Fill(pt,m,1./pt) ;
                        if(p1->DistToBad()>1 && p2->DistToBad()>1){
                          fhMi3     [index]->Fill(pt,m) ;
                          if(fMakeInvPtPlots)fhMiInvPt3[index]->Fill(pt,m,1./pt) ;
                        }
                      }
                    }// Fill bad dist histo
                  }//Asymmetry cut
                }// Asymmetry loop
              }//PID cut
            }//loop for histograms
            
            //-----------------------
            //Multi cuts analysis 
            //-----------------------            
            if(fMultiCutAna){
              //Several pt,ncell and asymmetry cuts
              for(Int_t ipt=0; ipt<fNPtCuts; ipt++){          
                for(Int_t icell=0; icell<fNCellNCuts; icell++){
                  for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
                    Int_t index = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
                    if(p1->Pt() >   fPtCuts[ipt]      && p2->Pt() > fPtCuts[ipt]        && 
                       a        <   fAsymCuts[iasym]                                    && 
                       p1->GetBtag() >=  fCellNCuts[icell] && p2->GetBtag() >= fCellNCuts[icell]){
                      fhMiPtNCellAsymCuts[index]->Fill(pt,m) ;
                      //printf("ipt %d, icell%d, iasym %d, name %s\n",ipt, icell, iasym,  fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym]->GetName());
                    }
                  }// pid bit cut loop
                }// icell loop
              }// pt cut loop
            } // Multi cut ana
            
            //Fill histograms with opening angle
            fhMixedOpeningAngle   ->Fill(pt,angle);
            fhMixedCosOpeningAngle->Fill(pt,TMath::Cos(angle));          
          }//ok
        }// second cluster loop
      }//first cluster loop
    }//loop on mixed events
    
    //--------------------------------------------------------
    //Add the current event to the list of events for mixing
    //--------------------------------------------------------
    TClonesArray *currentEvent = new TClonesArray(*GetInputAODBranch());
    //Add current event to buffer and Remove redundant events 
    if(currentEvent->GetEntriesFast()>0){
      evMixList->AddFirst(currentEvent) ;
      currentEvent=0 ; //Now list of particles belongs to buffer and it will be deleted with buffer
      if(evMixList->GetSize()>=fNmaxMixEv)
      {
        TClonesArray * tmp = (TClonesArray*) (evMixList->Last()) ;
        evMixList->RemoveLast() ;
        delete tmp ;
      }
    } 
    else{ //empty event
      delete currentEvent ;
      currentEvent=0 ; 
    }
  }// DoOwnMix
  
}	

//________________________________________________________________________
void AliAnaPi0::ReadHistograms(TList* outputList)
{
  // Needed when Terminate is executed in distributed environment
  // Refill analysis histograms of this class with corresponding histograms in output list. 
  
  // Histograms of this analsys are kept in the same list as other analysis, recover the position of
  // the first one and then add the next.
  Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"hRe_cen0_pid0_dist1"));
  
  if(!fhRe1) fhRe1 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhRe2) fhRe2 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhRe3) fhRe3 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhMi1) fhMi1 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhMi2) fhMi2 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhMi3) fhMi3 = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;	
  if(!fhReInvPt1) fhReInvPt1  = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhReInvPt2) fhReInvPt2  = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhReInvPt3) fhReInvPt3  = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhMiInvPt1) fhMiInvPt1  = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhMiInvPt2) fhMiInvPt2  = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;
  if(!fhMiInvPt3) fhMiInvPt3  = new TH2D*[fNCentrBin*fNPIDBits*fNAsymCuts] ;	
  if(!fhReMod)    fhReMod     = new TH2D*[fNModules]   ;	
  if(!fhReDiffPHOSMod)       fhReDiffPHOSMod        = new TH2D*[fNModules] ;	
  if(!fhReSameSectorEMCALMod)fhReSameSectorEMCALMod = new TH2D*[fNModules/2] ;	
  if(!fhReSameSideEMCALMod)  fhReSameSideEMCALMod   = new TH2D*[fNModules-2] ;	
  if(!fhMiMod)    fhMiMod     = new TH2D*[fNModules]   ;	
  if(!fhMiDiffPHOSMod)       fhMiDiffPHOSMod        = new TH2D*[fNModules] ;	
  if(!fhMiSameSectorEMCALMod)fhMiSameSectorEMCALMod = new TH2D*[fNModules/2] ;	
  if(!fhMiSameSideEMCALMod)  fhMiSameSideEMCALMod   = new TH2D*[fNModules-2] ;	
    
  fhReConv  = (TH2D*) outputList->At(index++);
  fhMiConv  = (TH2D*) outputList->At(index++);
  fhReConv2 = (TH2D*) outputList->At(index++);
  fhMiConv2 = (TH2D*) outputList->At(index++);

  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
      for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
        Int_t ihisto = ((ic*fNPIDBits)+ipid)*fNAsymCuts + iasym;

        fhRe1[ihisto] = (TH2D*) outputList->At(index++);
        fhRe2[ihisto] = (TH2D*) outputList->At(index++);
        fhRe3[ihisto] = (TH2D*) outputList->At(index++);
      
        fhReInvPt1[ihisto] = (TH2D*) outputList->At(index++);
        fhReInvPt2[ihisto] = (TH2D*) outputList->At(index++);
        fhReInvPt3[ihisto] = (TH2D*) outputList->At(index++);
      
        if(fDoOwnMix){
          fhMi1[ihisto] = (TH2D*) outputList->At(index++);
          fhMi2[ihisto] = (TH2D*) outputList->At(index++);
          fhMi3[ihisto] = (TH2D*) outputList->At(index++);
      
          fhMiInvPt1[ihisto] = (TH2D*) outputList->At(index++);
          fhMiInvPt2[ihisto] = (TH2D*) outputList->At(index++);
          fhMiInvPt3[ihisto] = (TH2D*) outputList->At(index++); 
        }//Own mix
      }//asymmetry loop
    }// pid loop
  }// centrality loop
  
  fhRePtAsym    = (TH2D*)outputList->At(index++);
  fhRePtAsymPi0 = (TH2D*)outputList->At(index++);
  fhRePtAsymEta = (TH2D*)outputList->At(index++);
  
  if(fMultiCutAna){
    
    if(!fhRePtNCellAsymCuts) fhRePtNCellAsymCuts = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
    if(!fhRePIDBits)         fhRePIDBits         = new TH2D*[fNPIDBits];

    for(Int_t ipid=0; ipid<fNPIDBits; ipid++){
      fhRePIDBits[ipid] = (TH2D*) outputList->At(index++);
    }// ipid loop
    
    for(Int_t ipt=0; ipt<fNPtCuts; ipt++){
      for(Int_t icell=0; icell<fNCellNCuts; icell++){
        for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
          fhRePtNCellAsymCuts[((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym] = (TH2D*) outputList->At(index++);
        }// iasym
      }// icell loop
    }// ipt loop
    
    if(!fhRePtMult) fhRePtMult  = new TH3D*[fNAsymCuts]  ;
    for(Int_t iasym = 0; iasym < fNAsymCuts; iasym++)
      fhRePtMult[iasym] = (TH3D*) outputList->At(index++);
  }// multi cut analysis 
  
  fhEvents     = (TH3D *) outputList->At(index++); 
  fhCentrality = (TH1D *) outputList->At(index++); 

  fhRealOpeningAngle     = (TH2D*)  outputList->At(index++);
  fhRealCosOpeningAngle  = (TH2D*)  outputList->At(index++);
  if(fDoOwnMix){
    fhMixedOpeningAngle     = (TH2D*)  outputList->At(index++);
    fhMixedCosOpeningAngle  = (TH2D*)  outputList->At(index++);
  }
  
  //Histograms filled only if MC data is requested 	
  if(IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC) ){
    fhPrimPi0Pt     = (TH1D*)  outputList->At(index++);
    fhPrimPi0AccPt  = (TH1D*)  outputList->At(index++);
    fhPrimPi0Y      = (TH2D*)  outputList->At(index++);
    fhPrimPi0AccY   = (TH2D*)  outputList->At(index++);
    fhPrimPi0Phi    = (TH2D*)  outputList->At(index++);
    fhPrimPi0AccPhi = (TH2D*)  outputList->At(index++);
    fhPrimEtaPt     = (TH1D*)  outputList->At(index++);
    fhPrimEtaAccPt  = (TH1D*)  outputList->At(index++);
    fhPrimEtaY      = (TH2D*)  outputList->At(index++);
    fhPrimEtaAccY   = (TH2D*)  outputList->At(index++);
    fhPrimEtaPhi    = (TH2D*)  outputList->At(index++);
    fhPrimEtaAccPhi = (TH2D*)  outputList->At(index++);
    for(Int_t i = 0; i<13; i++){
      fhMCOrgMass[i]     = (TH2D*)  outputList->At(index++);
      fhMCOrgAsym[i]     = (TH2D*)  outputList->At(index++);
      fhMCOrgDeltaEta[i] = (TH2D*)  outputList->At(index++);
      fhMCOrgDeltaPhi[i] = (TH2D*)  outputList->At(index++);
    }
    
    if(fMultiCutAnaSim){
      fhMCPi0MassPtTrue  = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCPi0MassPtRec   = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCPi0PtTruePtRec = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCEtaMassPtTrue  = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCEtaMassPtRec   = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      fhMCEtaPtTruePtRec = new TH2D*[fNPtCuts*fNAsymCuts*fNCellNCuts];
      for(Int_t ipt=0; ipt<fNPtCuts; ipt++){
        for(Int_t icell=0; icell<fNCellNCuts; icell++){
          for(Int_t iasym=0; iasym<fNAsymCuts; iasym++){
            Int_t in = ((ipt*fNCellNCuts)+icell)*fNAsymCuts + iasym;
            fhMCPi0MassPtTrue[in]  = (TH2D*)  outputList->At(index++);
            fhMCPi0PtTruePtRec[in] = (TH2D*)  outputList->At(index++);
            fhMCEtaMassPtTrue[in]  = (TH2D*)  outputList->At(index++);
            fhMCEtaPtTruePtRec[in] = (TH2D*)  outputList->At(index++);
          }
        }
      }
    }
    else{
      fhMCPi0MassPtTrue  = new TH2D*[1];
      fhMCPi0PtTruePtRec = new TH2D*[1];
      fhMCEtaMassPtTrue  = new TH2D*[1];
      fhMCEtaPtTruePtRec = new TH2D*[1];
      
      fhMCPi0MassPtTrue[0]  = (TH2D*)  outputList->At(index++);
      fhMCPi0PtTruePtRec[0] = (TH2D*)  outputList->At(index++);
      fhMCEtaMassPtTrue[0]  = (TH2D*)  outputList->At(index++);
      fhMCEtaPtTruePtRec[0] = (TH2D*)  outputList->At(index++);
    }
  }
  
  for(Int_t imod=0; imod < fNModules; imod++){
    fhReMod[imod]                = (TH2D*) outputList->At(index++);
    if(fCalorimeter=="EMCAL"){
      if(imod < fNModules/2) fhReSameSectorEMCALMod[imod] = (TH2D*) outputList->At(index++);
      if(imod < fNModules-2) fhReSameSideEMCALMod[imod]   = (TH2D*) outputList->At(index++);
    }
    else     fhReDiffPHOSMod[imod]        = (TH2D*) outputList->At(index++);

    if(fDoOwnMix){
      fhMiMod[imod]         = (TH2D*) outputList->At(index++);
      if(fCalorimeter=="EMCAL"){
        if(imod < fNModules/2) fhMiSameSectorEMCALMod[imod] = (TH2D*) outputList->At(index++);
        if(imod < fNModules-2) fhMiSameSideEMCALMod[imod]   = (TH2D*) outputList->At(index++);
      }
      else     fhMiDiffPHOSMod[imod]        = (TH2D*) outputList->At(index++);
    }
  }
  
}


//____________________________________________________________________________________________________________________________________________________
void AliAnaPi0::Terminate(TList* outputList) 
{
  //Do some calculations and plots from the final histograms.
  
  printf(" *** %s Terminate:\n", GetName()) ; 
  
  //Recover histograms from output histograms list, needed for distributed analysis.    
  ReadHistograms(outputList);
  
  if (!fhRe1) {
    printf("AliAnaPi0::Terminate() - Error: Remote output histograms not imported in AliAnaPi0 object");
    return;
  }
  
  printf("AliAnaPi0::Terminate()         Mgg Real        : %5.3f , RMS : %5.3f \n", fhRe1[0]->GetMean(),   fhRe1[0]->GetRMS() ) ;
    
  const Int_t buffersize = 255;

  char nameIM[buffersize];
  snprintf(nameIM, buffersize,"AliAnaPi0_%s_cPt",fCalorimeter.Data());
  TCanvas  * cIM = new TCanvas(nameIM, "", 400, 10, 600, 700) ;
  cIM->Divide(2, 2);
  
  cIM->cd(1) ; 
  //gPad->SetLogy();
  TH1D * hIMAllPt = (TH1D*) fhRe1[0]->ProjectionY(Form("IMPtAll_%s",fCalorimeter.Data()));
  hIMAllPt->SetLineColor(2);
  hIMAllPt->SetTitle("No cut on  p_{T, #gamma#gamma} ");
  hIMAllPt->Draw();

  cIM->cd(2) ; 
  TH1D * hIMPt5 = (TH1D*) fhRe1[0]->ProjectionY(Form("IMPt0-5_%s",fCalorimeter.Data()),0, fhRe1[0]->GetXaxis()->FindBin(5.));
//  hRe1Pt5->GetXaxis()->SetRangeUser(0,5);
//  TH1D * hIMPt5 = (TH1D*) hRe1Pt5->Project3D(Form("IMPt5_%s_pz",fCalorimeter.Data()));
  hIMPt5->SetLineColor(2);  
  hIMPt5->SetTitle("0 < p_{T, #gamma#gamma} < 5 GeV/c");
  hIMPt5->Draw();
  
  cIM->cd(3) ; 
  TH1D * hIMPt10 =  (TH1D*) fhRe1[0]->ProjectionY(Form("IMPt5-10_%s",fCalorimeter.Data()), fhRe1[0]->GetXaxis()->FindBin(5.),fhRe1[0]->GetXaxis()->FindBin(10.));
//  hRe1Pt10->GetXaxis()->SetRangeUser(5,10);
//  TH1D * hIMPt10 = (TH1D*) hRe1Pt10->Project3D(Form("IMPt10_%s_pz",fCalorimeter.Data()));
  hIMPt10->SetLineColor(2);  
  hIMPt10->SetTitle("5 < p_{T, #gamma#gamma} < 10 GeV/c");
  hIMPt10->Draw();
  
  cIM->cd(4) ; 
  TH1D * hIMPt20 =  (TH1D*) fhRe1[0]->ProjectionY(Form("IMPt10-20_%s",fCalorimeter.Data()), fhRe1[0]->GetXaxis()->FindBin(10.),fhRe1[0]->GetXaxis()->FindBin(20.));
 // TH3F * hRe1Pt20 =  (TH3F*)fhRe1[0]->Clone(Form("IMPt20_%s",fCalorimeter.Data()));
//  hRe1Pt20->GetXaxis()->SetRangeUser(10,20);
//  TH1D * hIMPt20 = (TH1D*) hRe1Pt20->Project3D(Form("IMPt20_%s_pz",fCalorimeter.Data()));
  hIMPt20->SetLineColor(2);  
  hIMPt20->SetTitle("10 < p_{T, #gamma#gamma} < 20 GeV/c");
  hIMPt20->Draw();
   
  char nameIMF[buffersize];
  snprintf(nameIMF,buffersize,"AliAnaPi0_%s_Mgg.eps",fCalorimeter.Data());
  cIM->Print(nameIMF);

  char namePt[buffersize];
  snprintf(namePt,buffersize,"AliAnaPi0_%s_cPt",fCalorimeter.Data());
  TCanvas  * cPt = new TCanvas(namePt, "", 400, 10, 600, 700) ;
  cPt->Divide(2, 2);

  cPt->cd(1) ; 
  //gPad->SetLogy();
  TH1D * hPt = (TH1D*) fhRe1[0]->ProjectionX(Form("Pt0_%s",fCalorimeter.Data()),-1,-1);
  hPt->SetLineColor(2);
  hPt->SetTitle("No cut on  M_{#gamma#gamma} ");
  hPt->Draw();

  cPt->cd(2) ; 
  TH1D * hPtIM1 = (TH1D*)fhRe1[0]->ProjectionX(Form("Pt1_%s",fCalorimeter.Data()), fhRe1[0]->GetZaxis()->FindBin(0.05),fhRe1[0]->GetZaxis()->FindBin(0.21)); 
//  TH3F * hRe1IM1 = (TH3F*)fhRe1[0]->Clone(Form("Pt1_%s",fCalorimeter.Data()));
//  hRe1IM1->GetZaxis()->SetRangeUser(0.05,0.21);
//  TH1D * hPtIM1 = (TH1D*) hRe1IM1->Project3D("x");
  hPtIM1->SetLineColor(2);  
  hPtIM1->SetTitle("0.05 < M_{#gamma#gamma} < 0.21 GeV/c^{2}");
  hPtIM1->Draw();
  
  cPt->cd(3) ; 
  TH1D * hPtIM2 = (TH1D*)fhRe1[0]->ProjectionX(Form("Pt2_%s",fCalorimeter.Data()), fhRe1[0]->GetZaxis()->FindBin(0.09),fhRe1[0]->GetZaxis()->FindBin(0.17)); 
//  TH3F * hRe1IM2 = (TH3F*)fhRe1[0]->Clone(Form("Pt2_%s",fCalorimeter.Data()));
//  hRe1IM2->GetZaxis()->SetRangeUser(0.09,0.17);
//  TH1D * hPtIM2 = (TH1D*) hRe1IM2->Project3D("x");
  hPtIM2->SetLineColor(2);  
  hPtIM2->SetTitle("0.09 < M_{#gamma#gamma} < 0.17 GeV/c^{2}");
  hPtIM2->Draw();

  cPt->cd(4) ; 
  TH1D * hPtIM3 = (TH1D*)fhRe1[0]->ProjectionX(Form("Pt3_%s",fCalorimeter.Data()), fhRe1[0]->GetZaxis()->FindBin(0.11),fhRe1[0]->GetZaxis()->FindBin(0.15)); 
//  TH3F * hRe1IM3 = (TH3F*)fhRe1[0]->Clone(Form("Pt3_%s",fCalorimeter.Data()));
//  hRe1IM3->GetZaxis()->SetRangeUser(0.11,0.15);
//  TH1D * hPtIM3 = (TH1D*) hRe1IM1->Project3D("x");
  hPtIM3->SetLineColor(2);  
  hPtIM3->SetTitle("0.11 < M_{#gamma#gamma} < 0.15 GeV/c^{2}");
  hPtIM3->Draw();
   
  char namePtF[buffersize];
  snprintf(namePtF,buffersize,"AliAnaPi0_%s_Pt.eps",fCalorimeter.Data());
  cPt->Print(namePtF);

  char line[buffersize] ; 
  snprintf(line,buffersize,".!tar -zcf %s_%s.tar.gz *.eps", GetName(),fCalorimeter.Data()) ; 
  gROOT->ProcessLine(line);
  snprintf(line, buffersize,".!rm -fR AliAnaPi0_%s*.eps",fCalorimeter.Data()); 
  gROOT->ProcessLine(line);
 
  printf(" AliAnaPi0::Terminate() - !! All the eps files are in %s_%s.tar.gz !!!\n", GetName(), fCalorimeter.Data());

}
  //____________________________________________________________________________________________________________________________________________________
Int_t AliAnaPi0::GetEventIndex(AliAODPWG4Particle * part, Double_t * vert)  
{
  // retieves the event index and checks the vertex
  //    in the mixed buffer returns -2 if vertex NOK
  //    for normal events   returns 0 if vertex OK and -1 if vertex NOK
  
  Int_t evtIndex = -1 ; 
  if(GetReader()->GetDataType()!=AliCaloTrackReader::kMC){
    
    if (GetMixedEvent()){
      
      evtIndex = GetMixedEvent()->EventIndexForCaloCluster(part->GetCaloLabel(0)) ;
      GetVertex(vert,evtIndex); 
      
      if(TMath::Abs(vert[2])> GetZvertexCut())
        evtIndex = -2 ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
    } else {// Single event
      
      GetVertex(vert);
      
      if(TMath::Abs(vert[2])> GetZvertexCut())
        evtIndex = -1 ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
      else 
        evtIndex = 0 ;
    }
  }//No MC reader
  else {
    evtIndex = 0;
    vert[0] = 0. ; 
    vert[1] = 0. ; 
    vert[2] = 0. ; 
  }
  
  return evtIndex ; 
}

