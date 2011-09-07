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
// Class to check results from simulations or reconstructed real data. 
// Fill few histograms and do some checking plots
//
//-- Author: Gustavo Conesa (INFN-LNF)
//_________________________________________________________________________


// --- ROOT system ---
//#include "Riostream.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TStyle.h"
#include <TObjString.h>

//---- AliRoot system ----
#include "AliAnaCalorimeterQA.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliAODTrack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODMCParticle.h"
#include "AliMCAnalysisUtils.h"
#include "AliAODPid.h"
#include "AliExternalTrackParam.h"

ClassImp(AliAnaCalorimeterQA)

//____________________________________________________________________________
AliAnaCalorimeterQA::AliAnaCalorimeterQA() : 
AliAnaPartCorrBaseClass(), fCalorimeter(""),           fStyleMacro(""), 
fFillAllPosHisto(kFALSE),  fFillAllPosHisto2(kTRUE), 
fFillAllTH12(kFALSE),      fFillAllTH3(kTRUE), 
fFillAllTMHisto(kTRUE),    fFillAllPi0Histo(kTRUE),
fCorrelate(kTRUE),         fNModules(12),              fNRCU(2),
fTimeCutMin(-1),           fTimeCutMax(9999999),
fEMCALCellAmpMin(0),       fPHOSCellAmpMin(0), 
fhE(0),                    fhPt(0),                    fhPhi(0),                fhEta(0),        fhEtaPhiE(0),
fhECharged(0),             fhPtCharged(0),             fhPhiCharged(0),         fhEtaCharged(0), fhEtaPhiECharged(0), 

//Invariant mass
fhIM(0 ),                  fhIMCellCut(0),             fhAsym(0), 
fhNCellsPerCluster(0),     fhNCellsPerClusterNoCut(0), fhNCellsPerClusterMIP(0),   fhNCellsPerClusterMIPCharged(0), 
fhNCellsvsClusterMaxCellDiffE0(0),    fhNCellsvsClusterMaxCellDiffE2(0),        fhNCellsvsClusterMaxCellDiffE6(0),
fhNClusters(0),    

//Timing
fhClusterTimeEnergy(0),               fhCellTimeSpreadRespectToCellMax(0),  
fhCellIdCellLargeTimeSpread(0),       fhClusterPairDiffTimeE(0),

fhClusterMaxCellCloseCellRatio(0),    fhClusterMaxCellDiff(0),                  fhClusterMaxCellDiffNoCut(0), 
//fhClusterMaxCellDiffDivLambda0(0),
fhLambda0vsClusterMaxCellDiffE0(0),   fhLambda0vsClusterMaxCellDiffE2(0),       fhLambda0vsClusterMaxCellDiffE6(0),

//
//bad cells
fhBadClusterEnergy(0),                fhBadClusterTimeEnergy(0),            fhBadClusterPairDiffTimeE(0),
fhBadClusterMaxCellCloseCellRatio(0), fhBadClusterMaxCellDiff(0),

//Position
fhRNCells(0),              fhXNCells(0),               fhYNCells(0),            fhZNCells(0),
fhRE(0),                   fhXE(0),                    fhYE(0),                 fhZE(0),    
fhXYZ(0),
fhRCellE(0),               fhXCellE(0),                fhYCellE(0),             fhZCellE(0),
fhXYZCell(0),
fhDeltaCellClusterRNCells(0),fhDeltaCellClusterXNCells(0),fhDeltaCellClusterYNCells(0),fhDeltaCellClusterZNCells(0),
fhDeltaCellClusterRE(0),     fhDeltaCellClusterXE(0),     fhDeltaCellClusterYE(0),     fhDeltaCellClusterZE(0),
// Cells
fhNCells(0),               fhAmplitude(0),             fhAmpId(0),              fhEtaPhiAmp(0), 
fhTime(0),                 fhTimeId(0),                fhTimeAmp(0), 
//fhT0Time(0),               fhT0TimeId(0),              fhT0TimeAmp(0), 
fhCaloCorrNClusters(0),    fhCaloCorrEClusters(0),     fhCaloCorrNCells(0),     fhCaloCorrECells(0),
fhCaloV0SCorrNClusters(0), fhCaloV0SCorrEClusters(0),  fhCaloV0SCorrNCells(0),  fhCaloV0SCorrECells(0),
fhCaloV0MCorrNClusters(0), fhCaloV0MCorrEClusters(0),  fhCaloV0MCorrNCells(0),  fhCaloV0MCorrECells(0),
fhCaloTrackMCorrNClusters(0), fhCaloTrackMCorrEClusters(0), fhCaloTrackMCorrNCells(0), fhCaloTrackMCorrECells(0),
//Super-Module dependent histgrams
fhEMod(0),                 fhNClustersMod(0),          fhNCellsPerClusterMod(0),  fhNCellsPerClusterModNoCut(0), fhNCellsMod(0),  
fhGridCellsMod(0),         fhGridCellsEMod(0),         fhGridCellsTimeMod(0), 
fhAmplitudeMod(0),         fhAmplitudeModFraction(0),  fhTimeAmpPerRCU(0), 
//fhT0TimeAmpPerRCU(0),      fhTimeCorrRCU(0),
fhIMMod(0),                fhIMCellCutMod(0),

// MC and reco
fhDeltaE(0),               fhDeltaPt(0),               fhDeltaPhi(0),           fhDeltaEta(0),   
fhRatioE(0),               fhRatioPt(0),               fhRatioPhi(0),           fhRatioEta(0),
fh2E(0),                   fh2Pt(0),                   fh2Phi(0),               fh2Eta(0),

// MC only
fhGenGamPt(0),             fhGenGamEta(0),             fhGenGamPhi(0),
fhGenPi0Pt(0),             fhGenPi0Eta(0),             fhGenPi0Phi(0),
fhGenEtaPt(0),             fhGenEtaEta(0),             fhGenEtaPhi(0),
fhGenOmegaPt(0),           fhGenOmegaEta(0),           fhGenOmegaPhi(0),
fhGenElePt(0),             fhGenEleEta(0),             fhGenElePhi(0), 
fhEMVxyz(0),               fhEMR(0),                   fhHaVxyz(0),             fhHaR(0),
fhGamE(0),                 fhGamPt(0),                 fhGamPhi(0),             fhGamEta(0), 
fhGamDeltaE(0),            fhGamDeltaPt(0),            fhGamDeltaPhi(0),        fhGamDeltaEta(0), 
fhGamRatioE(0),            fhGamRatioPt(0),            fhGamRatioPhi(0),        fhGamRatioEta(0),
fhEleE(0),                 fhElePt(0),                 fhElePhi(0),             fhEleEta(0),
fhPi0E(0),                 fhPi0Pt(0),                 fhPi0Phi(0),             fhPi0Eta(0), 
fhNeHadE(0),               fhNeHadPt(0),               fhNeHadPhi(0),           fhNeHadEta(0), 
fhChHadE(0),               fhChHadPt(0),               fhChHadPhi(0),           fhChHadEta(0),
fhGamECharged(0),          fhGamPtCharged(0),          fhGamPhiCharged(0),      fhGamEtaCharged(0), 
fhEleECharged(0),          fhElePtCharged(0),          fhElePhiCharged(0),      fhEleEtaCharged(0),
fhPi0ECharged(0),          fhPi0PtCharged(0),          fhPi0PhiCharged(0),      fhPi0EtaCharged(0), 
fhNeHadECharged(0),        fhNeHadPtCharged(0),        fhNeHadPhiCharged(0),    fhNeHadEtaCharged(0), 
fhChHadECharged(0),        fhChHadPtCharged(0),        fhChHadPhiCharged(0),    fhChHadEtaCharged(0),
fhGenGamAccE(0),           fhGenGamAccPt(0),           fhGenGamAccEta(0),       fhGenGamAccPhi(0),
fhGenPi0AccE(0),           fhGenPi0AccPt(0),           fhGenPi0AccEta(0),       fhGenPi0AccPhi(0),
fh1pOverE(0),              fh1dR(0),                   fh2EledEdx(0),           fh2MatchdEdx(0),
fhMCEle1pOverE(0),         fhMCEle1dR(0),              fhMCEle2MatchdEdx(0),
fhMCChHad1pOverE(0),       fhMCChHad1dR(0),            fhMCChHad2MatchdEdx(0),
fhMCNeutral1pOverE(0),     fhMCNeutral1dR(0),          fhMCNeutral2MatchdEdx(0),fh1pOverER02(0),           
fhMCEle1pOverER02(0),      fhMCChHad1pOverER02(0),     fhMCNeutral1pOverER02(0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
}

//________________________________________________________________________
TObjString *  AliAnaCalorimeterQA::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaCalorimeterQA ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Time Cut : %2.2f < T < %2.2f ns  \n",fTimeCutMin, fTimeCutMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"PHOS Cell Amplitude > %2.2f GeV, EMCAL Cell Amplitude > %2.2f GeV  \n",fPHOSCellAmpMin, fEMCALCellAmpMin) ;
  parList+=onePar ;
  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
	
  return new TObjString(parList) ;
}


//________________________________________________________________________
TList *  AliAnaCalorimeterQA::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("QAHistos") ; 
  
  //Histograms
  Int_t nptbins     = GetHistoPtBins(); 	        Float_t ptmax     = GetHistoPtMax();           Float_t ptmin     = GetHistoPtMin();
  Int_t nfineptbins = GetHistoFinePtBins(); 	    Float_t ptfinemax = GetHistoFinePtMax();       Float_t ptfinemin = GetHistoFinePtMin();
  Int_t nphibins    = GetHistoPhiBins();     	    Float_t phimax    = GetHistoPhiMax();          Float_t phimin    = GetHistoPhiMin();
  Int_t netabins    = GetHistoEtaBins();          Float_t etamax    = GetHistoEtaMax();          Float_t etamin    = GetHistoEtaMin();	
  Int_t nmassbins   = GetHistoMassBins();         Float_t massmax   = GetHistoMassMax(); 	       Float_t massmin   = GetHistoMassMin();
  Int_t nasymbins   = GetHistoAsymmetryBins();    Float_t asymmax   = GetHistoAsymmetryMax();    Float_t asymmin   = GetHistoAsymmetryMin();
  Int_t nPoverEbins = GetHistoPOverEBins();       Float_t pOverEmax = GetHistoPOverEMax();       Float_t pOverEmin = GetHistoPOverEMin();
  Int_t ndedxbins   = GetHistodEdxBins();         Float_t dedxmax   = GetHistodEdxMax();         Float_t dedxmin   = GetHistodEdxMin();
  Int_t ndRbins     = GetHistodRBins();           Float_t dRmax     = GetHistodRMax();           Float_t dRmin     = GetHistodRMin();
  Int_t ntimebins   = GetHistoTimeBins();         Float_t timemax   = GetHistoTimeMax();         Float_t timemin   = GetHistoTimeMin();       
  Int_t nbins       = GetHistoNClusterCellBins(); Int_t   nmax      = GetHistoNClusterCellMax(); Int_t   nmin      = GetHistoNClusterCellMin(); 
  Int_t nratiobins  = GetHistoRatioBins();        Float_t ratiomax  = GetHistoRatioMax();        Float_t ratiomin  = GetHistoRatioMin();
  Int_t nvdistbins  = GetHistoVertexDistBins();   Float_t vdistmax  = GetHistoVertexDistMax();   Float_t vdistmin  = GetHistoVertexDistMin();
  Int_t rbins       = GetHistoRBins();            Float_t rmax      = GetHistoRMax();            Float_t rmin      = GetHistoRMin(); 
  Int_t xbins       = GetHistoXBins();            Float_t xmax      = GetHistoXMax();            Float_t xmin      = GetHistoXMin(); 
  Int_t ybins       = GetHistoYBins();            Float_t ymax      = GetHistoYMax();            Float_t ymin      = GetHistoYMin(); 
  Int_t zbins       = GetHistoZBins();            Float_t zmax      = GetHistoZMax();            Float_t zmin      = GetHistoZMin(); 
  Int_t ssbins      = GetHistoShowerShapeBins();  Float_t ssmax     = GetHistoShowerShapeMax();  Float_t ssmin     = GetHistoShowerShapeMin();
  Int_t tdbins      = GetHistoDiffTimeBins() ;    Float_t tdmax     = GetHistoDiffTimeMax();     Float_t tdmin     = GetHistoDiffTimeMin();

  Int_t nv0sbins    = GetHistoV0SignalBins();          Int_t nv0smax = GetHistoV0SignalMax();          Int_t nv0smin = GetHistoV0SignalMin(); 
  Int_t nv0mbins    = GetHistoV0MultiplicityBins();    Int_t nv0mmax = GetHistoV0MultiplicityMax();    Int_t nv0mmin = GetHistoV0MultiplicityMin(); 
  Int_t ntrmbins    = GetHistoTrackMultiplicityBins(); Int_t ntrmmax = GetHistoTrackMultiplicityMax(); Int_t ntrmmin = GetHistoTrackMultiplicityMin(); 
  
  
  //EMCAL
  Int_t colmax = 48;
  Int_t rowmax = 24;
  fNRCU   = 2 ;
  //PHOS
  if(fCalorimeter=="PHOS"){
    colmax = 56;
    rowmax = 64;
    fNRCU   = 4 ;
  }
  
  
  fhE  = new TH1F ("hE","E reconstructed clusters ", nptbins*5,ptmin,ptmax*5);  
  fhE->SetXTitle("E (GeV)");
  outputContainer->Add(fhE);
  
  if(fFillAllTH12){
    fhPt  = new TH1F ("hPt","p_{T} reconstructed clusters", nptbins,ptmin,ptmax); 
    fhPt->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPt);
    
    fhPhi  = new TH1F ("hPhi","#phi reconstructed clusters ",nphibins,phimin,phimax); 
    fhPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhi);
    
    fhEta  = new TH1F ("hEta","#eta reconstructed clusters ",netabins,etamin,etamax); 
    fhEta->SetXTitle("#eta ");
    outputContainer->Add(fhEta);
  }
  
  fhEtaPhiE  = new TH3F ("hEtaPhiE","#eta vs #phi vs energy, reconstructed clusters",
                         netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
  fhEtaPhiE->SetXTitle("#eta ");
  fhEtaPhiE->SetYTitle("#phi (rad)");
  fhEtaPhiE->SetZTitle("E (GeV) ");
  outputContainer->Add(fhEtaPhiE);
  
  fhClusterTimeEnergy  = new TH2F ("hClusterTimeEnergy","energy vs TOF, reconstructed clusters",
                                   nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
  fhClusterTimeEnergy->SetXTitle("E (GeV) ");
  fhClusterTimeEnergy->SetYTitle("TOF (ns)");
  outputContainer->Add(fhClusterTimeEnergy);
    
  fhClusterPairDiffTimeE = new TH2F("hClusterPairDiffTimeE","cluster pair time difference vs E, only good clusters",
                                    nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhClusterPairDiffTimeE->SetXTitle("E_{cluster} (GeV)");
  fhClusterPairDiffTimeE->SetYTitle("#Delta t (ns)");
  outputContainer->Add(fhClusterPairDiffTimeE);  
  
  
  fhClusterMaxCellCloseCellRatio  = new TH2F ("hClusterMaxCellCloseCell","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                              nptbins,ptmin,ptmax, 100,0,1.); 
  fhClusterMaxCellCloseCellRatio->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellCloseCellRatio->SetYTitle("ratio");
  outputContainer->Add(fhClusterMaxCellCloseCellRatio);
  
  fhClusterMaxCellDiff  = new TH2F ("hClusterMaxCellDiff","energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",
                                       nptbins,ptmin,ptmax, 500,0,1.); 
  fhClusterMaxCellDiff->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellDiff->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  outputContainer->Add(fhClusterMaxCellDiff);  

  fhClusterMaxCellDiffNoCut  = new TH2F ("hClusterMaxCellDiffNoCut","energy vs difference of cluster energy - max cell energy / cluster energy",
                                    nptbins,ptmin,ptmax, 500,0,1.); 
  fhClusterMaxCellDiffNoCut->SetXTitle("E_{cluster} (GeV) ");
  fhClusterMaxCellDiffNoCut->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  outputContainer->Add(fhClusterMaxCellDiffNoCut);  
  
//  fhClusterMaxCellDiffDivLambda0  = new TH2F ("hClusterMaxCellDiffDivLambda0;","",
//                                         nptbins,ptmin,ptmax, 500,0,5.); 
//  fhClusterMaxCellDiffDivLambda0->SetXTitle("E_{cluster} (GeV) ");
//  fhClusterMaxCellDiffDivLambda0->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster} / #lambda_{0}");
//  outputContainer->Add(fhClusterMaxCellDiffDivLambda0);    
  
  fhLambda0vsClusterMaxCellDiffE0  = new TH2F ("hLambda0vsClusterMaxCellDiffE0","shower shape, #lambda^{2}_{0} vs fraction of energy carried by max cell, E < 2 GeV ",
                                               ssbins,ssmin,ssmax,500,0,1.); 
  fhLambda0vsClusterMaxCellDiffE0->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhLambda0vsClusterMaxCellDiffE0->SetXTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0vsClusterMaxCellDiffE0); 
  
  fhLambda0vsClusterMaxCellDiffE2  = new TH2F ("hLambda0vsClusterMaxCellDiffE2","shower shape, #lambda^{2}_{0} vs fraction of energy carried by max cell, 2 < E < 6 GeV ",
                                               ssbins,ssmin,ssmax,500,0,1.); 
  fhLambda0vsClusterMaxCellDiffE2->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhLambda0vsClusterMaxCellDiffE2->SetXTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0vsClusterMaxCellDiffE2); 
  
  fhLambda0vsClusterMaxCellDiffE6  = new TH2F ("hLambda0vsClusterMaxCellDiffE6","shower shape, #lambda^{2}_{0} vs fraction of energy carried by max cell, E > 6 ",
                                               ssbins,ssmin,ssmax,500,0,1.); 
  fhLambda0vsClusterMaxCellDiffE6->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhLambda0vsClusterMaxCellDiffE6->SetXTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0vsClusterMaxCellDiffE6); 

  fhNCellsvsClusterMaxCellDiffE0  = new TH2F ("hNCellsvsClusterMaxCellDiffE0","N cells per cluster vs fraction of energy carried by max cell, E < 2 GeV ",
                                               nbins/5,nmin,nmax/5,500,0,1.); 
  fhNCellsvsClusterMaxCellDiffE0->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhNCellsvsClusterMaxCellDiffE0->SetXTitle("N cells per cluster");
  outputContainer->Add(fhNCellsvsClusterMaxCellDiffE0); 
  
  fhNCellsvsClusterMaxCellDiffE2  = new TH2F ("hNCellsvsClusterMaxCellDiffE2","N cells per cluster vs fraction of energy carried by max cell, 2 < E < 6 GeV ",
                                               nbins/5,nmin,nmax/5,500,0,1.); 
  fhNCellsvsClusterMaxCellDiffE2->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhNCellsvsClusterMaxCellDiffE2->SetXTitle("N cells per cluster");
  outputContainer->Add(fhNCellsvsClusterMaxCellDiffE2); 
  
  fhNCellsvsClusterMaxCellDiffE6  = new TH2F ("hNCellsvsClusterMaxCellDiffE6","N cells per cluster vs fraction of energy carried by max cell, E > 6 ",
                                               nbins/5,nmin,nmax/5,500,0,1.); 
  fhNCellsvsClusterMaxCellDiffE6->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  fhNCellsvsClusterMaxCellDiffE6->SetXTitle("N cells per cluster");
  outputContainer->Add(fhNCellsvsClusterMaxCellDiffE6); 
  
  
  if(fCalorimeter=="EMCAL" && !GetCaloUtils()->GetEMCALRecoUtils()->IsRejectExoticCluster()){
    
    fhBadClusterEnergy  = new TH1F ("hBadClusterEnergy","Bad cluster energy", nptbins,ptmin,ptmax); 
    fhBadClusterEnergy->SetXTitle("E_{cluster} (GeV) ");
    outputContainer->Add(fhBadClusterEnergy);
    
    fhBadClusterMaxCellCloseCellRatio  = new TH2F ("hBadClusterMaxCellCloseCell","energy vs ratio of max cell / neighbour cell constributing cell, reconstructed bad clusters",
                                                   nptbins,ptmin,ptmax, 100,0,1.); 
    fhBadClusterMaxCellCloseCellRatio->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterMaxCellCloseCellRatio->SetYTitle("ratio");
    outputContainer->Add(fhBadClusterMaxCellCloseCellRatio);
        
    fhBadClusterMaxCellDiff  = new TH2F ("hBadClusterMaxCellDiff","energy vs difference of cluster energy - max cell energy / cluster energy for bad clusters",
                                                   nptbins,ptmin,ptmax, 500,0,1.); 
    fhBadClusterMaxCellDiff->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterMaxCellDiff->SetYTitle("(E_{cluster} - E_{cell max}) / E_{cluster}");
    outputContainer->Add(fhBadClusterMaxCellDiff);
    
    fhBadClusterTimeEnergy  = new TH2F ("hBadClusterTimeEnergy","energy vs TOF of reconstructed bad clusters",
                                               nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhBadClusterTimeEnergy->SetXTitle("E_{cluster} (GeV) ");
    fhBadClusterTimeEnergy->SetYTitle("TOF (ns)");
    outputContainer->Add(fhBadClusterTimeEnergy);    

    fhBadClusterPairDiffTimeE = new TH2F("hBadClusterPairDiffTimeE","cluster pair time difference (bad - good) vs E from bad cluster",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhBadClusterPairDiffTimeE->SetXTitle("E_{bad cluster} (GeV)");
    fhBadClusterPairDiffTimeE->SetYTitle("#Delta t (ns)");
    outputContainer->Add(fhBadClusterPairDiffTimeE);    
      
  }
  
  //Track Matching
  if(fFillAllTMHisto){
    if(fFillAllTH12){
      fhECharged  = new TH1F ("hECharged","E reconstructed clusters, matched with track", nptbins,ptmin,ptmax); 
      fhECharged->SetXTitle("E (GeV)");
      outputContainer->Add(fhECharged);
      
      fhPtCharged  = new TH1F ("hPtCharged","p_{T} reconstructed clusters, matched with track", nptbins,ptmin,ptmax); 
      fhPtCharged->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtCharged);
      
      fhPhiCharged  = new TH1F ("hPhiCharged","#phi reconstructed clusters, matched with track",nphibins,phimin,phimax); 
      fhPhiCharged->SetXTitle("#phi (rad)");
      outputContainer->Add(fhPhiCharged);
      
      fhEtaCharged  = new TH1F ("hEtaCharged","#eta reconstructed clusters, matched with track",netabins,etamin,etamax); 
      fhEtaCharged->SetXTitle("#eta ");
      outputContainer->Add(fhEtaCharged);
    }
    if(fFillAllTH3){
      fhEtaPhiECharged  = new TH3F ("hEtaPhiECharged","#eta vs #phi, reconstructed clusters, matched with track",
                                    netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
      fhEtaPhiECharged->SetXTitle("#eta ");
      fhEtaPhiECharged->SetYTitle("#phi ");
      fhEtaPhiECharged->SetZTitle("E (GeV) ");
      outputContainer->Add(fhEtaPhiECharged);	
    }
    
    fh1pOverE = new TH2F("h1pOverE","TRACK matches p/E",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fh1pOverE->SetYTitle("p/E");
    fh1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fh1pOverE);
    
    fh1dR = new TH1F("h1dR","TRACK matches dR",ndRbins,dRmin,dRmax);
    fh1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fh1dR) ;
    
    fh2MatchdEdx = new TH2F("h2MatchdEdx","dE/dx vs. p for all matches",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fh2MatchdEdx->SetXTitle("p (GeV/c)");
    fh2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fh2MatchdEdx);
    
    fh2EledEdx = new TH2F("h2EledEdx","dE/dx vs. p for electrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fh2EledEdx->SetXTitle("p (GeV/c)");
    fh2EledEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fh2EledEdx) ;
    
    fh1pOverER02 = new TH2F("h1pOverER02","TRACK matches p/E, all",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fh1pOverER02->SetYTitle("p/E");
    fh1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fh1pOverER02);	
  }
  
  if(fFillAllPi0Histo){
    fhIM  = new TH2F ("hIM","Cluster pairs Invariant mass vs reconstructed pair energy",nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
    fhIM->SetXTitle("p_{T, cluster pairs} (GeV) ");
    fhIM->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
    outputContainer->Add(fhIM);
    
    fhIMCellCut  = new TH2F ("hIMCellCut","Cluster (n cell > 1) pairs Invariant mass vs reconstructed pair energy",nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
    fhIMCellCut->SetXTitle("p_{T, cluster pairs} (GeV) ");
    fhIMCellCut->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
    outputContainer->Add(fhIMCellCut);
    
    fhAsym  = new TH2F ("hAssym","Cluster pairs Asymmetry vs reconstructed pair energy",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax); 
    fhAsym->SetXTitle("p_{T, cluster pairs} (GeV) ");
    fhAsym->SetYTitle("Asymmetry");
    outputContainer->Add(fhAsym);	
    
  }

  fhNCellsPerClusterNoCut  = new TH2F ("hNCellsPerClusterNoCut","# cells per cluster vs energy vs #eta, no bad clusters cut",nptbins,ptmin,ptmax, nbins/5,nmin,nmax/5); 
  fhNCellsPerClusterNoCut->SetXTitle("E (GeV)");
  fhNCellsPerClusterNoCut->SetYTitle("n cells");
  outputContainer->Add(fhNCellsPerClusterNoCut);
    
  fhNCellsPerCluster  = new TH2F ("hNCellsPerCluster","# cells per cluster vs energy vs #eta",nptbins,ptmin,ptmax, nbins/5,nmin,nmax/5); 
  fhNCellsPerCluster->SetXTitle("E (GeV)");
  fhNCellsPerCluster->SetYTitle("n cells");
  outputContainer->Add(fhNCellsPerCluster);
    
  if((fCalorimeter=="EMCAL" && GetReader()->GetEMCALPtMin() < 0.3) ||
     (fCalorimeter=="PHOS"  && GetReader()->GetPHOSPtMin()  < 0.3)) {
    fhNCellsPerClusterMIP  = new TH2F ("hNCellsPerClusterMIP","# cells per cluster vs energy vs #eta, smaller bin for MIP search", 
                                       40,0.,2., 11,0,10); 
    fhNCellsPerClusterMIP->SetXTitle("E (GeV)");
    fhNCellsPerClusterMIP->SetYTitle("n cells");
    outputContainer->Add(fhNCellsPerClusterMIP);
    
    
    if(fFillAllTMHisto){
      fhNCellsPerClusterMIPCharged  = new TH2F ("hNCellsPerClusterMIPCharged","# cells per track-matched cluster vs energy vs #eta, smaller bin for MIP search", 
                                                40,0.,2., 11,0,10); 
      fhNCellsPerClusterMIPCharged->SetXTitle("E (GeV)");
      fhNCellsPerClusterMIPCharged->SetYTitle("n cells");
      outputContainer->Add(fhNCellsPerClusterMIPCharged);
    }
	}
  
  fhNClusters  = new TH1F ("hNClusters","# clusters", nbins,nmin,nmax); 
  fhNClusters->SetXTitle("number of clusters");
  outputContainer->Add(fhNClusters);
  
  if(fFillAllPosHisto2){
    
    if(fFillAllTH3){
      fhXYZ  = new TH3F ("hXYZ","Cluster: x vs y vs z",xbins,xmin,xmax,ybins,ymin,ymax,zbins,zmin,zmax); 
      fhXYZ->SetXTitle("x (cm)");
      fhXYZ->SetYTitle("y (cm)");
      fhXYZ->SetZTitle("z (cm) ");
      outputContainer->Add(fhXYZ);  
    }
    
    fhXNCells  = new TH2F ("hXNCells","Cluster X position vs N Clusters per Cell",xbins,xmin,xmax,nbins,nmin,nmax); 
    fhXNCells->SetXTitle("x (cm)");
    fhXNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhXNCells);
    
    fhZNCells  = new TH2F ("hZNCells","Cluster Z position vs N Clusters per Cell",zbins,zmin,zmax,nbins,nmin,nmax); 
    fhZNCells->SetXTitle("z (cm)");
    fhZNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhZNCells);
    
    fhXE  = new TH2F ("hXE","Cluster X position vs cluster energy",xbins,xmin,xmax,nptbins,ptmin,ptmax); 
    fhXE->SetXTitle("x (cm)");
    fhXE->SetYTitle("E (GeV)");
    outputContainer->Add(fhXE);
    
    fhZE  = new TH2F ("hZE","Cluster Z position vs cluster energy",zbins,zmin,zmax,nptbins,ptmin,ptmax); 
    fhZE->SetXTitle("z (cm)");
    fhZE->SetYTitle("E (GeV)");
    outputContainer->Add(fhZE);    
    
    
    fhRNCells  = new TH2F ("hRNCells","Cluster R position vs N Clusters per Cell",rbins,rmin,rmax,nbins,nmin,nmax); 
    fhRNCells->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhRNCells);
    
    
    fhYNCells  = new TH2F ("hYNCells","Cluster Y position vs N Clusters per Cell",ybins,ymin,ymax,nbins,nmin,nmax); 
    fhYNCells->SetXTitle("y (cm)");
    fhYNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhYNCells);
    
    fhRE  = new TH2F ("hRE","Cluster R position vs cluster energy",rbins,rmin,rmax,nptbins,ptmin,ptmax); 
    fhRE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRE->SetYTitle("E (GeV)");
    outputContainer->Add(fhRE);
    
    fhYE  = new TH2F ("hYE","Cluster Y position vs cluster energy",ybins,ymin,ymax,nptbins,ptmin,ptmax); 
    fhYE->SetXTitle("y (cm)");
    fhYE->SetYTitle("E (GeV)");
    outputContainer->Add(fhYE);
  }
  if(fFillAllPosHisto){
    
    fhRCellE  = new TH2F ("hRCellE","Cell R position vs cell energy",rbins,rmin,rmax,nptbins,ptmin,ptmax); 
    fhRCellE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhRCellE);
    
    fhXCellE  = new TH2F ("hXCellE","Cell X position vs cell energy",xbins,xmin,xmax,nptbins,ptmin,ptmax); 
    fhXCellE->SetXTitle("x (cm)");
    fhXCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhXCellE);
    
    fhYCellE  = new TH2F ("hYCellE","Cell Y position vs cell energy",ybins,ymin,ymax,nptbins,ptmin,ptmax); 
    fhYCellE->SetXTitle("y (cm)");
    fhYCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhYCellE);
    
    fhZCellE  = new TH2F ("hZCellE","Cell Z position vs cell energy",zbins,zmin,zmax,nptbins,ptmin,ptmax); 
    fhZCellE->SetXTitle("z (cm)");
    fhZCellE->SetYTitle("E (GeV)");
    outputContainer->Add(fhZCellE);
    
    fhXYZCell  = new TH3F ("hXYZCell","Cell : x vs y vs z",xbins,xmin,xmax,ybins,ymin,ymax,zbins,zmin,zmax); 
    fhXYZCell->SetXTitle("x (cm)");
    fhXYZCell->SetYTitle("y (cm)");
    fhXYZCell->SetZTitle("z (cm)");
    outputContainer->Add(fhXYZCell);
    
    
    Float_t dx = TMath::Abs(xmin)+TMath::Abs(xmax);
    Float_t dy = TMath::Abs(ymin)+TMath::Abs(ymax);
    Float_t dz = TMath::Abs(zmin)+TMath::Abs(zmax);
    Float_t dr = TMath::Abs(rmin)+TMath::Abs(rmax);
    
    fhDeltaCellClusterRNCells  = new TH2F ("hDeltaCellClusterRNCells","Cluster-Cell R position vs N Clusters per Cell",rbins*2,-dr,dr,nbins,nmin,nmax); 
    fhDeltaCellClusterRNCells->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhDeltaCellClusterRNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterRNCells);
    
    fhDeltaCellClusterXNCells  = new TH2F ("hDeltaCellClusterXNCells","Cluster-Cell X position vs N Clusters per Cell",xbins*2,-dx,dx,nbins,nmin,nmax); 
    fhDeltaCellClusterXNCells->SetXTitle("x (cm)");
    fhDeltaCellClusterXNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterXNCells);
    
    fhDeltaCellClusterYNCells  = new TH2F ("hDeltaCellClusterYNCells","Cluster-Cell Y position vs N Clusters per Cell",ybins*2,-dy,dy,nbins,nmin,nmax); 
    fhDeltaCellClusterYNCells->SetXTitle("y (cm)");
    fhDeltaCellClusterYNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterYNCells);
    
    fhDeltaCellClusterZNCells  = new TH2F ("hDeltaCellClusterZNCells","Cluster-Cell Z position vs N Clusters per Cell",zbins*2,-dz,dz,nbins,nmin,nmax); 
    fhDeltaCellClusterZNCells->SetXTitle("z (cm)");
    fhDeltaCellClusterZNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterZNCells);
    
    fhDeltaCellClusterRE  = new TH2F ("hDeltaCellClusterRE","Cluster-Cell R position vs cluster energy",rbins*2,-dr,dr,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterRE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhDeltaCellClusterRE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterRE);		
    
    fhDeltaCellClusterXE  = new TH2F ("hDeltaCellClusterXE","Cluster-Cell X position vs cluster energy",xbins*2,-dx,dx,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterXE->SetXTitle("x (cm)");
    fhDeltaCellClusterXE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterXE);
    
    fhDeltaCellClusterYE  = new TH2F ("hDeltaCellClusterYE","Cluster-Cell Y position vs cluster energy",ybins*2,-dy,dy,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterYE->SetXTitle("y (cm)");
    fhDeltaCellClusterYE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterYE);
    
    fhDeltaCellClusterZE  = new TH2F ("hDeltaCellClusterZE","Cluster-Cell Z position vs cluster energy",zbins*2,-dz,dz,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterZE->SetXTitle("z (cm)");
    fhDeltaCellClusterZE->SetYTitle("E (GeV)");
    outputContainer->Add(fhDeltaCellClusterZE);
    
    fhEtaPhiAmp  = new TH3F ("hEtaPhiAmp","Cell #eta vs cell #phi vs cell energy",netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
    fhEtaPhiAmp->SetXTitle("#eta ");
    fhEtaPhiAmp->SetYTitle("#phi (rad)");
    fhEtaPhiAmp->SetZTitle("E (GeV) ");
    outputContainer->Add(fhEtaPhiAmp);		
    
  }
  
  //Calo cells
  fhNCells  = new TH1F ("hNCells","# cells", colmax*rowmax*fNModules,0,colmax*rowmax*fNModules); 
  fhNCells->SetXTitle("n cells");
  outputContainer->Add(fhNCells);
  
  fhAmplitude  = new TH1F ("hAmplitude","Cell Energy", nptbins*2,ptmin,ptmax); 
  fhAmplitude->SetXTitle("Cell Energy (GeV)");
  outputContainer->Add(fhAmplitude);
  
  fhAmpId  = new TH2F ("hAmpId","Cell Energy", nfineptbins,ptfinemin,ptfinemax,rowmax*colmax*fNModules,0,rowmax*colmax*fNModules); 
  fhAmpId->SetXTitle("Cell Energy (GeV)");
  outputContainer->Add(fhAmpId);
  
  
  //Cell Time histograms, time only available in ESDs
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
    
    fhCellTimeSpreadRespectToCellMax = new TH1F ("hCellTimeSpreadRespectToCellMax","t_{cell max}-t_{cell i} per cluster", 100,-200,200); 
    fhCellTimeSpreadRespectToCellMax->SetXTitle("#Delta t (ns)");
    outputContainer->Add(fhCellTimeSpreadRespectToCellMax);
    
    fhCellIdCellLargeTimeSpread= new TH1F ("hCellIdCellLargeTimeSpread","", colmax*rowmax*fNModules,0,colmax*rowmax*fNModules); 
    fhCellIdCellLargeTimeSpread->SetXTitle("Absolute Cell Id");
    outputContainer->Add(fhCellIdCellLargeTimeSpread);
    
    fhTime  = new TH1F ("hTime","Cell Time",ntimebins,timemin,timemax); 
    fhTime->SetXTitle("Cell Time (ns)");
    outputContainer->Add(fhTime);
    
    fhTimeId  = new TH2F ("hTimeId","Cell Time vs Absolute Id",ntimebins,timemin,timemax,rowmax*colmax*fNModules,0,rowmax*colmax*fNModules); 
    fhTimeId->SetXTitle("Cell Time (ns)");
    fhTimeId->SetYTitle("Cell Absolute Id");
    outputContainer->Add(fhTimeId);
    
    fhTimeAmp  = new TH2F ("hTimeAmp","Cell Time vs Cell Energy",nptbins*2,ptmin,ptmax,ntimebins,timemin,timemax); 
    fhTimeAmp->SetYTitle("Cell Time (ns)");
    fhTimeAmp->SetXTitle("Cell Energy (GeV)");
    outputContainer->Add(fhTimeAmp);
    
    //		fhT0Time  = new TH1F ("hT0Time","Cell Time",ntimebins,timemin,timemax); 
    //		fhT0Time->SetXTitle("T_{0} - T_{EMCal} (ns)");
    //		outputContainer->Add(fhT0Time);
    //		
    //		fhT0TimeId  = new TH2F ("hT0TimeId","Cell Time vs Absolute Id",ntimebins,timemin,timemax,rowmax*colmax*fNModules,0,rowmax*colmax*fNModules); 
    //		fhT0TimeId->SetXTitle("T_{0} - T_{EMCal} (ns)");
    //		fhT0TimeId->SetYTitle("Cell Absolute Id");
    //		outputContainer->Add(fhT0TimeId);
    //		
    //		fhT0TimeAmp  = new TH2F ("hT0TimeAmp","Cell Time vs Cell Energy",nptbins*2,ptmin,ptmax,ntimebins,timemin,timemax); 
    //		fhT0TimeAmp->SetYTitle("T_{0} - T_{EMCal} (ns)");
    //		fhT0TimeAmp->SetXTitle("Cell Energy (GeV)");
    //		outputContainer->Add(fhT0TimeAmp);
  }
	
  if(fCorrelate){
    //PHOS vs EMCAL
    fhCaloCorrNClusters  = new TH2F ("hCaloCorrNClusters","# clusters in EMCAL vs PHOS", nbins,nmin,nmax,nbins,nmin,nmax); 
    fhCaloCorrNClusters->SetXTitle("number of clusters in EMCAL");
    fhCaloCorrNClusters->SetYTitle("number of clusters in PHOS");
    outputContainer->Add(fhCaloCorrNClusters);
    
    fhCaloCorrEClusters  = new TH2F ("hCaloCorrEClusters","summed energy of clusters in EMCAL vs PHOS", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fhCaloCorrEClusters->SetXTitle("#Sigma E of clusters in EMCAL (GeV)");
    fhCaloCorrEClusters->SetYTitle("#Sigma E of clusters in PHOS (GeV)");
    outputContainer->Add(fhCaloCorrEClusters);
    
    fhCaloCorrNCells  = new TH2F ("hCaloCorrNCells","# Cells in EMCAL vs PHOS", nbins,nmin,nmax, nbins,nmin,nmax); 
    fhCaloCorrNCells->SetXTitle("number of Cells in EMCAL");
    fhCaloCorrNCells->SetYTitle("number of Cells in PHOS");
    outputContainer->Add(fhCaloCorrNCells);
    
    fhCaloCorrECells  = new TH2F ("hCaloCorrECells","summed energy of Cells in EMCAL vs PHOS", nptbins*2,ptmin,ptmax*2,nptbins*2,ptmin,ptmax*2); 
    fhCaloCorrECells->SetXTitle("#Sigma E of Cells in EMCAL (GeV)");
    fhCaloCorrECells->SetYTitle("#Sigma E of Cells in PHOS (GeV)");
    outputContainer->Add(fhCaloCorrECells);
    
    //Calorimeter VS V0 signal
    fhCaloV0SCorrNClusters  = new TH2F ("hCaloV0SNClusters",Form("# clusters in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax,nbins,nmin,nmax); 
    fhCaloV0SCorrNClusters->SetXTitle("V0 signal");
    fhCaloV0SCorrNClusters->SetYTitle(Form("number of clusters in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrNClusters);
    
    fhCaloV0SCorrEClusters  = new TH2F ("hCaloV0SEClusters",Form("summed energy of clusters in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax,nptbins,ptmin,ptmax); 
    fhCaloV0SCorrEClusters->SetXTitle("V0 signal");
    fhCaloV0SCorrEClusters->SetYTitle(Form("#Sigma E of clusters in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrEClusters);
    
    fhCaloV0SCorrNCells  = new TH2F ("hCaloV0SNCells",Form("# Cells in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax, nbins,nmin,nmax); 
    fhCaloV0SCorrNCells->SetXTitle("V0 signal");
    fhCaloV0SCorrNCells->SetYTitle(Form("number of Cells in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrNCells);
    
    fhCaloV0SCorrECells  = new TH2F ("hCaloV0SECells",Form("summed energy of Cells in %s vs V0 signal",fCalorimeter.Data()), nv0sbins,nv0smin,nv0smax,nptbins,ptmin,ptmax); 
    fhCaloV0SCorrECells->SetXTitle("V0 signal");
    fhCaloV0SCorrECells->SetYTitle(Form("#Sigma E of Cells in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0SCorrECells);    
    
    //Calorimeter VS V0 multiplicity
    fhCaloV0MCorrNClusters  = new TH2F ("hCaloV0MNClusters",Form("# clusters in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax,nbins,nmin,nmax); 
    fhCaloV0MCorrNClusters->SetXTitle("V0 signal");
    fhCaloV0MCorrNClusters->SetYTitle(Form("number of clusters in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrNClusters);
    
    fhCaloV0MCorrEClusters  = new TH2F ("hCaloV0MEClusters",Form("summed energy of clusters in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax,nptbins,ptmin,ptmax); 
    fhCaloV0MCorrEClusters->SetXTitle("V0 signal");
    fhCaloV0MCorrEClusters->SetYTitle(Form("#Sigma E of clusters in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrEClusters);
    
    fhCaloV0MCorrNCells  = new TH2F ("hCaloV0MNCells",Form("# Cells in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax, nbins,nmin,nmax); 
    fhCaloV0MCorrNCells->SetXTitle("V0 signal");
    fhCaloV0MCorrNCells->SetYTitle(Form("number of Cells in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrNCells);
    
    fhCaloV0MCorrECells  = new TH2F ("hCaloV0MECells",Form("summed energy of Cells in %s vs V0 signal",fCalorimeter.Data()), nv0mbins,nv0mmin,nv0mmax,nptbins,ptmin,ptmax); 
    fhCaloV0MCorrECells->SetXTitle("V0 signal");
    fhCaloV0MCorrECells->SetYTitle(Form("#Sigma E of Cells in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloV0MCorrECells);    
    
    //Calorimeter VS Track multiplicity
    fhCaloTrackMCorrNClusters  = new TH2F ("hCaloTrackMNClusters",Form("# clusters in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax,nbins,nmin,nmax); 
    fhCaloTrackMCorrNClusters->SetXTitle("# tracks");
    fhCaloTrackMCorrNClusters->SetYTitle(Form("number of clusters in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrNClusters);
    
    fhCaloTrackMCorrEClusters  = new TH2F ("hCaloTrackMEClusters",Form("summed energy of clusters in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax,nptbins,ptmin,ptmax); 
    fhCaloTrackMCorrEClusters->SetXTitle("# tracks");
    fhCaloTrackMCorrEClusters->SetYTitle(Form("#Sigma E of clusters in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrEClusters);
    
    fhCaloTrackMCorrNCells  = new TH2F ("hCaloTrackMNCells",Form("# Cells in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax, nbins,nmin,nmax); 
    fhCaloTrackMCorrNCells->SetXTitle("# tracks");
    fhCaloTrackMCorrNCells->SetYTitle(Form("number of Cells in %s",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrNCells);
    
    fhCaloTrackMCorrECells  = new TH2F ("hCaloTrackMECells",Form("summed energy of Cells in %s vs # tracks",fCalorimeter.Data()), ntrmbins,ntrmmin,ntrmmax,nptbins,ptmin,ptmax); 
    fhCaloTrackMCorrECells->SetXTitle("# tracks");
    fhCaloTrackMCorrECells->SetYTitle(Form("#Sigma E of Cells in %s (GeV)",fCalorimeter.Data()));
    outputContainer->Add(fhCaloTrackMCorrECells);    
    
    
  }//correlate calorimeters
  
  //Module histograms
  fhEMod                 = new TH1F*[fNModules];
  fhNClustersMod         = new TH1F*[fNModules];
  fhNCellsPerClusterMod  = new TH2F*[fNModules];
  fhNCellsPerClusterModNoCut  = new TH2F*[fNModules];
  fhNCellsMod            = new TH1F*[fNModules];
  fhGridCellsMod         = new TH2F*[fNModules];
  fhGridCellsEMod        = new TH2F*[fNModules];
  fhGridCellsTimeMod     = new TH2F*[fNModules];
  fhAmplitudeMod         = new TH1F*[fNModules];
  if(fCalorimeter=="EMCAL")
    fhAmplitudeModFraction = new TH1F*[fNModules*3];
  
  fhTimeAmpPerRCU        = new TH2F*[fNModules*fNRCU];
  //fhT0TimeAmpPerRCU      = new TH2F*[fNModules*fNRCU];
  //fhTimeCorrRCU          = new TH2F*[fNModules*fNRCU*fNModules*fNRCU];
  
  fhIMMod                = new TH2F*[fNModules];
  fhIMCellCutMod         = new TH2F*[fNModules];
  
  for(Int_t imod = 0; imod < fNModules; imod++){
    
    fhEMod[imod]  = new TH1F (Form("hE_Mod%d",imod),Form("Cluster reconstructed Energy in Module %d ",imod), nptbins,ptmin,ptmax); 
    fhEMod[imod]->SetXTitle("E (GeV)");
    outputContainer->Add(fhEMod[imod]);
    
    fhNClustersMod[imod]  = new TH1F (Form("hNClusters_Mod%d",imod),Form("# clusters in Module %d",imod), nbins,nmin,nmax); 
    fhNClustersMod[imod]->SetXTitle("number of clusters");
    outputContainer->Add(fhNClustersMod[imod]);
    
    fhNCellsPerClusterMod[imod]  = new TH2F (Form("hNCellsPerCluster_Mod%d",imod),
                                             Form("# cells per cluster vs cluster energy in Module %d",imod), 
                                             nptbins,ptmin,ptmax, nbins,nmin,nmax); 
    fhNCellsPerClusterMod[imod]->SetXTitle("E (GeV)");
    fhNCellsPerClusterMod[imod]->SetYTitle("n cells");
    outputContainer->Add(fhNCellsPerClusterMod[imod]);

    fhNCellsPerClusterModNoCut[imod]  = new TH2F (Form("hNCellsPerClusterNoCut_Mod%d",imod),
                                             Form("# cells per cluster vs cluster energy in Module %d, no cut",imod), 
                                             nptbins,ptmin,ptmax, nbins,nmin,nmax); 
    fhNCellsPerClusterModNoCut[imod]->SetXTitle("E (GeV)");
    fhNCellsPerClusterModNoCut[imod]->SetYTitle("n cells");
    outputContainer->Add(fhNCellsPerClusterModNoCut[imod]);
    
    
    fhNCellsMod[imod]  = new TH1F (Form("hNCells_Mod%d",imod),Form("# cells in Module %d",imod), colmax*rowmax,0,colmax*rowmax); 
    fhNCellsMod[imod]->SetXTitle("n cells");
    outputContainer->Add(fhNCellsMod[imod]);
    fhGridCellsMod[imod]  = new TH2F (Form("hGridCells_Mod%d",imod),Form("Entries in grid of cells in Module %d",imod), 
                                      colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhGridCellsMod[imod]->SetYTitle("row (phi direction)");
    fhGridCellsMod[imod]->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsMod[imod]);
    
    fhGridCellsEMod[imod]  = new TH2F (Form("hGridCellsE_Mod%d",imod),Form("Accumulated energy in grid of cells in Module %d",imod), 
                                       colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhGridCellsEMod[imod]->SetYTitle("row (phi direction)");
    fhGridCellsEMod[imod]->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsEMod[imod]);
    
    fhGridCellsTimeMod[imod]  = new TH2F (Form("hGridCellsTime_Mod%d",imod),Form("Accumulated time in grid of cells in Module %d, with E > 0.5 GeV",imod), 
                                          colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
    fhGridCellsTimeMod[imod]->SetYTitle("row (phi direction)");
    fhGridCellsTimeMod[imod]->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsTimeMod[imod]);
    
    fhAmplitudeMod[imod]  = new TH1F (Form("hAmplitude_Mod%d",imod),Form("Cell Energy in Module %d",imod), nptbins*2,ptmin,ptmax); 
    fhAmplitudeMod[imod]->SetXTitle("Cell Energy (GeV)");
    outputContainer->Add(fhAmplitudeMod[imod]);
    
    if(fCalorimeter == "EMCAL"){
      for(Int_t ifrac = 0; ifrac < 3; ifrac++){
        fhAmplitudeModFraction[imod*3+ifrac]  = new TH1F (Form("hAmplitude_Mod%d_Frac%d",imod,ifrac),Form("Cell reconstructed Energy in Module %d, Fraction %d ",imod,ifrac), nptbins,ptmin,ptmax); 
        fhAmplitudeModFraction[imod*3+ifrac]->SetXTitle("E (GeV)");
        outputContainer->Add(fhAmplitudeModFraction[imod*3+ifrac]);
      }
      
    }
    if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
      
      for(Int_t ircu = 0; ircu < fNRCU; ircu++){
        fhTimeAmpPerRCU[imod*fNRCU+ircu]  = new TH2F (Form("hTimeAmp_Mod%d_RCU%d",imod,ircu),
                                                      Form("Cell Energy vs Cell Time in Module %d, RCU %d ",imod,ircu), 
                                                      nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetXTitle("E (GeV)");
        fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetYTitle("time (ns)");
        outputContainer->Add(fhTimeAmpPerRCU[imod*fNRCU+ircu]);
        
        //				fhT0TimeAmpPerRCU[imod*fNRCU+ircu]  = new TH2F (Form("hT0TimeAmp_Mod%d_RCU%d",imod,ircu),
        //															  Form("Cell Energy vs T0-Cell Time in Module %d, RCU %d ",imod,ircu), 
        //															  nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        //				fhT0TimeAmpPerRCU[imod*fNRCU+ircu]->SetXTitle("E (GeV)");
        //				fhT0TimeAmpPerRCU[imod*fNRCU+ircu]->SetYTitle("T_{0} - T_{EMCal} (ns)");
        //				outputContainer->Add(fhT0TimeAmpPerRCU[imod*fNRCU+ircu]);
        //			
        
        //				for(Int_t imod2 = 0; imod2 < fNModules; imod2++){
        //						for(Int_t ircu2 = 0; ircu2 < fNModules; ircu2++){
        //							Int_t index =  (imod2*fNRCU+ircu2)+(fNModules*fNRCU)*(ircu+imod)+fNRCU*fNModules*imod; 
        //							fhTimeCorrRCU[index]  = new TH2F (Form("hTimeCorrRCU_Mod%d_RCU%d_CompareTo_Mod%d_RCU%d",imod, ircu,imod2, ircu2),
        //																			Form("Cell Energy > 0.3, Correlate cell Time in Module %d, RCU %d to Module %d, RCU %d",imod,ircu,imod2, ircu2),
        //																			ntimebins,timemin,timemax,ntimebins,timemin,timemax); 
        //							fhTimeCorrRCU[index]->SetXTitle("Trigger Cell Time (ns)");
        //							fhTimeCorrRCU[index]->SetYTitle("Cell Time (ns)");
        //							outputContainer->Add(fhTimeCorrRCU[index]);
        //						}
        //				}
      }
    }
    if(fFillAllPi0Histo){
      fhIMMod[imod]  = new TH2F (Form("hIM_Mod%d",imod),
                                 Form("Cluster pairs Invariant mass vs reconstructed pair energy in Module %d",imod),
                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMMod[imod]->SetXTitle("p_{T, cluster pairs} (GeV) ");
      fhIMMod[imod]->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
      outputContainer->Add(fhIMMod[imod]);
      
      fhIMCellCutMod[imod]  = new TH2F (Form("hIMCellCut_Mod%d",imod),
                                        Form("Cluster (n cells > 1) pairs Invariant mass vs reconstructed pair energy in Module %d",imod),
                                        nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMCellCutMod[imod]->SetXTitle("p_{T, cluster pairs} (GeV) ");
      fhIMCellCutMod[imod]->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
      outputContainer->Add(fhIMCellCutMod[imod]);
    }
  }
  
  
  //Monte Carlo Histograms
  if(IsDataMC()){
    
    fhDeltaE  = new TH1F ("hDeltaE","MC - Reco E ", nptbins*2,-ptmax,ptmax); 
    fhDeltaE->SetXTitle("#Delta E (GeV)");
    outputContainer->Add(fhDeltaE);
    
    fhDeltaPt  = new TH1F ("hDeltaPt","MC - Reco p_{T} ", nptbins*2,-ptmax,ptmax); 
    fhDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
    outputContainer->Add(fhDeltaPt);
    
    fhDeltaPhi  = new TH1F ("hDeltaPhi","MC - Reco #phi ",nphibins*2,-phimax,phimax); 
    fhDeltaPhi->SetXTitle("#Delta #phi (rad)");
    outputContainer->Add(fhDeltaPhi);
    
    fhDeltaEta  = new TH1F ("hDeltaEta","MC- Reco #eta",netabins*2,-etamax,etamax); 
    fhDeltaEta->SetXTitle("#Delta #eta ");
    outputContainer->Add(fhDeltaEta);
    
    fhRatioE  = new TH1F ("hRatioE","Reco/MC E ", nratiobins,ratiomin,ratiomax); 
    fhRatioE->SetXTitle("E_{reco}/E_{gen}");
    outputContainer->Add(fhRatioE);
    
    fhRatioPt  = new TH1F ("hRatioPt","Reco/MC p_{T} ", nratiobins,ratiomin,ratiomax); 
    fhRatioPt->SetXTitle("p_{T, reco}/p_{T, gen}");
    outputContainer->Add(fhRatioPt);
    
    fhRatioPhi  = new TH1F ("hRatioPhi","Reco/MC #phi ",nratiobins,ratiomin,ratiomax); 
    fhRatioPhi->SetXTitle("#phi_{reco}/#phi_{gen}");
    outputContainer->Add(fhRatioPhi);
    
    fhRatioEta  = new TH1F ("hRatioEta","Reco/MC #eta",nratiobins,ratiomin,ratiomax); 
    fhRatioEta->SetXTitle("#eta_{reco}/#eta_{gen} ");
    outputContainer->Add(fhRatioEta);
    
    fh2E  = new TH2F ("h2E","E distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fh2E->SetXTitle("E_{rec} (GeV)");
    fh2E->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fh2E);	  
    
    fh2Pt  = new TH2F ("h2Pt","p_T distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
    fh2Pt->SetXTitle("p_{T,rec} (GeV/c)");
    fh2Pt->SetYTitle("p_{T,gen} (GeV/c)");
    outputContainer->Add(fh2Pt);
    
    fh2Phi  = new TH2F ("h2Phi","#phi distribution, reconstructed vs generated", nphibins,phimin,phimax, nphibins,phimin,phimax); 
    fh2Phi->SetXTitle("#phi_{rec} (rad)");
    fh2Phi->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fh2Phi);
    
    fh2Eta  = new TH2F ("h2Eta","#eta distribution, reconstructed vs generated", netabins,etamin,etamax,netabins,etamin,etamax); 
    fh2Eta->SetXTitle("#eta_{rec} ");
    fh2Eta->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fh2Eta);
    
    //Fill histos depending on origin of cluster
    fhGamE  = new TH2F ("hGamE","E reconstructed vs E generated from #gamma", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhGamE->SetXTitle("E_{rec} (GeV)");
    fhGamE->SetXTitle("E_{gen} (GeV)");
    outputContainer->Add(fhGamE);
    
    fhGamPt  = new TH2F ("hGamPt","p_{T} reconstructed vs E generated from #gamma", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhGamPt->SetXTitle("p_{T rec} (GeV/c)");
    fhGamPt->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhGamPt);
    
    fhGamPhi  = new TH2F ("hGamPhi","#phi reconstructed vs E generated from #gamma",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhGamPhi->SetXTitle("#phi_{rec} (rad)");
    fhGamPhi->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhGamPhi);
    
    fhGamEta  = new TH2F ("hGamEta","#eta reconstructed vs E generated from #gamma",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhGamEta->SetXTitle("#eta_{rec} ");
    fhGamEta->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhGamEta);
    
    fhGamDeltaE  = new TH1F ("hGamDeltaE","#gamma MC - Reco E ", nptbins*2,-ptmax,ptmax); 
    fhGamDeltaE->SetXTitle("#Delta E (GeV)");
    outputContainer->Add(fhGamDeltaE);
    
    fhGamDeltaPt  = new TH1F ("hGamDeltaPt","#gamma MC - Reco p_{T} ", nptbins*2,-ptmax,ptmax); 
    fhGamDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
    outputContainer->Add(fhGamDeltaPt);
    
    fhGamDeltaPhi  = new TH1F ("hGamDeltaPhi","#gamma MC - Reco #phi ",nphibins*2,-phimax,phimax); 
    fhGamDeltaPhi->SetXTitle("#Delta #phi (rad)");
    outputContainer->Add(fhGamDeltaPhi);
    
    fhGamDeltaEta  = new TH1F ("hGamDeltaEta","#gamma MC- Reco #eta",netabins*2,-etamax,etamax); 
    fhGamDeltaEta->SetXTitle("#Delta #eta ");
    outputContainer->Add(fhGamDeltaEta);
    
    fhGamRatioE  = new TH1F ("hGamRatioE","#gamma Reco/MC E ", nratiobins,ratiomin,ratiomax); 
    fhGamRatioE->SetXTitle("E_{reco}/E_{gen}");
    outputContainer->Add(fhGamRatioE);
    
    fhGamRatioPt  = new TH1F ("hGamRatioPt","#gamma Reco/MC p_{T} ", nratiobins,ratiomin,ratiomax); 
    fhGamRatioPt->SetXTitle("p_{T, reco}/p_{T, gen}");
    outputContainer->Add(fhGamRatioPt);
    
    fhGamRatioPhi  = new TH1F ("hGamRatioPhi","#gamma Reco/MC #phi ",nratiobins,ratiomin,ratiomax); 
    fhGamRatioPhi->SetXTitle("#phi_{reco}/#phi_{gen}");
    outputContainer->Add(fhGamRatioPhi);
    
    fhGamRatioEta  = new TH1F ("hGamRatioEta","#gamma Reco/MC #eta",nratiobins,ratiomin,ratiomax); 
    fhGamRatioEta->SetXTitle("#eta_{reco}/#eta_{gen} ");
    outputContainer->Add(fhGamRatioEta);
    
    fhPi0E  = new TH2F ("hPi0E","E reconstructed vs E generated from #pi^{0}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhPi0E->SetXTitle("E_{rec} (GeV)");
    fhPi0E->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fhPi0E);
    
    fhPi0Pt  = new TH2F ("hPi0Pt","p_{T} reconstructed vs E generated from #pi^{0}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhPi0Pt->SetXTitle("p_{T rec} (GeV/c)");
    fhPi0Pt->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhPi0Pt);
    
    fhPi0Phi  = new TH2F ("hPi0Phi","#phi reconstructed vs E generated from #pi^{0}",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhPi0Phi->SetXTitle("#phi_{rec} (rad)");
    fhPi0Phi->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhPi0Phi);
    
    fhPi0Eta  = new TH2F ("hPi0Eta","#eta reconstructed vs E generated from #pi^{0}",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhPi0Eta->SetXTitle("#eta_{rec} ");
    fhPi0Eta->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhPi0Eta);
    
    fhEleE  = new TH2F ("hEleE","E reconstructed vs E generated from e^{#pm}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhEleE->SetXTitle("E_{rec} (GeV)");
    fhEleE->SetXTitle("E_{gen} (GeV)");		
    outputContainer->Add(fhEleE);		
    
    fhElePt  = new TH2F ("hElePt","p_{T} reconstructed vs E generated from e^{#pm}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhElePt->SetXTitle("p_{T rec} (GeV/c)");
    fhElePt->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhElePt);
    
    fhElePhi  = new TH2F ("hElePhi","#phi reconstructed vs E generated from e^{#pm}",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhElePhi->SetXTitle("#phi_{rec} (rad)");
    fhElePhi->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhElePhi);
    
    fhEleEta  = new TH2F ("hEleEta","#eta reconstructed vs E generated from e^{#pm}",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhEleEta->SetXTitle("#eta_{rec} ");
    fhEleEta->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhEleEta);
    
    fhNeHadE  = new TH2F ("hNeHadE","E reconstructed vs E generated from neutral hadron", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhNeHadE->SetXTitle("E_{rec} (GeV)");
    fhNeHadE->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fhNeHadE);
    
    fhNeHadPt  = new TH2F ("hNeHadPt","p_{T} reconstructed vs E generated from neutral hadron", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhNeHadPt->SetXTitle("p_{T rec} (GeV/c)");
    fhNeHadPt->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhNeHadPt);
    
    fhNeHadPhi  = new TH2F ("hNeHadPhi","#phi reconstructed vs E generated from neutral hadron",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhNeHadPhi->SetXTitle("#phi_{rec} (rad)");
    fhNeHadPhi->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhNeHadPhi);
    
    fhNeHadEta  = new TH2F ("hNeHadEta","#eta reconstructed vs E generated from neutral hadron",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhNeHadEta->SetXTitle("#eta_{rec} ");
    fhNeHadEta->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhNeHadEta);
    
    fhChHadE  = new TH2F ("hChHadE","E reconstructed vs E generated from charged hadron", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhChHadE->SetXTitle("E_{rec} (GeV)");
    fhChHadE->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fhChHadE);
    
    fhChHadPt  = new TH2F ("hChHadPt","p_{T} reconstructed vs E generated from charged hadron", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhChHadPt->SetXTitle("p_{T rec} (GeV/c)");
    fhChHadPt->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhChHadPt);
    
    fhChHadPhi  = new TH2F ("hChHadPhi","#phi reconstructed vs E generated from charged hadron",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhChHadPhi->SetXTitle("#phi_{rec} (rad)");
    fhChHadPhi->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhChHadPhi);
    
    fhChHadEta  = new TH2F ("hChHadEta","#eta reconstructed vs E generated from charged hadron",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhChHadEta->SetXTitle("#eta_{rec} ");
    fhChHadEta->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhChHadEta);
    
    //Charged clusters
    
    fhGamECharged  = new TH2F ("hGamECharged","E reconstructed vs E generated from #gamma, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhGamECharged->SetXTitle("E_{rec} (GeV)");
    fhGamECharged->SetXTitle("E_{gen} (GeV)");
    outputContainer->Add(fhGamECharged);
    
    fhGamPtCharged  = new TH2F ("hGamPtCharged","p_{T} reconstructed vs E generated from #gamma, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhGamPtCharged->SetXTitle("p_{T rec} (GeV/c)");
    fhGamPtCharged->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhGamPtCharged);
    
    fhGamPhiCharged  = new TH2F ("hGamPhiCharged","#phi reconstructed vs E generated from #gamma, track matched cluster",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhGamPhiCharged->SetXTitle("#phi_{rec} (rad)");
    fhGamPhiCharged->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhGamPhiCharged);
    
    fhGamEtaCharged  = new TH2F ("hGamEtaCharged","#eta reconstructed vs E generated from #gamma, track matched cluster",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhGamEtaCharged->SetXTitle("#eta_{rec} ");
    fhGamEtaCharged->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhGamEtaCharged);
    
    fhPi0ECharged  = new TH2F ("hPi0ECharged","E reconstructed vs E generated from #pi^{0}, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhPi0ECharged->SetXTitle("E_{rec} (GeV)");
    fhPi0ECharged->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fhPi0ECharged);
    
    fhPi0PtCharged  = new TH2F ("hPi0PtCharged","p_{T} reconstructed vs E generated from #pi^{0}, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhPi0PtCharged->SetXTitle("p_{T rec} (GeV/c)");
    fhPi0PtCharged->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhPi0PtCharged);
    
    fhPi0PhiCharged  = new TH2F ("hPi0PhiCharged","#phi reconstructed vs E generated from #pi^{0}, track matched cluster",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhPi0PhiCharged->SetXTitle("#phi_{rec} (rad)");
    fhPi0PhiCharged->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhPi0PhiCharged);
    
    fhPi0EtaCharged  = new TH2F ("hPi0EtaCharged","#eta reconstructed vs E generated from #pi^{0}, track matched cluster",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhPi0EtaCharged->SetXTitle("#eta_{rec} ");
    fhPi0EtaCharged->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhPi0EtaCharged);
    
    fhEleECharged  = new TH2F ("hEleECharged","E reconstructed vs E generated from e^{#pm}, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhEleECharged->SetXTitle("E_{rec} (GeV)");
    fhEleECharged->SetXTitle("E_{gen} (GeV)");		
    outputContainer->Add(fhEleECharged);		
    
    fhElePtCharged  = new TH2F ("hElePtCharged","p_{T} reconstructed vs E generated from e^{#pm}, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhElePtCharged->SetXTitle("p_{T rec} (GeV/c)");
    fhElePtCharged->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhElePtCharged);
    
    fhElePhiCharged  = new TH2F ("hElePhiCharged","#phi reconstructed vs E generated from e^{#pm}, track matched cluster",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhElePhiCharged->SetXTitle("#phi_{rec} (rad)");
    fhElePhiCharged->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhElePhiCharged);
    
    fhEleEtaCharged  = new TH2F ("hEleEtaCharged","#eta reconstructed vs E generated from e^{#pm}, track matched cluster",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhEleEtaCharged->SetXTitle("#eta_{rec} ");
    fhEleEtaCharged->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhEleEtaCharged);
    
    fhNeHadECharged  = new TH2F ("hNeHadECharged","E reconstructed vs E generated from neutral hadron, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhNeHadECharged->SetXTitle("E_{rec} (GeV)");
    fhNeHadECharged->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fhNeHadECharged);
    
    fhNeHadPtCharged  = new TH2F ("hNeHadPtCharged","p_{T} reconstructed vs E generated from neutral hadron, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhNeHadPtCharged->SetXTitle("p_{T rec} (GeV/c)");
    fhNeHadPtCharged->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhNeHadPtCharged);
    
    fhNeHadPhiCharged  = new TH2F ("hNeHadPhiCharged","#phi reconstructed vs E generated from neutral hadron, track matched cluster",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhNeHadPhiCharged->SetXTitle("#phi_{rec} (rad)");
    fhNeHadPhiCharged->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhNeHadPhiCharged);
    
    fhNeHadEtaCharged  = new TH2F ("hNeHadEtaCharged","#eta reconstructed vs E generated from neutral hadron, track matched cluster",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhNeHadEtaCharged->SetXTitle("#eta_{rec} ");
    fhNeHadEtaCharged->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhNeHadEtaCharged);
    
    fhChHadECharged  = new TH2F ("hChHadECharged","E reconstructed vs E generated from charged hadron, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhChHadECharged->SetXTitle("E_{rec} (GeV)");
    fhChHadECharged->SetYTitle("E_{gen} (GeV)");
    outputContainer->Add(fhChHadECharged);
    
    fhChHadPtCharged  = new TH2F ("hChHadPtCharged","p_{T} reconstructed vs E generated from charged hadron, track matched cluster", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
    fhChHadPtCharged->SetXTitle("p_{T rec} (GeV/c)");
    fhChHadPtCharged->SetYTitle("p_{T gen} (GeV/c)");
    outputContainer->Add(fhChHadPtCharged);
    
    fhChHadPhiCharged  = new TH2F ("hChHadPhiCharged","#phi reconstructed vs E generated from charged hadron, track matched cluster",nphibins,phimin,phimax,nphibins,phimin,phimax); 
    fhChHadPhiCharged->SetXTitle("#phi (rad)");
    fhChHadPhiCharged->SetXTitle("#phi_{rec} (rad)");
    fhChHadPhiCharged->SetYTitle("#phi_{gen} (rad)");
    outputContainer->Add(fhChHadPhiCharged);
    
    fhChHadEtaCharged  = new TH2F ("hChHadEtaCharged","#eta reconstructed vs E generated from charged hadron, track matched cluster",netabins,etamin,etamax,netabins,etamin,etamax); 
    fhChHadEtaCharged->SetXTitle("#eta_{rec} ");
    fhChHadEtaCharged->SetYTitle("#eta_{gen} ");
    outputContainer->Add(fhChHadEtaCharged);
    
    //Vertex of generated particles 
    
    fhEMVxyz  = new TH2F ("hEMVxyz","Production vertex of reconstructed ElectroMagnetic particles",nvdistbins,vdistmin,vdistmax,nvdistbins,vdistmin,vdistmax);//,100,0,500); 
    fhEMVxyz->SetXTitle("v_{x}");
    fhEMVxyz->SetYTitle("v_{y}");
    //fhEMVxyz->SetZTitle("v_{z}");
    outputContainer->Add(fhEMVxyz);
    
    fhHaVxyz  = new TH2F ("hHaVxyz","Production vertex of reconstructed hadrons",nvdistbins,vdistmin,vdistmax,nvdistbins,vdistmin,vdistmax);//,100,0,500); 
    fhHaVxyz->SetXTitle("v_{x}");
    fhHaVxyz->SetYTitle("v_{y}");
    //fhHaVxyz->SetZTitle("v_{z}");
    outputContainer->Add(fhHaVxyz);
    
    fhEMR  = new TH2F ("hEMR","Distance to production vertex of reconstructed ElectroMagnetic particles vs E rec",nptbins,ptmin,ptmax,nvdistbins,vdistmin,vdistmax); 
    fhEMR->SetXTitle("E (GeV)");
    fhEMR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
    outputContainer->Add(fhEMR);
    
    fhHaR  = new TH2F ("hHaR","Distance to production vertex of reconstructed Hadrons vs E rec",nptbins,ptmin,ptmax,nvdistbins,vdistmin,vdistmax); 
    fhHaR->SetXTitle("E (GeV)");
    fhHaR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
    outputContainer->Add(fhHaR);
    
    
    
    //Pure MC
    fhGenGamPt  = new TH1F("hGenGamPt" ,"p_{T} of generated #gamma",nptbins,ptmin,ptmax);
    fhGenGamEta = new TH1F("hGenGamEta","Y of generated #gamma",netabins,etamin,etamax);
    fhGenGamPhi = new TH1F("hGenGamPhi","#phi of generated #gamma",nphibins,phimin,phimax);
    
    fhGenPi0Pt  = new TH1F("hGenPi0Pt" ,"p_{T} of generated #pi^{0}",nptbins,ptmin,ptmax);
    fhGenPi0Eta = new TH1F("hGenPi0Eta","Y of generated #pi^{0}",netabins,etamin,etamax);
    fhGenPi0Phi = new TH1F("hGenPi0Phi","#phi of generated #pi^{0}",nphibins,phimin,phimax);
    
    fhGenEtaPt  = new TH1F("hGenEtaPt" ,"p_{T} of generated #eta",nptbins,ptmin,ptmax);
    fhGenEtaEta = new TH1F("hGenEtaEta","Y of generated #eta",netabins,etamin,etamax);
    fhGenEtaPhi = new TH1F("hGenEtaPhi","#phi of generated #eta",nphibins,phimin,phimax);
    
    fhGenOmegaPt  = new TH1F("hGenOmegaPt" ,"p_{T} of generated #omega",nptbins,ptmin,ptmax);
    fhGenOmegaEta = new TH1F("hGenOmegaEta","Y of generated #omega",netabins,etamin,etamax);
    fhGenOmegaPhi = new TH1F("hGenOmegaPhi","#phi of generated #omega",nphibins,phimin,phimax);		
    
    fhGenElePt  = new TH1F("hGenElePt" ,"p_{T} of generated e^{#pm}",nptbins,ptmin,ptmax);
    fhGenEleEta = new TH1F("hGenEleEta","Y of generated  e^{#pm}",netabins,etamin,etamax);
    fhGenElePhi = new TH1F("hGenElePhi","#phi of generated  e^{#pm}",nphibins,phimin,phimax);		
    
    fhGenGamPt->SetXTitle("p_{T} (GeV/c)");
    fhGenGamEta->SetXTitle("#eta");
    fhGenGamPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenGamPt);
    outputContainer->Add(fhGenGamEta);
    outputContainer->Add(fhGenGamPhi);
    
    fhGenPi0Pt->SetXTitle("p_{T} (GeV/c)");
    fhGenPi0Eta->SetXTitle("#eta");
    fhGenPi0Phi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenPi0Pt);
    outputContainer->Add(fhGenPi0Eta);
    outputContainer->Add(fhGenPi0Phi);
    
    fhGenEtaPt->SetXTitle("p_{T} (GeV/c)");
    fhGenEtaEta->SetXTitle("#eta");
    fhGenEtaPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenEtaPt);
    outputContainer->Add(fhGenEtaEta);
    outputContainer->Add(fhGenEtaPhi);
    
    fhGenOmegaPt->SetXTitle("p_{T} (GeV/c)");
    fhGenOmegaEta->SetXTitle("#eta");
    fhGenOmegaPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenOmegaPt);
    outputContainer->Add(fhGenOmegaEta);
    outputContainer->Add(fhGenOmegaPhi);
    
    fhGenElePt->SetXTitle("p_{T} (GeV/c)");
    fhGenEleEta->SetXTitle("#eta");
    fhGenElePhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenElePt);
    outputContainer->Add(fhGenEleEta);
    outputContainer->Add(fhGenElePhi);
    
    fhGenGamAccE   = new TH1F("hGenGamAccE" ,"E of generated #gamma in calorimeter acceptance",nptbins,ptmin,ptmax);
    fhGenGamAccPt  = new TH1F("hGenGamAccPt" ,"p_{T} of generated #gamma in calorimeter acceptance",nptbins,ptmin,ptmax);
    fhGenGamAccEta = new TH1F("hGenGamAccEta","Y of generated #gamma in calorimeter acceptance",netabins,etamin,etamax);
    fhGenGamAccPhi = new TH1F("hGenGamAccPhi","#phi of generated #gamma  in calorimeter acceptance",nphibins,phimin,phimax);
    
    fhGenPi0AccE   = new TH1F("hGenPi0AccE" ,"E of generated #pi^{0} in calorimeter acceptance",nptbins,ptmin,ptmax);
    fhGenPi0AccPt  = new TH1F("hGenPi0AccPt" ,"p_{T} of generated #pi^{0} in calorimeter acceptance",nptbins,ptmin,ptmax);
    fhGenPi0AccEta = new TH1F("hGenPi0AccEta","Y of generated #pi^{0} in calorimeter acceptance",netabins,etamin,etamax);
    fhGenPi0AccPhi = new TH1F("hGenPi0AccPhi","#phi of generated #pi^{0} in calorimeter acceptance",nphibins,phimin,phimax);
    
    fhGenGamAccE  ->SetXTitle("E (GeV)");
    fhGenGamAccPt ->SetXTitle("p_{T} (GeV/c)");
    fhGenGamAccEta->SetXTitle("#eta");
    fhGenGamAccPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenGamAccE);		
    outputContainer->Add(fhGenGamAccPt);
    outputContainer->Add(fhGenGamAccEta);
    outputContainer->Add(fhGenGamAccPhi);
    
    fhGenPi0AccE  ->SetXTitle("E (GeV)");		
    fhGenPi0AccPt ->SetXTitle("p_{T} (GeV/c)");
    fhGenPi0AccEta->SetXTitle("#eta");
    fhGenPi0AccPhi->SetXTitle("#phi (rad)");
    outputContainer->Add(fhGenPi0AccE);		
    outputContainer->Add(fhGenPi0AccPt);
    outputContainer->Add(fhGenPi0AccEta);
    outputContainer->Add(fhGenPi0AccPhi);
    
    //Track Matching 
    
    fhMCEle1pOverE = new TH2F("hMCEle1pOverE","TRACK matches p/E, MC electrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCEle1pOverE->SetYTitle("p/E");
    fhMCEle1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCEle1pOverE);
    
    fhMCEle1dR = new TH1F("hMCEle1dR","TRACK matches dR, MC electrons",ndRbins,dRmin,dRmax);
    fhMCEle1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCEle1dR) ;
    
    fhMCEle2MatchdEdx = new TH2F("hMCEle2MatchdEdx","dE/dx vs. p for all matches, MC electrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCEle2MatchdEdx->SetXTitle("p (GeV/c)");
    fhMCEle2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhMCEle2MatchdEdx);
    
    fhMCChHad1pOverE = new TH2F("hMCChHad1pOverE","TRACK matches p/E, MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCChHad1pOverE->SetYTitle("p/E");
    fhMCChHad1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCChHad1pOverE);
    
    fhMCChHad1dR = new TH1F("hMCChHad1dR","TRACK matches dR, MC charged hadrons",ndRbins,dRmin,dRmax);
    fhMCChHad1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCChHad1dR) ;
    
    fhMCChHad2MatchdEdx = new TH2F("hMCChHad2MatchdEdx","dE/dx vs. p for all matches, MC charged hadrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCChHad2MatchdEdx->SetXTitle("p (GeV/c)");
    fhMCChHad2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhMCChHad2MatchdEdx);
    
    fhMCNeutral1pOverE = new TH2F("hMCNeutral1pOverE","TRACK matches p/E, MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCNeutral1pOverE->SetYTitle("p/E");
    fhMCNeutral1pOverE->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCNeutral1pOverE);
    
    fhMCNeutral1dR = new TH1F("hMCNeutral1dR","TRACK matches dR, MC neutrals",ndRbins,dRmin,dRmax);
    fhMCNeutral1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCNeutral1dR) ;
    
    fhMCNeutral2MatchdEdx = new TH2F("hMCNeutral2MatchdEdx","dE/dx vs. p for all matches, MC neutrals",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCNeutral2MatchdEdx->SetXTitle("p (GeV/c)");
    fhMCNeutral2MatchdEdx->SetYTitle("<dE/dx>");
    outputContainer->Add(fhMCNeutral2MatchdEdx);
    
    fhMCEle1pOverER02 = new TH2F("hMCEle1pOverER02","TRACK matches p/E, MC electrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCEle1pOverER02->SetYTitle("p/E");
    fhMCEle1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCEle1pOverER02);
    
    fhMCChHad1pOverER02 = new TH2F("hMCChHad1pOverER02","TRACK matches p/E, MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCChHad1pOverER02->SetYTitle("p/E");
    fhMCChHad1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCChHad1pOverER02);
    
    fhMCNeutral1pOverER02 = new TH2F("hMCNeutral1pOverER02","TRACK matches p/E, MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,pOverEmin,pOverEmax);
    fhMCNeutral1pOverER02->SetYTitle("p/E");
    fhMCNeutral1pOverER02->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhMCNeutral1pOverER02);
  }
  
  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //    printf("i=%d, name= %s\n",i,outputContainer->At(i)->GetName());
  
  return outputContainer;
}

//_______________________________________________________________________________________________________________________________________
Int_t AliAnaCalorimeterQA::GetNewRebinForRePlotting(TH1D* histo, const Float_t newXmin, const Float_t newXmax,const Int_t newXnbins) const
{
  //Calculate the rebinning for the new requested bin size, only used when replotting executing the Terminte
 
  Float_t oldbinsize =  histo->GetBinWidth(0);
  Float_t newbinsize = TMath::Abs(newXmax-newXmin) / newXnbins;

  //printf("bin size, old %f, new %f\n",oldbinsize,newbinsize);
  if(newbinsize > oldbinsize) return (Int_t) (newbinsize/oldbinsize);
  else  return 1;

}

//__________________________________________________
void AliAnaCalorimeterQA::Init()
{ 
  //Check if the data or settings are ok
  
  if(fCalorimeter != "PHOS" && fCalorimeter !="EMCAL")
    AliFatal(Form("Wrong calorimeter name <%s>", fCalorimeter.Data()));
  
  if(GetReader()->GetDataType()== AliCaloTrackReader::kMC)
    AliFatal("Analysis of reconstructed data, MC reader not aplicable");
  
}


//__________________________________________________
void AliAnaCalorimeterQA::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaCaloQA_");
  
  fCalorimeter     = "EMCAL"; //or PHOS
  fStyleMacro      = "" ;
  fNModules        = 12; // set maximum to maximum number of EMCAL modules
  fNRCU            = 2;  // set maximum number of RCU in EMCAL per SM
  fTimeCutMin      = -1;
  fTimeCutMax      = 9999999;
  fEMCALCellAmpMin = 0.0;
  fPHOSCellAmpMin  = 0.0;
  	
}

//__________________________________________________________________
void AliAnaCalorimeterQA::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  
  printf("Select Calorimeter %s \n",fCalorimeter.Data());
  printf("Plots style macro  %s \n",fStyleMacro.Data()); 
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("EMCAL Min Amplitude   : %2.1f GeV/c\n", fEMCALCellAmpMin) ;
  printf("PHOS Min Amplitude    : %2.1f GeV/c\n", fPHOSCellAmpMin) ;

} 

//__________________________________________________________________
void  AliAnaCalorimeterQA::MakeAnalysisFillHistograms() 
{
  //Fill Calorimeter QA histograms
  TLorentzVector mom  ;
  TLorentzVector mom2 ;
  TObjArray * caloClusters = NULL;
  Int_t nLabel = 0;
  Int_t *labels=0x0;
  Int_t nCaloClusters = 0;
  Int_t nCaloClustersAccepted = 0;
  Int_t nCaloCellsPerCluster = 0;
  Int_t nTracksMatched = 0;
  Int_t trackIndex = 0;
  Int_t nModule = -1;
  
  //Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  if (TMath::Abs(v[2]) > GetZvertexCut()) return ;  
  
  //Play with the MC stack if available	
  //Get the MC arrays and do some checks
  if(IsDataMC()){
    if(GetReader()->ReadStack()){
      
      if(!GetMCStack()) 
        AliFatal("Stack not available, is the MC handler called?\n");
      
      //Fill some pure MC histograms, only primaries.
      for(Int_t i=0 ; i<GetMCStack()->GetNprimary(); i++){//Only primary particles, for all MC transport put GetNtrack()
        TParticle *primary = GetMCStack()->Particle(i) ;
        //printf("i %d, %s: status = %d, primary? %d\n",i, primary->GetName(), primary->GetStatusCode(), primary->IsPrimary());
        if (primary->GetStatusCode() > 11) continue; //Working for PYTHIA and simple generators, check for HERWIG 
        primary->Momentum(mom);
        MCHistograms(mom,TMath::Abs(primary->GetPdgCode()));
      } //primary loop
    }
    else if(GetReader()->ReadAODMCParticles()){
      
      if(!GetReader()->GetAODMCParticles(0)) 	
        AliFatal("AODMCParticles not available!");
      
      //Fill some pure MC histograms, only primaries.
      for(Int_t i=0 ; i < (GetReader()->GetAODMCParticles(0))->GetEntriesFast(); i++){
        AliAODMCParticle *aodprimary = (AliAODMCParticle*) (GetReader()->GetAODMCParticles(0))->At(i) ;
        //printf("i %d, %s: primary? %d physical primary? %d, flag %d\n",
        //	   i,(TDatabasePDG::Instance()->GetParticle(aodprimary->GetPdgCode()))->GetName(), 
        //	   aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary(), aodprimary->GetFlag());
        if (!aodprimary->IsPrimary()) continue; //accept all which is not MC transport generated. Don't know how to avoid partons
        //aodprimary->Momentum(mom);
        mom.SetPxPyPzE(aodprimary->Px(),aodprimary->Py(),aodprimary->Pz(),aodprimary->E());
        MCHistograms(mom,TMath::Abs(aodprimary->GetPdgCode()));
      } //primary loop
      
    }
  }// is data and MC	
  
  
  //Get List with CaloClusters  
  if      (fCalorimeter == "PHOS")  caloClusters = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL") caloClusters = GetEMCALClusters();
  else 
    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Wrong calorimeter name <%s>, END\n", fCalorimeter.Data()));
  
  //  if     (fCalorimeter == "EMCAL") GetReader()->GetInputEvent()->GetEMCALClusters(caloClusters);//GetEMCALClusters();
  //  else if(fCalorimeter == "PHOS")  GetReader()->GetInputEvent()->GetPHOSClusters (caloClusters);//GetPHOSClusters();
  //  else 
  //    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Wrong calorimeter name <%s>, END\n", fCalorimeter.Data()));
  
  if(!caloClusters) {
    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - No CaloClusters available\n"));
  }
  else{
    //----------------------------------------------------------
    //Correlate Calorimeters and V0 and track Multiplicity
    //----------------------------------------------------------
    if(fCorrelate)	Correlate();
    
    //----------------------------------------------------------
    // CALOCLUSTERS
    //----------------------------------------------------------
    
    nCaloClusters = caloClusters->GetEntriesFast() ; 
    Int_t *nClustersInModule = new Int_t[fNModules];
    for(Int_t imod = 0; imod < fNModules; imod++ ) nClustersInModule[imod] = 0;
    
    if(GetDebug() > 0)
      printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - In %s there are %d clusters \n", fCalorimeter.Data(), nCaloClusters);
    
    AliVTrack * track = 0x0;
    Float_t pos[3] ;
    Double_t tof = 0;
    //Loop over CaloClusters
    //if(nCaloClusters > 0)printf("QA  : Vertex Cut passed %f, cut %f, entries %d, %s\n",v[2], 40., nCaloClusters, fCalorimeter.Data());
    for(Int_t iclus = 0; iclus < nCaloClusters; iclus++){
      
      if(GetDebug() > 0) printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - cluster: %d/%d, data %d \n",
                                iclus+1,nCaloClusters,GetReader()->GetDataType());
      
      AliVCluster* clus =  (AliVCluster*)caloClusters->At(iclus);
      AliVCaloCells * cell = 0x0; 
      if(fCalorimeter == "PHOS") cell =  GetPHOSCells();
      else			                 cell =  GetEMCALCells();
      
      //Get cluster kinematics
      clus->GetPosition(pos);
      clus->GetMomentum(mom,v);
      tof = clus->GetTOF()*1e9;
      if(tof < fTimeCutMin || tof > fTimeCutMax) continue;
      
      //Check only certain regions
      Bool_t in = kTRUE;
      if(IsFiducialCutOn()) in =  GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(!in) continue;
      
      //MC labels
      nLabel = clus->GetNLabels();
      labels = clus->GetLabels();
      
      //Cells per cluster
      nCaloCellsPerCluster = clus->GetNCells();
      //if(mom.E() > 10 && nCaloCellsPerCluster == 1 ) printf("%s:************** E = %f ********** ncells = %d\n",fCalorimeter.Data(), mom.E(),nCaloCellsPerCluster);
      
      //matched cluster with tracks
      nTracksMatched = clus->GetNTracksMatched();
      if(GetReader()->GetDataType() == AliCaloTrackReader::kESD){
        trackIndex     = clus->GetTrackMatchedIndex();
        if(trackIndex >= 0){
          track = (AliVTrack*)GetReader()->GetInputEvent()->GetTrack(trackIndex);
        }
        else{
          if(nTracksMatched == 1) nTracksMatched = 0;
          track = 0;
        }
      }//kESD
      else{//AODs
        if(nTracksMatched > 0) track = (AliVTrack*)clus->GetTrackMatched(0);
      }
      
      //======================
      //Cells in cluster
      //======================
      
      //Get list of contributors
      UShort_t * indexList = clus->GetCellsAbsId() ;
      // check time of cells respect to max energy cell
      //Get maximum energy cell
      Int_t absId   = -1 ;
      //printf("nCaloCellsPerCluster %d\n",nCaloCellsPerCluster);
      if(fFillAllPosHisto){
        //Loop on cluster cells
        for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) {
          //	printf("Index %d\n",ipos);
          absId  = indexList[ipos]; 
          
          //Get position of cell compare to cluster
          
          if(fCalorimeter=="EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
            
            Double_t cellpos[] = {0, 0, 0};
            GetEMCALGeometry()->GetGlobal(absId, cellpos);
            
            fhDeltaCellClusterXNCells->Fill(pos[0]-cellpos[0],nCaloCellsPerCluster) ; 
            fhDeltaCellClusterYNCells->Fill(pos[1]-cellpos[1],nCaloCellsPerCluster) ; 
            fhDeltaCellClusterZNCells->Fill(pos[2]-cellpos[2],nCaloCellsPerCluster) ;
            
            fhDeltaCellClusterXE->Fill(pos[0]-cellpos[0],mom.E())  ; 
            fhDeltaCellClusterYE->Fill(pos[1]-cellpos[1],mom.E())  ; 
            fhDeltaCellClusterZE->Fill(pos[2]-cellpos[2],mom.E())  ; 
            
            Float_t r     = TMath::Sqrt(pos[0]*pos[0]        +pos[1]*pos[1]);//     +pos[2]*pos[2]);
            Float_t rcell = TMath::Sqrt(cellpos[0]*cellpos[0]+cellpos[1]*cellpos[1]);//+cellpos[2]*cellpos[2]);
            fhDeltaCellClusterRNCells->Fill(r-rcell, nCaloCellsPerCluster) ; 
            fhDeltaCellClusterRE     ->Fill(r-rcell, mom.E())  ; 
            
            //					Float_t celleta = 0, cellphi = 0;
            //					GetEMCALGeometry()->EtaPhiFromIndex(absId, celleta, cellphi); 
            //					Int_t imod = -1, iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1;
            //					GetEMCALGeometry()->GetCellIndex(absId,imod,iTower,iIphi,iIeta); 
            //					GetEMCALGeometry()->GetCellPhiEtaIndexInSModule(imod,iTower,
            //																				 iIphi, iIeta,iphi,ieta);
            //					printf("AbsId %d, SM %d, Index eta %d, phi %d\n", absId, imod, ieta, iphi);
            //					printf("Cluster E %f, eta %f, phi %f; Cell: Amp %f, eta %f, phi%f\n", mom.E(),mom.Eta(), mom.Phi()*TMath::RadToDeg(), cell->GetCellAmplitude(absId),celleta, cellphi*TMath::RadToDeg());
            //					printf("x cluster %f, x cell %f, cluster-cell %f\n",pos[0], cellpos[0],pos[0]-cellpos[0]);
            //					printf("y cluster %f, y cell %f, cluster-cell %f\n",pos[1], cellpos[1],pos[1]-cellpos[1]);
            //					printf("z cluster %f, z cell %f, cluster-cell %f\n",pos[2], cellpos[2],pos[2]-cellpos[2]);
            //					printf("r cluster %f, r cell %f, cluster-cell %f\n",r,      rcell,     r-rcell);
            //					
            
          }//EMCAL and its matrices are available
          else if(fCalorimeter=="PHOS" && GetCaloUtils()->IsPHOSGeoMatrixSet()){
            TVector3 xyz;
            Int_t relId[4], module;
            Float_t xCell, zCell;
            
            GetPHOSGeometry()->AbsToRelNumbering(absId,relId);
            module = relId[0];
            GetPHOSGeometry()->RelPosInModule(relId,xCell,zCell);
            GetPHOSGeometry()->Local2Global(module,xCell,zCell,xyz);
            
            fhDeltaCellClusterXNCells->Fill(pos[0]-xyz.X(),nCaloCellsPerCluster) ; 
            fhDeltaCellClusterYNCells->Fill(pos[1]-xyz.Y(),nCaloCellsPerCluster) ; 
            fhDeltaCellClusterZNCells->Fill(pos[2]-xyz.Z(),nCaloCellsPerCluster) ;
            
            fhDeltaCellClusterXE->Fill(pos[0]-xyz.X(),mom.E())  ; 
            fhDeltaCellClusterYE->Fill(pos[1]-xyz.Y(),mom.E())  ; 
            fhDeltaCellClusterZE->Fill(pos[2]-xyz.Z(),mom.E())  ; 
            
            Float_t r     = TMath::Sqrt(pos[0]*pos[0]  +pos[1]*pos[1]);//     +pos[2]*pos[2]);
            Float_t rcell = TMath::Sqrt(xyz.X()*xyz.X()+xyz.Y()*xyz.Y());//+xyz.Z()*xyz.Z());
            fhDeltaCellClusterRNCells->Fill(r-rcell, nCaloCellsPerCluster) ; 
            fhDeltaCellClusterRE     ->Fill(r-rcell, mom.E())  ; 
            
            //		    	  printf("x cluster %f, x cell %f, cluster-cell %f\n",pos[0], cellpos[0],pos[0]-cellpos[0]);
            //		        printf("y cluster %f, y cell %f, cluster-cell %f\n",pos[1], cellpos[1],pos[1]-cellpos[1]);
            //	       		printf("z cluster %f, z cell %f, cluster-cell %f\n",pos[2], cellpos[2],pos[2]-cellpos[2]);
            //     				printf("r cluster %f, r cell %f, cluster-cell %f\n",r,      rcell,     r-rcell);
          }//PHOS and its matrices are available
          
          
        }// cluster cell loop
      }//Fill all position histograms

      // Get the fraction of the cluster energy that carries the cell with highest energy
      Float_t maxCellFraction = 0.;
      Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cell, clus,maxCellFraction);
      Double_t tmax  = cell->GetCellTime(absIdMax)*1e9;

      if     (clus->E() < 2.){
        fhLambda0vsClusterMaxCellDiffE0->Fill(clus->GetM02(),      maxCellFraction);
        fhNCellsvsClusterMaxCellDiffE0 ->Fill(nCaloCellsPerCluster,maxCellFraction);
      }
      else if(clus->E() < 6.){
        fhLambda0vsClusterMaxCellDiffE2->Fill(clus->GetM02(),      maxCellFraction);
        fhNCellsvsClusterMaxCellDiffE2 ->Fill(nCaloCellsPerCluster,maxCellFraction);
      }
      else{
        fhLambda0vsClusterMaxCellDiffE6->Fill(clus->GetM02(),      maxCellFraction);  
        fhNCellsvsClusterMaxCellDiffE6 ->Fill(nCaloCellsPerCluster,maxCellFraction);
      }
      
      fhNCellsPerClusterNoCut  ->Fill(clus->E(), nCaloCellsPerCluster);
      nModule = GetModuleNumber(clus);
      if(nModule >=0 && nModule < fNModules) fhNCellsPerClusterModNoCut[nModule]->Fill(clus->E(), nCaloCellsPerCluster);

      fhClusterMaxCellDiffNoCut->Fill(clus->E(),maxCellFraction);
      //fhClusterMaxCellDiffDivLambda0->Fill(clus->E(),maxCellFraction / clus->GetNCells());

      //Check bad clusters if rejection was not on
      Bool_t badCluster = kFALSE;
      if(fCalorimeter=="EMCAL" && !GetCaloUtils()->GetEMCALRecoUtils()->IsRejectExoticCluster()){
        //Bad clusters histograms
        //Float_t minNCells = TMath::Max(1,TMath::Nint(1 + TMath::Log(clus->E() - 5 )*1.5 ));
        //if(nCaloCellsPerCluster <= minNCells) {
        if(clus->GetM02() < 0.05) {

          //if(clus->GetM02() > 0 || TMath::Abs(clus->GetM20()) > 0 || clus->GetDispersion() > 0)
          
          Int_t sm =0; Int_t ietaa=-1; Int_t iphii = 0; Int_t rcu = 0;
          sm = GetModuleNumberCellIndexes(absId,fCalorimeter, ietaa, iphii, rcu);
//          if(clus->GetNCells() > 3){
//            printf("Bad  : E %f, ncells %d, nclusters %d, dist to bad %f, l0 %f, l1 %f, d %f, cell max t %f, cluster TOF %f, sm %d, icol %d, irow %d, rcu %d\n", 
//                   clus->E(), clus->GetNCells(),nCaloClusters, clus->GetDistanceToBadChannel(), 
//                   clus->GetM02(), clus->GetM20(), clus->GetDispersion(),tmax, tof,sm,ietaa,iphii,rcu);
//          }
          
          badCluster = kTRUE;
          
          fhBadClusterEnergy     ->Fill(clus->E());
          fhBadClusterMaxCellDiff->Fill(clus->E(),maxCellFraction);
          fhBadClusterTimeEnergy ->Fill(clus->E(),tof);
          //printf("bad tof : %2.3f\n",tof);
          //if(clus->E() - emax < 0)printf("What?\n");

          //Clusters in event time difference
          
            for(Int_t iclus2 = 0; iclus2 < nCaloClusters; iclus2++ ){
              
              AliVCluster* clus2 =  (AliVCluster*)caloClusters->At(iclus2);
              
              if(clus->GetID()==clus2->GetID()) continue;
              
              if(clus->GetM02() > 0.01) {
                fhBadClusterPairDiffTimeE  ->Fill(clus->E(), tof-clus2->GetTOF()*1.e9);
//               if(clus->GetNCells()>3) printf("\t i %d, E %f, nCells %d, dist to bad %f, good tof %f, bad tof %f, diff %f \n",
//                       iclus2, clus2->E(), clus2->GetNCells(), clus->GetDistanceToBadChannel(), tof,clus2->GetTOF()*1.e9,tof-clus2->GetTOF()*1.e9);
              }
//              else{  
//                Int_t absId2 = clus2->GetCellsAbsId()[0];
//                Int_t sm2 =0; Int_t ietaa2=-1; Int_t iphii2 = 0; Int_t rcu2 = 0;
//                sm2 = GetModuleNumberCellIndexes(absId2,fCalorimeter, ietaa2, iphii2, rcu2);
//                 if(clus->GetNCells()>3) printf("i %d, E %f, nCells %d, bad tof %f, bad tof %f, diff %f, sm %d, icol %d, irow %d, rcu %d\n",
//                       iclus2, clus2->E(), clus2->GetNCells(), tof,clus2->GetTOF()*1.e9,tof-clus2->GetTOF()*1.e9,sm2,ietaa2,iphii2,rcu2);
//                
//              }
              
            }
          
          
          for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) {
            //	printf("Index %d\n",ipos);  
            absId  = indexList[ipos]; 
            if(absId!=absIdMax){
              Float_t frac = cell->GetCellAmplitude(absId)/cell->GetCellAmplitude(absIdMax);
              //printf("bad frac : %2.3f, e %2.2f, ncells %d, min %2.1f\n",frac,mom.E(),nCaloCellsPerCluster,minNCells);
              fhBadClusterMaxCellCloseCellRatio->Fill(mom.E(),frac);
            }
          }
        }//Bad cluster
      }
      
      if(!badCluster){
        
        //        if(TMath::Abs(clus->GetM20()) < 0.0001 && clus->GetNCells() > 3){
        //          Int_t sm =0; Int_t ietaa=-1; Int_t iphii = 0; Int_t rcu = 0;
        //          sm = GetModuleNumberCellIndexes(absId,fCalorimeter, ietaa, iphii, rcu);
        //          printf("Good : E %f, mcells %d, l0 %f, l1 %f, d %f, cell max t %f, cluster TOF %f, sm %d, icol %d, irow %d \n", 
        //                 clus->E(), clus->GetNCells(),clus->GetM02(), clus->GetM20(), clus->GetDispersion(),tmax, tof,sm,ietaa,iphii);
        //
        //        }
        
        fhClusterMaxCellDiff->Fill(clus->E(),maxCellFraction);
        fhClusterTimeEnergy ->Fill(mom.E(),tof);
        
        //Clusters in event time difference
        for(Int_t iclus2 = 0; iclus2 < nCaloClusters; iclus2++ ){
          
          AliVCluster* clus2 =  (AliVCluster*) caloClusters->At(iclus2);
          
          if(clus->GetID()==clus2->GetID()) continue;
          
          if(clus->GetM02() > 0.01) {
            fhClusterPairDiffTimeE  ->Fill(clus->E(), tof-clus2->GetTOF()*1.e9);
          }
        }        
        
        for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) {
          //	printf("Index %d\n",ipos);            
          absId  = indexList[ipos]; 
          if(absId!=absIdMax){
            Float_t frac = cell->GetCellAmplitude(absId)/cell->GetCellAmplitude(absIdMax);
            //printf("good frac : %2.3f\n",frac);
            fhClusterMaxCellCloseCellRatio->Fill(mom.E(),frac);
          }
        }
        
        // check time of cells respect to max energy cell
        if(nCaloCellsPerCluster > 1 &&  GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
          for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) {
            absId  = indexList[ipos];             
            if(absId == absIdMax) continue;
            Float_t diff = (tmax-cell->GetCellTime(absId))*1e9;
            fhCellTimeSpreadRespectToCellMax->Fill(diff);
            if(TMath::Abs(TMath::Abs(diff) > 100)) fhCellIdCellLargeTimeSpread->Fill(absId);
          }// fill cell-cluster histogram loop
        }//check time of cells respect to max energy cell
        
        
        //Get module of cluster
        nCaloClustersAccepted++;
        if(nModule >=0 && nModule < fNModules) nClustersInModule[nModule]++;
        
        //-----------------------------------------------------------
        //Fill histograms related to single cluster or track matching
        //-----------------------------------------------------------
        ClusterHistograms(mom, pos, nCaloCellsPerCluster, nModule, nTracksMatched, track, labels, nLabel);	
        
        
        //-----------------------------------------------------------
        //Invariant mass
        //-----------------------------------------------------------
        if(fFillAllPi0Histo){
          if(GetDebug()>1) printf("Invariant mass \n");
          
          //do not do for bad vertex
          // Float_t fZvtxCut = 40. ;	
          if(v[2]<-GetZvertexCut() || v[2]> GetZvertexCut()) continue ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
          
          Int_t nModule2 = -1;
          Int_t nCaloCellsPerCluster2=0;
          if (nCaloClusters > 1 ) {
            for(Int_t jclus = iclus + 1 ; jclus < nCaloClusters ; jclus++) {
              AliVCluster* clus2 =  (AliVCluster*)caloClusters->At(jclus);
              
              //Get cluster kinematics
              clus2->GetMomentum(mom2,v);
              //Check only certain regions
              Bool_t in2 = kTRUE;
              if(IsFiducialCutOn()) in2 =  GetFiducialCut()->IsInFiducialCut(mom2,fCalorimeter) ;
              if(!in2) continue;	
              //Get module of cluster
              nModule2 = GetModuleNumber(clus2);
              //Cells per cluster
              nCaloCellsPerCluster2 = clus2->GetNCells();
            }
            //Fill invariant mass histograms
            //All modules
            
            //printf("QA : Fill inv mass histo: pt1 %f, pt2 %f, pt12 %f, mass %f, calo %s \n",mom.Pt(),mom2.Pt(),(mom+mom2).Pt(),(mom+mom2).M(), fCalorimeter.Data());
            fhIM  ->Fill((mom+mom2).Pt(),(mom+mom2).M());
            //Single module
            if(nModule == nModule2 && nModule >=0 && nModule < fNModules)
              fhIMMod[nModule]->Fill((mom+mom2).Pt(),(mom+mom2).M());
            
            //Select only clusters with at least 2 cells
            if(nCaloCellsPerCluster > 1 && nCaloCellsPerCluster2 > 1) {
              //All modules
              fhIMCellCut  ->Fill((mom+mom2).Pt(),(mom+mom2).M());
              //Single modules
              if(nModule == nModule2 && nModule >=0 && nModule < fNModules)
                fhIMCellCutMod[nModule]->Fill((mom+mom2).Pt(),(mom+mom2).M());
            }
            
            //Asymetry histograms
            fhAsym->Fill((mom+mom2).Pt(),TMath::Abs((mom.E()-mom2.E())/(mom.E()+mom2.E())));
            
          }// 2nd cluster loop
        }//Fill Pi0
        
      }//good cluster
      
    }//cluster loop

    //Number of clusters histograms
    if(nCaloClustersAccepted > 0) fhNClusters->Fill(nCaloClustersAccepted);
    //  Number of clusters per module
    for(Int_t imod = 0; imod < fNModules; imod++ ){ 
      if(GetDebug() > 1) 
        printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - module %d calo %s clusters %d\n", imod, fCalorimeter.Data(), nClustersInModule[imod]); 
      fhNClustersMod[imod]->Fill(nClustersInModule[imod]);
    }
    delete [] nClustersInModule;
    //delete caloClusters;
  }// calo clusters array exists
  
  //----------------------------------------------------------
  // CALOCELLS
  //----------------------------------------------------------
  
  AliVCaloCells * cell = 0x0; 
  Int_t ncells = 0;
  if(fCalorimeter == "PHOS") 
    cell = GetPHOSCells();
  else		              
    cell = GetEMCALCells();
  
  if(!cell){ 
    AliFatal(Form("No %s CELLS available for analysis",fCalorimeter.Data()));
    return; // just to trick coverity
  }
  
  if(GetDebug() > 0) 
    printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - %s cell entries %d\n", fCalorimeter.Data(), cell->GetNumberOfCells());    
  
  //Init arrays and used variables
  Int_t *nCellsInModule = new Int_t[fNModules];
  for(Int_t imod = 0; imod < fNModules; imod++ ) nCellsInModule[imod] = 0;
  Int_t icol     = -1;
  Int_t irow     = -1;
  Int_t iRCU     = -1;
  Float_t amp    = 0.;
  Float_t time   = 0.;
  Int_t id       = -1;
  Float_t recalF = 1.;  
  
  for (Int_t iCell = 0; iCell < cell->GetNumberOfCells(); iCell++) {      
    if(GetDebug() > 2)  
      printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Cell : amp %f, absId %d \n", cell->GetAmplitude(iCell), cell->GetCellNumber(iCell));
    nModule = GetModuleNumberCellIndexes(cell->GetCellNumber(iCell),fCalorimeter, icol, irow, iRCU);
    if(GetDebug() > 2) 
      printf("\t module %d, column %d, row %d \n", nModule,icol,irow);
    
    if(nModule < fNModules) {	
      
      //Check if the cell is a bad channel
      if(GetCaloUtils()->IsBadChannelsRemovalSwitchedOn()){
        if(fCalorimeter=="EMCAL"){
          if(GetCaloUtils()->GetEMCALChannelStatus(nModule,icol,irow)) continue;
        }
        else {
          if(GetCaloUtils()->GetPHOSChannelStatus(nModule,icol,irow)) {
            printf("PHOS bad channel\n");
            continue;
          }
        }
      } // use bad channel map
      
      //Get Recalibration factor if set
      if (GetCaloUtils()->IsRecalibrationOn()) {
        if(fCalorimeter == "PHOS") recalF = GetCaloUtils()->GetPHOSChannelRecalibrationFactor(nModule,icol,irow);
        else		                   recalF = GetCaloUtils()->GetEMCALChannelRecalibrationFactor(nModule,icol,irow);
        //if(fCalorimeter == "PHOS")printf("Recalibration factor (sm,row,col)=(%d,%d,%d) -  %f\n",nModule,icol,irow,recalF);
      }
      
      amp     = cell->GetAmplitude(iCell)*recalF;
      time    = cell->GetTime(iCell)*1e9;//transform time to ns
      
      //Remove noisy channels, only possible in ESDs
      if(GetReader()->GetDataType() == AliCaloTrackReader::kESD){
        if(time < fTimeCutMin || time > fTimeCutMax) continue;
      }
      //if(amp > 3 && fCalorimeter=="EMCAL") printf("Amp = %f, time = %f, (mod, col, row)= (%d,%d,%d)\n",
      //										   amp,time,nModule,icol,irow);
      
      id      = cell->GetCellNumber(iCell);
      fhAmplitude->Fill(amp);
      fhAmpId    ->Fill(amp,id);
      
      fhAmplitudeMod[nModule]->Fill(amp);
      if(fCalorimeter=="EMCAL"){
        Int_t ifrac = 0;
        if(icol > 15 && icol < 32) ifrac = 1;
        else if(icol > 31) ifrac = 2;
        fhAmplitudeModFraction[nModule*3+ifrac]->Fill(amp);
      }
      
      nCellsInModule[nModule]++;
      fhGridCellsMod[nModule]    ->Fill(icol,irow);
      fhGridCellsEMod[nModule]   ->Fill(icol,irow,amp);
      
      if(GetReader()->GetDataType() == AliCaloTrackReader::kESD){
        //printf("%s: time %g\n",fCalorimeter.Data(), time);
        fhTime     ->Fill(time);
        fhTimeId   ->Fill(time,id);
        fhTimeAmp  ->Fill(amp,time);
        
        //Double_t t0 = GetReader()->GetInputEvent()->GetT0();
        //printf("---->>> Time EMCal %e, T0 %e, T0 vertex %e, T0 clock %e, T0 trig %d \n",time,t0, 
        //	   GetReader()->GetInputEvent()->GetT0zVertex(),
        //	   GetReader()->GetInputEvent()->GetT0clock(),
        //	   GetReader()->GetInputEvent()->GetT0Trig());
        //fhT0Time     ->Fill(time-t0);
        //fhT0TimeId   ->Fill(time-t0,id);
        //fhT0TimeAmp  ->Fill(amp,time-t0);
        
        //printf("id %d, nModule %d, iRCU %d: Histo Name %s\n",id, nModule,iRCU, fhTimeAmpPerRCU[nModule*fNRCU+iRCU]->GetName());
        //fhT0TimeAmpPerRCU[nModule*fNRCU+iRCU]->Fill(amp, time-t0);
        
        fhTimeAmpPerRCU  [nModule*fNRCU+iRCU]->Fill(amp, time);
        
        if(amp > 0.3){
          fhGridCellsTimeMod[nModule]->Fill(icol,irow,time);
          
          //					AliESDCaloCells * cell2 = 0x0; 
          //					if(fCalorimeter == "PHOS") cell2 =  GetReader()->GetInputEvent()->GetPHOSCells();
          //					else		           cell2 = GetReader()->GetInputEvent()->GetEMCALCells();
          //					Int_t icol2    = -1;
          //					Int_t irow2    = -1;
          //					Int_t iRCU2    = -1;
          //					Float_t amp2   =  0.;
          //					Float_t time2  =  0.;
          //					Int_t id2      = -1;
          //					Int_t nModule2 = -1;
          //					for (Int_t iCell2 = 0; iCell2 < ncells; iCell2++) {  
          //						amp2    = cell2->GetAmplitude(iCell2);
          //						if(amp2 < 0.3) continue;
          //						if(iCell2 == iCell) continue;
          //						time2    = cell2->GetTime(iCell2)*1e9;//transform time to ns
          //						//printf("%s: time %g\n",fCalorimeter.Data(), time);
          //						id2      = cell2->GetCellNumber(iCell2);
          //						nModule2 = GetModuleNumberCellIndexes(cell2->GetCellNumber(iCell2), fCalorimeter, icol2, irow2, iRCU2);
          //						Int_t index = (nModule2*fNRCU+iRCU2)+(fNModules*fNRCU)*(iRCU+fNRCU*nModule); 
          //						//printf("id %d, nModule %d, iRCU %d, id2 %d, nModule2 %d, iRCU2 %d, index %d: Histo Name %s\n",id, nModule,iRCU,cell2->GetCellNumber(iCell2),nModule2,iRCU2,index, fhTimeCorrRCU[index]->GetName());
          //						fhTimeCorrRCU[index]->Fill(time,time2);	
          //						
          //					}// second cell loop
          
        }// amplitude cut
      }
      
      
      //Get Eta-Phi position of Cell
      if(fFillAllPosHisto)
      {
        if(fCalorimeter=="EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
          Float_t celleta = 0.;
          Float_t cellphi = 0.;
          GetEMCALGeometry()->EtaPhiFromIndex(id, celleta, cellphi); 
          
          fhEtaPhiAmp->Fill(celleta,cellphi,amp);
          Double_t cellpos[] = {0, 0, 0};
          GetEMCALGeometry()->GetGlobal(id, cellpos);
          fhXCellE->Fill(cellpos[0],amp)  ; 
          fhYCellE->Fill(cellpos[1],amp)  ; 
          fhZCellE->Fill(cellpos[2],amp)  ;
          Float_t rcell = TMath::Sqrt(cellpos[0]*cellpos[0]+cellpos[1]*cellpos[1]);//+cellpos[2]*cellpos[2]);
          fhRCellE->Fill(rcell,amp)  ;
          fhXYZCell->Fill(cellpos[0],cellpos[1],cellpos[2])  ;
        }//EMCAL Cells
        else if(fCalorimeter=="PHOS" && GetCaloUtils()->IsPHOSGeoMatrixSet()){
          TVector3 xyz;
          Int_t relId[4], module;
          Float_t xCell, zCell;
          
          GetPHOSGeometry()->AbsToRelNumbering(id,relId);
          module = relId[0];
          GetPHOSGeometry()->RelPosInModule(relId,xCell,zCell);
          GetPHOSGeometry()->Local2Global(module,xCell,zCell,xyz);
          Float_t rcell = TMath::Sqrt(xyz.X()*xyz.X()+xyz.Y()*xyz.Y());
          fhXCellE ->Fill(xyz.X(),amp)  ; 
          fhYCellE ->Fill(xyz.Y(),amp)  ; 
          fhZCellE ->Fill(xyz.Z(),amp)  ;
          fhRCellE ->Fill(rcell  ,amp)  ;
          fhXYZCell->Fill(xyz.X(),xyz.Y(),xyz.Z())  ;
        }//PHOS cells
      }//fill cell position histograms
      
      if     (fCalorimeter=="EMCAL" && amp > fEMCALCellAmpMin) ncells ++ ;
      else if(fCalorimeter=="PHOS"  && amp > fPHOSCellAmpMin)  ncells ++ ;
      //else  
      //  printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - no %s CELLS passed the analysis cut\n",fCalorimeter.Data());    
    }//nmodules
  }//cell loop
  if(ncells > 0 )fhNCells->Fill(ncells) ; //fill the cells after the cut 
  
  //Number of cells per module
  for(Int_t imod = 0; imod < fNModules; imod++ ) {
    if(GetDebug() > 1) 
      printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - module %d calo %s cells %d\n", imod, fCalorimeter.Data(), nCellsInModule[imod]); 
    fhNCellsMod[imod]->Fill(nCellsInModule[imod]) ;
  }
  delete [] nCellsInModule;
  
  if(GetDebug() > 0)
    printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - End \n");
}


//_____________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterHistograms(const TLorentzVector mom, 
                                            Float_t *pos, const Int_t nCaloCellsPerCluster,const Int_t nModule,
                                            const Int_t nTracksMatched,  const AliVTrack * track,  
                                            const Int_t * labels, const Int_t nLabels){
  //Fill CaloCluster related histograms
	
  AliAODMCParticle * aodprimary  = 0x0;
  TParticle * primary = 0x0;
  Int_t tag = 0;	
  
  Float_t e   = mom.E();
  Float_t pt  = mom.Pt();
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  if(GetDebug() > 0) {
    printf("AliAnaCalorimeterQA::ClusterHistograms() - cluster: E %2.3f, pT %2.3f, eta %2.3f, phi %2.3f \n",e,pt,eta,phi*TMath::RadToDeg());
    if(IsDataMC()) {
      //printf("\t Primaries: nlabels %d, labels pointer %p\n",nLabels,labels);
      printf("\t Primaries: nlabels %d\n",nLabels);
      if(!nLabels || !labels) printf("\t Strange, no labels!!!\n");
    }
  }
  
  fhE     ->Fill(e);	
  if(nModule >=0 && nModule < fNModules) fhEMod[nModule]->Fill(e);
  if(fFillAllTH12){
    fhPt     ->Fill(pt);
    fhPhi    ->Fill(phi);
    fhEta    ->Fill(eta);
  }
  
  fhEtaPhiE->Fill(eta,phi,e);
  
  //Cells per cluster
  fhNCellsPerCluster   ->Fill(e, nCaloCellsPerCluster);
  if((fCalorimeter=="EMCAL" && GetReader()->GetEMCALPtMin() < 0.3) ||
     (fCalorimeter=="PHOS"  && GetReader()->GetPHOSPtMin()  < 0.3)) fhNCellsPerClusterMIP->Fill(e, nCaloCellsPerCluster);
  
  //Position
  if(fFillAllPosHisto2){
    fhXE     ->Fill(pos[0],e);
    fhYE     ->Fill(pos[1],e);
    fhZE     ->Fill(pos[2],e);
    if(fFillAllTH3)
      fhXYZ    ->Fill(pos[0], pos[1],pos[2]);
    
    fhXNCells->Fill(pos[0],nCaloCellsPerCluster);
    fhYNCells->Fill(pos[1],nCaloCellsPerCluster);
    fhZNCells->Fill(pos[2],nCaloCellsPerCluster);
    Float_t rxyz = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]);//+pos[2]*pos[2]);
    fhRE     ->Fill(rxyz,e);
    fhRNCells->Fill(rxyz  ,nCaloCellsPerCluster);
  }
  
  if(nModule >=0 && nModule < fNModules) fhNCellsPerClusterMod[nModule]->Fill(e, nCaloCellsPerCluster);
  
  //Fill histograms only possible when simulation
  if(IsDataMC() && nLabels > 0 && labels){
    
    //Play with the MC stack if available
    Int_t label = labels[0];
    
    if(label < 0) {
      if(GetDebug() >= 0) printf("AliAnaCalorimeterQA::ClusterHistograms() *** bad label ***:  label %d \n", label);
      return;
    }
    
    Int_t pdg  =-1; Int_t pdg0  =-1;Int_t status = -1; Int_t iMother = -1; Int_t iParent = -1;
    Float_t vxMC= 0; Float_t vyMC = 0;	
    Float_t eMC = 0; Float_t ptMC= 0; Float_t phiMC =0; Float_t etaMC = 0;
    Int_t charge = 0;	
    
    //Check the origin.
    tag = GetMCAnalysisUtils()->CheckOrigin(labels,nLabels, GetReader(),0);
    
    if(GetReader()->ReadStack() && !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown)){ //it MC stack and known tag
      
      if( label >= GetMCStack()->GetNtrack()) {
        if(GetDebug() >= 0) printf("AliAnaCalorimeterQA::ClusterHistograms() *** large label ***:  label %d, n tracks %d \n", label, GetMCStack()->GetNtrack());
        return ;
      }
      
      primary = GetMCStack()->Particle(label);
      iMother = label;
      pdg0    = TMath::Abs(primary->GetPdgCode());
      pdg     = pdg0;
      status  = primary->GetStatusCode();
      vxMC    = primary->Vx();
      vyMC    = primary->Vy();
      iParent = primary->GetFirstMother();
      
      if(GetDebug() > 1 ) {
        printf("AliAnaCalorimeterQA::ClusterHistograms() - Cluster most contributing mother: \n");
        printf("\t Mother label %d, pdg %d, %s, status %d, parent %d \n",iMother, pdg0, primary->GetName(),status, iParent);
      }
      
      //Get final particle, no conversion products
      if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)){
        //Get the parent
        primary = GetMCStack()->Particle(iParent);
        pdg = TMath::Abs(primary->GetPdgCode());
        if(GetDebug() > 1 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted cluster!. Find before conversion: \n");
        while((pdg == 22 || pdg == 11) && status != 1){
          iMother = iParent;
          primary = GetMCStack()->Particle(iMother);
          status  = primary->GetStatusCode();
          iParent = primary->GetFirstMother();
          pdg     = TMath::Abs(primary->GetPdgCode());
          if(GetDebug() > 1 )printf("\t pdg %d, index %d, %s, status %d \n",pdg, iMother,  primary->GetName(),status);	
        }	
        
        if(GetDebug() > 1 ) {
          printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted Cluster mother before conversion: \n");
          printf("\t Mother label %d, pdg %d, %s, status %d, parent %d \n",iMother, pdg, primary->GetName(), status, iParent);
        }
        
      }
      
      //Overlapped pi0 (or eta, there will be very few), get the meson
      if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
         GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)){
        if(GetDebug() > 1 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped Meson decay!, Find it: \n");
        while(pdg != 111 && pdg != 221){
          iMother = iParent;
          primary = GetMCStack()->Particle(iMother);
          status  = primary->GetStatusCode();
          iParent = primary->GetFirstMother();
          pdg     = TMath::Abs(primary->GetPdgCode());
          if(GetDebug() > 1 ) printf("\t pdg %d, %s, index %d\n",pdg,  primary->GetName(),iMother);
          if(iMother==-1) {
            printf("AliAnaCalorimeterQA::ClusterHistograms() - Tagged as Overlapped photon but meson not found, why?\n");
            //break;
          }
        }
        
        if(GetDebug() > 2 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped %s decay, label %d \n", 
                                   primary->GetName(),iMother);
      }
      
      eMC    = primary->Energy();
      ptMC   = primary->Pt();
      phiMC  = primary->Phi();
      etaMC  = primary->Eta();
      pdg    = TMath::Abs(primary->GetPdgCode());
      charge = (Int_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      
    }
    else if(GetReader()->ReadAODMCParticles() && !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown)){//it MC AOD and known tag
      //Get the list of MC particles
      if(!GetReader()->GetAODMCParticles(0)) 
        AliFatal("MCParticles not available!");
      
      aodprimary = (AliAODMCParticle*) (GetReader()->GetAODMCParticles(0))->At(label);
      iMother = label;
      pdg0    = TMath::Abs(aodprimary->GetPdgCode());
      pdg     = pdg0;
      status  = aodprimary->IsPrimary();
      vxMC    = aodprimary->Xv();
      vyMC    = aodprimary->Yv();
      iParent = aodprimary->GetMother();
      
      if(GetDebug() > 1 ) {
        printf("AliAnaCalorimeterQA::ClusterHistograms() - Cluster most contributing mother: \n");
        printf("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d, parent %d \n",
               iMother, pdg0, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary(), iParent);
      }
      
      //Get final particle, no conversion products
      if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)){
        if(GetDebug() > 1 ) 
          printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted cluster!. Find before conversion: \n");
        //Get the parent
        aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles(0))->At(iParent);
        pdg = TMath::Abs(aodprimary->GetPdgCode());
        while ((pdg == 22 || pdg == 11) && !aodprimary->IsPhysicalPrimary()) {
          iMother    = iParent;
          aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles(0))->At(iMother);
          status     = aodprimary->IsPrimary();
          iParent    = aodprimary->GetMother();
          pdg        = TMath::Abs(aodprimary->GetPdgCode());
          if(GetDebug() > 1 )
            printf("\t pdg %d, index %d, Primary? %d, Physical Primary? %d \n",
                   pdg, iMother, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary());	
        }	
        
        if(GetDebug() > 1 ) {
          printf("AliAnaCalorimeterQA::ClusterHistograms() - Converted Cluster mother before conversion: \n");
          printf("\t Mother label %d, pdg %d, parent %d, Primary? %d, Physical Primary? %d \n",
                 iMother, pdg, iParent, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary());
        }
        
      }
      
      //Overlapped pi0 (or eta, there will be very few), get the meson
      if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
         GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)){
        if(GetDebug() > 1 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped Meson decay!, Find it: PDG %d, mom %d \n",pdg, iMother);
        while(pdg != 111 && pdg != 221){
          
          iMother    = iParent;
          aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles(0))->At(iMother);
          status     = aodprimary->IsPrimary();
          iParent    = aodprimary->GetMother();
          pdg        = TMath::Abs(aodprimary->GetPdgCode());
          
          if(GetDebug() > 1 ) printf("\t pdg %d, index %d\n",pdg, iMother);
          
          if(iMother==-1) {
            printf("AliAnaCalorimeterQA::ClusterHistograms() - Tagged as Overlapped photon but meson not found, why?\n");
            //break;
          }
        }	
        
        if(GetDebug() > 2 ) printf("AliAnaCalorimeterQA::ClusterHistograms() - Overlapped %s decay, label %d \n", 
                                   aodprimary->GetName(),iMother);
      }	
      
      status = aodprimary->IsPrimary();
      eMC    = aodprimary->E();
      ptMC   = aodprimary->Pt();
      phiMC  = aodprimary->Phi();
      etaMC  = aodprimary->Eta();
      pdg    = TMath::Abs(aodprimary->GetPdgCode());
      charge = aodprimary->Charge();
      
    }
    
    //Float_t vz = primary->Vz();
    Float_t rVMC = TMath::Sqrt(vxMC*vxMC + vyMC*vyMC);
    if((pdg == 22 || TMath::Abs(pdg)==11) && status!=1) {
      fhEMVxyz   ->Fill(vxMC,vyMC);//,vz);
      fhEMR      ->Fill(e,rVMC);
    }
    
    //printf("reco e %f, pt %f, phi %f, eta %f \n", e, pt, phi, eta);
    //printf("prim e %f, pt %f, phi %f, eta %f \n", eMC,ptMC,phiMC ,etaMC );
    //printf("vertex: vx %f, vy %f, vz %f, r %f \n", vxMC, vyMC, vz, r);
    
    
    fh2E      ->Fill(e, eMC);
    fh2Pt     ->Fill(pt, ptMC);
    fh2Phi    ->Fill(phi, phiMC);
    fh2Eta    ->Fill(eta, etaMC);
    fhDeltaE  ->Fill(eMC-e);
    fhDeltaPt ->Fill(ptMC-pt);
    fhDeltaPhi->Fill(phiMC-phi);
    fhDeltaEta->Fill(etaMC-eta);
    if(eMC   > 0) fhRatioE  ->Fill(e/eMC);
    if(ptMC  > 0) fhRatioPt ->Fill(pt/ptMC);
    if(phiMC > 0) fhRatioPhi->Fill(phi/phiMC);
    if(etaMC > 0) fhRatioEta->Fill(eta/etaMC);			
    
    
    //Overlapped pi0 (or eta, there will be very few)
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
       GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)){
      fhPi0E     ->Fill(e,eMC);	
      fhPi0Pt    ->Fill(pt,ptMC);
      fhPi0Eta   ->Fill(eta,etaMC);	
      fhPi0Phi   ->Fill(phi,phiMC);
      if( nTracksMatched > 0){
        fhPi0ECharged     ->Fill(e,eMC);		
        fhPi0PtCharged    ->Fill(pt,ptMC);
        fhPi0PhiCharged   ->Fill(phi,phiMC);
        fhPi0EtaCharged   ->Fill(eta,etaMC);
      }
    }//Overlapped pizero decay
    else if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)){
      fhGamE     ->Fill(e,eMC);	
      fhGamPt    ->Fill(pt,ptMC);
      fhGamEta   ->Fill(eta,etaMC);	
      fhGamPhi   ->Fill(phi,phiMC);
      fhGamDeltaE  ->Fill(eMC-e);
      fhGamDeltaPt ->Fill(ptMC-pt);	
      fhGamDeltaPhi->Fill(phiMC-phi);
      fhGamDeltaEta->Fill(etaMC-eta);
      if(eMC > 0) fhGamRatioE  ->Fill(e/eMC);
      if(ptMC     > 0) fhGamRatioPt ->Fill(pt/ptMC);
      if(phiMC    > 0) fhGamRatioPhi->Fill(phi/phiMC);
      if(etaMC    > 0) fhGamRatioEta->Fill(eta/etaMC);
      if( nTracksMatched > 0){
        fhGamECharged     ->Fill(e,eMC);		
        fhGamPtCharged    ->Fill(pt,ptMC);
        fhGamPhiCharged   ->Fill(phi,phiMC);
        fhGamEtaCharged   ->Fill(eta,etaMC);
      }
    }//photon
    else if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron)){
      fhEleE     ->Fill(e,eMC);	
      fhElePt    ->Fill(pt,ptMC);
      fhEleEta   ->Fill(eta,etaMC);	
      fhElePhi   ->Fill(phi,phiMC);
      fhEMVxyz   ->Fill(vxMC,vyMC);//,vz);
      fhEMR      ->Fill(e,rVMC);
      if( nTracksMatched > 0){
        fhEleECharged     ->Fill(e,eMC);		
        fhElePtCharged    ->Fill(pt,ptMC);
        fhElePhiCharged   ->Fill(phi,phiMC);
        fhEleEtaCharged   ->Fill(eta,etaMC);
      }
    }
    else if(charge == 0){
      fhNeHadE     ->Fill(e,eMC);	
      fhNeHadPt    ->Fill(pt,ptMC);
      fhNeHadEta   ->Fill(eta,etaMC);	
      fhNeHadPhi   ->Fill(phi,phiMC);	
      fhHaVxyz     ->Fill(vxMC,vyMC);//,vz);
      fhHaR        ->Fill(e,rVMC);
      if( nTracksMatched > 0){
        fhNeHadECharged     ->Fill(e,eMC);		
        fhNeHadPtCharged    ->Fill(pt,ptMC);
        fhNeHadPhiCharged   ->Fill(phi,phiMC);
        fhNeHadEtaCharged   ->Fill(eta,etaMC);
      }
    }
    else if(charge!=0){
      fhChHadE     ->Fill(e,eMC);	
      fhChHadPt    ->Fill(pt,ptMC);
      fhChHadEta   ->Fill(eta,etaMC);	
      fhChHadPhi   ->Fill(phi,phiMC);	
      fhHaVxyz     ->Fill(vxMC,vyMC);//,vz);
      fhHaR        ->Fill(e,rVMC);
      if( nTracksMatched > 0){
        fhChHadECharged     ->Fill(e,eMC);		
        fhChHadPtCharged    ->Fill(pt,ptMC);
        fhChHadPhiCharged   ->Fill(phi,phiMC);
        fhChHadEtaCharged   ->Fill(eta,etaMC);
      }
    }
  }//Work with MC
  
	
  //Match tracks and clusters
  //To be Modified in case of AODs
  
  if( nTracksMatched > 0 &&  fFillAllTMHisto){
    if(fFillAllTH12 && fFillAllTMHisto){
      fhECharged      ->Fill(e);	
      fhPtCharged     ->Fill(pt);
      fhPhiCharged    ->Fill(phi);
      fhEtaCharged    ->Fill(eta);
    }
    
    if(fFillAllTMHisto){
      if(fFillAllTH3)fhEtaPhiECharged->Fill(eta,phi,e);		
      if((fCalorimeter=="EMCAL" && GetReader()->GetEMCALPtMin() < 0.3) ||
         (fCalorimeter=="PHOS"  && GetReader()->GetPHOSPtMin()  < 0.3))   fhNCellsPerClusterMIPCharged->Fill(e, nCaloCellsPerCluster);
    }
    //printf("track index %d ntracks %d\n", esd->GetNumberOfTracks());	
    //Study the track and matched cluster if track exists.
    if(!track) return;
    Double_t emcpos[3] = {0.,0.,0.};
    Double_t emcmom[3] = {0.,0.,0.};
    Double_t radius    = 441.0; //[cm] EMCAL radius +13cm
    Double_t bfield    = 0.;
    Double_t tphi      = 0;
    Double_t teta      = 0;
    Double_t tmom      = 0;
    Double_t tpt       = 0;
    Double_t tmom2     = 0;
    Double_t tpcSignal = 0;
    Bool_t okpos = kFALSE;
    Bool_t okmom = kFALSE;
    Bool_t okout = kFALSE;
    Int_t nITS   = 0;
    Int_t nTPC   = 0;
    
    //In case of ESDs get the parameters in this way
    if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
      if (track->GetOuterParam() ) {
        okout = kTRUE;
        
        bfield = GetReader()->GetInputEvent()->GetMagneticField();
        okpos = track->GetOuterParam()->GetXYZAt(radius,bfield,emcpos);
        okmom = track->GetOuterParam()->GetPxPyPzAt(radius,bfield,emcmom);
        if(!(okpos && okmom)) return;
        
        TVector3 position(emcpos[0],emcpos[1],emcpos[2]);
        TVector3 momentum(emcmom[0],emcmom[1],emcmom[2]);
        tphi = position.Phi();
        teta = position.Eta();
        tmom = momentum.Mag();
        
        //Double_t tphi  = track->GetOuterParam()->Phi();
        //Double_t teta  = track->GetOuterParam()->Eta();
        //Double_t tmom  = track->GetOuterParam()->P();
        tpt       = track->Pt();
        tmom2     = track->P();
        tpcSignal = track->GetTPCsignal();
        
        nITS = track->GetNcls(0);
        nTPC = track->GetNcls(1);
      }//Outer param available 
    }// ESDs
    else if(GetReader()->GetDataType()==AliCaloTrackReader::kAOD) {
      AliAODPid* pid = (AliAODPid*) ((AliAODTrack *) track)->GetDetPid();
      if (pid) {
        okout = kTRUE;
        pid->GetEMCALPosition(emcpos);
        pid->GetEMCALMomentum(emcmom);	
        
        TVector3 position(emcpos[0],emcpos[1],emcpos[2]);
        TVector3 momentum(emcmom[0],emcmom[1],emcmom[2]);
        tphi = position.Phi();
        teta = position.Eta();
        tmom = momentum.Mag();
        
        tpt       = track->Pt();
        tmom2     = track->P();
        tpcSignal = pid->GetTPCsignal();
        
        //nITS = ((AliAODTrack*)track)->GetNcls(0);
        //nTPC = ((AliAODTrack*)track)->GetNcls(1);
      }//pid 
    }//AODs
		
    if(okout){
      Double_t deta = teta - eta;
      Double_t dphi = tphi - phi;
      if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
      if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
      Double_t dR = sqrt(dphi*dphi + deta*deta);
			
      Double_t pOverE = tmom/e;
			
      fh1pOverE->Fill(tpt, pOverE);
      if(dR < 0.02) fh1pOverER02->Fill(tpt,pOverE);
			
      fh1dR->Fill(dR);
      fh2MatchdEdx->Fill(tmom2,tpcSignal);
			
      if(IsDataMC() && primary){ 
        Int_t pdg = primary->GetPdgCode();
        Double_t  charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
				
        if(TMath::Abs(pdg) == 11){
          fhMCEle1pOverE->Fill(tpt,pOverE);
          fhMCEle1dR->Fill(dR);
          fhMCEle2MatchdEdx->Fill(tmom2,tpcSignal);		
          if(dR < 0.02) fhMCEle1pOverER02->Fill(tpt,pOverE);
        }
        else if(charge!=0){
          fhMCChHad1pOverE->Fill(tpt,pOverE);
          fhMCChHad1dR->Fill(dR);
          fhMCChHad2MatchdEdx->Fill(tmom2,tpcSignal);	
          if(dR < 0.02) fhMCChHad1pOverER02->Fill(tpt,pOverE);
        }
        else if(charge == 0){
          fhMCNeutral1pOverE->Fill(tpt,pOverE);
          fhMCNeutral1dR->Fill(dR);
          fhMCNeutral2MatchdEdx->Fill(tmom2,tpcSignal);	
          if(dR < 0.02) fhMCNeutral1pOverER02->Fill(tpt,pOverE);
        }
      }//DataMC
      
      if(dR < 0.02 && pOverE > 0.5 && pOverE < 1.5
         && nCaloCellsPerCluster > 1 && nITS > 3 && nTPC > 20) {
        fh2EledEdx->Fill(tmom2,tpcSignal);
      }
    }
    else{//no ESD external param or AODPid
      //					ULong_t status=AliESDtrack::kTPCrefit;
      //				status|=AliESDtrack::kITSrefit;
      //printf("track status %d\n", track->GetStatus() );
      //				fhEChargedNoOut      ->Fill(e);		
      //				fhPtChargedNoOut     ->Fill(pt);
      //				fhPhiChargedNoOut    ->Fill(phi);
      //				fhEtaChargedNoOut    ->Fill(eta);
      //				fhEtaPhiChargedNoOut ->Fill(eta,phi);
      //				if(GetDebug() >= 0 && ((track->GetStatus() & status) == status)) printf("ITS+TPC\n");
      if(GetDebug() >= 0) printf("No ESD external param or AliAODPid \n");
      
    }//No out params
  }//matched clusters with tracks
  
}// Clusters


//__________________________________
void AliAnaCalorimeterQA::Correlate(){
  // Correlate information from PHOS and EMCAL and with V0 and track multiplicity
  
  //Clusters 
  TObjArray * caloClustersEMCAL = GetEMCALClusters();
  TObjArray * caloClustersPHOS  = GetPHOSClusters();
  
  Int_t nclEMCAL = caloClustersEMCAL->GetEntriesFast();
  Int_t nclPHOS  = caloClustersPHOS ->GetEntriesFast();
  
  Float_t sumClusterEnergyEMCAL = 0;
  Float_t sumClusterEnergyPHOS  = 0;
  Int_t iclus = 0;
  for(iclus = 0 ; iclus <  caloClustersEMCAL->GetEntriesFast() ; iclus++) 
    sumClusterEnergyEMCAL += ((AliVCluster*)caloClustersEMCAL->At(iclus))->E();
  for(iclus = 0 ; iclus <  caloClustersPHOS->GetEntriesFast(); iclus++) 
    sumClusterEnergyPHOS += ((AliVCluster*)caloClustersPHOS->At(iclus))->E();
  
  
  //Cells
  
  AliVCaloCells * cellsEMCAL = GetEMCALCells();
  AliVCaloCells * cellsPHOS  = GetPHOSCells();
  
  Int_t ncellsEMCAL = cellsEMCAL->GetNumberOfCells();
  Int_t ncellsPHOS  = cellsPHOS ->GetNumberOfCells();
  
  Float_t sumCellEnergyEMCAL = 0;
  Float_t sumCellEnergyPHOS  = 0;
  Int_t icell = 0;
  for(icell = 0 ; icell < cellsEMCAL->GetNumberOfCells()  ; icell++) 
    sumCellEnergyEMCAL += cellsEMCAL->GetAmplitude(icell);
  for(icell = 0 ; icell <  cellsPHOS->GetNumberOfCells(); icell++) 
    sumCellEnergyPHOS += cellsPHOS->GetAmplitude(icell);
  
  
  //Fill Histograms
  fhCaloCorrNClusters->Fill(nclEMCAL,nclPHOS);
  fhCaloCorrEClusters->Fill(sumClusterEnergyEMCAL,sumClusterEnergyPHOS);
  fhCaloCorrNCells   ->Fill(ncellsEMCAL,ncellsPHOS);
  fhCaloCorrECells   ->Fill(sumCellEnergyEMCAL,sumCellEnergyPHOS);
  
  Int_t v0S = GetV0Signal(0)+GetV0Signal(1);
  Int_t v0M = GetV0Multiplicity(0)+GetV0Multiplicity(1);
  Int_t trM = GetTrackMultiplicity();
  if(fCalorimeter=="PHOS"){
    fhCaloV0MCorrNClusters   ->Fill(v0M,nclPHOS);
    fhCaloV0MCorrEClusters   ->Fill(v0M,sumClusterEnergyPHOS);
    fhCaloV0MCorrNCells      ->Fill(v0M,ncellsPHOS);
    fhCaloV0MCorrECells      ->Fill(v0M,sumCellEnergyPHOS);
    
    fhCaloV0SCorrNClusters   ->Fill(v0S,nclPHOS);
    fhCaloV0SCorrEClusters   ->Fill(v0S,sumClusterEnergyPHOS);
    fhCaloV0SCorrNCells      ->Fill(v0S,ncellsPHOS);
    fhCaloV0SCorrECells      ->Fill(v0S,sumCellEnergyPHOS);
    
    fhCaloTrackMCorrNClusters->Fill(trM,nclPHOS);
    fhCaloTrackMCorrEClusters->Fill(trM,sumClusterEnergyPHOS);    
    fhCaloTrackMCorrNCells   ->Fill(trM,ncellsPHOS);
    fhCaloTrackMCorrECells   ->Fill(trM,sumCellEnergyPHOS);
  }
  else{
    fhCaloV0MCorrNClusters   ->Fill(v0M,nclEMCAL);
    fhCaloV0MCorrEClusters   ->Fill(v0M,sumClusterEnergyEMCAL);
    fhCaloV0MCorrNCells      ->Fill(v0M,ncellsEMCAL);
    fhCaloV0MCorrECells      ->Fill(v0M,sumCellEnergyEMCAL);
    
    fhCaloV0SCorrNClusters   ->Fill(v0S,nclEMCAL);
    fhCaloV0SCorrEClusters   ->Fill(v0S,sumClusterEnergyEMCAL);
    fhCaloV0SCorrNCells      ->Fill(v0S,ncellsEMCAL);
    fhCaloV0SCorrECells      ->Fill(v0S,sumCellEnergyEMCAL);
    
    fhCaloTrackMCorrNClusters->Fill(trM,nclEMCAL);
    fhCaloTrackMCorrEClusters->Fill(trM,sumClusterEnergyEMCAL);    
    fhCaloTrackMCorrNCells   ->Fill(trM,ncellsEMCAL);
    fhCaloTrackMCorrECells   ->Fill(trM,sumCellEnergyEMCAL);
  }
  
  if(GetDebug() > 0 )
  {
    printf("AliAnaCalorimeterQA::Correlate(): \n");
    printf("\t EMCAL: N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f \n",
           ncellsEMCAL,nclEMCAL, sumCellEnergyEMCAL,sumClusterEnergyEMCAL);
    printf("\t PHOS : N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f \n",
           ncellsPHOS,nclPHOS,sumCellEnergyPHOS,sumClusterEnergyPHOS);
    printf("\t V0 : Signal %d, Multiplicity  %d, Track Multiplicity %d \n", v0S,v0M,trM);
  }
}


//______________________________________________________________________________
void AliAnaCalorimeterQA::MCHistograms(const TLorentzVector mom, const Int_t pdg){
  //Fill pure monte carlo related histograms
	
  Float_t eMC    = mom.E();
  Float_t ptMC   = mom.Pt();
  Float_t phiMC  = mom.Phi();
  if(phiMC < 0) 
    phiMC  += TMath::TwoPi();
  Float_t etaMC  = mom.Eta();
  
  if (TMath::Abs(etaMC) > 1) return;
  
  Bool_t in = kTRUE;
  if(IsFiducialCutOn()) in =  GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
  
  if (pdg==22) {
    fhGenGamPt ->Fill(ptMC);
    fhGenGamEta->Fill(etaMC);
    fhGenGamPhi->Fill(phiMC);
    if(in){
      fhGenGamAccE  ->Fill(eMC);
      fhGenGamAccPt ->Fill(ptMC);
      fhGenGamAccEta->Fill(etaMC);
      fhGenGamAccPhi->Fill(phiMC);					
    }
  }
  else if (pdg==111) {
    fhGenPi0Pt ->Fill(ptMC);
    fhGenPi0Eta->Fill(etaMC);
    fhGenPi0Phi->Fill(phiMC);
    if(in){
      fhGenPi0AccE  ->Fill(eMC);					
      fhGenPi0AccPt ->Fill(ptMC);
      fhGenPi0AccEta->Fill(etaMC);
      fhGenPi0AccPhi->Fill(phiMC);					
    }
  }
  else if (pdg==221) {
    fhGenEtaPt ->Fill(ptMC);
    fhGenEtaEta->Fill(etaMC);
    fhGenEtaPhi->Fill(phiMC);
  }
  else if (pdg==223) {
    fhGenOmegaPt ->Fill(ptMC);
    fhGenOmegaEta->Fill(etaMC);
    fhGenOmegaPhi->Fill(phiMC);
  }
  else if (TMath::Abs(pdg)==11) {
    fhGenElePt ->Fill(ptMC);
    fhGenEleEta->Fill(etaMC);
    fhGenElePhi->Fill(phiMC);
  }	
  
}

//________________________________________________________________________
void AliAnaCalorimeterQA::ReadHistograms(TList* outputList)
{
  // Needed when Terminate is executed in distributed environment
  // Refill analysis histograms of this class with corresponding histograms in output list. 
	
  // Histograms of this analsys are kept in the same list as other analysis, recover the position of
  // the first one and then add the next 
  Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"hE"));
  //printf("Calo: %s, index: %d, nmodules %d\n",fCalorimeter.Data(),index,fNModules);
  
  //Read histograms, must be in the same order as in GetCreateOutputObject.
  fhE       = (TH1F *) outputList->At(index++); 	
  if(fFillAllTH12){
    fhPt      = (TH1F *) outputList->At(index++); 
    fhPhi     = (TH1F *) outputList->At(index++); 
    fhEta     = (TH1F *) outputList->At(index++);
  }
  
  fhEtaPhiE = (TH3F *) outputList->At(index++);
  
  fhClusterTimeEnergy = (TH2F*) outputList->At(index++);
  
  if(fFillAllTH12){
    fhECharged       = (TH1F *) outputList->At(index++); 	
    fhPtCharged      = (TH1F *) outputList->At(index++); 
    fhPhiCharged     = (TH1F *) outputList->At(index++); 
    fhEtaCharged     = (TH1F *) outputList->At(index++);
  }
  
  fhEtaPhiECharged = (TH3F *) outputList->At(index++);
  
  fh1pOverE =    (TH2F *) outputList->At(index++);
  fh1dR =        (TH1F *) outputList->At(index++);
  fh2MatchdEdx = (TH2F *) outputList->At(index++);
  fh2EledEdx =   (TH2F *) outputList->At(index++);
  fh1pOverER02 = (TH2F *) outputList->At(index++);
  
  fhIM        = (TH2F *) outputList->At(index++);
  fhIMCellCut = (TH2F *) outputList->At(index++);
  fhAsym      = (TH2F *) outputList->At(index++);
  
  fhNCellsPerCluster           = (TH2F *) outputList->At(index++);
  fhNCellsPerClusterMIP        = (TH2F *) outputList->At(index++);
  fhNCellsPerClusterMIPCharged = (TH2F *) outputList->At(index++);
  fhNClusters  = (TH1F *) outputList->At(index++); 
  
  fhRNCells = (TH2F *) outputList->At(index++);
  fhXNCells = (TH2F *) outputList->At(index++);
  fhYNCells = (TH2F *) outputList->At(index++);
  fhZNCells = (TH2F *) outputList->At(index++);
  fhRE	  = (TH2F *) outputList->At(index++);
  fhXE	  = (TH2F *) outputList->At(index++);
  fhYE	  = (TH2F *) outputList->At(index++);
  fhZE	  = (TH2F *) outputList->At(index++); 
  fhXYZ	  = (TH3F *) outputList->At(index++);
  if(fFillAllPosHisto){
    fhRCellE	  = (TH2F *) outputList->At(index++);
    fhXCellE	  = (TH2F *) outputList->At(index++);
    fhYCellE	  = (TH2F *) outputList->At(index++);
    fhZCellE	  = (TH2F *) outputList->At(index++); 
    fhXYZCell	  = (TH3F *) outputList->At(index++); 
    fhDeltaCellClusterRNCells = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterXNCells = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterYNCells = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterZNCells = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterRE	  = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterXE	  = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterYE	  = (TH2F *) outputList->At(index++);
    fhDeltaCellClusterZE	  = (TH2F *) outputList->At(index++); 
    fhEtaPhiAmp               = (TH3F *) outputList->At(index++); 
  }
  
  fhNCells     = (TH1F *) outputList->At(index++); 
  fhAmplitude  = (TH1F *) outputList->At(index++); 
  fhAmpId      = (TH2F *) outputList->At(index++); 
  
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
    
    fhCellTimeSpreadRespectToCellMax = (TH1F *) outputList->At(index++);
    fhCellIdCellLargeTimeSpread      = (TH1F *) outputList->At(index++);
    
    fhTime       = (TH1F *) outputList->At(index++); 
    fhTimeId     = (TH2F *) outputList->At(index++); 
    fhTimeAmp    = (TH2F *) outputList->At(index++); 
    
    //		fhT0Time       = (TH1F *) outputList->At(index++); 
    //		fhT0TimeId     = (TH2F *) outputList->At(index++); 
    //		fhT0TimeAmp    = (TH2F *) outputList->At(index++); 
    
  }
  
  
  if(fCorrelate){
    fhCaloCorrNClusters = (TH2F *) outputList->At(index++);
    fhCaloCorrEClusters = (TH2F *) outputList->At(index++); 
    fhCaloCorrNCells    = (TH2F *) outputList->At(index++); 
    fhCaloCorrECells    = (TH2F *) outputList->At(index++); 
    
    fhCaloV0SCorrNClusters = (TH2F *) outputList->At(index++);
    fhCaloV0SCorrEClusters = (TH2F *) outputList->At(index++); 
    fhCaloV0SCorrNCells    = (TH2F *) outputList->At(index++); 
    fhCaloV0SCorrECells    = (TH2F *) outputList->At(index++); 
    
    fhCaloV0MCorrNClusters = (TH2F *) outputList->At(index++);
    fhCaloV0MCorrEClusters = (TH2F *) outputList->At(index++); 
    fhCaloV0MCorrNCells    = (TH2F *) outputList->At(index++); 
    fhCaloV0MCorrECells    = (TH2F *) outputList->At(index++); 
    
    fhCaloTrackMCorrNClusters = (TH2F *) outputList->At(index++);
    fhCaloTrackMCorrEClusters = (TH2F *) outputList->At(index++); 
    fhCaloTrackMCorrNCells    = (TH2F *) outputList->At(index++); 
    fhCaloTrackMCorrECells    = (TH2F *) outputList->At(index++); 
  }
  
  //Module histograms
  fhEMod                 = new TH1F*[fNModules];
  fhNClustersMod         = new TH1F*[fNModules];
  fhNCellsPerClusterMod  = new TH2F*[fNModules];
  fhNCellsMod            = new TH1F*[fNModules];
  fhGridCellsMod         = new TH2F*[fNModules];
  fhGridCellsEMod        = new TH2F*[fNModules];
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) 
    fhGridCellsTimeMod     = new TH2F*[fNModules];
  fhAmplitudeMod         = new TH1F*[fNModules];
  if(fCalorimeter=="EMCAL")
    fhAmplitudeModFraction = new TH1F*[fNModules*3];
  
  //EMCAL
  fhTimeAmpPerRCU        = new TH2F*[fNModules*fNRCU];
  
  fhIMMod                = new TH2F*[fNModules];
  fhIMCellCutMod         = new TH2F*[fNModules];
  
  for(Int_t imod = 0 ; imod < fNModules; imod++){
    fhEMod[imod]                 = (TH1F *) outputList->At(index++);
    fhNClustersMod[imod]         = (TH1F *) outputList->At(index++); 
    fhNCellsPerClusterMod[imod]  = (TH2F *) outputList->At(index++); 
    fhNCellsMod[imod]            = (TH1F *) outputList->At(index++); 	
    fhGridCellsMod[imod]         = (TH2F *) outputList->At(index++);
    fhGridCellsEMod[imod]        = (TH2F *) outputList->At(index++); 
    if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) 
      fhGridCellsTimeMod[imod]     = (TH2F *) outputList->At(index++); 
    fhAmplitudeMod[imod]         = (TH1F *) outputList->At(index++);
    
    if(fCalorimeter=="EMCAL"){
      for(Int_t ifrac = 0; ifrac < 3; ifrac++){
        fhAmplitudeModFraction[imod*3+ifrac] = (TH1F *) outputList->At(index++); 
      }
    }
    
    for(Int_t ircu = 0; ircu < fNRCU; ircu++){
      fhTimeAmpPerRCU[imod*fNRCU+ircu] = (TH2F *) outputList->At(index++); 
      //fhT0TimeAmpPerRCU[imod*fNRCU+ircu] = (TH2F *) outputList->At(index++); 
      //			for(Int_t imod2 = 0; imod2 < fNModules; imod2++){
      //				for(Int_t ircu2 = 0; ircu2 < fNModules; ircu2++){
      //					fhTimeCorrRCU[imod*fNRCU+ircu+imod2*fNRCU+ircu2]  = (TH2F *) outputList->At(index++);
      //				}
      //			}
    }
    fhIMMod[imod]                = (TH2F *) outputList->At(index++); 
    fhIMCellCutMod[imod]         = (TH2F *) outputList->At(index++); 	
    
  }
  
  if(IsDataMC()){
    fhDeltaE   = (TH1F *) outputList->At(index++); 
    fhDeltaPt  = (TH1F *) outputList->At(index++); 
    fhDeltaPhi = (TH1F *) outputList->At(index++); 
    fhDeltaEta = (TH1F *) outputList->At(index++); 
    
    fhRatioE   = (TH1F *) outputList->At(index++); 
    fhRatioPt  = (TH1F *) outputList->At(index++); 
    fhRatioPhi = (TH1F *) outputList->At(index++); 
    fhRatioEta = (TH1F *) outputList->At(index++); 
    
    fh2E       = (TH2F *) outputList->At(index++); 
    fh2Pt      = (TH2F *) outputList->At(index++); 
    fh2Phi     = (TH2F *) outputList->At(index++); 
    fh2Eta     = (TH2F *) outputList->At(index++); 
    
    fhGamE     = (TH2F *) outputList->At(index++); 
    fhGamPt    = (TH2F *) outputList->At(index++); 
    fhGamPhi   = (TH2F *) outputList->At(index++); 
    fhGamEta   = (TH2F *) outputList->At(index++); 
    
    fhGamDeltaE   = (TH1F *) outputList->At(index++); 
    fhGamDeltaPt  = (TH1F *) outputList->At(index++); 
    fhGamDeltaPhi = (TH1F *) outputList->At(index++); 
    fhGamDeltaEta = (TH1F *) outputList->At(index++); 
    
    fhGamRatioE   = (TH1F *) outputList->At(index++); 
    fhGamRatioPt  = (TH1F *) outputList->At(index++); 
    fhGamRatioPhi = (TH1F *) outputList->At(index++); 
    fhGamRatioEta = (TH1F *) outputList->At(index++); 
    
    fhPi0E     = (TH2F *) outputList->At(index++); 
    fhPi0Pt    = (TH2F *) outputList->At(index++); 
    fhPi0Phi   = (TH2F *) outputList->At(index++); 
    fhPi0Eta   = (TH2F *) outputList->At(index++); 		
    
    fhEleE     = (TH2F *) outputList->At(index++); 
    fhElePt    = (TH2F *) outputList->At(index++); 
    fhElePhi   = (TH2F *) outputList->At(index++); 
    fhEleEta   = (TH2F *) outputList->At(index++); 		
    
    fhNeHadE     = (TH2F *) outputList->At(index++); 
    fhNeHadPt    = (TH2F *) outputList->At(index++); 
    fhNeHadPhi   = (TH2F *) outputList->At(index++); 
    fhNeHadEta   = (TH2F *) outputList->At(index++); 		
    
    fhChHadE     = (TH2F *) outputList->At(index++); 
    fhChHadPt    = (TH2F *) outputList->At(index++); 
    fhChHadPhi   = (TH2F *) outputList->At(index++); 
    fhChHadEta   = (TH2F *) outputList->At(index++); 				
    
    fhGamECharged     = (TH2F *) outputList->At(index++); 
    fhGamPtCharged    = (TH2F *) outputList->At(index++); 
    fhGamPhiCharged   = (TH2F *) outputList->At(index++); 
    fhGamEtaCharged   = (TH2F *) outputList->At(index++); 
    
    fhPi0ECharged     = (TH2F *) outputList->At(index++); 
    fhPi0PtCharged    = (TH2F *) outputList->At(index++); 
    fhPi0PhiCharged   = (TH2F *) outputList->At(index++); 
    fhPi0EtaCharged   = (TH2F *) outputList->At(index++); 		
    
    fhEleECharged     = (TH2F *) outputList->At(index++); 
    fhElePtCharged    = (TH2F *) outputList->At(index++); 
    fhElePhiCharged   = (TH2F *) outputList->At(index++); 
    fhEleEtaCharged   = (TH2F *) outputList->At(index++); 		
    
    fhNeHadECharged     = (TH2F *) outputList->At(index++); 
    fhNeHadPtCharged    = (TH2F *) outputList->At(index++); 
    fhNeHadPhiCharged   = (TH2F *) outputList->At(index++); 
    fhNeHadEtaCharged   = (TH2F *) outputList->At(index++); 		
    
    fhChHadECharged     = (TH2F *) outputList->At(index++); 
    fhChHadPtCharged    = (TH2F *) outputList->At(index++); 
    fhChHadPhiCharged   = (TH2F *) outputList->At(index++); 
    fhChHadEtaCharged   = (TH2F *) outputList->At(index++); 				
		
    //		fhEMVxyz     = (TH3F *) outputList->At(index++); 
    //		fhHaVxyz     = (TH3F *) outputList->At(index++); 
		
    fhEMVxyz     = (TH2F *) outputList->At(index++); 
    fhHaVxyz     = (TH2F *) outputList->At(index++); 
    fhEMR        = (TH2F *) outputList->At(index++); 
    fhHaR        = (TH2F *) outputList->At(index++); 
    
    fhGenGamPt    = (TH1F *) outputList->At(index++); 
    fhGenGamEta   = (TH1F *) outputList->At(index++); 
    fhGenGamPhi   = (TH1F *) outputList->At(index++); 
    
    fhGenPi0Pt    = (TH1F *) outputList->At(index++); 
    fhGenPi0Eta   = (TH1F *) outputList->At(index++); 
    fhGenPi0Phi   = (TH1F *) outputList->At(index++); 
    
    fhGenEtaPt    = (TH1F *) outputList->At(index++); 
    fhGenEtaEta   = (TH1F *) outputList->At(index++); 
    fhGenEtaPhi   = (TH1F *) outputList->At(index++); 
    
    fhGenOmegaPt  = (TH1F *) outputList->At(index++); 
    fhGenOmegaEta = (TH1F *) outputList->At(index++); 
    fhGenOmegaPhi = (TH1F *) outputList->At(index++); 
    
    fhGenElePt    = (TH1F *) outputList->At(index++); 
    fhGenEleEta   = (TH1F *) outputList->At(index++); 
    fhGenElePhi   = (TH1F *) outputList->At(index++); 
    
    fhGenGamAccE   = (TH1F *) outputList->At(index++); 		
    fhGenGamAccPt  = (TH1F *) outputList->At(index++); 
    fhGenGamAccEta = (TH1F *) outputList->At(index++); 
    fhGenGamAccPhi = (TH1F *) outputList->At(index++); 
    
    fhGenPi0AccE   = (TH1F *) outputList->At(index++); 		
    fhGenPi0AccPt  = (TH1F *) outputList->At(index++); 
    fhGenPi0AccEta = (TH1F *) outputList->At(index++); 
    fhGenPi0AccPhi = (TH1F *) outputList->At(index++); 
    
    fhMCEle1pOverE =    (TH2F *) outputList->At(index++);
    fhMCEle1dR =        (TH1F *) outputList->At(index++);
    fhMCEle2MatchdEdx = (TH2F *) outputList->At(index++);
    
    fhMCChHad1pOverE =    (TH2F *) outputList->At(index++);
    fhMCChHad1dR =        (TH1F *) outputList->At(index++);
    fhMCChHad2MatchdEdx = (TH2F *) outputList->At(index++);
    
    fhMCNeutral1pOverE    = (TH2F *) outputList->At(index++);
    fhMCNeutral1dR        = (TH1F *) outputList->At(index++);
    fhMCNeutral2MatchdEdx = (TH2F *) outputList->At(index++);
    
    fhMCEle1pOverER02     =    (TH2F *) outputList->At(index++);
    fhMCChHad1pOverER02   =    (TH2F *) outputList->At(index++);
    fhMCNeutral1pOverER02 =    (TH2F *) outputList->At(index++);
  }
}

//__________________________________________________________________
void  AliAnaCalorimeterQA::Terminate(TList* outputList) 
{
  //Do plots if requested	
  
  if(GetDebug() > 0) printf("AliAnaCalorimeterQA::Terminate() - Make plots for %s? %d\n",fCalorimeter.Data(), MakePlotsOn());
  
  //Do some plots to end
  if(fStyleMacro!="")gROOT->Macro(fStyleMacro); 
  //Recover histograms from output histograms list, needed for distributed analysis.	
  ReadHistograms(outputList);
  
  //printf(" AliAnaCalorimeterQA::Terminate()  *** %s Report:", GetName()) ; 
  //printf(" AliAnaCalorimeterQA::Terminate()        pt         : %5.3f , RMS : %5.3f \n", fhPt->GetMean(),   fhPt->GetRMS() ) ;
  
  const Int_t buffersize = 255;
  char name[buffersize];
  char cname[buffersize];
  
  //In case terminate is executed after the analysis, in a second step, and we want to rebin or to change the range of the histograms for plotting
  Int_t nptbins     = GetHistoPtBins(); 	        Float_t ptmax     = GetHistoPtMax();           Float_t ptmin     = GetHistoPtMin();
  Int_t nphibins    = GetHistoPhiBins();          Float_t phimax    = GetHistoPhiMax();          Float_t phimin    = GetHistoPhiMin();
  Int_t netabins    = GetHistoEtaBins();          Float_t etamax    = GetHistoEtaMax();          Float_t etamin    = GetHistoEtaMin();	
  //	Int_t nmassbins   = GetHistoMassBins();         Float_t massmax   = GetHistoMassMax(); 	       Float_t massmin   = GetHistoMassMin();
  //	Int_t nasymbins   = GetHistoAsymmetryBins();    Float_t asymmax   = GetHistoAsymmetryMax();    Float_t asymmin   = GetHistoAsymmetryMin();
  //	Int_t nPoverEbins = GetHistoPOverEBins();       Float_t pOverEmax = GetHistoPOverEMax();       Float_t pOverEmin = GetHistoPOverEMin();
  //	Int_t ndedxbins   = GetHistodEdxBins();         Float_t dedxmax   = GetHistodEdxMax();         Float_t dedxmin   = GetHistodEdxMin();
  //	Int_t ndRbins     = GetHistodRBins();           Float_t dRmax     = GetHistodRMax();           Float_t dRmin     = GetHistodRMin();
  Int_t ntimebins   = GetHistoTimeBins();         Float_t timemax   = GetHistoTimeMax();         Float_t timemin   = GetHistoTimeMin();       
  Int_t nbins       = GetHistoNClusterCellBins(); Int_t nmax        = GetHistoNClusterCellMax(); Int_t nmin        = GetHistoNClusterCellMin(); 
  //	Int_t nratiobins  = GetHistoRatioBins();        Float_t ratiomax  = GetHistoRatioMax();        Float_t ratiomin  = GetHistoRatioMin();
  //	Int_t nvdistbins  = GetHistoVertexDistBins();   Float_t vdistmax  = GetHistoVertexDistMax();   Float_t vdistmin  = GetHistoVertexDistMin();
  Int_t rbins       = GetHistoRBins();            Float_t rmax        = GetHistoRMax();          Float_t rmin      = GetHistoRMin(); 
  Int_t xbins       = GetHistoXBins();            Float_t xmax        = GetHistoXMax();          Float_t xmin      = GetHistoXMin(); 
  Int_t ybins       = GetHistoYBins();            Float_t ymax        = GetHistoYMax();          Float_t ymin      = GetHistoYMin(); 
  Int_t zbins       = GetHistoZBins();            Float_t zmax        = GetHistoZMax();          Float_t zmin      = GetHistoZMin(); 
  
  //Color code for the different modules
  Int_t modColorIndex[]={2,4,6,8};
  
  //--------------------------------------------------
  // Cluster energy distributions, module dependence
  //--------------------------------------------------
  snprintf(cname,buffersize,"QA_%s_ClusterEnergy",fCalorimeter.Data());
  TCanvas  * c = new TCanvas(cname, "Energy distributions", 800, 400) ;
  c->Divide(2, 1);
  Int_t rbE = GetNewRebinForRePlotting((TH1D*)fhE, ptmin, ptmax,nptbins) ;
  //printf("new E rb %d\n",rbE);
  fhE->Rebin(rbE);
  fhE->SetAxisRange(ptmin,ptmax,"X");
  c->cd(1) ; 
  if(fhE->GetEntries() > 0) gPad->SetLogy();
  TLegend pLegendE(0.7,0.6,0.9,0.8);
  pLegendE.SetTextSize(0.03);
  pLegendE.AddEntry(fhE,"all modules","L");
  pLegendE.SetFillColor(10);
  pLegendE.SetBorderSize(1);
  
  fhE->SetMinimum(1);	
  fhE->SetLineColor(1);
  fhE->Draw("HE");
  for(Int_t imod = 0; imod < fNModules; imod++){
    fhEMod[imod]->Rebin(rbE);
    fhEMod[imod]->SetLineColor(modColorIndex[imod]);
    fhEMod[imod]->Draw("HE same");
    pLegendE.AddEntry(fhEMod[imod],Form("module %d",imod),"L");
  }
  pLegendE.Draw();
  
  //Ratio of modules
  c->cd(2) ; 
  TLegend pLegendER(0.55,0.8,0.9,0.9);
  pLegendER.SetTextSize(0.03);
  pLegendER.SetFillColor(10);
  pLegendER.SetBorderSize(1);
  
  for(Int_t imod = 1; imod < fNModules; imod++){
    TH1D * htmp = (TH1D*)fhEMod[imod]->Clone(Form("hERat%d",imod));
    htmp->Divide(fhEMod[0]);
    htmp->SetLineColor(modColorIndex[imod]);
    if(imod==1){
      htmp->SetTitle("Ratio module X / module 0");
      htmp->SetAxisRange(ptmin,ptmax,"X");
      htmp->SetMaximum(5);
      htmp->SetMinimum(0);
      htmp->SetAxisRange(ptmin,ptmax,"X");
      htmp->Draw("HE");
    }
    else 
      htmp->Draw("same HE");
    
    pLegendER.AddEntry(fhEMod[imod],Form("module %d / module 0",imod),"L");
  }
  pLegendER.Draw();
  
  snprintf(name,buffersize,"QA_%s_ClusterEnergy.eps",fCalorimeter.Data());
  c->Print(name); printf("Plot: %s\n",name);
  
  //--------------------------------------------------
  // Cell energy distributions, module dependence
  //--------------------------------------------------
  snprintf(cname,buffersize,"%s_QA_CellEnergy",fCalorimeter.Data());
  TCanvas  * ca = new TCanvas(cname, "Cell Energy distributions", 800, 400) ;
  ca->Divide(2, 1);
  
  Int_t rbAmp = GetNewRebinForRePlotting((TH1D*)fhAmplitude, ptmin, ptmax,nptbins*2) ;
  //printf("new Amp rb %d\n",rbAmp);
  fhAmplitude->Rebin(rbAmp);
  fhAmplitude->SetAxisRange(ptmin,ptmax,"X");
  
  ca->cd(1) ; 
  if(fhAmplitude->GetEntries() > 0) gPad->SetLogy();
  TLegend pLegendA(0.7,0.6,0.9,0.8);
  pLegendA.SetTextSize(0.03);
  pLegendA.AddEntry(fhE,"all modules","L");
  pLegendA.SetFillColor(10);
  pLegendA.SetBorderSize(1);
  fhAmplitude->SetMinimum(0.1);
  fhAmplitude->SetLineColor(1);
  fhAmplitude->Draw("HE");
  
  for(Int_t imod = 0; imod < fNModules; imod++){
    fhAmplitudeMod[imod]->Rebin(rbAmp);
    fhAmplitudeMod[imod]->SetLineColor(modColorIndex[imod]);
    fhAmplitudeMod[imod]->Draw("HE same");
    pLegendA.AddEntry(fhAmplitudeMod[imod],Form("module %d",imod),"L");
  }
  pLegendA.Draw();
  
  
  ca->cd(2) ; 
  TLegend pLegendAR(0.55,0.8,0.9,0.9);
  pLegendAR.SetTextSize(0.03);
  pLegendAR.SetFillColor(10);
  pLegendAR.SetBorderSize(1);
  
  for(Int_t imod = 1; imod < fNModules; imod++){
    TH1D * htmp = (TH1D*)fhAmplitudeMod[imod]->Clone(Form("hAmpRat%d",imod));
    htmp->Divide(fhAmplitudeMod[0]);
    htmp->SetLineColor(modColorIndex[imod]);
    if(imod==1){
      htmp->SetTitle("Ratio cells energy in  module X / module 0");
      htmp->SetAxisRange(ptmin,ptmax,"X");
      htmp->SetMaximum(5);
      htmp->SetMinimum(0);
      htmp->Draw("HE");
    }
    else 
      htmp->Draw("same HE");
    pLegendAR.AddEntry(fhAmplitudeMod[imod],Form("module %d",imod),"L");
  }
  
  pLegendAR.Draw();
  snprintf(name,buffersize,"QA_%s_CellEnergy.eps",fCalorimeter.Data());
  ca->Print(name); printf("Plot: %s\n",name);	
  
  //----------------------------------------------------------
  // Cell energy distributions, FRACTION of module dependence
  // See Super Module calibration difference
  //---------------------------------------------------------	
  if(fCalorimeter=="EMCAL"){
    //Close To Eta 0 
    snprintf(cname,buffersize,"%s_QA_SMThirds",fCalorimeter.Data());
    TCanvas  * cfrac = new TCanvas(cname, "SM Thirds ratios", 800, 1200) ;
    cfrac->Divide(2, 3);
    cfrac->cd(1) ; 
    if(fhAmplitude->GetEntries() > 0) 
      gPad->SetLogy();
    TLegend pLegend1(0.6,0.6,0.9,0.8);
    pLegend1.SetTextSize(0.03);
    pLegend1.SetFillColor(10);
    pLegend1.SetBorderSize(1);
    pLegend1.SetHeader("Third close to Eta=0");
    fhAmplitudeModFraction[0]->SetTitle("Third close to Eta=0");
    fhAmplitudeModFraction[0]->SetAxisRange(ptmin,ptmax,"X");
    fhAmplitudeModFraction[0]->Draw("axis");
    TH1D * hAverageThird1 = (TH1D *)fhAmplitudeModFraction[3*0+2]->Clone("AverageThird1");
    for(Int_t imod = 0; imod < fNModules; imod++){
      Int_t ifrac = 0;
      if(imod%2==0) ifrac = 2;
      if(imod > 0) hAverageThird1->Add( fhAmplitudeModFraction[3*imod+ifrac]);
      fhAmplitudeModFraction[3*imod+ifrac]->SetLineColor(modColorIndex[imod]);
      fhAmplitudeModFraction[3*imod+ifrac]->Draw("HE same");
      pLegend1.AddEntry(fhAmplitudeModFraction[3*imod+ifrac],Form("super module %d",imod),"L");
    }
    hAverageThird1 ->Scale(1./fNModules);
    pLegend1.Draw();
    //Ratio
    cfrac->cd(2) ; 
    for(Int_t imod = 0; imod < fNModules; imod++){
      Int_t ifrac = 0;
      if(imod%2==0) ifrac = 2;
      TH1D * htmp =  (TH1D*)fhAmplitudeModFraction[3*imod+ifrac]->Clone(Form("ThirdFractionAverage_%d_%d",imod,ifrac));
      htmp->Divide(hAverageThird1);
      if(imod ==0) {
        htmp ->SetTitle("Close to eta = 0");
        htmp ->SetMaximum(5);
        htmp ->SetMinimum(0);
        htmp ->SetAxisRange(ptmin,ptmax,"X");
        htmp ->SetYTitle("ratio third to average");
        htmp -> Draw("HE");
      }
      else htmp -> Draw("same HE");
    }
    //pLegend1.Draw();
    
    //Middle Eta
    cfrac->cd(3) ; 
    if(fhAmplitude->GetEntries() > 0) 
      gPad->SetLogy();
    TLegend pLegend2(0.6,0.6,0.9,0.8);
    pLegend2.SetTextSize(0.03);
    pLegend2.SetFillColor(10);
    pLegend2.SetBorderSize(1);
    pLegend2.SetHeader("Middle Third");
    
    fhAmplitudeModFraction[0]->SetTitle("Middle Third");
    fhAmplitudeModFraction[0]->SetAxisRange(ptmin,ptmax,"X");
    fhAmplitudeModFraction[0]->Draw("axis");
    
    TH1D * hAverageThird2 = (TH1D *)fhAmplitudeModFraction[3*0+1]->Clone("AverageThird2");
    for(Int_t imod = 0; imod < fNModules; imod++){
      Int_t ifrac = 1;
      if(imod > 0) hAverageThird2->Add( fhAmplitudeModFraction[3*imod+ifrac]);
      fhAmplitudeModFraction[3*imod+ifrac]->SetLineColor(modColorIndex[imod]);
      fhAmplitudeModFraction[3*imod+ifrac]->Draw("HE same");
      pLegend2.AddEntry(fhAmplitudeModFraction[3*imod+ifrac],Form("super module %d",imod),"L");
    }
    hAverageThird2->Scale(1./fNModules);
    pLegend2.Draw();
    
    //Ratio
    cfrac->cd(4) ; 
    
    for(Int_t imod = 0; imod < fNModules; imod++){
      Int_t ifrac = 1;
      TH1D * htmp =  (TH1D*)fhAmplitudeModFraction[3*imod+ifrac]->Clone(Form("ThirdFractionAverage_%d_%d",imod,ifrac));
      htmp->Divide(hAverageThird2);
      if(imod ==0) {
        htmp ->SetTitle("Middle");
        htmp ->SetMaximum(5);
        htmp ->SetMinimum(0);
        htmp ->SetAxisRange(ptmin,ptmax,"X");
        htmp ->SetYTitle("ratio third to average");
        htmp -> Draw("HE");
      }
      else htmp -> Draw("same HE");
    }
    //pLegend2.Draw();
    
    //Close To Eta 0.7 
    cfrac->cd(5) ; 
    if(fhAmplitude->GetEntries() > 0) 
      gPad->SetLogy();
    TLegend pLegend3(0.6,0.6,0.9,0.8);
    pLegend3.SetTextSize(0.03);
    pLegend3.SetFillColor(10);
    pLegend3.SetBorderSize(1);
    pLegend3.SetHeader("Third close to Eta=0.7");
    
    fhAmplitudeModFraction[0]->SetTitle("Third close to Eta=0.7");
    fhAmplitudeModFraction[0]->SetAxisRange(ptmin,ptmax,"X");
    fhAmplitudeModFraction[0]->Draw("axis");
    
    TH1D * hAverageThird3 = (TH1D *)fhAmplitudeModFraction[3*0+0]->Clone("AverageThird3");
    for(Int_t imod = 0; imod < 4; imod++){
      Int_t ifrac = 2;
      if(imod%2==0) ifrac = 0;
      if(imod > 0) hAverageThird3->Add( fhAmplitudeModFraction[3*imod+ifrac]);
      fhAmplitudeModFraction[3*imod+ifrac]->SetLineColor(modColorIndex[imod]);
      fhAmplitudeModFraction[3*imod+ifrac]->Draw("HE same");
      pLegend3.AddEntry(fhAmplitudeModFraction[3*imod+ifrac],Form("super module %d",imod),"L");
    }
    hAverageThird3 ->Scale(1./fNModules);
    pLegend3.Draw();
    
    cfrac->cd(6) ; 
    
    for(Int_t imod = 0; imod < fNModules; imod++){
      Int_t ifrac = 2;
      if(imod%2==0) ifrac = 0;
      TH1D * htmp =  (TH1D*)fhAmplitudeModFraction[3*imod+ifrac]->Clone(Form("ThirdFractionAverage_%d_%d",imod,ifrac));
      htmp->Divide(hAverageThird3);
      if(imod ==0) {
        htmp ->SetTitle("Close to eta = 0.7");
        htmp ->SetMaximum(5);
        htmp ->SetMinimum(0);
        htmp ->SetAxisRange(ptmin,ptmax,"X");
        htmp ->SetYTitle("ratio third to average");
        htmp ->Draw("HE");
      }
      else htmp ->Draw("same HE");
    }
    //pLegend3.Draw();
    
    snprintf(name,buffersize,"QA_%s_CellEnergyModuleFraction.eps",fCalorimeter.Data());
    cfrac->Print(name); printf("Create plot %s\n",name);
  }//EMCAL	
  
  
  //----------------------------------------------------------
  // Cluster eta and phi distributions, energy cut dependence
  //---------------------------------------------------------	
  
  snprintf(cname,buffersize,"%s_QA_EtaPhiCluster",fCalorimeter.Data());
  TCanvas  * cetaphic = new TCanvas(cname, "Eta-Phi Reconstructed distributions", 1200, 400) ;
  cetaphic->Divide(3, 1);
  Int_t binmin = 0;
  Int_t rbPhi  = 1;
  Int_t rbEta  = 1;
  Int_t ncuts  = 7;
  Float_t ecut[]     = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3};
  Int_t   ecutcolor[]= {2, 4, 6, 7, 8, 9, 12};
  TH1D * hE = fhEtaPhiE->ProjectionZ();
  
  //PHI
  cetaphic->cd(1) ; 
  gPad->SetLogy();
  gPad->SetGridy();
  
  TLegend pLegendPhiCl(0.83,0.6,0.95,0.93);
  pLegendPhiCl.SetTextSize(0.03);
  pLegendPhiCl.SetFillColor(10);
  pLegendPhiCl.SetBorderSize(1);
  
  TH1D * htmp = fhEtaPhiE->ProjectionY("hphi_cluster_nocut",0,-1,0,-1);
  if(htmp){
    htmp->SetMinimum(1);
    rbPhi =  GetNewRebinForRePlotting(htmp, phimin, phimax,nphibins) ;
    //printf("new Phi rb %d\n",rbPhi);
    htmp->Rebin(rbPhi);
    htmp->SetTitle("#phi of clusters for energy in cluster > threshold");
    htmp->SetAxisRange(phimin,phimax,"X");
    htmp->Draw("HE");
    pLegendPhiCl.AddEntry(htmp,"No cut","L");
    
    for (Int_t i = 0; i < ncuts; i++) {
      binmin =  hE->FindBin(ecut[i]);
      //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
      htmp = fhEtaPhiE->ProjectionY(Form("hphi_cluster_cut%d",i),0,-1,binmin,-1);
      htmp->SetLineColor(ecutcolor[i]);
      htmp->Rebin(rbPhi);
      htmp->Draw("same HE");
      pLegendPhiCl.AddEntry(htmp,Form("E>%1.1f",ecut[i]),"L");
      
    }
  }
  pLegendPhiCl.Draw();
  
  //ETA
  cetaphic->cd(2) ; 
  gPad->SetLogy();
  gPad->SetGridy();
  
  delete htmp; 
  htmp = fhEtaPhiE->ProjectionX("heta_cluster_nocut",0,-1,0,-1);
  if(htmp){
    rbEta =  GetNewRebinForRePlotting(htmp,etamin, etamax,netabins) ;
    //printf("new Eta rb %d\n",rbEta);
    htmp->Rebin(rbEta);
    htmp->SetMinimum(1);
    htmp ->SetLineColor(1);
    htmp->SetTitle("#eta of clusters for energy in cluster > threshold");
    htmp->SetAxisRange(etamin,etamax,"X");
    htmp->Draw("HE");
    
    for (Int_t i = 0; i < ncuts; i++) {
      binmin =  hE->FindBin(ecut[i]);
      //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
      htmp = fhEtaPhiE->ProjectionX(Form("heta_cluster_cut%d",i),0,-1,binmin,-1);
      htmp->SetLineColor(ecutcolor[i]);
      htmp->Rebin(rbEta);
      htmp->Draw("same HE");	
    }
  }
  //ETA vs PHI	
  cetaphic->cd(3) ;
  TH2D* hEtaPhiCl = (TH2D*) fhEtaPhiE->Project3D("xy");
  hEtaPhiCl->SetAxisRange(etamin,etamax,"X");
  hEtaPhiCl->SetAxisRange(phimin,phimax,"Y");
  hEtaPhiCl->Draw("colz");
  
  snprintf(name,buffersize,"QA_%s_ClusterEtaPhi.eps",fCalorimeter.Data());
  cetaphic->Print(name); printf("Create plot %s\n",name);
  
  //----------------------------------------------------------
  // Cell eta and phi distributions, energy cut dependence
  //---------------------------------------------------------	
	
  snprintf(cname,buffersize,"%s_QA_EtaPhiCell",fCalorimeter.Data());
  TCanvas  * cetaphicell = new TCanvas(cname, "Eta-Phi Cells distributions", 1200, 400) ;
  cetaphicell->Divide(3, 1);
  
  //PHI
  cetaphicell->cd(1) ; 
  gPad->SetLogy();
  gPad->SetGridy();
  
  TLegend pLegendPhiCell(0.83,0.6,0.95,0.93);
  pLegendPhiCell.SetTextSize(0.03);
  pLegendPhiCell.SetFillColor(10);
  pLegendPhiCell.SetBorderSize(1);
  
  delete htmp; 
  htmp = fhEtaPhiAmp->ProjectionY("hphi_cell_nocut",0,-1,0,-1);
  if(htmp){
    htmp->SetMinimum(1);
    htmp->Rebin(rbPhi);
    htmp->SetTitle("#phi of cells for cell energy > threshold");
    htmp->SetAxisRange(phimin,phimax,"X");
    htmp->Draw("HE");
    pLegendPhiCell.AddEntry(htmp,"No cut","L");
    
    for (Int_t i = 0; i < ncuts; i++) {
      binmin =  hE->FindBin(ecut[i]);
      //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
      htmp = fhEtaPhiAmp->ProjectionY(Form("hphi_cell_cut%d",i),0,-1,binmin,-1);
      htmp->SetLineColor(ecutcolor[i]);
      htmp->Rebin(rbPhi);
      htmp->Draw("same HE");
      pLegendPhiCl.AddEntry(htmp,Form("E>%1.1f",ecut[i]),"L");
      
    }
  }
  pLegendPhiCell.Draw();
  
  //ETA
  cetaphicell->cd(2) ; 
  gPad->SetLogy();
  gPad->SetGridy();
  
  delete htmp; 
  htmp = fhEtaPhiAmp->ProjectionX("heta_cell_nocut",0,-1,0,-1);
  if(htmp){
    htmp ->SetLineColor(1);
    htmp->Rebin(rbEta);
    htmp->SetMinimum(1);
    htmp->SetTitle("#eta of cells for cell energy > threshold");
    htmp->SetAxisRange(etamin,etamax,"X");
    htmp->Draw("HE");
    
    for (Int_t i = 0; i < ncuts; i++) {
      binmin =  hE->FindBin(ecut[i]);
      //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
      htmp = fhEtaPhiAmp->ProjectionX(Form("heta_cell_cut%d",i),0,-1,binmin,-1);
      htmp->SetLineColor(ecutcolor[i]);
      htmp->Rebin(rbEta);
      htmp->Draw("same HE");
      
    }
  }
  //ETA vs PHI	
  cetaphicell->cd(3) ;
  TH2D* hEtaPhiCell = (TH2D*) fhEtaPhiAmp->Project3D("xy");
  hEtaPhiCell->SetAxisRange(etamin,etamax,"X");
  hEtaPhiCell->SetAxisRange(phimin,phimax,"Y");
  hEtaPhiCell->Draw("colz");
  
  snprintf(name,buffersize,"QA_%s_CellEtaPhi.eps",fCalorimeter.Data());
  cetaphicell->Print(name); printf("Create plot %s\n",name);
  
  
  ////////////////////////////////////////        
  ///////// Global Positions /////////////       
  ////////////////////////////////////////       
	
  //CLUSTERS
  Int_t rbX = 1;
  Int_t rbY = 1;
  Int_t rbZ = 1;
  if(fFillAllPosHisto)
  {
    snprintf(cname,buffersize,"%s_QA_ClusterXY",fCalorimeter.Data());
    TCanvas  * cxyz = new TCanvas(cname, "Cluster XY distributions", 1200, 400) ;
    cxyz->Divide(3, 1);
    
    cxyz->cd(1) ; 
    TH2D * hXY = (TH2D*) fhXYZ->Project3D("yx" );
    hXY->SetTitle("Cluster X vs Y");
    hXY->GetYaxis()->SetTitleOffset(1.6);
    hXY->Draw("colz");
    cxyz->cd(2) ; 
    TH2D * hYZ = (TH2D*) fhXYZ->Project3D("yz" );
    hYZ->SetTitle("Cluster Z vs Y");
    hYZ->GetYaxis()->SetTitleOffset(1.6);
    hYZ->Draw("colz");	
    cxyz->cd(3) ; 
    TH2D * hXZ = (TH2D*) fhXYZ->Project3D("zx" );
    hXZ->SetTitle("Cluster X vs Z");
    hXZ->GetYaxis()->SetTitleOffset(1.6);
    hXZ->Draw("colz");
    
    snprintf(name,buffersize,"QA_%s_ClusterXY_YZ_XZ.eps",fCalorimeter.Data());
    cxyz->Print(name); printf("Create plot %s\n",name);
    
    snprintf(cname,buffersize,"QA_%s_ClusterX",fCalorimeter.Data());
    TCanvas  * cx = new TCanvas(cname, "Cluster X distributions", 1200, 400) ;
    cx->Divide(3, 1);
    
    cx->cd(1) ; 
    TH1D * hX = (TH1D*) fhXYZ->Project3D("xe" );
    //gPad->SetLogy();
    gPad->SetGridy();
    hX->SetTitle("Cluster X ");
    hX->Draw("HE");
    rbX =  GetNewRebinForRePlotting(hX, xmin, xmax,xbins) ;
    //printf("new X rb %d\n",rbX);
    hX->Rebin(rbX);
    hX->SetMinimum(hX->GetMaximum()/2);
    hX->SetAxisRange(xmin,xmax);
    
    cx->cd(2) ; 
    TH1D * hY = (TH1D*) fhXYZ->Project3D("ye" );
    //gPad->SetLogy();
    hY->SetTitle("Cluster Y ");
    rbY =  GetNewRebinForRePlotting(hY, ymin, ymax, ybins) ;
    //printf("new Y rb %d\n",rbY);
    hY->Rebin(rbY);
    hY->SetMinimum(1);
    hY->SetAxisRange(ymin,ymax);
    hY->Draw("HE");	
    
    cx->cd(3) ; 
    TH1D * hZ = (TH1D*) fhXYZ->Project3D("ze" );
    //gPad->SetLogy();
    gPad->SetGridy();
    rbZ =  GetNewRebinForRePlotting(hZ,zmin, zmax,zbins) ;
    //printf("new Z rb %d\n",rbZ);
    hZ->Rebin(rbZ);	
    hZ->SetMinimum(hZ->GetMaximum()/2);
    hZ->SetAxisRange(zmin,zmax);
    hZ->Draw("HE");
    
    snprintf(name,buffersize,"QA_%s_ClusterX_Y_Z.eps",fCalorimeter.Data());
    cx->Print(name); printf("Create plot %s\n",name);
  }
  //CELLS
  if(fFillAllPosHisto)
  { 
    snprintf(cname,buffersize,"%s_QA_CellXY",fCalorimeter.Data());
    TCanvas  * cellxyz = new TCanvas(cname, "Cell XY distributions", 1200, 400) ;
    cellxyz->Divide(3, 1);
    
    cellxyz->cd(1) ; 
    TH2D * hXYCell = (TH2D*) fhXYZCell->Project3D("yx" );
    hXYCell->SetTitle("Cell X vs Y");
    hXYCell->GetYaxis()->SetTitleOffset(1.6);
    hXYCell->Draw("colz");
    cellxyz->cd(2) ; 
    TH2D * hYZCell = (TH2D*) fhXYZCell->Project3D("yz" );
    hYZCell->SetTitle("Cell Z vs Y");
    hYZCell->GetYaxis()->SetTitleOffset(1.6);
    hYZCell->Draw("colz");	
    cellxyz->cd(3) ; 
    TH2D * hXZCell = (TH2D*) fhXYZCell->Project3D("zx" );
    hXZCell->SetTitle("Cell X vs Z");
    hXZCell->GetYaxis()->SetTitleOffset(1.6);
    hXZCell->Draw("colz");
    
    snprintf(name,buffersize,"QA_%s_CellXY_YZ_XZ.eps",fCalorimeter.Data());
    cellxyz->Print(name); printf("Create plot %s\n",name);
    
    
    snprintf(cname,buffersize,"%s_QA_CellX",fCalorimeter.Data());
    TCanvas  * cellx = new TCanvas(cname, "Cell X distributions", 1200, 400) ;
    cellx->Divide(3, 1);
    
    cellx->cd(1) ; 
    TH1D * hXCell = (TH1D*) fhXYZCell->Project3D("xe" );
    //gPad->SetLogy();
    gPad->SetGridy();
    hXCell->SetTitle("Cell X ");
    hXCell->Rebin(rbX);
    hXCell->SetMinimum(hXCell->GetMaximum()/2);
    hXCell->SetAxisRange(xmin,xmax);
    hXCell->Draw("HE");
    
    cellx->cd(2) ; 
    TH1D * hYCell = (TH1D*) fhXYZCell->Project3D("ye" );
    //gPad->SetLogy();
    hYCell->SetTitle("Cell Y ");
    hYCell->Rebin(rbY);
    hYCell->SetAxisRange(ymin,ymax);
    hYCell->SetMinimum(1);
    hYCell->Draw("HE");	
    
    cellx->cd(3) ; 
    TH1D * hZCell = (TH1D*) fhXYZCell->Project3D("ze" );
    //gPad->SetLogy();
    gPad->SetGridy();
    hZCell->SetAxisRange(zmin,zmax);
    hZCell->SetTitle("Cell Z ");
    hZCell->Rebin(rbZ);
    hZCell->SetMinimum(hZCell->GetMaximum()/2);
    hZCell->Draw("HE");
    
    snprintf(name,buffersize,"QA_%s_CellX_Y_Z.eps",fCalorimeter.Data());
    cellx->Print(name); printf("Create plot %s\n",name);
    
    
    //----------------------------------------------------------
    // Cluster X, Y, Z, R, energy cut dependence
    //---------------------------------------------------------	
    
    snprintf(cname,buffersize,"%s_QA_ClusterX_Y_Z_R_ECut",fCalorimeter.Data());
    TCanvas  * cxe = new TCanvas(cname, "Cluster X Y Z R, E cut", 800, 800) ;
    cxe->Divide(2, 2);		
    //R
    cxe->cd(1) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    TLegend pLegendXCl(0.83,0.6,0.95,0.93);
    pLegendXCl.SetTextSize(0.03);
    pLegendXCl.SetFillColor(10);
    pLegendXCl.SetBorderSize(1);
    
    delete htmp; 
    htmp = fhRE->ProjectionX("hre_cluster_nocut",0,-1);
    Int_t rbR=1;
    if(htmp){
      htmp->SetMinimum(1);
      rbR =  GetNewRebinForRePlotting(htmp, rmin, rmax,rbins) ;
      //printf("new R rb %d\n",rbR);
      htmp->Rebin(rbR);
      htmp->SetTitle("r of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(rmin,rmax,"X");
      htmp->Draw("HE");
      pLegendXCl.AddEntry(htmp,"No cut","L");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        htmp = fhRE->ProjectionX(Form("hre_cluster_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbR);
        htmp->Draw("same HE");
        pLegendXCl.AddEntry(htmp,Form("E>%1.1f",ecut[i]),"L");
      }
    }
    pLegendXCl.Draw();
    
    //X
    cxe->cd(2) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhXE->ProjectionX("hxe_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbX);
      htmp->SetTitle("x of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(xmin,xmax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        htmp = fhXE->ProjectionX(Form("hxe_cluster_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbX);
        htmp->Draw("same HE");
      }
    }
    //Y
    cxe->cd(3) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    delete htmp; 
    htmp = fhYE->ProjectionX("hye_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbY);
      htmp->SetTitle("y of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(ymin,ymax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        htmp = fhYE->ProjectionX(Form("hye_cluster_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbY);
        htmp->Draw("same HE");
      }
    }
    //Z
    cxe->cd(4) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    delete htmp; 
    htmp = fhZE->ProjectionX("hze_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbZ);
      htmp->SetTitle("z of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(zmin,zmax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        htmp = fhZE->ProjectionX(Form("hze_cluster_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbZ);
        htmp->Draw("same HE"); 
      }
    }
    
    snprintf(name,buffersize,"QA_%s_ClusterX_Y_Z_R_ECut.eps",fCalorimeter.Data());
    cxe->Print(name); printf("Create plot %s\n",name);
    
    
    //----------------------------------------------------------
    // Cluster X, Y, Z, R, NCells in cluster dependence
    //---------------------------------------------------------	
    Int_t ncellcut[]={2, 3, 4};
    Int_t ncellcuts = 3;
    snprintf(cname,buffersize,"%s_QA_ClusterX_Y_Z_R_NCellsCut",fCalorimeter.Data());
    TCanvas  * cxn = new TCanvas(cname, "Cluster X Y Z R, NCells cut", 800, 800) ;
    cxn->Divide(2, 2);		
    //R
    cxn->cd(1) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    TLegend pLegendXClN(0.83,0.6,0.95,0.93);
    pLegendXClN.SetTextSize(0.03);
    pLegendXClN.SetFillColor(10);
    pLegendXClN.SetBorderSize(1);
    
    delete htmp; 
    htmp = fhRNCells->ProjectionX("hrn_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbR);
      htmp->SetTitle("r of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(rmin,rmax,"X");
      htmp->Draw("HE");
      pLegendXClN.AddEntry(htmp,"No cut","L");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        if(i < ncellcuts-1) htmp = fhRNCells->ProjectionX(Form("hrn_cluster_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhRNCells->ProjectionX(Form("hrn_cluster_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbR);
        htmp->Draw("same HE");
        if(i < ncellcuts-1) pLegendXClN.AddEntry(htmp,Form("n = %1.1d",ncellcut[i]-1),"L");
        else pLegendXClN.AddEntry(htmp,Form("n >= %1.1d",ncellcut[i]-1),"L");
        
      }
    }
    pLegendXClN.Draw();
    
    //X
    cxn->cd(2) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    delete htmp; 
    htmp = fhXNCells->ProjectionX("hxn_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbX);
      htmp->SetTitle("x of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(xmin,xmax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        if(i < ncellcuts-1)htmp = fhXNCells->ProjectionX(Form("hxn_cluster_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhXNCells->ProjectionX(Form("hxn_cluster_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbX);
        htmp->Draw("same HE");   
      }
    }
    //Y
    cxn->cd(3) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhYNCells->ProjectionX("hyn_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbY);
      htmp->SetTitle("y of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(ymin,ymax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        if(i < ncellcuts-1) htmp = fhYNCells->ProjectionX(Form("hyn_cluster_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhYNCells->ProjectionX(Form("hyn_cluster_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbY);
        htmp->Draw("same HE");  
      }
    }
    //Z
    cxn->cd(4) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    delete htmp; 
    htmp = fhZNCells->ProjectionX("hzn_cluster_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbZ);
      htmp->SetTitle("z of clusters for energy in cluster > threshold");
      htmp->SetAxisRange(zmin,zmax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        if(i < ncellcuts-1)htmp = fhZNCells->ProjectionX(Form("hzn_cluster_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhZNCells->ProjectionX(Form("hzn_cluster_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbZ);
        htmp->Draw("same HE");    
      }
    }
    
    snprintf(name,buffersize,"QA_%s_ClusterX_Y_Z_R_NCellsCut.eps",fCalorimeter.Data());
    cxn->Print(name); printf("Create plot %s\n",name);
    
    
    //----------------------------------------------------------
    // Cell X, Y, Z, R, energy cut dependence
    //---------------------------------------------------------	
    
    snprintf(cname,buffersize,"%s_QA_CellX_Y_Z_R_ECut",fCalorimeter.Data());
    TCanvas  * cxecell = new TCanvas(cname, "Cell X Y Z R, E cut", 800, 800) ;
    cxecell->Divide(2, 2);		
    //R
    cxecell->cd(1) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    TLegend pLegendXCell(0.83,0.6,0.95,0.93);
    pLegendXCell.SetTextSize(0.03);
    pLegendXCell.SetFillColor(10);
    pLegendXCell.SetBorderSize(1);
    
    delete htmp; 
    htmp = fhRCellE->ProjectionX("hre_cell_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbR);
      htmp->SetTitle("r of cells for energy in cluster > threshold");
      htmp->SetAxisRange(rmin,rmax,"X");
      htmp->Draw("HE");
      pLegendXCell.AddEntry(htmp,"No cut","L");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        htmp = fhRCellE->ProjectionX(Form("hre_celr_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbR);
        htmp->Draw("same HE");
        pLegendXCell.AddEntry(htmp,Form("E>%1.1f",ecut[i]),"L"); 
      }
    }
    pLegendXCell.Draw();
    
    //X
    cxecell->cd(2) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    delete htmp; 
    htmp = fhXCellE->ProjectionX("hxe_cells_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbX);
      htmp->SetTitle("x of cells for energy in cluster > threshold");
      htmp->SetAxisRange(xmin,xmax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        htmp = fhXCellE->ProjectionX(Form("hxe_cells_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbX);
        htmp->Draw("same HE");
      }
    }
    //Y
    cxecell->cd(3) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhYCellE->ProjectionX("hye_cells_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbY);
      htmp->SetTitle("y of cells for energy in cluster > threshold");
      htmp->SetAxisRange(ymin,ymax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        delete htmp; 
        htmp = fhYCellE->ProjectionX(Form("hye_cells_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbY);
        htmp->Draw("same HE");
      }
    }
    //Z
    cxecell->cd(4) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhZCellE->ProjectionX("hze_cells_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbZ);
      htmp->SetTitle("z of cells for energy in cluster > threshold");
      htmp->SetAxisRange(zmin,zmax,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        delete htmp; 
        htmp = fhZCellE->ProjectionX(Form("hze_cells_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbZ);
        htmp->Draw("same HE"); 
      }
    }
    snprintf(name,buffersize,"QA_%s_CellX_Y_Z_R_ECut.eps",fCalorimeter.Data());
    cxecell->Print(name); printf("Create plot %s\n",name);
    
    
    //----------------------------------------------------------
    // Cluster-Cell X, Y, Z, R, cluster energy cut dependence
    //---------------------------------------------------------	
    Int_t rbDR= 1;//rbR;
    Int_t rbDX= 1;//rbX;
    Int_t rbDY= 1;//rbY;
    Int_t rbDZ= 1;//rbZ;
    
    snprintf(cname,buffersize,"%s_QA_DeltaClusterCellX_Y_Z_R_ECut",fCalorimeter.Data());
    TCanvas  * cxde = new TCanvas(cname, "Cluster-Cell X, Y, Z, R, E cut", 800, 800) ;
    cxde->Divide(2, 2);		
    //R
    cxde->cd(1) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    TLegend pLegendXClD(0.83,0.6,0.95,0.93);
    pLegendXClD.SetTextSize(0.03);
    pLegendXClD.SetFillColor(10);
    pLegendXClD.SetBorderSize(1);
    
    delete htmp; 
    htmp = fhDeltaCellClusterRE->ProjectionX("hrde_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDR);
      htmp->SetTitle("r clusters - r cells for energy in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      pLegendXCl.AddEntry(htmp,"No cut","L");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        delete htmp; 
        htmp = fhDeltaCellClusterRE->ProjectionX(Form("hrde_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDR);
        htmp->Draw("same HE");
        pLegendXClD.AddEntry(htmp,Form("E>%1.1f",ecut[i]),"L");
      }
    }
    pLegendXClD.Draw();
    
    //X
    cxde->cd(2) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhDeltaCellClusterXE->ProjectionX("hxde_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDX);
      htmp->SetTitle("x clusters -x cells for energy in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        delete htmp; 
        htmp = fhDeltaCellClusterXE->ProjectionX(Form("hxde_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDX);
        htmp->Draw("same HE");
        
      }
    }
    //Y
    cxde->cd(3) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhDeltaCellClusterYE->ProjectionX("hyde_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDY);
      htmp->SetTitle("y clusters - ycells for energy in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        delete htmp; 
        htmp = fhDeltaCellClusterYE->ProjectionX(Form("hyde_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDY);
        htmp->Draw("same HE");
        
      }
    }
    //Z
    cxde->cd(4) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    delete htmp; 
    htmp = fhDeltaCellClusterZE->ProjectionX("hzde_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbZ);
      htmp->SetTitle("z clusters - z cells for energy in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncuts; i++) {
        binmin =  hE->FindBin(ecut[i]);
        //printf(" bins %d for e %f\n",binmin[i],ecut[i]);
        delete htmp; 
        htmp = fhDeltaCellClusterZE->ProjectionX(Form("hzde_cut%d",i),binmin,-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbZ);
        htmp->Draw("same HE");
        
      }
    }
    
    snprintf(name,buffersize,"QA_%s_DeltaClusterCellX_Y_Z_R_ECut.eps",fCalorimeter.Data());
    cxde->Print(name); printf("Create plot %s\n",name);
    
    
    //----------------------------------------------------------
    // Cluster-Cell X, Y, Z, R, NCells in cluster dependence
    //---------------------------------------------------------	
    snprintf(cname,buffersize,"%s_QA_DeltaClusterCellX_Y_Z_R_NCellsCut",fCalorimeter.Data());
    TCanvas  * cxdn = new TCanvas(cname, "Cluster-Cell X Y Z R, NCells cut", 800, 800) ;
    cxdn->Divide(2, 2);		
    //R
    cxdn->cd(1) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    
    TLegend pLegendXClDN(0.83,0.6,0.95,0.93);
    pLegendXClDN.SetTextSize(0.03);
    pLegendXClDN.SetFillColor(10);
    pLegendXClDN.SetBorderSize(1);
    delete htmp; 
    htmp = fhDeltaCellClusterRNCells->ProjectionX("hrdn_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDR);
      htmp->SetTitle("r clusters - r cells for n cells in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      pLegendXClDN.AddEntry(htmp,"No cut","L");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        delete htmp; 
        if(i < ncellcuts-1) htmp = fhDeltaCellClusterRNCells->ProjectionX(Form("hrdn_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhDeltaCellClusterRNCells->ProjectionX(Form("hrdn_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDR);
        htmp->Draw("same HE");
        if(i < ncellcuts-1) pLegendXClDN.AddEntry(htmp,Form("n = %1.1d",ncellcut[i]-1),"L");
        else pLegendXClDN.AddEntry(htmp,Form("n >= %1.1d",ncellcut[i]-1),"L");
        
      }
    }
    pLegendXClDN.Draw();
    
    //X
    cxdn->cd(2) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhDeltaCellClusterXNCells->ProjectionX("hxdn_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDX);
      htmp->SetTitle("x clusters - x cells for n cells in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        delete htmp; 
        if(i < ncellcuts-1)htmp = fhDeltaCellClusterXNCells->ProjectionX(Form("hxdn_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhDeltaCellClusterXNCells->ProjectionX(Form("hxdn_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDX);
        htmp->Draw("same HE");
        
      }
    }
    //Y
    cxdn->cd(3) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhDeltaCellClusterYNCells->ProjectionX("hydn_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDY);
      htmp->SetTitle("y clusters - y cells for n cells in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        delete htmp; 
        if(i < ncellcuts-1) htmp = fhDeltaCellClusterYNCells->ProjectionX(Form("hydn_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhDeltaCellClusterYNCells->ProjectionX(Form("hydn_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDY);
        htmp->Draw("same HE");
        
      }
    }
    //Z
    cxdn->cd(4) ; 
    gPad->SetLogy();
    gPad->SetGridy();
    delete htmp; 
    htmp = fhDeltaCellClusterZNCells->ProjectionX("hzdn_nocut",0,-1);
    if(htmp){
      htmp->SetMinimum(1);
      htmp->Rebin(rbDZ);
      htmp->SetTitle("z clusters - z cells for ncells in cluster > threshold");
      htmp->SetAxisRange(-50,50,"X");
      htmp->Draw("HE");
      
      for (Int_t i = 0; i < ncellcuts; i++) {
        delete htmp; 
        if(i < ncellcuts-1)htmp = fhDeltaCellClusterZNCells->ProjectionX(Form("hzdn_cut%d",i),ncellcut[i],ncellcut[i]);
        else htmp = fhDeltaCellClusterZNCells->ProjectionX(Form("hzdn_cut%d",i),ncellcut[i],-1);
        htmp->SetLineColor(ecutcolor[i]);
        htmp->Rebin(rbDZ);
        htmp->Draw("same HE");
        
      }
    }
    
    snprintf(name,buffersize,"QA_%s_DeltaClusterCellX_Y_Z_R_NCellsCut.eps",fCalorimeter.Data());
    cxdn->Print(name); printf("Create plot %s\n",name);
    
  }
  
  //----------------------------------------------------------
  //Reconstructed clusters energy-eta-phi distributions, matched with tracks
  //----------------------------------------------------------
  TH1F *	hEChargedClone   = 0 ;
  TH1F *	hPtChargedClone  = 0 ;
  TH1F *	hEtaChargedClone = 0 ;
  TH1F *	hPhiChargedClone = 0 ;
  if(fFillAllTH12){
    hEChargedClone   = (TH1F*)   fhECharged->Clone(Form("%sClone",fhECharged->GetName()));
    hPtChargedClone  = (TH1F*)   fhPtCharged->Clone(Form("%sClone",fhPtCharged->GetName()));
    hEtaChargedClone = (TH1F*)   fhEtaCharged->Clone(Form("%sClone",fhEtaCharged->GetName()));
    hPhiChargedClone = (TH1F*)   fhPhiCharged->Clone(Form("%sClone",fhPhiCharged->GetName()));
    
    snprintf(cname,buffersize,"QA_%s_rectrackmatch",fCalorimeter.Data());
    TCanvas  * ccltm = new TCanvas(cname, "Reconstructed clusters E-Phi-Eta, matched with tracks", 1200, 400) ;
    ccltm->Divide(3, 1);
    
    ccltm->cd(1) ; 
    if(fhECharged->GetEntries() > 0) gPad->SetLogy();
    fhECharged->Rebin(rbE);
    fhECharged->SetAxisRange(ptmin,ptmax,"X");
    fhECharged->SetMinimum(1);
    fhECharged->Draw();
    
    ccltm->cd(2) ; 
    if(fhPhiCharged->GetEntries() > 0) gPad->SetLogy();
    fhPhiCharged->Rebin(rbPhi);
    fhPhiCharged->SetAxisRange(phimin,phimax,"X");
    fhPhiCharged->Draw();
    fhPhiCharged->Draw();
    
    ccltm->cd(3) ;
    if(fhEtaCharged->GetEntries() > 0) gPad->SetLogy();
    fhEtaCharged->Rebin(rbEta);
    fhEtaCharged->SetAxisRange(etamin,etamax,"X");	
    fhEtaCharged->Draw();
    fhEtaCharged->Draw();
    
    snprintf(name,buffersize,"QA_%s_ClusterEnergyPhiEta_TrackMatched.eps",fCalorimeter.Data());
    ccltm->Print(name); printf("Plot: %s\n",name);
    
    //----------------------------------------------------------
    // Ratio  of reconstructed clusters energy-eta-phi distributions, matched with tracks over all
    //----------------------------------------------------------
    
    snprintf(cname,buffersize,"%s_QA_ChargedRatio",fCalorimeter.Data());
    TCanvas  * ccharge = new TCanvas(cname, "Charged clusters over all clusters", 1200, 400) ;
    ccharge->Divide(3, 1);
    
    ccharge->cd(1) ; 
    fhECharged->Sumw2();
    fhE->Sumw2();
    fhECharged->Divide(fhE);
    fhECharged->SetAxisRange(ptmin,ptmax,"X");
    fhECharged->SetMaximum(0.5);
    fhECharged->SetYTitle("track-matched clusters / all clusters");
    fhECharged->Draw("HE");
    
    ccharge->cd(2) ; 
    fhPhiCharged->Sumw2();
    fhPhi->Rebin(rbPhi);
    fhPhi->Sumw2();
    fhPhiCharged->Divide(fhPhi);
    fhPhiCharged->SetAxisRange(phimin,phimax,"X");
    fhPhiCharged->SetMaximum(0.5);
    fhPhiCharged->SetYTitle("track-matched clusters / all clusters");
    fhPhiCharged->Draw("HE");
    
    ccharge->cd(3) ; 
    fhEtaCharged->Sumw2();
    fhEta->Rebin(rbEta);
    fhEta->Sumw2();
    fhEtaCharged->Divide(fhEta);
    fhEtaCharged->SetAxisRange(etamin,etamax,"X");
    fhEtaCharged->SetMaximum(0.5);
    fhEtaCharged->SetYTitle("track-matched clusters / all clusters");
    fhEtaCharged->Draw("HE");
    
    snprintf(name,buffersize,"QA_%s_ClustersMatchedToAllRatios.eps",fCalorimeter.Data());
    ccharge->Print(name); printf("Create plot %s\n",name);
  }
  //-------------------------------------------	
  // N Cells - N Clusters - N Cells per cluster
  //-------------------------------------------
  snprintf(cname,buffersize,"QA_%s_nclustercells",fCalorimeter.Data());
  TCanvas  * cN = new TCanvas(cname, " Number of CaloClusters and CaloCells", 800, 1200) ;
  cN->Divide(2, 3);
  
  cN->cd(1) ; 
  
  TLegend pLegendN(0.7,0.6,0.9,0.8);
  pLegendN.SetTextSize(0.03);
  pLegendN.AddEntry(fhNClusters,"all modules","L");
  pLegendN.SetFillColor(10);
  pLegendN.SetBorderSize(1);
  
  if(fhNClusters->GetEntries() > 0) gPad->SetLogy();
  gPad->SetLogx();
  fhNClusters->SetLineColor(1);
  
  Int_t rbN = 1;
  if(fhNClusters->GetNbinsX()> nbins) rbN = fhNClusters->GetNbinsX()/nbins;
  
  fhNClusters->SetAxisRange(nmin,nmax,"X");
  fhNClusters->Draw("HE");
  for(Int_t imod = 0; imod < fNModules; imod++){
    fhNClustersMod[imod]->SetAxisRange(nmin,nmax,"X");
    fhNClustersMod[imod]->SetLineColor(modColorIndex[imod]);
    fhNClustersMod[imod]->Draw("same");
    pLegendN.AddEntry(fhNClustersMod[imod],Form("module %d",imod),"L");
  }
  pLegendN.Draw();
  
  cN->cd(2) ; 
  gPad->SetLogx();
  for(Int_t imod = 1; imod < fNModules; imod++){
    delete htmp; 
    htmp = (TH1D*)fhNClustersMod[imod]->Clone(Form("hNClustersRat%d",imod));
    htmp->Divide(fhNClustersMod[0]);
    htmp->SetLineColor(modColorIndex[imod]);
    if(imod==1){
      htmp->SetTitle("Ratio # clusters in  module X / module 0");
      htmp->SetMaximum(5);
      htmp->SetMinimum(0);
      htmp->Draw("HE");
    }
    else 
      htmp->Draw("same HE");
    
  }
  
  cN->cd(3) ; 
  if(fhNCells->GetEntries() > 0) gPad->SetLogy();
  gPad->SetLogx();
  fhNCells->SetLineColor(1);
  fhNCells->SetAxisRange(nmin,nmax,"X");
  fhNCells->Draw("HE");
  for(Int_t imod = 0; imod < fNModules; imod++){
    fhNCellsMod[imod]->SetAxisRange(nmin,nmax,"X");
    fhNCellsMod[imod]->SetLineColor(modColorIndex[imod]);
    fhNCellsMod[imod]->Draw("same HE");
  }
  
  
  cN->cd(4) ; 
  gPad->SetLogx();
  for(Int_t imod = 1; imod < fNModules; imod++){
    delete htmp; 
    htmp = (TH1D*)fhNCellsMod[imod]->Clone(Form("hNCellsRat%d",imod));
    htmp->Divide(fhNCellsMod[0]);
    htmp->SetLineColor(modColorIndex[imod]);
    if(imod==1){
      htmp->SetTitle("Ratio # cells in  module X / module 0");
      htmp->SetMaximum(5);
      htmp->SetMinimum(0);
      htmp->Draw("HE");
    }
    else 
      htmp->Draw("same HE");
    
  }
  
  cN->cd(5) ; 
  if(fhNCellsPerCluster->GetEntries() > 0) gPad->SetLogy();
  gPad->SetLogx();
  TH1D *cpc = fhNCellsPerCluster->ProjectionY("cpc",-1,-1);
  cpc->SetLineColor(1);
  cpc->SetTitle("# cells per cluster");
  cpc->Draw("HE"); 
  TH1D ** hNCellsCluster1D = new TH1D*[fNModules];
  
  for(Int_t imod = 0; imod < fNModules; imod++){
    hNCellsCluster1D[imod] = fhNCellsPerClusterMod[imod]->ProjectionY(Form("cpc_%d",imod),-1,-1);
    hNCellsCluster1D[imod]->SetLineColor(modColorIndex[imod]);
    hNCellsCluster1D[imod]->Draw("same HE");
  }
  
  
  cN->cd(6) ; 
  gPad->SetLogx();
  for(Int_t imod = 1; imod < fNModules; imod++){
    delete htmp; 
    htmp = (TH1D*)hNCellsCluster1D[imod]->Clone(Form("hNClustersCells1DRat%d",imod));
    htmp->Divide(hNCellsCluster1D[0]);
    htmp->SetLineColor(modColorIndex[imod]);
    if(imod==1){
      htmp->SetTitle("Ratio # cells per cluster in  module X / module 0");
      //htmp->SetAxisRange(ptmin,ptmax,"X");
      htmp->SetMaximum(3.5);
      htmp->SetMinimum(0);
      htmp->Draw("HE");
    }
    else 
      htmp->Draw("same HE");
  }
  delete [] hNCellsCluster1D;
  
  snprintf(name,buffersize,"QA_%s_NumberCaloClustersAndCaloCells.eps",fCalorimeter.Data());
  cN->Print(name); printf("Print plot %s\n",name);
  
  //----------------------------------------------------	
  // Cell Time histograms, time only available in ESDs
  //----------------------------------------------------
  if(GetReader()->GetDataType()==AliCaloTrackReader::kESD) {
    
    snprintf(cname,buffersize,"QA_%s_cellstime",fCalorimeter.Data());
    TCanvas  * ctime = new TCanvas(cname, " Cells time", 1200, 400) ;
    ctime->Divide(3, 1);
    
    Int_t rbTime = 1;
    if(fhTime->GetNbinsX()> ntimebins) rbTime = fhTime->GetNbinsX()/ntimebins;
    
    ctime->cd(1) ; 
    if(fhTime->GetEntries() > 0) gPad->SetLogy();
    fhTime->Rebin(rbTime);
    fhTime->SetAxisRange(timemin,timemax,"X");
    fhTime->Draw();
    
    ctime->cd(2) ; 
    fhTimeId->SetTitleOffset(1.8,"Y");
    fhTimeId->SetAxisRange(timemin,timemax,"X");
    fhTimeId->Draw("colz");
    
    ctime->cd(3) ; 
    fhTimeAmp->SetTitle("Cell Energy vs Cell Time");
    fhTimeAmp->SetTitleOffset(1.8,"Y");
    fhTimeAmp->SetAxisRange(timemin,timemax,"Y");
    fhTimeAmp->SetAxisRange(ptmin,ptmax,"X");		
    fhTimeAmp->Draw("colz");
    
    snprintf(name,buffersize,"QA_%s_CellsTime.eps",fCalorimeter.Data());
    ctime->Print(name); printf("Plot: %s\n",name);
  }
  
  
  //---------------------------------
  //Grid of cell per module plots 
  //---------------------------------
  {
    //Number of entries per cell
    gStyle->SetPadRightMargin(0.15);
    snprintf(cname,buffersize,"%s_QA_GridCellEntries",fCalorimeter.Data());
    TCanvas *cgrid   = new TCanvas("cgrid","Number of entries per cell", 12,12,800,400);
    if(fNModules%2 == 0)
      cgrid->Divide(fNModules/2,2); 
    else
      cgrid->Divide(fNModules/2+1,2); 
		
    for(Int_t imod = 0; imod < fNModules ; imod++){
      cgrid->cd(imod+1);
      gPad->SetLogz();
      gPad->SetGridy();
      gPad->SetGridx();
      //fhGridCellsMod[imod]->GetYAxis()->SetTitleColor(1);
      fhGridCellsMod[imod]->SetZTitle("Counts    ");
      fhGridCellsMod[imod]->SetYTitle("row (phi direction)    ");
      //fhGridCellsMod[imod]->SetLabelSize(0.025,"z");
      fhGridCellsMod[imod]->Draw("colz");
    }
    snprintf(name,buffersize,"QA_%s_GridCellsEntries.eps",fCalorimeter.Data());
    cgrid->Print(name); printf("Create plot %s\n",name);
    
    snprintf(cname,buffersize,"%s_QA_GridCellAccumEnergy",fCalorimeter.Data());
    TCanvas *cgridE   = new TCanvas("cgridE","Summed energy per cell", 12,12,800,400);
    if(fNModules%2 == 0)
      cgridE->Divide(fNModules/2,2); 
    else
      cgridE->Divide(fNModules/2+1,2); 
    for(Int_t imod = 0; imod < fNModules ; imod++){
      cgridE->cd(imod+1);
      gPad->SetLogz();
      gPad->SetGridy();
      gPad->SetGridx();
      //fhGridCellsEMod[imod]->SetLabelSize(0.025,"z");
      fhGridCellsEMod[imod]->SetZTitle("Accumulated Energy (GeV)    ");
      fhGridCellsEMod[imod]->SetYTitle("row (phi direction)    ");
      fhGridCellsEMod[imod]->Draw("colz");
    }
    snprintf(name,buffersize,"QA_%s_GridCellsAccumEnergy.eps",fCalorimeter.Data());
    cgridE->Print(name); printf("Create plot %s\n",name);
    
    //Accumulated energy per cell
    snprintf(cname,buffersize,"%s_QA_GridCellAverageEnergy",fCalorimeter.Data());
    TCanvas *cgridEA   = new TCanvas("cgridEA","Average energy per cell", 12,12,800,400);
    if(fNModules%2 == 0)	  
      cgridEA->Divide(fNModules/2,2);
    else
      cgridEA->Divide(fNModules/2+1,2);  
    for(Int_t imod = 0; imod < fNModules ; imod++){
      cgridEA->cd(imod+1);
      gPad->SetLogz();
      gPad->SetGridy();
      gPad->SetGridx();
      //fhGridCellsEMod[imod]->SetLabelSize(0.025,"z");
      fhGridCellsEMod[imod]->SetZTitle("Average Energy (GeV)    ");
      fhGridCellsEMod[imod]->Divide(fhGridCellsMod[imod]);
      fhGridCellsEMod[imod]->Draw("colz");
    }
    snprintf(name,buffersize,"QA_%s_GridCellsAverageEnergy.eps",fCalorimeter.Data());
    cgridEA->Print(name); printf("Create plot %s\n",name);
		
    //Accumulated Time per cell, E > 0.5 GeV
		
    snprintf(cname,buffersize,"%s_QA_GridCellAccumTime",fCalorimeter.Data());
    TCanvas *cgridT   = new TCanvas("cgridT","Summed time per cell", 12,12,800,400);
    if(fNModules%2 == 0)
      cgridT->Divide(fNModules/2,2); 
    else
      cgridE->Divide(fNModules/2+1,2); 
    for(Int_t imod = 0; imod < fNModules ; imod++){
      cgridT->cd(imod+1);
      gPad->SetLogz();
      gPad->SetGridy();
      gPad->SetGridx();
      //fhGridCellsTimeMod[imod]->SetLabelSize(0.025,"z");
      fhGridCellsTimeMod[imod]->SetZTitle("Accumulated Time (ns)    ");
      fhGridCellsTimeMod[imod]->SetYTitle("row (phi direction)    ");
      fhGridCellsTimeMod[imod]->Draw("colz");
    }
    snprintf(name,buffersize,"QA_%s_GridCellsAccumTime.eps",fCalorimeter.Data());
    cgridT->Print(name); printf("Create plot %s\n",name);
		
  }
  
  //---------------------------------------------
  //Calorimeter Correlation, PHOS vs EMCAL
  //---------------------------------------------
  if(fCorrelate){
    
    snprintf(cname,buffersize,"QA_%s_CaloCorr_EMCALvsPHOS",fCalorimeter.Data());
    TCanvas  * ccorr = new TCanvas(cname, " EMCAL vs PHOS", 400, 400) ;
    ccorr->Divide(2, 2);
    
    ccorr->cd(1) ; 
    //gPad->SetLogy();
    //gPad->SetLogx();
    fhCaloCorrNClusters->SetAxisRange(nmin,nmax,"X");
    fhCaloCorrNClusters->SetAxisRange(nmin,nmax,"Y");		
    fhCaloCorrNClusters ->Draw();
    
    ccorr->cd(2) ; 
    //gPad->SetLogy();
    //gPad->SetLogx();
    fhCaloCorrNCells->SetAxisRange(nmin,nmax,"X");
    fhCaloCorrNCells->SetAxisRange(nmin,nmax,"Y");		
    fhCaloCorrNCells->Draw();
    
    //gPad->SetLogy();
    //gPad->SetLogx();
    fhCaloCorrEClusters->SetAxisRange(ptmin,ptmax,"X");
    fhCaloCorrEClusters->SetAxisRange(ptmin,ptmax,"Y");		
    fhCaloCorrEClusters->Draw();
    
    ccorr->cd(4) ; 
    //gPad->SetLogy();
    //gPad->SetLogx();
    fhCaloCorrECells->SetAxisRange(ptmin,ptmax,"X");
    fhCaloCorrECells->SetAxisRange(ptmin,ptmax,"Y");		
    fhCaloCorrECells->Draw();
    
    snprintf(name,buffersize,"QA_%s_CaloCorr_EMCALvsPHOS.eps",fCalorimeter.Data());
    ccorr->Print(name); printf("Plot: %s\n",name);
  }
  
  //----------------------------
  //Invariant mass
  //-----------------------------
	
  Int_t imbinmin = -1;
  Int_t imbinmax = -1;
  
  if(fhIM->GetEntries() > 1){
    Int_t nebins  = fhIM->GetNbinsX();
    Int_t emax = (Int_t) fhIM->GetXaxis()->GetXmax();
    Int_t emin = (Int_t) fhIM->GetXaxis()->GetXmin();
    if (emin != 0 ) printf("emin != 0 \n");
    //printf("IM: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
    
    snprintf(cname,buffersize,"QA_%s_IM",fCalorimeter.Data());
    //	printf("c5\n");
    TCanvas  * c5 = new TCanvas(cname, "Invariant mass", 600, 400) ;
    c5->Divide(2, 3);
    
    c5->cd(1) ; 
    //fhIM->SetLineColor(4);
    //fhIM->Draw();
    imbinmin = 0;
    imbinmax =  (Int_t) (1-emin)*nebins/emax;
    TH1D *pyim1 = fhIM->ProjectionY(Form("%s_py1",fhIM->GetName()),imbinmin,imbinmax);
    pyim1->SetTitle("E_{pair} < 1 GeV");
    pyim1->SetLineColor(1);
    pyim1->Draw();
    TLegend pLegendIM(0.7,0.6,0.9,0.8);
    pLegendIM.SetTextSize(0.03);
    pLegendIM.AddEntry(pyim1,"all modules","L");
    pLegendIM.SetFillColor(10);
    pLegendIM.SetBorderSize(1);
    //FIXME
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyim1 = fhIMMod[imod]->ProjectionY(Form("%s_py1",fhIMMod[imod]->GetName()),imbinmin,imbinmax);
      pLegendIM.AddEntry(pyim1,Form("module %d",imod),"L");
      pyim1->SetLineColor(imod+1);
      pyim1->Draw("same");
    }
    pLegendIM.Draw();
    
    c5->cd(2) ; 
    imbinmin =  (Int_t) (1-emin)*nebins/emax;
    imbinmax =  (Int_t) (2-emin)*nebins/emax;
    TH1D *pyim2 = fhIM->ProjectionY(Form("%s_py2",fhIM->GetName()),imbinmin,imbinmax);
    pyim2->SetTitle("1 < E_{pair} < 2 GeV");
    pyim2->SetLineColor(1);
    pyim2->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyim2 = fhIMMod[imod]->ProjectionY(Form("%s_py2",fhIMMod[imod]->GetName()),imbinmin,imbinmax);
      pyim2->SetLineColor(imod+1);
      pyim2->Draw("same");
    }
    
    c5->cd(3) ; 
    imbinmin =  (Int_t) (2-emin)*nebins/emax;
    imbinmax =  (Int_t) (3-emin)*nebins/emax;
    TH1D *pyim3 = fhIM->ProjectionY(Form("%s_py3",fhIM->GetName()),imbinmin,imbinmax);
    pyim3->SetTitle("2 < E_{pair} < 3 GeV");
    pyim3->SetLineColor(1);
    pyim3->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyim3 = fhIMMod[imod]->ProjectionY(Form("%s_py3",fhIMMod[imod]->GetName()),imbinmin,imbinmax);
      pyim3->SetLineColor(imod+1);
      pyim3->Draw("same");
    }
    
    c5->cd(4) ;
    imbinmin =  (Int_t) (3-emin)*nebins/emax;
    imbinmax =  (Int_t) (4-emin)*nebins/emax;
    TH1D *pyim4 = fhIM->ProjectionY(Form("%s_py4",fhIM->GetName()),imbinmin,imbinmax);
    pyim4->SetTitle("3 < E_{pair} < 4 GeV");
    pyim4->SetLineColor(1);
    pyim4->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyim4 = fhIMMod[imod]->ProjectionY(Form("%s_py4",fhIMMod[imod]->GetName()),imbinmin,imbinmax);
      pyim4->SetLineColor(imod+1);
      pyim4->Draw("same");
    }
    
    c5->cd(5) ;
    imbinmin =  (Int_t) (4-emin)*nebins/emax;
    imbinmax =  (Int_t) (5-emin)*nebins/emax;
    TH1D *pyim5 = fhIM->ProjectionY(Form("%s_py5",fhIM->GetName()),imbinmin,imbinmax);
    pyim5->SetTitle("4< E_{pair} < 5 GeV");
    pyim5->SetLineColor(1);
    pyim5->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyim5 = fhIMMod[imod]->ProjectionY(Form("%s_py5",fhIMMod[imod]->GetName()),imbinmin,imbinmax);
      pyim5->SetLineColor(imod+1);
      pyim5->Draw("same");
    }
    
    c5->cd(6) ;
    imbinmin =  (Int_t) (5-emin)*nebins/emax;
    imbinmax =  -1;
    TH1D *pyim10 = fhIM->ProjectionY(Form("%s_py6",fhIM->GetName()),imbinmin,imbinmax);
    pyim10->SetTitle("E_{pair} > 5 GeV");
    pyim10->SetLineColor(1);
    pyim10->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyim10 = fhIMMod[imod]->ProjectionY(Form("%s_py6",fhIMMod[imod]->GetName()),imbinmin,imbinmax);
      pyim10->SetLineColor(imod+1);
      pyim10->Draw("same");
    }
    
    snprintf(name,buffersize,"QA_%s_InvariantMass.eps",fCalorimeter.Data());
    c5->Print(name); printf("Plot: %s\n",name);
  }
  
  //--------------------------------------------------
  //Invariant mass, clusters with more than one cell
  //-------------------------------------------------
  if(fhIMCellCut->GetEntries() > 1){
    Int_t nebins  = fhIMCellCut->GetNbinsX();
    Int_t emax = (Int_t) fhIMCellCut->GetXaxis()->GetXmax();
    Int_t emin = (Int_t) fhIMCellCut->GetXaxis()->GetXmin();
    if (emin != 0 ) printf("emin != 0 \n");
    //printf("IMCellCut: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
		
    snprintf(cname,buffersize,"QA_%s_IMCellCut",fCalorimeter.Data());
    //	printf("c5cc\n");
    TCanvas  * c5cc = new TCanvas(cname, "Invariant mass, Cell Cut", 600, 400) ;
    c5cc->Divide(2, 3);
    
    c5cc->cd(1) ; 
    //fhIMCellCut->SetLineColor(4);
    //fhIMCellCut->Draw();
    imbinmin = 0;
    imbinmax =  (Int_t) (1-emin)*nebins/emax;
    TH1D *pyimcc1 = fhIMCellCut->ProjectionY(Form("%s_py1",fhIMCellCut->GetName()),imbinmin,imbinmax);
    pyimcc1->SetTitle("E_{pair} < 1 GeV");
    pyimcc1->SetLineColor(1);
    pyimcc1->Draw();
    TLegend pLegendIMCellCut(0.7,0.6,0.9,0.8);
    pLegendIMCellCut.SetTextSize(0.03);
    pLegendIMCellCut.AddEntry(pyimcc1,"all modules","L");
    pLegendIMCellCut.SetFillColor(10);
    pLegendIMCellCut.SetBorderSize(1);
    
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyimcc1 = fhIMCellCutMod[imod]->ProjectionY(Form("%s_py1",fhIMCellCutMod[imod]->GetName()),imbinmin,imbinmax);
      pLegendIMCellCut.AddEntry(pyimcc1,Form("module %d",imod),"L");
      pyimcc1->SetLineColor(imod+1);
      pyimcc1->Draw("same");
    }
    pLegendIMCellCut.Draw();
    
    c5cc->cd(2) ; 
    imbinmin =  (Int_t) (1-emin)*nebins/emax;
    imbinmax =  (Int_t) (2-emin)*nebins/emax;
    TH1D *pyimcc2 = fhIMCellCut->ProjectionY(Form("%s_py2",fhIMCellCut->GetName()),imbinmin,imbinmax);
    pyimcc2->SetTitle("1 < E_{pair} < 2 GeV");
    pyimcc2->SetLineColor(1);
    pyimcc2->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyimcc2 = fhIMCellCutMod[imod]->ProjectionY(Form("%s_py1",fhIMCellCutMod[imod]->GetName()),imbinmin,imbinmax);
      pyimcc2->SetLineColor(imod+1);
      pyimcc2->Draw("same");
    }
    
    c5cc->cd(3) ; 
    imbinmin =  (Int_t) (2-emin)*nebins/emax;
    imbinmax =  (Int_t) (3-emin)*nebins/emax;
    TH1D *pyimcc3 = fhIMCellCut->ProjectionY(Form("%s_py3",fhIMCellCut->GetName()),imbinmin,imbinmax);
    pyimcc3->SetTitle("2 < E_{pair} < 3 GeV");
    pyimcc3->SetLineColor(1);
    pyimcc3->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyimcc3 = fhIMCellCutMod[imod]->ProjectionY(Form("%s_py1",fhIMCellCutMod[imod]->GetName()),imbinmin,imbinmax);
      pyimcc3->SetLineColor(imod+1);
      pyimcc3->Draw("same");
    }
    
    c5cc->cd(4) ;
    imbinmin =  (Int_t) (3-emin)*nebins/emax;
    imbinmax =  (Int_t) (4-emin)*nebins/emax;
    TH1D *pyimcc4 = fhIMCellCut->ProjectionY(Form("%s_py4",fhIMCellCut->GetName()),imbinmin,imbinmax);
    pyimcc4->SetTitle("3 < E_{pair} < 4 GeV");
    pyimcc4->SetLineColor(1);
    pyimcc4->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyimcc4 = fhIMCellCutMod[imod]->ProjectionY(Form("%s_py5",fhIMCellCutMod[imod]->GetName()),imbinmin,imbinmax);
      pyimcc4->SetLineColor(imod+1);
      pyimcc4->Draw("same");
    }
    
    c5cc->cd(5) ;
    imbinmin =  (Int_t) (4-emin)*nebins/emax;
    imbinmax =  (Int_t) (5-emin)*nebins/emax;
    TH1D *pyimcc5cc = fhIMCellCut->ProjectionY(Form("%s_py5",fhIMCellCut->GetName()),imbinmin,imbinmax);
    pyimcc5cc->SetTitle("4< E_{pair} < 5 GeV");
    pyimcc5cc->SetLineColor(1);
    pyimcc5cc->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyimcc5cc = fhIMCellCutMod[imod]->ProjectionY(Form("%s_py5",fhIMCellCutMod[imod]->GetName()),imbinmin,imbinmax);
      pyimcc5cc->SetLineColor(imod+1);
      pyimcc5cc->Draw("same");
    }
    
    c5cc->cd(6) ;
    imbinmin =  (Int_t) (5-emin)*nebins/emax;
    imbinmax =  -1;
    TH1D *pyimcc10 = fhIMCellCut->ProjectionY(Form("%s_py6",fhIMCellCut->GetName()),imbinmin,imbinmax);
    pyimcc10->SetTitle("E_{pair} > 5 GeV");
    pyimcc10->SetLineColor(1);
    pyimcc10->Draw();
    for(Int_t imod = 0; imod < fNModules; imod++){
      pyimcc10 = fhIMCellCutMod[imod]->ProjectionY(Form("%s_py1",fhIMCellCutMod[imod]->GetName()),imbinmin,imbinmax);
      pyimcc10->SetLineColor(imod+1);
      pyimcc10->Draw("same");
    }
    
    snprintf(name,buffersize,"QA_%s_InvariantMass_CellCut.eps",fCalorimeter.Data());
    c5cc->Print(name); printf("Plot: %s\n",name);
  }
  
  
  //Asymmetry
  if(fhAsym->GetEntries() > 1){
    Int_t nebins  = fhAsym->GetNbinsX();
    Int_t emax = (Int_t) fhAsym->GetXaxis()->GetXmax();
    Int_t emin = (Int_t) fhAsym->GetXaxis()->GetXmin();
    if (emin != 0 ) printf("emin != 0 \n");
    //printf("Asym: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
    
    snprintf(cname,buffersize,"QA_%s_Asym",fCalorimeter.Data());
    //	printf("c5\n");
    TCanvas  * c5b = new TCanvas(cname, "Asymmetry", 400, 400) ;
    c5b->Divide(2, 2);
    
    c5b->cd(1) ; 
    fhAsym->SetTitleOffset(1.6,"Y");
    fhAsym->SetLineColor(4);
    fhAsym->Draw();
    
    c5b->cd(2) ; 
    imbinmin = 0;
    imbinmax = (Int_t) (5-emin)*nebins/emax;
    TH1D *pyAsym5 = fhAsym->ProjectionY(Form("%s_py5",fhAsym->GetName()),imbinmin,imbinmax);
    pyAsym5->SetTitle("E_{pair} < 5 GeV");
    pyAsym5->SetLineColor(4);
    pyAsym5->Draw();
    
    c5b->cd(3) ; 
    imbinmin = (Int_t) (5-emin)*nebins/emax;
    imbinmax = (Int_t) (10-emin)*nebins/emax;
    TH1D *pyAsym510 = fhAsym->ProjectionY(Form("%s_py510",fhAsym->GetName()),imbinmin,imbinmax);
    pyAsym510->SetTitle("5 < E_{pair} < 10 GeV");
    pyAsym510->SetLineColor(4);
    pyAsym510->Draw();
    
    c5b->cd(4) ;
    imbinmin = (Int_t) (10-emin)*nebins/emax;
    imbinmax = -1;
    TH1D *pyAsym10 = fhAsym->ProjectionY(Form("%s_py10",fhAsym->GetName()),imbinmin,imbinmax);
    pyAsym10->SetTitle("E_{pair} > 10 GeV");
    pyAsym10->SetLineColor(4);
    pyAsym10->Draw();
    
    snprintf(name,buffersize,"QA_%s_Asymmetry.eps",fCalorimeter.Data());
    c5b->Print(name); printf("Plot: %s\n",name);
  }
  
  
  if(IsDataMC()){
    //Reconstructed vs MC distributions
    //printf("c6\n");
    snprintf(cname,buffersize,"QA_%s_recvsmc",fCalorimeter.Data());
    TCanvas  * c6 = new TCanvas(cname, "Reconstructed vs MC distributions", 400, 400) ;
    c6->Divide(2, 2);
    
    c6->cd(1) ; 
    fh2E->SetTitleOffset(1.6,"Y");
    fh2E->SetLineColor(4);
    fh2E->Draw();
    
    c6->cd(2) ; 
    fh2Pt->SetTitleOffset(1.6,"Y");
    fh2Pt->SetLineColor(4);
    fh2Pt->Draw();
    
    c6->cd(3) ; 
    fh2Phi->SetTitleOffset(1.6,"Y");
    fh2Phi->SetLineColor(4);
    fh2Phi->Draw();
    
    c6->cd(4) ; 
    fh2Eta->SetTitleOffset(1.6,"Y");
    fh2Eta->SetLineColor(4);
    fh2Eta->Draw();
    
    snprintf(name,buffersize,"QA_%s_ReconstructedVSMCDistributions.eps",fCalorimeter.Data());
    c6->Print(name); printf("Plot: %s\n",name);	
    
    //Reconstructed vs MC distributions
    //printf("c6\n");
    snprintf(cname,buffersize,"QA_%s_gamrecvsmc",fCalorimeter.Data());
    TCanvas  * c6Gam = new TCanvas(cname, "Reconstructed vs MC distributions", 400, 400) ;
    c6Gam->Divide(2, 2);
    
    c6Gam->cd(1) ; 
    fhGamE->Draw();
    
    c6Gam->cd(2) ; 
    fhGamPt->Draw();
    
    c6Gam->cd(3) ; 
    fhGamPhi->Draw();
    
    c6Gam->cd(4) ; 
    fhGamEta->Draw();
    
    snprintf(name,buffersize,"QA_%s_GammaReconstructedVSMCDistributions.eps",fCalorimeter.Data());
    c6->Print(name); printf("Plot: %s\n",name);	
    
    //Generated - reconstructed  
    //printf("c7\n");
    snprintf(cname,buffersize,"QA_%s_diffgenrec",fCalorimeter.Data());
    TCanvas  * c7 = new TCanvas(cname, "generated - reconstructed", 400, 400) ;
    c7->Divide(2, 2);
    
    c7->cd(1) ; 
    if(fhDeltaE->GetEntries() > 0) gPad->SetLogy();
    fhGamDeltaE->SetLineColor(4);
    fhDeltaE->Draw();
    fhGamDeltaE->Draw("same");
    
    TLegend pLegendd(0.65,0.55,0.9,0.8);
    pLegendd.SetTextSize(0.06);
    pLegendd.AddEntry(fhDeltaE,"all","L");
    pLegendd.AddEntry(fhGamDeltaE,"from  #gamma","L");
    pLegendd.SetFillColor(10);
    pLegendd.SetBorderSize(1);
    pLegendd.Draw();
    
    c7->cd(2) ; 
    if(fhDeltaPt->GetEntries() > 0) gPad->SetLogy();
    fhGamDeltaPt->SetLineColor(4);
    fhDeltaPt->Draw();
    fhGamDeltaPt->Draw("same");
    
    c7->cd(3) ; 
    fhGamDeltaPhi->SetLineColor(4);
    fhDeltaPhi->Draw();
    fhGamDeltaPhi->Draw("same");
    
    c7->cd(4) ; 
    fhGamDeltaEta->SetLineColor(4);
    fhDeltaEta->Draw();
    fhGamDeltaEta->Draw("same");
    
    snprintf(name,buffersize,"QA_%s_DiffGeneratedReconstructed.eps",fCalorimeter.Data());
    c7->Print(name); printf("Plot: %s\n",name);
    
    // Reconstructed / Generated 
    //printf("c8\n");
    snprintf(cname,buffersize,"QA_%s_ratiorecgen",fCalorimeter.Data());
    TCanvas  * c8 = new TCanvas(cname, " reconstructed / generated", 400, 400) ;
    c8->Divide(2, 2);
    
    c8->cd(1) ; 
    if(fhRatioE->GetEntries() > 0) gPad->SetLogy();
    fhGamRatioE->SetLineColor(4);
    fhRatioE->Draw();
    fhGamRatioE->Draw("same");
    
    TLegend pLegendr(0.65,0.55,0.9,0.8);
    pLegendr.SetTextSize(0.06);
    pLegendr.AddEntry(fhRatioE,"all","L");
    pLegendr.AddEntry(fhGamRatioE,"from  #gamma","L");
    pLegendr.SetFillColor(10);
    pLegendr.SetBorderSize(1);
    pLegendr.Draw();
    
    c8->cd(2) ; 
    if(fhRatioPt->GetEntries() > 0) gPad->SetLogy();
    fhGamRatioPt->SetLineColor(4);
    fhRatioPt->Draw();
    fhGamRatioPt->Draw("same");
    
    c8->cd(3) ; 
    fhGamRatioPhi->SetLineColor(4);
    fhRatioPhi->Draw();
    fhGamRatioPhi->Draw("same");
    
    c8->cd(4) ; 
    fhGamRatioEta->SetLineColor(4);
    fhRatioEta->Draw();
    fhGamRatioEta->Draw("same");
    
    snprintf(name,buffersize,"QA_%s_ReconstructedDivGenerated.eps",fCalorimeter.Data());
    c8->Print(name); printf("Plot: %s\n",name);
    
    //MC
    
    //Generated distributions
    //printf("c1\n");
    snprintf(cname,buffersize,"QA_%s_gen",fCalorimeter.Data());
    TCanvas  * c10 = new TCanvas(cname, "Generated distributions", 600, 200) ;
    c10->Divide(3, 1);
    
    c10->cd(1) ; 
    gPad->SetLogy();
    TH1F * haxispt  = (TH1F*) fhGenPi0Pt->Clone(Form("%s_axispt",fhGenPi0Pt->GetName()));  
    haxispt->SetTitle("Generated Particles p_{T}, |#eta| < 1");
    fhGenPi0Pt->SetLineColor(1);
    fhGenGamPt->SetLineColor(4);
    fhGenEtaPt->SetLineColor(2);
    fhGenOmegaPt->SetLineColor(7);
    fhGenElePt->SetLineColor(6);
    
    //Select the maximum of the histogram to show all lines.
    if(fhGenPi0Pt->GetMaximum() >= fhGenGamPt->GetMaximum() && fhGenPi0Pt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
       fhGenPi0Pt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenPi0Pt->GetMaximum() >= fhGenElePt->GetMaximum())
      haxispt->SetMaximum(fhGenPi0Pt->GetMaximum());
    else if(fhGenGamPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenGamPt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
            fhGenGamPt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenGamPt->GetMaximum() >= fhGenElePt->GetMaximum())
      haxispt->SetMaximum(fhGenGamPt->GetMaximum());
    else if(fhGenEtaPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenEtaPt->GetMaximum() >= fhGenGamPt->GetMaximum() && 
            fhGenEtaPt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenEtaPt->GetMaximum() >= fhGenElePt->GetMaximum())
      haxispt->SetMaximum(fhGenEtaPt->GetMaximum());	
    else if(fhGenOmegaPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenOmegaPt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
            fhGenOmegaPt->GetMaximum() >= fhGenGamPt->GetMaximum() && fhGenOmegaPt->GetMaximum() >= fhGenElePt->GetMaximum())
      haxispt->SetMaximum(fhGenOmegaPt->GetMaximum());
    else if(fhGenElePt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenElePt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
            fhGenElePt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenElePt->GetMaximum() >= fhGenGamPt->GetMaximum())
      haxispt->SetMaximum(fhGenElePt->GetMaximum());
    haxispt->SetMinimum(1);
    haxispt->Draw("axis");
    fhGenPi0Pt->Draw("same");
    fhGenGamPt->Draw("same");
    fhGenEtaPt->Draw("same");
    fhGenOmegaPt->Draw("same");
    fhGenElePt->Draw("same");
    
    TLegend pLegend(0.85,0.65,0.95,0.93);
    pLegend.SetTextSize(0.06);
    pLegend.AddEntry(fhGenPi0Pt,"  #pi^{0}","L");
    pLegend.AddEntry(fhGenGamPt,"  #gamma","L");
    pLegend.AddEntry(fhGenEtaPt,"  #eta","L");
    pLegend.AddEntry(fhGenOmegaPt,"  #omega","L");
    pLegend.AddEntry(fhGenElePt,"  e^{#pm}","L");
    pLegend.SetFillColor(10);
    pLegend.SetBorderSize(1);
    pLegend.Draw();
    
    c10->cd(2) ;
    gPad->SetLogy();
    TH1F * haxiseta  = (TH1F*) fhGenPi0Eta->Clone(Form("%s_axiseta",fhGenPi0Eta->GetName()));  
    haxiseta->SetTitle("Generated Particles #eta, |#eta| < 1");
    fhGenPi0Eta->SetLineColor(1);
    fhGenGamEta->SetLineColor(4);
    fhGenEtaEta->SetLineColor(2);
    fhGenOmegaEta->SetLineColor(7);
    fhGenEleEta->SetLineColor(6);
    //Select the maximum of the histogram to show all lines.
    if(fhGenPi0Eta->GetMaximum() >= fhGenGamEta->GetMaximum() && fhGenPi0Eta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
       fhGenPi0Eta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenPi0Eta->GetMaximum() >= fhGenEleEta->GetMaximum())
      haxiseta->SetMaximum(fhGenPi0Eta->GetMaximum());
    else if(fhGenGamEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenGamEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
            fhGenGamEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenGamEta->GetMaximum() >= fhGenEleEta->GetMaximum())
      haxiseta->SetMaximum(fhGenGamEta->GetMaximum());
    else if(fhGenEtaEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenEtaEta->GetMaximum() >= fhGenGamEta->GetMaximum() && 
            fhGenEtaEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenEtaEta->GetMaximum() >= fhGenEleEta->GetMaximum())
      haxiseta->SetMaximum(fhGenEtaEta->GetMaximum());	
    else if(fhGenOmegaEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenOmegaEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
            fhGenOmegaEta->GetMaximum() >= fhGenGamEta->GetMaximum() && fhGenOmegaEta->GetMaximum() >= fhGenEleEta->GetMaximum())
      haxiseta->SetMaximum(fhGenOmegaEta->GetMaximum());
    else if(fhGenEleEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenEleEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
            fhGenEleEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenEleEta->GetMaximum() >= fhGenGamEta->GetMaximum())
      haxiseta->SetMaximum(fhGenEleEta->GetMaximum());
    haxiseta->SetMinimum(100);
    haxiseta->Draw("axis");
    fhGenPi0Eta->Draw("same");
    fhGenGamEta->Draw("same");
    fhGenEtaEta->Draw("same");
    fhGenOmegaEta->Draw("same");
    fhGenEleEta->Draw("same");
    
    
    c10->cd(3) ; 
    gPad->SetLogy();
    TH1F * haxisphi  = (TH1F*) fhGenPi0Phi->Clone(Form("%s_axisphi",fhGenPi0Phi->GetName()));  
    haxisphi->SetTitle("Generated Particles #phi, |#eta| < 1");
    fhGenPi0Phi->SetLineColor(1);
    fhGenGamPhi->SetLineColor(4);
    fhGenEtaPhi->SetLineColor(2);
    fhGenOmegaPhi->SetLineColor(7);
    fhGenElePhi->SetLineColor(6);
    //Select the maximum of the histogram to show all lines.
    if(fhGenPi0Phi->GetMaximum() >= fhGenGamPhi->GetMaximum() && fhGenPi0Phi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
       fhGenPi0Phi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenPi0Phi->GetMaximum() >= fhGenElePhi->GetMaximum())
      haxisphi->SetMaximum(fhGenPi0Phi->GetMaximum());
    else if(fhGenGamPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenGamPhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
            fhGenGamPhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenGamPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
      haxisphi->SetMaximum(fhGenGamPhi->GetMaximum());
    else if(fhGenEtaPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenEtaPhi->GetMaximum() >= fhGenGamPhi->GetMaximum() && 
            fhGenEtaPhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenEtaPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
      haxisphi->SetMaximum(fhGenEtaPhi->GetMaximum());	
    else if(fhGenOmegaPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenOmegaPhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
            fhGenOmegaPhi->GetMaximum() >= fhGenGamPhi->GetMaximum() && fhGenOmegaPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
      haxisphi->SetMaximum(fhGenOmegaPhi->GetMaximum());
    else if(fhGenElePhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenElePhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
            fhGenElePhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenElePhi->GetMaximum() >= fhGenGamPhi->GetMaximum())
      haxisphi->SetMaximum(fhGenElePhi->GetMaximum());
    haxisphi->SetMinimum(100);
    haxisphi->Draw("axis");
    fhGenPi0Phi->Draw("same");
    fhGenGamPhi->Draw("same");
    fhGenEtaPhi->Draw("same");
    fhGenOmegaPhi->Draw("same");
    fhGenElePhi->Draw("same");
    
    snprintf(name,buffersize,"QA_%s_GeneratedDistributions.eps",fCalorimeter.Data());
    c10->Print(name); printf("Plot: %s\n",name);
    
    
    //Reconstructed clusters depending on its original particle.
    //printf("c1\n");
    snprintf(cname,buffersize,"QA_%s_recgenid",fCalorimeter.Data());
    TCanvas  * c11 = new TCanvas(cname, "Reconstructed particles, function of their original particle ID", 400, 400) ;
    c11->Divide(2, 2);
    
    
    c11->cd(1) ; 
    gPad->SetLogy();
    TH1F * hGamE   = (TH1F*) fhGamE->ProjectionX(Form("%s_px",fhGamE->GetName()),-1,-1);
    TH1F * hPi0E   = (TH1F*) fhPi0E->ProjectionX(Form("%s_px",fhPi0E->GetName()),-1,-1);
    TH1F * hEleE   = (TH1F*) fhEleE->ProjectionX(Form("%s_px",fhEleE->GetName()),-1,-1);
    TH1F * hNeHadE = (TH1F*) fhNeHadE->ProjectionX(Form("%s_px",fhNeHadE->GetName()),-1,-1);
    TH1F * hChHadE = (TH1F*) fhChHadE->ProjectionX(Form("%s_px",fhChHadE->GetName()),-1,-1);
    TH1F * haxisE  = (TH1F*) hPi0E->Clone(Form("%s_axisE",fhPi0E->GetName()));  
    haxisE->SetTitle("Reconstructed particles E, function of their original particle ID");
    hPi0E->SetLineColor(1);
    hGamE->SetLineColor(4);
    hNeHadE->SetLineColor(2);
    hChHadE->SetLineColor(7);
    hEleE->SetLineColor(6);
    
    //Select the maximum of the histogram to show all lines.
    if(hPi0E->GetMaximum() >= hGamE->GetMaximum() && hPi0E->GetMaximum() >= hNeHadE->GetMaximum() && 
       hPi0E->GetMaximum() >= hChHadE->GetMaximum() && hPi0E->GetMaximum() >= hEleE->GetMaximum())
      haxisE->SetMaximum(hPi0E->GetMaximum());
    else if(hGamE->GetMaximum() >= hPi0E->GetMaximum() && hGamE->GetMaximum() >= hNeHadE->GetMaximum() && 
            hGamE->GetMaximum() >= hChHadE->GetMaximum() && hGamE->GetMaximum() >= hEleE->GetMaximum())
      haxisE->SetMaximum(hGamE->GetMaximum());
    else if(hNeHadE->GetMaximum() >= hPi0E->GetMaximum() && hNeHadE->GetMaximum() >= hGamE->GetMaximum() && 
            hNeHadE->GetMaximum() >= hChHadE->GetMaximum() && hNeHadE->GetMaximum() >= hEleE->GetMaximum())
      haxisE->SetMaximum(hNeHadE->GetMaximum());	
    else if(hChHadE->GetMaximum() >= hPi0E->GetMaximum() && hChHadE->GetMaximum() >= hNeHadE->GetMaximum() && 
            hChHadE->GetMaximum() >= hGamE->GetMaximum() && hChHadE->GetMaximum() >= hEleE->GetMaximum())
      haxisE->SetMaximum(hChHadE->GetMaximum());
    else if(hEleE->GetMaximum() >= hPi0E->GetMaximum() && hEleE->GetMaximum() >= hNeHadE->GetMaximum() && 
            hEleE->GetMaximum() >= hChHadE->GetMaximum() && hEleE->GetMaximum() >= hGamE->GetMaximum())
      haxisE->SetMaximum(hEleE->GetMaximum());
    haxisE->SetXTitle("E (GeV)");
    haxisE->SetMinimum(1);
    haxisE->Draw("axis");
    hPi0E->Draw("same");
    hGamE->Draw("same");
    hNeHadE->Draw("same");
    hChHadE->Draw("same");
    hEleE->Draw("same");
    
    TLegend pLegend2(0.8,0.65,0.95,0.93);
    pLegend2.SetTextSize(0.06);
    pLegend2.AddEntry(hPi0E,"  #pi^{0}","L");
    pLegend2.AddEntry(hGamE,"  #gamma","L");
    pLegend2.AddEntry(hEleE,"  e^{#pm}","L");
    pLegend2.AddEntry(hChHadE,"  h^{#pm}","L");
    pLegend2.AddEntry(hNeHadE,"  h^{0}","L");
    pLegend2.SetFillColor(10);
    pLegend2.SetBorderSize(1);
    pLegend2.Draw();
    
    
    c11->cd(2) ; 
    gPad->SetLogy();
    //printf("%s, %s, %s, %s, %s\n",fhGamPt->GetName(),fhPi0Pt->GetName(),fhElePt->GetName(),fhNeHadPt->GetName(), fhChHadPt->GetName());
    TH1F * hGamPt   = (TH1F*) fhGamPt->ProjectionX(Form("%s_px",fhGamPt->GetName()),-1,-1);
    TH1F * hPi0Pt   = (TH1F*) fhPi0Pt->ProjectionX(Form("%s_px",fhPi0Pt->GetName()),-1,-1);
    TH1F * hElePt   = (TH1F*) fhElePt->ProjectionX(Form("%s_px",fhElePt->GetName()),-1,-1);
    TH1F * hNeHadPt = (TH1F*) fhNeHadPt->ProjectionX(Form("%s_px",fhNeHadPt->GetName()),-1,-1);
    TH1F * hChHadPt = (TH1F*) fhChHadPt->ProjectionX(Form("%s_px",fhChHadPt->GetName()),-1,-1);
    haxispt  = (TH1F*) hPi0Pt->Clone(Form("%s_axisPt",fhPi0Pt->GetName()));  
    haxispt->SetTitle("Reconstructed particles p_{T}, function of their original particle ID");
    hPi0Pt->SetLineColor(1);
    hGamPt->SetLineColor(4);
    hNeHadPt->SetLineColor(2);
    hChHadPt->SetLineColor(7);
    hElePt->SetLineColor(6);
    
    //Select the maximum of the histogram to show all lines.
    if(hPi0Pt->GetMaximum() >= hGamPt->GetMaximum() && hPi0Pt->GetMaximum() >= hNeHadPt->GetMaximum() && 
       hPi0Pt->GetMaximum() >= hChHadPt->GetMaximum() && hPi0Pt->GetMaximum() >= hElePt->GetMaximum())
      haxispt->SetMaximum(hPi0Pt->GetMaximum());
    else if(hGamPt->GetMaximum() >= hPi0Pt->GetMaximum() && hGamPt->GetMaximum() >= hNeHadPt->GetMaximum() && 
            hGamPt->GetMaximum() >= hChHadPt->GetMaximum() && hGamPt->GetMaximum() >= hElePt->GetMaximum())
      haxispt->SetMaximum(hGamPt->GetMaximum());
    else if(hNeHadPt->GetMaximum() >= hPi0Pt->GetMaximum() && hNeHadPt->GetMaximum() >= hGamPt->GetMaximum() && 
            hNeHadPt->GetMaximum() >= hChHadPt->GetMaximum() && hNeHadPt->GetMaximum() >= hElePt->GetMaximum())
      haxispt->SetMaximum(hNeHadPt->GetMaximum());	
    else if(hChHadPt->GetMaximum() >= hPi0Pt->GetMaximum() && hChHadPt->GetMaximum() >= hNeHadPt->GetMaximum() && 
            hChHadPt->GetMaximum() >= hGamPt->GetMaximum() && hChHadPt->GetMaximum() >= hElePt->GetMaximum())
      haxispt->SetMaximum(hChHadPt->GetMaximum());
    else if(hElePt->GetMaximum() >= hPi0Pt->GetMaximum() && hElePt->GetMaximum() >= hNeHadPt->GetMaximum() && 
            hElePt->GetMaximum() >= hChHadPt->GetMaximum() && hElePt->GetMaximum() >= hGamPt->GetMaximum())
      haxispt->SetMaximum(hElePt->GetMaximum());
    haxispt->SetXTitle("p_{T} (GeV/c)");
    haxispt->SetMinimum(1);
    haxispt->Draw("axis");
    hPi0Pt->Draw("same");
    hGamPt->Draw("same");
    hNeHadPt->Draw("same");
    hChHadPt->Draw("same");
    hElePt->Draw("same");
    
    c11->cd(3) ;
    gPad->SetLogy();
    
    TH1F * hGamEta   = (TH1F*) fhGamEta->ProjectionX(Form("%s_px",fhGamEta->GetName()),-1,-1);
    TH1F * hPi0Eta   = (TH1F*) fhPi0Eta->ProjectionX(Form("%s_px",fhPi0Eta->GetName()),-1,-1);
    TH1F * hEleEta   = (TH1F*) fhEleEta->ProjectionX(Form("%s_px",fhEleEta->GetName()),-1,-1);
    TH1F * hNeHadEta = (TH1F*) fhNeHadEta->ProjectionX(Form("%s_px",fhNeHadEta->GetName()),-1,-1);
    TH1F * hChHadEta = (TH1F*) fhChHadEta->ProjectionX(Form("%s_px",fhChHadEta->GetName()),-1,-1);
    haxiseta  = (TH1F*) hPi0Eta->Clone(Form("%s_axisEta",fhPi0Eta->GetName()));  
    haxiseta->SetTitle("Reconstructed particles #eta, function of their original particle ID");
    hPi0Eta->SetLineColor(1);
    hGamEta->SetLineColor(4);
    hNeHadEta->SetLineColor(2);
    hChHadEta->SetLineColor(7);
    hEleEta->SetLineColor(6);
    //Select the maximum of the histogram to show all lines.
    if(hPi0Eta->GetMaximum() >= hGamEta->GetMaximum() && hPi0Eta->GetMaximum() >= hNeHadEta->GetMaximum() && 
       hPi0Eta->GetMaximum() >= hChHadEta->GetMaximum() && hPi0Eta->GetMaximum() >= hEleEta->GetMaximum())
      haxiseta->SetMaximum(hPi0Eta->GetMaximum());
    else if(hGamEta->GetMaximum() >= hPi0Eta->GetMaximum() && hGamEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
            hGamEta->GetMaximum() >= hChHadEta->GetMaximum() && hGamEta->GetMaximum() >= hEleEta->GetMaximum())
      haxiseta->SetMaximum(hGamEta->GetMaximum());
    else if(hNeHadEta->GetMaximum() >= hPi0Eta->GetMaximum() && hNeHadEta->GetMaximum() >= hGamEta->GetMaximum() && 
            hNeHadEta->GetMaximum() >= hChHadEta->GetMaximum() && hNeHadEta->GetMaximum() >= hEleEta->GetMaximum())
      haxiseta->SetMaximum(hNeHadEta->GetMaximum());	
    else if(hChHadEta->GetMaximum() >= hPi0Eta->GetMaximum() && hChHadEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
            hChHadEta->GetMaximum() >= hGamEta->GetMaximum() && hChHadEta->GetMaximum() >= hEleEta->GetMaximum())
      haxiseta->SetMaximum(hChHadEta->GetMaximum());
    else if(hEleEta->GetMaximum() >= hPi0Eta->GetMaximum() && hEleEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
            hEleEta->GetMaximum() >= hChHadEta->GetMaximum() && hEleEta->GetMaximum() >= hGamEta->GetMaximum())
      haxiseta->SetMaximum(hEleEta->GetMaximum());
    
    haxiseta->SetXTitle("#eta");
    haxiseta->Draw("axis");
    hPi0Eta->Draw("same");
    hGamEta->Draw("same");
    hNeHadEta->Draw("same");
    hChHadEta->Draw("same");
    hEleEta->Draw("same");
    
    
    c11->cd(4) ; 
    gPad->SetLogy();
    TH1F * hGamPhi   = (TH1F*) fhGamPhi->ProjectionX(Form("%s_px",fhGamPhi->GetName()),-1,-1);
    TH1F * hPi0Phi   = (TH1F*) fhPi0Phi->ProjectionX(Form("%s_px",fhPi0Phi->GetName()),-1,-1);
    TH1F * hElePhi   = (TH1F*) fhElePhi->ProjectionX(Form("%s_px",fhElePhi->GetName()),-1,-1);
    TH1F * hNeHadPhi = (TH1F*) fhNeHadPhi->ProjectionX(Form("%s_px",fhNeHadPhi->GetName()),-1,-1);
    TH1F * hChHadPhi = (TH1F*) fhChHadPhi->ProjectionX(Form("%s_px",fhChHadPhi->GetName()),-1,-1);
    haxisphi  = (TH1F*) hPi0Phi->Clone(Form("%s_axisPhi",fhPi0Phi->GetName()));  
    haxisphi->SetTitle("Reconstructed particles #phi, function of their original particle ID");
    
    hPi0Phi->SetLineColor(1);
    hGamPhi->SetLineColor(4);
    hNeHadPhi->SetLineColor(2);
    hChHadPhi->SetLineColor(7);
    hElePhi->SetLineColor(6);
    //Select the maximum of the histogram to show all lines.
    if(hPi0Phi->GetMaximum() >= hGamPhi->GetMaximum() && hPi0Phi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
       hPi0Phi->GetMaximum() >= hChHadPhi->GetMaximum() && hPi0Phi->GetMaximum() >= hElePhi->GetMaximum())
      haxisphi->SetMaximum(hPi0Phi->GetMaximum());
    else if(hGamPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hGamPhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
            hGamPhi->GetMaximum() >= hChHadPhi->GetMaximum() && hGamPhi->GetMaximum() >= hElePhi->GetMaximum())
      haxisphi->SetMaximum(hGamPhi->GetMaximum());
    else if(hNeHadPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hNeHadPhi->GetMaximum() >= hGamPhi->GetMaximum() && 
            hNeHadPhi->GetMaximum() >= hChHadPhi->GetMaximum() && hNeHadPhi->GetMaximum() >= hElePhi->GetMaximum())
      haxisphi->SetMaximum(hNeHadPhi->GetMaximum());	
    else if(hChHadPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hChHadPhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
            hChHadPhi->GetMaximum() >= hGamPhi->GetMaximum() && hChHadPhi->GetMaximum() >= hElePhi->GetMaximum())
      haxisphi->SetMaximum(hChHadPhi->GetMaximum());
    else if(hElePhi->GetMaximum() >= hPi0Phi->GetMaximum() && hElePhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
            hElePhi->GetMaximum() >= hChHadPhi->GetMaximum() && hElePhi->GetMaximum() >= hGamPhi->GetMaximum())
      haxisphi->SetMaximum(hElePhi->GetMaximum());
    haxisphi->SetXTitle("#phi (rad)");
    haxisphi->Draw("axis");
    hPi0Phi->Draw("same");
    hGamPhi->Draw("same");
    hNeHadPhi->Draw("same");
    hChHadPhi->Draw("same");
    hElePhi->Draw("same");
    
    snprintf(name,buffersize,"QA_%s_RecDistributionsGenID.eps",fCalorimeter.Data());
    c11->Print(name); printf("Plot: %s\n",name);
    
    
    //Ratio reconstructed clusters / generated particles in acceptance, for different particle ID
    //printf("c1\n");
    
    TH1F *	hPi0EClone   = (TH1F*)   hPi0E  ->Clone(Form("%s_Clone",fhPi0E->GetName()));
    TH1F *	hGamEClone   = (TH1F*)   hGamE  ->Clone(Form("%s_Clone",fhGamE->GetName()));
    TH1F *	hPi0PtClone  = (TH1F*)   hPi0Pt ->Clone(Form("%s_Clone",fhPi0Pt->GetName()));
    TH1F *	hGamPtClone  = (TH1F*)   hGamPt ->Clone(Form("%s_Clone",fhGamPt->GetName()));	
    TH1F *	hPi0EtaClone = (TH1F*)   hPi0Eta->Clone(Form("%s_Clone",fhPi0Eta->GetName()));
    TH1F *	hGamEtaClone = (TH1F*)   hGamEta->Clone(Form("%s_Clone",fhGamEta->GetName()));	
    TH1F *	hPi0PhiClone = (TH1F*)   hPi0Phi->Clone(Form("%s_Clone",fhPi0Phi->GetName()));
    TH1F *	hGamPhiClone = (TH1F*)   hGamPhi->Clone(Form("%s_Clone",fhGamPhi->GetName()));	
    
    snprintf(cname,buffersize,"QA_%s_recgenidratio",fCalorimeter.Data());
    TCanvas  * c12 = new TCanvas(cname, "Ratio reconstructed clusters / generated particles in acceptance, for different particle ID", 400, 400) ;
    c12->Divide(2, 2);
    
    c12->cd(1) ; 
    gPad->SetLogy();
    haxisE->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
    hPi0EClone->Divide(fhGenPi0AccE);
    hGamEClone->Divide(fhGenGamAccE);
    haxisE->SetMaximum(5);
    haxisE->SetMinimum(1e-2);
    haxisE->SetXTitle("E (GeV)");
    haxisE->SetYTitle("ratio = rec/gen");
    haxisE->Draw("axis");
    hPi0E->Draw("same");
    hGamE->Draw("same");
    
    TLegend pLegend3(0.75,0.2,0.9,0.4);
    pLegend3.SetTextSize(0.06);
    pLegend3.AddEntry(hPi0EClone,"  #pi^{0}","L");
    pLegend3.AddEntry(hGamEClone,"  #gamma","L");
    pLegend3.SetFillColor(10);
    pLegend3.SetBorderSize(1);
    pLegend3.Draw();
    
    c12->cd(2) ; 
    gPad->SetLogy();
    haxispt->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
    hPi0PtClone->Divide(fhGenPi0AccPt);
    hGamPtClone->Divide(fhGenGamAccPt);
    haxispt->SetMaximum(5);
    haxispt->SetMinimum(1e-2);
    haxispt->SetXTitle("p_{T} (GeV/c)");
    haxispt->SetYTitle("ratio = rec/gen");
    haxispt->Draw("axis");
    hPi0PtClone->Draw("same");
    hGamPtClone->Draw("same");
    
    c12->cd(3) ;
    gPad->SetLogy();
    
    haxiseta->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
    hPi0EtaClone->Divide(fhGenPi0AccEta);
    hGamEtaClone->Divide(fhGenGamAccEta);
    haxiseta->SetMaximum(1.2);
    haxiseta->SetMinimum(1e-2);
    haxiseta->SetYTitle("ratio = rec/gen");
    haxiseta->SetXTitle("#eta");
    haxiseta->Draw("axis");
    hPi0EtaClone->Draw("same");
    hGamEtaClone->Draw("same");
    
    
    c12->cd(4) ; 
    gPad->SetLogy();
    haxisphi->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
    hPi0PhiClone->Divide(fhGenPi0AccPhi);
    hGamPhiClone->Divide(fhGenGamAccPhi);
    haxisphi->SetYTitle("ratio = rec/gen");
    haxisphi->SetXTitle("#phi (rad)");
    haxisphi->SetMaximum(1.2);
    haxisphi->SetMinimum(1e-2);
    haxisphi->Draw("axis");
    hPi0PhiClone->Draw("same");
    hGamPhiClone->Draw("same");
    
    snprintf(name,buffersize,"QA_%s_EfficiencyGenID.eps",fCalorimeter.Data());
    c12->Print(name); printf("Plot: %s\n",name);
    
    
    
    //Reconstructed distributions
    //printf("c1\n");
    snprintf(cname,buffersize,"QA_%s_vertex",fCalorimeter.Data());
    TCanvas  * c13 = new TCanvas(cname, "Particle vertex", 400, 400) ;
    c13->Divide(2, 2);
    
    c13->cd(1) ; 
    //gPad->SetLogy();
    fhEMVxyz->SetTitleOffset(1.6,"Y");
    fhEMVxyz->Draw();
    
    c13->cd(2) ; 
    //gPad->SetLogy();
    fhHaVxyz->SetTitleOffset(1.6,"Y");
    fhHaVxyz->Draw();
    
    c13->cd(3) ;
    gPad->SetLogy();
    TH1F * hEMR = (TH1F*) fhEMR->ProjectionY(Form("%s_py",fhEMR->GetName()),-1,-1); 
    hEMR->SetLineColor(4);
    hEMR->Draw();
    
    c13->cd(4) ; 
    gPad->SetLogy();
    TH1F * hHaR = (TH1F*) fhHaR->ProjectionY(Form("%s_py",fhHaR->GetName()),-1,-1); 
    hHaR->SetLineColor(4);
    hHaR->Draw();
    
    
    snprintf(name,buffersize,"QA_%s_ParticleVertex.eps",fCalorimeter.Data());
    c13->Print(name); printf("Plot: %s\n",name);
    
    
    //Track-matching distributions
    if(fFillAllTH12){
      //Reconstructed distributions, matched with tracks, generated particle dependence
      //printf("c2\n");
      snprintf(cname,buffersize,"QA_%s_rectrackmatchGenID",fCalorimeter.Data());
      TCanvas  * c22ch = new TCanvas(cname, "Reconstructed distributions, matched with tracks, for different particle ID", 400, 400) ;
      c22ch->Divide(2, 2);
      
      c22ch->cd(1) ; 
      
      TH1F * hGamECharged   = (TH1F*) fhGamECharged->ProjectionX(Form("%s_px",fhGamECharged->GetName()),-1,-1);
      TH1F * hPi0ECharged   = (TH1F*) fhPi0ECharged->ProjectionX(Form("%s_px",fhPi0ECharged->GetName()),-1,-1);
      TH1F * hEleECharged   = (TH1F*) fhEleECharged->ProjectionX(Form("%s_px",fhEleECharged->GetName()),-1,-1);
      TH1F * hNeHadECharged = (TH1F*) fhNeHadECharged->ProjectionX(Form("%s_px",fhNeHadECharged->GetName()),-1,-1);
      TH1F * hChHadECharged = (TH1F*) fhChHadECharged->ProjectionX(Form("%s_px",fhChHadECharged->GetName()),-1,-1);
      hPi0ECharged->SetLineColor(1);
      hGamECharged->SetLineColor(4);
      hNeHadECharged->SetLineColor(2);
      hChHadECharged->SetLineColor(7);
      hEleECharged->SetLineColor(6);	
      gPad->SetLogy();
      fhECharged->SetLineColor(3);
      fhECharged->SetMinimum(0.5);
      fhECharged->Draw();
      hPi0ECharged->Draw("same");
      hGamECharged->Draw("same");
      hNeHadECharged->Draw("same");
      hChHadECharged->Draw("same");
      hEleECharged->Draw("same");
      TLegend pLegend22(0.75,0.45,0.9,0.8);
      pLegend22.SetTextSize(0.06);
      pLegend22.AddEntry(fhECharged,"all","L");
      pLegend22.AddEntry(hPi0ECharged,"#pi^{0}","L");
      pLegend22.AddEntry(hGamECharged,"#gamma","L");
      pLegend22.AddEntry(hEleECharged,"e^{#pm}","L");
      pLegend22.AddEntry(hChHadECharged,"h^{#pm}","L");
      pLegend22.AddEntry(hNeHadECharged,"h^{0}","L");
      pLegend22.SetFillColor(10);
      pLegend22.SetBorderSize(1);
      pLegend22.Draw();
      
      c22ch->cd(2) ; 
      
      TH1F * hGamPtCharged   = (TH1F*) fhGamPtCharged->ProjectionX(Form("%s_px",fhGamPtCharged->GetName()),-1,-1);
      TH1F * hPi0PtCharged   = (TH1F*) fhPi0PtCharged->ProjectionX(Form("%s_px",fhPi0PtCharged->GetName()),-1,-1);
      TH1F * hElePtCharged   = (TH1F*) fhElePtCharged->ProjectionX(Form("%s_px",fhElePtCharged->GetName()),-1,-1);
      TH1F * hNeHadPtCharged = (TH1F*) fhNeHadPtCharged->ProjectionX(Form("%s_px",fhNeHadPtCharged->GetName()),-1,-1);
      TH1F * hChHadPtCharged = (TH1F*) fhChHadPtCharged->ProjectionX(Form("%s_px",fhChHadPtCharged->GetName()),-1,-1);
      hPi0PtCharged->SetLineColor(1);
      hGamPtCharged->SetLineColor(4);
      hNeHadPtCharged->SetLineColor(2);
      hChHadPtCharged->SetLineColor(7);
      hElePtCharged->SetLineColor(6);	
      gPad->SetLogy();
      fhPtCharged->SetLineColor(3);
      fhPtCharged->SetMinimum(0.5);
      fhPtCharged->Draw();
      hPi0PtCharged->Draw("same");
      hGamPtCharged->Draw("same");
      hNeHadPtCharged->Draw("same");
      hChHadPtCharged->Draw("same");
      hElePtCharged->Draw("same");	
      
      c22ch->cd(4) ; 
      
      TH1F * hGamEtaCharged   = (TH1F*) fhGamEtaCharged->ProjectionX(Form("%s_px",fhGamEtaCharged->GetName()),-1,-1);
      TH1F * hPi0EtaCharged   = (TH1F*) fhPi0EtaCharged->ProjectionX(Form("%s_px",fhPi0EtaCharged->GetName()),-1,-1);
      TH1F * hEleEtaCharged   = (TH1F*) fhEleEtaCharged->ProjectionX(Form("%s_px",fhEleEtaCharged->GetName()),-1,-1);
      TH1F * hNeHadEtaCharged = (TH1F*) fhNeHadEtaCharged->ProjectionX(Form("%s_px",fhNeHadEtaCharged->GetName()),-1,-1);
      TH1F * hChHadEtaCharged = (TH1F*) fhChHadEtaCharged->ProjectionX(Form("%s_px",fhChHadEtaCharged->GetName()),-1,-1);
      hPi0EtaCharged->SetLineColor(1);
      hGamEtaCharged->SetLineColor(4);
      hNeHadEtaCharged->SetLineColor(2);
      hChHadEtaCharged->SetLineColor(7);
      hEleEtaCharged->SetLineColor(6);	
      gPad->SetLogy();
      fhEtaCharged->SetLineColor(3);
      fhEtaCharged->SetMinimum(0.5);
      fhEtaCharged->Draw();
      hPi0EtaCharged->Draw("same");
      hGamEtaCharged->Draw("same");
      hNeHadEtaCharged->Draw("same");
      hChHadEtaCharged->Draw("same");
      hEleEtaCharged->Draw("same");
      
      c22ch->cd(3) ; 
      
      TH1F * hGamPhiCharged   = (TH1F*) fhGamPhiCharged->ProjectionX(Form("%s_px",fhGamPhiCharged->GetName()),-1,-1);
      TH1F * hPi0PhiCharged   = (TH1F*) fhPi0PhiCharged->ProjectionX(Form("%s_px",fhPi0PhiCharged->GetName()),-1,-1);
      TH1F * hElePhiCharged   = (TH1F*) fhElePhiCharged->ProjectionX(Form("%s_px",fhElePhiCharged->GetName()),-1,-1);
      TH1F * hNeHadPhiCharged = (TH1F*) fhNeHadPhiCharged->ProjectionX(Form("%s_px",fhNeHadPhiCharged->GetName()),-1,-1);
      TH1F * hChHadPhiCharged = (TH1F*) fhChHadPhiCharged->ProjectionX(Form("%s_px",fhChHadPhiCharged->GetName()),-1,-1);
      hPi0PhiCharged->SetLineColor(1);
      hGamPhiCharged->SetLineColor(4);
      hNeHadPhiCharged->SetLineColor(2);
      hChHadPhiCharged->SetLineColor(7);
      hElePhiCharged->SetLineColor(6);	
      gPad->SetLogy();
      fhPhiCharged->SetLineColor(3);
      fhPhiCharged->SetMinimum(0.5);
      fhPhiCharged->Draw();
      hPi0PhiCharged->Draw("same");
      hGamPhiCharged->Draw("same");
      hNeHadPhiCharged->Draw("same");
      hChHadPhiCharged->Draw("same");
      hElePhiCharged->Draw("same");
      
      
      snprintf(name,buffersize,"QA_%s_ReconstructedDistributions_TrackMatchedGenID.eps",fCalorimeter.Data());
      c22ch->Print(name); printf("Plot: %s\n",name);
      
      TH1F *	hGamEChargedClone   = (TH1F*)   hGamECharged->Clone(Form("%s_Clone",fhGamECharged->GetName()));
      TH1F *	hGamPtChargedClone  = (TH1F*)   hGamPtCharged->Clone(Form("%s_Clone",fhGamPtCharged->GetName()));
      TH1F *	hGamEtaChargedClone = (TH1F*)   hGamEtaCharged->Clone(Form("%s_Clone",fhGamEtaCharged->GetName()));
      TH1F *	hGamPhiChargedClone = (TH1F*)   hGamPhiCharged->Clone(Form("%s_Clone",fhGamPhiCharged->GetName()));
      
      TH1F *	hPi0EChargedClone   = (TH1F*)   hPi0ECharged->Clone(Form("%s_Clone",fhPi0ECharged->GetName()));
      TH1F *	hPi0PtChargedClone  = (TH1F*)   hPi0PtCharged->Clone(Form("%s_Clone",fhPi0PtCharged->GetName()));
      TH1F *	hPi0EtaChargedClone = (TH1F*)   hPi0EtaCharged->Clone(Form("%s_Clone",fhPi0EtaCharged->GetName()));
      TH1F *	hPi0PhiChargedClone = (TH1F*)   hPi0PhiCharged->Clone(Form("%s_Clone",fhPi0PhiCharged->GetName()));
      
      TH1F *	hEleEChargedClone   = (TH1F*)   hEleECharged->Clone(Form("%s_Clone",fhEleECharged->GetName()));
      TH1F *	hElePtChargedClone  = (TH1F*)   hElePtCharged->Clone(Form("%s_Clone",fhElePtCharged->GetName()));
      TH1F *	hEleEtaChargedClone = (TH1F*)   hEleEtaCharged->Clone(Form("%s_Clone",fhEleEtaCharged->GetName()));
      TH1F *	hElePhiChargedClone = (TH1F*)   hElePhiCharged->Clone(Form("%s_Clone",fhElePhiCharged->GetName()));	
      
      TH1F *	hNeHadEChargedClone   = (TH1F*)   hNeHadECharged->Clone(Form("%s_Clone",fhNeHadECharged->GetName()));
      TH1F *	hNeHadPtChargedClone  = (TH1F*)   hNeHadPtCharged->Clone(Form("%s_Clone",fhNeHadPtCharged->GetName()));
      TH1F *	hNeHadEtaChargedClone = (TH1F*)   hNeHadEtaCharged->Clone(Form("%s_Clone",fhNeHadEtaCharged->GetName()));
      TH1F *	hNeHadPhiChargedClone = (TH1F*)   hNeHadPhiCharged->Clone(Form("%s_Clone",fhNeHadPhiCharged->GetName()));
      
      TH1F *	hChHadEChargedClone   = (TH1F*)   hChHadECharged->Clone(Form("%s_Clone",fhChHadECharged->GetName()));
      TH1F *	hChHadPtChargedClone  = (TH1F*)   hChHadPtCharged->Clone(Form("%s_Clone",fhChHadPtCharged->GetName()));
      TH1F *	hChHadEtaChargedClone = (TH1F*)   hChHadEtaCharged->Clone(Form("%s_Clone",fhChHadEtaCharged->GetName()));
      TH1F *	hChHadPhiChargedClone = (TH1F*)   hChHadPhiCharged->Clone(Form("%s_Clone",fhChHadPhiCharged->GetName()));	
      
      //Ratio: reconstructed track matched/ all reconstructed
      //printf("c3\n");
      snprintf(cname,buffersize,"QA_%s_rectrackmatchratGenID",fCalorimeter.Data());
      TCanvas  * c3ch = new TCanvas(cname, "Ratio: reconstructed track matched/ all reconstructed, for different particle ID", 400, 400) ;
      c3ch->Divide(2, 2);
      
      c3ch->cd(1) ;
      hEChargedClone->SetMaximum(1.2);
      hEChargedClone->SetMinimum(0.001);	
      hEChargedClone->SetLineColor(3);
      hEChargedClone->SetYTitle("track matched / all");
      hPi0EChargedClone->Divide(hPi0E);
      hGamEChargedClone->Divide(hGamE);
      hEleEChargedClone->Divide(hEleE);
      hNeHadEChargedClone->Divide(hNeHadE);
      hChHadEChargedClone->Divide(hChHadE);
      hEChargedClone->Draw();
      hPi0EChargedClone->Draw("same");
      hGamEChargedClone->Draw("same");
      hEleEChargedClone->Draw("same");
      hNeHadEChargedClone->Draw("same");
      hChHadEChargedClone->Draw("same");
      
      TLegend pLegend3ch(0.75,0.45,0.9,0.8);
      pLegend3ch.SetTextSize(0.06);
      pLegend3ch.AddEntry(hEChargedClone,"all","L");
      pLegend3ch.AddEntry(hPi0EChargedClone,"#pi^{0}","L");
      pLegend3ch.AddEntry(hGamEChargedClone,"#gamma","L");
      pLegend3ch.AddEntry(hEleEChargedClone,"e^{#pm}","L");
      pLegend3ch.AddEntry(hChHadEChargedClone,"h^{#pm}","L");
      pLegend3ch.AddEntry(hNeHadEChargedClone,"h^{0}","L");
      pLegend3ch.SetFillColor(10);
      pLegend3ch.SetBorderSize(1);
      pLegend3ch.Draw();
      
      c3ch->cd(2) ;
      hPtChargedClone->SetMaximum(1.2);
      hPtChargedClone->SetMinimum(0.001);	
      hPtChargedClone->SetLineColor(3);
      hPtChargedClone->SetYTitle("track matched / all");
      hPi0PtChargedClone->Divide(hPi0Pt);
      hGamPtChargedClone->Divide(hGamPt);
      hElePtChargedClone->Divide(hElePt);
      hNeHadPtChargedClone->Divide(hNeHadPt);
      hChHadPtChargedClone->Divide(hChHadPt);
      hPtChargedClone->Draw();
      hPi0PtChargedClone->Draw("same");
      hGamPtChargedClone->Draw("same");
      hElePtChargedClone->Draw("same");
      hNeHadPtChargedClone->Draw("same");
      hChHadPtChargedClone->Draw("same");
      
      c3ch->cd(4) ;
      hEtaChargedClone->SetMaximum(1.2);
      hEtaChargedClone->SetMinimum(0.001);	
      hEtaChargedClone->SetLineColor(3);
      hEtaChargedClone->SetYTitle("track matched / all");
      hPi0EtaChargedClone->Divide(hPi0Eta);
      hGamEtaChargedClone->Divide(hGamEta);
      hEleEtaChargedClone->Divide(hEleEta);
      hNeHadEtaChargedClone->Divide(hNeHadEta);
      hChHadEtaChargedClone->Divide(hChHadEta);
      hEtaChargedClone->Draw();
      hPi0EtaChargedClone->Draw("same");
      hGamEtaChargedClone->Draw("same");
      hEleEtaChargedClone->Draw("same");
      hNeHadEtaChargedClone->Draw("same");
      hChHadEtaChargedClone->Draw("same");
      
      c3ch->cd(3) ;
      hPhiChargedClone->SetMaximum(1.2);
      hPhiChargedClone->SetMinimum(0.001);
      hPhiChargedClone->SetLineColor(3);
      hPhiChargedClone->SetYTitle("track matched / all");
      hPi0PhiChargedClone->Divide(hPi0Phi);
      hGamPhiChargedClone->Divide(hGamPhi);
      hElePhiChargedClone->Divide(hElePhi);
      hNeHadPhiChargedClone->Divide(hNeHadPhi);
      hChHadPhiChargedClone->Divide(hChHadPhi);
      hPhiChargedClone->Draw();
      hPi0PhiChargedClone->Draw("same");
      hGamPhiChargedClone->Draw("same");
      hElePhiChargedClone->Draw("same");
      hNeHadPhiChargedClone->Draw("same");
      hChHadPhiChargedClone->Draw("same");
      
      snprintf(name,buffersize,"QA_%s_RatioReconstructedMatchedDistributionsGenID.eps",fCalorimeter.Data());
      c3ch->Print(name); printf("Plot: %s\n",name);
      
    }	
  }
  //Track-matching distributions
  
  snprintf(cname,buffersize,"QA_%s_trkmatch",fCalorimeter.Data());
  TCanvas *cme = new TCanvas(cname,"Track-matching distributions", 400, 400);
  cme->Divide(2,2);
  
  TLegend pLegendpE0(0.6,0.55,0.9,0.8);
  pLegendpE0.SetTextSize(0.04);
  pLegendpE0.AddEntry(fh1pOverE,"all","L");
  pLegendpE0.AddEntry(fh1pOverER02,"dR < 0.02","L");		
  pLegendpE0.SetFillColor(10);
  pLegendpE0.SetBorderSize(1);
  //pLegendpE0.Draw();
  
  cme->cd(1);
  if(fh1pOverE->GetEntries() > 0) gPad->SetLogy();
  fh1pOverE->SetTitle("Track matches p/E");
  fh1pOverE->Draw();
  fh1pOverER02->SetLineColor(4);
  fh1pOverER02->Draw("same");
  pLegendpE0.Draw();
  
  cme->cd(2);
  if(fh1dR->GetEntries() > 0) gPad->SetLogy();
  fh1dR->Draw();
  
  cme->cd(3);
  fh2MatchdEdx->Draw();
  
  cme->cd(4);
  fh2EledEdx->Draw();
  
  snprintf(name,buffersize,"QA_%s_TrackMatchingEleDist.eps",fCalorimeter.Data());
  cme->Print(name); printf("Plot: %s\n",name);       
  
  if(IsDataMC()){
    snprintf(cname,buffersize,"QA_%s_trkmatchMCEle",fCalorimeter.Data());
    TCanvas *cmemc = new TCanvas(cname,"Track-matching distributions from MC electrons", 600, 200);
    cmemc->Divide(3,1);
    
    cmemc->cd(1);
    gPad->SetLogy();
    fhMCEle1pOverE->Draw();
    fhMCEle1pOverER02->SetLineColor(4);
    fhMCEle1pOverE->SetLineColor(1);
    fhMCEle1pOverER02->Draw("same");
    pLegendpE0.Draw();
		
    cmemc->cd(2);
    gPad->SetLogy();
    fhMCEle1dR->Draw();
		
    cmemc->cd(3);
    fhMCEle2MatchdEdx->Draw();
		
    snprintf(name,buffersize,"QA_%s_TrackMatchingDistMCEle.eps",fCalorimeter.Data());
    cmemc->Print(name); printf("Plot: %s\n",name);  
    
		
    snprintf(cname,buffersize,"QA_%s_trkmatchMCChHad",fCalorimeter.Data());
    TCanvas *cmemchad = new TCanvas(cname,"Track-matching distributions from MC charged hadrons", 600, 200);
    cmemchad->Divide(3,1);
		
    cmemchad->cd(1);
    gPad->SetLogy();
    fhMCChHad1pOverE->Draw();
    fhMCChHad1pOverER02->SetLineColor(4);
    fhMCChHad1pOverE->SetLineColor(1);
    fhMCChHad1pOverER02->Draw("same");
    pLegendpE0.Draw();
		
    cmemchad->cd(2);
    gPad->SetLogy();
    fhMCChHad1dR->Draw();
    
    cmemchad->cd(3);
    fhMCChHad2MatchdEdx->Draw();
		
    snprintf(name,buffersize,"QA_%s_TrackMatchingDistMCChHad.eps",fCalorimeter.Data());
    cmemchad->Print(name); printf("Plot: %s\n",name);       
    
    snprintf(cname,buffersize,"QA_%s_trkmatchMCNeutral",fCalorimeter.Data());
    TCanvas *cmemcn = new TCanvas(cname,"Track-matching distributions from MC neutrals", 600, 200);
    cmemcn->Divide(3,1);
		
    cmemcn->cd(1);
    gPad->SetLogy();
    fhMCNeutral1pOverE->Draw();
    fhMCNeutral1pOverE->SetLineColor(1);
    fhMCNeutral1pOverER02->SetLineColor(4);
    fhMCNeutral1pOverER02->Draw("same");
    pLegendpE0.Draw();
		
    cmemcn->cd(2);
    gPad->SetLogy();
    fhMCNeutral1dR->Draw();
		
    cmemcn->cd(3);
    fhMCNeutral2MatchdEdx->Draw();
		
    snprintf(name,buffersize,"QA_%s_TrackMatchingDistMCNeutral.eps",fCalorimeter.Data());
    cmemcn->Print(name); printf("Plot: %s\n",name);       
    
    snprintf(cname,buffersize,"QA_%s_trkmatchpE",fCalorimeter.Data());
    TCanvas *cmpoe = new TCanvas(cname,"Track-matching distributions, p/E", 400, 200);
    cmpoe->Divide(2,1);
		
    cmpoe->cd(1);
    gPad->SetLogy();
    fh1pOverE->SetLineColor(1);
    fhMCEle1pOverE->SetLineColor(4);
    fhMCChHad1pOverE->SetLineColor(2);
    fhMCNeutral1pOverE->SetLineColor(7);
    fh1pOverER02->SetMinimum(0.5);
    fh1pOverE->Draw();
    fhMCEle1pOverE->Draw("same");
    fhMCChHad1pOverE->Draw("same");
    fhMCNeutral1pOverE->Draw("same");
    TLegend pLegendpE(0.65,0.55,0.9,0.8);
    pLegendpE.SetTextSize(0.06);
    pLegendpE.AddEntry(fh1pOverE,"all","L");
    pLegendpE.AddEntry(fhMCEle1pOverE,"e^{#pm}","L");
    pLegendpE.AddEntry(fhMCChHad1pOverE,"h^{#pm}","L");
    pLegendpE.AddEntry(fhMCNeutral1pOverE,"neutrals","L");
    pLegendpE.SetFillColor(10);
    pLegendpE.SetBorderSize(1);
    pLegendpE.Draw();
    
    cmpoe->cd(2);
    gPad->SetLogy();
    fh1pOverER02->SetTitle("Track matches p/E, dR<0.2");
    fh1pOverER02->SetLineColor(1);
    fhMCEle1pOverER02->SetLineColor(4);
    fhMCChHad1pOverER02->SetLineColor(2);
    fhMCNeutral1pOverER02->SetLineColor(7);
    fh1pOverER02->SetMaximum(fh1pOverE->GetMaximum());
    fh1pOverER02->SetMinimum(0.5);
    fh1pOverER02->Draw();
    fhMCEle1pOverER02->Draw("same");
    fhMCChHad1pOverER02->Draw("same");
    fhMCNeutral1pOverER02->Draw("same");
    
    //		TLegend pLegendpE2(0.65,0.55,0.9,0.8);
    //		pLegendpE2.SetTextSize(0.06);
    //		pLegendpE2.SetHeader("dR < 0.02");
    //		pLegendpE2.SetFillColor(10);
    //		pLegendpE2.SetBorderSize(1);
    //		pLegendpE2.Draw();
    
    snprintf(name,buffersize,"QA_%s_TrackMatchingPOverE.eps",fCalorimeter.Data());
    cmpoe->Print(name); printf("Plot: %s\n",name);       			
  }
  
  char line[buffersize] ; 
  snprintf(line, buffersize,".!tar -zcf QA_%s_%s.tar.gz *%s*.eps", fCalorimeter.Data(), GetName(),fCalorimeter.Data()) ; 
  gROOT->ProcessLine(line);
  snprintf(line, buffersize,".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
  
  printf("AliAnaCalorimeterQA::Terminate() - !! All the eps files are in QA_%s_%s.tar.gz !!!\n",  fCalorimeter.Data(), GetName());
  
}
