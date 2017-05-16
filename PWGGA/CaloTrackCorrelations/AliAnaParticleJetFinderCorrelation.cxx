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

// --- ROOT system ---
#include "TH2F.h"
#include "TClonesArray.h"
#include "TClass.h"
//#include "Riostream.h"

// ---- AliRoot system ----
#include "AliCaloTrackReader.h"
#include "AliAODJet.h"
#include "AliAnaParticleJetFinderCorrelation.h" 
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliVTrack.h"
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

// ---- Jets ----
#include "AliAODJetEventBackground.h"
#include "TRandom2.h"
//MC access
#include "AliMCAnalysisUtils.h"
#include "AliVParticle.h"

// --- Detectors ---
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaParticleJetFinderCorrelation) ;
/// \endcond

//________________________________________________________________________
/// Default constructor. Initialize parameters.
//________________________________________________________________________
AliAnaParticleJetFinderCorrelation::AliAnaParticleJetFinderCorrelation() :
AliAnaCaloTrackCorrBaseClass(),  
  fDeltaPhiMaxCut(0.), fDeltaPhiMinCut(0.), fRatioMaxCut(0.),  fRatioMinCut(0.), 
  fConeSize(0.), fPtThresholdInCone(0.),fUseJetRefTracks(kTRUE),
  fMakeCorrelationInHistoMaker(kFALSE), fSelectIsolated(kTRUE),
  fJetConeSize(0.4),fJetMinPt(5),fJetMinPtBkgSub(-100.),fJetAreaFraction(0.6),
//fNonStandardJetFromReader(kTRUE), 
  fJetBranchName("jets"),
  fBackgroundJetFromReader(kTRUE),
  fBkgJetBranchName("jets"),
  fGammaConeSize(0.3),fUseBackgroundSubtractionGamma(kFALSE),fSaveGJTree(kTRUE),
  fMostEnergetic(kFALSE),fMostOpposite(kTRUE), fUseHistogramJetBkg(kTRUE),
  fUseHistogramTracks(kTRUE),fUseHistogramJetTracks(kTRUE),fMCStudies(kFALSE),fGenerator(0),
  fMomentum(),
  fhDeltaEta(0), /*fhDeltaPhi(0),*/fhDeltaPhiCorrect(0),fhDeltaPhi0PiCorrect(0), fhDeltaPt(0), fhPtRatio(0), fhPt(0),
  fhFFz(0),fhFFxi(0),fhFFpt(0),fhNTracksInCone(0),
  fhJetFFz(0),fhJetFFxi(0),fhJetFFpt(0),fhJetFFzCor(0),fhJetFFxiCor(0),
  fhGamPtPerTrig(0),fhPtGamPtJet(0),
  fhBkgFFz(),fhBkgFFxi(),fhBkgFFpt(),fhBkgNTracksInCone(),fhBkgSumPtInCone(),fhBkgSumPtnTracksInCone(),
  fhNjetsNgammas(0),fhCuts(0),
  fhDeltaEtaBefore(0),fhDeltaPhiBefore(0),fhDeltaPtBefore(0),fhPtRatioBefore(0),
  fhPtBefore(0),fhDeltaPhi0PiCorrectBefore(0),
  fhJetPtBefore(0),fhJetPtBeforeCut(0),fhJetPt(0),fhJetPtMostEne(0),fhJetPhi(0),fhJetEta(0),fhJetEtaVsPt(0),
  fhJetPhiVsEta(0),fhJetEtaVsNpartInJet(0),fhJetEtaVsNpartInJetBkg(0),fhJetChBkgEnergyVsPt(0),fhJetChAreaVsPt(0),/*fhJetNjet(0),*/
  fhTrackPhiVsEta(0),fhTrackAveTrackPt(0),fhJetNjetOverPtCut(),
/*fhJetChBkgEnergyVsPtEtaGt05(0),fhJetChBkgEnergyVsPtEtaLe05(0),fhJetChAreaVsPtEtaGt05(0),fhJetChAreaVsPtEtaLe05(0),*/
  fhJetChBkgEnergyVsArea(0),fhJetRhoVsPt(0),fhJetRhoVsCentrality(0),//fhJetBkgRho(0),
  fhJetNparticlesInJet(0),fhJetDeltaEtaDeltaPhi(0),fhJetDeltaEtaDeltaPhiAllTracks(0),
  fhJetAveTrackPt(0),fhJetNtracksInJetAboveThr(),fhJetRatioNTrkAboveToNTrk(),fhJetNtrackRatioMostEne(),
  fhJetNtrackRatioJet5GeV(),fhJetNtrackRatioLead5GeV(),
  fhBkgJetBackground(),fhBkgJetSigma(),fhBkgJetArea(),fhPhotonPtMostEne(0),
  fhPhotonAverageEnergy(0),fhPhotonRatioAveEneToMostEne(0),fhPhotonAverageEnergyMinus1(0),fhPhotonRatioAveEneMinus1ToMostEne(0),
  fhPhotonNgammaMoreAverageToNgamma(0),fhPhotonNgammaMoreAverageMinus1ToNgamma(0),fhPhotonNgammaOverPtCut(),
  fhPhotonBkgRhoVsNtracks(0),fhPhotonBkgRhoVsNclusters(0),fhPhotonBkgRhoVsCentrality(0),
  fhPhotonBkgRhoVsNcells(0),fhPhotonPt(0),fhPhotonPtCorrected(0),fhPhotonPtCorrectedZoom(0),fhPhotonPtDiff(0),
  fhPhotonPtDiffVsCentrality(0),fhPhotonPtDiffVsNcells(0),fhPhotonPtDiffVsNtracks(0),fhPhotonPtDiffVsNclusters(0),
  fhPhotonSumPtInCone(0),fhPhotonSumPtCorrectInCone(0),fhPhotonSumPtChargedInCone(0),
  fhSelectedJetPhiVsEta(0),fhSelectedJetChBkgEnergyVsPtJet(0),fhSelectedJetChAreaVsPtJet(0),fhSelectedJetNjet(0),fhSelectedNtracks(0),
  fhSelectedTrackPhiVsEta(0),fhCuts2(0),
  fhSelectedPhotonNLMVsPt(0),fhSelectedPhotonLambda0VsPt(0), fhRandomPhiEta(),
  fhMCPhotonCuts(0),fhMCPhotonPt(0),fhMCPhotonEtaPhi(0),fhMCJetOrigin(0),
  fhMCJetNPartVsPt(0),fhMCJetChNPartVsPt(0),fhMCJetNPart150VsPt(0),fhMCJetChNPart150VsPt(0),fhMCJetChNPart150ConeVsPt(0),
  fhMCJetRatioChFull(0),fhMCJetRatioCh150Ch(0),
  fhMCJetEtaPhi(0),fhMCJetChEtaPhi(0),fhMCJet150EtaPhi(0),fhMCJetCh150EtaPhi(0),fhMCJetCh150ConeEtaPhi(0),
fTreeGJ     (0),
fGamPt	    (0),
fGamLambda0 (0),
fGamNLM	    (0),
fGamSumPtCh (0),
fGamTime    (0),
fGamNcells  (0),
fGamEta	    (0),
fGamPhi	    (0),
fGamSumPtNeu(0),
fGamNtracks (0),
fGamNclusters(0),
fGamAvEne   (0),
fJetPhi	    (0),
fJetEta	    (0),
fJetPt	    (0),
fJetBkgChEne(0),
fJetArea    (0),
fJetNtracks (0),
fJetNtracks1(0),
fJetNtracks2(0),
fJetRho(0),
fEventNumber(0),
fNtracks    (0),
fZvertex    (0),
fCentrality (0),
fIso(0),
fGamRho(0),
fMCGamPt        (0),
fMCGamEta	(0),
fMCGamPhi	(0),
fMCPartonType	(0),
fMCJetPt	(0),
fMCJetChPt	(0),
fMCJet150Pt	(0),
fMCJetCh150Pt	(0),
fMCJetNPart	(0),
fMCJetChNPart	(0),
fMCJet150NPart	(0),
fMCJetCh150NPart(0),
fMCJetEta	(0),
fMCJetPhi	(0),
fMCJetChEta	(0),
fMCJetChPhi	(0),
fMCJet150Eta	(0),
fMCJet150Phi	(0),
fMCJetCh150Eta	(0),
fMCJetCh150Phi  (0),
fMCJetCh150ConePt(0),
fMCJetCh150ConeNPart(0),
fMCJetCh150ConeEta(0),
fMCJetCh150ConePhi(0)
{
  InitParameters();
  for(Int_t i=0;i<10;i++)
  {
    fhJetNjetOverPtCut     [i] = 0;
    fhPhotonNgammaOverPtCut[i] = 0;
  }
    
  fGenerator = new TRandom2();
  fGenerator->SetSeed(0);
}

//___________________________________________________________________
/// Destructor.
//___________________________________________________________________
AliAnaParticleJetFinderCorrelation::~AliAnaParticleJetFinderCorrelation()
{
  delete fGenerator;
}

//___________________________________________________________________
/// Create histograms to be saved in output file and
/// store them in fOutputContainer.
//___________________________________________________________________
TList *  AliAnaParticleJetFinderCorrelation::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ParticleJetFinderHistos") ; 
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();
  //	Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins();
  //	Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();
  //	Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();
  //	Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  //	Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
//	Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();	
  
//  fhDeltaPhi  = new TH2F("DeltaPhi","#phi_{trigger} - #phi_{jet} vs p_{T trigger}",nptbins,ptmin,ptmax,100,-6,4); 
//  fhDeltaPhi->SetYTitle("#Delta #phi");
//  fhDeltaPhi->SetXTitle("p_{T trigger} (GeV/c)");
//  outputContainer->Add(fhDeltaPhi);

  fhDeltaPhiCorrect  = new TH2F("DeltaPhiCorrect","#phi_{trigger} - #phi_{jet} vs p_{T trigger}",nptbins,ptmin,ptmax,100,0,6.5); 
  fhDeltaPhiCorrect->SetYTitle("#Delta #phi");
  fhDeltaPhiCorrect->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhDeltaPhiCorrect);

  fhDeltaPhi0PiCorrect  = new TH2F("DeltaPhi0PiCorrect","#phi_{trigger} - #phi_{jet} (0,#pi) vs p_{T trigger}",nptbins,ptmin,ptmax,100,0,3.5); 
  fhDeltaPhi0PiCorrect->SetYTitle("#Delta #phi");
  fhDeltaPhi0PiCorrect->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhDeltaPhi0PiCorrect);


  fhDeltaEta  = new TH2F("DeltaEta","#eta_{trigger} - #eta_{jet} vs p_{T trigger}",nptbins,ptmin,ptmax,100,-2,2); 
  fhDeltaEta->SetYTitle("#Delta #eta");
  fhDeltaEta->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhDeltaEta);
  
  fhDeltaPt  = new TH2F("DeltaPt","p_{T trigger} - p_{T jet} vs p_{T trigger}",nptbins,ptmin,ptmax,150,-50,100); 
  fhDeltaPt->SetYTitle("#Delta p_{T}");
  fhDeltaPt->SetXTitle("p_{T trigger} (GeV/c)"); 
  outputContainer->Add(fhDeltaPt);
  
  fhPtRatio  = new TH2F("PtRatio","p_{T jet} / p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,200,0,2.); 
  fhPtRatio->SetYTitle("ratio");
  fhPtRatio->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhPtRatio);
  
  fhPt  = new TH2F("Pt","p_{T jet} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
  fhPt->SetYTitle("p_{T jet}(GeV/c)");
  fhPt->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhPt);
  
  fhFFz  = new TH2F("FFz","z = p_{T i charged}/p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,200,0.,2);  
  fhFFz->SetYTitle("z");
  fhFFz->SetXTitle("p_{T trigger}");
  outputContainer->Add(fhFFz) ;
	
  fhFFxi  = new TH2F("FFxi","#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger}", nptbins,ptmin,ptmax,100,0.,10.); 
  fhFFxi->SetYTitle("#xi");
  fhFFxi->SetXTitle("p_{T trigger}");
  outputContainer->Add(fhFFxi) ;
  
  fhFFpt  = new TH2F("FFpt","p_{T i charged} vs p_{T trigger}", nptbins,ptmin,ptmax,100,0.,50.); 
  fhFFpt->SetYTitle("p_{T charged hadron}");
  fhFFpt->SetXTitle("p_{T trigger}");
  outputContainer->Add(fhFFpt) ;
  
  fhNTracksInCone  = new TH2F("NTracksInCone","Number of tracks in cone vs p_{T trigger}", nptbins,ptmin,ptmax,100,0.,150.); 
  fhNTracksInCone->SetYTitle("p_{T charged hadron}");
  fhNTracksInCone->SetXTitle("p_{T trigger}");
  outputContainer->Add(fhNTracksInCone) ;
  
  //FF according to jet axis
  fhJetFFz  = new TH2F("JetFFz","z = p_{T i charged}/p_{T jet} vs p_{T jet}",nptbins,ptmin,ptmax,200,0.,2);  
  fhJetFFz->SetYTitle("z");
  fhJetFFz->SetXTitle("p_{T jet}");
  outputContainer->Add(fhJetFFz) ;
	
  fhJetFFxi  = new TH2F("JetFFxi","#xi = ln(p_{T jet}/p_{T i charged}) vs p_{T jet}", nptbins,ptmin,ptmax,100,0.,10.); 
  fhJetFFxi->SetYTitle("#xi");
  fhJetFFxi->SetXTitle("p_{T jet}");
  outputContainer->Add(fhJetFFxi) ;
  
  fhJetFFpt  = new TH2F("JetFFpt","p_{T i charged} vs p_{T jet}", nptbins,ptmin,ptmax,100,0.,50.); 
  fhJetFFpt->SetYTitle("p_{T charged hadron}");
  fhJetFFpt->SetXTitle("p_{T jet}");
  outputContainer->Add(fhJetFFpt) ;

  fhJetFFzCor  = new TH2F("JetFFzCor","z = -cos(#alpha(jet,trig))*p_{T i charged}/p_{T jet} vs p_{T jet}",nptbins,ptmin,ptmax,200,0.,2);  
  fhJetFFzCor->SetYTitle("z");
  fhJetFFzCor->SetXTitle("p_{T jet}");
  outputContainer->Add(fhJetFFzCor) ;
	
  fhJetFFxiCor  = new TH2F("JetFFxiCor","#xi = ln(p_{T jet}/(-cos(#alpha(jet,trig))*p_{T i charged})) vs p_{T jet}", nptbins,ptmin,ptmax,100,0.,10.); 
  fhJetFFxiCor->SetYTitle("#xi");
  fhJetFFxiCor->SetXTitle("p_{T jet}");
  outputContainer->Add(fhJetFFxiCor) ;

  fhGamPtPerTrig  = new TH1F("GamPtPerTrig","GamPtPerTrig", nptbins,ptmin,ptmax); 
  fhGamPtPerTrig->SetYTitle("Counts");
  fhGamPtPerTrig->SetXTitle("p_{T, #gamma}");
  outputContainer->Add(fhGamPtPerTrig) ;
  
  fhPtGamPtJet  = new TH2F("PtGamPtJet","p_{T #gamma} vs p_{T jet}", nptbins,ptmin,ptmax,150,-50.,100.); 
  fhPtGamPtJet->SetXTitle("p_{T #gamma}");
  fhPtGamPtJet->SetYTitle("p_{T jet}");
  outputContainer->Add(fhPtGamPtJet) ;


  //background FF
  fhBkgFFz[0]  = new TH2F("BkgFFzRC",  "z = p_{T i charged}/p_{T trigger} vs p_{T trigger} Bkg RC"  ,nptbins,ptmin,ptmax,200,0.,2);  
  fhBkgFFz[1]  = new TH2F("BkgFFzPCG", "z = p_{T i charged}/p_{T trigger} vs p_{T trigger} Bkg PCG" ,nptbins,ptmin,ptmax,200,0.,2);  
  fhBkgFFz[2]  = new TH2F("BkgFFzPCJ", "z = p_{T i charged}/p_{T trigger} vs p_{T trigger} Bkg PCJ" ,nptbins,ptmin,ptmax,200,0.,2);  
  fhBkgFFz[3]  = new TH2F("BkgFFzMP",  "z = p_{T i charged}/p_{T trigger} vs p_{T trigger} Bkg MP"  ,nptbins,ptmin,ptmax,200,0.,2);  
  fhBkgFFz[4]  = new TH2F("BkgFFzTest","z = p_{T i charged}/p_{T trigger} vs p_{T trigger} Bkg Test",nptbins,ptmin,ptmax,200,0.,2);  
  for(Int_t i=0;i<5;i++){
    fhBkgFFz[i]->SetYTitle("z");
    fhBkgFFz[i]->SetXTitle("p_{T trigger}");
    outputContainer->Add(fhBkgFFz[i]) ;
  }

  fhBkgFFxi[0]  = new TH2F("BkgFFxiRC",  "#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger} Bkg RC",  nptbins,ptmin,ptmax,100,0.,10.); 
  fhBkgFFxi[1]  = new TH2F("BkgFFxiPCG", "#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger} Bkg PCG", nptbins,ptmin,ptmax,100,0.,10.); 
  fhBkgFFxi[2]  = new TH2F("BkgFFxiPCJ", "#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger} Bkg PCJ", nptbins,ptmin,ptmax,100,0.,10.); 
  fhBkgFFxi[3]  = new TH2F("BkgFFxiMP",  "#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger} Bkg MP",  nptbins,ptmin,ptmax,100,0.,10.); 
  fhBkgFFxi[4]  = new TH2F("BkgFFxiTest","#xi = ln(p_{T trigger}/p_{T i charged}) vs p_{T trigger} Bkg Test",nptbins,ptmin,ptmax,100,0.,10.); 
  for(Int_t i=0;i<5;i++){
    fhBkgFFxi[i]->SetYTitle("#xi");
    fhBkgFFxi[i]->SetXTitle("p_{T trigger}");
    outputContainer->Add(fhBkgFFxi[i]) ;
  }

  fhBkgFFpt[0]  = new TH2F("BkgFFptRC",  "p_{T i charged} vs p_{T trigger} Bkg RC",   nptbins,ptmin,ptmax,100,0.,50.); 
  fhBkgFFpt[1]  = new TH2F("BkgFFptPCG", "p_{T i charged} vs p_{T trigger} Bkg PCG",  nptbins,ptmin,ptmax,100,0.,50.); 
  fhBkgFFpt[2]  = new TH2F("BkgFFptPCJ", "p_{T i charged} vs p_{T trigger} Bkg PCJ",  nptbins,ptmin,ptmax,100,0.,50.); 
  fhBkgFFpt[3]  = new TH2F("BkgFFptMP",  "p_{T i charged} vs p_{T trigger} Bkg MP",   nptbins,ptmin,ptmax,100,0.,50.); 
  fhBkgFFpt[4]  = new TH2F("BkgFFptTest","p_{T i charged} vs p_{T trigger} Bkg Test", nptbins,ptmin,ptmax,100,0.,50.); 
  for(Int_t i=0;i<5;i++){
    fhBkgFFpt[i]->SetYTitle("p_{T charged hadron}");
    fhBkgFFpt[i]->SetXTitle("p_{T trigger}");
    outputContainer->Add(fhBkgFFpt[i]) ;
  }

  fhBkgNTracksInCone[0]  = new TH2F("BkgNTracksInConeRC",  "Number of tracks in cone vs p_{T trigger} Bkg RC",   nptbins,ptmin,ptmax,100,0.,150.); 
  fhBkgNTracksInCone[1]  = new TH2F("BkgNTracksInConePCG", "Number of tracks in cone vs p_{T trigger} Bkg PCG",  nptbins,ptmin,ptmax,100,0.,150.); 
  fhBkgNTracksInCone[2]  = new TH2F("BkgNTracksInConePCJ", "Number of tracks in cone vs p_{T trigger} Bkg PCJ",  nptbins,ptmin,ptmax,100,0.,150.); 
  fhBkgNTracksInCone[3]  = new TH2F("BkgNTracksInConeMP",  "Number of tracks in cone vs p_{T trigger} Bkg MP",   nptbins,ptmin,ptmax,100,0.,150.); 
  fhBkgNTracksInCone[4]  = new TH2F("BkgNTracksInConeTest","Number of tracks in cone vs p_{T trigger} Bkg Test", nptbins,ptmin,ptmax,100,0.,150.); 
  for(Int_t i=0;i<5;i++){
    fhBkgNTracksInCone[i]->SetYTitle("Number of tracks");
    fhBkgNTracksInCone[i]->SetXTitle("p_{T trigger}");
    outputContainer->Add(fhBkgNTracksInCone[i]) ;
  }

  fhBkgSumPtInCone[0]  = new TH2F("BkgSumPtInConeRC",  "Sum P_{T} in cone vs p_{T trigger} Bkg RC",   nptbins,ptmin,ptmax,100,0.,100.);
  fhBkgSumPtInCone[1]  = new TH2F("BkgSumPtInConePCG", "Sum P_{T} in cone vs p_{T trigger} Bkg PCG",  nptbins,ptmin,ptmax,100,0.,100.);
  fhBkgSumPtInCone[2]  = new TH2F("BkgSumPtInConePCJ", "Sum P_{T} in cone vs p_{T trigger} Bkg PCJ",  nptbins,ptmin,ptmax,100,0.,100.);
  fhBkgSumPtInCone[3]  = new TH2F("BkgSumPtInConeMP",  "Sum P_{T} in cone vs p_{T trigger} Bkg MP",   nptbins,ptmin,ptmax,100,0.,100.);
  fhBkgSumPtInCone[4]  = new TH2F("BkgSumPtInConeTest","Sum P_{T} in cone vs p_{T trigger} Bkg Test", nptbins,ptmin,ptmax,100,0.,100.);
  for(Int_t i=0;i<5;i++){
    fhBkgSumPtInCone[i]->SetYTitle("Sum P_{T}");
    fhBkgSumPtInCone[i]->SetXTitle("p_{T trigger}");
    outputContainer->Add(fhBkgSumPtInCone[i]) ;
  }

  fhBkgSumPtnTracksInCone[0]  = new TH2F("BkgSumPtnTracksInConeRC",  "Sum p_{T} / Number of tracks in cone vs p_{T trigger} Bkg RC",   nptbins,ptmin,ptmax,100,0.,20.);
  fhBkgSumPtnTracksInCone[1]  = new TH2F("BkgSumPtnTracksInConePCG", "Sum p_{T} / Number of tracks in cone vs p_{T trigger} Bkg PCG",  nptbins,ptmin,ptmax,100,0.,20.);
  fhBkgSumPtnTracksInCone[2]  = new TH2F("BkgSumPtnTracksInConePCJ", "Sum p_{T} / Number of tracks in cone vs p_{T trigger} Bkg PCJ",  nptbins,ptmin,ptmax,100,0.,20.);
  fhBkgSumPtnTracksInCone[3]  = new TH2F("BkgSumPtnTracksInConeMP",  "Sum p_{T} / Number of tracks in cone vs p_{T trigger} Bkg MP",   nptbins,ptmin,ptmax,100,0.,20.);
  fhBkgSumPtnTracksInCone[4]  = new TH2F("BkgSumPtnTracksInConeTest","Sum p_{T} / Number of tracks in cone vs p_{T trigger} Bkg Test", nptbins,ptmin,ptmax,100,0.,20.);
  for(Int_t i=0;i<5;i++){
    fhBkgSumPtnTracksInCone[i]->SetYTitle("Sum p_{T}/Number of tracks");
    fhBkgSumPtnTracksInCone[i]->SetXTitle("p_{T trigger}");
    outputContainer->Add(fhBkgSumPtnTracksInCone[i]) ;
  }


  //temporary histograms
  fhNjetsNgammas  = new TH2F("NjetsNgammas"," Number of jets vs number of gammas in event",20,0.,100.,10,0.,80.);  
  fhNjetsNgammas->SetYTitle("Number of gammas");
  fhNjetsNgammas->SetXTitle("Number of jets");
  outputContainer->Add(fhNjetsNgammas) ;

  fhCuts  = new TH1F("Cuts"," Cuts",10,0.,10.);  
  fhCuts->SetYTitle("Counts");
  fhCuts->SetXTitle("Cut number");
  outputContainer->Add(fhCuts) ;

  fhDeltaPhiBefore  = new TH2F("DeltaPhiBefore","#phi_{trigger} - #phi_{jet} vs p_{T trigger}",nptbins,ptmin,ptmax,100,0,6.5); 
  fhDeltaPhiBefore->SetYTitle("#Delta #phi");
  fhDeltaPhiBefore->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhDeltaPhiBefore);
  
  fhDeltaEtaBefore  = new TH2F("DeltaEtaBefore","#eta_{trigger} - #eta_{jet} vs p_{T trigger}",nptbins,ptmin,ptmax,100,-2,2); 
  fhDeltaEtaBefore->SetYTitle("#Delta #eta");
  fhDeltaEtaBefore->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhDeltaEtaBefore);
  
  fhDeltaPtBefore  = new TH2F("DeltaPtBefore","p_{T trigger} - p_{T jet} vs p_{T trigger}",nptbins,ptmin,ptmax,100,-50,50); 
  fhDeltaPtBefore->SetYTitle("#Delta p_{T}");
  fhDeltaPtBefore->SetXTitle("p_{T trigger} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtBefore);
  
  fhPtRatioBefore  = new TH2F("PtRatioBefore","p_{T jet} / p_{T trigger} vs p_{T trigger}",nptbins,ptmin,ptmax,200,0,2.); 
  fhPtRatioBefore->SetYTitle("ratio");
  fhPtRatioBefore->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhPtRatioBefore);
  
  fhPtBefore  = new TH2F("PtBefore","p_{T jet} vs p_{T trigger}",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
  fhPtBefore->SetYTitle("p_{T jet}(GeV/c)");
  fhPtBefore->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhPtBefore);

  fhDeltaPhi0PiCorrectBefore  = new TH2F("DeltaPhi0PiCorrectBefore","#phi_{trigger} - #phi_{jet} (0,#pi) vs p_{T trigger}",nptbins,ptmin,ptmax,100,0,3.5); 
  fhDeltaPhi0PiCorrectBefore->SetYTitle("#Delta #phi");
  fhDeltaPhi0PiCorrectBefore->SetXTitle("p_{T trigger} (GeV/c)");
  outputContainer->Add(fhDeltaPhi0PiCorrectBefore);

  //temporary jet histograms
  fhJetPtBefore            = new TH1F("JetPtBefore","JetPtBefore",150,-50,100); 
  fhJetPtBefore->SetYTitle("Counts");
  fhJetPtBefore->SetXTitle("p_{T jet}(GeV/c)");
  outputContainer->Add(fhJetPtBefore) ;

  fhJetPtBeforeCut            = new TH1F("JetPtBeforeCut","JetPtBeforeCut",150,-50,100); 
  fhJetPtBeforeCut->SetYTitle("Counts");
  fhJetPtBeforeCut->SetXTitle("p_{T jet}(GeV/c)");
  outputContainer->Add(fhJetPtBeforeCut) ;

  fhJetPt            = new TH1F("JetPt","JetPt",150,-50,100); 
  fhJetPt->SetYTitle("Counts");
  fhJetPt->SetXTitle("p_{T jet}(GeV/c)");
  outputContainer->Add(fhJetPt) ;

  fhJetPtMostEne            = new TH1F("JetPtMostEne","JetPtMostEne",150,0,150); 
  fhJetPtMostEne->SetYTitle("Counts");
  fhJetPtMostEne->SetXTitle("p_{T jet}(GeV/c)");
  outputContainer->Add(fhJetPtMostEne) ;

  fhJetPhi	     = new TH1F("JetPhi","JetPhi",130,0,6.5); 
  fhJetPhi->SetYTitle("Counts");
  fhJetPhi->SetXTitle("#phi_{jet}");
  outputContainer->Add(fhJetPhi) ;

  fhJetEta	     = new TH1F("JetEta","JetEta",100,-1,1); 
  fhJetEta->SetYTitle("Counts");
  fhJetEta->SetXTitle("#eta_{jet}");
  outputContainer->Add(fhJetEta) ;

  fhJetEtaVsPt      = new TH2F("JetEtaVsPt","JetEtaVsPt",100,0,100,50,-1,1);
  fhJetEtaVsPt->SetYTitle("#eta_{jet}");
  fhJetEtaVsPt->SetXTitle("p_{T,jet}(GeV/c)");
  outputContainer->Add(fhJetEtaVsPt) ;

  fhJetPhiVsEta      = new TH2F("JetPhiVsEta","JetPhiVsEta",65,0,6.5,50,-1,1); 
  fhJetPhiVsEta->SetYTitle("#eta_{jet}");
  fhJetPhiVsEta->SetXTitle("#phi_{jet}");
  outputContainer->Add(fhJetPhiVsEta) ;

  fhJetEtaVsNpartInJet= new TH2F("JetEtaVsNpartInJet","JetEtaVsNpartInJet",50,-1,1,100,0.,200.); 
  fhJetEtaVsNpartInJet->SetYTitle("N_{tracks-in-jet}");
  fhJetEtaVsNpartInJet->SetXTitle("#eta_{jet}");
  outputContainer->Add(fhJetEtaVsNpartInJet) ;

  fhJetEtaVsNpartInJetBkg= new TH2F("JetEtaVsNpartInJetBkg","JetEtaVsNpartInJetBkg",50,-1,1,100,0.,200.); 
  fhJetEtaVsNpartInJetBkg->SetYTitle("N_{tracks-in-jet}");
  fhJetEtaVsNpartInJetBkg->SetXTitle("#eta_{jet}");
  outputContainer->Add(fhJetEtaVsNpartInJetBkg) ;

  fhJetChBkgEnergyVsPt = new TH2F("JetBkgChEnergyVsPt","JetBkgChEnergyVsPt",100,0,100,90,0,90); 
  fhJetChBkgEnergyVsPt->SetYTitle("Jet Bkg Energy (GeV)");
  fhJetChBkgEnergyVsPt->SetXTitle("p_{T,jet} (GeV/c)");
  outputContainer->Add(fhJetChBkgEnergyVsPt);
  
  fhJetChAreaVsPt      = new TH2F("JetChAreaVsPt","JetChAreaVsPt",100,0,100,50,0,1); 
  fhJetChAreaVsPt->SetYTitle("Jet Area");
  fhJetChAreaVsPt->SetXTitle("p_{T,jet} (GeV/c)");
  outputContainer->Add(fhJetChAreaVsPt);
  
  if(IsHistogramTracks()){
    fhTrackPhiVsEta      = new TH2F("TrackPhiVsEta","TrackPhiVsEta",65,0,6.5,50,-1,1); 
    fhTrackPhiVsEta->SetYTitle("#eta_{track}");
    fhTrackPhiVsEta->SetXTitle("#phi_{track}");
    outputContainer->Add(fhTrackPhiVsEta) ;

    fhTrackAveTrackPt      = new TH1F("TrackAveTrackPt","TrackAveTrackPt",45,0,1.5);
    fhTrackAveTrackPt->SetYTitle("Counts");
    fhTrackAveTrackPt->SetXTitle("Average p_{T,track} (GeV/c)");
    outputContainer->Add(fhTrackAveTrackPt);
  
  }//end of IsHistogramTracks()

  for(Int_t i=0;i<10;i++){
    fhJetNjetOverPtCut[i]      = new TH1F(Form("JetNjetOverPtCut%d", i),Form("JetNjetOverPtCut%d", i),100,0,100);
    fhJetNjetOverPtCut[i]->SetYTitle("Counts");
    fhJetNjetOverPtCut[i]->SetXTitle("N_{jets} over threshold");
    outputContainer->Add(fhJetNjetOverPtCut[i]);
  }

  fhJetChBkgEnergyVsArea = new TH2F("JetBkgChEnergyVsArea","JetBkgChEnergyVsArea",100,0,100,70,0,0.7); 
  fhJetChBkgEnergyVsArea->SetXTitle("Jet Bkg Energy (GeV)");
  fhJetChBkgEnergyVsArea->SetYTitle("Area");
  outputContainer->Add(fhJetChBkgEnergyVsArea);

  fhJetRhoVsPt           = new TH2F("JetRhoVsPt","JetRhoVsPt",100,0,100,100,0,150); 
  fhJetRhoVsPt->SetYTitle("Rho");
  fhJetRhoVsPt->SetXTitle("p_{T,jet} (GeV/c)");
  outputContainer->Add(fhJetRhoVsPt);

  if(IsHistogramJetBkg()){
    fhJetRhoVsCentrality           = new TH2F("JetRhoVsCentrality","JetRhoVsCentrality",100,0,100,100,0,200);
    fhJetRhoVsCentrality->SetYTitle("Rho");
    fhJetRhoVsCentrality->SetXTitle("Centrality");
    outputContainer->Add(fhJetRhoVsCentrality);
  }

  fhJetNparticlesInJet           = new TH1F("JetNparticlesInJet","JetNparticlesInJet",100,0,200);
  fhJetNparticlesInJet->SetXTitle("N^{particles}");
  fhJetNparticlesInJet->SetYTitle("N^{jets}");
  outputContainer->Add(fhJetNparticlesInJet);

  fhJetDeltaEtaDeltaPhi      = new TH2F("JetDeltaEtaDeltaPhi","#Delta #eta^{jet-track} vs. #Delta #phi^{jet-track} for jet tracks",100,-0.8,0.8,100,-0.8,0.8);
  fhJetDeltaEtaDeltaPhi->SetXTitle("#Delta #eta^{jet-track}");
  fhJetDeltaEtaDeltaPhi->SetYTitle("#Delta #phi^{jet-track}");
  outputContainer->Add(fhJetDeltaEtaDeltaPhi );


  fhJetDeltaEtaDeltaPhiAllTracks      = new TH2F("JetDeltaEtaDeltaPhiAllTracks","#Delta #eta^{jet-track} vs. #Delta #phi^{jet-track} for all tracks",100,-3.2,3.2,100,-3.2,3.2);
  fhJetDeltaEtaDeltaPhiAllTracks->SetXTitle("#Delta #eta^{jet-track}");
  fhJetDeltaEtaDeltaPhiAllTracks->SetYTitle("#Delta #phi^{jet-track}");
  outputContainer->Add(fhJetDeltaEtaDeltaPhiAllTracks);


  if(IsHistogramJetTracks()){
    fhJetAveTrackPt           = new TH1F("JetAveTrackPt","JetAveTrackPt",45,0.,1.5);
    fhJetAveTrackPt->SetXTitle("Average p_{T,track} (GeV/c)");
    fhJetAveTrackPt->SetYTitle("Counts");
    outputContainer->Add(fhJetAveTrackPt);
    
    for(Int_t i=0;i<6;i++){
      if(i==0) fhJetNtracksInJetAboveThr[i]      = new TH2F(Form("JetNtracksInJetAboveThr%d", i),Form("JetNtracksInJetAboveThr%d", i),100,0,100,100,0,200);
      else fhJetNtracksInJetAboveThr[i]      = new TH2F(Form("JetNtracksInJetAboveThr%d", i),Form("JetNtracksInJetAboveThr%d", i),100,0,100,100,0,100);
      fhJetNtracksInJetAboveThr[i]->SetXTitle("p_{T,jet} (GeV/c)");
      fhJetNtracksInJetAboveThr[i]->SetYTitle("N_{tracks} over threshold");
      outputContainer->Add(fhJetNtracksInJetAboveThr[i]);
    }
    
    for(Int_t i=0;i<5;i++){
      fhJetRatioNTrkAboveToNTrk[i]      = new TH2F(Form("JetRatioNTrkAboveToNTrk%d", i),Form("JetRatioNTrkAboveToNTrk%d", i),100,0,100,40,0,1);
      fhJetRatioNTrkAboveToNTrk[i]->SetXTitle("p_{T,jet} (GeV/c)");
      fhJetRatioNTrkAboveToNTrk[i]->SetYTitle("Ratio N_{tracks} over threshold to N_{tracks}");
      outputContainer->Add(fhJetRatioNTrkAboveToNTrk[i]);
      
      fhJetNtrackRatioMostEne[i]      = new TH2F(Form("JetNtrackRatioMostEne%d", i),Form("JetNtrackRatioMostEne%d", i),100,0,100,40,0,1);
      fhJetNtrackRatioMostEne[i]->SetXTitle("p_{T,jet} (GeV/c)");
      fhJetNtrackRatioMostEne[i]->SetYTitle("Ratio N_{tracks} over threshold to N_{tracks}");
      outputContainer->Add(fhJetNtrackRatioMostEne[i]);
      
      fhJetNtrackRatioJet5GeV[i]      = new TH2F(Form("JetNtrackRatioJet5GeV%d", i),Form("JetNtrackRatioJet5GeV%d", i),100,0,100,40,0,1);
      fhJetNtrackRatioJet5GeV[i]->SetXTitle("p_{T,jet} (GeV/c)");
      fhJetNtrackRatioJet5GeV[i]->SetYTitle("Ratio N_{tracks} over threshold to N_{tracks}");
      outputContainer->Add(fhJetNtrackRatioJet5GeV[i]);
      
      fhJetNtrackRatioLead5GeV[i]      = new TH2F(Form("JetNtrackRatioLead5GeV%d", i),Form("JetNtrackRatioLead5GeV%d", i),100,0,100,40,0,1);
      fhJetNtrackRatioLead5GeV[i]->SetXTitle("p_{T,jet} (GeV/c)");
      fhJetNtrackRatioLead5GeV[i]->SetYTitle("Ratio N_{tracks} over threshold to N_{tracks}");
      outputContainer->Add(fhJetNtrackRatioLead5GeV[i]);
    }
  }//end of if IsHistogramJetTracks

  //temporary background jets histograms
  if(IsHistogramJetBkg()){
    for(Int_t i=0;i<4;i++){
      fhBkgJetBackground[i]      = new TH1F(Form("BkgJetBackground%d", i),Form("BkgJetBackground%d", i),100,0,200);
      fhBkgJetBackground[i]->SetXTitle("<#rho> (GeV/c)");
      fhBkgJetBackground[i]->SetYTitle("Counts");
      outputContainer->Add(fhBkgJetBackground[i]);
      
      fhBkgJetSigma[i]      = new TH1F(Form("BkgJetSigma%d", i),Form("BkgJetSigma%d", i),100,0,50);
      fhBkgJetSigma[i]->SetXTitle("#sigma (GeV/c)");
      fhBkgJetSigma[i]->SetYTitle("Counts");
      outputContainer->Add(fhBkgJetSigma[i]);
      
      fhBkgJetArea[i]      = new TH1F(Form("BkgJetArea%d", i),Form("BkgJetArea%d", i),100,0,1);
      fhBkgJetArea[i]->SetXTitle("<A>");
      fhBkgJetArea[i]->SetYTitle("Counts");
      outputContainer->Add(fhBkgJetArea[i]);
    }
  }

  //temporary photon histograms
  fhPhotonPtMostEne = new TH1F("PhotonPtMostEne","PhotonPtMostEne",100,0,100);
  fhPhotonPtMostEne->SetYTitle("Counts");
  fhPhotonPtMostEne->SetXTitle("p_{T,#gamma} (GeV/c)");
  outputContainer->Add(fhPhotonPtMostEne);

//  fhPhotonIndexMostEne = new TH1F("PhotonIndexMostEne","PhotonIndexMostEne",100,0,100);
//  fhPhotonIndexMostEne->SetYTitle("Counts");
//  fhPhotonIndexMostEne->SetXTitle("Index");
//  outputContainer->Add(fhPhotonIndexMostEne);

  fhPhotonAverageEnergy = new TH1F("PhotonAverageEnergy","PhotonAverageEnergy",100,0,10);
  fhPhotonAverageEnergy->SetYTitle("Counts");
  fhPhotonAverageEnergy->SetXTitle("p_{T,#gamma} (GeV/c)");
  outputContainer->Add(fhPhotonAverageEnergy);

  fhPhotonRatioAveEneToMostEne = new TH1F("PhotonRatioAveEneToMostEne","PhotonRatioAveEneToMostEne",100,0,1);
  fhPhotonRatioAveEneToMostEne->SetYTitle("Counts");
  fhPhotonRatioAveEneToMostEne->SetXTitle("Ratio");
  outputContainer->Add(fhPhotonRatioAveEneToMostEne);

  fhPhotonAverageEnergyMinus1 = new TH1F("PhotonAverageEnergyMinus1","PhotonAverageEnergyMinus1",100,0,10);
  fhPhotonAverageEnergyMinus1->SetYTitle("Counts");
  fhPhotonAverageEnergyMinus1->SetXTitle("p_{T,#gamma} (GeV/c)");
  outputContainer->Add(fhPhotonAverageEnergyMinus1);

  fhPhotonRatioAveEneMinus1ToMostEne = new TH1F("PhotonRatioAveEneMinus1ToMostEne","PhotonRatioAveEneMinus1ToMostEne",100,0,1);
  fhPhotonRatioAveEneMinus1ToMostEne->SetYTitle("Counts");
  fhPhotonRatioAveEneMinus1ToMostEne->SetXTitle("Ratio");
  outputContainer->Add(fhPhotonRatioAveEneMinus1ToMostEne);

  fhPhotonNgammaMoreAverageToNgamma = new TH1F("PhotonNgammaMoreAverageToNgamma","PhotonNgammaMoreAverageToNgamma",100,0,1);
  fhPhotonNgammaMoreAverageToNgamma->SetYTitle("Counts");
  fhPhotonNgammaMoreAverageToNgamma->SetXTitle("Ratio");
  outputContainer->Add(fhPhotonNgammaMoreAverageToNgamma);

  fhPhotonNgammaMoreAverageMinus1ToNgamma = new TH1F("PhotonNgammaMoreAverageMinus1ToNgamma","PhotonNgammaMoreAverageMinus1ToNgamma",100,0,1);
  fhPhotonNgammaMoreAverageMinus1ToNgamma->SetYTitle("Counts");
  fhPhotonNgammaMoreAverageMinus1ToNgamma->SetXTitle("Ratio");
  outputContainer->Add(fhPhotonNgammaMoreAverageMinus1ToNgamma);

  for(Int_t i=0;i<10;i++){
    fhPhotonNgammaOverPtCut[i]      = new TH1F(Form("PhotonNgammaOverPtCut%d",i),Form("PhotonNgammaOverPtCut%d",i),100,0,100);
    fhPhotonNgammaOverPtCut[i]->SetYTitle("Counts");
    fhPhotonNgammaOverPtCut[i]->SetXTitle("N_{#gamma} over threshold");
    outputContainer->Add(fhPhotonNgammaOverPtCut[i]);
  }

  fhPhotonBkgRhoVsNtracks = new TH2F("PhotonBkgRhoVsNtracks","PhotonBkgRhoVsNtracks",200,0,2500,75,0,1.5);
  //fhPhotonBkgRhoVsNtracks->SetXTitle("Counts");
  fhPhotonBkgRhoVsNtracks->SetXTitle("Ntracks");
  fhPhotonBkgRhoVsNtracks->SetYTitle("Rho");
  outputContainer->Add(fhPhotonBkgRhoVsNtracks);

  fhPhotonBkgRhoVsNclusters = new TH2F("PhotonBkgRhoVsNclusters","PhotonBkgRhoVsNclusters",50,0,100,75,0,1.5);
  fhPhotonBkgRhoVsNclusters->SetXTitle("Nclusters");
  fhPhotonBkgRhoVsNclusters->SetYTitle("Rho");
  outputContainer->Add(fhPhotonBkgRhoVsNclusters);

  fhPhotonBkgRhoVsCentrality = new TH2F("PhotonBkgRhoVsCentrality","PhotonBkgRhoVsCentrality",100,0,100,75,0,1.5);
  fhPhotonBkgRhoVsCentrality->SetXTitle("Centrality");
  fhPhotonBkgRhoVsCentrality->SetYTitle("Rho");
  outputContainer->Add(fhPhotonBkgRhoVsCentrality);

  fhPhotonBkgRhoVsNcells = new TH2F("PhotonBkgRhoVsNcells","PhotonBkgRhoVsNcells",100,0,200,75,0,1.5);
  fhPhotonBkgRhoVsNcells->SetXTitle("N_{cells}");
  fhPhotonBkgRhoVsNcells->SetYTitle("Rho");
  outputContainer->Add(fhPhotonBkgRhoVsNcells);

  fhPhotonPt = new TH1F("PhotonPt","PhotonPt",220,-10,100);
  fhPhotonPt->SetXTitle("p_{T,#gamma} (GeV/c)");
  fhPhotonPt->SetYTitle("Counts");
  outputContainer->Add(fhPhotonPt);

  fhPhotonPtCorrected = new TH1F("PhotonPtCorrected","PhotonPtCorrected",220,-10,100);
  fhPhotonPtCorrected->SetXTitle("p_{T,#gamma} (GeV/c)");
  fhPhotonPtCorrected->SetYTitle("Counts");
  outputContainer->Add(fhPhotonPtCorrected);

  fhPhotonPtDiff = new TH1F("PhotonPtDiff","PhotonPtDiff",50,0,10);
  fhPhotonPtDiff->SetXTitle("p_{T,#gamma} (GeV/c)");
  fhPhotonPtDiff->SetYTitle("Counts");
  outputContainer->Add(fhPhotonPtDiff);

  fhPhotonPtDiffVsNtracks = new TH2F("PhotonPtDiffVsNtracks","PhotonPtDiffVsNtracks",200,0,2500,50,0,5);
  fhPhotonPtDiffVsNtracks->SetXTitle("N_{tracks}");
  fhPhotonPtDiffVsNtracks->SetYTitle("<#rho^{#gamma}>*N_{cells}");
  outputContainer->Add(fhPhotonPtDiffVsNtracks);

  fhPhotonPtDiffVsNclusters = new TH2F("PhotonPtDiffVsNclusters","PhotonPtDiffVsNclusters",50,0,100,50,0,5);
  fhPhotonPtDiffVsNclusters->SetXTitle("N_{clusters}");
  fhPhotonPtDiffVsNclusters->SetYTitle("<#rho^{#gamma}>*N_{cells}");
  outputContainer->Add(fhPhotonPtDiffVsNclusters);

  fhPhotonPtDiffVsCentrality = new TH2F("PhotonPtDiffVsCentrality","PhotonPtDiffVsCentrality",100,0,100,50,0,5);
  fhPhotonPtDiffVsCentrality->SetXTitle("Centrality");
  fhPhotonPtDiffVsCentrality->SetYTitle("<#rho^{#gamma}>*N_{cells}");
  outputContainer->Add(fhPhotonPtDiffVsCentrality);

  fhPhotonPtDiffVsNcells = new TH2F("PhotonPtDiffVsNcells","PhotonPtDiffVsNcells",100,0,200,50,0,5);
  fhPhotonPtDiffVsNcells->SetXTitle("N_{cells}");
  fhPhotonPtDiffVsNcells->SetYTitle("<#rho^{#gamma}>*N_{cells}");
  outputContainer->Add(fhPhotonPtDiffVsNcells);


  fhPhotonPtCorrectedZoom = new TH1F("PhotonPtCorrectedZoom","PhotonPtCorrectedZoom",200,-5,5);
  fhPhotonPtCorrectedZoom->SetXTitle("p_{T,#gamma} (GeV/c)");
  fhPhotonPtCorrectedZoom->SetYTitle("Counts");
  outputContainer->Add(fhPhotonPtCorrectedZoom);

  fhPhotonSumPtInCone = new TH1F("PhotonSumPtInCone","PhotonSumPtInCone",100,0,100);
  fhPhotonSumPtInCone->SetXTitle("#Sigma p_{T,#gamma} (GeV/c)");
  fhPhotonSumPtInCone->SetYTitle("Counts");
  outputContainer->Add(fhPhotonSumPtInCone);

  fhPhotonSumPtCorrectInCone = new TH1F("PhotonSumPtCorrectInCone","PhotonSumPtCorrectInCone",100,-20,80);
  fhPhotonSumPtCorrectInCone->SetXTitle("#Sigma p_{T,#gamma} (GeV/c)");
  fhPhotonSumPtCorrectInCone->SetYTitle("Counts");
  outputContainer->Add(fhPhotonSumPtCorrectInCone);

  fhPhotonSumPtChargedInCone = new TH1F("PhotonSumPtChargedInCone","PhotonSumPtChargedInCone",100,0,100);
  fhPhotonSumPtChargedInCone->SetXTitle("#Sigma p_{T,#gamma}^{ch} (GeV/c)");
  fhPhotonSumPtChargedInCone->SetYTitle("Counts");
  outputContainer->Add(fhPhotonSumPtChargedInCone);


  //temporary jet histograms after selection
  fhSelectedJetPhiVsEta      = new TH2F("SelectedJetSelectedPhiVsEta","SelectedJetPhiVsEta",65,0,6.5,50,-1,1); 
  fhSelectedJetPhiVsEta->SetYTitle("#eta_{jet}");
  fhSelectedJetPhiVsEta->SetXTitle("#phi_{jet}");
  outputContainer->Add(fhSelectedJetPhiVsEta) ;

  fhSelectedJetChBkgEnergyVsPtJet = new TH2F("SelectedJetBkgChEnergyVsPtJet","SelectedJetBkgChEnergyVsPtJet",100,0,100,90,0,90); 
  fhSelectedJetChBkgEnergyVsPtJet->SetYTitle("Jet Bkg Energy (GeV)");
  fhSelectedJetChBkgEnergyVsPtJet->SetXTitle("p_{T,jet} (GeV/c)");
  outputContainer->Add(fhSelectedJetChBkgEnergyVsPtJet);
  
  fhSelectedJetChAreaVsPtJet      = new TH2F("SelectedJetChAreaVsPtJet","SelectedJetChAreaVsPtJet",100,0,100,50,0,1); 
  fhSelectedJetChAreaVsPtJet->SetYTitle("Jet Area");
  fhSelectedJetChAreaVsPtJet->SetXTitle("p_{T,jet} (GeV/c)");
  outputContainer->Add(fhSelectedJetChAreaVsPtJet);

  fhSelectedJetNjet      = new TH1F("SelectedJetNjet","SelectedJetNjet",100,0,100); 
  fhSelectedJetNjet->SetYTitle("Counts");
  fhSelectedJetNjet->SetXTitle("N_{jets} per event");
  outputContainer->Add(fhSelectedJetNjet);

  fhSelectedNtracks      = new TH1F("SelectedNtracks","SelectedNtracks",100,0,2000); 
  fhSelectedNtracks->SetYTitle("Counts");
  fhSelectedNtracks->SetXTitle("N_{tracks} per event");
  outputContainer->Add(fhSelectedNtracks);

  fhSelectedTrackPhiVsEta      = new TH2F("SelectedTrackPhiVsEta","SelectedTrackPhiVsEta",65,0,6.5,50,-1,1); 
  fhSelectedTrackPhiVsEta->SetYTitle("#eta_{track}");
  fhSelectedTrackPhiVsEta->SetXTitle("#phi_{track}");
  outputContainer->Add(fhSelectedTrackPhiVsEta) ;

  fhCuts2      = new TH1F("Cuts2","Cuts2",10,0.,10.); 
  fhCuts2->SetYTitle("Counts");
  fhCuts2->SetXTitle("Cut number");
  outputContainer->Add(fhCuts2);

  fhSelectedPhotonNLMVsPt      = new TH2F("SelectedPhotonNLMVsPt","SelectedPhotonNLMVsPt",100,0,100,10,0,10);
  fhSelectedPhotonNLMVsPt->SetYTitle("NLM");
  fhSelectedPhotonNLMVsPt->SetXTitle("p_{T,#gamma} (GeV/c)");
  outputContainer->Add(fhSelectedPhotonNLMVsPt);

  fhSelectedPhotonLambda0VsPt      = new TH2F("SelectedPhotonLambda0VsPt","SelectedPhotonLambda0VsPt",100,0,100,50,0,5);
  fhSelectedPhotonLambda0VsPt->SetYTitle("#lambda_{0}");
  fhSelectedPhotonLambda0VsPt->SetXTitle("p_{T,#gamma} (GeV/c)");
  outputContainer->Add(fhSelectedPhotonLambda0VsPt);

  //random
  fhRandomPhiEta[0]  = new TH2F("RandomPhiEtaRC",  "RandomPhiEtaRC",            100,0,6.5,100,-1.,1.);
  fhRandomPhiEta[1]  = new TH2F("RandomPhiEtaPCG", "RandomPhiEtaPerpConePhoton",100,0,6.5,100,-1.,1.);
  fhRandomPhiEta[2]  = new TH2F("RandomPhiEtaPCJ", "RandomPhiEtaPerpConeJet",   100,0,6.5,100,-1.,1.);
  fhRandomPhiEta[3]  = new TH2F("RandomPhiEtaMP",  "RandomPhiEtaMidPoint",      100,0,6.5,100,-1.,1.);
  fhRandomPhiEta[4]  = new TH2F("RandomPhiEtaTest","RandomPhiEtaTest",          100,0,6.5,100,-1.,1.);
  for(Int_t i=0;i<5;i++){
    fhRandomPhiEta[i]->SetXTitle("#phi");
    fhRandomPhiEta[i]->SetYTitle("#eta");
    outputContainer->Add(fhRandomPhiEta[i]);
  }

  //MC generated photons and jets information  
  if(fMCStudies) {
    fhMCPhotonCuts      = new TH1F("MCPhotonCuts","MCPhotonCuts",10,0.,10.); 
    fhMCPhotonCuts->SetYTitle("Counts");
    fhMCPhotonCuts->SetXTitle("Cut number");
    outputContainer->Add(fhMCPhotonCuts);
    
    fhMCPhotonPt      = new TH1F("MCPhotonPt","MCPhotonPt",100,0.,100.); 
    fhMCPhotonPt->SetYTitle("Counts");
    fhMCPhotonPt->SetXTitle("p_{T,#gamma}^{gen} (GeV/c)");
    outputContainer->Add(fhMCPhotonPt);
    
    fhMCPhotonEtaPhi      = new TH2F("MCPhotonEtaPhi","MCPhotonEtaPhi",65,0,6.5,50,-1,1); 
    fhMCPhotonEtaPhi->SetYTitle("#eta_{#gamma}");
    fhMCPhotonEtaPhi->SetXTitle("#phi_{#gamma}");
    outputContainer->Add(fhMCPhotonEtaPhi) ;
    
    fhMCJetOrigin      = new TH1F("MCJetOrigin","MCJetOrigin",35,-10.,25.); 
    fhMCJetOrigin->SetYTitle("Counts");
    fhMCJetOrigin->SetXTitle("PDG number");
    outputContainer->Add(fhMCJetOrigin);
    
    fhMCJetNPartVsPt      = new TH2F("MCJetNPartVsPt","MCJetNPartVsPt",100,0.,100.,100,0.,100.); 
    fhMCJetNPartVsPt->SetYTitle("N_{particles}");
    fhMCJetNPartVsPt->SetXTitle("p_{T,full-jet}^{gen} (GeV/c)");
    outputContainer->Add(fhMCJetNPartVsPt) ;
    
    fhMCJetChNPartVsPt      = new TH2F("MCJetChNPartVsPt","MCJetChNPartVsPt",100,0.,100.,100,0.,100.); 
    fhMCJetChNPartVsPt->SetYTitle("N_{particles}");
    fhMCJetChNPartVsPt->SetXTitle("p_{T,charged-jet}^{gen} (GeV/c)");
    outputContainer->Add(fhMCJetChNPartVsPt) ;
    
    fhMCJetNPart150VsPt      = new TH2F("MCJetNPart150VsPt","MCJetNPart150VsPt",100,0.,100.,100,0.,100.); 
    fhMCJetNPart150VsPt->SetYTitle("N_{particles (p_{T}>150 MeV/c)}");
    fhMCJetNPart150VsPt->SetXTitle("p_{T,full-jet}^{gen} (GeV/c)");
    outputContainer->Add(fhMCJetNPart150VsPt) ;
    
    fhMCJetChNPart150VsPt      = new TH2F("MCJetChNPart150VsPt","MCJetChNPart150VsPt",100,0.,100.,100,0.,100.); 
    fhMCJetChNPart150VsPt->SetYTitle("N_{particles (p_{T}>150 MeV/c)}");
    fhMCJetChNPart150VsPt->SetXTitle("p_{T,charged-jet}^{gen} (GeV/c)");
    outputContainer->Add(fhMCJetChNPart150VsPt) ;
    
    fhMCJetChNPart150ConeVsPt      = new TH2F("MCJetChNPart150ConeVsPt","MCJetChNPart150ConeVsPt",100,0.,100.,50,0.,50.); 
    fhMCJetChNPart150ConeVsPt->SetYTitle("N_{particles (p_{T}>150 MeV/c)}");
    fhMCJetChNPart150ConeVsPt->SetXTitle("p_{T,charged-jet,R=0.4}^{gen} (GeV/c)");
    outputContainer->Add(fhMCJetChNPart150ConeVsPt) ;

    fhMCJetRatioChFull      = new TH1F("MCJetRatioChFull","MCJetRatioChFull",100,0.,1.); 
    fhMCJetRatioChFull->SetYTitle("Counts");
    fhMCJetRatioChFull->SetXTitle("p_{T,charged-jet}^{gen}/p_{T,full-jet}^{gen}");
    outputContainer->Add(fhMCJetRatioChFull);
    
    fhMCJetRatioCh150Ch      = new TH1F("MCJetRatioCh150Ch","MCJetRatioCh150Ch",100,0.,1.); 
    fhMCJetRatioCh150Ch->SetYTitle("Counts");
    fhMCJetRatioCh150Ch->SetXTitle("p_{T,charged-jet,pT>150MeV/c}^{gen}/p_{T,charged-jet}^{gen}");
    outputContainer->Add(fhMCJetRatioCh150Ch);
    
    fhMCJetEtaPhi      = new TH2F("MCJetEtaPhi","MCJetEtaPhi",65,0,6.5,50,-1,1); 
    fhMCJetEtaPhi->SetYTitle("#eta_{jet}");
    fhMCJetEtaPhi->SetXTitle("#phi_{jet}");
    outputContainer->Add(fhMCJetEtaPhi) ;
    
    fhMCJetChEtaPhi      = new TH2F("MCJetChEtaPhi","MCJetChEtaPhi",65,0,6.5,50,-1,1); 
    fhMCJetChEtaPhi->SetYTitle("#eta_{jet}");
    fhMCJetChEtaPhi->SetXTitle("#phi_{jet}");
    outputContainer->Add(fhMCJetChEtaPhi) ;
    
    fhMCJet150EtaPhi      = new TH2F("MCJet150EtaPhi","MCJet150EtaPhi",65,0,6.5,50,-1,1); 
    fhMCJet150EtaPhi->SetYTitle("#eta_{jet}");
    fhMCJet150EtaPhi->SetXTitle("#phi_{jet}");
    outputContainer->Add(fhMCJet150EtaPhi) ;
    
    fhMCJetCh150EtaPhi      = new TH2F("MCJetCh150EtaPhi","MCJetCh150EtaPhi",65,0,6.5,50,-1,1); 
    fhMCJetCh150EtaPhi->SetYTitle("#eta_{jet}");
    fhMCJetCh150EtaPhi->SetXTitle("#phi_{jet}");
    outputContainer->Add(fhMCJetCh150EtaPhi) ;

    fhMCJetCh150ConeEtaPhi      = new TH2F("MCJetCh150ConeEtaPhi","MCJetCh150ConeEtaPhi",65,0,6.5,50,-1,1); 
    fhMCJetCh150ConeEtaPhi->SetYTitle("#eta_{jet}");
    fhMCJetCh150ConeEtaPhi->SetXTitle("#phi_{jet}");
    outputContainer->Add(fhMCJetCh150ConeEtaPhi) ;
  }

  //tree with data
  if(fSaveGJTree){
    fTreeGJ=new TTree("fTreeGJ","fTreeGJ");
    fTreeGJ->Branch("fGamPt"      ,&fGamPt    ,"fGamPt/D");//! fGamPt
    fTreeGJ->Branch("fGamLambda0" ,&fGamLambda0  ,"fGamLambda0/D");
    fTreeGJ->Branch("fGamNLM"     ,&fGamNLM      ,"fGamNLM/I");
    fTreeGJ->Branch("fGamSumPtCh" ,&fGamSumPtCh  ,"fGamSumPtCh/D");
    fTreeGJ->Branch("fGamNtracks" ,&fGamNtracks  ,"fGamNtracks/I");
    fTreeGJ->Branch("fGamTime"    ,&fGamTime     ,"fGamTime/D");
    fTreeGJ->Branch("fGamNcells"  ,&fGamNcells   ,"fGamNcells/I");
    fTreeGJ->Branch("fGamEta"     ,&fGamEta      ,"fGamEta/D");
    fTreeGJ->Branch("fGamPhi"     ,&fGamPhi      ,"fGamPhi/D");
    fTreeGJ->Branch("fGamSumPtNeu",&fGamSumPtNeu ,"fGamSumPtNeu/D");
    fTreeGJ->Branch("fGamNclusters" ,&fGamNclusters  ,"fGamNclusters/I");
    fTreeGJ->Branch("fGamAvEne"   ,&fGamAvEne    ,"fGamAvEne/D");
    fTreeGJ->Branch("fJetPhi"     ,&fJetPhi      ,"fJetPhi/D");
    fTreeGJ->Branch("fJetEta"     ,&fJetEta      ,"fJetEta/D");
    fTreeGJ->Branch("fJetPt"      ,&fJetPt       ,"fJetPt/D");
    fTreeGJ->Branch("fJetBkgChEne",&fJetBkgChEne ,"fJetBkgChEne/D");
    fTreeGJ->Branch("fJetArea"    ,&fJetArea     ,"fJetArea/D");
    fTreeGJ->Branch("fJetNtracks" ,&fJetNtracks  ,"fJetNtracks/I");
    fTreeGJ->Branch("fJetNtracks1" ,&fJetNtracks1  ,"fJetNtracks1/I");
    fTreeGJ->Branch("fJetNtracks2" ,&fJetNtracks2  ,"fJetNtracks2/I");
    fTreeGJ->Branch("fJetRho" ,&fJetRho  ,"fJetRho/D");
    fTreeGJ->Branch("fEventNumber",&fEventNumber ,"fEventNumber/I");
    fTreeGJ->Branch("fNtracks"    ,&fNtracks     ,"fNtracks/I");
    fTreeGJ->Branch("fZvertex"    ,&fZvertex     ,"fZvertex/D");
    fTreeGJ->Branch("fCentrality" ,&fCentrality  ,"fCentrality/D");
    fTreeGJ->Branch("fIso" ,&fIso  ,"fIso/O");
    fTreeGJ->Branch("fGamRho" ,&fGamRho  ,"fGamRho/D");

    fTreeGJ->Branch("fMCGamPt"         ,&fMCGamPt         ,"fMCGamPt/D"        );
    fTreeGJ->Branch("fMCGamEta"        ,&fMCGamEta        ,"fMCGamEta/D"       );
    fTreeGJ->Branch("fMCGamPhi"        ,&fMCGamPhi        ,"fMCGamPhi/D"       );
    fTreeGJ->Branch("fMCPartonType"    ,&fMCPartonType    ,"fMCPartonType/I"   );
    fTreeGJ->Branch("fMCJetPt"         ,&fMCJetPt         ,"fMCJetPt/D"        );
    fTreeGJ->Branch("fMCJetChPt"       ,&fMCJetChPt       ,"fMCJetChPt/D"      );
    fTreeGJ->Branch("fMCJet150Pt"      ,&fMCJet150Pt      ,"fMCJet150Pt/D"     );
    fTreeGJ->Branch("fMCJetCh150Pt"    ,&fMCJetCh150Pt    ,"fMCJetCh150Pt/D"   );
    fTreeGJ->Branch("fMCJetNPart"      ,&fMCJetNPart      ,"fMCJetNPart/I"     );
    fTreeGJ->Branch("fMCJetChNPart"    ,&fMCJetChNPart    ,"fMCJetChNPart/I"   );
    fTreeGJ->Branch("fMCJet150NPart"   ,&fMCJet150NPart   ,"fMCJet150NPart/I"  );
    fTreeGJ->Branch("fMCJetCh150NPart" ,&fMCJetCh150NPart ,"fMCJetCh150NPart/I");
    fTreeGJ->Branch("fMCJetEta"        ,&fMCJetEta        ,"fMCJetEta/D"       );
    fTreeGJ->Branch("fMCJetPhi"        ,&fMCJetPhi        ,"fMCJetPhi/D"       );
    fTreeGJ->Branch("fMCJetChEta"      ,&fMCJetChEta      ,"fMCJetChEta/D"     );
    fTreeGJ->Branch("fMCJetChPhi"      ,&fMCJetChPhi      ,"fMCJetChPhi/D"     );
    fTreeGJ->Branch("fMCJet150Eta"     ,&fMCJet150Eta     ,"fMCJet150Eta/D"    );
    fTreeGJ->Branch("fMCJet150Phi"     ,&fMCJet150Phi     ,"fMCJet150Phi/D"    );
    fTreeGJ->Branch("fMCJetCh150Eta"   ,&fMCJetCh150Eta   ,"fMCJetCh150Eta/D"  );
    fTreeGJ->Branch("fMCJetCh150Phi"   ,&fMCJetCh150Phi   ,"fMCJetCh150Phi/D"  );  

    fTreeGJ->Branch("fMCJetCh150ConePt"    ,&fMCJetCh150ConePt    ,"fMCJetCh150ConePt/D"  );    
    fTreeGJ->Branch("fMCJetCh150ConeNPart" ,&fMCJetCh150ConeNPart ,"fMCJetCh150ConeNPart/I");
    fTreeGJ->Branch("fMCJetCh150ConeEta"   ,&fMCJetCh150ConeEta   ,"fMCJetCh150ConeEta/D"  );   
    fTreeGJ->Branch("fMCJetCh150ConePhi"   ,&fMCJetCh150ConePhi   ,"fMCJetCh150ConePhi/D"  );   
    
    outputContainer->Add(fTreeGJ);
  }

  return outputContainer;
}

//_______________________________________________________
/// Initialize the parameters of the analysis.
//_______________________________________________________
void AliAnaParticleJetFinderCorrelation::InitParameters()
{
  SetInputAODName("PWG4Particle");
  AddToHistogramsName("AnaJetFinderCorr_");

  fDeltaPhiMinCut    = 2.6 ;
  fDeltaPhiMaxCut    = 3.7 ; 
  fRatioMaxCut       = 1.2 ; 
  fRatioMinCut       = 0.3 ; 
  fConeSize          = 0.4 ;
  fPtThresholdInCone = 0.5 ;
  fUseJetRefTracks   = kFALSE ;
  fMakeCorrelationInHistoMaker = kFALSE ;
  fSelectIsolated = kTRUE;
  fJetConeSize = 0.4 ;
  fJetMinPt = 15. ; //GeV/c
  fJetMinPtBkgSub = -100. ;//GeV/c
  fJetAreaFraction = 0.6 ;
  fJetBranchName = "jets";
  fBkgJetBranchName = "jets";
  fGammaConeSize = 0.4;
  fUseBackgroundSubtractionGamma = kFALSE;
  fSaveGJTree = kTRUE;
  fMostEnergetic = kFALSE;
  fMostOpposite = kTRUE;
  fUseHistogramJetBkg = kTRUE;
  fUseHistogramTracks = kTRUE;
  fUseHistogramJetTracks = kTRUE;
  fMCStudies = kFALSE;
}

//__________________________________________________________________
/// Input for jets is TClonesArray *aodRecJets.
/// \return the index of the jet that is opposite to the particle.
//__________________________________________________________________
Int_t  AliAnaParticleJetFinderCorrelation::SelectJet(AliAODPWG4Particle * particle, TClonesArray *aodRecJets)
{
  Double_t particlePt=particle->Pt();
  if(fUseBackgroundSubtractionGamma) {
      particlePt-=(fGamRho*particle->GetNCells());

//    Int_t clusterIDtmp = particle->GetCaloLabel(0) ;
//    Int_t nCells=0;
//    AliVCluster *cluster=0;
//    if(!(clusterIDtmp<0) ){
//      Int_t iclustmp = -1;
//      TObjArray* clusters = GetEMCALClusters();
//      cluster = FindCluster(clusters,clusterIDtmp,iclustmp);
//      nCells = cluster->GetNCells();
//    }
//    particlePt-=(fGamRho*nCells);
  }
  if(particlePt<=0) {
    //printf("Particle with negative  or 0 pt\n");
    return -2;
  }
  
  Int_t njets = aodRecJets->GetEntriesFast();
  AliAODJet * jet = 0 ;
  Int_t index = -1;
  Double_t dphiprev= 10000;
  Double_t particlePhi=particle->Phi();
  Double_t deltaPhi=-10000.;// in the range (0; 2*pi)
  Double_t jetPt=0.;
  
  Bool_t photonOnlyOnce=kTRUE;  

  for(Int_t ijet = 0; ijet < njets ; ijet++)
  {
    jet = dynamic_cast<AliAODJet*>(aodRecJets->At(ijet));
      
    if(!jet)
    {
      AliInfo("Jet not in container");
      continue;
    }
      
    fhCuts2->Fill(2., GetEventWeight());
      
    jetPt=jet->Pt();
    if(jetPt<fJetMinPt) continue;
      
    fhCuts2->Fill(3., GetEventWeight());
    //put jet eta requirement here |eta_jet|<0.9-jet_cone_size
    if(TMath::Abs(jet->Eta()) > (0.9 - fJetConeSize) ) continue;
      
    fhCuts2->Fill(4., GetEventWeight());
      
    if(jet->EffectiveAreaCharged()<fJetAreaFraction*TMath::Pi()*fJetConeSize*fJetConeSize) continue;
      
    fhCuts2->Fill(5., GetEventWeight());
      
    if(fBackgroundJetFromReader )
    {
      jetPt-= (fJetRho * jet->EffectiveAreaCharged() );
    }

    if(jetPt<fJetMinPtBkgSub) continue;
      
    fhCuts2->Fill(6., GetEventWeight());
      
    //printf("jet found\n");
    Double_t deltaPhi0pi  = TMath::Abs(particle->Phi()-jet->Phi());
    Double_t ratio = jetPt/particlePt;
    
    deltaPhi = particlePhi - jet->Phi() ;
    if ( deltaPhi0pi > TMath::Pi() ) deltaPhi0pi = 2. * TMath::Pi() - deltaPhi0pi ;
    if(deltaPhi<0) deltaPhi +=(TMath::Pi()*2.);
    
    // new histogram for Leticia x-check
    // isolated photon + jet(s)
    if(OnlyIsolated() && !particle->IsIsolated() && 
       (deltaPhi > fDeltaPhiMinCut) && (deltaPhi < fDeltaPhiMaxCut) )
    {
      //fill 1D photon + 2D photon+jets
      if(photonOnlyOnce)
      {
	    fhGamPtPerTrig->Fill(particlePt, GetEventWeight());
	    photonOnlyOnce=kFALSE;
      }
        
      fhPtGamPtJet->Fill(particlePt, jetPt, GetEventWeight());
    }
    

    fhDeltaPtBefore ->Fill(particlePt, particlePt - jetPt          , GetEventWeight());
    fhDeltaPhiBefore->Fill(particlePt, deltaPhi                    , GetEventWeight());
    fhDeltaEtaBefore->Fill(particlePt, particle->Eta() - jet->Eta(), GetEventWeight());
    fhPtRatioBefore ->Fill(particlePt, jetPt/particlePt            , GetEventWeight());
    fhPtBefore      ->Fill(particlePt, jetPt                       , GetEventWeight());
    
    fhDeltaPhi0PiCorrectBefore->Fill(particlePt, deltaPhi0pi , GetEventWeight());//good
    
    AliDebug(5,Form("Jet %d, Ratio pT %2.3f, Delta phi %2.3f",ijet,ratio,deltaPhi));
    
    if((deltaPhi > fDeltaPhiMinCut) && (deltaPhi < fDeltaPhiMaxCut) &&
       (ratio > fRatioMinCut) && (ratio < fRatioMaxCut)  &&
       (TMath::Abs(deltaPhi-TMath::Pi()) < TMath::Abs(dphiprev-TMath::Pi())  )){//In case there is more than one jet select the most opposite.
      dphiprev = deltaPhi;
      index = ijet ;
    }//Selection cuts
  }//AOD Jet loop
  
  return index ;
}

//__________________________________________________________________
/// Particle-Jet Correlation Analysis, fill AODs.
//__________________________________________________________________
void  AliAnaParticleJetFinderCorrelation::MakeAnalysisFillAOD()
{
  // Get the event, check if there are AODs available, if not, skip this analysis
  AliAODEvent * event = NULL;
  
  //  cout<<"GetReader()->GetOutputEvent() "<<GetReader()->GetOutputEvent()<<endl;
  //  cout<<"GetReader()->GetDataType() "<<GetReader()->GetDataType() <<endl;
  //  cout<<"AliCaloTrackReader::kAOD "<<AliCaloTrackReader::kAOD<<endl;
  //  cout<<"GetReader()->GetInputEvent() "<<GetReader()->GetInputEvent()<<endl;
  
  if(GetReader()->GetOutputEvent())
  {
    event = dynamic_cast<AliAODEvent*>(GetReader()->GetOutputEvent());
  }
  else if(GetReader()->GetDataType() == AliCaloTrackReader::kAOD)
  {
    event = dynamic_cast<AliAODEvent*>(GetReader()->GetInputEvent());
  }
  else
  {
    AliDebug(1,"There are no jets available for this analysis");
    return;
  }
  
  if(!GetInputAODBranch() || !event)
  {
    AliFatal(Form("No input particles in AOD with name branch < %s > \n",
                  GetInputAODName().Data()));
    return; // Trick coverity
  }
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation"))
  {
    AliFatal(Form("Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s>",
                  GetInputAODBranch()->GetClass()->GetName()));
    return; // Trick coverity
  }
  
  //
  // non-standard jets
  //
  Int_t nJets=-1;
  TClonesArray *aodRecJets = 0;
  //if(IsNonStandardJetFromReader()){//jet branch from reader
  AliDebug(3,Form("GetNonStandardJets function (from reader) is called"));
  aodRecJets = GetNonStandardJets();
  AliDebug(3,Form("aodRecJets %p",aodRecJets));
  if(aodRecJets==0x0)
    {
      if(GetDebug() > 3) event->Print();
      AliFatal("List of jets is null");
      return;
    }
  nJets=aodRecJets->GetEntries();
  AliDebug(3,Form("nJets %d",nJets));
  //}
  //end of new part
  
  if(nJets==0) {
    //printf("Why number of jets = 0? Check what type of collision it is. If PbPb -problem.\n");
    return;
  }
  
  //
  // Background jets
  //
  AliAODJetEventBackground* aodBkgJets = 0;
  if(IsBackgroundJetFromReader()){//jet branch from reader
   AliDebug(3,"GetBackgroundJets function is called");
    aodBkgJets = GetBackgroundJets();
    AliDebug(3,Form("aodBkgJets %p",aodBkgJets));
    if(aodBkgJets==0x0){
      if(GetDebug() > 3) event->Print();
      AliFatal("No jet background found\n");
      return; // Trick coverity        
    }
    if(GetDebug() > 3) aodBkgJets->Print("c");
  }
  
  Double_t rhoEvent=0.;
  if(aodBkgJets && IsBackgroundJetFromReader() )
  {
    rhoEvent = aodBkgJets->GetBackground(2);//hardcoded
  }
    
  fJetRho = rhoEvent;
  
  
  Int_t ntrig =  GetInputAODBranch()->GetEntriesFast() ;
  
  AliDebug(3,"Begin jet finder  correlation analysis, fill AODs");
  AliDebug(3,Form("In particle branch aod entries %d\n", ntrig));
  AliDebug(3,Form("In standard jet branch aod entries %d\n", event->GetNJets()));
  AliDebug(3,Form("In non standard jet branch aod entries %d\n", nJets));
  
  //if(nJets==0)   return;//to speed up
  //  cout<<"ntrig po return "<<ntrig<<endl;
  
  //
  //calculate average cell energy without most energetic photon
  //
  Double_t medianPhotonRho=0.;
  //TObjArray* clusters = GetEMCALClusters();
  //Int_t clusterIDtmp;
  //Int_t iclustmp = -1;
  //AliVCluster *cluster=0;
  
  if(IsBackgroundSubtractionGamma()){
    //
    // Find most energetic photon without background subtraction
    //
    Double_t maxPt=0.;
    Int_t maxIndex=-1;
    AliAODPWG4ParticleCorrelation* particlecorr =0;
    for(Int_t iaod = 0; iaod < ntrig ; iaod++){
      particlecorr =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
      if(particlecorr->Pt() > maxPt) {
        maxPt = particlecorr->Pt();
        maxIndex = iaod;
      }
    }
    
    //
    // calculate background energy per cell
    //
    Int_t numberOfcells=0;
    if(ntrig>1){
      Double_t *photonRhoArr=new Double_t[ntrig-1];
      Int_t photonRhoArrayIndex=0;
      for(Int_t iaod = 0; iaod < ntrig ; iaod++){
        particlecorr =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
        if(iaod==maxIndex) continue;
//        clusterIDtmp = particlecorr->GetCaloLabel(0) ;
//        if(clusterIDtmp < 0) continue;
//        cluster = FindCluster(clusters,clusterIDtmp,iclustmp);
//        photonRhoArr[photonRhoArrayIndex]=particlecorr->Pt()/ cluster->GetNCells();
//        numberOfcells+=cluster->GetNCells();
        photonRhoArr[photonRhoArrayIndex]=particlecorr->Pt()/ particlecorr->GetNCells();
        numberOfcells+=particlecorr->GetNCells();

        photonRhoArrayIndex++;
      }
      if(photonRhoArrayIndex>0) medianPhotonRho=TMath::Median(photonRhoArrayIndex,photonRhoArr);
      delete [] photonRhoArr;
    }
  }//end of if background calculation for gamma
  
  fGamRho = medianPhotonRho;
  
  //
  //take most energetic photon and most energetic jet and combine
  //
  if(fMostEnergetic){
    //
    // most energetic photon with background subtraction
    //
    Double_t mostEnePhotonPt=0.;
    Int_t indexMostEnePhoton=-1;
    AliAODPWG4ParticleCorrelation* particle =0;
    Double_t ptCorrect=0.;
//    Int_t nCells=0;
    for(Int_t iaod = 0; iaod < ntrig ; iaod++){
      particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
//      clusterIDtmp = particle->GetCaloLabel(0) ;
//      if(!(clusterIDtmp<0)){
//        cluster = FindCluster(clusters,clusterIDtmp,iclustmp);
//        nCells = cluster->GetNCells();
//      }
//      ptCorrect = particle->Pt() - medianPhotonRho * nCells;
      ptCorrect = particle->Pt() - medianPhotonRho * particle->GetNCells();
      
      if( ptCorrect > mostEnePhotonPt ){
        mostEnePhotonPt = ptCorrect;
        indexMostEnePhoton = iaod ;
      }
    }
    //    printf ("AliAnaParticleJetFinderCorrelation::MakeAnalysisFillAOD() - Photon with index %d selected \n",indexMostEnePhoton);
    //
    // most energetic jet with background subtraction
    //
    Double_t mostEneJetPt=0.;
    Int_t indexMostEneJet=-1;
    AliAODJet * jet = 0 ;
    //Double_t ptCorrect=0.;
    for(Int_t ijet = 0; ijet < nJets ; ijet++){
      jet = dynamic_cast<AliAODJet*>(aodRecJets->At(ijet));
      if(!jet) continue;
      if(jet->Pt()<fJetMinPt) continue;
      if(TMath::Abs(jet->Eta()) > (0.9 - fJetConeSize) ) continue;
      if(jet->EffectiveAreaCharged()<fJetAreaFraction*TMath::Pi()*fJetConeSize*fJetConeSize) continue;
      ptCorrect = jet->Pt() - rhoEvent * jet->EffectiveAreaCharged();
      if(ptCorrect > mostEneJetPt){
        mostEneJetPt = ptCorrect;
        indexMostEneJet = ijet;
      }
    }
    //    printf ("AliAnaParticleJetFinderCorrelation::MakeAnalysisFillAOD() - Jet with index %d selected \n",indexMostEneJet);
    
    //
    // assign jet to photon
    //
    if(indexMostEneJet>=0 && indexMostEnePhoton>=0){
      particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(indexMostEnePhoton));
      jet = dynamic_cast<AliAODJet*>(aodRecJets-> At(indexMostEneJet));
      if(jet)particle->SetRefJet(jet);
    }
  }//end of take most energetic photon and most ene. jet after bkg subtraction
  
  if(fMostOpposite){
    //Bool_t isJetFound=kFALSE;//new
    //Loop on stored AOD particles, trigger
    for(Int_t iaod = 0; iaod < ntrig ; iaod++){
      AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
      
      //Correlate with jets
      Int_t ijet = SelectJet(particle,aodRecJets);//input for jets is TClonesArray
      if(ijet > -1){
        //isJetFound=kTRUE;
        AliDebug(2,Form("Jet with index %d selected",ijet));
        AliAODJet *jet = dynamic_cast<AliAODJet*>(aodRecJets-> At(ijet));
        if(jet)particle->SetRefJet(jet);
        //printf("Most opposite found\n");
      }
    } // input aod loop
    //  if(GetReader()->WriteDeltaAODToFile() && isJetFound) WriteJetsToOutputBranch(aodRecJets);
  }//end of take most opposite photon and jet after bkg subtraction
  
  AliDebug(1," End fill AODs \n");
} 

//__________________________________________________________________
/// Particle-Jet Correlation Analysis, fill histograms.
//__________________________________________________________________
void  AliAnaParticleJetFinderCorrelation::MakeAnalysisFillHistograms()
{
  AliDebug(3,"I use MakeAnalysisFillHistograms");
  AliDebug(3,Form("ntracks before iso %d\n",GetCTSTracks()->GetEntriesFast()));

  // Get the event, check if there are AODs available, if not, skip this analysis
  AliAODEvent * event = NULL;
  
  //printf("GetReader()->GetOutputEvent() %d\n",GetReader()->GetOutputEvent() );
  //printf("GetReader()->GetDataType() %d\n",GetReader()->GetDataType() );
  //printf("AliCaloTrackReader::kAOD %d\n",AliCaloTrackReader::kAOD );
  //printf("GetReader()->GetInputEvent() %d\n",GetReader()->GetInputEvent() );
  
  if(GetReader()->GetOutputEvent())
  {
    event = dynamic_cast<AliAODEvent*>(GetReader()->GetOutputEvent());
  }
  else if(GetReader()->GetDataType() == AliCaloTrackReader::kAOD)
  {
    event = dynamic_cast<AliAODEvent*>(GetReader()->GetInputEvent());
  }
  else
  {
    AliDebug(3,"There are no jets available for this analysis");
    return;
  }
  
  if(!GetInputAODBranch() || !event)
  {
    AliFatal(Form("No input particles in AOD with name branch < %s >",
                  GetInputAODName().Data()));
    return; // Trick coverity        
  }
  
  Int_t nJets=-1;
  TClonesArray *aodRecJets = 0;
  //if(IsNonStandardJetFromReader()){//branch read in reader from reader
  AliDebug(3,"GetNonStandardJets function (from reader) is called");
  aodRecJets = GetNonStandardJets();
  if(aodRecJets==0x0)
  {
    if(GetDebug() > 3) event->Print();
    AliFatal("Jets container not found\n");
    return; // trick coverity
  }
  nJets=aodRecJets->GetEntries();
  //}
  if(nJets==0)
  {
    //    printf("Why number of jets = 0? Check what type of collision it is. If PbPb -problem.\n");
    GetReader()->FillInputNonStandardJets();
    aodRecJets = GetNonStandardJets();
    if(aodRecJets) nJets=aodRecJets->GetEntries();
    //    printf("nJets = %d\n",nJets);
    return;
  }
  
  //
  // Background jets
  //
  AliAODJetEventBackground* aodBkgJets = 0;
  if(IsBackgroundJetFromReader()){//jet branch from reader
    AliDebug(3,"GetBackgroundJets function is called");
    aodBkgJets = GetBackgroundJets();
    AliDebug(3,Form("aodBkgJets %p",aodBkgJets));
    if(aodBkgJets==0x0)
    {
      if(GetDebug() > 3) event->Print();
      AliFatal("No jet background container found");
      return; // trick coverity  
    }
    if(GetDebug() > 3) aodBkgJets->Print("c");
  }
  
  //
  // only background jets informations
  //
  //  Float_t pTback = 0;
  Double_t rhoEvent=0.;
  if(aodBkgJets) {
    if(IsBackgroundJetFromReader() ) rhoEvent = aodBkgJets->GetBackground(2);
    if(IsHistogramJetBkg())
    {
      fhJetRhoVsCentrality->Fill(GetEventCentrality(), rhoEvent, GetEventWeight());
        
      for(Int_t i=0;i<4;i++)
      {
        fhBkgJetBackground[i]->Fill(aodBkgJets->GetBackground(i), GetEventWeight());
        fhBkgJetSigma     [i]->Fill(aodBkgJets->GetSigma(i)     , GetEventWeight());
        fhBkgJetArea      [i]->Fill(aodBkgJets->GetMeanarea(i)  , GetEventWeight());
      }
    }//end of if fill HistogramJetBkg
  }//end if aodBkgJets exists
  
  //
  // only track information
  //
  Int_t nCTSTracks = GetCTSTracks()->GetEntriesFast();
  AliAODTrack *aodtrack;
  Int_t itrack = 0;
  if(IsHistogramTracks())
  {
    Double_t sumTrackPt=0;
    for(itrack = 0; itrack < nCTSTracks ; itrack++)
    {
      aodtrack = dynamic_cast <AliAODTrack*>(GetCTSTracks()->At(itrack));
      
      if(!aodtrack) continue;
        
      fhTrackPhiVsEta->Fill(aodtrack->Phi(), aodtrack->Eta(), GetEventWeight());
        
      sumTrackPt+=aodtrack->Pt();
    }
      
    if(nCTSTracks)
      fhTrackAveTrackPt->Fill(sumTrackPt/nCTSTracks, GetEventWeight());
  }//end if IsHistogramTracks
  
  //
  // only jet informations
  //
  AliAODJet * jettmp = 0 ;
  Double_t jetPttmp = 0.;
  Double_t var1tmp = 0.;
  Double_t var2tmp = 0.;
  //  fhJetNjet->Fill(nJets, GetEventWeight());
  Double_t ptMostEne=0;
  //  Int_t indexMostEne=-1;
  Int_t nJetsOverThreshold[10]={nJets,0,0,0,0,0,0,0,0,0};
  Int_t iCounter=0;
  Double_t sumJetTrackPt=0.;
  Int_t sumNJetTrack=0;
  Int_t nTracksInJet=0;
  Int_t itrk=0;
  for(Int_t ijet = 0; ijet < nJets ; ijet++)
  {
    jettmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ijet));
      
    if(!jettmp) continue;
      
    fhJetPtBefore->Fill(jettmp->Pt(), GetEventWeight());
      
    jetPttmp  = jettmp->Pt() - rhoEvent * jettmp->EffectiveAreaCharged();//<<---changed here
    
    //calculate number of tracks above 1,2,3,4,5 GeV in jet
    AliVTrack* jettrack = 0x0 ;
    
    Int_t nTrackThrGeV[5]={0,0,0,0,0};
    
    nTracksInJet=(jettmp->GetRefTracks())->GetEntriesFast();
      
    fhJetNparticlesInJet->Fill(nTracksInJet, GetEventWeight());
      
    if(nTracksInJet==0) continue;
      
    sumNJetTrack+=nTracksInJet;
      
    for(itrack=0;itrack<nTracksInJet;itrack++)
    {
      jettrack=(AliVTrack *) ((jettmp->GetRefTracks())->At(itrack));
        
      if(!jettrack) continue;
      
      fhJetDeltaEtaDeltaPhi->Fill(jettmp->Eta()-jettrack->Eta(), jettmp->Phi()-jettrack->Phi(), GetEventWeight());
        
      sumJetTrackPt+=jettrack->Pt();
      
      if(IsHistogramJetTracks())
      {
        if(jettrack->Pt()>1.) nTrackThrGeV[0]++;
        if(jettrack->Pt()>2.) nTrackThrGeV[1]++;
        if(jettrack->Pt()>3.) nTrackThrGeV[2]++;
        if(jettrack->Pt()>4.) nTrackThrGeV[3]++;
        if(jettrack->Pt()>5.) nTrackThrGeV[4]++;
      }//end of if IsHistogramJetTracks
    }//end of jet track loop
    //printf("jetPt %f ntrks %d ntrks>1,2,3,4,5GeV in jet %d, %d, %d, %d, %d\n",jetPttmp,nTracksInJet,nTrackThrGeV[0],nTrackThrGeV[1],nTrackThrGeV[2],nTrackThrGeV[3],nTrackThrGeV[4]);
    
    for(itrack = 0; itrack < nCTSTracks ; itrack++)
    {
      aodtrack = dynamic_cast <AliAODTrack*>(GetCTSTracks()->At(itrack));
        
      if(aodtrack) fhJetDeltaEtaDeltaPhiAllTracks->Fill(jettmp->Eta()-aodtrack->Eta(),
                                                        jettmp->Phi()-aodtrack->Phi(), GetEventWeight());
    }
    
    
    if(IsHistogramJetTracks()){
      fhJetNtracksInJetAboveThr[0]->Fill(jetPttmp, nTracksInJet, GetEventWeight());//all jets
      
      for(itrk=0;itrk<5;itrk++) {
        fhJetNtracksInJetAboveThr[itrk+1]->Fill(jetPttmp, nTrackThrGeV[itrk], GetEventWeight());//all jets
        fhJetRatioNTrkAboveToNTrk[itrk]  ->Fill(jetPttmp, (Double_t)nTrackThrGeV[itrk]/(Double_t)nTracksInJet, GetEventWeight());//all jets
      }
      if(ijet==0){//most ene jet
        for(itrk=0;itrk<5;itrk++)
          fhJetNtrackRatioMostEne[itrk]->Fill(jetPttmp, (Double_t)nTrackThrGeV[itrk]/(Double_t)nTracksInJet, GetEventWeight());
      }
      if(jetPttmp>5){//jet with pt>5GeV/c
        for(itrk=0;itrk<5;itrk++)
          fhJetNtrackRatioJet5GeV[itrk]->Fill(jetPttmp, (Double_t)nTrackThrGeV[itrk]/(Double_t)nTracksInJet, GetEventWeight());
      }
      if(nTrackThrGeV[4]>0){//jet with leading particle pt>5GeV/c
        for(itrk=0;itrk<5;itrk++)
          fhJetNtrackRatioLead5GeV[itrk]->Fill(jetPttmp, (Double_t)nTrackThrGeV[itrk]/(Double_t)nTracksInJet, GetEventWeight());
      }
    }//end of if IsHistogramJetTracks
    
    
    Double_t effectiveChargedBgEnergy=(IsBackgroundJetFromReader()?rhoEvent * jettmp->EffectiveAreaCharged():jettmp->ChargedBgEnergy());
    
    
    fhJetChBkgEnergyVsArea->Fill(effectiveChargedBgEnergy, jettmp->EffectiveAreaCharged(), GetEventWeight());
    //if(jettmp->EffectiveAreaCharged()>0)
    fhJetRhoVsPt->Fill(jetPttmp, jettmp->ChargedBgEnergy()*jettmp->EffectiveAreaCharged(), GetEventWeight());
    
    if(jetPttmp>ptMostEne)
    {
      ptMostEne = jetPttmp;
      //indexMostEne=ijet;
    }
    if(jettmp->Pt()>=fJetMinPt)
      fhJetPtBeforeCut->Fill(jetPttmp, GetEventWeight());

    fhJetPt->Fill(jetPttmp, GetEventWeight());
    fhJetChBkgEnergyVsPt->Fill(jetPttmp,effectiveChargedBgEnergy, GetEventWeight());
    fhJetChAreaVsPt->Fill(jetPttmp,jettmp->EffectiveAreaCharged(), GetEventWeight());
    
    AliDebug(5,Form("ChargedBgEnergy %f EffectiveAreaCharged %f\n", jettmp->ChargedBgEnergy(),jettmp->EffectiveAreaCharged()));
   
    for(iCounter=1;iCounter<10;iCounter++)
    {
      if(jetPttmp>iCounter) nJetsOverThreshold[iCounter]++;
    }
    
    var1tmp  = jettmp->Phi();
    var2tmp  = jettmp->Eta();
      
    fhJetPhi->Fill(var1tmp, GetEventWeight());
    fhJetEta->Fill(var2tmp, GetEventWeight());
      
    fhJetPhiVsEta->Fill(var1tmp, var2tmp, GetEventWeight());
    fhJetEtaVsPt->Fill(jetPttmp, var2tmp, GetEventWeight());
      
    fhJetEtaVsNpartInJet->Fill(var2tmp, nTracksInJet, GetEventWeight());
      
    if(jetPttmp>0)
      fhJetEtaVsNpartInJetBkg->Fill(var2tmp, nTracksInJet, GetEventWeight());
    
  }//end of jet loop
    
  if(IsHistogramJetTracks())
  {
    if(sumNJetTrack>0)
    {
      //printf("average track pt %f\n",sumJetTrackPt/sumNJetTrack);
      fhJetAveTrackPt->Fill(sumJetTrackPt/sumNJetTrack, GetEventWeight());
    }
  }//end of if IsHistogramJetTracks
  
  
  fhJetPtMostEne->Fill(ptMostEne, GetEventWeight());
  for(iCounter=0;iCounter<10;iCounter++)
  {
    fhJetNjetOverPtCut[iCounter]->Fill(nJetsOverThreshold[iCounter], GetEventWeight());
    nJetsOverThreshold[iCounter]=0;
  }
  
  //end of only jet part
  
  //
  // Photons
  //
  Int_t ntrig   =  GetInputAODBranch()->GetEntriesFast() ;
  
  AliDebug(1,"Begin jet finder  correlation analysis, fill histograms");
  AliDebug(1,Form("In particle branch aod entries %d\n", ntrig));
  AliDebug(1,Form("In jet output branch aod entries %d\n", event->GetNJets()));
  
  fhNjetsNgammas->Fill(nJets, ntrig, GetEventWeight());
  
  //if(nJets==0)   return;//to speed up
  //printf("ntrig %d, njets %d\n",ntrig,nJets);
  
  //
  //Fill only temporary photon histograms
  //
  Double_t maxPt=0.;
  Int_t maxIndex=-1;
  Double_t tmpPt=0.;
  Double_t sumPt=0.;
  nJetsOverThreshold[0]=ntrig;
    
  for(Int_t iaod = 0; iaod < ntrig ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particlecorr =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    tmpPt = particlecorr->Pt();
    if(tmpPt>maxPt)
    {
      maxPt = tmpPt;
      maxIndex = iaod;
    }
    
    for(iCounter=1;iCounter<10;iCounter++)
    {
      if(tmpPt>iCounter) nJetsOverThreshold[iCounter]++;
    }
    sumPt+=tmpPt;
  }
  
  for(iCounter=0;iCounter<10;iCounter++)
  {
    fhPhotonNgammaOverPtCut[iCounter]->Fill(nJetsOverThreshold[iCounter], GetEventWeight());
    //    nJetsOverThreshold[iCounter]=0;
  }
  
  fGamAvEne=0;
  //TObjArray* clusters = GetEMCALClusters();
  //printf("calculate median bkg energy for photons ");
  Double_t medianPhotonRho=0.;
  //Int_t clusterID;
  //Int_t iclustmp = -1;
  Int_t numberOfcells=0;
  //AliVCluster *cluster = 0;
  if(ntrig>1){
    Double_t *photonRhoArr=new Double_t[ntrig-1];
    fhPhotonPtMostEne->Fill(maxPt, GetEventWeight());
    //    fhPhotonIndexMostEne->Fill(indexMaxPt, GetEventWeight());
    fhPhotonAverageEnergy       ->Fill(sumPt/ntrig, GetEventWeight());
    fhPhotonRatioAveEneToMostEne->Fill(sumPt/(ntrig*maxPt), GetEventWeight());
    fhPhotonAverageEnergyMinus1 ->Fill((sumPt-maxPt)/(ntrig-1), GetEventWeight());
    fGamAvEne=(sumPt-maxPt)/(ntrig-1);
    fhPhotonRatioAveEneMinus1ToMostEne->Fill((sumPt-maxPt)/((ntrig-1)*maxPt), GetEventWeight());
    
    Int_t counterGamma=0;
    Int_t counterGammaMinus1=0;
    
    Int_t photonRhoArrayIndex=0;
    //TObjArray* clusterstmp = GetEMCALClusters();
    for(Int_t iaod = 0; iaod < ntrig ; iaod++){
      AliAODPWG4ParticleCorrelation* particlecorr =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
      if( particlecorr->Pt() > sumPt/ntrig ) counterGamma++;
      if( particlecorr->Pt() > (sumPt-maxPt)/(ntrig-1) ) counterGammaMinus1++;
      
      if(iaod==maxIndex) continue;
//      clusterID = particlecorr->GetCaloLabel(0) ;
//      if(clusterID < 0) continue;
//      cluster = FindCluster(clusters,clusterID,iclustmp);
//      photonRhoArr[photonRhoArrayIndex]=particlecorr->Pt()/ cluster->GetNCells();
//      numberOfcells+=cluster->GetNCells();
      photonRhoArr[photonRhoArrayIndex]=particlecorr->Pt()/ particlecorr->GetNCells();
      numberOfcells+=particlecorr->GetNCells();
      photonRhoArrayIndex++;
    }
    if(photonRhoArrayIndex>0) medianPhotonRho=TMath::Median(photonRhoArrayIndex,photonRhoArr);
    delete [] photonRhoArr;
    fhPhotonNgammaMoreAverageToNgamma      ->Fill((Double_t)counterGamma / (Double_t)ntrig, GetEventWeight());
    fhPhotonNgammaMoreAverageMinus1ToNgamma->Fill((Double_t)counterGammaMinus1 / (Double_t)ntrig, GetEventWeight());
  }
  
  //printf("median = %f\n",medianPhotonRho);
    
  fhPhotonBkgRhoVsNtracks   ->Fill(GetCTSTracks()->GetEntriesFast(), medianPhotonRho, GetEventWeight());
  fhPhotonBkgRhoVsNclusters ->Fill(ntrig, medianPhotonRho, GetEventWeight());
  fhPhotonBkgRhoVsCentrality->Fill(GetEventCentrality(), medianPhotonRho, GetEventWeight());
  fhPhotonBkgRhoVsNcells    ->Fill(numberOfcells, medianPhotonRho, GetEventWeight());
  
  //AliVCluster *cluster2 = 0;
  Double_t photon2Corrected=0;
  Double_t sumPtTmp=0.;
  Double_t sumPtCorrectTmp=0.;
  AliVTrack* trackTmp = 0x0 ;
  TVector3 p3Tmp;
  
  for(Int_t iaod = 0; iaod < ntrig ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particlecorr =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
//    clusterID = particlecorr->GetCaloLabel(0) ;
//    if(clusterID < 0) continue;
//    cluster = FindCluster(clusters,clusterID,iclustmp);
//  Int_t ncells = cluster->GetNCells();
    Int_t ncells = particlecorr->GetNCells();
    fhPhotonPt                ->Fill(particlecorr->Pt(), GetEventWeight());
    fhPhotonPtCorrected       ->Fill(particlecorr->Pt() - ncells * medianPhotonRho, GetEventWeight());
    fhPhotonPtDiff            ->Fill(ncells * medianPhotonRho, GetEventWeight());
    fhPhotonPtDiffVsCentrality->Fill(GetEventCentrality(),ncells * medianPhotonRho, GetEventWeight());
    fhPhotonPtDiffVsNcells    ->Fill(numberOfcells,ncells * medianPhotonRho, GetEventWeight());
    fhPhotonPtDiffVsNtracks   ->Fill(GetCTSTracks()->GetEntriesFast(),ncells * medianPhotonRho, GetEventWeight());
    fhPhotonPtDiffVsNclusters ->Fill(ntrig,ncells * medianPhotonRho, GetEventWeight());
    fhPhotonPtCorrectedZoom   ->Fill(particlecorr->Pt() - ncells * medianPhotonRho, GetEventWeight());
    
    //test: sum_pt in the cone 0.3 for each photon
    //should be: random fake gamma from MB
    //is: each gamma for EMCEGA
    sumPtTmp=0.;
    sumPtCorrectTmp=0.;
    
    for(Int_t iaod2 = 0; iaod2 < ntrig ; iaod2++)
    {
      if(iaod==iaod2) continue;
        
      AliAODPWG4ParticleCorrelation* particlecorr2 =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod2));
//      clusterID = particlecorr2->GetCaloLabel(0) ;
//      if(clusterID < 0) continue;
//      cluster2 = FindCluster(clusters,clusterID,iclustmp);
//      photon2Corrected = particlecorr2->Pt() - cluster2->GetNCells() * medianPhotonRho;
      photon2Corrected = particlecorr2->Pt() - particlecorr2->GetNCells() * medianPhotonRho;

      //if(Pt()<0.5) continue; //<<hardcoded here //FIXME
      if( TMath::Sqrt((particlecorr->Eta()-particlecorr2->Eta())*(particlecorr->Eta()-particlecorr2->Eta()) +
                      (particlecorr->Phi()-particlecorr2->Phi())*(particlecorr->Phi()-particlecorr2->Phi()) )<fGammaConeSize ){//if(/*cone is correct*/){
        sumPtTmp+= particlecorr2->Pt();
        sumPtCorrectTmp+=photon2Corrected;
      }
    }
      
    fhPhotonSumPtInCone->Fill(sumPtTmp, GetEventWeight());
    fhPhotonSumPtCorrectInCone->Fill(sumPtCorrectTmp, GetEventWeight());
    
    //test: sum_pt in the cone 0.3 for each track
    //should be: random fake gamma from MB
    //is: each gamma for EMCEGA
    sumPtTmp=0.;
    for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++){
      trackTmp = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
      p3Tmp.SetXYZ(trackTmp->Px(),trackTmp->Py(),trackTmp->Pz());
      if( TMath::Sqrt((particlecorr->Eta()-p3Tmp.Eta())*(particlecorr->Eta()-p3Tmp.Eta()) +
                      (particlecorr->Phi()-p3Tmp.Phi())*(particlecorr->Phi()-p3Tmp.Phi()) )<fGammaConeSize ){
        sumPtTmp+=p3Tmp.Pt();
      }
    }//end of loop over tracks
      
    fhPhotonSumPtChargedInCone->Fill(sumPtTmp, GetEventWeight());
  }
  
  //End of Fill temporary photon histograms
  
  //
  // Apply background subtraction for photons
  //
  fGamRho = medianPhotonRho;
  if(!IsBackgroundSubtractionGamma()) medianPhotonRho=0;
  
  
  //
  //Get vertex for cluster momentum calculation <<----new here
  //
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  fZvertex = vertex[2];
  
  //
  //Loop on stored AOD particles, trigger
  //
  for(Int_t iaod = 0; iaod < ntrig ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particlecorr =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    fhCuts ->Fill(0., GetEventWeight());
    fhCuts2->Fill(0., (Double_t)nJets * GetEventWeight());
    AliDebug(1,Form("OnlyIsolated %d  !particlecorr->IsIsolated() %d \n",OnlyIsolated(), !particlecorr->IsIsolated()));
    
    if(OnlyIsolated() && !particlecorr->IsIsolated()) continue;
      
    fhCuts ->Fill(1., GetEventWeight());
    fhCuts2->Fill(1., nJets * GetEventWeight());
    
    if(nJets>0)
    {
      fhCuts->Fill(2., GetEventWeight());
    }
    
    //Recover the jet correlated, found previously.
    AliAODJet	* jet = particlecorr->GetJet();
    //If correlation not made before, do it now.
    if(fMakeCorrelationInHistoMaker){
      //Correlate with jets
      Int_t ijet = SelectJet(particlecorr,aodRecJets);//input for jets is TClonesArray
      if(ijet > -1)
      {
        AliDebug(1,Form("Jet with index %d selected \n",ijet));
        //jet = event->GetJet(ijet);
        jet = dynamic_cast<AliAODJet*>(aodRecJets-> At(ijet));
        
       if(jet) particlecorr->SetRefJet(jet);
        
      }
    }
    
    if (!jet) continue ;
      
    fhCuts ->Fill(3., GetEventWeight());
    fhCuts2->Fill(7., GetEventWeight());

    //check MC genereted information
    if(fMCStudies) FindMCgenInfo();

    //
    //Fill Correlation Histograms
    //
//    clusterID = particlecorr->GetCaloLabel(0) ;
//    if(!(clusterID<0)){
//      cluster = FindCluster(clusters,clusterID,iclustmp);
//      //fill tree variables
//      fGamNcells = cluster->GetNCells();
//    }
    
    fGamNcells = particlecorr->GetNCells();
    
    Double_t ptTrig = particlecorr->Pt() - medianPhotonRho * fGamNcells;//<<---changed here
    Double_t ptJet = jet->Pt() - rhoEvent * jet->EffectiveAreaCharged();//<<---changed here
    Double_t phiTrig = particlecorr->Phi();
    Double_t phiJet = jet->Phi();
    Double_t etaTrig = particlecorr->Eta();
    Double_t etaJet = jet->Eta();
    Double_t deltaPhi=phiTrig-phiJet;
    if(deltaPhi<0)deltaPhi+=(TMath::Pi()*2.);
    //printf("pT trigger %2.3f, pT jet %2.3f, Delta phi %2.3f, Delta eta %2.3f, Delta pT %2.3f, ratio %2.3f \n",
    //	ptTrig,ptJet, phiJet-phiTrig, etaJet-etaTrig, ptTrig-ptJet, ptJet/ptTrig);
    fhDeltaPt ->Fill(ptTrig, ptTrig-ptJet, GetEventWeight());
    //    fhDeltaPhi->Fill(ptTrig, phiTrig-phiJet);//need to be shifted by 2pi
    
    fhDeltaPhiCorrect->Fill(ptTrig, deltaPhi, GetEventWeight());//correct
    
    Double_t deltaPhiCorrect = TMath::Abs( particlecorr->Phi() - jet->Phi() );
    if ( deltaPhiCorrect > TMath::Pi() ) deltaPhiCorrect = 2. * TMath::Pi() - deltaPhiCorrect ;
    fhDeltaPhi0PiCorrect->Fill(ptTrig, deltaPhiCorrect, GetEventWeight());
    
    fhDeltaEta->Fill(ptTrig, etaTrig-etaJet, GetEventWeight());
    fhPtRatio ->Fill(ptTrig, ptJet/ptTrig  , GetEventWeight());
    fhPt      ->Fill(ptTrig, ptJet         , GetEventWeight());
    
    fhSelectedJetPhiVsEta->Fill(phiJet, etaJet, GetEventWeight());
  
    fhSelectedJetChBkgEnergyVsPtJet->Fill(ptJet, (IsBackgroundJetFromReader()?rhoEvent * jet->EffectiveAreaCharged():jet->ChargedBgEnergy()), GetEventWeight());
      
    fhSelectedJetChAreaVsPtJet->Fill(ptJet, jet->EffectiveAreaCharged(), GetEventWeight());
    fhSelectedJetNjet         ->Fill(nJets, GetEventWeight());
    fhSelectedNtracks         ->Fill(GetCTSTracks()->GetEntriesFast(), GetEventWeight());//to be checked
    fhSelectedPhotonNLMVsPt   ->Fill(ptTrig, particlecorr->GetNLM(), GetEventWeight());
    
//    if(clusterID < 0 ){
//      fhSelectedPhotonLambda0VsPt->Fill(ptTrig, -1, GetEventWeight());
//      //fill tree variables
//      fGamLambda0  = -1;
//      fGamTime = -1;
//      fGamNcells = 0;
//      fGamSumPtNeu=0;
//    }
//    else
//    {
      //Int_t iclus = -1;
      TObjArray* clusters = GetEMCALClusters();
      //cluster = FindCluster(clusters,clusterID,iclustmp);
      Double_t lambda0=particlecorr->GetM02();
      fhSelectedPhotonLambda0VsPt->Fill(ptTrig, lambda0, GetEventWeight());
      //fill tree variables
      fGamLambda0  = lambda0;
      fGamTime = particlecorr->GetTime();
      //fGamNcells = cluster->GetNCells();
      
      fGamSumPtNeu=0;
      fGamNclusters=0;
      //TVector3 p3Tmp;
      //Double_t scalarProduct=0;
      //Double_t vectorLength=particlecorr->P();
      for(Int_t icalo=0; icalo <clusters->GetEntriesFast(); icalo++){
        AliVCluster* calo = (AliVCluster *) clusters->At(icalo);
        //if(clusterID==calo->GetID()) continue;//the same cluster as trigger
        calo->GetMomentum(fMomentum,vertex) ;//Assume that come from vertex in straight line
        //printf("min pt %f\n",GetMinPt());
        if(fMomentum.Pt()<GetMinPt()) continue; //<<hardcoded here //FIXME 0.5 check if correct
        p3Tmp.SetXYZ(fMomentum.Px(),fMomentum.Py(),fMomentum.Pz());
        //calculate sum pt in the cone
        if( TMath::Sqrt((particlecorr->Eta()-p3Tmp.Eta())*(particlecorr->Eta()-p3Tmp.Eta()) +
                        (particlecorr->Phi()-p3Tmp.Phi())*(particlecorr->Phi()-p3Tmp.Phi()) )<fGammaConeSize ){
          //scalarProduct = particlecorr->Px()*fMomentum.Px() + particlecorr->Py()*fMomentum.Py() + particlecorr->Pz()*fMomentum.Pz();
          //scalarProduct/=fMomentum.P();
          //scalarProduct/=vectorLength;
          //if(scalarProduct>TMath::Cos(0.3)) {//FIXME photon radius
          fGamSumPtNeu+=fMomentum.Pt();
          fGamNclusters++;
        }
      }
//    }
    
    //sum pt of charged tracks in the gamma isolation cone
    //starts here
    fGamSumPtCh=0;
    fGamNtracks=0;
    for(itrack = 0; itrack < nCTSTracks ; itrack++){
      aodtrack = dynamic_cast <AliAODTrack*>(GetCTSTracks()->At(itrack));
      if(!aodtrack) continue;
      fhSelectedTrackPhiVsEta->Fill(aodtrack->Phi(), aodtrack->Eta(), GetEventWeight());//fill histogram here
      //      if(aodtrack->Pt()<0.15) continue;//hardcoded
      if(aodtrack->Pt()<fPtThresholdInCone) continue;
      if(!aodtrack->IsHybridGlobalConstrainedGlobal()) continue;
      if(TMath::Sqrt((particlecorr->Phi() - aodtrack->Phi())*(particlecorr->Phi() - aodtrack->Phi()) +
                     (particlecorr->Eta() - aodtrack->Eta())*(particlecorr->Eta() - aodtrack->Eta()) ) <fGammaConeSize ) {
        fGamSumPtCh+=aodtrack->Pt();
        fGamNtracks++;
      }
    }
    //ends here
    
    //    for(Int_t itrack = 0; itrack < nCTSTracks ; itrack++){
    //      aodtrack = dynamic_cast <AliAODTrack*>(GetCTSTracks()->At(itrack));
    //      fhSelectedTrackPhiVsEta->Fill(aodtrack->Phi(), aodtrack->Eta(), GetEventWeight());
    //    }
    
    //
    // Background Fragmentation function
    //
    TVector3 gammaVector,jetVector;
    gammaVector.SetXYZ(particlecorr->Px(),particlecorr->Py(),particlecorr->Pz());
    jetVector.SetXYZ(jet->Px(),jet->Py(),jet->Pz());
    CalculateBkg(gammaVector,jetVector,vertex,1);//jet perp
    CalculateBkg(gammaVector,jetVector,vertex,2);//RC
    CalculateBkg(gammaVector,jetVector,vertex,3);//mid point
    CalculateBkg(gammaVector,jetVector,vertex,4);//gamma perp
    //CalculateBkg(gammaVector,jetVector,vertex,5);/test
    Double_t angleJetGam = gammaVector.Angle(jetVector);
    //printf("angleJetGam %f\n",angleJetGam*180/TMath::Pi());
    
    //
    // Fragmentation function
    //
    Float_t	 rad = 0, pt = 0, eta = 0, phi = 0;
    Int_t	 npartcone = 0;
    TVector3 p3;
    
    Int_t ntracks =  0;

    AliDebug(1,Form("fUseJetRefTracks %d"   ,fUseJetRefTracks   ));
    AliDebug(1,Form("jet->GetRefTracks() %p",jet->GetRefTracks()));
    AliDebug(1,Form("GetCTSTracks() %p"     ,GetCTSTracks()     ));
    
    if(!fUseJetRefTracks)
      ntracks =GetCTSTracks()->GetEntriesFast();
    else //If you want to use jet tracks from JETAN
      ntracks =  (jet->GetRefTracks())->GetEntriesFast();
    
    AliDebug(3,Form("ntracks %d\n",ntracks));
    AliVTrack* track = 0x0 ;
    for(Int_t ipr = 0;ipr < ntracks ; ipr ++ ){
      if(!fUseJetRefTracks)
        track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
      else //If you want to use jet tracks from JETAN
        track = (AliVTrack *) ((jet->GetRefTracks())->At(ipr));
      
      p3.SetXYZ(track->Px(),track->Py(),track->Pz());
      pt    = p3.Pt();
      eta  = p3.Eta();
      phi  = p3.Phi() ;
      if(phi < 0) phi+=TMath::TwoPi();
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt((eta-etaJet)*(eta-etaJet)+ (phi-phiJet)*(phi-phiJet));
      if(rad < fConeSize  && pt > fPtThresholdInCone)
      {
        //printf("charged in jet cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
        npartcone++;
        fhFFz ->Fill(ptTrig, pt/ptTrig, GetEventWeight());
        fhFFxi->Fill(ptTrig, TMath::Log(ptTrig/pt), GetEventWeight());
        fhFFpt->Fill(ptTrig, pt, GetEventWeight());
        
        //according to jet axis
        fhJetFFz ->Fill(ptJet, pt/ptJet, GetEventWeight());
        fhJetFFxi->Fill(ptJet, TMath::Log(ptJet/pt), GetEventWeight());
        fhJetFFpt->Fill(ptJet, pt, GetEventWeight());
        
        
        if(TMath::Cos(angleJetGam)<0 && ptJet!=0 && pt!=0 )
        {
          fhJetFFzCor ->Fill(ptJet, -pt*TMath::Cos(angleJetGam)/ptJet, GetEventWeight());
          fhJetFFxiCor->Fill(ptJet, TMath::Log(ptJet/(-pt*TMath::Cos(angleJetGam))), GetEventWeight());
        }
      }
    }//Tracks
    fhNTracksInCone->Fill(ptTrig, npartcone, GetEventWeight());
    //fill tree here for each photon-jet (isolation required usually)
    
    fGamPt      = ptTrig;
    //fGamLambda0  = ;//filled earlier
    fGamNLM      = particlecorr->GetNLM();
    //fGamSumPtCh  = ;//filled earlier
    //fGamTime     = particlecorr->GetTOF();//filled earlier
    //fGamNcells   = particlecorr->GetNCells();//filled earlier
    fGamEta      = etaTrig;
    fGamPhi      = phiTrig;
    //fGamSumPtNeu = ;//filled earlier
    //fGamNtracks  = ;//filled earlier
    //fGamNclusters= ;//filled earlier
    //fGamAvEne    = ;//filled earlier
    fJetPhi      = phiJet;
    fJetEta      = etaJet;
    fJetPt       = ptJet;
    fJetBkgChEne = (IsBackgroundJetFromReader()?rhoEvent * jet->EffectiveAreaCharged():jet->ChargedBgEnergy());
    fJetArea     = jet->EffectiveAreaCharged();
    fJetNtracks  = (jet->GetRefTracks())->GetEntriesFast();
    fEventNumber = 0;
    fNtracks     = GetCTSTracks()->GetEntriesFast();
    fCentrality  = GetEventCentrality();
    fIso         = particlecorr->IsIsolated();
    
    Int_t nTrk1GeV=0;
    Int_t nTrk2GeV=0;
    for(itrack=0;itrack < fJetNtracks;itrack++){
      track = (AliVTrack *) ((jet->GetRefTracks())->At(itrack));
      if(track->Pt()>1.) nTrk1GeV++;
      if(track->Pt()>2.) nTrk2GeV++;
    }
    
    fJetNtracks1 = nTrk1GeV;
    fJetNtracks2 = nTrk2GeV;
    
    if(fSaveGJTree) fTreeGJ->Fill();
  }//AOD trigger particle loop
  AliDebug(1,"End fill histograms");
}

//__________________________________________________________________
/// Print some relevant parameters set for the analysis
//__________________________________________________________________
void AliAnaParticleJetFinderCorrelation::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");

  printf("Phi trigger-jet        <     %3.2f\n", fDeltaPhiMaxCut) ; 
  printf("Phi trigger-jet        >     %3.2f\n", fDeltaPhiMinCut) ;
  printf("pT Ratio trigger/jet   <     %3.2f\n", fRatioMaxCut) ; 
  printf("pT Ratio trigger/jet   >     %3.2f\n", fRatioMinCut) ;
  printf("fConeSize              =     %3.2f\n", fConeSize) ; 
  printf("fPtThresholdInCone     =     %3.2f\n", fPtThresholdInCone) ;
  printf("fUseJetRefTracks	   =     %d\n",    fUseJetRefTracks) ;
  printf("fMakeCorrelationInHistoMaker	   =     %d\n",    fMakeCorrelationInHistoMaker) ;	
  printf("Isolated Trigger?  %d\n", fSelectIsolated) ;
  printf("Reconstructed jet cone size = %3.2f\n", fJetConeSize) ;
  printf("Reconstructed jet minimum pt before background subtraction = %3.2f\n", fJetMinPt) ;
  printf("Reconstructed jet minimum pt after background subtraction = %3.2f\n", fJetMinPtBkgSub) ;
  printf("Reconstructed jet minimum area fraction = %3.2f\n", fJetAreaFraction) ;

  //if(!fNonStandardJetFromReader){
  printf("fJetBranchName =   %s\n", fJetBranchName.Data()) ;
  //}
  if(!fBackgroundJetFromReader){
    printf("fBkgJetBranchName =   %s\n", fBkgJetBranchName.Data()) ;
  }

  printf("Isolation cone size = %3.2f\n", fGammaConeSize) ;
  printf("fUseBackgroundSubtractionGamma = %d\n",fUseBackgroundSubtractionGamma);
  printf("fSaveGJTree = %d\n",fSaveGJTree);
  printf("fMostEnergetic = %d\n",fMostEnergetic);
  printf("fMostOpposite = %d\n",fMostOpposite);

  printf("fUseHistogramJetBkg = %d\n",fUseHistogramJetBkg);
  printf("fUseHistogramTracks = %d\n",fUseHistogramTracks);
  printf("fUseHistogramJetTracks = %d\n",fUseHistogramJetTracks);
  printf("fMCStudies = %d\n",fMCStudies);
} 

//__________________________________________________________________
/// Calculate background for fragmentation function and fill histograms:
///  * 1. 90 degrees from jet axis in random place = perpendicular cone
///  * 2. Random cone not belonging to jet cone nor photon cone
///  * 3. In the middle point from jet and photon momentum vectors
///  * 4. 90 degrees from photon direction in random place = perpendicular cone 2
//__________________________________________________________________
void AliAnaParticleJetFinderCorrelation::CalculateBkg(TVector3 gamma, TVector3 jet,Double_t vertex[3],Int_t type=2)
{
  //
  // implementation of 2 works, 1 and 4 works
  //
  Double_t gammaPt  = gamma.Pt();
  Double_t gammaEta = gamma.Eta();
  Double_t gammaPhi = gamma.Phi();
  Double_t jetEta   = jet.Eta();
  Double_t jetPhi   = jet.Phi();

  //refference direction of background
  Double_t refEta=0.,refPhi=0.;
  Double_t rad = 0,rad2 = 0.;
  if(type==1){//perpendicular to jet axis
    //printf("vertex: %f %f %f \n",vertex[0],vertex[1],vertex[2]);

    Double_t xVar;
    Double_t yVar;
    Double_t newX=0.;
    Double_t newY=0.;
    Double_t newZ=0.;
    //Z axis vector: [jet.Px(),jet.Py(),jet.Pz()]
    Double_t jx=jet.Px();
    Double_t jy=jet.Py();
    Double_t jz=jet.Pz();
    //if(jz==0)printf("problem\n");
    //X axis 
    Double_t Xx=1.0-vertex[0];
    Double_t Xy=-1.0*vertex[1];
    Double_t Xz=jx/jz*(vertex[0]-1.)+vertex[1]*jy/jz;
    //Y axis
    Double_t Yx=jy*Xz-jz*Xy;
    Double_t Yy=jz*Xx-jx*Xz;
    Double_t Yz=jx*Xy-jy*Xx;
    //Determinant
    Double_t det = Xx*Yy*jz + Xy*Yz*jx + Xz*Yx*jy - Xx*Yz*jy - Xy*Yx*jz - Xz*Yy*jx;
    if(det==0)AliWarning("problem det==0\n");
    Double_t detX = 0.;
    Double_t detY = 0.;
    Double_t detZ = 0.;

//    Double_t tmpScalar = jx*Xx+jy*Xy+jz*Xz;
//    printf("scalar jet.P o X: %f\n",tmpScalar);
//    tmpScalar = jet.Px()*Yx+jet.Py()*Yy+jet.Pz()*Yz;
//    printf("scalar jet.P o Y: %f\n",tmpScalar);
//    tmpScalar = Xx*Yx+Xy*Yy+Xz*Yz;
//    printf("scalar X o Y: %f\n",tmpScalar);

    TVector3 perp;
    //randomise
    do
    {
      refPhi=fGenerator->Rndm()*TMath::Pi()*2.;
      //refPhi=fGenerator->Uniform(-1,1)*TMath::Pi();
      xVar=TMath::Cos(refPhi);
      yVar=TMath::Sin(refPhi);
      //yVar=TMath::Cos(-TMath::Pi()/2.+refPhi);
      //zVar=0 in new surface frame
      detX = xVar*Yy*jz + Xz*yVar*jy - xVar*Yz*jy - Xy*yVar*jz;
      detY = Xx*yVar*jz + xVar*Yz*jx - xVar*Yx*jz - Xz*yVar*jx;
      detZ = Xy*yVar*jx + xVar*Yx*jy - Xx*yVar*jy - xVar*Yy*jx;

      newX=detX/det;
      newY=detY/det;
      newZ=detZ/det;

      perp.SetXYZ(newX,newY,newZ);
      refEta = perp.Eta();
      refPhi = perp.Phi();//output in <-pi, pi> range
      if(refPhi<0)refPhi+=2*TMath::Pi();
      rad = TMath::Sqrt((gammaEta-refEta)*(gammaEta-refEta) + (gammaPhi-refPhi)*(gammaPhi-refPhi));
      rad2 = TMath::Sqrt((jetEta-refEta)*(jetEta-refEta) + (jetPhi-refPhi)*(jetPhi-refPhi));
      //printf("refEta,refPhi,rad,rad2: %f %f %f %f\n",refEta,refPhi,rad,rad2);
    } while (rad<fJetConeSize+fGammaConeSize || rad2<2.*fJetConeSize || TMath::Abs(refEta)>0.9-fJetConeSize);
      fhRandomPhiEta[2]->Fill(refPhi, refEta, GetEventWeight());

  }
  else if(type==2)
  {//random cone
    //randomise
    do
    {
      refPhi=fGenerator->Rndm()*TMath::Pi()*2.;
      refEta=fGenerator->Uniform(-(0.9-fJetConeSize),0.9-fJetConeSize);
      rad = TMath::Sqrt((gammaEta-refEta)*(gammaEta-refEta) + (gammaPhi-refPhi)*(gammaPhi-refPhi));
      rad2 = TMath::Sqrt((jetEta-refEta)*(jetEta-refEta) + (jetPhi-refPhi)*(jetPhi-refPhi));
      //check if reference is not within jet cone or gamma cone +0.4
      //example: FF fConeSize=0.4, fJetConeSize=0.4, fIsoGammaConeSize=0.3
    } while (rad<fJetConeSize+fGammaConeSize || rad2<2.*fJetConeSize);
    //photon:0.7=0.4+0.3; jets:0.8=0.4 +0.4 //rad<fConeSize+fJetConeSize rad2<fConeSize+0.3
    fhRandomPhiEta[0]->Fill(refPhi, refEta, GetEventWeight());
  }
  else if(type==4){//perpendicular to photon axis
    Double_t xVar;
    Double_t yVar;
    Double_t newX=0.;
    Double_t newY=0.;
    Double_t newZ=0.;
    //Z axis vector: [jet.Px(),jet.Py(),jet.Pz()]
    Double_t jx=gamma.Px();
    Double_t jy=gamma.Py();
    Double_t jz=gamma.Pz();
    //if(jz==0)printf("problem\n");
    //X axis 
    Double_t Xx=1.0-vertex[0];
    Double_t Xy=-1.0*vertex[1];
    Double_t Xz=jx/jz*(vertex[0]-1.)+vertex[1]*jy/jz;
    //Y axis
    Double_t Yx=jy*Xz-jz*Xy;
    Double_t Yy=jz*Xx-jx*Xz;
    Double_t Yz=jx*Xy-jy*Xx;
    //Determinant
    Double_t det = Xx*Yy*jz + Xy*Yz*jx + Xz*Yx*jy - Xx*Yz*jy - Xy*Yx*jz - Xz*Yy*jx;
    if(det==0)AliWarning("problem det==0");
    Double_t detX = 0.;
    Double_t detY = 0.;
    Double_t detZ = 0.;

//    Double_t tmpScalar = jx*Xx+jy*Xy+jz*Xz;
//    printf("scalar jet.P o X: %f\n",tmpScalar);
//    tmpScalar = jet.Px()*Yx+jet.Py()*Yy+jet.Pz()*Yz;
//    printf("scalar jet.P o Y: %f\n",tmpScalar);
//    tmpScalar = Xx*Yx+Xy*Yy+Xz*Yz;
//    printf("scalar X o Y: %f\n",tmpScalar);

    TVector3 perp;
    //randomise
    do{
      refPhi=fGenerator->Rndm()*TMath::Pi()*2.;
      //refPhi=fGenerator->Uniform(-1,1)*TMath::Pi();
      xVar=TMath::Cos(refPhi);
      yVar=TMath::Sin(refPhi);
      //yVar=TMath::Cos(-TMath::Pi()/2.+refPhi);
      //zVar=0 in new surface frame
      detX = xVar*Yy*jz + Xz*yVar*jy - xVar*Yz*jy - Xy*yVar*jz;
      detY = Xx*yVar*jz + xVar*Yz*jx - xVar*Yx*jz - Xz*yVar*jx;
      detZ = Xy*yVar*jx + xVar*Yx*jy - Xx*yVar*jy - xVar*Yy*jx;

      newX=detX/det;
      newY=detY/det;
      newZ=detZ/det;

      perp.SetXYZ(newX,newY,newZ);
      refEta = perp.Eta();
      refPhi = perp.Phi();//output in <-pi, pi> range
      if(refPhi<0)refPhi+=2*TMath::Pi();
      rad = TMath::Sqrt((gammaEta-refEta)*(gammaEta-refEta) + (gammaPhi-refPhi)*(gammaPhi-refPhi));
      rad2 = TMath::Sqrt((jetEta-refEta)*(jetEta-refEta) + (jetPhi-refPhi)*(jetPhi-refPhi));
      //printf("refEta,refPhi,rad,rad2: %f %f %f %f\n",refEta,refPhi,rad,rad2);
    } while (rad<fJetConeSize+fGammaConeSize || rad2<2.*fJetConeSize || TMath::Abs(refEta)>0.9-fJetConeSize);
    fhRandomPhiEta[1]->Fill(refPhi, refEta, GetEventWeight());

  }
  else if(type==3){//mid point

    Double_t jx=jet.Px();
    Double_t jy=jet.Py();
    Double_t jz=jet.Pz();
    //    if(jz==0)printf("problem\n");
    Double_t gx=gamma.Px();
    Double_t gy=gamma.Py();
    Double_t gz=gamma.Pz();

    Double_t cosAlpha=(jx*gx+jy*gy+jz*gz)/(jet.Mag()*gamma.Mag());
    Double_t cosinus=TMath::Sqrt((cosAlpha+1.)/2.);
    //perpendicular axis
    Double_t Zx=gy*jz-gz*jy;
    Double_t Zy=gz*jx-gx*jz;
    Double_t Zz=gx*jy-gy*jx;

    //Determinant
    Double_t det = Zx*gy*jz + Zy*gz*jx + Zz*gx*jy - Zz*gy*jx - Zy*gx*jz - Zx*gz*jy;

    Double_t newX=0.;
    Double_t newY=0.;
    Double_t newZ=0.;
    if(det!=0) {
      Double_t detX =            -Zy*gz*cosinus +Zz*cosinus*jy + Zz*gy*cosinus - Zy*cosinus*jz;
      Double_t detY = Zx*cosinus*jz          - Zz*gx*cosinus - Zz*cosinus*jx            + Zx*gz*cosinus;
      Double_t detZ = -Zx*gy*cosinus + Zy*cosinus*jx + Zy*gx*cosinus - Zx*cosinus*jy;

      newX=detX/det;
      newY=detY/det;
      newZ=detZ/det;
    }

    TVector3 perp;
    perp.SetXYZ(newX,newY,newZ);
    refEta = perp.Eta();
    refPhi = perp.Phi();//output in <-pi, pi> range
    if(refPhi<0)refPhi+=2*TMath::Pi();
    rad = TMath::Sqrt((gammaEta-refEta)*(gammaEta-refEta) + (gammaPhi-refPhi)*(gammaPhi-refPhi));
    rad2 = TMath::Sqrt((jetEta-refEta)*(jetEta-refEta) + (jetPhi-refPhi)*(jetPhi-refPhi));
      //printf("refEta,refPhi,rad,rad2: %f %f %f %f\n",refEta,refPhi,rad,rad2);

    if (rad<fJetConeSize+fGammaConeSize || rad2<2.*fJetConeSize || TMath::Abs(refEta)>0.9-fJetConeSize)
        fhRandomPhiEta[3]->Fill(refPhi, refEta, GetEventWeight());
  }
  else if(type==5){//tmp                                                                                                                                                   
    //printf("vertex: %f %f %f \n",vertex[0],vertex[1],vertex[2]);                                                                                                         

    Double_t xVar;
    Double_t newX=0.;
    Double_t newY=0.;
    Double_t newZ=0.;
    //Z axis vector: [jet.Px(),jet.Py(),jet.Pz()]                                                                                                                          
    Double_t jx=jet.Px();
    Double_t jy=jet.Py();
    Double_t jz=jet.Pz();
    //    if(jz==0)printf("problem\n");
    //X axis                                                                                                                                                               
    Double_t Xx=1.0-vertex[0];
    Double_t Xy=-1.0*vertex[1];
    Double_t Xz=jx/jz*(vertex[0]-1.)+vertex[1]*jy/jz;
    //Y axis                                                                                                                                                               
    Double_t Yx=jy*Xz-jz*Xy;
    Double_t Yy=jz*Xx-jx*Xz;
    Double_t Yz=jx*Xy-jy*Xx;

    // X and Y length                                                                                                                                                      
    Double_t Xlength=TMath::Sqrt(Xx*Xx+Xy*Xy+Xz*Xz);
    Double_t Ylength=TMath::Sqrt(Yx*Yx+Yy*Yy+Yz*Yz);
    Double_t ratio=Ylength/Xlength;

    TVector3 perp;
    //randomise                                                                                                                                                            
    do{
      refPhi=fGenerator->Rndm()*TMath::Pi()*2.;
      xVar=TMath::Tan(refPhi)/ratio;
      newX=xVar*Yx+Xx;
      newY=xVar*Yy+Xy;
      newZ=xVar*Yz+Xz;

      perp.SetXYZ(newX,newY,newZ);
      refEta = perp.Eta();
      refPhi = perp.Phi();//output in <-pi, pi> range                                                                                                                      
      if(refPhi<0)refPhi+=2*TMath::Pi();
      rad = TMath::Sqrt((gammaEta-refEta)*(gammaEta-refEta) + (gammaPhi-refPhi)*(gammaPhi-refPhi));
      rad2 = TMath::Sqrt((jetEta-refEta)*(jetEta-refEta) + (jetPhi-refPhi)*(jetPhi-refPhi));
      //printf("refEta,refPhi,rad,rad2: %f %f %f %f\n",refEta,refPhi,rad,rad2);
    } while (rad<fJetConeSize+fGammaConeSize || rad2<2.*fJetConeSize || TMath::Abs(refEta)>0.9-fJetConeSize);
    fhRandomPhiEta[4]->Fill(refPhi, refEta, GetEventWeight());
  }

  // calculate FF in background
  Int_t ntracks =  0;
  ntracks =GetCTSTracks()->GetEntriesFast();
  AliVTrack* track = 0x0 ;
  TVector3 p3;

  Double_t pt = 0, eta = 0, phi = 0;
  Int_t	 npartcone = 0;
  Double_t sumPt=0.;
  for(Int_t ipr = 0;ipr < ntracks ; ipr ++ ){
    track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    p3.SetXYZ(track->Px(),track->Py(),track->Pz());
    pt   = p3.Pt();
    if(pt<fPtThresholdInCone) {//0.150
      //printf("problem: track pt < %f MeV/c \n",fPtThresholdInCone);
      continue;
    }
    eta  = p3.Eta() ;
    phi  = p3.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    //Check if there is any particle inside cone with pt larger than  fPtThreshold
    rad = TMath::Sqrt((eta-refEta)*(eta-refEta) + (phi-refPhi)*(phi-refPhi));
    if(rad < fConeSize  && pt > fPtThresholdInCone){	
      //printf("charged in jet cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
        npartcone++;
	sumPt+=pt;
	if(type==1){//perp jet
	  fhBkgFFz[1] ->Fill(gammaPt, pt/gammaPt, GetEventWeight());
	  fhBkgFFxi[1]->Fill(gammaPt, TMath::Log(gammaPt/pt), GetEventWeight());
	  fhBkgFFpt[1]->Fill(gammaPt, pt, GetEventWeight());
	}
	else if(type==2){//RC
	  fhBkgFFz[0] ->Fill(gammaPt, pt/gammaPt, GetEventWeight());
	  fhBkgFFxi[0]->Fill(gammaPt, TMath::Log(gammaPt/pt), GetEventWeight());
	  fhBkgFFpt[0]->Fill(gammaPt, pt, GetEventWeight());
	}
	else if(type==3){//mid point
	  fhBkgFFz[3] ->Fill(gammaPt, pt/gammaPt, GetEventWeight());
	  fhBkgFFxi[3]->Fill(gammaPt, TMath::Log(gammaPt/pt), GetEventWeight());
	  fhBkgFFpt[3]->Fill(gammaPt, pt, GetEventWeight());
	}
	else if(type==4){//perp jet
	  fhBkgFFz[2] ->Fill(gammaPt, pt/gammaPt, GetEventWeight());
          fhBkgFFxi[2]->Fill(gammaPt, TMath::Log(gammaPt/pt), GetEventWeight());
          fhBkgFFpt[2]->Fill(gammaPt, pt, GetEventWeight());
        }
	else if(type==5){//test
	  fhBkgFFz[4] ->Fill(gammaPt, pt/gammaPt, GetEventWeight());
	  fhBkgFFxi[4]->Fill(gammaPt, TMath::Log(gammaPt/pt), GetEventWeight());
	  fhBkgFFpt[4]->Fill(gammaPt, pt, GetEventWeight());
	}
    }
  }//end of loop over tracks
    
  Double_t sumOverTracks=0.;
  if(npartcone!=0) sumOverTracks = sumPt/npartcone;
  if(type==1)
  {
    fhBkgNTracksInCone     [1]->Fill(gammaPt, npartcone    , GetEventWeight());
    fhBkgSumPtInCone       [1]->Fill(gammaPt, sumPt        , GetEventWeight());
    fhBkgSumPtnTracksInCone[1]->Fill(gammaPt, sumOverTracks, GetEventWeight());
  }
  else if(type==2)
  {
    fhBkgNTracksInCone     [0]->Fill(gammaPt, npartcone    , GetEventWeight());
    fhBkgSumPtInCone       [0]->Fill(gammaPt, sumPt        , GetEventWeight());
    fhBkgSumPtnTracksInCone[0]->Fill(gammaPt, sumOverTracks, GetEventWeight());
  }
  else if(type==3)
  {
    fhBkgNTracksInCone     [3]->Fill(gammaPt, npartcone    , GetEventWeight());
    fhBkgSumPtInCone       [3]->Fill(gammaPt, sumPt        , GetEventWeight());
    fhBkgSumPtnTracksInCone[3]->Fill(gammaPt, sumOverTracks, GetEventWeight());
  }
  else if(type==4)
  {
    fhBkgNTracksInCone     [2]->Fill(gammaPt, npartcone   , GetEventWeight());
    fhBkgSumPtInCone       [2]->Fill(gammaPt,sumPt        , GetEventWeight());
    fhBkgSumPtnTracksInCone[2]->Fill(gammaPt,sumOverTracks, GetEventWeight());
  }
  else if(type==5)
  {
    fhBkgNTracksInCone     [4]->Fill(gammaPt, npartcone    , GetEventWeight());
    fhBkgSumPtInCone       [4]->Fill(gammaPt, sumPt        , GetEventWeight());
    fhBkgSumPtnTracksInCone[4]->Fill(gammaPt, sumOverTracks, GetEventWeight());
  }
}

//__________________________________________________________________
/// Find information about photon and (quark or gluon) on generated level.
//__________________________________________________________________
void AliAnaParticleJetFinderCorrelation::FindMCgenInfo()
{
  if ( !GetMC() ) return ;
  
  // frequently used variables
  Int_t pdg    = 0 ;
  Int_t mother = -1 ; 
  Int_t absID  = 0 ;
  
  //Double_t photonY   = -100 ;
  //Double_t photonE   = -1 ;
  Double_t photonPt  = -1 ;
  Double_t photonPhi =  100 ;
  Double_t photonEta = -1 ;
  Bool_t   inacceptance = kFALSE;
  AliVParticle * primTmp = NULL;
  
  // jet counters
  Int_t nParticlesInJet=0;
  Int_t nChargedParticlesInJet=0;
  Int_t nParticlesInJet150=0;
  Int_t nChargedParticlesInJet150=0;
  Int_t nChargedParticlesInJet150Cone=0;
  
  Double_t eneParticlesInJet=0.;
  Double_t eneChargedParticlesInJet=0.;
  Double_t eneParticlesInJet150=0.;
  Double_t eneChargedParticlesInJet150=0.;
  Double_t eneChargedParticlesInJet150Cone=0.;
  
  Double_t pxParticlesInJet=0.;
  Double_t pxChargedParticlesInJet=0.;
  Double_t pxParticlesInJet150=0.;
  Double_t pxChargedParticlesInJet150=0.;
  Double_t pxChargedParticlesInJet150Cone=0.;
  
  Double_t pyParticlesInJet=0.;
  Double_t pyChargedParticlesInJet=0.;
  Double_t pyParticlesInJet150=0.;
  Double_t pyChargedParticlesInJet150=0.;
  Double_t pyChargedParticlesInJet150Cone=0.;
  
  Double_t etaParticlesInJet=0.;
  Double_t etaChargedParticlesInJet=0.;
  Double_t etaParticlesInJet150=0.;
  Double_t etaChargedParticlesInJet150=0.;
  Double_t etaChargedParticlesInJet150Cone=0.;
  
  Double_t phiParticlesInJet=0.;
  Double_t phiChargedParticlesInJet=0.;
  Double_t phiParticlesInJet150=0.;
  Double_t phiChargedParticlesInJet150=0.;
  Double_t phiChargedParticlesInJet150Cone=0.;
  
  Double_t ptParticlesInJet=0.;
  Double_t ptChargedParticlesInJet=0.;
  Double_t ptParticlesInJet150=0.;
  Double_t ptChargedParticlesInJet150=0.;
  Double_t ptChargedParticlesInJet150Cone=0.;
  
  Double_t coneJet=0.;
  Double_t coneChargedJet=0.;
  Double_t coneJet150=0.;
  Double_t coneChargedJet150=0.;
  
  std::vector<Int_t> jetParticleIndex;
  
  //jet origin
  //index =6 and 7 is hard scattering (jet-quark or photon)
  primTmp =  GetMC()->GetTrack(6);
  pdg=primTmp->PdgCode();
  AliDebug(3,Form("id 6 pdg %d, pt %f ",pdg,primTmp->Pt() ));
  
  if(TMath::Abs(pdg)<=6 ||pdg==21) 
  {
    fhMCJetOrigin->Fill(pdg, GetEventWeight());
    fMCPartonType=pdg;
  }
  
  primTmp =  GetMC()->GetTrack(7);
  pdg=primTmp->PdgCode();
  
  AliDebug(3,Form("id 7 pdg %d, pt %f",pdg,primTmp->Pt() ));
  
  if(TMath::Abs(pdg)<=6 ||pdg==21) 
  {
    fhMCJetOrigin->Fill(pdg, GetEventWeight());
    fMCPartonType=pdg;
  }
  //end of jet origin
  
  Int_t nprim = GetMC()->GetNumberOfTracks();
  for(Int_t i=0; i < nprim; i++) 
  {
    if ( !GetReader()->AcceptParticleMCLabel( i ) ) continue ;
    
    AliVParticle * prim =  GetMC()->GetTrack(i);
    
    pdg = prim->PdgCode();
    mother=prim->GetMother();
    
    //photon=22, gluon=21, quarks=(1,...,6), antiquarks=(-1,...,-6)
    if(pdg == 22)
    {//photon
      fhMCPhotonCuts->Fill(0., GetEventWeight());
      
      if(prim->MCStatusCode()!=1) continue;
      
      fhMCPhotonCuts->Fill(1., GetEventWeight());
      
      AliDebug(5,Form("id %d, prim %d, physPrim %d, status %d\n",i,prim->IsPrimary(),prim->IsPhysicalPrimary(),prim->MCStatusCode()));
      while(mother>7)
      {
        primTmp =  GetMC()->GetTrack(mother);
        mother=primTmp->GetMother();
      }
      
      if(mother<6)continue;
      
      fhMCPhotonCuts->Fill(2., GetEventWeight());
      
      primTmp =  GetMC()->GetTrack(mother);
      
      if(primTmp->PdgCode()!=22)continue;
      
      fhMCPhotonCuts->Fill(3., GetEventWeight());
      
      //Get photon kinematics
      photonPt  = prim->Pt() ;
      photonPhi = prim->Phi() ;
      if(photonPhi < 0) photonPhi+=TMath::TwoPi();
      photonEta = prim->Eta() ;
      
      fhMCPhotonPt    ->Fill(photonPt, GetEventWeight());
      fhMCPhotonEtaPhi->Fill(photonPhi, photonEta, GetEventWeight());
      
      //Check if photons hit the Calorimeter
      fMomentum.SetPxPyPzE(prim->Px(),prim->Py(),prim->Pz(),prim->E());
      inacceptance = kFALSE;
      if(GetCaloUtils()->IsEMCALGeoMatrixSet())
      {
        fhMCPhotonCuts->Fill(4, GetEventWeight());
        
        //check if in EMCAL
        if(GetEMCALGeometry())
        {
          GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(prim->Eta(),prim->Phi(),absID);
          if(absID >= 0) inacceptance = kTRUE;
          AliDebug(3,Form("In EMCAL Real acceptance? %d",inacceptance));
        }
        else{
          if(GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kEMCAL)) inacceptance = kTRUE ;
          AliDebug(1,Form("In EMCAL fiducial cut acceptance? %d",inacceptance));
        }
      }
      else
      {//no EMCAL nor EMCALGeoMatrixSet
        AliWarning("not EMCALGeoMatrix set");
      }//end of check if EMCAL
      
      if(inacceptance)fhMCPhotonCuts->Fill(5, GetEventWeight());
      
      AliDebug(5,Form("Photon Energy %f, Pt %f",prim->E(),prim->Pt()));
      fMCGamPt=photonPt;
      fMCGamEta=photonEta;
      fMCGamPhi=photonPhi;
    }//end of check if photon
    else
    {//not photon
      if(prim->MCStatusCode()!=1) continue;
      
      AliDebug(5,Form("id %d, prim %d, physPrim %d, status %d, pdg %d, E %f",
                      i,prim->IsPrimary(),prim->IsPhysicalPrimary(),prim->MCStatusCode(),prim->PdgCode(),prim->E()));
      
      while(mother>7)
      {
        primTmp =  GetMC()->GetTrack(mother);
        mother=primTmp->GetMother();
        AliDebug(5,Form("next mother %d",mother));
      }
      
      if(mother<6)continue;//soft part
      
      primTmp =  GetMC()->GetTrack(mother);
      pdg=primTmp->PdgCode();
      if( !(TMath::Abs(pdg)<=6 || pdg==21) ) continue;//origin not hard q or g
      
      //jetParticleIndex.Add(&i);
      jetParticleIndex.push_back(i);
      
      nParticlesInJet++;
      eneParticlesInJet+=prim->E();
      pxParticlesInJet+=prim->Px();
      pyParticlesInJet+=prim->Py();
      etaParticlesInJet+=(prim->E()*prim->Eta());
      photonPhi = prim->Phi() ;
      if(photonPhi < 0) photonPhi+=TMath::TwoPi();
      phiParticlesInJet+=(prim->E()*photonPhi);
      
      
      if(prim->Charge()!=0)
      {
        nChargedParticlesInJet++;
        eneChargedParticlesInJet+=prim->E();
        pxChargedParticlesInJet+=prim->Px();
        pyChargedParticlesInJet+=prim->Py();
        etaChargedParticlesInJet+=(prim->E()*prim->Eta());
        phiChargedParticlesInJet+=(prim->E()*photonPhi);
      }
      
      if(prim->Pt()>0.150 && TMath::Abs(prim->Eta())<0.9 ) 
      {//in TPC acceptance :)
        nParticlesInJet150++;
        eneParticlesInJet150+=prim->E();
        pxParticlesInJet150+=prim->Px();
        pyParticlesInJet150+=prim->Py();
        etaParticlesInJet150+=(prim->E()*prim->Eta());
        phiParticlesInJet150+=(prim->E()*photonPhi);
      }
      
      if(prim->Charge()!=0 && prim->Pt()>0.150 && TMath::Abs(prim->Eta())<0.9 )
      {//in TPC acceptance :)
        nChargedParticlesInJet150++;
        eneChargedParticlesInJet150+=prim->E();
        pxChargedParticlesInJet150+=prim->Px();
        pyChargedParticlesInJet150+=prim->Py();
        etaChargedParticlesInJet150+=(prim->E()*prim->Eta());
        phiChargedParticlesInJet150+=(prim->E()*photonPhi);
      }
      
    }//end of check pdg
  }//end of loop over primaries
  
  if(eneParticlesInJet != 0.) {
    etaParticlesInJet/=eneParticlesInJet ;
    phiParticlesInJet/=eneParticlesInJet ;
  }
  if(eneChargedParticlesInJet != 0) {
    etaChargedParticlesInJet/=eneChargedParticlesInJet;
    phiChargedParticlesInJet/=eneChargedParticlesInJet;
  }
  if(eneParticlesInJet150 != 0) {
    etaParticlesInJet150/=eneParticlesInJet150;
    phiParticlesInJet150/=eneParticlesInJet150;
  }
  if(eneChargedParticlesInJet150 != 0) {
    etaChargedParticlesInJet150/=eneChargedParticlesInJet150;
    phiChargedParticlesInJet150/=eneChargedParticlesInJet150;
  }
  
  ptParticlesInJet=TMath::Sqrt(pxParticlesInJet*pxParticlesInJet+pyParticlesInJet*pyParticlesInJet);
  ptChargedParticlesInJet=TMath::Sqrt(pxChargedParticlesInJet*pxChargedParticlesInJet+pyChargedParticlesInJet*pyChargedParticlesInJet);
  ptParticlesInJet150=TMath::Sqrt(pxParticlesInJet150*pxParticlesInJet150+pyParticlesInJet150*pyParticlesInJet150);
  ptChargedParticlesInJet150=TMath::Sqrt(pxChargedParticlesInJet150*pxChargedParticlesInJet150+pyChargedParticlesInJet150*pyChargedParticlesInJet150);
  
  Double_t distance=0.;
  Double_t eta=0.;
  Double_t phi=0.;
  Double_t mostPtCharged=0.;
  Int_t mostmostPtChargedId=-1;
  std::vector<Int_t>::iterator it;
  for( it=jetParticleIndex.begin(); it!=jetParticleIndex.end(); ++it )
  {
    AliVParticle * prim = GetMC()->GetTrack(*it);
    eta = prim->Eta();
    phi = prim->Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    //full jet
    distance=TMath::Sqrt((eta-etaParticlesInJet)*(eta-etaParticlesInJet)+(phi-phiParticlesInJet)*(phi-phiParticlesInJet));
    if(distance>coneJet) coneJet=distance;
    //charged jet
    distance=TMath::Sqrt((eta-etaChargedParticlesInJet)*(eta-etaChargedParticlesInJet)+(phi-phiChargedParticlesInJet)*(phi-phiChargedParticlesInJet));
    if(distance>coneChargedJet) coneChargedJet=distance;
    //
    distance=TMath::Sqrt((eta-etaParticlesInJet150)*(eta-etaParticlesInJet150)+(phi-phiParticlesInJet150)*(phi-phiParticlesInJet150));
    if(distance>coneJet150 && TMath::Abs(eta)<0.9 ) coneJet150=distance;
    //
    distance=TMath::Sqrt((eta-etaChargedParticlesInJet150)*(eta-etaChargedParticlesInJet150)+(phi-phiChargedParticlesInJet150)*(phi-phiChargedParticlesInJet150));
    if(distance>coneChargedJet150 && TMath::Abs(eta)<0.9) coneChargedJet150=distance;
    
    if(prim->Charge()!=0 && prim->Pt()>0.150 && TMath::Abs(eta)<0.9) {
      if(prim->Pt()>mostPtCharged) {
        mostPtCharged=prim->Pt();
        mostmostPtChargedId=(*it);
      }
    }
    
    if(distance<=0.4){
      if(prim->Charge()!=0 && prim->Pt()>0.150 && TMath::Abs(eta)<0.9) {
        nChargedParticlesInJet150Cone++;
        eneChargedParticlesInJet150Cone+=prim->E();
        pxChargedParticlesInJet150Cone+=prim->Px();
        pyChargedParticlesInJet150Cone+=prim->Py();
        etaChargedParticlesInJet150Cone+=(prim->E()*eta);
        phiChargedParticlesInJet150Cone+=(prim->E()*phi);
      }
      
    }
    
  }//end of loop over jet particle indexes
  if(eneChargedParticlesInJet150Cone != 0) {
    etaChargedParticlesInJet150Cone/=eneChargedParticlesInJet150Cone;
    phiChargedParticlesInJet150Cone/=eneChargedParticlesInJet150Cone;
  }
  
  ptChargedParticlesInJet150Cone=TMath::Sqrt(pxChargedParticlesInJet150Cone*pxChargedParticlesInJet150Cone+pyChargedParticlesInJet150Cone*pyChargedParticlesInJet150Cone);
  
  if(nChargedParticlesInJet150>0 && nChargedParticlesInJet150Cone<1)
  {//no particles in cone? take the most energetic one
    nChargedParticlesInJet150Cone=1;
    etaChargedParticlesInJet150Cone=(GetMC()->GetTrack(mostmostPtChargedId))->Eta();
    phiChargedParticlesInJet150Cone=(GetMC()->GetTrack(mostmostPtChargedId))->Phi();
    ptChargedParticlesInJet150Cone =(GetMC()->GetTrack(mostmostPtChargedId))->Pt();
  }
  
  jetParticleIndex.clear();
  
  //printouts
  
  AliDebug(3,Form("cone full %f, charged %f, full150 %f, charged150 %f",coneJet,coneChargedJet,coneJet150,coneChargedJet150));
  AliDebug(3,Form("Npart %d, NchPart %d, Npart(pt>150M) %d, NchPart(pt>150M) %d, NchPart(pt>150M)Cone %d\n",nParticlesInJet,nChargedParticlesInJet,nParticlesInJet150,nChargedParticlesInJet150,nChargedParticlesInJet150Cone));
  AliDebug(3,Form("Etot %f, Ech %f, E(pt>150M) %f, Ech(pt>150M) %f\n",eneParticlesInJet,eneChargedParticlesInJet,eneParticlesInJet150,eneChargedParticlesInJet150));
  AliDebug(3,Form("pt %f, ptch %f, pt(pt>150M) %f,ptch(pt>150M) %f,ptch(pt>150M)Cone %f\n",ptParticlesInJet,ptChargedParticlesInJet,ptParticlesInJet150,ptChargedParticlesInJet150,ptChargedParticlesInJet150Cone));
  AliDebug(3,Form("eta/phi tot %f/%f, ch %f/%f, tot150 %f/%f,  ch150 %f/%f, ch150cone %f/%f\n",etaParticlesInJet,phiParticlesInJet,etaChargedParticlesInJet,phiChargedParticlesInJet,etaParticlesInJet150,phiParticlesInJet150,etaChargedParticlesInJet150,phiChargedParticlesInJet150,etaChargedParticlesInJet150Cone,phiChargedParticlesInJet150Cone));
  
  //fill histograms
  if(ptParticlesInJet) fhMCJetRatioChFull->Fill(ptChargedParticlesInJet/ptParticlesInJet, GetEventWeight());
  if(ptChargedParticlesInJet) fhMCJetRatioCh150Ch->Fill(ptChargedParticlesInJet150/ptChargedParticlesInJet, GetEventWeight());
  
  fhMCJetNPartVsPt     ->Fill(ptParticlesInJet,nParticlesInJet, GetEventWeight());
  fhMCJetChNPartVsPt   ->Fill(ptChargedParticlesInJet,nChargedParticlesInJet, GetEventWeight());
  fhMCJetNPart150VsPt  ->Fill(ptParticlesInJet150,nParticlesInJet150, GetEventWeight());
  fhMCJetChNPart150VsPt->Fill(ptChargedParticlesInJet150,nChargedParticlesInJet150, GetEventWeight());
  fhMCJetChNPart150ConeVsPt->Fill(ptChargedParticlesInJet150Cone,nChargedParticlesInJet150Cone, GetEventWeight());
  
  fhMCJetEtaPhi->Fill(phiParticlesInJet,etaParticlesInJet, GetEventWeight());
  fhMCJetChEtaPhi->Fill(phiChargedParticlesInJet,etaChargedParticlesInJet, GetEventWeight());
  fhMCJet150EtaPhi->Fill(phiParticlesInJet150,etaParticlesInJet150, GetEventWeight());
  fhMCJetCh150EtaPhi->Fill(phiChargedParticlesInJet150,etaChargedParticlesInJet150, GetEventWeight());
  fhMCJetCh150ConeEtaPhi->Fill(phiChargedParticlesInJet150Cone,etaChargedParticlesInJet150Cone, GetEventWeight());
  
  //fill tree
  fMCJetPt      = ptParticlesInJet;
  fMCJetChPt    = ptChargedParticlesInJet;      
  fMCJet150Pt   = ptParticlesInJet150;     
  fMCJetCh150Pt = ptChargedParticlesInJet150;   
  fMCJetNPart      = nParticlesInJet;     
  fMCJetChNPart    = nChargedParticlesInJet;   
  fMCJet150NPart   = nParticlesInJet150;  
  fMCJetCh150NPart = nChargedParticlesInJet150;
  fMCJetEta      = etaParticlesInJet          ;
  fMCJetPhi      = phiParticlesInJet	      ;
  fMCJetChEta    = etaChargedParticlesInJet   ;
  fMCJetChPhi    = phiChargedParticlesInJet   ;
  fMCJet150Eta   = etaParticlesInJet150	      ;
  fMCJet150Phi   = phiParticlesInJet150	      ;
  fMCJetCh150Eta = etaChargedParticlesInJet150;
  fMCJetCh150Phi = phiChargedParticlesInJet150;
  
  fMCJetCh150ConePt    = ptChargedParticlesInJet150Cone;
  fMCJetCh150ConeNPart = nChargedParticlesInJet150Cone;
  fMCJetCh150ConeEta   = etaChargedParticlesInJet150Cone;
  fMCJetCh150ConePhi   = phiChargedParticlesInJet150Cone;
  
} // end of method FindMCgenInfo


