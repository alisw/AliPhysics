#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <AliHFOfflineCorrelator.h>

//
// On AliRoot shell, call the following before loading the macro:
//
// gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ANALYSIS/macros -I$ROOTSYS/include");
// gROOT->LoadMacro("AliHFOfflineCorrelator.cxx++g"); // set the right path! - if not already in AliPhysics!!
// .L CorrelateOffline.C++ //(needs to be campiled!!)
//
// Set:
// - the arguments of the macro
// - the names of the inputs in method SetInputNames
// - the pT ranges for D, ptAssoc and (if requiring 2D plots) the mass ranges in method SetPtRanges
// - the pool edges in method SetPools

void SetInputNames(AliHFOfflineCorrelator *correlator);
void SetPtRanges(AliHFOfflineCorrelator *correlator, Bool_t make2Dplots);
void SetPools(AliHFOfflineCorrelator *correlator);

void CorrelateOffline(
   AliHFOfflineCorrelator::DMesonSpecies specie=AliHFOfflineCorrelator::kD0toKpi, //the D-meson decay channel (check the enumerator for the options)
   AliHFOfflineCorrelator::AnalysisType analysisType=AliHFOfflineCorrelator::kME, //SE or ME correlation distributions are produced
   Int_t useEff=kTRUE, //use efficiency of trigger and tracks to weight correlation entries
   Int_t make2Dplots=kTRUE, //produce also 2D (dPhi,dEta) distributions in signal and sideband regions - ***REQUIRES DEFINITION OF SIGNAL AND SB REGIONS***
   Int_t maxTracks=-1, //in ME, uses a maximum number of tracks in the tree (for each input file)
   Int_t minD=-1, Int_t maxD=-1, //start and end of loop on D mesons
   Int_t numSelD=0, Int_t numSelTr=0, //# of selection for D and tracks (0=default selection; 1,2,3,... = alternate selections)
   Int_t wgtPeriods=kTRUE, //weight periods in ME analysis (keep enabled)
   TString nameOutputFile="OfflineCorrelations.root",
   Int_t firstBinNum=0, //start of numbering for the pTbins in input file
   Double_t mincent=0., Double_t maxcent=0., //centrality (or multiplicity) selection ***ACTIVE ONLY IF BOTH VALS ARE =! 0*** 
   Bool_t rejectSoftPi=kTRUE) //if active, removes 'fake' soft pions in ME (for SE, softpicut flag is in analysis task). No effect on D* and D+ analyses
{

  AliHFOfflineCorrelator *correlator = new AliHFOfflineCorrelator();
  Bool_t flagSpecie = correlator->SetDmesonSpecie(specie);
  correlator->SetAnalysisType(analysisType);
  correlator->SetUseEfficiency(useEff);
  correlator->SetOutputFileName(nameOutputFile);
  correlator->SetMaxTracksToCorrelate(maxTracks);
  correlator->SetDLoopRange(minD,maxD);
  correlator->SetDebugLevel(1);
  correlator->SetWeightPeriods(wgtPeriods);
  correlator->SetFirstBinNum(firstBinNum);
  correlator->SetNumSelD(numSelD);
  correlator->SetNumSelTr(numSelTr);
  correlator->SetCentralitySelection(mincent,maxcent);
  correlator->SetRejectSoftPion(rejectSoftPi);
 
  if(!flagSpecie) return;

  SetPtRanges(correlator,make2Dplots);
  SetInputNames(correlator);
  SetPools(correlator);

  Bool_t success = correlator->Correlate();

  if(success) std::cout << "***** Everything went fine! *****" << endl;
  else std::cout << "***** Errors in the extraction of the correlations! *****" << endl;

  return;

}

//________________________________________
void SetInputNames(AliHFOfflineCorrelator *correlator) {

  correlator->AddInputFile("./AnalysisResults_pp_10b.root");
  correlator->AddInputFile("./AnalysisResults_pp_10c.root");
  correlator->AddInputFile("./AnalysisResults_pp_10d.root");
  correlator->AddInputFile("./AnalysisResults_pp_10e.root");
  correlator->SetDirName("PWG3_D2H_D0InvMassOutputTTree");
  correlator->SetTreeNames("fTreeD","fTreeTr");
  correlator->SetEffMapNames("AssociatedTrkCuts","h_Eff","heff_rebin");

  return;
}

//________________________________________
void SetPtRanges(AliHFOfflineCorrelator *correlator, Bool_t make2Dplots) {

  const Int_t nptbins = 7;
  Double_t ptBinsD[nptbins+1] = {3.,4.,5.,6.,7.,8.,12.,16.};
  correlator->SetDPtBins(nptbins,ptBinsD);  
  correlator->AddAssocPtRange(0.3,99.);
  correlator->AddAssocPtRange(0.3,1.);
  correlator->AddAssocPtRange(1.,99.);

  if(make2Dplots) { //***Warning! Ranges have to be exactly those coming from the mass fits (for signal ragion) and from SB definition - otherwise results will be biased!***
    Double_t SignLowLim[nptbins] = {1.8328,1.8328,1.8248,1.8248,1.8248,1.8088,1.8248}; //to be filled looking at results from invariant mass fits!
    Double_t SignUppLim[nptbins] = {1.8968,1.8968,1.9088,1.9088,1.9168,1.9208,1.9208}; //current ranges are from pp fits of slowTrain output
    Double_t LSBLowLim[nptbins] = {1.7928,1.7768,1.7728,1.7648,1.7488,1.7448,1.7728};
    Double_t LSBUppLim[nptbins] = {1.8288,1.8208,1.8208,1.8128,1.8088,1.8048,1.8048};
    Double_t RSBLowLim[nptbins] = {1.9008,1.9088,1.9128,1.9168,1.9288,1.9288,1.9288};
    Double_t RSBUppLim[nptbins] = {1.9408,1.9528,1.9608,1.9688,1.9848,1.9888,1.9928}; 
    correlator->MakeDetaDphiPlots(SignLowLim,SignUppLim,LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  }

  return;
}

//________________________________________
void SetPools(AliHFOfflineCorrelator *correlator) {

  const Int_t nMultPools = 3;
 // Double_t multPools[nMultPools+1] = {0,40,65,500}; //pPb
  Double_t multPools[nMultPools+1] = {0,20,35,500}; //pp
  const Int_t nzVtxPools = 3;
  Double_t zVtxPools[nzVtxPools+1] = {-10,-2.5,2.5,10};

  correlator->SetPoolBins(nMultPools,multPools,nzVtxPools,zVtxPools);  

  return;
}

