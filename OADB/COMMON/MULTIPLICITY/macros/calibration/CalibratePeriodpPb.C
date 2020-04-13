#include "AliMultEstimator.h"
#include "AliMultSelectionCuts.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCalibrator.h"
#include <TString.h>
#include <TSystem.h>
#include <TF1.h>
#include <TFile.h>

void CalibratePeriodpPb(  TString lPeriodName         = "LHC16s",
                          TString inputFileName       = "input.root",
                          TString outputOADBFileName  = "OADB.root",
                          TString bufferFileName      = "buffer.root",
                          TString fileNameZVtxCorr    = "ZVtzCorrection.root"
) {

  //Load ALICE stuff
  TString gLibs[] = {"STEER", "ANALYSIS", "ANALYSISalice", "ANALYSIScalib","OADB"
                          };
  TString thislib = "lib";
  for(Int_t ilib = 0; ilib<5; ilib++) {
      thislib="lib";
      thislib.Append(gLibs[ilib].Data());
      cout<<"Will load "<<thislib.Data()<<endl;
      gSystem->Load(thislib.Data());
  }
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  cout<<"Alive! "<<endl;

  //All fine, let's try the calibrator
  AliMultSelectionCalibrator *lCalib = new AliMultSelectionCalibrator("lCalib");

  //============================================================
  // --- Definition of Boundaries ---
  //============================================================

  //Set Adaptive Percentile Boundaries, adjust if finer selection desired
  Double_t lDesiredBoundaries[1000];
  Long_t   lNDesiredBoundaries=0;
  lDesiredBoundaries[0] = 100;

  //Very simple 1%-wide bins all the way
  for( Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 1.0;
  }
  for( Int_t ib = 1; ib < 101; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.1;
  }

  lNDesiredBoundaries++;
  lDesiredBoundaries[lNDesiredBoundaries] = 0;

  //cout<< "Dump for debug: "<<endl;
  //for( Int_t ib=0;ib<101; ib++)
  //    cout<<"Boundary #"<<ib<<" at "<<lDesiredBoundaries[ib]<<endl;

  lCalib->SetBoundaries( lNDesiredBoundaries, lDesiredBoundaries );
  cout<<"Boundaries set. Will attempt calibration now... "<<endl;

  //============================================================
  // --- Set event selection criteria ---
  //============================================================

  if ( lPeriodName.Contains("LHC16q") || lPeriodName.Contains("LHC16t") || lPeriodName.Contains("LHC13b") || lPeriodName.Contains("LHC13c")) {
    cout<<"Setting event selection criteria for p-Pb..."<<endl;
    lCalib->GetEventCuts()->SetVzCut(10.0);
    lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
    lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
    lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
    lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
    lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
    lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
  }

  if ( lPeriodName.Contains("LHC16r") || lPeriodName.Contains("LHC16s")  || lPeriodName.Contains("LHC13d") || lPeriodName.Contains("LHC13e") || lPeriodName.Contains("LHC13f") ) {
    cout<<"Setting event selection criteria for HI p-Pb and Pb-p ..."<<endl;
    lCalib->GetEventCuts()->SetVzCut(10.0);
    lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
    lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
    lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
    lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kTRUE);
    lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
    lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
  }

  //============================================================
  // --- Definition of run ranges for triggered periods ---
  //============================================================

  // switch for run-wise or runRange calibration
  Bool_t doRunWiseCalib = kTRUE;
  if (lPeriodName.Contains("LHC13d") || lPeriodName.Contains("LHC13e") || lPeriodName.Contains("LHC13f"))
    doRunWiseCalib      = kFALSE;

  // for low MB stat periods set run ranges according to fills
  Int_t nMaxRunsLHC13d    = 6;
  Int_t runRangeLHC13d[7] = { 195679, 195720, 195760, 195783, 195826, 195867, 195874};
  Int_t nMaxRunsLHC13e    = 5;
  Int_t runRangeLHC13e[6] = { 195935, 195988, 196083, 196184, 196308, 196312 };
  // some fills were merged due to too little statistics
  Int_t nMaxRunsLHC13f    = 13;
  Int_t runRangeLHC13f[15]= { 196528, 196563, 196601, 196701, 196772, 196965, 197002, 197084, 197138, 197247,
                              197296, 197341, 197386, 197389};

  //============================================================
  // --- Definition of Input Variables ---
  //============================================================

  lCalib->SetupStandardInput();

  //Changes in new version: create AliMultSelection here
  AliMultSelection *lMultSel = new AliMultSelection();

  if ( lPeriodName.Contains("LHC13b") ) lCalib->SetRunToUseAsDefault( 195483 ); // LHC13b
  if ( lPeriodName.Contains("LHC13c") ) lCalib->SetRunToUseAsDefault( 195593 ); // LHC13c
  if ( lPeriodName.Contains("LHC13d") ) lCalib->SetRunToUseAsDefault( 195831 ); // LHC13d
  if ( lPeriodName.Contains("LHC13e") ) lCalib->SetRunToUseAsDefault( 196310 ); // LHC13e
  if ( lPeriodName.Contains("LHC13f") ) lCalib->SetRunToUseAsDefault( 196563 ); // LHC13f
  if ( lPeriodName.Contains("LHC16q") ) lCalib->SetRunToUseAsDefault( 265309 ); // LHC16q
  if ( lPeriodName.Contains("LHC16r") ) lCalib->SetRunToUseAsDefault( 266318 ); // LHC16r
  if ( lPeriodName.Contains("LHC16s") ) lCalib->SetRunToUseAsDefault( 266441 ); // LHC16s
  if ( lPeriodName.Contains("LHC16t") ) lCalib->SetRunToUseAsDefault( 267163 ); // LHC16t


  TString estimatorName[20]         = { "V0M", "V0A", "V0C", "V0MEq", "V0AEq", "V0CEq", "OnlineV0M", "OnlineV0A", "OnlineV0C", "ADM",
                                        "ADA", "ADC", "CL0", "CL1", "SPDClusters", "SPDTracklets", "RefMult05", "RefMult08", "ZNA", "ZNC" };
  TString estimatorNameFit[20]      = { "V0M", "V0A", "V0C", "V0MEq", "V0AEq", "V0CEq", "OnlineV0M", "OnlineV0A", "OnlineV0C", "ADM",
                                        "ADA", "ADC", "SPDCl0", "SPDCl1", "SPDCl", "NTracklets", "RefMultEta5", "RefMultEta8", "ZNA", "ZNC" };
  TString estimatorBasicString[20]  = { "(fAmplitude_V0A)+(fAmplitude_V0C)", "(fAmplitude_V0A)", "(fAmplitude_V0C)",
                                        "(fAmplitude_V0AEq)+(fAmplitude_V0CEq)", "(fAmplitude_V0AEq)", "(fAmplitude_V0CEq)",
                                        "(fAmplitude_OnlineV0A)+(fAmplitude_OnlineV0C)", "(fAmplitude_OnlineV0A)", "(fAmplitude_OnlineV0C)",
                                        "(fMultiplicity_ADA)+(fMultiplicity_ADC)", "(fMultiplicity_ADA)", "(fMultiplicity_ADC)",
                                        "(fnSPDClusters0)", "(fnSPDClusters1)", "(fnSPDClusters)", "(fnTracklets)",
                                        "(fRefMultEta5)", "(fRefMultEta8)",
                                        "(fZnaFired)*(fZnaTower)+!(fZnaFired)*(0)", "(fZncFired)*(fZncTower)+!(fZncFired)*(0)" };
  TString estimatorFullString[20];
  Bool_t isIntegerEst[20]           = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE, kFALSE, kFALSE};
  Bool_t enableEst[20]              = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kFALSE, kTRUE, kTRUE };
  Bool_t enableZVtxCorr[20]         = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kFALSE, kTRUE, kTRUE };

  //============================================================
  // --- Load Z vtx position correction functions
  //============================================================
  TFile* fileZVtxCorr     = new TFile(fileNameZVtxCorr.Data(), "OPEN");
  TF1* fitCorrZVtx[20]              = {NULL};
  for (Int_t i = 0; i < 20; i++){
    fitCorrZVtx[i]        = (TF1*)fileZVtxCorr->Get( Form("%s/ZVtxCorr_%s", lPeriodName.Data(),estimatorNameFit[i].Data() ));
  }
  //============================================================
  // --- Definition of Estimators ---
  //============================================================
  for (Int_t i = 0; i < 20; i++){
    if (fitCorrZVtx[i] && enableZVtxCorr[i]){
      TString currentFormula = fitCorrZVtx[i]->GetExpFormula();
      for (Int_t k = 0; k < fitCorrZVtx[i]->GetNpar(); k++){
        #ifndef __CLING__
          currentFormula.ReplaceAll(Form("[p%d]",k), Form("(%.10f)",fitCorrZVtx[i]->GetParameter(k)));
        #else
          currentFormula.ReplaceAll(Form("[%d]",k), Form("(%.10f)",fitCorrZVtx[i]->GetParameter(k)));
        #endif
      }
      currentFormula.ReplaceAll("x","(fEvSel_VtxZ)");
      estimatorFullString[i] = Form("((%s)/(%s))", estimatorBasicString[i].Data(),currentFormula.Data());
    } else {
      estimatorFullString[i] = estimatorBasicString[i];
    }
    cout << "calibrating: "<< estimatorName[i].Data() << "\twith :\t" << estimatorFullString[i].Data() << endl;;
    AliMultEstimator *fEst = new AliMultEstimator(estimatorName[i].Data(), "", estimatorFullString[i].Data());
    if (isIntegerEst[i])fEst->SetIsInteger(kTRUE);
    if (enableEst[i]) lMultSel-> AddEstimator( fEst );
  }

  //============================================================
  // --- Add run ranges where appropriate ---
  //============================================================
  if (doRunWiseCalib){
    lCalib->SetMultSelection(lMultSel);
  } else {
    if (!lPeriodName.CompareTo("LHC13d")){
      for (Int_t r = 0; r<nMaxRunsLHC13d; r++ ) lCalib->AddRunRange(runRangeLHC13d[r], runRangeLHC13d[r+1]-1, lMultSel);
    }
    if (!lPeriodName.CompareTo("LHC13e")){
      for (Int_t r = 0; r<nMaxRunsLHC13e; r++ ) lCalib->AddRunRange(runRangeLHC13e[r], runRangeLHC13e[r+1]-1, lMultSel);
    }
    if (!lPeriodName.CompareTo("LHC13f")){
      for (Int_t r = 0; r<nMaxRunsLHC13f; r++ ) lCalib->AddRunRange(runRangeLHC13f[r], runRangeLHC13f[r+1]-1, lMultSel);
    }
  }

  //============================================================
  // --- Set file names and start calibration ---
  //============================================================
  lCalib -> SetInputFile  ( inputFileName.Data() );
  lCalib -> SetBufferFile ( bufferFileName.Data() );
  lCalib -> SetOutputFile ( outputOADBFileName.Data() );
  lCalib -> Calibrate     ();
}
