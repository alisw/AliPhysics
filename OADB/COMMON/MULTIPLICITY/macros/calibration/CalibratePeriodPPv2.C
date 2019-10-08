#include "AliMultEstimator.h"
#include "AliMultSelectionCuts.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCalibrator.h"
#include <TString.h>
#include <TSystem.h>
#include <TF1.h>
#include <TFile.h>

////////////////////////////////////////////////////////////
//
// Default macro for calibrating minimum bias pp data.
//
////////////////////////////////////////////////////////////

void CalibratePeriodPPv2( TString inputFileName       = "",
                          TString lPeriodName         = "LHC16k",
                          TString lWhichTrigger       = "MB",
                          TString lAddName            = "",
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
  if (!lWhichTrigger.CompareTo("MB")) lCalib->SetSelectedTriggerClass(AliVEvent::kINT7);

  //============================================================
  // --- Definition of default runs per period               ---
  //============================================================
  if ( lPeriodName.Contains("LHC15n") ) lCalib->SetRunToUseAsDefault( 244531 ); // LHC15n
  if ( lPeriodName.Contains("LHC16k") ) lCalib->SetRunToUseAsDefault( 257630 ); // LHC16k
  if ( lPeriodName.Contains("LHC17p") ) lCalib->SetRunToUseAsDefault( 282306 ); // LHC17p
  if ( lPeriodName.Contains("LHC17q") ) lCalib->SetRunToUseAsDefault( 282366 ); // LHC17q

  //============================================================
  // --- Definition of run ranges for triggered periods ---
  //============================================================

  // switch for run-wise or runRange calibration
  Bool_t doRunWiseCalib = kTRUE;
//   if (lPeriodName.Contains("LHC13d") || lPeriodName.Contains("LHC13e") || lPeriodName.Contains("LHC13f"))
//     doRunWiseCalib      = kFALSE;


  //============================================================
  // --- Definition of Boundaries ---
  //============================================================

  //Set Adaptive Percentile Boundaries, adjust if finer selection desired
  Double_t lDesiredBoundaries[1000];
  Long_t   lNDesiredBoundaries=0;
  lDesiredBoundaries[0] = 100;
  //From Low To High Multiplicity
  for( Int_t ib = 1; ib < 91; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 1.0;
  }
  for( Int_t ib = 1; ib < 91; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.1;
  }
  for( Int_t ib = 1; ib < 91; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.01;
  }
  for( Int_t ib = 1; ib < 101; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.001;
  }
  lNDesiredBoundaries++;
  lDesiredBoundaries[lNDesiredBoundaries] = 0;

  lCalib->SetBoundaries( lNDesiredBoundaries, lDesiredBoundaries );
  cout<<"Boundaries set. Will attempt calibration now... "<<endl;

  //============================================================
  // --- Set event selection criteria ---
  //============================================================
  lCalib->GetEventCuts()->SetVzCut(10.0);
  lCalib->GetEventCuts()->SetTriggerCut                (kTRUE);
  lCalib->GetEventCuts()->SetINELgtZEROCut             (kTRUE);
  lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kTRUE);
  lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kTRUE);
  lCalib->GetEventCuts()->SetVertexConsistencyCut      (kTRUE);
  lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);

  //Additional selections for pp: incompleteDAQ
  lCalib->GetEventCuts()->SetIsNotIncompleteDAQ        (kTRUE);


  //============================================================
  // --- Definition of Input Variables ---
  //============================================================

  lCalib->SetupStandardInput();

  //Changes in new version: create AliMultSelection here
  AliMultSelection *lMultSel = new AliMultSelection();
  lCalib->SetMultSelection(lMultSel);

  //============================================================
  // --- Definition of Estimators ---
  //============================================================


  Double_t lDefaultV0MAnchor     = 0;
  Double_t lDefaultV0MPercentile = 0;

  Double_t lDefaultCL0Anchor     = 0;
  Double_t lDefaultCL0Percentile = 0;

  Double_t lDefaultCL1Anchor     = 0;
  Double_t lDefaultCL1Percentile = 0;

//   if ( lPeriodName.Contains("LHC15m") ) {
//       lDefaultV0MAnchor     = 115.0;
//       lDefaultV0MPercentile = 87.5;
//       lDefaultCL0Anchor     = 39.5;
//       lDefaultCL0Percentile = 88.9;
//       lDefaultCL1Anchor     = 40.5;
//       lDefaultCL1Percentile = 88.1;
//   }

  TString estimatorName[21]         = { "V0M", "V0A", "V0C", "V0MEq", "V0AEq", "V0CEq", "OnlineV0M", "OnlineV0A", "OnlineV0C", "ADM",
                                        "ADA", "ADC", "CL0", "CL1", "SPDClusters", "SPDTracklets", "RefMult05", "RefMult08", "ZNApp", "ZNCpp",
                                        "ZNACpp" };
  TString estimatorNameFit[21]      = { "V0M", "V0A", "V0C", "V0MEq", "V0AEq", "V0CEq", "OnlineV0M", "OnlineV0A", "OnlineV0C", "ADM",
                                        "ADA", "ADC", "SPDCl0", "SPDCl1", "SPDCl", "NTracklets", "RefMultEta5", "RefMultEta8", "ZNA", "ZNC",
                                        "ZNAC" };
  TString estimatorBasicString[21]  = { "(fAmplitude_V0A)+(fAmplitude_V0C)", "(fAmplitude_V0A)", "(fAmplitude_V0C)",
                                        "(fAmplitude_V0AEq)+(fAmplitude_V0CEq)", "(fAmplitude_V0AEq)", "(fAmplitude_V0CEq)",
                                        "(fAmplitude_OnlineV0A)+(fAmplitude_OnlineV0C)", "(fAmplitude_OnlineV0A)", "(fAmplitude_OnlineV0C)",
                                        "(fMultiplicity_ADA)+(fMultiplicity_ADC)", "(fMultiplicity_ADA)", "(fMultiplicity_ADC)",
                                        "(fnSPDClusters0)", "(fnSPDClusters1)", "(fnSPDClusters)", "(fnTracklets)",
                                        "(fRefMultEta5)", "(fRefMultEta8)",
                                        "-(fZnaFired) * (fZnaTower) + !(fZnaFired) * 1e6", "-(fZncFired) * (fZncTower) + !(fZncFired) * 1e6",
                                        "-0.89 * (fZnaFired) * (fZnaTower) - (fZncFired) * (fZncTower) + !(fZnaFired) * !(fZncFired) * 1e6" };
  TString estimatorFullString[21];
  Bool_t isIntegerEst[21]           = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE };
  Bool_t enableEst[21]              = { kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE,
                                        kTRUE, kTRUE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
                                        kTRUE  };
  Bool_t enableZVtxCorr[21]         = { kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE,
                                        kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
                                        kFALSE };

  //============================================================
  // --- Load Z vtx position correction functions
  //============================================================
  TFile* fileZVtxCorr     = new TFile(fileNameZVtxCorr.Data(), "OPEN");
  TF1* fitCorrZVtx[21]              = {NULL};
  for (Int_t i = 0; i < 21; i++){
    fitCorrZVtx[i]        = (TF1*)fileZVtxCorr->Get( Form("%s/ZVtxCorr_%s", lPeriodName.Data(),estimatorNameFit[i].Data() ));
    if (fitCorrZVtx[i]){
      cout << fitCorrZVtx[i]->GetExpFormula() << endl;
      for (Int_t k = 0; k < fitCorrZVtx[i]->GetNpar(); k++){
        cout << k << "\t" << fitCorrZVtx[i]->GetParameter(k) << endl;
      }
    }
  }
  //============================================================
  // --- Definition of Estimators ---
  //============================================================
  for (Int_t i = 0; i < 21; i++){
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
    cout << "no period yet which needs run-range calib defined" << endl;
  }


  //============================================================
  // --- Definition of Input/Output ---
  //============================================================

  lCalib -> SetInputFile  ( Form("%s", inputFileName.Data()) );
  lCalib -> SetBufferFile ( Form("buffer-%s-%s%s.root", lPeriodName.Data(), lWhichTrigger.Data(), lAddName.Data()) );
  lCalib -> SetOutputFile ( Form("OADB-%s-%s%s.root", lPeriodName.Data(), lWhichTrigger.Data(), lAddName.Data()) );
  lCalib -> Calibrate     ();

}
