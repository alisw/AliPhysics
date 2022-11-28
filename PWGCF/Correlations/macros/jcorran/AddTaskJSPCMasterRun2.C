#include <iostream>
#include <sstream>
#include <string>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <vector>
#include <TString.h>

AliAnalysisTask *AddTaskJSPCMasterRun2(TString taskName = "JSPCMaster", UInt_t period = 0, std::string Variations = "default",
                                  int whichNUAmap = 3, double ptMin = 0.2, double ptMax = 5.0, Bool_t saveCatalystQA = kFALSE,
                                  bool cutESDpileup = true, double ESDintercept = 15000,
                                  bool cutTPCpileup = false, bool saveQA_TPCpileup = false,
                                  Bool_t ComputeEtaGap = kFALSE, Float_t EtaGap = 0.8,
                                  Bool_t useWeightsNUE = kTRUE, Bool_t useWeightsNUA = kFALSE,
                                  Int_t doSPC = 0)
{
  double ESDslope = 3.38; bool saveQA_ESDpileup = false;
  bool removeBadArea = kFALSE; bool useTightCuts = kFALSE;
  int debug = 0;
  const int maxNrComb = 12;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //Explanation: 
  //doSPC   0: 2SPC, 4SPC and 5SPC
  //        1: 3SPC

  if(doSPC < 0 || doSPC > 1) return 0;

  //-------- Read in passed Variations -------- 
  std::istringstream iss(Variations);

  const int NPossibleVariations = 18;
  std::string PossibleVariations[NPossibleVariations] = { "default","hybrid", "SPD", "noPileup",
                                                          "pileup10", "zvtx6","zvtx10","zvtx7",
                                                          "NTPC80", "NTPC90", "NTPC100",
                                                          "chi2def", "chi2tight23", "DCAz1", "DCAz05",
                                                          "nqq", "pqq", "subA"}; //for reference, these variations are allowed
    
  int PassedVariations = 0;
  std::vector<TString> configNames;

  do {
    std::string subs;
    iss >> subs;
    
    // Check if valid variation
    bool exists = std::find(std::begin(PossibleVariations), std::end(PossibleVariations), subs) != std::end(PossibleVariations); 
    if(exists) {
      PassedVariations++;
      configNames.push_back(subs);
    }

  } while (iss);

  if (PassedVariations == 0) return 0; //Protection in case no valid Variation is passed

  //-------- JFlucWagons -------
  const int Nsets  = PassedVariations; // number of configurations // TBC: if this does not work, then do const int Nsets  = 10; //default max number of variations
 
  // Loading of the correction map.
  TString MAPfileNames[Nsets];
  TString MAPdirName = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask ("JCorrectionMapTask");

  enum { lhc15o = 0, lhc18q = 1, lhc18r = 2 };
  TString sPeriod[3] = {"15o", "18q", "18r"};   // Needed to load the correct map config.
  std::cout << "AddTaskJHOCFAMaster:: taskName = " << taskName << "\nPeriod = " << period 
    << std::endl;

  // Corrections fetching 
  if (period == lhc18q || period == lhc18r) {   // 2018 PbPb datasets.
    cMapTask->EnableCentFlattening(Form(
      "alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC%s_pass13.root",
      sPeriod[period].Data() ));
    cMapTask->EnableEffCorrection(Form(
      "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC18l8-0-Lists.root",
      sPeriod[period].Data() ));
  } else if (period == lhc15o) {                // 2015 PbPb dataset.
    cMapTask->EnableEffCorrection(Form(
      "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC16g-0-Lists.root",
      sPeriod[period].Data() ));
  }

  TString sCorrection[3] = { "15o", "18q", "18r" };

  for (int i = 0; i < Nsets; i++) {
    switch (whichNUAmap) {
    case 0:   // 0: Coarse binning, minPt = 0.2 for all.
      MAPfileNames[i] = Form("%sPhiWeights_LHC%s_Error_pt02_s_%s.root",
        MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      break;
    case 1:   // 1: Coarse binning, tuned minPt map.
      MAPfileNames[i] = Form("%sPhiWeights_LHC%s_Error_pt%02d_s_%s.root",
        MAPdirName.Data(), sCorrection[period].Data(), Int_t(ptMin * 10), configNames[i].Data());
      break;
    case 2:   // 2; Fine binning, minPt = 0.2 for all.
      if (strcmp(configNames[i].Data(), "default") == 0) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_Error_finerBins_Default_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      } else {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_Error_finerBins_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      }
      break;
    case 3:   // 3: Coarse binning, full PU cuts, minPt = 0.2 for all.
      if (strcmp(configNames[i].Data(), "default") == 0) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_fullPUcuts_Default_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      } else if ((strcmp(configNames[i].Data(), "zvtx6") == 0) || (strcmp(configNames[i].Data(), "zvtx7") == 0)) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_fullPUcuts_s_zvtx10.root",
          MAPdirName.Data(), sCorrection[period].Data());
      } else {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_fullPUcuts_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      }
      break;
    default:
      std::cout << "ERROR: Invalid configuration index. Skipping this element."
        << std::endl;   
    } // End switch.

    cMapTask->EnablePhiCorrection(i, MAPfileNames[i]);  // i: index for 'SetPhiCorrectionIndex(i)'.
  } // End loop over Nsets.

  // Setting of the general parameters.
  Int_t hybridCut = 768;
  Int_t globalCut = 96;

  /// Trigger - common to all configurations.
  UInt_t selEvt;
  if (period == lhc15o) {   // Minimum bias only.
    selEvt = AliVEvent::kINT7;
  } else if (period == lhc18q || period == lhc18r) {  // Minimum bias + central + semicentral.
    selEvt = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  }

  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (Int_t i = 0; i < PassedVariations; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", taskName.Data(), configNames[i].Data()));
    cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << endl;

    // Set the correct flags to use.
    if (strcmp(configNames[i].Data(), "noPileup") != 0) {     // Set flag only if we cut on pileup.
      fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);

      if (strcmp(configNames[i].Data(), "pileup10") == 0) {   // Vary the cut on the ESD pileup.
        fJCatalyst[i]->SetESDpileupCuts(true, ESDslope, 10000, saveQA_ESDpileup);
      } else {fJCatalyst[i]->SetESDpileupCuts(cutESDpileup, ESDslope, ESDintercept, saveQA_ESDpileup);}

      fJCatalyst[i]->SetTPCpileupCuts(cutTPCpileup, saveQA_TPCpileup); // Reject the TPC pileup.
    }
    if (period == lhc18q || period == lhc18r) {
      fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CENT_FLATTENING);
    }

    fJCatalyst[i]->SelectCollisionCandidates(selEvt);
   
    fJCatalyst[i]->SetSaveAllQA(saveCatalystQA);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();

    if (strcmp(configNames[i].Data(), "SPD") == 0) {
      fJCatalyst[i]->SetCentDetName("CL1");
    } else {  // Default: V0M.
      fJCatalyst[i]->SetCentDetName("V0M");
    } 

    if (strcmp(configNames[i].Data(), "zvtx10") == 0) {
      fJCatalyst[i]->SetZVertexCut(10.0);
    } else if (strcmp(configNames[i].Data(), "zvtx6") == 0) {
      fJCatalyst[i]->SetZVertexCut(6.0);
    } else if (strcmp(configNames[i].Data(), "zvtx7") == 0) {
      fJCatalyst[i]->SetZVertexCut(7.0);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetZVertexCut(8.0);
    }

    /// Filtering and detector cuts.
    if (strcmp(configNames[i].Data(), "hybrid") == 0) {
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
    } else {  // Default: global tracks.
      fJCatalyst[i]->SetTestFilterBit(globalCut);
    }

    if (strcmp(configNames[i].Data(), "NTPC80") == 0) {    
      fJCatalyst[i]->SetNumTPCClusters(80);
    } else if (strcmp(configNames[i].Data(), "NTPC90") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(90);
    } else if (strcmp(configNames[i].Data(), "NTPC100") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(100);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetNumTPCClusters(70);
    }

    if (strcmp(configNames[i].Data(), "chi2def") == 0) {    
      fJCatalyst[i]->SetChi2Cuts(0.0, 4.0);
    } else if (strcmp(configNames[i].Data(), "chi2tight23") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.1, 2.3);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetChi2Cuts(0.1, 4.0);
    }

    if (strcmp(configNames[i].Data(), "DCAz1") == 0) {    
      fJCatalyst[i]->SetDCAzCut(1.0);
    } else if (strcmp(configNames[i].Data(), "DCAz05") == 0) {
      fJCatalyst[i]->SetDCAzCut(0.5);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetDCAzCut(2.0);
    }

    if (strcmp(configNames[i].Data(), "nqq") == 0) {
      fJCatalyst[i]->SetParticleCharge(-1);
    } else if (strcmp(configNames[i].Data(), "pqq") == 0) {
     fJCatalyst[i]->SetParticleCharge(1);
    }   // Default: charge = 0 to accept all charges.

    /// Kinematic cuts and last fine tuning.
    fJCatalyst[i]->SetPtRange(ptMin, ptMax);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex(i);
    fJCatalyst[i]->SetRemoveBadArea(removeBadArea);
    fJCatalyst[i]->SetTightCuts(useTightCuts);
    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configuration of the analysis task itself.
  const int SPCCombination = 3; // We do 3SPC only if doSPC == 1 for time running sake.
  TString SPC[SPCCombination] = { "2SPC", "4SPC", "5SPC"};
  AliJSPCTaskRun2 *myTask[Nsets][SPCCombination];

  for (Int_t i = 0; i < PassedVariations; i++) {
    if (doSPC == 0) {
      for (Int_t j = 0; j < SPCCombination; j++) {
        myTask[i][j] = new AliJSPCTaskRun2(Form("%s_%s_%s", taskName.Data(), configNames[i].Data(), SPC[j].Data()));
        myTask[i][j]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
        myTask[i][j]->AliSPCRun2SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
        myTask[i][j]->AliSPCRun2SetSaveAllQA(kTRUE);
        myTask[i][j]->AliSPCRun2SetMinNuPar(14.);
        myTask[i][j]->AliSPCRun2SetUseWeights(useWeightsNUE, useWeightsNUA);
        myTask[i][j]->AliSPCRun2SetEtaGaps(ComputeEtaGap, EtaGap);

        if (j == 0) {
          Int_t harmonicArray1[maxNrComb][8] = {
                                          {4, 6,-2,-2,-2, 0, 0, 0},
                                          {3, 6,-3,-3, 0, 0, 0, 0},
                                          {3, 4,-2,-2, 0, 0, 0, 0},
                                          {3, 8,-4,-4, 0, 0, 0, 0},
                                          {5, 3, 3,-2,-2,-2, 0, 0},
                                          {5, 8,-2,-2,-2,-2, 0, 0},
                                          {6, 2, 2, 2, 2, -4, -4, 0},
                                          {0, 4, 4, 4, -3, -3, -3, -3},  // Too heavy for the moment.
                                          {0, 5, 5, -2, -2, -2, -2, -2}, // Too heavy for the moment.
                                          {0, 3, 3, 3, 3, -4, -4,-4}, // Too heavy for the moment.
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0}
                                        };

          for (int k = 0; k<maxNrComb; k++){
            myTask[i][j]->AliSPCRun2SetCorrSet(k,harmonicArray1[k]);
          }
        } // End j == 0, 2SPC.

        if (j == 1) {
          Int_t harmonicArray3[maxNrComb][8] = {       
                                          {4, 2,-3,-4, 5, 0, 0,0}, 
                                          {6, 2, 2, 2, 3,-4,-5,0}, 
                                          {5, 2, 2,-3, 4,-5, 0,0},
                                          {6, 2, 2, 3, 3,-4,-6,0},
                                          {4, 2, 7, -4, -5, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0}
                                        };

          for (int k = 0; k<maxNrComb; k++){
            myTask[i][j]->AliSPCRun2SetCorrSet(k,harmonicArray3[k]);
          }
        } // End j == 1, 4SPC.

        if (j == 2) {
          Int_t harmonicArray4[maxNrComb][8] = {
                                          {5, 2, 3, -4, 5, -6, 0,0},
                                          {6, 2, 3, 4, 4, -6, -7,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0}
                                        };
                                        
          for (int k = 0; k<maxNrComb; k++){
            myTask[i][j]->AliSPCRun2SetCorrSet(k,harmonicArray4[k]);
          }
        } // End j == 2, 5SPC.

        mgr->AddTask((AliAnalysisTask *) myTask[i][j]);
      } // End for (Int_t j = 0; j < SPCCombination; j++).
    } // End if (doSPC == 0)

    else if (doSPC == 1) {
      myTask[i][0] = new AliJSPCTaskRun2(Form("%s_%s_3SPC", taskName.Data(), configNames[i].Data()));
      myTask[i][0]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
      myTask[i][0]->AliSPCRun2SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
      myTask[i][0]->AliSPCRun2SetSaveAllQA(kTRUE);
      myTask[i][0]->AliSPCRun2SetMinNuPar(14.);
      myTask[i][0]->AliSPCRun2SetUseWeights(useWeightsNUE, useWeightsNUA);
      myTask[i][0]->AliSPCRun2SetEtaGaps(ComputeEtaGap, EtaGap);

      Int_t harmonicArray2[maxNrComb][8] = {
                                      {4, 3,-4,-4, 5, 0, 0,0},
                                      {3, 2, 4,-6, 0, 0, 0,0},
                                      {3, 2, 3,-5, 0, 0, 0,0},
                                      {4, 2,-3,-3, 4, 0, 0,0},
                                      {5, 2, 3, 3,-4,-4, 0,0},
                                      {6, 2, 2, 2, 2,-3,-5,0},
                                      {5, 3, 3, 3,-4,-5, 0,0},
                                      {4, 2, 2, 3, -7, 0, 0,0},
                                      {3, 3, 4, -7, 0, 0, 0,0},
                                      {3, 2, 5, -7, 0, 0, 0,0},
                                      {3, 3, 5, -8, 0, 0, 0,0},
                                      {4, 2, 2, 4, -8, 0, 0,0}
                                      // {0, 2, 2, 2, 2, 2, -4,-6} // Too heavy for the moment.
                                    };

      for (int k = 0; k<maxNrComb; k++){
        myTask[i][0]->AliSPCRun2SetCorrSet(k,harmonicArray2[k]);
      }

    } // End else if (doSPC == 1)
  } // End for (Int_t i = 0; i < PassedVariations; i++)

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[Nsets][SPCCombination+1];

  for (Int_t i = 0; i < PassedVariations; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);

    if (doSPC == 0) {
      for(Int_t j = 0; j < SPCCombination; j++) {
        mgr->ConnectInput(myTask[i][j], 0, cinput);
        jHist[i][j] = new AliAnalysisDataContainer();     
        jHist[i][j] = mgr->CreateContainer(Form ("%s", myTask[i][j]->GetName()),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s:outputAnalysis", AliAnalysisManager::GetCommonFileName()));
        mgr->ConnectOutput(myTask[i][j], 1, jHist[i][j]);
      }

      jHist[i][SPCCombination] = new AliAnalysisDataContainer();
      jHist[i][SPCCombination] = mgr->CreateContainer(Form ("%s", fJCatalyst[i]->GetName()),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
      mgr->ConnectOutput(fJCatalyst[i], 1, jHist[i][SPCCombination]);

    } else if (doSPC == 1) {
      mgr->ConnectInput(myTask[i][0], 0, cinput);
      jHist[i][0] = new AliAnalysisDataContainer();     
      jHist[i][0] = mgr->CreateContainer(Form ("%s", myTask[i][0]->GetName()),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s:outputAnalysis", AliAnalysisManager::GetCommonFileName()));
      mgr->ConnectOutput(myTask[i][0], 1, jHist[i][0]);

      jHist[i][SPCCombination] = new AliAnalysisDataContainer();
      jHist[i][SPCCombination] = mgr->CreateContainer(Form ("%s", fJCatalyst[i]->GetName()),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s", AliAnalysisManager::GetCommonFileName()));
      mgr->ConnectOutput(fJCatalyst[i], 1, jHist[i][SPCCombination]);
    }
  } // End for (Int_t i = 0; i < PassedVariations; i++)

  return myTask[0][0];
}
