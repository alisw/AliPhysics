#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TString.h"

AliAnalysisTask *AddTaskJHOCFAMaster(TString taskName = "JHOCFAMaster", UInt_t period = 0,
    std::string configArray = "0 1 3 4 6", int whichNUAmap = 3,
    double ptMin = 0.2, double ptMax = 5.0,
    bool saveFullQA = false,
    bool cutESDpileup = true, double ESDintercept = 15000,
    bool cutTPCpileup = false, bool saveQA_TPCpileup = false,
    bool useEtaGap = true, float etaGap = 1.0,
    bool useWeightsNUE = true, bool useWeightsNUA = false,
    bool useWeightsCent = false,
    bool getSC = true, bool getLower = true,
    bool Aside = false, bool Cside = false, bool saveQCNUA = false, bool frmBadArea18q = kFALSE,float eta = 0.8)
{
  // Configuration of the analysis.
  double ESDslope = 3.38; bool saveQA_ESDpileup = false;
  bool removeBadArea = kFALSE; bool useTightCuts = kFALSE;
  int debug = 0;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // Prepare the list of selection configurations for the wagons.
  enum { lhc15o = 0, lhc18q = 1, lhc18r = 2 };
  TString sPeriod[3] = {"15o", "18q", "18r"};   // Needed to load the correct map config.
  std::cout << "AddTaskJHOCFAMaster:: taskName = " << taskName << "\nPeriod = " << period 
    << std::endl;

  int iConfig = -1;
  int iOldConfig = -2;
  int index = 0;
  std::vector<TString> configNames;
  std::istringstream sConfig(configArray);

  do {
    sConfig >> iConfig;
    if (iOldConfig == iConfig) {break;}

    switch(iConfig) { // Hardcoded names to prevent typo in phi weights files.
    case 0 :    // Default selection.     // V0M + |zVtx < 8| + (pileup > 15000)
      configNames.push_back("default");   // + global tracks 96 + (NTPC < 70) + (chi2 in [0.1, 4]).
      break;
    case 1 :    // Syst: global changed to hybrid.
      configNames.push_back("hybrid");
      break;
    case 2 :    // Syst: V0M changed to SPD clusters.
      configNames.push_back("SPD");
      break;
    case 3 :    // Syst: (pileup > 15000) changed to (no pileup cut).
      configNames.push_back("noPileup");
      break;
    case 4 :    // Syst: (pileup > 15000) changed to (pileup > 10000).
      configNames.push_back("pileup10");
      break;
    case 5 :    // Syst: |zVtx < 8| changed to |zVtx < 9|.
      configNames.push_back("zvtx9");
      break;
    case 6 :    // Syst: |zVtx < 8| changed to |zVtx < 6|.
      configNames.push_back("zvtx6");
      break;
    case 7 :    // Syst: |zVtx < 8| changed to |zVtx < 7|.
      configNames.push_back("zvtx7");
      break;
    case 8 :    // Syst: (NTPC > 70) changed to (NTPC > 80).
      configNames.push_back("NTPC80");
      break;
    case 9 :    // Syst: (NTPC > 70) changed to (NTPC > 90).
      configNames.push_back("NTPC90");
      break;
    case 10 :    // Syst: (NTPC > 70) changed to (NTPC > 100).
      configNames.push_back("NTPC100");
      break;
    case 11 :    // Syst: (chi2 in [0.1, 4]) changed to (chi2 < 4).
      configNames.push_back("chi2def");
      break;
    case 12 :    // Syst: (chi2 in [0.1, 4]) changed to (chi2 < 2.5).
      configNames.push_back("chi2tight");
      break;
    case 13 :     // Syst: (DCAz < 2cm - default in global) changed to (DCAz < 1cm).
      configNames.push_back("DCAz1");
      break;
    case 14 :     // Syst: (DCAz < 2cm - default in global) changed to (DCAz < 0.5cm).
      configNames.push_back("DCAz05");
      break;
    case 15 :     // Syst: (all charges) changed to (negative charges only).
      configNames.push_back("nqq");
      break;
    case 16 :     // Syst: (all charges) changed to (positive charges only).
      configNames.push_back("pqq");
      break;
    case 17 :     // Syst: subA. TBI
      configNames.push_back("subA");
      break;
    case 18 :     // Syst: (chi2 in [0.1, 4]) changed to (chi2 in [0.1, 2.5]).
      configNames.push_back("chi2tight01");
      break;
    case 19 :     // Syst: (chi2 in [0.1, 4]) changed to (chi2 < 2.3).
      configNames.push_back("chi2tight23");
      break; 
    case 20 :     // Syst: (chi2 in [0, 1.36])
      configNames.push_back("chi2tight136");
      break;
    case 21 :     // Syst: (chi2 in [0.1, 1.36])
      configNames.push_back("chi2low136");
      break;  
    case 22 :     // Syst: (chi2 in [0.1, 2.3])
      configNames.push_back("chi2low23");
      break;
    case 23 :
      configNames.push_back("hybridBaseDCA");
      break;
    case 24 :
      configNames.push_back("DCAz15");
      break;
    case 25 :
      configNames.push_back("pileup500");
      break; 
    case 26 :
      configNames.push_back("NTPC65");
      break; 
    case 27 :    // Syst: |zVtx < 8| changed to |zVtx < 4|. In order to check if the tracking quality make a difference to our measurements
      configNames.push_back("zvtx4");
      break;
    default :
      std::cout << "ERROR: Invalid configuration index. Skipping this element."
        << std::endl;
    }
    index++;
    iOldConfig = iConfig;
  } while (sConfig);
  std::cout << "List of selection configuration has been prepared." << std::endl;

  // Get and load the NUE and NUA correction maps.
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask("JCorrectionMapTask");

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

  const int Nsets = index;
  TString MAPfileNames[Nsets];  // NUA correction maps.
  TString MAPdirName = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
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
    case 3:   // 3: Coarse binning, full PU cuts (15000), minPt = 0.2 for all.
      if (strcmp(configNames[i].Data(), "default") == 0) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_fullPUcuts_Default_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      } else if ((strcmp(configNames[i].Data(), "zvtx6") == 0)) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_fullPUcuts_s_zvtx10.root",
          MAPdirName.Data(), sCorrection[period].Data());
      } else {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_fullPUcuts_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      }
      break;
    case 4:   // Same as case 3 but with PU=500
      if (strcmp(configNames[i].Data(), "zvtx7") == 0) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_PUcuts500_s_zvtx9.root",
          MAPdirName.Data(), sCorrection[period].Data());
      } else {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_PUcuts500_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data());
      }
      break;
    case 5:   // Same as case 3 but for 18q
      if ((strcmp(configNames[i].Data(), "chi2low23") == 0)) {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_pt02_Chi2low23_s_chi2low23.root",
          MAPdirName.Data(), sCorrection[period].Data());
      } else {
        MAPfileNames[i] = Form("%sPhiWeights_LHC%s_pt02_%s_s_%s.root",
          MAPdirName.Data(), sCorrection[period].Data(), configNames[i].Data(),configNames[i].Data());
      }
      break;
    default:
      std::cout << "ERROR: Invalid configuration index. Skipping this element."
        << std::endl;   
    }

    cMapTask->EnablePhiCorrection(i, MAPfileNames[i]);  // i: index for 'SetPhiCorrectionIndex(i)'.
  }

  mgr->AddTask((AliAnalysisTask *) cMapTask);
  std::cout << "The NUA/NUE correction maps have been added to the manager." << std::endl;

  // Configure the catalyst task for each prepared configuration.
  // 'taskName' added in the name of the catalyst to prevent merging issues between wagons.
  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  int globalCut = 96;       // Global tracks - default.
  int hybridCut = 768;      // Global hybrid tracks.

  /// Trigger - common to all configurations.
  UInt_t selEvt;
  if (period == lhc15o) {   // Minimum bias only.
    selEvt = AliVEvent::kINT7;
  } else if (period == lhc18q || period == lhc18r) {  // Minimum bias + central + semicentral.
    selEvt = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  }

  if ((period == lhc18q || period == lhc18r ) && frmBadArea18q ) fJCatalyst[i]->SetRemoveBadArea18q(frmBadArea18q); 


  for (int i = 0; i < Nsets; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", 
      taskName.Data(), configNames[i].Data()));
    std::cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << std::endl;
    fJCatalyst[i]->SetSaveAllQA(saveFullQA);
    fJCatalyst[i]->SetSaveQCNUA(saveQCNUA);

    // Set the correct flags to use.
    if (strcmp(configNames[i].Data(), "noPileup") != 0) {     // Set flag only if we cut on pileup.
      fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);

      if (strcmp(configNames[i].Data(), "pileup10") == 0) {   // Vary the cut on the ESD pileup.
        fJCatalyst[i]->SetESDpileupCuts(true, ESDslope, 10000, saveQA_ESDpileup);
      } else if (strcmp(configNames[i].Data(), "pileup500") == 0) {   // Vary the cut on the ESD pileup.
        fJCatalyst[i]->SetESDpileupCuts(true, ESDslope, 500, saveQA_ESDpileup);
      } else {fJCatalyst[i]->SetESDpileupCuts(cutESDpileup, ESDslope, ESDintercept, saveQA_ESDpileup);}

      fJCatalyst[i]->SetTPCpileupCuts(cutTPCpileup, saveQA_TPCpileup); // Reject the TPC pileup.
    }
    if (period == lhc18q || period == lhc18r) {
      fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CENT_FLATTENING);
    }

    // Set the trigger, centrality and event selection.
    fJCatalyst[i]->SelectCollisionCandidates(selEvt);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();
    if (strcmp(configNames[i].Data(), "SPD") == 0) {
      fJCatalyst[i]->SetCentDetName("CL1");
    } else {  // Default: V0M.
      fJCatalyst[i]->SetCentDetName("V0M");
    } 

    if (strcmp(configNames[i].Data(), "zvtx9") == 0) {
      fJCatalyst[i]->SetZVertexCut(9.0);
    } else if (strcmp(configNames[i].Data(), "zvtx6") == 0) {
      fJCatalyst[i]->SetZVertexCut(6.0);
    } else if (strcmp(configNames[i].Data(), "zvtx7") == 0) {
      fJCatalyst[i]->SetZVertexCut(7.0);
    } else if (strcmp(configNames[i].Data(), "zvtx4") == 0) {
      fJCatalyst[i]->SetZVertexCut(4.0);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetZVertexCut(8.0);
    }

    /// Filtering and detector cuts.
    if (strcmp(configNames[i].Data(), "hybrid") == 0 || strcmp(configNames[i].Data(), "hybridBaseDCA") == 0) {
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
    }  else if (strcmp(configNames[i].Data(), "NTPC65") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(65);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetNumTPCClusters(70);
    }

    if (strcmp(configNames[i].Data(), "chi2def") == 0) {    
      fJCatalyst[i]->SetChi2Cuts(0.0, 4.0);
    } else if (strcmp(configNames[i].Data(), "chi2tight") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.0, 2.5);
    } else if (strcmp(configNames[i].Data(), "chi2tight01") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.1, 2.5);
    } else if (strcmp(configNames[i].Data(), "chi2tight23") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.0, 2.3);
    } else if (strcmp(configNames[i].Data(), "chi2tight136") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.0, 1.36);
    } else if (strcmp(configNames[i].Data(), "chi2low136") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.1, 1.36);
    } else if (strcmp(configNames[i].Data(), "chi2low23") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.1, 2.3);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetChi2Cuts(0.1, 4.0);
    }

    if (strcmp(configNames[i].Data(), "DCAz1") == 0) {    
      fJCatalyst[i]->SetDCAzCut(1.0);
    } else if (strcmp(configNames[i].Data(), "DCAz05") == 0) {
      fJCatalyst[i]->SetDCAzCut(0.5);
    } else if (strcmp(configNames[i].Data(), "DCAz15") == 0) {
      fJCatalyst[i]->SetDCAzCut(1.5);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetDCAzCut(2.0);
    }

    if (strcmp(configNames[i].Data(), "nqq") == 0) {
      fJCatalyst[i]->SetParticleCharge(-1);
    } else if (strcmp(configNames[i].Data(), "pqq") == 0) {
     fJCatalyst[i]->SetParticleCharge(1);
    }   // Default: charge = 0 to accept all charges.

    if (strcmp(configNames[i].Data(), "hybridBaseDCA") == 0) {
      fJCatalyst[i]->SetDCABaseCuts(true);
    }

    /// Kinematic cuts and last fine tuning.
    fJCatalyst[i]->SetPtRange(ptMin, ptMax);
    if (Aside){
      fJCatalyst[i]->SetEtaRange(0.0,eta);
    } else if (Cside){
      fJCatalyst[i]->SetEtaRange(-eta,0.0);
    } else {
      fJCatalyst[i]->SetEtaRange(-eta, eta);
    }
    fJCatalyst[i]->SetPhiCorrectionIndex(i);
    fJCatalyst[i]->SetRemoveBadArea(removeBadArea);
    fJCatalyst[i]->SetTightCuts(useTightCuts);
    if (period == lhc18q || period == lhc18r) {useWeightsCent = false;} // Security for 18qr.
    if (useWeightsCent) {   // Centrality weight correction for LHC15o.
      TString centWeightFile = Form(
        "alien:///alice/cern.ch/user/c/cimordas/CentWeights/CentralityWeights_LHC15oPass2_%s.root", configNames[i].Data());
      fJCatalyst[i]->SetInputCentralityWeight15o(true, centWeightFile);
      printf("Centrality weight will be used.\n");
    }

    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  } // Go to the next configuration.

  // Configure the analysis task wagons.
  AliJHOCFATask *myTask[Nsets];
  for (int i = 0; i < Nsets; i++) {
    myTask[i] = new AliJHOCFATask(Form("%s_s_%s", 
      taskName.Data(), configNames[i].Data()));
    myTask[i]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
    myTask[i]->HOCFASetDebugLevel(debug);

    myTask[i]->HOCFASetCentralityBinning(9);
    myTask[i]->HOCFASetCentralityArray("0. 5. 10. 20. 30. 40. 50. 60. 70. 80.");
    myTask[i]->HOCFASetMinMultiplicity(10);

    myTask[i]->HOCFASetPtRange(ptMin, ptMax);
    myTask[i]->HOCFASetEtaGap(useEtaGap, etaGap);
    myTask[i]->HOCFASetParticleWeights(useWeightsNUE, useWeightsNUA);
    myTask[i]->HOCFASetCentralityWeights(useWeightsCent);    

    myTask[i]->HOCFASetObservable(getSC, getLower);

    mgr->AddTask((AliAnalysisTask *)myTask[i]);
  }

  // Create the containers for input/output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[2*Nsets];

  for (int i = 0; i < Nsets; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);
    mgr->ConnectInput(myTask[i], 0, cinput);

    jHist[i] = new AliAnalysisDataContainer();
    jHist[i] = mgr->CreateContainer(Form("%s", myTask[i]->GetName()), 
      TList::Class(), AliAnalysisManager::kOutputContainer, 
      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(myTask[i], 1, jHist[i]);

    jHist[Nsets+i] = new AliAnalysisDataContainer();
    jHist[Nsets+i] = mgr->CreateContainer(Form("%s", fJCatalyst[i]->GetName()), 
      TList::Class(), AliAnalysisManager::kOutputContainer, 
      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(fJCatalyst[i], 1, jHist[Nsets+i]);
  }

  return myTask[0];
}
