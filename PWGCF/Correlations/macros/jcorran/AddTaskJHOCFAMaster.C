#include "TString.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

AliAnalysisTask *AddTaskJHOCFAMaster(TString taskName = "JHOCFAMaster", UInt_t period = 0, double ptMin = 0.2, double ptMax = 5.0, std::string configArray = "0 1 3 4 6", bool saveQA = kFALSE, bool removeBadArea = kFALSE, int debug = 0, bool useWeightsNUE = kTRUE, bool useWeightsNUA = kFALSE, bool getSC3h = kTRUE, bool getEtaGap = kFALSE, float etaGap = 1.0, int Ncombi = 6, TString combiArray = "2 3 4 2 3 5 2 3 6 2 4 5 2 4 6 3 4 5")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

// Prepare the configuration of the wagons.
  enum { lhc15o = 0, lhc18q = 1, lhc18r = 2 };
  TString speriod[3] = { "15o", "18q", "18r" };   // Needed to load correct map config.

  cout << "AddTaskJHOCFAMaster:: period =" << period << "\t pT range = (" << ptMin << "," << ptMax << ")." << endl;

  int iConfig = -1;
  int iOldConfig = -2;
  int index = 0;
  std::vector<TString> configNames;
  std::istringstream sConfig(configArray);

  do {
    sConfig >> iConfig;
    if (iOldConfig == iConfig) {break;}

    switch(iConfig) { // Hardcoded names to prevent typo in phi weights files.
    case 0 :  // Default: hybrid.
      configNames.push_back("hybrid");
      break;
    case 1 :  // Syst: global.
      configNames.push_back("global");
      break;
    case 2 :  // Syst: nqq. TBI
      configNames.push_back("nqq");
      break;
    case 3 :  // Syst: pileup.
      configNames.push_back("pileup");
      break;
    case 4 :  // Syst: SPD.
      configNames.push_back("SPD");
      break;
    case 5 :  // Syst: subA. TBI
      configNames.push_back("subA");
      break;
    case 6 :  // Syst: ztx < 10.
      configNames.push_back("zvtx");
      break;
    case 7 :  // Syst: pqq. TBI.
      configNames.push_back("pqq");
      break;
    default :
      std::cout << "ERROR: Invalid configuration index. Skipping this element."
        << std::endl;
    }
    index++;
    iOldConfig = iConfig;
  } while (sConfig);

// Load the correction maps.
  const int Nsets = index;
  TString MAPfilenames[Nsets];  // Azimuthal corrections.
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cmaptask = new AliJCorrectionMapTask("JCorrectionMapTask");
  TString sCorrection[3] = { "17i2a", "18q", "18r" }; // 17i2a for 15o.

  if (period == lhc18q || period == lhc18r) {   // 2018 PbPb datasets.
    cmaptask->EnableCentFlattening(Form("alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC%s_pass13.root", speriod[period].Data()));
    cmaptask->EnableEffCorrection(Form("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC18l8-0-Lists.root", speriod[period].Data()));
  } else if (period == lhc15o) {    // 2015 PbPb dataset.
    cmaptask->EnableEffCorrection(Form("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC16g-0-Lists.root", speriod[period].Data()));
  }

  for (int i = 0; i < Nsets; i++) {
    MAPfilenames[i] = Form("%sPhiWeights_LHC%s_Error_pt%02d_s_%s.root",
      MAPdirname.Data(), sCorrection[period].Data(), Int_t(ptMin * 10), configNames[i].Data());
    cmaptask->EnablePhiCorrection(i, MAPfilenames[i]);  // i: index for 'SetPhiCorrectionIndex(i)'.
  }
  mgr->AddTask((AliAnalysisTask *) cmaptask);

// Set the general variables.
  int hybridCut = 768;      // Global hybrid tracks.
  int globalCut = 96;       // Global tracks.
  UInt_t selEvt;            // Trigger.
  if (period == lhc15o) {   // Minimum bias.
    selEvt = AliVEvent::kINT7;
  } else if (period == lhc18q || period == lhc18r) {  // Minimum bias + central + semicentral.
    selEvt = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  }

// Configure the catalyst tasks for each config
// 'taskName' added in the name of the catalyst to prevent merging issues between wagons.
  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (int i = 0; i < Nsets; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", 
      taskName.Data(), configNames[i].Data()));
    std::cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << std::endl;
    fJCatalyst[i]->SetSaveAllQA(saveQA);

  // Set the correct flags to use.
    if (strcmp(configNames[i].Data(), "pileup") != 0) {
      fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);
    }
    if (period == lhc18q || period == lhc18r) {fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CENT_FLATTENING);}

  // Set the trigger and centrality selection.
    fJCatalyst[i]->SelectCollisionCandidates(selEvt);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();
    if (strcmp(configNames[i].Data(), "SPD") == 0) {
      fJCatalyst[i]->SetCentDetName("CL1");
    } else {  // Default: V0M.
      fJCatalyst[i]->SetCentDetName("V0M");
    }

  // Set the filtering and kinematic cuts.
    if (strcmp(configNames[i].Data(), "global") == 0) {
      fJCatalyst[i]->SetTestFilterBit(globalCut);
    } else {  // Default: Hybrid tracks.
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
    }
  
    if (strcmp(configNames[i].Data(), "zvtx") == 0) {    
      fJCatalyst[i]->SetZVertexCut(10.0);
    } // Else: do nothing, default is 8.

    if (strcmp(configNames[i].Data(), "nqq") == 0) {
      fJCatalyst[i]->SetParticleCharge(-1);
    } else if (strcmp(configNames[i].Data(), "pqq") == 0) {
     fJCatalyst[i]->SetParticleCharge(1);
    }   // Default: charge = 0 to accept all charges.

    fJCatalyst[i]->SetPtRange(ptMin, ptMax);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex(i);
    fJCatalyst[i]->SetRemoveBadArea(removeBadArea);
    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

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
    myTask[i]->HOCFASetParticleWeights(useWeightsNUE, useWeightsNUA);
    myTask[i]->HOCFASetObservable(getSC3h);
    myTask[i]->HOCFASetEtaGaps(getEtaGap, etaGap);
    myTask[i]->HOCFASetNumberCombi(Ncombi);
    myTask[i]->HOCFASetHarmoArray(Form("%s", combiArray.Data()));
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
