#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

AliAnalysisTask *AddTaskJFFlucJCMAP(TString taskName = "JFFlucJCMAP_Run2_pass2", UInt_t period = 0, std::string ptMinArray = "0.2 0.3 0.4 0.5", double ptMax = 5.0, int cutConfig = 0, bool saveQA = kFALSE, bool removeBadArea = kFALSE, int debug = 0, bool useWeightsNUE = kTRUE, bool useWeightsNUA = kFALSE, bool useTightCuts = kFALSE, bool ESDpileup = false, double slope = 3.38, double intercept = 15000, bool saveQApileup = false, bool getSC3h = kTRUE, bool getEtaGap = kFALSE, float etaGap = 1.0, int Ncombi = 6, TString combiArray = "2 3 4 2 3 5 2 3 6 2 4 5 2 4 6 3 4 5")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

// Prepare the configuration of the wagons.
  enum { lhc15o = 0, lhc18q = 1, lhc18r = 2 };
  TString speriod[3] = { "15o", "18q", "18r" };   // Needed to load correct map config.

  TString configName;   // Configuration name corresponding to the cutConfig.
  switch(cutConfig) { // Hardcoded names to prevent typo in phi weights files.
  case 0 :  // Default: global.
    configName = "global";
    break;
  case 1 :  // Syst: hybrid.
    configName = "hybrid";
    break;
  case 2 :  // Syst: SPD
    configName = "SPD";
    break;
  case 3 :  // Syst: no pileup.
    configName = "pileup";
    break;
  case 4 :  // Syst: pileup > 10000.
    configName = "pileup10";
    break;
  case 5 :  // Syst: zvtx < 10.
    configName = "zvtx";
    break;
  case 6 :  // Syst: NTPC > 80
    configName = "NTPC80";
    break;
  case 7 :  // Syst: NTPC > 90
    configName = "NTPC90";
    break;
  case 8 :  // Syst: NTPC > 100
    configName = "NTPC100";
    break;
  case 9 :  // Syst: chi2 TPC default.
    configName = "chi2def";
    break;
  case 10 :  // Syst: chi2 < 2.5.
    configName = "chi2tight";
    break;
  case 11 :  // Syst: DCAz < 1cm.
    configName = "DCAz1";
    break;
  case 12 :  // Syst: DCAz < 0.5cm.
    configName = "DCAz05";
    break;
  case 13 :  // Syst: zvtx < 7.
    configName = "zvtx7";
    break;
  case 14 :  // Syst: zvtx < 9.
    configName = "zvtx9";
    break;
  case 20 :  // Syst: nqq. TBI
    configName = "nqq";
    break;
  case 21 :  // Syst: pqq. TBI.
    configName = "pqq";
    break;
  case 22 :  // Syst: subA. TBI
    configName = "subA";
    break;
  default :
    std::cout << "ERROR: Invalid configuration index." << std::endl;
  }

  std::cout << "AddTaskJFFlucJCMAP:: period = " << period << "\t max pT = " << ptMax << std::endl;
  std::cout << "Config of the selection = " << configName.Data() << std::endl;

  // Prepare the array of min pT values from the provided string.
  int index = 0;
  double thisMinPt = -1.;
  double thisOldMinPt = -2.;
  std::vector<double> configMinPt;
  std::istringstream sConfig(ptMinArray);

  do {
    sConfig >> thisMinPt;
    if ( TMath::Abs(thisOldMinPt - thisMinPt) < 1e-6 ) {break;}

    configMinPt.push_back(thisMinPt); 
    index++;
    thisOldMinPt = thisMinPt;
  } while (sConfig);

// Load the correction maps.
// We assume the same maps for all minPt values.
/*
  TString MAPfilenames;  // Azimuthal corrections.
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cmaptask = new AliJCorrectionMapTask("JCorrectionMapTask");
  TString sCorrection[3] = { "15o", "18q", "18r" }; // 17i2a for 15o?

  if (period == lhc18q || period == lhc18r) {   // 2018 PbPb datasets.
    cmaptask->EnableCentFlattening(Form("alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC%s_pass13.root", speriod[period].Data()));
    cmaptask->EnableEffCorrection(Form("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC18l8-0-Lists.root", speriod[period].Data()));
  } else if (period == lhc15o) {    // 2015 PbPb dataset.
    cmaptask->EnableEffCorrection(Form("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC16g-0-Lists.root", speriod[period].Data()));
  }  

  MAPfilenames = Form("%sPhiWeights_LHC%s_Error_pt02_s_%s.root", MAPdirname.Data(), sCorrection[period].Data(), configName.Data());
  cmaptask->EnablePhiCorrection(0, MAPfilenames);  // i = 0: index for 'SetPhiCorrectionIndex(i)'.
  mgr->AddTask((AliAnalysisTask *) cmaptask);
 */

  // Set the general variables.
  int hybridCut = 768;      // Hybrid tracks.
  int globalCut = 96;       // Global tracks.
  UInt_t selEvt;            // Trigger.
  if (period == lhc15o) {   // Minimum bias.
    selEvt = AliVEvent::kINT7;
  } else if (period == lhc18q || period == lhc18r) {  // Minimum bias + central + semicentral.
    selEvt = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  }

  // Configure the catalyst tasks for each value of minPt.
  // taskName added in the name of the catalyst to prevent merging issues between wagons.
  const int Nsets = index;
  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (int i = 0; i < Nsets; i++) {
    std::cout << "Current min pT = " << configMinPt[i] << std::endl;
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s_minPt%02d", taskName.Data(), configName.Data(), Int_t(configMinPt[i] * 10)));
    std::cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << std::endl;
    fJCatalyst[i]->SetSaveAllQA(saveQA);

  // Set the correct flags to use.
    if (strcmp(configName.Data(), "pileup") != 0) {
      fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);
      if (strcmp(configName.Data(), "pileup10") == 0) {fJCatalyst[i]->SetESDpileupCuts(true, slope, 10000, saveQApileup);}
      else {fJCatalyst[i]->SetESDpileupCuts(ESDpileup, slope, intercept, saveQApileup);}
    }
    if (period == lhc18q || period == lhc18r) {fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CENT_FLATTENING);}

  // Set the trigger and centrality selection.
    fJCatalyst[i]->SelectCollisionCandidates(selEvt);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();
    if (strcmp(configName.Data(), "SPD") == 0) {
      fJCatalyst[i]->SetCentDetName("CL1");
    } else {  // Default: V0M.
      fJCatalyst[i]->SetCentDetName("V0M");
    }

  // Set the filtering and kinematic cuts.
    if (strcmp(configName.Data(), "hybrid") == 0) {
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
    } else {  // Default: global tracks.
      fJCatalyst[i]->SetTestFilterBit(globalCut);
    }
  
    if (strcmp(configName.Data(), "zvtx") == 0) {    
      fJCatalyst[i]->SetZVertexCut(10.0);
    } else if (strcmp(configName.Data(), "zvtx7") == 0) {
      fJCatalyst[i]->SetZVertexCut(7.0);
    } else if (strcmp(configName.Data(), "zvtx9") == 0) {
      fJCatalyst[i]->SetZVertexCut(9.0);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetZVertexCut(8.0);
    }

    if (strcmp(configName.Data(), "NTPC80") == 0) {    
      fJCatalyst[i]->SetNumTPCClusters(80);
    } else if (strcmp(configName.Data(), "NTPC90") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(90);
    } else if (strcmp(configName.Data(), "NTPC100") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(100);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetNumTPCClusters(70);
    }

    if (strcmp(configName.Data(), "chi2def") == 0) {    
      fJCatalyst[i]->SetChi2Cuts(0.0, 4.0);
    } else if (strcmp(configName.Data(), "chi2tight") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.0, 2.5);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetChi2Cuts(0.1, 4.0);
    }

    if (strcmp(configName.Data(), "DCAz1") == 0) {    
      fJCatalyst[i]->SetDCAzCut(1.0);
    } else if (strcmp(configName.Data(), "DCAz05") == 0) {
      fJCatalyst[i]->SetDCAzCut(0.5);
    } else {  // Default value for JCorran analyses in Run 2.
      fJCatalyst[i]->SetDCAzCut(2.0);
    }

    if (strcmp(configName.Data(), "nqq") == 0) {
      fJCatalyst[i]->SetParticleCharge(-1);
    } else if (strcmp(configName.Data(), "pqq") == 0) {
      fJCatalyst[i]->SetParticleCharge(1);
    }   // Default: charge = 0 to accept all charges.

    fJCatalyst[i]->SetPtRange(configMinPt[i], ptMax);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex(0);  // Instead of i.
    fJCatalyst[i]->SetRemoveBadArea(removeBadArea);
    fJCatalyst[i]->SetTightCuts(useTightCuts);
    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configure the analysis task wagons.
  AliJFFlucJCTask *myTask[Nsets];
  for (int i = 0; i < Nsets; i++) {
    myTask[i] = new AliJFFlucJCTask(Form("%s_s_%s_minPt%02d", 
      taskName.Data(), configName.Data(), Int_t(configMinPt[i] * 10)));
    myTask[i]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
   
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
