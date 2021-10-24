#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <TString.h>

AliAnalysisTask *AddTaskJHOCFAMaster10h (TString taskName = "JHOCFAMaster10h", double ptMin = 0.2, double ptMax = 5.0, Bool_t removebadarea = kFALSE,
  Bool_t applyHMOcut = kTRUE, Bool_t saveCatalystQA = kFALSE, Bool_t saveHMOQA = kFALSE, std::string configArray = "0 1 6 8 9 10 11 12 13 14 15 16",
  Int_t Ncombi = 6, TString combiArray = "2 3 4 2 3 5 2 3 6 2 4 5 2 4 6 3 4 5", Bool_t newWeightNaming = kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

// Configuration of the wagons.
  int iConfig = -1;
  int iOldConfig = -2;
  int index = 0;
  std::vector<TString> configNames;
  std::istringstream sConfig(configArray);

  do {
    sConfig >> iConfig;
    if (iOldConfig == iConfig) {break;}

    switch(iConfig) { // Hardcoded names to prevent typo in file names.
      case 0 :  // Default: tpconly.
        configNames.push_back("tpconly");
        break;
      case 1 :  // Syst: hybrid.
        configNames.push_back("hybrid");
        break;
      case 2 :  // Syst: nqq. TBI
        configNames.push_back("nqq");
        break;
      case 3 :  // Syst: pileup. TBI
        configNames.push_back("pileup");
        break;
      case 4 :  // Syst: SPD. Not included in 9904.
        configNames.push_back("SPD");
        break;
      case 5 :  // Syst: subA. TBI
        configNames.push_back("subA");
        break;
      case 6 :  // Syst: V0M.
        configNames.push_back("V0M");
        break;
      case 7 :  // Syst: ztx < 8. Not included in 9904.
        configNames.push_back("zvtx");
        break;
      case 8 :  // Syst: ztx < 6.
        configNames.push_back("zvtx6");
        break;
      case 9 :  // Syst: ztx < 12.
        configNames.push_back("zvtx12");
        break;
      case 10 :  // Syst: chi^2 in [0.3, 4.0].
        configNames.push_back("chi03");
        break;
      case 11 :  // Syst: chi^2 in [0.1, 3.5].
        configNames.push_back("chi35");
        break;
      case 12 :  // Syst: dcaxy < 1cm
        configNames.push_back("dcaxy1");
        break;
      case 13 :  // Syst: dcaz < 2cm.
        configNames.push_back("dcaz2");
        break;
      case 14 :  // Syst: Ntpc > 80.
        configNames.push_back("ntpc80");
        break;
      case 15 :  // Syst: Ntpc > 90.
        configNames.push_back("ntpc90");
        break;
      case 16 :  // Syst: Ntpc > 100.
        configNames.push_back("ntpc100");
        break;
      default :
        printf("ERROR: Invalid configuration index. Skipping this element.\n");
    }
    index++;
    iOldConfig = iConfig;
  } while (sConfig);
  const int Nsets  = index;

// Loading of the correction map.
  TString MAPfilenames[Nsets];
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask ("JCorrectionMapTask");

  cMapTask->EnableEffCorrection ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC10h-LHC11a10a_bis-0-Lists.root"); // Efficiency correction.
  for (Int_t i = 0; i < Nsets; i++) {
    if (newWeightNaming) {
      MAPfilenames[i] = Form("%sPhiWeights_LHC10h_%s_pt%02d_9904.root", MAPdirname.Data(), configNames[i].Data(), Int_t (ptMin * 10));  // Azimuthal correction.
    }
    else {
      MAPfilenames[i] = Form("%sPhiWeights_LHC10h_Error_pt%02d_s_%s.root", MAPdirname.Data(), Int_t (ptMin * 10), configNames[i].Data());  // Azimuthal correction.
    }
    cMapTask->EnablePhiCorrection(i, MAPfilenames[i]); // i is index for set file correction ->SetPhiCorrectionIndex(i);
  }
  mgr->AddTask((AliAnalysisTask *) cMapTask); // Loaded correction map added to the analysis manager.

// Setting of the general parameters.
  UInt_t configTrigger = AliVEvent::kMB; // Minimum bias trigger for LHC10h.
  Int_t TPConlyFilter = 128;  // Filterbit for TPConly tracks.
  Int_t hybridCut = 768;

// Configuration of the catalyst tasks for the data selection.
  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (Int_t i = 0; i < Nsets; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", taskName.Data(), configNames[i].Data()));
    cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << endl;

    // Ensure that the outlier cut is not on for the pileup systematics.
    if (strcmp(configNames[i].Data(), "pileup") == 0) {applyHMOcut = kFALSE;}
    if (applyHMOcut) {fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);}
    
    fJCatalyst[i]->SelectCollisionCandidates(configTrigger);
    fJCatalyst[i]->SetSaveAllQA(saveCatalystQA);
    fJCatalyst[i]->SetSaveHMOhist(saveHMOQA);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();

    if (strcmp(configNames[i].Data(), "V0M") == 0) {  
      fJCatalyst[i]->SetCentDetName("V0M"); // V0M for syst
    }
    else {
      fJCatalyst[i]->SetCentDetName("CL1"); // SPD clusters in the default analysis.
    }

    if (strcmp(configNames[i].Data(), "hybrid") == 0) {
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
    }
    else {
      fJCatalyst[i]->SetTestFilterBit(TPConlyFilter); // default
    }

    if (strcmp(configNames[i].Data(), "zvtx6") == 0) {
      fJCatalyst[i]->SetZVertexCut(6.0);
    }
    else if (strcmp(configNames[i].Data(), "zvtx12") == 0) {
      fJCatalyst[i]->SetZVertexCut(12.0);
    }
    else if (strcmp(configNames[i].Data(), "zvtx") == 0) {
      fJCatalyst[i]->SetZVertexCut(8.0);  // Not in 9904.
    } // Else: do nothing, default is 10.

    if (strcmp(configNames[i].Data(), "nqq") == 0) {
      fJCatalyst[i]->SetParticleCharge(-1);
    } // Default: charge = 0 to accept all charges.
  // TBA: Systematics for positive charges.
    
/*
    if (strcmp(configNames[i].Data(), "subA") == 0) {
      fJCatalyst[i]->
    } // Default: charge = 0 to accept all charges.
  // TBA: Systematics for positive charges.
*/

    if (strcmp(configNames[i].Data(), "chi03") == 0) {  
      fJCatalyst[i]->SetChi2Cuts(0.3, 4.0); // chi2 in [0.3, 4.0] for syst.
    }
    else if (strcmp(configNames[i].Data(), "chi35") == 0) {
      fJCatalyst[i]->SetChi2Cuts(0.1, 3.5); // chi2 in [0.1, 3.5] for syst.
    } // Else: do nothing, default is [0.1, 4.0].

    if (strcmp(configNames[i].Data(), "dcaxy1") == 0) {  
      fJCatalyst[i]->SetDCAxyCut(1.0); // DCAxy < 1cm for syst.
    } // Else: do nothing, default is DCAxy < 2.4cm.

    if (strcmp(configNames[i].Data(), "dcaz2") == 0) {  
      fJCatalyst[i]->SetDCAzCut(2.0); // DCAz < 2cm for syst.
    } // Else: do nothing, default is DCAz < 3.2cm.

    if (strcmp(configNames[i].Data(), "ntpc80") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(80);
    }
    else if (strcmp(configNames[i].Data(), "ntpc90") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(90);
    }
    else if (strcmp(configNames[i].Data(), "ntpc100") == 0) {
      fJCatalyst[i]->SetNumTPCClusters(100);
    } // Else: do nothing, default is 70.


    fJCatalyst[i]->SetPtRange(ptMin, ptMax);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex (i);
    fJCatalyst[i]->SetRemoveBadArea(removebadarea);

    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configuration of the analysis task itself.
  AliJHOCFATask *myTask[Nsets];
  for (Int_t i = 0; i < Nsets; i++){
    myTask[i] = new AliJHOCFATask(Form("%s_s_%s", taskName.Data(), configNames[i].Data()));
    myTask[i]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
    myTask[i]->SetDebugLevel(0);
    myTask[i]->HOCFASetCentralityBinning(9);
    myTask[i]->HOCFASetCentralityArray("0. 5. 10. 20. 30. 40. 50. 60. 70. 80.");
    myTask[i]->HOCFASetMinMultiplicity(10);
    myTask[i]->HOCFASetParticleWeights(kTRUE);
    myTask[i]->HOCFASetNumberCombi(Ncombi);
    myTask[i]->HOCFASetHarmoArray(Form("%s", combiArray.Data())); // 6 combis.
    mgr->AddTask((AliAnalysisTask *) myTask[i]);
  }

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[2*Nsets];

  for (Int_t i = 0; i < Nsets; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);
    mgr->ConnectInput(myTask[i], 0, cinput);

    jHist[i] = new AliAnalysisDataContainer();
    jHist[Nsets+i] = new AliAnalysisDataContainer();
    jHist[i] = mgr->CreateContainer(Form ("%s", myTask[i]->GetName()),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(myTask[i], 1, jHist[i]);

    jHist[Nsets+i] = mgr->CreateContainer(Form ("%s", fJCatalyst[i]->GetName()),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(fJCatalyst[i], 1, jHist[Nsets+i]);
  }

  return myTask[0];
}
