#include <iostream>
#include <sstream>
#include <string>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <vector>
#include <TString.h>

AliAnalysisTask *AddTaskJSPCMasterRun1(TString taskName = "JSPCMasterRun1", std::string Variations = "tpconly",
                                  double ptMin = 0.2, Bool_t saveCatalystQA = kFALSE,
                                  Bool_t applyHMOcut = kTRUE, 
                                  bool UseITS = false, double DCAzHybrid = 3.2, double DCAxyHybrid = 2.4,
                                  Bool_t ComputeEtaGap = kFALSE, Float_t EtaGap = 0.8,
                                  Bool_t useWeightsNUE = kTRUE, Bool_t useWeightsNUA = kFALSE,
                                  Bool_t newWeightNaming = kTRUE, Bool_t UseAlternativeWeights = kFALSE,
                                  TString AlternativeWeightFile = "alien:///alice/cern.ch/user/m/mlesch/Weights/WeightsLHC10h_Filter768_OnlyPrimaries_vAN-20201209.root",
                                  Int_t doSPC = 0)
{
  Bool_t saveHMOQA = kFALSE;
  int debug = 0;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //Explanation: 
  //doSPC  	0: normal SPC
  //		1: test v2 and SC
  // 		2: test rho 

  if(doSPC < 0 || doSPC > 2) return 0;

  //-------- Read in passed Variations -------- 
  std::istringstream iss(Variations);

  const int NPossibleVariations = 14;
  std::string PossibleVariations[NPossibleVariations] = { "tpconly","hybrid", "V0M","zvtx6","zvtx12","zvtx", "chi03",
                                                          "chi35", "dcaxy1", "dcaz2", "ntpc80", "ntpc90", "ntpc100", "global"}; 
                                                          //for reference, these variations are allowed
    
  int PassedVariations = 0;
  std::vector<TString> configNames;

  do {
        std::string subs;
        iss >> subs;
  
        // Check if valid variation
        bool exists = std::find(std::begin(PossibleVariations), std::end(PossibleVariations), subs) != std::end(PossibleVariations); 
	if(exists)
	{
		PassedVariations++;
		configNames.push_back(subs);
	}

    } while (iss);

  if(PassedVariations == 0) return 0; //Protection in case no valid Variation is passed

  //-------- JFlucWagons -------
  const int Nsets  = PassedVariations; // number of configurations // TBC: if this does not work, then do const int Nsets  = 10; //default max number of variations
 
  // Loading of the correction map.
  TString MAPfilenames[Nsets];
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask ("JCorrectionMapTask");

  // Corrections fetching 
  cMapTask->EnableEffCorrection ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC10h-LHC11a10a_bis-0-Lists.root"); // Efficiency correction.

  for (Int_t i = 0; i < PassedVariations; i++) {

    if (strcmp(configNames[i].Data(), "global") == 0) { continue; } //no map for global

    if (newWeightNaming) {
      MAPfilenames[i] = Form("%sPhiWeights_LHC10h_%s_pt%02d_9904.root", MAPdirname.Data(), configNames[i].Data(), Int_t (ptMin * 10));  // Azimuthal correction.
    }
    else {
      MAPfilenames[i] = Form("%sPhiWeights_LHC10h_Error_pt%02d_s_%s.root", MAPdirname.Data(), Int_t (ptMin * 10), configNames[i].Data());  // Azimuthal correction.
    }
    cMapTask->EnablePhiCorrection (i, MAPfilenames[i]); // i is index for set file correction ->SetPhiCorrectionIndex(i);
  }
  
  mgr->AddTask((AliAnalysisTask *) cMapTask); // Loaded correction map added to the analysis manager.


// Setting of the general parameters.
  Int_t tpconlyCut = 128;
  Int_t hybridCut = 768;
  Int_t globalCut = 96;

  /// Trigger - common to all configurations.
  UInt_t selEvt;
  selEvt = AliVEvent::kMB; // Minimum bias trigger for LHC10h

  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (Int_t i = 0; i < PassedVariations; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", taskName.Data(), configNames[i].Data()));
    cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << endl;
    if (applyHMOcut) {fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);}
    fJCatalyst[i]->SelectCollisionCandidates(selEvt);
   
    fJCatalyst[i]->SetSaveAllQA(saveCatalystQA);
    fJCatalyst[i]->SetSaveHMOhist(saveHMOQA);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();

    if (strcmp(configNames[i].Data(), "V0M") == 0) {    
      fJCatalyst[i]->SetCentDetName("V0M"); // V0M for syst
    } else {
      fJCatalyst[i]->SetCentDetName("CL1"); // SPD clusters in the default analysis.
    }
    if (strcmp(configNames[i].Data(), "hybrid") == 0) {
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
      if(UseAlternativeWeights){
        fJCatalyst[i]->SetInputAlternativeNUAWeights10h(kTRUE, AlternativeWeightFile);
      }
      fJCatalyst[i]->SetDCAxyCut(DCAxyHybrid); 
      fJCatalyst[i]->SetDCAzCut(DCAzHybrid); 
      fJCatalyst[i]->SetITSCuts(UseITS,2);
    } else if (strcmp(configNames[i].Data(), "global") == 0) {
      fJCatalyst[i]->SetTestFilterBit(globalCut);
      if(UseAlternativeWeights){
        fJCatalyst[i]->SetInputAlternativeNUAWeights10h(kTRUE, AlternativeWeightFile);
      }
    } else {
      fJCatalyst[i]->SetTestFilterBit(tpconlyCut); // default
    }

    if (strcmp(configNames[i].Data(), "zvtx6") == 0) {
      fJCatalyst[i]->SetZVertexCut(6.);
    }
    if (strcmp(configNames[i].Data(), "zvtx") == 0) {
      fJCatalyst[i]->SetZVertexCut(8.);
    }
    if (strcmp(configNames[i].Data(), "zvtx12") == 0) {
      fJCatalyst[i]->SetZVertexCut(12.);
    }

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

    fJCatalyst[i]->SetPtRange(ptMin, 5.0);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    if (strcmp(configNames[i].Data(), "global") != 0) { fJCatalyst[i]->SetPhiCorrectionIndex(i); } //enable map only if we don't have global
   
    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configuration of the analysis task itself.
  const int SPCCombination = 3;
  TString SPC[SPCCombination] = { "2SPC", "3SPC", "4SPC" };
  AliJSPCTaskRun2 *myTask[Nsets][SPCCombination];

  for (Int_t i = 0; i < PassedVariations; i++){
    if (doSPC == 0) {  
      for(Int_t j = 0; j < SPCCombination; j++){
        myTask[i][j] = new AliJSPCTaskRun2(Form("%s_%s_%s", taskName.Data(), configNames[i].Data(), SPC[j].Data()));
      	myTask[i][j]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());

      	myTask[i][j]->AliSPCRun2SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
      	myTask[i][j]->AliSPCRun2SetSaveAllQA(kTRUE);
      	myTask[i][j]->AliSPCRun2SetMinNuPar(14.);
      	myTask[i][j]->AliSPCRun2SetUseWeights(useWeightsNUE, useWeightsNUA);
        myTask[i][j]->AliSPCRun2SetEtaGaps(ComputeEtaGap, EtaGap);

      	if(j==0){

      	  Int_t harmonicArray1[10][8] = {
                                          {4, 6,-2,-2,-2, 0, 0, 0},
                                          {3, 6,-3,-3, 0, 0, 0, 0},
                                          {3, 4,-2,-2, 0, 0, 0, 0},
                                          {3, 8,-4,-4, 0, 0, 0, 0},
                                          {5, 3, 3,-2,-2,-2, 0, 0},
                                          {5, 8,-2,-2,-2,-2, 0, 0},
                                          {6, 2, 2, 2, 2, -4, -4, 0},
                                          {0, 4, 4, 4, -3, -3, -3, -3},  // Too heavy for the moment.
                                          {0, 5, 5, -2, -2, -2, -2, -2}, // Too heavy for the moment.
                                          {0, 3, 3, 3, 3, -4, -4,-4}     // Too heavy for the moment.
                                        };

            for (int k = 0; k<10; k++){
              myTask[i][j]->AliSPCRun2SetCorrSet(k,harmonicArray1[k]);
            }
      	}

      	if(j==1){

      	  Int_t harmonicArray2[10][8] = {
                                          {4, 3,-4,-4, 5, 0, 0,0},
                                          {3, 2, 4,-6, 0, 0, 0,0},
                                          {3, 2, 3,-5, 0, 0, 0,0},
                                          {4, 2,-3,-3, 4, 0, 0,0},
                                          {5, 2, 3, 3,-4,-4, 0,0},
                                          {6, 2, 2, 2, 2,-3,-5,0},
                                          {5, 3, 3, 3,-4,-5, 0,0},
                                          {0, 2, 2, 2, 2, 2, -4,-6}, // Too heavy for the moment.
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0}
                                        };


            for (int k = 0; k<10; k++){
              myTask[i][j]->AliSPCRun2SetCorrSet(k,harmonicArray2[k]);
            }
      	}

      	if(j==2){

          Int_t harmonicArray3[10][8] = {       
                                          {4, 2,-3,-4, 5, 0, 0,0}, 
                                          {6, 2, 2, 2, 3,-4,-5,0}, 
                                          {5, 2, 2,-3, 4,-5, 0,0},
                                          {6, 2, 2, 3, 3,-4,-6,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0},
                                          {0, 0, 0,0, 0, 0, 0,0}
                                        };

            for (int k = 0; k<10; k++){
              myTask[i][j]->AliSPCRun2SetCorrSet(k,harmonicArray3[k]);
            }

      	}

        mgr->AddTask((AliAnalysisTask *) myTask[i][j]);
      }
    } // if (doSPC)
  }

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[Nsets][SPCCombination+1];

  for (Int_t i = 0; i < PassedVariations; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);

    if (doSPC == 0) {
      for(Int_t j = 0; j < SPCCombination; j++){
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
    }
    else if (doSPC == 1) {

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
    else if (doSPC == 2) {
      for(Int_t j = 0; j < 2; j++){
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
    }


  }

  return myTask[0][0];
}

