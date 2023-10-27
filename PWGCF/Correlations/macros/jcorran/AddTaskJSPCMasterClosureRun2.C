#include <iostream>
#include <sstream>
#include <string>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <vector>
#include <TString.h>

AliAnalysisTask *AddTaskJSPCMasterClosureRun2(TString taskName = "JSPCMaster", double ptMin = 0.2, bool newWeightNaming = true, bool IsLHC10h = true, Int_t doSPC = 0)
{
  const int maxNrComb = 12;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //-------- JFlucWagons -------
  const int Nsets  = 3; // number of configurations // MC=0, hybrid=1 

  // Loading of the correction map. //ADDED + CHECK changed maps to 10h maps to run on 10h AMPT
  TString MAPfilenames[Nsets];
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask ("JCorrectionMapTask");

  if(IsLHC10h){
   cMapTask->EnableEffCorrection ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC10h-LHC11a10a_bis-0-Lists.root"); // Efficiency correction. 
  }
  else {  
   cMapTask->EnableEffCorrection ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC15o-LHC16g-0-Lists.root"); // Efficiency correction.
  }
  for (Int_t i = 0; i < Nsets; i++) {
    if (newWeightNaming) {
      if(IsLHC10h){
        MAPfilenames[i] = Form("%sPhiWeights_LHC10h_tpconly_pt%02d_9904.root", MAPdirname.Data(), Int_t (ptMin * 10));  // Azimuthal correction.  
      }
      else {  
        MAPfilenames[i] = Form("%sPhiWeights_LHC17i2_default_s_default.root", MAPdirname.Data());
        // MAPfilenames[i] = Form("/home/maxim/Downloads/PhiWeights_LHC17i2_default.root");
      }      
    }
    else {
      MAPfilenames[i] = Form("%sPhiWeights_LHC10h_Error_pt%02d_s_tpconly.root", MAPdirname.Data(), Int_t (ptMin * 10));  // Azimuthal correction.
    }
    cMapTask->EnablePhiCorrection(i, MAPfilenames[i]); // i is index for set file correction ->SetPhiCorrectionIndex(i);
  }
  mgr->AddTask((AliAnalysisTask *) cMapTask); // Loaded correction map added to the analysis manager.


  //=====================================================================================================================
  // Loading of the correction map.
  /*TString MAPfilenames;
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask ("JCorrectionMapTask");

  cMapTask->EnableEffCorrection ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC15o-LHC16g-0-.root"); // Efficiency correction.
  MAPfilenames = "PhiWeights_LHC17i2a_Error_pt02_s_hybrid.root";  // Azimuthal correction.
  cMapTask->EnablePhiCorrection (0, MAPfilenames); // i is index for set file correction ->SetPhiCorrectionIndex(i);
  mgr->AddTask((AliAnalysisTask *) cMapTask); // Loaded correction map added to the analysis manager.*/

// Setting of the general parameters.

  Int_t tpconlyCut = 128;
  Int_t hybridCut = 768;
  Int_t defaultCut = 96;

  Int_t UsedCut=0;
  if(IsLHC10h){
    UsedCut = tpconlyCut;
  }
  else {  
    UsedCut = defaultCut;  
  } 

  UInt_t selEvt;
  selEvt = AliVEvent::kINT7;// Minimum bias trigger for LHC15o.

  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.

  //MC True
  fJCatalyst[0] = new AliJCatalystTask("JCatalystTaskMC");
  cout << "Setting the catalyst: " << fJCatalyst[0]->GetJCatalystTaskName() << endl;
  fJCatalyst[0]->AddFlags(AliJCatalystTask::FLUC_ALICE_IPINFO|AliJCatalystTask::FLUC_MC|AliJCatalystTask::FLUC_EXCLUDEWDECAY); 
  //fJCatalyst[0]->AddFlags(AliJCatalystTask::FLUC_ALICE_IPINFO|AliJCatalystTask::FLUC_KINEONLY|AliJCatalystTask::FLUC_EXCLUDEWDECAY);
  //fJCatalyst[0]->AddFlags(AliJCatalystTask::FLUC_ALICE_IPINFO|AliJCatalystTask::FLUC_MC|AliJCatalystTask::FLUC_EXCLUDEWDECAY);

  fJCatalyst[0]->SelectCollisionCandidates(selEvt);
  fJCatalyst[0]->SetSaveAllQA(kTRUE);
  fJCatalyst[0]->SetSaveHMOhist(kFALSE);
  fJCatalyst[0]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
  fJCatalyst[0]->SetInitializeCentralityArray();
  fJCatalyst[0]->SetCentDetName("V0M"); //ADDED + CHECK changed to V0M as in default analysis 
  fJCatalyst[0]->SetTestFilterBit(UsedCut); //ADDED + CHECK changed to tpc only like in default analysis
  fJCatalyst[0]->SetZVertexCut(8.);
  fJCatalyst[0]->SetPtRange(ptMin, 5.0);
  fJCatalyst[0]->SetEtaRange(-0.8, 0.8);
  fJCatalyst[0]->SetNumTPCClusters(70);
  fJCatalyst[0]->SetDCAzCut(2.0);
  fJCatalyst[0]->SetChi2Cuts(0.1, 4.0);
  mgr->AddTask((AliAnalysisTask *)fJCatalyst[0]); 

  //MC Reco
  fJCatalyst[1] = new AliJCatalystTask("JCatalystTaskRecoWithWeights");
  cout << "Setting the catalyst: " << fJCatalyst[1]->GetJCatalystTaskName() << endl;
  //fJCatalyst[1]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS); //Add FLUC_KINEONLY???
  //fJCatalyst[1]->AddFlags(AliJCatalystTask::FLUC_ALICE_IPINFO);
  fJCatalyst[1]->SelectCollisionCandidates(selEvt); 
  fJCatalyst[1]->SetSaveAllQA(kTRUE);
  fJCatalyst[1]->SetSaveHMOhist(kFALSE);
  fJCatalyst[1]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
  fJCatalyst[1]->SetInitializeCentralityArray();
  fJCatalyst[1]->SetCentDetName("V0M"); //ADDED + CHECK changed to V0M as in default analysis 
  fJCatalyst[1]->SetTestFilterBit(UsedCut); //ADDED + CHECK changed to tpc only like in default analysis
  fJCatalyst[1]->SetZVertexCut(8.);
  fJCatalyst[1]->SetPtRange(ptMin, 5.0);
  fJCatalyst[1]->SetEtaRange(-0.8, 0.8);
  fJCatalyst[1]->SetNumTPCClusters(70);
  fJCatalyst[1]->SetDCAzCut(2.0);
  fJCatalyst[1]->SetChi2Cuts(0.1, 4.0);
  fJCatalyst[1]->SetPhiCorrectionIndex(0); //ADDED + CHECK
  mgr->AddTask((AliAnalysisTask *)fJCatalyst[1]); //ADDED + CHECK

  //MC Reco without weights
  fJCatalyst[2] = new AliJCatalystTask("JCatalystTaskRecoNoWeights");
  cout << "Setting the catalyst: " << fJCatalyst[2]->GetJCatalystTaskName() << endl;
  //fJCatalyst[2]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);
  //fJCatalyst[2]->AddFlags(AliJCatalystTask::FLUC_ALICE_IPINFO);
  fJCatalyst[2]->SelectCollisionCandidates(selEvt);
  fJCatalyst[2]->SetSaveAllQA(kTRUE);
  fJCatalyst[2]->SetSaveHMOhist(kFALSE);
  fJCatalyst[2]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
  fJCatalyst[2]->SetInitializeCentralityArray();
  fJCatalyst[2]->SetCentDetName("V0M"); //ADDED + CHECK changed to V0M as in default analysis 
  fJCatalyst[2]->SetTestFilterBit(UsedCut); //ADDED + CHECK changed to tpc only like in default analysis
  fJCatalyst[2]->SetZVertexCut(8.);
  fJCatalyst[2]->SetPtRange(ptMin, 5.0);
  fJCatalyst[2]->SetEtaRange(-0.8, 0.8);
  fJCatalyst[2]->SetNumTPCClusters(70);
  fJCatalyst[2]->SetDCAzCut(2.0);
  fJCatalyst[2]->SetChi2Cuts(0.1, 4.0);
  mgr->AddTask((AliAnalysisTask *)fJCatalyst[2]); //ADDED + CHECK

  //Configuration of the analysis task itself.
  const int SPCCombination = 3;
  TString SPC[SPCCombination] = { "2SPC", "4SPC", "5SPC" };
  AliJSPCTaskRun2 *myTask[Nsets][SPCCombination];
  
  for(Int_t i = 0; i < Nsets; i++){
    if (doSPC == 0) {
      for (Int_t j = 0; j < SPCCombination; j++) {
        myTask[i][j] = new AliJSPCTaskRun2(Form("%s_%s_%s", taskName.Data(), (fJCatalyst[i]->GetJCatalystTaskName()).Data(), SPC[j].Data()));
        myTask[i][j]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
        myTask[i][j]->AliSPCRun2SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
        myTask[i][j]->AliSPCRun2SetSaveAllQA(kTRUE);
        myTask[i][j]->AliSPCRun2SetMinNuPar(14.);
        if(i==1){
          myTask[i][j]->AliSPCRun2SetUseWeights(kTRUE, kTRUE);
        } else {
          myTask[i][j]->AliSPCRun2SetUseWeights(kTRUE, kFALSE);
        }
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
      myTask[i][0] = new AliJSPCTaskRun2(Form("%s_%s_3SPC", taskName.Data(),(fJCatalyst[i]->GetJCatalystTaskName()).Data()));
      myTask[i][0]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
      myTask[i][0]->AliSPCRun2SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
      myTask[i][0]->AliSPCRun2SetSaveAllQA(kTRUE);
      myTask[i][0]->AliSPCRun2SetMinNuPar(14.);
        if(i==1){
          myTask[i][0]->AliSPCRun2SetUseWeights(kTRUE, kTRUE);
        } else {
          myTask[i][0]->AliSPCRun2SetUseWeights(kTRUE, kFALSE);
        }

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
  } //for(Int_t i = 0; i < Nsets; i++)

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[Nsets][SPCCombination+1];

  for (Int_t i = 0; i < Nsets; i++) {
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
  }//for (Int_t i = 0; i < Nsets; i++)

  return myTask[0][0];
}

