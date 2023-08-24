#include <iostream>
#include <sstream>
#include <string>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <vector>
#include <TString.h>

AliAnalysisTask *AddTaskJSPCMasterClosure(TString taskName = "JSPCMaster10h", double ptMin = 0.2, bool newWeightNaming = true, bool IsLHC10h = true)
{

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
        MAPfilenames[i] = Form("%sPhiWeights_LHC17i2a_Error_pt02_s_hybrid.root", MAPdirname.Data());
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

  Int_t UsedCut=0;
  if(IsLHC10h){
    UsedCut = tpconlyCut;
  }
  else {  
    UsedCut = hybridCut;  
  } 

  UInt_t selEvt;
  selEvt = AliVEvent::kAny;// Minimum bias trigger for LHC10h.

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
  fJCatalyst[0]->SetCentDetName("CL1"); //ADDED + CHECK changed to CL1 as in default analysis 
  fJCatalyst[0]->SetTestFilterBit(UsedCut); //ADDED + CHECK changed to tpc only like in default analysis
  fJCatalyst[0]->SetZVertexCut(10.);
  fJCatalyst[0]->SetPtRange(ptMin, 5.0);
  fJCatalyst[0]->SetEtaRange(-0.8, 0.8);
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
  fJCatalyst[1]->SetCentDetName("CL1"); //ADDED + CHECK changed to CL1 as in default analysis 
  fJCatalyst[1]->SetTestFilterBit(UsedCut); //ADDED + CHECK changed to tpc only like in default analysis
  fJCatalyst[1]->SetZVertexCut(10.);
  fJCatalyst[1]->SetPtRange(ptMin, 5.0);
  fJCatalyst[1]->SetEtaRange(-0.8, 0.8);
  fJCatalyst[1]->SetPhiCorrectionIndex(1); //ADDED + CHECK
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
  fJCatalyst[2]->SetCentDetName("CL1"); //ADDED + CHECK changed to CL1 as in default analysis 
  fJCatalyst[2]->SetTestFilterBit(UsedCut); //ADDED + CHECK changed to tpc only like in default analysis
  fJCatalyst[2]->SetZVertexCut(10.);
  fJCatalyst[2]->SetPtRange(ptMin, 5.0);
  fJCatalyst[2]->SetEtaRange(-0.8, 0.8);
  fJCatalyst[2]->SetPhiCorrectionIndex(2); //ADDED + CHECK
  mgr->AddTask((AliAnalysisTask *)fJCatalyst[2]); //ADDED + CHECK

  //Configuration of the analysis task itself.
  const int SPCCombination = 3;
  TString SPC[SPCCombination] = { "2SPC", "3SPC", "4SPC" };
  AliJSPCTask *myTask[Nsets][SPCCombination];
  
  for(Int_t i = 0; i < Nsets; i++){
      for(Int_t j = 0; j < SPCCombination; j++){

        std::cout<<Form("%s",(fJCatalyst[i]->GetJCatalystTaskName()).Data())<<std::endl;

        myTask[i][j] = new AliJSPCTask(Form("%s_%s_%s", taskName.Data(),(fJCatalyst[i]->GetJCatalystTaskName()).Data(), SPC[j].Data()));
      	myTask[i][j]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());

      	myTask[i][j]->JSPCSetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
      	myTask[i][j]->JSPCSetSaveAllQA(kTRUE);
      	myTask[i][j]->JSPCSetMinNuPar(14.);
      	myTask[i][j]->JSPCSetFisherYates(kFALSE, 1.); 

        if(i==0){ myTask[i][j]->JSPCSetUseWeights(kFALSE, kFALSE); } 
      	else if(i==1){ 
          if(IsLHC10h){
             myTask[i][j]->JSPCSetUseWeights(kTRUE, kFALSE); 
          }
          else {  
             myTask[i][j]->JSPCSetUseWeights(kTRUE, kTRUE); 
         }     
        } 
        else if(i==2){ myTask[i][j]->JSPCSetUseWeights(kFALSE, kFALSE); } 

        myTask[i][j]->JSPCSetEtaGaps(kFALSE, 0.8);

      	if(j==0){
      	  myTask[i][j]->JSPCSetCorrSet1(4., 6.,-2.,-2.,-2., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet2(3., 6.,-3.,-3., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet3(3., 4.,-2.,-2., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet4(3., 8.,-4.,-4., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet5(5., 3., 3.,-2.,-2.,-2., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet6(5., 8.,-2.,-2.,-2.,-2., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet7(0., 0., 0., 0., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
      	}

      	if(j==1){
      	  myTask[i][j]->JSPCSetCorrSet1(4., 3.,-4.,-4., 5., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet2(3., 2., 4.,-6., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet3(3., 2., 3.,-5., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet4(4., 2.,-3.,-3., 4., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet5(5., 2., 3., 3.,-4.,-4., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet6(6., 2., 2., 2., 2.,-3.,-5.,0.);
      	  myTask[i][j]->JSPCSetCorrSet7(5., 3., 3., 3.,-4.,-5., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
      	}

      	if(j==2){
      	  myTask[i][j]->JSPCSetCorrSet1(4., 2.,-3.,-4., 5., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet2(6., 2., 2., 2., 3.,-4.,-5.,0.);
      	  myTask[i][j]->JSPCSetCorrSet3(5., 2., 2.,-3., 4.,-5., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet4(6., 2., 2., 3., 3.,-4.,-6.,0.);
      	  myTask[i][j]->JSPCSetCorrSet5(0., 0., 0., 0., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet6(0., 0., 0., 0., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet7(0., 0., 0., 0., 0., 0., 0.,0.);
      	  myTask[i][j]->JSPCSetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
      	}

        myTask[i][j]->JSPCSetMixed(kFALSE,2., kFALSE, kTRUE);

        mgr->AddTask((AliAnalysisTask *) myTask[i][j]);

      }//for(Int_t j = 0; j < SPCCombination; j++)
  }//for(Int_t i = 0; i < Nsets; i++)

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[Nsets][SPCCombination+1];

  for (Int_t i = 0; i < Nsets; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);

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

  }//for (Int_t i = 0; i < Nsets; i++)

  return myTask[0][0];
}

