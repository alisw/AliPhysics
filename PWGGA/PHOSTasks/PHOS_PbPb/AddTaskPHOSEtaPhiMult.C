class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;

AliAnalysisTaskEtaPhigg* AddTaskPHOSEtaPhiMult(const char* name = "PHOSEtaPhiMult",
					   UInt_t offlineTriggerMask = AliVEvent::kCentral,
                               Int_t trainConfig = 1,  //change different set of cuts
                               //Bool_t isMC   = kFALSE, //run MC
                               Int_t enableQAMesonTask = 0, //enable QA in AddTask_GammaConvFlow_PbPb2
                               Int_t enableQAPhotonTask = 0, // enable additional QA task
                               //TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                               Bool_t doWeighting = kFALSE,  //enable Weighting
                               TString cutnumberAODBranch = "1000000060084000001500000",
                               Bool_t BasicHistoSP = kFALSE,
                               Bool_t debug = kFALSE,
                               Bool_t UseMassSel = kFALSE,
                               Float_t MinMass = 0,
                               Float_t MaxMass = 0.2,
                               Bool_t UseKappaSel = kFALSE,
                               Float_t MinKappa = 0,
                               Float_t MaxKappa = 10
                               )
{
  //Add a task AliAnalysisTaskGammaFlow to the analysis train
  //Author: Dmitri Peressounko
  /* $Id$ */

 
  // ================= Load Librariers =================================
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGGAGammaConv");
  gSystem->Load("libPWGGAPHOSTasks");
  
  Int_t isHeavyIon = 1;


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskEtaPhigg", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisTaskEtaPhigg", "This task requires an input event handler");
    return NULL;
  }
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

//========= Add PID Reponse to ANALYSIS manager ====
   if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AddTaskPIDResponse(kFALSE);
   }//========= Add PID Reponse to ANALYSIS manager ====
  

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
  
    //=========  Set Cutnumber for V0Reader ================================
    TString cutnumberPhoton = "00000008400100001500000000";
    TString cutnumberEvent = "10000003";
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");

      
      fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
      fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
      fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
      fV0ReaderV1->SetUseMassToZero(kFALSE);
      if (!mgr) {
          Error("AddTask_V0ReaderV1", "No analysis manager found.");
          return;
      }
      
      AliConvEventCuts *fEventCuts=NULL;
      if(cutnumberEvent!=""){
          fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
          fEventCuts->SetPreSelectionCutFlag(kTRUE);
          if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
              fV0ReaderV1->SetEventCuts(fEventCuts);
              fEventCuts->SetFillCutHistograms("",kTRUE);
          }
      }
      
      // Set AnalysisCut Number
      AliConversionPhotonCuts *fCuts=NULL;
      if(cutnumberPhoton!=""){
          fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
          fCuts->SetPreSelectionCutFlag(kTRUE);
          fCuts->SetIsHeavyIon(isHeavyIon);
          if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
              fV0ReaderV1->SetConversionCuts(fCuts);
              fCuts->SetFillCutHistograms("",kTRUE);
          }
      }
      
      if(inputHandler->IsA()==AliAODInputHandler::Class()){
          // AOD mode
          fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
      }
      fV0ReaderV1->Init();
      
      AliLog::SetGlobalLogLevel(AliLog::kInfo);
      
      //connect input V0Reader
      mgr->AddTask(fV0ReaderV1);
      mgr->ConnectInput(fV0ReaderV1,0,cinput);        
  }
    

  
  AliAnalysisTaskEtaPhiMultgg * task = new AliAnalysisTaskEtaPhiMultgg(Form("%sEtaPhiMult", name));
  task->SelectCollisionCandidates(offlineTriggerMask);

  if(UseMassSel==kTRUE)  task->SetMassWindow(MinMass,MaxMass);
  if(UseKappaSel==kTRUE) task->SetKappaWindow(MinKappa,MaxKappa);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  TString trigger;
  if(offlineTriggerMask==AliVEvent::kCentral) trigger = "Central" ;
  else if(offlineTriggerMask==AliVEvent::kSemiCentral) trigger = "SemiCentral" ;
  else trigger = "Other" ;
  TString cname(Form("%s%s", name,trigger.Data()));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 1;
//     if(trainConfig == 37){
// 		numberOfCuts = 6;
// 		cout << "number of cuts: " << numberOfCuts << endl;
// 	}
  TString *eventCutArray = new TString[numberOfCuts];
  TString *photonCutArray = new TString[numberOfCuts];
  TString *mesonCutArray = new TString[numberOfCuts];

  if (trainConfig == 1){
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";
  } else if (trainConfig == 2) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";
  } else if (trainConfig == 3) {
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";
  } else if (trainConfig == 4) {
    eventCutArray[ 0] = "50200013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";
  } else if (trainConfig == 5) {
    eventCutArray[ 0] = "51200013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";
  } else if (trainConfig == 6) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000";
  } else if (trainConfig == 7) {
    eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 8) {
    eventCutArray[ 0] = "54800013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 9) {
    eventCutArray[ 0] = "54500013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 10) {
    eventCutArray[ 0] = "55600013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 11) {
    eventCutArray[ 0] = "56800013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 12) {
    eventCutArray[ 0] = "56700013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 13) {
    eventCutArray[ 0] = "57800013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 14) {
    eventCutArray[ 0] = "46900013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else if (trainConfig == 15) {
    eventCutArray[ 0] = "58900013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152206500900000";
  } else  if (trainConfig == 16){
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 17) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 18) {
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 19) {
    eventCutArray[ 0] = "50200013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 20) {
    eventCutArray[ 0] = "51200013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 21) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 22) {
    eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 23) {
    eventCutArray[ 0] = "54800013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 24) {
    eventCutArray[ 0] = "54500013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 25) {
    eventCutArray[ 0] = "55600013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 26) {
    eventCutArray[ 0] = "56800013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 27) {
    eventCutArray[ 0] = "56700013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 28) {
    eventCutArray[ 0] = "57800013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 29) {
    eventCutArray[ 0] = "46900013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 30) {
    eventCutArray[ 0] = "58900013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 31) {
    eventCutArray[ 0] = "50800013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 32) {
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 33) {
    eventCutArray[ 0] = "53500013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 34) {
    eventCutArray[ 0] = "54500013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 35) {
    eventCutArray[ 0] = "53400013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 36) {
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 37){                           
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 38) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 39) {
    eventCutArray[ 0] = "51200013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 40) {
    eventCutArray[ 0] = "54600013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 41) {
    eventCutArray[ 0] = "56800013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 42) {
    eventCutArray[ 0] = "53400013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 43) {
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 44) {
    eventCutArray[ 0] = "50200013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 45) {
    eventCutArray[ 0] = "50400013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 46) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "00200009697004000500000000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 47) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "00200009697005000500000000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 48) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "00200009797004000500000000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 49) {
    eventCutArray[ 0] = "61200013"; photonCutArray[ 0] = "00200009797005000500000000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 50) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 51) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 52) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 53) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 54) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 55) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 56) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 57) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009007000008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else if (trainConfig == 58) {
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
  } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return;
  }

      
  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();
  
  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("BOX");
  HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  
  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
//    task->SetMesonCutList(numberOfCuts,MesonCutList);
//    task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kFALSE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA


  return task;
}
