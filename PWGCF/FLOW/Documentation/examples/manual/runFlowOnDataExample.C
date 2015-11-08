void runFlowOnDataExample() {
    // example macro of running a flow analysis on local 
    // data
    //
    // all steps are explained in detail in
    // chapter 3.4.3 of the flow package manual
    // $ALICE_PHYSICS/PWGCF/FLOW/Documentation/FlowPackageManual.pdf
    // which is also available on the twiki page
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FlowPAGFlowPackageManual
    //
    // author: Redmer Alexander Bertens, Utrecht University
    // rbertens@cern.ch , rbertens@nikhef.nl , r.a.bertens@uu.nl

  // load libraries
  gSystem->Load("libPWGflowTasks");

    // create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("MyManager");

  // create a tchain which will point to an aod tree
  TChain* chain = new TChain("aodTree");
  // add a few files to the chain (change this so that your local files are added)
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0003/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0003/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0004/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0005/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0006/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0007/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0008/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0009/AliAOD.root");
  chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0010/AliAOD.root");
  // create an input handler
  AliVEventHandler* inputH = new AliAODInputHandler();
  // and connect it to the manager
  mgr->SetInputEventHandler(inputH);

   // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      printf("No analysis manager to connect to!\n");
      return NULL;
  }
        
  // just to see if all went well, check if the input event handler has been connected
  if (!mgr->GetInputEventHandler()) {
      printf("This task requires an input event handler!\n");
      return NULL;
    }

  // create instance of the class. because possible qa plots are added in a second ouptut slot,
  // the flow analysis task must know if you want to save qa plots at the time of class construction
  Bool_t doQA = kTRUE;
  // craete instance of the class
  AliAnalysisTaskFlowEvent* taskFE = new AliAnalysisTaskFlowEvent("FlowEventTask", "", doQA);
  // add the task to the manager
  mgr->AddTask(taskFE);
  // set the trigger selection
  taskFE->SelectCollisionCandidates(AliVEvent::kMB);

    // define the event cuts object
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("EventCuts");
  // configure some event cuts, starting with centrality
  cutsEvent->SetCentralityPercentileRange(20., 30.);
  // method used for centrality determination
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  // vertex-z cut
  cutsEvent->SetPrimaryVertexZrange(-10.,10.);
  // enable the qa plots
  cutsEvent->SetQA(doQA);
  // explicit multiplicity outlier cut
  cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
  cutsEvent->SetLHC10h(kTRUE);
  
  
  // and, last but not least, pass these cuts to your flow event task
  taskFE->SetCutsEvent(cutsEvent);

  //create the track cuts object using a static function of AliFlowTrackCuts
  AliFlowTrackCuts* cutsRP = AliFlowTrackCuts::GetAODTrackCutsForFilterBit(1, "RP_cuts");
  // specify the pt range
  cutsRP->SetPtRange(0.2, 5.);
  // specify eta range
  cutsRP->SetEtaRange(-0.8, 0.8);
  // specify track type
  cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  // enable saving qa histograms
  cutsRP->SetQA(kTRUE);

    //create the track cuts object using a static function of AliFlowTrackCuts
  AliFlowTrackCuts* cutsPOI = AliFlowTrackCuts::GetAODTrackCutsForFilterBit(1, "pion selection");
  // specify the pt range
  cutsPOI->SetPtRange(0.2, 5.);
  // specify eta range
  cutsPOI->SetEtaRange(-0.8, 0.8);
  // specify the track type
  cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  // enable saving qa histograms
  cutsPOI->SetQA(kTRUE);

    // which particle do we want to identify ?
  AliPID::EParticleType particleType=AliPID::kPion;
  // specify the pid method that we want to use  
  AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian;
  // define the probability (between 0 and 1) 
  Double_t probability = .9;
  // pass these variables to the track cut object
  cutsPOI->SetPID(particleType, sourcePID, probability);
  // the bayesian pid routine uses priors tuned to an average centrality
  cutsPOI->SetPriors(35.);

    // connect the RP's to the flow event task
  taskFE->SetCutsRP(cutsRP);
  // connect the POI's to the flow event task
  taskFE->SetCutsPOI(cutsPOI);

    // get the default name of the output file ("AnalysisResults.root")
  TString file = AliAnalysisManager::GetCommonFileName();
  // get the common input container from the analysis manager
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   // create a data container for the output of the flow event task  
   // the output of the task is the AliFlowEventSimle class which will
   // be passed to the flow analysis tasks. note that we use a kExchangeContainer here,
   // which exchanges data between classes of the analysis chain, but is not
   // written to the output file
  AliAnalysisDataContainer *coutputFE = mgr->CreateContainer(
      "FlowEventContainer",
      AliFlowEventSimple::Class(),
      AliAnalysisManager::kExchangeContainer);
  // connect the input data to the flow event task
  mgr->ConnectInput(taskFE,0,cinput);
  // and connect the output to the flow event task
  mgr->ConnectOutput(taskFE,1,coutputFE);
  // create an additional container for the QA output of the flow event task
  // the QA histograms will be stored in a sub-folder of the output file called 'QA'
  TString taskFEQAname = file;
  taskFEQAname += ":QA";
  AliAnalysisDataContainer* coutputFEQA = mgr->CreateContainer(
      "FlowEventContainerQA",
       TList::Class(),
       AliAnalysisManager::kOutputContainer,
       taskFEQAname.Data()       
       );
  // and connect the qa output container to the flow event. 
  // this container will be written to the output file
  mgr->ConnectOutput(taskFE,2,coutputFEQA);

    // declare necessary pointers
  AliAnalysisDataContainer *coutputQC[3];
  AliAnalysisTaskQCumulants *taskQC[3];

  // the tasks will be creaated and added to the manager in a loop
  for(Int_t i = 0; i < 3; i++) {
      // create the flow analysis tasks
      taskQC[i] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_n=%i", i+2));
      // set the triggers 
      taskQC[i]->SelectCollisionCandidates(AliVEvent::kMB);
      // and set the correct harmonic n
      taskQC[i]->SetHarmonic(i+2);

      // connect the task to the analysis manager
      mgr->AddTask(taskQC[i]);

      // create and connect the output containers
      TString outputQC = file;
      // create a sub-folder in the output file for each flow analysis task's output
      outputQC += Form(":QC_output_for_n=%i", i+2);
      /// create the output containers
      coutputQC[i] = mgr->CreateContainer(
          outputQC.Data(),
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          outputQC);
      // connect the output of the flow event task to the flow analysis task
      mgr->ConnectInput(taskQC[i], 0, coutputFE);
      // and connect the output of the flow analysis task to the output container
      // which will be written to the output file
      mgr->ConnectOutput(taskQC[i], 1, coutputQC[i]);
  }
  
  // check if we can initialize the manager
  if(!mgr->InitAnalysis()) return;   
  // print the status of the manager to screen 
  mgr->PrintStatus();
  // print to screen how the analysis is progressing
  mgr->SetUseProgressBar(1, 25);
  // start the analysis locally, reading the events from the tchain
  mgr->StartAnalysis("local", chain);

  
}
