/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTask* macro for flow analysis
// Creates a Flow Event task and adds it to the analysis manager.
// Also creates Flow Analysis tasks and connects them to the output of the flow event task.
//
/////////////////////////////////////////////////////////////////////////////////////////////

AliAnalysisTaskFlowEvent *AddTaskAMPTTest(TString fileNameBase="AnalysisResults",
					  Double_t etaMin = -.9,
                                          Double_t etaMax = .9,
					  Double_t ptMin  = 0.05,
					  Double_t ptMax  = 20.0,
					  Int_t chargePOI = 0,
					  Int_t harmonic  = 2,
                                          Int_t pdg_pid   = 0,
					  Double_t gImpactParameterMin = 0.0,
					  Double_t gImpactParameterMax = 100.0,
					  Int_t gRefMultMin = 0,
					  Int_t gRefMultMax = 100,
					  Bool_t MCEP     = kFALSE,
					  Bool_t SP       = kTRUE,
					  Bool_t GFC      = kFALSE,
					  Bool_t QC       = kTRUE,
					  Bool_t FQD      = kFALSE,
					  Bool_t LYZ1SUM  = kFALSE,
					  Bool_t LYZ1PROD = kFALSE,
					  Bool_t LYZ2SUM  = kFALSE,
					  Bool_t LYZ2PROD = kFALSE,
					  Bool_t LYZEP    = kFALSE,
					  Bool_t MH       = kFALSE,
					  Bool_t NL       = kFALSE,
					  Int_t side      = 0) {	  
  // Define the range for eta subevents (for SP method)
  Double_t minA = etaMin;
  Double_t maxA = 0.0;
  Double_t minB = 0.0;
  Double_t maxB = etaMax;

  // AFTERBURNER
  Bool_t useAfterBurner=kFALSE;
  Double_t v1=0.0;
  Double_t v2=0.0;
  Double_t v3=0.0;
  Double_t v4=0.0;
  Int_t numberOfTrackClones=0; //non-flow

  // Define a range of the detector to exclude
  Bool_t ExcludeRegion = kFALSE;
  Double_t excludeEtaMin = -0.;
  Double_t excludeEtaMax = 0.;
  Double_t excludePhiMin = 0.;
  Double_t excludePhiMax = 0.;

  // use physics selection class
  Bool_t  UsePhysicsSelection = kFALSE;

  // QA
  Bool_t runQAtask=kFALSE;
  Bool_t FillQAntuple=kFALSE;
  Bool_t DoQAcorrelations=kFALSE;

  //Methods
  Bool_t METHODS[] = {SP,QC};

  // Boolean to use/not use weights for the Q vector
  Bool_t WEIGHTS[] = {kFALSE,kFALSE,kFALSE}; //Phi, v'(pt), v'(eta)

  // SETTING THE CUTS
  //---------Data selection----------
  //kMC, kGlobal, kTPCstandalone, kESD_SPDtracklet
  AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kMC;
  AliFlowTrackCuts::trackParameterType poitype = AliFlowTrackCuts::kMC;

  //---------Parameter mixing--------
  //kPure - no mixing, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt
  AliFlowTrackCuts::trackParameterMix rpmix = AliFlowTrackCuts::kTrackWithMCkine;
  AliFlowTrackCuts::trackParameterMix poimix = AliFlowTrackCuts::kTrackWithMCkine;


  const char* rptypestr = AliFlowTrackCuts::GetParamTypeName(rptype);
  const char* poitypestr = AliFlowTrackCuts::GetParamTypeName(poitype);

  //===========================================================================
  // EVENTS CUTS:
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("event cuts");
  
  // Ref mult TRACK CUTS:
  AliFlowTrackCuts* cutsRefMult = new AliFlowTrackCuts("MCRefMult");
  cutsRefMult->SetParamType(rptype);
  cutsRefMult->SetParamMix(rpmix);
  cutsRefMult->SetPtRange(ptMin,ptMax);
  cutsRefMult->SetEtaRange(etaMin,etaMax);
  cutsRefMult->SetRequireCharge(kTRUE);
  cutsRefMult->SetQA(kFALSE);

  // RP TRACK CUTS:
  AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("MCRP");
  cutsRP->SetParamType(rptype);
  cutsRP->SetParamMix(rpmix);
  cutsRP->SetPtRange(ptMin,ptMax);
  cutsRP->SetEtaRange(etaMin,etaMax);
  cutsRP->SetQA(kFALSE);

  // POI TRACK CUTS:
  AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("MCPOI");
  cutsPOI->SetParamType(poitype);
  cutsPOI->SetParamMix(poimix);
  cutsPOI->SetPtRange(ptMin,ptMax);
  if(pdg_pid != 0) {
      cutsPOI->SetMCPID(pdg_pid);
      cutsPOI->SetIgnoreSignInMCPID(kTRUE);
  }

  // side A
  if(side < 0)
    cutsPOI->SetEtaRange(etaMin,-0.0);

  // side C
  else if(side > 0)
  cutsPOI->SetEtaRange(0.0,etaMax);

  // both sides
  else
    cutsPOI->SetEtaRange(etaMin,etaMax);
  if(chargePOI != 0)
    cutsPOI->SetCharge(chargePOI);
  cutsPOI->SetQA(kFALSE);

  TString outputSlotName("");
  outputSlotName+=Form("V%i_",harmonic);
  outputSlotName+=cutsRP->GetName();
  outputSlotName+="_";
  outputSlotName+=cutsPOI->GetName();
  outputSlotName+="_";
  outputSlotName+=Form("PDG_%i", pdg_pid);


  TString fileName(fileNameBase);
  fileName.Append(".root");

  Bool_t useWeights  = WEIGHTS[0] || WEIGHTS[1] || WEIGHTS[2];
  if (useWeights) cout<<"Weights are used"<<endl;
  else cout<<"Weights are not used"<<endl;
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowEvent", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis
  // manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlowEvent", "This task requires an input event handler");
    return NULL;
  }  

  // Open external input files
  //===========================================================================
  //weights: 
  TFile *weightsFile = NULL;
  TList *weightsList = NULL;

  if(useWeights) {
    //open the file with the weights:
    weightsFile = TFile::Open("weights.root","READ");
    if(weightsFile) {
      //access the list which holds the histos with weigths:
      weightsList = (TList*)weightsFile->Get("weights");
    }
    else {
      cout<<" WARNING: the file <weights.root> with weights from the previous run was not available."<<endl;
      return NULL;
    } 
  }
  
  //LYZ2
  if (LYZ2SUM || LYZ2PROD) {
    //read the outputfile of the first run
    TString outputFileName = "AnalysisResults1.root";
    TString pwd(gSystem->pwd());
    pwd+="/";
    pwd+=outputFileName.Data();
    TFile *outputFile = NULL;
    if(gSystem->AccessPathName(pwd.Data(),kFileExists)) {
      cout<<"WARNING: You do not have an output file:"<<endl;
      cout<<"         "<<pwd.Data()<<endl;
      exit(0);
    } else { outputFile = TFile::Open(pwd.Data(),"READ");}
    
    if (LYZ2SUM){  
      // read the output directory from LYZ1SUM 
      TString inputFileNameLYZ2SUM = "outputLYZ1SUManalysis" ;
      inputFileNameLYZ2SUM += rptypestr;
      cout<<"The input directory is "<<inputFileNameLYZ2SUM.Data()<<endl;
      TFile* fInputFileLYZ2SUM = (TFile*)outputFile->FindObjectAny(inputFileNameLYZ2SUM.Data());
      if(!fInputFileLYZ2SUM || fInputFileLYZ2SUM->IsZombie()) { 
	cerr << " ERROR: To run LYZ2SUM you need the output file from LYZ1SUM. This file is not there! Please run LYZ1SUM first." << endl ; 
	return NULL;
      }
      else {
	TList* fInputListLYZ2SUM = (TList*)fInputFileLYZ2SUM->Get("LYZ1SUM");
	if (!fInputListLYZ2SUM) {cout<<"list is NULL pointer!"<<endl;}
      }
      cout<<"LYZ2SUM input file/list read..."<<endl;
    }

    if (LYZ2PROD){  
      // read the output directory from LYZ1PROD 
      TString inputFileNameLYZ2PROD = "outputLYZ1PRODanalysis" ;
      inputFileNameLYZ2PROD += rptypestr;
      cout<<"The input directory is "<<inputFileNameLYZ2PROD.Data()<<endl;
      TFile* fInputFileLYZ2PROD = (TFile*)outputFile->FindObjectAny(inputFileNameLYZ2PROD.Data());
      if(!fInputFileLYZ2PROD || fInputFileLYZ2PROD->IsZombie()) { 
	cerr << " ERROR: To run LYZ2PROD you need the output file from LYZ1PROD. This file is not there! Please run LYZ1PROD first." << endl ; 
	return NULL;
      }
      else {
	TList* fInputListLYZ2PROD = (TList*)fInputFileLYZ2PROD->Get("LYZ1PROD");
	if (!fInputListLYZ2PROD) {cout<<"list is NULL pointer!"<<endl;}
      }
      cout<<"LYZ2PROD input file/list read..."<<endl;
    }
  }

  if (LYZEP) {
    //read the outputfile of the second run
    TString outputFileName = "AnalysisResults2.root";
    TString pwd(gSystem->pwd());
    pwd+="/";
    pwd+=outputFileName.Data();
    TFile *outputFile = NULL;
    if(gSystem->AccessPathName(pwd.Data(),kFileExists)) {
      cout<<"WARNING: You do not have an output file:"<<endl;
      cout<<"         "<<pwd.Data()<<endl;
      exit(0);
    } else {
      outputFile = TFile::Open(pwd.Data(),"READ");
    }
    
    // read the output file from LYZ2SUM
    TString inputFileNameLYZEP = "outputLYZ2SUManalysis" ;
    inputFileNameLYZEP += rptypestr;
    cout<<"The input file is "<<inputFileNameLYZEP.Data()<<endl;
    TFile* fInputFileLYZEP = (TFile*)outputFile->FindObjectAny(inputFileNameLYZEP.Data());
    if(!fInputFileLYZEP || fInputFileLYZEP->IsZombie()) { 
      cerr << " ERROR: To run LYZEP you need the output file from LYZ2SUM. This file is not there! Please run LYZ2SUM first." << endl ; 
      return NULL;
    }
    else {
      TList* fInputListLYZEP = (TList*)fInputFileLYZEP->Get("LYZ2SUM");
      if (!fInputListLYZEP) {cout<<"list is NULL pointer!"<<endl;}
    }
    cout<<"LYZEP input file/list read..."<<endl;
  }
  
  
  //===========================================================================
  AliAnalysisTaskFlowEvent *taskFE = NULL;

  if(useAfterBurner)
    { 
      taskFE = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent %s",outputSlotName.Data()),"",kFALSE,1);
      taskFE->SetFlow(v1,v2,v3,v4); 
      taskFE->SetNonFlowNumberOfTrackClones(numberOfTrackClones);
      taskFE->SetAfterburnerOn();
    }
  else {taskFE = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent %s",outputSlotName.Data()),"",kFALSE); }
  if (ExcludeRegion) {
    taskFE->DefineDeadZone(excludeEtaMin, excludeEtaMax, excludePhiMin, excludePhiMax); 
  }
  taskFE->SetSubeventEtaRange(minA, maxA, minB, maxB);
  if (UsePhysicsSelection) {
    taskFE->SelectCollisionCandidates(AliVEvent::kMB);
    cout<<"Using Physics Selection"<<endl;
  }
  mgr->AddTask(taskFE);
  
  // Pass cuts for RPs and POIs to the task:
  taskFE->SetCutsEvent(cutsEvent);
  taskFE->SetCutsRP(cutsRP);
  taskFE->SetCutsPOI(cutsPOI);

  // Create the analysis tasks, add them to the manager.
  //===========================================================================
    AliAnalysisTaskScalarProduct *taskSP;
    AliAnalysisTaskQCumulants *taskQC;

  if (SP){
    taskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct %s",outputSlotName.Data()),WEIGHTS[0]);
    taskSP->SetRelDiffMsub(1.0);
    taskSP->SetApplyCorrectionForNUA(kTRUE);
    taskSP->SetHarmonic(harmonic);
    if(side < 0)   taskSP->SetTotalQvector("Qb");
    else     taskSP->SetTotalQvector("Qa");
    mgr->AddTask(taskSP);
  }
  if (QC){
    taskQC = new AliAnalysisTaskQCumulants(Form("TaskQCumulants %s",outputSlotName.Data()),useWeights);
    //taskQC->SetMultiplicityIs(AliFlowCommonConstants::kRP); 
    taskQC->SetMultiplicityIs(AliFlowCommonConstants::kExternal); 
    //taskQC->SetFillProfilesVsMUsingWeights(kFALSE);
    //taskQC->SetUseQvectorTerms(kTRUE);
    taskQC->SetUsePhiWeights(WEIGHTS[0]); 
    taskQC->SetUsePtWeights(WEIGHTS[1]);
    taskQC->SetUseEtaWeights(WEIGHTS[2]); 
    taskQC->SetCalculateCumulantsVsM(kTRUE);
    taskQC->SetnBinsMult(10000);
    taskQC->SetMinMult(0.);
    taskQC->SetMaxMult(10000.);
    taskQC->SetHarmonic(harmonic);
    taskQC->SetApplyCorrectionForNUA(kFALSE);
    taskQC->SetFillMultipleControlHistograms(kFALSE);     
    mgr->AddTask(taskQC);
  }

  // Create the output container for the data produced by the task
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutputFE = 
  mgr->CreateContainer(Form("FlowEventSimple %s",outputSlotName.Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput(taskFE,0,cinput1); 
  mgr->ConnectOutput(taskFE,1,coutputFE);
 
  if (taskFE->GetQAOn())
  {
    TString outputQA = fileName;
    outputQA += ":QA";
    AliAnalysisDataContainer* coutputFEQA = 
    mgr->CreateContainer(Form("QA %s",outputSlotName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA);
    mgr->ConnectOutput(taskFE,2,coutputFEQA);
  }

  // Create the output containers for the data produced by the analysis tasks
  // Connect to the input and output containers
  //===========================================================================
    AliAnalysisDataContainer *cinputWeights;
   if( useWeights) cinputWeights = mgr->CreateContainer(Form("Weights %s",outputSlotName.Data()),
								   TList::Class(),AliAnalysisManager::kInputContainer); 

  if(SP) {
    TString outputSP = fileName;
    outputSP += ":outputSPanalysis";
    outputSP+= rptypestr;
    AliAnalysisDataContainer *coutputSP = mgr->CreateContainer(Form("SP%s",outputSlotName.Data()), 
							       TList::Class(),AliAnalysisManager::kOutputContainer,outputSP); 
    mgr->ConnectInput(taskSP,0,coutputFE); 
    mgr->ConnectOutput(taskSP,1,coutputSP); 
    if (WEIGHTS[0]) {
      mgr->ConnectInput(taskSP,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    }
  }
  if(QC) {
    TString outputQC = fileName;
    outputQC += ":outputQCanalysis";
    outputQC+= rptypestr;

    AliAnalysisDataContainer *coutputQC = mgr->CreateContainer(Form("QC%s",outputSlotName.Data()), 
							       TList::Class(),AliAnalysisManager::kOutputContainer,outputQC); 
    mgr->ConnectInput(taskQC,0,coutputFE); 
    mgr->ConnectOutput(taskQC,1,coutputQC);
    if (useWeights) {
      mgr->ConnectInput(taskQC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  if (runQAtask)
  {
    AliAnalysisTaskQAflow* taskQAflow = new AliAnalysisTaskQAflow(Form("TaskQAflow %s",outputSlotName.Data()));
    taskQAflow->SetEventCuts(cutsEvent);
    taskQAflow->SetTrackCuts(cutsRP);
    taskQAflow->SetFillNTuple(FillQAntuple);
    taskQAflow->SetDoCorrelations(DoQAcorrelations);
    mgr->AddTask(taskQAflow);
    
    Printf("outputSlotName %s",outputSlotName.Data());
    TString taskQAoutputFileName(fileNameBase);
    taskQAoutputFileName.Append("_QA.root");
    AliAnalysisDataContainer* coutputQAtask = mgr->CreateContainer(Form("flowQA %s",outputSlotName.Data()),
                                              TObjArray::Class(),
                                              AliAnalysisManager::kOutputContainer,
                                              taskQAoutputFileName);
    AliAnalysisDataContainer* coutputQAtaskTree = mgr->CreateContainer(Form("flowQAntuple %s",outputSlotName.Data()),
                                              TNtuple::Class(),
                                              AliAnalysisManager::kOutputContainer,
                                              taskQAoutputFileName);
    mgr->ConnectInput(taskQAflow,0,mgr->GetCommonInputContainer());
    mgr->ConnectInput(taskQAflow,1,coutputFE);
    mgr->ConnectOutput(taskQAflow,1,coutputQAtask);
    if (FillQAntuple) mgr->ConnectOutput(taskQAflow,2,coutputQAtaskTree);
  }

  return taskFE;
}






