///////////////////////////////////////////////////////////////////
//                                                               //
// AddTaskFlowModes                                              //
// Author: Naghmeh Mohammadi (nmohamma@cern.ch), Nikhef, 2017    //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskFlowModes* AddTaskFlowModes(TString name = "FlowHarmonics", 
					   AliAnalysisTaskFlowModes::ColSystem collisionSystem = AliAnalysisTaskFlowModes::kPbPb,
					   Int_t PVtxZMax = 10,
					   Int_t NumTPCclsMin = 70,
					   Int_t TrackFilterBit = 96,
                       			   Double_t MaxChi2perTPC = 4,
					   Double_t DCAzMax = 0.,
                       			   Double_t DCAxyMax = 0.,
					   TString MultEstimator = "V0M",
					   Bool_t PileUp=kFALSE,
					   Bool_t DoOnlyMixedFlow = kTRUE,
	                                   Bool_t FillWeights = kFALSE,
					   Bool_t AntiProtonOnly = kFALSE,
					   Int_t PIDComb = 2,					   
					   Bool_t PIDbayesian = kFALSE,
                       			   Double_t PIDprob = 0.9)
{
  cout<<"=========================================================="<<endl;
  cout<<"====================AddTaskFlowModes======================"<<endl;
  cout<<"=========================================================="<<endl;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;

  TString fileName = AliAnalysisManager::GetCommonFileName();   // by default, a file is open for writing. here, we get the filename
  fileName += Form(":%s",name.Data());      // create a subfolder in the file

  AliAnalysisTaskFlowModes* task1 = new AliAnalysisTaskFlowModes(name.Data()); // now we create an instance of your task
  if(!task1) return 0x0;

 // Analysis
    task1->SetRunMode(AliAnalysisTaskFlowModes::kFull);//kFull
    task1->SetColisionSystem(collisionSystem);//AliAnalysisTaskFlowModes::kPbPb);// collisionSystem  kPP, kPbPb
    task1->SetNumEventsAnalyse(50);//In case of fRunMode == kTest it only analyses up to 50 events.  
    task1->SetSampling(kTRUE);
    task1->SetAnalysisType(AliAnalysisTaskFlowModes::kAOD);
    task1->SetFillQAhistos(kTRUE);
    task1->SetProcessCharged(kTRUE);
    // PID selection
    task1->SetProcessPID(kTRUE,PIDbayesian);
    task1->SetPIDnsigmaCombination(PIDComb); // applies PIDnsigma combination 2 out of the 3 possible combinations.
    task1->SetPIDUseAntiProtonOnly(AntiProtonOnly);
    task1->SetPIDNumSigmasPionMax(3);
    task1->SetPIDNumSigmasKaonMax(3);
    task1->SetPIDNumSigmasProtonMax(3);
    task1->SetPIDNumSigmasCombinedNoTOFrejection(kTRUE);
    if(PIDbayesian) task1->SetBayesianProbability(PIDprob);
    


    // Flow
    task1->SetFlowRFPsPtMin(0.2);
    task1->SetFlowRFPsPtMax(5.);
    // task1->SetFlowDoFourCorrelations(kFALSE);
    task1->SetFlowDoOnlyMixedCorrelations(DoOnlyMixedFlow);
    task1->SetFlowFillWeights(FillWeights);
    //task1->SetUseNUAWeigthsFile("alice/cern.ch/user/n/nmohamma/CorrectionMaps/hadronPID/PIDnSigma/fb96/NUACorrectionMap.root"); 
    //task1->SetUseNUEWeigthsFile("alice/cern.ch/user/n/nmohamma/CorrectionMaps/fb768/NUECorrectionMap_fb768_nSigmaComb2.root"); 

    //task1->SetPositivelyChargedRef(kFALSE);//reference particles both positively charged
    //task1->SetNegativelyChargedRef(kFALSE);//reference particles both negatively charged
    //task1->SetPositivelyChargedPOI(kFALSE);//POIs positively charged
    //task1->SetNegativelyChargedPOI(kFALSE);//POIs negatively charged

    // Events selection
    task1->SetTrigger(0); // case 0: kINT7, case 1: kHighMultV0, case 2: kHighMultSPD.
    task1->SetMultEstimator(MultEstimator);//CL1, CL0
    task1->SetCentralityRange(kTRUE); // kTRUE: run from 0-100% kFALSE: run from 50-100%
    task1->SetPVtxZMax(PVtxZMax);
    // Charged selection
    task1->SetChargedEtaMax(0.8);
    // task1->SetChargedPtMin(0.2);
    // task1->SetChargedPtMax(5.);
    // In PID vn paper: DCAz < 3.2 cm and DCAxy < 2.4 cm
    // if DCAxy and DCAz set to 0. then the DCA cuts only come from the filterbit
    task1->SetChargedDCAzMax(DCAzMax); //DCAz max is set to 2 in filterbit 32 and 96 
    task1->SetChargedDCAxyMax(DCAxyMax); // in filterbit 32 and 96 is a pt dependant tight cut and in 768 it is not set at all
    if(PileUp) task1->SetExtraPileUpCut(); 
    task1->SetMaxChi2perTPCcls(MaxChi2perTPC);   
    task1->SetChargedNumTPCclsMin(NumTPCclsMin);
    task1->SetChargedTrackFilterBit(TrackFilterBit);

  mgr->AddTask(task1); // add your task to the manager
  // Creating containers
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(Form("Flow_Refs_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput2 = mgr->CreateContainer(Form("Flow_Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput3 = mgr->CreateContainer(Form("Flow_PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput4 = mgr->CreateContainer(Form("QA_Events_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput5 = mgr->CreateContainer(Form("QA_Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput6 = mgr->CreateContainer(Form("QA_PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput7 = mgr->CreateContainer(Form("Flow_Weights_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));

  // Connecting containers to task
  mgr->ConnectInput(task1,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task1,1,cOutput1);
  mgr->ConnectOutput(task1,2,cOutput2);
  mgr->ConnectOutput(task1,3,cOutput3);
  mgr->ConnectOutput(task1,4,cOutput4);
  mgr->ConnectOutput(task1,5,cOutput5);
  mgr->ConnectOutput(task1,6,cOutput6);
  mgr->ConnectOutput(task1,7,cOutput7);

  return task1;
}
