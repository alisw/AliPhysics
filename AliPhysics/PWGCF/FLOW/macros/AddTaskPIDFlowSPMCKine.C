class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;

void AddTaskPIDFlowSPMCKine(Int_t uptoWhichHarmonics = 5, 
			    Float_t etamin=-0.8,
			    Float_t etamax=0.8,
			    Float_t gEtaGap=0.0,
			    TString fileNameBase="AnalysisResults",
			    TString uniqueStr="Pion_00",
			    TString gQVector ="Qa",
			    Int_t charge=0,
			    Double_t bmin = 0.0,
			    Double_t bmax = 12.0,
			    Bool_t isPID = kTRUE,
			    Bool_t isVZERO = kFALSE, // use vzero sp method
			    Bool_t useAfterBurner=kFALSE,
			    AliPID::EParticleType particleType=AliPID::kPion) {
  //macro to run over the kine tree for higher harmonics (pid with pdg)
  TF1 *gV2Param = 0x0;
  Bool_t doQA=kTRUE;
    
  // AFTERBURNER
  Double_t v1=0.0;
  Double_t v2=0.0;
  Double_t v3=0.0;
  Double_t v4=0.0;
  Int_t numberOfTrackClones=0; //non-flow
    
  //Define the range for eta subevents (for SP method) with TPC
  Double_t minA = etamin;//
  Double_t maxA = -0.5*gEtaGap;//
  Double_t minB = +0.5*gEtaGap;//
  Double_t maxB = etamax;//
    
  //---------Data selection---------- ESD only!!!
  //kMC, kGlobal, kESD_TPConly, kESD_SPDtracklet
  AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kMC;
  AliFlowTrackCuts::trackParameterType poitype = AliFlowTrackCuts::kMC;
  
  //---------Parameter mixing--------
  //kPure - no mixing, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt
  AliFlowTrackCuts::trackParameterMix rpmix = AliFlowTrackCuts::kPure;
  AliFlowTrackCuts::trackParameterMix poimix = AliFlowTrackCuts::kPure;
  
  const char* rptypestr = AliFlowTrackCuts::GetParamTypeName(rptype);  //ESD
  const char* poitypestr = AliFlowTrackCuts::GetParamTypeName(poitype); //ESD
  
  //general parameters    
  const int nharmonics = uptoWhichHarmonics-1;
  TString outputSlotName[nharmonics];
  TString suffixName;
  
  //========================================//
  //Define event cuts
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("eventcuts");
  if(bmin > = 0) cutsEvent->SetImpactParameterRange(bmin,bmax);
  cutsEvent->SetPrimaryVertexZrange(-10.,10.);
  cutsEvent->SetQA(doQA);
  
  //========================================//
  //Definition of reference particles
  AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts("RP");
  cutsRP->SetParamType(rptype);
  cutsRP->SetParamMix(rpmix);
  cutsRP->SetPtRange(0.2,10.);
  cutsRP->SetEtaRange(etamin,etamax);
  cutsRP->SetMCisPrimary(kTRUE);
  cutsRP->SetRequireCharge(kTRUE);
  cutsRP->SetQA(doQA);
  
  //========================================//
  //Definition of POIs
  AliFlowTrackCuts* SP_POI = new AliFlowTrackCuts("POI");
  SP_POI->SetParamType(poitype);
  SP_POI->SetParamMix(poimix);
  SP_POI->SetMCisPrimary(kTRUE);
  SP_POI->SetPtRange(0.2,10.);
  SP_POI->SetIgnoreSignInMCPID(kTRUE);
  if(isPID){
    if(particleType==AliPID::kPion) SP_POI->SetMCPID(211);
    if(particleType==AliPID::kKaon) SP_POI->SetMCPID(321);
    if(particleType==AliPID::kProton) SP_POI->SetMCPID(2212);
  }        
  if(charge != 0) SP_POI->SetCharge(charge);
  
  //=====================================================================
  //Define the sub events for the POIs
  if(!isVZERO && gQVector=="Qa") {
    SP_POI->SetEtaRange( +0.5*gEtaGap, etamax );
    printf(" > NOTE: Using half TPC (Qb) as POI selection u < \n");
    suffixName = "Qa";
  }
  if(!isVZERO && gQVector=="Qb") {
    SP_POI->SetEtaRange( etamin,-0.5*gEtaGap );
    printf(" > NOTE: Using half TPC (Qa) as POI selection u < \n");
    suffixName = "Qb";
  }
  if(isVZERO) {
    SP_POI->SetEtaRange( etamin,etamax );
    printf(" > NOTE: Using full TPC as POI selection u < \n");
    suffixName = "vzero";
  }
  SP_POI->SetQA(doQA);
  //=====================================================================
    
  suffixName += "_flow_";
  if(bmin >= 0) {
    suffixName += Form("%.2f_", bmin);
    suffixName += Form("%.2f_", bmax);
  }
  suffixName += Form("%.f_", gEtaGap*10);
    
  if(isPID) 
    suffixName+=uniqueStr;
  else
    suffixName+="AllCharged";
  
  if (charge<0) suffixName+="-";
  if (charge>0) suffixName+="+";
  
  //for v2,v3,v4 and v5
  for(int harmonic = 2; harmonic < uptoWhichHarmonics + 1; harmonic++){  
    outputSlotName[harmonic-2] = "";
    outputSlotName[harmonic-2]+=uniqueStr;
    outputSlotName[harmonic-2]+=Form("_v%i_",harmonic);
    if(bmin >= 0) {
      outputSlotName[harmonic-2]+=Form("%.2f-",bmin);
      outputSlotName[harmonic-2]+=Form("%.2f_",bmax);
    }
    if(isPID)
      outputSlotName[harmonic-2]+=AliPID::ParticleName(particleType);
    else
      outputSlotName[harmonic-2]+="AllCharged";
    if (charge<0) outputSlotName[harmonic-2]+="-";
    if (charge>0) outputSlotName[harmonic-2]+="+";
  }    
    
  TString fileName(fileNameBase);
  fileName.Append(".root");
        
  //========================FLOWPACKAGE TASKS=========================//
  AliAnalysisDataContainer *cinput1;
  AliAnalysisDataContainer *coutputFE;
  AliAnalysisDataContainer* coutputFEQA;
  AliAnalysisTaskFlowEvent *taskFE;
  
  AliAnalysisDataContainer *flowEvent[nharmonics];
  AliAnalysisTaskFilterFE *tskFilter[nharmonics];
  
  AliAnalysisDataContainer *coutputSP[nharmonics];
  AliAnalysisTaskScalarProduct *taskSP[nharmonics];
  
  TString outputQA;
  TString myNameSP[nharmonics];
  TString slot[nharmonics];
  
  //======================================================================
  //Get the pointer to the existing analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowEvent", "No analysis manager to connect to.");
    return NULL;
  }
  
  //Check the analysis type using the event handlers connected to the 
  //analysis manager. The availability of MC handler can also be checked here.
  //=====================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlowEvent", "This task requires an input event handler");
    return NULL;
  }
  
  taskFE = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName.Data()),"",doQA); 
  
  //Use after burner: rotate phi for each particle
  if(useAfterBurner) {
    if(gV2Param) 
      taskFE->SetPtDifferentialV2(gV2Param);
    else 
      taskFE->SetFlow(v1,v2,v3,v4);
    taskFE->SetNonFlowNumberOfTrackClones(numberOfTrackClones);
    taskFE->SetAfterburnerOn();
  }
  
  if(!isVZERO) taskFE->SetSubeventEtaRange(minA, maxA, minB, maxB);
  if(isVZERO)  taskFE->SetSubeventEtaRange(-5,-1.5,+1.5,5);
  mgr->AddTask(taskFE);
  
  // Pass cuts for RPs and POIs to the task:
  taskFE->SetCutsEvent(cutsEvent);
  taskFE->SetCutsRP(cutsRP);
  taskFE->SetCutsPOI(SP_POI);
  if (cutsRP->GetParamType()==AliFlowTrackCuts::kVZERO) {
    taskFE->SetHistWeightvsPhiMin(0.);
    taskFE->SetHistWeightvsPhiMax(200.);
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  cinput1 = mgr->GetCommonInputContainer();
  
  coutputFE = mgr->CreateContainer(Form("FlowEvent_%s",suffixName.Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  
  mgr->ConnectInput(taskFE,0,cinput1);
  mgr->ConnectOutput(taskFE,1,coutputFE);
  //==========================================================
  TString Species = "";
  if(isPID) 
    Species += AliPID::ParticleName(particleType);
  else
    Species += "Allcharged";
          
  for(int harm = 2; harm < uptoWhichHarmonics + 1; harm++){
    myNameSP[harm-2] = "SP_";
    myNameSP[harm-2] += gQVector;
    myNameSP[harm-2] += Form("_v%i_%s_%.f",harm,outputSlotName[harm-2].Data(),gEtaGap*10);
    
    flowEvent[harm-2] = mgr->CreateContainer( Form("Filter_%s", myNameSP[harm-2].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer );
    
    tskFilter[harm-2] = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP[harm-2].Data()),cutsRP, NULL);//SP_POI
    if(!isVZERO)
      tskFilter[harm-2]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) tskFilter[harm-2]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(tskFilter[harm-2]);
    mgr->ConnectInput( tskFilter[harm-2],0,coutputFE);
    mgr->ConnectOutput(tskFilter[harm-2],1,flowEvent[harm-2]);
            
    taskSP[harm-2] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotName[harm-2].Data()),kFALSE);
    taskSP[harm-2]->SetHarmonic(harm);
    taskSP[harm-2]->SetRelDiffMsub(1.0);
    taskSP[harm-2]->SetTotalQvector(gQVector);
    taskSP[harm-2]->SetApplyCorrectionForNUA(kTRUE);
    
    TString outputSP = fileName;
    outputSP += ":outputSPanalysis";
    outputSP+= rptypestr;
    slot[harm-2] = "SP_";
    slot[harm-2] += outputSlotName[harm-2];
    slot[harm-2] += "_";
    slot[harm-2] += gQVector;
    coutputSP[harm-2] = mgr->CreateContainer(Form("%s_%.f",slot[harm-2].Data(),gEtaGap*10),TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
    mgr->AddTask(taskSP[harm-2]);
    mgr->ConnectInput(taskSP[harm-2],0,flowEvent[harm-2]);
    mgr->ConnectInput(taskSP[harm-2],0,coutputFE);
    mgr->ConnectOutput(taskSP[harm-2],1,coutputSP[harm-2]);
  }
  
  if (taskFE->GetQAOn()) {
    outputQA = fileName;
    outputQA += ":QA";
    coutputFEQA = mgr->CreateContainer(Form("QA_%s",suffixName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA);
    mgr->ConnectOutput(taskFE,2,coutputFEQA);
  }
}

