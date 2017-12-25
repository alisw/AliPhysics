AliAnalysisTaskPHOSPi0EtaToGammaGamma* AddTaskPHOSPi0EtaToGammaGamma_PbPb_5TeV(
    const char* name     = "Pi0EtaToGammaGamma",
    const UInt_t trigger = AliVEvent::kINT7,
    const TString CollisionSystem = "PbPb",
    const Bool_t isMC = kFALSE,
    const Int_t L1input = -1,//L1H,L1M,L1L
    const Int_t L0input = -1,//L0
    const Float_t CenMin = 0.,
    const Float_t CenMax = 90.,
    const Int_t NMixed   = 10,
    const Bool_t FlowTask = kFALSE,
    const Int_t harmonics = -1,
    const Int_t FlowMethod = -1,
    const Int_t QnDetector = -1,
    const Bool_t useCoreE = kFALSE,
    const Bool_t useCoreDisp = kFALSE,
    const Double_t NsigmaCPV  = 2.5,
    const Double_t NsigmaDisp = 2.5,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t TOFcorrection = kTRUE,
    const Bool_t Trgcorrection = kFALSE,
    const Bool_t NonLinStudy = kFALSE,
    const Double_t bs = 100.,//bunch space in ns.
    const Double_t distBC = -1,//minimum distance to bad channel.
    const Bool_t isJJMC = kFALSE,
    const TString MCtype = "MBMC",
    const Bool_t ForceActiveTRU = kFALSE,
    const Bool_t ApplyTOFTrigger = kFALSE,
    const AliPHOSEventCuts::PileupFinder pf = AliPHOSEventCuts::kMultiVertexer
    )
{
  //Add a task AliAnalysisTaskPHOSPi0EtaToGammaGamma to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPi0EtaToGammaGamma", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPi0EtaToGammaGamma", "This task requires an input event handler");
    return NULL;
  }

	TString TriggerName="";
	if     (trigger == (UInt_t)AliVEvent::kAny)  TriggerName = "kAny";
	else if(trigger == (UInt_t)AliVEvent::kINT7) TriggerName = "kINT7";
	else if(trigger == (UInt_t)AliVEvent::kPHI7) TriggerName = "kPHI7";

  //if(trigger == (UInt_t)AliVEvent::kPHI7){
  //  if(triggerinput.Contains("L1") || triggerinput.Contains("L0")){
  //    TriggerName = TriggerName + "_" + triggerinput;

  //  }
  //  else{
  //    ::Error("AddTaskPHOSPi0EtaToGammaGamma", "PHOS trigger analysis requires at least trigger input (L0 or L1[H,M,L]).");
  //    return NULL;
  //  }
  //}
  if(trigger == (UInt_t)AliVEvent::kPHI7){
    if(L1input > 0){
      if(L1input == 7)      TriggerName = TriggerName + "_" + "L1H";
      else if(L1input == 6) TriggerName = TriggerName + "_" + "L1M";
      else if(L1input == 5) TriggerName = TriggerName + "_" + "L1L";
    }
    else if(L0input > 0)    TriggerName = TriggerName + "_" + "L0";
    else{
      ::Error("AddTaskPHOSPi0EtaToGammaGamma", "PHOS trigger analysis requires at least 1 trigger input (L0 or L1[H,M,L]).");
      return NULL;
    }
  }

  Int_t systemID = -1;
  if(CollisionSystem=="pp")                                 systemID = 0;
  else if(CollisionSystem=="PbPb")                          systemID = 1;
  else if(CollisionSystem=="pPb" || CollisionSystem=="Pbp") systemID = 2;

  TString PIDname="";
  if(NsigmaCPV > 0) PIDname += Form("_CPV%d",(Int_t)(NsigmaCPV*10));
  if(NsigmaDisp > 0){
    if(useCoreDisp) PIDname += Form("_CoreDisp%d",(Int_t)(NsigmaDisp*10));
    else            PIDname += Form("_FullDisp%d",(Int_t)(NsigmaDisp*10));
  }
  if(useCoreE) PIDname += "_CoreE";
  else         PIDname += "_FullE";

  TString taskname = "";
  if(FlowTask){
    if(harmonics <= 0){
      ::Error("AddTaskPHOSPi0EtaToGammaGamma", "Qn flow vector correction is ON, but you do not set harmonics.");
      return NULL;
    }

    TString FMname = "";
    if(FlowMethod == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      FMname = "EP";
    else if(FlowMethod == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) FMname = "SP";

    TString detname = "";
    if(QnDetector == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullTPC)        detname = "FullTPC";
    else if(QnDetector == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCNegEta) detname = "TPCNegEta";
    else if(QnDetector == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCPosEta) detname = "TPCPosEta";
    else if(QnDetector == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullV0)    detname = "FullV0";
    else if(QnDetector == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A)       detname = "V0A";
    else if(QnDetector == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C)       detname = "V0C";
    else{
      ::Error("AddTaskPHOSPi0EtaToGammaGamma", "detector to measure Qn vector does not exist.");
      return NULL;
    }

    taskname = Form("%s_%s_%s_Cen%d_%d%s_Harmonics%d_%s_%s_BS%dns_DBC%dcell",name,CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),harmonics,FMname.Data(),detname.Data(),(Int_t)bs,(Int_t)(distBC));

  }
  else taskname = Form("%s_%s_%s_Cen%d_%d%s_BS%dns_DBC%dcell",name,CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),(Int_t)bs,(Int_t)(distBC));

  if(trigger == (UInt_t)AliVEvent::kPHI7 && ApplyTOFTrigger) taskname += "_TOFTrigger";

  if(ForceActiveTRU) taskname += "_ForceActiveTRU";

  AliAnalysisTaskPHOSPi0EtaToGammaGamma* task = new AliAnalysisTaskPHOSPi0EtaToGammaGamma(taskname);

  Double_t Ethre = 0.0;
  if(L1input == 7)       Ethre = 0.0;
  else if(L1input == 6)  Ethre = 0.0;
  else if(L1input == 5)  Ethre = 0.0;
  else if(L0input == 9)  Ethre = 0.0;//LHC15n
  else if(L0input == 17) Ethre = 0.0;//LHC17p
  if(trigger == (UInt_t)AliVEvent::kPHI7) task->SetPHOSTriggerAnalysis(L1input,L0input,Ethre,isMC,ApplyTOFTrigger);
  if(kMC && trigger == (UInt_t)AliVEvent::kPHI7) trigger = AliVEvent::kINT7;//change trigger selection in MC when you do PHOS trigger analysis.
  if(ForceActiveTRU) task->SetForceActiveTRU(L1input,L0input,Ethre,isMC);//this is to measure rejection factor from cluster energy kPHI7/kINT7 with same acceptance.
  task->SelectCollisionCandidates(trigger);

  task->SetCollisionSystem(systemID);//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
  task->SetJetJetMC(isJJMC);
  task->SetMCType(MCtype);
 
  task->SetNonLinearityStudy(NonLinStudy);
 
  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
  task->SetCoreEnergyFlag(useCoreE);

  task->SetEventCuts(isMC,pf);
  task->SetClusterCuts(useCoreDisp,NsigmaCPV,NsigmaDisp,distBC);

  task->SetCentralityMin(CenMin);
  task->SetCentralityMax(CenMax);
  task->SetDepthNMixed(NMixed);
  task->SetQnVectorTask(FlowTask);
  task->SetHarmonics(harmonics);
  task->SetFlowMethod(FlowMethod);
  task->SetQnDetector(QnDetector);


  //centrality setting
  task->SetCentralityEstimator("V0M");

  //setting esd track selection for hybrid track
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
  AliESDtrackCuts *cutsG = CreateTrackCutsPWGJE(10001008);//for good global tracks
  task->SetESDtrackCutsForGlobal(cutsG);
  AliESDtrackCuts *cutsGC = CreateTrackCutsPWGJE(10011008);//for good global-constrained tracks
  task->SetESDtrackCutsForGlobalConstrained(cutsGC);

  //bunch space for TOF cut
  task->SetBunchSpace(bs);//in unit of ns.
  if(!isMC && TOFcorrection){
    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) - ( 0 + [3]/(exp( -(x-[4]) / [5] ) + 1)  )",0,100);
    f1tof->SetNpx(1000);
    f1tof->SetParameters(0.996,5.61,-0.146,0.036,7.39,0.054);
    task->SetTOFCutEfficiencyFunction(f1tof);
  }

  if(isMC){
    //const Int_t centrality_Pi0[] = {0,10,30,50,90};
    //const Int_t Ncen_Pi0 = sizeof(centrality_Pi0)/sizeof(centrality_Pi0[0]);
    //TF1 *f1Pi0Weight = new TF1("f1Pi0Weight","1.",0,100);
    //task->SetAdditionalPi0PtWeightFunction(f1Pi0Weight);

    const Int_t Ncen_K0S = 7;
    const Double_t centrality_K0S[Ncen_K0S] = {0,5,10,20,40,60,100};
    TArrayD *centarray_K0S = new TArrayD(Ncen_K0S,centrality_K0S);

    TObjArray *farray_K0S = new TObjArray(Ncen_K0S-1);
    TF1 *f1weightK0S[Ncen_K0S-1];
    const Double_t p0[Ncen_K0S-1] = {  1.81,   1.83,   1.81,   1.70,   1.91,   1.85};
    const Double_t p1[Ncen_K0S-1] = {  1.38,   1.31,   1.27,   1.30,   1.61,   1.06};
    const Double_t p2[Ncen_K0S-1] = {-0.373, -0.398, -0.424, -0.565, -0.640, -0.714};
    const Double_t p3[Ncen_K0S-1] = {  3.93,   3.70,   4.14,   4.88,  0.542,   1.30};
    const Double_t p4[Ncen_K0S-1] = {-0.235, -0.279, -0.435,  -2.56, -0.356, -0.543};

    for(Int_t icen=0;icen<Ncen_K0S-1;icen++){
      f1weightK0S[icen] = new TF1(Form("f1weightK0S_%d",icen),"[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1) )",0,100);
      f1weightK0S[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
      farray_K0S->Add(f1weightK0S[icen]);
    }

    task->SetAdditionalK0SPtWeightFunction(centarray_K0S,farray_K0S);
  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString prefix = Form("hist_%s",taskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

