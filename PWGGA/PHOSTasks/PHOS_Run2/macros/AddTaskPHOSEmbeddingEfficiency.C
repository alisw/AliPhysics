AliAnalysisTaskPHOSEmbeddingEfficiency* AddTaskPHOSEmbeddingEfficiency(
    const char* name     = "EmbeddingEfficiency",
    const TString parname = "Pi0",
    const UInt_t trigger = AliVEvent::kINT7,
    const TString CollisionSystem = "PbPb",
    const Bool_t isMC = kFALSE,
    const TString triggerinput = "",//L1H,L1M,L1L,L0
    const Float_t CenMin = 0.,
    const Float_t CenMax = 90.,
    const Int_t NMixed   = 10,
    const Bool_t FlowTask = kFALSE,
    const Int_t harmonics = -1,
    const Bool_t useCoreE = kTRUE,
    const Bool_t useCoreDisp = kTRUE,
    const Double_t NsigmaCPV  = 3.0,
    const Double_t NsigmaDisp = 3.0,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t TOFcorrection = kTRUE,
    const Bool_t NonLinStudy = kFALSE,
    const Double_t bs = 100.,//bunch space in ns.
    const Double_t distBC = 0.0,//minimum distance to bad channel.
    const Double_t Emin = 0.2,
    const Bool_t isJJMC = kFALSE,
    const TString MCtype = "MBMC"
    )
{
  //Add a task AliAnalysisTaskPHOSEmbeddingEfficiency to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSEmbeddingEfficiency", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSEmbeddingEfficiency", "This task requires an input event handler");
    return NULL;
  }

	TString TriggerName="";
	if     (trigger == (UInt_t)AliVEvent::kAny)  TriggerName = "kAny";
	else if(trigger == (UInt_t)AliVEvent::kINT7) TriggerName = "kINT7";
	else if(trigger == (UInt_t)AliVEvent::kPHI7) TriggerName = "kPHI7";

  if(trigger == (UInt_t)AliVEvent::kPHI7){
    if(triggerinput.Contains("L1") || triggerinput.Contains("L0")){
      TriggerName = TriggerName + "_" + triggerinput;

    }
    else{
      ::Error("AddTaskPHOSEmbeddingEfficiency", "PHOS trigger analysis requires at least trigger input (L0 or L1[H,M,L]).");
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
     if(harmonics > 0) taskname = Form("%s_%s_%s_%s_Cen%d_%d%s_Harmonics%d_BS%dns_DBC%dcell_Emin%dMeV",name,parname.Data(),CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),harmonics,(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));
      else{
        ::Error("AddTaskPHOSEmbeddingEfficiency", "Qn flow vector correction is ON, but you do not set harmonics.");
        return NULL;
      }
  }
  else taskname = Form("%s_%s_%s_%s_Cen%d_%d%s_BS%dns_DBC%dcell_Emin%dMeV",name,parname.Data(),CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));

  AliAnalysisTaskPHOSEmbeddingEfficiency* task = new AliAnalysisTaskPHOSEmbeddingEfficiency(taskname);
  task->SelectCollisionCandidates(trigger);

  if(trigger == (UInt_t)AliVEvent::kPHI7) task->SetPHOSTriggerAnalysis(triggerinput);

  task->SetEmbeddedParticle(parname);
  task->SetCollisionSystem(systemID);//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
  task->SetJetJetMC(isJJMC);
  task->SetMCType(MCtype);
 
  task->SetNonLinearityStudy(NonLinStudy);
 
  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
  task->SetCoreEnergyFlag(useCoreE);

  task->SetEventCuts(isMC);
  task->SetClusterCuts(useCoreDisp,NsigmaCPV,NsigmaDisp,distBC);

  task->SetCentralityMin(CenMin);
  task->SetCentralityMax(CenMax);
  task->SetDepthNMixed(NMixed);
  task->SetQnVectorTask(FlowTask);
  task->SetHarmonics(harmonics);

  //set minimum energy
  task->SetEmin(Emin);

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
    f1tof->SetParameters(0.996,5.61,-0.146,0.036,7.39,0.054);
    task->SetTOFCutEfficiencyFunction(f1tof);
    printf("TOF cut efficiency as a function of E is %s\n",f1tof->GetTitle());
  }

  if(isMC){
    if(parname.Contains("Pi0",TString::kIgnoreCase)){
      //for pi0 pT weighting
      const Int_t Ncen_Pi0 = 7;
      const Double_t centrality_Pi0[Ncen_Pi0] = {0,5,10,20,40,60,100};
      TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0,centrality_Pi0);

      TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
      TF1 *f1weightPi0[Ncen_Pi0-1];

      //fWeightCen0005 = new TF1("fWeightCen0005" ,"64.41*TMath::Power(x,-(5.88 + -92.9/(TMath::Power(x,4.12) + 54.1)))" ,0.,100.);
      //fWeightCen0510 = new TF1("fWeightCen0510" ,"64.88*TMath::Power(x,-(5.90 + -430/ (TMath::Power(x,5.49) + 411)))"   ,0.,100.);
      //fWeightCen1020 = new TF1("fWeightCen1020" ,"43.92*TMath::Power(x,-(5.80 + -396/ (TMath::Power(x,5.25) + 382)))"   ,0.,100.);
      //fWeightCen2040 = new TF1("fWeightCen2040" ,"24.84*TMath::Power(x,-(5.74 + -83.6/(TMath::Power(x,4.01) + 72.3)))" ,0.,100.);
      //fWeightCen4060 = new TF1("fWeightCen4060" ," 8.33*TMath::Power(x,-(5.61 + -25.1/(TMath::Power(x,2.93) + 21.1)))" ,0.,100.);
      //fWeightCen6080 = new TF1("fWeightCen6080" ," 2.08*TMath::Power(x,-(5.66 + -6.97/(TMath::Power(x,1.66) + 4.06)))" ,0.,100.);

      const Double_t p0[Ncen_Pi0-1] = {64.41, 64.88, 43.92, 24.84,  8.33,  2.08};
      const Double_t p1[Ncen_Pi0-1] = { 5.88,  5.90,  5.80,  5.74,  5.61,  5.66};
      const Double_t p2[Ncen_Pi0-1] = {-92.9,  -430,  -396, -83.6, -25.1, -6.97};
      const Double_t p3[Ncen_Pi0-1] = { 4.12,  5.49,  5.25,  4.01,  2.93,  1.66};
      const Double_t p4[Ncen_Pi0-1] = { 54.1,   411,   382,  72.3,  21.1,  4.06};

      for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
        f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"[0]*TMath::Power(x,-([1] + [2]/(TMath::Power(x,[3]) + [4])))",0,100);//this is iterative procedure.
        f1weightPi0[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightPi0[icen]->SetNpx(1000);
        farray_Pi0->Add(f1weightPi0[icen]);
      }
      task->SetAdditionalPi0PtWeightFunction(centarray_Pi0,farray_Pi0);
    }

    else if(parname.Contains("Eta",TString::kIgnoreCase)){
      //for eta pT weighting
      const Int_t Ncen_Eta = 7;
      const Double_t centrality_Eta[Ncen_Eta] = {0,5,10,20,40,60,100};
      TArrayD *centarray_Eta = new TArrayD(Ncen_Eta,centrality_Eta);

      TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      TF1 *f1weightEta[Ncen_Eta-1];

      //fWeightCen0005 = new TF1("fWeightCen0005" ,"28.96*0.8*TMath::Power(x,-(5.85 + -199.17/(TMath::Power(x,4.64)+95.30)))",0.,100.);
      //fWeightCen0510 = new TF1("fWeightCen0510" ,"21.97*0.8*TMath::Power(x,-(5.79 + -33.54/(TMath::Power(x,2.96) +10.84)))",0.,100.);
      //fWeightCen1020 = new TF1("fWeightCen1020" ,"18.91*0.8*TMath::Power(x,-(5.71 + -44.76/(TMath::Power(x,3.37) +19.66)))",0.,100.);
      //fWeightCen2040 = new TF1("fWeightCen2040" ,"11.54*0.8*TMath::Power(x,-(5.74 + -18.43/(TMath::Power(x,2.62) +7.37)))" ,0.,100.);
      //fWeightCen4060 = new TF1("fWeightCen4060" ," 4.18*0.8*TMath::Power(x,-(5.67 + -9.43/(TMath::Power(x,2.00)  +3.39)))" ,0.,100.);
      //fWeightCen6080 = new TF1("fWeightCen6080" ," 1.03*0.8*TMath::Power(x,-(2.01 + 2.32/(TMath::Power(x,-1.67)  +0.661)))",0.,100.);

      const Double_t p0[Ncen_Eta-1] = {64.41, 64.88, 43.92, 24.84,  8.33,  2.08};
      const Double_t p1[Ncen_Eta-1] = { 5.88,  5.90,  5.80,  5.74,  5.61,  5.66};
      const Double_t p2[Ncen_Eta-1] = {-92.9,  -430,  -396, -83.6, -25.1, -6.97};
      const Double_t p3[Ncen_Eta-1] = { 4.12,  5.49,  5.25,  4.01,  2.93,  1.66};
      const Double_t p4[Ncen_Eta-1] = { 54.1,   411,   382,  72.3,  21.1,  4.06};

      for(Int_t icen=0;icen<Ncen_Eta-1;icen++){
        f1weightEta[icen] = new TF1(Form("f1weightEta_%d",icen),"[0]*TMath::Power(x,-([1] + [2]/(TMath::Power(x,[3]) + [4])))",0,100);//this is iterative procedure.
        f1weightEta[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightEta[icen]->SetNpx(1000);
        farray_Eta->Add(f1weightEta[icen]);
      }
      task->SetAdditionalEtaPtWeightFunction(centarray_Eta,farray_Eta);
    }

    else if(parname.Contains("Gamma",TString::kIgnoreCase)){
      //for gamma pT weighting
      const Int_t Ncen_Gamma = 7;
      const Double_t centrality_Gamma[Ncen_Gamma] = {0,5,10,20,40,60,100};
      TArrayD *centarray_Gamma = new TArrayD(Ncen_Gamma,centrality_Gamma);

      TObjArray *farray_Gamma = new TObjArray(Ncen_Gamma-1);
      TF1 *f1weightGamma[Ncen_Gamma-1];

      //fWeightCen0005 = new TF1("fWeightCen0005" ,"28.96*0.8*TMath::Power(x,-(5.85 + -199.17/(TMath::Power(x,4.64)+95.30)))",0.,100.);
      //fWeightCen0510 = new TF1("fWeightCen0510" ,"21.97*0.8*TMath::Power(x,-(5.79 + -33.54/(TMath::Power(x,2.96) +10.84)))",0.,100.);
      //fWeightCen1020 = new TF1("fWeightCen1020" ,"18.91*0.8*TMath::Power(x,-(5.71 + -44.76/(TMath::Power(x,3.37) +19.66)))",0.,100.);
      //fWeightCen2040 = new TF1("fWeightCen2040" ,"11.54*0.8*TMath::Power(x,-(5.74 + -18.43/(TMath::Power(x,2.62) +7.37)))" ,0.,100.);
      //fWeightCen4060 = new TF1("fWeightCen4060" ," 4.18*0.8*TMath::Power(x,-(5.67 + -9.43/(TMath::Power(x,2.00)  +3.39)))" ,0.,100.);
      //fWeightCen6080 = new TF1("fWeightCen6080" ," 1.03*0.8*TMath::Power(x,-(2.01 + 2.32/(TMath::Power(x,-1.67)  +0.661)))",0.,100.);

      const Double_t p0[Ncen_Gamma-1] = {64.41, 64.88, 43.92, 24.84,  8.33,  2.08};
      const Double_t p1[Ncen_Gamma-1] = { 5.88,  5.90,  5.80,  5.74,  5.61,  5.66};
      const Double_t p2[Ncen_Gamma-1] = {-92.9,  -430,  -396, -83.6, -25.1, -6.97};
      const Double_t p3[Ncen_Gamma-1] = { 4.12,  5.49,  5.25,  4.01,  2.93,  1.66};
      const Double_t p4[Ncen_Gamma-1] = { 54.1,   411,   382,  72.3,  21.1,  4.06};

      for(Int_t icen=0;icen<Ncen_Gamma-1;icen++){
        f1weightGamma[icen] = new TF1(Form("f1weightGamma_%d",icen),"[0]*TMath::Power(x,-([1] + [2]/(TMath::Power(x,[3]) + [4])))",0,100);//this is iterative procedure.
        f1weightGamma[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightGamma[icen]->SetNpx(1000);
        farray_Gamma->Add(f1weightGamma[icen]);
      }
      task->SetAdditionalGammaPtWeightFunction(centarray_Gamma,farray_Gamma);
    }

    //const Int_t Ncen_K0S = 7;
    //const Double_t centrality_K0S[Ncen_K0S] = {0,5,10,20,40,60,100};
    //TArrayD *centarray_K0S = new TArrayD(Ncen_K0S,centrality_K0S);
    //TObjArray *farray_K0S = new TObjArray(Ncen_K0S-1);
    //TF1 *f1weightK0S[Ncen_K0S-1];
    //const Double_t p0[Ncen_K0S-1] = {  1.81,   1.83,   1.81,   1.70,   1.91,   1.85};
    //const Double_t p1[Ncen_K0S-1] = {  1.38,   1.31,   1.27,   1.30,   1.61,   1.06};
    //const Double_t p2[Ncen_K0S-1] = {-0.373, -0.398, -0.424, -0.565, -0.640, -0.714};
    //const Double_t p3[Ncen_K0S-1] = {  3.93,   3.70,   4.14,   4.88,  0.542,   1.30};
    //const Double_t p4[Ncen_K0S-1] = {-0.235, -0.279, -0.435,  -2.56, -0.356, -0.543};
    //for(Int_t icen=0;icen<Ncen_K0S-1;icen++){
    //  f1weightK0S[icen] = new TF1(Form("f1weightK0S_%d",icen),"[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1) )",0,100);
    //  f1weightK0S[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
    //  farray_K0S->Add(f1weightK0S[icen]);
    //}
    //task->SetAdditionalK0SPtWeightFunction(centarray_K0S,farray_K0S);
  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString prefix = Form("hist_%s",taskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

