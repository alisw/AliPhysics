AliAnalysisTaskPHOSPi0EtaToGammaGamma* AddTaskPHOSPi0EtaToGammaGamma_PbPb_5TeV(
    const char* name     = "Pi0EtaToGammaGamma",
    UInt_t trigger = AliVEvent::kINT7,
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
    const Double_t Emin = 0.2,//minimum energy for photon selection in GeV
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
  else if(trigger & AliVEvent::kPHI7)          TriggerName = "kPHI7";

  if(trigger & AliVEvent::kPHI7){
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

    taskname = Form("%s_%s_%s_Cen%d_%d%s_Harmonics%d_%s_%s_BS%dns_DBC%dcell_Emin%dMeV",name,CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),harmonics,FMname.Data(),detname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));

  }
  else taskname = Form("%s_%s_%s_Cen%d_%d%s_BS%dns_DBC%dcell_Emin%dMeV",name,CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));

  if((trigger & AliVEvent::kPHI7) && ApplyTOFTrigger) taskname += "_TOFTrigger";
  if(ForceActiveTRU) taskname += "_ForceActiveTRU";

  AliAnalysisTaskPHOSPi0EtaToGammaGamma* task = new AliAnalysisTaskPHOSPi0EtaToGammaGamma(taskname);

  Double_t Ethre = 0.0;
  if(L1input == 7)       Ethre = 0.0;
  else if(L1input == 6)  Ethre = 0.0;
  else if(L1input == 5)  Ethre = 0.0;
  else if(L0input == 9)  Ethre = 0.0;//LHC15n
  else if(L0input == 17) Ethre = 0.0;//LHC17p
  if(trigger & AliVEvent::kPHI7)      task->SetPHOSTriggerAnalysis(L1input,L0input,Ethre,isMC,ApplyTOFTrigger,-1);
  else if(trigger & AliVEvent::kINT7) task->SetPHOSTriggerAnalysisMB(L1input,L0input,Ethre,isMC,ApplyTOFTrigger,-1);
  if(isMC && (trigger & AliVEvent::kPHI7)) trigger = AliVEvent::kINT7;//change trigger selection in MC when you do PHOS trigger analysis.
  if(ForceActiveTRU) task->SetForceActiveTRU(L1input,L0input,Ethre,isMC);//this is to measure rejection factor from cluster energy kPHI7/kINT7 with same acceptance.
  task->SelectCollisionCandidates(trigger);
  task->SetTriggerThreshold(Ethre);

  task->SetCollisionSystem(systemID);//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
  task->SetJetJetMC(isJJMC);
  task->SetMCType(MCtype);

  task->SetNonLinearityStudy(NonLinStudy);

  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
  //  task->SetCoreEnergyFlag(useCoreE);

  task->SetEventCuts(isMC,pf);
  task->SetClusterCuts(useCoreDisp,NsigmaCPV,NsigmaDisp,useCoreE,distBC);

  task->SetCentralityMin(CenMin);
  task->SetCentralityMax(CenMax);
  task->SetDepthNMixed(NMixed);
  task->SetQnVectorTask(FlowTask);
  task->SetHarmonics(harmonics);
  task->SetFlowMethod(FlowMethod);
  task->SetQnDetector(QnDetector);

  //set minimum energy
  task->SetEmin(Emin);

  //centrality setting
  task->SetCentralityEstimator("V0M");

  AliESDtrackCuts *cutsG = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);//standard cuts with very loose DCA
  cutsG->SetMaxDCAToVertexXY(2.4);
  cutsG->SetMaxDCAToVertexZ(3.2);
  cutsG->SetDCAToVertex2D(kTRUE);
  task->SetESDtrackCutsForGlobal(cutsG);

  AliESDtrackCuts *cutsGC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();//standard cuts with tight DCA cut
  task->SetESDtrackCutsForGlobalConstrained(cutsGC);

  //bunch space for TOF cut
  task->SetBunchSpace(bs);//in unit of ns.
  if(!isMC && TOFcorrection){
    //TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) - ( 0 + [3]/(exp( -(x-[4]) / [5] ) + 1)  )",0,100);
    //f1tof->SetParameters(0.996,5.61,-0.146,0.036,7.39,0.054);

    //TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) * ( 1 + [3]/(TMath::TwoPi() * [5]) * exp(-(x-[4]) * (x-[4])/(2 * [5]*[5] )) *( 1 + TMath::Erf([6]*((x-[4])/[5]))) )",0,100);
    //f1tof->SetNpx(1000);
    //f1tof->FixParameter(0, 9.97378e-01);
    //f1tof->FixParameter(1, 6.66818e+00);
    //f1tof->FixParameter(2,-5.65437e-02);
    //f1tof->FixParameter(3,-7.38995e-01);
    //f1tof->FixParameter(4, 9.72815e+00);
    //f1tof->FixParameter(5, 1.17920e+00);
    //f1tof->FixParameter(6, 5.55735e-02);

    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0]/(exp(-(x-[1])/[2]) + 1) - ([3]/(exp( -(x-[4]) / [5] ) + 1))",0,100);
    f1tof->SetNpx(1000);
    f1tof->FixParameter(0,9.99550e-01);
    f1tof->FixParameter(1,7.61897e-03);
    f1tof->FixParameter(2,1.42936e-01);
    f1tof->FixParameter(3,3.70000e-02);
    f1tof->FixParameter(4,7.17525e+00);
    f1tof->FixParameter(5,4.66735e-01);
    task->SetTOFCutEfficiencyFunction(f1tof);
  }
  if(!isMC && Trgcorrection){
    if(L1input == 7){
      printf("L1H is selected\n");
      //TF1 *f1trg = new TF1("f1TriggerEfficiency","[0]/(TMath::Exp(-(x-[1])/[2]) + 1)",0,100);
      //f1trg->SetParameters(0.431,8.83,0.79);//from MB //acc x trigger efficiency 6-50GeV
      TF1 *f1trg = new TF1("f1TriggerEfficiency","[0]/(TMath::Exp(-(x-[1])/[2]) + 1) + [3]/(TMath::Exp(-(x-[4])/[5]) + 1)",0,100);
      //f1trg->SetParameters(0.218,8.02,0.588,0.236,10.3,0.883);//from MB //acc x trigger efficiency 5-40GeV
      f1trg->SetParameters(0.180,7.89,0.544,0.270,9.84,0.875);//from MB //acc x trigger efficiency 5-40GeV
      f1trg->SetNpx(1000);
      task->SetTriggerEfficiency(f1trg);
    }
    else if(L1input==6){
      printf("L1M is selected\n");
      //TF1 *f1trg = new TF1("f1TriggerEfficiency","[0]/(TMath::Exp(-(x-[1])/[2]) + 1)",0,100);
      //f1trg->SetParameters(0.445,4.43,0.72);//from MB //acc x trigger efficiency 6-50GeV
      TF1 *f1trg = new TF1("f1TriggerEfficiency","[0]/(TMath::Exp(-(x-[1])/[2]) + 1) + [3]/(TMath::Exp(-(x-[4])/[5]) + 1)",0,100);
      //f1trg->SetParameters(0.220,3.87,0.334,0.230,5.32,0.523);//from MB //acc x trigger efficiency 2-40GeV
      f1trg->SetParameters(0.230,3.87,0.339,0.220,5.22,0.496);//from MB //acc x trigger efficiency 2-40GeV
      f1trg->SetNpx(1000);
      task->SetTriggerEfficiency(f1trg);
    }
  }

  if(isMC){
    //for pi0
    const Int_t Ncen_Pi0 = 11;
    const Double_t centrality_Pi0[Ncen_Pi0] = {0,5,10,20,30,40,50,60,70,80,100};
    TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0,centrality_Pi0);

    TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
    TF1 *f1weightPi0[Ncen_Pi0-1];
    const Double_t p0_Pi0[Ncen_Pi0-1] = {2. * 8.52796e-02,2. * 9.57970e-02,2. * 1.09042e-01,2. * 1.28762e-01 ,2. * 1.51087e-01 ,2. * 1.82705e-01 ,2. * 2.16360e-01 ,2. * 2.37666e-01 ,2. * 2.52706e-01 ,2. * 3.34001e-01};
    const Double_t p1_Pi0[Ncen_Pi0-1] = {2. * 8.34243e-01,2. * 8.11715e-01,2. * 7.73274e-01,2. * 7.28962e-01 ,2. * 6.77506e-01 ,2. * 6.06502e-01 ,2. * 5.31093e-01 ,2. * 4.52193e-01 ,2. * 3.86976e-01 ,2. * 3.22488e-01};
    const Double_t p2_Pi0[Ncen_Pi0-1] = {9.27577e-01,9.53380e-01,9.52280e-01,9.78872e-01 ,9.82192e-01 ,1.01124e+00 ,1.08236e+00 ,1.14572e+00 ,1.12243e+00 ,2.16920e+00};
    const Double_t p3_Pi0[Ncen_Pi0-1] = {2.13453e-01,2.09818e-01,2.03573e-01,2.00238e-01 ,1.94211e-01 ,1.87993e-01 ,1.94509e-01 ,1.95069e-01 ,1.75698e-01 ,3.62140e-01};

    for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
      f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"[0] + [1]*pow(x,[2])*exp(-[3]*x*x)",0,100);
      f1weightPi0[icen]->SetParameters(p0_Pi0[icen],p1_Pi0[icen],p2_Pi0[icen],p3_Pi0[icen]);
      farray_Pi0->Add(f1weightPi0[icen]);
    }

    task->SetAdditionalPi0PtWeightFunction(centarray_Pi0,farray_Pi0);//do not change pi0 spectra in MC

    //for K0S
    const Int_t Ncen_K0S = 7;
    const Double_t centrality_K0S[Ncen_K0S] = {0,5,10,20,40,60,100};
    TArrayD *centarray_K0S = new TArrayD(Ncen_K0S,centrality_K0S);

    TObjArray *farray_K0S = new TObjArray(Ncen_K0S-1);
    TF1 *f1weightK0S[Ncen_K0S-1];
    const Double_t p0_K0S[Ncen_K0S-1] = {  1.81,   1.83,   1.81,   1.70,   1.91,   1.85};
    const Double_t p1_K0S[Ncen_K0S-1] = {  1.38,   1.31,   1.27,   1.30,   1.61,   1.06};
    const Double_t p2_K0S[Ncen_K0S-1] = {-0.373, -0.398, -0.424, -0.565, -0.640, -0.714};
    const Double_t p3_K0S[Ncen_K0S-1] = {  3.93,   3.70,   4.14,   4.88,  0.542,   1.30};
    const Double_t p4_K0S[Ncen_K0S-1] = {-0.235, -0.279, -0.435,  -2.56, -0.356, -0.543};

    for(Int_t icen=0;icen<Ncen_K0S-1;icen++){
      f1weightK0S[icen] = new TF1(Form("f1weightK0S_%d",icen),"[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1) )",0,100);
      f1weightK0S[icen]->SetParameters(p0_K0S[icen],p1_K0S[icen],p2_K0S[icen],p3_K0S[icen],p4_K0S[icen]);
      farray_K0S->Add(f1weightK0S[icen]);
    }

    task->SetAdditionalK0SPtWeightFunction(centarray_K0S,farray_K0S);

    //for L0
    const Int_t Ncen_L0 = 7;
    const Double_t centrality_L0[Ncen_L0] = {0,5,10,20,40,60,100};
    TArrayD *centarray_L0 = new TArrayD(Ncen_L0,centrality_L0);

    TObjArray *farray_L0 = new TObjArray(Ncen_L0-1);
    TF1 *f1weightL0[Ncen_L0-1];
    const Double_t p0_L0[Ncen_L0-1] = {27.2   ,272    ,548    , 496    ,1880   ,2220  };
    const Double_t p1_L0[Ncen_L0-1] = {0.0562 ,0.0420 ,0.0385 , 0.0420 ,0.0442 ,0.0562};
    const Double_t p2_L0[Ncen_L0-1] = {-7.79  ,-11.6  ,-12.8  , -12.2  ,-13.1  ,-11.7 };
    const Double_t p3_L0[Ncen_L0-1] = {0.231  ,0.245  ,0.286  , 0.371  ,0.470  ,0.570 };

    for(Int_t icen=0;icen<Ncen_L0-1;icen++){
      f1weightL0[icen] = new TF1(Form("f1weightL0_%d",icen),"[0] * TMath::Power(x,4) * TMath::Exp(-[1] * TMath::Power(x-[2],2)) + [3]",0,100);
      f1weightL0[icen]->SetParameters(p0_L0[icen],p1_L0[icen],p2_L0[icen],p3_L0[icen]);
      farray_L0->Add(f1weightL0[icen]);
    }

    //task->SetAdditionalL0PtWeightFunction(centarray_L0,farray_L0);

  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

  TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString prefix = Form("hist_%s",taskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

