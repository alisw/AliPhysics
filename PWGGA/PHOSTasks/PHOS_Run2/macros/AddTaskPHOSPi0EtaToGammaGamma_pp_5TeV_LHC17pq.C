AliAnalysisTaskPHOSPi0EtaToGammaGamma* AddTaskPHOSPi0EtaToGammaGamma_pp_5TeV_LHC17pq(
    const char* name     = "Pi0EtaToGammaGamma",
    UInt_t trigger = AliVEvent::kINT7,
    const TString CollisionSystem = "pp",
    const Bool_t isMC = kFALSE,
    const Int_t L1input = -1,//L1H,L1M,L1L
    const Int_t L0input = -1,//L0
    const Float_t CenMin = 0.,
    const Float_t CenMax = 90.,
    const Int_t NMixed   = 10,
    const Bool_t FlowTask = kFALSE,
    const Int_t harmonics = -1,
    const Bool_t useCoreE = kFALSE,
    const Bool_t useCoreDisp = kFALSE,
    const Double_t NsigmaCPV  = 2.5,
    const Double_t NsigmaDisp = 2.5,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t TOFcorrection = kTRUE,
    const Bool_t Trgcorrection = kFALSE,
    const Bool_t NonLinStudy = kFALSE,
    const Double_t bs = 25.,//bunch space in ns.
    const Double_t distBC = 0,//minimum distance to bad channel.
    const Double_t Emin = 0.2,//minimum energy for photon selection in GeV
    const Bool_t isJJMC = kFALSE,
    const TString MCtype = "MBMC",
    const Bool_t ForceActiveTRU = kFALSE,
    const Bool_t ApplyTOFTrigger = kFALSE,
    const AliPHOSEventCuts::PileupFinder pf = AliPHOSEventCuts::kSPDInMultBins
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
	else if(trigger & AliVEvent::kCaloOnly)      TriggerName = "kCaloOnly";

  if(trigger & AliVEvent::kPHI7 || trigger & AliVEvent::kCaloOnly){
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
     if(harmonics > 0) taskname = Form("%s_%s_%s_Cen%d_%d%s_Harmonics%d_BS%dns_DBC%dcell_Emin%dMeV",name,CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),harmonics,(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));
      else{
        ::Error("AddTaskPHOSPi0EtaToGammaGamma", "Qn flow vector correction is ON, but you do not set harmonics.");
        return NULL;
      }
  }
  else taskname = Form("%s_%s_%s_Cen%d_%d%s_BS%dns_DBC%dcell_Emin%dMeV",name,CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));
  if((trigger & AliVEvent::kPHI7) && ApplyTOFTrigger) taskname += "_TOFTrigger";

  if(ForceActiveTRU) taskname += "_ForceActiveTRU";

  AliAnalysisTaskPHOSPi0EtaToGammaGamma* task = new AliAnalysisTaskPHOSPi0EtaToGammaGamma(taskname);

  Double_t Ethre = 0.0;
  if(L1input == 7)       Ethre = 8.0;
  else if(L1input == 6)  Ethre = 6.0;
  else if(L1input == 5)  Ethre = 4.0;
  else if(L0input == 9)  Ethre = 0.0;//LHC15n//threshold was s et to 3 GeV, but efficiency can be measured down to 2 GeV
  else if(L0input == 17) Ethre = 0.0;//LHC17p//threshold was s et to 4 GeV, but efficiency can be measured down to 2 GeV

  if(trigger & AliVEvent::kPHI7 || trigger & AliVEvent::kCaloOnly) task->SetPHOSTriggerAnalysis(L1input,L0input,Ethre,isMC,ApplyTOFTrigger,-1);
  else if(trigger & AliVEvent::kINT7) task->SetPHOSTriggerAnalysisMB(L1input,L0input,Ethre,isMC,ApplyTOFTrigger,-1);
  if(isMC && (trigger & AliVEvent::kPHI7 || trigger & AliVEvent::kCaloOnly)) trigger = AliVEvent::kINT7;//change trigger selection in MC when you do PHOS trigger analysis.
  if(ForceActiveTRU) task->SetForceActiveTRU(L1input,L0input,Ethre,isMC);//this is to measure rejection factor from cluster energy kPHI7/kINT7 with same acceptance.

  task->SelectCollisionCandidates(trigger);

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

  //set minimum energy
  task->SetEmin(Emin);

  //centrality setting
  task->SetCentralityEstimator("HybridTrack");

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
    //TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) - ( 0 + [3]/(exp( -(x-[4]) / [5] ) + 1)  )",0,100);//this was not sufficient in LHC17pq
    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - [9]) - (  [3]/(exp( -(x-[4]) / [5] ) + 1) * [6]/(exp( -(x-[7]) / [8] ) + 1) )",0,100);
    f1tof->SetNpx(1000);
    f1tof->SetParameters( 3.04101e+00, 2.05831e+00, -6.79500e-01, 7.71770e-01, 6.43927e+00, 4.26270e-01, 7.20653e-01, 7.75214e+00, 1.07235e+00, 1.67197e+00);//20180314

    //f1tof->SetParameters(0.995,2.44,5.49e-2,0.770,6.19,0.298,0.719,7.90,0.917);//20180305
    //f1tof->SetParameters(0.992,2.47,5.80e-2,0.750,6.23,0.280,0.731,7.82,0.857);//20180302
    //f1tof->SetParameters(0.992,2.47,5.81e-2,0.523,7.71,0.586);//2018023
    //f1tof->SetParameters(0.990,2.52,6.10e-2,0.513,7.62,0.546);
    //f1tof->SetParameters(0.991,2.38,8.12e-3,0.425,7.63,0.677);
    //f1tof->SetParameters(0.996,2.33,1.97e-3,0.332,7.56,0.774);
    //f1tof->SetParameters(0.993,2.34,2.95e-3,0.367,7.27,0.556);
    task->SetTOFCutEfficiencyFunction(f1tof);
    //printf("TOF cut efficiency as a function of E is %s\n",f1tof->GetTitle());
  }
  if(!isMC && Trgcorrection){
    TF1 *f1trg = new TF1("f1TriggerEfficiency","[0]/(TMath::Exp(-(x-[1])/[2]) + 1)",0,100);
    f1trg->SetNpx(1000);
    f1trg->SetParameters(0.597,3.72,0.276);//from MB //acc x trigger efficiency 2-30GeV
    //f1trg->SetParameters(0.603,3.73,0.292);//from MB //acc x trigger efficiency 3-30GeV
    //f1trg->SetParameters(0.616,3.72,0.298);//from MB //acc x trigger efficiency 3-30GeV
    //f1trg->SetParameters(0.609,3.73,0.301);//from MB //acc x trigger efficiency 3-30GeV//old
    task->SetTriggerEfficiency(f1trg);
  }

  if(isMC){
    //for pi0 pT weighting
    const Int_t Ncen_Pi0 = 2;
    const Double_t centrality_Pi0[Ncen_Pi0] = {0,9999};
    TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0,centrality_Pi0);

    TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
    TF1 *f1weightPi0[Ncen_Pi0-1];

    if(MCtype.Contains("P8") || MCtype.Contains("Pythia8")){
      printf("Pythia8 is selected.\n");
      const Double_t p0_pi0[Ncen_Pi0-1] = {-0.58};
      const Double_t p1_pi0[Ncen_Pi0-1] = { 0.76};
      const Double_t p2_pi0[Ncen_Pi0-1] = { 1.14};

      for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
        f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"[2]*(1.+[0]/(1. + TMath::Power(x/[1],2)))",0,100);//this is iterative procedure.
        f1weightPi0[icen]->SetParameters(p0_pi0[icen],p1_pi0[icen],p2_pi0[icen]);
        farray_Pi0->Add(f1weightPi0[icen]);
      }
      task->SetAdditionalPi0PtWeightFunction(centarray_Pi0,farray_Pi0);

      //for eta/pi ratio
      const Int_t Ncen_Eta = 2;
      const Double_t centrality_Eta[Ncen_Eta] = {0,9999};
      TArrayD *centarray_Eta = new TArrayD(Ncen_Eta,centrality_Eta);

      TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      TF1 *f1weightEta[Ncen_Eta-1];
      const Double_t p0_eta[Ncen_Eta-1] = { 1.21};
      const Double_t p1_eta[Ncen_Eta-1] = {-0.499};
      const Double_t p2_eta[Ncen_Eta-1] = {1.44};

      for(Int_t icen=0;icen<Ncen_Eta-1;icen++){
        f1weightEta[icen] = new TF1(Form("f1weightEta_%d",icen),"[0]/(exp(-(x-[1])/[2]) + 1)",0,100);
        f1weightEta[icen]->SetParameters(p0_eta[icen],p1_eta[icen],p2_eta[icen]);
        farray_Eta->Add(f1weightEta[icen]);
      }

      task->SetAdditionalEtaPtWeightFunction(centarray_Eta,farray_Eta);

      //for K/pi ratio
      const Int_t Ncen_K0S = 2;
      const Double_t centrality_K0S[Ncen_K0S] = {0,9999};
      TArrayD *centarray_K0S = new TArrayD(Ncen_K0S,centrality_K0S);

      TObjArray *farray_K0S = new TObjArray(Ncen_K0S-1);
      TF1 *f1weightK0S[Ncen_K0S-1];
      const Double_t p0_K0S[Ncen_K0S-1] = { 1.37};
      const Double_t p1_K0S[Ncen_K0S-1] = { 4.98};
      const Double_t p2_K0S[Ncen_K0S-1] = {0.156};
      const Double_t p3_K0S[Ncen_K0S-1] = { 2.79};
      const Double_t p4_K0S[Ncen_K0S-1] = {0.239};

      for(Int_t icen=0;icen<Ncen_K0S-1;icen++){
        f1weightK0S[icen] = new TF1(Form("f1weightK0S_%d",icen),"[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1) )",0,100);
        f1weightK0S[icen]->SetParameters(p0_K0S[icen],p1_K0S[icen],p2_K0S[icen],p3_K0S[icen],p4_K0S[icen]);
        farray_K0S->Add(f1weightK0S[icen]);
      }

      task->SetAdditionalK0SPtWeightFunction(centarray_K0S,farray_K0S);

    }
    else if(MCtype.Contains("P6") || MCtype.Contains("Pythia6")){
      printf("Pythia6 is selected.\n");
      const Double_t p0_pi0[Ncen_Pi0-1] = {0.22};
      const Double_t p1_pi0[Ncen_Pi0-1] = { 2.2};
      const Double_t p2_pi0[Ncen_Pi0-1] = { 1.2};
      const Double_t p3_pi0[Ncen_Pi0-1] = { 0.8};

      for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
        f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"[0]*exp(-pow((x-[1]),2)/[2])+[3]",0,100);
        f1weightPi0[icen]->SetParameters(p0_pi0[icen],p1_pi0[icen],p2_pi0[icen],p3_pi0[icen]);
        farray_Pi0->Add(f1weightPi0[icen]);
      }
      task->SetAdditionalPi0PtWeightFunction(centarray_Pi0,farray_Pi0);


      //for eta/pi ratio
      const Int_t Ncen_Eta = 2;
      const Double_t centrality_Eta[Ncen_Eta] = {0,9999};
      TArrayD *centarray_Eta = new TArrayD(Ncen_Eta,centrality_Eta);

      TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      TF1 *f1weightEta[Ncen_Eta-1];
      const Double_t p0_eta[Ncen_Eta-1] = { 1.38};
      const Double_t p1_eta[Ncen_Eta-1] = {-0.320};
      const Double_t p2_eta[Ncen_Eta-1] = {1.52};

      for(Int_t icen=0;icen<Ncen_Eta-1;icen++){
        f1weightEta[icen] = new TF1(Form("f1weightEta_%d",icen),"[0]/(exp(-(x-[1])/[2]) + 1)",0,100);
        f1weightEta[icen]->SetParameters(p0_eta[icen],p1_eta[icen],p2_eta[icen]);
        farray_Eta->Add(f1weightEta[icen]);
      }

      task->SetAdditionalEtaPtWeightFunction(centarray_Eta,farray_Eta);


      //for K/pi ratio
      const Int_t Ncen_K0S = 2;
      const Double_t centrality_K0S[Ncen_K0S] = {0,9999};
      TArrayD *centarray_K0S = new TArrayD(Ncen_K0S,centrality_K0S);

      TObjArray *farray_K0S = new TObjArray(Ncen_K0S-1);
      TF1 *f1weightK0S[Ncen_K0S-1];
      const Double_t p0_K0S[Ncen_K0S-1] = { 1.44};
      const Double_t p1_K0S[Ncen_K0S-1] = { 5.76};

      for(Int_t icen=0;icen<Ncen_K0S-1;icen++){
        f1weightK0S[icen] = new TF1(Form("f1weightK0S_%d",icen),"[0] * (2/(1+exp(-[1]*x)) - 1)",0,100);
        f1weightK0S[icen]->SetParameters(p0_K0S[icen],p1_K0S[icen]);
        farray_K0S->Add(f1weightK0S[icen]);
      }

      task->SetAdditionalK0SPtWeightFunction(centarray_K0S,farray_K0S);
    }

  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString prefix = Form("hist_%s",taskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

