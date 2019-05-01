AliAnalysisTaskPHOSEmbeddingEfficiency* AddTaskPHOSEmbeddingEfficiency(
    const char* name     = "EmbeddingEfficiency",
    const TString parname = "Pi0",
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
    const TString MCtype = "MBMC",
    const Bool_t ForceActiveTRU = kFALSE,
    const Bool_t ApplyTOFTrigger = kFALSE,
    const AliPHOSEventCuts::PileupFinder pf = AliPHOSEventCuts::kMultiVertexer
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

  if(trigger & AliVEvent::kPHI7){
    if(L1input > 0){
      if(L1input == 7)      TriggerName = TriggerName + "_" + "L1H";
      else if(L1input == 6) TriggerName = TriggerName + "_" + "L1M";
      else if(L1input == 5) TriggerName = TriggerName + "_" + "L1L";
    }
    else if(L0input > 0)    TriggerName = TriggerName + "_" + "L0";
    else{
      ::Error("AddTaskPHOSEmbeddingEfficiency", "PHOS trigger analysis requires at least 1 trigger input (L0 or L1[H,M,L]).");
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
      ::Error("AddTaskPHOSEmbeddingEfficiency", "Qn flow vector correction is ON, but you do not set harmonics.");
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

    taskname = Form("%s_%s_%s_%s_Cen%d_%d%s_Harmonics%d_%s_%s_BS%dns_DBC%dcell_Emin%dMeV",name,parname.Data(),CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),harmonics,FMname.Data(),detname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));

  }
  else taskname = Form("%s_%s_%s_%s_Cen%d_%d%s_BS%dns_DBC%dcell_Emin%dMeV",name,parname.Data(),CollisionSystem.Data(),TriggerName.Data(),(Int_t)CenMin,(Int_t)CenMax,PIDname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));
  if((trigger & AliVEvent::kPHI7) && ApplyTOFTrigger) taskname += "_TOFTrigger";
  if(ForceActiveTRU) taskname += "_ForceActiveTRU";

  AliAnalysisTaskPHOSEmbeddingEfficiency* task = new AliAnalysisTaskPHOSEmbeddingEfficiency(taskname);

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

  task->SetEmbeddedParticle(parname);
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
    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) - ( 0 + [3]/(exp( -(x-[4]) / [5] ) + 1)  )",0,100);
    f1tof->SetParameters(0.996,5.61,-0.146,0.036,7.39,0.054);
    task->SetTOFCutEfficiencyFunction(f1tof);
    printf("TOF cut efficiency as a function of E is %s\n",f1tof->GetTitle());
  }

  if(isMC){
    if(parname.Contains("Pi0",TString::kIgnoreCase)){
      //for pi0 pT weighting
      const Int_t Ncen_Pi0 = 11;
      const Double_t centrality_Pi0[Ncen_Pi0] = {0,5,10,20,30,40,50,60,70,80,100};
      TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0,centrality_Pi0);

      TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
      TF1 *f1weightPi0[Ncen_Pi0-1];

      const Double_t p0[Ncen_Pi0-1] = {2.39991e+02, 1.78111e+02, 1.21109e+02, 6.91786e+01, 4.07880e+01, 2.05577e+01, 1.08079e+01, 4.98463e+00, 2.23119e+00, 1.16590e+00};//Ae
      const Double_t p1[Ncen_Pi0-1] = {3.84238e-01, 3.90424e-01, 3.96450e-01, 4.05978e-01, 4.11701e-01, 4.22538e-01, 4.23961e-01, 4.31410e-01, 4.28510e-01, 3.85667e-01};//Te
      const Double_t p2[Ncen_Pi0-1] = {1.16561e+03, 9.25084e+02, 7.15782e+02, 5.01254e+02, 3.55730e+02, 2.36175e+02, 1.49593e+02, 8.91887e+01, 4.56253e+01, 1.97146e+01};//A
      const Double_t p3[Ncen_Pi0-1] = {3.06202e-01, 3.12875e-01, 3.13599e-01, 3.14247e-01, 3.05660e-01, 3.01051e-01, 2.89748e-01, 2.75489e-01, 2.64481e-01, 2.50078e-01};//T
      const Double_t p4[Ncen_Pi0-1] = {2.73068e+00, 2.72595e+00, 2.70812e+00, 2.68939e+00, 2.64558e+00, 2.61370e+00, 2.56506e+00, 2.53088e+00, 2.51428e+00, 2.43636e+00};//n

      //printf("reading...alien:///alice/cern.ch/user/d/dsekihat/InputPtSpectra/InputPtSpectra_Embedding_PbPb_5.02TeV.root\n");
      //TFile *rootfile_pi0 = TFile::Open("alien:///alice/cern.ch/user/d/dsekihat/InputPtSpectra/InputPtSpectra_Embedding_PbPb_5.02TeV.root","READ");

      for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
        f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"[0] * TMath::Exp(-(TMath::Sqrt(x*x + 0.139*0.139) - 0.139) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4])",0,100);//TCM fit to PbPb
        f1weightPi0[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightPi0[icen]->SetNpx(1000);
        farray_Pi0->Add(f1weightPi0[icen]);

        //Int_t cen1 = centrality_Pi0[icen];
        //Int_t cen2 = centrality_Pi0[icen+1];
        //TF1 *f1weight = (TF1*)rootfile_pi0->Get(Form("f1weightPi0_Cen%d_%d",cen1,cen2));
        //farray_Pi0->Add(f1weight);
      }
      task->SetAdditionalPi0PtWeightFunction(centarray_Pi0,farray_Pi0);
      //rootfile_pi0->Close();
    }

    else if(parname.Contains("Eta",TString::kIgnoreCase)){
      //for eta pT weighting
      const Int_t Ncen_Eta = 11;
      const Double_t centrality_Eta[Ncen_Eta] = {0,5,10,20,30,40,50,60,70,80,100};
      TArrayD *centarray_Eta = new TArrayD(Ncen_Eta,centrality_Eta);

      TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      TF1 *f1weightEta[Ncen_Eta-1];

      const Double_t p0[Ncen_Eta-1] = {2.39991e+02, 1.78111e+02, 1.21109e+02, 6.91786e+01, 4.07880e+01, 2.05577e+01, 1.08079e+01, 4.98463e+00, 2.23119e+00, 1.16590e+00};//Ae
      const Double_t p1[Ncen_Eta-1] = {3.84238e-01, 3.90424e-01, 3.96450e-01, 4.05978e-01, 4.11701e-01, 4.22538e-01, 4.23961e-01, 4.31410e-01, 4.28510e-01, 3.85667e-01};//Te
      const Double_t p2[Ncen_Eta-1] = {1.16561e+03, 9.25084e+02, 7.15782e+02, 5.01254e+02, 3.55730e+02, 2.36175e+02, 1.49593e+02, 8.91887e+01, 4.56253e+01, 1.97146e+01};//A
      const Double_t p3[Ncen_Eta-1] = {3.06202e-01, 3.12875e-01, 3.13599e-01, 3.14247e-01, 3.05660e-01, 3.01051e-01, 2.89748e-01, 2.75489e-01, 2.64481e-01, 2.50078e-01};//T
      const Double_t p4[Ncen_Eta-1] = {2.73068e+00, 2.72595e+00, 2.70812e+00, 2.68939e+00, 2.64558e+00, 2.61370e+00, 2.56506e+00, 2.53088e+00, 2.51428e+00, 2.43636e+00};//n

      for(Int_t icen=0;icen<Ncen_Eta-1;icen++){
        f1weightEta[icen] = new TF1(Form("f1weightEta_%d",icen),"0.48 * ([0] * TMath::Exp(-(TMath::Sqrt(x*x + 0.547*0.547) - 0.547) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4]))",0,100);//TCM fit to PbPb
        f1weightEta[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightEta[icen]->SetNpx(1000);
        farray_Eta->Add(f1weightEta[icen]);
      }
      task->SetAdditionalEtaPtWeightFunction(centarray_Eta,farray_Eta);
    }

    else if(parname.Contains("Gamma",TString::kIgnoreCase)){
      //for gamma pT weighting
      const Int_t Ncen_Gamma = 11;
      const Double_t centrality_Gamma[Ncen_Gamma] = {0,5,10,20,30,40,50,60,70,80,100};
      TArrayD *centarray_Gamma = new TArrayD(Ncen_Gamma,centrality_Gamma);

      TObjArray *farray_Gamma = new TObjArray(Ncen_Gamma-1);
      TF1 *f1weightGamma[Ncen_Gamma-1];

      const Double_t p0[Ncen_Gamma-1] = {4.12097e+02 , 2.72720e+02, 1.85512e+02, 6.40594e+01, 6.40594e+01, 6.71162e+00, 6.71162e+00, 5.94487e+02, 5.94487e+02, 5.94487e+02};//Ae
      const Double_t p1[Ncen_Gamma-1] = {3.19655e-01 , 3.28095e-01, 3.31559e-01, 3.42014e-01, 3.42014e-01, 4.27787e-02, 4.27787e-02, 1.16145e-01, 1.16145e-01, 1.16145e-01};//Te
      const Double_t p2[Ncen_Gamma-1] = {5.26532e+03 , 3.86198e+03, 2.22781e+03, 8.59703e+02, 8.59703e+02, 1.74914e+02, 1.74914e+02, 1.32598e+01, 1.32598e+01, 1.32598e+01};//A
      const Double_t p3[Ncen_Gamma-1] = {2.25277e-01 , 2.38582e-01, 2.56551e-01, 2.93980e-01, 2.93980e-01, 3.75784e-01, 3.75784e-01, 4.84576e-01, 4.84576e-01, 4.84576e-01};//T
      const Double_t p4[Ncen_Gamma-1] = {2.76231e+00 , 2.79069e+00, 2.80401e+00, 2.86491e+00, 2.86491e+00, 3.02951e+00, 3.02951e+00, 3.09311e+00, 3.09311e+00, 3.09311e+00};//n

      for(Int_t icen=0;icen<Ncen_Gamma-1;icen++){
        f1weightGamma[icen] = new TF1(Form("f1weightGamma_%d",icen),"[0] * TMath::Exp(-(TMath::Sqrt(x*x + 0*0) - 0) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4])",0,100);//TCM fit to PbPb
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

