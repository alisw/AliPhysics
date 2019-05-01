AliAnalysisTaskPHOSSingleSim* AddTaskPHOSSingleSim(
    const char* name     = "SingleSim",
    const TString parname = "Pi0",
    UInt_t trigger = AliVEvent::kAny,
    const Bool_t isMC = kTRUE,
    const Int_t L1input = -1,//L1H,L1M,L1L
    const Int_t L0input = -1,//L0
    const Bool_t useCoreE = kFALSE,
    const Bool_t useCoreDisp = kFALSE,
    const Double_t NsigmaCPV  = 2.5,
    const Double_t NsigmaDisp = 2.5,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t NonLinStudy = kFALSE,
    const Double_t bs = 25.,//bunch space in ns.
    const Double_t distBC = 0,//minimum distance to bad channel.
    const Double_t Emin = 0.2,//minimum energy for photon selection in GeV
    const Int_t dummy = -1
    )
{
  //Add a task AliAnalysisTaskPHOSSingleSim to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSSingleSim", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSSingleSim", "This task requires an input event handler");
    return NULL;
  }

	TString TriggerName="";
	if     (trigger == (UInt_t)AliVEvent::kAny)  TriggerName = "kAny";
	else if(trigger == (UInt_t)AliVEvent::kINT7) TriggerName = "kINT7";
	else if(trigger == (UInt_t)AliVEvent::kPHI7) TriggerName = "kPHI7";

  if(trigger == (UInt_t)AliVEvent::kPHI7){
    if(L1input > 0){
      if(L1input == 7)      TriggerName = TriggerName + "_" + "L1H";
      else if(L1input == 6) TriggerName = TriggerName + "_" + "L1M";
      else if(L1input == 5) TriggerName = TriggerName + "_" + "L1L";
    }
    else if(L0input > 0)    TriggerName = TriggerName + "_" + "L0";
    else{
      ::Error("AddTaskPHOSSingleSim", "PHOS trigger analysis requires at least 1 trigger input (L0 or L1[H,M,L]).");
      return NULL;
    }
  }



  TString PIDname="";
  if(NsigmaCPV > 0) PIDname += Form("_CPV%d",(Int_t)(NsigmaCPV*10));
  if(NsigmaDisp > 0){
    if(useCoreDisp) PIDname += Form("_CoreDisp%d",(Int_t)(NsigmaDisp*10));
    else            PIDname += Form("_FullDisp%d",(Int_t)(NsigmaDisp*10));
  }
  if(useCoreE) PIDname += "_CoreE";
  else         PIDname += "_FullE";

  TString taskname = "";
  taskname = Form("%s_pp_%s%s_BS%dns_DBC%dcell_Emin%dMeV",name,TriggerName.Data(),PIDname.Data(),(Int_t)bs,(Int_t)(distBC),(Int_t)(Emin*1e+3));

  AliAnalysisTaskPHOSSingleSim* task = new AliAnalysisTaskPHOSSingleSim(taskname);

  Double_t Ethre = 0.0;

  //if(trigger == (UInt_t)AliVEvent::kPHI7) task->SetPHOSTriggerAnalysis(L1input,L0input,Ethre,isMC,kFALSE,dummy);
  if(trigger == (UInt_t)AliVEvent::kPHI7)       task->SetPHOSTriggerAnalysis(L1input,L0input,Ethre,isMC,kFALSE,dummy);
  else if(trigger  == (UInt_t) AliVEvent::kAny) task->SetPHOSTriggerAnalysisMB(L1input,L0input,Ethre,isMC,kFALSE,dummy);
  if(isMC && trigger == (UInt_t)AliVEvent::kPHI7) trigger = AliVEvent::kAny;//change trigger selection in MC when you do PHOS trigger analysis.

  //task->SelectCollisionCandidates(trigger);

  task->SetCollisionSystem(0);//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
  task->SetJetJetMC(kFALSE);
  task->SetMCType("MBMC");
  task->SetNonLinearityStudy(NonLinStudy,1.018);
  task->SetParticle(parname); 
  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
//  task->SetCoreEnergyFlag(useCoreE);

  const AliPHOSEventCuts::PileupFinder pf = AliPHOSEventCuts::kSPDInMultBins;
  task->SetEventCuts(isMC,pf);
  task->SetClusterCuts(useCoreDisp,NsigmaCPV,NsigmaDisp,useCoreE,distBC);

  task->SetCentralityMin(0);
  task->SetCentralityMax(9999);
  task->SetDepthNMixed(10);
  task->SetQnVectorTask(kFALSE);
  task->SetHarmonics(-1);

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

  printf("Your particle is %s\n",parname.Data());

  if(isMC){
    if(parname=="Pi0"){
      //for pi0 pT weighting
      const Int_t Ncen_Pi0 = 2;
      const Double_t centrality_Pi0[Ncen_Pi0] = {0,9999};
      TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0,centrality_Pi0);

      TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
      TF1 *f1weightPi0[Ncen_Pi0-1];

      const Double_t p0[Ncen_Pi0-1] = {2.70};
      const Double_t p1[Ncen_Pi0-1] = {0.132};
      const Double_t p2[Ncen_Pi0-1] = {6.64};

      for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
        f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"1.0 * ([0]/TMath::TwoPi() * ([2]-1)*([2]-2)/([2]*[1]*([2]*[1] + 0.139*([2]-2) )) * TMath::Power(1+(TMath::Sqrt(x*x+0.139*0.139) - 0.139)/([2]*[1]),-[2]))",0,100);//1/2pi x 1/Nev x 1/pT x d2N/dpTdy
        f1weightPi0[icen]->SetNpx(1000);
        f1weightPi0[icen]->SetParameters(p0[icen],p1[icen],p2[icen]);
        farray_Pi0->Add(f1weightPi0[icen]);
      }
      task->SetAdditionalPi0PtWeightFunction(centarray_Pi0,farray_Pi0);
    }
    else if(parname=="Eta"){
      //for eta pT weighting
      const Int_t Ncen_Eta = 2;
      const Double_t centrality_Eta[Ncen_Eta] = {0,9999};
      TArrayD *centarray_Eta = new TArrayD(Ncen_Eta,centrality_Eta);

      TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      TF1 *f1weightEta[Ncen_Eta-1];

      const Double_t p0[Ncen_Eta-1] = {0.25};
      const Double_t p1[Ncen_Eta-1] = {0.218};
      const Double_t p2[Ncen_Eta-1] = {6.60};

      for(Int_t icen=0;icen<Ncen_Eta-1;icen++){
        f1weightEta[icen] = new TF1(Form("f1weightEta_%d",icen),"[0]/TMath::TwoPi() * ([2]-1)*([2]-2)/([2]*[1]*([2]*[1] + 0.547*([2]-2) )) * TMath::Power(1+(TMath::Sqrt(x*x+0.547*0.547) - 0.547)/([2]*[1]),-[2])",0,100);//1/2pi x 1/Nev x 1/pT x d2N/dpTdy
        f1weightEta[icen]->SetNpx(1000);
        f1weightEta[icen]->SetParameters(p0[icen],p1[icen],p2[icen]);
        farray_Eta->Add(f1weightEta[icen]);
      }
      task->SetAdditionalEtaPtWeightFunction(centarray_Eta,farray_Eta);
    }
    else if(parname=="Gamma"){
      //for gamma pT weighting
      const Int_t Ncen_Gamma = 2;
      const Double_t centrality_Gamma[Ncen_Gamma] = {0,9999};
      TArrayD *centarray_Gamma = new TArrayD(Ncen_Gamma,centrality_Gamma);

      TObjArray *farray_Gamma = new TObjArray(Ncen_Gamma-1);
      TF1 *f1weightGamma[Ncen_Gamma-1];

      const Double_t p0[Ncen_Gamma-1] = {4.77586e+01};//Ae
      const Double_t p1[Ncen_Gamma-1] = {1.15413e-01};//Te
      const Double_t p2[Ncen_Gamma-1] = {9.08380e-01};//A
      const Double_t p3[Ncen_Gamma-1] = {4.76998e-01};//T
      const Double_t p4[Ncen_Gamma-1] = {3.04367e+00};//n

      for(Int_t icen=0;icen<Ncen_Gamma-1;icen++){
        f1weightGamma[icen] = new TF1(Form("f1weightGamma_%d",icen),"[0] * TMath::Exp(-(TMath::Sqrt(x*x + 0*0) - 0) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4])",0,100);//TCM fit to pp inclusive photon
        f1weightGamma[icen]->SetNpx(1000);
        f1weightGamma[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        farray_Gamma->Add(f1weightGamma[icen]);
      }
      task->SetAdditionalGammaPtWeightFunction(centarray_Gamma,farray_Gamma);
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

