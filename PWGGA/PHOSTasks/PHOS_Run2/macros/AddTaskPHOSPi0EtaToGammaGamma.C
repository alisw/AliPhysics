AliAnalysisTaskPHOSPi0EtaToGammaGamma* AddTaskPHOSPi0EtaToGammaGamma(
    const char* name     = "Pi0EtaToGammaGamma",
    const UInt_t trigger = AliVEvent::kINT7,
    const TString CollisionSystem = "pp",
    const Bool_t isMC = kFALSE,
    const Int_t L1input = -1,
    const Int_t L0input = -1,
    const Bool_t useCoreE = kFALSE,
    const Bool_t useCoreDisp = kFALSE,
    const Double_t NsigmaCPV  = 2.5,
    const Double_t NsigmaDisp = 2.5,
    const Bool_t usePHOSTender = kTRUE,
    const Double_t bs = 25.,//bunch space in ns.
    const Double_t distBC = -1,//minimum distance to bad channel.
    const Int_t pThardbin = -1
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

  const Bool_t rejectPileup = kTRUE;
  const Bool_t rejectDAQincomplete = kTRUE;
  const Double_t MaxAbsZ = 10.;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_Run2/macros/CreatePHOSEventCuts.C");
  AliPHOSTriggerHelper *helper = 0x0;

  if(trigger == (UInt_t)AliVEvent::kPHI7){
    if(L1input > -1 && L0input > -1){
      ::Error("AddTaskPHOSPi0EtaToGammaGamma", "Please select L0 or L1 trigger input.");
      return NULL;
    }

    if(L1input > -1){
      if     (L1input==7) TriggerName += "_L1H";
      else if(L1input==6) TriggerName += "_L1M";
      else if(L1input==5) TriggerName += "_L1L";
    }
    else if(L0input > -1){
      if(L0input==9)      TriggerName += "_L0";
    }
    else{
      ::Error("AddTaskPHOSPi0EtaToGammaGamma", "PHOS trigger analysis requires at least trigger input (L0 or L1).");
      return NULL;
    }

    helper = new AliPHOSTriggerHelper(L1input,L0input);
    helper->SetMatchingDistance(-3,-3,0,0);

    const TString pathBCMTRU = "alien:///alice/cern.ch/user/d/dsekihat/BadMap/TRUBadMap_LHC15n.root";
    TFile *rootfile_TRUmap = TFile::Open(pathBCMTRU,"READ");
    for(Int_t imod=1;imod<5;imod++){
      TH2I *h2 = (TH2I*)rootfile_TRUmap->Get(Form("hTRUBadMap_M%d",imod));
      helper->SetPHOSTRUBadMap(imod,h2);
    }


  }
  AliPHOSEventCuts *eventcuts = CreatePHOSEventCuts(kMC,helper);

  Int_t systemID = -1;
  if(CollisionSystem=="pp")                                 systemID = 0;
  else if(CollisionSystem=="PbPb")                          systemID = 1;
  else if(CollisionSystem=="pPb" || CollisionSystem=="Pbp") systemID = 2;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_Run2/macros/CreatePHOSClusterCuts.C");
  AliPHOSClusterCuts *clustercuts = CreatePHOSClusterCuts(useCoreDisp,NsigmaCPV,NsigmaDisp);

  TString PIDname="";
  if(NsigmaCPV > 0) PIDname += Form("_CPV%d",(Int_t)(NsigmaCPV*10));
  if(NsigmaDisp > 0){
    if(useCoreDisp) PIDname += Form("_CoreDisp%d",(Int_t)(NsigmaDisp*10));
    else            PIDname += Form("_FullDisp%d",(Int_t)(NsigmaDisp*10));
  }
  if(useCoreE) PIDname += "_CoreE";
  else         PIDname += "_FullE";

  TString taskname = Form("%s_%s_%s%s_BS%dns_DBC%02dmm",name,CollisionSystem.Data(),TriggerName.Data(),PIDname.Data(),(Int_t)bs,(Int_t)(distBC*10));

  AliAnalysisTaskPHOSPi0EtaToGammaGamma* task = new AliAnalysisTaskPHOSPi0EtaToGammaGamma(taskname);
  task->SelectCollisionCandidates(trigger);
  task->SetCollisionSystem(systemID);//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
  if(pThardbin>0) task->SetPtHardBin(pThardbin);

  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
  task->SetCoreEnergyFlag(useCoreE);

  task->SetEventCuts(eventcuts);
  task->SetClusterCuts(clustercuts);

  //centrality setting
  if(CollisionSystem=="pp"){//for pp
    const Int_t nbins = 0;
    Double_t cbin_data[nbins+1] = {0};//added lastbin + 1 as "Ntrack is greater than cbin[lastbin] upto infinity"
    Double_t cbin_MBMC[nbins+1] = {0};//added lastbin + 1 as "Ntrack is greater than cbin[lastbin] upto infinity"
    Double_t cbin_JJMC[nbins+1] = {0};//added lastbin + 1 as "Ntrack is greater than cbin[lastbin] upto infinity"

    TArrayD tbin;
    if(!isMC)           tbin = TArrayD(nbins+1, cbin_data);//data
    else{//for MC
      if(pThardbin < 0) tbin = TArrayD(nbins+1, cbin_MBMC);//MB MC to get same cluster occupance on PHOS
      else              tbin = TArrayD(nbins+1, cbin_JJMC);//MB MC to get same cluster occupance on PHOS
    }

    Int_t nMixed[nbins+1] = {100};
    TArrayI tNMixed(nbins+1, nMixed);
    task->SetCentralityBinningPP(tbin, tNMixed);
  }
  else if(CollisionSystem=="PbPb"){//for PbPb
    task->SetCentralityEstimator("V0M");

    const Int_t nbins = 9;
    Double_t cbin_data[nbins+1] = {0,10,20,30,40,50,60,70,80,90};//for data
    Double_t cbin_MBMC[nbins+1] = {1,11,21,31,41,51,61,71,81,91};//for MBMC to get same cluster occupance on PHOS
    Double_t cbin_JJMC[nbins+1] = {2,12,22,32,42,52,62,72,82,92};//for JJMC to get same cluster occupance on PHOS

    TArrayD tbin;
    if(!isMC)           tbin = TArrayD(nbins+1, cbin_data);//data
    else{//for MC
      if(pThardbin < 0) tbin = TArrayD(nbins+1, cbin_MBMC);//MB MC
      else              tbin = TArrayD(nbins+1, cbin_JJMC);//MB MC
    }

    Int_t nMixed[nbins] = {10,10,20,20,50,50,100,100,200};
    TArrayI tNMixed(nbins, nMixed);
    task->SetCentralityBinningPbPb(tbin, tNMixed);
  }
  else if(CollisionSystem=="pPb" || CollisionSystem=="Pbp"){//for pPb
    if(CollisionSystem == "pPb")      task->SetCentralityEstimator("V0A");//Pb-going side C->A
    else if(CollisionSystem == "Pbp") task->SetCentralityEstimator("V0C");//Pb-going side A->C

    const Int_t nbins = 5;
    Double_t cbin_data[nbins+1] = {0,20,40,60,80,100};//for data
    Double_t cbin_MBMC[nbins+1] = {0,20,40,60,80,100};//for MBMC to get same cluster occupance on PHOS
    Double_t cbin_JJMC[nbins+1] = {0,20,40,60,80,100};//for JJMC to get same cluster occupance on PHOS

    TArrayD tbin;
    if(!isMC)           tbin = TArrayD(nbins+1, cbin_data);//data
    else{//for MC
      if(pThardbin < 0) tbin = TArrayD(nbins+1, cbin_MBMC);//MB MC
      else              tbin = TArrayD(nbins+1, cbin_JJMC);//MB MC
    }

    Int_t nMixed[nbins] = {20,20,50,100,100};
    TArrayI tNMixed(nbins, nMixed);
    task->SetCentralityBinningPbPb(tbin, tNMixed);
  }

  //setting esd track selection for hybrid track
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
  AliESDtrackCuts *cutsG = CreateTrackCutsPWGJE(10001008);//for good global tracks
  task->SetESDtrackCutsForGlobal(cutsG);
  AliESDtrackCuts *cutsGC = CreateTrackCutsPWGJE(10011008);//for good global-constrained tracks
  task->SetESDtrackCutsForGlobalConstrained(cutsGC);

  //bunch space for TOF cut
  if(CollisionSystem=="pp"){//pp
    task->SetBunchSpace(25.);//in unit of ns.
    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) - ( 0 + [3]/(exp( -(x-[4]) / [5] ) + 1)  )",0,100);
    f1tof->SetParameters(0.996,2.45,0.039,0.31,7.16,0.620);

    //TF1 *f1tof = new TF1("f1TOFCutEfficiency","1.",0,100);

    task->SetTOFCutEfficiencyFunction(f1tof);
  }
  else if(CollisionSystem=="PbPb"){//PbPb
    task->SetBunchSpace(100.);//in unit of ns.
    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1)  )",0,100);
    f1tof->SetParameters(0.997,7.38,0.029,7.6,0.5);
    task->SetTOFCutEfficiencyFunction(f1tof);
  }
  else if(CollisionSystem=="pPb" || CollisionSystem=="pPb"){//pPb
    task->SetBunchSpace(100.);//in unit of ns.
    //TF1 *f1tof = new TF1("f1TOFCutEfficiency","1.",0,100);
    TF1 *f1tof = new TF1("f1TOFCutEfficiency","[0] * (2/(1+exp(-[1]*(x-[2]))) - 1) - ( 0 + [3]/(exp( -(x-[4]) / [5] ) + 1)  )",0,100);
    f1tof->SetParameters(0.996,6.79,0.0035,0.114,7.67,0.455);
    task->SetTOFCutEfficiencyFunction(f1tof);
  }

  if(isMC){
    //TF1 *f1Pi0Weight = new TF1("f1Pi0Weight","1+[0]*exp(-[1]*x)",0,100);
    //f1Pi0Weight->SetParameters(0.77,0.57);
    TF1 *f1Pi0Weight = new TF1("f1Pi0Weight","1.",0,100);
    task->SetAdditionalPi0PtWeightFunction(f1Pi0Weight);

    TF1 *f1K0SWeight = new TF1("f1K0SWeight","[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1)  )",0,100);
    f1K0SWeight->SetParameters(1.37,4.98,0.156,2.79,0.238);
    task->SetAdditionalK0SPtWeightFunction(f1K0SWeight);
  }


  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  TString outputFile = AliAnalysisManager::GetCommonFileName();

  TString prefix = Form("hist_%s",taskname.Data());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());//event characerization
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

