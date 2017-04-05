//
// Particle cuts
//
const Double_t etamin = -0.8;
const Double_t etamax =  0.8;
const Double_t ptmin = 0.0;
const Double_t ptmax = 30.0;
const Double_t phimin = 0.;
const Double_t phimax = 2*TMath::Pi();
const Double_t thetamin = 0;
const Double_t thetamax = TMath::Pi();
const Double_t zvtxmin = -10.0;
const Double_t zvtxmax =  10.0;
//
const Int_t mintrackrefsTPC = 5;
const Int_t mintrackrefsITS = 4;
const Int_t mintrackrefsTOF = 0;
const Int_t mintrackrefsMUON = 0;
const Bool_t ischarged = kTRUE;

//PID Threshold
const Float_t thresholdPID = 0.8;
AliCFSingleTrackEfficiencyTask *AddSingleTrackEfficiencyTaskPbPb(const Bool_t readAOD = 0, // Flag to read AOD:1 or ESD:0
								 TString suffix="default", // suffix for the output directory
								 AliPID::EParticleType specie=AliPID::kPion, 
								 Int_t pdgcode=0, //particle specie
								 Bool_t useMCtruthForKine=kFALSE,
								 ULong64_t triggerMask=AliVEvent::kAnyINT,
								 TString centralityEstimator = "V0M",
								 Int_t fBit=0,
								 Bool_t TPCRefit = kTRUE,
								 Int_t minclustersTPC = 0,
								 Bool_t ITSRefit = kTRUE,
								 Int_t spdHits=AliESDtrackCuts::kAny,
								 Int_t minclustersITS = 0,
								 Int_t configuration=AliCFSingleTrackEfficiencyTask::kFast,
								 Int_t usageOfBayesianPID=AliSingleTrackEffCuts::kNoBayesianPID)
{

  Info("AliCFSingleTrackEfficiencyTask","SETUP CONTAINER");

  //
  // Setting up the container
  // 
  // Variables
  const Int_t nvar = 7; // number of variables on the grid: pt, y, phi, theta, zvtx, multiplicity, centrality
  UInt_t nstep = 8;     // number of container steps
  const UInt_t ipt = 0;
  const UInt_t iy  = 1;
  const UInt_t iphi = 2;
  const UInt_t itheta = 3;
  const UInt_t izvtx = 4;
  const UInt_t imult = 5;
  const UInt_t icent = 6;
  //
  // Containter bining

  const Int_t nbinpt=24;
  //   A2. Bins variation by hand for other variables
  const Int_t nbin2 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 9 : 2; //bins in eta
  const Int_t nbin3 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 9 : 9; //bins in phi
  const Int_t nbin4 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 9 : 1; //bins in theta
  const Int_t nbin5 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 10 : 1; //bins in zvtx

  const Int_t nbinmult=8;
  const Int_t nbincent=12;

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbinpt;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
  iBin[3]=nbin4;
  iBin[4]=nbin5;
  iBin[5]=nbinmult;
  iBin[6]=nbincent;

  //arrays for lower bounds :
  Double_t binLimpT[nbinpt+1] = {0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,
				 2.,2.5,3.0,3.5,4.0,4.5,5.0,5.5,
				 6.0,7.0,8.0,10.,12.,14.,16.,20.,30.};
  Double_t *binLim2 = new Double_t[iBin[1]+1];
  Double_t *binLim3 = new Double_t[iBin[2]+1];
  Double_t *binLim4 = new Double_t[iBin[3]+1];
  Double_t *binLim5 = new Double_t[iBin[4]+1];
  Double_t binLimmult[nbinmult+1] = {0.,100.,500.,1000.,2000.,3000.,4000.,5000.,10000.};
  Double_t binLimcent[nbincent+1] = {0.,2.5,5.0,7.5,10.,15.,20.,30.,40.,50.,60.,80.,100.};

  // Other Variables
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)etamin + (etamax-etamin)/nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)phimin + (phimax-phimin)/nbin3*(Double_t)i ;
  for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)thetamin + (thetamax-thetamin)/nbin4*(Double_t)i ;
  for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)zvtxmin + (zvtxmax-zvtxmin)/nbin5*(Double_t)i ;

  // Container  
  AliCFContainer* container = new AliCFContainer(Form("container%s",suffix.Data()),"container for tracks",nstep,nvar,iBin);
  container -> SetBinLimits(ipt,binLimpT);    // pt
  container -> SetBinLimits(iy,binLim2);      // eta
  container -> SetBinLimits(iphi,binLim3);    // phi
  container -> SetBinLimits(itheta,binLim4);  // theta
  container -> SetBinLimits(izvtx,binLim5);   // Zvtx
  container -> SetBinLimits(imult,binLimmult);// multiplicity
  container -> SetBinLimits(icent,binLimcent);// centrality

  // Variable Titles
  container -> SetVarTitle(ipt,"pt");
  container -> SetVarTitle(iy, "y");
  container -> SetVarTitle(iphi,"phi");
  container -> SetVarTitle(itheta, "theta");
  container -> SetVarTitle(izvtx, "Zvtx");
  container -> SetVarTitle(imult, "Multiplicity");
  container -> SetVarTitle(icent, "Centrality");

  // Step Titles
  container -> SetStepTitle(0, " MC Particle with Generated Cuts");
  container -> SetStepTitle(1, " MC Particle with Kine Acceptance Cuts");
  container -> SetStepTitle(2, " MC Particle with Track Ref Acceptance Cuts");
  container -> SetStepTitle(3, " Total Reconstructed  Particle ");
  container -> SetStepTitle(4, " Reco Particle With Kine Acceptance Cuts");
  container -> SetStepTitle(5, " Reco Particle to MC True pt particles ");
  container -> SetStepTitle(6, " Reco Particle With Quality Cuts");
  container -> SetStepTitle(7, " Reco PID With Quality Cuts");


  // SET TLIST FOR QA HISTOS
  TList* qaList = new TList();
  TObjArray* emptyList = new TObjArray(0);

  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager();

  man->SetNStepEvent(2);
  man->SetEventContainer(container);
  man->SetEventCutsList(0,emptyList);//evtmcList);
  man->SetEventCutsList(1,emptyList);//evtrecoList);
  
  man->SetParticleContainer(container);
  man->SetParticleCutsList(0,emptyList);//mcGenList);
  man->SetParticleCutsList(1,emptyList);//mcKineList);
  man->SetParticleCutsList(2,emptyList);//mcaccList);
  man->SetParticleCutsList(3,emptyList);//evtrecoPureList);
  man->SetParticleCutsList(4,emptyList);//recKineList);
  man->SetParticleCutsList(5,emptyList);
  man->SetParticleCutsList(6,emptyList);
  man->SetParticleCutsList(7,emptyList);
  
  // Simulated particle & event cuts
  AliSingleTrackEffCuts* cuts = new AliSingleTrackEffCuts();
  cuts->SetPtRange(ptmin,ptmax);
  cuts->SetEtaRange(etamin,etamax);
  cuts->SetIsCharged(ischarged);
  cuts->SetMinVtxContr(1);
  cuts->SetMaxVtxZ(zvtxmax);
  cuts->SetNumberOfClusters(mintrackrefsITS,mintrackrefsTPC,mintrackrefsTOF,mintrackrefsMUON);
  cuts->SetTriggerMask(triggerMask);
  cuts->SetIsAOD(readAOD);
  //
  // Pid selection here
  //
  if(pdgcode>0){
    cuts->SetUsePid(true);
    cuts->SetParticleSpecie(specie);
    cuts->SetPdgCode(pdgcode);
    // 
    const Int_t nlims=1;
    Float_t plims[nlims+1]={0.,999.}; //TPC limits in momentum [GeV/c]
    Float_t sigmas[nlims]={3.};
    cuts->SetUseTPCPid();
    cuts->SetTPCSigmaPtBins(nlims,plims,sigmas);
    cuts->SetMaximumPTPC(4.);
    // 
    const Int_t nlims2=1;
    Float_t plims2[nlims2+1]={0.,999.}; //TPC limits in momentum [GeV/c]
    Float_t sigmas2[nlims2]={3.};
    cuts->SetUseTOFPid();
    cuts->SetTOFSigmaPtBins(nlims2,plims2,sigmas2);
    cuts->SetMaximumPTOF(4.);

    if(usageOfBayesianPID>0) {
      cuts->SetUseCombinPID(usageOfBayesianPID);
      if(usageOfBayesianPID==AliSingleTrackEffCuts::kThresholdBayesianProb)
	cuts->SetPIDThreshold(thresholdPID);
    }

  }

  //
  //  Track Quality cuts via ESD track cuts
  //
  AliESDtrackCuts* QualityCuts = new AliESDtrackCuts();
  QualityCuts->SetRequireSigmaToVertex(kFALSE);
  QualityCuts->SetMinNClustersTPC(minclustersTPC);
  QualityCuts->SetMinNClustersITS(minclustersITS);
  QualityCuts->SetRequireTPCRefit(TPCRefit);
  QualityCuts->SetRequireITSRefit(ITSRefit);
  QualityCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,spdHits);
  QualityCuts->SetMinDCAToVertexXY(0.);
  QualityCuts->SetEtaRange(etamin,etamax);
  QualityCuts->SetPtRange(ptmin,ptmax);


  //CREATE THE TASK
  printf("CREATE CF Single track task\n");

  AliCFSingleTrackEfficiencyTask *task = new AliCFSingleTrackEfficiencyTask("AliCFSingleTrackEfficiencyTask",QualityCuts,cuts);
  if(readAOD && fBit>=0){
    task->SetFilterBit(kTRUE);
    task->SetFilterType(fBit);
  }else{
    task->SetFilterBit(kFALSE);
  }
  //  task->SelectCollisionCandidates(triggerMask);//AliVEvent::kMB);
  if(centralityEstimator != "") task->SetUseCentrality(kTRUE,centralityEstimator);
  task->SetConfiguration(configuration);
  task->SetUseGeneratedKine(useMCtruthForKine);
  task->SetCFManager(man); //here is set the CF manager

  //
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }

  // This task requires an ESD or AOD input handler and an AOD output handler.
  // Check this using the analysis manager.
  //===============================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddSingleTrackEfficiencyTaskPbPb", "AliCFSingleTrackEfficiency task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }
  
  mgr->AddTask(task);
  printf(" Create the output container\n");

  //
  // Create and connect containers for input/output
  //
  // ----- output data -----
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString input1name="cchain0";
  TString output2name="HistEventsProcessed", output3name="container",output4name="list",output5name="ESDtrackCuts",output6name="MCtrackCuts";
  outputfile += ":PWGPP_CFSingleTrack";
//  outputfile += suffix;
  output2name += suffix;
  output3name += suffix;
  output4name += suffix;
  output5name += suffix;
  output6name += suffix;


  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  // ----- output data -----
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output2name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output3name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output QA histograms
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output4name, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output ESD track cuts for book keeping
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output5name, AliESDtrackCuts::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output event and particle selection cuts for book keeping
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(output6name, AliSingleTrackEffCuts::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);
  mgr->ConnectOutput(task,5,coutput5);

  return task;
}
