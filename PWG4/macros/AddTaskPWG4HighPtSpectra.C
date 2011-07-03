//DEFINITION OF A FEW CONSTANTS
const Float_t phimin = 0.;
const Float_t phimax = 2.*TMath::Pi();
const Float_t etamin = -0.9;
const Float_t etamax = 0.9;

const Int_t   mintrackrefsTPC = 1;
const Int_t   mintrackrefsITS = 1;

void AddTaskPWG4HighPtSpectraAll(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0)
{
  int cent = 10;

  AliPWG4HighPtSpectra *taskSpectra00cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,0);
  AliPWG4HighPtSpectra *taskSpectra01cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,1);
  AliPWG4HighPtSpectra *taskSpectra02cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,2);
  AliPWG4HighPtSpectra *taskSpectra03cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,3);
  //  AliPWG4HighPtSpectra *taskSpectra10cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,1,0);
  AliPWG4HighPtSpectra *taskSpectra20cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,2,0);

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtSpectra *taskSpectra00 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,0);
      AliPWG4HighPtSpectra *taskSpectra01 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,1);
      AliPWG4HighPtSpectra *taskSpectra02 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,2);
      AliPWG4HighPtSpectra *taskSpectra03 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,3);
      //      AliPWG4HighPtSpectra *taskSpectra10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,1,0);
      AliPWG4HighPtSpectra *taskSpectra20 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,2,0);
    }
  }

}


AliPWG4HighPtSpectra* AddTaskPWG4HighPtSpectra(char *prodType = "LHC10e14", Bool_t isPbPb=kTRUE,Int_t centClass = 0, Int_t trackType = 0, Int_t cuts = 0)
{

  /*
    trackType: 0 = global
               1 = TPC stand alone
               2 = TPC stand alone constrained to SPD vertex
    cuts:      0 (global) = standard ITSTPC2010 a la RAA analysis
               1 (global) = ITSrefit, no SPD requirements -> standard for jet analysis
               2 (global) = ITSrefit + no hits in SPD
	       3 (global) = standard ITS tight cuts with nCrossed rows cut for hybrid tracks
               0 (TPC)    = standard TPC + NClusters>70
               1 (TPC)    = standard TPC + NClusters>0 --> to study new TPC QA recommendations
   */

  // Creates HighPtSpectra analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtSpectra", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPWG4HighPtSpectra", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  const char *analysisType = "ESD";//"TPC"

  // C. Create the task, add it to manager.
  //===========================================================================
 //CONTAINER DEFINITION
  Info("AliPWG4HighPtSpectra","SETUP CONTAINER");
  //the sensitive variables, their indices
  UInt_t ipt   = 0;
  UInt_t iphi  = 1;
  UInt_t ieta  = 2;

  //Setting up the container grid... 
  UInt_t nstep = 7; //Steps/Modes for containers
  Int_t kStepReconstructed          = 0;
  Int_t kStepSecondaries            = 1;
  Int_t kStepReconstructedMC        = 2;
  Int_t kStepMCAcceptance           = 3;

  //redefine pt ranges in case of Jet-Jet production
  Float_t ptBinEdges[2][2];
  Float_t ptmin =  2.0 ;
  Float_t ptmax =  50.0 ;
  Float_t binWidth3 = 5.;
  if(!strcmp(prodType, "LHC10e14")) { 
    ptmin =  0.0 ;
    ptmax =  500.0 ;

    ptBinEdges[0][0] = 100.;
    ptBinEdges[0][1] = 5.;
    ptBinEdges[1][0] = 300.;
    ptBinEdges[1][1] = 10.;
    binWidth3 = 20.;
  } else { 
    ptmin =  2.0 ;
    ptmax =  50.0 ;
    
    ptBinEdges[0][0] = 10.;
    ptBinEdges[0][1] = 0.5;
    ptBinEdges[1][0] = 60.;
    ptBinEdges[1][1] = 2.;
    binWidth3 = 5.;
  }
    
  const Int_t nvar   = 3; //number of variables on the grid: pt:phi:eta
  const Int_t nbin11 = (int)((ptBinEdges[0][0]-ptmin)/ptBinEdges[0][1]);
  const Int_t nbin12 = (int)((ptBinEdges[1][0]-ptBinEdges[0][0])/ptBinEdges[1][1])+nbin11;
  const Int_t nbin13 = (int)((ptmax-ptBinEdges[1][0])/binWidth3)+nbin12;
  const Int_t nbin1  = nbin13; //bins in pt 
  const Int_t nbin2  =  18;    //bins in phi
  const Int_t nbin3  =  2;     //bins in eta

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
   
  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];
  Double_t *binLim3=new Double_t[nbin3+1];
  
  //values for bin lower bounds 
  for(Int_t i=0; i<=nbin1; i++) {
    if(i<=nbin11) binLim1[i]=(Double_t)ptmin + (ptBinEdges[0][0]-ptmin)/nbin11*(Double_t)i ;  
    if(i<=nbin12 && i>nbin11) binLim1[i]=(Double_t)ptBinEdges[0][0] + (ptBinEdges[1][0]-ptBinEdges[0][0])/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;  
    if(i<=nbin13 && i>nbin12) binLim1[i]=(Double_t)ptBinEdges[1][0] + (ptmax-ptBinEdges[1][0])/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;  
  }
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)phimin + (phimax-phimin)/nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)etamin + (etamax-etamin)/nbin3*(Double_t)i ;  
  

  AliCFContainer* containerPos = new AliCFContainer("containerPos","container for positive tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerPos -> SetBinLimits(ipt,binLim1);
  containerPos -> SetBinLimits(iphi,binLim2);
  containerPos -> SetBinLimits(ieta,binLim3);

  AliCFContainer* containerNeg = new AliCFContainer("containerNeg","container for negative tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerNeg -> SetBinLimits(ipt,binLim1);
  containerNeg -> SetBinLimits(iphi,binLim2);
  containerNeg -> SetBinLimits(ieta,binLim3);
  
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
 AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0 && cuts==0) {
    // tight global tracks
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);// essentially swittches it off
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
  }
  if(trackType==0 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement for jet analysis
    // TPC  
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);    
    
    trackCuts->SetRequireSigmaToVertex(kFALSE);
    
  }
  if(trackType==0 && cuts==2) {
    //Cuts global tracks with ITSrefit requirement but without SPD
    // TPC  
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);    
    trackCuts->SetRequireSigmaToVertex(kFALSE);
    //no SPD points
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  }
  if(trackType==0 && cuts==3) {
    // tight global tracks
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,1);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);// essentially swittches it off
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
  }
  if(trackType==1 && cuts==0) {
    //Set track cuts for TPConly tracks
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts(); 
    trackCuts->SetMinNClustersTPC(70);
  }
  if(trackType==2 && cuts==0) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts();
    trackCuts->SetMinNClustersTPC(70);
  }
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(0.15,1e10);
  mcKineCuts->SetRapidityRange(-0.9,0.9);//-0.5,0.5);
  mcKineCuts->SetRequireIsCharged(kTRUE);

  //Acceptance Cuts
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  // mcAccCuts->SetMinNHitITS(mintrackrefsITS);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);

  TObjArray* recList = new TObjArray(0);
  TObjArray* secList = new TObjArray(0) ;
  TObjArray* recMCList = new TObjArray(0);

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcAccCuts);

  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* manPos = new AliCFManager("manPos","Manager for Positive tracks") ;
  manPos->SetParticleContainer(containerPos);
  manPos->SetParticleCutsList(kStepReconstructed,recList);
  manPos->SetParticleCutsList(kStepSecondaries,secList);
  manPos->SetParticleCutsList(kStepReconstructedMC,recMCList);
  manPos->SetParticleCutsList(kStepMCAcceptance,mcList);

  AliCFManager* manNeg = new AliCFManager("manNeg","Manager for Negative tracks") ;
  manNeg->SetParticleContainer(containerNeg);
  manNeg->SetParticleCutsList(kStepReconstructed,recList);
  manNeg->SetParticleCutsList(kStepSecondaries,secList);
  manNeg->SetParticleCutsList(kStepReconstructedMC,recMCList);
  manNeg->SetParticleCutsList(kStepMCAcceptance,mcList);


  printf("Create task AliPWG4HighPtSpectra\n");
  AliPWG4HighPtSpectra *taskPWG4HighPtSpectra = new AliPWG4HighPtSpectra(Form("AliPWG4HighPtSpectraCent%dTrackType%dCuts%d",centClass,trackType,cuts));
  taskPWG4HighPtSpectra->SetTrackType(trackType);
  taskPWG4HighPtSpectra->SetCuts(trackCuts);
  taskPWG4HighPtSpectra->SetCFManagerPos(manPos); //here is set the CF manager +
  taskPWG4HighPtSpectra->SetCFManagerNeg(manNeg); //here is set the CF manager -

  if(isPbPb) {
    taskPWG4HighPtSpectra->SetIsPbPb(kTRUE);
    taskPWG4HighPtSpectra->SetCentralityClass(centClass);
  }
  taskPWG4HighPtSpectra->SetSigmaConstrainedMax(5.);


  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ output containers ------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtSpectraCent%dTrackType%dCuts%d",centClass,trackType,cuts);

  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer(Form("chist0HighPtSpectraCent%dTrackType%dCuts%d",centClass,trackType,cuts), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("ccontainer0HighPtSpectraCent%dTrackType%dCuts%d",centClass,trackType,cuts), AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("ccontainer1HighPtSpectraCent%dTrackType%dCuts%d",centClass,trackType,cuts), AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *cout_cuts0 = mgr->CreateContainer(Form("qa_SpectraTrackCutsCent%dTrackType%dCuts%d",centClass,trackType,cuts), AliESDtrackCuts::Class(), AliAnalysisManager::kParamContainer,outputfile);

  mgr->AddTask(taskPWG4HighPtSpectra);

  mgr->ConnectInput(taskPWG4HighPtSpectra,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4HighPtSpectra,0,coutput0);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,1,coutput1);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,2,coutput2);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,3,cout_cuts0);

  // Return task pointer at the end
  return taskPWG4HighPtSpectra;
}
