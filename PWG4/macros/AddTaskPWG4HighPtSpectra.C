//DEFINITION OF A FEW CONSTANTS
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

const Float_t phimin = 0.;
const Float_t phimax = 2.*TMath::Pi();
const Float_t etamin = -0.9;
const Float_t etamax = 0.9;

const Int_t   mintrackrefsTPC = 0;
const Int_t   mintrackrefsITS = 0;

void AddTaskPWG4HighPtSpectraAll(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0)
{
  int cent = 10;

  AliPWG4HighPtSpectra *taskSpectra00cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,0);
  AliPWG4HighPtSpectra *taskSpectra01cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,1);
  //  AliPWG4HighPtSpectra *taskSpectra02cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,2);
  //  AliPWG4HighPtSpectra *taskSpectra10cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,1,0);
  //  AliPWG4HighPtSpectra *taskSpectra20cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,2,0);
  AliPWG4HighPtSpectra *taskSpectra70cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,7,0);
  AliPWG4HighPtSpectra *taskSpectra71cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,7,1);
  AliPWG4HighPtSpectra *taskSpectra72cent10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,7,2);

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtSpectra *taskSpectra00 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,0);
      AliPWG4HighPtSpectra *taskSpectra01 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,1);
      //      AliPWG4HighPtSpectra *taskSpectra02 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,0,2);
      //      AliPWG4HighPtSpectra *taskSpectra10 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,1,0);
      //      AliPWG4HighPtSpectra *taskSpectra20 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,2,0);
      AliPWG4HighPtSpectra *taskSpectra70 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,7,0);
      AliPWG4HighPtSpectra *taskSpectra71 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,7,1);
      AliPWG4HighPtSpectra *taskSpectra72 = AddTaskPWG4HighPtSpectra(prodType,isPbPb,cent,7,2);
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
  
  //Load common track cut class
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/CreateTrackCutsPWG4.C");

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
  UInt_t incls = 3;

  //Setting up the container grid... 
  UInt_t nstep = 4; //Steps/Modes for containers
  Int_t kStepReconstructed          = 0;
  Int_t kStepSecondaries            = 1;
  Int_t kStepReconstructedMC        = 2;
  Int_t kStepMCAcceptance           = 3;

  //redefine pt ranges in case of Jet-Jet production
  Float_t ptBinEdges[4][4];
  Float_t ptmin =  0.0 ;
  Float_t ptmax =  100.0 ;

  ptBinEdges[0][0] = 2.;
  ptBinEdges[0][1] = 0.2;
  ptBinEdges[1][0] = 6.;
  ptBinEdges[1][1] = 0.5;
  ptBinEdges[2][0] = 20.;
  ptBinEdges[2][1] = 2.;
  ptBinEdges[3][0] = 100.;
  ptBinEdges[3][1] = 5.;
  
  const Int_t nvar   = 4; //number of variables on the grid: pt:phi:eta:NClsIter1

  const Int_t nbin11 = round((ptBinEdges[0][0]-ptmin)/ptBinEdges[0][1]);
  const Int_t nbin12 = round((ptBinEdges[1][0]-ptBinEdges[0][0])/ptBinEdges[1][1])+nbin11;
  const Int_t nbin13 = round((ptBinEdges[2][0]-ptBinEdges[1][0])/ptBinEdges[2][1])+nbin12;
  const Int_t nbin14 = round((ptBinEdges[3][0]-ptBinEdges[2][0])/ptBinEdges[3][1])+nbin13;

  const Int_t nbin1  = nbin14; //bins in pt 
  const Int_t nbin2  =  18;    //bins in phi
  const Int_t nbin3  =  2;     //bins in eta
  const Int_t nbin4  =  6;     //bins in NClsIter1: 0 70 80 90 100 120


  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
  iBin[3]=nbin4;
   
  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];
  Double_t *binLim3=new Double_t[nbin3+1];
  Double_t *binLim4=new Double_t[nbin4+1];
  
  //values for bin lower bounds 
  for(Int_t i=0; i<=nbin1; i++) {
    if(i<=nbin11) binLim1[i]=(Double_t)ptmin + (ptBinEdges[0][0]-ptmin)/nbin11*(Double_t)i ;  
    if(i<=nbin12 && i>nbin11) binLim1[i]=(Double_t)ptBinEdges[0][0] + (ptBinEdges[1][0]-ptBinEdges[0][0])/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ; 
    if(i<=nbin13 && i>nbin12) binLim1[i]=(Double_t)ptBinEdges[1][0] + (ptBinEdges[2][0]-ptBinEdges[1][0])/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ; 
    if(i<=nbin14 && i>nbin13) binLim1[i]=(Double_t)ptBinEdges[2][0] + (ptBinEdges[3][0]-ptBinEdges[2][0])/(nbin14-nbin13)*((Double_t)i-(Double_t)nbin13) ; 
  }
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)phimin + (phimax-phimin)/nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)etamin + (etamax-etamin)/nbin3*(Double_t)i ;  
  binLim4[0] = 0.;
  binLim4[1] = 70.;
  binLim4[2] = 80.;
  binLim4[3] = 90.;
  binLim4[4] = 100.;
  binLim4[5] = 120.;
  binLim4[6] = 160.;


  AliCFContainer* containerPos = new AliCFContainer("containerPos","container for positive tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerPos -> SetBinLimits(ipt,binLim1);
  containerPos -> SetBinLimits(iphi,binLim2);
  containerPos -> SetBinLimits(ieta,binLim3);
  containerPos -> SetBinLimits(incls,binLim4);

  AliCFContainer* containerNeg = new AliCFContainer("containerNeg","container for negative tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerNeg -> SetBinLimits(ipt,binLim1);
  containerNeg -> SetBinLimits(iphi,binLim2);
  containerNeg -> SetBinLimits(ieta,binLim3);
  containerNeg -> SetBinLimits(incls,binLim4);
  
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  AliESDtrackCuts *trackCutsReject = 0x0;
  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0 && cuts==0) {
    // tight global tracks - RAA analysis
    trackCuts = CreateTrackCutsPWG4(1000);
  }
  if(trackType==0 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis
    trackCuts = CreateTrackCutsPWG4(10001005);
  }
  if(trackType==0 && cuts==2) {
    //Cuts global tracks with ITSrefit requirement but without SPD
    trackCuts = CreateTrackCutsPWG4(10011005);
  }
  if(trackType==7 && cuts==0) {
    // no requirements on SPD and ITSrefit failed
    trackCuts = CreateTrackCutsPWG4(10041005);   //no ITSrefit requirement
    trackCutsReject = CreateTrackCutsPWG4(1005); //ITSrefit requirement
    trackCutsReject->SetEtaRange(etamin,etamax);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }
  if(trackType==7 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement but without SPD
    trackCuts = CreateTrackCutsPWG4(10011005);
  }
  if(trackType==7 && cuts==2) {
    // no requirements on SPD and ITSrefit failed
    trackCuts = CreateTrackCutsPWG4(10041005);       //no ITSrefit requirement filter 256
    trackCutsReject = CreateTrackCutsPWG4(10001005); //ITSrefit requirement filter 16
    trackCutsReject->SetEtaRange(etamin,etamax);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }

  if(trackType==1 && cuts==0) {
    //Set track cuts for TPConly tracks
    trackCuts = CreateTrackCutsPWG4(2001);
  }
  if(trackType==2 && cuts==0) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWG4(2001);
  }
  trackCuts->SetEtaRange(etamin,etamax);
  trackCuts->SetPtRange(0.15, 1e10);

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(0.15,1e10);
  mcKineCuts->SetEtaRange(etamin,etamax);//-0.5,0.5);
  mcKineCuts->SetRequireIsCharged(kTRUE);

  //Acceptance Cuts
  //AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  // mcAccCuts->SetMinNHitITS(mintrackrefsITS);
  //mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);

  TObjArray* recList = new TObjArray(0);
  TObjArray* secList = new TObjArray(0) ;
  TObjArray* recMCList = new TObjArray(0);

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  //mcList->AddLast(mcAccCuts);

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


  AliPWG4HighPtSpectra *taskPWG4HighPtSpectra = new AliPWG4HighPtSpectra(Form("AliPWG4HighPtSpectraCent%dTrackType%dCuts%d",centClass,trackType,cuts));
  taskPWG4HighPtSpectra->SetTrackType(trackType);
  taskPWG4HighPtSpectra->SetCuts(trackCuts);
  taskPWG4HighPtSpectra->SetCutsReject(trackCutsReject);
  taskPWG4HighPtSpectra->SetCFManagerPos(manPos); //here is set the CF manager +
  taskPWG4HighPtSpectra->SetCFManagerNeg(manNeg); //here is set the CF manager -

  if(isPbPb) {
    taskPWG4HighPtSpectra->SetIsPbPb(kTRUE);
    taskPWG4HighPtSpectra->SetCentralityClass(centClass);
  }
  //  taskPWG4HighPtSpectra->SetSigmaConstrainedMax(5.);


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
