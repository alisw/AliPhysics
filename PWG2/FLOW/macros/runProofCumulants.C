//RUN SETTINGS
//analysis type can be ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";


//SETTING THE CUTS

//for integrated flow
const Double_t ptmin1 = 0.0;
const Double_t ptmax1 = 1000.0;
const Double_t ymin1  = -2.;
const Double_t ymax1  = 2.;
const Int_t mintrackrefsTPC1 = 2;
const Int_t mintrackrefsITS1 = 3;
const Int_t charge1 = 1;
const Int_t PDG1 = 211;
const Int_t minclustersTPC1 = 50;
const Int_t maxnsigmatovertex1 = 3;

//for differential flow
const Double_t ptmin2 = 0.0;
const Double_t ptmax2 = 1000.0;
const Double_t ymin2  = -2.;
const Double_t ymax2  = 2.;
const Int_t mintrackrefsTPC2 = 2;
const Int_t mintrackrefsITS2 = 3;
const Int_t charge2 = 1;
const Int_t PDG2 = 321;
const Int_t minclustersTPC2 = 50;
const Int_t maxnsigmatovertex2 = 3;


  void runProofCumulants(const Char_t* data="/PWG2/akisiel/LHC500C2030", Int_t nRuns=4, Int_t offset=0) {
  //void runProofCumulants(const Char_t* data="/PWG2/akisiel/LHC500C0005", Int_t nRuns=-1, Int_t offset=0){ 

  //void runProofCumulants(const Char_t* data="/PWG2/pganoti/Pythia6At10TeV_05T_b", Int_t nRuns=1000, Int_t offset=0){ 
  //void runProofCumulants(const Char_t* data="/PWG2/jgrosseo/sim_1600XX_esd", Int_t nRuns=-1, Int_t offset=0){ 
  
  //void runProofCumulants(const Char_t* data="/COMMON/COMMON/run15035_PbPb", Int_t nRuns=-1, Int_t offset=0){

  //void runProofCumulants(const Char_t* data="/PWG2/hricaud/LHC07f_160033DataSet", Int_t nRuns=-1, Int_t offset=0){

  //void runProofCumulants(const Char_t* data="/PWG2/hricaud/LHC07f_160038_root_archiveDataSet", Int_t nRuns=-1, Int_t offset=0){

 //void runProofCumulants(const Char_t* data="/PWG2/belikov/40825", Int_t nRuns=-1, Int_t offset=0){
 


  TStopwatch timer;
  timer.Start();

  printf("*** Connect to PROOF ***\n");
  TProof::Open("abilandz@lxb6046.cern.ch");
  //TProof::Open("snelling@localhost");

 gProof->UploadPackage("STEERBase.par");
 gProof->EnablePackage("STEERBase");
 gProof->UploadPackage("ESD.par");
 gProof->EnablePackage("ESD");
 gProof->UploadPackage("AOD.par");
 gProof->EnablePackage("AOD");
 gProof->UploadPackage("ANALYSIS.par");
 gProof->EnablePackage("ANALYSIS");
 gProof->UploadPackage("ANALYSISalice.par");
 gProof->EnablePackage("ANALYSISalice");
 gProof->UploadPackage("PWG2AOD.par");
 gProof->EnablePackage("PWG2AOD");
 gProof->UploadPackage("CORRFW.par");
 gProof->EnablePackage("CORRFW");
 gProof->ClearPackage("PWG2flow");
 gProof->UploadPackage("PWG2flow.par");
 gProof->EnablePackage("PWG2flow");



//____________________________________________//
 //Create cuts using correction framework

 //############# cuts on MC
 AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
 mcKineCuts1->SetPtRange(ptmin1,ptmax1);
 mcKineCuts1->SetRapidityRange(ymin1,ymax1);
 mcKineCuts1->SetChargeMC(charge1);

 AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
 mcKineCuts2->SetPtRange(ptmin2,ptmax2);
 mcKineCuts2->SetRapidityRange(ymin2,ymax2);
 mcKineCuts2->SetChargeMC(charge2);

 AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts");
 mcGenCuts1->SetRequireIsPrimary();
 mcGenCuts1->SetRequirePdgCode(PDG1);

 AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts");
 mcGenCuts2->SetRequireIsPrimary();
 mcGenCuts2->SetRequirePdgCode(PDG2);

 //############# Acceptance Cuts  ????????
 AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
 mcAccCuts->SetMinNHitITS(mintrackrefsITS1);
 mcAccCuts->SetMinNHitTPC(mintrackrefsTPC1);

 //############# Rec-Level kinematic cuts
 AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
 recKineCuts1->SetPtRange(ptmin1,ptmax1);
 recKineCuts1->SetRapidityRange(ymin1,ymax1);
 recKineCuts1->SetChargeRec(charge1);

 AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
 recKineCuts2->SetPtRange(ptmin2,ptmax2);
 recKineCuts2->SetRapidityRange(ymin2,ymax2);
 recKineCuts2->SetChargeRec(charge2);

 AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts","rec-level quality cuts");
 recQualityCuts->SetMinNClusterTPC(minclustersTPC1);
 recQualityCuts->SetRequireITSRefit(kTRUE);

 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
 recIsPrimaryCuts->SetMaxNSigmaToVertex(maxnsigmatovertex1);

 AliCFTrackCutPid* cutPID1 = new AliCFTrackCutPid("cutPID1","ESD_PID") ;
 AliCFTrackCutPid* cutPID2 = new AliCFTrackCutPid("cutPID2","ESD_PID") ;
 int n_species = AliPID::kSPECIES ;
 Double_t* prior = new Double_t[n_species];
 
 prior[0] = 0.0244519 ;
 prior[1] = 0.0143988 ;
 prior[2] = 0.805747  ;
 prior[3] = 0.0928785 ;
 prior[4] = 0.0625243 ;

 cutPID1->SetPriors(prior);
 cutPID1->SetProbabilityCut(0.0);
 cutPID1->SetDetectors("TPC TOF");
 switch(TMath::Abs(PDG1)) {
 case 11   : cutPID1->SetParticleType(AliPID::kElectron, kTRUE); break;
 case 13   : cutPID1->SetParticleType(AliPID::kMuon    , kTRUE); break;
 case 211  : cutPID1->SetParticleType(AliPID::kPion    , kTRUE); break;
 case 321  : cutPID1->SetParticleType(AliPID::kKaon    , kTRUE); break;
 case 2212 : cutPID1->SetParticleType(AliPID::kProton  , kTRUE); break;
 default   : printf("UNDEFINED PID\n"); break;
 }

 cutPID2->SetPriors(prior);
 cutPID2->SetProbabilityCut(0.0);
 cutPID2->SetDetectors("TPC TOF");
 switch(TMath::Abs(PDG2)) {
 case 11   : cutPID2->SetParticleType(AliPID::kElectron, kTRUE); break;
 case 13   : cutPID2->SetParticleType(AliPID::kMuon    , kTRUE); break;
 case 211  : cutPID2->SetParticleType(AliPID::kPion    , kTRUE); break;
 case 321  : cutPID2->SetParticleType(AliPID::kKaon    , kTRUE); break;
 case 2212 : cutPID2->SetParticleType(AliPID::kProton  , kTRUE); break;
 default   : printf("UNDEFINED PID\n"); break;
 }

 printf("CREATE MC KINE CUTS\n");
 TObjArray* mcList1 = new TObjArray(0);
 mcList1->AddLast(mcKineCuts1);
 mcList1->AddLast(mcGenCuts1);

 TObjArray* mcList2 = new TObjArray(0);
 mcList2->AddLast(mcKineCuts2);
 mcList2->AddLast(mcGenCuts2);

 printf("CREATE ACCEPTANCE CUTS\n");
 TObjArray* accList = new TObjArray(0) ;
 accList->AddLast(mcAccCuts);

 printf("CREATE RECONSTRUCTION CUTS\n");
 TObjArray* recList1 = new TObjArray(0) ;
 recList1->AddLast(recKineCuts1);
 recList1->AddLast(recQualityCuts);
 recList1->AddLast(recIsPrimaryCuts);

 TObjArray* recList2 = new TObjArray(0) ;
 recList2->AddLast(recKineCuts2);
 recList2->AddLast(recQualityCuts);
 recList2->AddLast(recIsPrimaryCuts);

 printf("CREATE PID CUTS\n");
 TObjArray* fPIDCutList1 = new TObjArray(0) ;
 fPIDCutList1->AddLast(cutPID1);

 TObjArray* fPIDCutList2 = new TObjArray(0) ;
 fPIDCutList2->AddLast(cutPID2);

printf("CREATE INTERFACE AND CUTS\n");
 AliCFManager* cfmgr1 = new AliCFManager();
 cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);
 //cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);

 AliCFManager* cfmgr2 = new AliCFManager();
 cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
 //cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);
 
 
 
 
 //____________________________________________//
 // Make the analysis manager
 AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

 if (type == "ESD"){
   AliVEventHandler* esdH = new AliESDInputHandler;
   //esdH->SetInactiveBranches("FMD CaloCluster");   //needed?
   mgr->SetInputEventHandler(esdH); }

 if (type == "AOD"){
   AliVEventHandler* aodH = new AliAODInputHandler;
   mgr->SetInputEventHandler(aodH); }

 if (type == "MC" || type == "ESDMC0" || type == "ESDMC1"){
   AliVEventHandler* esdH = new AliESDInputHandler;
   mgr->SetInputEventHandler(esdH);

   AliMCEventHandler *mc = new AliMCEventHandler();
   mgr->SetMCtruthEventHandler(mc); }

  
  
  
  
  
  
  
  
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskCumulants *task1 = new AliAnalysisTaskCumulants("TaskCumulants");
  task1->SetAnalysisType(type);
  task1->SetCFManager1(cfmgr1);
  task1->SetCFManager2(cfmgr2);
  mgr->AddTask(task1);

 
  // Create containers for input/output
 AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
 TString outputName = "outputFromCumulantAnalysis";
 outputName+=type;
 outputName+=".root";
 AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,outputName);
  
 //added according to Andrei:
 //coutput1->SetSpecialOutput();
  
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // old way with a chain
  //  mgr->StartAnalysis("proof",chain);
  mgr->StartAnalysis("proof",data,nRuns,offset);

  timer.Stop();
  timer.Print();
}
