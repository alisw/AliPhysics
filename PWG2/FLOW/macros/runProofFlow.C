//RUN SETTINGS

//Flow analysis method can be:
// SP    = Scalar Product
// LYZ1  = Lee Yang Zeroes first run
// LYZ2  = Lee Yang Zeroes second run
// LYZEP = Lee Yang Zeroes Event Plane
// CUM   = Cumulants
// MCEP  = Flow calculated from the real MC event plane (only for simulated data)
const TString method = "LYZ1";


//analysis type can be ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";

//Bolean to fill/not fill the QA histograms
Bool_t QA = kFALSE;   

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
const Int_t PDG2 = 211;
const Int_t minclustersTPC2 = 50;
const Int_t maxnsigmatovertex2 = 3;

//ESD (pp)
//void runProofFlow(const Char_t* data="/COMMON/COMMON/LHC08c11_10TeV_0.5T", Int_t nRuns=-1, Int_t offset=0) {

//AOD
//void runProofFlow(const Char_t* data="/PWG2/pulvir/PDC08_pythia10TeV_ESD", Int_t nRuns=100, Int_t offset=0) {
//void runProofFlow(const Char_t* data="/PWG2/mvala/pp_09_run82xxT_ESD", Int_t nRuns=100, Int_t offset=0) {
//void runProofFlow(const Char_t* data="/PWG2/pulvir/PDC08_pythia10TeV_MC", Int_t nRuns=100, Int_t offset=0) {
void runProofFlow(const Char_t* data="/PWG2/akisiel/Therminator_c2030", Int_t nRuns=-1, Int_t offset=0) {
//  void runProofFlow(const Char_t* data="/PWG2/nkolk/myDataSet", Int_t nRuns=-1, Int_t offset=0) {


 TStopwatch timer;
 timer.Start();

 printf("*** Connect to PROOF ***\n");
 // TProof::Open("snelling@alicecaf.cern.ch");
 TProof::Open("snelling@localhost");
 

 gProof->UploadPackage("AF-v4-15"); 
 gProof->EnablePackage("AF-v4-15");
 // gProof->UploadPackage("STEERBase.par");
 // gProof->EnablePackage("STEERBase");
 // gProof->UploadPackage("ESD.par");
 // gProof->EnablePackage("ESD");
 // gProof->UploadPackage("AOD.par");
 // gProof->EnablePackage("AOD");
 // gProof->UploadPackage("ANALYSIS.par");
 // gProof->EnablePackage("ANALYSIS");
 // gProof->UploadPackage("ANALYSISalice.par");
 // gProof->EnablePackage("ANALYSISalice");
 // gProof->UploadPackage("PWG2AOD.par");
 // gProof->EnablePackage("PWG2AOD");
 gProof->UploadPackage("CORRFW.par");
 gProof->EnablePackage("CORRFW");
 //gProof->ClearPackage("PWG2flow");
 gProof->UploadPackage("PWG2flow.par");
 gProof->EnablePackage("PWG2flow");


//____________________________________________//
 //Create cuts using correction framework

 //Set TList for the QA histograms
 if (QA) {
   TList* qaInt = new TList();
   TList* qaDiff = new TList();
 }

 //############# cuts on MC
 AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
 mcKineCuts1->SetPtRange(ptmin1,ptmax1);
 mcKineCuts1->SetRapidityRange(ymin1,ymax1);
 mcKineCuts1->SetChargeMC(charge1);
 if (QA) { mcKineCuts1->SetQAOn(qaInt); }
 
 AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
 mcKineCuts2->SetPtRange(ptmin2,ptmax2);
 mcKineCuts2->SetRapidityRange(ymin2,ymax2);
 mcKineCuts2->SetChargeMC(charge2);
 if (QA) { mcKineCuts2->SetQAOn(qaDiff); }
 
 AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts");
 mcGenCuts1->SetRequireIsPrimary();
 mcGenCuts1->SetRequirePdgCode(PDG1);
 if (QA) { mcGenCuts1->SetQAOn(qaInt); }
 
 AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts");
 mcGenCuts2->SetRequireIsPrimary();
 mcGenCuts2->SetRequirePdgCode(PDG2);
 if (QA) { mcGenCuts2->SetQAOn(qaDiff); }
 
 //############# Acceptance Cuts  
 AliCFAcceptanceCuts *mcAccCuts1 = new AliCFAcceptanceCuts("mcAccCuts1","MC acceptance cuts");
 mcAccCuts1->SetMinNHitITS(mintrackrefsITS1);
 mcAccCuts1->SetMinNHitTPC(mintrackrefsTPC1);
 if (QA) { mcAccCuts1->SetQAOn(qaInt); }
 
 AliCFAcceptanceCuts *mcAccCuts2 = new AliCFAcceptanceCuts("mcAccCuts2","MC acceptance cuts");
 mcAccCuts2->SetMinNHitITS(mintrackrefsITS2);
 mcAccCuts2->SetMinNHitTPC(mintrackrefsTPC2);
 if (QA) { mcAccCuts2->SetQAOn(qaDiff); }
 
 //############# Rec-Level kinematic cuts
 AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
 recKineCuts1->SetPtRange(ptmin1,ptmax1);
 recKineCuts1->SetRapidityRange(ymin1,ymax1);
 recKineCuts1->SetChargeRec(charge1);
 if (QA) { recKineCuts1->SetQAOn(qaInt); }
 
 AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
 recKineCuts2->SetPtRange(ptmin2,ptmax2);
 recKineCuts2->SetRapidityRange(ymin2,ymax2);
 recKineCuts2->SetChargeRec(charge2);
 if (QA) { recKineCuts2->SetQAOn(qaDiff); }
 
 AliCFTrackQualityCuts *recQualityCuts1 = new AliCFTrackQualityCuts("recQualityCuts1","rec-level quality cuts");
 recQualityCuts1->SetMinNClusterTPC(minclustersTPC1);
 // recQualityCuts1->SetRequireITSRefit(kTRUE);
 if (QA) { recQualityCuts1->SetQAOn(qaInt); }
 
 AliCFTrackQualityCuts *recQualityCuts2 = new AliCFTrackQualityCuts("recQualityCuts2","rec-level quality cuts");
 recQualityCuts2->SetMinNClusterTPC(minclustersTPC2);
 // recQualityCuts2->SetRequireITSRefit(kTRUE);
 if (QA) { recQualityCuts2->SetQAOn(qaDiff); }
 
 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts1 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts1","rec-level isPrimary cuts");
 recIsPrimaryCuts1->SetMaxNSigmaToVertex(maxnsigmatovertex1);
 if (QA) { recIsPrimaryCuts1->SetQAOn(qaInt); }
 
 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts2 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts2","rec-level isPrimary cuts");
 recIsPrimaryCuts2->SetMaxNSigmaToVertex(maxnsigmatovertex2);
 if (QA) { recIsPrimaryCuts2->SetQAOn(qaDiff); }
 
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
 if (QA) { cutPID1->SetQAOn(qaInt);
   cutPID2->SetQAOn(qaDiff); }
 
 printf("CREATE MC KINE CUTS\n");
 TObjArray* mcList1 = new TObjArray(0);
 mcList1->AddLast(mcKineCuts1);
 mcList1->AddLast(mcGenCuts1);
 
 TObjArray* mcList2 = new TObjArray(0);
 mcList2->AddLast(mcKineCuts2);
 mcList2->AddLast(mcGenCuts2);

 printf("CREATE ACCEPTANCE CUTS\n");
 TObjArray* accList1 = new TObjArray(0) ;
 accList1->AddLast(mcAccCuts1);

 TObjArray* accList2 = new TObjArray(0) ;
 accList2->AddLast(mcAccCuts2);

 printf("CREATE RECONSTRUCTION CUTS\n");
 TObjArray* recList1 = new TObjArray(0) ;
 recList1->AddLast(recKineCuts1);
 recList1->AddLast(recQualityCuts1);
 recList1->AddLast(recIsPrimaryCuts1);

 TObjArray* recList2 = new TObjArray(0) ;
 recList2->AddLast(recKineCuts2);
 recList2->AddLast(recQualityCuts2);
 recList2->AddLast(recIsPrimaryCuts2);

 printf("CREATE PID CUTS\n");
 TObjArray* fPIDCutList1 = new TObjArray(0) ;
 fPIDCutList1->AddLast(cutPID1);

 TObjArray* fPIDCutList2 = new TObjArray(0) ;
 fPIDCutList2->AddLast(cutPID2);

 printf("CREATE INTERFACE AND CUTS\n");
 AliCFManager* cfmgr1 = new AliCFManager();
 cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);

 AliCFManager* cfmgr2 = new AliCFManager();
 cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);

 if (method == "LYZ2"){  
   // read the input file from the first run 
   TString inputFileName = "outputLYZ1analysis" ;
   inputFileName += type;
   inputFileName += "_firstrun.root";
   cout<<"The input file is "<<inputFileName.Data()<<endl;
   TFile* fInputFile = new TFile(inputFileName.Data(),"READ");
   if(!fInputFile || fInputFile->IsZombie()) { 
     cerr << " ERROR: NO First Run file... " << endl ; }
   else {
     TList* fInputList = (TList*)fInputFile->Get("cobj1");
     if (!fInputList) {cout<<"list is NULL pointer!"<<endl;}
   }
   cout<<"input file/list read..."<<endl;
 }

 if (method == "LYZEP") {
   // read the input file from the second LYZ run
   TString inputFileName = "outputLYZ2analysis" ;
   inputFileName += type;
   inputFileName += "_secondrun.root";
   cout<<"The input file is "<<inputFileName.Data()<<endl;
   TFile* fInputFile = new TFile(inputFileName.Data(),"READ");
   if(!fInputFile || fInputFile->IsZombie()) { 
     cerr << " ERROR: NO First Run file... " << endl ; }
   else {
     TList* fInputList = (TList*)fInputFile->Get("cobj1");
     if (!fInputList) {cout<<"list is NULL pointer!"<<endl;}
   }
   cout<<"input file/list read..."<<endl;
 }
    
 //____________________________________________//
 // Make the analysis manager
 AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

 if (type == "ESD"){
   AliVEventHandler* esdH = new AliESDInputHandler;
   mgr->SetInputEventHandler(esdH);
   if (method == "MCEP") {
     AliMCEventHandler *mc = new AliMCEventHandler();
     mgr->SetMCtruthEventHandler(mc);}  }

 if (type == "AOD"){
   AliVEventHandler* aodH = new AliAODInputHandler;
   mgr->SetInputEventHandler(aodH); 
   if (method == "MCEP") {
     AliMCEventHandler *mc = new AliMCEventHandler();
     mgr->SetMCtruthEventHandler(mc);}  }

 if (type == "MC" || type == "ESDMC0" || type == "ESDMC1"){
   AliVEventHandler* esdH = new AliESDInputHandler;
   mgr->SetInputEventHandler(esdH);

   AliMCEventHandler *mc = new AliMCEventHandler();
   mgr->SetMCtruthEventHandler(mc); }

 //____________________________________________//
 // tasks
 if (method == "SP"){
   if (QA) { AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct",kTRUE); }
   else { AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct",kFALSE); }
   taskSP->SetAnalysisType(type);
   taskSP->SetCFManager1(cfmgr1);
   taskSP->SetCFManager2(cfmgr2);
   if (QA) { 
     taskSP->SetQAList1(qaInt);
     taskSP->SetQAList2(qaDiff); }
   mgr->AddTask(taskSP);
  }
  else if (method == "LYZ1"){
    if (QA) { AliAnalysisTaskLeeYangZeros *taskLYZ1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE,kTRUE);}
    else { AliAnalysisTaskLeeYangZeros *taskLYZ1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE,kFALSE);}
    taskLYZ1->SetAnalysisType(type);
    taskLYZ1->SetFirstRunLYZ(kTRUE);
    taskLYZ1->SetUseSumLYZ(kTRUE);
    taskLYZ1->SetCFManager1(cfmgr1);
    taskLYZ1->SetCFManager2(cfmgr2);
    if (QA) { 
      taskLYZ1->SetQAList1(qaInt);
      taskLYZ1->SetQAList2(qaDiff);}
    mgr->AddTask(taskLYZ1);
  }
  else if (method == "LYZ2"){
    if (QA) { AliAnalysisTaskLeeYangZeros *taskLYZ2 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE,kTRUE);}
    else { AliAnalysisTaskLeeYangZeros *taskLYZ2 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE,kFALSE); }
    taskLYZ2->SetAnalysisType(type);
    taskLYZ2->SetFirstRunLYZ(kFALSE);
    taskLYZ2->SetUseSumLYZ(kTRUE);
    taskLYZ2->SetCFManager1(cfmgr1);
    taskLYZ2->SetCFManager2(cfmgr2);
    if (QA) { 
      taskLYZ2->SetQAList1(qaInt);
      taskLYZ2->SetQAList2(qaDiff); }
    mgr->AddTask(taskLYZ2);
  }
  else if (method == "LYZEP"){
    if (QA) { AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane",kTRUE); }
    else { AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane",kFALSE); }
    taskLYZEP->SetAnalysisType(type);
    taskLYZEP->SetCFManager1(cfmgr1);
    taskLYZEP->SetCFManager2(cfmgr2);
    if (QA) { 
      taskLYZEP->SetQAList1(qaInt);
      taskLYZEP->SetQAList2(qaDiff); }
    mgr->AddTask(taskLYZEP);
  }
  else if (method == "CUM"){
    if (QA) { AliAnalysisTaskCumulants *taskCUM = new AliAnalysisTaskCumulants("TaskCumulants",kTRUE);}
    else { AliAnalysisTaskCumulants *taskCUM = new AliAnalysisTaskCumulants("TaskCumulants",kFALSE);}
    taskCUM->SetAnalysisType(type);
    taskCUM->SetCFManager1(cfmgr1);
    taskCUM->SetCFManager2(cfmgr2);
    if (QA) { 
      taskCUM->SetQAList1(qaInt);
      taskCUM->SetQAList2(qaDiff); }
    mgr->AddTask(taskCUM);
  }
  else if (method == "MCEP"){
    if (QA) { AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane",kTRUE);}
    else { AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane",kFALSE);}
    taskMCEP->SetAnalysisType(type);
    taskMCEP->SetCFManager1(cfmgr1);
    taskMCEP->SetCFManager2(cfmgr2);
    if (QA) { 
      taskMCEP->SetQAList1(qaInt);
      taskMCEP->SetQAList2(qaDiff); }
    mgr->AddTask(taskMCEP);
  }


 // Create containers for input/output
 AliAnalysisDataContainer *cinput1 = 
    mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);

  if (method == "LYZ2" || method == "LYZEP"){ 
    AliAnalysisDataContainer *cinput2 = 
		    mgr->CreateContainer("cobj2",TList::Class(),AliAnalysisManager::kInputContainer); } 


  TString outputName = "output";
  outputName+= method;
  outputName+= "analysis";
  outputName+= type;
  if (method == "LYZ1") outputName+= "_firstrun";
  if (method == "LYZ2") outputName+= "_secondrun";
  outputName+= ".root";
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("cobj1", TList::Class(),AliAnalysisManager::kOutputContainer,outputName);

  if (QA) { 
    TString qaNameInt = "QAforInt_";
    qaNameInt += method;
    qaNameInt += "_";
    qaNameInt += type;
    qaNameInt += ".root";
    AliAnalysisDataContainer *coutput2 = 
      mgr->CreateContainer("QAint", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameInt);

    TString qaNameDiff = "QAforDiff_";
    qaNameDiff += method;
    qaNameDiff += "_";
    qaNameDiff += type;
    qaNameDiff += ".root";
    AliAnalysisDataContainer *coutput3 = 
      mgr->CreateContainer("QAdiff", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiff);
  }
 
 //____________________________________________//
  
  if (method == "SP")   { 
    mgr->ConnectInput(taskSP,0,cinput1); 
    mgr->ConnectOutput(taskSP,0,coutput1);
    if (QA) { mgr->ConnectOutput(taskSP,1,coutput2);
    mgr->ConnectOutput(taskSP,2,coutput3); }
  } 
  else if (method == "LYZ1") { 
    mgr->ConnectInput(taskLYZ1,0,cinput1); 
    mgr->ConnectOutput(taskLYZ1,0,coutput1);
    if (QA) { mgr->ConnectOutput(taskLYZ1,1,coutput2);
    mgr->ConnectOutput(taskLYZ1,2,coutput3); }
  }  
  else if (method == "LYZ2") { 
    mgr->ConnectInput(taskLYZ2,0,cinput1); 
    mgr->ConnectInput(taskLYZ2,1,cinput2);
    mgr->ConnectOutput(taskLYZ2,0,coutput1);
    if (QA) { mgr->ConnectOutput(taskLYZ2,1,coutput2);
    mgr->ConnectOutput(taskLYZ2,2,coutput3); }
    cinput2->SetData(fInputList);
  }  
  else if (method == "LYZEP") { 
    mgr->ConnectInput(taskLYZEP,0,cinput1); 
    mgr->ConnectInput(taskLYZEP,1,cinput2);
    mgr->ConnectOutput(taskLYZEP,0,coutput1);
    if (QA) { mgr->ConnectOutput(taskLYZEP,1,coutput2);
    mgr->ConnectOutput(taskLYZEP,2,coutput3); }
    cinput2->SetData(fInputList);
  }
  else if (method == "CUM")   { 
    mgr->ConnectInput(taskCUM,0,cinput1); 
    mgr->ConnectOutput(taskCUM,0,coutput1);
    if (QA) { mgr->ConnectOutput(taskCUM,1,coutput2);
    mgr->ConnectOutput(taskCUM,2,coutput3); }
  } 
  else if (method == "MCEP")  { 
    mgr->ConnectInput(taskMCEP,0,cinput1); 
    mgr->ConnectOutput(taskMCEP,0,coutput1);
    if (QA) { mgr->ConnectOutput(taskMCEP,1,coutput2);
    mgr->ConnectOutput(taskMCEP,2,coutput3); }
  }   

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // old way with a chain
  //  mgr->StartAnalysis("proof",chain);
  mgr->StartAnalysis("proof",data,nRuns,offset);

  timer.Stop();
  timer.Print();
}
