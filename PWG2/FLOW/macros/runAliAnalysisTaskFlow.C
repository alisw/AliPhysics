/////////////////////////////////////////////////////////////////////////////////
//
// HOW TO USE THIS MACRO:
//
// With this macro several flow analysis can be run.
// SP    = Scalar Product                (for PbPb or pp)
// LYZ1  = Lee Yang Zeroes first run     (for PbPb)
// LYZ2  = Lee Yang Zeroes second run    (for PbPb)
// LYZEP = Lee Yang Zeroes Event Plane   (for PbPb)
// GFC   = Cumulants                     (for PbPb)
// QC    = Q-cumulants                   (for PbPb or pp)
// FQD   = Fitting q-distribution        (for PbPb)
// MCEP  = Flow calculated from the real MC event plane (for PbPb only)
//
// The LYZ analysis should be done in the following order;
// LYZ1 -> LYZ2 -> LYZEP,
// because LYZ2 depends on the outputfile of LYZ1 and LYZEP on the outputfile
// of LYZ2.
//
// The MCEP method is a reference method. 
// It can only be run when MC information (kinematics.root & galice.root file) 
// is available in which the reaction plane is stored.
//
// One can run on ESD, AOD or MC.
// Additional options are ESDMC0, ESDMC1. In these options the ESD and MC 
// information is combined. Tracks are selected in the ESD, the PID information 
// is taken from the MC (perfect PID). For ESDMC0 the track kinematics is taken 
// from the ESD and for ESDMC1 it is taken from the MC information.
//
// the macro can be used to run local in aliroot or root, on the grid and on caf
///////////////////////////////////////////////////////////////////////////////////

enum anaModes {mLocal,mLocalPAR,mPROOF,mGRID};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF

// RUN SETTINGS

// Flow analysis method can be:(set to kTRUE or kFALSE)
Bool_t SP    = kTRUE;
Bool_t LYZ1  = kTRUE;
Bool_t LYZ2  = kFALSE;
Bool_t LYZEP = kFALSE;
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kTRUE;
Bool_t MCEP  = kTRUE;

// Analysis type can be ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";

// Boolean to fill/not fill the QA histograms
Bool_t QA = kFALSE;   

// Weights 
//Use weights for Q vector
Bool_t usePhiWeights = kFALSE; //Phi
Bool_t usePtWeights  = kFALSE; //v'(pt)
Bool_t useEtaWeights = kFALSE; //v'(eta)
Bool_t useWeights = usePhiWeights||usePtWeights||useEtaWeights;


// SETTING THE CUTS

// For integrated flow
const Double_t ptmin1 = 0.0;
const Double_t ptmax1 = 10.0;
const Double_t ymin1  = -1.;
const Double_t ymax1  = 1.;
const Int_t mintrackrefsTPC1 = 2;
const Int_t mintrackrefsITS1 = 3;
const Int_t charge1 = 1;
Bool_t UsePIDIntegratedFlow = kFALSE;
const Int_t PDG1 = 211;
const Int_t minclustersTPC1 = 50;
const Int_t maxnsigmatovertex1 = 3;

// For differential flow
const Double_t ptmin2 = 0.0;
const Double_t ptmax2 = 10.0;
const Double_t ymin2  = -1.;
const Double_t ymax2  = 1.;
const Int_t mintrackrefsTPC2 = 2;
const Int_t mintrackrefsITS2 = 3;
const Int_t charge2 = 1;
Bool_t UsePIDDifferentialFlow = kFALSE;
const Int_t PDG2 = 211;
const Int_t minclustersTPC2 = 50;
const Int_t maxnsigmatovertex2 = 3;

// ESD data on PROOF cluster
//void runAliAnalysisTaskFlow(Int_t mode=mPROOF, const Char_t* data="/PWG2/akisiel/Therminator_c2030", Int_t nRuns=-1, Int_t offset=0)

// AOD data on PROOF cluster
//void runAliAnalysisTaskFlow(Int_t mode=mPROOF, const Char_t* data="/PWG2/nkolk/myDataSet", Int_t nRuns=-1, Int_t offset=0)
//void runAliAnalysisTaskFlow(Int_t mode=mPROOF, const Char_t* data="/PWG2/akisiel/Therminator_midcentral_AOD", Int_t nRuns=44, Int_t offset=0)

// Data at Nikhef
void runAliAnalysisTaskFlow(Int_t mode=mLocal, Int_t nRuns = 64, const Char_t* dataDir="/data/alice2/kolk/Therminator_midcentral", Int_t offset = 0) 
//void runAliAnalysisTaskFlow(Int_t mode=mLocalPAR, Int_t nRuns = 55, const Char_t* dataDir="/data/alice2/kolk/Therminator_midcentral", Int_t offset = 0) 

// Data on my mac
//void runAliAnalysisTaskFlow(Int_t mode=mLocal, Int_t nRuns = -1, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0) 
//void runAliAnalysisTaskFlow(Int_t mode=mLocalPAR, Int_t nRuns = 55, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0) 

{

 TStopwatch timer;
 timer.Start();

 if (LYZ1 && LYZ2) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(); }
 if (LYZ2 && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
 if (LYZ1 && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }

 LoadLibraries(mode);

if (mode==mLocal || mode == mLocalPAR || mode == mGRID) {
  if (type!="AOD") { TChain* chain = CreateESDChain(dataDir, nRuns, offset);}
  else { TChain* chain = CreateAODChain(dataDir, nRuns, offset);}
}

 //Create cuts using correction framework

 //Set TList for the QA histograms
 if (QA) {
   TList* qaIntFE = new TList(); TList* qaDiffFE = new TList();
 }

 //############# cuts on MC
 AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
 mcKineCuts1->SetPtRange(ptmin1,ptmax1);
 mcKineCuts1->SetRapidityRange(ymin1,ymax1);
 mcKineCuts1->SetChargeMC(charge1);
 if (QA) { 
   mcKineCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
 mcKineCuts2->SetPtRange(ptmin2,ptmax2);
 mcKineCuts2->SetRapidityRange(ymin2,ymax2);
 mcKineCuts2->SetChargeMC(charge2);
 if (QA) { 
   mcKineCuts2->SetQAOn(qaDiffFE);
 }

 AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts for integrated flow");
 mcGenCuts1->SetRequireIsPrimary();
 if (UsePIDIntegratedFlow) {mcGenCuts1->SetRequirePdgCode(PDG1);}
 if (QA) { 
   mcGenCuts1->SetQAOn(qaIntFE);
 }

 AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts for differential flow");
 mcGenCuts2->SetRequireIsPrimary();
 if (UsePIDDifferentialFlow) {mcGenCuts2->SetRequirePdgCode(PDG2);}
 if (QA) { 
   mcGenCuts2->SetQAOn(qaDiffFE);
 }

 //############# Acceptance Cuts  
 AliCFAcceptanceCuts *mcAccCuts1 = new AliCFAcceptanceCuts("mcAccCuts1","MC acceptance cuts");
 mcAccCuts1->SetMinNHitITS(mintrackrefsITS1);
 mcAccCuts1->SetMinNHitTPC(mintrackrefsTPC1);
 if (QA) { 
   mcAccCuts1->SetQAOn(qaIntFE);
 }

 AliCFAcceptanceCuts *mcAccCuts2 = new AliCFAcceptanceCuts("mcAccCuts2","MC acceptance cuts");
 mcAccCuts2->SetMinNHitITS(mintrackrefsITS2);
 mcAccCuts2->SetMinNHitTPC(mintrackrefsTPC2);
 if (QA) { 
   mcAccCuts2->SetQAOn(qaDiffFE);
 }

 //############# Rec-Level kinematic cuts
 AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
 recKineCuts1->SetPtRange(ptmin1,ptmax1);
 recKineCuts1->SetRapidityRange(ymin1,ymax1);
 recKineCuts1->SetChargeRec(charge1);
 if (QA) { 
   recKineCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
 recKineCuts2->SetPtRange(ptmin2,ptmax2);
 recKineCuts2->SetRapidityRange(ymin2,ymax2);
 recKineCuts2->SetChargeRec(charge2);
 if (QA) { 
   recKineCuts2->SetQAOn(qaDiffFE);
 }

 AliCFTrackQualityCuts *recQualityCuts1 = new AliCFTrackQualityCuts("recQualityCuts1","rec-level quality cuts");
 recQualityCuts1->SetMinNClusterTPC(minclustersTPC1);
 recQualityCuts1->SetStatus(AliESDtrack::kITSrefit);
 if (QA) { 
   recQualityCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackQualityCuts *recQualityCuts2 = new AliCFTrackQualityCuts("recQualityCuts2","rec-level quality cuts");
 recQualityCuts2->SetMinNClusterTPC(minclustersTPC2);
 recQualityCuts2->SetStatus(AliESDtrack::kITSrefit);
 if (QA) { 
   recQualityCuts2->SetQAOn(qaDiffFE);
 }

 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts1 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts1","rec-level isPrimary cuts");
 recIsPrimaryCuts1->SetMaxNSigmaToVertex(maxnsigmatovertex1);
 if (QA) { 
   recIsPrimaryCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts2 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts2","rec-level isPrimary cuts");
 recIsPrimaryCuts2->SetMaxNSigmaToVertex(maxnsigmatovertex2);
 if (QA) { 
   recIsPrimaryCuts2->SetQAOn(qaDiffFE);
 }

 int n_species = AliPID::kSPECIES ;
 Double_t* prior = new Double_t[n_species];

 prior[0] = 0.0244519 ;
 prior[1] = 0.0143988 ;
 prior[2] = 0.805747  ;
 prior[3] = 0.0928785 ;
 prior[4] = 0.0625243 ;

 AliCFTrackCutPid* cutPID1 = NULL;
 if(UsePIDIntegratedFlow) {
   AliCFTrackCutPid* cutPID1 = new AliCFTrackCutPid("cutPID1","ESD_PID for integrated flow") ;
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
   if (QA) { 
     cutPID1->SetQAOn(qaIntFE); 
   }
 }
		  
 AliCFTrackCutPid* cutPID2 = NULL;
 if (UsePIDDifferentialFlow) {
   AliCFTrackCutPid* cutPID2 = new AliCFTrackCutPid("cutPID2","ESD_PID for differential flow") ;
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
   if (QA) { 
     cutPID2->SetQAOn(qaIntFE);
   }
 }

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
 if(UsePIDIntegratedFlow) {fPIDCutList1->AddLast(cutPID1);}

 TObjArray* fPIDCutList2 = new TObjArray(0) ;
 if (UsePIDDifferentialFlow)  {fPIDCutList2->AddLast(cutPID2);}

 printf("CREATE INTERFACE AND CUTS\n");
 AliCFManager* cfmgr1 = new AliCFManager();
 cfmgr1->SetNStepParticle(4); 
 cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);

 AliCFManager* cfmgr2 = new AliCFManager();
 cfmgr2->SetNStepParticle(4); 
 cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);

 //weights: 
 TFile *weightsFile = NULL;
 TList *weightsList = NULL;

 if(useWeights) {
   //open the file with the weights:
   weightsFile = TFile::Open("weights.root","READ");
   if(weightsFile) {
     //access the list which holds the histos with weigths:
     weightsList = (TList*)weightsFile->Get("weights");
   }
   else {
     cout<<" WARNING: the file <weights.root> with weights from the previous run was not available."<<endl;
     break;
   } 
 }


 if (LYZ2){  
   // read the input file from the first run 
   TString inputFileNameLYZ2 = "outputLYZ1analysis" ;
   inputFileNameLYZ2 += type;
   inputFileNameLYZ2 += ".root";
   cout<<"The input file is "<<inputFileNameLYZ2.Data()<<endl;
   TFile* fInputFileLYZ2 = new TFile(inputFileNameLYZ2.Data(),"READ");
   if(!fInputFileLYZ2 || fInputFileLYZ2->IsZombie()) { 
     cerr << " ERROR: NO First Run file... " << endl ; 
     break;
   }
   else {
     TList* fInputListLYZ2 = (TList*)fInputFileLYZ2->Get("cobjLYZ1");
     if (!fInputListLYZ2) {cout<<"list is NULL pointer!"<<endl;}
   }
   cout<<"LYZ2 input file/list read..."<<endl;
 }

 if (LYZEP) {
   // read the input file from the second LYZ run
   TString inputFileNameLYZEP = "outputLYZ2analysis" ;
   inputFileNameLYZEP += type;
   inputFileNameLYZEP += ".root";
   cout<<"The input file is "<<inputFileNameLYZEP.Data()<<endl;
   TFile* fInputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
   if(!fInputFileLYZEP || fInputFileLYZEP->IsZombie()) { 
     cerr << " ERROR: NO First Run file... " << endl ; 
     break;
   }
   else {
     TList* fInputListLYZEP = (TList*)fInputFileLYZEP->Get("cobjLYZ2");
     if (!fInputListLYZEP) {cout<<"list is NULL pointer!"<<endl;}
   }
   cout<<"LYZEP input file/list read..."<<endl;
 }

 //____________________________________________//
 // Make the analysis manager
 AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");

 if (type == "ESD"){
   AliVEventHandler* esdH = new AliESDInputHandler;
   mgr->SetInputEventHandler(esdH);
   if (MCEP) {
     AliMCEventHandler *mc = new AliMCEventHandler();
     mgr->SetMCtruthEventHandler(mc);}  }

 if (type == "AOD"){
   AliVEventHandler* aodH = new AliAODInputHandler;
   mgr->SetInputEventHandler(aodH); 
   if (MCEP) {
     AliMCEventHandler *mc = new AliMCEventHandler();
     mgr->SetMCtruthEventHandler(mc);}  }

 if (type == "MC" || type == "ESDMC0" || type == "ESDMC1"){
   AliVEventHandler* esdH = new AliESDInputHandler;
   mgr->SetInputEventHandler(esdH);

   AliMCEventHandler *mc = new AliMCEventHandler();
   mgr->SetMCtruthEventHandler(mc); }
 
 //____________________________________________//
 // tasks
 AliAnalysisTaskFlowEvent *taskFE = NULL;
 if (QA) { 
   taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",kTRUE); 
   taskFE->SetQAList1(qaIntFE);
   taskFE->SetQAList2(qaDiffFE);
   taskFE->SetAnalysisType(type);
   taskFE->SetCFManager1(cfmgr1);
   taskFE->SetCFManager2(cfmgr2);
   mgr->AddTask(taskFE);
 }
 else { 
   taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",kFALSE); 
   taskFE->SetAnalysisType(type);
   taskFE->SetCFManager1(cfmgr1);
   taskFE->SetCFManager2(cfmgr2);
   mgr->AddTask(taskFE);
 }
 if (FQD){
   AliAnalysisTaskFittingQDistribution *taskFQD = new AliAnalysisTaskFittingQDistribution("TaskFittingQDistribution",useWeights);
   taskFQD->SetUsePhiWeights(usePhiWeights); 
   mgr->AddTask(taskFQD);
 }
 if (SP){
   AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct");
   mgr->AddTask(taskSP);
 }
 if (LYZ1){
   AliAnalysisTaskLeeYangZeros *taskLYZ1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE,kFALSE);
   taskLYZ1->SetAnalysisType(type);
   taskLYZ1->SetFirstRunLYZ(kTRUE);
   taskLYZ1->SetUseSumLYZ(kTRUE);
   taskLYZ1->SetCFManager1(cfmgr1);
   taskLYZ1->SetCFManager2(cfmgr2);
   mgr->AddTask(taskLYZ1);
 }
 if (LYZ2){
   AliAnalysisTaskLeeYangZeros *taskLYZ2 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE,kFALSE);
   taskLYZ2->SetAnalysisType(type);
   taskLYZ2->SetFirstRunLYZ(kFALSE);
   taskLYZ2->SetUseSumLYZ(kTRUE);
   taskLYZ2->SetCFManager1(cfmgr1);
   taskLYZ2->SetCFManager2(cfmgr2);
   mgr->AddTask(taskLYZ2);
 }
 if (LYZEP){
   AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane",kFALSE);
   taskLYZEP->SetAnalysisType(type);
   taskLYZEP->SetCFManager1(cfmgr1);
   taskLYZEP->SetCFManager2(cfmgr2);
   mgr->AddTask(taskLYZEP);
 }
 if (GFC){
   AliAnalysisTaskCumulants *taskGFC = new AliAnalysisTaskCumulants("TaskCumulants",useWeights);
   taskGFC->SetUsePhiWeights(usePhiWeights); 
   taskGFC->SetUsePtWeights(usePtWeights);
   taskGFC->SetUseEtaWeights(useEtaWeights); 
   mgr->AddTask(taskGFC);
 }
 if (QC){
   AliAnalysisTaskQCumulants *taskQC = new AliAnalysisTaskQCumulants("TaskQCumulants",useWeights);
   taskQC->SetUsePhiWeights(usePhiWeights); 
   taskQC->SetUsePtWeights(usePtWeights);
   taskQC->SetUseEtaWeights(useEtaWeights); 
   mgr->AddTask(taskQC);
 }
 if (MCEP){
   AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane",kFALSE);
   taskMCEP->SetAnalysisType(type);
   taskMCEP->SetCFManager1(cfmgr1);
   taskMCEP->SetCFManager2(cfmgr2);
   mgr->AddTask(taskMCEP);
 }


 // Create containers for input/output
 AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
 // TString outputFE = "outputFlowEvent";
 // outputFE+= type;
 // outputFE+= ".root";
 AliAnalysisDataContainer *coutputFE = mgr->CreateContainer("cobjFlowEventSimple",  AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
 
 if (useWeights) {    
   AliAnalysisDataContainer *cinputWeights = mgr->CreateContainer("cobjWeights",TList::Class(),AliAnalysisManager::kInputContainer); 
 }

 if (LYZ2){ 
   AliAnalysisDataContainer *cinputLYZ2 = mgr->CreateContainer("cobjLYZ2in",TList::Class(),AliAnalysisManager::kInputContainer); } 
 if (LYZEP){ 
   AliAnalysisDataContainer *cinputLYZEP = mgr->CreateContainer("cobjLYZEPin",TList::Class(),AliAnalysisManager::kInputContainer); } 

 if(SP) {
   TString outputSP = "outputSPanalysis";
   outputSP+= type;
   outputSP+= ".root";
   AliAnalysisDataContainer *coutputSP = mgr->CreateContainer("cobjSP", TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
 }

 if(LYZ1) {
   TString outputLYZ1 = "outputLYZ1analysis";
   outputLYZ1+= type;
   outputLYZ1+= ".root";
   AliAnalysisDataContainer *coutputLYZ1 = mgr->CreateContainer("cobjLYZ1", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1);
 }

 if(LYZ2) {
   TString outputLYZ2 = "outputLYZ2analysis";
   outputLYZ2+= type;
   outputLYZ2+= ".root";
   AliAnalysisDataContainer *coutputLYZ2 = mgr->CreateContainer("cobjLYZ2", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2);
 }

 if(LYZEP) {
   TString outputLYZEP = "outputLYZEPanalysis";
   outputLYZEP+= type;
   outputLYZEP+= ".root";
   AliAnalysisDataContainer *coutputLYZEP = mgr->CreateContainer("cobjLYZEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZEP);
 }

 if(GFC) {
   TString outputGFC = "outputGFCanalysis";
   outputGFC+= type;
   outputGFC+= ".root";
   AliAnalysisDataContainer *coutputGFC = mgr->CreateContainer("cobjGFC", TList::Class(),AliAnalysisManager::kOutputContainer,outputGFC);
 }

 if(QC) {
   TString outputQC = "outputQCanalysis";
   outputQC+= type;
   outputQC+= ".root";
   AliAnalysisDataContainer *coutputQC = mgr->CreateContainer("cobjQC", TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
 }

 if(FQD) {
   TString outputFQD = "outputFQDanalysis";
   outputFQD+= type;
   outputFQD+= ".root";
   AliAnalysisDataContainer *coutputFQD = mgr->CreateContainer("cobjFQD", TList::Class(),AliAnalysisManager::kOutputContainer,outputFQD);
 }

 if(MCEP) {
   TString outputMCEP = "outputMCEPanalysis";
   outputMCEP+= type;
   outputMCEP+= ".root";
   AliAnalysisDataContainer *coutputMCEP = mgr->CreateContainer("cobjMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputMCEP);
 }


 if (QA) { 
   TString qaNameIntFE = "QAforInt_FE_";
   qaNameIntFE += type;
   qaNameIntFE += ".root";
   AliAnalysisDataContainer *coutputQA1FE = 
     mgr->CreateContainer("QAintFE", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntFE);
   
   TString qaNameDiffFE = "QAforDiff_FE_";
   qaNameDiffFE += type;
   qaNameDiffFE += ".root";
   AliAnalysisDataContainer *coutputQA2FE = 
     mgr->CreateContainer("QAdiffFE", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffFE);
 }


 //____________________________________________//


 // the flow event simple is produced here
 mgr->ConnectInput(taskFE,0,cinput1); 
 mgr->ConnectOutput(taskFE,0,coutputFE);
 if (QA) { 
   mgr->ConnectOutput(taskFE,1,coutputQA1FE);
   mgr->ConnectOutput(taskFE,2,coutputQA2FE); 
 }

 if (FQD) { 
   mgr->ConnectInput(taskFQD,0,coutputFE); 
   mgr->ConnectOutput(taskFQD,0,coutputFQD);
   if(useWeights) {
     mgr->ConnectInput(taskFQD,1,cinputWeights);
     cinputWeights->SetData(weightsList);
   } 
 }    
 if (SP) { 
   mgr->ConnectInput(taskSP,0,coutputFE); 
   mgr->ConnectOutput(taskSP,0,coutputSP);
 } 
 if (LYZ1) { 
   mgr->ConnectInput(taskLYZ1,0,cinput1); 
   mgr->ConnectOutput(taskLYZ1,0,coutputLYZ1);
 }  
 if (LYZ2) { 
   mgr->ConnectInput(taskLYZ2,0,cinput1); 
   mgr->ConnectInput(taskLYZ2,1,cinputLYZ2);
   mgr->ConnectOutput(taskLYZ2,0,coutputLYZ2);
   cinputLYZ2->SetData(fInputListLYZ2);
 }  
 if (LYZEP) { 
   mgr->ConnectInput(taskLYZEP,0,cinput1); 
   mgr->ConnectInput(taskLYZEP,1,cinputLYZEP);
   mgr->ConnectOutput(taskLYZEP,0,coutputLYZEP);
   cinputLYZEP->SetData(fInputListLYZEP);
 }
 if (GFC) { 
   mgr->ConnectInput(taskGFC,0,coutputFE); 
   mgr->ConnectOutput(taskGFC,0,coutputGFC);
   if (useWeights) {
     mgr->ConnectInput(taskGFC,1,cinputWeights);
     cinputWeights->SetData(weightsList);
   } 
 }  
 if (QC) { 
   mgr->ConnectInput(taskQC,0,coutputFE); 
   mgr->ConnectOutput(taskQC,0,coutputQC);
   if (useWeights) {
     mgr->ConnectInput(taskQC,1,cinputWeights);
     cinputWeights->SetData(weightsList);
   } 
 }
 if (MCEP) { 
   mgr->ConnectInput(taskMCEP,0,cinput1); 
   mgr->ConnectOutput(taskMCEP,0,coutputMCEP);
 }  



 //----------------------------------------------------------
 // Run the analysis
 //----------------------------------------------------------

 if (!mgr->InitAnalysis()) return;
 mgr->PrintStatus();

 if (mode==mLocal || mode == mLocalPAR) {
   mgr->StartAnalysis("local",chain);
 }
 else if (mode==mPROOF) {
   //  mgr->StartAnalysis("proof",chain);
   mgr->StartAnalysis("proof",data,nRuns,offset);
 }
 else if (mode==mGRID) { 
   mgr->StartAnalysis("local",chain);
 }

 timer.Stop();
 timer.Print();
}


void LoadLibraries(const anaModes mode) {

 //--------------------------------------
 // Load the needed libraries most of them already loaded by aliroot
 //--------------------------------------
 gSystem->Load("libTree.so");
 gSystem->Load("libGeom.so");
 gSystem->Load("libVMC.so");
 gSystem->Load("libXMLIO.so");
 gSystem->Load("libPhysics.so");

 //----------------------------------------------------------
 // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
 //----------------------------------------------------------
 if (mode==mLocal) {
   //--------------------------------------------------------
   // If you want to use already compiled libraries 
   // in the aliroot distribution
   //--------------------------------------------------------
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libCORRFW.so");
   cerr<<"libCORRFW.so loaded..."<<endl;
   gSystem->Load("libPWG2flowCommon.so");
   cerr<<"libPWG2flowCommon.so loaded..."<<endl;
   gSystem->Load("libPWG2flowTasks.so");
   cerr<<"libPWG2flowTasks.so loaded..."<<endl;
 }

 else if (mode == mLocalPAR || mode == mGRID) {
   //--------------------------------------------------------
   //If you want to use root and par files from aliroot
   //--------------------------------------------------------  
   SetupPar("STEERBase");
   SetupPar("ESD");
   SetupPar("AOD");
   SetupPar("ANALYSIS");
   SetupPar("ANALYSISalice");
   SetupPar("PWG2AOD");
   SetupPar("CORRFW");
   SetupPar("PWG2flowCommon");
   cerr<<"PWG2flowCommon.par loaded..."<<endl;
   SetupPar("PWG2flowTasks");
   cerr<<"PWG2flowTasks.par loaded..."<<endl;
 }

 //---------------------------------------------------------
 // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
 //---------------------------------------------------------
 else if (mode==mPROOF) {
   //

   //  set to debug root versus if needed
   //  TProof::Mgr("alicecaf")->SetROOTVersion("v5-21-01-alice_dbg");
   //  TProof::Mgr("alicecaf")->SetROOTVersion("v5-21-01-alice");

   // Connect to proof
   // Put appropriate username here
   // TProof::Reset("proof://snelling@alicecaf.cern.ch"); 
   printf("*** Connect to PROOF ***\n");
   //  TProof::Open("abilandz@alicecaf.cern.ch");
   TProof::Open("snelling@localhost");

   // Enable the STEERBase Package
   gProof->ClearPackage("STEERBase.par");
   gProof->UploadPackage("STEERBase.par");
   gProof->EnablePackage("STEERBase");
   // Enable the ESD Package
   gProof->ClearPackage("ESD.par");
   gProof->UploadPackage("ESD.par");
   gProof->EnablePackage("ESD");
   // Enable the AOD Package
   gProof->ClearPackage("AOD.par");
   gProof->UploadPackage("AOD.par");
   gProof->EnablePackage("AOD");
   // Enable the Analysis Package
   gProof->ClearPackage("ANALYSIS.par");
   gProof->UploadPackage("ANALYSIS.par");
   gProof->EnablePackage("ANALYSIS");
   // Enable the Analysis Package alice
   gProof->ClearPackage("ANALYSISalice.par");
   gProof->UploadPackage("ANALYSISalice.par");
   gProof->EnablePackage("ANALYSISalice");
   // Load the PWG2 AOD
   gProof->ClearPackage("PWG2AOD.par");
   gProof->UploadPackage("PWG2AOD.par");
   gProof->EnablePackage("PWG2AOD");
   // Enable the Correction Framework
   gProof->ClearPackage("CORRFW.par");
   gProof->UploadPackage("CORRFW.par");
   gProof->EnablePackage("CORRFW");
   // Enable Flow Analysis
   gProof->ClearPackage("PWG2flowCommon");
   gProof->UploadPackage("PWG2flowCommon.par");
   gProof->EnablePackage("PWG2flowCommon");
   gProof->ClearPackage("PWG2flowTasks");
   gProof->UploadPackage("PWG2flowTasks.par");
   gProof->EnablePackage("PWG2flowTasks");
   //
   gProof->ShowEnabledPackages();
 }  

}

void SetupPar(char* pararchivename) {
 //Load par files, create analysis libraries
 //For testing, if par file already decompressed and modified
 //classes then do not decompress.

 TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
 TString parpar(Form("%s.par", pararchivename)) ; 
 if ( gSystem->AccessPathName(parpar.Data()) ) {
   gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
   TString processline(Form(".! make %s", parpar.Data())) ; 
   gROOT->ProcessLine(processline.Data()) ;
   gSystem->ChangeDirectory(cdir) ; 
   processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
   gROOT->ProcessLine(processline.Data()) ;
 } 
 if ( gSystem->AccessPathName(pararchivename) ) {  
   TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
   gROOT->ProcessLine(processline.Data());
 }

 TString ocwd = gSystem->WorkingDirectory();
 gSystem->ChangeDirectory(pararchivename);

 // check for BUILD.sh and execute
 if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
   printf("*******************************\n");
   printf("*** Building PAR archive    ***\n");
   cout<<pararchivename<<endl;
   printf("*******************************\n");

   if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
     Error("runProcess","Cannot Build the PAR Archive! - Abort!");
     return -1;
   }
 }
 // check for SETUP.C and execute
 if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
   printf("*******************************\n");
   printf("*** Setup PAR archive       ***\n");
   cout<<pararchivename<<endl;
   printf("*******************************\n");
   gROOT->Macro("PROOF-INF/SETUP.C");
 }

 gSystem->ChangeDirectory(ocwd.Data());
 printf("Current dir: %s\n", ocwd.Data());
}


// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
 // creates chain of files in a given directory or file containing a list.
 // In case of directory the structure is expected as:
 // <aDataDir>/<dir0>/AliESDs.root
 // <aDataDir>/<dir1>/AliESDs.root
 // ...

 if (!aDataDir)
   return 0;

 Long_t id, size, flags, modtime;
 if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
   {
     printf("%s not found.\n", aDataDir);
     return 0;
   }

 TChain* chain = new TChain("esdTree");
 TChain* chaingAlice = 0;

 if (flags & 2)
   {
     TString execDir(gSystem->pwd());
     TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
     TList* dirList            = baseDir->GetListOfFiles();
     Int_t nDirs               = dirList->GetEntries();
     gSystem->cd(execDir);

     Int_t count = 0;

     for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliESDs.root/esdTree");
	  //  cerr<<presentDirName<<endl;
	}

   }
 else
   {
     // Open the input stream
     ifstream in;
     in.open(aDataDir);

     Int_t count = 0;

     // Read the input list of files and add them to the chain
     TString esdfile;
     while(in.good()) {
	in >> esdfile;
	if (!esdfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
       // add esd file
	chain->Add(esdfile);
     }

     in.close();
   }

 return chain;
}

// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
 // creates chain of files in a given directory or file containing a list.
 // In case of directory the structure is expected as:
 // <aDataDir>/<dir0>/AliAOD.root
 // <aDataDir>/<dir1>/AliAOD.root
 // ...

 if (!aDataDir)
   return 0;

 Long_t id, size, flags, modtime;
 if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
   {
     printf("%s not found.\n", aDataDir);
     return 0;
   }

 TChain* chain = new TChain("aodTree");
 TChain* chaingAlice = 0;

 if (flags & 2)
   {
     TString execDir(gSystem->pwd());
     TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
     TList* dirList            = baseDir->GetListOfFiles();
     Int_t nDirs               = dirList->GetEntries();
     gSystem->cd(execDir);

     Int_t count = 0;

     for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliAOD.root/aodTree");
	  // cerr<<presentDirName<<endl;
	}

   }
 else
   {
     // Open the input stream
     ifstream in;
     in.open(aDataDir);

     Int_t count = 0;

     // Read the input list of files and add them to the chain
     TString aodfile;
     while(in.good()) {
	in >> aodfile;
	if (!aodfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
       // add aod file
	chain->Add(aodfile);
     }

     in.close();
   }

 return chain;
}


