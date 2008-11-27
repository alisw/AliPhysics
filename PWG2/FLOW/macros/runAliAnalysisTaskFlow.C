// from CreateESDChain.C - instead of  #include "CreateESDChain.C"
TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 10, Int_t offset = 0) ;
TChain* CreateAODChain(const char* aDataDir = "AODfiles.txt", Int_t aRuns = 10, Int_t offset = 0) ;
void LookupWrite(TChain* chain, const char* target) ;

 ///////////////////////////////////////////////////////////////////////////////// 	 
// 	 
// HOW TO USE THIS MACRO: 	 
// 	 
// With this macro several flow analysis can be run. 	 
// SP    = Scalar Product                (for PbPb or pp) 	 
// LYZ1  = Lee Yang Zeroes first run     (for PbPb) 	 
// LYZ2  = Lee Yang Zeroes second run    (for PbPb) 	 
// LYZEP = Lee Yang Zeroes Event Plane   (for PbPb) 	 
// GFC   = Generating Function Cumulants (for PbPb)  
// QC    = Q-cumulants                   (for PbPb)  
// FQD   = Fitting q-distribution        (for PbPb) 
// MCEP  = Flow calculated from the real MC event plane (for PbPb only!) 	 
// 	 
// The LYZ analysis should be done in the following order; 	 
// LYZ1 -> LYZ2 -> LYZEP, 	 
// because LYZ2 depends on the outputfile of LYZ1 and LYZEP on the outputfile 	 
// of LYZ2. 	 
// 	 
// The MCEP method is a reference method. 	 
// It can only be run when MC information (kinematics.root & galice.root file) is available 	 
// in which the reaction plane is stored. 	 
// 	 
// One can run on ESD, AOD or MC. 	 
// Additional options are ESDMC0, ESDMC1. In these options the ESD and MC 	 
// information is combined. Tracks are selected in the ESD, the PID information 	 
// is taken from the MC (perfect PID). For ESDMC0 the track kinematics is taken 	 
// from the ESD and for ESDMC1 it is taken from the MC information. 	 
// 	 
/////////////////////////////////////////////////////////////////////////////////// 	 
	 

//SETTING THE ANALYSIS

//Flow analysis methods can be: (set to kTRUE or kFALSE)
Bool_t SP    = kFALSE;
Bool_t LYZ1  = kFALSE;
Bool_t LYZ2  = kFALSE;
Bool_t LYZEP = kFALSE;
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kTRUE;
Bool_t MCEP  = kFALSE;


//Type of analysis can be:
// ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";

//Bolean to fill/not fill the QA histograms
Bool_t QA = kTRUE;    

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


void runAliAnalysisTaskFlow(Int_t nRuns = 44, const Char_t* dataDir="/data/alice2/ante/ab2", Int_t offset = 0) 

//void runAliAnalysisTaskFlow(Int_t nRuns = -1, const Char_t* dataDir="/data/alice2/kolk/Therminator_midcentral", Int_t offset = 0) 


//void runAliAnalysisTaskFlow(Int_t nRuns = -1, const Char_t* dataDir="/data/alice2/LHyquid3_rot", Int_t offset = 0) 

//void runAliAnalysisTaskFlow(Int_t nRuns = 4, const Char_t* dataDir="/data/alice2/ante/AOD", Int_t offset = 0) 
  //void runAliAnalysisTaskFlowMore(Int_t nRuns = 3, const Char_t* dataDir="/data/alice1/kolk/TherminatorNov08", Int_t offset = 0) 
  //void runAliAnalysisTaskFlowMore(Int_t nRuns = 1, const Char_t* dataDir="/data/alice1/kolk/AOD", Int_t offset = 0)
  //void runAliAnalysisTaskFlowMore(Int_t nRuns = 3, const Char_t* dataDir="/data/alice1/kolk/testCrashIfNoKine", Int_t offset = 0)
  //void runAliAnalysisTaskFlowMore(Int_t nRuns = 3, const Char_t* dataDir="/data/alice1/kolk/testCrashIfMissingKine", Int_t offset = 0)
{
  TStopwatch timer;
  timer.Start();
  
  if (LYZ1 && LYZ2) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(); }
  
  if (LYZ2 && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
  
  if (LYZ1 && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
  
  
  // include path (to find the .h files when compiling)
  gSystem->AddIncludePath("-I$ALICE_ROOT/include") ;
  gSystem->AddIncludePath("-I$ROOTSYS/include") ;
  
  // load needed libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libESD.so");
  cerr<<"libESD loaded..."<<endl;
  gSystem->Load("libANALYSIS.so");
  cerr<<"libANALYSIS.so loaded..."<<endl;
  gSystem->Load("libANALYSISRL.so");
  cerr<<"libANALYSISRL.so loaded..."<<endl;
  gSystem->Load("libCORRFW.so");
  cerr<<"libCORRFW.so loaded..."<<endl;
  gSystem->Load("libPWG2flow.so");
  cerr<<"libPWG2flow.so loaded..."<<endl;

  // create the TChain. CreateESDChain() is defined in CreateESDChain.C
  if (type!="AOD") { TChain* chain = CreateESDChain(dataDir, nRuns, offset);
  cout<<"chain ("<<chain<<")"<<endl; }
  else { TChain* chain = CreateAODChain(dataDir, nRuns, offset);
  cout<<"chain ("<<chain<<")"<<endl; }
 
  //____________________________________________//
  //Create cuts using correction framework

  //Set TList for the QA histograms
  if (QA) {
    if (SP){
      TList* qaIntSP = new TList();
      TList* qaDiffSP = new TList(); }
    if (LYZ1) {
      TList* qaIntLYZ1 = new TList();
      TList* qaDiffLYZ1 = new TList(); }
    if (LYZ2) {
      TList* qaIntLYZ2 = new TList();
      TList* qaDiffLYZ2 = new TList(); }
    if (LYZEP) {
      TList* qaIntLYZEP = new TList();
      TList* qaDiffLYZEP = new TList(); }
    if (GFC) {
      TList* qaIntGFC = new TList();
      TList* qaDiffGFC = new TList(); }
    if (QC) {
      TList* qaIntQC = new TList();
      TList* qaDiffQC = new TList(); }
    if (FQD) {
      TList* qaIntFQD = new TList();
      TList* qaDiffFQD = new TList(); }
    if (MCEP) {
      TList* qaIntMCEP = new TList();
      TList* qaDiffMCEP = new TList(); }
  }
  
  //############# cuts on MC
  AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
  mcKineCuts1->SetPtRange(ptmin1,ptmax1);
  mcKineCuts1->SetRapidityRange(ymin1,ymax1);
  mcKineCuts1->SetChargeMC(charge1);
  if (QA) { 
    if (SP)   { mcKineCuts1->SetQAOn(qaIntSP); }
    if (LYZ1) { mcKineCuts1->SetQAOn(qaIntLYZ1); }
    if (LYZ2) { mcKineCuts1->SetQAOn(qaIntLYZ2); }
    if (LYZEP){ mcKineCuts1->SetQAOn(qaIntLYZEP); }
    if (GFC)  { mcKineCuts1->SetQAOn(qaIntGFC); }
    if (QC)   { mcKineCuts1->SetQAOn(qaIntQC); }
    if (FQD)  { mcKineCuts1->SetQAOn(qaIntFQD); }
    if (MCEP) { mcKineCuts1->SetQAOn(qaIntMCEP); }
  }
  
  AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
  mcKineCuts2->SetPtRange(ptmin2,ptmax2);
  mcKineCuts2->SetRapidityRange(ymin2,ymax2);
  mcKineCuts2->SetChargeMC(charge2);
  if (QA) { 
    if (SP)   { mcKineCuts2->SetQAOn(qaDiffSP); }
    if (LYZ1) { mcKineCuts2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { mcKineCuts2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ mcKineCuts2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { mcKineCuts2->SetQAOn(qaDiffGFC); }
    if (QC)   { mcKineCuts2->SetQAOn(qaDiffQC); } 
    if (FQD)  { mcKineCuts2->SetQAOn(qaDiffFQD); }    
    if (MCEP) { mcKineCuts2->SetQAOn(qaDiffMCEP); }
  }

  AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts");
  mcGenCuts1->SetRequireIsPrimary();
  mcGenCuts1->SetRequirePdgCode(PDG1);
  if (QA) { 
    if (SP)   { mcGenCuts1->SetQAOn(qaIntSP); }
    if (LYZ1) { mcGenCuts1->SetQAOn(qaIntLYZ1); }
    if (LYZ2) { mcGenCuts1->SetQAOn(qaIntLYZ2); }
    if (LYZEP){ mcGenCuts1->SetQAOn(qaIntLYZEP); }
    if (GFC)  { mcGenCuts1->SetQAOn(qaIntGFC); }
    if (QC)   { mcGenCuts1->SetQAOn(qaIntQC); }
    if (FQD)  { mcGenCuts1->SetQAOn(qaIntFQD); }
    if (MCEP) { mcGenCuts1->SetQAOn(qaIntMCEP); }
  }

  AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts");
  mcGenCuts2->SetRequireIsPrimary();
  mcGenCuts2->SetRequirePdgCode(PDG2);
  if (QA) { 
    if (SP)   { mcGenCuts2->SetQAOn(qaDiffSP); }
    if (LYZ1) { mcGenCuts2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { mcGenCuts2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ mcGenCuts2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { mcGenCuts2->SetQAOn(qaDiffGFC); }
    if (QC)   { mcGenCuts2->SetQAOn(qaDiffQC); }
    if (FQD)  { mcGenCuts2->SetQAOn(qaDiffFQD); }
    if (MCEP) { mcGenCuts2->SetQAOn(qaDiffMCEP); }
  }
  
  //############# Acceptance Cuts
  AliCFAcceptanceCuts *mcAccCuts1 = new AliCFAcceptanceCuts("mcAccCuts1","MC acceptance cuts");
  mcAccCuts1->SetMinNHitITS(mintrackrefsITS1);
  mcAccCuts1->SetMinNHitTPC(mintrackrefsTPC1);
  if (QA) { 
    if (SP)   { mcAccCuts1->SetQAOn(qaIntSP); }
    if (LYZ1) { mcAccCuts1->SetQAOn(qaIntLYZ1); }
    if (LYZ2) { mcAccCuts1->SetQAOn(qaIntLYZ2); }
    if (LYZEP){ mcAccCuts1->SetQAOn(qaIntLYZEP); }
    if (GFC)  { mcAccCuts1->SetQAOn(qaIntGFC); }
    if (QC)   { mcAccCuts1->SetQAOn(qaIntQC); }
    if (FQD)  { mcAccCuts1->SetQAOn(qaIntFQD); }
    if (MCEP) { mcAccCuts1->SetQAOn(qaIntMCEP); }
  }

  AliCFAcceptanceCuts *mcAccCuts2 = new AliCFAcceptanceCuts("mcAccCuts2","MC acceptance cuts");
  mcAccCuts2->SetMinNHitITS(mintrackrefsITS2);
  mcAccCuts2->SetMinNHitTPC(mintrackrefsTPC2);
  if (QA) { 
    if (SP)   { mcAccCuts2->SetQAOn(qaDiffSP); }
    if (LYZ1) { mcAccCuts2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { mcAccCuts2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ mcAccCuts2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { mcAccCuts2->SetQAOn(qaDiffGFC); }
    if (QC)   { mcAccCuts2->SetQAOn(qaDiffQC); }
    if (FQD)  { mcAccCuts2->SetQAOn(qaDiffFQD); }    
    if (MCEP) { mcAccCuts2->SetQAOn(qaDiffMCEP); }
  }
  
  //############# Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
  recKineCuts1->SetPtRange(ptmin1,ptmax1);
  recKineCuts1->SetRapidityRange(ymin1,ymax1);
  recKineCuts1->SetChargeRec(charge1);
  if (QA) { 
    if (SP)   { recKineCuts1->SetQAOn(qaIntSP); }
    if (LYZ1) { recKineCuts1->SetQAOn(qaIntLYZ1); }
    if (LYZ2) { recKineCuts1->SetQAOn(qaIntLYZ2); }
    if (LYZEP){ recKineCuts1->SetQAOn(qaIntLYZEP); }
    if (GFC)  { recKineCuts1->SetQAOn(qaIntGFC); }
    if (QC)   { recKineCuts1->SetQAOn(qaIntQC); }
    if (FQD)  { recKineCuts1->SetQAOn(qaIntFQD); }
    if (MCEP) { recKineCuts1->SetQAOn(qaIntMCEP); }
  }

  AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
  recKineCuts2->SetPtRange(ptmin2,ptmax2);
  recKineCuts2->SetRapidityRange(ymin2,ymax2);
  recKineCuts2->SetChargeRec(charge2);
  if (QA) { 
    if (SP)   { recKineCuts2->SetQAOn(qaDiffSP); }
    if (LYZ1) { recKineCuts2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { recKineCuts2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ recKineCuts2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { recKineCuts2->SetQAOn(qaDiffGFC); }
    if (QC)   { recKineCuts2->SetQAOn(qaDiffQC); }
    if (FQD)  { recKineCuts2->SetQAOn(qaDiffFQD); }
    if (MCEP) { recKineCuts2->SetQAOn(qaDiffMCEP); }
  }
  
  AliCFTrackQualityCuts *recQualityCuts1 = new AliCFTrackQualityCuts("recQualityCuts1","rec-level quality cuts");
  recQualityCuts1->SetMinNClusterTPC(minclustersTPC1);
  recQualityCuts1->SetStatus(AliESDtrack::kITSrefit);
  if (QA) { 
    if (SP)   { recQualityCuts1->SetQAOn(qaIntSP); }
    if (LYZ1) { recQualityCuts1->SetQAOn(qaIntLYZ1); }
    if (LYZ2) { recQualityCuts1->SetQAOn(qaIntLYZ2); }
    if (LYZEP){ recQualityCuts1->SetQAOn(qaIntLYZEP); }
    if (GFC)  { recQualityCuts1->SetQAOn(qaIntGFC); }
    if (QC)   { recQualityCuts1->SetQAOn(qaIntQC); }
    if (FQD)  { recQualityCuts1->SetQAOn(qaIntFQD); }
    if (MCEP) { recQualityCuts1->SetQAOn(qaIntMCEP); }
  }

  AliCFTrackQualityCuts *recQualityCuts2 = new AliCFTrackQualityCuts("recQualityCuts2","rec-level quality cuts");
  recQualityCuts2->SetMinNClusterTPC(minclustersTPC2);
  recQualityCuts2->SetStatus(AliESDtrack::kITSrefit);
  if (QA) { 
    if (SP)   { recQualityCuts2->SetQAOn(qaDiffSP); }
    if (LYZ1) { recQualityCuts2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { recQualityCuts2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ recQualityCuts2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { recQualityCuts2->SetQAOn(qaDiffGFC); }
    if (QC)   { recQualityCuts2->SetQAOn(qaDiffQC); }
    if (FQD)  { recQualityCuts2->SetQAOn(qaDiffFQD); }
    if (MCEP) { recQualityCuts2->SetQAOn(qaDiffMCEP); }
  }

  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts1 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts1","rec-level isPrimary cuts");
  recIsPrimaryCuts1->SetMaxNSigmaToVertex(maxnsigmatovertex1);
  if (QA) { 
    if (SP)   { recIsPrimaryCuts1->SetQAOn(qaIntSP); }
    if (LYZ1) { recIsPrimaryCuts1->SetQAOn(qaIntLYZ1); }
    if (LYZ2) { recIsPrimaryCuts1->SetQAOn(qaIntLYZ2); }
    if (LYZEP){ recIsPrimaryCuts1->SetQAOn(qaIntLYZEP); }
    if (GFC)  { recIsPrimaryCuts1->SetQAOn(qaIntGFC); }
    if (QC)   { recIsPrimaryCuts1->SetQAOn(qaIntQC); }
    if (FQD)  { recIsPrimaryCuts1->SetQAOn(qaIntFQD); }
    if (MCEP) { recIsPrimaryCuts1->SetQAOn(qaIntMCEP); }
  }

  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts2 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts2","rec-level isPrimary cuts");
  recIsPrimaryCuts2->SetMaxNSigmaToVertex(maxnsigmatovertex2);
  if (QA) { 
    if (SP)   { recIsPrimaryCuts2->SetQAOn(qaDiffSP); }
    if (LYZ1) { recIsPrimaryCuts2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { recIsPrimaryCuts2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ recIsPrimaryCuts2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { recIsPrimaryCuts2->SetQAOn(qaDiffGFC); }
    if (QC)   { recIsPrimaryCuts2->SetQAOn(qaDiffQC); }
    if (FQD)  { recIsPrimaryCuts2->SetQAOn(qaDiffFQD); }
    if (MCEP) { recIsPrimaryCuts2->SetQAOn(qaDiffMCEP); }
  }
  
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
  if (QA) { 
    if (SP)   { 
      cutPID1->SetQAOn(qaIntSP); 
      cutPID2->SetQAOn(qaDiffSP); }
    if (LYZ1) { 
      cutPID1->SetQAOn(qaIntLYZ1); 
      cutPID2->SetQAOn(qaDiffLYZ1); }
    if (LYZ2) { 
      cutPID1->SetQAOn(qaIntLYZ2); 
      cutPID2->SetQAOn(qaDiffLYZ2); }
    if (LYZEP){ 
      cutPID1->SetQAOn(qaIntLYZEP); 
      cutPID2->SetQAOn(qaDiffLYZEP); }
    if (GFC)  { 
      cutPID1->SetQAOn(qaIntGFC); 
      cutPID2->SetQAOn(qaDiffGFC); }
    if (QC)  { 
      cutPID1->SetQAOn(qaIntQC); 
      cutPID2->SetQAOn(qaDiffQC); }
    if (FQD)  { 
      cutPID1->SetQAOn(qaIntFQD); 
      cutPID2->SetQAOn(qaDiffFQD); }
    if (MCEP) { 
      cutPID1->SetQAOn(qaIntMCEP); 
      cutPID2->SetQAOn(qaDiffMCEP); }
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
  fPIDCutList1->AddLast(cutPID1);
  
  TObjArray* fPIDCutList2 = new TObjArray(0) ;
  fPIDCutList2->AddLast(cutPID2);
  
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* cfmgr1 = new AliCFManager();
  cfmgr1->SetNStepParticle(4); //05nov08
  cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);       //on MC
  cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList1);      //on MC
  cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);      //on ESD
  cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);  //on ESD
  
  AliCFManager* cfmgr2 = new AliCFManager();
  cfmgr2->SetNStepParticle(4); //05nov08
  cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);
  
  
  if (LYZ2){  
    // read the input file from the first run 
    TString inputFileNameLYZ2 = "outputLYZ1analysis" ;
    inputFileNameLYZ2 += type;
    inputFileNameLYZ2 += "_firstrun.root";
    cout<<"The input file is "<<inputFileNameLYZ2.Data()<<endl;
    TFile* fInputFileLYZ2 = new TFile(inputFileNameLYZ2.Data(),"READ");
    if(!fInputFileLYZ2 || fInputFileLYZ2->IsZombie()) { 
      cerr << " ERROR: NO First Run file... " << endl ; }
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
    inputFileNameLYZEP += "_secondrun.root";
    cout<<"The input file is "<<inputFileNameLYZEP.Data()<<endl;
    TFile* fInputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
    if(!fInputFileLYZEP || fInputFileLYZEP->IsZombie()) { cerr << " ERROR: NO First Run file... " << endl ; }
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
  // Task
  
  if (SP){
    if (QA) { AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct",kTRUE); }
    else { AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct",kFALSE); }
    taskSP->SetAnalysisType(type);
    taskSP->SetCFManager1(cfmgr1);
    taskSP->SetCFManager2(cfmgr2);
    if (QA) { 
      taskSP->SetQAList1(qaIntSP);
      taskSP->SetQAList2(qaDiffSP); }
    mgr->AddTask(taskSP);
  }
  if (LYZ1){
    if (QA) { AliAnalysisTaskLeeYangZeros *taskLYZ1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE,kTRUE);}
    else { AliAnalysisTaskLeeYangZeros *taskLYZ1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE,kFALSE);}
    taskLYZ1->SetAnalysisType(type);
    taskLYZ1->SetFirstRunLYZ(kTRUE);
    taskLYZ1->SetUseSumLYZ(kTRUE);
    taskLYZ1->SetCFManager1(cfmgr1);
    taskLYZ1->SetCFManager2(cfmgr2);
    if (QA) { 
      taskLYZ1->SetQAList1(qaIntLYZ1);
      taskLYZ1->SetQAList2(qaDiffLYZ1);}
    mgr->AddTask(taskLYZ1);
  }
  if (LYZ2){
    if (QA) { AliAnalysisTaskLeeYangZeros *taskLYZ2 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE,kTRUE);}
    else { AliAnalysisTaskLeeYangZeros *taskLYZ2 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE,kFALSE); }
    taskLYZ2->SetAnalysisType(type);
    taskLYZ2->SetFirstRunLYZ(kFALSE);
    taskLYZ2->SetUseSumLYZ(kTRUE);
    taskLYZ2->SetCFManager1(cfmgr1);
    taskLYZ2->SetCFManager2(cfmgr2);
    if (QA) { 
      taskLYZ2->SetQAList1(qaIntLYZ2);
      taskLYZ2->SetQAList2(qaDiffLYZ2); }
    mgr->AddTask(taskLYZ2);
  }
  if (LYZEP){
    if (QA) { AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane",kTRUE); }
    else { AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane",kFALSE); }
    taskLYZEP->SetAnalysisType(type);
    taskLYZEP->SetCFManager1(cfmgr1);
    taskLYZEP->SetCFManager2(cfmgr2);
    if (QA) { 
      taskLYZEP->SetQAList1(qaIntLYZEP);
      taskLYZEP->SetQAList2(qaDiffLYZEP); }
    mgr->AddTask(taskLYZEP);
  }
  if (GFC){
    if (QA) { AliAnalysisTaskCumulants *taskGFC = new AliAnalysisTaskCumulants("TaskCumulants",kTRUE);}
    else { AliAnalysisTaskCumulants *taskGFC = new AliAnalysisTaskCumulants("TaskCumulants",kFALSE);}
    taskGFC->SetAnalysisType(type);
    taskGFC->SetCFManager1(cfmgr1);
    taskGFC->SetCFManager2(cfmgr2);
    if (QA) { 
      taskGFC->SetQAList1(qaIntGFC);
      taskGFC->SetQAList2(qaDiffGFC); }
    mgr->AddTask(taskGFC);
  }
  if (QC){
    if (QA) { AliAnalysisTaskQCumulants *taskQC = new AliAnalysisTaskQCumulants("TaskQCumulants",kTRUE);}
    else { AliAnalysisTaskQCumulants *taskQC = new AliAnalysisTaskQCumulants("TaskQCumulants",kFALSE);}
    taskQC->SetAnalysisType(type);
    taskQC->SetCFManager1(cfmgr1);
    taskQC->SetCFManager2(cfmgr2);
    if (QA) { 
      taskQC->SetQAList1(qaIntQC);
      taskQC->SetQAList2(qaDiffQC); }
    mgr->AddTask(taskQC);
  }
  if (FQD){
    if (QA) { AliAnalysisTaskFittingQDistribution *taskFQD = new AliAnalysisTaskFittingQDistribution("TaskFittingQDistribution",kTRUE);}
    else { AliAnalysisTaskFittingQDistribution *taskFQD = new AliAnalysisTaskFittingQDistribution("TaskFittingQDistribution",kFALSE);}
    taskFQD->SetAnalysisType(type);
    taskFQD->SetCFManager1(cfmgr1);
    taskFQD->SetCFManager2(cfmgr2);
    if (QA) { 
      taskFQD->SetQAList1(qaIntFQD);
      taskFQD->SetQAList2(qaDiffFQD); }
    mgr->AddTask(taskFQD);
  }
  if (MCEP){
    if (QA) { AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane",kTRUE);}
    else { AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane",kFALSE);}
    taskMCEP->SetAnalysisType(type);
    taskMCEP->SetCFManager1(cfmgr1);
    taskMCEP->SetCFManager2(cfmgr2);
    if (QA) { 
      taskMCEP->SetQAList1(qaIntMCEP);
      taskMCEP->SetQAList2(qaDiffMCEP); }
    mgr->AddTask(taskMCEP);
  }
  
  
  // Create containers for input/output
  
  AliAnalysisDataContainer *cinput1 = 
    mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  
  if (LYZ2){ 
    AliAnalysisDataContainer *cinputLYZ2 = 
      mgr->CreateContainer("cobjLYZ2in",TList::Class(),AliAnalysisManager::kInputContainer); } 
  if (LYZEP){ 
    AliAnalysisDataContainer *cinputLYZEP = 
      mgr->CreateContainer("cobjLYZEPin",TList::Class(),AliAnalysisManager::kInputContainer); } 
  
  if(SP) {
    TString outputSP = "outputSPanalysis";
    outputSP+= type;
    outputSP+= ".root";
    AliAnalysisDataContainer *coutputSP = 
      mgr->CreateContainer("cobjSP", TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
  }
  
  if(LYZ1) {
    TString outputLYZ1 = "outputLYZ1analysis";
    outputLYZ1+= type;
    outputLYZ1+= "_firstrun.root";
    AliAnalysisDataContainer *coutputLYZ1 = 
      mgr->CreateContainer("cobjLYZ1", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1);
  }
  
  if(LYZ2) {
    TString outputLYZ2 = "outputLYZ2analysis";
    outputLYZ2+= type;
    outputLYZ2+= "_secondrun.root";
    AliAnalysisDataContainer *coutputLYZ2 = 
      mgr->CreateContainer("cobjLYZ2", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2);
  }
  
  if(LYZEP) {
    TString outputLYZEP = "outputLYZEPanalysis";
    outputLYZEP+= type;
    outputLYZEP+= ".root";
    AliAnalysisDataContainer *coutputLYZEP = 
      mgr->CreateContainer("cobjLYZEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZEP);
  }
  
  if(GFC) {
    TString outputGFC = "outputGFCanalysis";
    outputGFC+= type;
    outputGFC+= ".root";
    AliAnalysisDataContainer *coutputGFC = 
      mgr->CreateContainer("cobjGFC", TList::Class(),AliAnalysisManager::kOutputContainer,outputGFC);
  }
  
  if(QC) {
    TString outputQC = "outputQCanalysis";
    outputQC+= type;
    outputQC+= ".root";
    AliAnalysisDataContainer *coutputQC = 
      mgr->CreateContainer("cobjQC", TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
  }
  
  if(FQD) {
    TString outputFQD = "outputFQDanalysis";
    outputFQD+= type;
    outputFQD+= ".root";
    AliAnalysisDataContainer *coutputFQD = 
      mgr->CreateContainer("cobjFQD", TList::Class(),AliAnalysisManager::kOutputContainer,outputFQD);
  } 

  if(MCEP) {
    TString outputMCEP = "outputMCEPanalysis";
    outputMCEP+= type;
    outputMCEP+= ".root";
    AliAnalysisDataContainer *coutputMCEP = 
      mgr->CreateContainer("cobjMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputMCEP);
  }
  
  if (QA) { 
    if(SP) {
      TString qaNameIntSP = "QAforInt_SP_";
      qaNameIntSP += type;
      qaNameIntSP += ".root";
      AliAnalysisDataContainer *coutputQA1SP = 
	mgr->CreateContainer("QAintSP", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntSP);
      
      TString qaNameDiffSP = "QAforDiff_SP_";
      qaNameDiffSP += type;
      qaNameDiffSP += ".root";
      AliAnalysisDataContainer *coutputQA2SP = 
	mgr->CreateContainer("QAdiffSP", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffSP);
    }
    if(LYZ1) {
      TString qaNameIntLYZ1 = "QAforInt_LYZ1_";
      qaNameIntLYZ1 += type;
      qaNameIntLYZ1 += ".root";
      AliAnalysisDataContainer *coutputQA1LYZ1 = 
	mgr->CreateContainer("QAintLYZ1", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntLYZ1);
      
      TString qaNameDiffLYZ1 = "QAforDiff_LYZ1_";
      qaNameDiffLYZ1 += type;
      qaNameDiffLYZ1 += ".root";
      AliAnalysisDataContainer *coutputQA2LYZ1 = 
	mgr->CreateContainer("QAdiffLYZ1", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffLYZ1);
    }
    if(LYZ2) {
      TString qaNameIntLYZ2 = "QAforInt_LYZ2_";
      qaNameIntLYZ2 += type;
      qaNameIntLYZ2 += ".root";
      AliAnalysisDataContainer *coutputQA1LYZ2 = 
	mgr->CreateContainer("QAintLYZ2", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntLYZ2);
      
      TString qaNameDiffLYZ2 = "QAforDiff_LYZ2_";
      qaNameDiffLYZ2 += type;
      qaNameDiffLYZ2 += ".root";
      AliAnalysisDataContainer *coutputQA2LYZ2 = 
	mgr->CreateContainer("QAdiffLYZ2", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffLYZ2);
    }
    if(LYZEP) {
      TString qaNameIntLYZEP = "QAforInt_LYZEP_";
      qaNameIntLYZEP += type;
      qaNameIntLYZEP += ".root";
      AliAnalysisDataContainer *coutputQA1LYZEP = 
	mgr->CreateContainer("QAintLYZEP", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntLYZEP);
      
      TString qaNameDiffLYZEP = "QAforDiff_LYZEP_";
      qaNameDiffLYZEP += type;
      qaNameDiffLYZEP += ".root";
      AliAnalysisDataContainer *coutputQA2LYZEP = 
	mgr->CreateContainer("QAdiffLYZEP", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffLYZEP);
    }
    if(GFC) { 
      TString qaNameIntGFC = "QAforInt_GFC_";
      qaNameIntGFC += type;
      qaNameIntGFC += ".root";
      AliAnalysisDataContainer *coutputQA1GFC = 
	mgr->CreateContainer("QAintGFC", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntGFC);
      
      TString qaNameDiffGFC = "QAforDiff_GFC_";
      qaNameDiffGFC += type;
      qaNameDiffGFC += ".root";
      AliAnalysisDataContainer *coutputQA2GFC = 
	mgr->CreateContainer("QAdiffGFC", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffGFC);
    }
    if(QC) { 
      TString qaNameIntQC = "QAforInt_QC_";
      qaNameIntQC += type;
      qaNameIntQC += ".root";
      AliAnalysisDataContainer *coutputQA1QC = 
	mgr->CreateContainer("QAintQC", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntQC);
      
      TString qaNameDiffQC = "QAforDiff_QC_";
      qaNameDiffQC += type;
      qaNameDiffQC += ".root";
      AliAnalysisDataContainer *coutputQA2QC = 
	mgr->CreateContainer("QAdiffQC", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffQC);
    }
    if(FQD) { 
      TString qaNameIntFQD = "QAforInt_FQD_";
      qaNameIntFQD += type;
      qaNameIntFQD += ".root";
      AliAnalysisDataContainer *coutputQA1FQD = 
	mgr->CreateContainer("QAintFQD", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntFQD);
      
      TString qaNameDiffFQD = "QAforDiff_FQD_";
      qaNameDiffFQD += type;
      qaNameDiffFQD += ".root";
      AliAnalysisDataContainer *coutputQA2FQD = 
	mgr->CreateContainer("QAdiffFQD", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffFQD);
    }
    if(MCEP) {
      TString qaNameIntMCEP = "QAforInt_MCEP_";
      qaNameIntMCEP += type;
      qaNameIntMCEP += ".root";
      AliAnalysisDataContainer *coutputQA1MCEP = 
	mgr->CreateContainer("QAintMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameIntMCEP);
      
      TString qaNameDiffMCEP = "QAforDiff_MCEP_";
      qaNameDiffMCEP += type;
      qaNameDiffMCEP += ".root";
      AliAnalysisDataContainer *coutputQA2MCEP = 
	mgr->CreateContainer("QAdiffMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameDiffMCEP);
    }
  }
  
  //____________________________________________//
  
  if (SP)   { 
    mgr->ConnectInput(taskSP,0,cinput1); 
    mgr->ConnectOutput(taskSP,0,coutputSP);
    if (QA) { mgr->ConnectOutput(taskSP,1,coutputQA1SP);
    mgr->ConnectOutput(taskSP,2,coutputQA2SP); }
  } 
  if (LYZ1) { 
    mgr->ConnectInput(taskLYZ1,0,cinput1); 
    mgr->ConnectOutput(taskLYZ1,0,coutputLYZ1);
    if (QA) { mgr->ConnectOutput(taskLYZ1,1,coutputQA1LYZ1);
    mgr->ConnectOutput(taskLYZ1,2,coutputQA2LYZ1); }
  }  
  if (LYZ2) { 
    mgr->ConnectInput(taskLYZ2,0,cinput1); 
    mgr->ConnectInput(taskLYZ2,1,cinputLYZ2);
    mgr->ConnectOutput(taskLYZ2,0,coutputLYZ2);
    if (QA) { mgr->ConnectOutput(taskLYZ2,1,coutputQA1LYZ2);
    mgr->ConnectOutput(taskLYZ2,2,coutputQA2LYZ2); }
    cinputLYZ2->SetData(fInputListLYZ2);
  }  
  if (LYZEP) { 
    mgr->ConnectInput(taskLYZEP,0,cinput1); 
    mgr->ConnectInput(taskLYZEP,1,cinputLYZEP);
    mgr->ConnectOutput(taskLYZEP,0,coutputLYZEP);
    if (QA) { mgr->ConnectOutput(taskLYZEP,1,coutputQA1LYZEP);
    mgr->ConnectOutput(taskLYZEP,2,coutputQA2LYZEP); }
    cinputLYZEP->SetData(fInputListLYZEP);
  }
  if (GFC)   { 
    mgr->ConnectInput(taskGFC,0,cinput1); 
    mgr->ConnectOutput(taskGFC,0,coutputGFC);
    if (QA) { mgr->ConnectOutput(taskGFC,1,coutputQA1GFC);
    mgr->ConnectOutput(taskGFC,2,coutputQA2GFC); }
  }  
  if (QC)   { 
    mgr->ConnectInput(taskQC,0,cinput1); 
    mgr->ConnectOutput(taskQC,0,coutputQC);
    if (QA) { mgr->ConnectOutput(taskQC,1,coutputQA1QC);
    mgr->ConnectOutput(taskQC,2,coutputQA2QC); }
  }
  if (FQD)   { 
    mgr->ConnectInput(taskFQD,0,cinput1); 
    mgr->ConnectOutput(taskFQD,0,coutputFQD);
    if (QA) { mgr->ConnectOutput(taskFQD,1,coutputQA1FQD);
    mgr->ConnectOutput(taskFQD,2,coutputQA2FQD); }
  }    
  if (MCEP)  { 
    mgr->ConnectInput(taskMCEP,0,cinput1); 
    mgr->ConnectOutput(taskMCEP,0,coutputMCEP);
    if (QA) { mgr->ConnectOutput(taskMCEP,1,coutputQA1MCEP);
    mgr->ConnectOutput(taskMCEP,2,coutputQA2MCEP); }
  }  
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
  
  timer.Stop();
  timer.Print();
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
	  cerr<<presentDirName<<endl;
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
	  cerr<<presentDirName<<endl;
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



void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target
  
  chain->Lookup();
  
  TObjArray* list = chain->GetListOfFiles();
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  ofstream outfile;
  outfile.open(target);
  
  while ((obj = iter->Next()))
    outfile << obj->GetTitle() << "#AliESDs.root" << endl;
  
  outfile.close();
  
  delete iter;
}
