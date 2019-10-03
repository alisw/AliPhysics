void runAnalysisTaskPtMothFromPtDaug(TH1F* histPtDaug=0x0, 
         Bool_t readKineFromNtupla=kTRUE,const char* galiceList="galice.txt")
 {
/////////////////////////////////////////////////////////////////////////////
//  Run-macro to extract pt-spectra (and ptMin-spectra) for mothers        //
//  input: 1) pT histogram of daughter particles                           //
//         2) boolean flag:kFALSE -> read Kinematics.root to evaluate      //
//                                 correction factors, create a TNtuple    //
//                                 with kinematic informations of mothers  //
//                                 and daughters and store it in the file  //
//                                 "DecayKine.root"                        //
//                         kTRUE  -> read the TNtupla from the file        //
//                                 "DecayKine.root" (after it is created)  //
//                                 to evaluate correction factors          //
//         3) name of file with the list of "galice.root" files to read    //
//            kinematics (not needed after the TNtupla is created)         //
//                                                                         //
//  output: 1) file Mothers.root which contains pt-spectra and ptMin       //
//             spectra of mothers particles                                //  
//          2) TNtupla with kinematic informations (optional)              //
//                                                                         //
//      Origin:  Giuseppe.Bruno@ba.infn.it, Fiorella.Fionda@ba.infn.it     //
//									   //	
/////////////////////////////////////////////////////////////////////////////
  
  char *ntuplaFileName = "DecayKine.root"; // default name of the Ntupla
  char *mode = "local"; // analysis mode (select local or proof) 
  char *dataset = "/COMMON/COMMON/LHC09a14_0.9TeV_0.5T"; // define dataset for proof

  if(mode == "proof") loadLib();   
  else{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGHFbase");
  }
  // Set AliPtMothFromPtDaugh object 
  AliPtMothFromPtDaugh *ptExtr = new AliPtMothFromPtDaugh();
  ptExtr->SetDefaultAnalysis(AliPtMothFromPtDaugh::kBtoJPSI);
  ptExtr->SetBinsPtMoth(0.,10,20,1);
  ptExtr->SetBinsPtMinMoth(0.,10,20,1);
  ptExtr->SetEtaMothers(-1.5,1.5);
  ptExtr->SetEtaDaughter(-1.,1.);
  if(!ptExtr->ReadHistoPtDaught(histPtDaug))
     { printf("Daughter pt-Histogram is not defined \n"); return; }
  if(!ptExtr->CreateWeights()) return;
  //
  // create Analysis manager with MC and Input handlers
  //
  AliAnalysisManager *mgr  = new AliAnalysisManager("mgr", "Analysis Manager");
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  AliInputEventHandler *inputHandler =0x0; 
  mgr->SetInputEventHandler(inputHandler);

  TChain *chain=0x0;
  if(mode=="local"){
  if(!readKineFromNtupla) chain = CreateChain(galiceList); // create a chain with galice.root files
  else  mgr->SetExternalLoop(kTRUE); // to abort event loop when Ntupla is already created
  }
  //
  // create task and add it to the manager
  //
  AliAnalysisTaskPtMothFromPtDaugh *task = new AliAnalysisTaskPtMothFromPtDaugh(readKineFromNtupla);
  task->SetPtMothFromPtDaugh(ptExtr); // set AliPtMothFromPtDaugh object to the task
  task->SetNtuplaFileName(ntuplaFileName);
  mgr->AddTask(task);  
   //
   // create input / output containers
   //
  AliAnalysisDataContainer *cOutput = mgr->CreateContainer("Mothers", TList::Class(), AliAnalysisManager::kOutputContainer,"Mothers.root");
  mgr->ConnectOutput(task, 1, cOutput);
  // optional output for TNtupla
  AliAnalysisDataContainer *cOutput1 = 0x0;
  if(!readKineFromNtupla){
  cOutput1 = mgr->CreateContainer("DecayKine", TNtuple::Class(), AliAnalysisManager::kOutputContainer,ntuplaFileName);
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,2,cOutput1);
  }
  //
  // run analysis locally
  //
  Int_t result=mgr->InitAnalysis();
  if (!result){
    Error("RunAnalysis","Could not initialise Analysis.");
    return 0;
  }
  mgr->PrintStatus();
  if(mode=="local") mgr->StartAnalysis(mode,chain);
  else if(mode=="proof") mgr->StartAnalysis(mode,dataset);
  if(readKineFromNtupla) mgr->Terminate();
  return;
 }

TChain *CreateChain(const char* galiceName)
 {
  //
  // Create a chain from a list of galice.root
  //
  ifstream in;
  in.open(galiceName);
  if (!in.is_open()) return 0x0;
  TChain *chain=new TChain("TE");
  TString line;
   while(in.good()) {
     in >> line;
     if (!line.IsNull()) chain->AddFile(line);
     }
  return chain;
 }

void loadLib(){
printf("****** Connect to PROOF *******\n");
gEnv->SetValue("XSec.GSI.DelegProxy","2");
TProof::Open("alicecaf.cern.ch");
gProof->UploadPackage("STEERBase.par");
gProof->EnablePackage("STEERBase.par");
gProof->UploadPackage("ESD.par");
gProof->EnablePackage("ESD.par");
gProof->UploadPackage("AOD.par");
gProof->EnablePackage("AOD.par");
gProof->UploadPackage("ANALYSIS.par");
gProof->EnablePackage("ANALYSIS.par");
gProof->UploadPackage("ANALYSISalice.par");
gProof->EnablePackage("ANALYSISalice.par");
gProof->UploadPackage("PWG3base.par");
gProof->EnablePackage("PWG3base.par");
}
