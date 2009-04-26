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

// Analysis type can be ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";

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

// Boolean to fill/not fill the QA histograms
Bool_t QA = kFALSE;   

// Weights 
//Use weights for Q vector
Bool_t usePhiWeights = kFALSE; //Phi
Bool_t usePtWeights  = kFALSE; //v'(pt)
Bool_t useEtaWeights = kFALSE; //v'(eta)
Bool_t useWeights = usePhiWeights||usePtWeights||useEtaWeights;


void runFlowTask(Int_t mode=mLocal, Int_t nRuns = -1, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset=0)
{

  TStopwatch timer;
  timer.Start();

  if (LYZ1 && LYZ2) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(); }
  if (LYZ2 && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
  if (LYZ1 && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }

  LoadLibraries(mode);

  if (mode==mLocal || mode == mLocalPAR || mode == mGRID) {
    if (type!="AOD") {TChain* chain = CreateESDChain(dataDir, nRuns, offset);}
    else { TChain* chain = CreateAODChain(dataDir, nRuns, offset);}
  }

  //____________________________________________//
  //external input files
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
    taskFE->SetAnalysisType(type);
    mgr->AddTask(taskFE);
  }
  else { 
    taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",kFALSE);
    taskFE->SetAnalysisType(type);
    mgr->AddTask(taskFE);
  }
  if (FQD){
    AliAnalysisTaskFittingQDistribution *taskFQD = new AliAnalysisTaskFittingQDistribution("TaskFittingQDistribution",kFALSE);
    taskFQD->SetUsePhiWeights(usePhiWeights); 
    mgr->AddTask(taskFQD);
  }
  if (SP){
    AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct");
    mgr->AddTask(taskSP);
  }
  if (LYZ1){
    AliAnalysisTaskLeeYangZeros *taskLYZ1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE);
    taskLYZ1->SetFirstRunLYZ(kTRUE);
    taskLYZ1->SetUseSumLYZ(kTRUE);
    mgr->AddTask(taskLYZ1);
  }
  if (LYZ2){
    AliAnalysisTaskLeeYangZeros *taskLYZ2 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE);
    taskLYZ2->SetFirstRunLYZ(kFALSE);
    taskLYZ2->SetUseSumLYZ(kTRUE);
    mgr->AddTask(taskLYZ2);
  }
  if (LYZEP){
    AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane");
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
    AliAnalysisTaskQCumulants *taskQC = new AliAnalysisTaskQCumulants("TaskQCumulants",kFALSE);
    taskQC->SetUsePhiWeights(usePhiWeights); 
    taskQC->SetUsePtWeights(usePtWeights);
    taskQC->SetUseEtaWeights(useEtaWeights); 
    mgr->AddTask(taskQC);
  }
  if (MCEP){
    AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane");
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
    mgr->ConnectInput(taskLYZ1,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1,0,coutputLYZ1);
  }  
  if (LYZ2) { 
    mgr->ConnectInput(taskLYZ2,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2,1,cinputLYZ2);
    mgr->ConnectOutput(taskLYZ2,0,coutputLYZ2);
    cinputLYZ2->SetData(fInputListLYZ2);
  }  
  if (LYZEP) { 
    mgr->ConnectInput(taskLYZEP,0,coutputFE); 
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
    mgr->ConnectInput(taskMCEP,0,coutputFE); 
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
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
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

   gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/PWG2/FLOW/AliFlowTasks");
   gROOT->LoadMacro("ConfigFlowAnalysis.C++");
   

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

   gSystem->AddIncludePath("-I./PWG2FlowTasks -I./CORRFW");
   gROOT->LoadMacro("ConfigFlowAnalysis.C++");

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
   //TProof::Open("abilandz@alicecaf.cern.ch");
   //TProof::Open("nkolk@alicecaf.cern.ch");
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

   gSystem->AddIncludePath("-I./PWG2flowTasks -I./CORRFW");
   gProof->Load("ConfigFlowAnalysis.C++");
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


