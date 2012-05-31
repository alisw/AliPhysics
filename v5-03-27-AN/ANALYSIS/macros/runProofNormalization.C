#define TENDER

void runProofNormalization(const char * dataset = "LHC09b12_7TeV_0.5T", TString dataSetPath ="/PWG0/jgrosseo/",const char * filename = "LHC09b12_7TeV_0.5T_norm.root", Bool_t isMC = 1,Int_t nev =123456789) {

#ifdef TENDER
  TGrid::Connect("alien://");
#endif
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("alice-caf","workers=1x");// limit the number of workers
  //  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
  //  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/");
  //  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/");
    // Use a precompiled tag
    TString alirootMode="";    // STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice (default aliroot mode)
    //alirootMode="ALIROOT";     // $ALICE_ROOT/macros/loadlibs.C
    //  alirootMode="REC";     // $ALICE_ROOT/macros/loadlibsrec.C
    //  alirootMode="SIM";     // $ALICE_ROOT/macros/loadlibssim.C
    //  alirootMode="TRAIN";   // $ALICE_ROOT/macros/loadlibstrain.C (not working yet)
    //  alirootMode="CUSTOM";  // nothing is loaded, but aliroot variables are set (not working yet)
 
    TString extraLibs;
    extraLibs= ""; // not needed in default aliroot mode
    extraLibs+="CDB:RAWDatabase:STEER:TENDER:TRDbase:STAT:TRDrec:VZERObase:VZEROsim:VZEROrec:RAWDatarec:TPCbase:TPCrec:TPCcalib:TENDERSupplies:RAWDatabase:RAWDatarec:RAWDatasim:TOFbase:TOFrec";
    TList *list = new TList();
    // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
    list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
    // sets $ALIROOT_EXTRA_LIBS on each worker to let proof to know to load extra libs
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
#ifdef TENDER
    list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
#endif
    // connect to proof
    gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-22-AN", list);
    //    gProof->Exec("TGrid::Connect(\"alien://\");");
  // gProof->UploadPackage("$ALICE_ROOT/obj/STEERBase");
  // gProof->EnablePackage("$ALICE_ROOT/obj/STEERBase");
  // gProof->UploadPackage("$ALICE_ROOT/obj/ESD");
  // gProof->EnablePackage("$ALICE_ROOT/obj/ESD");
  // gProof->UploadPackage("$ALICE_ROOT/obj/AOD");
  // gProof->EnablePackage("$ALICE_ROOT/obj/AOD");
  // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSIS");
  // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSIS");
  // gProof->UploadPackage("$ALICE_ROOT/obj/OADB");
  // gProof->EnablePackage("$ALICE_ROOT/obj/OADB");
  // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSISalice");
  // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSISalice");
//   gProof->UploadPackage("STEERBase.par");
//   gProof->EnablePackage("STEERBase");
//   gProof->UploadPackage("ESD.par");
//   gProof->EnablePackage("ESD");
//   gProof->UploadPackage("AOD.par");
//   gProof->EnablePackage("AOD");
//   gProof->UploadPackage("ANALYSIS.par");
//   gProof->EnablePackage("ANALYSIS");
//   gProof->UploadPackage("ANALYSISalice.par");
//   gProof->EnablePackage("ANALYSISalice");
//   gProof->UploadPackage("CORRFW.par");
//   gProof->EnablePackage("CORRFW"); 

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //  mgr->SetDebugLevel(3);
  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler; 

  mgr->SetInputEventHandler(esdH);
	
  if(isMC) {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mc->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mc);
  }
  //____________________________________________//

#ifdef TENDER
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/TenderSupplies/AddTaskTender.C");
  AliAnalysisTask* tender=0x0;
  if(!isMC)
    {
      tender = AddTaskTender(kTRUE);
      // tender->SetDebugLevel(10);
    }
  else
    {
      tender = AddTaskTender(kFALSE);
      // tender->SetDebugLevel(10);
    }
#endif

  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,1,!isMC); // Use Physics Selection. Enable computation of BG if is not MC
  //  task->SelectCollisionCandidates(); /// This should be disabled, at least for MC: we need all the events
  physSelTask->GetPhysicsSelection()->SetBin0Callback("TaskNormalization");

  // assign simple task
  AliCollisionNormalizationTask * task = new AliCollisionNormalizationTask("TaskNormalization");
  //  task->SetMC();
  task->SetMC(isMC);
  mgr->AddTask(task);



  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();	
  mgr->ConnectInput(task,0,cinput1);


  
  // Attach output
  cOutput = mgr->CreateContainer("Norm", TList::Class(), AliAnalysisManager::kOutputContainer,filename);
  mgr->ConnectOutput(task, 1, cOutput);      
	
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",dataSetPath+dataset+"#esdTree",nev);

}
