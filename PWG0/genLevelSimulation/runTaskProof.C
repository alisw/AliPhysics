void runTaskProof(const char * dataset, const char * datasetpath="/COMMON/COMMON/", const char * outdir="") {
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("alice-caf","workers=1x");
  gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-33-AN");
  
  //  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/");
  // gSystem->AddIncludePath("-I${ALICE_ROOT}/include/");
  // gProof->UploadPackage("$ALICE_ROOT/STEERBase");
  // gProof->EnablePackage("$ALICE_ROOT/STEERBase");
  // gProof->UploadPackage("$ALICE_ROOT/ESD");
  // gProof->EnablePackage("$ALICE_ROOT/ESD");
  // gProof->UploadPackage("$ALICE_ROOT/AOD");
  // gProof->EnablePackage("$ALICE_ROOT/AOD");
  // gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
  // gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
  // gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
  // gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
 

    // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //AliVEventHandler *esdH = new AliESDInputHandler();
  AliESDInputHandler* esdH = new AliESDInputHandler();
  //esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);
	
  AliMCEventHandler *mc = new AliMCEventHandler();
  mc->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(mc);

  // assign simple task
  gProof->Load(gSystem->ExpandPathName("$(ALICE_ROOT)/PWG0/genLevelSimulation/AliAnalysisTaskdNdetaMC.cxx+g"));
  AliAnalysisTask *task = new AliAnalysisTaskdNdetaMC("TaskdNdeta");
  mgr->AddTask(task);

  char outfname[1000];
  sprintf(outfname, "dndeta_%s.root",dataset);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1",TList::Class(),
							    AliAnalysisManager::kOutputContainer,outfname);

  mgr->ConnectInput(task,0,cinput1);
  mgr->ConnectOutput(task,1,coutput1);
	
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",Form("%s%s#TE",datasetpath,dataset),5000);
  // mgr->StartAnalysis("proof",Form("%s%s#TE",datasetpath,dataset),5000);
  //  mgr->StartAnalysis("proof","/COMMON/COMMON/LHC09b14_7TeV_0.5T_Phojet#esdTree");

  if (!strcmp(outdir,"")){
    gSystem->Exec(Form("mv %s %s", outfname, outdir));
  }

}
