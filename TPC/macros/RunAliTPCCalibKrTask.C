/// \file RunAliTPCCalibKrTask.C
///
/// Example usage:
///
/// 0. Load neccessary libraries
///
/// ~~~{.cpp}
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libANALYSISalice");
/// gSystem->Load("libTPCcalib");
/// gSystem->Load("libXrdClient");
/// gSystem->Load("libNetx");
/// TGrid::Connect("alien://",0,0,"t");
/// gSystem->Load("$ROOTSYS/lib/libXrdClient");
/// ~~~
///
/// 1. Make list of the files
///
/// ~~~{.cpp}
/// .L $ALICE_ROOT/TPC/macros/testTPC/AlienToolkit.cxx+
/// gSystem->Load("libXrdClient");
/// gSystem->Load("libNetx");
/// AlienToolkit toolkit;
/// char *path = "/alice/cern.ch/user/a/amatyja/alice/data/"
/// toolkit.MakeCollection(path,"32129*Krypton.root");
/// toolkit.MakeCollection(path,"32231*Krypton.root");
/// toolkit.MakeCollection(path,"32284*Krypton.root");
/// toolkit.PrintPFN(); > list.txt
/// ~~~
///
/// 2. Initialization of proof
///
/// ~~~{.cpp}
/// TProofMgr * proofmgr = TProof::Mgr("lxgrid5.gsi.de");
/// TProof * proof = proofmgr->CreateSession();
/// proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)1000);
/// .L /u/miranov/macros/ProofEnableAliRoot.C
/// ProofEnableAliRoot("/usr/local/grid/AliRoot/HEAD0108");
/// gProof->Exec("gSystem->Load(\"libANALYSIS\")",kTRUE);
/// gProof->Exec("gSystem->Load(\"libSTAT\")",kTRUE);
/// gProof->Exec("gSystem->Load(\"libTPCcalib\")",kTRUE);
/// gProof->Exec(".x $ALICE_ROOT/TPC/macros/ConfigOCDB.C");
/// ~~~
///
/// 3. Run analysis on PROOF
///
/// ~~~{.cpp}
/// // Create chain of input files
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
/// gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
/// .L $ALICE_ROOT/TPC/macros/RunAliTPCCalibKrTask.C
/// RunAliTPCCalibKrTask(kTRUE);
///
/// // Check the cuts for clusters
/// AliXRDPROOFtoolkit tool;
/// TChain * chain = tool.MakeChain("list.txt","Kr",0,1000);
/// chain->Lookup();
/// chain->SetProof(kTRUE);
///
/// TCut cutR0("cutR0","fADCcluster/fSize<200");        // adjust it according v seetings -
/// TCut cutR1("cutR1","fADCcluster/fSize>7");          // cosmic tracks and noise removal
/// TCut cutR2("cutR2","fMax.fAdc/fADCcluster<0.4");    // digital noise removal
/// TCut cutR3("cutR3","fMax.fAdc/fADCcluster>0.01");   // noise removal
/// TCut cutR4("cutR4","fMax.fTime>200");   // noise removal
/// TCut cutR5("cutR5","fMax.fTime<600");   // noise removal
/// TCut cutS1("cutS1","fSize<200");    // adjust it according v seetings - cosmic tracks
/// TCut cutAll = cutR0+cutR1+cutR2+cutR3+cutR4+cutR5+cutS1;
///
/// // Example usage
/// TFile f("KrHisto.root");
/// AliTPCCalibKr *kr = f.Get("AliTPCCalibKr");
///
/// kr->ProjectHisto(kr->GetHistoKr(71),"aaa",30,36,30,40)->Draw()
///
/// MakeTree();
///
/// //default cuts
/// TCut cutKr("cutKr","entries.fElements<5000&&fitNormChi2.fElements<3&&fitNormChi2.fElements>0.2&&abs(fitRMS.fElements/fitMean.fElements-0.06)<0.025");
///
/// TObjArray * array = AliTPCCalibViewerGUI::ShowGUI("kryptonTree.root");
/// AliTPCCalibViewerGUI * viewer = (AliTPCCalibViewerGUI*)array->At(0);
/// TTree * tree = viewer->GetViewer()->GetTree();
///
/// tree->SetAlias("cutAll","abs(fitNormChi2.fElements-2.)<1.8&&entries.fElements/entries_Median.fElements<4&&entries.fElements/entries_Median.fElements>0.4&&fitRMS.fElements/fitMean.fElements<0.09&&fitRMS.fElements/fitMean.fElements>0.02")
/// ~~~

TChain * chain = 0;

void RunAliTPCCalibKrTask(Bool_t bProof = kFALSE)
{
  ///

  AliXRDPROOFtoolkit tool;
  chain = tool.MakeChain("list.txt","Kr","",20000,0);
  chain->Lookup();
  chain->SetBranchStatus("Cl.fCluster",kFALSE);

  //
  // Create the analysis manager
  //
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Calibration component
  AliTPCCalibKr *calibObj = new AliTPCCalibKr;
  //calibObj->SetASide(kFALSE);

  // Add task
  AliTPCCalibKrTask *task = new AliTPCCalibKrTask;
  task->SetInputChain(chain);
  task->SetTPCCalibKr(calibObj);
  mgr->AddTask(task);

  // Attach input
  cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer,"outHistFile.root");
  mgr->ConnectOutput(task, 0, cOutput);
  //
  cOutput->SetSpecialOutput(kTRUE);
  cOutput->SetFileName("CalibObjectFile.root");
  //
  // Run analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->SetDebugLevel(1);

  if(bProof) {
    mgr->StartAnalysis("proof", chain);
  }
  else mgr->StartAnalysis("local", chain);
}



void MakeTree(){

  TFile fpad("calibKr.root");
  AliTPCPreprocessorOnline * preprocesor = new AliTPCPreprocessorOnline;
  preprocesor->AddComponent(spectrMean->Clone());
  preprocesor->AddComponent(spectrRMS->Clone());
  preprocesor->AddComponent(fitMean->Clone());
  preprocesor->AddComponent(fitRMS->Clone());
  preprocesor->AddComponent(fitNormChi2->Clone());
  preprocesor->AddComponent(entries->Clone());
  preprocesor->DumpToFile("kryptonTree.root");

}
