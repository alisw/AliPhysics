void RunIPTask(const char* mode)
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT "
			  "-I$ALICE_ROOT/include -I$ALICE_ROOT/ITS "
			  "-I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS "
			  "-I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros "
			  "-I$ALICE_ROOT/ANALYSIS -g"); 
  //
  // Load analysis libraries
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWG1.so");
  //
  TChain *chainESD = 0;
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(mode,-1);
  //chain->SetBranchStatus("*ESDfriend*",1);
  //
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  AliESDInputHandler *esdH = new AliESDInputHandler();
  esdH->SetActiveBranches("ESDfriend");
  //
  mgr->SetInputEventHandler(esdH);
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskIntSpotESD.C");
  AliAnalysisTaskIPInfo* iptask = AddTaskIntSpotESD();

  if(!mgr->InitAnalysis()) return;
  //
  mgr->StartAnalysis("local",chain);
}
