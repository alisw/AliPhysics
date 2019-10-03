void RunIPTask(const char* mode)
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS "
			  "-I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS "
			  "-I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS "
			  "-I$ALICE_PHYSICS/STEER -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros "
			  "-I$ALICE_PHYSICS/ANALYSIS -g"); 
  //
  // Load analysis libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  //
  TChain *chainESD = 0;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(mode,-1);
  //chain->SetBranchStatus("*ESDfriend*",1);
  //
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  AliESDInputHandler *esdH = new AliESDInputHandler();
  esdH->SetActiveBranches("ESDfriend");
  //
  mgr->SetInputEventHandler(esdH);
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskIntSpotESD.C");
  AliAnalysisTaskIPInfo* iptask = AddTaskIntSpotESD();

  if(!mgr->InitAnalysis()) return;
  //
  mgr->StartAnalysis("local",chain);
}
