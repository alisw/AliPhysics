/*
  Macros to run PWGPP train:
  
  void runPWGPPTrain(const char *macros = "AddTask*.C", const char *fname="AliESDs.root");
  //Parameters:
  //    macros: run the train for selected tasks
  //            tasks are expected to be in the working  directory
  //    fname : name of the input file or input list


  .L $ALICE_PHYSICS/PWGPP/PWGPPmacros/runPWGPPTrain.C


*/

void LoadTrainLibs(){
  //
  // load libraries needed for train
  //
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/PWGPP/TRD/macros/ ");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWG2");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libPWGmuondep");
  gSystem->Load("libPWG2forward");
  gSystem->Load("libPWG4PartCorrBase");
  gSystem->Load("libPWG4PartCorrDep");
}

void TestConfig(){

}

void AddMacros(const char *addmacros){
  //
  // add tasks from selected macros - * convention can be used
  // Macro expected to be in the working directory
  // Macros has to be without arguments
  //
  //
  TString  macroList = gSystem->GetFromPipe(Form("cat ConfigTask.txt |grep  %s",addmacros));  
  TObjArray * array  = macroList.Tokenize("\n");

  if (!array) { printf("No task specified"); return;}
  if (array->GetEntries()==0) {
    printf("Empty list of tasks"); 
    return;
  }
  for (Int_t itask=0; itask<array->GetEntries(); itask++){
    gROOT->Macro(array->At(itask)->GetName());    
  }
}

void AddToChain(TChain * chain, TString inputList){
  //
  // add the files form inputList into chain
  //
  ifstream in;
  in.open(inputList.Data());
  Int_t counter=0; 
  TString currentFile;
  while(in.good()) {
    in >> currentFile;
    if (!currentFile.Contains(".root")) continue;
    chain->Add(currentFile.Data());
    printf("%d\t%s\n",counter,currentFile.Data()); 
    counter++;
  }
}

void PrintSysInfo(){
  //
  // print sysinfo
  //
  TF1 f1("f1","pol1");
  TTree * tree = AliSysInfo::MakeTree("syswatch.log");
  Int_t entries=0;
  TGraph *gr = 0;
  entries=tree->Draw("VM:id0","id0>10","goff");
  gr=new TGraph(entries, tree->GetV2(),tree->GetV1());
  gr->Fit(&f1);
  gr->Draw("alp");
  printf("SysInfoMem:\t%f\n",f1->GetParameter(1));
  //
  entries=tree->Draw("T:id0","id0>10","goff");
  gr=new TGraph(entries, tree->GetV2(),tree->GetV1());
  gr->Fit(&f1);
  gr->Draw("alp");
  printf("SysInfoTime:\t%f\n",f1->GetParameter(1));
}

void runPWGPPTrain(const char *macros = "AddTask*.C", TString inputList ="esd.list", Int_t debugLevel=0) {
  //Parameters:
  //    macros: run the train for selected tasks
  //            tasks are expected to be in the working  directory
  //    fname : name of the input file or input list
  TStopwatch timer;
  timer.Start();
  LoadTrainLibs();

  //____________________________________________
  // Make the analysis manager
  //
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  mgr->SetDebugLevel(debugLevel);
  
  AliInputEventHandler* esdH = new AliESDInputHandlerRP();
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  
  //
  // make chain
  //
  TChain* chain = new TChain("esdTree");
  if (inputList.Contains(".root")){
    chain->AddFile(inputList);
  }
  if (inputList.Contains(".list")){
    AddToChain(chain, inputList);
  }
  chain->Lookup();
  //
  //
  // Init
  AddMacros(macros);
  mgr->SetNSysInfo(1000);
  if (!mgr->InitAnalysis()) 
    mgr->PrintStatus();
  mgr->PrintStatus();
  // Run on dataset
  //
  mgr->StartAnalysis("local", chain);
  PrintSysInfo();

  timer.Stop();
  timer.Print();
}

