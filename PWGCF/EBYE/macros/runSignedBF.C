void runSignedBF(const char* directoryFile = "/dcache/alice/panosch/alice/data/2015/LHC15o/000244917/pass2_lowIR/",  const char* runtype = "local" // local, proof or grid
		 ) {
  //Author: Panos.Christakoglou@cern.ch
  gROOT->ProcessLine(".include $ROOTSYS/include"); 
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    
  gSystem->Load("libCore");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libEventMixing");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisMyTask"); 
  
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  TChain *chain = new TChain("aodTree");  
  TSystemDirectory dir(directoryFile ,directoryFile);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString filename;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      filename = directoryFile;
      filename += file->GetName();
      //filename = directoryFile;
      filename += "/AliAOD.root";
      cout << "Adding: "<<filename.Data() << endl;
      chain->Add(filename.Data());
    }
  }

  //=========Add the Qn vector framework task==================//
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrections.C");
  AddTaskFlowQnVectorCorrections();

  
  //=========On the fly compilation (testing mode)==================//
  gROOT->LoadMacro("AliAnalysisTaskSignedBF.cxx++g"); 
  
  //=========Add the signed BF task==================//
  gROOT->LoadMacro("AddTaskSignedBF.C"); 
  AliAnalysisTaskSignedBF *taskBF = AddTaskSignedBF();
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  mgr->StartAnalysis("local", chain); 
  
}
