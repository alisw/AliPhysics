#include "AliPSQA.h"

void PSQA() {
  //
  // Produce QA  output of physics selection.
  // Opens selected files for QA, then processes the historgrams necessary, 
  // prints results to a root file and/or text files specified by user.
  //
  //  Please read the README in case the documentation provided in this macro, .h, and .cxx files is not enough!
  //
  //gROOT->LoadMacro("AliPSQA.cxx++g");
  AliPSQA * QAobject = new AliPSQA;  // Initialize relevant object
  //
  // configure relevant parameters
  //
  Bool_t runFromLocalInputFiles = kTRUE; // if kTRUE, run over either cache or test, otherwise use alien files
  Bool_t saveTxtFiles           = kTRUE; 
  Bool_t saveRootFile           = kTRUE;
  Bool_t drawGraphs             = kTRUE;
  //
  // settings
  //
  QAobject->InitializeRuns("../QAoutputPerPeriod/12d_Pass1/GoodRunList.list"); 
  QAobject->SetLocalPath("../InputFilesFromGridPerPeriod/12d_Pass1/"); // Location of files to be processed locally (whether test or cache)
  QAobject->SetTxtOutDirectory("../QAoutputPerPeriod/12d_Pass1/"); // Path for macro output text files
  QAobject->SetRootOutDirectory("../QAoutputPerPeriod/12d_Pass1/"); // Path for macro output root files
  QAobject->SetOutRootName("LHC12d_PSQA.root");
  //
  QAobject->SetPSInput(""); // Subdirectory for files to be processed, which pass to use
  QAobject->SetROOTInput("event_stat.root"); //root file to be processed
  QAobject->InitializeTriggers("triggerNames.list","triggerMaskChannels.list"); // trig names to be used, plus trig channels (fast vs regular)
  QAobject->InitializePlots("plots.list");     //Initialize macro with QA plots to be used
  QAobject->SetTxtFileName("PSQA_"); // This is appended with the run name and .txt
  //
  // enable settings
  //
  QAobject->SetLocalMode(runFromLocalInputFiles); 
  QAobject->SetSaveQAValToTxtFile(saveTxtFiles);
  QAobject->SetSaveQAToROOTFile(saveRootFile);
  //
  // run and extract    
  //
  Bool_t goodQACheck = QAobject->ComputeQAPerRun();
  //
  // Print a list of the runs that failed, with some hint as to why they failed. Implementation file has comments.
  //
  QAobject->GetBadFiles()->Print();
  //
  if (goodQACheck == kFALSE){return;} // No plots saved at all, skip altogether the option to plot all canvases.
  //
  QAobject->SetDrawAllTGraphs(drawGraphs); // If kTRUE, will draw all TGraphErrors saved in the ROOT file.
  if (QAobject->GetDrawAllTGraphs()){
    QAobject->DrawAllTGraphs();
  }

};

