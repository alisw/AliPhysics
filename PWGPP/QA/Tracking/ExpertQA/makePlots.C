/*
  .L $ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/makePlots.C

*/

void makePlots(const char *inputFile, Int_t mode=0)
{
  if (mode==1) {
    return  testMakeDCArPullFitsMI(inputFile);
  }
  
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+");
  if (gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++")<0){
    gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+");
  }
  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
  a->InitAnalysis(inputFile );
  a->Loop();
  
  
}



void testMakeDCArPullFitsMI()(const char *inputFile){
  //
  //
  //
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+");
  if (gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++")<0){
    gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+");
  }
  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
  a->InitAnalysis("");
  a->MakeDCArPullFitsMI();
  (*(a->GetStreamer()))<<"trending"<<"\n";

}




void testQPt(const char *inputFile){
  //
  //
  //
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+");
  if (gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++")<0){
    gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+");
  }
  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
  a->InitAnalysis("");
  a->MakePowerFit(-1);
}

