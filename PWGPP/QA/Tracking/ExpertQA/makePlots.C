/*
  .L $ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/makePlots.C

*/

void makePlots(const char *inputFile)
{
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  if (gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++")<0){
    gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  }
  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
  a->InitAnalysis(inputFile );
  a->Loop();
  
  
}

void testQPt(const char *inputFile){
  //
  //
  //
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  if (gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++")<0){
    gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  }
  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
  a->InitAnalysis("");
  MakePowerFit(-1);
  
}

