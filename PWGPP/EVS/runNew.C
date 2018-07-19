#ifndef __CINT__
#include "Includes.h"
#endif

void runNew(TString fileList="runlist.txt"){
  TString runlist = gSystem->GetFromPipe("ls all");
  ifstream f(fileList.Data());
  Int_t run;
  gSystem->Unlink("trending_merged.root");
  while (!f.eof()){
    f >> run;
    if (runlist.Contains(Form("%i",run))) continue;
    gSystem->Exec(Form("aliroot -l -b -q 'runLevelEventStatQA.C\(\"%s\",%i,\"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2017/OCDB\"\)'","",run));
    gSystem->mkdir(Form("all/%i",run));
//    gSystem->Exec(Form("mv trending.root all/%i/trending.root",run)); // Old solution
    gSystem->Rename("trending.root",Form("all/%i/trending.root",run));
  }
  gSystem->Exec("hadd -f trending_merged.root all/*/trending.root");
  gSystem->Exec("aliroot -l -b -q periodLevelQA.C");
  gSystem->Exec("root -l -b -q integrated_lumi_pp.C");
}
