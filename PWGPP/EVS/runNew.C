#include "fstream"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
void runNew(TString fileList="runlist.txt"){
  TString runlist = gSystem->GetFromPipe("ls all");
  ifstream f(fileList.Data());
  Int_t run;
  gSystem->Exec("rm trending.root");
  while (!f.eof()){
    f >> run;
    if (runlist.Contains(Form("%i",run))) continue;
    gSystem->Exec(Form("aliroot -l -b -q 'runLevelEventStatQA.C\(\"%s\",%i,\"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB\"\)'","",run));
    gSystem->mkdir(Form("all/%i",run));
    gSystem->Exec(Form("mv trending.root all/%i/trending.root",run));
  }
  gSystem->Exec("hadd -f trending.root all/*/trending.root");
  gSystem->Exec("aliroot -l -b -q periodLevelQA.C");
}
