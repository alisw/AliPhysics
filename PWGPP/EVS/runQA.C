#include "fstream"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
void runQA(TString fileList="all.list"){
  TString outDir = fileList;
  outDir.ReplaceAll(".list","");
  gSystem->mkdir(outDir.Data());
  ifstream f;
  f.open(fileList);
  Char_t buffer[1000];
  gSystem->Exec("rm trending.root");
  while (!f.eof()){
    f >> buffer;
    Int_t run = TString(buffer).Atoi();
    printf("%i\n",run);
    TString statfile = "";
    gSystem->Exec(Form("aliroot -b -q 'runLevelEventStatQA.C\(\"%s\",%i\)'",statfile.Data(),run));
    gSystem->mkdir(Form("%s/%i",outDir.Data(),run));
    gSystem->Exec(Form("mv trending.root %s/%i/trending.root",outDir.Data(),run));
  }
  f.close();
  
  gSystem->Exec(Form("hadd -f trending.root %s/*/trending.root",outDir.Data()));
  gSystem->Exec("aliroot -b -q periodLevelQA.C");
}
