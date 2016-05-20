#include "fstream"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
void runQA(TString path="/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16h/manual/"){
  gSystem->Exec(Form("ls %s > run.list",path.Data()));
  ifstream f;
  f.open("run.list");
  Char_t buffer[1000];
  while (!f.eof()){
    f >> buffer;
    Int_t run = TString(buffer).Atoi();
    if (run<99999) continue;
    printf("%i\n",run);
    TString statfile = "event_stat.root";
    gSystem->Exec(Form("cp triggerInfo.C %s/000%i/",path.Data(),run));
    gSystem->Exec(Form("cp runLevelEventStatQA.C %s/000%i/",path.Data(),run));
//    gSystem->Exec(Form("cd %s/000%i",path.Data(),run));
    gSystem->Exec(Form("cd %s/000%i; aliroot -b -q 'runLevelEventStatQA.C\(\"%s\",%i\,\"raw://\")'",path.Data(),run,statfile.Data(),run));
  }
  f.close();
  gSystem->Exec(Form("cd %s; hadd -f trending.root */trending.root",path.Data()));
  gSystem->Exec(Form("cp periodLevelQA.C %s/",path.Data()));
  gSystem->Exec(Form("cd %s; aliroot -b -q periodLevelQA.C",path.Data()));
}
