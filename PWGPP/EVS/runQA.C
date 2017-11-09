#include "fstream"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
void runQA(){
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2012/LHC12h/manual/");

//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16d/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16e/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16f/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16g/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16h/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16i/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16j/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16k/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15f/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15h/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15i/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15j/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15l/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15l/manual/");
//  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2015/LHC15o/manual/");
  runQA_period("/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16r/16r_old/");
}
void runQA_period(TString path="/afs/cern.ch/work/a/aliqaevs/www/data/2016/LHC16k/manual/"){
  gSystem->Exec(Form("ls %s > run.list",path.Data()));
  ifstream f;

  f.open("run.list");
  Char_t buffer[1000];
  while (!f.eof()){
    f >> buffer;
    Int_t run = TString(buffer).Atoi();
    if (run<99999) continue;
    printf("%i\n",run);
    TString statfile = "EventStat_temp.root";
    gSystem->Exec(Form("cp triggerInfo.C %s/000%i/",path.Data(),run));
    gSystem->Exec(Form("cp runLevelEventStatQA.C %s/000%i/",path.Data(),run));
    gSystem->Exec(Form("cd %s/000%i; aliroot -l -b -q 'runLevelEventStatQA.C\(\"%s\",%i\,\"local:///cvmfs/alice.cern.ch/calibration/data/2016/OCDB\")'",path.Data(),run,statfile.Data(),run));
  }
  f.close();

  gSystem->Exec(Form("cp triggerInfo.C %s/",path.Data()));
  gSystem->Exec(Form("cp periodLevelQA.C %s/",path.Data()));
  gSystem->Exec(Form("cp runLevelEventStatQA.C %s/",path.Data()));
  gSystem->Exec(Form("cd %s; hadd -f EventStat_temp.root */EventStat_temp.root",path.Data()));
  gSystem->Exec(Form("cd %s; aliroot -l -b -q runLevelEventStatQA.C",path.Data()));
  gSystem->Exec(Form("cd %s; hadd -f trending.root */trending.root",path.Data()));
  gSystem->Exec(Form("cd %s; aliroot -l -b -q periodLevelQA.C",path.Data()));
}
