/*
  merge output calib objects on Alien
  using AliFileMerger functionality

  Directory with runCalibTrain output: outputDir
  pattern: AliESDfriends_v1.root
  Output file name: CalibObjects.root

  Example:
  .L $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/merge.C
  merge("alien:///alice/cern.ch/user/j/jotwinow/CalibTrain/output","AliESDfriends_v1.root");
*/

void merge(const char* outputDir, const char* pattern)
{
  //
  // load libraries
  //
  gROOT->Macro("LoadLibraries.C");
  //
  TH1::AddDirectory(0);

  //
  AliFileMerger merger;
  merger.AddReject("esdFriend"); // do not merge
  cpTimeOut(outputDir, pattern,10);

  // local
   merger.IterTXT("calib.list","CalibObjects.root",kFALSE);

  // alien
  //merger.IterAlien(outputDir, "CalibObjects.root", pattern);

return;
}

void cpTimeOut(const char * searchdir, const char* pattern, Int_t timeOut=10)
{

  gSystem->Setenv("XRDCLIENTMAXWAIT",Form("%d",timeOut));
  gEnv->SetValue("XNet.RequestTimeout", timeOut);
  gEnv->SetValue("XNet.ConnectTimeout", timeOut);
  gEnv->SetValue("XNet.TransactionTimeout", timeOut);
  TFile::SetOpenTimeout(timeOut);

  TGrid::Connect("alien");

  TString filelist;
  TString command;
  command = Form("find %s/ %s", searchdir, pattern);
  cerr<<"command: "<<command<<endl;
  TGridResult *res = gGrid->Command(command);
  if (!res) return;
  TIter nextmap(res);
  TMap *map = 0;

  ofstream outputFile;
  outputFile.open(Form("calib.list"));
  Int_t counter=0;

  while((map=(TMap*)nextmap()))
  {
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length())
    {
      delete res;
      break;
    }

    TString src=Form("%s",objs->GetString().Data());
    TString dst=src;
    dst.ReplaceAll("alien:///","");
    dst.ReplaceAll("/","_");
    TTimeStamp s1;
    Bool_t result = TFile::Cp(src.Data(),dst.Data(),kTRUE);
    TTimeStamp s2;
    AliSysInfo::AddStamp(dst.Data(),counter, result);
    if (result)
    {
      counter++;
      outputFile << dst.Data()<< endl;
    }
  }
  cout<<counter<<" files copied!"<<endl;

  outputFile.close();
  gSystem->Exec("mv syswatch.log syswatch_copy.log");
  return;
}

