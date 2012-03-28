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
void mergeInChunksTXT(const char* mlist, const char* dest, int maxFiles=700);


void merge(const char* outputDir, const char* pattern, Bool_t copyLocal=kFALSE)
{
  //
  // load libraries
  //
  printf("Merging with chunks copying turned %s\n",copyLocal ? "ON":"OFF");
  gROOT->Macro("LoadLibraries.C");
  //
  cpTimeOut(outputDir, pattern,10, copyLocal);
  //
  // local
  mergeInChunksTXT("calib.list","CalibObjects.root");
  //  AliFileMerger merger;
  //  merger.AddReject("esdFriend"); // do not merge
  //  merger.SetMaxFilesOpen(700);
  //  merger.IterTXT("calib.list","CalibObjects.root",kFALSE);

  // alien
  //merger.IterAlien(outputDir, "CalibObjects.root", pattern);

  return;
}

void cpTimeOut(const char * searchdir, const char* pattern, Int_t timeOut=10, Bool_t copyLocal)
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
  res->Print();
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
    Bool_t result = kTRUE;
    if (copyLocal) {
      dst.ReplaceAll("alien:///","");
      dst.ReplaceAll("/","_");
      TTimeStamp s1;
      result = TFile::Cp(src.Data(),dst.Data(),kTRUE);
      TTimeStamp s2;
      AliSysInfo::AddStamp(dst.Data(),counter, result);
    }
    if (result) {
      counter++;
      outputFile << dst.Data()<< endl;
    }
  }
  if (copyLocal) cout<<counter<<" files copied!"<<endl;
  else           cout<<counter<<" files registerd!"<<endl;
  outputFile.close();
  gSystem->Exec("mv syswatch.log syswatch_copy.log");
  return;
}

void mergeInChunksTXT(const char* mlist, const char* dest, int maxFiles)
{
  TH1::AddDirectory(0);
  AliFileMerger merger;
  //  merger.SetMaxFilesOpen(999);
  merger.AddReject("esdFriend"); // do not merge
  //
  if (maxFiles<2) maxFiles = 2;
  TString filesToMerge = mlist, fileDest = dest;
  if (filesToMerge.IsNull()) {printf("List to merge is not provided\n"); return;}
  if (fileDest.IsNull())     {printf("Merging destination is not provided\n"); return;}
  const char* tmpMerge[3]={"__merge_part0.root","__merge_part1.root","__part_of_calib.list"};
  //
  gSystem->ExpandPathName(filesToMerge);
  ofstream outfile;
  ifstream infile(filesToMerge.Data());
  if (!infile) {printf("No %s file\n",filesToMerge.Data()); return;}
  //
  int currTmp = 0, nfiles = 0, nparts = 0; // counter for number of merging operations
  string line;
  TString lineS;
  while ( !infile.eof() ) {
    getline(infile, line); 
    lineS = line;
    if (lineS.IsNull() || lineS.BeginsWith("#")) continue;
    int st = nfiles%maxFiles;
    if (st==0) { // new chunk should be started
      if (nfiles) { // merge prev. chunk
	outfile.close();
	merger.IterTXT(tmpMerge[2], tmpMerge[currTmp] ,kFALSE);
	printf("Merging to %s | %d files at step %d\n",tmpMerge[currTmp], nfiles,nparts);
      }
      outfile.open(tmpMerge[2], ios::out); // start list for new chunk  
      if (nparts++) {
	printf("Adding previous result %s | %d files %d at part\n",tmpMerge[currTmp], nfiles,nparts);
	outfile << tmpMerge[currTmp] << endl; // result of previous merge goes 1st
      }
      currTmp = (currTmp==0) ? 1:0;         // swap tmp files
    }
    outfile << line << endl;
    nfiles++;
  }
  // merge the rest
  merger.IterTXT(tmpMerge[2], dest ,kFALSE);
  outfile.close();
  infile.close();
  for (int i=0;i<3;i++) gSystem->Exec(Form("if [ -e %s ]; then \nrm %s\nfi",tmpMerge[i],tmpMerge[i]));
  printf("Merged %d files in %d steps\n",nfiles, nparts);
  //
}

