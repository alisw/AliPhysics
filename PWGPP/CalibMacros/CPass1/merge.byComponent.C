//___________________________________________________________________

void merge(TString component, const Char_t *inputstring)
{

  /* load libs */
  gROOT->Macro("LoadLibraries.C");
  TH1::AddDirectory(0);

  /* copy only */
  if (component == "COPY") {
    CopyCPass(inputstring, "AliESDfriends_v1.root");
    return;
  }
  
  /* merge component */
  MergeCPass(inputstring, component);
  
}

//___________________________________________________________________

void MergeCPass(const Char_t *list, TString component)
{
  AliFileMerger merger;
  /* select what to merge */
  if (component == "ALL")
    merger.AddReject("esdFriend");
  else
    merger.AddAccept(component.Data());
  /* merge */
  merger.IterTXT(list, "CalibObjects.root", kFALSE);
  /* notify */
  gSystem->Exec(Form("touch %s_merge_done", component.Data()));
  return;
}

//___________________________________________________________________

void CopyCPass(const char *searchdir, const char *pattern, Int_t timeOut=10)
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
  gSystem->Exec("touch copy_done");
  return;
}

