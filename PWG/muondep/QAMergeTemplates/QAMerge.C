#include "TSystem.h"
#include "TList.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TRegexp.h"
#include "AliMUONVTrackerData.h"
#include "TFile.h"
#include "TKey.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TUrl.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TGridCollection.h"
#include "TROOT.h"
#include "TError.h"
#include "TH1.h"
#include "TMethodCall.h"

void CopyFromRemote(Int_t runNumber, const char* period, const char* pass);
void CopyFromRemote(const TList& files);
TList* GetFileList(Int_t runNumber, const char* period, const char* pass);
void QAMerge(Int_t runNumber, const char* period, const char* pass);
void QAMerge(const TList& files, const char* outputfile);

////______________________________________________________________________________
//AliMUONVTrackerData* GetTrackerData(TDirectory* dir, const char* trackerDataName)
//{
//  AliMUONVTrackerData* trackerData(0x0);
//
//  TIter next(dir->GetListOfKeys());
//  TKey* key;
//  
//  while ( ( key = static_cast<TKey*>(next()) ) && !trackerData )
//  {
//    if ( key->IsFolder() )
//    {
//      trackerData = GetTrackerData(static_cast<TDirectory*>(key->ReadObj()),trackerDataName);
//    }
//    else
//    {
//      TString name(key->GetName());
//      
//      if ( name.Contains(trackerDataName) )
//      {
//        trackerData = dynamic_cast<AliMUONVTrackerData*>(key->ReadObj());
//      }
//    }
//  }
//  
//  return trackerData;
//  
//}
//
////______________________________________________________________________________
//AliMUONVTrackerData* GetTrackerData(const char* filename, const char* trackerDataName)
//{
//  TFile* file = TFile::Open(filename);
//  
//  cout << filename << endl;
//  
//  if (!file) return 0x0;
//
//  AliMUONVTrackerData* trackerData = GetTrackerData(file,trackerDataName);  
//
//  delete file;
//  
//  return trackerData;
//  
//}

////______________________________________________________________________________
//AliMUONVTrackerData* MergeTrackerData(const TList& filelist, const char* trackerDataName)
//{
//  TIter next(&filelist);
//  TObjString* s;
//  
//  s = static_cast<TObjString*>(next());
//  
//  AliMUONVTrackerData* d = GetTrackerData(s->String().Data(),trackerDataName);
//  
//  if (!d)
//  {
//    cout << Form("ERROR : no tracker data found in %s",s->String().Data()) << endl;
//    return 0x0;
//  }
//  
//  AliMUONVTrackerData* merged = static_cast<AliMUONVTrackerData*>(d->Clone(Form("%s-MERGED",trackerDataName)));
//  
//  while ( ( s = static_cast<TObjString*>(next()) ) )
//  {
//    d = GetTrackerData(s->String().Data(),trackerDataName);
//    if (d)
//    {
//      TList list;
//      list.Add(d);
//      cout << "Merging " << s->String().Data() << " " << d->GetName() << " " << d->NumberOfEvents(-1) << endl;
//      
//      merged->Merge(&list);
//      delete d;
//    }
//  }
//  
//  merged->SetName(Form("%s-MERGED-%d",trackerDataName,runNumber));
//  
//  merged->SetTitle(Form("%s merged for run %d",trackerDataName,runNumber));
//  
//
//  return merged;
//}

//______________________________________________________________________________
void GenerateFileList(const char* runlist="philippe.runlist")
{
  if (!gGrid)
  {
    TGrid::Connect("alien://");    
  }
  if (!gGrid) return;
  
  ifstream in(runlist);
  Int_t run;
  
  Int_t n(0);
  
  while ( in >> run )
  {
    TString basedir("/alice/data/2009/LHC09d/");
    
    TString dir(Form("%s%09d/ESDs/pass1/*/Merged.QA.Data.root",basedir.Data(),run));
    
    TGridResult* r = gGrid->Ls(dir.Data());

    Int_t i(0);
    
    while ( r->GetFileName(i) ) 
    {
      cout << "alien://" << r->GetFileNamePath(i) << endl;
      ++i;      
    }
    
    delete r;
    
    ++n;
  }
}

//______________________________________________________________________________
void CopyFromRemote(Int_t runNumber, const char* period, const char* pass)
{
  TList* list = GetFileList(runNumber,period,pass);
  if ( list )
  {
    CopyFromRemote(*list);
  }
  delete list;
}

//______________________________________________________________________________
void CopyFromRemote(const TList& files)
{
  TStopwatch timer;
  
  TIter next(&files);
  TObjString* line;
  
  while ( ( line = static_cast<TObjString*>(next()) ) )
  {
    TUrl url(line->String());

    if ( TString(url.GetProtocol()) == "alien" )
    {
      if (!gGrid)
      {
        TGrid::Connect("alien://");
        if (!gGrid)
        {
          cout << "Cannot get alien connection" << endl;
          return;
        }
      }
    }

    TString file(url.GetFile());

    TString dir(gSystem->DirName(file));

    gSystem->mkdir(dir.Data(),kTRUE);

    if ( gSystem->AccessPathName(file.Data())==kFALSE)
    {
      cout << "Skipping copy of " << file.Data() << " as it already exists" << endl;
    }
    else
    {
      TFile::Cp(line->String().Data(),file.Data());
      if ( line->String().Contains("root_archive.zip") )
      {
        gSystem->Exec(Form("unzip %s -d %s",file.Data(),gSystem->DirName(file.Data())));
        gSystem->Exec(Form("rm %s",file.Data()));
      }
    }
  }
  
  timer.Print();
}

//______________________________________________________________________________
TList* GetFileList(Int_t runNumber, const char* period, const char* pass)
{
  if (!gGrid) TGrid::Connect("alien://");
  if (!gGrid)
  {
    return 0x0;
  }

  Int_t year;
  
  TString speriod(period);
  
  if (!speriod.BeginsWith("LHC"))
  {
    std::cerr << "Period not starting with LHC !" << std::endl;
    return 0x0;
  }
  
  year = 2000 + TString(speriod(3,3)).Atoi();
  
  TString sbasedir(Form("/alice/data/%d/%s/%09d/ESDs/%s",year,period,runNumber,pass));
  
  std::cout << sbasedir.Data() << std::endl;
  
  TString search("Merged.QA.Data.root");
  
  TGridResult* res = gGrid->Query(sbasedir.Data(),search.Data());
  
  Int_t nFiles = res->GetEntries();
  
  nFiles = res->GetEntries();
  
  Long64_t size(0);
  Int_t count(0);
  TList* files(0x0);
  
  for (Int_t i = 0; i < nFiles; ++i)
  {
    Long64_t thisSize = TString(res->GetKey(i,"size")).Atoll();
    
//    TUrl url(res->GetKey(i,"turl"));
    
    if (!files)
    {
      files = new TList;
      files->SetOwner(kTRUE);
    }
    
    files->Add(new TObjString(res->GetKey(i,"turl")));
    
    ++count;
    
    size += thisSize;
  }
  
  std::cout << Form("%d files for a total of %d MB",count,TMath::Nint(size/1024.0/1024)) << std::endl;

  delete res;
  
  return files;
}

//______________________________________________________________________________
void QAMerge(Int_t runNumber, const char* period, const char* pass)
{
  TList* list = GetFileList(runNumber,period,pass);
  if ( list )
  {
    QAMerge(*list,Form("QA.%s.%s.%d.root",period,pass,runNumber));
  }
  delete list;

}

//______________________________________________________________________________
TObject* GetObject(const char* filepath, const char* objectPath)
{
  TFile* f = TFile::Open(filepath);
  
  if (!f)
  {
    std::cerr << "Could not open file " << filepath << std::endl;
    return 0x0;    
  }
  
  TObject* o = f->Get(objectPath);
  
  if (!o)
  {
    std::cout << "Cannot get object " << objectPath << " from " << filepath << std::endl;
    return 0x0;
  }
  
  if ( o->InheritsFrom("TH1") )
  {
    TH1* h = static_cast<TH1*>(o);
    h->SetDirectory(0);
  }
  
  delete f;
  
  return o;
}

//______________________________________________________________________________
TObject* Merge(const TList& files, const char* objectPath)
{
  TIter next(&files);
  TObjString* filename;
  
  filename = static_cast<TObjString*>(next());
  
  TObject* d = GetObject(filename->String().Data(),objectPath);
  
  if (!d)
  {
    cout << Form("ERROR : no object found in %s",filename->String().Data()) << endl;
    return 0x0;
  }
  
  TString oname(gSystem->BaseName(objectPath));
  
  TObject* merged = d->Clone();
  
  TMethodCall callEnv;

  if ( merged->IsA() )
  {
    callEnv.InitWithPrototype(merged->IsA(), "Merge", "TCollection*");
  }
  
  if (!callEnv.IsValid())
  {
    std::cout << "callEnv invalid" << std::endl;
    delete merged;
    return 0x0;
  }

  while ( ( filename = static_cast<TObjString*>(next()) ) )
  {
//    std::cout << filename->String().Data() << " : " << objectPath << std::endl;
    
    d = GetObject(filename->String().Data(),objectPath);
    
    if (!d)
    {
      continue;
    }
    
    TList list;
    list.Add(d);
    callEnv.SetParam((Long_t) (&list));
    callEnv.Execute(merged);
    
    delete d;
  }
  
  return merged;

}

//______________________________________________________________________________
void QAMerge(const TList& files, const char* outputfile)
{
  const char* eventSpecie = "HighMultiplicity";
  
  TList objectsToMerge;
  objectsToMerge.SetOwner(kTRUE);

  objectsToMerge.Add(new TObjString(Form("MUON/Raws/%s/Expert/%s_hTrackerDDLNofEventsSeen",eventSpecie,eventSpecie)));

  objectsToMerge.Add(new TObjString(Form("MUON/RecPoints/%s/Expert/%s_RecPoints",eventSpecie,eventSpecie)));
  objectsToMerge.Add(new TObjString(Form("MUON/RecPoints/%s/Expert/%s_RecCharges",eventSpecie,eventSpecie)));

  objectsToMerge.Add(new TObjString(Form("MUON/Raws/%s/Expert/%s_RawCharges",eventSpecie,eventSpecie)));

  for ( Int_t i = 1; i <= 10; ++i )
  {
    objectsToMerge.Add(new TObjString(Form("MUON/RecPoints/%s/Expert/%s_hTrackerClusterHitMapForChamber%d",eventSpecie,eventSpecie,i)));
    objectsToMerge.Add(new TObjString(Form("MUON/ESDs/%s/Expert/%s_hESDClusterHitMap%d",eventSpecie,eventSpecie,i)));
  }

  TObjArray output;
  output.SetOwner(kTRUE);

  TIter next(&objectsToMerge);
  TObjString* s;
  
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    output.Add(Merge(files,s->String().Data()));
  }
  
  TFile* f = TFile::Open(outputfile,"recreate");
  
  next.Reset();

  Int_t i(0);
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TObject* o = output.At(i);
    if (o)
    {
      TString dir = gSystem->DirName(s->String());
      TString file = gSystem->BaseName(s->String());
      TDirectory* tdir = static_cast<TDirectory*>(f->Get(dir.Data()));
      if (!tdir)
      {
        f->mkdir(dir.Data()," ");
        tdir = static_cast<TDirectory*>(f->Get(dir.Data()));
      }
      tdir->cd();
      o->Write(file.Data());
    }
    ++i;
  }

  f->Write();
  
  f->Close();
  
  delete f;
}

//______________________________________________________________________________
void QAMerge(const char* inputfilelist, const char* mergedFile)
{
  if (!gGrid)
  {
    TGrid::Connect("alien://");    
  }
  if (!gGrid) return;
  
  TString sinput(inputfilelist);
  TList filelist;
  filelist.SetOwner(kTRUE);
  
  if (sinput.Contains(".xml"))
  {
    TGridCollection *coll = gGrid->OpenCollection(inputfilelist);
    if (!coll)
    {
      ::Error("MergeOutput", "Input XML collection empty.");
      return;
    }
    // Iterate grid collection
    while (coll->Next())
    {
      filelist.Add(new TObjString(coll->GetTURL()));
    }
  }
  else
  {
    ifstream in(inputfilelist);
    char filename[1024];

    while (in.getline(filename,1024,'\n'))
    {
      filelist.Add(new TObjString(filename));
    }
  }
    
  if ( filelist.GetEntries()==0 )
  {
    std::cout << "No file to be merged !" << std::endl;
  }
  else
  {
    QAMerge(filelist,mergedFile);
  }
  
  ofstream out;
  out.open("outputs_valid_merge", ios::out);
  out.close();
}

