// $Id:$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <vector>
#include <map>
#include <Riostream.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TError.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TObjectTable.h>
#include "TreeClasses.C"
#endif

#define DUPLICATECHECK
#ifdef DUPLICATECHECK
#define DUPLICATEVERBOSE
#endif

struct EvInfo
{
  Int_t         fRun;
  Short_t       fArrId;
  UInt_t        fEntry; 
};


void mergeCorr(Int_t run_from=-1, 
               Int_t run_to=-1,
               const char *outPrefix = "merged",
               const char *inFileNames = "res/output_*.root");

void mergeCorr(Int_t nEvents,
               const char *outFileName,
               const char *inFileNames);

//-----------------------------------------------------------------------------------------------------

void mergeCorr(Int_t run_from,
               Int_t run_to,
               const char *outPrefix,
               const char *inFileNames)
{
  TChain *c = new TChain("MyTree");
  c->Add(inFileNames);

  std::map<Int_t,Int_t> goodruns;

  TObjArray *filenames = c->GetListOfFiles();
  Int_t nfiles = filenames->GetEntries();
  for(Int_t i=0; i<nfiles; ++i) {
    TString fname(filenames->At(i)->GetTitle());
    TFile *file = TFile::Open(fname);
    if (!file)
      continue;
    TList *output = dynamic_cast<TList*>(file->Get("output"));
    output->SetOwner(1);
    TTree *tree = dynamic_cast<TTree*>(output->FindObject("MyTree"));
    Int_t nent = tree->GetEntries();
    if (nent<1) {
      delete output;
      file->Close();
      delete file;
      continue;
    }

    MyHeader *header = 0;
    TBranch *br = tree->GetBranch("header");
    br->SetAddress(&header);

    std::vector<Int_t> multruns;

    Int_t runno = -1;
    for (Int_t ev=0;ev<nent;++ev) {
      br->GetEntry(ev);
      if (header->fRun==0) {
        continue;
      }
      if (runno<0) 
        runno = header->fRun;
      else if (runno!=header->fRun) {
        Bool_t found = 0;
        for (UInt_t j=0;j<multruns.size();++j) {
          if (multruns.at(j)==header->fRun) {
            found = 1;
            break;
          }
        }
        if (!found)
          multruns.push_back(header->fRun);
      }
    }
    if (multruns.size()<=0) {
      if ((runno>run_from) && (run_to<0 || runno<run_to)) {
        cout << "Found run " << runno << " associated to " << fname << endl;
        goodruns.insert(pair<Int_t,Int_t>(runno,0));
        TString dir(gSystem->DirName(inFileNames));
        dir+="/runs";
        TString link(Form("%s/%s",gSystem->pwd(),fname.Data()));
        TString rundir(Form("%s/%d", dir.Data(), runno));
        TString cmd(Form("mkdir -p %s && cd %s && ln -sf %s", 
                         rundir.Data(), rundir.Data(), link.Data()));
        gSystem->Exec(cmd);
      }
    } else {
      cout << "Multiple runs (" << runno;
      for (UInt_t j=0;j<multruns.size();++j)
        cout << ", " << multruns.at(j);
      cout << ") found in " << fname << endl;
      cout << " Cannot deal with this case yet!!!" << endl;
    }
    delete output;
    file->Close();
    delete file;
  }

  for (map<Int_t, Int_t>::iterator it = goodruns.begin();it!=goodruns.end();++it) {
    TString infiles(Form("%s/runs/%d/*.root",gSystem->DirName(inFileNames),it->first));
    TString outdir(Form("%s/mergedruns",gSystem->DirName(inFileNames)));
    TString outFileName(Form("%s/%s_run%d.root",outdir.Data(),outPrefix,it->first));
    gSystem->Exec(Form("mkdir -p %s",outdir.Data()));
    gSystem->Exec(Form("root -b -q mergeCorr.C+'(-1,\"%s\",\"%s\")'",
                       outFileName.Data(),infiles.Data()));
  }
}

//-----------------------------------------------------------------------------------------------------

void mergeCorr(Int_t nEvents,
               const char *outFileName,
               const char *inFileNames)
{
  TChain *c = new TChain("MyTree");
  c->Add(inFileNames);

  map<ULong64_t, EvInfo*> einfos;
#ifdef DUPLICATEVERBOSE
  map<ULong64_t, MyHeader*> eheaders;
#endif

  Int_t totnent = 0;
  Int_t totnsus = 0;
#ifdef DUPLICATECHECK
  Int_t totncan = 0;
#ifdef DUPLICATEVERBOSE
  Int_t totndup = 0;
#endif
#endif

  TObjArray *maps = new TObjArray;
  maps->SetOwner(1);
  Int_t lrun = 0;
  TObjArray *filenames = c->GetListOfFiles();
  Int_t nfiles = filenames->GetEntries();
  for(Int_t i=0; i<nfiles; ++i) {
    TString fname(filenames->At(i)->GetTitle());
    TFile *file = TFile::Open(fname);
    if (!file)
      continue;
    TList *output = dynamic_cast<TList*>(file->Get("output"));
    output->SetOwner(1);
    TTree *tree = dynamic_cast<TTree*>(output->FindObject("MyTree"));
    Int_t nent = tree->GetEntries();
    if (nent<1) {
      delete output;
      file->Close();
      delete file;
      continue;
    }

    totnent += nent;
    cout << "Found "<< nent << " entries in " << fname << endl;
    TObjArray *filesnames = dynamic_cast<TObjArray*>(output->FindObject("filesnames"));

    MyHeader *header = 0;
    TBranch *br = tree->GetBranch("header");
    br->SetAddress(&header);

    for (Int_t ev=0;ev<nent;++ev) {
      br->GetEntry(ev);
      if (header->fRun==0) {
#ifdef DUPLICATEVERBOSE
        cout << "--> Suspicious found for " << header->fRun << " " 
             << header->fTime << " " << header->fPeriod << " " << header->fBx << endl;
#endif
        ++totnsus;
        continue;
      }

      //ULong64_t evtid = header->fTime;
      ULong64_t evtid = header->GetEventId();
#ifdef DUPLICATECHECK
      if (einfos.find(evtid)!=einfos.end()) {
#ifdef DUPLICATEVERBOSE        
        cout << "--> Potential duplicate found for " << header->fRun << " " 
             << header->fTime << " " << header->fPeriod << " " << header->fBx << endl;
#endif
      ++totncan;
#ifdef DUPLICATEVERBOSE        
        map<ULong64_t, MyHeader*>::iterator it = eheaders.find(evtid);
        if ( (it->second->fRun==header->fRun) && 
             (it->second->fNTracks==header->fNTracks) &&
             (it->second->fNSelTracks==header->fNSelTracks) &&
             (it->second->fNTracklets==header->fNTracklets) &&
             (it->second->fVz==header->fVz)) {
          cout << "--> Duplicate found:" << endl;
          cout << " Run: " << it->second->fRun << " vs " << header->fRun << endl; 
          cout << " Ntracks: " << it->second->fNTracks << " vs " << header->fNTracks << endl; 
          cout << " Nseltracks: " << it->second->fNSelTracks << " vs " << header->fNSelTracks << endl; 
          cout << " Ntracklets: " << it->second->fNTracklets << " vs " << header->fNTracklets << endl; 
          cout << " Vz: " << it->second->fVz << " vs " << header->fVz << endl;
          cout << " Duplicate event " << header->fEvNumberInFile << " in file " 
               << filesnames->At(header->fFileId)->GetName() << endl;
          ++totndup;
        }
#endif
        continue;
      }
#ifdef DUPLICATEVERBOSE        
      MyHeader *cheader = new MyHeader(*header);
      eheaders.insert( pair<ULong64_t, MyHeader*>(evtid,cheader) );
#endif
#endif
      EvInfo *cinfo = new EvInfo;
      cinfo->fRun = header->fRun;
      cinfo->fArrId = i;
      cinfo->fEntry = ev;
      einfos.insert( pair<ULong64_t, EvInfo*>(evtid,cinfo) );
      if (lrun != header->fRun) {
        for (Int_t lx=0;lx<4;++lx) {
          lrun = header->fRun;
          TString cname(Form("L%d_run%d_map",lx,lrun));
          if (lx==3) 
            cname=Form("TrigClass_run%d_map",lrun);
          TMap *map = dynamic_cast<TMap*>(output->FindObject(cname));
          if (map) {
            maps->Add(map);
            output->Remove(map);
          }
        }
      }
    }
    delete output;
    file->Close();
    delete file;
    if ((nEvents>0) && (einfos.size()>=(UInt_t)nEvents))
      break;
  }
  cout << "Found total: " << totnent << endl;
  cout << "Found suspicious: " << totnsus << endl;
#ifdef DUPLICATECHECK
  cout << "Found potential duplicates: " << totncan << endl;
#ifdef DUPLICATEVERBOSE        
  cout << "Found duplicates: " << totndup << endl;
#endif
#endif

  // sort output tree
  std::vector<TFile*> files(nfiles);
  std::vector<TTree*> trees(nfiles);
  std::vector<MyHeader*> headers(nfiles);
  std::vector<TClonesArray*> partas(nfiles);
  std::vector<TClonesArray*> trklas(nfiles);

  MyHeader *newHeader = new MyHeader;
  if (TClass::GetClass("MyPart"))
    TClass::GetClass("MyPart")->IgnoreTObjectStreamer();
  TClonesArray *newParts     = new TClonesArray("MyPart",10000);
  if (TClass::GetClass("MyTracklet"))
    TClass::GetClass("MyTracklet")->IgnoreTObjectStreamer();
  TClonesArray *newTracklets = new TClonesArray("MyTracklet",10000);

  TFile *newFile = TFile::Open(outFileName,"recreate","Merged/Sorted MyTree created by mergeCorr.C",5);
  TDirectory::TContext cd(newFile);
  for (Int_t m=0;m<maps->GetEntries();++m) {
    TMap *map = (TMap*)maps->At(m);
    //map->Print();
    map->Write(map->GetName(),TObject::kSingleKey);
  }
  TTree *newTree = new TTree("MyTree", "MyTree created by MyAliTask.cxx");
  newTree->Branch("header",&newHeader, 32*1024, 99);
  newTree->Branch("parts",&newParts, 256*1024, 99);
  newTree->Branch("tracklets",&newTracklets, 256*1024, 99);
  Int_t newEntries = 0;
  for (map<ULong64_t, EvInfo*>::iterator it = einfos.begin();it!=einfos.end();++it) {
    Short_t id = it->second->fArrId;
    UInt_t entry = it->second->fEntry;
    if (!trees.at(id)) { // connect new tree
      TString fname(filenames->At(id)->GetTitle());
      TFile *file = TFile::Open(fname);
      TDirectory::TContext cd2(newFile,file);
      files[id] = file;
      TList *output = dynamic_cast<TList*>(file->Get("output"));
      TTree *tree = dynamic_cast<TTree*>(output->FindObject("MyTree"));
      trees[id] = tree;
      output->Remove(tree);
      output->SetOwner(1);
      delete output;
      tree->SetBranchAddress("header",&headers[id]);
      tree->SetBranchAddress("parts",&partas[id]);
      tree->SetBranchAddress("tracklets",&trklas[id]);
    }
    files[id]->cd();
    trees.at(id)->GetEntry(entry);
    *newHeader = *(headers[id]);
    newParts->Clear();
    Int_t nparts = partas[id]->GetEntries();
    for (Int_t p = 0; p < nparts; ++p) {
      MyPart *orig = dynamic_cast<MyPart*>(partas[id]->At(p));
      new((*newParts)[p]) MyPart(*orig);
    }
    Int_t ntrkls = trklas[id]->GetEntries();
    for (Int_t p = 0; p < ntrkls; ++p) {
      MyTracklet *orig = dynamic_cast<MyTracklet*>(trklas[id]->At(p));
      new((*newTracklets)[p]) MyTracklet(*orig);
    }
    newFile->cd();
    newTree->Fill();
    ++newEntries;
    delete it->second;
    if ((nEvents>0) && (newEntries>=nEvents))
      break;
  }

  newTree->Write();
  newFile->Close();
  delete newFile;

  for(Int_t i=0;i<nfiles;++i) {
    delete trees.at(i);
    delete files.at(i);
  }
}
