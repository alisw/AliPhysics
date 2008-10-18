#include "TGrid.h"
#include "TGridCollection.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TMap.h"
#include "TGridResult.h"
#include "TAlien.h"
#include "TAlienResult.h"
#include <fstream>

/*
  Simple toolkit for Alien
  Use ROOT C++ interfaces classes to access information from ALIEN

  Important functionality:
  MakeJobList - prepared the list which is later used by agent.sh
  to analyze data on batch system 



  See example bellow and class heeader for details.
  !!! The functionality not expelicitly mentioned in example can be not working
  Code is still under development

*/


 
/*
  gSystem->Load("libXrdClient.so");
  gSystem->Load("libNetx.so");
  //Raw data example
  char *mask = "20225";
  char *path = "/alice/data/2008/LHC08d"

  .L $ALICE_ROOT/TPC/macros/testTPC/AlienToolkit.cxx+
  //
  AlienToolkit toolkit;
  toolkit.MakeCollection(path,mask); // make a list of the registerd data
  //toolkit.StageCastor();    // stage files on castor
  //
  toolkit.MakeJobList("job.list","", "", "FindKrClustersRaw");

*/



class AlienToolkit: public TObject{
public:
  AlienToolkit();
  TGridCollection* MakeCollection(const char *path,  char *mask);
  void             Stage();
  void             StageCastor();
  void             LocalCopy(const char* destination);
  void             RemoteCopy(const char* destination="root://gsiaf.gsi.de:1094/", Int_t maxfiles=20);
  void             RemoteCopyAlien(const char* destination="root://gsiaf.gsi.de:1094/", Int_t maxfiles=20);

  void             PrintPFN();
  void             PrintLFN();
  void             MakeJobList(const  char * outname, const char *outputPrefix,  const char *action, const char *suffix);
  static Bool_t    IsDir(const char * name);
  static Bool_t    IsFile(const char * name);
  static Bool_t    ResubmitJobs();

public:
  TGridCollection *fCollection;
  TObjArray        fInfoArray;
  ClassDef(AlienToolkit,1)
};

ClassImp(AlienToolkit)

AlienToolkit::AlienToolkit():
  fCollection(0)
{
  if (!gGrid)
    TGrid::Connect("alien://",0,0,"t");
}




Bool_t AlienToolkit::IsDir(const char * name){
  //
  // Check if it is directory
  // 
  void *dir = gSystem->OpenDirectory(name);
  return (dir!=0);
}


Bool_t AlienToolkit::IsFile(const  char *name ){
  //
  // Check if is it file
  //
  Long_t id, size, flags, modtime ;
  Int_t ret = gSystem->GetPathInfo(name,&id, &size, &flags, &modtime);
  if (ret==0) return kFALSE;
  if ((flags&2)>0) return kFALSE;
  return kTRUE;
}


TGridCollection* AlienToolkit::MakeCollection(const char *path,  char *mask){
  //
  /*char *mask = "root_archive.zip";
    char *path = "/alice/cern.ch/user/m/miranov/test2007/"
  */
  //
  TGridResult* query = gGrid->Query(path,mask);  
  TGridCollection* collection=gGrid->OpenCollectionQuery(query);
  //  collection->SelectFile(mask,1,20); // select files 0-100
  collection->LookupSUrls();
  //collection->CheckIfOnline();
  collection->Reset();
  fCollection = collection;
  //
  collection->Reset();
  Int_t counter=0;
  Int_t counterm=0;
  Int_t counterf=0;
  TMap *filemap;
  TObjString  lastLFN;

  while ( (filemap = collection->Next())) { 
    TIterator *nextfile = filemap->MakeIterator();
    TMap *attributes;
    while ((attributes = (TMap *) nextfile->Next())) {
      printf("%d\t%d\t%d\n",counter, counterm,counterf); 
      TMap * map = new TMap;

      TObjString *surl = new TObjString(collection->GetSURL(attributes->GetName()));
      TObjString *turl =  new TObjString(collection->GetTURL(attributes->GetName()));
      TObjString *lfn =  new TObjString(collection->GetLFN(attributes->GetName()));
      map->Add(new TObjString("alienLFN"),lfn);
      map->Add(new TObjString("alienSURL"),surl);
      map->Add(new TObjString("alienTURL"),turl);
      
      if (lastLFN.String().CompareTo(lfn->String())==0) continue;
      lastLFN = *lfn;
      Int_t   isOnline = collection->IsOnline(attributes->GetName());
      printf("Base Name:\t%s\n",attributes->GetName());
      printf("Size:\t%d\n",(Int_t)collection->GetSize());
      printf("IsOnline:\t%d\n",collection->IsOnline());
      printf("IsSelected:\t%d\n",collection->IsSelected());
      printf("IsOnline:\t%d\n",isOnline);
      printf("LFN  Name:\t%s\n",lfn->String().Data());
      printf("TURL Name:\t%s\n",turl->String().Data());
      printf("SURL Name:\t%s\n",surl->String().Data());
      counter++;
      counterf++;
      fInfoArray.AddLast(map);
      
    }
    counterf=0;
    counterm++;
  }
}


void AlienToolkit::Stage(){
  //
  // Stage selected alien files
  //
  Int_t entries = fInfoArray.GetEntries();
  ofstream aout("stage.txt");
  for (Int_t i=0; i<entries;i++){
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *lfn = (TObjString*)map("alienLFN");
    if (!lfn) continue;
    printf("Staging submitfor\t%s\n",lfn->String().Data());
    aout<<"stage "<<lfn->String().Data()<<endl;
  }
  aout.close();
  gSystem->Exec("aliensh file:stage.txt");
}

void AlienToolkit::StageCastor(){
  //
  // Stage selected alien files
  //
  Int_t entries = fInfoArray.GetEntries();
  ofstream aout("stage.sh");
  aout<<"#!/usr/local/bin/bash"<<endl;
  for (Int_t i=0; i<entries;i++){
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *pfn = (TObjString*)map("alienSURL");
    if (!pfn) continue;
    if (!pfn->String().Contains("//castor")) continue;
    Char_t * cstr = &( pfn->String()[pfn->String().Index("/castor")]);
    if (!cstr) continue;
    printf("Staging submitfor\t%s\n",cstr);
    char command[1000];
    sprintf(command,"stager_qry -M %s | grep /castor | gawk  \'{print $3;}\'", cstr);
    gSystem->Exec(command);    
    sprintf(command,"stager_get -M %s", cstr);
    aout<<command<<endl;
    //gSystem->Exec(command);
  }
  aout.close();
  //gSystem->Exec("source stage.sh &");
}

void AlienToolkit::LocalCopy(const char *destination){
  //
  // Copy selected files to the destination directory
  // the LFN path name translated to the directory name replacing
  // separtor - the flat structure is created 

  Int_t entries = fInfoArray.GetEntries();
  ofstream aout("stage.txt");
  for (Int_t i=0; i<entries;i++){
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *lfn = (TObjString*)map("alienLFN");
    if (!lfn) continue;
    printf("Staging submitfor\t%s\n",lfn->String().Data());
    TString dnames=lfn->String().Data();
    dnames.ReplaceAll(".root_dir","");
    dnames.ReplaceAll("/","_");
    TString dname=destination;
    dname+=dnames;
    aout<<"cp  "<<lfn->String().Data()<<" "<<dname.Data()<<endl;
  }
  aout.close();
  gSystem->Exec("aliensh file:stage.txt");
}


void              AlienToolkit::PrintPFN(){
  //
  //
  // 
  Int_t entries = fInfoArray.GetEntries();
  for (Int_t i=0; i<entries;i++){
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *lfn = (TObjString*)map("alienLFN");
    TObjString *pfn = (TObjString*)map("alienSURL");
    if (!lfn) continue;
    if (!pfn) continue;
    printf("%s\n",pfn->String().Data());
  }
}

void              AlienToolkit::PrintLFN(){
  //
  //
  // 
  Int_t entries = fInfoArray.GetEntries();
  for (Int_t i=0; i<entries;i++){
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *lfn = (TObjString*)map("alienLFN");
    TObjString *pfn = (TObjString*)map("alienSURL");
    if (!lfn) continue;
    if (!pfn) continue;
    printf("%s\n",lfn->String().Data());
  }
}


void   AlienToolkit::MakeJobList(const  char * outname, const char *outputPrefix, const char *action, const char *suffix){
  //
  //
  //
   Int_t entries = fInfoArray.GetEntries();
  ofstream aout(outname);
  for (Int_t i=0; i<entries;i++){
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *lfn = (TObjString*)map("alienLFN");
    TObjString *pfn = (TObjString*)map("alienSURL");
    if (!lfn) continue;
    if (!pfn) continue;
     if (lfn->String().Contains(".tag.")) continue;
    printf("Job info\t%s\n",lfn->String().Data());
    TString jobID=lfn->String().Data();
    jobID.ReplaceAll("/","_");
    //
    //
    //
    //
    TString outputDir=outputPrefix;
    outputDir+=lfn->String().Data();
    
    aout<<jobID<<" "<<pfn->String().Data()<<"  "<<outputDir.Data()<<" "<<action;
    if (suffix)  aout<<" "<<suffix<<"\n";
    aout<<endl;
  }
  aout.close();
  //
}


void AlienToolkit::RemoteCopy(const char *destination,Int_t maxfiles){
  //
  // Copy selected files to the destination directory
  // the LFN path name translated to the directory name replacing
  // separtor - the flat structure is created 

  Int_t entries = fInfoArray.GetEntries();
  ofstream *aout=0;
  for (Int_t i=0; i<entries;i++){
    if (i%maxfiles==0){
      if (aout) aout->close();
      aout = new ofstream(Form("stage_%d.sh",i/maxfiles));
      (*aout)<<"!/bin/bash\n";
      (*aout)<<"source ~/.balice\n";
    }
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *pfn = (TObjString*)map("alienSURL");
    TObjString *lfn = (TObjString*)map("alienLFN");
    if (!pfn) continue;
    if (!lfn) continue;
    TString dnames=lfn->String().Data();
    TString dname=destination;
    dname+=dnames;
    
    (*aout)<< "mkdirhier "<<gSystem->DirName(dname.Data())<<endl;
    (*aout)<<"xrdcp -d 1 "<<pfn->String().Data()<<" "<<dname.Data()<<endl;
    if (dnames.Contains(".zip")){
      (*aout)<<"unzip "<<pfn->String().Data()<<endl;
    }
  }
  aout->close();
  //gSystem->Exec("aliensh file:stage.txt");
}

void AlienToolkit::RemoteCopyAlien(const char *destination,Int_t maxfiles){
  //
  // Copy selected files to the destination directory
  // the LFN path name translated to the directory name replacing
  // separtor - the flat structure is created 
 
  Int_t entries = fInfoArray.GetEntries();
  ofstream *aout=0;
  for (Int_t i=0; i<entries;i++){
    if (i%maxfiles==0){
      if (aout) aout->close();
      aout = new ofstream(Form("stage_%d.sh",i/maxfiles));
      (*aout)<<"!/bin/bash\n";
      (*aout)<<"source ~/.balice\n";
      (*aout)<<"source ~/.aliensetup\n";      
    }
    TMap &map = *((TMap*)fInfoArray.At(i));
    TObjString *pfn = (TObjString*)map("alienSURL");
    TObjString *lfn = (TObjString*)map("alienLFN");
    if (!pfn) continue;
    if (!lfn) continue;
    TString dnames=lfn->String().Data();
    TString dname=destination;
    dname+=dnames;
    (*aout)<< " echo Copy"<<gSystem->DirName(dname.Data())<<endl;
    (*aout)<< "mkdirhier "<<gSystem->DirName(dname.Data())<<endl;
    (*aout)<<"alien_cp  alien://"<<lfn->String().Data()<<" "<<dname.Data()<<endl;
    if (dnames.Contains(".zip")){
      (*aout)<< "cd  "<<gSystem->DirName(dname.Data())<<endl;
      (*aout)<<"unzip "<<gSystem->BaseName(dname.Data())<<endl;
    }
  }
  aout->close();
  //gSystem->Exec("aliensh file:stage.txt");
}

Bool_t AlienToolkit::ResubmitJobs(){
  //
  // Resubmit the processes finished in error state
  //
  //gSystem->Exec("alien_ps | grep EE > ps.txt")
  gSystem->Exec("alien_ps | grep EE |gawk '{print $2}' > psEE.txt");
  gSystem->Exec("alien_ps | grep EXPIRED |gawk '{print $2}' > psEXPIRED.txt");
  //
  ifstream in;
  in.open("psEXPIRED.txt");
  TString line;
  TString result;
  while(in.good()) {
    in >> line;
    //    line.ReplaceAll("-","");
    printf("Resubmiting %s\n",line.Data());
    Int_t pos = line.First("-");
    printf("%s\n",Form("alien_resubmit %s",&line.Data()[pos+1]));
    gSystem->Exec(Form("alien_resubmit %s",&line.Data()[pos+1]));
  }
  in.close();
}





/*
  Problems:
1. ????  fCollection.IsOnline("")

  #12 0x04e9a03f in XrdClientAdmin::IsFileOnline (this=0x94a91b0, vs=@0xbfe07f60, vb=@0xbfe07f80) at ../../src/XrdClient/XrdClientVector.hh:55
#13 0x0423a975 in TXNetSystem::IsOnline (this=0x9652800,
    path=0x9661488 "root://voalice04.cern.ch:1094//castor/cern.ch/alice/2005_castor2/15/23612/9e787bec-b010-11dc-a80c-000e0c3e6d91")
    at netx/src/TXNetSystem.cxx:476
#14 0x04237c6e in TXNetFileStager::IsStaged (this=0x96526b8,
    path=0x9661400 "root://voalice04.cern.ch:1094//castor/cern.ch/alice/2005_castor2/15/23612/9e787bec-b010-11dc-a80c-000e0c3e6d91")
    at netx/src/TXNetFileStager.cxx:63
#15 0x00551117 in TFileStager::GetStaged (this=0x96526b8, pathlist=0x9649820) at net/src/TFileStager.cxx:58
#16 0x010e5d96 in TAlienCollection::CheckIfOnline (this=0x95da8d8, bulk=true) at alien/src/TAlienCollection.cxx:1143


2. id xrootd on voalice04.cern.ch ?
   "root://voalice04.cern.ch:1094//castor/cern.ch/alice/2005_castor2/15/23612/9e787bec-b010-11dc-a80c-000e0c3e6d91"

3. How to check if the file exist in redirector ?
   
4. Ho

*/

