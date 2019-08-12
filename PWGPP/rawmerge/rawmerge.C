//
// Macro to create the "raw" data file with selected events
// source is 
//     $ALICE_PHYSICS/../src/PWGPP/rawmerge/rawmerge.C+
// 
// Original code: 
//    marian.ivanov@cern.ch
//    modifications:
//       mikolaj.krzewicki@cern.ch
//       mesut.arslandok@cern.ch
//

#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TGridCollection.h"
#include "TGrid.h"
#include "TStopwatch.h"
#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"

#include "AliSysInfo.h"

Bool_t makeAlienInputEventList( TString outputFileName="wn.list",
                                TString referenceFileName="filteredEvents.list",
                                TString inputCollectionName="wn.xml" );


void rawmerge( TString inputFileName="wn.xml",
               TString fullEventListFileName="event.list",
               TString outputFileNameTMP="filtered.root", Int_t timeOut=30)
{
   // Create the filtered raw data files using a file with list of raw chunks with event numbers.
   // inputFileName         - either the text file with chunkname+event number or xml collection
   // fullEventListFileName - if 1st arg is an xml collection, this is the full list of
   //                         chunks and events to be filtered acoording the xml contents
  //
   // if the file list is an xml collection (for running on alien),
   // first extract the available chunks and event numbers from the
   // reference file, and use that as input
  
//   TGrid::Connect("alien://") 
  gSystem->Setenv("XRDCLIENTMAXWAIT",Form("%d",timeOut));
  gEnv->SetValue("XNet.RequestTimeout", timeOut);
  gEnv->SetValue("XNet.ConnectTimeout", timeOut);
  gEnv->SetValue("XNet.TransactionTimeout", timeOut);
  gEnv->SetValue("XNet.FirstConnectMaxCnt", 2);
  TGrid * alien = TGrid::Connect("alien://",0,0,"t");
  TFile::SetOpenTimeout(timeOut);
  printf(" ------ TFile open timeout limit = %ds ------ \n",TFile::GetOpenTimeout()); 
  
  if (inputFileName.Contains(".xml")){
    TString tmp=inputFileName+".list";
    makeAlienInputEventList(tmp,fullEventListFileName,inputFileName);
    inputFileName=tmp;
  }

  Int_t eventNumber;
  ifstream files;
  files.open(inputFileName.Data());
  if (!files.is_open()){
    fprintf(stderr,"error: could not read event list file \"%s\". Exiting.\n",inputFileName.Data());
    return;
  }
   
  //TSeqCollection* listOfFiles = (TSeqCollection*)gROOT->GetListOfFiles();
  TList* listOfFiles = new TList();
  TString outputFileName;
  TString line;
  TString iURI;
  TString triggerType;
  TString  iURIold;
  TFile *ifile=0;
  TFile *ofile=0;
  TTree *itree=0;
  TTree *otree=0;
  Long64_t ievent=0;
  ULong64_t igid=0;
  Int_t ofilenumber=0;
  Int_t lineNumber=0;
  Long64_t eventold=0;
  ULong64_t igidOld=0;
  Double_t realTime=0;
  //
  AliRawVEvent*    fEvent=0;              // event
  TBranch *fEventBranch = 0;              // branch for event header   

  while (files.good()){
    TStopwatch timer; timer.Start();
    ++lineNumber;
      //read the line, do some checks
    line.ReadLine(files);
    TObjArray* entries = line.Tokenize(" ");
    TObjString* iURIobjstr = (TObjString*)entries->At(0);
    iURI.Clear();
    if (iURIobjstr) iURI=iURIobjstr->String();
    TObjString* ieventobjstr = (TObjString*)entries->At(1);
    if (ieventobjstr) ievent = ieventobjstr->String().Atoi();
    if (iURI.IsNull() || !ieventobjstr){
      printf("malformed line: %s, skipping...\n",line.Data());
      continue;
    }
    TObjString* triggerTypeobjstr = (TObjString*)entries->At(2);
    triggerType.Clear();
    if (triggerTypeobjstr) triggerType=triggerTypeobjstr->String();

    printf("> processing \"%s\" event %llu trigger \"%s\"...\n",iURI.Data(),ievent,triggerType.Data());
    if (ievent==eventold && iURI==iURIold){
      printf("duplicated continue\n");
      continue;
    }
    //
    if (!iURIold.Contains(iURI.Data())){
      //if a new file - open it
      printf("new file: %s\n",iURI.Data());
      delete ifile;
      ifile=0;
      AliSysInfo::AddStamp((iURI+"_OpenBegin").Data(),11,lineNumber); // open + file counter
      ifile=TFile::Open(iURI.Data());
      AliSysInfo::AddStamp((iURI+"_OpenEnd").Data(),10,lineNumber); // open + file counter
      if (!ifile){
        fprintf(stderr,"warning: could not open file for event \"%s\", skipping it...\n",iURI.Data());
        continue;
      }
    }
    else{
      //if same file, reuse it
      printf("using already open file: %s\n",iURI.Data());
    }
    iURIold=iURI;
    AliSysInfo::AddStamp((iURI+"_Begin").Data(),100,lineNumber); // dump file   + event   counter within file
    //
    TTree *itree=dynamic_cast<TTree*>(ifile->Get("RAW"));
    if (!itree){
      fprintf(stderr,"warning: could not find RAW tree for event \"%s\", skipping it...\n",iURI.Data());
      continue;
    }
    
    // manage output files
    if (!triggerType.IsNull()) triggerType.Prepend("_");
    outputFileName=outputFileNameTMP;
    outputFileName.ReplaceAll(".root","");
    outputFileName+=triggerType;
    outputFileName+=".root";
      
    ofile=dynamic_cast<TFile*>(listOfFiles->FindObject(outputFileName.Data()));
    if (!ofile) {
      printf("< creating output file \"%s\"\n",outputFileName.Data());
      ofile=TFile::Open(outputFileName,"recreate");
      if (!ofile){
        fprintf(stderr,"error: could not create output file: \"%s\" Exiting.\n",outputFileName.Data());
        return;
      }
      listOfFiles->Add(ofile);
    }
    ofile->cd();
    otree=dynamic_cast<TTree*>(ofile->Get("RAW"));
    if (!otree) {
      otree=itree->CloneTree(0);
    }

    fEventBranch = itree->GetBranch("rawevent");  // as in AliRawReaderRoot::AliRawReaderRoot  
    fEventBranch->SetAddress(&fEvent);           // access event header

    AliSysInfo::AddStamp((iURI+"_GetBegin").Data(),1001,lineNumber); // get entry + file counter
    Int_t size= itree->GetEntry(ievent);
    Int_t readEntry=itree->GetReadEntry();   
    otree->CopyAddresses(itree);

    Bool_t isOK=kTRUE;    
    AliSysInfo::AddStamp((iURI+"_GetEnd").Data(),1000,lineNumber); // get entry + file counter
    ULong64_t  timeStamp;
    {
      const UInt_t *id  = fEvent->GetHeader()->GetP("Id");                             // copy of AliRawReaderRoot::GetEventId()
      ULong64_t  period = id ? (((id)[0]>>4)&0x0fffffff): 0;                           // AliRawReader::Get<>
      ULong64_t  orbit  = id ? ((((id)[0]<<20)&0xf00000)|(((id)[1]>>12)&0xfffff)) : 0; // AliRawReader::Get<>
      ULong64_t  bcID   =  id ? ((id)[1]&0x00000fff) : 0;                              // AliRawReader::Get<>
      igid    = (((ULong64_t)period << 36) | ((ULong64_t)orbit << 12) |(ULong64_t)bcID); // AliRawReader::GetEventIdAsLong() 
      timeStamp=fEvent->GetHeader()->Get("Timestamp");  
      if (igid==igidOld){  // check debugger
	isOK=kFALSE;
	Int_t sizeNew=itree->GetEntry(ievent);
	Int_t readEntryNew=itree->GetReadEntry(); 
	fEvent->GetHeader()->Dump();      
      }      
    }
    printf("filling event (%llu)  gid (%llu) timestamp (%llu) in file %s\n",ievent,igid, timeStamp, ofile->GetName());
    AliSysInfo::AddStamp((iURI+"_FillBegin").Data(),2001,lineNumber); // get entry + file counter
    otree->Fill();
    AliSysInfo::AddStamp((iURI+"_FillEnd").Data(),2000,lineNumber); // get entry + file counter
    eventold=ievent;
    igidOld=igid;  
    //otree->CopyEntries(itree,Form("Entry$==%d",ievent),1);

      // reset input
    itree->ResetBranchAddresses();
    
    realTime+=timer.RealTime();
    printf(" ================================================= \n");
    printf(" ----- Merging time of event %d ----- and total time elapsed so far %f ----- \n", lineNumber, realTime);
    timer.Stop(); timer.Print();
    printf(" ================================================= \n");
    AliSysInfo::AddStamp((iURI+"End").Data(),3,lineNumber); // dump file   + event   counter within file
  }

   //close the files
  for (Int_t i=0; i<listOfFiles->GetEntries(); i++)
  {
    ofile=dynamic_cast<TFile*>(listOfFiles->At(i));
    if (!ofile) {continue;}
    otree=dynamic_cast<TTree*>(ofile->Get("RAW"));
    Long64_t nEntries=0;
    if (otree) { nEntries = otree->GetEntries(); }

     //write the file and close
    ofile->Write();
    printf("closing file: %s with %llu entries\n",ofile->GetName(),nEntries);
    delete ofile;

     //remove empty files
    if (nEntries==0)
    {
      gSystem->Unlink(ofile->GetName());
    }
  }
}

Bool_t makeAlienInputEventList( TString outputFileName,
                                TString referenceFileName,
                                TString inputCollectionName )
{
  TGridCollection *coll = gGrid->OpenCollection(inputCollectionName);
  if (!coll)
  {
    ::Error("makeAlienInputEventList", "Cannot open collection from %s", inputCollectionName.Data());
      return NULL;
  }
  TString configLine;
  TString chunkPath;
  TString chunkName;
  ifstream referenceFile;
  referenceFile.open(referenceFileName.Data());
  if (!referenceFile.is_open())
  {
    printf("could not open %s\n",referenceFileName.Data());
    return kFALSE;
  }
  gSystem->Unlink(outputFileName);
  ofstream outputFile;
  outputFile.open(outputFileName.Data());
  if (!outputFile.is_open())
  {
    printf("could not open %s\n",referenceFileName.Data());
    return kFALSE;
  }

  while (coll->Next())
  {
    chunkPath = coll->GetTURL();
    TObjArray* a = chunkPath.Tokenize("/");
    TObjString* objstr = dynamic_cast<TObjString*>(a->At(a->GetEntries()-1));
    if (!objstr)
    {
      printf("empty chunkPath from collection!\n");
      return kFALSE;
    }
    TString chunkName = objstr->String();
    chunkName.ReplaceAll("alien://","");
    while (referenceFile.good())
    {
      configLine.ReadLine(referenceFile);
      if (configLine.Contains(chunkName))
      {
        configLine=configLine.Strip(TString::kBoth);
        if (!configLine.BeginsWith("alien://"))
        {
          configLine.Prepend("alien://");
        }
        outputFile << configLine << endl;
        cout << configLine << endl;
      }
    }
      //jump to beginning and clear flags
    referenceFile.clear();
    referenceFile.seekg(0,ios::beg);
  }
  outputFile.close();
  referenceFile.close();
  return kTRUE;
}

