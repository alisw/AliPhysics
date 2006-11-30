// @(#) $Id$
// Original: AliHLTFileHandler.cxx,v 1.49 2005/06/23 17:46:55 hristov 

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>, Anders Vestbo <mailto:vestbo$fi.uib.no>, C. Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TClonesArray.h>
#include <TSystem.h>
#include <TMath.h>

#ifdef use_newio
#include <AliRunLoader.h>
#endif
#include <AliTPCParamSR.h>
#include <AliTPCDigitsArray.h>
#include <AliTPCClustersArray.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersRow.h>
#include <AliSimDigits.h>

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCFileHandler.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/** \class AliHLTTPCFileHandler
<pre>
//_____________________________________________________________
// AliHLTTPCFileHandler
//
// The HLT ROOT <-> binary files handling class
//
// This class provides the interface between AliROOT files,
// and HLT binary files. It should be used for converting 
// TPC data stored in AliROOT format (outputfile from a simulation),
// into the data format currently used by in the HLT framework. 
// This enables the possibility to always use the same data format, 
// whether you are using a binary file as an input, or a AliROOT file.
//
// For example on how to create binary files from a AliROOT simulation,
// see example macro exa/Binary.C.
//
// For reading a AliROOT file into HLT format in memory, do the following:
//
// AliHLTTPCFileHandler file;
// file.Init(slice,patch);
// file.SetAliInput("galice.root");
// AliHLTTPCDigitRowData *dataPt = (AliHLTTPCDigitRowData*)file.AliDigits2Memory(nrows,eventnr);
// 
// All the data are then stored in memory and accessible via the pointer dataPt.
// Accesing the data is then identical to the example 1) showed in AliHLTTPCMemHandler class.
//
// For converting the data back, and writing it to a new AliROOT file do:
//
// AliHLTTPCFileHandler file;
// file.Init(slice,patch);
// file.SetAliInput("galice.root");
// file.Init(slice,patch,NumberOfRowsInPatch);
// file.AliDigits2RootFile(dataPt,"new_galice.root");
// file.CloseAliInput();
</pre>
*/

ClassImp(AliHLTTPCFileHandler)

AliHLTTPCFileHandler::AliHLTTPCFileHandler(Bool_t b)
  :
  fInAli(NULL),
#ifdef use_newio
  fUseRunLoader(kFALSE),
#endif
  fParam(NULL),
  fMC(NULL),
  fDigits(NULL),
  fDigitsTree(NULL),
  fIndexCreated(kFALSE),
  fUseStaticIndex(b)
{
  //Default constructor

  for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++)
    for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++) 
      fIndex[i][j]=-1;

  if(fUseStaticIndex&&!fgStaticIndexCreated) CleanStaticIndex();
}

AliHLTTPCFileHandler::AliHLTTPCFileHandler(const AliHLTTPCFileHandler& ref)
  :
  fInAli(NULL),
#ifdef use_newio
  fUseRunLoader(kFALSE),
#endif
  fParam(NULL),
  fMC(NULL),
  fDigits(NULL),
  fDigitsTree(NULL),
  fIndexCreated(kFALSE),
  fUseStaticIndex(ref.fUseStaticIndex)
{
  HLTFatal("copy constructor untested");
}

AliHLTTPCFileHandler& AliHLTTPCFileHandler::operator=(const AliHLTTPCFileHandler&)
{ 
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTTPCFileHandler::~AliHLTTPCFileHandler()
{
  //Destructor
  if(fMC) CloseMCOutput();
  FreeDigitsTree();
  if(fInAli) CloseAliInput();
}

// of course on start up the index is not created
Bool_t AliHLTTPCFileHandler::fgStaticIndexCreated=kFALSE;
Int_t  AliHLTTPCFileHandler::fgStaticIndex[36][159]; 

void AliHLTTPCFileHandler::CleanStaticIndex() 
{ 
  // use this static call to clean static index after
  // running over one event
  for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++){
    for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++)
      fgStaticIndex[i][j]=-1;
  }
  fgStaticIndexCreated=kFALSE;
}

Int_t AliHLTTPCFileHandler::SaveStaticIndex(Char_t *prefix,Int_t event) 
{ 
  // use this static call to store static index after
  if(!fgStaticIndexCreated) return -1;

  Char_t fname[1024];
  if(prefix)
    sprintf(fname,"%s-%d.txt",prefix,event);
  else
    sprintf(fname,"TPC.Digits.staticindex-%d.txt",event);

  ofstream file(fname,ios::trunc);
  if(!file.good()) return -1;

  for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++){
    for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++)
      file << fgStaticIndex[i][j] << " ";
    file << endl;
  }
  file.close();
  return 0;
}

Int_t AliHLTTPCFileHandler::LoadStaticIndex(Char_t *prefix,Int_t event) 
{ 
  // use this static call to store static index after
  if(fgStaticIndexCreated){
      LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::LoadStaticIndex","Inxed")
	<<"Static index already created, will overwrite"<<ENDLOG;
      CleanStaticIndex();
  }

  Char_t fname[1024];
  if(prefix)
    sprintf(fname,"%s-%d.txt",prefix,event);
  else
    sprintf(fname,"TPC.Digits.staticindex-%d.txt",event);

  ifstream file(fname);
  if(!file.good()) return -1;

  for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++){
    for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++)
      file >> fgStaticIndex[i][j];
  }
  file.close();

  fgStaticIndexCreated=kTRUE;
  return 0;
}

void AliHLTTPCFileHandler::FreeDigitsTree()
{ 
  //free digits tree
  if(!fDigitsTree)
    {
      LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::FreeDigitsTree()","Pointer")
	<<"Cannot free digitstree, it is not present"<<ENDLOG;
      return;
    }
  delete fDigits;
  fDigits=0;
#ifndef use_newio
  fDigitsTree->Delete();
#endif
  fDigitsTree=0;

  for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++){
    for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++)
      fIndex[i][j]=-1;
  }
  fIndexCreated=kFALSE;
}

Bool_t AliHLTTPCFileHandler::SetMCOutput(Char_t *name)
{ 
  //set mc input
  fMC = fopen(name,"w");
  if(!fMC){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetMCOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTTPCFileHandler::SetMCOutput(FILE *file)
{ 
  //set mc output
  fMC = file;
  if(!fMC){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetMCOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliHLTTPCFileHandler::CloseMCOutput()
{ 
  //close mc output
  if(!fMC){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::CloseMCOutPut","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fMC);
  fMC =0;
}

Bool_t AliHLTTPCFileHandler::SetAliInput()
{ 
  //set ali input
#ifdef use_newio
  fInAli->CdGAFile();
  fParam = (AliTPCParam*)gFile->Get("75x40_100x60_150x60");
  if(!fParam){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File")
      <<"No TPC parameters found in \""<<gFile->GetName()
      <<"\", creating standard parameters "
      <<"which might not be what you want!"<<ENDLOG;
    fParam = new AliTPCParamSR;
  }
  if(!fParam){ 
    LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::SetAliInput","File Open")
      <<"No AliTPCParam "<<AliHLTTPCTransform::GetParamName()<<" in File "<<gFile->GetName()<<ENDLOG;
    return kFALSE;
  }
#else
  if(!fInAli->IsOpen()){
    LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::SetAliInput","File Open")
      <<"Ali File "<<fInAli->GetName()<<" does not exist"<<ENDLOG;
    return kFALSE;
  }
  fParam = (AliTPCParam*)fInAli->Get(AliHLTTPCTransform::GetParamName());
  if(!fParam){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File")
      <<"No TPC parameters found in \""<<fInAli->GetName()
      <<"\", creating standard parameters "
      <<"which might not be what you want!"<<ENDLOG;
    fParam = new AliTPCParamSR;
  }
  if(!fParam){ 
    LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::SetAliInput","File Open")
      <<"No AliTPCParam "<<AliHLTTPCTransform::GetParamName()<<" in File "<<fInAli->GetName()<<ENDLOG;
    return kFALSE;
  }
#endif

  return kTRUE;
}

Bool_t AliHLTTPCFileHandler::SetAliInput(Char_t *name)
{ 
  //Open the AliROOT file with name.
#ifdef use_newio
  fInAli= AliRunLoader::Open(name);
#else
  fInAli= new TFile(name,"READ");
#endif
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File Open")
    <<"Pointer to fInAli = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

#ifdef use_newio
Bool_t AliHLTTPCFileHandler::SetAliInput(AliRunLoader *runLoader)
{ 
  //set ali input as runloader
  fInAli=runLoader;
  fUseRunLoader = kTRUE;
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File Open")
    <<"Pointer to AliRunLoader = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}
#endif

#ifdef use_newio
Bool_t AliHLTTPCFileHandler::SetAliInput(TFile */*file*/)
{
  //Specify already opened AliROOT file to use as an input.
  LOG(AliHLTTPCLog::kFatal,"AliHLTTPCFileHandler::SetAliInput","File Open")
    <<"This function is not supported for NEWIO, check ALIHLT_USENEWIO settings in Makefile.conf"<<ENDLOG;
  return kFALSE;
}
#else
Bool_t AliHLTTPCFileHandler::SetAliInput(TFile *file)
{ 
  //Specify already opened AliROOT file to use as an input.
  fInAli=file;
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File Open")
    <<"Pointer to fInAli = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}
#endif

void AliHLTTPCFileHandler::CloseAliInput()
{ 
  //close ali input
#ifdef use_newio
  if(fUseRunLoader) return;
#endif
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::CloseAliInput","RunLoader")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
#ifndef use_newio
  if(fInAli->IsOpen()) fInAli->Close();
#endif

  delete fInAli;
  fInAli = 0;
}

Bool_t AliHLTTPCFileHandler::IsDigit(Int_t event)
{
  //Check if there is a TPC digit tree in the current file.
  //Return kTRUE if tree was found, and kFALSE if not found.
  
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::IsDigit","File")
    <<"Pointer to fInAli = 0x0 "<<ENDLOG;
    return kTRUE;  //maybe you are using binary input which is Digits!
  }
#ifdef use_newio
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandlerNewIO::IsDigit","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(event);
  tpcLoader->LoadDigits();
  TTree *t=tpcLoader->TreeD();
#else
  Char_t name[1024];
  sprintf(name,"TreeD_%s_%d",AliHLTTPCTransform::GetParamName(),event);
  TTree *t=(TTree*)fInAli->Get(name);
#endif
  if(t){
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandlerNewIO::IsDigit","File Type")
    <<"Found Digit Tree -> Use Fast Cluster Finder"<<ENDLOG;
    return kTRUE;
  }
  else{
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandlerNewIO::IsDigit","File Type")
    <<"No Digit Tree -> Use Cluster Tree"<<ENDLOG;
    return kFALSE;
  }
}

///////////////////////////////////////// Digit IO  
Bool_t AliHLTTPCFileHandler::AliDigits2Binary(Int_t event,Bool_t altro)
{
  //save alidigits as binary
  Bool_t out = kTRUE;
  UInt_t nrow;
  AliHLTTPCDigitRowData* data = 0;
  if(altro)
    data = AliAltroDigits2Memory(nrow,event);
  else
    data = AliDigits2Memory(nrow,event);
  out = Memory2Binary(nrow,data);
  Free();
  return out;
}

Bool_t AliHLTTPCFileHandler::AliDigits2CompBinary(Int_t event,Bool_t altro)
{
  //Convert AliROOT TPC data, into HLT data format.
  //event specifies the event you want in the aliroot file.
  
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliHLTTPCDigitRowData *digits=0;
  if(altro)
    digits = AliAltroDigits2Memory(ndigits,event);
  else
    digits = AliDigits2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}

Bool_t AliHLTTPCFileHandler::CreateIndex()
{
  //create the access index or copy from static index
  fIndexCreated=kFALSE;

  if(!fgStaticIndexCreated || !fUseStaticIndex) { //we have to create index 
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::CreateIndex","Index")
      <<"Starting to create index, this can take a while."<<ENDLOG;

    for(Int_t n=0; n<fDigitsTree->GetEntries(); n++) {
      Int_t sector, row;
      Int_t lslice,lrow;
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      if(!AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row)){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::CreateIndex","Slice/Row")
	  <<AliHLTTPCLog::kDec<<"Index could not be created. Wrong values "
	  <<sector<<" "<<row<<ENDLOG;
	return kFALSE;
      }
      if(fIndex[lslice][lrow]==-1) {
	fIndex[lslice][lrow]=n;
      }
    }
    if(fUseStaticIndex) { // create static index
      for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++){
	for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++)
	  fgStaticIndex[i][j]=fIndex[i][j];
      }
      fgStaticIndexCreated=kTRUE; //remember that index has been created
    }

  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::CreateIndex","Index")
    <<"Index successfully created."<<ENDLOG;

  } else if(fUseStaticIndex) { //simply copy static index
    for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++){
      for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++)
	fIndex[i][j]=fgStaticIndex[i][j];
    }

  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::CreateIndex","Index")
    <<"Index successfully taken from static copy."<<ENDLOG;
  }
  fIndexCreated=kTRUE;
  return kTRUE;
}

AliHLTTPCDigitRowData * AliHLTTPCFileHandler::AliDigits2Memory(UInt_t & nrow,Int_t event)
{
  //Read data from AliROOT file into memory, and store it in the HLT data format.
  //Returns a pointer to the data.

  AliHLTTPCDigitRowData *data = 0;
  nrow=0;
  
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2Memory","File")
    <<"No Input avalible: Pointer to fInAli == NULL"<<ENDLOG;
    return 0; 
  }

#ifndef use_newio
  if(!fInAli->IsOpen()){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2Memory","File")
    <<"No Input avalible: TFile not opened"<<ENDLOG;
    return 0;
  }
#endif

  if(!fDigitsTree)
    if(!GetDigitsTree(event)) return 0;

  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t lslice,lrow;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)fDigitsTree->GetEntries();
  if(entries==0) {
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2Memory","ndigits")
      <<"No TPC digits (entries==0)!"<<ENDLOG;
    nrow = (UInt_t)(fRowMax-fRowMin+1);
    Int_t size = nrow*sizeof(AliHLTTPCDigitRowData);
    data=(AliHLTTPCDigitRowData*) Allocate(size);
    AliHLTTPCDigitRowData *tempPt = data;
    for(Int_t r=fRowMin;r<=fRowMax;r++){
      tempPt->fRow = r;
      tempPt->fNDigit = 0;
      tempPt++;
    }
    return data;
  }

  Int_t * ndigits = new Int_t[fRowMax+1];
  Float_t xyz[3];

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];
    if(n!=-1){ //data on that row
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);

      if(lrow!=r){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Rows in slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      ndigits[lrow] = 0;
      fDigits->First();
      do {
	time=fDigits->CurrentRow();
	pad=fDigits->CurrentColumn();
	dig = fDigits->GetDigit(time,pad);
	if(dig <= fParam->GetZeroSup()) continue;
	if(dig >= AliHLTTPCTransform::GetADCSat())
	  dig = AliHLTTPCTransform::GetADCSat();
      
	AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
	//	if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	//	  continue; // why 230???

	ndigits[lrow]++; //for this row only
	ndigitcount++;   //total number of digits to be published

      } while (fDigits->Next());
      //cout << lrow << " " << ndigits[lrow] << " - " << ndigitcount << endl;
    }
    nrows++;
  }

  Int_t size = sizeof(AliHLTTPCDigitData)*ndigitcount
    + nrows*sizeof(AliHLTTPCDigitRowData);

  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCFileHandler::AliDigits2Memory","Digits")
    <<AliHLTTPCLog::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliHLTTPCDigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliHLTTPCDigitRowData *tempPt = data;

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];
    tempPt->fRow = r;
    tempPt->fNDigit = 0;

    if(n!=-1){//data on that row
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
      if(lrow!=r){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Rows on slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      tempPt->fNDigit = ndigits[lrow];

      Int_t localcount=0;
      fDigits->First();
      do {
	time=fDigits->CurrentRow();
	pad=fDigits->CurrentColumn();
	dig = fDigits->GetDigit(time,pad);
	if (dig <= fParam->GetZeroSup()) continue;
	if(dig >= AliHLTTPCTransform::GetADCSat())
	  dig = AliHLTTPCTransform::GetADCSat();

	//Exclude data outside cone:
	AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
	//	if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	//	  continue; // why 230???

	if(localcount >= ndigits[lrow])
	  LOG(AliHLTTPCLog::kFatal,"AliHLTTPCFileHandler::AliDigits2Binary","Memory")
	    <<AliHLTTPCLog::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	    <<ndigits[lrow]<<ENDLOG;
	
	tempPt->fDigitData[localcount].fCharge=dig;
	tempPt->fDigitData[localcount].fPad=pad;
	tempPt->fDigitData[localcount].fTime=time;
#ifdef do_mc
	tempPt->fDigitData[localcount].fTrackID[0] = fDigits->GetTrackID(time,pad,0);
	tempPt->fDigitData[localcount].fTrackID[1] = fDigits->GetTrackID(time,pad,1);
	tempPt->fDigitData[localcount].fTrackID[2] = fDigits->GetTrackID(time,pad,2);
#endif
	localcount++;
      } while (fDigits->Next());
    }
    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliHLTTPCDigitRowData)
                                      + ndigits[lrow]*sizeof(AliHLTTPCDigitData);
    tmp += size;
    tempPt = (AliHLTTPCDigitRowData*)tmp;
  }
  delete [] ndigits;
  return data;
}

AliHLTTPCDigitRowData * AliHLTTPCFileHandler::AliAltroDigits2Memory(UInt_t & nrow,Int_t event,Bool_t eventmerge)
{
  //Read data from AliROOT file into memory, and store it in the HLT data format.
  //Returns a pointer to the data.
  //This functions filter out single timebins, which is noise. The timebins which
  //are removed are timebins which have the 4 zero neighbours; 
  //(pad-1,time),(pad+1,time),(pad,time-1),(pad,time+1).
  
  AliHLTTPCDigitRowData *data = 0;
  nrow=0;
  
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliAltroDigits2Memory","File")
    <<"No Input avalible: Pointer to TFile == NULL"<<ENDLOG;
    return 0; 
  }
#ifndef use_newio
  if(!fInAli->IsOpen()){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliAltroDigits2Memory","File")
    <<"No Input avalible: TFile not opened"<<ENDLOG;
    return 0;
  }
#endif
  if(eventmerge == kTRUE && event >= 1024)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliAltroDigits2Memory","TrackIDs")
	<<"Too many events if you want to merge!"<<ENDLOG;
      return 0;
    }
  delete fDigits;
  fDigits=0;
#ifdef use_newio 
  /* Dont understand why we have to do 
     reload the tree, but otherwise the code crashes */
  fDigitsTree=0;
  if(!GetDigitsTree(event)) return 0;
#else
  if(!fDigitsTree){
    if(!GetDigitsTree(event)) return 0;
  }
#endif

  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)fDigitsTree->GetEntries();
  if(entries==0) {
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliAltroDigits2Memory","ndigits")
      <<"No TPC digits (entries==0)!"<<ENDLOG;
    nrow = (UInt_t)(fRowMax-fRowMin+1);
    Int_t size = nrow*sizeof(AliHLTTPCDigitRowData);
    data=(AliHLTTPCDigitRowData*) Allocate(size);
    AliHLTTPCDigitRowData *tempPt = data;
    for(Int_t r=fRowMin;r<=fRowMax;r++){
      tempPt->fRow = r;
      tempPt->fNDigit = 0;
      tempPt++;
    }
    return data;
  }
  Int_t * ndigits = new Int_t[fRowMax+1];
  Int_t lslice,lrow;
  Int_t zerosupval=AliHLTTPCTransform::GetZeroSup();
  Float_t xyz[3];

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];

    ndigits[r] = 0;

    if(n!=-1){//data on that row
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
      //cout << lslice << " " << fSlice << " " << lrow << " " << r << " " << sector << " " << row << endl;
      if((lslice!=fSlice)||(lrow!=r)){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliAltroDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Rows on slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      fDigits->ExpandBuffer();
      fDigits->ExpandTrackBuffer();
      for(Int_t i=0; i<fDigits->GetNCols(); i++){
	for(Int_t j=0; j<fDigits->GetNRows(); j++){
	  pad=i;
	  time=j;
	  dig = fDigits->GetDigitFast(time,pad);
	  if(dig <= zerosupval) continue;
	  if(dig >= AliHLTTPCTransform::GetADCSat())
	    dig = AliHLTTPCTransform::GetADCSat();

	  //Check for single timebins, and remove them because they are noise for sure.
	  if(i>0 && i<fDigits->GetNCols()-1 && j>0 && j<fDigits->GetNRows()-1)
	    if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
	       fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
	       fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
	       fDigits->GetDigitFast(time,pad+1)<=zerosupval)
	      continue;
	      
	  //Boundaries:
	  if(i==0) //pad==0
	    {
	      if(j < fDigits->GetNRows()-1 && j > 0) 
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval)
		    continue;
		}
	      else if(j > 0)
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval)
		    continue;
		}
	    }
	  if(j==0)
	    {
	      if(i < fDigits->GetNCols()-1 && i > 0)
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval)
		    continue;
		}
	      else if(i > 0)
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval)
		    continue;
		}
	    }

	  if(i==fDigits->GetNCols()-1)
	    {
	      if(j>0 && j<fDigits->GetNRows()-1)
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad-1)<=zerosupval)
		    continue;
		}
	      else if(j==0 && j<fDigits->GetNRows()-1)
		{
		  if(fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad-1)<=zerosupval)
		    continue;
		}
	      else 
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad-1)<=zerosupval)
		    continue;
		}
	    }
	
	  if(j==fDigits->GetNRows()-1)
	    {
	      if(i>0 && i<fDigits->GetNCols()-1)
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval &&
		     fDigits->GetDigitFast(time-1,pad)<=zerosupval)
		    continue;
		}
	      else if(i==0 && fDigits->GetNCols()-1)
		{
		  if(fDigits->GetDigitFast(time,pad+1)<=zerosupval &&
		     fDigits->GetDigitFast(time-1,pad)<=zerosupval)
		    continue;
		}
	      else 
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time-1,pad)<=zerosupval)
		    continue;
		}
	    }

	  AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
	  //	  if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	  //	  continue; 
	      
	  ndigits[lrow]++; //for this row only
	  ndigitcount++;   //total number of digits to be published
	}
      }
    }
    nrows++;
  }
  
  Int_t size = sizeof(AliHLTTPCDigitData)*ndigitcount
    + nrows*sizeof(AliHLTTPCDigitRowData);

  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCFileHandler::AliAltroDigits2Memory","Digits")
    <<AliHLTTPCLog::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliHLTTPCDigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliHLTTPCDigitRowData *tempPt = data;
 
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];
    tempPt->fRow = r;
    tempPt->fNDigit = 0;
    if(n!=-1){ //no data on that row
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);

      if(lrow!=r){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliAltroDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Rows on slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      tempPt->fNDigit = ndigits[lrow];

      Int_t localcount=0;
      fDigits->ExpandBuffer();
      fDigits->ExpandTrackBuffer();
      for(Int_t i=0; i<fDigits->GetNCols(); i++){
	for(Int_t j=0; j<fDigits->GetNRows(); j++){
	  pad=i;
	  time=j;
	  dig = fDigits->GetDigitFast(time,pad);
	  if(dig <= zerosupval) continue;
	  if(dig >= AliHLTTPCTransform::GetADCSat())
	    dig = AliHLTTPCTransform::GetADCSat();
	      
	  //Check for single timebins, and remove them because they are noise for sure.
	  if(i>0 && i<fDigits->GetNCols()-1 && j>0 && j<fDigits->GetNRows()-1)
	    if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
	       fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
	       fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
	       fDigits->GetDigitFast(time,pad+1)<=zerosupval)
	      continue;
	  
	  //Boundaries:
	  if(i==0) //pad ==0
	    {
	      if(j < fDigits->GetNRows()-1 && j > 0) 
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval)
		    continue;
		}
	      else if(j > 0)
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval)
		    continue;
		}
	    }
	  if(j==0)
	    {
	      if(i < fDigits->GetNCols()-1 && i > 0)
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval)
		    continue;
		}
	      else if(i > 0)
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval)
		    continue;
		}
	    }
	
	  if(i == fDigits->GetNCols()-1)
	    {
	      if(j>0 && j<fDigits->GetNRows()-1)
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad-1)<=zerosupval)
		    continue;
		}
	      else if(j==0 && j<fDigits->GetNRows()-1)
		{
		  if(fDigits->GetDigitFast(time+1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad-1)<=zerosupval)
		    continue;
		}
	      else 
		{
		  if(fDigits->GetDigitFast(time-1,pad)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad-1)<=zerosupval)
		    continue;
		}
	    }
	  if(j==fDigits->GetNRows()-1)
	    {
	      if(i>0 && i<fDigits->GetNCols()-1)
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time,pad+1)<=zerosupval &&
		     fDigits->GetDigitFast(time-1,pad)<=zerosupval)
		    continue;
		}
	      else if(i==0 && fDigits->GetNCols()-1)
		{
		  if(fDigits->GetDigitFast(time,pad+1)<=zerosupval &&
		     fDigits->GetDigitFast(time-1,pad)<=zerosupval)
		    continue;
		}
	      else 
		{
		  if(fDigits->GetDigitFast(time,pad-1)<=zerosupval &&
		     fDigits->GetDigitFast(time-1,pad)<=zerosupval)
		    continue;
		}
	    }
	
	  AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
	  //	  if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	  //	    continue;
	  
	  if(localcount >= ndigits[lrow])
	    LOG(AliHLTTPCLog::kFatal,"AliHLTTPCFileHandler::AliAltroDigits2Binary","Memory")
	      <<AliHLTTPCLog::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	      <<ndigits[lrow]<<ENDLOG;
	
	  tempPt->fDigitData[localcount].fCharge=dig;
	  tempPt->fDigitData[localcount].fPad=pad;
	  tempPt->fDigitData[localcount].fTime=time;
#ifdef do_mc
	  tempPt->fDigitData[localcount].fTrackID[0] = (fDigits->GetTrackIDFast(time,pad,0)-2);
	  tempPt->fDigitData[localcount].fTrackID[1] = (fDigits->GetTrackIDFast(time,pad,1)-2);
	  tempPt->fDigitData[localcount].fTrackID[2] = (fDigits->GetTrackIDFast(time,pad,2)-2);
	  if(eventmerge == kTRUE) //careful track mc info will be touched
	    {//Event are going to be merged, so event number is stored in the upper 10 bits.
	      tempPt->fDigitData[localcount].fTrackID[0] += 128; //leave some room
	      tempPt->fDigitData[localcount].fTrackID[1] += 128; //for neg. numbers
	      tempPt->fDigitData[localcount].fTrackID[2] += 128;
	      tempPt->fDigitData[localcount].fTrackID[0] += ((event&0x3ff)<<22);
	      tempPt->fDigitData[localcount].fTrackID[1] += ((event&0x3ff)<<22);
	      tempPt->fDigitData[localcount].fTrackID[2] += ((event&0x3ff)<<22);
	    }
#endif
	  localcount++;
	}
      }
    }
    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliHLTTPCDigitRowData)
      + ndigits[r]*sizeof(AliHLTTPCDigitData);
    tmp += size;
    tempPt = (AliHLTTPCDigitRowData*)tmp;
  }
  delete [] ndigits;
  return data;
}
 
Bool_t AliHLTTPCFileHandler::GetDigitsTree(Int_t event)
{
  //Connects to the TPC digit tree in the AliROOT file.
#ifdef use_newio
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::GetDigitsTree","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(event);
  tpcLoader->LoadDigits();
  fDigitsTree = tpcLoader->TreeD();
#else  
  fInAli->cd();
  Char_t dname[100];
  sprintf(dname,"TreeD_%s_%d",AliHLTTPCTransform::GetParamName(),event);
  fDigitsTree = (TTree*)fInAli->Get(dname);
#endif
  if(!fDigitsTree) 
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::GetDigitsTree","Digits Tree")
	<<AliHLTTPCLog::kHex<<"Error getting digitstree "<<(void*)fDigitsTree<<ENDLOG;
      return kFALSE;
    }
  fDigitsTree->GetBranch("Segment")->SetAddress(&fDigits);

  if(!fIndexCreated) return CreateIndex();
  else return kTRUE;
}

void AliHLTTPCFileHandler::AliDigits2RootFile(AliHLTTPCDigitRowData *rowPt,Char_t *new_digitsfile)
{
  //Write the data stored in rowPt, into a new AliROOT file.
  //The data is stored in the AliROOT format 
  //This is specially a nice thing if you have modified data, and wants to run it  
  //through the offline reconstruction chain.
  //The arguments is a pointer to the data, and the name of the new AliROOT file.
  //Remember to pass the original AliROOT file (the one that contains the original
  //simulated data) to this object, in order to retrieve the MC id's of the digits.

  if(!fInAli)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2RootFile","File")
	<<"No rootfile "<<ENDLOG;
      return;
    }
  if(!fParam)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2RootFile","File")
	<<"No parameter object. Run on rootfile "<<ENDLOG;
      return;
    }

#ifdef use_newio
  //Get the original digitstree:
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2RootFile","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return;
  }
  tpcLoader->LoadDigits();
  TTree *t=tpcLoader->TreeD();

  AliTPCDigitsArray *old_array = new AliTPCDigitsArray();
  old_array->Setup(fParam);
  old_array->SetClass("AliSimDigits");

  Bool_t ok = old_array->ConnectTree(t);
  if(!ok)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2RootFile","File")
	<< "No digits tree object" << ENDLOG;
      return;
    }

  tpcLoader->SetDigitsFileName(new_digitsfile);
  tpcLoader->MakeDigitsContainer();
    
  //setup a new one, or connect it to the existing one:
  AliTPCDigitsArray *arr = new AliTPCDigitsArray(); 
  arr->SetClass("AliSimDigits");
  arr->Setup(fParam);
  arr->MakeTree(tpcLoader->TreeD());
#else
  
  //Get the original digitstree:
  Char_t dname[100];
  sprintf(dname,"TreeD_%s_0",AliHLTTPCTransform::GetParamName());

  fInAli->cd();
  AliTPCDigitsArray *old_array = new AliTPCDigitsArray();
  old_array->Setup(fParam);
  old_array->SetClass("AliSimDigits");

  Bool_t ok = old_array->ConnectTree(dname);
  if(!ok)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2RootFile","File")
	<<"No digits tree object." <<ENDLOG;
      return;
    }

  Bool_t create=kFALSE;
  TFile *digFile;
  
  if(gSystem->AccessPathName(new_digitsfile))
    {
      LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::AliDigits2RootFile","File")
	<<"Creating new file "<<new_digitsfile<<ENDLOG;
      create = kTRUE;
      digFile = TFile::Open(new_digitsfile,"RECREATE");
      fParam->Write(fParam->GetTitle());
    }
  else
    {
      create = kFALSE;
      digFile = TFile::Open(new_digitsfile,"UPDATE");
      
    }
  if(!digFile->IsOpen())
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2RootFile","Rootfile")
	<<"Error opening rootfile "<<new_digitsfile<<ENDLOG;
      return;
    }
  
  digFile->cd();
    
  //setup a new one, or connect it to the existing one:
  AliTPCDigitsArray *arr = new AliTPCDigitsArray(); 
  arr->SetClass("AliSimDigits");
  arr->Setup(fParam);
  if(create)
    arr->MakeTree();
  else
    {
      Bool_t ok = arr->ConnectTree(dname);
      if(!ok)
	{
	  LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2RootFile","Rootfile")
	    <<"No digits tree object in existing file"<<ENDLOG;
	  return;
	}
    }
#endif

  Int_t digcounter=0,trackID[3];

  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      
      if((Int_t)rowPt->fRow != i) 
	LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2RootFile","Data")
	  <<"Mismatching row numbering "<<(Int_t)rowPt->fRow<<" "<<i<<ENDLOG;
            
      Int_t sector,row;
      AliHLTTPCTransform::Slice2Sector(fSlice,i,sector,row);
      
      AliSimDigits *old_dig = (AliSimDigits*)old_array->LoadRow(sector,row);
      AliSimDigits * dig = (AliSimDigits*)arr->CreateRow(sector,row);
      old_dig->ExpandBuffer();
      old_dig->ExpandTrackBuffer();
      dig->ExpandBuffer();
      dig->ExpandTrackBuffer();
      
      if(!old_dig)
	LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2RootFile","Data")
	  <<"No padrow " << sector << " " << row <<ENDLOG;

      AliHLTTPCDigitData *digPt = rowPt->fDigitData;
      digcounter=0;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  Short_t charge = (Short_t)digPt[j].fCharge;
	  Int_t pad = (Int_t)digPt[j].fPad;
	  Int_t time = (Int_t)digPt[j].fTime;
	  
	  if(charge == 0) //Only write the digits that has not been removed
	    {
	      LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2RootFile","Data")
		<<"Zero charge" <<ENDLOG;
	      continue;
	    }

	  digcounter++;
	  
	  //Tricks to get and set the correct track id's. 
	  for(Int_t t=0; t<3; t++)
	    {
	      Int_t label = old_dig->GetTrackIDFast(time,pad,t);
	      if(label > 1)
		trackID[t] = label - 2;
	      else if(label==0)
		trackID[t] = -2;
	      else
		trackID[t] = -1;
	    }
	  
	  dig->SetDigitFast(charge,time,pad);
	  
	  for(Int_t t=0; t<3; t++)
	    ((AliSimDigits*)dig)->SetTrackIDFast(trackID[t],time,pad,t);
	  
	}
      //cout<<"Wrote "<<digcounter<<" on row "<<i<<endl;
      UpdateRowPointer(rowPt);
      arr->StoreRow(sector,row);
      arr->ClearRow(sector,row);  
      old_array->ClearRow(sector,row);
    }

  char treeName[100];
  sprintf(treeName,"TreeD_%s_0",fParam->GetTitle());
  
#ifdef use_newio
  arr->GetTree()->SetName(treeName);
  arr->GetTree()->AutoSave();
  tpcLoader->WriteDigits("OVERWRITE");
#else
  digFile->cd();
  arr->GetTree()->SetName(treeName);
  arr->GetTree()->AutoSave();
  digFile->Close();
#endif
  delete arr;
  delete old_array;
}

///////////////////////////////////////// Point IO  
Bool_t AliHLTTPCFileHandler::AliPoints2Binary(Int_t eventn)
{
  //points to binary
  Bool_t out = kTRUE;
  UInt_t npoint;
  AliHLTTPCSpacePointData *data = AliPoints2Memory(npoint,eventn);
  out = Memory2Binary(npoint,data);
  Free();
  return out;
}

AliHLTTPCSpacePointData * AliHLTTPCFileHandler::AliPoints2Memory(UInt_t & npoint,Int_t eventn)
{
  //points to memory
  AliHLTTPCSpacePointData *data = 0;
  npoint=0;
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliPoints2Memory","File")
    <<"No Input avalible: no object fInAli"<<ENDLOG;
    return 0;
  }
#ifndef use_newio
  if(!fInAli->IsOpen()){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliPoints2Memory","File")
    <<"No Input avalible: TFile not opend"<<ENDLOG;
    return 0;
  }
#endif

  TDirectory *savedir = gDirectory;
#ifdef use_newio
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliPoints2Memory","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return 0;
  }
  fInAli->GetEvent(eventn);
  tpcLoader->LoadRecPoints();

  AliTPCClustersArray carray;
  carray.Setup(fParam);
  carray.SetClusterType("AliTPCcluster");
  Bool_t clusterok = carray.ConnectTree(tpcLoader->TreeR());
#else
  fInAli->cd();
  
  Char_t cname[100];
  sprintf(cname,"TreeC_TPC_%d",eventn);
  AliTPCClustersArray carray;
  carray.Setup(fParam);
  carray.SetClusterType("AliTPCcluster");
  Bool_t clusterok = carray.ConnectTree(cname);
#endif

  if(!clusterok) return 0;

  AliTPCClustersRow ** clusterrow = 
               new AliTPCClustersRow*[ (int)carray.GetTree()->GetEntries()];
  Int_t *rows = new int[ (int)carray.GetTree()->GetEntries()];
  Int_t *sects = new int[  (int)carray.GetTree()->GetEntries()];
  Int_t sum=0;

  Int_t lslice,lrow;
  for(Int_t i=0; i<carray.GetTree()->GetEntries(); i++){
    AliSegmentID *s = carray.LoadEntry(i);
    Int_t sector,row;
    fParam->AdjustSectorRow(s->GetID(),sector,row);
    rows[i] = row;
    sects[i] = sector;
    clusterrow[i] = 0;
    AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
    if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
    clusterrow[i] = carray.GetRow(sector,row);
    if(clusterrow[i])
      sum+=clusterrow[i]->GetArray()->GetEntriesFast();
  }
  UInt_t size = sum*sizeof(AliHLTTPCSpacePointData);

  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCFileHandler::AliPoints2Memory","File")
  <<AliHLTTPCLog::kDec<<"Found "<<sum<<" SpacePoints"<<ENDLOG;

  data = (AliHLTTPCSpacePointData *) Allocate(size);
  npoint = sum;
  UInt_t n=0; 
  Int_t pat=fPatch;
  if(fPatch==-1)
    pat=0;
  for(Int_t i=0; i<carray.GetTree()->GetEntries(); i++){
    if(!clusterrow[i]) continue;
    Int_t row = rows[i];
    Int_t sector = sects[i];
    AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
    Int_t entries_in_row = clusterrow[i]->GetArray()->GetEntriesFast();
    for(Int_t j = 0;j<entries_in_row;j++){
      AliTPCcluster *c = (AliTPCcluster*)(*clusterrow[i])[j];
      data[n].fZ = c->GetZ();
      data[n].fY = c->GetY();
      data[n].fX = fParam->GetPadRowRadii(sector,row);
      data[n].fCharge = (UInt_t)c->GetQ();
      data[n].fID = n+((fSlice&0x7f)<<25)+((pat&0x7)<<22);//uli
      data[n].fPadRow = lrow;
      data[n].fSigmaY2 = c->GetSigmaY2();
      data[n].fSigmaZ2 = c->GetSigmaZ2();
#ifdef do_mc
      data[n].fTrackID[0] = c->GetLabel(0);
      data[n].fTrackID[1] = c->GetLabel(1);
      data[n].fTrackID[2] = c->GetLabel(2);
#endif
      if(fMC) fprintf(fMC,"%d %d\n",data[n].fID,c->GetLabel(0));
      n++;
    }
  }
  for(Int_t i=0;i<carray.GetTree()->GetEntries();i++){
    Int_t row = rows[i];
    Int_t sector = sects[i];
    if(carray.GetRow(sector,row))
      carray.ClearRow(sector,row);
  }

  delete [] clusterrow;
  delete [] rows;
  delete [] sects;
  savedir->cd();   

  return data;
}

