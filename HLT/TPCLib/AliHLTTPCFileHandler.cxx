// $Id$
// Original: AliHLTFileHandler.cxx,v 1.49 2005/06/23 17:46:55 hristov 

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: U. Frankenfeld, A. Vestbo, C. Loizides                *
//*                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTPCFileHandler.cxx
/// @author U. Frankenfeld, A. Vestbo, C. Loizides, maintained by
///         Matthias Richter
/// @date   
/// @brief  file input for the TPC tracking code before migration to the
///         HLT component framework

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include <TClonesArray.h>
#include <TSystem.h>
#include <TMath.h>

#include "AliLoader.h"
#include <AliRunLoader.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <AliTPCParamSR.h>
#include <AliTPCDigitsArray.h>
#include <AliTPCclusterMI.h>
#include <AliTPCClustersRow.h>
#include <AliSimDigits.h>
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDigitData.h"
//#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCSpacePointData.h"
//#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCFileHandler.h"
#include "AliHLTTPCMapping.h"
#include "AliHLTAltroEncoder.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/**
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
  fUseRunLoader(kFALSE),
  fParam(NULL),
  fDigits(NULL),
  fDigitsTree(NULL),
  fMC(NULL),
  fIndexCreated(kFALSE),
  fUseStaticIndex(b)
{
  //Default constructor

  for(Int_t i=0;i<AliHLTTPCTransform::GetNSlice();i++)
    for(Int_t j=0;j<AliHLTTPCTransform::GetNRows();j++) 
      fIndex[i][j]=-1;

  if(fUseStaticIndex&&!fgStaticIndexCreated) CleanStaticIndex();
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

  const int fnamelen=1024;
  Char_t fname[fnamelen];
  if(prefix)
    snprintf(fname,fnamelen, "%s-%d.txt",prefix,event);
  else
    snprintf(fname,fnamelen, "TPC.Digits.staticindex-%d.txt",event);

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

  const int fnamelen=1024;
  Char_t fname[fnamelen];
  if(prefix)
    snprintf(fname,fnamelen,"%s-%d.txt",prefix,event);
  else
    snprintf(fname,fnamelen,"TPC.Digits.staticindex-%d.txt",event);

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
  if (fDigits) {
    delete fDigits;
  }
  fDigits=0;
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

  // fParam is in all cases an external object
  const char* cdbEntry="TPC/Calib/Parameters";
  AliCDBManager* pMan=AliCDBManager::Instance();
  if (pMan) {
    AliCDBEntry *pEntry = pMan->Get(cdbEntry);
    if (pEntry==NULL || 
	pEntry->GetObject()==NULL ||
	(fParam=dynamic_cast<AliTPCParam*>(pEntry->GetObject()))==NULL) {
      LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File")
	<<"can not load AliTPCParam object from OCDB entry " << cdbEntry <<ENDLOG;
    }
  }
  if (!fParam) {
  // the old solution until Nov 2008
  fInAli->CdGAFile();
  fParam = (AliTPCParam*)gFile->Get("75x40_100x60_150x60");
  if(!fParam){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File")
      <<"No TPC parameters found in \""<<gFile->GetName()
      <<"\", creating standard parameters "
      <<"which might not be what you want!"<<ENDLOG;
    fParam = new AliTPCParamSR;
  }
  }
  if(!fParam){ 
    LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::SetAliInput","File Open")
      <<"No AliTPCParam "<<AliHLTTPCTransform::GetParamName()<<" in File "<<gFile->GetName()<<ENDLOG;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliHLTTPCFileHandler::SetAliInput(Char_t *name)
{ 
  //Open the AliROOT file with name.
  fInAli= AliRunLoader::Open(name);
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File Open")
    <<"Pointer to fInAli = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

Bool_t AliHLTTPCFileHandler::SetAliInput(AliRunLoader *runLoader)
{ 
  //set ali input as runloader
  if(!runLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::SetAliInput","File Open")
      <<"invalid agument: pointer to AliRunLoader NULL "<<ENDLOG;
    return kFALSE;
  }
  if (fInAli!=NULL && fInAli!=runLoader) {
    LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::SetAliInput","File Open")
    <<"Pointer to AliRunLoader already set"<<ENDLOG;
    return kFALSE;
  }
  fInAli=runLoader;
  fUseRunLoader = kTRUE;
  return SetAliInput();
}

void AliHLTTPCFileHandler::CloseAliInput()
{ 
  //close ali input
  if(fUseRunLoader) return;
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::CloseAliInput","RunLoader")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }

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
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandlerNewIO::IsDigit","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(event);
  tpcLoader->LoadDigits();
  TTree *t=tpcLoader->TreeD();
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
Bool_t AliHLTTPCFileHandler::AliDigits2BinaryFile(Int_t event,Bool_t altro)
{
  //save alidigits as binary
  Bool_t out = kTRUE;
  UInt_t nrow;
  AliHLTTPCDigitRowData* data = 0;
  if(altro)
    data = AliAltroDigits2Memory(nrow,event);
  else
    data = AliDigits2Memory(nrow,event);
  out = Memory2BinaryFile(nrow,data);
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

    //Int_t lslice,lrow;
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
      // this is just to make sure the same array dimensions are
      // used as in the AliHLTTPCTransform class. The check for
      // correct bounds is done in AliHLTTPCTransform::Sector2Slice
      assert(lslice>=0 && lslice<fgkNSlice);
      assert(lrow>=0 && lrow<fgkNRow);
      if(fIndex[lslice][lrow]==-1) {
	fIndex[lslice][lrow]=n;
      }
    }
    assert(AliHLTTPCTransform::GetNSlice()==fgkNSlice);
    assert(AliHLTTPCTransform::GetNRows()==fgkNRow);
    if(fUseStaticIndex) { // create static index
      for(Int_t islice=0;islice<AliHLTTPCTransform::GetNSlice() && islice<fgkNSlice;islice++){
	for(Int_t irow=0;irow<AliHLTTPCTransform::GetNRows() && irow<fgkNRow;irow++)
	  fgStaticIndex[islice][irow]=fIndex[islice][irow];
      }
      fgStaticIndexCreated=kTRUE; //remember that index has been created
    }

  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::CreateIndex","Index")
    <<"Index successfully created."<<ENDLOG;

  } else if(fUseStaticIndex) { //simply copy static index
    for(Int_t islice=0;islice<AliHLTTPCTransform::GetNSlice() && islice<fgkNSlice;islice++){
      for(Int_t irow=0;irow<AliHLTTPCTransform::GetNRows() && irow<fgkNRow;irow++)
	fIndex[islice][irow]=fgStaticIndex[islice][irow];
    }

  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::CreateIndex","Index")
    <<"Index successfully taken from static copy."<<ENDLOG;
  }
  fIndexCreated=kTRUE;
  return kTRUE;
}

AliHLTTPCDigitRowData * AliHLTTPCFileHandler::AliDigits2Memory(UInt_t & nrow,Int_t event, Byte_t* tgtBuffer, UInt_t *pTgtSize)
{
  //Read data from AliROOT file into memory, and store it in the HLT data format
  //in the provided buffer or an allocated buffer.
  //Returns a pointer to the data.

  AliHLTTPCDigitRowData *data = 0;
  nrow=0;
  
  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2Memory","File")
    <<"No Input avalible: Pointer to fInAli == NULL"<<ENDLOG;
    return 0; 
  }

  if(!fDigitsTree)
    if(!GetDigitsTree(event)) return 0;

  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)fDigitsTree->GetEntries();
  if(entries==0) {
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2Memory","ndigits")
      <<"No TPC digits (entries==0)!"<<ENDLOG;
    nrow = (UInt_t)(fRowMax-fRowMin+1);
    UInt_t size = nrow*sizeof(AliHLTTPCDigitRowData);
    if (tgtBuffer!=NULL && pTgtSize!=NULL && *pTgtSize>0) {
      if (size<=*pTgtSize) {
	data=reinterpret_cast<AliHLTTPCDigitRowData*>(tgtBuffer);
      } else {
      }
    } else {
      data=reinterpret_cast<AliHLTTPCDigitRowData*>(Allocate(size));
    }
    AliHLTTPCDigitRowData *tempPt = data;
    if (data) {
    if (pTgtSize) *pTgtSize=size;
    for(Int_t r=fRowMin;r<=fRowMax;r++){
      tempPt->fRow = r;
      tempPt->fNDigit = 0;
      tempPt++;
    }
    }
    return data;
  }

  Int_t * ndigits = new Int_t[fRowMax+1];
  memset(ndigits, 0, (fRowMax+1)*sizeof(Int_t));
  Float_t xyz[3];

  // The digits of the current event have been indexed: all digits are organized in
  // rows, all digits of one row are stored in a AliSimDigits object (fDigit) which
  // are stored in the digit tree.
  // The index map relates the AliSimDigits objects in the tree to dedicated pad rows
  // in the TPC
  // This loop filters the pad rows according to the slice no set via Init
  assert(fRowMax<fgkNRow);
  for(Int_t r=fRowMin;r<=fRowMax && r<fgkNRow;r++){
    Int_t n=fIndex[fSlice][r];
    if(n!=-1){ // there is data on that row available
      Int_t lslice,lrow;
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
//       LOG(AliHLTTPCLog::kInformational,"AliHLTTPCFileHandler::AliDigits2Memory","Digits")
// 	<< "Sector "<<sector<<" Row " << row << " Slice " << lslice << " lrow " << lrow<<ENDLOG;

      if(lrow!=r){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Rows in slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      ndigits[lrow] = 0;
      for (bool bHaveData=fDigits->First();
	   bHaveData;
	   bHaveData=fDigits->Next()) {
	time=fDigits->CurrentRow();
	pad=fDigits->CurrentColumn();
	dig = fDigits->GetDigit(time,pad);
	if(dig <= fParam->GetZeroSup()) continue;
	if(dig >= AliHLTTPCTransform::GetADCSat())
	  dig = AliHLTTPCTransform::GetADCSat();

	// if we switch to AliTPCTransform, this maybe needs to be 
	// adjusted as well
	AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
	//	if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	//	  continue; // why 230???

	ndigits[lrow]++; //for this row only
	ndigitcount++;   //total number of digits to be published

      }
      //cout << lrow << " " << ndigits[lrow] << " - " << ndigitcount << endl;
    }
    //see comment below//nrows++;
  }
  // Matthias 05.11.2007
  // The question is whether we should always return a AliHLTTPCDigitRowData
  // for each row, even the empty ones or only for the ones filled with data.
  // the AliHLTTPCDigitReaderUnpacked as the counnterpart so far assumes 
  // empty RawData structs for empty rows. But some of the code here implies
  // the latter approach, e.g. the count of nrows in the loop above (now
  // commented). At least the two loops were not consistent, it's fixed now.
  nrows=fRowMax-fRowMin+1;

  UInt_t bufferSize = sizeof(AliHLTTPCDigitData)*ndigitcount
    + nrows*sizeof(AliHLTTPCDigitRowData);

  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCFileHandler::AliDigits2Memory","Digits")
    << "Found "<<ndigitcount<<" Digits in " << nrows << " rows out of [" << fRowMin << "," << fRowMax <<"]"<<ENDLOG;

  if (tgtBuffer!=NULL && pTgtSize!=NULL && *pTgtSize>0) {
    if (bufferSize<=*pTgtSize) {
      data=reinterpret_cast<AliHLTTPCDigitRowData*>(tgtBuffer);
    } else {
    }
  } else if (bufferSize>0) {
    data=reinterpret_cast<AliHLTTPCDigitRowData*>(Allocate(bufferSize));
  }
  if (pTgtSize) *pTgtSize=bufferSize;
  if (data==NULL) {
    delete [] ndigits;
    return NULL;
  }
  nrow = (UInt_t)nrows;
  AliHLTTPCDigitRowData *tempPt = data;
  memset(data, 0, bufferSize);

  for(Int_t r=fRowMin;r<=fRowMax && r<fgkNRow;r++){
    Int_t n=fIndex[fSlice][r];

    AliHLTTPCTransform::Slice2Sector(fSlice,r,sector,row);
    tempPt->fRow = row;
    tempPt->fNDigit = 0;

    if(n!=-1){//data on that row
      Int_t lslice,lrow;
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
      if(lrow!=r){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Rows on slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      // set the correct row no and digit count
      tempPt->fRow = row;
      tempPt->fNDigit = ndigits[lrow];

      Int_t localcount=0;
      for (bool bHaveData=fDigits->First();
	   bHaveData;
	   bHaveData=fDigits->Next()) {
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
	tempPt->fDigitData[localcount].fTrackID[0] = fDigits->GetTrackID(time,pad,0);
	tempPt->fDigitData[localcount].fTrackID[1] = fDigits->GetTrackID(time,pad,1);
	tempPt->fDigitData[localcount].fTrackID[2] = fDigits->GetTrackID(time,pad,2);
	localcount++;
      }
    }

    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t blockSize = sizeof(AliHLTTPCDigitRowData)
      + tempPt->fNDigit*sizeof(AliHLTTPCDigitData);
    tmp += blockSize;
    tempPt = (AliHLTTPCDigitRowData*)tmp;
  }
  assert((Byte_t*)tempPt==((Byte_t*)data)+bufferSize);
  delete [] ndigits;
  return data;
}

int AliHLTTPCFileHandler::AliDigits2Altro(Int_t event, Byte_t* tgtBuffer, UInt_t tgtSize)
{
  //Read data from AliROOT file into memory, and store it in the ALTRO data format
  //in the provided buffer 
  //Returns: size of the encoded data in bytes
  int iResult=0;

  if(!fInAli){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::AliDigits2Altro","File")
    <<"No Input avalible: Pointer to fInAli == NULL"<<ENDLOG;
    return 0; 
  }

  if (!tgtBuffer) {
    return -EINVAL;
  }

  if(!fDigitsTree)
    if(!GetDigitsTree(event)) return 0;

  UShort_t dig;
  Int_t time=0;
  Int_t pad=0;
  Int_t sector=0;

  AliHLTTPCMapping mapper(fPatch);
  AliHLTAltroEncoder encoder;
  encoder.SetBuffer(tgtBuffer, tgtSize);

  // The digits of the current event have been indexed: all digits are organized in
  // rows, all digits of one row are stored in a AliSimDigits object (fDigit) which
  // are stored in the digit tree.
  // The index map relates the AliSimDigits objects in the tree to dedicated pad rows
  // in the TPC
  // This loop filters the pad rows according to the slice no set via Init

  assert(fRowMax<fgkNRow);
  for(Int_t r=fRowMin;r<=fRowMax && r<fgkNRow && iResult>=0;r++){
    Int_t n=fIndex[fSlice][r];

    Int_t row=0;
    Int_t rowOffset=0;
    AliHLTTPCTransform::Slice2Sector(fSlice,r,sector,row);
    AliHLTTPCTransform::Slice2Sector(fSlice,AliHLTTPCTransform::GetFirstRow(fPatch),sector,rowOffset);

    if(n!=-1){//data on that row
      Int_t lslice,lrow;
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row);
      if(lrow!=r){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliDigits2Altro","Row")
	  <<AliHLTTPCLog::kDec<<"Rows on slice " << fSlice << " dont match "<<lrow<<" "<<r<<ENDLOG;
	continue;
      }

      Int_t channelAddress=-1;
      fDigits->First();
      do {
	time=fDigits->CurrentRow();
	pad=fDigits->CurrentColumn();
	dig = fDigits->GetDigit(time,pad);
	if (dig <= fParam->GetZeroSup()) continue;
	if(dig >= AliHLTTPCTransform::GetADCSat())
	  dig = AliHLTTPCTransform::GetADCSat();
	
	channelAddress=mapper.GetHwAddress(row-rowOffset, pad);
	iResult=encoder.AddChannelSignal(dig, time, channelAddress);
      } while (fDigits->Next() && iResult>=0);
      if (iResult>=0 && channelAddress>=0) {
	iResult=encoder.SetChannel(channelAddress);
      }
    }
  }

  if (iResult>=0) {
    iResult=(encoder.GetTotal40bitWords()*5)/4;
  }

  return iResult;
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
  if(eventmerge == kTRUE && event >= 1024)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFileHandler::AliAltroDigits2Memory","TrackIDs")
	<<"Too many events if you want to merge!"<<ENDLOG;
      return 0;
    }
  delete fDigits;
  fDigits=0;
  /* Dont understand why we have to do 
     reload the tree, but otherwise the code crashes */
  fDigitsTree=0;
  if(!GetDigitsTree(event)) return 0;

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

  for(Int_t r=fRowMin;r<=fRowMax && r<fgkNRow;r++){
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
	  localcount++;
	}
      }
    }
    Byte_t *tmp = (Byte_t*)tempPt;
    tmp += sizeof(AliHLTTPCDigitRowData)
      + ndigits[r]*sizeof(AliHLTTPCDigitData);
    tempPt = (AliHLTTPCDigitRowData*)tmp;
  }
  delete [] ndigits;
  return data;
}
 
Bool_t AliHLTTPCFileHandler::GetDigitsTree(Int_t event)
{
  //Connects to the TPC digit tree in the AliROOT file.
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCFileHandler::GetDigitsTree","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(event);
  tpcLoader->LoadDigits();
  fDigitsTree = tpcLoader->TreeD();
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

#if 0
// Is this really still used : JMT 2013-03-03
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

  TDirectory *savedir = gDirectory;
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
  carray.SetClusterType("AliTPCclusterMI");
  Bool_t clusterok = carray.ConnectTree(tpcLoader->TreeR());

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
    Int_t entriesInRow = clusterrow[i]->GetArray()->GetEntriesFast();
    for(Int_t j = 0;j<entriesInRow;j++){
      AliTPCclusterMI *c = (AliTPCclusterMI*)(*clusterrow[i])[j];
      data[n].fZ = c->GetZ();
      data[n].fY = c->GetY();
      data[n].fX = fParam->GetPadRowRadii(sector,row);
      data[n].fCharge = (UInt_t)c->GetQ();
      data[n].SetID( fSlice, pat, n );
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
#endif
