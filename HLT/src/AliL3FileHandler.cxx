// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>, Anders Vestbo <mailto:vestbo$fi.uib.no>, C. Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group 

#include "AliL3StandardIncludes.h"
#include <TClonesArray.h>
#include <TSystem.h>
#include <AliTPC.h>

#include <AliTPCParamSR.h>
#include <AliTPCDigitsArray.h>
#include <AliTPCClustersArray.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersRow.h>
#include <AliSimDigits.h>

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3FileHandler.h"
#include "AliL3DigitData.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3SpacePointData.h"
#include "AliL3TrackArray.h"

#if GCCVERSION == 3
using namespace std;
#endif

/** \class AliL3FileHandler
<pre>
//_____________________________________________________________
// AliL3FileHandler
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
// AliL3FileHandler file;
// file.Init(slice,patch);
// file.SetAliInput("galice.root");
// AliL3DigitRowData *dataPt = (AliL3DigitRowData*)file.AliDigits2Memory(nrows,eventnr);
// 
// All the data are then stored in memory and accessible via the pointer dataPt.
// Accesing the data is then identical to the example 1) showed in AliL3MemHandler class.
//
// For converting the data back, and writing it to a new AliROOT file do:
//
// AliL3FileHandler file;
// file.Init(slice,patch);
// file.SetAliInput("galice.root");
// file.Init(slice,patch,NumberOfRowsInPatch);
// file.AliDigits2RootFile(dataPt,"new_galice.root");
// file.CloseAliInput();
</pre>
*/

ClassImp(AliL3FileHandler)

AliL3FileHandler::AliL3FileHandler()
{
  //Default constructor
  fInAli = 0;
  fParam = 0;
  fMC =0;
  fLastIndex=0;
  fDigits=0;
  fDigitsTree=0;
  for(Int_t i=0;i<AliL3Transform::GetNSlice();i++){
    for(Int_t j=0;j<AliL3Transform::GetNRows();j++)
      fIndex[i][j]=-1;
  }
  fIndexCreated=kFALSE;
}

AliL3FileHandler::~AliL3FileHandler()
{
  //Destructor
  if(fMC) CloseMCOutput();
  FreeDigitsTree();
  if(fInAli) CloseAliInput();
}

void AliL3FileHandler::FreeDigitsTree()
{
  if(!fDigitsTree)
    {
      LOG(AliL3Log::kWarning,"AliL3FileHandler::FreeDigitsTree()","Pointer")
	<<"Cannot free digitstree, it is not present"<<ENDLOG;
      return;
    }
  fDigits=0;
  fDigitsTree->Delete();
  fDigitsTree=0;
  fLastIndex=0;
  for(Int_t i=0;i<AliL3Transform::GetNSlice();i++){
    for(Int_t j=0;j<AliL3Transform::GetNRows();j++)
      fIndex[i][j]=-1;
  }
  fIndexCreated=kFALSE;
}

Bool_t AliL3FileHandler::SetMCOutput(Char_t *name)
{
  fMC = fopen(name,"w");
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetMCOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliL3FileHandler::SetMCOutput(FILE *file)
{
  fMC = file;
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetMCOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliL3FileHandler::CloseMCOutput()
{
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::CloseMCOutPut","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fMC);
  fMC =0;
}

Bool_t AliL3FileHandler::SetAliInput()
{
  fInAli->CdGAFile();
  fParam = AliTPC::LoadTPCParam(gFile);
  if(!fParam){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File")
      <<"No TPC parameters found in \""<<gFile->GetName()
      <<"\", creating standard parameters "
      <<"which might not be what you want!"<<ENDLOG;
    fParam = new AliTPCParamSR;
  }
  if(!fParam){ 
    LOG(AliL3Log::kError,"AliL3FileHandler::SetAliInput","File Open")
      <<"No AliTPCParam "<<AliL3Transform::GetParamName()<<" in File "<<gFile->GetName()<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliL3FileHandler::SetAliInput(Char_t *name)
{
  //Open the AliROOT file with name.
  
  fInAli= AliRunLoader::Open(name);
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File Open")
    <<"Pointer to AliRunLoader = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

Bool_t AliL3FileHandler::SetAliInput(AliRunLoader *runLoader)
{
  //Specify already opened AliROOT file to use as an input.
  
  fInAli=runLoader;
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File Open")
    <<"Pointer to AliRunLoader = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

void AliL3FileHandler::CloseAliInput()
{
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::CloseAliInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  delete fInAli;
  fInAli = 0;
}

Bool_t AliL3FileHandler::IsDigit(Int_t event)
{
  //Check if there is a TPC digit tree in the current file.
  //Return kTRUE if tree was found, and kFALSE if not found.
  
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::IsDigit","File")
    <<"Pointer to AliRunLoader = 0x0 "<<ENDLOG;
    return kTRUE;  //maybe you are using binary input which is Digits!!  
  }
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::IsDigit","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(event);
  tpcLoader->LoadDigits();
  TTree *t=tpcLoader->TreeD();
  if(t){
    LOG(AliL3Log::kInformational,"AliL3FileHandler::IsDigit","File Type")
    <<"Found Digit Tree -> Use Fast Cluster Finder"<<ENDLOG;
    return kTRUE;
  }
  else{
    LOG(AliL3Log::kInformational,"AliL3FileHandler::IsDigit","File Type")
    <<"No Digit Tree -> Use Cluster Tree"<<ENDLOG;
    return kFALSE;
  }
}

///////////////////////////////////////// Digit IO  
Bool_t AliL3FileHandler::AliDigits2Binary(Int_t event,Bool_t altro)
{
  Bool_t out = kTRUE;
  UInt_t nrow;
  AliL3DigitRowData* data = 0;
  if(altro)
    data = AliAltroDigits2Memory(nrow,event);
  else
    data = AliDigits2Memory(nrow,event);
  out = Memory2Binary(nrow,data);
  Free();
  return out;
}

Bool_t AliL3FileHandler::AliDigits2CompBinary(Int_t event,Bool_t altro)
{
  //Convert AliROOT TPC data, into HLT data format.
  //event specifies the event you want in the aliroot file.
  
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliL3DigitRowData *digits=0;
  if(altro)
    digits = AliAltroDigits2Memory(ndigits,event);
  else
    digits = AliDigits2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}

Bool_t AliL3FileHandler::CreateIndex()
{
  //create the access index

  LOG(AliL3Log::kInformational,"AliL3FileHandler::CreateIndex","Index")
    <<"Starting to create index, this can take a while."<<ENDLOG;

  fIndexCreated=kFALSE;
  for(Int_t i=0;i<AliL3Transform::GetNSlice();i++){
    for(Int_t j=0;j<AliL3Transform::GetNRows();j++)
      fIndex[i][j]=-1;
  }

  for(Int_t n=0; n<fDigitsTree->GetEntries(); n++){
    Int_t sector, row;
    Int_t lslice,lrow;
    fDigitsTree->GetEvent(n);
    fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
    if(!AliL3Transform::Sector2Slice(lslice,lrow,sector,row)){
      LOG(AliL3Log::kError,"AliL3FileHandler::CreateIndex","Slice/Row")
	<<AliL3Log::kDec<<"Index could not be created. Wrong values "
	<<sector<<" "<<row<<ENDLOG;
      return kFALSE;
    }
    //cout << lslice << " " << lrow << " " << sector << " " << row << endl;
    if(fIndex[lslice][lrow]==-1) fIndex[lslice][lrow]=n;
  }

  //for(Int_t i=0;i<AliL3Transform::GetNSlice();i++) cout << i << " " << fIndex[i][0] << endl;

  LOG(AliL3Log::kInformational,"AliL3FileHandler::CreateIndex","Index")
    <<"Index successfully created."<<ENDLOG;

  fIndexCreated=kTRUE;
  return fIndexCreated;
}

AliL3DigitRowData * AliL3FileHandler::AliDigits2Memory(UInt_t & nrow,Int_t event)
{
  //Read data from AliROOT file into memory, and store it in the HLT data format.
  //Returns a pointer to the data.

  AliL3DigitRowData *data = 0;
  nrow=0;
  
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliDigits2Memory","File")
    <<"No Input avalible: no object AliRunLoader"<<ENDLOG;
    return 0; 
  }
  
  if(!fDigitsTree)
    if(!GetDigitsTree(event)) return 0;
  
  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t lslice,lrow;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)fDigitsTree->GetEntries();
  Int_t ndigits[entries];
  Float_t xyz[3];

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];
    if(n==-1) continue; //no data on that row
    
    fDigitsTree->GetEvent(n);
    fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);

    if(lrow!=r){
      LOG(AliL3Log::kError,"AliL3FileHandler::AliDigits2Memory","Row")
	<<AliL3Log::kDec<<"Rows dont match "<<lrow<<" "<<r<<ENDLOG;
    }

    ndigits[lrow] = 0;
    fDigits->First();
    do {
      time=fDigits->CurrentRow();
      pad=fDigits->CurrentColumn();
      dig = fDigits->GetDigit(time,pad);
      if(dig <= fParam->GetZeroSup()) continue;
      if(dig >= AliL3Transform::GetADCSat())
	dig = AliL3Transform::GetADCSat();
      
      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
      if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	continue; // why 230???

      ndigits[lrow]++; //for this row only
      ndigitcount++;   //total number of digits to be published

    } while (fDigits->Next());
    //cout << lrow << " " << ndigits[lrow] << " - " << ndigitcount << endl;
    nrows++;
  }

  Int_t size = sizeof(AliL3DigitData)*ndigitcount
    + nrows*sizeof(AliL3DigitRowData);

  LOG(AliL3Log::kDebug,"AliL3FileHandler::AliDigits2Memory","Digits")
    <<AliL3Log::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliL3DigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliL3DigitRowData *tempPt = data;

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];

    if(n==-1) continue; //no data on that row
    
    fDigitsTree->GetEvent(n);
    fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
    if(lrow!=r){
      LOG(AliL3Log::kError,"AliL3FileHandler::AliDigits2Memory","Row")
	<<AliL3Log::kDec<<"Rows dont match "<<lrow<<" "<<r<<ENDLOG;
    }
    tempPt->fRow = lrow;
    tempPt->fNDigit = ndigits[lrow];

    Int_t localcount=0;
    fDigits->First();
    do {
      time=fDigits->CurrentRow();
      pad=fDigits->CurrentColumn();
      dig = fDigits->GetDigit(time,pad);
      if (dig <= fParam->GetZeroSup()) continue;
      if(dig >= AliL3Transform::GetADCSat())
	dig = AliL3Transform::GetADCSat();

      //Exclude data outside cone:
      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
      if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	continue; // why 230???

      if(localcount >= ndigits[lrow])
	LOG(AliL3Log::kFatal,"AliL3FileHandler::AliDigits2Binary","Memory")
	  <<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
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

    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliL3DigitRowData)
                                      + ndigits[lrow]*sizeof(AliL3DigitData);
    tmp += size;
    tempPt = (AliL3DigitRowData*)tmp;
  }
  return data;
}

AliL3DigitRowData * AliL3FileHandler::AliAltroDigits2Memory(UInt_t & nrow,Int_t event,Bool_t eventmerge)
{
  //Read data from AliROOT file into memory, and store it in the HLT data format.
  //Returns a pointer to the data.
  //This functions filter out single timebins, which is noise. The timebins which
  //are removed are timebins which have the 4 zero neighbours; 
  //(pad-1,time),(pad+1,time),(pad,time-1),(pad,time+1).
  
  AliL3DigitRowData *data = 0;
  nrow=0;
  
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliAltroDigits2Memory","File")
    <<"No Input avalible: no object AliRunLoader"<<ENDLOG;
    return 0; 
  }
  if(eventmerge == kTRUE && event >= 1024)
    {
      LOG(AliL3Log::kError,"AliL3FileHandler::AliAltroDigits2Memory","TrackIDs")
	<<"Too many events if you want to merge!!"<<ENDLOG;
      return 0;
    }
  
  if(!fDigitsTree)
    if(!GetDigitsTree(event)) return 0;
  
  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)fDigitsTree->GetEntries();
  Int_t ndigits[entries];
  Int_t lslice,lrow;
  Int_t zerosupval=AliL3Transform::GetZeroSup();
  Float_t xyz[3];
  
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];
    if(n==-1) continue; //no data on that row
    
    fDigitsTree->GetEvent(n);
    fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);

    if(lrow!=r){
      LOG(AliL3Log::kError,"AliL3FileHandler::AliAltroDigits2Memory","Row")
	<<AliL3Log::kDec<<"Rows dont match "<<lrow<<" "<<r<<ENDLOG;
    }

    ndigits[lrow] = 0;
    fDigits->ExpandBuffer();
    fDigits->ExpandTrackBuffer();
    for(Int_t i=0; i<fDigits->GetNCols(); i++){
      for(Int_t j=0; j<fDigits->GetNRows(); j++){
	pad=i;
	time=j;
	dig = fDigits->GetDigitFast(time,pad);
	if(dig <= zerosupval) continue;
	if(dig >= AliL3Transform::GetADCSat())
	  dig = AliL3Transform::GetADCSat();

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

	AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	  continue; 
	      
	ndigits[lrow]++; //for this row only
	ndigitcount++;   //total number of digits to be published
      }
    }
    
    nrows++;
  }
  
  Int_t size = sizeof(AliL3DigitData)*ndigitcount
    + nrows*sizeof(AliL3DigitRowData);

  LOG(AliL3Log::kDebug,"AliL3FileHandler::AliAltroDigits2Memory","Digits")
    <<AliL3Log::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliL3DigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliL3DigitRowData *tempPt = data;
 
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t n=fIndex[fSlice][r];
    if(n==-1) continue; //no data on that row
    
    fDigitsTree->GetEvent(n);
    fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);

    if(lrow!=r){
      LOG(AliL3Log::kError,"AliL3FileHandler::AliAltroDigits2Memory","Row")
	<<AliL3Log::kDec<<"Rows dont match "<<lrow<<" "<<r<<ENDLOG;
    }

    tempPt->fRow = lrow;
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
	if(dig >= AliL3Transform::GetADCSat())
	  dig = AliL3Transform::GetADCSat();
	      
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
	
	AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
	  continue;
	    
	if(localcount >= ndigits[lrow])
	  LOG(AliL3Log::kFatal,"AliL3FileHandler::AliAltroDigits2Binary","Memory")
	    <<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
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
      
    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliL3DigitRowData)
      + ndigits[lrow]*sizeof(AliL3DigitData);
    tmp += size;
    tempPt = (AliL3DigitRowData*)tmp;
  }
  return data;
}
 
Bool_t AliL3FileHandler::GetDigitsTree(Int_t event)
{
  //Connects to the TPC digit tree in the AliROOT file.
  
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::GetDigitsTree","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(event);
  tpcLoader->LoadDigits();
  fDigitsTree = tpcLoader->TreeD();
  if(!fDigitsTree) 
    {
      LOG(AliL3Log::kError,"AliL3FileHandler::GetDigitsTree","Digits Tree")
	<<AliL3Log::kHex<<"Error getting digitstree "<<(Int_t)fDigitsTree<<ENDLOG;
      return kFALSE;
    }
  fDigitsTree->GetBranch("Segment")->SetAddress(&fDigits);

  return CreateIndex();
}

void AliL3FileHandler::AliDigits2RootFile(AliL3DigitRowData *rowPt,Char_t *new_digitsfile)
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
      printf("AliL3FileHandler::AliDigits2RootFile : No rootfile\n");
      return;
    }
  if(!fParam)
    {
      printf("AliL3FileHandler::AliDigits2RootFile : No parameter object. Run on rootfile\n");
      return;
    }
  
  //Get the original digitstree:
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliDigits2RootFile","File")
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
      printf("AliL3FileHandler::AliDigits2RootFile : No digits tree object\n");
      return;
    }

  tpcLoader->SetDigitsFileName(new_digitsfile);
  tpcLoader->MakeDigitsContainer();
    
  //setup a new one, or connect it to the existing one:
  AliTPCDigitsArray *arr = new AliTPCDigitsArray(); 
  arr->SetClass("AliSimDigits");
  arr->Setup(fParam);
  arr->MakeTree(tpcLoader->TreeD());

  Int_t digcounter=0,trackID[3];

  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      
      if((Int_t)rowPt->fRow != i) 
	cerr<<"AliL3FileHandler::AliDigits2RootFile : Mismatching row numbering!!!"<<endl;
            
      Int_t sector,row;
      AliL3Transform::Slice2Sector(fSlice,i,sector,row);
      
      AliSimDigits *old_dig = (AliSimDigits*)old_array->LoadRow(sector,row);
      AliSimDigits * dig = (AliSimDigits*)arr->CreateRow(sector,row);
      old_dig->ExpandBuffer();
      old_dig->ExpandTrackBuffer();
      dig->ExpandBuffer();
      dig->ExpandTrackBuffer();
      
      if(!old_dig)
      	printf("AliL3FileHandler::AliDigits2RootFile : No padrow %d %d\n",sector,row);
      
      AliL3DigitData *digPt = rowPt->fDigitData;
      digcounter=0;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  Short_t charge = (Short_t)digPt[j].fCharge;
	  Int_t pad = (Int_t)digPt[j].fPad;
	  Int_t time = (Int_t)digPt[j].fTime;
	  
	  if(charge == 0) //Only write the digits that has not been removed
	    {
	      cerr<<"AliL3FileHandler::AliDigits2RootFile : Zero charge!!! "<<endl;
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

  arr->GetTree()->SetName(treeName);
  arr->GetTree()->AutoSave();
  tpcLoader->WriteDigits("OVERWRITE");
  delete arr;
  delete old_array;
  
}

///////////////////////////////////////// Point IO  
Bool_t AliL3FileHandler::AliPoints2Binary(Int_t eventn){
  Bool_t out = kTRUE;
  UInt_t npoint;
  AliL3SpacePointData *data = AliPoints2Memory(npoint,eventn);
  out = Memory2Binary(npoint,data);
  Free();
  return out;
}

AliL3SpacePointData * AliL3FileHandler::AliPoints2Memory(UInt_t & npoint,Int_t eventn){
  AliL3SpacePointData *data = 0;
  npoint=0;
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliPoints2Memory","File")
    <<"No Input avalible: no object AliRunLoader"<<ENDLOG;
    return 0;
  }

  TDirectory *savedir = gDirectory;
  AliLoader* tpcLoader = fInAli->GetLoader("TPCLoader");
  if(!tpcLoader){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliPoints2Memory","File")
    <<"Pointer to AliLoader for TPC = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  fInAli->GetEvent(eventn);
  tpcLoader->LoadRecPoints();

  AliTPCClustersArray carray;
  carray.Setup(fParam);
  carray.SetClusterType("AliTPCcluster");
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
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
    if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
    clusterrow[i] = carray.GetRow(sector,row);
    if(clusterrow[i])
      sum+=clusterrow[i]->GetArray()->GetEntriesFast();
  }
  UInt_t size = sum*sizeof(AliL3SpacePointData);

  LOG(AliL3Log::kDebug,"AliL3FileHandler::AliPoints2Memory","File")
  <<AliL3Log::kDec<<"Found "<<sum<<" SpacePoints"<<ENDLOG;

  data = (AliL3SpacePointData *) Allocate(size);
  npoint = sum;
  UInt_t n=0; 
  Int_t pat=fPatch;
  if(fPatch==-1)
    pat=0;
  for(Int_t i=0; i<carray.GetTree()->GetEntries(); i++){
    if(!clusterrow[i]) continue;
    Int_t row = rows[i];
    Int_t sector = sects[i];
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
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

