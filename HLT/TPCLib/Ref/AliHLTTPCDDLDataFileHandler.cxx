// @(#) $Id$

// Author: C. Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDigitData.h"
#ifdef use_newio
#include "../RAW/AliRawReaderRoot.h"
#include "../RAW/AliRawReaderDate.h"
#else
#include "AliHLTTPCDDLTPCRawStream.h"
#include "AliHLTTPCDDLRawReaderFile.h"
#endif
#include "AliHLTTPCDDLDataFileHandler.h"

#if GCCVERSION == 3
using namespace std;
#endif

/** \class AliHLTTPCDDLDataFileHandler 
<pre>
//_____________________________________________________________
// AliHLTTPCDDLDataFileHandler
//
//  This class does converts from the DDL format of offline
//  into the memory I/O handling of the HLT binary files.
//  
//  Examples: see ddl2binary in exa and the general 
//            AliHLTTPCMemHandler class description
//
</pre>
*/

ClassImp(AliHLTTPCDDLDataFileHandler)

AliHLTTPCDDLDataFileHandler::AliHLTTPCDDLDataFileHandler()
{
  fReader=0;
  fTPCStream=0;
}

AliHLTTPCDDLDataFileHandler::~AliHLTTPCDDLDataFileHandler()
{
  FreeAll();
}

void AliHLTTPCDDLDataFileHandler::FreeAll()
{
  if(fReader) delete fReader;
  if(fTPCStream) delete fTPCStream;
  fReader = 0;
  fTPCStream = 0;
}


#ifdef use_newio
Bool_t AliHLTTPCDDLDataFileHandler::SetReaderInput(Char_t *name,Int_t event)
{
  fEvent=event;
  if(fReader) delete fReader;
  if(fTPCStream) delete fTPCStream;
  if(event>=0){
    fFilename=name;
    fReader=new AliRawReaderRoot(name,event);
  } else {
    fFilename="";
    fReader=new AliRawReaderDate(name);
  }
  fTPCStream=new AliTPCRawStream(fReader);

  return kTRUE;
}
#else
Bool_t AliHLTTPCDDLDataFileHandler::SetReaderInput(Char_t *name, Bool_t add)
{
  if(fReader){
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLDataFileHandler::SetReaderInput","File Open")
      <<"Reader ptr is already in use"<<ENDLOG;
    return kFALSE;
  }
 
  fReader=new AliHLTTPCDDLRawReaderFile(name,add);
  fTPCStream=new AliHLTTPCDDLTPCRawStream(fReader);
  
  return kTRUE;
}
Bool_t AliHLTTPCDDLDataFileHandler::SetReaderInput(AliHLTTPCDDLRawReaderFile *rf)
{
  if(fReader){
    LOG(AliHLTTPCLog::kError,"AliHLTTPCRawDataFileHandler::SetReaderInput","File Open")
      <<"Reader ptr is already in use, delete it first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with given file.
  fReader = rf;
  fTPCStream=new AliHLTTPCDDLTPCRawStream(fReader);

  return kTRUE;
}
#endif

void AliHLTTPCDDLDataFileHandler::CloseReaderInput()
{
  if(!fReader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCRawDataFileHandler::CloseReaderInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }

  delete fReader;
  delete fTPCStream;
  fReader = 0;
  fTPCStream = 0;
}

#ifdef use_newio
Bool_t AliHLTTPCDDLDataFileHandler::IsDigit(Int_t i)
{
  return kTRUE;
}
#endif

AliHLTTPCDigitRowData * AliHLTTPCDDLDataFileHandler::DDLData2Memory(UInt_t &nrow,Int_t event)
{
#ifdef use_newio
  if((fEvent>=0)&&(event!=fEvent)){
    fEvent=event;
    if(fReader) delete fReader;
    if(fTPCStream) delete fTPCStream;
    fReader=new AliRawReaderRoot(fFilename,event);
    fTPCStream=new AliTPCRawStream(fReader);
  }
#endif
  AliHLTTPCDigitRowData *data = 0;
  nrow=0;

  if(!fReader){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCDDLDataFileHandler::DDLData2Memory","File")
    <<"No Input avalible: no object AliHLTTPCDDLRawReaderFile"<<ENDLOG;
    return 0; 
  }
  
  Int_t nrows=fRowMax-fRowMin+1;
  Int_t ndigitcount=0;
  Int_t ndigits[nrows];
  UShort_t ***charges=new UShort_t**[nrows];
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t lrow=r-fRowMin;
    charges[lrow]=new UShort_t*[AliHLTTPCTransform::GetNPads(r)];
    for(Int_t k=0;k<AliHLTTPCTransform::GetNPads(r);k++){
      charges[lrow][k]=new UShort_t[AliHLTTPCTransform::GetNTimeBins()];
      for(Int_t j=0;j<AliHLTTPCTransform::GetNTimeBins();j++) charges[lrow][k][j]=0;
    }
  }
  
  Int_t ddls_to_search=0;
  Int_t ddls[9]={-1,-1,-1,-1,-1,-1,-1,-1,-1};
  Int_t ddlid=-1,lddlid=-1;
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    ndigits[r-fRowMin] = 0; //now digits on row

    Int_t patch=AliHLTTPCTransform::GetPatch(r);
    Int_t sector,row;
    AliHLTTPCTransform::Slice2Sector(fSlice,r,sector,row);

    if(sector<36) //taken from AliTPCBuffer160.cxx
      ddlid=sector*2+patch;
    else
      ddlid=70+(sector-36)*4+patch;

    if((lddlid!=ddlid-1)&&(r==30)){ //dont forget the split row on the last ddl
      ddls[ddls_to_search++]=ddlid-1;	
      lddlid=ddlid-1;
    }
    if((lddlid==-1)||(ddlid!=lddlid)){
      ddls[ddls_to_search++]=ddlid;
      lddlid=ddlid;
    }
    if((r==90)||(r==139)){ //dont forget the split row on the next ddl
      ddls[ddls_to_search++]=ddlid+1;	
      lddlid=ddlid+1;
    }
  }
  //for(Int_t i=0;i<ddls_to_search;i++) cout << ddls[i] <<endl;

  for(Int_t i=0;i<ddls_to_search;i++){
#ifdef use_newio
    fReader->Reset();
    fReader->Select(0,ddls[i],ddls[i]+1);
#else
    fTPCStream->SetDDLID(ddls[i]); //ddl to read out
#endif

    while (fTPCStream->Next()){
  
      UShort_t dig=fTPCStream->GetSignal();
      if(dig <= AliHLTTPCTransform::GetZeroSup()) continue;
      if(dig >= AliHLTTPCTransform::GetADCSat())
	dig = AliHLTTPCTransform::GetADCSat();

      Int_t time=fTPCStream->GetTime();
      Int_t pad=fTPCStream->GetPad();
      Int_t sector=fTPCStream->GetSector();
      Int_t row=fTPCStream->GetRow();
      Int_t slice,srow;

      //test row criteria (patch boundaries)
      AliHLTTPCTransform::Sector2Slice(slice,srow,sector,row);
      if((srow<fRowMin)||(srow>fRowMax))continue;
      if(slice!=fSlice){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLDataFileHandler::DDLDigits2Memory","Slice")
	  <<AliHLTTPCLog::kDec<<"Found slice "<<slice<<", expected "<<fSlice<<ENDLOG;
	continue;
      }

      //cut out the inner cone
      Float_t xyz[3];
      AliHLTTPCTransform::Raw2Local(xyz,sector,row,pad,time);
      if(AliHLTTPCTransform::Row2X(srow)<230./250.*fabs(xyz[2]))
	continue; // why 230???

      Int_t lrow=srow-fRowMin;
      if((lrow<0)||lrow>=nrows){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLDataFileHandler::DDLDigits2Memory","Row")
	  <<AliHLTTPCLog::kDec<<"Row value out of bounds "<<lrow<<" "<<nrows<<ENDLOG;
	continue;
      }
      if((pad<0)||(pad>=AliHLTTPCTransform::GetNPads(srow))){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLDataFileHandler::DDLDigits2Memory","Pad")
	  <<AliHLTTPCLog::kDec<<"Pad value out of bounds "<<pad<<" "
	  <<AliHLTTPCTransform::GetNPads(srow)<<ENDLOG;
	continue;
      }
      if((time<0)||(time>=AliHLTTPCTransform::GetNTimeBins())){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLDataFileHandler::DDLDigits2Memory","Time")
	  <<AliHLTTPCLog::kDec<<"Time out of bounds "<<time<<" "
	  <<AliHLTTPCTransform::GetNTimeBins()<<ENDLOG;
	continue;
      }

      //store digit
      ndigits[lrow]++; //for this row only
      ndigitcount++;   //total number of digits to be published

      charges[lrow][pad][time]=dig;
    }
  }
  
  Int_t size = sizeof(AliHLTTPCDigitData)*ndigitcount
    + nrows*sizeof(AliHLTTPCDigitRowData);

  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCDDLDataFileHandler::DDLDigits2Memory","Digits")
    <<AliHLTTPCLog::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliHLTTPCDigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliHLTTPCDigitRowData *tempPt = data;

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t lrow=r-fRowMin;
    tempPt->fRow = r;
    tempPt->fNDigit = ndigits[lrow];
  
    Int_t localcount=0;
    for(Int_t pad=0;pad<AliHLTTPCTransform::GetNPads(r);pad++){
      for(Int_t time=0;time<AliHLTTPCTransform::GetNTimeBins();time++){
	UShort_t dig=charges[lrow][pad][time];
	if(!dig) continue;

	if(localcount >= ndigits[lrow])
	  LOG(AliHLTTPCLog::kFatal,"AliHLTTPCDDLDataFileHandler::DDLDigits2Binary","Memory")
	    <<AliHLTTPCLog::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	    <<ndigits[lrow]<<ENDLOG;
	

	tempPt->fDigitData[localcount].fCharge=dig;
	tempPt->fDigitData[localcount].fPad=pad;
	tempPt->fDigitData[localcount].fTime=time;
#ifdef do_mc
	tempPt->fDigitData[localcount].fTrackID[0] = 0;
	tempPt->fDigitData[localcount].fTrackID[1] = 0;
	tempPt->fDigitData[localcount].fTrackID[2] = 0;
#endif
	localcount++;
      }
    }

    if(localcount != ndigits[lrow])
      LOG(AliHLTTPCLog::kFatal,"AliHLTTPCDDLDataFileHandler::DDLDigits2Binary","Memory")
	<<AliHLTTPCLog::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	<<ndigits[lrow]<<ENDLOG;


    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliHLTTPCDigitRowData)
                                      + ndigits[lrow]*sizeof(AliHLTTPCDigitData);
    tmp += size;
    tempPt = (AliHLTTPCDigitRowData*)tmp;
  }

  //delete charge array
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t lrow=r-fRowMin;
    for(Int_t k=0;k<AliHLTTPCTransform::GetNPads(r);k++)
	delete charges[lrow][k];
    delete charges[lrow];
  }
  delete charges;

  return data;
}


Bool_t AliHLTTPCDDLDataFileHandler::DDLData2CompBinary(Int_t event)
{
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliHLTTPCDigitRowData *digits=0;
  digits = DDLData2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}
