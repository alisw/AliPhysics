// @(#) $Id$

// Author: C. Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3RootTypes.h"
#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3DigitData.h"
#ifdef use_newio
#include "../RAW/AliRawReaderRoot.h"
#include "../RAW/AliRawReaderDate.h"
#else
#include "AliL3DDLTPCRawStream.h"
#include "AliL3DDLRawReaderFile.h"
#endif
#include "AliL3DDLDataFileHandler.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** \class AliL3DDLDataFileHandler 
<pre>
//_____________________________________________________________
// AliL3DDLDataFileHandler
//
//  This class does converts from the DDL format of offline
//  into the memory I/O handling of the HLT binary files.
//  
//  Examples: see ddl2binary in exa and the general 
//            AliL3MemHandler class description
//
</pre>
*/

ClassImp(AliL3DDLDataFileHandler)

AliL3DDLDataFileHandler::AliL3DDLDataFileHandler()
{
  fReader=0;
  fTPCStream=0;
}

AliL3DDLDataFileHandler::~AliL3DDLDataFileHandler()
{
  FreeAll();
}

void AliL3DDLDataFileHandler::FreeAll()
{
  if(fReader) delete fReader;
  if(fTPCStream) delete fTPCStream;
  fReader = 0;
  fTPCStream = 0;
}


#ifdef use_newio
Bool_t AliL3DDLDataFileHandler::SetReaderInput(Char_t *name,Int_t event)
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
Bool_t AliL3DDLDataFileHandler::SetReaderInput(Char_t *name, Bool_t add)
{
  if(fReader){
    LOG(AliL3Log::kError,"AliL3DDLDataFileHandler::SetReaderInput","File Open")
      <<"Reader ptr is already in use"<<ENDLOG;
    return kFALSE;
  }
 
  fReader=new AliL3DDLRawReaderFile(name,add);
  fTPCStream=new AliL3DDLTPCRawStream(fReader);
  
  return kTRUE;
}
Bool_t AliL3DDLDataFileHandler::SetReaderInput(AliL3DDLRawReaderFile *rf)
{
  if(fReader){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetReaderInput","File Open")
      <<"Reader ptr is already in use, delete it first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with given file.
  fReader = rf;
  fTPCStream=new AliL3DDLTPCRawStream(fReader);

  return kTRUE;
}
#endif

void AliL3DDLDataFileHandler::CloseReaderInput()
{
  if(!fReader){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::CloseReaderInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }

  delete fReader;
  delete fTPCStream;
  fReader = 0;
  fTPCStream = 0;
}

#ifdef use_newio
Bool_t AliL3DDLDataFileHandler::IsDigit(Int_t i)
{
  return kTRUE;
}
#endif

AliL3DigitRowData * AliL3DDLDataFileHandler::DDLData2Memory(UInt_t &nrow,Int_t event)
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
  AliL3DigitRowData *data = 0;
  nrow=0;

  if(!fReader){
    LOG(AliL3Log::kWarning,"AliL3DDLDataFileHandler::DDLData2Memory","File")
    <<"No Input avalible: no object AliL3DDLRawReaderFile"<<ENDLOG;
    return 0; 
  }
  
  Int_t nrows=fRowMax-fRowMin+1;
  Int_t ndigitcount=0;
  //  Int_t ndigits[nrows];
  Int_t * ndigits = new Int_t[nrows];
  UShort_t ***charges=new UShort_t**[nrows];
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t lrow=r-fRowMin;
    charges[lrow]=new UShort_t*[AliL3Transform::GetNPads(r)];
    for(Int_t k=0;k<AliL3Transform::GetNPads(r);k++){
      charges[lrow][k]=new UShort_t[AliL3Transform::GetNTimeBins()];
      for(Int_t j=0;j<AliL3Transform::GetNTimeBins();j++) charges[lrow][k][j]=0;
    }
  }
  
  Int_t ddls_to_search=0;
  Int_t ddls[9]={-1,-1,-1,-1,-1,-1,-1,-1,-1};
  Int_t ddlid=-1,lddlid=-1;
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    ndigits[r-fRowMin] = 0; //now digits on row

    Int_t patch=AliL3Transform::GetPatch(r);
    Int_t sector,row;
    AliL3Transform::Slice2Sector(fSlice,r,sector,row);

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
      if(dig <= AliL3Transform::GetZeroSup()) continue;
      if(dig >= AliL3Transform::GetADCSat())
	dig = AliL3Transform::GetADCSat();

      Int_t time=fTPCStream->GetTime();
      Int_t pad=fTPCStream->GetPad();
      Int_t sector=fTPCStream->GetSector();
      Int_t row=fTPCStream->GetRow();
      Int_t slice,srow;

      //test row criteria (patch boundaries)
      AliL3Transform::Sector2Slice(slice,srow,sector,row);
      if((srow<fRowMin)||(srow>fRowMax))continue;
      if(slice!=fSlice){
	LOG(AliL3Log::kError,"AliL3DDLDataFileHandler::DDLDigits2Memory","Slice")
	  <<AliL3Log::kDec<<"Found slice "<<slice<<", expected "<<fSlice<<ENDLOG;
	continue;
      }

      //cut out the inner cone
      Float_t xyz[3];
      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
      if(AliL3Transform::Row2X(srow)<230./250.*fabs(xyz[2]))
	continue; // why 230???

      Int_t lrow=srow-fRowMin;
      if((lrow<0)||lrow>=nrows){
	LOG(AliL3Log::kError,"AliL3DDLDataFileHandler::DDLDigits2Memory","Row")
	  <<AliL3Log::kDec<<"Row value out of bounds "<<lrow<<" "<<nrows<<ENDLOG;
	continue;
      }
      if((pad<0)||(pad>=AliL3Transform::GetNPads(srow))){
	LOG(AliL3Log::kError,"AliL3DDLDataFileHandler::DDLDigits2Memory","Pad")
	  <<AliL3Log::kDec<<"Pad value out of bounds "<<pad<<" "
	  <<AliL3Transform::GetNPads(srow)<<ENDLOG;
	continue;
      }
      if((time<0)||(time>=AliL3Transform::GetNTimeBins())){
	LOG(AliL3Log::kError,"AliL3DDLDataFileHandler::DDLDigits2Memory","Time")
	  <<AliL3Log::kDec<<"Time out of bounds "<<time<<" "
	  <<AliL3Transform::GetNTimeBins()<<ENDLOG;
	continue;
      }

      //store digit
      ndigits[lrow]++; //for this row only
      ndigitcount++;   //total number of digits to be published

      charges[lrow][pad][time]=dig;
    }
  }
  
  Int_t size = sizeof(AliL3DigitData)*ndigitcount
    + nrows*sizeof(AliL3DigitRowData);

  LOG(AliL3Log::kDebug,"AliL3DDLDataFileHandler::DDLDigits2Memory","Digits")
    <<AliL3Log::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliL3DigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliL3DigitRowData *tempPt = data;

  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t lrow=r-fRowMin;
    tempPt->fRow = r;
    tempPt->fNDigit = ndigits[lrow];
  
    Int_t localcount=0;
    for(Int_t pad=0;pad<AliL3Transform::GetNPads(r);pad++){
      for(Int_t time=0;time<AliL3Transform::GetNTimeBins();time++){
	UShort_t dig=charges[lrow][pad][time];
	if(!dig) continue;

	if(localcount >= ndigits[lrow])
	  LOG(AliL3Log::kFatal,"AliL3DDLDataFileHandler::DDLDigits2Binary","Memory")
	    <<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
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
      LOG(AliL3Log::kFatal,"AliL3DDLDataFileHandler::DDLDigits2Binary","Memory")
	<<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	<<ndigits[lrow]<<ENDLOG;


    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliL3DigitRowData)
                                      + ndigits[lrow]*sizeof(AliL3DigitData);
    tmp += size;
    tempPt = (AliL3DigitRowData*)tmp;
  }

  //delete charge array
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t lrow=r-fRowMin;
    for(Int_t k=0;k<AliL3Transform::GetNPads(r);k++)
	delete charges[lrow][k];
    delete charges[lrow];
  }
  delete charges;
  delete [] ndigits;

  return data;
}


Bool_t AliL3DDLDataFileHandler::DDLData2CompBinary(Int_t event)
{
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliL3DigitRowData *digits=0;
  digits = DDLData2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}
