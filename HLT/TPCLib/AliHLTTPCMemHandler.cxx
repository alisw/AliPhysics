// @(#) $Id$
// Original: AliHLTMemHandler.cxx,v 1.52 2005/06/14 10:55:21 cvetan 

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

//  @file   AliHLTTPCMemHandler.cxx
//  @author U. Frankenfeld, A. Vestbo, C. Loizides, maintained by
//          Matthias Richter
//  @date   
//  @brief  input interface base class for the TPC tracking code before
//          migration to the HLT component framework

#include <cassert>
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCMemHandler.h"
#include "TMath.h"

#if __GNUC__ >= 3
using namespace std;
#endif
  
ClassImp(AliHLTTPCMemHandler)
  
AliHLTTPCMemHandler::AliHLTTPCMemHandler()
  :
  fRowMin(0),
  fRowMax(0),
  fSlice(0),
  fPatch(0),
  fInBinary(NULL),
  fOutBinary(NULL),
  fPt(NULL),
  fSize(0),
  fIsRandom(kFALSE),
  fNRandom(0),
  fNGenerate(0),
  fNUsed(0),
  fNDigits(0),
  fDPt(NULL),
  fRandomDigits(NULL),
  fDummy(0)
{ 
  //Constructor
  Init(0,0);
  ResetROI();
}

AliHLTTPCMemHandler::~AliHLTTPCMemHandler()
{
  //Destructor
  if(fPt) delete[] fPt;
  if(fRandomDigits) delete [] fRandomDigits;
  if(fDPt) delete [] fDPt;
}

void AliHLTTPCMemHandler::Init(Int_t s,Int_t p, Int_t *r)
{
  //init handler
  assert(s<fgkNSlice);
  if (s>fgkNSlice) {
    fSlice=0;
    fPatch=0;
    fRowMin=0;
    fRowMax=0;
    if (r) *r=0;
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Init","sector coordinates")
      <<"Invalid slice no " << s <<ENDLOG;
    return;
  }
  fSlice=s;fPatch=p;
  if(r) {
    fRowMin=r[0];
    fRowMax=r[1];
  }else{
     fRowMin=AliHLTTPCTransform::GetFirstRow(p);
     fRowMax=AliHLTTPCTransform::GetLastRow(p); 
  }
  ResetROI();
}

void AliHLTTPCMemHandler::ResetROI()
{
  //Resets the Look-up table for Region of Interest mode.
  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      fEtaMinTimeBin[i] = 0;
      fEtaMaxTimeBin[i] = AliHLTTPCTransform::GetNTimeBins()-1;
    }
}

void AliHLTTPCMemHandler::SetROI(const Float_t *eta,Int_t */*slice*/)
{
  // Init the Look-up table for the Region of Interest mode.
  //   Here you can specify a certain etaregion, - all data
  //   outside this region will be discarded:
  //   eta[0] = mimium eta
  //   eta[1] = maximum eta
  //   slice[0] = mimumum slice
  //   slice[1] = maximum slice


  if(TMath::Abs(eta[1])<.00001)
    {
      LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::SetROI","Eta Values")
	<<"Bad ROI parameters."<<ENDLOG;
      for(Int_t i=fRowMin; i<=fRowMax; i++)
	{
	  fEtaMinTimeBin[i]=0;
	  fEtaMaxTimeBin[i]=0;
	}
      return;
    }
  
  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      Int_t sector,row;
      Float_t xyz[3];
      
      Float_t thetamax = 2*atan(exp(-1.*eta[1]));
      
      xyz[0] = AliHLTTPCTransform::Row2X(i);
      xyz[1]=0;
      xyz[2] = xyz[0]/tan(thetamax);
      AliHLTTPCTransform::Slice2Sector(fSlice,i,sector,row);
      AliHLTTPCTransform::Local2Raw(xyz,sector,row);
      
      fEtaMinTimeBin[i] = (Int_t)xyz[2];
      
      if(TMath::Abs(eta[0])<.00001)
	fEtaMaxTimeBin[i] = 445;
      else
	{
	  Float_t thetamin = 2*atan(exp(-1.*eta[0]));
	  xyz[0] = AliHLTTPCTransform::Row2X(i);
	  xyz[1] = AliHLTTPCTransform::GetMaxY(i);
	  Float_t radii = sqrt(pow(xyz[0],2) + pow(xyz[1],2));
	  xyz[2] = radii/tan(thetamin);
	  AliHLTTPCTransform::Local2Raw(xyz,sector,row);
	  fEtaMaxTimeBin[i] = (Int_t)xyz[2];
	}
    }
  
}

Bool_t AliHLTTPCMemHandler::SetBinaryInput(char *name)
{
  //Set the input binary file.
  fInBinary = fopen(name,"r");
  if(!fInBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::SetBinaryInput","File Open")
      <<"Error opening file "<<name<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::SetBinaryInput(FILE *file)
{
  //Set the input binary file.
  fInBinary = file;
  if(!fInBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::SetBinaryInput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliHLTTPCMemHandler::CloseBinaryInput()
{
  //Close the input file.
  if(!fInBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::CloseBinaryInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fInBinary);
  fInBinary =0;
}

Bool_t AliHLTTPCMemHandler::SetBinaryOutput(char *name)
{
  //Set the binary output file.
    fOutBinary = fopen(name,"w");
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::SetBinaryOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::SetBinaryOutput(FILE *file)
{
  //Set the binary output file.
    fOutBinary = file;
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::SetBinaryOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliHLTTPCMemHandler::CloseBinaryOutput()
{
  //close binary  
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::CloseBinaryOutPut","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fOutBinary);
  fOutBinary =0;
}

UInt_t AliHLTTPCMemHandler::GetFileSize()
{
  //Returns the file size in bytes of the input file.
  if(!fInBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::GetFileSize","File")
      <<"No Input File"<<ENDLOG;
    return 0;
  }
  fseek(fInBinary,0,SEEK_END);
  long size=ftell(fInBinary);
  rewind(fInBinary);
  if (size<0) return 0;
  return (UInt_t)size;
}

Byte_t *AliHLTTPCMemHandler::Allocate()
{
  //Allocate
  return Allocate(GetFileSize()); 
}

Byte_t *AliHLTTPCMemHandler::Allocate(AliHLTTPCTrackArray *array)
{
  //Allocate memory for tracks in memory. Used by TrackArray2Binary()
  if(!array){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Allocate","Memory")
      <<"Pointer to AliHLTTPCTrackArray = 0x0 "<<ENDLOG;
    return 0;
  }
  return Allocate(array->GetOutSize()); 
}

Byte_t *AliHLTTPCMemHandler::Allocate(UInt_t size)
{
  //Allocate memory of size in bytes.
  if(fPt){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Allocate","Memory")
      <<"Delete Memory"<<ENDLOG;
    Free();
  } 
  fPt = new Byte_t[size];
  fSize = size;
  memset(fPt,0,fSize);
  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCMemHandler::Allocate","Memory")
  <<AliHLTTPCLog::kDec<<"Allocate "<<size<<" Bytes of Memory"<<ENDLOG;
  return fPt;
}

void AliHLTTPCMemHandler::Free()
{
  //Clear the memory, if allocated.
  if(!fPt){
    //    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCMemHandler::Free","Memory")
    //      <<"No Memory allocated - can't Free"<<ENDLOG;
    return;
  }  
  delete[] fPt;
  fPt = 0;
  fSize =0;
}

///////////////////////////////////////// Random
void AliHLTTPCMemHandler::SetRandomSeed()
{
  //If you are adding random data to the original data.
  time_t *tp=0;
  SetRandomSeed(time(tp));
}

void AliHLTTPCMemHandler::SetRandomCluster(Int_t maxnumber)
{
  //If you are adding random data to the original data.
  
  fIsRandom = kTRUE;
  fNRandom = maxnumber;
  fNDigits = 0;
  if(fRandomDigits) delete [] fRandomDigits;
  fRandomDigits = new AliHLTTPCRandomDigitData[fNRandom*9];
  if(fDPt) delete [] fDPt;
  fDPt = new AliHLTTPCRandomDigitData *[fNRandom*9];
}

void AliHLTTPCMemHandler::QSort(AliHLTTPCRandomDigitData **a, Int_t first, Int_t last)
{

   // Sort array of AliHLTTPCRandomDigitData pointers using a quicksort algorithm.
   // Uses CompareDigits() to compare objects.
   // Thanks to Root!

   static AliHLTTPCRandomDigitData *tmp;
   static int i;           // "static" to save stack space
   int j;

   while (last - first > 1) {
      i = first;
      j = last;
      for (;;) {
         while (++i < last && CompareDigits(a[i], a[first]) < 0)
            ;
         while (--j > first && CompareDigits(a[j], a[first]) > 0)
            ;
         if (i >= j)
            break;

         tmp  = a[i];
         a[i] = a[j];
         a[j] = tmp;
      }
      if (j == first) {
         ++first;
         continue;
      }
      tmp = a[first];
      a[first] = a[j];
      a[j] = tmp;

      if (j - first < last - (j + 1)) {
	QSort(a, first, j);
	first = j + 1;   // QSort(j + 1, last);
      } else {
	QSort(a, j + 1, last);
	last = j;        // QSort(first, j);
      }
   }
}

UInt_t AliHLTTPCMemHandler::GetRandomSize() const
{
  //get random size
  Int_t nrandom = 0;
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t npad=AliHLTTPCTransform::GetNPads(r);
    nrandom  += Int_t (fNGenerate * ((Double_t) npad/141.));
  }
  return 9 * nrandom * sizeof(AliHLTTPCDigitData);
}

void AliHLTTPCMemHandler::DigitizePoint(Int_t row, Int_t pad, 
				    Int_t time,Int_t charge)
{
  //Making one single random cluster.
  for(Int_t j=-1;j<2;j++){
    for(Int_t k=-1;k<2;k++){
      Int_t dcharge = charge;
      if(j) dcharge /=2;
      if(k) dcharge /=2;
      if(dcharge<10) continue;
      Int_t dpad  = j + pad;
      Int_t dtime = k + time;
      
      if(dpad<0||dpad>=AliHLTTPCTransform::GetNPads(row))  continue;
      if(dtime<0||dtime>=AliHLTTPCTransform::GetNTimeBins()) continue;
      
      fRandomDigits[fNDigits].fCharge = dcharge;
      fRandomDigits[fNDigits].fRow = row;
      fRandomDigits[fNDigits].fPad = dpad;
      fRandomDigits[fNDigits].fTime = dtime;
      fDPt[fNDigits] = &fRandomDigits[fNDigits];
      fNDigits++;
    }
  }
}

///////////////////////////////////////// Digit IO  
Bool_t AliHLTTPCMemHandler::Memory2BinaryFile(UInt_t nrow,AliHLTTPCDigitRowData *data)
{
  //Write data to the outputfile as is. No run-length encoding is done.

  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2Binary","File")
      <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2Binary","Memory")
      <<"Pointer to AliHLTTPCDigitRowData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  
  AliHLTTPCDigitRowData *rowPt = data; 
  Int_t outsize = 0;
  for(UInt_t i=0;i<nrow;i++){
    Int_t size = sizeof(AliHLTTPCDigitData) * rowPt->fNDigit 
      + sizeof(AliHLTTPCDigitRowData);
    outsize += size;
    fwrite(rowPt,size,1,fOutBinary);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    rowPt = (AliHLTTPCDigitRowData *) bytePt;
  }
  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCMemHandler::Memory2Binary","Memory")
    <<AliHLTTPCLog::kDec<<"Wrote "<<outsize<<" Bytes to Memory ("
    <<nrow<<" Rows)"<<ENDLOG;
  return kTRUE;
}

void AliHLTTPCMemHandler::AddData(AliHLTTPCDigitData *data,UInt_t & ndata,
			      UInt_t /*row*/,UShort_t pad,UShort_t time,UShort_t charge) const
{
  //add some data
  data[ndata].fPad = pad;
  data[ndata].fTime = time;
  data[ndata].fCharge = charge;
  ndata++;
}

void AliHLTTPCMemHandler::AddRandom(AliHLTTPCDigitData *data, UInt_t & ndata)
{
  //add some random data
  data[ndata].fPad = fDPt[fNUsed]->fPad;
  data[ndata].fTime = fDPt[fNUsed]->fTime;
  data[ndata].fCharge = fDPt[fNUsed]->fCharge;
  ndata++;
  fNUsed++;
}

void AliHLTTPCMemHandler::MergeDataRandom(AliHLTTPCDigitData *data, UInt_t & ndata,
				      UInt_t row, UShort_t pad, UShort_t time, UShort_t charge)
{
  //merge random data
  data[ndata].fPad = pad;
  data[ndata].fTime = time;
  data[ndata].fCharge = charge;
  while(ComparePoints(row,pad,time)==0){
    Int_t ch = data[ndata].fCharge + fDPt[fNUsed]->fCharge;
    if(charge>=AliHLTTPCTransform::GetADCSat()) ch = AliHLTTPCTransform::GetADCSat();
    data[ndata].fCharge = ch;
    fNUsed++;
  }
  ndata++;
}

void AliHLTTPCMemHandler::AddDataRandom(AliHLTTPCDigitData *data, UInt_t & ndata,
                   UInt_t row, UShort_t pad, UShort_t time, UShort_t charge)
{
  //add data random
  Int_t action;
  while((action=ComparePoints(row,pad,time))==1){
    AddRandom(data,ndata);
  }
  if(action==0){
    MergeDataRandom(data,ndata,row,pad,time,charge);
  }
  if(action<0){
    AddData(data,ndata,row,pad,time,charge);
  }  
}

void AliHLTTPCMemHandler::Write(UInt_t *comp, UInt_t & index, 
			    UInt_t & subindex, UShort_t value) const
{
  //write compressed data
  UInt_t shift[3] = {0,10,20};
  if(subindex==0) comp[index] =0; //clean up memory
  comp[index] |= (value&0x03ff)<<shift[subindex];
  if(subindex == 2){
    subindex = 0;
    index++;
  }
  else subindex++;
}

UShort_t AliHLTTPCMemHandler::Read(UInt_t *comp, UInt_t & index, UInt_t & subindex) const
{ 
  //read compressed data
  UInt_t shift[3] = {0,10,20};
  UShort_t value = (comp[index]>>shift[subindex])&0x03ff;
  if(subindex == 2){
    subindex = 0;
    index++;
  }
  else subindex++;
  
  return value;
}

UShort_t AliHLTTPCMemHandler::Test(const UInt_t *comp, 
			       UInt_t index, UInt_t  subindex) const
{
  //supi dupi test
  UInt_t shift[3] = {0,10,20};
  return (comp[index]>>shift[subindex])&0x03ff;
}

Int_t AliHLTTPCMemHandler::Memory2CompMemory(UInt_t nrow,
					 AliHLTTPCDigitRowData *data,UInt_t *comp)
{
  //Performs run-length encoding on data stored in memory pointed to by data.
  //The compressed data is written to comp.
  if(!comp){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2CompMemory","Memory")
      <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2CompMemory","Memory")
      <<"Pointer to AliHLTTPCDigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  AliHLTTPCDigitRowData *rowPt = data;
  UInt_t index=0;
  UInt_t subindex=0;
  
  for(UInt_t i=0;i<nrow;i++){
    UShort_t value = rowPt->fRow;
    Write(comp,index,subindex,value);
    UShort_t maxpad=0; 
    UShort_t npad=0;
    Int_t ddd[1000];
    for(Int_t d=0;d<200;d++) ddd[d]=0;
    for(UInt_t dig=0;dig<rowPt->fNDigit;dig++){
      if(rowPt->fDigitData[dig].fPad <200){ 
        ddd[rowPt->fDigitData[dig].fPad]++;
      }
    }
    for(Int_t d=0;d<200;d++){ 
      if(ddd[d]){
        npad++;
        maxpad =d;
      }
    }
    Write(comp,index,subindex,npad);
    UInt_t digit=0;
    for(UShort_t pad=0;pad <= maxpad;pad++){
      if(digit>=rowPt->fNDigit || rowPt->fDigitData[digit].fPad !=  pad)
        continue;
      Write(comp,index,subindex,pad);
//    write zero if time != 0
      if(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        if(rowPt->fDigitData[digit].fTime>0){
          Write(comp,index,subindex,0);
          Write(comp,index,subindex,rowPt->fDigitData[digit].fTime);
        }
      }
      while(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        UShort_t charge = rowPt->fDigitData[digit].fCharge;
        if(charge>=1023){
          charge=1023;
        }
        Write(comp,index,subindex,charge);
        if(digit+1<rowPt->fNDigit&&rowPt->fDigitData[digit+1].fPad == pad){
          if(rowPt->fDigitData[digit].fTime +1 !=
                     rowPt->fDigitData[digit+1].fTime){
            Write(comp,index,subindex,0);
            UShort_t nzero = rowPt->fDigitData[digit+1].fTime - 
                             (rowPt->fDigitData[digit].fTime +1);
            Write(comp,index,subindex,nzero);
          }  
        }
        digit++;
      }
      Write(comp,index,subindex,0);
      Write(comp,index,subindex,0);
    }
    
    Int_t size = sizeof(AliHLTTPCDigitData) * rowPt->fNDigit+
                                            sizeof(AliHLTTPCDigitRowData);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    rowPt = (AliHLTTPCDigitRowData *) bytePt;
  }
  while(subindex)
    Write(comp,index,subindex,0);
  return index * sizeof(UInt_t);
}

UInt_t AliHLTTPCMemHandler::GetCompMemorySize(UInt_t nrow,AliHLTTPCDigitRowData *data) const
{
  //Return the size of RLE data, after compressing data.
  
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::GetCompMemorySize","Memory")
      <<"Pointer to AliHLTTPCDigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  AliHLTTPCDigitRowData *rowPt = data;
  UInt_t index=0;
  
  for(UInt_t i=0;i<nrow;i++){
    index++;
    UShort_t maxpad=0; 
    UShort_t npad=0;
    Int_t ddd[1000];
    for(Int_t d=0;d<200;d++) ddd[d]=0;
    for(UInt_t dig=0;dig<rowPt->fNDigit;dig++){
      if(rowPt->fDigitData[dig].fPad <200){ 
        ddd[rowPt->fDigitData[dig].fPad]++;
      }
    }
    for(Int_t d=0;d<200;d++){ 
      if(ddd[d]){
        npad++;
        maxpad =d;
      }
    }
    index++;
    UInt_t digit=0;
    for(UShort_t pad=0;pad <= maxpad;pad++){
      if(digit>=rowPt->fNDigit || rowPt->fDigitData[digit].fPad !=  pad)
        continue;
      index++;
      //    write zero if time != 0
      if(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        if(rowPt->fDigitData[digit].fTime>0){
          index++;
          index++;
        }
      }
      while(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        index++;
        if(digit+1<rowPt->fNDigit&&rowPt->fDigitData[digit+1].fPad == pad){
          if(rowPt->fDigitData[digit].fTime +1 !=
                     rowPt->fDigitData[digit+1].fTime){
            index++;
            index++;
          }  
        }
        digit++;
      }
      index++;
      index++;
    }

    Int_t size = sizeof(AliHLTTPCDigitData) * rowPt->fNDigit+
                                            sizeof(AliHLTTPCDigitRowData);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    rowPt = (AliHLTTPCDigitRowData *) bytePt;
  }
  while(index%3)
    index++;
  return (index/3) * sizeof(UInt_t);
}

UInt_t AliHLTTPCMemHandler::GetMemorySize(UInt_t nrow,UInt_t *comp) const
{
  //get memory size
  if(!comp){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::GetMemorySize","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  Int_t outsize=0;

  UInt_t index=0;
  UInt_t subindex=0;

  for(UInt_t i=0;i<nrow;i++){
    UInt_t ndigit=0;
    Read(comp,index,subindex);
    UShort_t npad = Read(comp,index,subindex);
    for(UShort_t p=0;p<npad;p++){
      Read(comp,index,subindex);
      if(Test(comp,index,subindex)==0){
        Read(comp,index,subindex);
        if(Read(comp,index,subindex)== 0) continue;
      }
      for(;;){
        while(Read(comp,index,subindex)!=0) ndigit++;
        if(Read(comp,index,subindex)==0) break;
      }
    }
    Int_t size = sizeof(AliHLTTPCDigitData) * ndigit+
                                        sizeof(AliHLTTPCDigitRowData);
    outsize += size;
  }
   
  return outsize;
}

UInt_t AliHLTTPCMemHandler::GetNRow(UInt_t *comp,UInt_t size)
{
  //get number of rows
  if(!comp){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::GetNRow","Memory")
      <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  size = size /4;
  UInt_t nrow=0;
  UInt_t index=0;
  UInt_t subindex=0;
  while(index<size-1){ //don't start with last word
    nrow++;
    UInt_t ndigit=0;
    Read(comp,index,subindex);
    UShort_t npad = Read(comp,index,subindex);
    for(UShort_t p=0;p<npad;p++){
      Read(comp,index,subindex);
      if(Test(comp,index,subindex)==0){
        Read(comp,index,subindex);
        if(Read(comp,index,subindex)==0)continue;
      }
      for(;;){
        while(Read(comp,index,subindex)!=0) ndigit++;
        if(Read(comp,index,subindex)==0) break;
      }
    }
  }
  if(index==size-1){  //last word
    if(subindex<2){
      if(Read(comp,index,subindex)!=0) nrow++;
    }
  }
  return nrow;
}

Bool_t AliHLTTPCMemHandler::CompMemory2CompBinary(UInt_t nrow,UInt_t *comp,
                                             UInt_t size)
{
  //Write the RLE data in comp to the output file.
  
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::CompMemory2CompBinary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!comp){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::CompMemory2CompBinary","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(size==0)
    size=GetMemorySize(nrow,comp);
  if(!size){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::CompMemory2CompBinary","Memory")
    <<"Memory size = 0 "<<ENDLOG;
    return kFALSE;
  }
  UInt_t length = size/sizeof(UInt_t);
  fwrite(&length,sizeof(UInt_t),1,fOutBinary);  
  fwrite(comp,size,1,fOutBinary);
  return kTRUE;
}


Bool_t AliHLTTPCMemHandler::Memory2CompBinary(UInt_t nrow,AliHLTTPCDigitRowData *data)
{
  //Perform RLE on the data, and write it to the output file.
  Bool_t out = kTRUE;
  AliHLTTPCMemHandler * handler = new AliHLTTPCMemHandler();
  UInt_t size = GetCompMemorySize(nrow,data);
  UInt_t *comp =(UInt_t *)handler->Allocate(size);
  Memory2CompMemory(nrow,data,comp);
  CompMemory2CompBinary(nrow,comp,size);
  handler->Free();
  delete handler;
  return out;
}


///////////////////////////////////////// Point IO  
Bool_t AliHLTTPCMemHandler::Memory2Binary(UInt_t npoint,AliHLTTPCSpacePointData *data)
{
  //Writing spacepoints stored in data to the outputfile.
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2Binary","File")
      <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2Binary","Memory")
      <<"Pointer to AliHLTTPCSpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  UInt_t size = npoint*sizeof(AliHLTTPCSpacePointData);
  fwrite(data,size,1,fOutBinary);
  
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::Transform(UInt_t npoint,AliHLTTPCSpacePointData *data,Int_t slice)
{
  //Transform the space points in data, to global coordinates in slice.
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Transform","Memory")
    <<"Pointer to AliHLTTPCSpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  
  for(UInt_t i=0;i<npoint;i++){
    Float_t xyz[3];
    xyz[0] = data[i].fX;
    xyz[1] = data[i].fY;
    xyz[2] = data[i].fZ;
    AliHLTTPCTransform::Local2Global(xyz,slice);
    data[i].fX = xyz[0];
    data[i].fY = xyz[1];
    data[i].fZ = xyz[2];
  }
  return kTRUE;
}

///////////////////////////////////////// Track IO  
Bool_t AliHLTTPCMemHandler::Memory2Binary(UInt_t ntrack,AliHLTTPCTrackSegmentData *data)
{
  //Write the tracks stored in data, to outputfile.
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2Binary","Memory")
    <<"Pointer to AliHLTTPCTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliHLTTPCTrackSegmentData *trackPt = data;
  for(UInt_t i=0;i<ntrack;i++){
    Int_t size=sizeof(AliHLTTPCTrackSegmentData)+trackPt->fNPoints*sizeof(UInt_t); 
    fwrite(trackPt,size,1,fOutBinary);
    Byte_t *bytePt = (Byte_t*) trackPt;
    bytePt += size; 
    trackPt = (AliHLTTPCTrackSegmentData*) bytePt;
  }
  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCMemHandler::Memory2Binary","File")
  <<AliHLTTPCLog::kDec<<"Wrote  "<<ntrack<<" Tracks to File"<<ENDLOG;
  
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::TrackArray2Binary(AliHLTTPCTrackArray *array)
{
  //Write the trackarray to the outputfile.
  if(!fOutBinary){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::TrackArray2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::TrackArray2Binary","Memory")
    <<"Pointer to AliHLTTPCTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliHLTTPCTrackSegmentData *data = (AliHLTTPCTrackSegmentData *)Allocate(array);

  UInt_t ntrack;
  TrackArray2Memory(ntrack,data,array);
  Memory2Binary(ntrack,data);
  Free();
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::TrackArray2Memory(UInt_t & ntrack,AliHLTTPCTrackSegmentData *data,AliHLTTPCTrackArray *array) const
{
  //Fill the trackarray into the AliTrackSegmentData structures before writing to outputfile.
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::TrackArray2Memory","Memory")
    <<"Pointer to AliHLTTPCTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::TrackArray2Memory","Memory")
    <<"Pointer to AliHLTTPCTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  array->WriteTracks(ntrack,data);
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::Memory2TrackArray(UInt_t ntrack,AliHLTTPCTrackSegmentData *data,AliHLTTPCTrackArray *array) const
{
  //Fill the tracks in data into trackarray.
  
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTPCTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTPCTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->FillTracks(ntrack,data);
  return kTRUE;
}

Bool_t AliHLTTPCMemHandler::Memory2TrackArray(UInt_t ntrack,AliHLTTPCTrackSegmentData *data,AliHLTTPCTrackArray *array,Int_t slice) const
{
  //Fill the tracks in data into trackarray, and rotate the tracks to global coordinates.
    
  if(!data){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTPCTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTPCTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->FillTracks(ntrack,data,slice);
  return kTRUE;
}

void AliHLTTPCMemHandler::UpdateRowPointer(AliHLTTPCDigitRowData *&tempPt)
{
  //Update the data pointer to the next padrow in memory.
  
  Byte_t *tmp = (Byte_t*)tempPt;
  Int_t size = sizeof(AliHLTTPCDigitRowData) + tempPt->fNDigit*sizeof(AliHLTTPCDigitData);
  tmp += size;
  tempPt = (AliHLTTPCDigitRowData*)tmp;
}

Int_t  AliHLTTPCMemHandler::ComparePoints(UInt_t /*row*/,UShort_t pad,UShort_t time) const
{
  //compare two points
  if(fNUsed>=fNDigits) return -2;

  if(pad==fDPt[fNUsed]->fPad&&time==fDPt[fNUsed]->fTime) return 0;

  if(pad<fDPt[fNUsed]->fPad) return -1;
  if(pad==fDPt[fNUsed]->fPad&&time<fDPt[fNUsed]->fTime)  return -1;

  return 1;
}

Int_t AliHLTTPCMemHandler::CompareDigits(const AliHLTTPCRandomDigitData *a,const AliHLTTPCRandomDigitData *b) const
{
  //compare two digits
  if(a->fPad==b->fPad && a->fTime == b->fTime) return 0;

  if(a->fPad<b->fPad) return -1;
  if(a->fPad==b->fPad && a->fTime<b->fTime) return -1;
  
  return 1;
}
