//Author:        Uli Frankenfeld
//Last Modified: 17.12.2000

#include <math.h>
#include <time.h>
#include <iostream.h>
#include <stdio.h>

#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3Logging.h"

#include "AliL3DigitData.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3SpacePointData.h"
#include "AliL3TrackArray.h"
//_____________________________________________________________
//
// The L3 Binary File handler 
//

ClassImp(AliL3MemHandler)

AliL3MemHandler::AliL3MemHandler(){
  //Default constructor
  fPt = 0;
  fSize =0;
  fInBinary = 0;
  fOutBinary = 0;
  fNRandom = 0;
  Int_t r[2]={0,0};
  Init(0,0,r);
  IsRandom = kFALSE;
  fDigits = 0;
  fDPt =0;
  fNGenerate = 0;
  fNUsed = 0;
  fNDigits = 0;
  
  Int_t row[2] = {0,173};
  Init(0,0,row);
  ResetROI();
}


AliL3MemHandler::~AliL3MemHandler(){
  //Destructor
  if(fPt) delete[] fPt;
  if(fDigits) delete [] fDigits;
  if(fDPt) delete [] fDPt;
}

void AliL3MemHandler::ResetROI(){
    for(Int_t i=fRowMin; i<=fRowMax; i++)
      {
      fEtaMinTimeBin[i] = 0;
      fEtaMaxTimeBin[i] = 445;
      }
}


void AliL3MemHandler::SetROI(Float_t *eta,Int_t *slice){
 if(eta[1]==0)
     {
     LOG(AliL3Log::kWarning,"AliL3MemHandler::SetROI","Eta Values")
	 <<"Bad ROI parameters. IDIOT! "<<ENDLOG;
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
     
     xyz[0] = fTransformer->Row2X(i);
     xyz[1]=0;
     xyz[2] = xyz[0]/tan(thetamax);
     fTransformer->Slice2Sector(fSlice,i,sector,row);
     fTransformer->Local2Raw(xyz,sector,row);
     
     fEtaMinTimeBin[i] = (Int_t)xyz[2];
     
     if(eta[0]==0)
	 fEtaMaxTimeBin[i] = 445;
     else
	 {
	 Float_t thetamin = 2*atan(exp(-1.*eta[0]));
	 xyz[0] = fTransformer->Row2X(i);
	 xyz[1] = fTransformer->GetMaxY(i);
	 Float_t radii = sqrt(pow(xyz[0],2) + pow(xyz[1],2));
	 xyz[2] = radii/tan(thetamin);
	 fTransformer->Local2Raw(xyz,sector,row);
	 fEtaMaxTimeBin[i] = (Int_t)xyz[2];
	 }
     }
 
}

Bool_t AliL3MemHandler::SetBinaryInput(char *name){
  fInBinary = fopen(name,"r");
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::SetBinaryInput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliL3MemHandler::SetBinaryInput(FILE *file){
  fInBinary = file;
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::SetBinaryInput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliL3MemHandler::CloseBinaryInput(){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CloseBinaryInput","File Close")
    <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fInBinary);
  fInBinary =0;
}

Bool_t AliL3MemHandler::SetBinaryOutput(char *name){
  fOutBinary = fopen(name,"w");
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::SetBinaryOutput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliL3MemHandler::SetBinaryOutput(FILE *file){
  fOutBinary = file;
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::SetBinaryOutput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliL3MemHandler::CloseBinaryOutput(){
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CloseBinaryOutPut","File Close")
    <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fOutBinary);
  fOutBinary =0;
}


UInt_t AliL3MemHandler::GetFileSize(){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::GetFileSize","File")
    <<"No Input File"<<ENDLOG;
    return 0;
  }
  fseek(fInBinary,0,SEEK_END);
  UInt_t size = (UInt_t) ftell(fInBinary);
  rewind(fInBinary);
  return size; 
}

Byte_t *AliL3MemHandler::Allocate(){
  return Allocate(GetFileSize()); 
}

Byte_t *AliL3MemHandler::Allocate(AliL3TrackArray *array){
  if(!array){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Allocate","Memory")
    <<"Pointer to AliL3TrackArray = 0x0 "<<ENDLOG;
    return 0;
  }
  return Allocate(array->GetOutSize()); 
}

Byte_t *AliL3MemHandler::Allocate(UInt_t size){
  if(fPt){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Allocate","Memory")
    <<"Delete Memory"<<ENDLOG;
    Free();
  } 
  fPt = new Byte_t[size];
  fSize = size;
  LOG(AliL3Log::kDebug,"AliL3MemHandler::Allocate","Memory")
  <<AliL3Log::kDec<<"Allocate "<<size<<" Bytes of Memory"<<ENDLOG;
  return fPt;
}

void AliL3MemHandler::Free(){
  if(!fPt){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Free","Memory")
    <<"No Memory allocated - can't Free"<<ENDLOG;
    return;
  }  
  delete[] fPt;
  fPt = 0;
  fSize =0;
}

///////////////////////////////////////// Random
void AliL3MemHandler::SetRandomSeed(){
  time_t *tp=0;
  SetRandomSeed(time(tp));
}

void AliL3MemHandler::SetRandomCluster(Int_t maxnumber){
  IsRandom = kTRUE;
  fNRandom = maxnumber;
  fNDigits = 0;
  if(fDigits) delete [] fDigits;
  fDigits = new AliL3RandomDigitData[fNRandom*9];
  if(fDPt) delete [] fDPt;
  fDPt = (AliL3RandomDigitData **)new Int_t[fNRandom*9];
}

void AliL3MemHandler::QSort(AliL3RandomDigitData **a, Int_t first, Int_t last){

   // Sort array of AliL3RandomDigitData pointers using a quicksort algorithm.
   // Uses CompareDigits() to compare objects.
   // Thanks to Root!

   static AliL3RandomDigitData *tmp;
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

UInt_t AliL3MemHandler::GetRandomSize(){
  Int_t nrandom = 0;
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t npad=fTransformer->GetNPads(r);
    nrandom  += Int_t (fNGenerate * ((Double_t) npad/141.));
  }
  return 9 * nrandom * sizeof(AliL3DigitData);
}

void AliL3MemHandler::Generate(Int_t row){
  if(!IsRandom) return;
  ResetRandom();
  fNDigits = 0;
  Int_t npad=fTransformer->GetNPads(row);
  Int_t ntime = fEtaMaxTimeBin[row] - fEtaMinTimeBin[row];
  Int_t nrandom  = Int_t (fNGenerate * ((Double_t) npad/141.) * 
                   (Double_t) ntime/(Double_t) fTransformer->GetNTimeBins() );
 
  for(Int_t n=0;n<nrandom;n++){
    Int_t pad = (int)((float)rand()/RAND_MAX*npad);
    Int_t time =(int)((float)rand()/RAND_MAX*ntime+fEtaMinTimeBin[row] );
    Int_t charge = (int)((float)rand()/RAND_MAX*1023);
    DigitizePoint(row,pad,time,charge);
  }
  QSort(fDPt,0,fNDigits);
//  for(Int_t d=0;d<fNDigits;d++)
//    fprintf(stderr,"%d %d %d %d\n",fDPt[d]->fRow,fDPt[d]->fPad,
//                             fDPt[d]->fTime,fDPt[d]->fCharge);
}


void AliL3MemHandler::DigitizePoint(Int_t row, Int_t pad, 
                                       Int_t time,Int_t charge){
  for(Int_t j=-1;j<2;j++){
    for(Int_t k=-1;k<2;k++){
      Int_t dcharge = charge;
      if(j) dcharge /=2;
      if(k) dcharge /=2;
      if(dcharge<10) continue;
      Int_t dpad  = j + pad;
      Int_t dtime = k + time;

      if(dpad<0||dpad>=fTransformer->GetNPads(row))  continue;
      if(dtime<0||dtime>=fTransformer->GetNTimeBins()) continue;

      fDigits[fNDigits].fCharge = dcharge;
      fDigits[fNDigits].fRow = row;
      fDigits[fNDigits].fPad = dpad;
      fDigits[fNDigits].fTime = dtime;
      fDPt[fNDigits] = &fDigits[fNDigits];
      fNDigits++;
    }
  }
}

///////////////////////////////////////// Digit IO  
Bool_t AliL3MemHandler::Memory2Binary(UInt_t nrow,AliL3DigitRowData *data){
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2Binary","Memory")
    <<"Pointer to AliL3DigitRowData = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  AliL3DigitRowData *row_pt = data; 
  Int_t outsize = 0;
  for(UInt_t i=0;i<nrow;i++){
    Int_t size = sizeof(AliL3DigitData) * row_pt->fNDigit 
                                       + sizeof(AliL3DigitRowData);
    outsize += size;
    fwrite(row_pt,size,1,fOutBinary);
    Byte_t  *byte_pt =(Byte_t *) row_pt;
    byte_pt += size;
    row_pt = (AliL3DigitRowData *) byte_pt;
  }
  LOG(AliL3Log::kDebug,"AliL3MemHandler::Memory2Binary","Memory")
  <<AliL3Log::kDec<<"Wrote "<<outsize<<" Bytes to Memory ("
  <<nrow<<" Rows)"<<ENDLOG;
  return kTRUE;
}

Bool_t AliL3MemHandler::Binary2Memory(UInt_t & nrow,AliL3DigitRowData *data){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2Memory","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2Memory","Memory")
    <<"Pointer to AliL3DigitRowData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  rewind(fInBinary);
  AliL3DigitRowData *row_pt = data;
  UInt_t rowcount = 0;
  Int_t outsize =0;
  while(!feof(fInBinary)){
    Byte_t  *byte_pt =(Byte_t *) row_pt;

    if(fread(row_pt,sizeof(AliL3DigitRowData),1,fInBinary)!=1) break;
    byte_pt += sizeof(AliL3DigitRowData);
    outsize += sizeof(AliL3DigitRowData);

    Int_t size = sizeof(AliL3DigitData) * row_pt->fNDigit;
   
    if(fread(byte_pt,size,1,fInBinary)!=1) break;
    byte_pt += size;
    outsize += size;

    row_pt = (AliL3DigitRowData *) byte_pt;
    rowcount++;
  }  
  nrow= rowcount;
    LOG(AliL3Log::kDebug,"AliL3MemHandler::Binary2Memory","Memory")
    <<AliL3Log::kDec<<"Wrote "<<outsize<<" Bytes to Memory ("
    <<rowcount<<" Rows)"<<ENDLOG;
  return kTRUE;
}

void AliL3MemHandler::AddData(AliL3DigitData *data,UInt_t & ndata,
                   UInt_t row,UShort_t pad,UShort_t time,UShort_t charge){
  data[ndata].fPad = pad;
  data[ndata].fTime = time;
  data[ndata].fCharge = charge;
  ndata++;
}

void AliL3MemHandler::AddRandom(AliL3DigitData *data, UInt_t & ndata){
  data[ndata].fPad = fDPt[fNUsed]->fPad;
  data[ndata].fTime = fDPt[fNUsed]->fTime;
  data[ndata].fCharge = fDPt[fNUsed]->fCharge;
  ndata++;
  fNUsed++;
}

void AliL3MemHandler::MergeDataRandom(AliL3DigitData *data, UInt_t & ndata,
                   UInt_t row, UShort_t pad, UShort_t time, UShort_t charge){
//  AddData(data,ndata,row,pad,time,charge);
  data[ndata].fPad = pad;
  data[ndata].fTime = time;
  data[ndata].fCharge = charge;
  while(ComparePoints(row,pad,time)==0){
    Int_t ch = data[ndata].fCharge + fDPt[fNUsed]->fCharge;
    if(charge>1023) ch = 1023;
    data[ndata].fCharge = ch;
    fNUsed++;
  }
  ndata++;
}

void AliL3MemHandler::AddDataRandom(AliL3DigitData *data, UInt_t & ndata,
                   UInt_t row, UShort_t pad, UShort_t time, UShort_t charge){
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

void AliL3MemHandler::Write(UInt_t *comp, UInt_t & index, 
                                    UInt_t & subindex, UShort_t value){
  UInt_t shift[3] = {0,10,20};
  comp[index] |= (value&0x03ff)<<shift[subindex];
  if(subindex == 2){
    subindex = 0;
    index++;
  }
  else subindex++;
}

UShort_t AliL3MemHandler::Read(UInt_t *comp, UInt_t & index, UInt_t & subindex){
  UInt_t shift[3] = {0,10,20};
  UShort_t value = (comp[index]>>shift[subindex])&0x03ff;
  if(subindex == 2){
    subindex = 0;
    index++;
  }
  else subindex++;
  
  return value;
}

UShort_t AliL3MemHandler::Test(UInt_t *comp, 
                         UInt_t index, UInt_t  subindex){
  UInt_t shift[3] = {0,10,20};
  return (comp[index]>>shift[subindex])&0x03ff;
}

Int_t AliL3MemHandler::Memory2CompMemory(UInt_t nrow,
                                  AliL3DigitRowData *data,UInt_t *comp){
  if(!comp){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2CompMemory","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2CompMemory","Memory")
    <<"Pointer to AliL3DigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  AliL3DigitRowData *row_pt = data;
  UInt_t index=0;
  UInt_t subindex=0;
  
  for(UInt_t i=0;i<nrow;i++){
    UShort_t value = row_pt->fRow;
    Write(comp,index,subindex,value);
    UShort_t maxpad=0; 
    UShort_t npad=0;
    Int_t ddd[1000];
    for(Int_t d=0;d<200;d++) ddd[d]=0;
    for(UInt_t dig=0;dig<row_pt->fNDigit;dig++){
      if(row_pt->fDigitData[dig].fPad <200){ 
        ddd[row_pt->fDigitData[dig].fPad]++;
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
      if(digit>=row_pt->fNDigit || row_pt->fDigitData[digit].fPad !=  pad)
        continue;
      Write(comp,index,subindex,pad);
//    write zero if time != 0
      if(digit<row_pt->fNDigit && row_pt->fDigitData[digit].fPad == pad){
        if(row_pt->fDigitData[digit].fTime>0){
          Write(comp,index,subindex,0);
          Write(comp,index,subindex,row_pt->fDigitData[digit].fTime);
        }
      }
      while(digit<row_pt->fNDigit && row_pt->fDigitData[digit].fPad == pad){
        UShort_t charge = row_pt->fDigitData[digit].fCharge;
        if(charge>=1024){
          charge=1023;
        }
        Write(comp,index,subindex,charge);
        if(digit+1<row_pt->fNDigit&&row_pt->fDigitData[digit+1].fPad == pad){
          if(row_pt->fDigitData[digit].fTime +1 !=
                     row_pt->fDigitData[digit+1].fTime){
            Write(comp,index,subindex,0);
            UShort_t nzero = row_pt->fDigitData[digit+1].fTime - 
                             (row_pt->fDigitData[digit].fTime +1);
            Write(comp,index,subindex,nzero);
          }  
        }
        digit++;
      }
      Write(comp,index,subindex,0);
      Write(comp,index,subindex,0);
    }

    Int_t size = sizeof(AliL3DigitData) * row_pt->fNDigit+
                                            sizeof(AliL3DigitRowData);
    Byte_t  *byte_pt =(Byte_t *) row_pt;
    byte_pt += size;
    row_pt = (AliL3DigitRowData *) byte_pt;
  }
  while(subindex)
    Write(comp,index,subindex,0);
  return index * sizeof(UInt_t);
}

Int_t AliL3MemHandler::CompMemory2Memory(UInt_t  nrow,
                                   AliL3DigitRowData *data,UInt_t *comp){
  if(!comp){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompMemory2Memory","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompMemory2Memory","Memory")
    <<"Pointer to AliL3DigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  Int_t outsize=0;

  AliL3DigitRowData *row_pt = data;
  UInt_t index=0;
  UInt_t subindex=0;

  for(UInt_t i=0;i<nrow;i++){
    UInt_t ndigit=0;
    UInt_t row =Read(comp,index,subindex);
    row_pt->fRow=row;
    Generate(row);
    UShort_t npad = Read(comp,index,subindex);
    for(UShort_t p=0;p<npad;p++){
      UShort_t charge;
      UShort_t time =0;
      UShort_t pad = Read(comp,index,subindex);
      if(Test(comp,index,subindex)==0){
        Read(comp,index,subindex);
        if( (time = Read(comp,index,subindex)) == 0 ){
          continue;
        }
      }
      for(;;){
        while( (charge=Read(comp,index,subindex)) != 0){
          if(time>=fEtaMinTimeBin[row]&&time<=fEtaMaxTimeBin[row])
//          AddData(row_pt->fDigitData,ndigit,row,pad,time,charge);
            AddDataRandom(row_pt->fDigitData,ndigit,row,pad,time,charge);
          time++;
        }
        UShort_t tshift = Read(comp,index,subindex);
        if(tshift ==0) break;
        time += tshift;
      }
    }
    row_pt->fNDigit = ndigit;
    Int_t size = sizeof(AliL3DigitData) * row_pt->fNDigit+
                                        sizeof(AliL3DigitRowData);
    Byte_t  *byte_pt =(Byte_t *) row_pt;
    byte_pt += size;
    outsize += size;
    row_pt = (AliL3DigitRowData *) byte_pt;
  }
  return outsize;
}

UInt_t AliL3MemHandler::GetCompMemorySize(UInt_t nrow,AliL3DigitRowData *data){
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::GetCompMemorySize","Memory")
    <<"Pointer to AliL3DigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  AliL3DigitRowData *row_pt = data;
  UInt_t index=0;
  
  for(UInt_t i=0;i<nrow;i++){
    index++;
    UShort_t maxpad=0; 
    UShort_t npad=0;
    Int_t ddd[1000];
    for(Int_t d=0;d<200;d++) ddd[d]=0;
    for(UInt_t dig=0;dig<row_pt->fNDigit;dig++){
      if(row_pt->fDigitData[dig].fPad <200){ 
        ddd[row_pt->fDigitData[dig].fPad]++;
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
      if(digit>=row_pt->fNDigit || row_pt->fDigitData[digit].fPad !=  pad)
        continue;
      index++;
//    write zero if time != 0
      if(digit<row_pt->fNDigit && row_pt->fDigitData[digit].fPad == pad){
        if(row_pt->fDigitData[digit].fTime>0){
          index++;
          index++;
        }
      }
      while(digit<row_pt->fNDigit && row_pt->fDigitData[digit].fPad == pad){
        index++;
        if(digit+1<row_pt->fNDigit&&row_pt->fDigitData[digit+1].fPad == pad){
          if(row_pt->fDigitData[digit].fTime +1 !=
                     row_pt->fDigitData[digit+1].fTime){
            index++;
            index++;
          }  
        }
        digit++;
      }
      index++;
      index++;
    }

    Int_t size = sizeof(AliL3DigitData) * row_pt->fNDigit+
                                            sizeof(AliL3DigitRowData);
    Byte_t  *byte_pt =(Byte_t *) row_pt;
    byte_pt += size;
    row_pt = (AliL3DigitRowData *) byte_pt;
  }
  while(index%3)
    index++;
  return (index/3) * sizeof(UInt_t);
}

UInt_t AliL3MemHandler::GetMemorySize(UInt_t nrow,UInt_t *comp){
  if(!comp){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::GetMemorySize","Memory")
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
    Int_t size = sizeof(AliL3DigitData) * ndigit+
                                        sizeof(AliL3DigitRowData);
    outsize += size;
  }
   
  return outsize;
}

UInt_t AliL3MemHandler::GetNRow(UInt_t *comp,UInt_t size){
  if(!comp){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::GetNRow","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  size = size /4;
  UInt_t nrow=0;
  UInt_t index=0;
  UInt_t subindex=0;
  while(index+1<size){
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
  if(index+1==size)
  if(subindex<2)
  if(Read(comp,index,subindex)!=0) nrow++;

  return nrow;
}

Bool_t AliL3MemHandler::CompMemory2CompBinary(UInt_t nrow,UInt_t *comp,
                                                               UInt_t size){
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompMemory2CompBinary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!comp){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompMemory2CompBinary","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(size==0)
    size=GetMemorySize(nrow,comp);
  if(!size){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompMemory2CompBinary","Memory")
    <<"Memory size = 0 "<<ENDLOG;
    return kFALSE;
  }
  UInt_t length = size/sizeof(UInt_t);
  fwrite(&length,sizeof(UInt_t),1,fOutBinary);  
  fwrite(comp,size,1,fOutBinary);
  return kTRUE;
}

Bool_t AliL3MemHandler::CompBinary2CompMemory(UInt_t & nrow,UInt_t *comp){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompBinary2CompMemory","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!comp){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::CompBinary2CompMemory","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  rewind(fInBinary);
  UInt_t length;
  if(fread(&length,sizeof(UInt_t),1,fInBinary)!=1) return kFALSE;
  UInt_t size = length*sizeof(UInt_t);
  if(fread(comp,size,1,fInBinary)!=1) return kFALSE;
  // now find the number of dig
  nrow =  GetNRow(comp,size);
  return kTRUE;
}

AliL3DigitRowData *AliL3MemHandler::CompBinary2Memory(UInt_t & nrow){
  AliL3MemHandler * handler = new AliL3MemHandler();
  handler->SetBinaryInput(fInBinary);
  UInt_t *comp =(UInt_t *)handler->Allocate();
  handler->CompBinary2CompMemory(nrow,comp);
  UInt_t size = GetMemorySize(nrow,comp);
  AliL3DigitRowData *data = (AliL3DigitRowData *)Allocate(size);
  CompMemory2Memory(nrow,data,comp);
  handler->Free();
  delete handler;
  return data;  
}

Bool_t AliL3MemHandler::Memory2CompBinary(UInt_t nrow,AliL3DigitRowData *data){
  Bool_t out = kTRUE;
  AliL3MemHandler * handler = new AliL3MemHandler();
  UInt_t size = GetCompMemorySize(nrow,data);
  UInt_t *comp =(UInt_t *)handler->Allocate(size);
  Memory2CompMemory(nrow,data,comp);
  CompMemory2CompBinary(nrow,comp,size);
  handler->Free();
  delete handler;
  return out;
}


///////////////////////////////////////// Point IO  
Bool_t AliL3MemHandler::Memory2Binary(UInt_t npoint,AliL3SpacePointData *data){
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2Binary","Memory")
    <<"Pointer to AliL3SpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  UInt_t size = npoint*sizeof(AliL3SpacePointData);
  fwrite(data,size,1,fOutBinary);

  return kTRUE;
}

Bool_t AliL3MemHandler::Transform(UInt_t npoint,AliL3SpacePointData *data,
                                 Int_t slice, AliL3Transform* trans){
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Transform","Memory")
    <<"Pointer to AliL3SpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!trans){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Transform","Object")
    <<"Pointer to AliL3Transform = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  for(UInt_t i=0;i<npoint;i++){
    Float_t xyz[3];
    xyz[0] = data[i].fX;
    xyz[1] = data[i].fY;
    xyz[2] = data[i].fZ;
    trans->Local2Global(xyz,slice);
    data[i].fX = xyz[0];
    data[i].fY = xyz[1];
    data[i].fZ = xyz[2];
  }
  return kTRUE;
}

Bool_t AliL3MemHandler::Binary2Memory(UInt_t & npoint,AliL3SpacePointData *data){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2Memory","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2Memory","Memory")
    <<"Pointer to AliL3SpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  Int_t size = GetFileSize(); 
/*
  UInt_t  size,slice,patch,row[2];
  AliL3EventDataTypeRoot datatype;
  UInt_t node;
  if(fread(&datatype,sizeof(AliL3EventDataTypeRoot),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&node,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&size,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&slice,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&patch,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(row,2*sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
*/
  npoint = size/sizeof(AliL3SpacePointData);
  if(fread(data,size,1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(size%sizeof(AliL3SpacePointData)){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File Size")
    <<"File Size wrong "<<ENDLOG;
    return kFALSE; 
  }
  LOG(AliL3Log::kDebug,"AliL3MemHandler::Binary2Memory","File")
  <<AliL3Log::kDec<<"Wrote  "<<size<<" Bytes to Memory"<<ENDLOG;
  return kTRUE;
}

///////////////////////////////////////// Track IO  
Bool_t AliL3MemHandler::Memory2Binary(UInt_t ntrack,AliL3TrackSegmentData *data){
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2Binary","Memory")
    <<"Pointer to AliL3TrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliL3TrackSegmentData *track_pt = data;
  for(UInt_t i=0;i<ntrack;i++){
    Int_t size=sizeof(AliL3TrackSegmentData)+track_pt->fNPoints*sizeof(UInt_t); 
    fwrite(track_pt,size,1,fOutBinary);
    Byte_t *byte_pt = (Byte_t*) track_pt;
    byte_pt += size; 
    track_pt = (AliL3TrackSegmentData*) byte_pt;
  }
  LOG(AliL3Log::kDebug,"AliL3MemHandler::Memory2Binary","File")
  <<AliL3Log::kDec<<"Wrote  "<<ntrack<<" Tracks to File"<<ENDLOG;
  
  return kTRUE;
}

Bool_t AliL3MemHandler::Binary2Memory(UInt_t & ntrack,AliL3TrackSegmentData *data){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2Memory","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2Memory","Memory")
    <<"Pointer to AliL3TrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  ntrack=0;
  AliL3TrackSegmentData *track_pt = data;
  rewind(fInBinary);
/*
  UInt_t  size,slice,patch,row[2];
  AliL3EventDataTypeRoot datatype;
  UInt_t node;
  if(fread(&datatype,sizeof(AliL3EventDataTypeRoot),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&node,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&size,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&slice,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(&patch,sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(fread(row,2*sizeof(UInt_t),1,fInBinary)!=1){
    LOG(AliL3Log::kFatal,"AliL3MemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
*/
  while(!feof(fInBinary)){
    if(fread(track_pt,sizeof(AliL3TrackSegmentData),1,fInBinary)!=1) break;
    Int_t size=track_pt->fNPoints*sizeof(UInt_t);
    if(fread(track_pt->fPointIDs,size,1,fInBinary)!=1) break;
    Byte_t *byte_pt = (Byte_t*) track_pt;
    byte_pt += sizeof(AliL3TrackSegmentData)+size;
    track_pt = (AliL3TrackSegmentData*) byte_pt;
    ntrack++; 
  }
  LOG(AliL3Log::kDebug,"AliL3MemHandler::Binary2Memory","File")
  <<AliL3Log::kDec<<"Wrote  "<<ntrack<<" Tracks to Memory"<<ENDLOG;
  return kTRUE;
}

Bool_t AliL3MemHandler::TrackArray2Binary(AliL3TrackArray *array){
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::TrackArray2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::TrackArray2Binary","Memory")
    <<"Pointer to AliL3TrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliL3TrackSegmentData *data = (AliL3TrackSegmentData *)Allocate(array);
  UInt_t ntrack;
  TrackArray2Memory(ntrack,data,array);
  Memory2Binary(ntrack,data);
  Free();
  return kTRUE;
}

Bool_t AliL3MemHandler::Binary2TrackArray(AliL3TrackArray *array){
  if(!fInBinary){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2TrackArray","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Binary2TrackArray","Memory")
    <<"Pointer to AliL3TrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliL3TrackSegmentData *data = (AliL3TrackSegmentData *)Allocate();
  UInt_t ntrack;
  Binary2Memory(ntrack,data);
  Memory2TrackArray(ntrack,data,array);  
  Free();
  return kTRUE;
}

Bool_t AliL3MemHandler::TrackArray2Memory(UInt_t & ntrack,AliL3TrackSegmentData *data,AliL3TrackArray *array){
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::TrackArray2Memory","Memory")
    <<"Pointer to AliL3TrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::TrackArray2Memory","Memory")
    <<"Pointer to AliL3TrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->WriteTracks(ntrack,data);
  return kTRUE;
}

Bool_t AliL3MemHandler::Memory2TrackArray(UInt_t ntrack,AliL3TrackSegmentData *data,AliL3TrackArray *array){
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliL3TrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliL3TrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->FillTracks(ntrack,data);
  return kTRUE;
}

Bool_t AliL3MemHandler::Memory2TrackArray(UInt_t ntrack,AliL3TrackSegmentData *data,AliL3TrackArray *array,Int_t slice, AliL3Transform* trans){
  if(!data){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliL3TrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliL3TrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!trans){
    LOG(AliL3Log::kWarning,"AliL3MemHandler::Memory2TrackArray","Object")
    <<"Pointer to AliL3Transform = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->FillTracks(ntrack,data,slice,trans);
  return kTRUE;
}

