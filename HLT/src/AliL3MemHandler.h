#ifndef ALIL3MEMHANDLER_H
#define ALIL3MEMHANDLER_H

#include "AliL3RootTypes.h"
#include <stdio.h>
#include <stdlib.h>
#include "AliL3DigitData.h"
#include "AliL3Transform.h"
class AliL3SpacePointData;
class AliL3DigitRowData;
class AliL3TrackSegmentData;
class AliL3TrackArray;
class AliL3RandomPointData;

class AliL3MemHandler{
 private:
  AliL3Transform *fTransformer;//!
  FILE *fInBinary;//!
  FILE *fOutBinary;//!
  Byte_t *fPt;//!
  UInt_t fSize;

  AliL3RandomDigitData **fDPt;//!
  AliL3RandomDigitData *fDigits;//!
  Bool_t IsRandom;
  Int_t fNRandom;
  Int_t fNGenerate;
  Int_t fNUsed;
  Int_t fNDigits;

  void Write(UInt_t *comp, UInt_t & index, UInt_t & subindex, UShort_t value);
  UShort_t Read(UInt_t *comp, UInt_t & index, UInt_t & subindex);
  UShort_t Test(UInt_t *comp, UInt_t index, UInt_t subindex); 
  
  void DigitizePoint(Int_t row,Int_t pad, Int_t time,Int_t charge);
  void QSort(AliL3RandomDigitData **a, Int_t first, Int_t last);
  Int_t ComparePoints(UInt_t row,UShort_t pad,UShort_t time);
  Int_t CompareDigits(AliL3RandomDigitData *a,AliL3RandomDigitData *b);
  void AddData(AliL3DigitData *data,UInt_t & ndata,
                      UInt_t row,UShort_t pad,UShort_t time,UShort_t charge);
  void AddRandom(AliL3DigitData *data,UInt_t & ndata);
  void MergeDataRandom(AliL3DigitData *data,UInt_t & ndata,
                      UInt_t row,UShort_t pad,UShort_t time,UShort_t charge);
  void AddDataRandom(AliL3DigitData *data,UInt_t & ndata,
                      UInt_t row,UShort_t pad,UShort_t time,UShort_t charge);

 protected:
  Int_t fRowMin;
  Int_t fRowMax;
  Int_t fSlice;
  Int_t fPatch;

  Int_t fEtaMinTimeBin[174];
  Int_t fEtaMaxTimeBin[174];
  

 public:
  AliL3MemHandler();
  virtual ~AliL3MemHandler();
  
  void Reset(){CloseBinaryInput();CloseBinaryOutput();Free();}  
  void Init(AliL3Transform *t){fTransformer = t;}
  void Init(Int_t s,Int_t p,Int_t* row){fSlice=s;fPatch=p;fRowMin=row[0];fRowMax=row[1]; ResetROI();}

  Bool_t SetBinaryInput(char *name);
  Bool_t SetBinaryInput(FILE *file);
  void CloseBinaryInput();
  
  Bool_t SetBinaryOutput(char *name);
  Bool_t SetBinaryOutput(FILE *file);
  void CloseBinaryOutput();

  //Random cluster
  void SetRandomCluster(Int_t maxnumber);
  void SetRandomSeed(UInt_t seed){srand(seed);}
  void SetRandomSeed();

  void ResetRandom(){fNDigits = 0; fNUsed = 0;}
  void Generate(Int_t row);
  void SetNGenerate(Int_t number){(number>fNRandom)?fNGenerate=fNRandom:fNGenerate = number;}

  void SetROI(Float_t *eta,Int_t *slice);
  void ResetROI();

  //Digit IO
  Bool_t Memory2Binary(UInt_t nrow,AliL3DigitRowData *data);
  Bool_t Binary2Memory(UInt_t & nrow,AliL3DigitRowData *data);

  Int_t Memory2CompMemory(UInt_t nrow,AliL3DigitRowData *data,UInt_t *comp);
  Int_t CompMemory2Memory(UInt_t nrow,AliL3DigitRowData *data,UInt_t *comp);
  Bool_t CompMemory2CompBinary(UInt_t nrow,UInt_t *comp, UInt_t size=0);
  Bool_t CompBinary2CompMemory(UInt_t & nrow,UInt_t *comp);

  AliL3DigitRowData *CompBinary2Memory(UInt_t & nrow);
  Bool_t Memory2CompBinary(UInt_t nrow,AliL3DigitRowData *data);
  UInt_t GetNRow(UInt_t *comp,UInt_t size);

  //Point IO
  Bool_t Memory2Binary(UInt_t npoint,AliL3SpacePointData *data);
  Bool_t Binary2Memory(UInt_t & npoint,AliL3SpacePointData *data);
  Bool_t Transform(UInt_t npoint,AliL3SpacePointData *data,
                              Int_t slice, AliL3Transform* trans);

  //Track IO
  Bool_t Memory2Binary(UInt_t ntrack,AliL3TrackSegmentData *data);
  Bool_t Binary2Memory(UInt_t & ntrack,AliL3TrackSegmentData *data);
  Bool_t TrackArray2Binary(AliL3TrackArray *array);
  Bool_t Binary2TrackArray(AliL3TrackArray *array);
  Bool_t TrackArray2Memory(UInt_t & ntrack,AliL3TrackSegmentData *data,
                                             AliL3TrackArray *array);
  Bool_t Memory2TrackArray(UInt_t ntrack,AliL3TrackSegmentData *data,
                                             AliL3TrackArray *array);
  Bool_t Memory2TrackArray(UInt_t ntrack,AliL3TrackSegmentData *data,
                   AliL3TrackArray *array,Int_t slice, AliL3Transform* trans);

  //Memory Allocation
  UInt_t GetAllocatedSize(){return fSize;}  
  UInt_t GetFileSize();
  UInt_t GetMemorySize(UInt_t nrow,UInt_t *comp);
  UInt_t GetCompMemorySize(UInt_t nrow,AliL3DigitRowData *data);
  UInt_t GetRandomSize();

  Byte_t *Allocate(UInt_t size);
  Byte_t *Allocate();  // allocate size of Binary Input File
  Byte_t *Allocate(AliL3TrackArray *array);
  void   Free();
  
  ClassDef(AliL3MemHandler,1) 
};

inline Int_t  AliL3MemHandler::ComparePoints(UInt_t row,UShort_t pad,UShort_t time){
  if(fNUsed>=fNDigits) return -2;

  if(pad==fDPt[fNUsed]->fPad&&time==fDPt[fNUsed]->fTime) return 0;

  if(pad<fDPt[fNUsed]->fPad) return -1;
  if(pad==fDPt[fNUsed]->fPad&&time<fDPt[fNUsed]->fTime)  return -1;

  return 1;
}

inline Int_t AliL3MemHandler::CompareDigits(AliL3RandomDigitData *a,AliL3RandomDigitData *b){
  if(a->fPad==b->fPad && a->fTime == b->fTime) return 0;

  if(a->fPad<b->fPad) return -1;
  if(a->fPad==b->fPad && a->fTime<b->fTime) return -1;
  
  return 1;
}

#endif
