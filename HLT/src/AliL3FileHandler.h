#ifndef ALIL3FILEHANDLER_H
#define ALIL3FILEHANDLER_H

#include "AliTPCParam.h"
#include "AliL3MemHandler.h"
#include <TObject.h>
#include <TFile.h>
#include <stdio.h>
class AliL3SpacePointData;
class AliL3DigitRowData;
class AliL3TrackSegmentData;
class AliL3TrackArray;

class AliL3FileHandler:public AliL3MemHandler{
 private:
  TFile *fInAli;
  AliTPCParam *fParam;
  AliL3Transform *fTransformer;//!
  Bool_t SetAliInput();

  FILE *fMC;//!

 public:
  AliL3FileHandler();
  virtual ~AliL3FileHandler();
//  void Init(Int_t s,Int_t p,Int_t* row){fSlice=s;fPatch=p;fRowMin=row[0];fRowMax=row[1];}

  Int_t GetRowMin(){return fRowMin;}
  Int_t GetRowMax(){return fRowMax;}
  Int_t GetSlice(){return fSlice;}
  Int_t GetPatch(){return fPatch;}

  Bool_t SetAliInput(char *name);
  Bool_t SetAliInput(TFile *file);
  void CloseAliInput(); 
  Bool_t IsDigit();
 
  Bool_t SetMCOutput(char *name);
  Bool_t SetMCOutput(FILE *file);
  void CloseMCOutput();

  //Digit IO
  Bool_t AliDigits2Binary();
  AliL3DigitRowData *AliDigits2Memory(UInt_t & nrow); //Allocates Memory
  Bool_t AliDigits2CompBinary();  

  //Point IO
  Bool_t AliPoints2Binary();
  AliL3SpacePointData *AliPoints2Memory(UInt_t & npoint);//Allocates Memory

  ClassDef(AliL3FileHandler,1)   // Level3 
};

#endif
