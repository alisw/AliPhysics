#ifndef ALIL3FILEHANDLER_H
#define ALIL3FILEHANDLER_H

#include "AliTPCParam.h"
#include "AliL3MemHandler.h"
#include "AliSimDigits.h"
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <stdio.h>

class AliL3SpacePointData;
class AliL3DigitRowData;
class AliL3TrackSegmentData;
class AliL3TrackArray;

class AliL3FileHandler:public AliL3MemHandler{
 private:
  TFile *fInAli;
  AliTPCParam *fParam;
  Bool_t SetAliInput();
  Int_t fLastIndex;
  AliSimDigits *fDigits;
  TTree *fDigitsTree;
  FILE *fMC;//!
  
  Bool_t GetDigitsTree(Int_t event);
  
 public:
  AliL3FileHandler();
  virtual ~AliL3FileHandler();

  void FreeDigitsTree();
  Bool_t SetAliInput(char *name);
  Bool_t SetAliInput(TFile *file);
  void CloseAliInput(); 
  Bool_t IsDigit(Int_t event);
  
  Bool_t SetMCOutput(char *name);
  Bool_t SetMCOutput(FILE *file);
  void CloseMCOutput();

  //Digit IO
  Bool_t AliDigits2Binary(Int_t event=0);
  AliL3DigitRowData *AliDigits2Memory(UInt_t & nrow,Int_t event=0); //Allocates Memory
  Bool_t AliDigits2CompBinary(Int_t event=0);  
  void AliDigits2RootFile(AliL3DigitRowData *rowPt,Char_t *new_digitsfile);

  //Point IO
  Bool_t AliPoints2Binary();
  AliL3SpacePointData *AliPoints2Memory(UInt_t & npoint);//Allocates Memory

  ClassDef(AliL3FileHandler,1)   //Filehandler class
};

#endif
