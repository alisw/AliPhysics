// @(#) $Id$

#ifndef ALIL3FILEHANDLER_H
#define ALIL3FILEHANDLER_H

#include "AliL3MemHandler.h"

class TClonesArray;

#include <AliSimDigits.h>
#include <AliTPCParam.h>
#include <AliRunLoader.h>

#include <TObject.h>
#include <TTree.h>

class AliL3SpacePointData;
class AliL3DigitRowData;
class AliL3TrackSegmentData;
class AliL3TrackArray;

class AliL3FileHandler:public AliL3MemHandler{
 private:
  AliRunLoader *fInAli;
  AliTPCParam *fParam;
  Bool_t SetAliInput();
  Int_t fLastIndex;
  AliSimDigits *fDigits;
  TTree *fDigitsTree;
  FILE *fMC;//!
  
  Bool_t fIndexCreated;   //is index created
  Int_t  fIndex[36][159]; //stores index over digitstree 
                          //for faster access w/o ASVVERSION

  Bool_t GetDigitsTree(Int_t event);
  Bool_t CreateIndex();  //create the index

 public:
  AliL3FileHandler();
  virtual ~AliL3FileHandler();

  void FreeDigitsTree();
  Bool_t SetAliInput(Char_t *name);
  Bool_t SetAliInput(AliRunLoader *runLoader);
  void CloseAliInput(); 
  Bool_t IsDigit(Int_t event);
  
  Bool_t SetMCOutput(Char_t *name);
  Bool_t SetMCOutput(FILE *file);
  void CloseMCOutput();

  //Digit IO
  Bool_t AliDigits2Binary(Int_t event=0,Bool_t altro=kFALSE);
  AliL3DigitRowData *AliDigits2Memory(UInt_t & nrow,Int_t event=0); //Allocates Memory
  AliL3DigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t eventmerge=kFALSE); //Allocates Memory
  Bool_t AliDigits2CompBinary(Int_t event=0,Bool_t altro=kFALSE);  
  void AliDigits2RootFile(AliL3DigitRowData *rowPt,Char_t *new_digitsfile);

  //Point IO
  Bool_t AliPoints2Binary(Int_t eventn=0);
  AliL3SpacePointData *AliPoints2Memory(UInt_t & npoint,Int_t eventn=0);//Allocates Memory

  ClassDef(AliL3FileHandler,1)   //Filehandler class
};

#endif
