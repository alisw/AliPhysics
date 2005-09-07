// @(#) $Id$

#ifndef ALIHLTTPCFILEHANDLER_H
#define ALIHLTTPCFILEHANDLER_H

#include "AliHLTTPCMemHandler.h"

class TClonesArray;

#include <AliSimDigits.h>
#include <AliTPCParam.h>
#ifdef use_newio
#include <AliRunLoader.h>
#endif

#include <TObject.h>
#include <TFile.h>
#include <TTree.h>

class AliHLTTPCSpacePointData;
class AliHLTTPCDigitRowData;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCTrackArray;

class AliHLTTPCFileHandler:public AliHLTTPCMemHandler{

 protected:
#ifdef use_newio
  AliRunLoader *fInAli;//!
  Bool_t fUseRunLoader; //use runloader
#else
  TFile *fInAli;//!
#endif

  AliTPCParam *fParam;//!
  AliSimDigits *fDigits;//!

  TTree *fDigitsTree;//!
  FILE *fMC;//!
  
  Bool_t fIndexCreated;   //is index created
  Int_t  fIndex[36][159]; //stores index over digitstree 
                          //for faster access w/o ASVVERSION
  Bool_t fUseStaticIndex; //take static index
  static Bool_t fgStaticIndexCreated;   //global index created
  static Int_t  fgStaticIndex[36][159]; //global index

  virtual Bool_t SetAliInput();
  Bool_t GetDigitsTree(Int_t event);
  Bool_t CreateIndex();  //create the index

 public:
  AliHLTTPCFileHandler(Bool_t b=kFALSE);
  ~AliHLTTPCFileHandler();

  void FreeDigitsTree();
  static void CleanStaticIndex();
  static Int_t SaveStaticIndex(Char_t *prefix=0,Int_t event=0);
  static Int_t LoadStaticIndex(Char_t *prefix=0,Int_t event=0);

  Bool_t SetAliInput(Char_t *name);
  Bool_t SetAliInput(TFile *file);
#ifdef use_newio
  Bool_t SetAliInput(AliRunLoader *runLoader);
#else
#endif
  void CloseAliInput(); 
  Bool_t IsDigit(Int_t event);
  
  Bool_t SetMCOutput(Char_t *name);
  Bool_t SetMCOutput(FILE *file);
  void CloseMCOutput();

  //Digit IO
  Bool_t AliDigits2Binary(Int_t event=0,Bool_t altro=kFALSE);
  AliHLTTPCDigitRowData *AliDigits2Memory(UInt_t & nrow,Int_t event=0); //Allocates Memory
  AliHLTTPCDigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t eventmerge=kFALSE); 
  //Allocates Memory
  Bool_t AliDigits2CompBinary(Int_t event=0,Bool_t altro=kFALSE);  
  void AliDigits2RootFile(AliHLTTPCDigitRowData *rowPt,Char_t *new_digitsfile);

  //Point IO
  Bool_t AliPoints2Binary(Int_t eventn=0);
  AliHLTTPCSpacePointData *AliPoints2Memory(UInt_t & npoint,Int_t eventn=0);//Allocates Memory

  ClassDef(AliHLTTPCFileHandler,1)   //Filehandler class
};

#endif
