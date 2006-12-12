// @(#) $Id$
// Original: AliHLTFileHandler.h,v 1.19 2004/06/11 16:06:33 loizides 

#ifndef ALIHLTTPCFILEHANDLER_H
#define ALIHLTTPCFILEHANDLER_H

#include "AliHLTTPCMemHandler.h"

class TClonesArray;

class AliSimDigits;
class AliTPCParam;
#include <AliRunLoader.h>

class TObject;
class TFile;
class TTree;
#include "AliHLTLogging.h"

class AliHLTTPCSpacePointData;
class AliHLTTPCDigitRowData;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCTrackArray;

/**
 * class AliHLTTPCFileHandler
 * This is the input interface class for the TPC tracking code before conversion to
 * the HLT component framework.
 * 
 */
class AliHLTTPCFileHandler:public AliHLTTPCMemHandler, public AliHLTLogging {

 protected:
  AliRunLoader *fInAli;//!
  Bool_t fUseRunLoader; //use runloader

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
  /** standard constructor */
  AliHLTTPCFileHandler(Bool_t b=kFALSE);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCFileHandler(const AliHLTTPCFileHandler&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCFileHandler& operator=(const AliHLTTPCFileHandler&);
  /** destructor */
  virtual ~AliHLTTPCFileHandler();

  void FreeDigitsTree();
  static void CleanStaticIndex();
  static Int_t SaveStaticIndex(Char_t *prefix=0,Int_t event=0);
  static Int_t LoadStaticIndex(Char_t *prefix=0,Int_t event=0);

  Bool_t SetAliInput(Char_t *name);
  Bool_t SetAliInput(AliRunLoader *runLoader);
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
  void AliDigits2RootFile(AliHLTTPCDigitRowData *rowPt,Char_t *newDigitsfile);

  //Point IO
  Bool_t AliPoints2Binary(Int_t eventn=0);
  AliHLTTPCSpacePointData *AliPoints2Memory(UInt_t & npoint,Int_t eventn=0);//Allocates Memory

  ClassDef(AliHLTTPCFileHandler,1)   //Filehandler class
};

#endif
