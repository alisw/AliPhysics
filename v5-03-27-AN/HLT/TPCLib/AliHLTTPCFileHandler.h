// $Id$
// Original: AliHLTFileHandler.h,v 1.19 2004/06/11 16:06:33 loizides 

#ifndef ALIHLTTPCFILEHANDLER_H
#define ALIHLTTPCFILEHANDLER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCFileHandler.h
/// @author U. Frankenfeld, A. Vestbo, C. Loizides, maintained by
///         Matthias Richter
/// @date   
/// @brief  file input for the TPC tracking code before migration to the
///         HLT component framework
///

#include "AliHLTTPCMemHandler.h"

class TClonesArray;

class AliSimDigits;
class AliTPCParam;
#include <AliRunLoader.h>

class TObject;
class TFile;
class TTree;

struct AliHLTTPCSpacePointData;
struct AliHLTTPCDigitRowData;
struct AliHLTTPCTrackSegmentData;
class AliHLTTPCTrackArray;

/**
 * @class AliHLTTPCFileHandler
 * This is the input interface class for the TPC tracking code before conversion to
 * the HLT component framework.
 * 
 * @ingroup alihlt_tpc
 */
class AliHLTTPCFileHandler:public AliHLTTPCMemHandler {

 public:
  /** standard constructor */
  AliHLTTPCFileHandler(Bool_t b=kFALSE);
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

  /**
   * Write AliDigits from AliRoot file to binary file.
   * @param event      event no
   * @param altro      use @ref AliDigits2Memory if kFALSE and @ref
   *                   AliDigits2Memory if kTRUE
   *
   * Calls the @ref AliHLTTPCMemHandler::Memory2BinaryFile to write the file.
   */
  Bool_t AliDigits2BinaryFile(Int_t event=0,Bool_t altro=kFALSE);

  /**
   * Convert AliDigits from AliRoot file to HLT Digit data in memory.
   * Read and convert/write digits to memory buffer. If no target buffer available,
   * an appropriate buffer is allocated.<br>
   * If the variable pTgtSize is prvided, the total size of the result array is
   * returned. \b Note: the total size differs as the @ref AliHLTTPCDigitRowData
   * structs are variable in size depending on the no of digits for that particular
   * row.
   * @param nrow       [OUT] number of rows
   * @param event      the event no
   * @param tgtBuffer  target buffer (optional)
   * @param pTgtSize   size of target buffer (optional)
   * @return pointer to array, size in nrow <br>
   *         NULL in case of failure, required size in pTgtSize
   */
  //TODO: Check that the following change works. It should, but just double check.
  AliHLTTPCDigitRowData *AliDigits2Memory(UInt_t & nrow,Int_t event, Byte_t* tgtBuffer, UInt_t* pTgtSize=NULL);

  AliHLTTPCDigitRowData *AliDigits2Memory(UInt_t & nrow,Int_t event=0)
  {
    return AliDigits2Memory(nrow, event, NULL, NULL);
  }

  /**
   * Convert and filter AliDigits from AliRoot file to HLT Digit data in memory.
   * This functions is the same as @ref AliDigits2Memory but in addition it
   * filters out single timebins, which is noise. The timebins which
   * are removed are timebins which have the 4 zero neighbours; 
   * (pad-1,time),(pad+1,time),(pad,time-1),(pad,time+1).
   *
   * This is legacy code, the two functions contain big portions of identical code
   * will be merged.
   * See @ref AliDigits2Memory for detailed description.
   */
  AliHLTTPCDigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t eventmerge=kFALSE); 

  /**
   * Convert AliDigits from AliRoot file to Altro data format in memory.
   */
  int AliDigits2Altro(Int_t event, Byte_t* tgtBuffer, UInt_t size);

  /**
   * Write AliDigits from AliRoot file to binary file.
   * @param event      event no
   * @param altro      use @ref AliDigits2Memory if kFALSE and @ref
   *                   AliDigits2Memory if kTRUE
   *
   * \b Note: pretty much the same as @ref AliDigits2BinaryFile.
   * Calls the @ref AliHLTTPCMemHandler::Memory2CompBinary to write the file.
   */
  Bool_t AliDigits2CompBinary(Int_t event=0,Bool_t altro=kFALSE);  

  //Point IO
  Bool_t AliPoints2Binary(Int_t eventn=0);
  AliHLTTPCSpacePointData *AliPoints2Memory(UInt_t & npoint,Int_t eventn=0);//Allocates Memory

 protected:
  AliRunLoader *fInAli;//!
  Bool_t fUseRunLoader; //use runloader

  AliTPCParam *fParam;//!
  AliSimDigits *fDigits;//!

  TTree *fDigitsTree;//!
  FILE *fMC;//!

  Bool_t fIndexCreated;   //is index created
  Int_t  fIndex[fgkNSlice][fgkNRow]; //stores index over digitstree 
                          //for faster access w/o ASVVERSION
  Bool_t fUseStaticIndex; //take static index
  static Bool_t fgStaticIndexCreated;   //global index created
  static Int_t  fgStaticIndex[fgkNSlice][fgkNRow]; //global index

  virtual Bool_t SetAliInput();
  Bool_t GetDigitsTree(Int_t event);
  Bool_t CreateIndex();  //create the index

 private:
  /** copy constructor prohibited */
  AliHLTTPCFileHandler(const AliHLTTPCFileHandler&);
  /** assignment operator prohibited */
  AliHLTTPCFileHandler& operator=(const AliHLTTPCFileHandler&);

  ClassDef(AliHLTTPCFileHandler,0)   //HLT TPC Filehandler IO class
};

#endif
