#ifndef ALITOFMERGER_H
#define ALITOFMERGER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// #include "AliMerger.h"
// #include "AliMergable.h"
#include "TRandom.h"
#include "AliDetector.h"
#include "AliTOF.h"
//Piotr.Skowronski@cern.ch :
//Corrections applied in order to compile (only) with new I/O and folder structure
//To be implemented correctly by responsible
class AliRunLoader;

typedef enum {kDigitize=0, kMerge = 1} MergeMode_t;

class AliTOFMerger {
 public:
  
  AliTOFMerger();
  virtual ~AliTOFMerger();
  
  
  // Initialize merging and digitisation
  virtual void Init();
  
  // Do the main work
  void Digitise() ;
  TClonesArray *SDigits() const {return fSDigits;}
 
  void ReadWriteDigit(Int_t);
  
  // Setters -> Later Communication with gAlice 
  void SetSignalEventNumber(Int_t i)     {fEvNrSig = i;}
  void SetBackgroundEventNumber(Int_t i) {fEvNrBgr = i;}    
  void SetBackgroundFileName(char* file) {fFnBgr = file;}        
  void SetSignalFileName(char* file)     {fFnSig = file;}        
  void SetMode(MergeMode_t mode) {fMerge = mode;}

    enum {kBgTag = -1};
      
 private:    
    // Open the bgr file
    TFile *InitBgr();

   
 private:
    TClonesArray *fDigits;               // ! array with digits
    TClonesArray *fSDigits      ; // List of summable digits
    Int_t fNDigits;                 // number of digits
    Int_t fEvNrSig;                 // signal     event number
    Int_t fEvNrBgr;                 // background event number    
    MergeMode_t fMerge;             // merging type kDigitize, kMerge
    char  *fFnBgr;                  // background file name
    char  *fFnSig;                  // signal file name
    TFile *fBgrFile;                // Pointer to background file
    
    AliRunLoader * fRunLoader;        //! Run Loader 
    ClassDef(AliTOFMerger,0)
};    
#endif

