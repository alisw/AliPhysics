#ifndef ALIFMDSDigitizer_H
#define ALIFMDSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Task Class for making SDigits in FMD      
//                  
//-- Author: Alla Maevskaia(INR)


// --- ROOT system ---
#include "TTask.h"
#include "TString.h"
#include "AliFMD.h"
#include "AliDetector.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliRunLoader;

class AliFMDSDigitizer: public TTask {

public:
  AliFMDSDigitizer() ;          // ctor
  AliFMDSDigitizer(const char* HeaderFile,char *SdigitsFile = 0) ; 

  virtual ~AliFMDSDigitizer() ; // dtor
  // Int_t    Digitize(Float_t Energy);

  char *GetSDigitsFile()const{return (char*) fSDigitsFile.Data();}  
  virtual void  Exec(Option_t *option); 
  void SetNEvents(Int_t Nevents){fNevents = Nevents;}
  Stat_t GetNEvents(){return fNevents;}
  void SetSDigitsFile(char * file ) ;
  virtual void Print(Option_t* option) const ;
  TClonesArray *SDigits() const {return fSDigits;}
  TClonesArray *Hits() const {return fHits;}
  // Granularity
   virtual void SetRingsSi1(Int_t ringsSi1);
   virtual void SetSectorsSi1(Int_t sectorsSi1);
   virtual void SetRingsSi2(Int_t ringsSi2);
   virtual void SetSectorsSi2(Int_t sectorsSi2);


private:
  Int_t   fNevents ;        // Number of events to digitize
  TString fSDigitsFile ;    //output file 
  TClonesArray *fSDigits      ; // List of summable digits
  TClonesArray *fHits      ; // List of summable digits
  TString fHeadersFile ;    //input file

  AliRunLoader *fRunLoader;//!Run Loader

 protected:
  //Granularity
   Int_t fRingsSi1;       // Number of rings
   Int_t fSectorsSi1;    // Number of sectors
   Int_t fRingsSi2;       // Number of rings
  Int_t fSectorsSi2;    // Number of sectors
  

  ClassDef(AliFMDSDigitizer,1)  // description 

};

#endif // AliFMDSDigitizer_H
