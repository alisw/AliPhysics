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

class AliFMDSDigitizer: public TTask {

public:
  AliFMDSDigitizer() ;          // ctor
  AliFMDSDigitizer(char* HeaderFile,char *SdigitsFile = 0) ; 

  virtual ~AliFMDSDigitizer() ; // dtor
  // Int_t    Digitize(Float_t Energy);

  char *GetSDigitsFile()const{return (char*) fSDigitsFile.Data();}  
  virtual void  Exec(Option_t *option); 
  void SetNEvents(Int_t Nevents){fNevents = Nevents;}
  Stat_t GetNEvents(){return fNevents;}
  void SetSDigitsFile(char * file ) ;
  virtual void Print(Option_t* option) const ;

private:
  Int_t   fNevents ;        // Number of events to digitize
  TString fSDigitsFile ;    //output file 
  TString fHeadersFile ;    //input file


  ClassDef(AliFMDSDigitizer,2)  // description 

};

#endif // AliFMDSDigitizer_H
