#ifndef ALITOFSDigitizer_H
#define ALITOFSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
//  Task Class for making SDigits in TOF      
//                  
//-- Author: F. Pierella


#include "TTask.h"
#include "TString.h"
#include "AliTOF.h"
#include "AliDetector.h"

class AliTOFSDigitizer: public TTask {

public:
  AliTOFSDigitizer() ;          // ctor
  AliTOFSDigitizer(char* HeaderFile, char *SdigitsFile = 0) ; 

  virtual ~AliTOFSDigitizer() ; // dtor
  // Int_t    Digitize(Float_t Energy);

//  char *GetSDigitsFile() const {return const_cast<char*>(fSDigitsFile.Data());}  
  const char *GetSDigitsFile() const {return fSDigitsFile.Data();}  
  virtual void  Exec(Option_t *option); 
  void  SetNEvents(Int_t Nevents) {fNevents = Nevents;}
  Int_t GetNEvents() const {return fNevents;}
  void SetSDigitsFile(char * file ) ;
  virtual void Print(Option_t* option) const ;
  TClonesArray *SDigits() const {return fSDigits;}
  TClonesArray *Hits() const {return fHits;}



private:
  Int_t   fNevents;         // Number of events to digitize
  TString fSDigitsFile;     // output file 
  TClonesArray *fSDigits;   // array of summable digits
  TClonesArray *fHits;      // array of summable digits
  TString fHeadersFile;     // input file

 protected:


  ClassDef(AliTOFSDigitizer,1)  // creates TOF SDigits

};

#endif // AliTOFSDigitizer_H
