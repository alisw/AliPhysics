#ifndef ALIVZERODigitizer_H
#define ALIVZERODigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
//_________________________________________________________________________
//
//  Class for making Digits in VZERO 
//_________________________________________________________________________   


#include "AliDigitizer.h"
#include "TString.h"

class TClonesArray;
class TFile;
class TMath;
class TObjArray;
class TParticle;
class TTree;
class TNtuple;

class AliLoader;
class AliRunLoader;
class AliRun;
class AliDetector;
class AliVZEROhit;
class AliHit;
class AliHeader;
class AliRunDigitizer;

class AliVZEROcell;
class AliVZEROsdigit;
class AliVZEROdigit;

// --- Standard library ---

// --- AliRoot header files ---

class AliRunLoader;

class AliVZERODigitizer: public AliDigitizer {

public:

  AliVZERODigitizer() ;                       // constructor
  AliVZERODigitizer(AliRunDigitizer *manager);// constructor
  virtual ~AliVZERODigitizer() ;              // destructor
  
  void OpengAliceFile(const char *file);  
  char *GetDigitsFile()const{return (char*) fDigitsFile.Data();}  
  virtual void  Exec();
  void AddDigit(Int_t /* eventnumber */, Int_t /* cellnumber */, Int_t /* adc */);		
  void SetNEvents(Int_t Nevents){fNevents = Nevents;}
  void ResetDigit();
  Stat_t GetNEvents(){return fNevents;}

 private:
 
  Int_t   fNevents;         // Number of events to digitize
  Int_t   fNdigits;         // Number of digits
  TString fDigitsFile ;     // output file   
  TString fHeadersFile;     // input file
  
  Float_t fPhotoCathodeEfficiency; // Photocathode efficiency
  Float_t fPMVoltage ;             // Photomultiplier voltage
  Float_t fPMGain;                 // Photomultiplier gain

 protected:
 
  AliRunLoader *fRunLoader;  // Pointer to Run Loader
  AliVZEROhit  *fVZEROHit;   // Pointer to specific detector hits
  AliDetector  *fVZERO;      // Get pointers to Alice detectors 
                             // and Hit containers 
  AliLoader    *fVZEROLoader;  // Pointer to specific detector loader

  TClonesArray *fHits;       // Pointer to hit array
  TParticle    *fParticle;   // Pointer to a given particle

  TTree        *fTreeH;      // Hits tree
  TTree        *fTreeD;      // Digits tree

  TClonesArray *fDigits;     // List of digits

  ClassDef(AliVZERODigitizer,1) 

};

#endif // AliVZERODigitizer_H
