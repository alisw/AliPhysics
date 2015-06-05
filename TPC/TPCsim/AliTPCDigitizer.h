#ifndef ALITPCDIGITIZER_H
#define ALITPCDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigitizer.h"
class TTreeSRedirector;

class AliDigitizationInput;
class AliTPCSAMPAEmulator;

class AliTPCDigitizer : public AliDigitizer {
 public:    
  enum EStreamFlags{
    kStreamCrosstalk        =0x1,     // flag: stream crosstalk signal ()
    kStreamSignal           =0x2,     // flag: stream signal per pad
    kStreamSignalAll        =0x4      // flag: stream signal per pad dump all signal (without 0 suppression)
  };
  AliTPCDigitizer();
  AliTPCDigitizer(AliDigitizationInput * digInput);
  virtual ~AliTPCDigitizer();
    // Initialize merging and digitization
    virtual Bool_t Init();
    // Do the main work
    virtual void Digitize(Option_t* option=0);    
    Int_t GetDebug() const {return fDebug;}       // get debug level
    void SetDebug(Int_t level){fDebug = level;}   // set debug level     
  static AliTPCSAMPAEmulator *GetEmulator(){return fgSAMPAEmulator;}
  static void SetEmulator( AliTPCSAMPAEmulator *emulator){fgSAMPAEmulator=emulator;}
 private: 
    void DigitizeFast(Option_t* option=0); //digitize - using row pointers
    void DigitizeSave(Option_t* option=0); // digitize using controlled arrays   
    void DigitizeWithTailAndCrossTalk(Option_t* option=0); 
    Int_t fDebug;                         //
    TTreeSRedirector *fDebugStreamer;     //!debug streamer
  static AliTPCSAMPAEmulator *fgSAMPAEmulator; 
 private:
    AliTPCDigitizer& operator=(const AliTPCDigitizer&);
    AliTPCDigitizer(const AliTPCDigitizer&);
    ClassDef(AliTPCDigitizer,2)  // MUON merging/digitization
};    
#endif

