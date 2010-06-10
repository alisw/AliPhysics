#ifndef AliVZEROBUFFER_H
#define AliVZEROBUFFER_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
// Class used for storing VZERO digits according to the DDLs format//
/////////////////////////////////////////////////////////////////////

#ifdef __CINT__
class fstream;
#else
#include "Riostream.h"
#endif

#include "AliFstream.h"

class AliVZEROBuffer:public TObject{

public:
  AliVZEROBuffer();
  AliVZEROBuffer(const char* fileName); //constructor
  virtual ~AliVZEROBuffer(); //destructor
  void    WriteTriggerInfo(UInt_t trigger);
  void    WriteTriggerScalers();
  void    WriteBunchNumbers();  
  void    WriteChannel(Int_t channel, Short_t *adc, Bool_t integrator);
  void    WriteBeamFlags(Bool_t *bbFlag, Bool_t *bgFlag);
  void    WriteMBInfo();
  void    WriteMBFlags();  
  void    WriteBeamScalers();
  void    WriteTiming(Float_t time, Float_t width);

private:
  AliVZEROBuffer(const AliVZEROBuffer &source); // copy constructor
  AliVZEROBuffer& operator=(const AliVZEROBuffer &source); // ass. op.

  UInt_t      fRemainingWord; // Remaining data word between even and odd channel's data
  AliFstream* f;      //The IO file name
  ClassDef(AliVZEROBuffer,2)
};

#endif
