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
  AliVZEROBuffer(const AliVZEROBuffer &source); // copy constructor
  AliVZEROBuffer& operator=(const AliVZEROBuffer &source); // ass. op.
  void    WriteTriggerInfo(UInt_t trigger);
  void    WriteTriggerScalers();
  void    WriteBunchNumbers();  
  void    WriteChannel(Int_t cell,UInt_t ADC, UInt_t Time);
  void    WriteBeamFlags();
  void    WriteMBInfo();
  void    WriteMBFlags();  
  void    WriteBeamScalers();
  void    WriteTiming(Int_t cell,UInt_t ADC, UInt_t Time);

  void    SetVerbose(Int_t val){fVerbose=val;}
  Int_t   GetVerbose() const{return  fVerbose;} 
  
private:
  Int_t fVerbose; //Verbosity level: 0-silent, 1:cout msg, 2: txt files for checking
  AliFstream* f;      //The IO file name
  ClassDef(AliVZEROBuffer,1)
};

#endif
