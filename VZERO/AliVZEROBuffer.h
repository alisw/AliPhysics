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
  void    WriteChannel(Int_t cell,UInt_t ADC, Float_t Time, Bool_t integrator);
  void    WriteBeamFlags();
  void    WriteMBInfo();
  void    WriteMBFlags();  
  void    WriteBeamScalers();
  void    WriteTiming(Int_t cell,UInt_t ADC, Float_t Time);

  void    SetVerbose(Int_t val){fVerbose=val;}
  Int_t   GetVerbose() const{return  fVerbose;} 
  
// Getter of Online channel number (as defined in the FEE readout) from the 
// Offline channel (as defined by aliroot numbering convention):
  Int_t              GetOnlineChannel(Int_t channel)  const
       { Int_t  fOnlineChannel[64] = {39, 38, 37, 36, 35, 34, 33, 32, 
                                      47, 46, 45, 44, 43, 42, 41, 40, 
			              55, 54, 53, 52, 51, 50, 49, 48, 
			              63, 62, 61, 60, 59, 58, 57, 56,
			               7,  6,  5,  4,  3,  2,  1,  0, 
			              15, 14, 13, 12, 11, 10,  9,  8,
			              23, 22, 21, 20, 19, 18, 17, 16, 
			              31, 30, 29, 28, 27, 26, 25, 24};
               return fOnlineChannel[channel]; }
	       
private:
  Int_t fVerbose; //Verbosity level: 0-silent, 1:cout msg, 2: txt files for checking
  AliFstream* f;      //The IO file name
  ClassDef(AliVZEROBuffer,1)
};

#endif
