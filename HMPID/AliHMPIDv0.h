#ifndef AliHMPIDv0_h
#define AliHMPIDv0_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHMPID.h"

class AliHMPIDv0 : public AliHMPID 
{    
public:
                 AliHMPIDv0():AliHMPID()                                              {;}       //default ctor
                 AliHMPIDv0(const char *name, const char *title):AliHMPID(name,title) {;}       //named ctor
  virtual       ~AliHMPIDv0()                                         {;}                      //dtor
  virtual void   Init()                                              {;}                      //interface from AliHMPID
  virtual Int_t  IsVersion()                                    const{return 0;}              //interface from AliHMPID
  virtual void   StepManager();                                                               //interface from AliHMPID
protected:
  ClassDef(AliHMPIDv0,1)  //HMPID coarse version for material budget study and debuging
};
	
#endif
