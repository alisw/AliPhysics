#ifndef CPVv0_H
#define CPVv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:CPV version 0    //
//  Coarse geometry                                   //
//                                                    //
//  Author: Yuri Kharlov, IHEP, Protvino              //
//  e-mail: Yuri.Kharlov@cern.ch                      //
//  Last modified: 17 September 1999                  //
////////////////////////////////////////////////////////

// --- galice header files ---
#include "AliCPV.h"

class AliCPVv0 : public AliCPV
{
  
public:
           AliCPVv0();
           AliCPVv0(const char *name, const char *title);
  virtual ~AliCPVv0(){}
  virtual void   CreateGeometry();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  
  ClassDef(AliCPVv0,1)  //Hits manager for set:CPV version 0
};
 
#endif

