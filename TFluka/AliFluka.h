#ifndef ALIFLUKA_H
#define ALIFLUKA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// FLUKA implementation of the AliMC Interface                               //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TFluka.h"

class AliFluka : public TFluka
{
 public:
    AliFluka(const char *name, const char *title) 
	{TFluka(name, title);}
    AliFluka() {;}
    virtual ~AliFluka() {;}
    virtual void ProcessRun(Int_t nevent);
    virtual void ProcessEvent();
  private:
  AliFluka(const AliFluka &) {}
  AliFluka & operator=(const AliFluka &) {return (*this);}

  ClassDef(AliFluka,1)  //Fluka implementation of AliMC including Fluka specific methods
};

#endif 

