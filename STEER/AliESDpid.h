#ifndef ALIESDPID_H
#define ALIESDPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    Combined PID class
//           for the Event Summary Data class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------
#include <Rtypes.h>

class AliESDEvent;

class AliESDpid {
public:
  AliESDpid(){}
  virtual ~AliESDpid() {}
  static Int_t MakePID(AliESDEvent *event);
private:
  ClassDef(AliESDpid,2)   // TPC PID class
};

#endif


