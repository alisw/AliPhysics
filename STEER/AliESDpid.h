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

class AliESD;

class AliESDpid {
public:
  AliESDpid(){}
  static Int_t MakePID(AliESD *event);
private:
  ClassDef(AliESDpid,1)   // TPC PID class
};

#endif


