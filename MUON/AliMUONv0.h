#ifndef ALIMUONV0_H
#define ALIMUONV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:
//  MUON version 0: Coarse geometry, no hits 
/////////////////////////////////////////////////////////
 
#include "AliMUON.h"

class AliMUONv0 : public AliMUON 
{
public:
   AliMUONv0();
   AliMUONv0(const char *name, const char *title);
   virtual  ~AliMUONv0() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 0;}  
private:
   ClassDef(AliMUONv0,1)  // MUON Detector class Version 0 
};
#endif







