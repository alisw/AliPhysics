#ifndef AliRICHReconstructor_h
#define AliRICHReconstructor_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliReconstructor.h>
#include <AliLog.h>
#include "AliRICHTracker.h"

class AliRICH;

class AliRICHReconstructor: public AliReconstructor 
{
public:
           AliRICHReconstructor(): AliReconstructor()          {AliDebug(1,"Start.");}
  virtual ~AliRICHReconstructor()                              {AliDebug(1,"Start.");}  
  virtual AliTracker*  CreateTracker(AliRunLoader*)const       {AliDebug(1,"Start.");return new AliRICHTracker;}    //virtual from AliReconstructor
  virtual void         Reconstruct(AliRunLoader* pAL) const;              //virtual from AliReconstructor
  virtual void         FillESD(AliRunLoader* pAL, AliESD* pESD) const;    //virtual from AliReconstructor
private:
  AliRICH*             GetRICH(AliRunLoader* pAL) const;

  ClassDef(AliRICHReconstructor, 0)   //class for the RICH reconstruction
};

#endif
