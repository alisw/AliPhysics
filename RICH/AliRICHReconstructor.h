#ifndef AliRICHReconstructor_h
#define AliRICHReconstructor_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliReconstructor.h>
#include "AliRICHTracker.h"
#include "AliRICHClusterFinder.h"

class AliRICHReconstructor: public AliReconstructor 
{
public:
           AliRICHReconstructor(): AliReconstructor()              {}//default ctor
  virtual ~AliRICHReconstructor()                                  {}//dtor  
  virtual AliTracker*  CreateTracker(AliRunLoader*           )const{return new AliRICHTracker;}                         //interface from AliReconstructor
  virtual void         Reconstruct  (AliRunLoader* pRunLoader)const{AliRICHClusterFinder clus(pRunLoader); clus.Exec();}//interface from AliReconstructor
  using AliReconstructor::Reconstruct;                                                                                  //to get rid of virtual hidden warning 
//  virtual void         FillESD(AliRunLoader* pAL, AliESD* pESD) const;    //virtual from AliReconstructor
//  using AliReconstructor::FillESD;                                        //to get rid of virtual hidden warning 
protected:
  ClassDef(AliRICHReconstructor, 0)   //class for the RICH reconstruction
};

#endif
