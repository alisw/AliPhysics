//   Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
//  See cxx source for full Copyright notice                               
//  AliFMDReconstructor.h 
//  Task Class for making TreeR in FMD                        
//-- Authors: Evgeny Karpechev (INR) and Alla Maevskaia (INR)
/*
    Reconstruct nember of particles 
    in given group of pads for given FMDvolume
    determine by numberOfVolume , 
    numberOfMinSector,numberOfMaxSector,
    numberOfMinRing, numberOfMaxRing
    Reconstruction method choose dependence on number of empty pads  
  */
/* $Id$ */


#ifndef ALIFMDRECONSTRUCTOR_H
#define ALIFMDRECONSTRUCTOR_H

#include "AliReconstructor.h"


class AliFMDReconstructor: public AliReconstructor 
{
 public:
  AliFMDReconstructor(): AliReconstructor() {}; 
  virtual ~AliFMDReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  
  ClassDef(AliFMDReconstructor, 0)  // class for the FMD reconstruction


}; 
#endif









