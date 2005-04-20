#ifndef ALITOFRECONSTRUCTOR_H
#define ALITOFRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"

class AliTOFGeometry;


class AliTOFReconstructor: public AliReconstructor {
public:
  AliTOFReconstructor(): AliReconstructor() {};
  virtual ~AliTOFReconstructor() {};

  virtual void         Reconstruct(AliRunLoader*, AliRawReader*) const { };
  virtual void         Reconstruct(AliRawReader*, TTree*) const { };
  virtual void         Reconstruct(TTree*, TTree*) const { };
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader*, AliRawReader*, AliESD*) const { };
  virtual void         FillESD(AliRawReader*, TTree*, AliESD*) const { };
  virtual void         FillESD(TTree*, TTree*, AliESD*) const { };
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  AliTOFGeometry*      GetTOFGeometry(AliRunLoader* runLoader) const;

  ClassDef(AliTOFReconstructor, 0)   // class for the TOF reconstruction
};

#endif
