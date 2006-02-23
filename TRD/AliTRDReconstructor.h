#ifndef ALITRDRECONSTRUCTOR_H
#define ALITRDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TRD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include "AliReconstructor.h"

class AliTRDparameter;
class AliRawReader;

class AliTRDReconstructor: public AliReconstructor {
public:
  AliTRDReconstructor(): AliReconstructor() {};
  virtual ~AliTRDReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader, AliRawReader* rawReader) const;
  virtual void         Reconstruct(AliRawReader*, TTree*) const { };
  virtual void         Reconstruct(TTree*, TTree*) const { };
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader*, AliRawReader*, AliESD*) const { };
  virtual void         FillESD(AliRawReader*, TTree*, AliESD*) const { };
  virtual void         FillESD(TTree*, TTree*, AliESD*) const { };
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  AliTRDparameter*     GetTRDparameter(AliRunLoader* runLoader) const;

  ClassDef(AliTRDReconstructor, 0)   // class for the TRD reconstruction
};

#endif
