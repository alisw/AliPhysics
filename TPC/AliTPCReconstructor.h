#ifndef ALITPCRECONSTRUCTOR_H
#define ALITPCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"

class AliTPCParam;


class AliTPCReconstructor: public AliReconstructor {
public:
  AliTPCReconstructor(): AliReconstructor() {};
  virtual ~AliTPCReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

  static void SetCtgRange(Double_t ctgRange = 1.05) {fgCtgRange = ctgRange;}
  static Double_t GetCtgRange(){ return fgCtgRange;}

private:
  AliTPCParam*         GetTPCParam(AliRunLoader* runLoader) const;

  static Double_t fgCtgRange; //! +-fCtgRange is the ctg(Theta) window used for clusterization and tracking (MI) 

  ClassDef(AliTPCReconstructor, 0)   // class for the TPC reconstruction
};

#endif
