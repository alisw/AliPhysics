#ifndef ALITOFRECONSTRUCTOR_H
#define ALITOFRECONSTRUCTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"

class TTree;

class AliESDEvent;
class AliRawReader;

class AliTOFcalib;

class AliTOFReconstructor: public AliReconstructor {
public:
  AliTOFReconstructor();
  AliTOFReconstructor(const AliTOFReconstructor &source); // copy constructor
  AliTOFReconstructor& operator=(const AliTOFReconstructor &source); // ass. op.
  virtual ~AliTOFReconstructor();

  virtual void         Reconstruct(AliRawReader* rawReader,
				   TTree* clusterTree) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clusterTree) const;

  virtual void         ConvertDigits(AliRawReader* reader, TTree* digitsTree) const;

  virtual AliTracker*  CreateTracker() const;

  virtual void         FillESD(AliRawReader*, TTree*clustersTree, AliESDEvent*esd) const
  {FillESD((TTree*)NULL,clustersTree,esd);}
  virtual void         FillESD(TTree*, TTree*, AliESDEvent*) const {}

private:
  AliTOFcalib    *fTOFcalib;    // pointer to TOF calib class

  ClassDef(AliTOFReconstructor, 3)   // class for the TOF reconstruction
};

#endif
