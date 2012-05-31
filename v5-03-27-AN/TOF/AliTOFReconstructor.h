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
#include "AliTOFRecoParam.h"

class TTree;

class AliESDEvent;
class AliRawReader;
class AliTOFcalib;
//class AliTOFT0maker;

class AliTOFReconstructor: public AliReconstructor {
public:
  AliTOFReconstructor();
  virtual ~AliTOFReconstructor();

  virtual void         Reconstruct(AliRawReader* rawReader,
				   TTree* clusterTree) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clusterTree) const;

  virtual void         ConvertDigits(AliRawReader* reader, TTree* digitsTree) const;

  virtual AliTracker*  CreateTracker() const;

  virtual void         FillESD(AliRawReader*, TTree*clustersTree, AliESDEvent* esd) const
  {FillESD((TTree*)NULL,clustersTree,esd);}
  virtual void         FillESD(TTree *, TTree *, AliESDEvent * /*esdEvent*/) const;

  static const AliTOFRecoParam* GetRecoParam() { return dynamic_cast<const AliTOFRecoParam*>(AliReconstructor::GetRecoParam(3)); } // getting RecoParam obj

  virtual void FillEventTimeWithTOF(AliESDEvent *event, AliESDpid *esdPID);

private:
  AliTOFReconstructor(const AliTOFReconstructor &); //Not implemented
  AliTOFReconstructor& operator=(const AliTOFReconstructor &); //Not implemented

  AliTOFcalib    *fTOFcalib;    // pointer to TOF calib class
  //AliTOFT0maker  *fTOFT0maker;  // pointer to TOF T0 maker class

  ClassDef(AliTOFReconstructor, 3)   // class for the TOF reconstruction
};

#endif
