#ifndef ALIITSURECONSTRUCTOR_H
#define ALIITSURECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSUReconstructor.h 54053 2012-01-22 22:12:15Z masera $ */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ITS reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliITSURecoParam.h"

class AliTracker;
class AliTrackleter;
class AliITSUGeomTGeo;

class AliITSUReconstructor: public AliReconstructor {
public:
  AliITSUReconstructor();
  virtual ~AliITSUReconstructor();
  virtual void         Init();
  virtual void         Reconstruct(TTree *digitsTree, TTree *clustersTree) const;
  virtual void         Reconstruct(AliRawReader*, TTree*) const {};

  virtual AliTracker*    CreateTracker() const;
  virtual AliVertexer*   CreateVertexer() const;
  virtual AliTrackleter* CreateMultFinder() const;
  virtual AliTracker*    CreateTrackleter() const;

  virtual const char*    GetDetectorName() const {return "ITS";}

  static const AliITSURecoParam* GetRecoParam() { 
    return dynamic_cast<const AliITSURecoParam*>(AliReconstructor::GetRecoParam(0)); }

private:
  AliITSUReconstructor(const AliITSUReconstructor &); //Not implemented
  AliITSUReconstructor& operator=(const AliITSUReconstructor &); //Not implemented

  AliITSUGeomTGeo* fGM;          // geometry wrapper
  TObjArray        fClusterFinders; // array of clusterfinders per layer
  TObjArray        fRecPoints;      // container for recpoints TClonesArrays
  //
  ClassDef(AliITSUReconstructor, 0)   // class for the ITSU reconstruction
};

#endif
