#ifndef ALIITSRECONSTRUCTOR_H
#define ALIITSRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ITS reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliITSRecoParam.h"
#include "AliITSDetTypeRec.h"
class AliESDpid;
class AliITSgeom;
class AliTracker;
class AliITStrackerMI;

class AliITSReconstructor: public AliReconstructor {
public:
  AliITSReconstructor();
  virtual ~AliITSReconstructor();
  virtual void         Init();
  virtual void         GetPidSettings(AliESDpid *ESDpid);
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual AliTracker*  CreateTracker() const;
  virtual AliVertexer* CreateVertexer() const;
  virtual AliTracker*  CreateTrackleter() const;

  virtual void         FillESD(TTree* /*digitsTree*/, TTree* clustersTree, 
			       AliESDEvent* /* esd */) const; 
  virtual void         FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, 
			       AliESDEvent* esd) const
  {FillESD((TTree*)NULL, clustersTree, esd);}
  void ResetRecPoints() {fDetTypeRec->ResetRecPoints();}

  static const AliITSRecoParam* GetRecoParam() { return dynamic_cast<const AliITSRecoParam*>(AliReconstructor::GetRecoParam(0)); }

private:
  AliITSReconstructor(const AliITSReconstructor &); //Not implemented
  AliITSReconstructor& operator=(const AliITSReconstructor &); //Not implemented

  //data
  AliITSDetTypeRec       *fDetTypeRec;  // reconstructor

  ClassDef(AliITSReconstructor, 6)   // class for the ITS reconstruction
};

#endif
