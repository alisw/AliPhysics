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

class AliITSgeom;
class AliTracker;
class AliITStrackerMI;
class AliITSpidESD;
class AliITSDetTypeRec;

class AliITSReconstructor: public AliReconstructor {
public:
  AliITSReconstructor();
  virtual ~AliITSReconstructor();
  AliITSReconstructor(const AliITSReconstructor &ob); // copy constructor
  AliITSReconstructor& operator=(const AliITSReconstructor & ob); // ass. op.
  virtual void         Init();
  
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual AliTracker*  CreateTracker() const;
  virtual AliVertexer* CreateVertexer() const;

  virtual void         FillESD(TTree* /*digitsTree*/, TTree* clustersTree, 
			       AliESDEvent* esd) const; 
  virtual void         FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, 
			       AliESDEvent* esd) const
  {FillESD((TTree*)NULL, clustersTree, esd);}

  void SetRecoParam(AliITSRecoParam * param){ fgkRecoParam = param;}
  static const AliITSRecoParam* GetRecoParam(){ return fgkRecoParam;}

private:
  //data
  static AliITSRecoParam *fgkRecoParam; // reconstruction parameters
  AliITSpidESD           *fItsPID;      // Pid for ITS
  AliITSDetTypeRec       *fDetTypeRec;  // reconstructor

  ClassDef(AliITSReconstructor, 3)   // class for the ITS reconstruction
};

#endif
