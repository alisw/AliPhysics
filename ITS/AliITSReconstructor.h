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
  virtual void Init(AliRunLoader* runLoader);
  
  virtual Bool_t       HasLocalReconstruction() const {return kTRUE;};

  virtual void         Reconstruct(AliRunLoader* runLoader) const
    {AliReconstructor::Reconstruct(runLoader);}
  virtual void         Reconstruct(AliRunLoader* runLoader,
				   AliRawReader* rawReader) const
    {AliReconstructor::Reconstruct(runLoader,rawReader);}
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual AliVertexer* CreateVertexer(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESDEvent* esd) const;

  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
				 AliESDEvent* esd) const 
    {AliReconstructor::FillESD(digitsTree, clustersTree, esd);}
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESDEvent* esd) const
    {AliReconstructor::FillESD(rawReader, clustersTree, esd);}
  virtual void         FillESD(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESDEvent* esd) const
    {AliReconstructor::FillESD(runLoader,rawReader, esd);}

  void SetRecoParam(AliITSRecoParam * param){ fgkRecoParam = param;}
  static const AliITSRecoParam* GetRecoParam(){ return fgkRecoParam;}

private:
  //data
  static AliITSRecoParam *fgkRecoParam; // reconstruction parameters
  AliITSpidESD           *fItsPID;      // Pid for ITS
  AliITSDetTypeRec       *fDetTypeRec;  // reconstructor
  ClassDef(AliITSReconstructor, 2)   // class for the ITS reconstruction
};

#endif
