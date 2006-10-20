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

class AliITSgeom;
class AliLoader;
class AliTracker;
class AliITStrackerMI;
class AliITSpidESD;

class AliITSReconstructor: public AliReconstructor {
public:
  AliITSReconstructor();
  virtual ~AliITSReconstructor();
  AliITSReconstructor(const AliITSReconstructor &ob); // copy constructor
  AliITSReconstructor& operator=(const AliITSReconstructor & ob); // ass. op.
  virtual void Init(AliRunLoader* runLoader) const;
  
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         Reconstruct(AliRunLoader* runLoader,
				   AliRawReader* rawReader) const;
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const 
    {AliReconstructor::Reconstruct(rawReader,clustersTree);}
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const 
    {AliReconstructor::Reconstruct(digitsTree, clustersTree);}

  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual AliVertexer* CreateVertexer(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
				 AliESD* esd) const 
    {AliReconstructor::FillESD(digitsTree, clustersTree, esd);}
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESD* esd) const
    {AliReconstructor::FillESD(rawReader, clustersTree, esd);}
  virtual void         FillESD(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESD* esd) const
    {AliReconstructor::FillESD(runLoader,rawReader, esd);}

private:
  // methods
  AliITSgeom*          GetITSgeom(AliRunLoader* runLoader) const;
  //data
  AliITSpidESD *fItsPID; //Pid for ITS
  ClassDef(AliITSReconstructor, 1)   // class for the ITS reconstruction
};

#endif
