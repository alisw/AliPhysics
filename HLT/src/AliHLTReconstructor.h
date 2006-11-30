// @(#) $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifdef use_reconstruction
#include "AliReconstructor.h"

class AliITSgeom;

class AliHLTReconstructor: public AliReconstructor {
public:
  AliHLTReconstructor();
  AliHLTReconstructor(Bool_t doTracker, Bool_t doHough);
  virtual ~AliHLTReconstructor();

  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const{
    AliReconstructor::Reconstruct(digitsTree,clustersTree);
  }
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const {
    AliReconstructor::Reconstruct(rawReader,clustersTree);
  }
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         Reconstruct(AliRunLoader* runLoader, 
				   AliRawReader* rawReader) const {
    AliReconstructor::Reconstruct(runLoader,rawReader);
  }
  virtual AliTracker*  CreateTracker(AliRunLoader*) const;
  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
			       AliESD* esd) const {
    AliReconstructor::FillESD(digitsTree,clustersTree,esd);
  }
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESD* esd) const {
    AliReconstructor::FillESD(rawReader,clustersTree,esd);
  }
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  virtual void         FillESD(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESD* esd) const {
    AliReconstructor:: FillESD(runLoader,rawReader,esd);
  }
  void SetDoBench(Bool_t b){fDoBench=b;}
  void SetDoCleanup(Bool_t b){fDoCleanUp=b;}

private:
  void ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const;
  void ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const;
  void FillESDforConformalMapping(AliESD* esd,Int_t iEvent) const;
  void FillESDforHoughTransform(AliESD* esd,Int_t iEvent) const;

  AliITSgeom*          GetITSgeom(AliRunLoader* runLoader) const;

  Bool_t fDoHough;   //do the hough transform
  Bool_t fDoTracker; //do the standard conformal tracker
  Bool_t fDoBench;   //store the benchmark results
  Bool_t fDoCleanUp; //delete tmp tracking files

  ClassDef(AliHLTReconstructor, 0)   // class for the TPC reconstruction
};
#endif

typedef AliHLTReconstructor AliL3Reconstructor; // for backward compatibility

#endif
