// @(#) $Id$

#ifndef ALIHLTRECONSTRUCTOR_H
#define ALIHLTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifdef use_reconstruction
#include "AliReconstructor.h"

class AliHLTSystem;

class AliHLTReconstructor: public AliReconstructor {
public:
  AliHLTReconstructor();
  AliHLTReconstructor(Bool_t doTracker, Bool_t doHough);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTReconstructor(const AliHLTReconstructor& src);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTReconstructor& operator=(const AliHLTReconstructor& src);
  /** destructor */
  virtual ~AliHLTReconstructor();

  /** init the reconstructor */
  void Init(AliRunLoader* runLoader);

  /** reconstruct simulated MC data */
  void Reconstruct(AliRunLoader* runLoader) const;
  /** reconstruct data from RawReader */
  void Reconstruct(AliRunLoader* runLoader, AliRawReader* rawReader) const;

  /** create a tracker */
  AliTracker*  CreateTracker(AliRunLoader*) const;

  /** fill esd for one event */
  void FillESD(AliRunLoader* runLoader, AliESD* esd) const;

  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const{
    AliReconstructor::Reconstruct(digitsTree,clustersTree);
  }
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const {
    AliReconstructor::Reconstruct(rawReader,clustersTree);
  }

  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
			       AliESD* esd) const {
    AliReconstructor::FillESD(digitsTree,clustersTree,esd);
  }
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESD* esd) const {
    AliReconstructor::FillESD(rawReader,clustersTree,esd);
  }
  virtual void         FillESD(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESD* esd) const {
    AliReconstructor:: FillESD(runLoader,rawReader,esd);
  }
  void SetDoBench(Bool_t b){fDoBench=b;}
  void SetDoCleanup(Bool_t b){fDoCleanUp=b;}
  virtual void         FillDHLTRecPoint(AliRawReader* rawReader, Int_t nofEvent, Int_t dcCut) const;
private:
  void ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const;
  void ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const;
  void FillESDforConformalMapping(AliESD* esd,Int_t iEvent) const;
  void FillESDforHoughTransform(AliESD* esd,Int_t iEvent) const;

  Bool_t fDoHough;   //do the hough transform
  Bool_t fDoTracker; //do the standard conformal tracker
  Bool_t fDoBench;   //store the benchmark results
  Bool_t fDoCleanUp; //delete tmp tracking files

  AliHLTSystem* fpSystem; //! HLT steering object
  Int_t  fRecEvents;      //! number of reconstructed events
  Int_t  fFilled;         //! number of event filled to ESD

  ClassDef(AliHLTReconstructor, 1)   // class for the TPC reconstruction
};
#endif

typedef AliHLTReconstructor AliL3Reconstructor; // for backward compatibility

#endif
