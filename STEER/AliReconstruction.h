#ifndef ALIRECONSTRUCTION_H
#define ALIRECONSTRUCTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>
#include "AliReconstructor.h"
#include "AliDetector.h"

class AliRunLoader;
class AliLoader;
class AliVertexer;
class AliTracker;
class AliESD;
class TFile;


class AliReconstruction: public TNamed {
public:
  AliReconstruction(const char* gAliceFilename = "galice.root",
		    const char* name = "AliReconstruction", 
		    const char* title = "reconstruction");
  AliReconstruction(const AliReconstruction& rec);
  AliReconstruction& operator = (const AliReconstruction& rec);
  virtual ~AliReconstruction();

  void           SetGAliceFile(const char* fileName);

  void           SetRunLocalReconstruction(const char* detectors) {
    fRunLocalReconstruction = detectors;};
  void           SetRunVertexFinder(Bool_t run) {fRunVertexFinder = run;};
  void           SetRunTracking(Bool_t run) {fRunTracking = run;};
  void           SetFillESD(const char* detectors) {fFillESD = detectors;};

  void           SetStopOnError(Bool_t stopOnError) 
    {fStopOnError = stopOnError;}
  void           SetCheckPointLevel(Int_t checkPointLevel)
    {fCheckPointLevel = checkPointLevel;}

  virtual Bool_t Run();

private:
  class AliDummyReconstructor: public AliReconstructor {
  public:
    AliDummyReconstructor(AliDetector* detector) {fDetector = detector;};
    virtual ~AliDummyReconstructor() {};

    virtual void         Reconstruct(AliRunLoader* /*runLoader*/) const
      {fDetector->Reconstruct();};
    virtual AliVertexer* CreateVertexer(AliRunLoader* /*runLoader*/) const 
      {return fDetector->CreateVertexer();}
    virtual AliTracker*  CreateTracker(AliRunLoader* /*runLoader*/) const 
      {return fDetector->CreateTracker();}
    virtual void         FillESD(AliRunLoader* /*runLoader*/, AliESD* esd) const
      {fDetector->FillESD(esd);};

    virtual const char*  GetDetectorName() const
      {return fDetector->GetName();};
  private:
    AliDetector*         fDetector;   // detector object
  };

  Bool_t         RunLocalReconstruction(const TString& detectors);
  Bool_t         RunVertexFinder(AliESD*& esd);
  Bool_t         RunTracking(AliESD*& esd);
  Bool_t         FillESD(AliESD*& esd, const TString& detectors);

  Bool_t         IsSelected(TString detName, TString& detectors) const;
  AliReconstructor* GetReconstructor(const char* detName) const;
  Bool_t         CreateVertexer();
  Bool_t         CreateTrackers();
  void           CleanUp(TFile* file = NULL);

  Bool_t         ReadESD(AliESD*& esd, const char* recStep) const;
  void           WriteESD(AliESD* esd, const char* recStep) const;

  TString        fRunLocalReconstruction; // run the local reconstruction for these detectors
  Bool_t         fRunVertexFinder;    // run the vertex finder
  Bool_t         fRunTracking;        // run the barrel tracking
  TString        fFillESD;            // fill ESD for these detectors
  TString        fGAliceFileName;     // name of the galice file
  Bool_t         fStopOnError;        // stop or continue on errors
  Int_t          fCheckPointLevel;    // level of ESD check points

  AliRunLoader*  fRunLoader;          //! current run loader object
  AliLoader*     fITSLoader;          //! loader for ITS
  AliVertexer*   fITSVertexer;        //! vertexer for ITS
  AliTracker*    fITSTracker;         //! tracker for ITS
  AliLoader*     fTPCLoader;          //! loader for TPC
  AliTracker*    fTPCTracker;         //! tracker for TPC
  AliLoader*     fTRDLoader;          //! loader for TRD
  AliTracker*    fTRDTracker;         //! tracker for TRD
  AliLoader*     fTOFLoader;          //! loader for TOF
  AliTracker*    fTOFTracker;         //! tracker for TOF

  static const Int_t fgkNDetectors = 14;   //! number of detectors
  static const char* fgkDetectorName[fgkNDetectors]; //! names of detectors
  TObjArray      fReconstructors;     //! array of reconstructor objects

  ClassDef(AliReconstruction, 1)      // class for running the reconstruction
};

#endif
