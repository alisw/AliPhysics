#ifndef ALIRECONSTRUCTION_H
#define ALIRECONSTRUCTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the reconstruction                                      //
// Clusters and tracks are created for all detectors and all events by       //
// typing:                                                                   //
//                                                                           //
//   AliReconstruction rec;                                                  //
//   rec.Run();                                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>

class AliReconstructor;
class AliRunLoader;
class AliRawReader;
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
  void           SetInput(const char* input) {fInput = input;};
  void           SetOption(const char* detector, const char* option);

  void           SetRunLocalReconstruction(const char* detectors) {
    fRunLocalReconstruction = detectors;};
  void           SetRunVertexFinder(Bool_t run) {fRunVertexFinder = run;};
  void           SetRunTracking(const char* detectors) {
    fRunTracking = detectors;};
  void           SetFillESD(const char* detectors) {fFillESD = detectors;};
  void           SetRunReconstruction(const char* detectors) {
    SetRunLocalReconstruction(detectors); 
    SetRunTracking(detectors);
    SetFillESD(detectors);};

  void           SetStopOnError(Bool_t stopOnError) 
    {fStopOnError = stopOnError;}
  void           SetCheckPointLevel(Int_t checkPointLevel)
    {fCheckPointLevel = checkPointLevel;}

  virtual Bool_t Run(const char* input = NULL);

private:
  Bool_t         RunLocalReconstruction(const TString& detectors);
  Bool_t         RunVertexFinder(AliESD*& esd);
  Bool_t         RunTracking(AliESD*& esd);
  Bool_t         FillESD(AliESD*& esd, const TString& detectors);

  Bool_t         IsSelected(TString detName, TString& detectors) const;
  AliReconstructor* GetReconstructor(Int_t iDet);
  Bool_t         CreateVertexer();
  Bool_t         CreateTrackers(const TString& detectors);
  void           CleanUp(TFile* file = NULL);

  Bool_t         ReadESD(AliESD*& esd, const char* recStep) const;
  void           WriteESD(AliESD* esd, const char* recStep) const;

  TString        fRunLocalReconstruction; // run the local reconstruction for these detectors
  Bool_t         fRunVertexFinder;    // run the vertex finder
  TString        fRunTracking;        // run the tracking for these detectors
  TString        fFillESD;            // fill ESD for these detectors
  TString        fGAliceFileName;     // name of the galice file
  TString        fInput;              // name of input file or directory
  Bool_t         fStopOnError;        // stop or continue on errors
  Int_t          fCheckPointLevel;    // level of ESD check points
  TObjArray      fOptions;            // options for reconstructor objects

  AliRunLoader*  fRunLoader;          //! current run loader object
  AliRawReader*  fRawReader;          //! current raw data reader

  static const Int_t fgkNDetectors = 15;   //! number of detectors
  static const char* fgkDetectorName[fgkNDetectors]; //! names of detectors
  AliReconstructor*  fReconstructor[fgkNDetectors];  //! array of reconstructor objects
  AliLoader*     fLoader[fgkNDetectors];   //! detector loaders
  AliVertexer*   fVertexer;                //! vertexer for ITS
  AliTracker*    fTracker[fgkNDetectors];  //! trackers

  ClassDef(AliReconstruction, 3)      // class for running the reconstruction
};

#endif
