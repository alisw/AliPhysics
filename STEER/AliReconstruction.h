#ifndef ALIRECONSTRUCTION_H
#define ALIRECONSTRUCTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TString.h>

class AliRunLoader;
class AliLoader;
class AliTracker;
class AliESD;


class AliReconstruction: public TNamed {
public:
  AliReconstruction(const char* name = "AliReconstruction", 
		    const char* title = "reconstruction");
  AliReconstruction(const AliReconstruction& rec);
  AliReconstruction& operator = (const AliReconstruction& rec);
  virtual ~AliReconstruction();

  void           SetGAliceFile(const char* fileName);

  void           SetRunReconstruction(const char* detectors) {
    fRunReconstruction = detectors;};
  void           SetRunTracking(Bool_t run) {fRunTracking = run;};
  void           SetFillESD(const char* detectors) {fFillESD = detectors;};

  virtual Bool_t Run();

private:
  void           Init();

  Bool_t         IsSelected(TString detName, TString& detectors) const;

  Bool_t         RunReconstruction(const TString& detectors);
  Bool_t         RunTracking(AliESD* esd);
  Bool_t         FillESD(AliESD* esd, const TString& detectors);

  TString        fRunReconstruction;  // run the reconstr. for these detectors
  Bool_t         fRunTracking;        // run the barrel tracking
  TString        fFillESD;            // fill ESD for these detectors
  Bool_t         fStopOnError;        // stop or continue on errors

  TString        fGAliceFileName;     // name of the galice file

  AliRunLoader*  fRunLoader;          //! current run loader object
  AliLoader*     fITSLoader;          //! loader for ITS
  AliTracker*    fITSTracker;         //! tracker for ITS
  AliLoader*     fTPCLoader;          //! loader for TPC
  AliTracker*    fTPCTracker;         //! tracker for TPC
  AliLoader*     fTRDLoader;          //! loader for TRD
  AliTracker*    fTRDTracker;         //! tracker for TRD
  AliLoader*     fTOFLoader;          //! loader for TOF
  AliTracker*    fTOFTracker;         //! tracker for TOF

  ClassDef(AliReconstruction, 1)      // class for running the reconstruction
};

#endif
