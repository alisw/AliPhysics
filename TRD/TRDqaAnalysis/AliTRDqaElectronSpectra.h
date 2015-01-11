#ifndef ALITRDQAELECTRONSPECTRA_H
#define ALITRDQAELECTRONSPECTRA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliTRDqaElectronSpectra.h  $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
//
// The transverse momentum spectrum is analyzed stack-by-stack
// for all tracks, and for electron tracks. 
// Tracks have to pass quality cuts. 
// Electrons are waighted with the PID LQ
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliAnalysisTask.h"  

class TTree; 
class AliESDEvent; 
class TH1D; 
class TH2D;
class AliExternalTrackParam;

class AliTRDqaElectronSpectra : public AliAnalysisTask {

 public:

  AliTRDqaElectronSpectra();
  AliTRDqaElectronSpectra(const char *name);
  AliTRDqaElectronSpectra(const AliTRDqaElectronSpectra & trd);
  AliTRDqaElectronSpectra &operator=(const AliTRDqaElectronSpectra & /*g*/) { return *this; };
  virtual ~AliTRDqaElectronSpectra() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

 private:
 
  TTree        * fChain;             //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD;               //! Declaration of leave types

  TObjArray * fOutputContainer; //! output data container
  
  // histograms
  
  TH1D *fStatus;            // track status
  TH1D *fSector;            // sector
  TH1D *fTheta;             // theta to decide on stack
  TH1D *fStack;             // stack ID

  TH1D *fnTracks;           // number of tracks in a stack
  TH1D *fnElTracks;         // number of electrons tracks in a stack
  TH1D *fTracksRatio;       // fraction of electron tracks in a stack
  
  TH1D *fPt;                // transverse momentum distribution
  TH1D *fPtElectron;        // transverse momentum of electrons

  TH1D *fMeanPt;            // all tracks
  TH1D *fMeanPtElectron;    // electrons

  TH2D *fPtStack;           // pt distribution per stack
  TH2D *fPtStackElectron;   // for electrons

  TH1D *fElectronLQ;        // electron likehood
  

  ClassDef(AliTRDqaElectronSpectra, 0); // a TRD analysis task 
};
#endif // ALITRDQAELECTRONSPECTRA_H
