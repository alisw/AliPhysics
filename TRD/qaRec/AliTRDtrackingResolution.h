#ifndef ALITRDTRACKINGRESOLUTION_H
#define ALITRDTRACKINGRESOLUTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingResolution.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"

class TObjArray;
class TList;
class TH2I;
class TProfile;
class TTreeSRedirector;
class AliTRDReconstructor;
class AliTRDrecoParam;
class AliTRDtrackingResolution : public AliAnalysisTask{
public:
  AliTRDtrackingResolution(const char *name = "TRD Tracking Resolution");
  ~AliTRDtrackingResolution(){};
  
  void ConnectInputData(Option_t *);
  void CreateOutputObjects();
  Int_t GetDebugLevel() const { return fDebugLevel;}
  void Exec(Option_t *);
  void SetDebugLevel(Int_t level);
  void SetRecoParam(AliTRDrecoParam *r);
  void Terminate(Option_t *);
  
private:
  AliTRDtrackingResolution(const AliTRDtrackingResolution&);
  AliTRDtrackingResolution& operator=(const AliTRDtrackingResolution&);

private:
  enum{
    kNLayers = 6
  };
  TObjArray *fTrackObjects;     // Input Track Info container
  TList     *fOutputHistograms; // Container for the output histograms
  TH2I      *fYRes;
  TH2I      *fPhiRes;
  
  AliTRDReconstructor   *fReconstructor;  //! local reconstructor
  Int_t fDebugLevel;											// Debug Level
  TTreeSRedirector *fDebugStream; 				// Debug stream
  
  ClassDef(AliTRDtrackingResolution, 1) // tracking resolution task
};
#endif
