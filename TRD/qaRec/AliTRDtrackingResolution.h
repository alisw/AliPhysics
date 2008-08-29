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
class TProfile;
class TTreeSRedirector;
class AliTRDReconstructor;
class AliTRDrecoParam;
class AliTRDseedV1;
class AliTRDtrackInfo;
class AliTRDtrackingResolution : public AliAnalysisTask{
public:
  AliTRDtrackingResolution(const char *name = "TRD Tracking Resolution");
  ~AliTRDtrackingResolution(){};
  
  void    ConnectInputData(Option_t *);
  void    CreateOutputObjects();
  Int_t   GetDebugLevel() const { return fDebugLevel;}
  void    Exec(Option_t *);
  void    SetDebugLevel(Int_t level);
  void    SetRecoParam(AliTRDrecoParam *r);
  void    Terminate(Option_t *);
  
private:
  AliTRDtrackingResolution(const AliTRDtrackingResolution&);
  AliTRDtrackingResolution& operator=(const AliTRDtrackingResolution&);
  Bool_t  Resolution(AliTRDseedV1 *tracklet, AliTRDtrackInfo *info, Float_t &phi);

private:
  enum{
    kNLayers = 6
  };
  TObjArray *fTracks;     // Input Track Info container
  TList     *fHistos;           // Container for the output histograms
  
  AliTRDReconstructor   *fReconstructor;  //! local reconstructor
  Int_t fDebugLevel;											// Debug Level
  TTreeSRedirector *fDebugStream; 				// Debug stream
  
  ClassDef(AliTRDtrackingResolution, 1) // tracking resolution task
};
#endif
