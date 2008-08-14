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

class AliTRDtrackingResolution : public AliAnalysisTask{
public:
  AliTRDtrackingResolution(const char *name = "TRD Tracking Resolution");
  ~AliTRDtrackingResolution(){};
  
  void ConnectInputData(Option_t *);
  void CreateOutputObjects();
  void Exec(Option_t *);
  void Terminate(Option_t *);
  
  void SetDebugLevel(Int_t level);
  Int_t GetDebugLevel() const { return fDebugLevel;}
  
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
  //TGraph      *fSigmaY;             // y-Resolution
/*		TH1F *fZres;														// z-Resolution
    TProfile *fYresAngle;										// y-Resolution dependent on Angle*/
/*		TH1F *fYresLayer[kNLayers];							// y-Resolution histograms for each Layer
    TH1F *fZresLayer[kNLayers];							// z-Resolution histograms for each Layer
    TProfile *fYresLayerAngle[kNLayers];		// y-Resolution histograms for each Layer - Angular Dependence*/
/*		TH1F *fPhiRes;													// Angular resolution in Phi-Direction
    TProfile *fPhiResAngle;									// Phi-resolution dependent on angle*/
/*		TH1F *fPhiResLayer[kNLayers];						// Phi-Resolution histograms for each Layer
    TProfile *fPhiResLayerAngle[kNLayers];	// Phi-resolution histograms for each Layer - Angular Dependence */
    Int_t fDebugLevel;											// Debug Level
    TTreeSRedirector *fDebugStream; 				// Debug stream
  
  ClassDef(AliTRDtrackingResolution, 1) // tracking resolution task
};
#endif
