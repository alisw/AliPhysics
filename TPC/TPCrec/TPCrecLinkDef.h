#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliTPCclusterer+;      // The TPC clusterer

#pragma link C++ class AliTPCtrack+;          // Derived from AliTrack base class for TPC tracks
#pragma link C++ class AliTPCpolyTrack+;      // Polynomial description of track (used in AliTPCtracker::MakeSeeds2)
                                              //  working in global coordinate frame
                                              // --- docu to be added
#pragma link C++ class AliTPCseed+;           // Derived from AliTPCtrack - the track seed

#pragma link C++ class AliTPCtrackerRow+;     // Container for info (cluster) on padrow level, method FindNearest ...
#pragma link C++ class AliTPCtrackerSector+;  // Container for info (cluster) on sector level (array rows)
#pragma link C++ class AliTPCtracker+;        // The TPC tracker

#pragma link C++ class AliTPCReconstructor+;  // The TPC reconstructor steering TPC reconstruction
#pragma link C++ class AliTPCTracklet+;       // Used inside calbration for global fitting
                                              // --- should be removed at a later point after calib reassessment

// Used in Krypton --- Update documentation for all 4 classes
#pragma link C++ class AliTPCvtpr+;           // Helper class for clusterer --- Rename such that is clear that it is used in Kr
#pragma link C++ class AliPadMax+;            // Helper class for clusterer --- Rename such that is clear that it is used in Kr
#pragma link C++ class AliTPCclusterKr+;      // Krypton cluster
#pragma link C++ class AliTPCclustererKr+;    // The Krypton clusterer

// Used for Cosmics
#pragma link C++ class AliTPCCosmicUtils+;    // Helper class for cosmic tracker
#pragma link C++ class AliTPCCosmicTrackfit+; // Helper class for cosmic tracker
#pragma link C++ class AliCosmicTracker+;     // Tracker for cosmics (combined fit for upper and lower half)

#pragma read \
  sourceClass="AliTPCseed" \
  targetClass="AliTPCseed" \
  source="Bool_t  fInDead;Bool_t fIsSeeding;Bool_t fBSigned; AliTPCTrackerPoint  fTrackPoints[160]" \
  version="[-8]"	   \
  target="fTrackPointsArr" \
  targetType="AliTPCTrackerPoints" \
  code="{\
  newObj->SetInDead(onfile.fInDead);		\
  newObj->SetIsSeeding(onfile.fIsSeeding);	\
  newObj->SetBSigned(onfile.fBSigned);           \
  for (int i=159;i--;) fTrackPointsArr.SetPoint(i,&onfile.fTrackPoints[i]);}"

#endif
