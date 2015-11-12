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
  sourceClass="AliTPCtrack" \
  targetClass="AliTPCtrack" \
  source="Int_t fIndex[200]" \
  version="[-4]" \
  target="fIndex" \
  targetType="Int_t[159]" \
  code="{ for (int i=159;i--;) fIndex[i]=onfile.fIndex[i]; }" 


#pragma read \
  sourceClass="AliTPCseed" \
  targetClass="AliTPCseed" \
  source="Bool_t fBSigned" \
  version="[-8]"	   \
  target="" \
  targetType="UInt_t" \
  code="{ newObj->SetBSigned(onfile.fBSigned); }"

#pragma read \
  sourceClass="AliTPCseed" \
  targetClass="AliTPCseed" \
  source="AliTPCTrackerPoint  fTrackPoints[160];" \
  version="[-8]"	   \
  target="fTrackPointsArr" \
  targetType="AliTPCTrackerPoints" \
  code="{for (int i=159;i--;) fTrackPointsArr.SetPoint(i,&onfile.fTrackPoints[i]);}"

#pragma read \
  sourceClass="AliTPCseed" \
  targetClass="AliTPCseed" \
  source="AliTPCclusterMI* fClusterPointer[160]" \
  version="[-8]"	   \
  target="fClusterPointer" \
  targetType="AliTPCclusterMI**" \
  code="{fClusterPointer = new AliTPCclusterMI*[159]; \
  for (int i=159;i--;) fClusterPointer[i] = onfile.fClusterPointer[i] ? new AliTPCclusterMI(*onfile.fClusterPointer[i]) : 0;}"

#pragma read \
  sourceClass="AliTPCseed" \
  targetClass="AliTPCseed" \
  source="" \
  version="[-8]"	   \
  target="fNClStore" \
  targetType="Int_t" \
  code="{fNClStore = 159;}"

#endif
