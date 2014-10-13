#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliTPCclusterMI+;      // Derived from Cluster (shape in addition)
                                              // --- ask Peter what can happen if renamed (add some pragma !??!)
#pragma link C++ class AliTPCclusterInfo+;    // additional info attach to cluster (add digit map)
                                              // currently not used by default

#pragma link C++ class AliComplexCluster+;    // Used to store additional cluster and tracklet information along track
                                              // Following classes are derived
                                              // --- Documentation to be added - classes to be cleaned
                                              // --- Marian to investigate
#pragma link C++ class AliTPCTrackerPoint+;   // defined in AliComplexCluster.h
#pragma link C++ class AliTPCClusterPoint+;   // defined in AliComplexCluster.h
#pragma link C++ class AliTPCExactPoint+;     // defined in AliComplexCluster.h
#pragma link C++ class AliTPCTrackPoint+;     // defined in AliComplexCluster.h
#pragma link C++ class AliTPCTrackPoint2+;    // defined in AliComplexCluster.h

#pragma link C++ class AliClusters+;          // Generic container for clusters derived from segmentID - all clusters
                                              //   contains 1 AliTPCClustersRow per segment (1 segment = 1 padrow)
#pragma link C++ class AliTPCClustersRow+;    // TPC Container array of "cluster" inside 1 padrow

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
#pragma link C++ class AliTPCRecoParam+;      // Config parameters for reconstruction
#pragma link C++ class AliTPCClusterParam+;   // Cluster parametrization
#pragma link C++ class AliTPCTracklet+;       // Used inside calbration for global fitting
                                              // --- should be removed at a later point after calib reassessment
#pragma link C++ class AliTPCQADataMakerRec+; // Offline QA

// Used in Krypton --- Update documentation for all 4 classes
#pragma link C++ class AliTPCvtpr+;           // Helper class for clusterer --- Rename such that is clear that it is used in Kr
#pragma link C++ class AliPadMax+;            // Helper class for clusterer --- Rename such that is clear that it is used in Kr
#pragma link C++ class AliTPCclusterKr+;      // Krypton cluster
#pragma link C++ class AliTPCclustererKr+;    // The Krypton clusterer

// Used for Cosmics
#pragma link C++ class AliTPCCosmicUtils+;    // Helper class for cosmic tracker
#pragma link C++ class AliTPCCosmicTrackfit+; // Helper class for cosmic tracker
#pragma link C++ class AliCosmicTracker+;     // Tracker for cosmics (combined fit for upper and lower half)

#endif
