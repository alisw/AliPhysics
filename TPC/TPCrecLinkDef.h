#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliTPCclusterMI+;      // derived from AliTPCCluster (shape in addition)
                                              // -- ask peter what can happen if renamed (add some pragma !??!)
             // if no used  ... remove non MI classes

#pragma link C++ class AliTPCclusterInfo+;    // additional info attach to cluster (add digit map)
                                              // currently not used by default

#pragma link C++ class AliComplexCluster+;    // Used to store additional cluster and tracklet information along track
                                              // -> the following classes 
                                              // --- Documentation to be added - classes to be cleaned
                                              // ... Marian to investigate

#pragma link C++ class AliTPCTrackerPoint+;   // defined in AliComplexCluster.h
#pragma link C++ class AliTPCClusterPoint+;   // defined in AliComplexCluster.h
#pragma link C++ class AliTPCExactPoint+;     // defined in AliComplexCluster.h
#pragma link C++ class AliTPCTrackPoint+;     // defined in AliComplexCluster.h
#pragma link C++ class AliTPCTrackPoint2+;    // defined in AliComplexCluster.h

#pragma link C++ class AliClusters+;          // Generic container for clusters derived from segmentID - all clusters
                                              //  contains 1 AliTPCClustersRow per segment (1 segment = 1 padrow)
#pragma link C++ class AliTPCClustersRow+;    // TPC Container array of "cluster" inside 1 padrow

#pragma link C++ class AliClustersArray+;     // Container of clusters ?!?!
                                              // --- docu to be added
                                              // --- remove if not needed 
#pragma link C++ class AliTPCClustersArray+;  // Container of clusters ?!?!
                                              // --- docu to be added
                                              // --- remove if not needed 

//#pragma link C++ class AliTPCclusterer+;  // --- obsolete remove
#pragma link C++ class AliTPCclustererMI+;    // The TPC clusterer -> rename -"MI"

#pragma link C++ class AliTPCtrack+;          // derived from AliTrack base class for TPC tracks

#pragma link C++ class AliTPCpolyTrack+;      // Polynomial description of track (used in AliTPCtrackerMI::MakeSeeds2)
                                              //  working in global coordinate frame
                                              // --- docu to be added
#pragma link C++ class AliTPCseed+;           // derived from AliTPCtrack - the track seed
#pragma link C++ class AliTPCtrackerRow+;     // container for info (cluster) on padrow level
                                              // functions: FindNearest ...
#pragma link C++ class AliTPCtrackerSector+;  // container for info (cluster) on sector level (array rows)
#pragma link C++ class AliTPCtrackerMI+;      // The TPC tracker -> rename -"MI"

#pragma link C++ class AliTPCReconstructor+;  // the TPC reconstructor steering TPC reconstruction
#pragma link C++ class AliTPCRecoParam+;      // config parameters for reconstruction
#pragma link C++ class AliTPCClusterParam+;   // cluster parametrization
#pragma link C++ class AliTPCTracklet+;       // used inside calbration for global fitting
                                              // -- should be removed at a later point after calib reassessment
#pragma link C++ class AliTPCQADataMakerRec+; // offline QA

// -- Used in Krypton .--- Update docuemntation for all 4
#pragma link C++ class AliTPCvtpr+;           // HElperClass for clyusterer - Rename such that is clear taht it is used in Kr
#pragma link C++ class AliPadMax+;            // Helper clas Rename such that is clear taht it is used in Kr
#pragma link C++ class AliTPCclusterKr+;      // Krypton Cluster
#pragma link C++ class AliTPCclustererKr+;    // The Krypron clsuterer

// -- Used for Cosmics
#pragma link C++ class AliTPCCosmicUtils+;    // HelperClass for cosmics tracker
#pragma link C++ class AliTPCCosmicTrackfit+; // dito
#pragma link C++ class AliCosmicTracker+;     // Tracker fro cosmics (combined fit for upper and lower half)

#endif

		       //#pragma link C++ class AliTPCPid+;            // --- move to attic
		       //#pragma link C++ class AliTPCtrackPid+;       // --- move to attic
