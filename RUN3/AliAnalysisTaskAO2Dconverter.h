/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskAO2Dconverter_H
#define AliAnalysisTaskAO2Dconverter_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAnalysisFilter.h"

#include <TString.h>

#include <Rtypes.h>

class AliESDEvent;

class AliAnalysisTaskAO2Dconverter : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskAO2Dconverter() = default;
  AliAnalysisTaskAO2Dconverter(const char *name);
  virtual ~AliAnalysisTaskAO2Dconverter();

  AliAnalysisTaskAO2Dconverter(const AliAnalysisTaskAO2Dconverter &) = default;
  AliAnalysisTaskAO2Dconverter &operator=(const AliAnalysisTaskAO2Dconverter &) = delete;

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetNumberOfEventsPerCluster(int n) { fNumberOfEventsPerCluster = n; }

  static AliAnalysisTaskAO2Dconverter* AddTask(TString suffix = "");
  enum TreeIndex {
    kEvents = 0,
    kTracks,
    kCalo,
    kTOF,
    kTrees
  };
  TTree* CreateTree(TreeIndex t);
  void EnableTree(TreeIndex t) { fTreeStatus[t] = kTRUE; };
  void DisableTree(TreeIndex t) { fTreeStatus[t] = kFALSE; };
  static const TString TreeName[kTrees];  //! Names of the TTree containers
  static const TString TreeTitle[kTrees]; //! Titles of the TTree containers

  AliAnalysisFilter fTrackFilter; // Standard track filter object
private:
  AliEventCuts fEventCuts;      //! Standard event cuts
  AliESDEvent *fESD = nullptr;  //! input event
  TList *fOutputList = nullptr; //! output list

  Bool_t fTreeStatus[kTrees];         // Status of the trees i.e. kTRUE (enabled) or kFALSE (disabled)
  TTree* fTree[kTrees] = { nullptr }; //! Array with all the output trees
  void FillTree(TreeIndex t);

  int fNumberOfEventsPerCluster = 1000;

  // fEventTree variables  

  // Event variables
  ULong64_t fEventId = 0u;      /// Event unique id
  Float_t fVtxX = -999.f;       /// Primary vertex x coordinate
  Float_t fVtxY = -999.f;       /// Primary vertex y coordinate
  Float_t fVtxZ = -999.f;       /// Primary vertex z coordinate
  Float_t fCentFwd = -1.f;      /// Centrality/Multiplicity percentile estimated with forward detectors
  Float_t fCentBarrel = -1.f;   /// Centrality/Multiplicity percentile estimated with barrel detectors

  Float_t fEventTime[10] = { -999.f };    /// Event time (t0) obtained with different methods (best, T0, T0-TOF, ...) for the whole event as a function of momentum
  Float_t fEventTimeRes[10] = { -999.f }; /// Resolution on the event time (t0) obtained with different methods (best, T0, T0-TOF, ...) for the whole event as a function of momentum
  UChar_t fEventTimeMask[10] = { 0u };    /// Mask with the method used to compute the event time (0x1=T0-TOF,0x2=T0A,0x3=TOC) for each momentum bins

  // fTrackTree variables

  // Identifier to associate tracks to collisions.
  // The inversed association (collisions to tracks) is generated on the fly

  // Coordinate system parameters
  Float_t fX = -999.f;     /// X coordinate for the point of parametrisation
  Float_t fAlpha = -999.f; /// Local <--> global coor.system rotation angle

  // Track parameters
  Float_t fY = -999.f;          /// fP[0] local Y-coordinate of a track (cm)
  Float_t fZ = -999.f;          /// fP[1] local Z-coordinate of a track (cm)
  Float_t fSnp = -999.f;        /// fP[2] local sine of the track momentum azimuthal angle
  Float_t fTgl = -999.f;        /// fP[3] tangent of the track momentum dip angle
  Float_t fSigned1Pt = -999.f;  /// fP[4] 1/pt (1/(GeV/c))

  // Covariance matrix
  Float_t fCYY = -999.f;       /// fC[0]
  Float_t fCZY = -999.f;       /// fC[1]
  Float_t fCZZ = -999.f;       /// fC[2]
  Float_t fCSnpY = -999.f;     /// fC[3]
  Float_t fCSnpZ = -999.f;     /// fC[4]
  Float_t fCSnpSnp = -999.f;   /// fC[5]
  Float_t fCTglY = -999.f;     /// fC[6]
  Float_t fCTglZ = -999.f;     /// fC[7]
  Float_t fCTglSnp = -999.f;   /// fC[8]
  Float_t fCTglTgl = -999.f;   /// fC[9]
  Float_t fC1PtY = -999.f;     /// fC[10]
  Float_t fC1PtZ = -999.f;     /// fC[11]
  Float_t fC1PtSnp = -999.f;   /// fC[12]
  Float_t fC1PtTgl = -999.f;   /// fC[13]
  Float_t fC1Pt21Pt2 = -999.f; /// fC[14]

  // Additional track parameters
  Float_t fTPCinnerP = -999.f; /// Full momentum at the inner wall of TPC for dE/dx PID

  // Track quality parameters
  ULong64_t fFlags = 0u;       /// Reconstruction status flags

  // Clusters
  UChar_t fITSClusterMap = 0u; /// ITS map of clusters, one bit per a layer
  UShort_t fTPCncls = 0u;      /// number of clusters assigned in the TPC
  UChar_t fTRDntracklets = 0u; /// number of TRD tracklets used for tracking/PID (TRD/TOF pattern)

  // Chi2
  Float_t fITSchi2Ncl = -999.f; /// chi2/Ncl ITS
  Float_t fTPCchi2Ncl = -999.f; /// chi2/Ncl TPC
  Float_t fTRDchi2 = -999.f;    /// chi2 TRD match (?)
  Float_t fTOFchi2 = -999.f;    /// chi2 TOF match (?)

  // PID
  Float_t fTPCsignal = -999.f; /// dE/dX TPC
  Float_t fTRDsignal = -999.f; /// dE/dX TRD
  Float_t fTOFsignal = -999.f; /// TOFsignal
  Float_t fLength = -999.f;    /// Int.Lenght @ TOF

  // TOF
  Int_t fTOFChannel = -1;    /// Index of the matched channel
  Short_t fTOFncls = -1;     /// Number of matchable clusters of the track
  Float_t fDx = -1;          /// Residual along x
  Float_t fDz = -1;          /// Residual along z
  Float_t fToT = -1;         /// ToT
  Float_t fLengthRatio = -1; /// Ratio of the integrated track length @ TOF to the cluster with respect to the matched cluster

  // fCaloTree variables
  Short_t fCellNumber = -1;     /// Cell absolute Id. number
  Float_t fAmplitude = -999.f;  /// Cell amplitude (= energy!)
  Float_t fTime = -999.f;       /// Cell time
  Char_t fType = -1;            /// Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)

  ClassDef(AliAnalysisTaskAO2Dconverter, 1);
};

#endif
