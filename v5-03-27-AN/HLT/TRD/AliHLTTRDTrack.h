//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDTRACK_H
#define ALIHLTTRDTRACK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliTRDtrackV1.h"
#include "AliHLTLogging.h"

class AliHLTTRDTracklet;

class AliHLTTRDTrack
{
 public:
  AliHLTTRDTrack();
  AliHLTTRDTrack(const AliTRDtrackV1* const inTrack);
  ~AliHLTTRDTrack();

  void CopyDataMembers(const AliTRDtrackV1* const inTrack);
  void ExportTRDTrack(AliTRDtrackV1* const outTrack) const;
  AliHLTUInt8_t *GetEndPointer() const // Returns pointer to the end of the track
    { return ((AliHLTUInt8_t *) this + fSize); };
  AliHLTUInt32_t GetSize() const {return fSize;};
  void Print(Bool_t printTracklets = kTRUE) const;
  static AliHLTUInt32_t SaveAt(AliHLTUInt8_t *const block, const AliTRDtrackV1* const inTrack);
  static AliHLTUInt32_t LoadFrom(AliTRDtrackV1 *const outTrack, const AliHLTUInt8_t *const block);

 private:
  AliHLTTRDTrack(const AliHLTTRDTrack& inTrack);
  AliHLTTRDTrack& operator=(const AliHLTTRDTrack& inTrack);
  void InitArrays();

  /* Probably need */
  Float_t      fPID[AliPID::kSPECIES];//  PID probabilities
  Float_t      fBudget[3];            //  Integrated material budget
  Float_t      fDE;                   //  Integrated delta energy

  /* ======== From AliKalmanTrack ======== */
  
  /* Defenitely need */
  Float_t      fFakeRatio;            // fake ratio
  Float_t      fChi2;                 // total chi2 value for this track
  // Float_t      fMass;                 // mass hypothesis
  // Int_t        fLab;                  // track label

  /* Probably need */
  Int_t        fN;                    // number of associated clusters
  Float_t      fIntegratedLength;     // integrated length  // variables for time integration (S.Radomski@gsi.de)

  /* ======= From AliExternalTrackParam ======== */

  /* Defenitely need */
  Float_t      fX;                    // X coordinate for the point of parametrisation
  Float_t      fAlpha;                // Local <-->global coor.system rotation angle
  Float_t      fP[5];                 // The track parameters
  Float_t      fC[15];                // The track parameter covariance matrix

  /* Not need */
  //  static Float_t    fgMostProbablePt; // "Most probable" pt (to be used if Bz=0)

  AliHLTUInt32_t fSize;               // Size of the track with tracklets and clusters in the memory
  UChar_t      fBits;
  Bool_t       fTrackletAtPlane[AliTRDtrackV1::kNplane];   // Used positions in the original array of tracklets


};

#endif
