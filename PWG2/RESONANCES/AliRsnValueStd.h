#ifndef ALIRSNVALUESTD_H
#define ALIRSNVALUESTD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Collection of all values which can be computed within the package
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnValue.h"

class AliRsnValueStd : public AliRsnValue {
public:

   enum EValueType {
      kTrackP,               // single track total momentum
      kTrackPt,              // single track transverse momentum
      kTrackPtpc,            // single track total momentum in the TPC inner wall
      kTrackEta,             // single track pseudo-rapidity
      kTrackY,               // single track rapidity
      kTrackITSsignal,       // single track ITS signal
      kTrackTPCsignal,       // single track TPC signal
      kTrackTOFsignal,       // single track TOF signal
      kTrackTOFbeta,         // single track beta from TOF
      kTrackLength,          // single track integrated length
      kTrackValues,          // --- limit for track values -----------------------------------------
                             
      kPairP1,               // total momentum of 1st daughter of a pair
      kPairP2,               // total momentum of 2nd daughter of a pair
      kPairP1t,              // transverse momentum of 1st daughter of a pair
      kPairP2t,              // transverse momentum of 2nd daughter of a pair
      kPairP1z,              // longitudinal momentum of 1st daughter of a pair
      kPairP2z,              // longitudinal momentum of 2nd daughter of a pair
      kPairInvMass,          // pair invariant mass (with reconstructed momenta)
      kPairInvMassMC,        // pair invariant mass (with MC momenta)
      kPairInvMassRes,       // pair invariant mass resolution
      kPairPt,               // pair transverse momentum
      kPairPz,               // pair longitudinal momentum
      kPairEta,              // pair pseudo-rapidity
      kPairMt,               // pair transverse mass (need a reference mass)
      kPairY,                // pair rapidity (need a reference mass)
      kPairPhi,              // pair azimuthal angle (with reconstructed momenta)
      kPairPhiMC,            // pair azimuthal angle (with MC momenta)
      kPairPtRatio,          // ratio |pt1 - pt2|/(pt1 + pt2) of daughter transverse momenta
      kPairDipAngle,         // inverse cosine of the angle between daughter vector momenta
      kPairCosThetaStar,     // polarization angle
      kPairQInv,             // invariant relative momentum of the two daughters
      kPairAngleToLeading,   // angle between pair momentum and leading particle
      kPairValues,           // --- limit for pair values ------------------------------------------
                             
      kEventLeadingPt,       // transverse momentum of the event leading particle
      kEventMult,            // multiplicity computed as the number of tracks
      kEventMultMC,          // multiplicity from MC
      kEventMultESDCuts,     // multiplicity of good quality tracks
      kEventMultSPD,         // multiplicity from SPD
      kEventVz,              // Z position of event primary vertex
      kEventCentralityV0,    // event centrality (V0 method)
      kEventCentralityTrack, // event centrality (tracks method)
      kEventCentralityCL1,   // event centrality (CL1 method)
      kValueTypes            // --- limit for event values (and global) ----------------------------
   };

   AliRsnValueStd();
   AliRsnValueStd(const char *name, EValueType type, Int_t nbins = 0, Double_t min = 0.0, Double_t max = 0.0);
   AliRsnValueStd(const char *name, EValueType type, Double_t min, Double_t max, Double_t step);
   AliRsnValueStd(const char *name, EValueType type, Int_t nbins, Double_t *array);
   AliRsnValueStd(const AliRsnValueStd& copy);
   AliRsnValueStd& operator=(const AliRsnValueStd& copy);
   virtual ~AliRsnValueStd() { /*does nothing, since pointers are not owned by this object*/ }

   EValueType        GetValueType() const           {return fValueType;}
   const char*       GetValueTypeName() const;
   TObject*          GetSupportObject()             {return fSupportObject;}
   void              SetSupportObject(TObject *obj) {fSupportObject = obj;}
   void              SetValueType(EValueType type)  {fValueType = type; fTargetType = TargetType(type);}

   virtual Bool_t    Eval(TObject *object, Bool_t useMC = kFALSE);
   virtual void      Print(Option_t *option = "") const;
   static  RSNTARGET TargetType(EValueType type);

protected:

   EValueType   fValueType;      // value type
   TObject     *fSupportObject;  // support object needed for computing some of the values

   ClassDef(AliRsnValueStd, 1)   // AliRsnValueStd class
};

#endif
