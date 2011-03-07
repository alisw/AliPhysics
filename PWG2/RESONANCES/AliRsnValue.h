#ifndef ALIRSNVALUE_H
#define ALIRSNVALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Collection of all values which can be computed within the package
//
////////////////////////////////////////////////////////////////////////////////

#include "TArrayD.h"
#include "AliRsnTarget.h"

class AliRsnValue : public AliRsnTarget {
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
      kTrackValues,          // --- limita for track values ----------------------------------------
                             
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

   AliRsnValue();
   AliRsnValue(const char *name, EValueType type, Int_t nbins = 0, Double_t min = 0.0, Double_t max = 0.0);
   AliRsnValue(const char *name, EValueType type, Double_t min, Double_t max, Double_t step);
   AliRsnValue(const char *name, EValueType type, Int_t nbins, Double_t *array);
   AliRsnValue(const AliRsnValue& copy);
   AliRsnValue& operator=(const AliRsnValue& copy);
   virtual ~AliRsnValue() { /*does nothing, since pointers are not owned by this object*/ }

   TArrayD     GetArray() const               {return fBinArray;}
   Double_t    GetComputedValue() const       {return fComputedValue;}
   EValueType  GetValueType() const           {return fValueType;}
   const char* GetValueTypeName() const;
   TObject*    GetSupportObject()             {return fSupportObject;}
   void        SetSupportObject(TObject *obj) {fSupportObject = obj;}
   void        SetValueType(EValueType type)  {fValueType = type; fTargetType = TargetType(type);}

   void        SetBins(Int_t n, Double_t min, Double_t max);
   void        SetBins(Int_t n, Double_t *array);
   void        SetBins(Double_t min, Double_t max, Double_t step);

   void        Set(EValueType type, Int_t n, Double_t min, Double_t max)       {SetValueType(type); SetBins(n, min, max);}
   void        Set(EValueType type, Int_t n, Double_t *array)                  {SetValueType(type); SetBins(n, array);}
   void        Set(EValueType type, Double_t min, Double_t max, Double_t step) {SetValueType(type); SetBins(min, max, step);}

   virtual Bool_t    Eval(TObject *object, Bool_t useMC = kFALSE);
   virtual void      Print(Option_t *option = "") const;
   static  RSNTARGET TargetType(EValueType type);

protected:

   Double_t     fComputedValue;  // computed value
   EValueType   fValueType;      // value type
   TArrayD      fBinArray;       // array of bins (when used for a histogram axis)
   TObject     *fSupportObject;  // support object needed for computing some of the values

   // ROOT dictionary
   ClassDef(AliRsnValue, 2)
};

#endif
