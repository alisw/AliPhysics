#ifndef ALIRSNMINIVALUE_H
#define ALIRSNMINIVALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
//
//  Values which depend on 4-momentum of the pair.
//
////////////////////////////////////////////////////////////////////////////////

class AliRsnMiniPair;
class AliRsnMiniEvent;

class AliRsnMiniValue : public TNamed {
public:

   enum EType {
      kVz,            // event Z position of primary vertex
      kSpherocity,    // Spherocity
      kMult,          // event multiplicity or centrality (depends on task settings)
      kRefMult,       // event reference multiplicity (depends on task settings) - may differ from centrality estimator
      kTracklets,     // event tracklets
      kPlaneAngle,    // event reaction plane angle
      kLeadingPt,     // event leading particle momentum
      kEventCuts,     // -- limit of event cuts ----------------------------------------------------
      kPt,            // pair transverse momentum
      kPz,            // pair longitudinal momentum
      kInvMass,       // pair invariant mass (with reconstructed momenta)
      kInvMassMother, // pair invariant mass, always returns mass of mother
      kInvMassRes,    // pair invariant mass resolution
      kInvMassDiff,   // pair invariant mass difference (MC - reconstructed)
      kEta,           // pair pseudo-rapidity
      kMt,            // pair transverse mass (need a reference mass)
      kY,             // pair rapidity (need a reference mass)
      kPtRatio,       // ratio |pt1 - pt2|/(pt1 + pt2) of daughter transverse momenta
      kDipAngle,      // inverse cosine of the angle between daughter vector momenta
      kCosThetaStar,  // polarization angle
      kCosThetaStarAbs,  // polarization angle
      kCosThetaJackson,  // polarization angle in Jackson frame
      kCosThetaTransversity, // polarization angle in transversity frame
      kCosThetaToEventPlane, // polarization angle with respect to Event Plane
      kAngleLeading,  // angle to leading particle
      kFirstDaughterPt,  //pt of the first daughter of the pair
      kSecondDaughterPt, //pt of the second daughter of the pair
      kFirstDaughterP,   //p of the first daughter of the pair
      kSecondDaughterP,  //p of the second daughter of the pair
      kDCAproduct,    // product of the daughter's dca to PV (same in AliRsnValuePair)
      kFirstDaughterDCA,  //DCA to PV of the first daughter of the pair
      kSecondDaughterDCA, //DCA to PV of the second daughter of the pair
      kNSisters,    // number of daughters (only for MC)
      kPairPtRes,       // pair pT resolution
      kPairYRes,        // pair rapidity resolution
      kPhiV,   // PhiV calculation
      kAsym,   // pair asymmetry
      kTypes          // -- general limit ----------------------------------------------------------
   };

   AliRsnMiniValue(EType type = kTypes, Bool_t useMC = kFALSE);
   AliRsnMiniValue(const AliRsnMiniValue &copy);
   AliRsnMiniValue &operator=(const AliRsnMiniValue &copy);
   virtual ~AliRsnMiniValue() { }

   void               SetType(EType type)   {fType = type;}
   EType              GetType()      const  {return fType;}
   const char        *GetTypeName()  const  {return TypeName(fType);}
   Bool_t             IsEventValue() const  {return (fType < kEventCuts);}

   Float_t            Eval(AliRsnMiniPair *pair, AliRsnMiniEvent *event = 0x0);

   static const char *TypeName(EType type);
   static const char *ValueName(EType type, Bool_t useMC);

protected:

   EType            fType;            //  type from enumeration
   Bool_t           fUseMCInfo;       //  switch to use rec/sim momentum

   ClassDef(AliRsnMiniValue, 1)       //  AliRsnMiniValue class
};

inline const char *AliRsnMiniValue::ValueName(EType type, Bool_t useMC)
{
//
// Define a criterion to name these object.
// They are not managed by the user, since each object is a singleton
//

   if (useMC)
      return Form("MC_%s", TypeName(type));
   else
      return TypeName(type);
}

#endif
