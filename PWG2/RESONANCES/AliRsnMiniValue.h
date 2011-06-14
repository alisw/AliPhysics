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
      kMult,          // event multiplicity or centrality (depends on task settings)
      kPlaneAngle,    // event reaction plane angle
      kLeadingPt,     // event leading particle momentum
      kEventCuts,     // -- limit of event cuts ----------------------------------------------------
      kPt,            // pair transverse momentum
      kPz,            // pair longitudinal momentum
      kInvMass,       // pair invariant mass (with reconstructed momenta)
      kInvMassRes,    // pair invariant mass resolution
      kEta,           // pair pseudo-rapidity
      kMt,            // pair transverse mass (need a reference mass)
      kY,             // pair rapidity (need a reference mass)
      kPtRatio,       // ratio |pt1 - pt2|/(pt1 + pt2) of daughter transverse momenta
      kDipAngle,      // inverse cosine of the angle between daughter vector momenta
      kCosThetaStar,  // polarization angle
      kAngleLeading,  // angle to leading particle
      kTypes          // -- general limit ----------------------------------------------------------
   };

   AliRsnMiniValue(EType type = kTypes, Bool_t useMC = kFALSE);
   AliRsnMiniValue(const AliRsnMiniValue& copy);
   AliRsnMiniValue& operator=(const AliRsnMiniValue& copy);
   virtual ~AliRsnMiniValue() { }

   void               SetType(EType type)   {fType = type;}
   EType              GetType()      const  {return fType;}
   const char*        GetTypeName()  const  {return TypeName(fType);}
   Bool_t             IsEventValue() const  {return (fType < kEventCuts);}
   
   Float_t            Eval(AliRsnMiniPair *pair, AliRsnMiniEvent *event = 0x0);
   
   static const char* TypeName(EType type);
   static const char* ValueName(EType type, Bool_t useMC);

protected:

   EType            fType;            //  type from enumeration
   Bool_t           fUseMCInfo;       //  switch to use rec/sim momentum
                                       
   ClassDef(AliRsnMiniValue, 1)       //  AliRsnMiniValue class
};

inline const char* AliRsnMiniValue::ValueName(EType type, Bool_t useMC)
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
