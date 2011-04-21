#ifndef ALIRSNVALUEPAIR_H
#define ALIRSNVALUEPAIR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Values which depend on 4-momentum of the pair.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnValue.h"

class AliRsnValuePair : public AliRsnValue {
public:

   enum EType {
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
      kTypes
   };

   AliRsnValuePair(const char *name = "valPair", EType type = kTypes);
   AliRsnValuePair(const AliRsnValuePair& copy);
   AliRsnValuePair& operator=(const AliRsnValuePair& copy);
   virtual ~AliRsnValuePair() { }

   void             SetType(EType type)  {fType = type;}
   EType            GetType()     const  {return fType;}
   const char*      GetTypeName() const;

   virtual Bool_t   Eval(TObject *object);

protected:

   EType           fType;                //  type from enumeration

   ClassDef(AliRsnValuePair, 1)  // AliRsnValuePair class
};

#endif
