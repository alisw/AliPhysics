#ifndef ALIRSNVALUEDAUGHTER_H
#define ALIRSNVALUEDAUGHTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Values which depend on 4-momentum of the daughters.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnValue.h"

class AliRsnValueDaughter : public AliRsnValue {
public:

   enum EType {
      kP,          // total momentum
      kPt,         // transverse momentum
      kPtpc,       // total momentum in the TPC inner wall
      kEta,        // pseudo-rapidity
      kITSsignal,  // ITS signal
      kTPCsignal,  // TPC signal
      kTOFsignal,  // TOF signal
      kTPCnsigmaPi,// TPC number of sigmas pion
      kTPCnsigmaK, // TPC number of sigmas kaon
      kTPCnsigmaP, // TPC number of sigmas proton
      kTOFnsigmaPi,// TOF number of sigmas pion
      kTOFnsigmaK, // TOF number of sigmas kaon
      kTOFnsigmaP, // TOF number of sigmas proton
      kTypes
   };

   AliRsnValueDaughter(const char *name = "valDaughter", EType type = kTypes);
   AliRsnValueDaughter(const AliRsnValueDaughter& copy);
   AliRsnValueDaughter& operator=(const AliRsnValueDaughter& copy);
   virtual ~AliRsnValueDaughter() { }

   void             SetType(EType type)  {fType = type;}
   EType            GetType()     const  {return fType;}
   const char*      GetTypeName() const;

   virtual Bool_t   Eval(TObject *object);

protected:

   EType           fType;                //  type from enumeration

   ClassDef(AliRsnValueDaughter, 1)  // AliRsnValueDaughter class
};

#endif
