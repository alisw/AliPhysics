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
      kY,          // rapidity
      kMass,       // mass
      kITSsignal,  // ITS signal
      kTPCsignal,  // TPC signal
      kTOFsignal,  // TOF signal
      kTPCnsigmaPi,// TPC number of sigmas pion
      kTPCnsigmaK, // TPC number of sigmas kaon
      kTPCnsigmaP, // TPC number of sigmas proton
      kTOFnsigmaPi,// TOF number of sigmas pion
      kTOFnsigmaK, // TOF number of sigmas kaon
      kTOFnsigmaP, // TOF number of sigmas proton
      kTOFdeltaPi, // TOF t-texp for pion hypothesis
      kTOFdeltaK, // TOF t-texp for pion hypothesis
      kTOFdeltaP, // TOF t-texp for pion hypothesis
      kNITSclusters,  // n ITS clusters
      kNTPCclusters,  // n TPC clusters
      kNTPCcrossedRows,  // n TPC crossed rows
      kNTPCcrossedRowsFclusters,  // n TPC crossed rows over findable clusters
      kITSchi2,     // ITS chi^2
      kTPCchi2,     // TPC chi^2
      kDCAXY,       // DCA xy
      kDCAZ,        // DCA z
      kCharge,     // charge
      kPhi,        // azimuthal angle at vertex
      kPhiOuterTPC,// azimuthal angle at TPC outer radius
      kV0DCA,       // V0 DCA
      kV0Radius,       // V0 radius
      kV0Mass,       // V0 mass
      kV0P,       // V0 momentum
      kV0Pt,       // V0 transverse momentum
      kV0NPt,       // transverse momentum of negative V0 daughter
      kV0PPt,       // transverse momentum of positive V0 daughter
      kV0DCAXY,     // DCA for Secondary Tracks to Primary Vertex
      kV0Lifetime,   //Lifetime for V0 particles
      kDaughterDCA, // DCA of V0 Daughters
      kCosPointAng, // V0 Cosine of Pointing Angle
      kLambdaProtonPIDCut,         //V0 - Lambda number of sigmas proton
      kAntiLambdaAntiProtonPIDCut, //V0 - AntiLambda number of sigmas antiproton
      kLambdaPionPIDCut,	          //V0 - Lambda number of sigmas pion
      kAntiLambdaAntiPionPIDCut,   //V0 - AntiLambda number of sigmas pion
      kK0SMass,               //V0 - mass for K0S hypothesis
      kLambdaMass,            //V0 - mass for Lambda hypothesis
      kAntiLambdaMass,        //V0 - mass for anti-Lambda hypothesis
      kXiMass,                //Cascade - mass for Xi hypothesis
      kOmegaMass,             //Cascade - mass for Omega hypothesis
      kCascadeP,              //Cascade - momentum
      kCascadePt,             //Cascade - pT
      kCascadeDCA,            //Cascade - DCA to primary vertex
      kCascadeRadius,         //Cascade - radius (cylindrical coordinates)
      kCascadeDaughterDCA,    //Cascade - DCA of daughters
      kCascadeCosPointAng,    //Cascade - Cosine of Pointing Angle of Cascade to pimary vertex
      kCascadeV0CosPointAng,  //Cascade - Cosine of Pointing Angle of V0 to Cascade vertex
      kCascadeV0Lifetime,     //Cascade - Lifetime of the V0 from the Cascade
      kCascadeV0Pt,           //Cascade - pT of V0
      kBachelorPt,            //Cascade - pT of bachelor track
      kBachelorPionTPCnsigma, //Cascade - TPC nsigma of bachelor track for pion hypothesis
      kBachelorKaonTPCnsigma, //Cascade - TPC nsigma of bachelor track for kaon hypothesis
      kTypes
   };

   AliRsnValueDaughter(const char *name = "valDaughter", EType type = kTypes);
   AliRsnValueDaughter(const AliRsnValueDaughter &copy);
   AliRsnValueDaughter &operator=(const AliRsnValueDaughter &copy);
   virtual ~AliRsnValueDaughter() { }

   void             SetType(EType type)  {fType = type;}
   EType            GetType()     const  {return fType;}
   const char      *GetTypeName() const;

   virtual Bool_t   Eval(TObject *object);

protected:

   EType           fType;                //  type from enumeration

   ClassDef(AliRsnValueDaughter, 1)  // AliRsnValueDaughter class
};

#endif
