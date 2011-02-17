//
// Class AliRsnDaughterDef
//
// Defines some requirements to select a good candidate daughter for a resonance
// among all available objects presend in each event (tracks, V0s, etc. ...).
// This includes:
// -- track type (charged track, V0, cascade)
// -- track charge
// -- track PID type and related mass (only for charged)
//
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNDAUGHTERDEF_H
#define ALIRSNDAUGHTERDEF_H

#include "AliPID.h"
#include "AliRsnDaughter.h"

class AliRsnDaughterDef : public TObject {
public:

   AliRsnDaughterDef();
   AliRsnDaughterDef(AliPID::EParticleType type, Char_t charge = 0);
   AliRsnDaughterDef(AliRsnDaughter::ERefType refType, Char_t charge = 0);
   AliRsnDaughterDef(const AliRsnDaughterDef &copy);
   const AliRsnDaughterDef& operator= (const AliRsnDaughterDef &copy);
   virtual ~AliRsnDaughterDef() { }

   // getters
   Double_t                 GetMass()         const {return fMass;}
   Char_t                   GetCharge()       const {return fCharge;}
   Short_t                  GetChargeShort()  const {if (fCharge == '+') return 1; else if (fCharge == '-') return -1; else return 0;}
   AliPID::EParticleType    GetPID()          const {return fPID;}
   AliRsnDaughter::ERefType GetDaughterType() const {return fDaughterType;}
   const char*              GetName()         const {return Form("%s%c", AliPID::ParticleShortName(fPID), fCharge);}

   // setters
   Bool_t SetDaughter(AliPID::EParticleType pid, Char_t charge);
   void   SetDaughterType(AliRsnDaughter::ERefType type) {fDaughterType = type;}
   
   // checker
   Bool_t MatchesDaughter(AliRsnDaughter *daughter, Bool_t truePID = kFALSE);

private:

   Double_t                  fMass;         // mass of particles
   Char_t                    fCharge;       // charge of particles
   AliPID::EParticleType     fPID;          // PID of particles
   AliRsnDaughter::ERefType  fDaughterType; // object type (track/V0)

   // ROOT dictionary
   ClassDef(AliRsnDaughterDef, 1)
};

#endif
