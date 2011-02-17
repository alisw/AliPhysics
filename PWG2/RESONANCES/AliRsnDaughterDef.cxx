//
// Class AliRsnDaughterDef
//
// Defines a decay channel for a resonance,
// resulting in a specified PDG code for the mother,
// and the particle type for the daughters, defined
// according to the internal PID format of the package
//
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"
#include "AliPID.h"
#include "AliRsnDaughterDef.h"

ClassImp(AliRsnDaughterDef)

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef() : 
   fMass(0.0),
   fCharge(0),
   fPID(AliPID::kUnknown),
   fDaughterType(AliRsnDaughter::kNoType)
{
//
// Default constructor.
// Initializes the data members to default (meaningless) values.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef
(AliPID::EParticleType type, Char_t sign) :
   fMass(0.0),
   fCharge(0),
   fPID(AliPID::kUnknown),
   fDaughterType(AliRsnDaughter::kNoType)
{
//
// Constructor with arguments.
// This constructor allows to define all the working parameters.
//

   SetDaughter(type, sign);
}


//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(const AliRsnDaughterDef &copy) :
   TObject(copy),
   fMass(copy.fMass),
   fCharge(copy.fCharge),
   fPID(copy.fPID),
   fDaughterType(copy.fDaughterType)
{
//
// Copy constructor with standard behavior
//
}

//_____________________________________________________________________________
const AliRsnDaughterDef& AliRsnDaughterDef::operator=(const AliRsnDaughterDef &copy)
{
//
// Assignment operator with standard behavior.
//

   fMass = copy.fMass;
   fCharge = copy.fCharge;
   fPID = copy.fPID;
   fDaughterType = copy.fDaughterType;

   return (*this);
}

//_____________________________________________________________________________
Bool_t AliRsnDaughterDef::SetDaughter(AliPID::EParticleType type, Char_t charge)
{
//
// Set one element of the pair
// and returns warnings if the type is not valid.
//

   AliPID pid;

   if (charge != '+' && charge != '-' && charge != '0') {
      AliError(Form("Character '%c' not recognized as charge sign", charge));
      return kFALSE;
   }
   if (type < 0 && type > (Int_t)AliPID::kSPECIESN) {
      AliError("Type index out of enumeration range");
      return kFALSE;
   }

   fCharge = charge;
   fPID    = type;
   fMass   = pid.ParticleMass(type);
   if ((Int_t)type < AliPID::kSPECIES) fDaughterType = AliRsnDaughter::kTrack;
   else if (type == AliPID::kKaon0) {
      fDaughterType = AliRsnDaughter::kV0;
      fCharge = '0';
   } else return kFALSE;

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnDaughterDef::MatchesDaughter(AliRsnDaughter *checked, Bool_t truePID)
{
//
// Checks if the argument matches the definitions here.
// If second argument is kTRUE and checked daughter has MC,
// check also that the particle species is correct.
//

   if (checked->RefType() != fDaughterType) return kFALSE;
   if (checked->ChargeChar() != fCharge) return kFALSE;
   if (truePID && checked->GetRef()) if (TMath::Abs(checked->GetPDG()) != AliPID::ParticleCode(fPID)) return kFALSE;

   return kTRUE;
}
