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
// Constructor.
// This version of constructor leaves all undefined,
// and all daughters will be accepted.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(AliPID::EParticleType type, Char_t sign) :
   fMass(0.0),
   fCharge(sign),
   fPID(type),
   fDaughterType(AliRsnDaughter::kNoType)
{
//
// Constructor.
// This version of constructor initializes the PID type
// and the charge (optional, leave its default to include both),
// and calls 'SetDaughter()' to assign the object type accordingly.
//

   SetDaughter(type, sign);
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(AliRsnDaughter::ERefType refType, Char_t sign) :
   fMass(0.0),
   fCharge(sign),
   fPID(AliPID::kUnknown),
   fDaughterType(refType)
{
//
// Constructor.
// This version of constructor initialized the object type
// and the charge (optiona, leave its defaul to include both),
// and sets the PID type undefined.
//
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

   // charge and type come from arguments
   fCharge = charge;
   fPID    = type;
   
   // mass is assigned by AliPID, if possible
   fMass = 0.0;
   if ((Int_t)type >= 0 && (Int_t)type < AliPID::kSPECIESN)
      fMass = pid.ParticleMass(type);
   
   // object type is determined by type itself
   if ((Int_t)type >= 0 && (Int_t)type < AliPID::kSPECIES) 
      fDaughterType = AliRsnDaughter::kTrack;
   else if (type == AliPID::kKaon0) {
      fDaughterType = AliRsnDaughter::kV0;
      fCharge = '0';
   }
   else
      fDaughterType = AliRsnDaughter::kNoType;

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

   // charge matching:
   // if charge was set to '+', '-' or '0'
   // the checked daughter charge must match that defined in the settings
   // otherwise no specific charge requirement was done
   Bool_t chargeMatch = kTRUE;
   switch (fCharge) {
      case '+': chargeMatch = checked->IsPos();
      case '-': chargeMatch = checked->IsNeg();
      case '0': chargeMatch = checked->IsNeutral();
      default : chargeMatch = kTRUE;
   }
   
   // object type matching
   Bool_t objMatch = kTRUE;
   if (fDaughterType != AliRsnDaughter::kNoType)
      objMatch = (checked->RefType() == fDaughterType);
      
   // particle type matching (only if MC is available and second arg is true)
   Bool_t pidMatch = kTRUE;
   if (truePID && fPID != AliPID::kUnknown && checked->GetRefMC())
      pidMatch = (AliPID::ParticleCode(fPID) == checked->GetPDG());
      
   // return the AND of all
   return (chargeMatch && objMatch && pidMatch);
}
