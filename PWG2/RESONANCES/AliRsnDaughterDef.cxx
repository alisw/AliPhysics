//
//  Definition for a single candidate daughter.
//  They can be chosen among the following possibilities:
//
//  -- particle species, chosen in the enum AliRsnDaughter::ESpecies;
//     this option does two things:
//     -- when possible (needs MC) and required, allow to check is a daughter
//        is of the defined species (e.g.: for selecting true daughters of a resonance)
//     -- defines the mass to be assigned to this object, for any purpose
//        (e.g.: to compute its 4-momentum to be used for a resonance mass)
//
//  -- track charge, which can be '+', '-', '0', 
//     -- any other char leaves this undefined (use only in daughter monitor loops)
//     -- when doing resonance analysis, or when RSN Input handler needs to be used,
//        this must always be defined
//     
//  -- object type (track/V0/cascade)
//     -- could be needed to select tracks when particle species is not specified
//     -- works only in single daughter loops
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//

#include "AliLog.h"
#include "AliRsnDaughterDef.h"

ClassImp(AliRsnDaughterDef)


//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef() :
   fPID(AliRsnDaughter::kUnknown),
   fMass(0.0),
   fCharge(0),
   fRefType(AliRsnDaughter::kNoType)
{
//
// This version of constructor leaves everything undefined:
// this will cause all daughters to be accepted.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(AliRsnDaughter::ESpecies type, Char_t sign) :
   fPID(type),
   fMass(AliRsnDaughter::SpeciesMass(type)),
   fCharge(sign),
   fRefType(AliRsnDaughter::RefType(type))
{
//
// This version of constructor initializes the PID type (and then the mass)
// and the charge (optional, leave 2nd argument to default to include both).
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(EPARTYPE type, Char_t sign) :
   fPID(AliRsnDaughter::FromAliPID(type)),
   fMass(AliRsnDaughter::SpeciesMass(AliRsnDaughter::FromAliPID(type))),
   fCharge(sign),
   fRefType(AliRsnDaughter::RefType(AliRsnDaughter::FromAliPID(type)))
{
//
// This version of constructor initializes the PID type (and then the mass)
// and the charge (optional, leave 2nd argument to default to include both).
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(AliRsnDaughter::ERefType refType, Char_t sign) :
   fPID(AliRsnDaughter::kUnknown),
   fMass(0.0),
   fCharge(sign),
   fRefType(refType)
{
//
// This version of constructor initializes the object type
// and the charge (optional, leave 2nd argument to default to include both),
// and leaves the PID type undefined.
// This is useful when one is interested in all tracks/V0s/cascades without
// requiring them to be identified as a certain species, but if one then requires
// an object linked to this definition to compute a rapidity or a transverse mass,
// this will not work.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(const AliRsnDaughterDef &copy) :
   TObject(copy),
   fPID(copy.fPID),
   fMass(copy.fMass),
   fCharge(copy.fCharge),
   fRefType(copy.fRefType)
{
//
// Copy constructor has standard behavior.
//
}

//_____________________________________________________________________________
const AliRsnDaughterDef& AliRsnDaughterDef::operator=(const AliRsnDaughterDef &copy)
{
//
// Assignment operator has standard behavior.
//

   fMass = copy.fMass;
   fCharge = copy.fCharge;
   fPID = copy.fPID;
   fRefType = copy.fRefType;

   return (*this);
}
