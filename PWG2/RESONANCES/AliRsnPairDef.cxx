//
// Class AliRsnPairDef
//
// Defines a decay channel for a resonance,
// resulting in a specified PDG code for the mother,
// and the particle type for the daughters, defined
// according to the internal PID format of the package
//
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"

ClassImp(AliRsnPairDef)

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef() : fMotherMass(0.0), fMotherPDG(0), fDef1(), fDef2()
{
//
// Empty constructor.
// Initializes the data members to default values:
//  - no definition of particles in the pair;
//  - histogram binning undefined.
// When using this constructor, all analysis elements (particles, histogram)
// must be defined before starting event processing.
//
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef
(AliPID::EParticleType type1, Char_t sign1, AliPID::EParticleType type2, Char_t sign2, Int_t motherPDG, Double_t motherMass) :
   fMotherMass(motherMass),
   fMotherPDG(motherPDG),
   fDef1(type1, sign1),
   fDef2(type2, sign2)
{
//
// Constructor with arguments.
// This constructor allows to define all the working parameters.
//
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef(const AliRsnPairDef &copy) :
   TObject(copy),
   fMotherMass(copy.fMotherMass),
   fMotherPDG(copy.fMotherPDG),
   fDef1(copy.fDef1),
   fDef2(copy.fDef2)
{
//
// Copy constructor with standard behavior
//
}

//_____________________________________________________________________________
const AliRsnPairDef& AliRsnPairDef::operator=(const AliRsnPairDef &copy)
{
//
// Assignment operator with standard behavior.
//

   fMotherMass = copy.fMotherMass;
   fMotherPDG = copy.fMotherPDG;
   fDef1 = copy.fDef1;
   fDef2 = copy.fDef2;

   return (*this);
}

//_____________________________________________________________________________
const char* AliRsnPairDef::GetPairName() const
{
//
// Returns a compact string with the name of the pair,
// to be used for naming objects related to it.
//

   return Form("%s%s", fDef1.GetName(), fDef2.GetName());
}
