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
#include "AliRsnPairDef.h"

ClassImp(AliRsnPairDef)

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef() : fMotherMass(0.0), fMotherPDG(0)
{
//
// Empty constructor.
// Initializes the data members to default values:
//  - no definition of particles in the pair;
//  - histogram binning undefined.
// When using this constructor, all analysis elements (particles, histogram)
// must be defined before starting event processing.
//

  Int_t i;
  for (i = 0; i < 2; i++) {
    fCharge[i] = '0';
    fMass[i] = 0.0;
    fType[i] = AliPID::kUnknown;
  }
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef
(Char_t sign1, AliPID::EParticleType type1, Char_t sign2, AliPID::EParticleType type2, Int_t motherPDG, Double_t motherMass) :
  fMotherMass(motherMass),
  fMotherPDG(motherPDG)
{
//
// Constructor with arguments.
// This constructor allows to define all the working parameters.
//

  SetPair(sign1, type1, sign2, type2);
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef
(AliPID::EParticleType type1, Char_t sign1, AliPID::EParticleType type2, Char_t sign2, Int_t motherPDG, Double_t motherMass) :
  fMotherMass(motherMass),
  fMotherPDG(motherPDG)
{
//
// Constructor with arguments.
// This constructor allows to define all the working parameters.
//

  SetPair(sign1, type1, sign2, type2);
}


//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef(const AliRsnPairDef &copy) :
    TObject(copy),
    fMotherMass(copy.fMotherMass),
    fMotherPDG(copy.fMotherPDG)
{
//
// Copy constructor with standard behavior
//

  SetPair(copy.fCharge[0], copy.fType[0], copy.fCharge[1], copy.fType[1]);
}

//_____________________________________________________________________________
const AliRsnPairDef& AliRsnPairDef::operator=(const AliRsnPairDef &copy)
{
//
// Assignment operator with standard behavior.
//

  fMotherMass = copy.fMotherMass;
  fMotherPDG = copy.fMotherPDG;
  SetPair(copy.fCharge[0], copy.fType[0], copy.fCharge[1], copy.fType[1]);

  return (*this);
}

//_____________________________________________________________________________
Bool_t AliRsnPairDef::SetPairElement(Int_t i, Char_t charge, AliPID::EParticleType type)
{
//
// Set one element of the pair
// and returns warnings if the type is not valid.
//

  AliPID pid;

  if (i < 0 || i > 1) {
    AliError("Index out of range");
    return kFALSE;
  }
  if (charge != '+' && charge != '-') {
    AliError(Form("Character '%c' not recognized as charge sign", charge));
    return kFALSE;
  }
  if (type < 0 && type > (Int_t)AliPID::kSPECIES) {
    AliError("Type index out of enumeration range");
    return kFALSE;
  }
  fCharge[i] = charge;
  fType[i] = type;
  fMass[i] = pid.ParticleMass(type);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnPairDef::SetPair
(Char_t charge1, AliPID::EParticleType type1, Char_t charge2, AliPID::EParticleType type2)
{
//
// Set both elements of the pair,
// returning logical AND of check for each one.
//
  Bool_t part1 = SetPairElement(0, charge1, type1);
  Bool_t part2 = SetPairElement(1, charge2, type2);

  return (part1 && part2);
}

//_____________________________________________________________________________
TString AliRsnPairDef::GetPairName() const
{
//
// Returns a compact string with the name of the pair,
// to be used for naming objects related to it.
//

  TString sName;
  sName += AliPID::ParticleShortName(fType[0]);
  sName += fCharge[0];
  sName += AliPID::ParticleShortName(fType[1]);
  sName += fCharge[1];

  return sName;
}
