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
#include "AliPID.h"
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
    fPID[i] = AliPID::kUnknown;
    fDaughterType[i] = AliRsnDaughter::kTrack;
  }
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

  SetDaughters(type1, sign1, type2, sign2);
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

  SetDaughters(copy.fPID[0], copy.fCharge[0], copy.fPID[1], copy.fCharge[1]);
}

//_____________________________________________________________________________
const AliRsnPairDef& AliRsnPairDef::operator=(const AliRsnPairDef &copy)
{
//
// Assignment operator with standard behavior.
//

  fMotherMass = copy.fMotherMass;
  fMotherPDG = copy.fMotherPDG;
  SetDaughters(copy.fPID[0], copy.fCharge[0], copy.fPID[1], copy.fCharge[1]);

  return (*this);
}

//_____________________________________________________________________________
Bool_t AliRsnPairDef::SetDaughter(Int_t i, AliPID::EParticleType type, Char_t charge)
{
//
// Set one element of the pair
// and returns warnings if the type is not valid.
//

  AliPID pid;

  if (i < 0 || i > 1) 
  {
    AliError("Index out of range");
    return kFALSE;
  }
  if (charge != '+' && charge != '-' && charge != '0')
  {
    AliError(Form("Character '%c' not recognized as charge sign", charge));
    return kFALSE;
  }
  if (type < 0 && type > (Int_t)AliPID::kSPECIESN) 
  {
    AliError("Type index out of enumeration range");
    return kFALSE;
  }
  
  fCharge[i] = charge;
  fPID[i] = type;
  fMass[i] = pid.ParticleMass(type);
  if ((Int_t)type < AliPID::kSPECIES) fDaughterType[i] = AliRsnDaughter::kTrack;
  else if (type == AliPID::kKaon0) 
  {
    fDaughterType[i] = AliRsnDaughter::kV0;
    fCharge[i] = '0';
  }
  else return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnPairDef::SetDaughters
(AliPID::EParticleType type1, Char_t charge1, AliPID::EParticleType type2, Char_t charge2)
{
//
// Set both elements of the pair,
// returning logical AND of check for each one.
//
  Bool_t part1 = SetDaughter(0, type1, charge1);
  Bool_t part2 = SetDaughter(1, type2, charge2);

  return (part1 && part2);
}

//_____________________________________________________________________________
const char* AliRsnPairDef::GetPairName() const
{
//
// Returns a compact string with the name of the pair,
// to be used for naming objects related to it.
//

  return Form("%s%c%s%c", AliPID::ParticleShortName(fPID[0]), fCharge[0], AliPID::ParticleShortName(fPID[1]), fCharge[1]);
}
