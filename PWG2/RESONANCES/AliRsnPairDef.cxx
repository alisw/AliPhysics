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
#include "AliRsnDaughter.h"
#include "AliRsnPairDef.h"

ClassImp(AliRsnPairDef)

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef() : fMotherPDG(0)
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
  for (i = 0; i < 2; i++)
  {
    fCharge[i] = '0';
    fMass[i] = 0.0;
    fType[i] = AliRsnPID::kUnknown;
  }
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef
(Char_t sign1, AliRsnPID::EType type1, Char_t sign2, AliRsnPID::EType type2, Int_t motherPDG) :
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
(AliRsnPID::EType type1, Char_t sign1, AliRsnPID::EType type2, Char_t sign2, Int_t motherPDG) :
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

  fMotherPDG = copy.fMotherPDG;
  SetPair(copy.fCharge[0], copy.fType[0], copy.fCharge[1], copy.fType[1]);

  return (*this);
}

//_____________________________________________________________________________
Bool_t AliRsnPairDef::SetPairElement(Int_t i, Char_t charge, AliRsnPID::EType type)
{
//
// Set one element of the pair
// and returns warnings if the type is not valid.
//
  if (i < 0 || i > 1)
  {
    AliError("Index out of range");
    return kFALSE;
  }
  if (charge != '+' && charge != '-')
  {
    AliError(Form("Character '%c' not recognized as charge sign"));
    return kFALSE;
  }
  if (type < AliRsnPID::kElectron && type > AliRsnPID::kUnknown)
  {
    AliError("Type index out of enumeration range");
    return kFALSE;
  }
  fCharge[i] = charge;
  fType[i] = type;
  fMass[i] = AliRsnPID::ParticleMass(type);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnPairDef::SetPair
(Char_t charge1, AliRsnPID::EType type1, Char_t charge2, AliRsnPID::EType type2)
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
Double_t AliRsnPairDef::ComputeWeight(AliRsnDaughter *d0, AliRsnDaughter *d1)
{
//
// Compute a weight for filling the histograms:
// probability of first track to be identified as 'type[0]' times
// the probability of second track to be identified as 'type[1]',
// according to the order of appearance in argument list.
//

  Double_t prob0 = d0->PIDProb()[fType[0]];
  Double_t prob1 = d1->PIDProb()[fType[1]];

  return prob0*prob1;
}

//_____________________________________________________________________________
TString AliRsnPairDef::GetPairName()
{
//
// Returns a compact string with the name of the pair,
// to be used for naming objects related to it.
//

  TString sName;
  sName += AliRsnPID::ParticleName(fType[0]);
  sName += fCharge[0];
  sName += AliRsnPID::ParticleName(fType[1]);
  sName += fCharge[1];

  return sName;
}
