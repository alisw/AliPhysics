#include "AliLog.h"

#include "AliRsnMCInfo.h"

ClassImp(AliRsnMCInfo)

//_____________________________________________________________________________
AliRsnMCInfo::AliRsnMCInfo() : TObject(),
    fEnergy(0),
    fPDG(0),
    fMother(-1),
    fMotherPDG(0)
{
//
// Default constructor.
// Initializes all data-members with meaningless values.
//
  for (Int_t i = 0; i < 3; i++) fP[i] = 0.0;
}

//_____________________________________________________________________________
AliRsnMCInfo::AliRsnMCInfo(const AliRsnMCInfo & copy) : TObject(copy),
    fEnergy(copy.fEnergy),
    fPDG(copy.fPDG),
    fMother(copy.fMother),
    fMotherPDG(copy.fMotherPDG)
{
//
// Copy constructor.
// Initializes all data-members with meaningless values.
//

  for (Int_t i = 0; i < 3; i++) fP[i] = copy.fP[i];

}

//_____________________________________________________________________________
AliRsnMCInfo::~AliRsnMCInfo()
{
}

//_____________________________________________________________________________
void AliRsnMCInfo::Adopt(TParticle * particle)
{
//
// Copies data from a TParticle into "this":
//  - PDG code
//  - GEANT label of mother (if any, otherwise -1)
// If the argument is NULL, nothing is done, and an alert
// is given by the method.
//

  if (!particle)
  {
    AliError("NULL argument passed. Nothing done.");
    return;
  }

  fP[0] = particle->Px();
  fP[1] = particle->Py();
  fP[2] = particle->Pz();

  fEnergy = particle->Energy();

  fPDG    = particle->GetPdgCode();
  fMother = particle->GetFirstMother();
  fMotherPDG = (Short_t) 0;
}
