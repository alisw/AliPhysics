//
// Class AliRsnMCInfo
//
// Contains informations from the MonteCarlo particle is associated to a track.
// It is used when looking at "perfect" PID and at "true" pairs, but the user
// does not need to access its methods.
//

#include "AliLog.h"

#include "AliRsnMCInfo.h"

ClassImp(AliRsnMCInfo)

//_____________________________________________________________________________
AliRsnMCInfo::AliRsnMCInfo() :
  TObject(), fEnergy(0), fPDG(0), fMother(-1), fMotherPDG(0)
{
//
// Default constructor.
// Initializes all data-members to meaningless values.
//

  Int_t i;
  for (i = 0; i < 3; i++) fP[i] = fV[i] = 0.0;
}

//_____________________________________________________________________________
AliRsnMCInfo::AliRsnMCInfo(const AliRsnMCInfo & copy) :
  TObject(copy),
  fEnergy(copy.fEnergy),
  fPDG(copy.fPDG),
  fMother(copy.fMother),
  fMotherPDG(copy.fMotherPDG)
{
//
// Copy constructor.
//

  Int_t i;
  for (i = 0; i < 3; i++) {
    fP[i] = copy.fP[i];
    fV[i] = copy.fV[i];
  }
}

//_____________________________________________________________________________
AliRsnMCInfo::~AliRsnMCInfo()
{
//
// Destructor.
// Does nothing because there are no pointers to clean.
//
}

//_____________________________________________________________________________
void AliRsnMCInfo::Adopt(TParticle * particle)
{
//
// Uses a TParticle to initialize its data members.
// If the argument is NULL, nothing is done and an error message is returned.
//

  if (!particle) return;

  fP[0] = particle->Px();
  fP[1] = particle->Py();
  fP[2] = particle->Pz();

  fV[0] = particle->Vx();
  fV[1] = particle->Vy();
  fV[2] = particle->Vz();

  fEnergy = particle->Energy();

  fPDG    = particle->GetPdgCode();
  fMother = particle->GetFirstMother();
  fMotherPDG = (Short_t) 0;
}
