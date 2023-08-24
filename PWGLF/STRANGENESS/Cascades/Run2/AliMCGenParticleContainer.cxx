#include "AliMCGenParticleContainer.h"

class AliMCGenParticleContainer;    // your analysis class

ClassImp(AliMCGenParticleContainer) // classimp: necessary for root

//_____________________________________________________________________________
AliMCGenParticleContainer::AliMCGenParticleContainer() :
TObject                     (), 

fTreeMCVarIsPhysicalPrimary (kFALSE),

fTreeMCVarPdgCode           (-9999),
fTreeMCVarMotherPdgCode     (-9999),
fTreeMCVarLabel             (-999),

fTreeMCVarCharge            (-999),
fTreeMCVarRap               (-999),
fTreeMCVarEta               (-999),
fTreeMCVarPtot              (-100),
fTreeMCVarPt                (-100),
fTreeMCVarPx                (-100),
fTreeMCVarPy                (-100),
fTreeMCVarPz                (-100),
fTreeMCVarPhi               (-1000)

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliMCGenParticleContainer::~AliMCGenParticleContainer() // destructor
{
    // destructor
}
//_____________________________________________________________________________
void AliMCGenParticleContainer::Reset()
{
    fTreeMCVarIsPhysicalPrimary = kFALSE;

    fTreeMCVarPdgCode           = -9999;
    fTreeMCVarMotherPdgCode     = -9999;
    fTreeMCVarLabel             = -999;

    fTreeMCVarCharge            = -999;
    fTreeMCVarRap               = -999;
    fTreeMCVarEta               = -999;
    fTreeMCVarPtot              = -100;
    fTreeMCVarPt                = -100;
    fTreeMCVarPx                = -100;
    fTreeMCVarPy                = -100;
    fTreeMCVarPz                = -100;
    fTreeMCVarPhi               = -1000;
}
//_____________________________________________________________________________
void AliMCGenParticleContainer::Fill(AliMCParticle* particle, Int_t label, Int_t lMotherPDG, Bool_t IsPrimary) 
{
    fTreeMCVarIsPhysicalPrimary = IsPrimary;
    
    fTreeMCVarPdgCode           = particle->PdgCode();
    fTreeMCVarMotherPdgCode     = lMotherPDG;
    fTreeMCVarLabel             = label;
    
    fTreeMCVarCharge            = particle->Charge();
    
    fTreeMCVarRap               = particle->Y();
    fTreeMCVarEta               = particle->Eta();
    fTreeMCVarPtot              = particle->P();
    fTreeMCVarPt                = particle->Pt();
    fTreeMCVarPx                = particle->Px();
    fTreeMCVarPy                = particle->Py();
    fTreeMCVarPz                = particle->Pz();
    fTreeMCVarPhi               = particle->Phi();
    
}

//___________________________________________________________

Int_t AliMCGenParticleContainer::GetPID() const
{
  /*
   * get PID
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
  {
    if (TMath::Abs(fTreeMCVarPdgCode) == AliPID::ParticleCode(ipart))
        return ipart;
  }
  return -1;

}

//___________________________________________________________

Double_t AliMCGenParticleContainer::GetMass() const
{
  /*
   * get mass
   */

  if (GetPID() == -1)
      return 0.;
  return AliPID::ParticleMass(GetPID());
}
