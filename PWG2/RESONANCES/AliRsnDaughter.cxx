//
// Class AliRsnDaughter
//
// Interface to candidate daughters of a resonance (tracks).
// Points to the source of information, which is generally an AliVParticle-derived object
// and contains few internal data-members to store "on fly" some important information
// for the computations required during resonance analysis.
// It contains a useful TLorentzVector data-member which, provided that a meaningful mass was assigned,
// eases a lot the computation of invariant masses from summing up several of these objects.
//
// authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          M. Vala (martin.vala@cern.ch)
//

#include <TParticle.h>
#include "AliAODVertex.h"
#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter() :
  fOK(kFALSE),
  fLabel(-1),
  fMotherPDG(0),
  fPrec(0.0, 0.0, 0.0, 0.0),
  fPsim(0.0, 0.0, 0.0, 0.0),
  fRef(0x0),
  fRefMC(0x0)
{
//
// Default constructor.
//
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(const AliRsnDaughter &copy) :
  TObject(copy),
  fOK(copy.fOK),
  fLabel(copy.fLabel),
  fMotherPDG(copy.fMotherPDG),
  fPrec(copy.fPrec),
  fPsim(copy.fPsim),
  fRef(copy.fRef),
  fRefMC(copy.fRefMC)
{
//
// Copy constructor.
// Pointers are NOT duplicated, since they don't come from a 'new'
// statement, but from just referencing something in the data source.
//
}

//_____________________________________________________________________________
AliRsnDaughter& AliRsnDaughter::operator=(const AliRsnDaughter &copy)
{
//
// Assignment operator.
// Pointers are NOT duplicated, since they don't come from a 'new'
// statement, but from just referencing something in the data source.
//

  fOK        = copy.fOK;
  fLabel     = copy.fLabel;
  fMotherPDG = copy.fMotherPDG;
  fPrec      = copy.fPrec;
  fPsim      = copy.fPsim;
  fRef       = copy.fRef;
  fRefMC     = copy.fRefMC;

  return (*this);
}

//_____________________________________________________________________________
void AliRsnDaughter::Reset()
{
//
// Reset this track to meaningless values and to a 'bad' status.
// After this has been done, this object should not be used
// for analysis unless initialized properly.
//

  fOK        = kFALSE;
  fLabel     = -1;
  fMotherPDG = 0;
  fRef       = 0x0;
  fRefMC     = 0x0;
  
  fPrec.SetXYZM(0.0, 0.0, 0.0, 0.0);
  fPsim.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetPDG(Bool_t abs)
{
//
// Return the PDG code of the particle from MC ref (if any).
// If argument is kTRUE, returns its absolute value.
//

  Int_t pdg = 0;

  // ESD
  AliMCParticle *esd = GetRefMCESD();
  if (esd) pdg = esd->Particle()->GetPdgCode();
  
  // AOD
  AliAODMCParticle *aod = GetRefMCAOD();
  if (aod) pdg = aod->GetPdgCode();
  
  // abs value if required
  if (abs) pdg = TMath::Abs(pdg);
  return pdg;
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetID()
{
//
// Return reference index, using the "GetID" method
// of the possible source object.
// In case of V0s, since this method is unsuccessful, return the label.
//

  // ESD tracks
  AliESDtrack *esd = GetRefESDtrack();
  if (esd) return esd->GetID();

  // AOD tracks
  AliAODTrack *aod = GetRefAODtrack();
  if (aod) return aod->GetID();

  // whatever else
  return GetLabel();
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::HasFlag(ULong_t flag)
{
//
// Checks that the 'status' flag of the source object has one or 
// a combination of status flags specified in argument.
// Works only with track-like objects, irrespectively if they
// are ESD or AOD tracks, since it refers to their AliVTrack base class.
//

  AliVTrack *track  = dynamic_cast<AliVTrack*>(fRef);
  if (!track) return kFALSE;
  
  ULong_t status = (ULong_t)track->GetStatus();
  
  return ((status & flag) != 0);
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::SetMass(Double_t mass)
{
//
// Assign a mass hypothesis to the track.
// This causes the 4-momentum data members to be initialized
// using the momenta of referenced tracks/v0s and this mass.
// This step is fundamental for the following of the analysis.
// Of course, this operation is successful only if the mass is
// a meaningful positive number, and the pointers are properly initialized.
//

  if (mass < 0.) return kFALSE;
  
  if (fRef)   fPrec.SetXYZM(fRef  ->Px(), fRef  ->Py(), fRef  ->Pz(), mass);
  if (fRefMC) fPsim.SetXYZM(fRefMC->Px(), fRefMC->Py(), fRefMC->Pz(), mass);
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::IsKinkDaughter()
{
//
// Checks if this track is a kink daughter.
// this information is important for some cuts, in some cases
// and it is retrieved differently from ESDs and AODs, so
// this is done here in order to have a unique outcome.
//

  AliESDtrack *etrack = GetRefESDtrack();
  AliAODTrack *atrack = GetRefAODtrack();
  
  if (etrack)
  {
    return (etrack->GetKinkIndex(0) > 0);
  }
  else if (atrack)
  {
    AliAODVertex *vertex = atrack->GetProdVertex();
    if (vertex) if (vertex->GetType() == AliAODVertex::kKink) return kTRUE;
  }
  
  return kFALSE;
}
