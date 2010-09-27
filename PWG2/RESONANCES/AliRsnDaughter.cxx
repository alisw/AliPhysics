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

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"

#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter() :
  fOK(kFALSE),
  fLabel(-1),
  fMotherPDG(0),
  fP(0.0, 0.0, 0.0, 0.0),
  fPMC(0.0, 0.0, 0.0, 0.0),
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
  fP(copy.fP),
  fPMC(copy.fPMC),
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
  fP         = copy.fP;
  fPMC       = copy.fPMC;
  fRef       = copy.fRef;
  fRefMC     = copy.fRefMC;

  return (*this);
}

//_____________________________________________________________________________
AliRsnDaughter::~AliRsnDaughter()
{
//
// Destructor.
// Since pointers do not allocate new objects, no 'delete' is done.
//
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
  
  fP  .SetXYZM(0.0, 0.0, 0.0, 0.0);
  fPMC.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

//_____________________________________________________________________________
void AliRsnDaughter::Print(Option_t * const /*option*/) const
{
//
// Prints the values of data members, using the options:
// - P --> momentum
// - V --> DCA vertex
// - C --> electric charge
// - F --> flags
// - I --> identification (PID, probability and mass)
// - W --> PID weights
// - M --> Montecarlo
// - L --> index & label
// - A --> angles
// - ALL --> All oprions switched on
//
// Index and label are printed by default.
//
  
  Char_t type[50], source[50];
  if (IsTrack()) snprintf(type, 5, "track");
  else if (IsV0()) snprintf(type, 2, "V0");
  else snprintf(type, 4, "none");
  
  if (IsESD()) snprintf(source, 3, "ESD");
  else if (IsAOD()) snprintf(source, 3, "AOD");
  else if (fRefMC != 0x0) snprintf(source, 7, "MC only");
  else snprintf(source, 4, "none");
  
  AliInfo("===== ALIRSNDAUGHTER INFORMATION ==============================================");
  AliInfo(Form(".......Index                  : %d", GetID()));
  AliInfo(Form(".......Label                  : %d", GetLabel()));
  AliInfo(Form(".......Type, source           : %s %s", type, source));
  AliInfo(Form(".......Charge                 : %c", ChargeChar()));
  
  if (fRef)
  {
    AliInfo(Form(".......Px, Py, Pz, Pt (ref)   : %f %f %f %f", fP.X(), fP.Y(), fP.Z(), fP.Perp()));
  } else AliInfo("....... absent REF");
  
  if (fRefMC) 
  {
    AliInfo(Form(".......Px, Py, Pz, Pt (ref MC): %f %f %f %f", fP.X(), fP.Y(), fP.Z(), fP.Perp())); 
    AliInfo(Form(".......PDG code               : %d", fRefMC->Particle()->GetPdgCode()));
    AliInfo(Form(".......Mother (label)         : %d", fRefMC->Particle()->GetFirstMother()));
    AliInfo(Form(".......Mother (PDG code)      : %d", fMotherPDG));
  } else AliInfo("....... absent REF MC");
  
  AliInfo("===== END ALIRSNDAUGHTER INFORMATION ==========================================");
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetID() const
{
//
// Return reference index, using the "GetID" method
// of the possible source object.
// In case of V0s, since this method is unsuccessful, return the label.
//

  const AliESDtrack *esd = GetRefESDtrack();
  if (esd) return esd->GetID();

  const AliAODTrack *aod = GetRefAODtrack();
  if (aod) return aod->GetID();

  return GetLabel();
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::HasFlag(ULong_t flag) const
{
//
// Checks that the 'status' flag of the source object has one or 
// a combination of status flags combined with the bitwise OR.
// Works only with track-like objects.
//

  if (IsTrack())
  {
    AliVTrack *track  = dynamic_cast<AliVTrack*>(fRef);
    ULong_t    status = (ULong_t)track->GetStatus();
    return ((status & flag) != 0);
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::SetMass(Double_t mass)
{
//
// The argument defines a mass to be assigned to the 4-momentum.
// This value is not stored as a data-member, but serves just to modify
// the energy of the 4-momentum both in the MC and non-MC objects,
// while the vector momentum is taken from the respective references.
// The method returns a BOOL outcome which is kFALSE in the case that
// mass is meaningless or when both references are NULL or the goodness 
// flag is not positive.
//

  if (!fOK) return kFALSE;
  if (mass < 0.) return kFALSE;
  if (!fRef && !fRefMC) return kFALSE;
  
  if (fRef)   fP  .SetXYZM(fRef  ->Px(), fRef  ->Py(), fRef  ->Pz(), mass);
  if (fRefMC) fPMC.SetXYZM(fRefMC->Px(), fRefMC->Py(), fRefMC->Pz(), mass);
  
  return kTRUE;
}
