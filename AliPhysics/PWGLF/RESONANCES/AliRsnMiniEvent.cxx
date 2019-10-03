//
// Mini-Event
// Contains only the useful quantities computed on the event
// which can be used for event mixing, or for direct output
// when doing analysis w.r. to multiplicity or event plane, for example.
//
// Author: A. Pulvirenti
// Developers: F. Bellini (fbellini@cern.ch)
//

#include "AliRsnMiniParticle.h"
#include "AliRsnMiniEvent.h"

ClassImp(AliRsnMiniEvent)

//__________________________________________________________________________________________________
AliRsnMiniEvent::AliRsnMiniEvent(const AliRsnMiniEvent &copy) :
   TObject(copy),
   fID(copy.fID),
   fVz(copy.fVz),
   fSpherocity(copy.fSpherocity),
   fMult(copy.fMult),
   fRefMult(copy.fRefMult),
   fTracklets(copy.fTracklets),
   fAngle(copy.fAngle),
   fQnVector(copy.fQnVector),
   fLeading(copy.fLeading),
   fParticles(copy.fParticles),
   fRef(copy.fRef),
   fRefMC(copy.fRefMC)
{
//
// Copy constructor.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//
}

//__________________________________________________________________________________________________
AliRsnMiniEvent &AliRsnMiniEvent::operator=(const AliRsnMiniEvent &copy)
{
//
// Assignment operator.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//

   TObject::operator=(copy);
   if (this == &copy)
      return *this;
   fID = copy.fID;
   fVz = copy.fVz;
   fSpherocity=copy.fSpherocity;
   fMult = copy.fMult;
   fRefMult = copy.fRefMult;
   fTracklets = copy.fTracklets;
   fAngle = copy.fAngle;
   fQnVector = copy.fQnVector;
   fLeading = copy.fLeading;
   fParticles = copy.fParticles;
   fRef = copy.fRef;
   fRefMC = copy.fRefMC;
   return (*this);
}

//__________________________________________________________________________________________________
void AliRsnMiniEvent::Clear(Option_t *)
{
//
// Clears event
//
    fID = 0;
    fVz = 0;
    fSpherocity=0;
    fMult = 0;
    fRefMult = 0;
    fTracklets = 0;
    fAngle = 0;
    fQnVector = 0;
    fLeading = -1;
    fRef = 0;
    fRefMC = 0;

    fParticles.Clear("C");
}

//__________________________________________________________________________________________________
AliRsnMiniParticle *AliRsnMiniEvent::AddParticle()
{
//
// Add a new particle to the list and returns a pointer to it,
// in order to allow to se its parameters.
//

   Int_t n = fParticles.GetEntries();
   return (AliRsnMiniParticle *)fParticles.ConstructedAt(n);
}

//__________________________________________________________________________________________________
void AliRsnMiniEvent::AddParticle(AliRsnMiniParticle copy)
{
//
// Add a new particle to the list and returns a pointer to it,
// in order to allow to se its parameters.
//

   Int_t n = fParticles.GetEntries();
   new (fParticles[n]) AliRsnMiniParticle(copy);
}

//__________________________________________________________________________________________________
AliRsnMiniParticle *AliRsnMiniEvent::GetParticle(Int_t i)
{
//
// Return the particle
//

   if (i < 0 || i > fParticles.GetEntriesFast()) return 0x0;

   return (AliRsnMiniParticle *)fParticles[i];
}

//__________________________________________________________________________________________________
AliRsnMiniParticle *AliRsnMiniEvent::LeadingParticle(Bool_t mc)
{
//
// Return the leading particle
//

   if (fLeading == -1 ) SelectLeadingParticle(mc);
   if (fLeading == -2) return 0x0;
   if (fLeading >= fParticles.GetEntriesFast()) return 0x0;

   return (AliRsnMiniParticle *)fParticles[fLeading];
}


//_____________________________________________________________________________
void AliRsnMiniEvent::SelectLeadingParticle(Bool_t mc)
{
//
// Searches for leading particle (particle with maximum Pt in event)
//

   Double_t ptMax = 0.0;
   Double_t pt;
   Int_t i;
   AliRsnMiniParticle *part = 0x0;
   fLeading = -2;
   for (i = 0; i < fParticles.GetEntriesFast(); i++) {
      part = (AliRsnMiniParticle *)fParticles[i];
      pt = TMath::Sqrt(TMath::Power(part->Px(mc),2)+TMath::Power(part->Py(mc),2));
      if (pt > ptMax) {
         ptMax = pt;
         fLeading = i;
      }
   }
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniEvent::CountParticles(TArrayI &found, Char_t charge, Int_t cutID)
{
//
// Counts how many particles have the specified charge and cut bit
// if charge is not '+', '-' or '0', all charges are considered
// if cut bit is < 0, it is not checked
//

   Int_t i, npart = fParticles.GetEntriesFast();
   Int_t    count = 0;
   AliRsnMiniParticle *part = 0x0;

   found.Set(npart);

   for (i = 0; i < npart; i++) {
      part = (AliRsnMiniParticle *)fParticles[i];
      if (charge == '+' || charge == '-' || charge == '0') {
         if (part->Charge() != charge) continue;
      }
      if (cutID >= 0) {
         if (!part->HasCutBit(cutID)) continue;
      }
      found[count] = i;
      count++;
   }

   found.Set(count);
   return count;
}
