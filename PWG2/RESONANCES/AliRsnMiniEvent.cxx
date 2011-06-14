//
// Mini-Event
// Contains only the useful quantities computed on the event
// which can be used for event mixing, or for direct output
// when doing analysis w.r. to multiplicity or event plane, for example.
//
// Author: A. Pulvirenti
//

#include "AliRsnMiniParticle.h"
#include "AliRsnMiniEvent.h"

ClassImp(AliRsnMiniEvent)

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
AliRsnMiniParticle* AliRsnMiniEvent::LeadingParticle()
{
//
// Return the leading particle
//

   if (fLeading < 0) return 0x0;
   if (fLeading >= fParticles.GetEntriesFast()) return 0x0;
   
   return (AliRsnMiniParticle*)fParticles[fLeading];
}
