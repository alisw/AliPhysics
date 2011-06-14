//
// This object is used as lightweight temporary container
// of all information needed from any input object and 
// useful for resonance analysis.
// Lists of such objects are stored in a buffer, in order
// to allow an event mixing.
//

#include "AliRsnDaughter.h"
#include "AliRsnMiniParticle.h"

ClassImp(AliRsnMiniParticle)

//__________________________________________________________________________________________________
void AliRsnMiniParticle::CopyDaughter(AliRsnDaughter *daughter)
{
//
// Sets data members from the passed object
//

   // reset what could not be initialized
   fPDG = 0;
   fMother = -1;
   fMotherPDG = 0;
   fCutBits = 0x0;
   fPsim[0] = fPrec[0] = fPsim[1] = fPrec[1] = fPsim[2] = fPrec[2] = 0.0;

   // charge
   if (daughter->IsPos())
      fCharge = '+';
   else if (daughter->IsNeg())
      fCharge = '-';
   else
      fCharge = '0';
   
   // rec info
   if (daughter->GetRef()) {
      fPrec[0] = daughter->GetRef()->Px();
      fPrec[1] = daughter->GetRef()->Py();
      fPrec[2] = daughter->GetRef()->Pz();
   }
   
   // MC info
   if (daughter->GetRefMC()) {
      fPsim[0] = daughter->GetRefMC()->Px();
      fPsim[1] = daughter->GetRefMC()->Py();
      fPsim[2] = daughter->GetRefMC()->Pz();
      fPDG = daughter->GetPDG();
      fMother = daughter->GetMother();
      fMotherPDG = daughter->GetMotherPDG();
   }
}
