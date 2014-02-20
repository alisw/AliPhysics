//
// This object is used as lightweight temporary container
// of all information needed from any input object and
// useful for resonance analysis.
// Lists of such objects are stored in a buffer, in order
// to allow an event mixing.
//

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "AliAODEvent.h"
#include "AliRsnEvent.h"
#include "AliRsnMiniParticle.h"

ClassImp(AliRsnMiniParticle)

//__________________________________________________________________________________________________
void AliRsnMiniParticle::CopyDaughter(AliRsnDaughter *daughter)
{
//
// Sets data members from the passed object
//

   // reset what could not be initialized
   fDCA = 0.0;  //typically used for D0 analysis
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
   
   AliRsnEvent *event = (AliRsnEvent *) daughter->GetOwnerEvent();
   if (event && event->IsAOD()){
     AliAODTrack *track = (AliAODTrack*) daughter->Ref2AODtrack();   
     AliAODEvent *aodEvent = (AliAODEvent*) event->GetRefAOD();
     if (track && aodEvent) {
       AliVVertex *vertex = (AliVVertex*) aodEvent->GetPrimaryVertex();
       Double_t b[2], cov[3]; 
       if (vertex) {
	 track->PropagateToDCA(vertex, aodEvent->GetMagneticField(), kVeryBig, b, cov); 
	 fDCA = b[0];
       }
     }
   } else {
     AliWarning("DCA not implemented for ESDs");
   }
   
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniParticle::Mass()
{
   //
   // return mass of particle
   //

   TDatabasePDG *db   = TDatabasePDG::Instance();
   TParticlePDG *part = db->GetParticle(PDG());
   return part->Mass();
}

//__________________________________________________________________________________________________
void AliRsnMiniParticle::Set4Vector(TLorentzVector &v, Float_t mass, Bool_t mc)
{
   //
   // return 4 vector of particle
   //

   if (mass<0.0) mass = Mass();
   v.SetXYZM(Px(mc), Py(mc), Pz(mc),mass);
}
