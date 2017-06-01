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
#include "AliMCEvent.h"
#include "AliRsnEvent.h"
#include "AliRsnMiniEvent.h"
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
   fNTotSisters = -1;
   fIsFromB = kFALSE;
   fIsQuarkFound = kFALSE;
   fCutBits = 0x0;
   fPsim[0] = fPrec[0] = fPmother[0] = fPsim[1] = fPrec[1] = fPmother[1] = fPsim[2] = fPrec[2] = fPmother[2] = 0.0;
   fIndexV0[0] = fIndexV0[1] = -1;

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
   if (!event) {
     AliWarning("Invalid reference event: cannot copy DCA nor Nsisters.");
     return;
   }
   if (event->IsAOD()){
     // DCA to Primary Vertex for AOD
     AliAODTrack *track = daughter->Ref2AODtrack();  
     AliAODv0 *v0 = daughter->Ref2AODv0();     
     AliAODEvent *aodEvent = (AliAODEvent*) event->GetRefAOD();
     if (track && aodEvent) {
       AliVVertex *vertex = (AliVVertex*) aodEvent->GetPrimaryVertex();
       Double_t b[2], cov[3]; 
       if (vertex) {
	 if ( !((track->GetStatus() & AliESDtrack::kTPCin) == 0) && !((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) && !((track->GetStatus() & AliESDtrack::kITSrefit) == 0) ){
	   if (track->PropagateToDCA(vertex, aodEvent->GetMagneticField(), kVeryBig, b, cov))
	     fDCA = b[0];
	 }
       } 
     }
      if (v0 && aodEvent) {
          fIndexV0[0] = v0->GetPosID();
          fIndexV0[1] = v0->GetNegID();
          // Printf("!!!!!!!! RSN index=%d v0Pos=%d v0Neg=%d", fIndex, fIndexV0[0], fIndexV0[1]);
       }
     // Number of Daughters from MC and Momentum of the Mother
     if (event->GetRefMC()) {
       TClonesArray * list = event->GetAODList();
       AliAODMCParticle *part = (AliAODMCParticle *)list->At(fMother);
       if (part) {
	 fNTotSisters = part->GetNDaughters();
	 fPmother[0]  = part->Px();
	 fPmother[1]  = part->Py();
	 fPmother[2]  = part->Pz();
	 Int_t istep = 0;
	 Int_t pdgGranma = 0;
	 Int_t abspdgGranma =0;
	 Int_t mother_temp = daughter->GetMother();
	 while (mother_temp >=0 ){
	       istep++;
	       AliDebug(2,Form("mother at step %d = %d", istep, mother_temp));
	       AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(list->At(mother_temp));
	       if (mcGranma){
	               pdgGranma = mcGranma->GetPdgCode();
	               AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
	               abspdgGranma = TMath::Abs(pdgGranma);
	               if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	        	 fIsFromB=kTRUE;
	               }
	               if(abspdgGranma==4 || abspdgGranma==5) fIsQuarkFound=kTRUE;
	               mother_temp = mcGranma->GetMother();
	       }else{
	               AliError("Failed casting the mother particle!");
	               break;
	       }
	 }
       }
     }
   } else {
     if (event->IsESD()){
       //DCA to Primary Vertex for ESD
       AliESDtrack *track = daughter->Ref2ESDtrack(); 
       AliESDv0 *v0 = daughter->Ref2ESDv0();    
       AliESDEvent *esdEvent = (AliESDEvent*) event->GetRefESD();
       if (track && esdEvent) {
	 AliVVertex *vertex = (AliVVertex*) esdEvent->GetPrimaryVertex();
	 Double_t b[2], cov[3]; 
	 if (vertex) {
	   if ( !((track->GetStatus() & AliESDtrack::kTPCin) == 0) && !((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) && !((track->GetStatus() & AliESDtrack::kITSrefit) == 0) ){
	     if (track->PropagateToDCA(vertex, esdEvent->GetMagneticField(), kVeryBig, b, cov))
	       fDCA = b[0];
	   }
	 }      
       }
       if (v0 && esdEvent) {
          fIndexV0[0] = v0->GetPindex();
          fIndexV0[1] = v0->GetNindex();
       }
       // Number of Daughters from MC and Momentum of the Mother
       if (event->GetRefMC()) {
	 AliMCParticle *part = (AliMCParticle *)event->GetRefMC()->GetTrack(fMother);
	 AliMCEvent * MCEvent = event->GetRefMCESD();
	 if(part){
	   fNTotSisters = part->Particle()->GetNDaughters();
	   fPmother[0]  = part->Px();
	   fPmother[1]  = part->Py();
	   fPmother[2]  = part->Pz();
	   Int_t istep = 0;
	   Int_t pdgGranma = 0;
	   Int_t abspdgGranma =0;
	   Int_t mother_temp = daughter->GetMother();
	   while (mother_temp >=0 ){
		 istep++;
		 AliDebug(2,Form("mother at step %d = %d", istep, mother_temp));
		 AliMCParticle* mcGranma = dynamic_cast<AliMCParticle*>(MCEvent->GetTrack(mother_temp));
		 if (mcGranma){
		 	 pdgGranma = mcGranma->PdgCode();
		 	 AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
		 	 abspdgGranma = TMath::Abs(pdgGranma);
		 	 if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
		 	   fIsFromB=kTRUE;
		 	 }
		 	 if(abspdgGranma==4 || abspdgGranma==5) fIsQuarkFound=kTRUE;
		 	 mother_temp = mcGranma->GetMother();
		 }else{
		 	 AliError("Failed casting the mother particle!");
		 	 break;
		 }
	   }
	 }
       }
     }
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
