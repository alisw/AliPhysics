//
// Class to identify resonance candidate pairs, which are then
// added back into the mini event to serve as daughters for
// another particle.
//

#include "Riostream.h"

#include "TMath.h"
#include "TString.h"

#include "AliRsnMiniParticle.h"
#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniAnalysisTask.h"

#include "AliLog.h"
#include "AliRsnCutSet.h"

ClassImp(AliRsnMiniResonanceFinder)

//__________________________________________________________________________________________________
AliRsnMiniResonanceFinder::AliRsnMiniResonanceFinder(const char *name) :
   TNamed(name, ""),
   fCutIDrsn(0),
   fMotherMass(0.0),
   fPairCuts(0x0),
   fPair(),
   fSel1(0),
   fSel2(0),
   fEvent(0x0)
{
//
// Constructor
//

   fCutID[0] = fCutID[1] = -1;
   fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown;
   fCharge[0] = fCharge[1] = 0;
}

//__________________________________________________________________________________________________
AliRsnMiniResonanceFinder::AliRsnMiniResonanceFinder(const AliRsnMiniResonanceFinder& copy) :
   TNamed(copy),
   fCutIDrsn(copy.fCutIDrsn),
   fMotherMass(copy.fMotherMass),
   fPairCuts(copy.fPairCuts),
   fPair(),
   fSel1(0),
   fSel2(0),
   fEvent(copy.fEvent)
{
//
// Copy constructor
//
  
   Int_t i;
   for (i = 0; i < 2; i++) {
      fCutID[i] = copy.fCutID[i];
      fDaughter[i] = copy.fDaughter[i];
      fCharge[i] = copy.fCharge[i];
   }
}

//__________________________________________________________________________________________________
AliRsnMiniResonanceFinder &AliRsnMiniResonanceFinder::operator=(const AliRsnMiniResonanceFinder &copy)
{
//
// Assignment operator
//
   if (this == &copy)
      return *this;
   fCutIDrsn = copy.fCutIDrsn;
   fMotherMass = copy.fMotherMass;
   fPairCuts = copy.fPairCuts;
   fEvent = copy.fEvent;

   Int_t i;
   for (i = 0; i < 2; i++) {
      fCutID[i] = copy.fCutID[i];
      fDaughter[i] = copy.fDaughter[i];
      fCharge[i] = copy.fCharge[i];
   }

   fSel1.Set(0);
   fSel2.Set(0);

   return (*this);
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniResonanceFinder::ConnectTo(AliRsnMiniAnalysisTask* task)
{
  //
  // IMPORTANT: Do not assign any additional track cuts to the AliRsnMiniAnalysisTask
  // after calling this ConnectTo function.
  //
  task->SetResonanceFinder(this);
  fCutIDrsn = task->GetNumberOfTrackCuts();
  return fCutIDrsn;// plug this value into AliRsnMiniOutput objects
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniResonanceFinder::RunResonanceFinder(AliRsnMiniEvent* event)
{
//
// Loops on the passed mini-event, and for each pair of particles
// which satisfy the charge and cut requirements defined here, add an entry.
// Returns the number of successful fillings.
//

   // loop variables
   Int_t i1, i2, start, nadded = 0;
   AliRsnMiniParticle *mother, *p1, *p2;

   // it is necessary to know if criteria for the two daughters are the same
   Bool_t sameCriteria = ((fCharge[0] == fCharge[1]) && (fDaughter[0] == fDaughter[1]));

   TString selList1  = "";
   TString selList2  = "";
   Int_t   n1 = event->CountParticles(fSel1, fCharge[0], fCutID[0]);
   Int_t   n2 = event->CountParticles(fSel2, fCharge[1], fCutID[1]);
   for (i1 = 0; i1 < n1; i1++) selList1.Append(Form("%d ", fSel1[i1]));
   for (i2 = 0; i2 < n2; i2++) selList2.Append(Form("%d ", fSel2[i2]));
   AliDebugClass(1, Form("[%10s] Part #1: -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", GetName(), event->ID(), fCharge[0], fCutID[0], n1, selList1.Data()));
   AliDebugClass(1, Form("[%10s] Part #2: -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", GetName(), event->ID(), fCharge[1], fCutID[1], n2, selList2.Data()));
   if (!n1 || !n2) {
      AliDebugClass(1, "No pairs to mix");
      return 0;
   }

   // external loop
   for (i1 = 0; i1 < n1; i1++) {
      p1 = event->GetParticle(fSel1[i1]);
      // define starting point for inner loop
      // if daughter selection criteria (charge, cuts) are the same
      // and the two events coincide, internal loop must start from
      // the first track *after* current one;
      // otherwise it starts from the beginning
      start = (sameCriteria ? i1 + 1 : 0);
      AliDebugClass(2, Form("Start point = %d", start));
      // internal loop
      for (i2 = start; i2 < n2; i2++) {
         p2 = event->GetParticle(fSel2[i2]);
         // avoid to mix a particle with itself
         if (p1->Index() == p2->Index()) {
            AliDebugClass(2, "Skipping same index");
            continue;
         }
         // sum momenta
         fPair.Fill(p1, p2, GetMass(0), GetMass(1), fMotherMass);

         // check pair against cuts
         if (fPairCuts) {
            if (!fPairCuts->IsSelected(&fPair)) continue;
         }
         // package the pair as a mini particle

         mother = event->AddParticle();
	     mother->Clear();
	     mother->Index() = -1;
	     mother->IndexV0Pos() = p1->Index();
	     mother->IndexV0Neg() = p2->Index();
	     mother->Charge() = '0';
	     mother->SetCutBit(fCutIDrsn);

	     mother->PrecX() = fPair.Sum(0).X();
	     mother->PrecY() = fPair.Sum(0).Y();
	     mother->PrecZ() = fPair.Sum(0).Z();
          
         if(p1->Mother() == p2->Mother()){
            mother->PsimX() = p1->PmotherX();
            mother->PsimY() = p1->PmotherY();
            mother->PsimZ() = p1->PmotherZ();
            mother->PDG() = p1->MotherPDG();
         }else{
            mother->PsimX() = fPair.Sum(1).X();
            mother->PsimY() = fPair.Sum(1).Y();
            mother->PsimZ() = fPair.Sum(1).Z();
         }

         nadded++;
      } // end internal loop
   } // end external loop

   AliDebugClass(1, Form("Pairs added in total = %4d", nadded));
   return nadded;
}
