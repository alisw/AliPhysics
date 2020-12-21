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
   fResMass(0.0),
   fResPDG(0),
   fPairCuts(0x0),
   fPair(),
   fSelRes(0),
   fSel1(0),
   fSel2(0),
   fPairMode(0),
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
   fResMass(copy.fResMass),
   fResPDG(copy.fResPDG),
   fPairCuts(copy.fPairCuts),
   fPair(),
   fSelRes(0),
   fSel1(0),
   fSel2(0),
   fPairMode(copy.fPairMode),
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
   fResMass = copy.fResMass;
   fResPDG = copy.fResPDG;
   fPairCuts = copy.fPairCuts;
   fPairMode = copy.fPairMode;
   fEvent = copy.fEvent;

   Int_t i;
   for (i = 0; i < 2; i++) {
      fCutID[i] = copy.fCutID[i];
      fDaughter[i] = copy.fDaughter[i];
      fCharge[i] = copy.fCharge[i];
   }

   fSelRes.Set(0);
   fSel1.Set(0);
   fSel2.Set(0);

   return (*this);
}

//__________________________________________________________________________________________________
AliRsnMiniResonanceFinder::~AliRsnMiniResonanceFinder()
{
//
// Destructor.
//

  return;
}

/*
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
 */

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
   AliRsnMiniParticle *r, *p1, *p2;

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
         fPair.Fill(p1, p2, GetMass(0), GetMass(1), fResMass);

         // check pair against cuts
         if (fPairCuts && !fPairCuts->IsSelected(&fPair)) continue;
         // package the pair as a mini particle

         if(fPairMode==1 && (p1->Mother()<0 || p1->Mother()!=p2->Mother() || !AliRsnDaughter::IsEquivalentPDGCode(p1->MotherPDG(),GetResonancePDG()) ) ) continue;// use only true pairs (MC)
         if(fPairMode==2 && p1->Mother()>=0 && p1->Mother()==p2->Mother()) continue;// use only false pairs (MC)

         r = event->AddParticle();
         r->Clear();
         r->SetResonance();
         r->Charge() = '0';
         r->SetCutBit(fCutIDrsn);

         if(p1->Charge()=='+'){
            if(p2->Charge()=='-' || p2->Charge()=='+'){
               r->IndexV0Pos() = p1->Index();
               r->IndexV0Neg() = p2->Index();
               r->Charge() = '0'; // does not account for double charges, e.g., Delta++
            }else if(p2->Charge()=='0'){
               r->IndexV0Pos() = p2->IndexV0Pos();
               r->IndexV0Neg() = p2->IndexV0Neg();
               r->IndexBachelor() = p1->Index();
               r->Charge() = '+';
            }
         }else if(p1->Charge()=='-'){
            if(p2->Charge()=='-' || p2->Charge()=='+'){
               r->IndexV0Pos() = p2->Index();
               r->IndexV0Neg() = p1->Index();
               r->Charge() = '0'; // does not account for double charges, e.g., anti-Delta--
            }else if(p2->Charge()=='0'){
               r->IndexV0Pos() = p2->IndexV0Pos();
               r->IndexV0Neg() = p2->IndexV0Neg();
               r->IndexBachelor() = p1->Index();
               r->Charge() = '-';
            }
         }else if(p1->Charge()=='0'){
            if(p2->Charge()=='-' || p2->Charge()=='+'){
               r->IndexV0Pos() = p1->IndexV0Pos();
               r->IndexV0Neg() = p1->IndexV0Neg();
               r->IndexBachelor() = p2->Index();
               r->Charge() = p2->Charge();
            }
         }

         r->PrecX() = fPair.Sum(0).X();
         r->PrecY() = fPair.Sum(0).Y();
         r->PrecZ() = fPair.Sum(0).Z();
         r->StoredMass(0) = fPair.InvMass(0);

         if(p1->Mother()>=0 && p1->Mother()==p2->Mother()){
            r->Index() = p1->Mother();
            r->PsimX() = p1->PmotherX();
            r->PsimY() = p1->PmotherY();
            r->PsimZ() = p1->PmotherZ();
            r->PDG() = p1->MotherPDG();
         }else{
            r->PsimX() = fPair.Sum(1).X();
            r->PsimY() = fPair.Sum(1).Y();
            r->PsimZ() = fPair.Sum(1).Z();
         }
         r->StoredMass(1) = fPair.InvMass(1);

         r->PmotherX() = r->PmotherY() = r->PmotherZ() = 0.0;// filled later

         nadded++;
      } // end internal loop
   } // end external loop

   event->CountParticles(fSelRes, '0', fCutIDrsn);

   AliDebugClass(1, Form("Pairs added in total = %4d", nadded));
   return nadded;
}

//__________________________________________________________________________________________________
void AliRsnMiniResonanceFinder::FillMother(AliMCEvent* event, AliRsnMiniEvent* mini)
{
//
// Loops over MC tracks to find matches to resonances that have been reconatructed and selected.
// Adds information about their mothers to the mini event.
// ESD version

   if (!event) {
      AliError("Missing AliMCEvent");
      return;
   }

   if (!mini) {
      AliError("Missing AliRsnMiniEvent");
      return;
   }

   AliMCParticle *p,*mother;
   AliRsnMiniParticle* r;
   Int_t j,k,imother, npart = event->GetNumberOfTracks();

   for (j=0; j<npart; j++) {
      p = (AliMCParticle*) event->GetTrack(j);
      if (!p) continue;
      if (!AliRsnDaughter::IsEquivalentPDGCode(p->Particle()->GetPdgCode() , GetResonancePDG())) continue;

      r = 0;
      imother = -2;
      mother = 0;

      for (k=0; k<fSelRes.GetSize(); k++){
         r = mini->GetParticle(fSelRes[k]);
         if(!r || r->Index()!=j) continue;

         imother=p->GetMother();
         if(imother<0) continue;
         mother = dynamic_cast<AliMCParticle*>(event->GetTrack(imother));
         if(!mother){
            AliError("Failed casting the mother particle!");
            continue;
         }

         r->Mother() = imother;
         r->PmotherX() = mother->Px();
         r->PmotherY() = mother->Py();
         r->PmotherZ() = mother->Pz();
         r->MotherPDG() = mother->PdgCode();
         break;
      }
   }
    
   return;
}

//__________________________________________________________________________________________________
void AliRsnMiniResonanceFinder::FillMother(TClonesArray* event, AliRsnMiniEvent* mini)
{
//
// Loops over MC tracks to find matches to resonances that have been reconatructed and selected.
// Adds information about their mothers to the mini event.
// AOD version

   if (!event) {
      AliError("Missing AOD Event");
      return;
   }

   if (!mini) {
      AliError("Missing AliRsnMiniEvent");
      return;
   }

   AliAODMCParticle *p,*mother;
   AliRsnMiniParticle* r;
   Int_t j,k,imother, npart = event->GetEntries();

   for (j=0; j<npart; j++) {
      p = (AliAODMCParticle*) event->At(j);
      if (!p) continue;
      if (!AliRsnDaughter::IsEquivalentPDGCode(p->GetPdgCode() , GetResonancePDG())) continue;

      r = 0;
      imother = -2;
      mother = 0;

      for (k=0; k<fSelRes.GetSize(); k++){
         r = mini->GetParticle(fSelRes[k]);
         if(!r || r->Index()!=j) continue;

         imother=p->GetMother();
         if(imother<0) continue;
         mother = dynamic_cast<AliAODMCParticle*>(event->At(imother));
         if(!mother){
            AliError("Failed casting the mother particle!");
            continue;
         }
          
         r->Mother() = imother;
         r->PmotherX() = mother->Px();
         r->PmotherY() = mother->Py();
         r->PmotherZ() = mother->Pz();
         r->MotherPDG() = mother->PdgCode();
         break;
      }
   }
    
   return;
}
