//
// Computator for single daughters.
// Implements a simple loop on tracks from one of the entry lists
// filled by the task AliRsnInputHandler, adding a check on their
// definition specified in the daughter def.
// Author: A. Pulvirenti
//

#include <Riostream.h>
#include <TEntryList.h>

#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughterDef.h"
#include "AliRsnDaughterSelector.h"

#include "AliRsnLoopDaughter.h"

ClassImp(AliRsnLoopDaughter)

//_____________________________________________________________________________
AliRsnLoopDaughter::AliRsnLoopDaughter(const char *name, Int_t listID, AliRsnDaughterDef *def) :
   AliRsnLoop(name),
   fTrueMC(kFALSE),
   fOnlyTrue(kFALSE),
   fUseMCRef(kFALSE),
   fListID(listID),
   fDef(def),
   fDaughter()
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnLoopDaughter::AliRsnLoopDaughter(const AliRsnLoopDaughter& copy) :
   AliRsnLoop(copy),
   fTrueMC(copy.fTrueMC),
   fOnlyTrue(copy.fOnlyTrue),
   fUseMCRef(copy.fUseMCRef),
   fListID(copy.fListID),
   fDef(copy.fDef),
   fDaughter(copy.fDaughter)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnLoopDaughter& AliRsnLoopDaughter::operator=(const AliRsnLoopDaughter& copy)
{
//
// Assignment operator
//

   AliRsnLoop::operator=(copy);
   fTrueMC = copy.fTrueMC;
   fOnlyTrue = copy.fOnlyTrue;
   fUseMCRef = copy.fUseMCRef;
   fListID = copy.fListID;
   fDaughter = copy.fDaughter;
   fDef = copy.fDef;

   return (*this);
}

//_____________________________________________________________________________
AliRsnLoopDaughter::~AliRsnLoopDaughter()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnLoopDaughter::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//
}

//_____________________________________________________________________________
Bool_t AliRsnLoopDaughter::Init(const char *prefix, TList *list)
{
//
// Initialization function.
// Loops on all functions and eventual the ntuple, to initialize output objects.
//

   return AliRsnLoop::Init(Form("%s_%s", prefix, GetName()), list);
}

//_____________________________________________________________________________
Int_t AliRsnLoopDaughter::DoLoop
(AliRsnEvent *evMain, AliRsnDaughterSelector *selMain, AliRsnEvent *, AliRsnDaughterSelector *)
{
//
// Loop function.
// Computes what is needed from passed events.
// Returns the number of pairs successfully processed.
//

   if (!OkEvent(evMain)) return 0;

   Int_t i, il, nadd = 0, nlist = 0;
   TEntryList *list[2] = {0, 0};
   
   if (fDef->IsChargeDefined()) {
      list[0] = selMain->GetSelected(fListID, fDef->GetChargeC());
      list[1] = 0x0;
      nlist = 1;
   } else {
      list[0] = selMain->GetSelected(fListID, '+');
      if (list[0]) {
         list[1] = selMain->GetSelected(fListID, '-');
         nlist = 2;
      } else {
         list[0] = selMain->GetSelected(fListID, '0');
         list[1] = 0x0;
         nlist = 1;
      }
   }
   
   // if it is required to loop over True MC, do this here and skip the rest of the method
   if (fTrueMC) return LoopTrueMC(evMain);
   
   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   
   for (il = 0; il < nlist; il++) {
      if (!list[il]) {
         AliError(Form("List #%d is null", il));
         continue;
      }
      for (i = 0; i < list[il]->GetN(); i++) {
         evMain->SetDaughter(fDaughter, (Int_t)list[il]->GetEntry(i));
         // check matching
         if (fOnlyTrue && !fDef->MatchesPID(&fDaughter)) continue;
         if (!fDef->MatchesCharge(&fDaughter)) continue;
         if (!fDef->MatchesRefType(&fDaughter)) continue;
         fDaughter.FillP(fDef->GetMass());
         // fill outputs
         nadd++;
         next.Reset();
         while ( (out = (AliRsnListOutput*)next()) ) {
            out->Fill(&fDaughter);
         }
      }
   }
   
   return nadd;
}

//_____________________________________________________________________________
Int_t AliRsnLoopDaughter::LoopTrueMC(AliRsnEvent *rsn)
{
//
// Loop on event and fill containers
//

   // check presence of MC reference
   if (!rsn->GetRefMC()) {
      AliError("Need a MC to compute efficiency");
      return 0;
   }
   
   // check event type:
   // must be ESD or AOD, and then use a bool to know in the rest
   if (!rsn->IsESD() && !rsn->IsAOD()) {
      AliError("Need to process ESD or AOD input");
      return 0;
   }
   
   // retrieve the MC primary vertex position
   // and do some additional coherence checks
   Int_t npart = 0;
   TClonesArray *listAOD = 0x0;
   if (rsn->IsESD()) {
      npart = rsn->GetRefMCESD()->GetNumberOfTracks();
   } else {
      AliAODEvent *aod = rsn->GetRefMCAOD();
      listAOD = (TClonesArray*)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
      if (listAOD) npart = listAOD->GetEntries();
   }
   
   // check number of particles
   if (!npart) {
      AliInfo("Empty event");
      return 0;
   }
   
   // utility variables
   Int_t ipart, count = 0;
   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   Int_t pdg = AliRsnDaughter::SpeciesPDG(fDef->GetPID());
   
   
   // loop over particles
   for (ipart = 0; ipart < npart; ipart++) {
      // check i-th particle
      if (rsn->IsESD()) {
         if (!rsn->GetRefMCESD()->Stack()->IsPhysicalPrimary(ipart)) continue;
         AliMCParticle *part = (AliMCParticle*)rsn->GetRefMCESD()->GetTrack(ipart);
         if (TMath::Abs(part->Particle()->GetPdgCode()) != pdg) continue;
         fDaughter.SetRef  (rsn->GetRefMCESD()->GetTrack(ipart));
         fDaughter.SetRefMC(rsn->GetRefMCESD()->GetTrack(ipart));
      } else {
         AliAODMCParticle *part = (AliAODMCParticle*)listAOD->At(ipart);
         if (!part->IsPhysicalPrimary()) continue;
         if (TMath::Abs(part->GetPdgCode()) != pdg) continue;
         fDaughter.SetRef  ((AliAODMCParticle*)listAOD->At(ipart));
         fDaughter.SetRefMC((AliAODMCParticle*)listAOD->At(ipart));
      }
      //if (fDaughter.GetPDG() != AliRsnDaughter::SpeciesPDG(fDef->GetPID())) continue;
      fDaughter.FillP(fDef->GetMass());
      // fill outputs
      count++;
      next.Reset();
      while ( (out = (AliRsnListOutput*)next()) ) {
         out->Fill(&fDaughter);
      }
   }
   
   return count;
}
