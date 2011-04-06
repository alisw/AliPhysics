#include <TEntryList.h>

#include "AliLog.h"
#include "AliPID.h"

#include "AliRsnCutSet.h"

#include "AliRsnPIDManager.h"

ClassImp(AliRsnPIDManager)

//_____________________________________________________________________________
AliRsnPIDManager::AliRsnPIDManager(const char *name, const char *title) : TNamed(name, title)
{
//
// Default constructor.
//
   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
AliRsnPIDManager::~AliRsnPIDManager()
{
   //
   // Destructor
   //
   AliDebug(AliLog::kDebug + 10, "<-");
   for (Int_t iCharge = 0; iCharge < 2; iCharge++) fIdParticles[iCharge].Delete();
   AliDebug(AliLog::kDebug + 10, "->");
}

void AliRsnPIDManager::Print(Option_t* option) const
{
   TNamed::Print(option);

   TEntryList *el;
   AliRsnCutSet *cut;
   for (Int_t iCharge = 0; iCharge < 2; iCharge++) {
      for (Int_t iParticle = 0; iParticle < AliPID::kSPECIES; iParticle++) {
         cut = (AliRsnCutSet*) fCuts.At(iParticle);
         if (!cut) continue;
         el = (TEntryList*)fIdParticles[iCharge].At(iParticle);
         if (el) AliInfo(Form("%d %d %lld %s", iCharge, iParticle, el->GetN(), cut->GetName()));
      }
   }
}

void AliRsnPIDManager::Init()
{
   for (Int_t iCharge = 0; iCharge < 2; iCharge++) {
      for (Int_t iParticle = 0; iParticle < AliPID::kSPECIES; iParticle++) {
         fIdParticles[iCharge].Add(new TEntryList());
//          fCuts.Add(new AliRsnCutSet("", AliRsnTarget::kDaughter));
      }
   }
}


void AliRsnPIDManager::Reset()
{
   TEntryList *el;
   for (Int_t iCharge = 0; iCharge < 2; iCharge++) {
      for (Int_t iParticle = 0; iParticle < AliPID::kSPECIES; iParticle++) {
         el = (TEntryList*)fIdParticles[iCharge].At(iParticle);
         el->Reset();
      }
   }
}


TEntryList* AliRsnPIDManager::GetParticles(const Int_t& charge, const Int_t& pidId)
{
   return (TEntryList*)fIdParticles[charge].At(pidId);
}

void AliRsnPIDManager::SetCut(AliRsnCutSet*cut, const Int_t& pidId)
{
   if (!cut) return;
   if (!fCuts.GetEntriesFast()) Init();
   fCuts.AddAt(cut, pidId);
}

void AliRsnPIDManager::ApplyCuts(AliRsnEvent* ev)
{

   // Loop over event and add entruies to entry list
   Int_t numTracks = ev->GetRefESD()->GetNumberOfTracks();
   Int_t iTrack;
   AliRsnDaughter daughter;
   for (iTrack = 0; iTrack < numTracks; iTrack++) {
      ev->SetDaughter(daughter, iTrack, AliRsnDaughter::kTrack);
      daughter.SetRsnID(iTrack);
      CheckTrack(&daughter, iTrack);
   }

}

void AliRsnPIDManager::CheckTrack(AliRsnDaughter* d, const Int_t & id)
{
   TEntryList *el;
   AliRsnCutSet *cut;
   TString cutName;
   for (Int_t iParticle = 0; iParticle < AliPID::kSPECIES; iParticle++) {

      cut = (AliRsnCutSet*) fCuts.At(iParticle);
      if (!cut) continue;
      cutName = cut->GetName();
      if (!cutName.IsNull()) {
         if (!cut->IsSelected(d)) continue;
         if (d->IsPos()) {
            el = (TEntryList*)fIdParticles[0].At(iParticle);
            el->Enter(id);
         } else if (d->IsNeg()) {
            el = (TEntryList*)fIdParticles[1].At(iParticle);
            el->Enter(id);
         }
      }
   }
}
