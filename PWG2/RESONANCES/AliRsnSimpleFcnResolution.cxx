//
// Class AliRsnSimpleFcnResolution
//
// This is the most fundamental AliRsnSimpleFunction,
// which computes the invariant mass spectrum of a resonance,
// by correlating pairs of tracks from an event (or mixing, for BG)
//

#include <TH1.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPID.h"
#include "AliRsnSimpleFcnResolution.h"

ClassImp(AliRsnSimpleFcnResolution)

//________________________________________________________________________________________
AliRsnSimpleFcnResolution::AliRsnSimpleFcnResolution() :
    AliRsnSimpleFunction()
{
//
// Constructor.
// Only default initializations.
//

  fTrueFlag = kTRUE;
  fMixFlag = kFALSE;
}

//________________________________________________________________________________________
AliRsnSimpleFcnResolution::AliRsnSimpleFcnResolution
(const char *name, AliRsnPairDef *pd,
 AliRsnHistoDef *hd, AliRsnCutMgr *cuts, Option_t *option) :
    AliRsnSimpleFunction(name, pd, hd, cuts, option)
{
//
// Constructor.
// Only default initializations.
//

  fTrueFlag = kTRUE;
  fMixFlag = kFALSE;
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFcnResolution::ProcessOne(AliRsnEvent* event)
{
//
// Process a single event to build the invariant mass resolution histogram
// from the correlation of track pairs, according to the AliRsnPairDef settings.
// This class uses always the perfect PID.
//

  if (!event)
  {
    AliError("Argument cannot be NULL.");
    return kFALSE;
  }

  // check event cut
  if (!CutPass(event)) return kFALSE;

  // assign pointers to the list of indexes to be used
  AliRsnDaughter::EPIDMethod mtd = AliRsnDaughter::kPerfect;
  TArrayI *plist1 = event->GetTracksArray(mtd, fPairDef->GetCharge(0), fPairDef->GetType(0));
  TArrayI *plist2 = event->GetTracksArray(mtd, fPairDef->GetCharge(1), fPairDef->GetType(1));
  if (!plist1 || !plist2) return 0.0;
  TArrayI &list1 = *plist1;
  TArrayI &list2 = *plist2;

  Stat_t count = 0;
  Int_t i1, i2, start2;
  AliRsnDaughter *track1 = 0x0, *track2 = 0x0;

  // loop on particle of type 1
  for (i1 = 0; i1 < list1.GetSize(); i1++)
  {
    track1 = event->GetTrack(list1[i1]);
    if (!CutPass(track1)) continue;
    // loop on particle of type 2
    // in case we are building a like-sign histogram with particles
    // of the same type, we must avoid that each pair is used twice
    start2 = 0;
    if (plist1 == plist2) start2 = i1 + 1;
    for (i2 = start2; i2 < list2.GetSize(); i2++)
    {
      track2 = event->GetTrack(list2[i2]);
      if (!CutPass(track2)) continue;
      if (Add(track1, track2)) count++;
    }
  }

  return (count > (Stat_t)0);
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFcnResolution::Add
(AliRsnDaughter *track1, AliRsnDaughter *track2)
{
//
// Add a histogram entry from two tracks, if they pass the cuts.
// The order matters, because track1 is processed using first element
// in the AliRsnPairDef definition, and track2 is processed using second ones.
// In case the study is done only for true pairs, this is checked automatically.
//

  // setup pair and check pair cuts
  fPair.SetPair(track1, track2);
  if (!CutPass(&fPair)) return kFALSE;

  // computation variables
  Double_t mass1 = fPairDef->GetMass(0);
  Double_t mass2 = fPairDef->GetMass(1);

  // make computation
  Double_t invmass = fPair.GetInvMass(mass1, mass2);
  Double_t invmassMC = fPair.GetInvMassMC(mass1, mass2);
  Double_t res = (invmassMC - invmass) / invmassMC;

  // fill
  fHisto1D->Fill(res);

  return kTRUE;
}
