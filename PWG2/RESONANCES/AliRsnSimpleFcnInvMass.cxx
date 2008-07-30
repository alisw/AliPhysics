//
// Class AliRsnSimpleFcnInvMass
//
// This is the most fundamental AliRsnSimpleFunction,
// which computes the invariant mass spectrum of a resonance,
// by correlating pairs of tracks from an event (or mixing, for BG)
//

#include <TH1.h>
#include <TArrayI.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnSimpleFcnInvMass.h"

ClassImp(AliRsnSimpleFcnInvMass)

//________________________________________________________________________________________
AliRsnSimpleFcnInvMass::AliRsnSimpleFcnInvMass() :
  AliRsnSimpleFunction(),
  fUseMCValues(kFALSE),
  fPIDMethod(AliRsnDaughter::kNoPID)
{
//
// Constructor.
// Only default initializations.
//
}

//________________________________________________________________________________________
AliRsnSimpleFcnInvMass::AliRsnSimpleFcnInvMass
(const char *name, AliRsnDaughter::EPIDMethod method, AliRsnPairDef *pd,
 AliRsnHistoDef *hd, AliRsnCutMgr *cuts, Option_t *option) :
  AliRsnSimpleFunction(name, pd, hd, cuts, option),
  fUseMCValues(kFALSE),
  fPIDMethod(method)
{
//
// Constructor.
// Only default initializations.
//
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFcnInvMass::ProcessOne(AliRsnEvent* event)
{
//
// Process a single event to build the invariant mass histogram
// from the correlation of track pairs, according to the AliRsnPairDef settings.
// The collection of used tracks depends on the PID method chosen.
//

    if (!event) {
        AliError("Argument cannot be NULL.");
        return kFALSE;
    }
    
    // check PID method
    if (fPIDMethod < AliRsnDaughter::kNoPID || fPIDMethod >= AliRsnDaughter::kMethods) {
        AliError("PID method not properly initialized");
        return kFALSE;
    }
    
    // check event cut
    if (!CutPass(event)) return kFALSE;

    // assign pointers to the list of indexes to be used
    TArrayI *listCharged1 = 0x0, *listCharged2 = 0x0;
    if (fPIDMethod == AliRsnDaughter::kNoPID || fPIDMethod == AliRsnDaughter::kWeighted) {
        listCharged1 = event->GetCharged(fPairDef->GetCharge(0));
        listCharged2 = event->GetCharged(fPairDef->GetCharge(1));
    }
    else {
        listCharged1 = event->GetTracksArray(fPIDMethod, fPairDef->GetCharge(0), fPairDef->GetType(0));
        listCharged2 = event->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), fPairDef->GetType(1));
    }
    if (!listCharged1 || !listCharged2) return 0.0;
    TArrayI &list1 = *listCharged1;
    TArrayI &list2 = *listCharged2;

    Stat_t count = 0;
    Int_t i1, i2, start2;
    AliRsnDaughter *track1 = 0x0, *track2 = 0x0;

    // loop on particle of type 1
    for (i1 = 0; i1 < list1.GetSize(); i1++) {
        track1 = event->GetTrack(list1[i1]);
        if (!CutPass(track1)) continue;
        // loop on particle of type 2
        // in case we are building a like-sign histogram with particles
        // of the same type, we must avoid that each pair is used twice
        start2 = 0;
        if (listCharged1 == listCharged2) start2 = i1 + 1;
        for (i2 = start2; i2 < list2.GetSize(); i2++) {
            track2 = event->GetTrack(list2[i2]);
            if (!CutPass(track2)) continue;
            if (Add(track1, track2)) count++;
        }
    }

    return (count > (Stat_t)0);
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFcnInvMass::ProcessTwo
(AliRsnEvent* event1, AliRsnEvent* event2)
{
//
// Process a single event to build the invariant mass histogram
// from the correlation of track pairs, according to the AliRsnPairDef settings.
// The collection of used tracks depends on the PID method chosen.
//

    if (!event1 || !event2) {
        AliError("Arguments cannot be NULL.");
        return kFALSE;
    }
    
    // check event cut
    if (!CutPass(event1)) return kFALSE;
    if (!CutPass(event2)) return kFALSE;

    // assign pointers to the list of indexes to be used
    TArrayI *listCharged1 = 0x0, *listCharged2 = 0x0;
    if (fPIDMethod == AliRsnDaughter::kNoPID || fPIDMethod == AliRsnDaughter::kWeighted) {
        listCharged1 = event1->GetCharged(fPairDef->GetCharge(0));
        listCharged2 = event2->GetCharged(fPairDef->GetCharge(1));
    }
    else {
        listCharged1 = event1->GetTracksArray(fPIDMethod, fPairDef->GetCharge(0), fPairDef->GetType(0));
        listCharged2 = event2->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), fPairDef->GetType(1));
    }
    if (!listCharged1 || !listCharged2) return 0.0;
    TArrayI &list1 = *listCharged1;
    TArrayI &list2 = *listCharged2;

    Stat_t count = 0;
    Int_t i1, i2;
    AliRsnDaughter *track1 = 0x0, *track2 = 0x0;

    // loop on particle of type 1
    for (i1 = 0; i1 < list1.GetSize(); i1++) {
        track1 = event1->GetTrack(list1[i1]);
        if (!CutPass(track1)) continue;
        // loop on particle of type 2
        for (i2 = 0; i2 < list2.GetSize(); i2++) {
            track2 = event2->GetTrack(list2[i2]);
            if (!CutPass(track2)) continue;
            if (Add(track1, track2)) count++;
        }
    }

    return (count > (Stat_t)0);
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFcnInvMass::Add
(AliRsnDaughter *track1, AliRsnDaughter *track2)
{
//
// Add a histogram entry from two tracks, if they pass the cuts.
// The order matters, because track1 is processed using first element
// in the AliRsnPairDef definition, and track2 is processed using second ones.
// In case the study is done only for true pairs, this is checked automatically.
//

    // set the static variable pointing to PID method
    // to the current one in use by this object
    AliRsnDaughter::SetPIDMethod(fPIDMethod);

    // setup pair and check pair cuts
    fPair.SetPair(track1, track2);
    if (!CutPass(&fPair)) return kFALSE;

    // computation variables
    Double_t mass1[2] = {fPairDef->GetMass(0), fPairDef->GetMass(1)};
    Double_t mass2[2] = {fPairDef->GetMass(1), fPairDef->GetMass(0)};
    
    // define the value of the entry weight according to PID method selected:
    // - in case of "weighted" computations, which is done with all tracks
    //   using the PID computed probabilities, the weight will be the product
    //   of the ones related to the used PID for the tracks
    Bool_t   useWeight = kFALSE;
    Double_t entryWeight[2];
    entryWeight[0] = fPairDef->ComputeWeight(track1, track2);
    entryWeight[1] = fPairDef->ComputeWeight(track2, track1);
    if (fPIDMethod != AliRsnDaughter::kWeighted) {
        useWeight = kTRUE;
        entryWeight[0] = 1.0;
        entryWeight[0] = 1.0;
    }
    
    // when filling the histogram, if the defined pair is a like-sign one,
    // and if the two PID involved are different, it is necessary
    // to fill it twice, exchanging the tracks with respect to pair definition
    // in all other cases (unlike sign, like sign with same PID) this is not done
    Int_t i, lastAdd = 0;
    if (fPairDef->IsLikeSign() && !fPairDef->HasEqualTypes()) lastAdd = 1;
    
    // make computation
    Double_t invmass;
    for (i = 0; i <= lastAdd; i++) {
        if (fUseMCValues) invmass = fPair.GetInvMassMC(mass1[i], mass2[i]);
        else invmass = fPair.GetInvMass(mass1[i], mass2[i]);
        
        if (useWeight) fHisto1D->Fill(invmass, entryWeight[i]);
        else fHisto1D->Fill(invmass);
    }

    return kTRUE;
}
