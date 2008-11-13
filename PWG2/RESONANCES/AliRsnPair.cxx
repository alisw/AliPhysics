//
// *** Class AliRsnPair ***
//
// "Core" method for defining the work on a pari of particles.
// For one analysis, one must setup one of this for each pair he wants to analyze,
// adding to it all analysis which he desires to do.
// Here he defines the cuts, and the particle types and charges, and can add
// functions which do different operations on the same pair, and some binning
// with respect to some kinematic variables (eta, momentum)
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <Riostream.h>
#include "TObjArray.h"

#include "AliLog.h"

#include "AliRsnFunction.h"
#include "AliRsnPairParticle.h"

#include "AliRsnPair.h"

ClassImp(AliRsnPair)

//_____________________________________________________________________________
AliRsnPair::AliRsnPair
(EPairType type, AliRsnPairDef *def, Int_t mixNum) :
  TObject(),
  fIsMixed(kFALSE),
  fUseMC(kFALSE),
  fIsLikeSign(kFALSE),
  fMixNum(mixNum),
  fMixingCut(0x0),
  fMatched(mixNum),
  fPairDef(def),
  fPairType(type),
  fTypePID(AliRsnDaughter::kRealistic),
  fCutMgr(0),
  fFunctions("AliRsnFunction", 0)
{
//
// Default constructor
//

  SetUp(type);
  if (!fIsMixed) fMatched.Set(1);
}
//_____________________________________________________________________________
AliRsnPair::~AliRsnPair()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnPair::SetUp(EPairType type)
{
//
// Sets up flag values by the pair types
//

  switch (type)
  {
    case kNoPID:
      SetAllFlags(AliRsnDaughter::kNoPID, kFALSE, kFALSE);
      break;
    case kNoPIDMix:
      SetAllFlags(AliRsnDaughter::kNoPID, kTRUE, kFALSE);
      break;
    case kRealisticPID:
      SetAllFlags(AliRsnDaughter::kRealistic, kFALSE, kFALSE);
      break;
    case kRealisticPIDMix:
      SetAllFlags(AliRsnDaughter::kRealistic, kTRUE, kFALSE);
      break;
    case kPerfectPID:
      // SetAllFlags (AliRsnDaughter::kPerfect, kFALSE, kFALSE);
      SetAllFlags(AliRsnDaughter::kPerfect, kFALSE, kTRUE);
      break;
    case kPerfectPIDMix:
      // SetAllFlags (AliRsnDaughter::kPerfect, kTRUE, kFALSE);
      SetAllFlags(AliRsnDaughter::kPerfect, kTRUE, kTRUE);
      break;
    default :
      AliWarning("Wrong type selected: setting up for kRealisticPID.");
      SetAllFlags(AliRsnDaughter::kRealistic, kFALSE, kFALSE);
      break;
  }
}

//_____________________________________________________________________________
void AliRsnPair::SetAllFlags(AliRsnDaughter::EPIDMethod pidType, Bool_t isMix, Bool_t useMC)
{
//
// Sets up all flags values
//

  fTypePID = pidType;
  fIsMixed = isMix;
  fUseMC = useMC;
}

//_____________________________________________________________________________
void AliRsnPair::Init()
{
//
// Init pair
//

  fIsLikeSign = fPairDef->IsLikeSign();
  Print();
}

//_____________________________________________________________________________
void AliRsnPair::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//
  AliInfo(Form("%s", GetPairHistTitle(0x0).Data()));
  AliInfo(Form("PDG %d %d", AliRsnPID::PDGCode(fPairDef->GetType(0)),
               AliRsnPID::PDGCode(fPairDef->GetType(1))));
  AliInfo(Form("Masses %f %f", fPairDef->GetMass(0), fPairDef->GetMass(1)));
  AliInfo(Form("Number of functions %d", fFunctions.GetEntries()));
  switch(fTypePID) {
    case AliRsnDaughter::kNoPID:
      AliInfo("PID method: none");
      break;
    case AliRsnDaughter::kRealistic:
      AliInfo("PID method: realistic");
      break;
    case AliRsnDaughter::kPerfect:
      AliInfo("PID method: perfect");
      break;
    default:
      AliInfo("PID method: undefined");
    }
}

//_____________________________________________________________________________
void AliRsnPair::ProcessPair(AliRsnEventBuffer *buf)
{
//
// Process current event in the passed buffer.
// If this pair is a single-event pair, only that event is processed,
// otherwise, if it is a mixing pair, all good matches are found and mixed.
//

  if (fIsMixed && buf->GetEventsBufferIndex() >= fMixNum) ProcessPairMix(buf);
  else ProcessPairSingle(buf);

  /*
  // track type/charge #0 in pairDef is taken from this event,
  // track tipe/charge #1 is taken from the matched event
  AliRsnEvent *e1 = buf->GetCurrentEvent();
  if (!e1) return;
  TArrayI* array1 = e1->GetTracksArray(fTypePID, fPairDef->GetCharge(0), fPairDef->GetType(0));

  // reset array of matched events
  Int_t i;
  for (i = 0; i < fMixNum; i++) fMatched[i] = -1;

  // here searches for all matched events:
  // in case of single-event pairs, there will be only one, equal to current,
  // otherwise there will be many
  Int_t lastOkEvent = buf->IndexOf(e1) - 1;
  Double_t nMatched = 0.0;
  if (fIsMixed) {
    for (i = 0; i < fMixNum; i++) {
      // find other event by event cut
      AliRsnEvent *e2 = FindEventByEventCut(buf, lastOkEvent);
      if (!e2) break;
      fMatched[i] = lastOkEvent;
      lastOkEvent--;
      nMatched++;
    }
    if (fMatched.GetSize() < 0) {
      AliWarning(Form("Event #%d: found no events to match", buf->IndexOf(e1)));
      return;
    }
  }
  else {
    fMatched[0] = buf->IndexOf(e1);
    nMatched = 1;
  }

  // in order to balance the problem that different events could be matched
  // a different number of times (between 1 and the maximum allowd quantity),
  // for each event, the histograms are filled with a variable weight depending
  // on the number of matches found:
  Double_t weight = (Double_t)fMixNum / (Double_t)nMatched;

  // now that matched events are found, they are used to fill the histograms
  TArrayI* array2 = 0;
  for (i = 0; i < fMixNum; i++) {
    if (fMatched[i] < 0) break;
    AliRsnEvent *e2 = buf->GetEvent(fMatched[i]);
    array2 = e2->GetTracksArray(fTypePID, fPairDef->GetCharge(1), fPairDef->GetType(1));
    if (fIsMixed) LoopPair(e1, array1, e2, array2, weight);
    else LoopPair(e1, array1, e2, array2);
  }
  */
}

//_____________________________________________________________________________
void AliRsnPair::ProcessPairSingle(AliRsnEventBuffer *buf)
{
//
// SINGLE EVENT PROCESSING
// This function fills the functions' histograms using tracks
// from the same event for both components of the defined pair.
// This function is used for signal, like-sign and rotated background.
//

  AliRsnEvent *event = buf->GetCurrentEvent();
  if (!event) return;

  TArrayI* array1 = event->GetTracksArray(fTypePID, fPairDef->GetCharge(0), fPairDef->GetType(0));
  TArrayI* array2 = event->GetTracksArray(fTypePID, fPairDef->GetCharge(1), fPairDef->GetType(1));

  LoopPair(event, array1, event, array2);
}

//_____________________________________________________________________________
void AliRsnPair::ProcessPairMix(AliRsnEventBuffer *buf)
{
//
// MIXED EVENT PROCESSING
// This function fills the functions' histograms using
// tracks from different events.
// Used for event mixing.
//

  // track type/charge #0 in pairDef is taken from this event,
  // track tipe/charge #1 is taken from the matched events
  AliRsnEvent *e1 = buf->GetCurrentEvent();
  if (!e1) return;
  TArrayI* array10 = e1->GetTracksArray(fTypePID, fPairDef->GetCharge(0), fPairDef->GetType(0));
  TArrayI* array11 = e1->GetTracksArray(fTypePID, fPairDef->GetCharge(1), fPairDef->GetType(1));

  // find matched events
  Int_t i, iev = buf->GetEventsBufferIndex();
  Int_t nMatched = FindMatchedEvents(iev, buf);
  if (!nMatched) {
    AliWarning(Form("Event #%d: found no events to match", iev));
    return;
  }
  else if (nMatched < fMixNum) {
    AliWarning(Form("Event #%d: found only %d events to match", iev, nMatched));
    return;
  }

  /*
  else {
    TString str("Matched events: ");
    for (i = 0; i < nMatched; i++) str.Append(Form("%d ", fMatched[i]));
    AliInfo(Form("Event #%d: %s", iev, str.Data()));
  }
  */

  // in order to balance the problem that different events could be matched
  // a different number of times (between 1 and the maximum allowd quantity),
  // for each event, the histograms are filled with a variable weight depending
  // on the number of matches found:
  // Double_t weight = 0.5 * (Double_t)fMixNum / (Double_t)nMatched;

  // now that matched events are found, they are used to fill the histograms
  TArrayI *array20 = 0, *array21 = 0;
  for (i = 0; i < fMixNum; i++) {
    if (fMatched[i] < 0) break;
    AliRsnEvent *e2 = buf->GetEvent(fMatched[i]);
    array20 = e2->GetTracksArray(fTypePID, fPairDef->GetCharge(0), fPairDef->GetType(0));
    array21 = e2->GetTracksArray(fTypePID, fPairDef->GetCharge(1), fPairDef->GetType(1));
    // track type/charge #0 from event1, track type/charge #1 from event2
    LoopPair(e1, array10, e2, array21);
    // track type/charge #1 from event1, track type/charge #0 from event2
    LoopPair(e1, array11, e2, array20);
  }
}

//_____________________________________________________________________________
Int_t AliRsnPair::FindMatchedEvents(Int_t evIndex, AliRsnEventBuffer *buf)
{
//
// Resets the 'fMatched' data member and stores into it
// a number of well-matched events found in the buffer,
// (maximum amount = fMixNum).
// The argument is the index of the event to be matched,
// in the event buffer.
//

  // reset array of matched events
  Int_t i;
  for (i = 0; i < fMixNum; i++) fMatched[i] = -1;

  // starts from the position behind the one
  // of the event to be matched and goes backward
  // if it reaches the value 0, stops
  AliRsnEvent *eventToBeMatched = buf->GetCurrentEvent();
  AliRsnEvent *matchEvent = 0x0;
  Int_t checkIndex = evIndex - 1, curIndex = 0;
  for (;;checkIndex--) {
    if (checkIndex < 0) checkIndex = buf->GetEventsBufferSize() - 1;
    if (checkIndex == evIndex) break;
    matchEvent = buf->GetEvent(checkIndex);
    if (!matchEvent) continue;
    if (fMixingCut) {
      if (fMixingCut->IsSelected(AliRsnCut::kMixEvent, eventToBeMatched, matchEvent)) continue;
    }
    // assign to current array slot the matched event
    // and increment current slot and stops if it exceeds array size
    fMatched[curIndex++] = checkIndex;
    if (curIndex >= fMixNum) break;
  }

  // returns the current index value,
  // which is also the number of matched events found
  return curIndex;
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnPair::FindEventByEventCut(AliRsnEventBuffer *buf, Int_t& num)
{
//
// For now it just returns num events before current event
// in buffer (buf)
// TODO event cut selection
//

  AliRsnEvent *returnEvent = 0x0;

  if (fIsMixed)
  {
    //returnEvent = buf->GetEvent(buf->GetEventsBufferIndex() - num);
    returnEvent = buf->GetNextGoodEvent(num, fMixingCut);
  }
  else
  {
    returnEvent = buf->GetCurrentEvent();
  }

  return returnEvent;
}

//_____________________________________________________________________________
void AliRsnPair::LoopPair
(AliRsnEvent * ev1, TArrayI * a1, AliRsnEvent * ev2, TArrayI * a2, Double_t weight)
{
//
// Loop on all pairs of tracks of the defined types/charges,
// using the arrays of indexes and the events containing them.
//

  if (!a1) {AliDebug(4, "No TArrayI 1 from currentEvent->GetTracksArray(...)"); return;}
  if (!a2) {AliDebug(4, "No TArrayI 2 from currentEvent->GetTracksArray(...)"); return;}

  AliRsnDaughter::SetPIDMethod(fTypePID);
  AliRsnDaughter *daughter1 = 0;
  AliRsnDaughter *daughter2 = 0;
  AliRsnFunction *fcn = 0;
  Int_t j, startj = 0;
    for (Int_t i = 0; i < a1->GetSize(); i++)
    {
        // get track #1
        daughter1 = (AliRsnDaughter *) ev1->GetTrack(a1->At(i));
        if (!daughter1) continue;
        // cuts on track #1
        if (!CutPass(daughter1)) continue;
        // get track #2
        daughter2 = 0;
        // check starting index for searching the event:
        // for like-sign pairs we avoid duplicating the pairs
        if (fIsLikeSign) startj = i+1; else startj = 0;
        // AliInfo(Form("%d",startj));
        // loop on event for all track #2 to be combined with the found track #1
        for (j = startj; j < a2->GetSize(); j++)
        {
            daughter2 = (AliRsnDaughter *) ev2->GetTrack(a2->At(j));
            if (!daughter2) continue;
            // cuts on track #2
            if (!CutPass(daughter2)) continue;
            // make pair
            AliRsnPairParticle pairParticle;
            pairParticle.SetPair(daughter1, daughter2);
            // cuts on pair
            if (!CutPass(&pairParticle)) continue;
            // fill all histograms
            TObjArrayIter nextFcn(&fFunctions);
            while ( (fcn = (AliRsnFunction*)nextFcn()) ) {
                fcn->Fill(&pairParticle, fPairDef, weight);
            }
        }
    }
}

//_____________________________________________________________________________
TList * AliRsnPair::GenerateHistograms(TString prefix)
{
//
// Generates needed histograms
//

  TList *list = new TList();
  list->SetName(GetPairHistName(0x0).Data());

  Char_t hName[255], hTitle[255];
  AliRsnFunction *fcn = 0;
  for (Int_t i=0;i< fFunctions.GetEntries();i++)
  {
    fcn = (AliRsnFunction*)fFunctions.At(i);
    sprintf(hName, "%s_%s", prefix.Data(), GetPairHistName(fcn).Data());
    sprintf(hTitle, "%s", GetPairHistTitle(fcn).Data());
    TList *histos = fcn->Init(hName, hTitle);
    histos->Print();
    list->Add(histos);
  }

  return list;
}

//_____________________________________________________________________________
void AliRsnPair::GenerateHistograms(TString prefix, TList *tgt)
{
//
// Generates needed histograms
//

  if (!tgt) {
    AliError("NULL target list!");
    return;
  }

  Char_t hName[255], hTitle[255];
  AliRsnFunction *fcn = 0;
  for (Int_t i=0;i< fFunctions.GetEntries();i++)
  {
    fcn = (AliRsnFunction*)fFunctions.At(i);
    sprintf(hName, "%s_%s", prefix.Data(), GetPairHistName(fcn).Data());
    sprintf(hTitle, "%s", GetPairHistTitle(fcn).Data());
    fcn->Init(hName, hTitle, tgt);
  }
}

//_____________________________________________________________________________
TString AliRsnPair::GetPairTypeName(EPairType type) const
{
//
// Returns type name, made with particle names ant chosen PID
//
  switch (type)
  {
    case kNoPID : return ("NOPID_");break;
    case kNoPIDMix : return ("NOPIDMIX_");break;
    case kRealisticPID : return ("REALISTIC_");break;
    case kRealisticPIDMix : return ("REALISTICMIX_");break;
    case kPerfectPID : return ("PERFECT_");break;
    case kPerfectPIDMix : return ("PERFECTMIX_");break;
    case kTruePairs : return ("TRUEPAIRS_"); break;
    default:
      AliWarning("Unrecognized value of EPairTypeName argument");
      break;
  }

  return "NOTYPE";
}

//_____________________________________________________________________________
TString AliRsnPair::GetPairName() const
{
//
// Retruns pair name
//
  TString sName;
  sName += GetPairTypeName(fPairType);
  sName += fPairDef->GetPairName();

  return sName;
}

//_____________________________________________________________________________
TString AliRsnPair::GetPairHistName(AliRsnFunction *fcn, TString text) const
{
//
// Returns eff. mass histogram name
//

  TString sName;
  if (fcn)
  {
    sName = fcn->GetFcnName();
    sName += "_";
  }
  sName += GetPairName();
  sName += "_";
  if (fCutMgr) sName += fCutMgr->GetName();
  sName += text;

  return sName;
}

//_____________________________________________________________________________
TString AliRsnPair::GetPairHistTitle(AliRsnFunction *fcn, TString text) const
{
//
// Returns eff. mass histogram title
//

  TString sTitle;
  if (fcn)
  {
    sTitle = fcn->GetFcnTitle();
    sTitle += " ";
  }
  sTitle += GetPairName();
  sTitle +=" ";
  if (fCutMgr) sTitle += fCutMgr->GetTitle();
  sTitle += text;

  return sTitle;
}

//_____________________________________________________________________________
void AliRsnPair::AddFunction(AliRsnFunction *fcn)
{
//
// Adds a new computing function
//
  Int_t size = fFunctions.GetEntries();
  new(fFunctions[size]) AliRsnFunction(*fcn);
}

//________________________________________________________________________________________
Bool_t AliRsnPair::CutPass(AliRsnDaughter *d)
{
//
// Check if the AliRsnDaughter argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
//
  if (!fCutMgr) return kTRUE;
  else return fCutMgr->IsSelected(AliRsnCut::kParticle, d);
}

//________________________________________________________________________________________
Bool_t AliRsnPair::CutPass(AliRsnPairParticle *p)
{
//
// Check if the AliRsnPairParticle argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
// In this case, another separate check which could be necessary
// concerns the possibility that the two tracks are a "true pair" of
// daughters of the same resonance. If the corresponding flag is set,
// this further check is done, and the method returns kTRUE only
// when also this check is passed.
//

  if (!fCutMgr) return kTRUE;
  else return fCutMgr->IsSelected(AliRsnCut::kPair, p);
}

//________________________________________________________________________________________
Bool_t AliRsnPair::CutPass(AliRsnEvent *e)
{
//
// Check if the AliRsnEvent argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
//
  if (!fCutMgr) return kTRUE;
  else return fCutMgr->IsSelected(AliRsnCut::kEvent, e);
}
