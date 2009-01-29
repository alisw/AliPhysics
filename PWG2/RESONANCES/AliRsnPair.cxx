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
#include <TObjArray.h>

#include "AliLog.h"

#include "AliRsnFunction.h"
#include "AliRsnPairParticle.h"

#include "AliRsnPair.h"

ClassImp(AliRsnPair)

//_____________________________________________________________________________
AliRsnPair::AliRsnPair(EPairType type, AliRsnPairDef *def) :
  TObject(),
  fIsMixed(kFALSE),
  fPairType(type),
  fPIDMethod(AliRsnDaughter::kRealistic),
  fPairDef(def),
  fCutMgr(0),
  fFunctions("AliRsnFunction", 0)
{
//
// Default constructor
//

  SetUp(type);
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
      SetAllFlags(AliRsnDaughter::kNoPID, kFALSE);
      break;
    case kNoPIDMix:
      SetAllFlags(AliRsnDaughter::kNoPID, kTRUE);
      break;
    case kRealisticPID:
      SetAllFlags(AliRsnDaughter::kRealistic, kFALSE);
      break;
    case kRealisticPIDMix:
      SetAllFlags(AliRsnDaughter::kRealistic, kTRUE);
      break;
    case kPerfectPID:
      SetAllFlags (AliRsnDaughter::kPerfect, kFALSE);
      break;
    case kPerfectPIDMix:
      SetAllFlags (AliRsnDaughter::kPerfect, kTRUE);
      break;
    default :
      AliWarning("Wrong type selected: setting up for realistic PID - no mixing.");
      SetAllFlags(AliRsnDaughter::kRealistic, kFALSE);
      break;
  }
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

  switch(fPIDMethod)
  {
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
void AliRsnPair::ProcessPair(AliRsnEvent *ev1, AliRsnEvent *ev2)
{
//
// Fills the functions' histograms using tracks from passed events.
// What tracks are taken in each event depend from the order of
// track types defined in the AliRsnPairDef for this object:
//  - tracks of type 0 are the ones stored as pairDef data members with index [0]
//    ---> taken from first argument (ev1)
//  - tracks of type 1 are the ones stored as pairDef data members with index [1]
//    ---> taken from second argument (ev2)
//
// When doing single-event analysis (e.g. signal, like-sign), second argument
// can be NULL, and it will be forced to point to the same object of first one.
//

  if (!ev2) ev2 = ev1;

  TArrayI *array1 = ev1->GetTracksArray(fPIDMethod, fPairDef->GetCharge(0), fPairDef->GetType(0));
  TArrayI *array2 = ev2->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), fPairDef->GetType(1));

  LoopPair(ev1, array1, ev2, array2);
}

//_____________________________________________________________________________
void AliRsnPair::LoopPair
(AliRsnEvent * ev1, TArrayI * a1, AliRsnEvent * ev2, TArrayI * a2)
{
//
// Loop on all pairs of tracks of the defined types/charges,
// using the arrays of indexes and the events containing them.
// This method is private, for safety reasons.
//

  if (!a1) {AliDebug(4, "No TArrayI 1 from currentEvent->GetTracksArray(...)"); return;}
  if (!a2) {AliDebug(4, "No TArrayI 2 from currentEvent->GetTracksArray(...)"); return;}

  // cuts on events
  if (!CutPass(ev1) || !CutPass(ev2)) return;

  AliRsnDaughter::SetPIDMethod(fPIDMethod);
  AliRsnDaughter *daughter1 = 0;
  AliRsnDaughter *daughter2 = 0;
  AliRsnFunction *fcn = 0;

  Bool_t isLikeSign = fPairDef->IsLikeSign();
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
        if (isLikeSign) startj = i+1; else startj = 0;
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
                fcn->Fill(&pairParticle, fPairDef);
            }
        }
    }
}

//_____________________________________________________________________________
TList * AliRsnPair::GenerateHistograms(TString prefix)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the output TList.
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
    //histos->Print();
    list->Add(histos);
  }

  return list;
}

//_____________________________________________________________________________
void AliRsnPair::GenerateHistograms(TString prefix, TList *tgt)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the TList passed as argument
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
// Returns definitive histogram name
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
// Returns definitive histogram title
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
  new (fFunctions[size]) AliRsnFunction(*fcn);
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
