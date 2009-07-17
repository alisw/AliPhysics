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

#include <TList.h>

#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnFunction.h"
#include "AliRsnPIDIndex.h"
#include "AliRsnCutMgr.h"
#include "AliRsnCutStd.h"

#include "AliRsnPair.h"

ClassImp(AliRsnPair)

//_____________________________________________________________________________
AliRsnPair::AliRsnPair(EPairType type, AliRsnPairDef *def) :
    TObject(),
    fOnlyTrue(kFALSE),
    fIsMixed(kFALSE),
    fPairType(type),
    fPIDMethod(AliRsnDaughter::kRealistic),
    fPairDef(def),
    fCutMgr(0),
    fFunctions("AliRsnFunction", 0),
    fTrack1(),
    fTrack2(),
    fPairParticle()
{
//
// Default constructor
//
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
  SetUp(type);
}
//_____________________________________________________________________________
AliRsnPair::~AliRsnPair()
{
//
// Destructor
//
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPair::SetUp(EPairType type)
{
//
// Sets up flag values by the pair types
//
  AliDebug(AliLog::kDebug+2,"<-");
  switch (type) {
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
    SetAllFlags(AliRsnDaughter::kPerfect, kFALSE);
    break;
  case kPerfectPIDMix:
    SetAllFlags(AliRsnDaughter::kPerfect, kTRUE);
    break;
  default :
    AliWarning("Wrong type selected: setting up for realistic PID - no mixing.");
    SetAllFlags(AliRsnDaughter::kRealistic, kFALSE);
    break;
  }
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPair::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//
  AliDebug(AliLog::kDebug+2,"<-");
  AliInfo(Form("%s", GetPairHistTitle(0x0).Data()));
  AliInfo(Form("PDG %d %d", AliPID::ParticleCode(fPairDef->GetType(0)),
               AliPID::ParticleCode(fPairDef->GetType(1))));
  AliInfo(Form("Masses %f %f", fPairDef->GetMass(0), fPairDef->GetMass(1)));
  AliInfo(Form("Number of functions %d", fFunctions.GetEntries()));

  switch (fPIDMethod) {
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
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPair::LoopPair(AliRsnPIDIndex *const pidIndex1, AliRsnEvent *const ev1, AliRsnPIDIndex *const pidIndex2, AliRsnEvent *const ev2)
{
//
// Prepare the loop for computation of functions.
// Each PIDIndex is used to retrieve the appropriate array of indexes
// of the tracks to be used in each event.
// In case of single-event analysis, only the first two arguments are used
// and both arrays are taken from the same PIDIndex and will loop on the same event
// In case of mixing, all arguments are used, and first set of tracks will be found
// in the first event with the first PIDIndex, and the second set of tracks will
// be found in second event with second PIDIndex.
//

  AliDebug(AliLog::kDebug+2,"<-");

  TArrayI *a1 = 0;
  TArrayI *a2 = 0;

  if (fPIDMethod == AliRsnDaughter::kNoPID) {
    AliDebug(AliLog::kDebug+2, Form("Returning indexes of with NO PID (%d) ...", fPIDMethod));
    a1 = pidIndex1->GetTracksArray(fPIDMethod, fPairDef->GetCharge(0), AliPID::kUnknown);
    if (pidIndex2 && ev2)
      a2 = pidIndex2->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), AliPID::kUnknown);
    else
      a2 = pidIndex1->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), AliPID::kUnknown);
  } else {
    AliDebug(AliLog::kDebug+2, Form("Returning indexes of with PID (%d) ...", fPIDMethod));
    a1 = pidIndex1->GetTracksArray(fPIDMethod,fPairDef->GetCharge(0), (AliPID::EParticleType)fPairDef->GetType(0));
    if (pidIndex2 && ev2)
      a2 = pidIndex2->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), (AliPID::EParticleType)fPairDef->GetType(1));
    else
      a2 = pidIndex1->GetTracksArray(fPIDMethod, fPairDef->GetCharge(1), (AliPID::EParticleType)fPairDef->GetType(1));
  }

  LoopPair(a1, a2, ev1, ev2);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPair::LoopPair(TArrayI *a1, TArrayI *a2, AliRsnEvent *ev1, AliRsnEvent *ev2)
{
//
// Loop on all pairs of tracks of the defined types/charges,
// using the arrays of indexes and the events containing them.
// This method is private, for safety reasons.
//
  AliDebug(AliLog::kDebug+2,"<-");
  if (!ev1) {AliError(Form("ev1 is %p. skipping LoopPair() ...",ev1))return;}
  AliDebug(AliLog::kDebug+1,"ev1 is OK ...");


  if (!ev2) {
    if (fIsMixed) {
      AliDebug(AliLog::kDebug+1, "ev2 is null and fIsMixed is true. Skipping ...");
      return;
    }
    ev2 = ev1;
  }

  if (!a1) {AliDebug(AliLog::kDebug+1, "No TArrayI 1 from currentEvent->GetTracksArray(...)"); return;}
  if (!a2) {AliDebug(AliLog::kDebug+1, "No TArrayI 2 from currentEvent->GetTracksArray(...)"); return;}

  AliDebug(AliLog::kDebug+1,Form("a1=%d a2=%d",a1->GetSize(),a2->GetSize()));
  if (a1->GetSize()<=0) {AliDebug(AliLog::kDebug+1, "Size of TArrayI 1 is 0 or less ..."); return;}
  if (a2->GetSize()<=0) {AliDebug(AliLog::kDebug+1, "Size of TArrayI 2 is 0 or less ..."); return;}

  AliDebug(AliLog::kDebug,Form("Indexes_a1=%d Indexes_a2=%d",a1->GetSize(),a2->GetSize()));

  // cuts on events
  if (!CutPass(ev1) || !CutPass(ev2)) return;
  AliDebug(AliLog::kDebug+1,"Event cut passed...");
  AliRsnDaughter::SetPIDMethod(fPIDMethod);
  AliRsnFunction *fcn = 0;

  Bool_t isLikeSign = fPairDef->IsLikeSign();
  Int_t j, startj = 0;

  for (Int_t i = 0; i < a1->GetSize(); i++) {
    // get track #1
    ev1->SetDaughter(fTrack1, a1->At(i));
    if (!fTrack1.IsOK()) continue;
    // assign the required PID type to track #1
    fTrack1.SetRequiredPID(fPairDef->GetType(0));
    AliDebug(AliLog::kDebug+1,"daughter1 is OK ...");
    // cuts on track #1
    if (!CutPass(&fTrack1)) continue;
    AliDebug(AliLog::kDebug+1,"daughter1 cut passed ...");
    // check starting index for searching the event:
    // for like-sign pairs we avoid duplicating the pairs
    if (isLikeSign) startj = i+1; else startj = 0;
    // loop on event for all track #2 to be combined with the found track #1
    for (j = startj; j < a2->GetSize(); j++) {
      ev2->SetDaughter(fTrack2, a2->At(j));
      if (!fTrack2.IsOK()) continue;
      AliDebug(AliLog::kDebug+1,"daughter2 is OK ...");
      // assign the required PID type to track #2
      fTrack2.SetRequiredPID(fPairDef->GetType(1));
      // cuts on track #2
      if (!CutPass(&fTrack2)) continue;
      AliDebug(AliLog::kDebug+1,"daughter2 cut passed ...");
      // make pair
      fPairParticle.SetPair(&fTrack1, &fTrack2);
      // in case of request, check that it is true pair
      if (fOnlyTrue)
      {
        if (fPairParticle.CommonMother() != fPairDef->GetMotherPDG()) continue;
        if (fPairParticle.GetDaughter(0)->PerfectPID() != fPairDef->GetType(0)) continue;
        if (fPairParticle.GetDaughter(1)->PerfectPID() != fPairDef->GetType(1)) continue;
      }
      // cuts on pair
      if (!CutPass(&fPairParticle)) continue;
      AliDebug(AliLog::kDebug+1, "pairParticle cut passed");

      if (AliLog::GetDebugLevel("",GetName()) == AliLog::kDebug)
        fPairParticle.PrintInfo();

      // fill all histograms
      TObjArrayIter nextFcn(&fFunctions);
      while ((fcn = (AliRsnFunction*)nextFcn())) {
        fcn->SetPairDef(fPairDef);
        fcn->SetPair(&fPairParticle);
        fcn->SetEvent(ev1);
        fcn->Fill();
      }
    }
  }
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
TList * AliRsnPair::GenerateHistograms(TString prefix,TList *list)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the output TList.
//
  AliDebug(AliLog::kDebug+2,"<-");
//   if (!list){
//     TList *list = new TList();
//     list->SetName(GetPairHistName(0x0).Data());
//   }
  Char_t hName[255], hTitle[255];
  //AliRsnFunction *fcn = 0;
  AliRsnFunction *fcn = 0;
  for (Int_t i=0;i< fFunctions.GetEntries();i++) {
    fcn = (AliRsnFunction*)fFunctions.At(i);
    sprintf(hName, "%s_%s", prefix.Data(), GetPairHistName(fcn).Data());
    sprintf(hTitle, "%s", GetPairHistTitle(fcn).Data());
    //TList *histos = fcn->Init(hName, hTitle);
    list->Add(fcn->CreateHistogram(hName, hTitle));
    //histos->Print();
    //list->Add(histos);
  }
//   cout << "PRINTING LIST" << endl;
//   list->Print();
  AliDebug(AliLog::kDebug+2,"->");
  return list;
}



//_____________________________________________________________________________
TString AliRsnPair::GetPairTypeName(EPairType type) const
{
//
// Returns type name, made with particle names ant chosen PID
//
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
  switch (type) {
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
  AliDebug(AliLog::kDebug+2,"<-");
  TString sName;
  sName += GetPairTypeName(fPairType);
  sName += fPairDef->GetPairName();
  AliDebug(AliLog::kDebug+2,"->");
  return sName;
}

//_____________________________________________________________________________
TString AliRsnPair::GetPairHistName(AliRsnFunction *const fcn, TString text) const
{
//
// Returns definitive histogram name
//
  AliDebug(AliLog::kDebug+2,"<-");
  TString sName;
  if (fcn) {
    sName = fcn->GetName();
    sName += "_";
  }
  sName += GetPairName();
  sName += "_";
  if (fCutMgr) sName += fCutMgr->GetName();
  sName += text;
  AliDebug(AliLog::kDebug+2,"->");
  return sName;
}

//_____________________________________________________________________________
TString AliRsnPair::GetPairHistTitle(AliRsnFunction *const fcn, TString text) const
{
//
// Returns definitive histogram title
//
  AliDebug(AliLog::kDebug+2,"<-");
  TString sTitle;
  if (fcn) {
    sTitle = fcn->GetTitle();
    sTitle += " ";
  }
  sTitle += GetPairName();
  sTitle +=" ";
  if (fCutMgr) sTitle += fCutMgr->GetTitle();
  sTitle += text;
  AliDebug(AliLog::kDebug+2,"->");
  return sTitle;
}

//_____________________________________________________________________________
void AliRsnPair::AddFunction(AliRsnFunction *const fcn)
{
//
// Adds a new computing function
//
  AliDebug(AliLog::kDebug+2,"<-");
  Int_t size = fFunctions.GetEntries();
  new(fFunctions[size]) AliRsnFunction(*fcn);
  AliDebug(AliLog::kDebug+2,"->");
}

//________________________________________________________________________________________
Bool_t AliRsnPair::CutPass(AliRsnDaughter *d)
{
//
// Check if the AliRsnDaughter argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
//
  AliDebug(AliLog::kDebug+2,"<-AliRsnDaughter");
  AliDebug(AliLog::kDebug+2,"->");
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
  AliDebug(AliLog::kDebug+2,"<-AliRsnPairParticle");
  AliDebug(AliLog::kDebug+2,"->");
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
  AliDebug(AliLog::kDebug+2,"<-AliRsnEvent");
  AliDebug(AliLog::kDebug+2,"->");
  if (!fCutMgr) return kTRUE;
  else return fCutMgr->IsSelected(AliRsnCut::kEvent, e);

}
