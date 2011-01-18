//
// Class AliRsnAnalysisPhiKK
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>
#include <TList.h>

#include "AliRsnFunction.h"
#include "AliRsnAnalysisPhiKK.h"

ClassImp(AliRsnAnalysisPhiKK)

//_____________________________________________________________________________
AliRsnAnalysisPhiKK::AliRsnAnalysisPhiKK(const char *name, Bool_t useKine) :
  AliRsnVAnalysisTaskSE(name, useKine),
  
  fPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455),
  
  fCutEvent      (Form("%s_cutEvent"   , name), AliRsnTarget::kEvent),
  fCutTrackCommon(Form("%s_cutTrackCom", name), AliRsnTarget::kDaughter),
  fCutTrackPos   (Form("%s_cutTrackPos", name), AliRsnTarget::kDaughter),
  fCutTrackNeg   (Form("%s_cutTrackNeg", name), AliRsnTarget::kDaughter),
  fCutPair       (Form("%s_cutPair", name), AliRsnTarget::kMother),
  
  fFuncPM  ("AliRsnFunction", 0),
  fFuncPP  ("AliRsnFunction", 0),
  fFuncMM  ("AliRsnFunction", 0),
  fFuncTrue("AliRsnFunction", 0),
  
  fOutList(0x0)
{
//
// Default constructor.
// Defines another output slot for histograms/ntuples
//

  DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliRsnAnalysisPhiKK::AliRsnAnalysisPhiKK(const AliRsnAnalysisPhiKK& copy) :
  AliRsnVAnalysisTaskSE(copy),
  
  fPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455),
  
  fCutEvent      (copy.fCutEvent),
  fCutTrackCommon(copy.fCutTrackCommon),
  fCutTrackPos   (copy.fCutTrackPos),
  fCutTrackNeg   (copy.fCutTrackNeg),
  fCutPair       (copy.fCutPair),
  
  fFuncPM  (copy.fFuncPM),
  fFuncPP  (copy.fFuncPP),
  fFuncMM  (copy.fFuncMM),
  fFuncTrue(copy.fFuncTrue),
  
  fOutList(0x0)
{
//
// Copy constructor.
//
}

//_____________________________________________________________________________
AliRsnAnalysisPhiKK& AliRsnAnalysisPhiKK::operator=(const AliRsnAnalysisPhiKK& copy)
{
//
// Assigment operator.
//

  AliRsnVAnalysisTaskSE::operator=(copy);
  
  fFuncPM   = copy.fFuncPM;
  fFuncPP   = copy.fFuncPP;
  fFuncMM   = copy.fFuncMM;
  fFuncTrue = copy.fFuncTrue;
  
  fCutEvent       = copy.fCutEvent;
  fCutTrackCommon = copy.fCutTrackCommon;
  fCutTrackPos    = copy.fCutTrackPos;
  fCutTrackNeg    = copy.fCutTrackNeg;
  fCutPair        = copy.fCutPair;
  
  if (fOutList) fOutList->Clear();
  
  return (*this);
}

//_____________________________________________________________________________
void AliRsnAnalysisPhiKK::AddFunction(AliRsnFunction* const fcn)
{
//
// Adds a new computing function to each collection,
// in order to have exactly the sames for each kind of pair.
//
  
  Int_t size = fFuncPM.GetEntries();
  
  new (fFuncPM  [size]) AliRsnFunction(*fcn);
  new (fFuncPP  [size]) AliRsnFunction(*fcn);
  new (fFuncMM  [size]) AliRsnFunction(*fcn);
  new (fFuncTrue[size]) AliRsnFunction(*fcn);
}

//_____________________________________________________________________________
void AliRsnAnalysisPhiKK::RsnUserCreateOutputObjects()
{
//
// Creation of output objects.
// These are created through the utility methods in the analysis manager,
// which asks all the AliRsnPair objects to initialize their output which
// is then linked to the TList data member of this, which will contain all the output.
//

  if (!fOutList) fOutList = new TList;
  fOutList->Clear();
  
  Int_t   i, j, nFunc = fFuncPM.GetEntries();
  TString hName(""), suf[4] = {"PM", "PP", "MM", "True"};
  AliRsnFunction *fcn[4] = {0, 0, 0, 0};
  
  for (i = 0; i < nFunc; i++)
  {
    fcn[0] = (AliRsnFunction*)fFuncPM.At(i);
    fcn[1] = (AliRsnFunction*)fFuncPP.At(i);
    fcn[2] = (AliRsnFunction*)fFuncMM.At(i);
    fcn[3] = (AliRsnFunction*)fFuncTrue.At(i);
    for (j = 0; j < 4; j++)
    {
      hName  = GetName();
      hName += '_';
      hName += suf[j];
      hName += '_';
      hName += fcn[j]->GetName();
      if (fcn[j]->IsUsingTH1()) fOutList->Add(fcn[j]->CreateHistogram(hName.Data(), ""));
      else fOutList->Add(fcn[j]->CreateHistogramSparse(hName.Data(), ""));
    }
  }

  PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisPhiKK::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects,
// using 'reconstructed' or 'MonteCarlo' functions depending on MC-only flag.
//

  // allocate statically all class objects used here
  static TArrayI                   good(0);
  static AliRsnDaughter            kaon[2], temp;
  static AliRsnMother              phi;
  static AliRsnDaughter::ERefType  type;
  static AliRsnFunction           *fcn = 0x0;
  static TClonesArray             *ref = 0x0;
  
  // define constants used for kinematics
  static const Double_t kaonMass = 0.493677;
  
  // simpler variables are declared non static
  Int_t  i, j, k, index, ngood = 0;
  Int_t  tot = AliRsnTarget::GetCurrentEvent()->GetAbsoluteSum();
  Bool_t assignOK, truePair;

  // point to first event in the target 
  AliRsnTarget::SwitchToFirst();
  if (!AliRsnTarget::GetCurrentEvent()) return;
  
  // initially, set the array of good indexes 
  // to the full number of tracks and reset the counter
  good.Set(tot);
  ngood = 0;
  
  // loop on tracks and get those which satisfy cuts
  for (i = 0; i < tot; i++)
  {
    // assign track and skip all that are not charged tracks
    assignOK = AliRsnTarget::GetCurrentEvent()->ConvertAbsoluteIndex(i, index, type);
    if (!assignOK) continue;
    if (type != AliRsnDaughter::kTrack) continue;
    AliRsnTarget::GetCurrentEvent()->SetDaughter(temp, index, AliRsnDaughter::kTrack);
    
    // skip tracks which don't pass common cuts
    if (!fCutTrackCommon.IsSelected(&temp)) continue;
    
    // accept tracks which pass also charge-related cuts
    if ( (temp.Charge() > 0) && (fCutTrackPos.IsSelected(&temp)) )
    {
      good[ngood] = index;
      ++ngood;
    }
    else if ( (temp.Charge() < 0) && (fCutTrackNeg.IsSelected(&temp)) )
    {
      good[ngood] = index;
      ++ngood;
    }
  }
  
  // rese the arrays to the real counts
  good.Set(ngood);
  
  // now that the 'tot' value is useless, set it to
  // the total number of functions, which by construction 
  // is THE SAME for all collections
  tot = fFuncPM.GetEntries();
  
  // fill histograms: do a unique loop on all good indexes
  // and choose the histogram to fill from track charges
  for (i = 0; i < ngood; i++)
  {
    AliRsnTarget::GetCurrentEvent()->SetDaughter(kaon[0], good[i], AliRsnDaughter::kTrack);
    
    for (j = 0; j < ngood; j++)
    {
      // reject equal indexes
      if (good[i] == good[j]) continue;
      AliRsnTarget::GetCurrentEvent()->SetDaughter(kaon[1], good[j], AliRsnDaughter::kTrack);
  
      // adjust charges of pair def
      fPairDef.SetDaughters(AliPID::kKaon, kaon[0].ChargeChar(), AliPID::kKaon, kaon[1].ChargeChar());
    
      // fill the pair using the kaon masses and the passed daughters
      phi.SetDaughters(&kaon[0], kaonMass, &kaon[1], kaonMass);
      
      // check pair cuts
      if (!fCutPair.IsSelected(&phi)) continue;
      
      // choose the functions to fill according to charges
      if (fPairDef.IsLikeSign())
      {
        if (kaon[0].IsPos()) ref = &fFuncPP; else ref = &fFuncMM;
        truePair = kFALSE;
      }
      else
      {
        ref = &fFuncPM;
        truePair = IsTruePair(&kaon[0], &kaon[1]);
      }
        
      // loop on functions in chosen collection and fill
      for (k = 0; k < tot; k++)
      {
        // fill standard histogram
        fcn = (AliRsnFunction*)fFuncPP[k];
        fcn->SetPairDef(&fPairDef);
        fcn->SetPair(&phi);
        fcn->Fill();
        
        // in case of true pair, fill its histogram
        if (truePair)
        {
          fcn = (AliRsnFunction*)fFuncTrue[k];
          fcn->Fill();
        }
      }
    } // end internal loop
  } // end external loop
  
  PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisPhiKK::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhiKK::EventProcess()
{
//
// Customized event pre-processing.
// First checks if the current event passes all cuts,
// and if it does, updates the informations and then
// call the operations which are already defined in the
// omonyme function in mother class
//

  // initially, an event is expected to be bad
  fTaskInfo.SetEventUsed(kFALSE);
  
  // check #1: number of tracks in event (reject empty events)
  Int_t ntracks = fRsnEvent.GetMultiplicity();
  if (ntracks < 1) 
  {
    // empty events are rejected by default
    fTaskInfo.SetEventUsed(kFALSE);
    AliDebug(AliLog::kDebug, "Empty event. Skipping...");
    return kFALSE;
  }

  // check the event cuts and update the info data accordingly
  // events not passing the cuts must be rejected
  if (!fCutEvent.IsSelected(&fRsnEvent))
  {
    fTaskInfo.SetEventUsed(kFALSE);
    return kFALSE;
  }
  
  // if we reach this point, cuts were passed;
  // then additional operations can be done
  
  // find leading particle (without any PID/momentum restriction)
  fRsnEvent.SelectLeadingParticle(0);
  
  // final return value is positive
  // but call the mother class method which updates info object
  fTaskInfo.SetEventUsed(kTRUE);
  return AliRsnVAnalysisTaskSE::EventProcess();
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhiKK::IsTruePair(AliRsnDaughter *d1, AliRsnDaughter *d2)
{
//
// Checks if the two daughters in argument come from the same phi resonance
// and, if they do, check also that they are both kaons
//

  // constants related to PDG
  static const Int_t phiPDG  = 333;
  static const Int_t kaonPDG = 321;

  // check #1: is MC present?
  if (!d1->GetRefMC() || !d2->GetRefMC()) return kFALSE;

  // check #2: same mother?
  Int_t m1 = -1;
  Int_t m2 = -2;
  if (d1->IsESD() && d2->IsESD() )
  {
    if (d1->GetRefMCESD() && d2->GetRefMCESD())
    {
      m1 = d1->GetRefMCESD()->Particle()->GetFirstMother();
      m2 = d2->GetRefMCESD()->Particle()->GetFirstMother();
    }
  }
  if (d1->IsAOD() && d2->IsAOD())
  {
    if (d1->GetRefMCAOD() && d2->GetRefMCAOD())
    {
      m1 = d1->GetRefMCAOD()->GetMother();
      m2 = d2->GetRefMCAOD()->GetMother();
    }
  }
  if (m1 < 0 || m2 < 0 || (m1 > 0 && m2 > 0 && m1 != m2)) return kFALSE;
  
  // check #3: is the common mother a phi (PDG = 333)?
  if (d1->GetMotherPDG() != phiPDG) return kFALSE;
  
  // check #4: are the two particles a K+K- pair?
  m1 = d1->GetPDG();
  m2 = d2->GetPDG();
  if (m1 == kaonPDG && m2 == -kaonPDG) 
    return kTRUE;
  else if (m1 == -kaonPDG && m2 == kaonPDG)
    return kTRUE;
  else
    return kFALSE;
}
