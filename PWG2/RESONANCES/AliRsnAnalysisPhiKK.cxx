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
  fGood(0),
  fMother(),
  fPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455),
  fCutEvent(Form("%s_cutEvent", name), AliRsnTarget::kEvent),
  fCutTrackCommon(Form("%s_cutTrackCom", name), AliRsnTarget::kDaughter),
  fCutTrackPos(Form("%s_cutTrackPos", name), AliRsnTarget::kDaughter),
  fCutTrackNeg(Form("%s_cutTrackNeg", name), AliRsnTarget::kDaughter),
  fCutPair(Form("%s_cutPair", name), AliRsnTarget::kMother),
  fFuncPM("AliRsnFunction", 0),
  fFuncPP("AliRsnFunction", 0),
  fFuncMM("AliRsnFunction", 0),
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
  fGood(0),
  fMother(),
  fPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455),
  fCutEvent(Form("%s_cutEvent", copy.GetName()), AliRsnTarget::kEvent),
  fCutTrackCommon(Form("%s_cutTrackCom", copy.GetName()), AliRsnTarget::kDaughter),
  fCutTrackPos(Form("%s_cutTrackPos", copy.GetName()), AliRsnTarget::kDaughter),
  fCutTrackNeg(Form("%s_cutTrackNeg", copy.GetName()), AliRsnTarget::kDaughter),
  fCutPair(Form("%s_cutPair", copy.GetName()), AliRsnTarget::kMother),
  fFuncPM(copy.fFuncPM),
  fFuncPP(copy.fFuncPP),
  fFuncMM(copy.fFuncMM),
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
      fcn[j]->SetPairDef(&fPairDef);
      fcn[j]->SetPair(&fMother);
      
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

  // skip if the global event pointers are NULL
  // and point to first event in the target 
  AliRsnEvent *event = AliRsnEvent::GetCurrentEvent1();
  if (!event) return;
  AliRsnTarget::SwitchToFirst();
  
  // select good kaons, applying all cuts
  Int_t i, index, tot = event->GetAbsoluteSum(), ngood = 0;
  Bool_t assignOK;
  AliRsnDaughter::ERefType type;
  fGood.Set(tot);
  for (i = 0; i < tot; i++)
  {
    // assign track and skip all that are not charged tracks
    assignOK = event->ConvertAbsoluteIndex(i, index, type);
    if (!assignOK || type != AliRsnDaughter::kTrack) continue;
    event->SetDaughter(fDaughter[0], index, AliRsnDaughter::kTrack);
    
    // skip tracks which don't pass common cuts
    if (!fCutTrackCommon.IsSelected(fDaughter)) continue;
    
    // accept tracks which pass also charge-related cuts
    if ((fDaughter[0].Charge() > 0) && (fCutTrackPos.IsSelected(&fDaughter[0])))
    {
      fGood[ngood] = index;
      ++ngood;
    }
    else if ((fDaughter[0].Charge() < 0) && (fCutTrackNeg.IsSelected(&fDaughter[0])))
    {
      fGood[ngood] = index;
      ++ngood;
    }
  }
  fGood.Set(ngood);
  
  // now that the 'tot' value is useless, set it to
  // the total number of functions, which by construction is THE SAME
  // for all collections
  tot = fFuncPM.GetEntries();
  
  // loop on selected tracks and fill histograms
  Int_t i0, i1, m0, m1, motherPDG;
  const Double_t  kaonMass = 0.493677;
  const Int_t     phiPDG   = 333;
  const Int_t     kaonPDG  = 313;
  AliRsnFunction *fcn;
  for (i0 = 0; i0 < ngood; i0++)
  {
    index = fGood[i0];
    event->SetDaughter(fDaughter[0], index, AliRsnDaughter::kTrack);
    
    for (i1 = i0 + 1; i1 < ngood; i1++)
    {
      index = fGood[i1];
      event->SetDaughter(fDaughter[1], index, AliRsnDaughter::kTrack);
      
      // skip in case the two indexes match
      if (fDaughter[0].GetID() == fDaughter[1].GetID()) continue;
      
      // adjust charges of pair def
      fPairDef.SetDaughters(AliPID::kKaon, fDaughter[0].ChargeChar(), AliPID::kKaon, fDaughter[1].ChargeChar());
    
      // fill the pair using the kaon masses and the passed daughters
      fMother.SetDaughters(&fDaughter[0], kaonMass, &fDaughter[1], kaonMass);
      
      // check pair cuts
      if (!fCutPair.IsSelected(&fMother)) continue;
      
      // if pair is like-sign, fill appropriate histos
      if (fPairDef.IsLikeSign())
      {
        if (fDaughter[0].Charge() > 0)
        {
          for (i = 0; i < tot; i++)
          {
            fcn = (AliRsnFunction*)fFuncPP[i];
            fcn->Fill();
          }
        }
        else
        {
          for (i = 0; i < tot; i++)
          {
            fcn = (AliRsnFunction*)fFuncMM[i];
            fcn->Fill();
          }
        }
      }
      else
      {
        // if pair is unlike-sign, check that it is true
        motherPDG = fMother.CommonMother(m0, m1);
        if (motherPDG == phiPDG)
        {
          if (m0 < 0 || m1 < 0) motherPDG = 0;
          if (TMath::Abs(fDaughter[0].GetPDG()) != kaonPDG) motherPDG = 0;
          if (TMath::Abs(fDaughter[1].GetPDG()) != kaonPDG) motherPDG = 0;
        }
        
        // fill unlike-sign histo
        for (i = 0; i < tot; i++)
        {
          fcn = (AliRsnFunction*)fFuncPM[i];
          fcn->Fill();
          if (motherPDG == phiPDG)
          {
            fcn = (AliRsnFunction*)fFuncTrue[i];
            fcn->Fill();
          }
        }
      }
    } // for (i1)
  } // for (i0)
  
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
