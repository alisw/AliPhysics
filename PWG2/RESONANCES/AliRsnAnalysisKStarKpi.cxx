//
// Class AliRsnAnalysisKStarKpi
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
#include "AliRsnAnalysisKStarKpi.h"

ClassImp(AliRsnAnalysisKStarKpi)

//_____________________________________________________________________________
AliRsnAnalysisKStarKpi::AliRsnAnalysisKStarKpi(const char *name, Bool_t useKine) :
  AliRsnVAnalysisTaskSE(name, useKine),
  fGoodK(0),
  fGoodPi(0),
  fKaon(),
  fPion(),
  fMother(),
  fPairDef(AliPID::kKaon, '+', AliPID::kPion, '-', 313, 0.896),
  fCutEvent(Form("%s_cutEvent", name), AliRsnTarget::kEvent),
  fCutTrackCommon(Form("%s_cutTrackCom", name), AliRsnTarget::kDaughter),
  fCutTrackKaon(Form("%s_cutTrackKaon", name), AliRsnTarget::kDaughter),
  fCutTrackPion(Form("%s_cutTrackPion", name), AliRsnTarget::kDaughter),
  fCutPair(Form("%s_cutPair", name), AliRsnTarget::kMother),
  fFuncPM("AliRsnFunction", 0),
  fFuncMP("AliRsnFunction", 0),
  fFuncPP("AliRsnFunction", 0),
  fFuncMM("AliRsnFunction", 0),
  fFuncTruePM("AliRsnFunction", 0),
  fFuncTrueMP("AliRsnFunction", 0),
  fOutList(0x0)
{
//
// Default constructor.
// Defines another output slot for histograms/ntuples
//

  DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliRsnAnalysisKStarKpi::AliRsnAnalysisKStarKpi(const AliRsnAnalysisKStarKpi& copy) :
  AliRsnVAnalysisTaskSE(copy),
  fGoodK(0),
  fGoodPi(0),
  fKaon(),
  fPion(),
  fMother(),
  fPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455),
  fCutEvent(copy.fCutEvent),
  fCutTrackCommon(copy.fCutTrackCommon),
  fCutTrackKaon(copy.fCutTrackKaon),
  fCutTrackPion(copy.fCutTrackPion),
  fCutPair(copy.fCutPair),
  fFuncPM(copy.fFuncPM),
  fFuncMP(copy.fFuncMP),
  fFuncPP(copy.fFuncPP),
  fFuncMM(copy.fFuncMM),
  fFuncTruePM(copy.fFuncTruePM),
  fFuncTrueMP(copy.fFuncTrueMP),
  fOutList(0x0)
{
//
// Copy constructor.
//
}

//_____________________________________________________________________________
AliRsnAnalysisKStarKpi& AliRsnAnalysisKStarKpi::operator=(const AliRsnAnalysisKStarKpi& copy)
{
//
// Assigment operator.
//

  AliRsnVAnalysisTaskSE::operator=(copy);
  
  fFuncPM     = copy.fFuncPM;
  fFuncMP     = copy.fFuncMP;
  fFuncPP     = copy.fFuncPP;
  fFuncMM     = copy.fFuncMM;
  fFuncTruePM = copy.fFuncTruePM;
  fFuncTrueMP = copy.fFuncTrueMP;
  
  fCutEvent       = copy.fCutEvent;
  fCutTrackCommon = copy.fCutTrackCommon;
  fCutTrackPion   = copy.fCutTrackPion;
  fCutTrackKaon   = copy.fCutTrackKaon;
  fCutPair        = copy.fCutPair;
  
  if (fOutList) fOutList->Clear();
  
  return (*this);
}

//_____________________________________________________________________________
void AliRsnAnalysisKStarKpi::AddFunction(AliRsnFunction* const fcn)
{
//
// Adds a new computing function to each collection,
// in order to have exactly the sames for each kind of pair.
//
  
  Int_t size = fFuncPM.GetEntries();
  
  new (fFuncPM    [size]) AliRsnFunction(*fcn);
  new (fFuncMP    [size]) AliRsnFunction(*fcn);
  new (fFuncPP    [size]) AliRsnFunction(*fcn);
  new (fFuncMM    [size]) AliRsnFunction(*fcn);
  new (fFuncTruePM[size]) AliRsnFunction(*fcn);
  new (fFuncTrueMP[size]) AliRsnFunction(*fcn);
}

//_____________________________________________________________________________
void AliRsnAnalysisKStarKpi::RsnUserCreateOutputObjects()
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
  TString hName(""), suf[6] = {"PM", "MP", "PP", "MM", "TruePM", "TrueMP"};
  AliRsnFunction *fcn[6] = {0, 0, 0, 0};
  
  for (i = 0; i < nFunc; i++)
  {
    fcn[0] = (AliRsnFunction*)fFuncPM.At(i);
    fcn[1] = (AliRsnFunction*)fFuncMP.At(i);
    fcn[2] = (AliRsnFunction*)fFuncPP.At(i);
    fcn[3] = (AliRsnFunction*)fFuncMM.At(i);
    fcn[4] = (AliRsnFunction*)fFuncTruePM.At(i);
    fcn[5] = (AliRsnFunction*)fFuncTrueMP.At(i);
    for (j = 0; j < 6; j++)
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
void AliRsnAnalysisKStarKpi::RsnUserExec(Option_t*)
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
  
  // select good tracks, applying all cuts
  // in this point we use only the 'fPion' data member
  // but it holds for both kinds of particles
  Int_t i, index, tot = event->GetAbsoluteSum(), ngoodK = 0, ngoodPi = 0;
  Bool_t assignOK;
  AliRsnDaughter::ERefType type;
  fGoodK.Set(tot);
  fGoodPi.Set(tot);
  for (i = 0; i < tot; i++)
  {
    // assign track and skip all that are not charged tracks
    assignOK = event->ConvertAbsoluteIndex(i, index, type);
    if (!assignOK || type != AliRsnDaughter::kTrack) continue;
    event->SetDaughter(fPion, index, AliRsnDaughter::kTrack);
    
    // skip tracks which don't pass common cuts
    if (!fCutTrackCommon.IsSelected(&fPion)) continue;
    
    // accept tracks which pass also charge-related cuts
    if (fCutTrackPion.IsSelected(&fPion))
    {
      fGoodPi[ngoodPi] = index;
      ++ngoodPi;
    }
    else if (fCutTrackKaon.IsSelected(&fPion))
    {
      fGoodK[ngoodK] = index;
      ++ngoodK;
    }
  }
  fGoodK.Set(ngoodK);
  fGoodPi.Set(ngoodPi);
  
  // now that the 'tot' value is useless, set it to
  // the total number of functions, which by construction is THE SAME
  // for all collections
  tot = fFuncPM.GetEntries();
  
  // loop on selected tracks and fill histograms
  // external loop is on kaons (first element in pairDef)
  // internal loop is on pions (second element in pairDef)
  Int_t i0, i1, m0, m1, motherPDG;
  const Double_t  kaonMass = 0.493677;
  const Double_t  pionMass = 0.13957;
  const Int_t     kstarPDG = 313;
  const Int_t     kaonPDG  = 313;
  const Int_t     pionPDG  = 211;
  AliRsnFunction *fcn;
  for (i0 = 0; i0 < ngoodK; i0++)
  {
    index = fGoodK[i0];
    event->SetDaughter(fKaon, index, AliRsnDaughter::kTrack);
    
    for (i1 = 0; i1 < ngoodPi; i1++)
    {      
      index = fGoodPi[i1];
      event->SetDaughter(fPion, index, AliRsnDaughter::kTrack);
      
      // skip in case the two indexes match
      if (fPion.GetID() == fKaon.GetID()) continue;
      
      // adjust charges of pair def
      fPairDef.SetDaughters(AliPID::kKaon, fKaon.ChargeChar(), AliPID::kPion, fPion.ChargeChar());
    
      // fill the pair using the kaon masses and the passed daughters
      fMother.SetDaughters(&fKaon, kaonMass, &fPion, pionMass);
      
      // check pair cuts
      if (!fCutPair.IsSelected(&fMother)) continue;
      
      // if pair is like-sign, fill appropriate histos
      if (fPairDef.IsLikeSign())
      {
        if (fKaon.Charge() > 0)
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
        if (motherPDG == kstarPDG)
        {
          if (m0 < 0 || m1 < 0) motherPDG = 0;
          if (TMath::Abs(fKaon.GetPDG()) != kaonPDG) motherPDG = 0;
          if (TMath::Abs(fPion.GetPDG()) != pionPDG) motherPDG = 0;
        }
        
        // fill unlike-sign histo (the appropriate one)
        if (fKaon.Charge() > 0)
        {
          for (i = 0; i < tot; i++)
          {
            fcn = (AliRsnFunction*)fFuncPM[i];
            fcn->Fill();
            if (motherPDG == kstarPDG)
            {
              fcn = (AliRsnFunction*)fFuncTruePM[i];
              fcn->Fill();
            }
          }
        }
        else
        {
          for (i = 0; i < tot; i++)
          {
            fcn = (AliRsnFunction*)fFuncMP[i];
            fcn->Fill();
            if (motherPDG == kstarPDG)
            {
              fcn = (AliRsnFunction*)fFuncTrueMP[i];
              fcn->Fill();
            }
          }
        }
      }
    } // for (i1)
  } // for (i0)
  
  PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisKStarKpi::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisKStarKpi::EventProcess()
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
