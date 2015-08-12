#include <iostream>

#include "TParticle.h"

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisTaskHMTFMCMultEst.h"
#include "AliMultiplicityEstimators.h"

using namespace std;

ClassImp(AliAnalysisTaskHMTFMCMultEst)

Bool_t IsPi0PhysicalPrimary(Int_t index, AliStack *stack)
{
  // pi0 are considered unstable and thus not primary in AliStack.
  // If the given track index is a pi0 and if that pi0 has no ancestors, it is a primary
  // Return kFALSE if the index points to anything but a pi0

  TParticle* p = stack->Particle(index);
  Int_t ist = p->GetStatusCode();

  // Initial state particle
  //if (ist > 1) return kFALSE; // pi0 are not initial state (?)
  // pi0 are unstable and thus this returned false in the original implementation.
  //if (!stack->IsStable(pdg)) return kFALSE;

  Int_t pdg = TMath::Abs(p->GetPdgCode());

  // The function is only for pi0's so I'm out of here if its not a pi'0 we are looking at
  if (pdg != 111) return kFALSE;

  if (index < stack->GetNprimary()) {
    // Particle produced by generator
    return kTRUE;
  }
  else {
    // Particle produced during transport
    Int_t imo =  p->GetFirstMother();
    TParticle* pm  = stack->Particle(imo);
    Int_t mpdg = TMath::Abs(pm->GetPdgCode());
    // Check for Sigma0
    if ((mpdg == 3212) &&  (imo <  stack->GetNprimary())) return kTRUE;

    // Check if it comes from a pi0 decay
    if ((mpdg == 111) && (imo < stack->GetNprimary()))   return kTRUE;

    // Check if this is a heavy flavor decay product
    Int_t mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

    // Light hadron
    if (mfl < 4) return kFALSE;

    // Heavy flavor hadron produced by generator
    if (imo <  stack->GetNprimary()) {
      return kTRUE;
    }

    // To be sure that heavy flavor has not been produced in a secondary interaction
    // Loop back to the generated mother
    while (imo >=  stack->GetNprimary()) {
      imo = pm->GetFirstMother();
      pm  =  stack->Particle(imo);
    }
    mpdg = TMath::Abs(pm->GetPdgCode());
    mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

    if (mfl < 4) {
      return kFALSE;
    } else {
      return kTRUE;
    }
    } // produced by generator ?
}


AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst()
  : AliAnalysisTaskSE(), fMyOut(0), fEstimatorsList(0), fEstimatorNames(0),
    festimators(0), fRequireINELgt0(0), fRunconditions(0), fParticleCounter(0)
{
}

//________________________________________________________________________
AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst(const char *name)
  : AliAnalysisTaskSE(name), fMyOut(0), fEstimatorsList(0), fEstimatorNames(0),
    festimators(0), fRequireINELgt0(0), fRunconditions(0), fParticleCounter(0)
{
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

void AliAnalysisTaskHMTFMCMultEst::AddEstimator(const char* n)
{
  if (!fEstimatorNames.IsNull()) fEstimatorNames.Append(",");
  fEstimatorNames.Append(n);
}

void AliAnalysisTaskHMTFMCMultEst::InitEstimators()
{
  fEstimatorsList = new TList;
  fEstimatorsList->SetOwner();
  fEstimatorsList->SetName("estimators");

  TObjArray* arr = fEstimatorNames.Tokenize(",");
  TObject*   obj = 0;
  TIter      next(arr);
  std::cout << "Init estimators... " << std::endl;
  while ((obj = next())) {
    MultiplicityEstimatorBase* e = MakeEstimator(obj->GetName());
    fEstimatorsList->Add(e);
  }
}
//________________________________________________________________________
MultiplicityEstimatorBase*
AliAnalysisTaskHMTFMCMultEst::MakeEstimator(const TString& name)
{
  if (name.BeginsWith("Total")) return new EtaBase("Total", "full #eta coverage ");
  if (name.BeginsWith("EtaLt05")) return new EtaBase("EtaLt05", "| #eta| #leq 0.5", -0.5, 0.0, 0.0, 0.5);
  if (name.BeginsWith("EtaLt08")) return new EtaBase("EtaLt08", "| #eta| #leq 0.8", -0.8, 0.0, 0.0, 0.8);
  if (name.BeginsWith("EtaLt15")) return new EtaBase("EtaLt15", "| #eta| #leq 1.5", -1.5, 0.0, 0.0, 1.5);
  if (name.BeginsWith("Eta08_15")) return new EtaBase("Eta08_15", "0.8 #leq | #eta| #leq 1.5",
						      -1.5, -0.8, 0.8, 1.5);
  if (name.BeginsWith("V0A")) return new EtaBase("V0A", "2.8 #leq #eta #leq 5.1",
						 0.0, 0.0, 2.8, 5.1);
  if (name.BeginsWith("V0C")) return new EtaBase("V0C", "-3.7 #leq #eta #leq -1.7",
						 -3.7, -1.7, 0.0, 0.0);
  if (name.BeginsWith("V0M")) return new EtaBase("V0M", "-3.7 #leq #eta #leq -1.7 || 2.8 #leq #eta #leq 5.1",
						 -3.7, -1.7, 2.8, 5.1);

  return 0;
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::UserCreateOutputObjects()
{
  fMyOut = new TList();
  fMyOut->SetOwner();

  InitEstimators();

  TIter next(fEstimatorsList);
  MultiplicityEstimatorBase* e = 0;
  while ((e = static_cast<MultiplicityEstimatorBase*>(next()))) {
    e->RegisterHistograms(fMyOut);
    // putting estimators into a vector for easier looping in UserExec.
    // it is only available on the slaves
    festimators.push_back(e);
  }

  // Suppress annoying printout
  AliLog::SetGlobalLogLevel(AliLog::kError);

  PostData(1, fMyOut);
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::UserExec(Option_t *)
{
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }
  std::vector<MultiplicityEstimatorBase*>::iterator iter, end;
  for (iter =festimators.begin(), end = festimators.end(); iter != end; ++iter) {//estimator loop
    (*iter)->PreEvent(mcEvent);
  }//estimator loop

  // Track loop for establishing multiplicity and checking for INEL > 0
  Bool_t isINEL_gt_0(kFALSE);
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    // Pass the particle on to the estimators if it is a primary. Extra check for pi0's is needed since they are unstable
    if (mcEvent->Stack()->IsPhysicalPrimary(iTrack) ||
	IsPi0PhysicalPrimary(iTrack, mcEvent->Stack())){
      if (TMath::Abs(track->Eta()) < 1) isINEL_gt_0 = kTRUE;
      for (iter =festimators.begin(), end = festimators.end(); iter != end; ++iter) { //estimator loop
	(*iter)->ProcessTrackForMultiplicityEstimation(track);
      }//estimator loop
    }
  }  //track loop

  // Track loop with known multiplicity in each estimator
  // Skip this if event is not INEL>0 and it is required to be so
  if (!fRequireINELgt0 || isINEL_gt_0){
    for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
      AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
      if (!track) {
	Printf("ERROR: Could not receive track %d", iTrack);
	continue;
      }
      if (mcEvent->Stack()->IsPhysicalPrimary(iTrack) ||
	  IsPi0PhysicalPrimary(iTrack, mcEvent->Stack())){
	for (iter =festimators.begin(), end = festimators.end(); iter != end; ++iter) { //estimator loop
	  (*iter)->ProcessTrackWithKnownMultiplicity(track);
	}//estimator loop
      }
    }  //track loop

    // Increment eventcounters etc.
    for (iter =festimators.begin(), end = festimators.end(); iter != end; ++iter) { //estimator loop
      (*iter)->PostEvent();
    }//estimator loop
  }
  // Post output data.
  PostData(1, fMyOut);
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::Terminate(Option_t *)
{
  // recreates the fEstimatorsList
  InitEstimators();

  // This list is associated to a read only file
  fMyOut = static_cast<TList*> (GetOutputData(1));
  if (!fMyOut) {
    Error("Terminate", "Didn't get sum container");
    return;
  }
  TList* results = new TList;
  results->SetName("terminateResults");

  TIter nextEst(fEstimatorsList);
  MultiplicityEstimatorBase* e = 0;
  std::cout << "Terminating estimators..." << std::endl;
  while ((e = static_cast<MultiplicityEstimatorBase*>(nextEst()))) {
    e->Terminate(fMyOut);
  }
  PostData(1, fMyOut);

  fRunconditions = new TList;
  fRunconditions->SetName("rclist");
  fRunconditions->SetOwner(0);
  TObjString *rcinfo = new TObjString;
  if (fRequireINELgt0) rcinfo->SetString("INELgt0_true");
  else rcinfo->SetString("INELgt0_false");
  fRunconditions->Add(rcinfo);
  PostData(2, fRunconditions);
}
