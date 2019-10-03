#include <iostream>
#include <new>

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisTaskHMTFMCMultEst.h"

// Classifiers:
#include "AliEventClassifierMult.h"
#include "AliEventClassifierMPI.h"
#include "AliEventClassifierQ2.h"
#include "AliEventClassifierSpherocity.h"
#include "AliEventClassifierSphericity.h"

// Observables (aka Plots)
#include "AliObservableEtaNch.h"
#include "AliObservableClassifierpTPID.h"
#include "AliObservableCorrelationsOfClassifiers.h"

using namespace std;

ClassImp(AliAnalysisTaskHMTFMCMultEst)

//________________________________________________________________________
AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst()
: AliAnalysisTaskSE(), fMyOut(0), fClassifiers(0), fObservables(0), fGlobalTrigger(0), fGlobalSystem(0),
  fGlobalTriggerClassifiers(0)
{
}

//________________________________________________________________________
AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst(const char *name)
  : AliAnalysisTaskSE(name), fMyOut(0), fClassifiers(0), fObservables(0), fGlobalTrigger(0), fGlobalSystem(0),
    fGlobalTriggerClassifiers(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::UserCreateOutputObjects()
{
  fMyOut = new TList();
  fMyOut->SetOwner();

  //////////////////////////////////////////
  // Multiplicity based classifiers       //
  //////////////////////////////////////////
  // Some classifiers are used as "reference" classifiers. Ie. We do not compute N! correlation
  // plots between N classifiers. We restrict ourselves to calculating only the correlations of every
  // classifier to a small number of reference classifiers.
  Float_t region2[] = {-0.5, 0.5};  // Use this array to pass two pass two edges to the constructor
  AliEventClassifierMult* refClassifierEtaLt05 =
    new AliEventClassifierMult("EtaLt05", "|#eta| < 0.5", region2, 2, true, true, fMyOut, fGlobalSystem);
  fClassifiers.push_back(refClassifierEtaLt05);

  region2[0] = -0.8; region2[1] = 0.8;
  fClassifiers.push_back(new AliEventClassifierMult("EtaLt08", "|#eta| < 0.8", region2, 2, true, true, fMyOut, fGlobalSystem));

  region2[0] = -1.0; region2[1] = 1.0;
  AliEventClassifierMult* etaLt1 = new AliEventClassifierMult("EtaLt1", "|#eta| < 1.0", region2, 2, true, true, fMyOut, fGlobalSystem);
  fClassifiers.push_back(etaLt1);

  region2[0] = -1.5; region2[1] = 1.5;
  fClassifiers.push_back(new AliEventClassifierMult("EtaLt15", "|#eta| < 1.5", region2, 2, true, true, fMyOut, fGlobalSystem));

  region2[0] = 2.8; region2[1] = 5.1;
  AliEventClassifierMult* v0a = new AliEventClassifierMult("V0A", "2.8 < #eta < 5.1", region2, 2, true, true, fMyOut, fGlobalSystem);
  fClassifiers.push_back(v0a);

  region2[0] = -3.7; region2[1] = -1.7;
  AliEventClassifierMult* v0c = new AliEventClassifierMult("V0C", "-3.7 < #eta < -1.7", region2, 2, true, true, fMyOut, fGlobalSystem);
  fClassifiers.push_back(v0c);

  region2[0] = -8.7; region2[1] = 8.7; // not inclusive region, not charged!
  fClassifiers.push_back(new AliEventClassifierMult("ZDC", "|#eta| > 8.7", region2, 2, false, false, fMyOut, fGlobalSystem));

  Float_t region4[] = {-3.7, -1.7, 2.8, 5.1};
  fClassifiers.push_back(new AliEventClassifierMult("V0M", "V0A + V0C", region4, 4, true, true, fMyOut, fGlobalSystem));
  region4[0] = -1.5;   region4[1] = -0.8;   region4[2] = 0.8;   region4[3] = 1.5; 
  fClassifiers.push_back(new AliEventClassifierMult("Eta08_15", "0.8 < |#eta| < 1.5",
   						    region4, 4, true, true, fMyOut, fGlobalSystem));

  ////////////////////////////////////
  // nMPI and Q^2 based classifiers //
  ////////////////////////////////////
  AliEventClassifierMPI* refClassifierMPI = new AliEventClassifierMPI("nMPI", "nMPI", fMyOut);
  AliEventClassifierQ2*  refClassifierQ2  = new AliEventClassifierQ2("Q2", "Q2", fMyOut);
  fClassifiers.push_back(refClassifierMPI);
  fClassifiers.push_back(refClassifierQ2);

  ///////////////////////////////
  // Spherocity and Sphericity //
  ///////////////////////////////
  // AliEventClassifierSpherocity* refClassifierSpherocity =
  //   new AliEventClassifierSpherocity("spherocity", "spherocity", fMyOut);
  //fClassifiers.push_back(refClassifierSpherocity);
  AliEventClassifierSphericity* refClassifierSphericity =
    new AliEventClassifierSphericity("sphericity", "sphericity", fMyOut);
  fClassifiers.push_back(refClassifierSphericity);

  // Set the global trigger if one was defined for this task
  // Remember to reset the classifiers used here as well (if new ones are defined just for the trigger)
  if (fGlobalTrigger == kINEL) {
    SetupInelAsGlobalTrigger();
    cout << "INEL global trigger set" << endl;
  }
  else if (fGlobalTrigger == kINELGT0) {
    SetupInelGt0AsGlobalTrigger(etaLt1);
    cout << "Inel > 0 global trigger set" << "\n";
  }
  else if (fGlobalTrigger == kV0AND) {
    SetupV0ANDAsGlobalTrigger(v0a, v0c);
    cout << "V0AND global trigger set" << "\n";
  }
  else
    AliError("Invalid global trigger name given");

  //////////////////////////////////////////
  // Connect classifiers with observables //
  //////////////////////////////////////////
  for (Int_t i=0; i < fClassifiers.size(); i++) {
     fObservables.push_back(new AliObservableEtaNch(fClassifiers.at(i)));
     fObservables.push_back(new AliObservableClassifierpTPID(fClassifiers.at(i)));
     fObservables.push_back(new AliObservableCorrelationsOfClassifiers(fClassifiers.at(i), refClassifierEtaLt05));
     fObservables.push_back(new AliObservableCorrelationsOfClassifiers(fClassifiers.at(i), refClassifierMPI));
     fObservables.push_back(new AliObservableCorrelationsOfClassifiers(fClassifiers.at(i), refClassifierQ2));
     //fObservables.push_back(new AliObservableCorrelationsOfClassifiers(fClassifiers.at(i), refClassifierSpherocity));
     fObservables.push_back(new AliObservableCorrelationsOfClassifiers(fClassifiers.at(i), refClassifierSphericity));
  }
  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fMyOut);
}


//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::UserExec(Option_t *)
{
  // Reset classifiers before doing anything with this new event
  for (Int_t i = 0; i < fClassifiers.size(); i++) {
    fClassifiers[i]->ResetClassifier();
  }

  // Load event
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     AliError("ERROR: Could not retrieve MC event");
     return;
  }
  AliStack  *stack = mcEvent->Stack();

  // do we have the right trigger?
  if (((fGlobalTrigger == kINEL) && IsInel(mcEvent, stack)) ||
      ((fGlobalTrigger == kINELGT0) && IsInelGt0(mcEvent, stack)) ||
      ((fGlobalTrigger == kV0AND) && IsV0AND(mcEvent, stack))) {
    for (Int_t i = 0; i < fObservables.size(); i++) {
      fObservables[i]->Fill(mcEvent, stack);
    }
  }
  
  // Post output data.
  PostData(1, fMyOut);
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::Terminate(Option_t *)
{
  //PostData(2, fRunconditions);
}



/////////////////////////////////////////////////////////////
// Setup the global triggers and their respose when called //
/////////////////////////////////////////////////////////////
void AliAnalysisTaskHMTFMCMultEst::SetupInelAsGlobalTrigger() {
  return;
}

void AliAnalysisTaskHMTFMCMultEst::SetupInelGt0AsGlobalTrigger(AliEventClassifierBase* etaLt1) {
  fGlobalTriggerClassifiers.push_back(etaLt1);
}

void AliAnalysisTaskHMTFMCMultEst::SetupV0ANDAsGlobalTrigger(AliEventClassifierBase* V0A, AliEventClassifierBase* V0C) {
  fGlobalTriggerClassifiers.push_back(V0A);
  fGlobalTriggerClassifiers.push_back(V0C);
}

/*
  Return true if the current event fulfills the trigger requiremtn
*/
Bool_t AliAnalysisTaskHMTFMCMultEst::IsInel(AliMCEvent *event, AliStack *stack) {
  return kTRUE;
}

Bool_t AliAnalysisTaskHMTFMCMultEst::IsInelGt0(AliMCEvent *event, AliStack *stack) {
  if (fGlobalTriggerClassifiers[0]->GetClassifierValue(event, stack) > 0)
    return kTRUE;
  else
    return kFALSE;
}

Bool_t AliAnalysisTaskHMTFMCMultEst::IsV0AND(AliMCEvent *event, AliStack *stack) {
  // The tow estimators in the vector are V0A and V0B
  if ((fGlobalTriggerClassifiers[0]->GetClassifierValue(event, stack) > 0)
      && (fGlobalTriggerClassifiers[1]->GetClassifierValue(event, stack) > 0))
    return kTRUE;
  else
    return kFALSE;
}

