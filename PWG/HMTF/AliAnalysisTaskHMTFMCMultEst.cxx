#include <iostream>

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
: AliAnalysisTaskSE(), fMyOut(0), fClassifiers(0), fObservables(0), fGlobalTriggerName(0),
  fGlobalTriggerMinValue(0), fGlobalTrigger(0)
{
}

//________________________________________________________________________
AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst(const char *name)
  : AliAnalysisTaskSE(name), fMyOut(0), fClassifiers(0), fObservables(0), fGlobalTriggerName(0),
    fGlobalTriggerMinValue(0), fGlobalTrigger(0)
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
  // Some classifieres are used as "reference" classifiers. Ie. We do not compute N! correlation
  // plots between N classifiers. We restrict oursefls to calculating only the correlations of every
  // classifier to a small number of reference classifiers.
  Float_t region2[] = {-0.5, 0.5};  // Use this array to pass two pass two edges to the constructor
  AliEventClassifierMult* refClassifierEtaLt05 = new AliEventClassifierMult("EtaLt05", "|#eta| < 0.5",
									    region2, 2, true, true, fMyOut);
  fClassifiers.push_back(refClassifierEtaLt05);
  region2[0] = -0.8; region2[1] = 0.8; 
  fClassifiers.push_back(new AliEventClassifierMult("EtaLt08", "|#eta| < 0.8", region2, 2, true, true, fMyOut));
  region2[0] = -1.5; region2[1] = 1.5; 
  fClassifiers.push_back(new AliEventClassifierMult("EtaLt15", "|#eta| < 0.8", region2, 2, true, true, fMyOut));
  region2[0] = 2.8; region2[1] = 5.1; 
  fClassifiers.push_back(new AliEventClassifierMult("V0A", "2.8 < #eta < 5.1", region2, 2, true, true, fMyOut));
  region2[0] = -3.7; region2[1] = -1.7; 
  fClassifiers.push_back(new AliEventClassifierMult("V0C", "-3.7 < #eta < -1.7", region2, 2, true, true, fMyOut));

  region2[0] = -8.7; region2[1] = 8.7; // not inclusive region, not charged!
  fClassifiers.push_back(new AliEventClassifierMult("ZDC", "|#eta| > 8.7", region2, 2, false, false, fMyOut));

  Float_t region4[] = {-3.7, -1.7, 2.8, 5.1};
  fClassifiers.push_back(new AliEventClassifierMult("V0M", "V0A + V0C", region4, 4, true, true, fMyOut));
  region4[0] = -1.5;   region4[1] = -0.8;   region4[2] = 0.8;   region4[3] = 1.5; 
  fClassifiers.push_back(new AliEventClassifierMult("Eta08_15", "0.8 < |#eta| < 1.5",
						    region4, 4, true, true, fMyOut));

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
  AliEventClassifierSphericity* refClassifierSphericity =
    new AliEventClassifierSphericity("sphericity", "sphericity", fMyOut);
  //fClassifiers.push_back(refClassifierSpherocity);
  fClassifiers.push_back(refClassifierSphericity);

  // Set the global trigger if one was defined for this task
  if (strcmp(fGlobalTriggerName.Data(), "InelGt0") == 0) {
    region2[0] = -1.0; region2[1] = 1.0;
    fGlobalTrigger = new AliEventClassifierMult("EtaLt1", "|#eta| < 1.0",
						region2, 2, true, true, fMyOut);
    fGlobalTriggerMinValue = 0.0;
    cout << "Inel > 0 global trigger set" << "\n";
  }
  else if (strcmp(fGlobalTriggerName.Data(), "SphericityGt09") == 0) {
    fGlobalTrigger = refClassifierSphericity;
    fGlobalTriggerMinValue = 0.9;
    cout << "Sphericity global trigger set" << "\n";
  }
  else
    cout << "No global trigger set" << endl;
  
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
    if (fGlobalTrigger)
      fGlobalTrigger->ResetClassifier();
  }

  // Load event
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     AliError("ERROR: Could not retrieve MC event");
     return;
  }
  AliStack  *stack = mcEvent->Stack();

  // Either no global trigger is defined or we require it to be fullfilled
  if (!fGlobalTrigger
      || (fGlobalTrigger->GetClassifierValue(mcEvent, stack) > fGlobalTriggerMinValue)) {
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
