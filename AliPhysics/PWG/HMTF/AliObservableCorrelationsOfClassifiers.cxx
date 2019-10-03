#include <iostream>

#include "TString.h"
#include "TH3F.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

#include "AliObservableBase.h"
#include "AliObservableCorrelationsOfClassifiers.h"
#include "AliEventClassifierBase.h"

#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliObservableCorrelationsOfClassifiers)

AliObservableCorrelationsOfClassifiers::AliObservableCorrelationsOfClassifiers()
  : AliObservableBase()
{
}


AliObservableCorrelationsOfClassifiers::AliObservableCorrelationsOfClassifiers(AliEventClassifierBase *classifier0,
									       AliEventClassifierBase *classifier1)
  :  AliObservableBase("correlationsOfClassifiers", "Correlation between two classifiers")
{
  // Correlations between the values of classifier0 and 1. The results are written to the folder of the former.
  // Hence, classifier1 can be regarded as the reference classifier, with which each other needs to be correlated.
  
  fclassifier0 = classifier0;
  fclassifier1 = classifier1;
  const Int_t classifier_bins  = 250;

  fhistogram = new TH2F("corr_this_with_" + TString(classifier1->GetName()),
			"#eta vs. classifier vs. N_{ch}",
			classifier_bins, fclassifier0->GetExpectedMinValue(), fclassifier0->GetExpectedMaxValue(),
			classifier_bins, fclassifier1->GetExpectedMinValue(), fclassifier1->GetExpectedMaxValue());
  fhistogram->GetXaxis()->SetTitle(classifier0->GetName());
  fhistogram->GetYaxis()->SetTitle(classifier1->GetName());
  fhistogram->GetZaxis()->SetTitle("Event count");
  fhistogram->Sumw2();
  fhistogram->SetDirectory(0);
  
  // Construct array for classifier values
  classifier0->GetClassifierOutputList()->Add(fhistogram);
}

void AliObservableCorrelationsOfClassifiers::Fill(AliMCEvent *event, AliStack *stack) {
  Float_t cls0 = fclassifier0->GetClassifierValue(event, stack);
  Float_t cls1 = fclassifier1->GetClassifierValue(event, stack);
  Float_t weight = event->GenEventHeader()->EventWeight();
  fhistogram->Fill(cls0, cls1, weight);
}
