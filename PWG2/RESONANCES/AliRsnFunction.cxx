//
// Class AliRsnFunction
//
// This class defines a base classe to implement a function
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from
// this class and override the virtual functions defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#include <TString.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnPairParticle.h"
#include "AliRsnFunctionAxis.h"

#include "AliRsnFunction.h"

ClassImp(AliRsnFunction)

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction() :
    TNamed(),
    fPairDef(0x0),
    fAxisList("AliRsnFunctionAxis", 0),
    fTrack(0x0),
    fPair(0x0),
    fEvent(0x0),
    fHistogram(0x0)
{
//
// Constructor.
//
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction(const AliRsnFunction &copy) :
    TNamed(copy),
    fPairDef(copy.fPairDef),
    fAxisList(copy.fAxisList),
    fTrack(copy.fTrack),
    fPair(copy.fPair),
    fEvent(copy.fEvent),
    fHistogram(0x0)
{
//
// Copy constructor.
//
}

//________________________________________________________________________________________
const AliRsnFunction& AliRsnFunction::operator=(const AliRsnFunction& copy)
{
//
// Assignment operator.
//

  SetName(copy.GetName());
  SetTitle(copy.GetTitle());

  fPairDef = copy.fPairDef;
  fTrack = copy.fTrack;
  fPair = copy.fPair;
  fEvent = copy.fEvent;

  if (fHistogram) delete fHistogram;
  fHistogram = 0x0;

  return (*this);
}

//________________________________________________________________________________________
const char* AliRsnFunction::GetName() const
{
//
// Defines the name of this object according to
// the function type and binning
//

  TString name("");

  TObjArrayIter next(&fAxisList);
  AliRsnFunctionAxis *axis = 0;

  while ((axis = (AliRsnFunctionAxis*)next())) {
    if (name.Length() > 1) name += '_';
    name += axis->GetName();
  }

  return name.Data();
}

//________________________________________________________________________________________
void AliRsnFunction::AddAxis(AliRsnFunctionAxis *const axis)
{
  Int_t size = fAxisList.GetEntries();
  new(fAxisList[size]) AliRsnFunctionAxis(*axis);
}

//________________________________________________________________________________________
THnSparseD* AliRsnFunction::CreateHistogram(const char *histoName, const char *histoTitle)
{
//
// Creates and returns the histogram defined using
// arguments fo name and title, and the first histoDef for binning.
// Variable-sized histogram binning is always called, due to use of histoDef,
// even if the bins are equal, since they are defined in this class.
// Eventually present histoDef's in other slots of array (1, 2) are ignored.
//

  Int_t size = fAxisList.GetEntries();
  if (!size) {
    AliError("No axes defined");
    return 0x0;
  }

  Int_t    *nbins = new Int_t[size];
  Double_t *min   = new Double_t[size];
  Double_t *max   = new Double_t[size];

  // retrieve binnings for main and secondary axes
  AliRsnFunctionAxis *fcnAxis = 0;
  for (Int_t i = 0; i < size; i++) {
    fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(i);
    if (!fcnAxis) {
      nbins[i] = 0;
      min[i] = 0.0;
      max[i] = 0.0;
      AliError("Empty axis");
      continue;
    }
    nbins[i] = fcnAxis->GetNBins();
    min[i] = fcnAxis->GetMin();
    max[i] = fcnAxis->GetMax();
  }

  // create histogram
  fHistogram = new THnSparseD(histoName, histoTitle, size, nbins, min, max);
  fHistogram->Sumw2();

  return fHistogram;
}

//________________________________________________________________________________________
Bool_t AliRsnFunction::Fill()
{
//
// Fill function histogram with values computed from given input object.
//

  AliDebug(AliLog::kDebug +2,"->");

  Int_t  i, nAxes = fAxisList.GetEntries();
  Double_t *values = new Double_t[nAxes];

  AliRsnFunctionAxis *fcnAxis = 0;
  for (i = 0; i < nAxes; i++) {
    fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(i);
    if (!fcnAxis) {
      values[i] = 0.0;
      continue;
    }
    switch (fcnAxis->GetAxisObject()) {
    case AliRsnFunctionAxis::kParticle:
      values[i] = fcnAxis->Eval(fTrack);
      break;
    case AliRsnFunctionAxis::kPair:
      values[i] = fcnAxis->Eval(fPair, fPairDef);
      break;
    case AliRsnFunctionAxis::kEvent:
      values[i] = fcnAxis->Eval(fEvent);
      break;
    default:
      values[i] = 0.0;
    }
  }

  // check presence of output histogram
  if (!fHistogram) {
    AliError("Histogram is not yet initialized");
    return kFALSE;
  }
  //TArrayD val(values); val->Print();
  fHistogram->Fill(values);

  AliDebug(AliLog::kDebug +2,"->");
  return kTRUE;
}
