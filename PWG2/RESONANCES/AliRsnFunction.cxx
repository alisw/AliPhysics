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
AliRsnFunction::AliRsnFunction(Bool_t useTH1) :
    TNamed(),
    fPairDef(0x0),
    fAxisList("AliRsnFunctionAxis", 0),
    fTrack(0x0),
    fPair(0x0),
    fEvent(0x0),
    fUseTH1(useTH1),
    fSize(0),
    fH1(0x0),
    fHSparse(0x0)
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
    fUseTH1(copy.fUseTH1),
    fSize(copy.fSize),
    fH1(0x0),
    fHSparse(0x0)
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
  fUseTH1 = copy.fUseTH1;
  fSize = copy.fSize;

  if (fH1) delete fH1;
  fH1 = 0x0;
  
  if (fHSparse) delete fHSparse;
  fHSparse = 0x0;

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
  if (size >= 3 && fUseTH1)
  {
    AliWarning("A TH1-type output cannot add more than 3 axes: switching to THnSparse -- THIS COULD CAUSE VERY LARGE FILES!!!");
    fUseTH1 = kFALSE;
  }
  new(fAxisList[size]) AliRsnFunctionAxis(*axis);
}

//________________________________________________________________________________________
TH1* AliRsnFunction::CreateHistogram(const char *histoName, const char *histoTitle)
{
//
// Creates and returns the histogram defined using
// arguments fo name and title, and the first histoDef for binning.
// Variable-sized histogram binning is always called, due to use of histoDef,
// even if the bins are equal, since they are defined in this class.
// Eventually present histoDef's in other slots of array (1, 2) are ignored.
//
// This version produces a THnSparseD.
//

  fSize = fAxisList.GetEntries();
  if (!fSize) {
    AliError("No axes defined");
    return 0x0;
  }
  else if (fSize < 1 || fSize > 3)
  {
    AliError("Too few or too many axes defined");
    return 0x0;
  }

  Int_t    *nbins = new Int_t   [fSize];
  Double_t *min   = new Double_t[fSize];
  Double_t *max   = new Double_t[fSize];

  // retrieve binnings for main and secondary axes
  AliRsnFunctionAxis *fcnAxis = 0;
  for (Int_t i = 0; i < fSize; i++) {
    fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(i);
    if (!fcnAxis) {
      nbins[i] = 0;
      min[i]   = 0.0;
      max[i]   = 0.0;
      AliError("Empty axis");
      continue;
    }
    nbins[i] = fcnAxis->GetNBins();
    min[i]   = fcnAxis->GetMin();
    max[i]   = fcnAxis->GetMax();
  }

  // create histogram depending on the number of axes
  switch (fSize)
  {
    case 1:
      fH1 = new TH1D(histoName, histoTitle, nbins[0], min[0], max[0]);
      break;
    case 2:
      fH1 = new TH2D(histoName, histoTitle, nbins[0], min[0], max[0], nbins[1], min[1], max[1]);
      break;
    case 3:
      fH1 = new TH3D(histoName, histoTitle, nbins[0], min[0], max[0], nbins[1], min[1], max[1], nbins[2], min[2], max[2]);
      break;
  }
  fH1->Sumw2();

  return fH1;
}

//________________________________________________________________________________________
THnSparseD* AliRsnFunction::CreateHistogramSparse(const char *histoName, const char *histoTitle)
{
//
// Creates and returns the histogram defined using
// arguments fo name and title, and the first histoDef for binning.
// Variable-sized histogram binning is always called, due to use of histoDef,
// even if the bins are equal, since they are defined in this class.
// Eventually present histoDef's in other slots of array (1, 2) are ignored.
//
// This version produces a THnSparseD.
//

  fSize = fAxisList.GetEntries();
  if (!fSize) {
    AliError("No axes defined");
    return 0x0;
  }

  Int_t    *nbins = new Int_t   [fSize];
  Double_t *min   = new Double_t[fSize];
  Double_t *max   = new Double_t[fSize];

  // retrieve binnings for main and secondary axes
  AliRsnFunctionAxis *fcnAxis = 0;
  for (Int_t i = 0; i < fSize; i++) {
    fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(i);
    if (!fcnAxis) {
      nbins[i] = 0;
      min[i]   = 0.0;
      max[i]   = 0.0;
      AliError("Empty axis");
      continue;
    }
    nbins[i] = fcnAxis->GetNBins();
    min[i]   = fcnAxis->GetMin();
    max[i]   = fcnAxis->GetMax();
  }

  Int_t size = fAxisList.GetEntries();
  if (!size) {
    AliError("No axes defined");
    return 0x0;
  }

  // create histogram
  fHSparse = new THnSparseD(histoName, histoTitle, size, nbins, min, max);
  fHSparse->Sumw2();
  
  // clean heap
  delete [] nbins;
  delete [] min;
  delete [] max;

  return fHSparse;
}


//________________________________________________________________________________________
Bool_t AliRsnFunction::Fill()
{
//
// Fill function histogram with values computed from given input object.
//

  AliDebug(AliLog::kDebug +2,"->");

  Int_t  i;
  Double_t *values = new Double_t[fSize];

  AliRsnFunctionAxis *fcnAxis = 0;
  for (i = 0; i < fSize; i++) {
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
  
  // fill histogram
  if (fUseTH1)
  {
    // check presence of output histogram
    if (!fH1) {
      AliError("Required a TH1 whish is not initialized");
      return kFALSE;
    }
    
    // fill according to dimensions
    switch (fSize)
    {
      case 1:
        {
          TH1D *h1 = (TH1D*)fH1;
          h1->Fill(values[0]);
        }
        break;
      case 2:
        {
          TH2D *h2 = (TH2D*)fH1;
          h2->Fill(values[0], values[1]);
        }
        break;
      case 3:
        {
          TH3D *h3 = (TH3D*)fH1;
          h3->Fill(values[0], values[1], values[2]);
        }
        break;
      default:
        AliError(Form("Wrong size : %d", fSize));
        return kFALSE;
    }
  }
  else
  {
    // check presence of output histogram
    if (!fHSparse) {
      AliError("Required a THnSparseD whish is not initialized");
      return kFALSE;
    }
    
    fHSparse->Fill(values);
  }
  
  delete [] values;

  AliDebug(AliLog::kDebug +2,"->");
  return kTRUE;
}
