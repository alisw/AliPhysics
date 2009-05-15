//
// Class AliRsnHistoDef
//
// Definition for a histogram type.
// Since one could do an analysis which is not an invariant mass
// the histogram definition should be more flexible, and it is stored
// separately in a new class.
// This class considers the possibility of a 1D or 2D histograms
// with its related binning, and can create a new histo from his definitions
//

#include <TObject.h>

//#include "AliLog.h"
#include "AliRsnHistoDef.h"

ClassImp(AliRsnHistoDef)

//_____________________________________________________________________________
AliRsnHistoDef::AliRsnHistoDef() :
    fNBins(0),
    fMin(0.0),
    fMax(0.0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnHistoDef::AliRsnHistoDef
(Int_t nbins, Double_t min, Double_t max) :
    fNBins(0),
    fMin(0.0),
    fMax(0.0)
{
//
// 1D histo definition.
//
  SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnHistoDef::AliRsnHistoDef
(Double_t min, Double_t max, Double_t step) :
    fNBins(0),
    fMin(0.0),
    fMax(0.0)
{
//
// 1D histo definition.
//
  SetBins(min, max, step);
}

//_____________________________________________________________________________
void AliRsnHistoDef::SetBins(Int_t n, Double_t min, Double_t max)
{
//
// Binning for histogram.
//

  fNBins = n;

  if (min < max)
  {
    fMin = min;
    fMax = max;
  }
  else
  {
    fMin = max;
    fMax = min;
  }
}

//_____________________________________________________________________________
void AliRsnHistoDef::SetBins(Double_t min, Double_t max, Double_t step)
{
//
// Binning for histogram.
//

  if (min < max)
  {
    fMin = min;
    fMax = max;
  }
  else
  {
    fMin = max;
    fMax = min;
  }

  fNBins = (Int_t)((fMax - fMin) / (step)) + 1;
}
