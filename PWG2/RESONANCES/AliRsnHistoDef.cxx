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

#include <TH1.h>

#include "AliLog.h"
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
void AliRsnHistoDef::SetBins(Int_t n, Double_t min, Double_t max)
{
//
// Binning for histogram.
//
    fNBins = n;
    fMin = min;
    fMax = max;
    
    CheckEdges();
}

//_____________________________________________________________________________
void AliRsnHistoDef::CheckEdges()
{
//
// Checks that histogram edges are appropriate,
// otherwise swaps them.
//

    if (fMin > fMax) {
        AliWarning(Form("min = %f -- max = %f --> swapping", fMin, fMax));
        Double_t temp = fMin;
        fMin = fMax;
        fMax = temp;
    }
}

//_____________________________________________________________________________
TH1D* AliRsnHistoDef::CreateHistogram(const char *name, const char *title)
{
//
// Create a histogram with given name and title,
// whose binning is specified according to this class members.
//
    TH1D *histo = new TH1D(name, title, fNBins, fMin, fMax);
    return histo;
}
