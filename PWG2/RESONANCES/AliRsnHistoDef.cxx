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
#include <TH2.h>

#include "AliLog.h"

#include "AliRsnHistoDef.h"

ClassImp(AliRsnHistoDef)

//_____________________________________________________________________________
AliRsnHistoDef::AliRsnHistoDef() :
  fNDim(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnHistoDef::AliRsnHistoDef
(Int_t nbins, Double_t min, Double_t max) :
  fNDim(1)
{
//
// 1D histo definition.
//

    SetBins(nbins, min, max);
}
//_____________________________________________________________________________
AliRsnHistoDef::AliRsnHistoDef
(Int_t nbinsX, Double_t minX, Double_t maxX, Int_t nbinsY, Double_t minY, Double_t maxY):
  fNDim(2)
{
//
// 1D histo definition.
//

    SetBins(nbinsX, minX, maxX, nbinsY, minY, maxY);
}

//_____________________________________________________________________________
void AliRsnHistoDef::SetBins(Int_t n, Double_t min, Double_t max)
{
//
// Binning for 1D histograms.
//

    fNDim = 1;
    
    fNBins[0] = n;
    fMin[0] = min;
    fMax[0] = max;
    
    CheckEdges();
}

//_____________________________________________________________________________
void AliRsnHistoDef::SetBins
(Int_t nx, Double_t minx, Double_t maxx, Int_t ny, Double_t miny, Double_t maxy)
{
//
// Binning for 1D histograms.
//

    fNDim = 2;
    
    fNBins[0] = nx;
    fMin[0] = minx;
    fMax[0] = maxx;
    
    fNBins[1] = ny;
    fMin[1] = miny;
    fMax[1] = maxy;
    
    CheckEdges();
}

//_____________________________________________________________________________
void AliRsnHistoDef::CheckEdges()
{
//
// Checks that histogram edges are appropriate,
// otherwise swaps them.
//

    Int_t i;
    for (i = 0; i < fNDim; i++) {
        if (fMin[i] > fMax[i]) {
            AliWarning(Form("min = %f -- max = %f --> swapping", fMin, fMax));
            Double_t temp = fMin[i];
            fMin[i] = fMax[i];
            fMax[i] = temp;
        }
    }
}

//_____________________________________________________________________________
TH1D* AliRsnHistoDef::Create1DHistogram(const char *name, const char *title)
{
//
// Create 1D histogram, if the configuration is appropriate for this.
//

    if (fNDim != 1) {
        AliError("Number of dimension not set to 1!");
        return 0x0;
    }

    TH1D *histo = new TH1D(name, title, fNBins[0], fMin[0], fMax[0]);
    return histo;
}

//_____________________________________________________________________________
TH2D* AliRsnHistoDef::Create2DHistogram(const char *name, const char *title)
{
//
// Create 2D histogram, if the configuration is appropriate for this.
//

    if (fNDim != 2) {
        AliError("Number of dimension not set to 2!");
        return 0x0;
    }
    
    TH2D *histo = new TH2D(name, title, fNBins[0], fMin[0], fMax[0], fNBins[1], fMin[1], fMax[1]);
    return histo;
}
