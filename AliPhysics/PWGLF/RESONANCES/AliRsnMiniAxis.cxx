//
// All implementations related to definition of an axis
// which is used in the output histogams.
// Simpler than TAxis, it defines an array of edges
// which is then ported to the output histogram definition.
// currently ported only in mini-package, but it could
// become a default also for general package.
//

#include "TMath.h"
#include "AliRsnMiniAxis.h"

ClassImp(AliRsnMiniAxis)

//_____________________________________________________________________________
void AliRsnMiniAxis::Set(Int_t nbins, Double_t min, Double_t max)
{
//
// Set binning for the axis in equally spaced bins
// where the number of bins, minimum and maximum are given.
//

   if (!nbins) {
      fBins.Set(0);
      return;
   }

   fBins.Set(nbins + 1);

   Double_t mymax = TMath::Max(min, max);
   Double_t mymin = TMath::Min(min, max);

   Int_t    k = 0;
   Double_t binSize = (mymax - mymin) / ((Double_t)nbins);

   fBins[0] = mymin;
   for (k = 1; k <= nbins; k++) fBins[k] = fBins[k - 1] + binSize;
}

//_____________________________________________________________________________
void AliRsnMiniAxis::Set(Double_t min, Double_t max, Double_t step)
{
//
// Set binning for the axis in equally spaced bins
// where the bin size, minimum and maximum are given.
//

   Double_t dblNbins = TMath::Abs(max - min) / step;
   Int_t    intNbins = ((Int_t)dblNbins) + 1;

   Set(intNbins, min, max);
}

//_____________________________________________________________________________
void AliRsnMiniAxis::Set(Int_t nbins, Double_t *array)
{
//
// Set binning for the axis in unequally spaced bins
// using the same way it is done in TAxis
//

   if (!nbins) {
      fBins.Set(0);
      return;
   }

   Int_t i;
   fBins.Set(nbins);
   for (i = 0; i < nbins; i++) fBins[i] = array[i];
}
