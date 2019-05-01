#ifndef UtilHisto_h
#define UtilHisto_h

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Utilities for histogram management et similia                         //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

class TH1;
class TH1F;
class TH2;
class TObjArray;

//_________________________________________________________________________________________________
//Function to check if two histograms have the same binning.
Bool_t SameBinning(const TH1* h1, const TH1* h2, const Bool_t check = kFALSE);

//_________________________________________________________________________________________________
//Function to get the starting point of the histogram with a minimum request on the count number
Int_t GetHistoLowRange(TH1F* h, const Int_t nmin, const Int_t ncounts = 3, const Int_t secondmin = 1);

//_________________________________________________________________________________________________
//Function to get the first histogram of a series
Int_t GetFirstHistogram(const TObjArray* templates, const Int_t nmin = 10);

//_________________________________________________________________________________________________
//Function to get the starting point after a certain value
Int_t GetHistoLowRangeAfter(TH1F* h, const Int_t nmin, const Int_t after = 1);

//_________________________________________________________________________________________________
//Function to get the starting point after a certain value, without any holes in the histogram
Int_t GetHistoNoHolesAfter(TH1F* h, const Int_t after = 1, const Int_t check = 20);

//_________________________________________________________________________________________________
//Function to get the starting point before a certain value, without any holes in the histogram
Int_t GetHistoNoHolesBefore(TH1F* h, const Double_t before = 0);

//_________________________________________________________________________________________________
//Function to get the last histogram of a series
Int_t GetHistoNoHolesBefore(TObjArray* h, const Double_t before = 0);

//_________________________________________________________________________________________________
//Function to check if the histogram is in range
Bool_t IsHistogramInRange(TH1* h, const Double_t rangelow, const Double_t rangehigh, const Double_t threshold = 1., const Bool_t verbose = kFALSE);

//_________________________________________________________________________________________________
//Function to get the residual yield from a fit, within a certain range
Double_t GetResidualYield(TH1* hdata, TH1* hfit, const Double_t rangelow, const Double_t rangehigh);

//_________________________________________________________________________________________________
//Function to get the overlap fraction of two histograms
Double_t GetOverlapFraction(const TH1F* h, const TH1F* bkg, Double_t& error, const Bool_t show = kFALSE);

//_________________________________________________________________________________________________
//Function to clean a histogram, i.e. remove the holes
TH1* Clean(TH1* h, const Float_t center, const Bool_t add, const Int_t nchain = 1, const Float_t thr = 1);

#endif
