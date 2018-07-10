#ifndef ALITOFTEMPLATEFITTER_H
#define ALITOFTEMPLATEFITTER_H
//Flags to set the modes to be used
// #define USECDECONVOLUTION//Cholesky-like
#define USEFITFUNCTIONS //Fit functions
class TH1;
class TH1F;
class TArrayD;
class TObjArray;
class TF1;
#include "Rtypes.h"

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Utilities for template fitting e.g. yield extraction                  //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

//_________________________________________________________________________________________________
Bool_t PerformFitWithTFF(TH1F* hData, TObjArray* mc, Double_t* range, Double_t* fitrange, TArrayD& fraction, TArrayD& fractionErr, TObjArray*& prediction);

//_________________________________________________________________________________________________
Double_t CHI2 = -1;
Bool_t PerformFitWithRooFit(TH1F* hData, TObjArray* mc, Double_t* range, Double_t* fitrange, TArrayD& fraction, TArrayD& fractionErr, TObjArray*& prediction, Double_t& chi2 = CHI2);

//_________________________________________________________________________________________________
Bool_t PerformFitWithFunctions(TH1F* hData, TObjArray* func, TF1* funcsum, Double_t* range, Double_t* fitrange, TArrayD& fraction, TArrayD& fractionErr, TObjArray*& prediction);

//_________________________________________________________________________________________________
Bool_t PerformFitWithCD(TH1F* hData, TObjArray* mc, Double_t* range, Double_t* fitrange, TArrayD& fraction, TArrayD& fractionErr, TObjArray*& prediction);

//_________________________________________________________________________________________________
Bool_t UseBinCounting(TH1F* hData, TObjArray* mc, Double_t* rangelol, TArrayD& fraction, TArrayD& fractionErr);

#endif
