#ifndef __ALICALORAWANALYZERGAMMA2_H__
#define __ALICALORAWANALYZERGAMMA2_H__


/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerGamma2
/// \ingroup EMCALraw
/// \brief  Raw data fitting: Gamma-2 function
///
/// Evaluation of amplitude and peak position using gamma-2 function.
/// Derivatives calculated analytically. 
/// Newton's method used for solving the set of non-linear equations.
///
/// \author Martin.Poghosyan@cern.ch, ORNL. 
//_________________________________________________________________________

#include "AliCaloRawAnalyzerFitter.h"



class AliCaloRawAnalyzerGamma2 :public AliCaloRawAnalyzerFitter {

public:
AliCaloRawAnalyzerGamma2();
~AliCaloRawAnalyzerGamma2(){;}


void SetNiterationsMax(Int_t n) {fNiterationsMax=n;}
Int_t GetNiterations()    {return fNiter;}
Int_t GetNiterationsMax() {return fNiterationsMax;}

virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector,
                                       UInt_t altrocfg1,
                                       UInt_t altrocfg2 );
  





private:

Int_t fNiter=0; ///< number of iteraions
Int_t fNiterationsMax=15; ///< max number of iteraions

Bool_t DoFit_1peak(Int_t ifirst, Int_t nS, Float_t &A0, Float_t &t0, Float_t &chi2);
void DoParabolaFit(Int_t x, Float_t &Amp, Float_t &time);


ClassDef(AliCaloRawAnalyzerGamma2,1)   // tbd

};
#endif
