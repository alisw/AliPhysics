// -*- mode: c++ -*-
#ifndef ALICALOFITRESULTS_H
#define ALICALOFITRESULTS_H
/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "Rtypes.h"

#include "AliCaloFitSubarray.h"

// Container class to hold results from fitting 
// as well as other methods for
// raw data signals extraction
class  AliCaloFitResults
{
 public:
  enum kReturnCode {kDummy=-1, kCrude=-9, kNoFit=-99, kInvalid=-9999};// possible return values

  explicit AliCaloFitResults( const Int_t maxSig, 
			      const Float_t ped, 
			      const Short_t fitStatus, 
			      const Float_t  amp, 
			      const Float_t t0,
			      const Float_t chi, 
			      const Int_t ndf, 
			      const Int_t minSig, 
			      const AliCaloFitSubarray fitSubarray); 

  explicit AliCaloFitResults( const Int_t maxSig, 
			      const Float_t ped, 
			      const Short_t fitStatus, 
			      const Float_t  amp, 
			      const Float_t t0,
			      const Float_t chi, 
			      const Int_t ndf, 
			      const Int_t minSig = kDummy); 

  explicit AliCaloFitResults( const Int_t maxSig, const Int_t minSig );
  //AliCaloFitResults( const Int_t maxSig, const Int_t minSig );


  virtual  ~AliCaloFitResults();
  Int_t  GetMaxSig() const  { return fMaxSig;};
  Float_t   GetPed() const { return fPed;};
  Int_t  GetMinSig() const { return fMinSig;};
  Int_t  GetStatus() const  { return fStatus;};
  Float_t   GetAmp() const {  return fAmpSig; };
  Float_t   GetTof() const {  return fT0; }; 
  Float_t   GetChisSquare() const { return fChi2Sig;};
  Int_t  GetNdf() const { return fNdfSig; };
  AliCaloFitSubarray  GetFitSubarray() const { return fFitSubarray; };
  
 private:
  AliCaloFitResults();
  Int_t   fMaxSig;   //Maximum sample value ( 0 - 1023 )
  Float_t    fPed;      //Pedestal 
  Int_t   fStatus;   //Sucess or failure of fitting pocedure
  Float_t    fAmpSig;   //Amplitude in entities of ADC counts
  Float_t    fT0;       //Start time of signal in entities of sample intervals 
  Float_t    fChi2Sig;  //Chi Square of fit 
  Int_t   fNdfSig;   //Number of degrees of freedom of fit
  Int_t   fMinSig;   //Pedestal 
  AliCaloFitSubarray fFitSubarray; // info on time-bin array used for the fitting
};

#endif
