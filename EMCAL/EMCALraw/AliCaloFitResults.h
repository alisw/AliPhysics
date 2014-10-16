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
#include "AliCaloConstants.h"

class  AliCaloFitResults
{
 public:
  //  enum kReturnCode {kFitPar=1, kDummy=-1, kCrude=-9, kNoFit=-99, kInvalid=-9999};// possible return values
  // kFitPar: method fit or parametrization was used
  // kDummy: just a filler parameter, if e.g. chi2 not available
  // kCrude: maximum was used
  // kNoFit: maximum was used, exception handling for fit invoked
  // kInvalid: could not even look for maximum

  
  explicit AliCaloFitResults( const Int_t maxSig, 
			      const Float_t ped, 
			      const Short_t fitStatus, 
			      const Float_t  amp, 
			      const double time,
			      const Int_t maxTimebin,
			      //   const Float_t chi, 
			      const Float_t chi,  
			      const Int_t ndf, 
			      const Int_t minSig, 
			      const AliCaloFitSubarray fitSubarray  ); 
  

  explicit AliCaloFitResults( const Int_t maxSig, 
			      const Float_t ped, 
			      const Short_t fitStatus, 
			      const Float_t  amp, 
			      const double time,
			      const Int_t maxTimebin,
			      //   const Float_t chi, 
			      const Float_t chi, 
			      const Int_t ndf, 
			      const Int_t minSig = Ret::kDummy);  
  //			      const Int_t minSig = CaloConstants::ReturnCodes::kDummy); 


  // shorter interface when no fit is done

  
  explicit AliCaloFitResults( const Int_t maxSig, 
			      const Float_t ped, 
			      const Short_t fitStatus, 
			      const Float_t  amp, 
			      const Int_t maxTimebin); 
  

  // minimum interface
  explicit AliCaloFitResults( const Int_t maxSig, const Int_t minSig );

  AliCaloFitResults();
  virtual  ~AliCaloFitResults();
  UShort_t  GetMaxSig() const  { return fMaxSig;};
  Float_t   GetPed() const { return fPed;};
  UShort_t  GetMinSig() const { return fMinSig;};
  Int_t  GetStatus() const  { return fStatus;};
  Float_t   GetAmp() const {  return fAmpSig; };
  Float_t   GetTof() const {  return fTime; }; 
  double   GetTime() const {  return fTime; };
  Int_t   GetMaxTimebin() const {  return fMaxTimebin; };
  Float_t   GetChi2() const { return fChi2Sig;};
  UShort_t  GetNdf() const { return fNdfSig; };
  AliCaloFitSubarray  GetFitSubarray() const { return fFitSubarray; };
  void SetTime(Float_t time ) { fTime = time; };
  void SetAmp(Float_t amp ) { fAmpSig = amp; };

 private:
  // AliCaloFitResults();
  UShort_t   fMaxSig;      //Maximum sample value ( 0 - 1023 )
  Float_t    fPed;      //Pedestal 
  Int_t   fStatus;      //Sucess or failure of fitting pocedure
  Float_t    fAmpSig;   //Amplitude in entities of ADC counts
  double    fTime;     //peak/max time of signal in entities of sample intervals 
  Int_t    fMaxTimebin; //timebin with maximum ADC value
  Float_t    fChi2Sig;  //Chi Square of fit 
  UShort_t   fNdfSig;      //Number of degrees of freedom of fit
  UShort_t   fMinSig;      //Pedestal 
  AliCaloFitSubarray fFitSubarray; // info on time-bin array used for the fitting
};

#endif
