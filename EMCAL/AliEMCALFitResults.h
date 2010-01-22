// -*- mode: c++ -*-
#ifndef ALIEMCALFITRESULTS_H
#define ALIEMCALFITRESULTS_H
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

// Container class to hold results from fitting 
// as well as other methods for
// raw data signals extraction
class  AliEMCALFitResults
{
 public:
  explicit AliEMCALFitResults( const UShort_t maxSig, 
 			       const Float_t ped, 
			       const Short_t fitStatus, 
			       const Float_t  amp, 
			       const Float_t t0,
			       const Float_t chi, 
			       const UShort_t ndf, 
			       const UShort_t minSig = -99); 

  explicit AliEMCALFitResults( const UShort_t maxSig, const UShort_t minSig );
  //AliEMCALFitResults( const UShort_t maxSig, const UShort_t minSig );


  virtual  ~AliEMCALFitResults();
  UShort_t  GetMaxSig() const  { return fMaxSig;};
  Float_t   GetPed() const { return fPed;};
  UShort_t  GetMinSig() const { return fMinSig;};
  UShort_t  GetStatus() const  { return fStatus;};
  Float_t   GetAmp() const {  return fAmpSig; };
  Float_t   GetTof() const {  return fT0; }; 
  Float_t   GetChisSquare() const { return fChi2Sig;};
  UShort_t  GetNdf() const { return fNdfSig; };
  
 private:
  AliEMCALFitResults();
  UShort_t   fMaxSig;   //Maximum sample value ( 0 - 1023 )
  Float_t    fPed;      //Pedestal 
  UShort_t   fStatus;   //Sucess or failure of fitting pocedure
  Float_t    fAmpSig;   //Amplitude in entities of ADC counts
  Float_t    fT0;       //Start time of signal in entities of sample intervals 
  Float_t    fChi2Sig;  //Chi Square of fit 
  UShort_t   fNdfSig;   //Number of degrees of freedom of fit
  UShort_t   fMinSig;   //Pedestal 
};

#endif
