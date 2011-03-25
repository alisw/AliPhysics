// -*- mode: c++ -*-
/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <perthomas.hille@yale.edu>            *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthomas.hille@yale.edu                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliCaloFitResults.h"
#include "AliCaloConstants.h"


AliCaloFitResults::AliCaloFitResults() : fMaxSig(0),
					 fPed(-1), 
					 fStatus(-1),
					 fAmpSig(-1),
					 fTime(-1),
					 fMaxTimebin(-1),
					 fChi2Sig(-1),
					 fNdfSig(0),
					 fMinSig(0),
					 fFitSubarray(10) 
{
  
  
}



// Container class to hold results from fitting 
// as well as other methods for
// raw data signals extraction. The class memebers
// fChi2Sig, fNdfSig is only relevant if a fitting procedure is
// Applied. fStatus holds information on wether or not 
// The signal was fitted sucessfully. fStatus might have a different meaning If other 
// procedures than  A different meaning Fitting is applied 


AliCaloFitResults::AliCaloFitResults(const Int_t maxSig, 
				     const Float_t ped, 
				     const Short_t fitstatus, 
				     const Float_t  amp,  
				     const double time,  
				     const Int_t maxTimebin, 
				     const Float_t chi,  
				     const Int_t ndf, 
				     Int_t minSig,
				     const AliCaloFitSubarray fitSubarray ) : 
  fMaxSig(maxSig),
  fPed(ped), 
  fStatus(fitstatus),
  fAmpSig(amp),
  fTime(time),
  fMaxTimebin(maxTimebin),
  fChi2Sig(chi),
  fNdfSig(ndf),
  fMinSig(minSig),
  fFitSubarray(fitSubarray) 
{

}



AliCaloFitResults::AliCaloFitResults(const Int_t maxSig, const Float_t ped, 
				     const Short_t fitstatus, const Float_t  amp,  
				     const double time, const Int_t maxTimebin, const Float_t chi,  
				     const Int_t ndf, Int_t minSig ) : 
  fMaxSig(maxSig),
  fPed(ped), 
  fStatus(fitstatus),
  fAmpSig(amp),
  fTime(time),
  fMaxTimebin(maxTimebin),
  fChi2Sig(chi),
  fNdfSig(ndf),
  fMinSig(minSig),
  fFitSubarray(Ret::kDummy)  
{
  
}



AliCaloFitResults::AliCaloFitResults(const Int_t maxSig, const Float_t ped, 
				     const Short_t fitstatus, const Float_t  amp,  
				     const Int_t maxTimebin) : 
  fMaxSig(maxSig),
  fPed(ped), 
  fStatus(fitstatus),
  fAmpSig(amp),
  fTime(maxTimebin),
  fMaxTimebin(maxTimebin),
  fChi2Sig( Ret::kNoFit ),
  fNdfSig( 0 ),
  fMinSig( 0 ),
  fFitSubarray( Ret::kNoFit ) 
{

}




AliCaloFitResults::AliCaloFitResults(const Int_t maxSig, const Int_t minSig) : 
  fMaxSig(maxSig),
  fPed( Ret::kInvalid ),
  fStatus( Ret::kInvalid ),
  fAmpSig( Ret::kInvalid ), 
  fTime( Ret::kInvalid ),
  fMaxTimebin( Ret::kInvalid ),
  fChi2Sig( Ret::kInvalid ),
  fNdfSig( 0 ),
  fMinSig (minSig),
  fFitSubarray( Ret::kInvalid )  
{

}
 


AliCaloFitResults::~AliCaloFitResults()
{

}

