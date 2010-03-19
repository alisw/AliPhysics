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


#include "AliCaloFitResults.h"


// Container class to hold results from fitting 
// as well as other methods for
// raw data signals extraction. The class memebers
// fChi2Sig, fNdfSig is only relevant if a fitting procedure is
// Applied. fStatus holds information on wether or not 
// The signal was fitted sucessfully. fStatus might have a different meaning If other 
// procedures than  A different meaning Fitting is applied 

AliCaloFitResults::AliCaloFitResults(const Int_t maxSig, const Float_t ped, 
				     const Short_t fitstatus, const Float_t  amp,  
				     const Float_t time,  const Int_t maxTimebin, const Float_t chi,  
				     const Int_t ndf, Int_t minSig,
				     const AliCaloFitSubarray fitSubarray) : 
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
				     const Float_t time, const Int_t maxTimebin, const Float_t chi,  
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
  fFitSubarray(kDummy)  
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
  fChi2Sig( kNoFit ),
  fNdfSig( kNoFit ),
  fMinSig( kNoFit ),
  fFitSubarray( kNoFit ) 
{
}


AliCaloFitResults::AliCaloFitResults(const Int_t maxSig, const Int_t minSig) : 
  fMaxSig(maxSig),
  fPed( kInvalid ),
  fStatus( kInvalid ),
  fAmpSig( kInvalid ), 
  fTime( kInvalid ),
  fMaxTimebin( kInvalid ),
  fChi2Sig( kInvalid ),
  fNdfSig( kInvalid),
  fMinSig (minSig),
  fFitSubarray(kInvalid)  
{

}



AliCaloFitResults::~AliCaloFitResults()
{

}

