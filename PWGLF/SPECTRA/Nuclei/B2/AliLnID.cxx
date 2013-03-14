/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// class for light nuclei identification
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <AliExternalTrackParam.h>
#include <TParticle.h>
#include <AliESDtrack.h>
#include <AliPID.h>
#include <AliITSPIDResponse.h>
#include <AliTPCPIDResponse.h>
#include <TPDGCode.h>
#include "AliLnID.h"

ClassImp(AliLnID)

AliLnID::AliLnID()
: TObject()
, fPidProcedure(kBayes)
, fRange(5.)
, fZexp(2)
, fITSpid(0)
, fTPCpid(0)
{
//
// Default constructor
//
	fSpecies[0] = AliPID::kElectron;
	fSpecies[1] = AliPID::kMuon;
	fSpecies[2] = AliPID::kPion;
	fSpecies[3] = AliPID::kKaon;
	fSpecies[4] = AliPID::kProton;
	fSpecies[5] = AliPID::kDeuteron;
	fSpecies[6] = AliPID::kTriton;
	fSpecies[7] = AliPID::kHe3;
	fSpecies[8] = AliPID::kAlpha;
	
	// equal prior probabilities for all particle species
	for(Int_t i=0; i<kSPECIES; ++i) fPrior[i]=1./kSPECIES;
	
	// ALEPH Bethe Bloch parameters for ITS
	Double_t param[] = { 0.13, 15.77/0.95, 4.95, 0.312, 2.14, 0.82};
	fITSpid = new AliITSPIDResponse(&param[0]);
	
	fTPCpid = new AliTPCPIDResponse();
	fTPCpid->SetSigma(0.06,0);
}

AliLnID::AliLnID(const AliLnID& other)
: TObject(other)
, fPidProcedure(other.fPidProcedure)
, fRange(other.fRange)
, fZexp(other.fZexp)
, fITSpid(0)
, fTPCpid(0)
{
//
// Copy constructor
//
	for(Int_t i=0; i < kSPECIES; ++i)
	{
		fSpecies[i] = other.fSpecies[i];
		fPrior[i]   = other.fPrior[i];
	}
	
	fITSpid = new AliITSPIDResponse(*(other.fITSpid));
	fTPCpid = new AliTPCPIDResponse(*(other.fTPCpid));
}

AliLnID& AliLnID::operator=(const AliLnID& other)
{
//
// Assignment operator
//
	if(&other == this) return *this; // check for self-assignment
	
	fPidProcedure = other.fPidProcedure;
	fRange = other.fRange;
	fZexp = other.fZexp;
	
	// deallocate memory
	delete fITSpid;
	delete fTPCpid;
	
	// copy
	TObject::operator=(other);
	
	for(Int_t i=0; i < kSPECIES; ++i)
	{
		fSpecies[i] = other.fSpecies[i];
		fPrior[i]   = other.fPrior[i];
	}
	
	fITSpid = new AliITSPIDResponse(*(other.fITSpid));
	fTPCpid = new AliTPCPIDResponse(*(other.fTPCpid));
	
	return *this;
}

void AliLnID::SetPriorProbabilities(Double_t e, Double_t mu, Double_t pi, Double_t k, Double_t p, Double_t d, Double_t t, Double_t he3, Double_t alpha)
{
//
// Set prior probabilities
//
	fPrior[0] = e; // electron
	fPrior[1] = mu; // muon
	fPrior[2] = pi; // pion
	fPrior[3] = k; // kaon
	fPrior[4] = p; // proton
	fPrior[5] = d; // deuteron
	fPrior[6] = t; // triton
	fPrior[7] = he3; // he3
	fPrior[8] = alpha; // alpha
}

void AliLnID::SetPriorProbabilities(const Double_t* prob)
{
//
// Set prior probabilities
//
	for(Int_t i=0; i < kSPECIES; ++i) fPrior[i] = prob[i];
}

void AliLnID::SetITSBetheBlochParams(const Double_t* par, Double_t res)
{
//
// Set ALEPH Bethe Bloch parameters for ITS
//
	delete fITSpid;
	Double_t param[] = { res, par[0], par[1], par[2], par[3], par[4]};
	fITSpid = new AliITSPIDResponse(&param[0]);
}

void AliLnID::SetTPCBetheBlochParams(const Double_t* par, Double_t mip, Double_t res)
{
//
// Set ALEPH Bethe Bloch parameters for TPC
//
	fTPCpid->SetMip(mip);
	fTPCpid->SetSigma(res,0);
	fTPCpid->SetBetheBlochParameters(par[0], par[1], par[2], par[3], par[4]);
}

void AliLnID::SetTPCBetheBlochParams(Double_t c0, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t mip, Double_t res)
{
//
// Set ALEPH Bethe Bloch parameters for TPC
//
	fTPCpid->SetMip(mip);
	fTPCpid->SetSigma(res,0);
	fTPCpid->SetBetheBlochParameters(c0, c1, c2, c3, c4);
}

AliLnID::~AliLnID()
{
//
// Default destructor
//
	delete fITSpid;
	delete fTPCpid;
}

Int_t AliLnID::GetPID(const TParticle* p) const
{
//
// Montecarlo PID
//
	enum { kDeuteron=1000010020, kTriton=1000010030, kHelium3=1000020030, kAlpha=1000020040 };
	
	switch(TMath::Abs(p->GetPdgCode()))
	{
		case kElectron:  return AliPID::kElectron;
		case kMuonMinus: return AliPID::kMuon;
		case kPiPlus:    return AliPID::kPion;
		case kKPlus:     return AliPID::kKaon;
		case kProton:    return AliPID::kProton;
		case kDeuteron:  return AliPID::kDeuteron;
		case kTriton:    return AliPID::kTriton;
		case kHelium3:   return AliPID::kHe3;
		case kAlpha:     return AliPID::kAlpha;
	}
	
	return -1;
}

Int_t AliLnID::GetPID(Int_t pidCode, Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta, Double_t nSigITS, Double_t nSigTPC, Double_t nSigTOF) const
{
//
// PID according to the selected procedure
//
	if(fPidProcedure == kBayes)
	{
		return this->GetBayesPID(pITS, dEdxITS, nPointsITS, pTPC, dEdxTPC, nPointsTPC, pTOF, beta);
	}
	else if(fPidProcedure == kMaxLikelihood)
	{
		return this->GetMaxLikelihoodPID(pITS, dEdxITS, nPointsITS, pTPC, dEdxTPC, nPointsTPC, pTOF, beta);
	}
	else if(fPidProcedure == kTPC)
	{
		return this->GetTPCpid(pidCode, pTPC, dEdxTPC, nPointsTPC, nSigTPC);
	}
	else if(fPidProcedure == kITSTPC)
	{
		return this->GetITSTPCpid(pidCode, pITS, dEdxITS, nPointsITS, nSigITS, pTPC, dEdxTPC, nPointsTPC, nSigTPC);
	}
	else if(fPidProcedure == kTPCTOF)
	{
		return this->GetTPCTOFpid(pidCode, pTPC, dEdxTPC, nPointsTPC, nSigTPC, pTOF, beta, nSigTOF);
	}
	
	return -1;
}

Int_t AliLnID::GetTPCpid(Int_t pidCode, Double_t pTPC, Double_t dEdxTPC, Double_t nPointsTPC, Double_t nSigmaTPC) const
{
//
// Check if particle with the given pid code is within
// +/- nSigma around the expected dEdx in the TPC
//
	if( this->GetTPCmatch(pidCode, pTPC, dEdxTPC, nPointsTPC, nSigmaTPC)) return pidCode;
	
	return -1;
}

Int_t AliLnID::GetTPCTOFpid(Int_t pidCode, Double_t pTPC, Double_t dEdxTPC, Double_t nPointsTPC, Double_t nSigmaTPC, Double_t pTOF, Double_t beta, Double_t nSigmaTOF) const
{
//
// Check if particle with the given pid code is within
// +/- nSigmaTPC around the expected dEdx in the TPC
// AND +/- nSigmaTOF around the expected beta in the TOF (when available)
//
	if(  beta > 0 )
	{
		if( this->GetTPCmatch(pidCode, pTPC, dEdxTPC, nPointsTPC, nSigmaTPC)
		 && this->GetTOFmatch(pidCode, pTOF, beta, nSigmaTOF))
		{
			return pidCode;
		}
	}
	else
	{
		if( this->GetTPCmatch(pidCode, pTPC, dEdxTPC, nPointsTPC, nSigmaTPC))
		{
			return pidCode;
		}
	}
	
	return -1;
}

Int_t AliLnID::GetITSTPCpid(Int_t pidCode, Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t nSigmaITS, Double_t pTPC, Double_t dEdxTPC, Double_t nPointsTPC, Double_t nSigmaTPC) const
{
//
// Check if particle with the given pid code is within
// +/- nSigmaITS around the expected dEdx in the ITS
// AND +/- nSigmaTPC around the expected dEdx in the TPC
//
	if( this->GetITSmatch(pidCode, pITS, dEdxITS, nPointsITS, nSigmaITS)
	  && this->GetTPCmatch(pidCode, pTPC, dEdxTPC, nPointsTPC, nSigmaTPC))
	{
		return pidCode;
	}
	
	return -1;
}

Int_t AliLnID::GetIndexOfMaxValue(const Double_t* w) const
{
//
// Index with maximum value in the array of size kSPECIES
//
	Double_t wmax= 0;
	Int_t imax = -1;
	for(Int_t i=0; i < kSPECIES; ++i)
	{
		if(w[i] > wmax)
		{
			wmax = w[i];
			imax = i;
		}
	}
	
	return imax;
}

Bool_t AliLnID::GetLikelihood(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta, Double_t* r) const
{
//
// Fill r array with the combined likelihood for ITS, TPC and TOF when possible
// return 1 if success
//
	Double_t its[kSPECIES];
	Double_t tpc[kSPECIES];
	Double_t tof[kSPECIES];
	
	Bool_t itsPid = this->GetITSlikelihood(pITS, dEdxITS, nPointsITS, its);
	Bool_t tpcPid = this->GetTPClikelihood(pTPC, dEdxTPC, nPointsTPC, tpc);
	Bool_t tofPid = this->GetTOFlikelihood(pTOF, beta, tof);
	
	if(!itsPid && !tpcPid && !tofPid) return kFALSE;
	
	// Combine detector responses
	
	for(Int_t i=0; i<kSPECIES; ++i) r[i]=1.;
	
	if(itsPid)
	{
		for(Int_t i=0; i < kSPECIES; ++i) r[i] *= its[i];
	}
	if(tpcPid)
	{
		for(Int_t i=0; i < kSPECIES; ++i) r[i] *= tpc[i];
	}
	if(tofPid)
	{
		for(Int_t i=0; i < kSPECIES; ++i) r[i] *= tof[i];
	}
	
	return kTRUE;
}

Int_t AliLnID::GetMaxLikelihoodPID(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta) const
{
//
// Maximum likelihood principle
//
	Double_t r[kSPECIES];
	
	if(!this->GetLikelihood(pITS, dEdxITS, nPointsITS, pTPC, dEdxTPC, nPointsTPC, pTOF, beta, r)) return -1;
	
	Int_t imax = this->GetIndexOfMaxValue(r);
	
	if(imax < 0) return -1;
	
	return fSpecies[imax];
}

Int_t AliLnID::GetBayesPID(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta) const
{
//
// Bayesian inference
//
	Double_t r[kSPECIES];
	
	if(!this->GetLikelihood(pITS, dEdxITS, nPointsITS, pTPC, dEdxTPC, nPointsTPC, pTOF, beta, r)) return -1;
	
	// Bayes' rule
	
	Double_t w[kSPECIES] = {0};
	for(Int_t i = 0; i < kSPECIES; ++i) w[i] = r[i]*fPrior[i];  // no need to normalize
	
	Int_t imax = this->GetIndexOfMaxValue(w);
	
	if(imax < 0) return -1;
	
	return fSpecies[imax];
}

Bool_t AliLnID::GetITSlikelihood(Double_t pITS, Double_t dEdx, Int_t nPointsITS, Double_t* r) const
{
//
// Probability of dEdx in the ITS for each particle species.
// Adapted from STEER/ESD/AliESDpid.cxx
// (truncated mean method)
//
	if( pITS <= 0) return kFALSE;
	if( dEdx < 1.) return kFALSE;
	if( nPointsITS < 3) return kFALSE;
	
	Bool_t mismatch = kTRUE;
	
	Double_t sum = 0;
	for (Int_t i=0; i<kSPECIES; ++i)
	{
		Double_t mass = AliPID::ParticleMass(fSpecies[i]);
		Double_t p = (fSpecies[i] > AliPID::kTriton) ? 2.*pITS : pITS; // correct by Z
		Double_t expDedx = (fSpecies[i] > AliPID::kTriton) ? 4.*fITSpid->BetheAleph(p, mass) : fITSpid->BetheAleph(p, mass); // correct by Z^2
		Double_t sigma = fITSpid->GetResolution(expDedx, nPointsITS);
		
		if (TMath::Abs(dEdx - expDedx) > fRange*sigma)
		{
			r[i]=TMath::Exp(-0.5*fRange*fRange)/sigma;
		}
		else
		{
			r[i]=TMath::Exp(-0.5*(dEdx-expDedx)*(dEdx-expDedx)/(sigma*sigma))/sigma;
			mismatch = kFALSE;
		}
		
		sum += r[i];
	}
	
	if(sum <= 0. || mismatch)
	{
		return kFALSE;
	}
	
	for(Int_t i=0; i < kSPECIES; ++i) r[i] /= sum;
	
	return kTRUE;
}

Bool_t AliLnID::GetTPClikelihood(Double_t pTPC, Double_t dEdx, Int_t /*nPointsTPC*/, Double_t* r) const
{
//
// Probability of dEdx for each particle species in the TPC.
// Adapted from STEER/ESD/AliESDpid.cxx
//
	if( pTPC <= 0) return kFALSE;
	if( dEdx < 1.) return kFALSE;
	
	Bool_t mismatch = kTRUE;
	
	Double_t sum = 0;
	for(Int_t i=0; i < kSPECIES; ++i)
	{
		Double_t m = AliPID::ParticleMass(fSpecies[i]);
		Double_t p = (fSpecies[i] > AliPID::kTriton) ? 2.*pTPC : pTPC; // correct by Z
		Double_t expDedx = (fSpecies[i] > AliPID::kTriton) ? TMath::Power(2.,fZexp)*fTPCpid->Bethe(p/m) : fTPCpid->Bethe(p/m); // correct by Z^X (as in AliTPCPIDResponse)
		Double_t sigma = fTPCpid->GetRes0()*expDedx;
		
		if(TMath::Abs(dEdx - expDedx) > fRange*sigma)
		{
			r[i] = TMath::Exp(-0.5*fRange*fRange)/sigma;
		}
		else
		{
			r[i] = TMath::Exp(-0.5*(dEdx-expDedx)*(dEdx-expDedx)/(sigma*sigma))/sigma;
			mismatch=kFALSE;
		}
		
		sum += r[i];
	}
	
	if(sum <= 0. || mismatch)
	{
		return kFALSE;
	}
	
	for(Int_t i=0; i < kSPECIES; ++i) r[i] /= sum;
	
	return kTRUE;
}

Bool_t AliLnID::GetTOFlikelihood(Double_t pTOF, Double_t beta, Double_t* r) const
{
//
// Probability of beta for each particle species
//
	if( pTOF <= 0) return kFALSE;
	if( beta <= 0) return kFALSE;
	
	Double_t mismatch = kTRUE;
	
	Double_t sum = 0;
	for(Int_t i=0; i < kSPECIES; ++i)
	{
		Double_t mass = AliPID::ParticleMass(fSpecies[i]);
		Double_t p = (fSpecies[i] > AliPID::kTriton) ? 2.*pTOF : pTOF; // correct by Z
		Double_t expBeta = this->Beta(p,mass);
		Double_t sigma = this->GetBetaExpectedSigma(p,mass);
		
		if(TMath::Abs(beta-expBeta) > fRange*sigma)
		{
			r[i] = TMath::Exp(-0.5*fRange*fRange)/sigma;
		}
		else
		{
			r[i] = TMath::Exp(-0.5*(beta-expBeta)*(beta-expBeta)/(sigma*sigma))/sigma;
			mismatch = kFALSE;
		}
		
		sum += r[i];
	}
	
	if(sum <= 0. || mismatch)
	{
		return kFALSE;
	}
	
	for(Int_t i=0; i < kSPECIES; ++i) r[i] /= sum;
	
	return kTRUE;
}

Double_t AliLnID::Beta(Double_t p, Double_t m) const
{
//
// Expected beta for mass hypothesis m
//
	return p/TMath::Sqrt(p*p+m*m);
}

Double_t AliLnID::GetBetaExpectedSigma(Double_t p, Double_t m) const
{
//
// Expected sigma for the given mass hypothesis
//
	const Double_t kC0 = 0.0131203;
	const Double_t kC1 = 0.00670148;
	
	Double_t kappa = kC0*m*m/(p*(p*p+m*m)) + kC1;
	return kappa*Beta(p,m);
}

Bool_t AliLnID::GetITSmatch(Int_t pid, Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t nSigma) const
{
//
// Check if the signal is less than nSigma from the expected ITS dEdx
//
	Double_t mass = AliPID::ParticleMass(pid);
	Double_t p = (pid > AliPID::kTriton) ? 2.*pITS : pITS; // correct by Z
	Double_t expDedx = (pid > AliPID::kTriton) ? 4.*fITSpid->BetheAleph(p, mass) : fITSpid->BetheAleph(p, mass); // correct by Z^2
	Double_t sigma = fITSpid->GetResolution(expDedx, nPointsITS);
	
	if(TMath::Abs(dEdxITS-expDedx) < nSigma*sigma)
	{
		return kTRUE;
	}
	
	return kFALSE;
}

Bool_t AliLnID::GetTPCmatch(Int_t pid, Double_t pTPC, Double_t dEdxTPC, Double_t /*nPointsTPC*/, Double_t nSigma) const
{
//
// Check if the signal is less than nSigma from the expected TPC dEdx
//
	Double_t m = AliPID::ParticleMass(pid);
	Double_t p = (pid > AliPID::kTriton) ? 2.*pTPC : pTPC; // correct by Z
	Double_t expDedx = (pid > AliPID::kTriton) ? TMath::Power(2.,fZexp)*fTPCpid->Bethe(p/m) : fTPCpid->Bethe(p/m); // correct by Z^X (as in AliTPCPIDResponse)
	Double_t sigma = fTPCpid->GetRes0()*expDedx;
	
	if(TMath::Abs(dEdxTPC-expDedx) < nSigma*sigma)
	{
		return kTRUE;
	}
	
	return kFALSE;
}

Bool_t AliLnID::GetTOFmatch(Int_t pid, Double_t pTOF, Double_t beta, Double_t nSigma) const
{
//
// Check if the signal is less than nSigma from the expected velocity
//
	Double_t mass = AliPID::ParticleMass(pid);
	Double_t p = (pid > AliPID::kTriton) ? 2.*pTOF : pTOF;
	Double_t expBeta = this->Beta(p, mass);
	Double_t sigma = this->GetBetaExpectedSigma(p, mass);
	
	if(TMath::Abs(beta-expBeta) < nSigma*sigma)
	{
		return kTRUE;
	}
	
	return kFALSE;
}

Bool_t AliLnID::IsITSTPCmismatch(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t nSigma) const
{
//
// Check track TPC mismatch with ITS
//
	for(Int_t i=0; i<kSPECIES; ++i)
	{
		if(this->GetTPCmatch(fSpecies[i], pTPC, dEdxTPC, nPointsTPC, nSigma) &&
		   this->GetITSmatch(fSpecies[i], pITS, dEdxITS, nPointsITS, 5.))
		{
			return kFALSE;
		}
	}
	
	return kTRUE;
}

Bool_t AliLnID::IsITSTOFmismatch(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTOF, Double_t beta, Double_t nSigma) const
{
//
// Check track TOF mismatch with ITS
//
	for(Int_t i=0; i<kSPECIES; ++i)
	{
		if(this->GetTOFmatch(fSpecies[i], pTOF, beta, nSigma) &&
		   this->GetITSmatch(fSpecies[i], pITS, dEdxITS, nPointsITS, 3.))
		{
			return kFALSE;
		}
	}
	
	return kTRUE;
}

Bool_t AliLnID::IsTPCTOFmismatch(Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta, Double_t nSigma) const
{
//
// Check track TOF mismatch with TPC
//
	for(Int_t i=0; i<kSPECIES; ++i)
	{
		if(this->GetTOFmatch(fSpecies[i], pTOF, beta, nSigma) &&
		   this->GetTPCmatch(fSpecies[i], pTPC, dEdxTPC, nPointsTPC, 3.))
		{
			  return kFALSE;
		}
	}
	
	return kTRUE;
}
