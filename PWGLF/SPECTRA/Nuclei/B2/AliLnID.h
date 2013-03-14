#ifndef ALILNID_H
#define ALILNID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// class for light nuclei identification
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <AliPID.h>

class TParticle;

class AliLnID: public TObject
{
  public:
	AliLnID();
	AliLnID(const AliLnID& other);
	AliLnID& operator=(const AliLnID& other);
	
	virtual ~AliLnID();
	
	Int_t GetPID(const TParticle* p) const;
	
	Int_t GetPID(Int_t partCode, Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta, Double_t nSigITS=3, Double_t nSigTPC=3, Double_t nSigTOF=3) const;
	
	Int_t GetTPCpid(Int_t partCode, Double_t pTPC, Double_t dEdx, Double_t nPoints, Double_t nSigma=3) const;
	
	Int_t GetITSTPCpid(Int_t partCode, Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t nSigmaITS, Double_t pTPC, Double_t dEdxTPC, Double_t nPointsTPC, Double_t nSigmaTPC) const;
	
	Int_t GetTPCTOFpid(Int_t partCode, Double_t pTPC, Double_t dEdx, Double_t nPoints, Double_t nSigmaTPC, Double_t pTOF, Double_t beta, Double_t nSigmaTOF) const;
	
	Int_t GetMaxLikelihoodPID(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta) const;
	
	Int_t GetBayesPID(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta) const;
	
	Bool_t GetITSlikelihood(Double_t p, Double_t dEdx, Int_t nPoints, Double_t* r) const;
	Bool_t GetTPClikelihood(Double_t p, Double_t dEdx, Int_t nPoints, Double_t* r) const;
	Bool_t GetTOFlikelihood(Double_t p, Double_t beta, Double_t* r) const;
	
	Double_t GetBetaExpectedSigma(Double_t p, Double_t mass) const;
	
	Bool_t GetITSmatch(Int_t pid, Double_t p, Double_t dEdx, Int_t nPoints, Double_t nSigma=3.) const;
	Bool_t GetTPCmatch(Int_t pid, Double_t pTPC, Double_t dEdx, Double_t nPoints, Double_t nSigma=3.) const;
	Bool_t GetTOFmatch(Int_t pid, Double_t pTOF, Double_t beta, Double_t nSigma=3.) const;
	
	Int_t GetPidProcedure() const { return fPidProcedure; }
	
	Bool_t IsITSTPCmismatch(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t nSigma=5.) const;
	Bool_t IsITSTOFmismatch(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTOF, Double_t beta, Double_t nSigma=5.) const;
	Bool_t IsTPCTOFmismatch(Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta, Double_t nSigma=5.) const;
	
	void SetPriorProbabilities(const Double_t* prob);
	void SetPriorProbabilities(Double_t e, Double_t mu, Double_t pi, Double_t k, Double_t p, Double_t d, Double_t t, Double_t he3, Double_t alpha);
	
	void SetPidProcedure(Int_t proc) { fPidProcedure = proc; }
	
	void SetITSBetheBlochParams(const Double_t* param, Double_t res=0.13);
	
	void SetTPCBetheBlochParams(const Double_t* param, Double_t mip=1., Double_t res=0.06);
	void SetTPCBetheBlochParams(Double_t c0, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t mip=1., Double_t res=0.06);
	
	void SetTPCChargeCorrection(Double_t zexp) { fZexp = zexp; }
	
	enum { kSPECIES = 9 };
	enum { kBayes, kMaxLikelihood, kTPC, kITSTPC, kTPCTOF };
	
  private:
 
	Double_t Beta(Double_t p, Double_t m) const;
	Bool_t GetLikelihood(Double_t pITS, Double_t dEdxITS, Int_t nPointsITS, Double_t pTPC, Double_t dEdxTPC, Int_t nPointsTPC, Double_t pTOF, Double_t beta, Double_t* r) const;
	Int_t GetIndexOfMaxValue(const Double_t* w) const;
	
  private:
 
	Int_t fPidProcedure; // PID procedure code
	AliPID::EParticleType fSpecies[kSPECIES]; // particle species known by the pid
	
	Double_t fPrior[kSPECIES]; // prior probabilities
	Double_t fRange; // number of sigmas to the expected values
	Double_t fZexp; // TPC BB charge dependence
	
	class AliITSPIDResponse* fITSpid; // ITS likelihood
	class AliTPCPIDResponse* fTPCpid; // TPC likelihood
	
	ClassDef(AliLnID, 1)
};

#endif // ALILNID_H
