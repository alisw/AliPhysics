#ifndef NONLINEARITYSCANSELECTION_H
#define NONLINEARITYSCANSELECTION_H

// --- Custom header files ---
#include "PhysPhotonSelection.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

class NonlinearityScanSelection : public PhysPhotonSelection
{
	enum ScanSize {kNbinsA = 11, kNbinsSigma = 11};
public:
	NonlinearityScanSelection(): PhysPhotonSelection() {}
	NonlinearityScanSelection(const char * name, const char * title, ClusterCuts cuts,
	                           Float_t nona = 0., Float_t nonsigma = 1., Float_t genergy = 1.):
		PhysPhotonSelection(name, title, cuts),
		fInvariantMass(),
		fMixInvariantMass(),
		fNonA(nona),
		fNonSigma(nonsigma),
		fGlobalEnergyScale(genergy),
		fPrecisionA(0.01),
		fPrecisionSigma(0.1)
	{
		fCuts.fTimingCut = 99999; // No timing cut in MC
	}

	virtual void InitSelectionHistograms();
protected:

	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags, Int_t ia = -1, Int_t ib = -1) const;
	virtual Float_t Nonlinearity(Float_t x, Int_t ia = -1, Int_t ib = -1) const;

	virtual Float_t GetA(Int_t ia) const
	{
		if(ia < 0)
			return fNonA;

		return fNonA - fPrecisionA * kNbinsA / 2 + ia * fPrecisionA;
	}

	virtual Float_t GetSigma(Int_t ib) const
	{
		if(ib < 0)
			return fNonSigma;

		return fNonSigma - fPrecisionSigma * kNbinsSigma / 2 + ib * fPrecisionSigma;
	}

	NonlinearityScanSelection(const NonlinearityScanSelection &);
	NonlinearityScanSelection & operator = (const NonlinearityScanSelection &);

private:

	// Parameters of nonlinearity parametrization
	TH1 * fInvariantMass[kNbinsA][kNbinsSigma];    //!
	TH1 * fMixInvariantMass[kNbinsA][kNbinsSigma]; //!

	Float_t fNonA;
	Float_t fNonSigma;
	Float_t fGlobalEnergyScale;

	Float_t fPrecisionA;
	Float_t fPrecisionSigma;

	ClassDef(NonlinearityScanSelection, 2)
};
#endif