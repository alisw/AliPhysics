#ifndef ALIPP13NONLINEARITYSCANSELECTION_H
#define ALIPP13NONLINEARITYSCANSELECTION_H

// --- Custom header files ---
#include "AliPP13PhysPhotonSelection.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

class AliPP13NonlinearityScanSelection : public AliPP13PhysPhotonSelection
{
	enum ScanSize {kNbinsA = 11, kNbinsSigma = 11};
public:
	AliPP13NonlinearityScanSelection(): AliPP13PhysPhotonSelection() {}
	AliPP13NonlinearityScanSelection(const char * name, const char * title, AliPP13ClusterCuts cuts,
	                           Float_t nona = 0., Float_t nonsigma = 1., Float_t genergy = 1.):
		AliPP13PhysPhotonSelection(name, title, cuts),
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

	AliPP13NonlinearityScanSelection(const AliPP13NonlinearityScanSelection &);
	AliPP13NonlinearityScanSelection & operator = (const AliPP13NonlinearityScanSelection &);

private:

	// Parameters of nonlinearity parametrization
	TH1 * fInvariantMass[kNbinsA][kNbinsSigma];    //!
	TH1 * fMixInvariantMass[kNbinsA][kNbinsSigma]; //!

	Float_t fNonA;
	Float_t fNonSigma;
	Float_t fGlobalEnergyScale;

	Float_t fPrecisionA;
	Float_t fPrecisionSigma;

	ClassDef(AliPP13NonlinearityScanSelection, 2)
};
#endif