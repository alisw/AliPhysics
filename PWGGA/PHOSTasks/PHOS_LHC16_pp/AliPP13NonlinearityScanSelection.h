#ifndef ALIPP13NONLINEARITYSCANSELECTION_H
#define ALIPP13NONLINEARITYSCANSELECTION_H

// --- Custom header files ---
#include "AliPP13SelectionWeights.h"
#include "AliPP13SpectrumSelectionMC.h"

// --- AliRoot header files ---
#include <AliVCluster.h>

// TODO: Fix logic for pointers
//

// NB: Don't use pointers in the array, this way will be easier
//
class AliPP13NonlinearityScanSelection : public AliPP13SpectrumSelectionMC
{
	enum ScanSize {kNbinsA = 9, kNbinsSigma = 9};
public:
	AliPP13NonlinearityScanSelection(): AliPP13SpectrumSelectionMC() {}
	AliPP13NonlinearityScanSelection(const char * name, const char * title, AliPP13ClusterCuts cuts, 
			AliPP13SelectionWeightsMC * sw, Float_t precA = 0.01, Float_t precSigma = 0.1):
		AliPP13SpectrumSelectionMC(name, title, cuts, sw),
		fInvariantMass(),
		fMixInvariantMass(),
		fPtPrimaryPi0(),
		fPrecisionA(precA),
		fPrecisionSigma(precSigma)
	{

		Float_t nona = sw->fNonA;
		Float_t nonSigma = sw->fNonSigma;

		for(Int_t ia = 0; ia < kNbinsA; ++ia)
		{
			for(Int_t ib = 0; ib < kNbinsSigma; ++ib)
			{
				AliPP13SelectionWeightsMC & swi = fWeightsScan[ia][ib];
				swi.fNonGlobal = sw->fNonGlobal;
				swi.fNonA = nona - fPrecisionA * kNbinsA / 2 + ia * fPrecisionA;
				swi.fNonSigma = nonSigma - fPrecisionSigma * kNbinsSigma / 2 + ib * fPrecisionSigma;
			}
		}
	}

	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);
	virtual void InitSelectionHistograms();
protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags, Int_t ia, Int_t ib) const;

	AliPP13NonlinearityScanSelection(const AliPP13NonlinearityScanSelection &);
	AliPP13NonlinearityScanSelection & operator = (const AliPP13NonlinearityScanSelection &);

private:

	// Set of weights that are need for NonlinearityScan
	AliPP13SelectionWeightsMC fWeightsScan[kNbinsA][kNbinsSigma];

	// Parameters of nonlinearity parametrization
	TH1 * fInvariantMass[kNbinsA][kNbinsSigma];    //!
	TH1 * fMixInvariantMass[kNbinsA][kNbinsSigma]; //!
	TH1 * fPtPrimaryPi0; //!

	Float_t fPrecisionA;
	Float_t fPrecisionSigma;

	ClassDef(AliPP13NonlinearityScanSelection, 2)
};
#endif
