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
	enum ScanSize {kNbinsA = 9, kNbinsB = 9};
public:
	AliPP13NonlinearityScanSelection(): AliPP13SpectrumSelectionMC() {}
	AliPP13NonlinearityScanSelection(const char * name, const char * title, AliPP13ClusterCuts cuts, 
			AliPP13SelectionWeightsScan * sw, Float_t precA = 0.01, Float_t precB = 0.1):
		AliPP13SpectrumSelectionMC(name, title, cuts, sw),
		fInvariantMass(),
		fMixInvariantMass(),
		fPtPrimaryPi0(),
		fPrecisionA(precA),
		fPrecisionB(precB)
	{

		Float_t nonA = sw->fE;
		Float_t nonB = sw->fD;

		for(Int_t ia = 0; ia < kNbinsA; ++ia)
		{
			for(Int_t ib = 0; ib < kNbinsB; ++ib)
			{
				AliPP13SelectionWeightsScan & swi = fWeightsScan[ia][ib];
				swi.fE = nonA + fPrecisionA * (kNbinsA / 2 - ia) / (kNbinsA / 2);
				swi.fD = nonB + fPrecisionB * (kNbinsB / 2 - ib) / (kNbinsB / 2);
			}
		}
	}

	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);
	virtual void InitSelectionHistograms();
protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	virtual TLorentzVector ClusterMomentumBinned(const AliVCluster * c1, const EventFlags & eflags, Int_t ia, Int_t ib) const;

	AliPP13NonlinearityScanSelection(const AliPP13NonlinearityScanSelection &);
	AliPP13NonlinearityScanSelection & operator = (const AliPP13NonlinearityScanSelection &);

private:

	// Set of weights that are need for NonlinearityScan
	AliPP13SelectionWeightsScan fWeightsScan[kNbinsA][kNbinsB];

	// Parameters of nonlinearity parametrization
	TH1 * fInvariantMass[kNbinsA][kNbinsB];    //!
	TH1 * fMixInvariantMass[kNbinsA][kNbinsB]; //!
	TH1 * fPtPrimaryPi0; //!

	Float_t fPrecisionA;
	Float_t fPrecisionB;

	ClassDef(AliPP13NonlinearityScanSelection, 2)
};
#endif
