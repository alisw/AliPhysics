#ifndef ALIPP13EFFICIENCYSELECTIONMC_H
#define ALIPP13EFFICIENCYSELECTIONMC_H


#include <map>

// --- Custom header files ---
#include "AliPP13SpectrumSelectionMC.h"
#include "AliPP13SelectionWeights.h"
#include "AliPP13MesonSelectionMC.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>
#include <TF1.h>

// --- AliRoot header files ---
#include <AliAODMCParticle.h>
#include <AliVCluster.h>
#include <AliStack.h>
#include <AliLog.h>


class AliPP13EfficiencySelectionMC: public AliPP13SpectrumSelectionMC
{
public:
	AliPP13EfficiencySelectionMC():
		AliPP13SpectrumSelectionMC(),
		fInvMass()
	{
		fPartNames[kGamma] = "#gamma";
		fPartNames[kPi0] = "#pi^{0}";
		fPartNames[kEta] = "#eta";
	}

	AliPP13EfficiencySelectionMC(const char * name, const char * title, AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13SpectrumSelectionMC(name, title, cuts, w),
		fInvMass()
	{
		// Force no timing cut for MC,
		// as there is no photons from different bunches
		fCuts.fTimingCut = 9999;

		// Don't use c++11 here, as it might fail at some nodes
		fPartNames[kGamma] = "#gamma";
		fPartNames[kPi0] = "#pi^{0}";
		fPartNames[kEta] = "#eta";

	}

	virtual void InitSelectionHistograms();
	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);

	virtual ~AliPP13EfficiencySelectionMC()
	{
		for (ParticleSpectrums::iterator i = fSpectrums.begin(); i != fSpectrums.end(); ++i)
			delete i->second;
	}

protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	AliPP13EfficiencySelectionMC(const AliPP13EfficiencySelectionMC &);
	AliPP13EfficiencySelectionMC & operator = (const AliPP13EfficiencySelectionMC &);
	// NB: This data structure contains all necesary histograms
	//     for the particles we want to get
	typedef std::map<Int_t, ParticleSpectrum * > ParticleSpectrums;
	ParticleSpectrums fSpectrums;

	typedef std::map<Int_t, TString> EnumNames;
	EnumNames fPartNames;

	TH1 * fInvMass[2]; //!
	ClassDef(AliPP13EfficiencySelectionMC, 2)
};
#endif
