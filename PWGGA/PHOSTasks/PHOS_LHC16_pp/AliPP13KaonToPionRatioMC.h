#ifndef ALIPP13KAONTOPIONRATIOMC_H
#define ALIPP13KAONTOPIONRATIOMC_H


#include <map>

// --- Custom header files ---
#include "AliPP13SpectrumSelectionMC.h"
#include "AliPP13ParticlesHistogram.h"
#include "AliPP13SelectionWeights.h"

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


// NB: This will simplify the code
//

// TODO: Split this class to have separate efficiency and contamination estimators
//

class AliPP13KaonToPionRatioMC: public AliPP13SpectrumSelectionMC
{
public:
	enum Modes {kGenerated = 0, kReconstructed = 1, kNhists = 2};
	enum Particles
	{
		kPiplus = 211,
		kPiminus = -211,
		kKplus = 321,
		kKminus = -321
	};


	AliPP13KaonToPionRatioMC():
		AliPP13SpectrumSelectionMC(),
		fPrimary(),
		fAll()
	{
		fPartNames[kPiplus] = "#pi^{+}";
		fPartNames[kPiminus] = "#pi^{-}";
		fPartNames[kKplus] = "K^{+}";
		fPartNames[kKminus] = "K^{-}";
	}

	AliPP13KaonToPionRatioMC(const char * name, const char * title,
			AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13SpectrumSelectionMC(name, title, cuts, w),
		fPrimary(),
		fAll()
	{
		// Force no timing cut for MC,
		// as there is no photons from different bunches
		fCuts.fTimingCut = 9999;
		fPartNames[kPiplus] = "#pi^{+}";
		fPartNames[kPiminus] = "#pi^{-}";
		fPartNames[kKplus] = "K^{+}";
		fPartNames[kKminus] = "K^{-}";
	}



	virtual void InitSelectionHistograms();
	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);

protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
	{
		(void) c1;
		(void) c2;
		(void) eflags;
	}

	AliPP13KaonToPionRatioMC(const AliPP13KaonToPionRatioMC &);
	AliPP13KaonToPionRatioMC & operator = (const AliPP13KaonToPionRatioMC &);

	typedef std::map<Int_t, TString> EnumNames;
	EnumNames fPartNames;

	typedef std::map<Int_t, TH1 *> EnumHists;
	EnumHists fPrimary;     //!
	EnumHists fAll;          //!

	// Parameters of weighed MC parametrization
	ClassDef(AliPP13KaonToPionRatioMC, 2)
};
#endif
