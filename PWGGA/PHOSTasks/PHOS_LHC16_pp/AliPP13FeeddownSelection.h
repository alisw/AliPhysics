#ifndef ALIPP13FEEDDOWNSELECTION_H
#define ALIPP13FEEDDOWNSELECTION_H


// --- Custom header files ---
#include "AliPP13SpectrumSelectionMC.h"
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



class AliPP13FeeddownSelection: public AliPP13SpectrumSelectionMC
{
public:
	enum Particles
	{
		kGamma = 22, kPi0 = 111, kEta = 221, kK0s = 310,
		kOmega = 223, kLambda = 3122, kPPion = 211, kNPion = -211,
		kPRho = 213, kNRho = -213,
		kKStarP = 323, kKStarN = -323, kKStar0 = 313, kBarKstar0 = -313,
		kKplus = 321, kKminus = -321, kSigmaZero = 3212
	};

	AliPP13FeeddownSelection():
		AliPP13SpectrumSelectionMC(),
		fInvMass(),
		fFeedownK0s()
	{
	}

	AliPP13FeeddownSelection(const char * name, const char * title,
	                        AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13SpectrumSelectionMC(name, title, cuts, w),
		fInvMass(),
		fFeedownK0s()
	{
	}

	virtual void InitSelectionHistograms();
	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);

protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);
	virtual AliAODMCParticle * GetMother(const AliAODMCParticle * particle, TClonesArray * particles) const;
	virtual AliAODMCParticle * GetMother(const AliVCluster * c1, TClonesArray * particles) const
	{
		Int_t label = c1->GetLabelAt(0);
		if (label <= -1)
			return 0;
		AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle * >(particles->At(label));
		return GetMother(particle, particles);
	}
	AliPP13FeeddownSelection(const AliPP13FeeddownSelection &);
	AliPP13FeeddownSelection & operator = (const AliPP13FeeddownSelection &);

	TH1 * fInvMass[2];     //!
	TH1 * fFeedownK0s[2];  //!

	// Parameters of weighed MC parametrization
	ClassDef(AliPP13FeeddownSelection, 2)
};
#endif
