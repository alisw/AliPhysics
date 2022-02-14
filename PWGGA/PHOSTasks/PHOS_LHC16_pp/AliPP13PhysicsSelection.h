#ifndef ALIPP13PHYSICSSELECTION_H
#define ALIPP13PHYSICSSELECTION_H

// --- Custom libraries ---
#include "AliPP13ClusterCuts.h"
#include "AliPP13SelectionWeights.h"


// --- ROOT system ---
#include <TObjArray.h>
#include "TNamed.h"
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

// --- AliRoot header files ---
#include <AliPIDResponse.h>
#include <AliPHOSGeometry.h>
#include <AliVCluster.h>
#include <AliLog.h>


struct SelectionLimits
{
	SelectionLimits()
	{
		nM = 750;
		mMin = 0.0;
		mMax = 1.5;
		nPt = 800;
		ptMin = 0;
		ptMax = 40;

	}

	Int_t nM;
	Double_t mMin;
	Double_t mMax;
	Int_t nPt;
	Double_t ptMin;
	Double_t ptMax;
};


class AliPP13PhysicsSelection : public TNamed
{
public:

	AliPP13PhysicsSelection():
		TNamed(),
		fListOfHistos(0),
		fCuts(),
		fWeights(),
		fEventCounter(0),
		fLimits()
	{}

	AliPP13PhysicsSelection(const char * name, const char * title, AliPP13ClusterCuts cuts,
	                        AliPP13SelectionWeights * sw):
		TNamed(name, title),
		fListOfHistos(0),
		fCuts(cuts),
		fWeights(dynamic_cast<AliPP13SelectionWeights *>(sw->Clone())),
		fEventCounter(0),
		fLimits()
	{}

	virtual ~AliPP13PhysicsSelection();

	virtual void InitSummaryHistograms();
	virtual void InitSelectionHistograms() = 0;
	virtual Bool_t SelectEvent(const EventFlags & flgs);

	// This is a dummy method to count number of Triggered Events.
	virtual void CountMBEvent();
	virtual void FillHistograms(TObjArray * clusArray, TList * pool, const EventFlags & eflags); // implements algorithm
	virtual void ConsiderGeneratedParticles(const EventFlags & eflags)
	{
		(void) eflags;
	}

	virtual void MixPhotons(TObjArray & photons, TList * pool, const EventFlags & eflags);
	virtual TList * GetListOfHistos() { return fListOfHistos; }

protected:
	virtual Bool_t AcceptPhotonCandidate(AliVCluster * clus, const EventFlags & eflags);
	virtual void SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags);
	virtual void SelectTwoParticleCombinations(const TObjArray & photonCandidates, const EventFlags & flags);
	virtual void FillClusterHistograms(const AliVCluster * c, const EventFlags & eflags)
	{
		(void) c;
		(void) eflags;
	}

	virtual Int_t CheckClusterGetSM(const AliVCluster * clus, Int_t & x, Int_t & z) const;
	virtual TLorentzVector ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const;

	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
	{
		(void) c1;
		(void) c2;
		(void) eflags;
	}

	AliPP13PhysicsSelection(const AliPP13PhysicsSelection &);
	AliPP13PhysicsSelection & operator = (const AliPP13PhysicsSelection &);

	TList  * fListOfHistos;  //! list of histograms
	AliPP13ClusterCuts fCuts;

	AliPP13SelectionWeights * fWeights;
	TH1 * fEventCounter;  //!

	SelectionLimits fLimits;
private:
	ClassDef(AliPP13PhysicsSelection, 2)
};
#endif
