#ifndef ALIPP13PYTHIAINFOSELECTION_H
#define ALIPP13PYTHIAINFOSELECTION_H

// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13SelectionWeights.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliAODMCParticle.h>
#include <AliVCluster.h>
#include <AliStack.h>
#include <AliLog.h>

class AliPP13PythiaInfoSelection: public AliPP13PhysicsSelection
{
public:
	AliPP13PythiaInfoSelection():
		AliPP13PhysicsSelection(),
		fXsec(0),
		fTrials(0)
	{

	}

	AliPP13PythiaInfoSelection(const char * name, const char * title):
		// NB: We don't need to set a pointer here 
		//
		AliPP13PhysicsSelection(name, title, AliPP13ClusterCuts(), 0), 
		fXsec(0),
		fTrials(0)
	{

	}

	virtual void InitSelectionHistograms();

	// Fetch all pythia info here
	virtual void CountMBEvent();

	// Make these methods empty
	virtual void FillHistograms(TObjArray * clusArray, TList * pool, const EventFlags & eflags);

protected:
	AliPP13PythiaInfoSelection(const AliPP13PythiaInfoSelection &);
	AliPP13PythiaInfoSelection & operator = (const AliPP13PythiaInfoSelection &);

private:
	// NB: Don't delete these pointers
	TH1F * fXsec;    //!
	TH1F * fTrials;  //!

	ClassDef(AliPP13PythiaInfoSelection, 2)
};
#endif
