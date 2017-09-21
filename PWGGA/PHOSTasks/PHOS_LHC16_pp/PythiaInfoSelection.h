#ifndef PYTHIAINFOSELECTION_H
#define PYTHIAINFOSELECTION_H

// --- Custom header files ---
#include "PhotonSelection.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliAODMCParticle.h>
#include <AliVCluster.h>
#include <AliStack.h>
#include <AliLog.h>

class PythiaInfoSelection: public PhotonSelection
{
public:
	PythiaInfoSelection():
		PhotonSelection(),
		fXsec(0),
		fTrials(0)
	{

	}

	PythiaInfoSelection(const char * name, const char * title):
		PhotonSelection(name, title, ClusterCuts()),
		fXsec(0),
		fTrials(0)
	{

	}

	virtual void InitSelectionHistograms();

	// Fetch all pythia info here
	virtual void CountMBEvent();

	// Make these methods empty
	virtual void FillPi0Mass(TObjArray * clusArray, TList * pool, const EventFlags & eflags);

protected:
	PythiaInfoSelection(const PythiaInfoSelection &);
	PythiaInfoSelection & operator = (const PythiaInfoSelection &);

private:
	// NB: Don't delete these pointers
	TH1F * fXsec;    //!
	TH1F * fTrials;  //!

	ClassDef(PythiaInfoSelection, 2)
};
#endif