#ifndef ALIANALYSISTASKPP13_H
#define ALIANALYSISTASKPP13_H

// --- ROOT system ---
#include <TString.h>
#include <TClonesArray.h>

// --- AliRoot header files ---
#include <AliAnalysisTaskSE.h>

// --- Custom libraries ---
#include "AliPP13PhysicsSelection.h"
#include "AliPP13MixingSample.h"

class AliAnalysisTaskPP13 : public AliAnalysisTaskSE
{

public:
	enum {kMinModule = 1, kMaxModule=4};
	AliAnalysisTaskPP13();
	AliAnalysisTaskPP13(const char * name, TList * selections, Int_t nmix = 100);
	virtual ~AliAnalysisTaskPP13();

	void   UserCreateOutputObjects();
	void   UserExec(Option_t *);
	void   Terminate(Option_t *);

	TList * GetSelections() { return fSelections; }

protected:
	TClonesArray * GetMCParticles(const AliVEvent * event) const;
	Bool_t EventSelected(const AliVEvent * event, EventFlags & eprops) const;

	AliAnalysisTaskPP13 & operator =(const AliAnalysisTaskPP13 &);
	AliAnalysisTaskPP13(const AliAnalysisTaskPP13 & c);

private:
	AliPP13MixingSample * fPreviousEvents;
	TList * fSelections;     // analysis instance
	Int_t fNMixedEvents;     // number of events used for mixing

	ClassDef(AliAnalysisTaskPP13, 2);
};

#endif
