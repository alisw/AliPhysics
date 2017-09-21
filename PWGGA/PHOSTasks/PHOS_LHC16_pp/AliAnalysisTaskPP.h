#ifndef ALIANALYSISTASKPROMPT_H
#define ALIANALYSISTASKPROMPT_H

// --- ROOT system ---
#include <TString.h>
#include <TClonesArray.h>

// --- AliRoot header files ---
#include <AliAnalysisTaskSE.h>

// --- Custom libraries ---
#include "PhotonSelection.h"
#include "MixingSample.h"

class AliAnalysisTaskPP : public AliAnalysisTaskSE
{

public:
	enum {kMinModule = 1, kMaxModule=4};
	AliAnalysisTaskPP();
	AliAnalysisTaskPP(const char * name, TList * selections, Int_t nmix = 100);
	virtual ~AliAnalysisTaskPP();

	void   UserCreateOutputObjects();
	void   UserExec(Option_t *);
	void   Terminate(Option_t *);

	void SetBadCells(Int_t badcells[], Int_t nbad);
	void SetBadMap(const char * filename);
	TList * GetSelections() { return fSelections; }

protected:
	TClonesArray * GetMCParticles(const AliVEvent * event) const;
	Bool_t EventSelected(const AliVEvent * event, EventFlags & eprops) const;
	Bool_t IsClusterBad(AliVCluster * clus) const;
	Bool_t CellInPhos(Int_t absId, Int_t & sm, Int_t & ix, Int_t & iz) const;

	AliAnalysisTaskPP & operator =(const AliAnalysisTaskPP &);
	AliAnalysisTaskPP(const AliAnalysisTaskPP & c);

private:
	MixingSample * fPreviousEvents;
	TList * fSelections;     // analysis instance
	TH2I  * fPHOSBadMap[kMaxModule - kMinModule + 1];  // bad channel maps
	Int_t fNMixedEvents;     // number of events used for mixing
	Int_t fNBad;             // number of entries in fBadCells
	Int_t * fBadCells;       //[fNBad] bad cells array

	ClassDef(AliAnalysisTaskPP, 2);
};

#endif