#ifndef ALIAODEVENTCUTSDIHADRONPID_H 
#define ALIAODEVENTCUTSDIHADRONPID_H

#include "TString.h"
#include "TH1F.h"
#include "TList.h"

class AliAODEventCutsDiHadronPID : public TNamed 

{

public:
	AliAODEventCutsDiHadronPID();								// Default Constructor
	AliAODEventCutsDiHadronPID(const char* name);				// Named Constructor
	AliAODEventCutsDiHadronPID(const AliAODEventCutsDiHadronPID &source);	// Copy Constructor
	virtual ~AliAODEventCutsDiHadronPID();						// Destructor
	virtual Long64_t Merge(TCollection* list);					// Merger

	void CreateHistos();					// Create QA histograms

public:
	Bool_t IsSelected(AliAODEvent* event);
	void PrintCuts();

// Setters
	void SetIsPbPb(Bool_t ispbpb = kTRUE) {fIsPbPb = ispbpb;}
	void SetIsMC(Bool_t ismc = kTRUE) {fIsMC = ismc;}

	void SetTrigger(UInt_t trigger) {
		fTrigger = trigger;
		fTestTrigger = kTRUE;
	}

	// Note that minCentrality is expected to be the biggest number.
	void SetCentrality(Float_t maxCentrality, Float_t minCentrality) {
		if (minCentrality > maxCentrality) {
			fMinCentrality = minCentrality;
			fMaxCentrality = maxCentrality;
		} else {
			fMinCentrality = maxCentrality;
			fMaxCentrality = minCentrality;
		}
		fTestCentrality = kTRUE;
	} 
	
	void SetCentralityEstimator(const char* centralityestimator) {
		fCentralityEstimator = centralityestimator;
	}
	void SetMaxVertexZ(Float_t maxVertexZ) {
		fMaxVertexZ = maxVertexZ;
		fTestVertexZ = kTRUE;
	}
	void SetMinReferenceMultiplicity(Int_t minrefmult) {
		fMinRefMult = minrefmult;
		fTestMinRefMult = kTRUE;
	}

	void SetDebugLevel(Int_t debuglevel) {fDebug = debuglevel;}

// Getters
	Bool_t GetIsPbPb() const {return fIsPbPb;} 
	Bool_t GetIsMC() const {return fIsMC;}
	UInt_t GetTrigger() const {return fTrigger;}
	Float_t GetMinCentrality() const {return fMinCentrality;}
	Float_t GetMaxCentrality() const {return fMaxCentrality;}
	TString GetCentralityEstimator() const {return fCentralityEstimator;}
	Float_t GetMaxVertexZ() const {return fMaxVertexZ;}

	// Functions returning pointer data members aren't very safe.
	TList* GetListOfSelectedEventQAHistos() {
		if (fSelectedEventQAHistos) return fSelectedEventQAHistos;
		else return 0x0;
	}
	TList* GetListOfAllEventQAHistos() {
		if (fAllEventQAHistos) return fAllEventQAHistos;
		else return 0x0;
	}

	TObject* GetHistSelectedEvents(const char* name) {return fSelectedEventQAHistos->FindObject(name);}
	TObject* GetHistAllEvents(const char* name) {return fAllEventQAHistos->FindObject(name);}

	// Cannot be made const because GetHistSelectedEvents() isn't safe.
	Int_t GetNAcceptedEvents() {return ((TH1F*)GetHistSelectedEvents("fHistTriggerSelected"))->GetEntries();}

	Int_t GetDebugLevel() const {return fDebug;}

private:
// Expected Event Details.
	Bool_t fIsPbPb;
	Bool_t fIsMC;

// Event Cuts.
	UInt_t fTrigger;
	Float_t fMinCentrality;
	Float_t fMaxCentrality;
	TString fCentralityEstimator;
	Float_t fMaxVertexZ;
	Int_t fMinRefMult;

// Which cuts to be checked.
	Bool_t fTestTrigger;
	Bool_t fTestCentrality;
	Bool_t fTestVertexZ;
	Bool_t fTestMinRefMult;

// QA histograms (don't stream)
	TList* fSelectedEventQAHistos;
	TList* fAllEventQAHistos;
	TH1F** fHistTrigger;						//! Trigger
	TH1F** fHistRefMultiplicity;				//! Number of tracks
	TH1F** fHistCentrality;						//! Centrality
	TH1F** fHistCentralityQuality;				//! Centrality Quality
	TH1F** fHistVertexZ;						//! VertexZ

	Int_t fDebug;								// Debug flag.

	ClassDef(AliAODEventCutsDiHadronPID,2);

};

#endif
