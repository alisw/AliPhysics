#ifndef AliLMREvent_H
#define AliLMREvent_H


#include "TObject.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"

class AliLMRMuon;

class AliLMREvent : public TObject {
public:
	AliLMREvent();
	AliLMREvent(Int_t run, Double_t evtPlane, Double_t Vert[3],Double_t *Activity);
	AliLMREvent (const AliLMREvent& evt); // copy constructor
	virtual ~AliLMREvent();
	AliLMREvent& operator=(const AliLMREvent&);
	void Clear(Option_t*);

	void SetEventPlane(Double_t eventPlane) {fEventPlane = eventPlane;};
	void SetRunNumber(Int_t RunNb) {fRunNb = RunNb;};
	void SetXVertex(Double_t x) {fVertex[0] = x;};
	void SetYVertex(Double_t y) {fVertex[1] = y;};
	void SetZVertex(Double_t z) {fVertex[2] = z;};
	void SetVertex(Double_t V[3]);
	void SetVtxContributors(Int_t nContributors){fNVtxContributors=nContributors;};
	void SetTriggerString(TString val) { fTriggerString = val;};
	void SetMultiplicity(TString method,Double_t multiplicity);
	void SetPhysicsSelectionMask(UShort_t mask){fPhysicsSelectionMask=mask;};
	Double_t GetEventPlane() {return fEventPlane ;};
	Int_t GetRunNumber() {return fRunNb ;};
	Double_t GetXVertex() {return fVertex[0] ;};
	Double_t GetYVertex() {return fVertex[1] ;};
	Double_t GetZVertex() {return fVertex[2] ;};
	Double_t*GetVertex()  {return fVertex    ;};
	TString GetTriggerString() {return fTriggerString ;};
	Int_t GetNVtxContributors(){return fNVtxContributors;};
	UShort_t GetPhysicsSelectionMask(){return fPhysicsSelectionMask;};
	Double_t GetMultiplicity(TString method);
	
	AliLMRMuon *AddMuon();
	AliLMRMuon *GetMuon(Int_t imu) { return ((imu < fMuons->GetEntriesFast()) ? (AliLMRMuon *)fMuons->At(imu) : 0);};
	Int_t GetNMuons() const {return (fMuons ? fMuons->GetEntriesFast() : -1);}
	
private:
	TClonesArray * fMuons;
	Double_t fEventPlane;
	Double_t fMultiplicity_V0M;
	Double_t fMultiplicity_ADM;
	Double_t fMultiplicity_SPDTracklets;
	Double_t fMultiplicity_SPDClusters;
	Double_t fMultiplicity_RefMult05;
	Double_t fMultiplicity_RefMult08;
	Int_t fRunNb;
	Short_t fNMuons;
	UShort_t fPhysicsSelectionMask;
	Double_t fVertex[3];
	Int_t fNVtxContributors;
	TString fTriggerString;

	ClassDef(AliLMREvent, 2)  //The class title

};

#endif

