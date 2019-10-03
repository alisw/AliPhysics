#ifndef AliLMRMuon_H
#define AliLMRMuon_H

#include "TObject.h"
#include "TMath.h"
#include "TLorentzVector.h"

class AliLMRMuon : public TObject {
public:
	AliLMRMuon();
	AliLMRMuon(const AliLMRMuon& part); // copy constructor
	virtual ~AliLMRMuon(){};

	void SetTriggerMatch(Short_t  triggerMatch) {fTriggerMatch = triggerMatch ;};
	void SetpDCA(Double_t pdca) {fpDCA = pdca ;};
	void SetRabs(Double_t rAbs) {fRabs =  rAbs;};
	void SetMomentum(Double_t x, Double_t y, Double_t z) {fPx=x;fPy=y;fPz=z;};
	void SetChi2(Double_t  chi2) {fChi2 = chi2 ;};
	void SetChi2Match(Double_t chi2Match) {fChi2Match = chi2Match;};
	void SetCharge(Short_t charge) {fCharge=charge;};
	void SetSelectionMask(UInt_t mask) {fSelectionMask=mask;};
	void SetLocalBoard(UShort_t LB) {fLocalBoard=LB;};

	UInt_t GetSekectionMask(){return fSelectionMask;};
	Short_t GetTriggerMatch() {return fTriggerMatch;};
	Double_t GetpDCA() {return fpDCA ;};
	Double_t GetRabs() {return fRabs ;};
	Double_t GetChi2() {return fChi2 ;};
	Double_t GetChi2Match() {return fChi2Match ;};
	Double_t P() const {return TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz) ;};
	Double_t Px() const {return fPx ;};
	Double_t Py() const {return fPy ;};
	Double_t Pz() const {return fPz ;};
	Double_t Energy() const {return (TMath::Sqrt(P()*P() + kMuonMass*kMuonMass)) ;};
	Short_t Charge() const {return fCharge;};
	TLorentzVector P4() const;
	UShort_t GetLocalBoard() {return fLocalBoard;};

private:
	Short_t fTriggerMatch;
	UShort_t fLocalBoard;
	UInt_t fSelectionMask;	
	Double_t fpDCA;
	Double_t fRabs;
	Short_t fCharge;
	Double_t fPx;
	Double_t fPy;
	Double_t fPz;
	Double_t fChi2;
	Double_t fChi2Match;
	static const Double_t kMuonMass;

	ClassDef(AliLMRMuon, 2)  //The class title
};

#endif

