// --- Standard libraries ---
#include <Riostream.h>

#include "TObject.h"

#include "TClonesArray.h"
#include "AliLMRMuon.h"

ClassImp(AliLMRMuon)

const Double_t AliLMRMuon::kMuonMass = 0.105658357;

AliLMRMuon::AliLMRMuon() : TObject() , 
fTriggerMatch(0),fpDCA(0.), fRabs(0.), fCharge(0), fPx(0.), 
fPy(0.), fPz(0.), fChi2(0.), fChi2Match(0.)
{
	
}

AliLMRMuon::AliLMRMuon(const AliLMRMuon& muon) : TObject(), 
	fTriggerMatch(muon.fTriggerMatch),
	fpDCA(muon.fpDCA),
	fRabs(muon.fRabs),
	fCharge(muon.fCharge),
	fPx(muon.fPx),
	fPy(muon.fPy),
	fPz(muon.fPz),
	fChi2(muon.fChi2),
	fChi2Match(muon.fChi2Match)
{
  ///copy constructor
    
}

TLorentzVector AliLMRMuon::P4() const
{
  TLorentzVector v(fPx,fPy,fPz,Energy());
  return v;
}
