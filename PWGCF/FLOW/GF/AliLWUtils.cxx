#include "AliLWUtils.h"
AliLWTPCTrack::AliLWTPCTrack():fPt(-1),fPhi(-999),fEta(-999),fTrFlag(1) {};
AliLWTPCTrack::AliLWTPCTrack(Float_t pt, Float_t phi, Float_t eta, Short_t trFlag):fPt(pt),fPhi(phi),fEta(eta),fTrFlag(trFlag) {};
AliLWTPCTrack::~AliLWTPCTrack() {};
AliLWFMDTrack::AliLWFMDTrack():fPhi(-999),fEta(-999),fMult(-1) {};
AliLWFMDTrack::AliLWFMDTrack(Float_t phi, Float_t eta, Short_t mult):fPhi(phi),fEta(eta),fMult(mult) {};
AliLWFMDTrack::~AliLWFMDTrack() {};
AliLWEvent::AliLWEvent():fRunNo(0),fVz(-999),fCent(-1),fEvFlag(1) {};
AliLWEvent::AliLWEvent(UInt_t runNo, Float_t vz, Float_t cent, Short_t evFlag):fRunNo(runNo),fVz(vz),fCent(cent),fEvFlag(evFlag) {};
AliLWEvent::~AliLWEvent() {};
void AliLWEvent::Setup(UInt_t runNo, Float_t vz, Float_t cent, Short_t evFlag) {
  fRunNo = runNo;
  fVz = vz;
  fCent = cent;
  fEvFlag = evFlag;
}
