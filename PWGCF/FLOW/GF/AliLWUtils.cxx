#include "AliLWUtils.h"
AliLWTPCTrack::AliLWTPCTrack():fPt(-1),fPhi(-999),fEta(-999),fTrFlag(1) {};
AliLWTPCTrack::AliLWTPCTrack(Float_t pt, Float_t phi, Float_t eta, Short_t trFlag):fPt(pt),fPhi(phi),fEta(eta),fTrFlag(trFlag) {};
AliLWTPCTrack::~AliLWTPCTrack() {};
Bool_t AliLWTPCTrack::IsEqual(TObject *obj) const
{
  AliLWTPCTrack *l_Tr = (AliLWTPCTrack*)obj;
  if(!l_Tr) return kFALSE;
  if(fPhi == l_Tr->fPhi &&
     fEta == l_Tr->fEta &&
     fPt  == l_Tr->fPt  &&
     fTrFlag== l_Tr->fTrFlag)
      return kTRUE;
  return kFALSE;
};
Int_t AliLWTPCTrack::Compare(TObject *obj) const
{
  AliLWTPCTrack *l_Tr = (AliLWTPCTrack*)obj;
  if(fPt < l_Tr->fPt) return -1;
  else if(fPt > l_Tr->fPt) return 1;
  else return 0;
};
Bool_t AliLWTPCTrack::IsSortable() const {return kTRUE; };
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
