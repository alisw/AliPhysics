/*
 * AliFemtoDreamTrack.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMTRACK_H_
#define ALIFEMTODREAMTRACK_H_

#include <vector>
#include "AliESDtrack.h"
#include "AliFemtoDreamBasePart.h"
#include "AliPIDResponse.h"
class AliFemtoDreamTrack : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamTrack();
  virtual ~AliFemtoDreamTrack();
  void SetTrack(AliAODTrack *track);
  void SetTrack(AliESDtrack *track);
  UInt_t GetilterMap() const {return fFilterMap;};
  bool TestFilterBit(UInt_t filterBit)
  {return (bool) ((filterBit & fFilterMap) != 0);}

  float GetDCAXY() const {  return fdcaXY;};
  float GetDCAXYProp() const {return fdcaXYProp;};
  float GetDCAZ() const {return fdcaZ;};
  float GetDCAZProp() const {return fdcaZProp;};
  float GetChiSquare() const { return fChi2; }

  //Quality Varaibles of the track
  float GetNClsTPC() const {return fNClsTPC;};
  float GetTPCCrossedRows() const {return fTPCCrossedRows;};
  float GetRatioCr() const {return fRatioCR;};
  bool isnoSharedClst() const {return fnoSharedClst;};
  float GetTPCClsC() const {return fTPCClsS;};
  bool GetSharedClusterITS(int i)const{return fSharedClsITSLayer.at(i);};
  bool GetHasSharedClsITS()const{return fHasSharedClsITSLayer;};
  bool GetHasITSHit() const {return fHasITSHit;};
  bool GetITSHit(int i) const {return fITSHit.at(i);};
  bool GetTOFTimingReuqirement() const {return fTOFTiming;};
  bool GetHasTPCRefit()const{return fTPCRefit;};
  //PID Getters
  AliPIDResponse::EDetPidStatus   GetstatusTOF() const {return fstatusTOF;};
  AliPIDResponse::EDetPidStatus   GetstatusTPC() const {return fstatusTPC;};
  float GetdEdxTPC() const {return fdEdxTPC;};
  float GetbetaTOF() const {return fbetaTOF;};
  float GetnSigmaTPC(Int_t i) const {return fnSigmaTPC[i];};
  float GetnSigmaTOF(Int_t i) const {return fnSigmaTOF[i];};
  TString ClassName(){return "TrackCuts";};
 private:
  void Reset();
  float GetBeta(AliAODTrack *track);
  bool CheckGlobalTrack(const Int_t TrackID);
  void SetAODTrackingInformation();
  void SetESDTrackingInformation();
  void SetPhiAtRadii();
  void SetPIDInformation();
  void SetMCInformation();
  AliPIDResponse *fPIDResponse;
  AliPIDResponse::EDetPidStatus fstatusTPC;
  AliPIDResponse::EDetPidStatus fstatusTOF;
  UInt_t fFilterMap;
  float fdcaXY;
  float fdcaZ;
  float fdcaXYProp;
  float fdcaZProp;
  float fNClsTPC;
  float fTPCCrossedRows;
  float fRatioCR;
  bool fnoSharedClst;
  float fTPCClsS;
  float fChi2;
  std::vector<bool> fSharedClsITSLayer;
  bool fHasSharedClsITSLayer;
  float fdEdxTPC;
  float fbetaTOF;
  bool fHasITSHit;
  std::vector<bool> fITSHit;
  bool fTOFTiming;
  bool fTPCRefit;
  AliESDtrack *fESDTrack;
  AliAODTrack *fAODTrack;
  AliAODTrack *fAODGlobalTrack;
  float fnSigmaTPC[5];
  float fnSigmaTOF[5];
  UInt_t fESDStatus;
  int fESDnClusterITS;
  int fESDnClustersTPC;
  float fESDnCrossedRowsTPC;
  ClassDef(AliFemtoDreamTrack,3)
};

#endif /* ALIFEMTODREAMTRACK_H_ */
