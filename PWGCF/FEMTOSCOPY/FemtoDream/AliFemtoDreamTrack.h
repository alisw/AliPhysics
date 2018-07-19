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
#include "AliESDtrackCuts.h"
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
  bool PassESDFiltering() const {return fPassFiltering;};
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
  AliFemtoDreamTrack &operator=(const AliFemtoDreamTrack &obj);
  AliFemtoDreamTrack(const AliFemtoDreamTrack&);
  void Reset();
  float GetBeta(AliAODTrack *track);
  float GetBeta(AliESDtrack *track);
  bool CheckGlobalTrack(const Int_t TrackID);
  void SetAODTrackingInformation();
  void ApplyESDtoAODFilter();
  void SetESDTrackingInformation();
  void SetPhiAtRadii();
  void SetAODPIDInformation();
  void SetESDPIDInformation();
  void SetMCInformation();
  AliPIDResponse *fPIDResponse;
  AliPIDResponse::EDetPidStatus fstatusTPC;
  AliPIDResponse::EDetPidStatus fstatusTOF;
  UInt_t fFilterMap;
  bool fPassFiltering;
  float fdcaXY;
  float fdcaZ;
  float fdcaXYProp;
  float fdcaZProp;
  float fNClsTPC;
  float fTPCCrossedRows;
  float fRatioCR;
  bool fnoSharedClst;
  float fTPCClsS;
  float fTPCClsSRatio;
  float fChi2;
  float fChi2TPC;
  float fChi2ITS;
  std::vector<bool> fSharedClsITSLayer;
  bool fHasSharedClsITSLayer;
  float fdEdxTPC;
  float fbetaTOF;
  bool fHasITSHit;
  std::vector<bool> fITSHit;
  bool fTOFTiming;
  bool fTPCRefit;
  float fnSigmaTPC[5];
  float fnSigmaTOF[5];
  ULong_t fESDStatus;
  int fESDnClusterITS;
  int fESDnClusterTPC;
  AliESDtrack *fESDTrack;
  AliESDtrack *fESDTPCOnlyTrack;
  AliESDtrackCuts *fESDTrackCuts;
  AliAODTrack *fAODTrack;
  AliAODTrack *fAODGlobalTrack;
  ClassDef(AliFemtoDreamTrack,3)
};

#endif /* ALIFEMTODREAMTRACK_H_ */
