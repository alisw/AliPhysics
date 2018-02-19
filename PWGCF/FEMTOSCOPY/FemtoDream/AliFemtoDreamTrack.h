/*
 * AliFemtoDreamTrack.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMTRACK_H_
#define ALIFEMTODREAMTRACK_H_

#include <vector>
#include "AliFemtoDreamBasePart.h"
#include "AliPIDResponse.h"
class AliFemtoDreamTrack : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamTrack();
  virtual ~AliFemtoDreamTrack();
  void SetTrack(AliAODTrack *track);
  UInt_t GetilterMap() const {return fFilterMap;};
  bool TestFilterBit(UInt_t filterBit)
  {return (bool) ((filterBit & fFilterMap) != 0);}

  double GetDCAXY() const {  return fdcaXY;};
  double GetDCAXYProp() const {return fdcaXYProp;};
  double GetDCAZ() const {return fdcaZ;};
  double GetDCAZProp() const {return fdcaZProp;};

  //Quality Varaibles of the track
  double GetNClsTPC() const {return fNClsTPC;};
  float GetTPCCrossedRows() const {return fTPCCrossedRows;};
  float GetRatioCr() const {return fRatioCR;};
  bool isnoSharedClst() const {return fnoSharedClst;};
  double GetTPCClsC() const {return fTPCClsS;};
  bool GetSharedClusterITS(int i)const{return fSharedClsITSLayer.at(i);};
  bool GetHasSharedClsITS()const{return fHasSharedClsITSLayer;};
  bool GetHasITSHit() const {return fHasITSHit;};
  bool GetITSHit(int i) const {return fITSHit.at(i);};
  bool GetTOFTimingReuqirement() const {return fTOFTiming;};
  bool GetHasTPCRefit()const{return fTPCRefit;};
  //PID Getters
  AliPIDResponse::EDetPidStatus   GetstatusTOF() const {return fstatusTOF;};
  AliPIDResponse::EDetPidStatus   GetstatusTPC() const {return fstatusTPC;};
  double GetdEdxTPC() const {return fdEdxTPC;};
  double GetbetaTOF() const {return fbetaTOF;};
  double GetnSigmaTPC(Int_t i) const {return fnSigmaTPC[i];};
  double GetnSigmaTOF(Int_t i) const {return fnSigmaTOF[i];};
  TString ClassName(){return "TrackCuts";};
 private:
  void Reset();
  float GetBeta(AliAODTrack *track);
  bool CheckGlobalTrack(const Int_t TrackID);
  void SetTrackingInformation();
  void SetPIDInformation();
  void SetMCInformation();
  AliPIDResponse *fPIDResponse;
  AliPIDResponse::EDetPidStatus fstatusTPC;
  AliPIDResponse::EDetPidStatus fstatusTOF;
  UInt_t fFilterMap;
  double fdcaXY;
  double fdcaZ;
  double fdcaXYProp;
  double fdcaZProp;
  double fNClsTPC;
  double fTPCCrossedRows;
  double fRatioCR;
  double fnoSharedClst;
  double fTPCClsS;
  std::vector<bool> fSharedClsITSLayer;
  bool fHasSharedClsITSLayer;
  double fdEdxTPC;
  double fbetaTOF;
  bool fHasITSHit;
  std::vector<bool> fITSHit;
  bool fTOFTiming;
  bool fTPCRefit;
  AliAODTrack *fTrack;
  AliAODTrack *fGlobalTrack;
  double fnSigmaTPC[5];
  double fnSigmaTOF[5];
  ClassDef(AliFemtoDreamTrack, 1)
};

#endif /* ALIFEMTODREAMTRACK_H_ */
