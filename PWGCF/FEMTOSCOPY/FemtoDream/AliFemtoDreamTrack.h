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
#include "AliNanoAODTrack.h"
class AliFemtoDreamTrack : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamTrack();
  virtual ~AliFemtoDreamTrack();
  void SetTrack(AliAODTrack *track, const int multiplicity = -1);
  void SetTrack(AliVTrack *track, AliVEvent *event,
                const int multiplicity = -1);
  void SetTrack(AliESDtrack *track, AliMCEvent *mcEvent = nullptr,
                const int multiplicity = -1, const bool TPCOnlyTrack = true,
                const bool IsOmegaTrack = false);
  UInt_t GetFilterMap() const {
    return fFilterMap;
  }
  ;
  bool TestFilterBit(UInt_t filterBit) {
    return (bool) ((filterBit & fFilterMap) != 0);
  }
  bool PassESDFiltering() const {
    return fPassFiltering;
  }
  ;
  float GetDCAXY() const {
    return fdcaXY;
  }
  ;
  float GetDCAXYProp() const {
    return fdcaXYProp;
  }
  ;
  float GetDCAZ() const {
    return fdcaZ;
  }
  ;
  float GetDCAZProp() const {
    return fdcaZProp;
  }
  ;
  float GetChiSquare() const {
    return fChi2;
  }

  //Quality Varaibles of the track
  float GetNClsTPC() const {
    return fNClsTPC;
  }
  ;
  float GetTPCCrossedRows() const {
    return fTPCCrossedRows;
  }
  ;
  float GetRatioCr() const {
    return fRatioCR;
  }
  ;
  bool isnoSharedClst() const {
    return fnoSharedClst;
  }
  ;
  float GetTPCClsC() const {
    return fTPCClsS;
  }
  ;
  std::vector<bool> GetSharedClusterITS() const {
    return fSharedClsITSLayer;
  }
  bool GetSharedClusterITS(int i) const {
    return fSharedClsITSLayer.at(i);
  }
  ;
  bool GetHasSharedClsITS() const {
    return fHasSharedClsITSLayer;
  }
  ;
  bool GetHasSPDHit() const {
    return (fITSHit.at(0) || fITSHit.at(1)) ? true : false;
  }
  ;
  bool GetHasITSHit() const {
    return fHasITSHit;
  }
  ;
  std::vector<bool> GetITSHits() const { return fITSHit; }
  bool GetITSHit(int i) const {
    return fITSHit.at(i);
  }
  ;
  bool GetTOFTimingReuqirement() const {
    return fTOFTiming;
  }
  ;
  bool GetHasTPCRefit() const {
    return fTPCRefit;
  }
  ;
  //PID Getters
  AliPIDResponse::EDetPidStatus GetstatusTOF() const {
    return fstatusTOF;
  }
  ;
  AliPIDResponse::EDetPidStatus GetstatusTPC() const {
    return fstatusTPC;
  }
  ;
  AliPIDResponse::EDetPidStatus GetstatusITS() const {
    return fstatusITS;
  }
  ;
  float GetdEdxTPC() const {
    return fdEdxTPC;
  }
  ;
  float GetbetaTOF() const {
    return fbetaTOF;
  }
  ;
  float GetnSigmaITS(Int_t i) const {
    return fnSigmaITS[i];
  }
  ;
  float GetnSigmaTPC(Int_t i) const {
    return fnSigmaTPC[i];
  }
  ;
  float GetnSigmaTOF(Int_t i) const {
    return fnSigmaTOF[i];
  }
  ;
  TString ClassName() {
    return "AliFemtoDreamTrack";
  }
  ;
 private:
  AliFemtoDreamTrack &operator=(const AliFemtoDreamTrack &obj);
  AliFemtoDreamTrack(const AliFemtoDreamTrack&);
  void Reset();
  float GetBeta(AliNanoAODTrack *track) const;
  float GetBeta(AliAODTrack *track) const;
  float GetBeta(AliESDtrack *track) const;
  bool CheckGlobalTrack(const Int_t TrackID);
  bool CheckGlobalVTrack(const Int_t TrackID);
  void SetAODTrackingInformation();
  void SetVInformation(AliVEvent *event);
  void ApplyESDtoAODFilter(const bool TPCOnlyTrack = true);
  void SetESDTrackingInformation(const bool TPCOnlyTrack = true);
  void SetESDTrackingInformationOmega();
  void SetPhiAtRadii(const float bfield);
  void SetGlobalCoordAtRadii(const float bfield);
  void SetAODPIDInformation();
  void SetESDPIDInformation();
  void SetMCInformation();
  void SetMCInformation(AliVEvent* event);
  void SetMCInformation(TClonesArray* array, int trackID);
  void SetMCInformation(AliMCEvent *mcEvent);
  AliPIDResponse *fPIDResponse;
  AliPIDResponse::EDetPidStatus fstatusITS;
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
  float fnSigmaITS[6];
  float fnSigmaTPC[6];
  float fnSigmaTOF[6];
  ULong_t fESDStatus;
  int fESDnClusterITS;
  int fESDnClusterTPC;
  AliESDtrack *fESDTrack;
  AliESDtrack *fESDTPCOnlyTrack;
  AliESDtrackCuts *fESDTrackCuts;
  AliVTrack *fVTrack;
  AliVTrack *fVGlobalTrack;
  AliAODTrack *fAODTrack;
  AliAODTrack *fAODGlobalTrack;
  ClassDef(AliFemtoDreamTrack,5)
};

#endif /* ALIFEMTODREAMTRACK_H_ */
