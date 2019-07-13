#ifndef ALILIGHTNTRACK_H
#define ALILIGHTNTRACK_H

/*
 * AliLightNTrack.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include <vector>
#include "AliLightNBasePart.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include <iostream>
#include "AliLog.h"
#include "TClonesArray.h"
#include "Rtypes.h"
class AliLightNTrack : public AliLightNBasePart {
public:
    AliLightNTrack();
    virtual ~AliLightNTrack();
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
    double GetNClsITS() const {return fNClsITS;};
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
    AliPIDResponse::EDetPidStatus   GetstatusITS() const {return fstatusITS;};
    double GetdEdxTPC() const {return fdEdxTPC;};
    double GetbetaTOF() const {return fbetaTOF;};
    double GetdEdxITS() const {return fdEdxITS;};
    double GetMassSquare() const {return fmass2sq;};
    double GetnSigmaTPC(Int_t i) const {return fnSigmaTPC[i];};
    double GetnSigmaTOF(Int_t i) const {return fnSigmaTOF[i];};
    double GetnSigmaITS(Int_t i) const {return fnSigmaITS[i];};
    double GetRapidity(AliPID::EParticleType pid) const {return (fTrack->Y(pid));};			//dy=0.465 for p-Pb collisions... return (fTrack->Y(pid)-0.465);
    TString ClassName(){return "TrackCuts";};
private:
    void Reset();
    float GetBeta(AliAODTrack *track);
    float GetMass2sq(Double_t beta,AliAODTrack *track);
    bool CheckGlobalTrack(const Int_t TrackID);
    void SetTrackingInformation();
    void SetPIDInformation();
    void SetMCInformation();
    AliPIDResponse *fPIDResponse;
    AliPIDResponse::EDetPidStatus fstatusTPC;
    AliPIDResponse::EDetPidStatus fstatusTOF;
    AliPIDResponse::EDetPidStatus fstatusITS;
    UInt_t fFilterMap;
    double fdcaXY;
    double fdcaZ;
    double fdcaXYProp;
    double fdcaZProp;
    double fNClsTPC;
    double fNClsITS;
    double fTPCCrossedRows;
    double fRatioCR;
    double fnoSharedClst;
    double fTPCClsS;
    std::vector<bool> fSharedClsITSLayer;
    bool fHasSharedClsITSLayer;
    double fdEdxTPC;
    double fbetaTOF;
    double fdEdxITS;
    double fmass2sq;
    bool fHasITSHit;
    std::vector<bool> fITSHit;
    bool fTOFTiming;
    bool fTPCRefit;
    AliAODTrack *fTrack;
    AliAODTrack *fGlobalTrack;
    double fnSigmaTPC[6];
    double fnSigmaTOF[6];
    double fnSigmaITS[6];
    ClassDef(AliLightNTrack, 1)
};

#endif /* ALILIGHTNTRACK_H */
