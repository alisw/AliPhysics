#ifndef ALIPIDRESPONSE_H
#define ALIPIDRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------//
//        Base class for handling the pid response               //
//        functions of all detectors                             //
//        and give access to the nsigmas                         //
//                                                               //
//   Origin: Jens Wiechula, Uni Tuebingen, jens.wiechula@cern.ch //
//---------------------------------------------------------------//

#include "AliITSPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliTRDPIDResponse.h"
#include "AliTOFPIDResponse.h"

#include "AliVParticle.h"
#include "AliVTrack.h"

#include "TNamed.h"

class AliVEvent;

class AliPIDResponse : public TNamed {
public:
  AliPIDResponse(Bool_t isMC=kFALSE) : fITSResponse(isMC), fTPCResponse(), fTRDResponse(), fTOFResponse() {;}
  virtual ~AliPIDResponse() {;}

  enum EStartTimeType_t {kFILL_T0,kTOF_T0, kT0_T0, kBest_T0};
  
  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}
  AliTRDPIDResponse &GetTRDResponse() {return fTRDResponse;}
  
  virtual Float_t NumberOfSigmasITS(const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC(const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTOF(const AliVParticle *track, AliPID::EParticleType type) const = 0;

  void SetTOFResponse(AliVEvent */*event*/,EStartTimeType_t /*option*/) {;}
  
protected:
  AliITSPIDResponse fITSResponse;    //PID response function of the ITS
  AliTPCPIDResponse fTPCResponse;    //PID response function of the TPC
  AliTRDPIDResponse fTRDResponse;    //PID response function of the TRD
  AliTOFPIDResponse fTOFResponse;    //PID response function of the TOF
  
  ClassDef(AliPIDResponse,1);  //PID response handling
};

inline Float_t AliPIDResponse::NumberOfSigmasTPC(const AliVParticle *vtrack, AliPID::EParticleType type) const {
  AliVTrack *track=(AliVTrack*)vtrack;
  Double_t mom  = track->GetTPCmomentum();
  Double_t sig  = track->GetTPCsignal();
  UInt_t   sigN = track->GetTPCsignalN();

  Double_t nSigma = -999.;
  if (sigN>0) nSigma=fTPCResponse.GetNumberOfSigmas(mom,sig,sigN,type);

  return nSigma;
}

inline Float_t AliPIDResponse::NumberOfSigmasITS(const AliVParticle *vtrack, AliPID::EParticleType type) const {
  AliVTrack *track=(AliVTrack*)vtrack;
  Float_t dEdx=track->GetITSsignal();
  if (dEdx<=0) return -999.;
  
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  Float_t mom=track->P();
  Bool_t isSA=kTRUE;
  if(track->GetTPCNcls()>0) isSA=kFALSE;
  return fITSResponse.GetNumberOfSigmas(mom,dEdx,type,nPointsForPid,isSA);
}

#endif
