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
#include "AliEMCALPIDResponse.h"

#include "AliVParticle.h"
#include "AliVTrack.h"

#include "TNamed.h"

class AliVEvent;
class TF1;

class AliPIDResponse : public TNamed {
public:
  AliPIDResponse(Bool_t isMC=kFALSE);
  virtual ~AliPIDResponse();

  enum EDetCode {
    kDetITS = 0x1,
    kDetTPC = 0x2,
    kDetTRD = 0x4,
    kDetTOF = 0x8,
    kDetHMPID = 0x10,
    kDetEMCAL = 0x20,
    kDetPHOS = 0x40
  };

  enum EStartTimeType_t {kFILL_T0,kTOF_T0, kT0_T0, kBest_T0};

  enum ITSPIDmethod { kITSTruncMean, kITSLikelihood };

  enum EDetPidStatus {
    kDetNoSignal=0,
    kDetPidOk=1,
    kDetMismatch=2
  };
  
  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}
  AliTRDPIDResponse &GetTRDResponse() {return fTRDResponse;}
  AliEMCALPIDResponse &GetEMCALResponse() {return fEMCALResponse;}

  Float_t NumberOfSigmas(EDetCode detCode, const AliVParticle *track, AliPID::EParticleType type) const;
  
  virtual Float_t NumberOfSigmasITS(const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC(const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasEMCAL(const AliVTrack *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTOF(const AliVParticle *track, AliPID::EParticleType type) const = 0;
  virtual Bool_t IdentifiedAsElectronTRD(const AliVTrack *track, Double_t efficiencyLevel) const;

  EDetPidStatus ComputePIDProbability  (EDetCode detCode, const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  
  EDetPidStatus ComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t  p[]) const;
  EDetPidStatus ComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeEMCALProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputePHOSProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;


  void SetITSPIDmethod(ITSPIDmethod pmeth) { fITSPIDmethod = pmeth; }
  virtual void SetTOFResponse(AliVEvent */*event*/,EStartTimeType_t /*option*/) {;}
  void SetTRDslicesForPID(UInt_t slice1, UInt_t slice2) {fTRDslicesForPID[0]=slice1;fTRDslicesForPID[1]=slice2;}
  
  void SetOADBPath(const char* path) {fOADBPath=path;}
  void InitialiseEvent(AliVEvent *event, Int_t pass);
  void SetCurrentFile(const char* file) { fCurrentFile=file; }

  AliVEvent * GetCurrentEvent() const {return fCurrentEvent;}

  // User settings for the MC period and reco pass
  void SetMCperiod(const char *mcPeriod) {fMCperiodUser=mcPeriod;}
  void SetRecoPass(Int_t recoPass)       {fRecoPassUser=recoPass;}
  
  AliPIDResponse(const AliPIDResponse &other);
  AliPIDResponse& operator=(const AliPIDResponse &other);
  
protected:
  AliITSPIDResponse fITSResponse;    //PID response function of the ITS
  AliTPCPIDResponse fTPCResponse;    //PID response function of the TPC
  AliTRDPIDResponse fTRDResponse;    //PID response function of the TRD
  AliTOFPIDResponse fTOFResponse;    //PID response function of the TOF
  AliEMCALPIDResponse fEMCALResponse;  //PID response function of the EMCAL

  Float_t           fRange;          // nSigma max in likelihood
  ITSPIDmethod      fITSPIDmethod;   // 0 = trunc mean; 1 = likelihood

private:
  Bool_t fIsMC;                        //  If we run on MC data

  TString fOADBPath;                   // OADB path to use
  
  TString fBeamType;                   //! beam type (PP) or (PBPB)
  TString fLHCperiod;                  //! LHC period
  TString fMCperiodTPC;                //! corresponding MC period to use for the TPC splines
  TString fMCperiodUser;               //  MC prodution requested by the user
  TString fCurrentFile;                //! name of currently processed file
  Int_t   fRecoPass;                   //! reconstruction pass
  Int_t   fRecoPassUser;               //  reconstruction pass explicitly set by the user
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  
  TObjArray *fArrPidResponseMaster;    //!  TPC pid splines
  TF1       *fResolutionCorrection;    //! TPC resolution correction

  AliTRDPIDParams *fTRDPIDParams;       //! TRD PID Params
  AliTRDPIDReference *fTRDPIDReference; //! TRD PID References
  UInt_t fTRDslicesForPID[2];           //! TRD PID slices

  Int_t   fTOFTimeZeroType;            //! default start time type for tof (ESD)
  Float_t fTOFres;                     //! TOF resolution

  AliVEvent *fCurrentEvent;            //! event currently being processed
  
  void ExecNewRun();
  
  //
  //setup parametrisations
  //

  //ITS
  void SetITSParametrisation();
  
  //TPC
  void SetTPCPidResponseMaster();
  void SetTPCParametrisation();
  Double_t GetTPCMultiplicityBin(const AliVEvent * const event);

  //TRD
  void SetTRDPidResponseMaster();
  void InitializeTRDResponse();

  //TOF
  
  //
  void SetRecoInfo();
  
  ClassDef(AliPIDResponse,2);  //PID response handling
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
