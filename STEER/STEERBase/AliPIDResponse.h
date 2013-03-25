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

#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliITSPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliTRDPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliHMPIDPIDResponse.h"
#include "AliEMCALPIDResponse.h"
#include "AliPID.h"

#include "TNamed.h"

class TF1;
class TObjArray;
class TLinearFitter;

class AliVEvent;
class AliTRDPIDResponseObject;
class AliTOFPIDParams;
class AliHMPIDPIDParams;
class AliOADBContainer;

class AliPIDResponse : public TNamed {
public:
  AliPIDResponse(Bool_t isMC=kFALSE);
  virtual ~AliPIDResponse();

  enum EDetector {
    kITS=0,
    kTPC=1,
    kTRD=2,
    kTOF=3,
    kHMPID=4,
    kEMCAL=5,
    kPHOS=6,
    kNdetectors=7
  };
  
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
    kDetMismatch=2,
    kDetNoParams=3
  };

  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}
  AliTRDPIDResponse &GetTRDResponse() {return fTRDResponse;}
  AliEMCALPIDResponse &GetEMCALResponse() {return fEMCALResponse;}

  // -----------------------------------------
  // buffered getters
  //
  
  // Number of sigmas
  EDetPidStatus NumberOfSigmas(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type, Double_t &val) const;
  
  Float_t NumberOfSigmas(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type) const;
  
  virtual Float_t NumberOfSigmasITS  (const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC  (const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC  (const AliVParticle *track, AliPID::EParticleType type, AliTPCPIDResponse::ETPCdEdxSource dedxSource) const;
  virtual Float_t NumberOfSigmasEMCAL(const AliVParticle *track, AliPID::EParticleType type, Double_t &eop, Double_t showershape[4]) const;
  virtual Float_t NumberOfSigmasTOF  (const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTOF  (const AliVParticle *track, AliPID::EParticleType type, const Float_t /*timeZeroTOF*/) const { return NumberOfSigmasTOF(track,type); }
  virtual Float_t NumberOfSigmasHMPID(const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasEMCAL(const AliVParticle *track, AliPID::EParticleType type) const;

  Bool_t IdentifiedAsElectronTRD(const AliVTrack *track, Double_t efficiencyLevel,Double_t centrality=-1,AliTRDPIDResponse::ETRDPIDMethod PIDmethod=AliTRDPIDResponse::kLQ1D) const;

  // Signal delta
  EDetPidStatus GetSignalDelta(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type, Double_t &val, Bool_t ratio=kFALSE) const;
  Double_t GetSignalDelta(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type, Bool_t ratio=kFALSE) const;
  
  // Probabilities
  EDetPidStatus ComputePIDProbability  (EDetCode  detCode, const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputePIDProbability  (EDetector detCode, const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  
  virtual EDetPidStatus ComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  virtual EDetPidStatus ComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  virtual EDetPidStatus ComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  virtual EDetPidStatus ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  virtual EDetPidStatus ComputeEMCALProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  virtual EDetPidStatus ComputePHOSProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  virtual EDetPidStatus ComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;

  virtual EDetPidStatus ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[],AliTRDPIDResponse::ETRDPIDMethod PIDmethod) const;
  
  // pid status
  EDetPidStatus CheckPIDStatus(EDetector detCode, const AliVTrack *track)  const;

  Float_t GetTOFMismatchProbability(const AliVTrack *track) const;
  
  void SetITSPIDmethod(ITSPIDmethod pmeth) { fITSPIDmethod = pmeth; }
  
  void SetOADBPath(const char* path) {fOADBPath=path;}
  const char *GetOADBPath() const {return fOADBPath.Data();}

  void SetCustomTPCpidResponse(const char* tpcpid) { fCustomTPCpidResponse = tpcpid; }
  const char* GetCustomTPCpidResponse() const { return fCustomTPCpidResponse.Data(); }
  
  void InitialiseEvent(AliVEvent *event, Int_t pass, Int_t run=-1);
  void SetCurrentFile(const char* file) { fCurrentFile=file; }

  // cache PID in the track
  void SetCachePID(Bool_t cache)    { fCachePID=cache;  }
  Bool_t GetCachePID() const { return fCachePID; }
  void FillTrackDetectorPID(const AliVTrack *track, EDetector detector) const;
  void FillTrackDetectorPID();

  AliVEvent * GetCurrentEvent() const {return fCurrentEvent;}

  // User settings for the MC period and reco pass
  void SetMCperiod(const char *mcPeriod) {fMCperiodUser=mcPeriod;}
  void SetRecoPass(Int_t recoPass)       {fRecoPassUser=recoPass;}

  // event info
  Float_t GetCurrentCentrality() const {return fCurrCentrality;};

  // TPC setting
  void SetUseTPCEtaCorrection(Bool_t useEtaCorrection = kTRUE) { fUseTPCEtaCorrection = useEtaCorrection; };
  Bool_t UseTPCEtaCorrection() const { return fUseTPCEtaCorrection; };
  
  // TOF setting
  void SetTOFtail(Float_t tail=1.1){if(tail > 0) fTOFtail=tail; else printf("TOF tail should be greater than 0 (nothing done)\n");};
  void SetTOFResponse(AliVEvent *vevent,EStartTimeType_t option);

  virtual Float_t GetTPCsignalTunedOnData(const AliVTrack *t) const {return t->GetTPCsignal();};
  Bool_t IsTunedOnData() const {return fTuneMConData;};
  void SetTunedOnData(Bool_t flag=kTRUE,Int_t recoPass=0){fTuneMConData = flag; if(recoPass>0) fRecoPassUser = recoPass;};

  AliPIDResponse(const AliPIDResponse &other);
  AliPIDResponse& operator=(const AliPIDResponse &other);

  
protected:
  AliITSPIDResponse   fITSResponse;    //PID response function of the ITS
  AliTPCPIDResponse   fTPCResponse;    //PID response function of the TPC
  AliTRDPIDResponse   fTRDResponse;    //PID response function of the TRD
  AliTOFPIDResponse   fTOFResponse;    //PID response function of the TOF
  AliHMPIDPIDResponse fHMPIDResponse;  //PID response function of the HMPID
  AliEMCALPIDResponse fEMCALResponse;  //PID response function of the EMCAL

  Float_t           fRange;          // nSigma max in likelihood
  ITSPIDmethod      fITSPIDmethod;   // 0 = trunc mean; 1 = likelihood

  //unbuffered PID calculation
  virtual Float_t GetNumberOfSigmasTOFold  (const AliVParticle */*track*/, AliPID::EParticleType /*type*/) const {return 0;}
  virtual Float_t GetSignalDeltaTOFold(const AliVParticle */*track*/, AliPID::EParticleType /*type*/, Bool_t /*ratio*/=kFALSE) const {return -9999.;}
  
  EDetPidStatus GetComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[],AliTRDPIDResponse::ETRDPIDMethod PIDmethod=AliTRDPIDResponse::kLQ1D) const;
  EDetPidStatus GetTOFPIDStatus(const AliVTrack *track) const;

private:
  Bool_t fIsMC;                        //  If we run on MC data
  Bool_t fCachePID;

  TString fOADBPath;                   // OADB path to use
  TString fCustomTPCpidResponse;       // Custom TPC Pid Response file for debugging purposes
  
  TString fBeamType;                   //! beam type (PP) or (PBPB)
  TString fLHCperiod;                  //! LHC period
  TString fMCperiodTPC;                //! corresponding MC period to use for the TPC splines
  TString fMCperiodUser;               //  MC prodution requested by the user
  TString fCurrentFile;                //! name of currently processed file
  Int_t   fRecoPass;                   //! reconstruction pass
  Int_t   fRecoPassUser;               //  reconstruction pass explicitly set by the user
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  Float_t fResT0A;                     //! T0A resolution in current run
  Float_t fResT0C;                     //! T0C resolution in current run
  Float_t fResT0AC;                    //! T0A.and.T0C resolution in current run
  
  TObjArray *fArrPidResponseMaster;    //!  TPC pid splines
  TF1       *fResolutionCorrection;    //! TPC resolution correction
  AliOADBContainer* fOADBvoltageMaps;  //! container with the voltage maps
  Bool_t fUseTPCEtaCorrection;         // Use TPC eta correction

  AliTRDPIDResponseObject *fTRDPIDResponseObject; //! TRD PID Response Object

  Float_t fTOFtail;                    //! TOF tail effect used in TOF probability
  AliTOFPIDParams *fTOFPIDParams;      //! TOF PID Params - period depending (OADB loaded)
  
  AliHMPIDPIDParams *fHMPIDPIDParams;  //! HMPID PID Params (OADB loaded)

  TObjArray *fEMCALPIDParams;             //! EMCAL PID Params

  AliVEvent *fCurrentEvent;            //! event currently being processed

  Float_t fCurrCentrality;             //! current centrality
  
  Bool_t fTuneMConData;                // switch to force the MC to be similar to data (dE/dx)
  
  void ExecNewRun();
  
  //
  //setup parametrisations
  //

  //ITS
  void SetITSParametrisation();
  
  //TPC
  void SetTPCEtaMaps(Double_t refineFactorMapX = 6.0, Double_t refineFactorMapY = 6.0, Double_t refineFactorSigmaMapX = 6.0,
                     Double_t refineFactorSigmaMapY = 6.0);
  void SetTPCPidResponseMaster();
  void SetTPCParametrisation();
  Double_t GetTPCMultiplicityBin(const AliVEvent * const event);
  
  // TPC helpers for the eta maps
  void AddPointToHyperplane(TH2D* h, TLinearFitter* linExtrapolation, Int_t binX, Int_t binY);
  TH2D* RefineHistoViaLinearInterpolation(TH2D* h, Double_t refineFactorX = 6.0, Double_t refineFactorY = 6.0);

  //TRD
  void SetTRDPidResponseMaster();
  void InitializeTRDResponse();
  void SetTRDSlices(UInt_t TRDslicesForPID[2],AliTRDPIDResponse::ETRDPIDMethod method) const;

  //TOF
  void SetTOFPidResponseMaster();
  void InitializeTOFResponse();

  //HMPID
  void SetHMPIDPidResponseMaster();
  void InitializeHMPIDResponse();

  //EMCAL
  void SetEMCALPidResponseMaster();
  void InitializeEMCALResponse();

  //
  void SetRecoInfo();

  //-------------------------------------------------
  //unbuffered PID calculation
  //
  
  // Number of sigmas
  Float_t GetNumberOfSigmas(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type) const;
  Float_t GetNumberOfSigmasITS  (const AliVParticle *track, AliPID::EParticleType type) const;
  Float_t GetNumberOfSigmasTPC  (const AliVParticle *track, AliPID::EParticleType type) const;
  Float_t GetNumberOfSigmasTOF  (const AliVParticle *track, AliPID::EParticleType type) const;
  Float_t GetNumberOfSigmasHMPID(const AliVParticle *track, AliPID::EParticleType type) const;
  Float_t GetNumberOfSigmasEMCAL(const AliVParticle *track, AliPID::EParticleType type, Double_t &eop, Double_t showershape[4]) const;
  Float_t GetNumberOfSigmasEMCAL(const AliVParticle *track, AliPID::EParticleType type) const;

  Float_t GetBufferedNumberOfSigmas(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type) const;

  // Signal deltas
  EDetPidStatus GetSignalDeltaITS(const AliVParticle *track, AliPID::EParticleType type, Double_t &val, Bool_t ratio=kFALSE) const;
  EDetPidStatus GetSignalDeltaTPC(const AliVParticle *track, AliPID::EParticleType type, Double_t &val, Bool_t ratio=kFALSE) const;
  EDetPidStatus GetSignalDeltaTOF(const AliVParticle *track, AliPID::EParticleType type, Double_t &val, Bool_t ratio=kFALSE) const;
  EDetPidStatus GetSignalDeltaHMPID(const AliVParticle *vtrack, AliPID::EParticleType type, Double_t &val, Bool_t ratio=kFALSE) const;
  
  // Probabilities
  EDetPidStatus GetComputePIDProbability  (EDetector detCode,  const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus GetComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus GetComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus GetComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus GetComputeEMCALProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus GetComputePHOSProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus GetComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;

  // pid status
  EDetPidStatus GetPIDStatus(EDetector det, const AliVTrack *track) const;
  EDetPidStatus GetITSPIDStatus(const AliVTrack *track) const;
  EDetPidStatus GetTPCPIDStatus(const AliVTrack *track) const;
  EDetPidStatus GetTRDPIDStatus(const AliVTrack *track) const;
  EDetPidStatus GetHMPIDPIDStatus(const AliVTrack *track) const;
  EDetPidStatus GetPHOSPIDStatus(const AliVTrack *track) const;
  EDetPidStatus GetEMCALPIDStatus(const AliVTrack *track) const;

  ClassDef(AliPIDResponse, 11);  //PID response handling
};

#endif
