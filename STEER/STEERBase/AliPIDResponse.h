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

class TF1;
class AliOADBContainer;
class TObjArray;

class AliVEvent;
class AliTRDPIDResponseObject;
class AliTOFPIDParams;

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
    kDetMismatch=2
  };

  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}
  AliTRDPIDResponse &GetTRDResponse() {return fTRDResponse;}
  AliEMCALPIDResponse &GetEMCALResponse() {return fEMCALResponse;}

  Float_t NumberOfSigmas(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type) const;
  Float_t NumberOfSigmas(EDetCode  detCode, const AliVParticle *track, AliPID::EParticleType type) const;
  
  virtual Float_t NumberOfSigmasITS  (const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC  (const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t NumberOfSigmasTPC  (const AliVParticle *track, AliPID::EParticleType type, AliTPCPIDResponse::ETPCdEdxSource dedxSource);
  virtual Float_t NumberOfSigmasEMCAL(const AliVParticle *track, AliPID::EParticleType type, Double_t &eop, Double_t showershape[4]) const;
  virtual Float_t NumberOfSigmasTOF  (const AliVParticle *track, AliPID::EParticleType type) const = 0;
  virtual Float_t NumberOfSigmasEMCAL(const AliVParticle *track, AliPID::EParticleType type) const;

  virtual Bool_t IdentifiedAsElectronTRD(const AliVTrack *track, Double_t efficiencyLevel) const;

  EDetPidStatus ComputePIDProbability  (EDetector detCode, const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputePIDProbability  (EDetCode  detCode, const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  
  EDetPidStatus ComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeEMCALProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputePHOSProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  EDetPidStatus ComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;

  void SetTRDPIDmethod(AliTRDPIDResponse::ETRDPIDMethod method=AliTRDPIDResponse::kLQ1D);
  
  void SetITSPIDmethod(ITSPIDmethod pmeth) { fITSPIDmethod = pmeth; }
  void SetTRDslicesForPID(UInt_t slice1, UInt_t slice2) {fTRDslicesForPID[0]=slice1;fTRDslicesForPID[1]=slice2;}
  
  void SetOADBPath(const char* path) {fOADBPath=path;}
  const char *GetOADBPath() const {return fOADBPath.Data();}

  void SetCustomTPCpidResponse(const char* tpcpid) { fCustomTPCpidResponse = tpcpid; }
  const char* GetCustomTPCpidResponse() const { return fCustomTPCpidResponse.Data(); }
  
  void InitialiseEvent(AliVEvent *event, Int_t pass, Int_t run=-1);
  void SetCurrentFile(const char* file) { fCurrentFile=file; }

  // cache PID in the track
  void FillTrackDetectorPID();

  AliVEvent * GetCurrentEvent() const {return fCurrentEvent;}

  // User settings for the MC period and reco pass
  void SetMCperiod(const char *mcPeriod) {fMCperiodUser=mcPeriod;}
  void SetRecoPass(Int_t recoPass)       {fRecoPassUser=recoPass;}

  // event info
  Float_t GetCurrentCentrality() const {return fCurrCentrality;};

  // TOF setting
  void SetTOFtail(Float_t tail=1.1){if(tail > 0) fTOFtail=tail; else printf("TOF tail should be greater than 0 (nothing done)\n");};
  void SetTOFResponse(AliVEvent *vevent,EStartTimeType_t option);

  virtual Float_t GetTPCsignalTunedOnData(const AliVTrack *t) const {return t->GetTPCsignal();};
  Bool_t IsTunedOnData() const {return fTuneMConData;};
  void SetTunedOnData(Bool_t flag=kTRUE,Int_t recoPass=0){fTuneMConData = flag; if(recoPass>0) fRecoPassUser = recoPass;};

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
  
  TObjArray *fArrPidResponseMaster;    //!  TPC pid splines
  TF1       *fResolutionCorrection;    //! TPC resolution correction
  AliOADBContainer* fOADBvoltageMaps;  //! container with the voltage maps

  AliTRDPIDResponseObject *fTRDPIDResponseObject; //! TRD PID Response Object
  UInt_t fTRDslicesForPID[2];           //! TRD PID slices

  Float_t fTOFtail;                    //! TOF tail effect used in TOF probability
  AliTOFPIDParams *fTOFPIDParams;      //! TOF PID Params - period depending (OADB loaded)

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
  void SetTPCPidResponseMaster();
  void SetTPCParametrisation();
  Double_t GetTPCMultiplicityBin(const AliVEvent * const event);

  //TRD
  void SetTRDPidResponseMaster();
  void InitializeTRDResponse();

  //TOF
  void SetTOFPidResponseMaster();
  void InitializeTOFResponse();

  //EMCAL
  void SetEMCALPidResponseMaster();
  void InitializeEMCALResponse();

  //
  void SetRecoInfo();
  
  ClassDef(AliPIDResponse, 9);  //PID response handling
};

#endif
