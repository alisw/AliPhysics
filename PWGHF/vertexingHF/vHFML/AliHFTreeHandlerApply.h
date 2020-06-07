#ifndef ALIHFTREEHANDLERAPPLY_H
#define ALIHFTREEHANDLERAPPLY_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApply
// \brief helper class to handle a tree for cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODMCParticle.h"
#include "AliAODPidHF.h"

class AliHFTreeHandlerApply : public TObject
{
public:
  
  enum candtype {
    kSelected        = BIT(0),
    kSignal          = BIT(1),
    kBkg             = BIT(2),
    kPrompt          = BIT(3),
    kFD              = BIT(4),
    kRefl            = BIT(5),
    kSelectedTopo    = BIT(6),
    kSelectedPID     = BIT(7),
    kSelectedTracks  = BIT(8) //up to BIT(10) included for general flags, following BITS particle-specific
  };
  
  enum optpid {
    kNoPID,
    kNsigmaPID,
    kNsigmaPIDint,
    kNsigmaPIDfloatandint, //--> to test
    kNsigmaCombPID,
    kNsigmaCombPIDint,
    kNsigmaCombPIDfloatandint, //--> to test
    kRawPID,
    kRawAndNsigmaPID,
    kNsigmaDetAndCombPID
  };
  
  enum piddet {
    kTPC,
    kTOF,
    kCombTPCTOF // must be the last element in the enum
  };
  
  enum optsingletrack {
    kNoSingleTrackVars, // single-track vars off
    kRedSingleTrackVars, // only pT, p, eta, phi
    kAllSingleTrackVars // all single-track vars
  };
  
  AliHFTreeHandlerApply();
  AliHFTreeHandlerApply(int PIDopt);
  virtual ~AliHFTreeHandlerApply();
  
  AliHFTreeHandlerApply(const AliHFTreeHandlerApply &source) = delete;
  AliHFTreeHandlerApply& operator=(const AliHFTreeHandlerApply &source) = delete;
  
  //core methods --> implemented in each derived class
  virtual TTree* BuildTree(TString name, TString title) = 0;
  virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, float mlprob, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse* pidrespo) = 0;

  //for MC gen --> common implementation
  TTree* BuildTreeMCGen(TString name, TString title);
  bool SetMCGenVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, AliAODMCParticle* mcpart);

  //to be called for each candidate
  void SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected);
  void FillTree() { //to be called for each candidate!
    if(fFillOnlySignal && !(fCandType&kSignal) && !(fCandType&kRefl)) { //if fill only signal and not signal/reflection candidate, do not store
      fCandType=0;
    }
    else {
      fTreeVar->Fill();
      fCandType=0;
      fRunNumberPrevCand = fRunNumber;
    }
  }
  
  //common methods
  void SetOptPID(int PIDopt) {fPidOpt=PIDopt;}
  void SetOptSingleTrackVars(int opt) {fSingleTrackOpt=opt;}
  void SetFillOnlySignal(bool fillopt=true) {fFillOnlySignal=fillopt;}
  void SetDauInAcceptance(bool dauinacc = true) {fDauInAcceptance=dauinacc;}

  void SetIsSelectedStd(bool isselected, bool isselectedTopo, bool isselectedPID, bool isselectedTracks) {
    if(isselected) fCandType |= kSelected;
    else fCandType &= ~kSelected;
    if(isselectedTopo) fCandType |= kSelectedTopo;
    else fCandType &= ~kSelectedTopo;
    if(isselectedPID) fCandType |= kSelectedPID;
    else fCandType &= ~kSelectedPID;
    if(isselectedTracks) fCandType |= kSelectedTracks;
    else fCandType &= ~kSelectedTracks;
  }
  
  static bool IsSelectedStd(int candtype) {
    if(candtype&1) return true;
    return false;
  }
  static bool IsSignal(int candtype) {
    if(candtype>>1&1) return true;
    return false;
  }
  static bool IsBkg(int candtype) {
    if(candtype>>2&1) return true;
    return false;
  }
  static bool IsPrompt(int candtype) {
    if(candtype>>3&1) return true;
    return false;
  }
  static bool IsFD(int candtype) {
    if(candtype>>4&1) return true;
    return false;
  }
  static bool IsRefl(int candtype) {
    if(candtype>>5&1) return true;
    return false;
  }
  static bool IsSelectedStdTopo(int candtype) {
    if(candtype>>6&1) return true;
    return false;
  }
  static bool IsSelectedStdPID(int candtype) {
    if(candtype>>7&1) return true;
    return false;
  }
  static bool IsSelectedStdTracks(int candtype) {
    if(candtype>>8&1) return true;
    return false;
  }
  
  void EnableNsigmaTPCDataDrivenCorrection(int syst) {
    fApplyNsigmaTPCDataCorr=true;
    fSystNsigmaTPCDataCorr=syst;
  }
  
protected:
  //constant variables
  static const unsigned int knMaxProngs   = 4;
  static const unsigned int knMaxDet4Pid  = 2;
  static const unsigned int knMaxHypo4Pid = 3;
  
  const float kCSPEED = 2.99792457999999984e-02; // cm / ps
  
  //helper methods for derived clases (to be used in BuildTree and SetVariables functions)
  void AddCommonDmesonVarBranches(Bool_t HasSecVtx = kTRUE);
  void AddSingleTrackBranches();
  void AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
  bool SetSingleTrackVars(AliAODTrack* prongtracks[]);
  bool SetPidVars(AliAODTrack* prongtracks[], AliPIDResponse* pidrespo, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
  
  //utils methods
  double CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF);
  int RoundFloatToInt(double num);
  float ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* cand, float bfield);
  float GetTOFmomentum(AliAODTrack* track, AliPIDResponse* pidrespo);
  
  void GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC, float eta);
  
  TTree* fTreeVar;                                                     //!<! tree with variables
  unsigned int fNProngs;                                               /// number of prongs
  
  int fCandType;                                                       /// flag for candidate type (bit map above)
  float fInvMass;                                                      /// candidate invariant mass
  float fPt;                                                           /// candidate pt
  float fPtGen;                                                        /// generated candidate pt
  float fY;                                                            /// candidate rapidity
  float fEta;                                                          /// candidate pseudorapidity
  float fPhi;                                                          /// candidate azimuthal angle
  float fMLProb;                                                       /// candidate's applied ML probability
  float fDecayLength;                                                  /// candidate decay length
  float fDecayLengthXY;                                                /// candidate decay length in the transverse plane
  float fNormDecayLengthXY;                                            /// candidate normalised decay length in the transverse plane
  float fCosP;                                                         /// candidate cosine of pointing angle
  float fCosPXY;                                                       /// candidate cosine of pointing angle in the transcverse plane
  float fImpParXY;                                                     /// candidate impact parameter in the transverse plane
  float fDCA;                                                          /// DCA of candidates prongs
  
  float fPProng[knMaxProngs];                                          /// prong momentum
  int fSPDhitsProng[knMaxProngs];                                      /// prong hits in the SPD
  float fTPCPProng[knMaxProngs];                                       /// prong TPC momentum
  float fTOFPProng[knMaxProngs];                                       /// prong TOF momentum
  float fPtProng[knMaxProngs];                                         /// prong pt
  float fEtaProng[knMaxProngs];                                        /// prong pseudorapidity
  float fPhiProng[knMaxProngs];                                        /// prong azimuthal angle
  int fNTPCclsProng[knMaxProngs];                                      /// prong track number of clusters in TPC
  int fNTPCclsPidProng[knMaxProngs];                                   /// prong track number of clusters in TPC used for PID
  float fNTPCCrossedRowProng[knMaxProngs];                             /// prong track crossed row in TPC
  float fChi2perNDFProng[knMaxProngs];                                 /// prong track chi2/ndf
  int fNITSclsProng[knMaxProngs];                                      /// prong track number of clusters in ITS
  int fITSclsMapProng[knMaxProngs];                                    /// prong track ITS cluster map
  float fTrackIntegratedLengthProng[knMaxProngs];                      /// prong track integrated lengths
  float fStartTimeResProng[knMaxProngs];                               /// prong track start time resolutions (for TOF)
  float fPIDNsigmaVector[knMaxProngs][knMaxDet4Pid+1][knMaxHypo4Pid];  /// PID nsigma variables
  int fPIDNsigmaIntVector[knMaxProngs][knMaxDet4Pid+1][knMaxHypo4Pid]; /// PID nsigma variables (integers)
  float fPIDrawVector[knMaxProngs][knMaxDet4Pid];                      /// raw PID variables
  
  int fPidOpt;                                                         /// option for PID variables
  int fSingleTrackOpt;                                                 /// option for single-track variables
  bool fFillOnlySignal;                                                /// flag to enable only signal filling
  bool fIsMCGenTree;                                                   /// flag to know if is a tree for MC generated particles
  bool fDauInAcceptance;                                               /// flag to know if the daughter are in acceptance in case of MC gen
  
  int fEvID;                                                           /// event ID corresponding to the one set in fTreeEvChar, first 32 bit of fEvIDLong
  int fEvIDExt;                                                        /// event ID corresponding to the one set in fTreeEvChar, second 32 bit of fEvIDLong
  Long64_t fEvIDLong;                                                  /// event ID corresponding to the one set in fTreeEvChar, full fEvIDLong
  int fRunNumber;                                                      /// run number
  int fRunNumberPrevCand;                                              /// run number of previous candidate

  bool fApplyNsigmaTPCDataCorr;                                        /// flag to enable data-driven NsigmaTPC correction
  int fSystNsigmaTPCDataCorr;                                          /// system for data-driven NsigmaTPC correction
  vector<vector<float> > fMeanNsigmaTPCPionData;                       /// array of NsigmaTPC pion mean in data
  vector<vector<float> > fMeanNsigmaTPCKaonData;                       /// array of NsigmaTPC kaon mean in data
  vector<vector<float> > fMeanNsigmaTPCProtonData;                     /// array of NsigmaTPC proton mean in data
  vector<vector<float> > fSigmaNsigmaTPCPionData;                      /// array of NsigmaTPC pion mean in data
  vector<vector<float> > fSigmaNsigmaTPCKaonData;                      /// array of NsigmaTPC kaon mean in data
  vector<vector<float> > fSigmaNsigmaTPCProtonData;                    /// array of NsigmaTPC proton mean in data
  float fPlimitsNsigmaTPCDataCorr[AliAODPidHF::kMaxPBins+1];           /// array of p limits for data-driven NsigmaTPC correction
  int fNPbinsNsigmaTPCDataCorr;                                        /// number of p bins for data-driven NsigmaTPC correction
  float fEtalimitsNsigmaTPCDataCorr[AliAODPidHF::kMaxEtaBins+1];       /// vector of eta limits for data-driven NsigmaTPC correction
  int fNEtabinsNsigmaTPCDataCorr;                                      /// number of eta bins for data-driven NsigmaTPC correction
  
  /// \cond CLASSIMP
  ClassDef(AliHFTreeHandlerApply,2); ///
  /// \endcond
};
#endif
