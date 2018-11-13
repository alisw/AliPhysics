#ifndef ALIHFTREEHANDLER_H
#define ALIHFTREEHANDLER_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandler
// \brief helper class to handle a tree for cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include "vector"
#include <TTree.h>
#include "AliAODTrack.h"
#include "AliAODPidHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODMCParticle.h"

using std::vector;

class AliHFTreeHandler : public TObject
{
  public:
  
    enum candtype {
      kSelected = BIT(0),
      kSignal   = BIT(1),
      kBkg      = BIT(2),
      kPrompt   = BIT(3),
      kFD       = BIT(4),
      kRefl     = BIT(5),
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
      kRawAndNsigmaPID
    };

    enum piddet {
      kTPC,
      kTOF
    };

    AliHFTreeHandler();
    AliHFTreeHandler(int PIDopt);

    virtual ~AliHFTreeHandler();

    //core methods --> implemented in each derived class
    virtual TTree* BuildTree(TString name, TString title) = 0;
    virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF* pidHF) = 0;
    //for MC gen --> common implementation
    TTree* BuildTreeMCGen(TString name, TString title);
    bool SetMCGenVariables(AliAODMCParticle* mcpart);

    virtual void FillTree() = 0; //to be called for each event, not each candidate!
    
    //common methods
    void SetOptPID(int PIDopt) {fPidOpt=PIDopt;}
    void SetFillOnlySignal(bool fillopt=true) {fFillOnlySignal=fillopt;}

    void SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected);
    void SetIsSelectedStd(bool isselected) {
      if(isselected) fCandTypeMap |= kSelected;
      else fCandTypeMap &= ~kSelected;
    }

    void SetDauInAcceptance(bool dauinacc = true) {fDauInAccFlag=dauinacc;}
  
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

  protected:  
    //constant variables
    static const unsigned int knMaxProngs   = 3;
    static const unsigned int knMaxDet4Pid  = 2;
    static const unsigned int knMaxHypo4Pid = 3;

    const float kCSPEED = 2.99792457999999984e-02; // cm / ps

    //helper methods for derived clases (to be used in BuildTree and SetVariables functions)
    void AddCommonDmesonVarBranches();
    void AddSingleTrackBranches();
    void AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
    bool SetSingleTrackVars(AliAODTrack* prongtracks[]);
    bool SetPidVars(AliAODTrack* prongtracks[], AliAODPidHF* pidHF, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
    void ResetDmesonCommonVarVectors();
    void ResetSingleTrackVarVectors();
    void ResetPidVarVectors();
    void ResetMCGenVectors();
  
    //utils methods
    double CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF);
    int RoundFloatToInt(double num);
    float ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* cand, float bfield);
    float GetTOFmomentum(AliAODTrack* track, AliPIDResponse* pidrespo);
  
    TTree* fTreeVar; /// tree with variables
    unsigned int fNProngs; /// number of prongs
    unsigned int fNCandidates; /// number of candidates in one fill (event)
    int fCandTypeMap; ///flag for candidate type (bit map above)
    vector<int> fCandType; ///vector of flag for candidate type (bit map above)
    vector<float> fInvMass; ///vector of candidate invariant mass
    vector<float> fPt; ///vector of candidate pt
    vector<float> fY; ///vector of candidate rapidity
    vector<float> fEta; ///vector of candidate pseudorapidity
    vector<float> fPhi; ///vector of candidate azimuthal angle
    vector<float> fDecayLength; ///vector of candidate decay length
    vector<float> fDecayLengthXY; ///vector of candidate decay length in the transverse plane
    vector<float> fNormDecayLengthXY; ///vector of candidate normalised decay length in the transverse plane
    vector<float> fCosP; ///vector of candidate cosine of pointing angle
    vector<float> fCosPXY; ///vector of candidate cosine of pointing angle in the transcverse plane
    vector<float> fImpParXY; ///vector of candidate impact parameter in the transverse plane
    vector<float> fPProng[knMaxProngs]; ///vectors of prong momentum
    vector<float> fTPCPProng[knMaxProngs]; ///vectors of prong TPC momentum
    vector<float> fTOFPProng[knMaxProngs]; ///vectors of prong TOF momentum
    vector<float> fPtProng[knMaxProngs]; ///vectors of prong pt
    vector<float> fEtaProng[knMaxProngs]; ///vectors of prong pseudorapidity
    vector<float> fPhiProng[knMaxProngs]; ///vectors of prong azimuthal angle
    vector<int> fNTPCclsProng[knMaxProngs]; ///vectors of prong track number of clusters in TPC
    vector<int> fNTPCclsPidProng[knMaxProngs]; ///vectors of prong track number of clusters in TPC used for PID
    vector<float> fNTPCCrossedRowProng[knMaxProngs]; ///vectors of prong track crossed row in TPC
    vector<float> fChi2perNDFProng[knMaxProngs]; ///vectors of prong track chi2/ndf
    vector<int> fNITSclsProng[knMaxProngs]; ///vectors of prong track number of clusters in ITS
    vector<int> fITSclsMapProng[knMaxProngs];///vectors of prong track ITS cluster map
    vector<float> fTrackIntegratedLengthProng[knMaxProngs]; /// vectors of prong track integrated lengths
    vector<float> fStartTimeResProng[knMaxProngs]; /// vectors of prong track start time resolutions (for TOF)
    vector<float> fPIDNsigmaVector[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid]; ///vectors of PID nsigma variables
    vector<int> fPIDNsigmaIntVector[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid]; ///vectors of PID nsigma variables (integers)
    vector<float> fPIDrawVector[knMaxProngs][knMaxDet4Pid]; ///vectors of raw PID variables
    int fPidOpt; /// option for PID variables
    bool fFillOnlySignal; ///flag to enable only signal filling
    bool fIsMCGenTree; ///flag to know if is a tree for MC generated particles
    bool fDauInAccFlag; ///flag to know if the daughter are in acceptance in case of MC gen
    vector<bool> fDauInAcceptance; ///vector of flags to know if the daughter are in acceptance in case of MC gen
  
  /// \cond CLASSIMP
  ClassDef(AliHFTreeHandler,1); /// 
  /// \endcond
};

#endif
