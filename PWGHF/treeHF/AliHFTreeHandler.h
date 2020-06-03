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
// L. van Doremalen, lennart.van.doremalen@cern.ch
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODMCParticle.h"
#include "AliAODPidHF.h"
#include "AliHFJet.h"

#ifdef HAVE_FASTJET
#include "AliHFJetFinder.h"
#endif

class AliHFTreeHandler : public TObject
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

    AliHFTreeHandler();
    AliHFTreeHandler(int PIDopt);

    virtual ~AliHFTreeHandler();

    //core methods --> implemented in each derived class
    virtual TTree* BuildTree(TString name, TString title) = 0;
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse* pidrespo) = 0;
    //for MC gen --> common implementation
    TTree* BuildTreeMCGen(TString name, TString title);
    bool SetMCGenVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, AliAODMCParticle* mcpart);

    void SetJetVars(TClonesArray *array, AliAODRecoDecayHF* cand, Double_t invmass, TClonesArray *mcarray, AliAODMCParticle* mcPart);
    void SetAndFillInclusiveJetVars(TClonesArray *array,TClonesArray *mcarray);
    void SetGenJetVars(TClonesArray *array, AliAODMCParticle* mcPart);
    void SetAndFillInclusiveGenJetVars(TClonesArray *array);
#ifdef HAVE_FASTJET
    void SetJetParameters(AliHFJetFinder& hfjetfinder);
#endif
    void SetJetTreeVars(AliHFJet hfjet);
    void SetGenJetTreeVars(AliHFJet hfjet);


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
    void SetFillJets(bool FillJets) {fFillJets=FillJets;}
    void SetDoJetSubstructure(bool DoJetSubstructure) {fDoJetSubstructure=DoJetSubstructure;}
    void SetTrackingEfficiency(Double_t TrackingEfficiency) {fTrackingEfficiency=TrackingEfficiency;}
    void SetJetProperties(Double_t JetRadius,Int_t JetAlgorithm,Double_t MinJetPt) {fJetRadius=JetRadius;fJetAlgorithm=JetAlgorithm;fMinJetPt=MinJetPt;}
    void SetSubJetProperties(Double_t SubJetRadius,Int_t SubJetAlgorithm,Double_t SoftDropZCut,Double_t SoftDropBeta) {fSubJetRadius=SubJetRadius;fSubJetAlgorithm=SubJetAlgorithm;fSoftDropZCut=SoftDropZCut;fSoftDropBeta=SoftDropBeta;}
    void SetOptPID(int PIDopt) {fPidOpt=PIDopt;}
    void SetOptSingleTrackVars(int opt) {fSingleTrackOpt=opt;}
    void SetFillOnlySignal(bool fillopt=true) {fFillOnlySignal=fillopt;}

    void SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected);
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

    void SetDauInAcceptance(bool dauinacc = true) {fDauInAcceptance=dauinacc;}
  
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
    void AddJetBranches();
    void AddGenJetBranches();
    void AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
    bool SetSingleTrackVars(AliAODTrack* prongtracks[]);
    bool SetPidVars(AliAODTrack* prongtracks[], AliPIDResponse* pidrespo, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
  
    //utils methods
    double CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF);
    int RoundFloatToInt(double num);
    float ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* cand, float bfield);
    float GetTOFmomentum(AliAODTrack* track, AliPIDResponse* pidrespo);
  
    void GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC, float eta);

    TTree* fTreeVar; /// tree with variables
    unsigned int fNProngs; /// number of prongs
    unsigned int fNCandidates; /// number of candidates in one fill (event)
    int fCandType; ///flag for candidate type (bit map above)
    float fInvMass; ///candidate invariant mass
    float fPt; ///candidate pt
    float fPtGen; ///generated candidate pt
    float fY; ///candidate rapidity
    float fEta; ///candidate pseudorapidity
    float fPhi; ///candidate azimuthal angle
    float fDecayLength; ///candidate decay length
    float fDecayLengthXY; ///candidate decay length in the transverse plane
    float fNormDecayLengthXY; ///candidate normalised decay length in the transverse plane
    float fCosP; ///candidate cosine of pointing angle
    float fCosPXY; ///candidate cosine of pointing angle in the transcverse plane
    float fImpParXY; ///candidate impact parameter in the transverse plane
    float fDCA; ///DCA of candidates prongs
    float fPProng[knMaxProngs]; ///prong momentum
    int fSPDhitsProng[knMaxProngs]; ///prong hits in the SPD
    float fTPCPProng[knMaxProngs]; ///prong TPC momentum
    float fTOFPProng[knMaxProngs]; ///prong TOF momentum
    float fPtProng[knMaxProngs]; ///prong pt
    float fEtaProng[knMaxProngs]; ///prong pseudorapidity
    float fPhiProng[knMaxProngs]; ///prong azimuthal angle
    int fNTPCclsProng[knMaxProngs]; ///prong track number of clusters in TPC
    int fNTPCclsPidProng[knMaxProngs]; ///prong track number of clusters in TPC used for PID
    float fNTPCCrossedRowProng[knMaxProngs]; ///prong track crossed row in TPC
    float fChi2perNDFProng[knMaxProngs]; ///prong track chi2/ndf
    int fNITSclsProng[knMaxProngs]; ///prong track number of clusters in ITS
    int fITSclsMapProng[knMaxProngs];///prong track ITS cluster map
    float fTrackIntegratedLengthProng[knMaxProngs]; /// prong track integrated lengths
    float fStartTimeResProng[knMaxProngs]; /// prong track start time resolutions (for TOF)
    float fPIDNsigmaVector[knMaxProngs][knMaxDet4Pid+1][knMaxHypo4Pid]; ///PID nsigma variables
    int fPIDNsigmaIntVector[knMaxProngs][knMaxDet4Pid+1][knMaxHypo4Pid]; ///PID nsigma variables (integers)
    float fPIDrawVector[knMaxProngs][knMaxDet4Pid]; ///raw PID variables
    int fPidOpt; ///option for PID variables
    int fSingleTrackOpt; ///option for single-track variables
    bool fFillOnlySignal; ///flag to enable only signal filling
    bool fIsMCGenTree; ///flag to know if is a tree for MC generated particles
    bool fDauInAcceptance; ///flag to know if the daughter are in acceptance in case of MC gen
    int fEvID; ///event ID corresponding to the one set in fTreeEvChar, first 32 bit of fEvIDLong
    int fEvIDExt; ///event ID corresponding to the one set in fTreeEvChar, second 32 bit of fEvIDLong
    Long64_t fEvIDLong; ///event ID corresponding to the one set in fTreeEvChar, full fEvIDLong
    int fRunNumber; ///run number
    int fRunNumberPrevCand; ///run number of previous candidate
    bool fApplyNsigmaTPCDataCorr; /// flag to enable data-driven NsigmaTPC correction
    int fSystNsigmaTPCDataCorr; /// system for data-driven NsigmaTPC correction
    vector<vector<float> > fMeanNsigmaTPCPionData; /// array of NsigmaTPC pion mean in data 
    vector<vector<float> > fMeanNsigmaTPCKaonData; /// array of NsigmaTPC kaon mean in data 
    vector<vector<float> > fMeanNsigmaTPCProtonData; /// array of NsigmaTPC proton mean in data 
    vector<vector<float> > fSigmaNsigmaTPCPionData; /// array of NsigmaTPC pion mean in data 
    vector<vector<float> > fSigmaNsigmaTPCKaonData; /// array of NsigmaTPC kaon mean in data 
    vector<vector<float> > fSigmaNsigmaTPCProtonData; /// array of NsigmaTPC proton mean in data 
    float fPlimitsNsigmaTPCDataCorr[AliAODPidHF::kMaxPBins+1]; /// array of p limits for data-driven NsigmaTPC correction
    int fNPbinsNsigmaTPCDataCorr;/// number of p bins for data-driven NsigmaTPC correction
    float fEtalimitsNsigmaTPCDataCorr[AliAODPidHF::kMaxEtaBins+1]; /// vector of eta limits for data-driven NsigmaTPC correction
    int fNEtabinsNsigmaTPCDataCorr; /// number of eta bins for data-driven NsigmaTPC correction

    float fPtJet; ///jet pt
    float fPtGenJet; ///gen jet pt
    float fEtaJet; ///jet pseudorapidity
    float fEtaGenJet; ///gen jet pseudorapidity
    float fPhiJet; ///jet azimuthal angle
    float fPhiGenJet; ///gen jet azimuthal angle
    float fLeadingPtJet; //jet leading track pT
    float fLeadingPtGenJet; //genjet leading track pT
    float fDeltaEtaJetHadron; ///jet hadron pseudorapidity
    float fDeltaEtaGenJetHadron; ///gen jet hadron pseudorapidity
    float fDeltaPhiJetHadron; ///jet hadron azimuthal angle
    float fDeltaPhiGenJetHadron; ///jet hadron azimuthal angle
    float fDeltaRJetHadron; ///jet hadron distance
    float fDeltaRGenJetHadron; ///gen jet hadron distance
    float fNTracksJet;  //number of tracks in the jet
    float fNTracksGenJet;  //number of tracks in the gen jet
    float fZJet; // fragmentation function in jet
    float fZGenJet; //fragmentation function in gen jet
    float fAngularityk1B1Jet; // Angularity with kappa = 1 and beta =1 in jet
    float fAngularityk1B1GenJet; // Angularity with kappa = 1 and beta =1 in gen jet  
    float fpTDispersionJet; // pT dispersion in jet
    float fpTDispersionGenJet; // pT dispersion in gen jet
    float fChargek03Jet; // jet charge with kappa = 0.3 in jet
    float fChargek03GenJet; // jet charge with kappa = 0.3 in gen jet
    float fChargek05Jet; // jet charge with kappa = 0.5 in jet
    float fChargek05GenJet; // jet charge with kappa = 0.5 in gen jet
    float fChargek07Jet; // jet charge with kappa = 0.7 in jet
    float fChargek07GenJet; // jet charge with kappa = 0.7 in gen jet
    float fZgJet; //zg
    float fZgGenJet; //gen zg
    float fRgJet; //Rg
    float fRgGenJet; //gen Rg
    float fNsdJet; //Nsd
    float fNsdGenJet; //gen Nsd
    float fPt_splittingJet; //Pt_splitting
    float fPt_splittingGenJet; //gen Pt_splitting
    float fk0Jet; //k0
    float fk0GenJet; //gen k0
    float fZk0Jet; //Zk0
    float fZk0GenJet; //gen Zk0
    float fRk0Jet; //Rk0
    float fRk0GenJet; //gen Rk0
    float fk1Jet; //k1
    float fk1GenJet; //gen k1
    float fZk1Jet; //Zk1
    float fZk1GenJet; //gen Zk1
    float fRk1Jet; //Rk1
    float fRk1GenJet; //gen Rk1
    float fk2Jet; //k2
    float fk2GenJet; //gen k2
    float fZk2Jet; //Zk2
    float fZk2GenJet; //gen Zk2
    float fRk2Jet; //Rk2
    float fRk2GenJet; //gen Rk2
    float fkTJet; //kT
    float fkTGenJet; //gen kT
    float fZkTJet; //ZkT
    float fZkTGenJet; //gen ZkT
    float fRkTJet; //RkT
    float fRkTGenJet; //gen RkT
    bool  fFillJets; //fill jets
    bool  fDoJetSubstructure; //fill jet substructure
    Double_t fJetRadius; //Jet finding radius
    Double_t fSubJetRadius; //Subjet finding radius
    Int_t fJetAlgorithm; //Jet finding algorithm
    Int_t fSubJetAlgorithm; //SubJet finding algorithm
    Double_t fMinJetPt; //Jet finding mimimum Jet pT
    Double_t fSoftDropZCut; //soft drop z parameter
    Double_t fSoftDropBeta; //soft drop beta  parameter
    Double_t fTrackingEfficiency;

  /// \cond CLASSIMP
  ClassDef(AliHFTreeHandler,9); ///
  /// \endcond
};
#endif
