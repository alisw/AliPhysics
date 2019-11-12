#ifndef ALIHFMLVARHANDLER_H
#define ALIHFMLVARHANDLER_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandler
// \brief helper class to handle a tree and variables for ML analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TString.h>

#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODPidHF.h"

class AliHFMLVarHandler : public TObject
{
    public:
    
        enum candtype {
            kSignal = BIT(0),
            kBkg    = BIT(1),
            kPrompt = BIT(2),
            kFD     = BIT(3),
            kRefl   = BIT(4) //up to BIT(8) included for general flags, following BITS particle-specific
        };
    
        enum optpid {
            kNoPID,
            kRawPID,
            kNsigmaPID,
            kNsigmaCombPID,
            kNsigmaDetAndCombPID,
            kRawAndNsigmaPID
        };

        enum piddet {
            kTPC,
            kTOF,
            kCombTPCTOF // must be the last element in the enum
        };

        AliHFMLVarHandler();
        AliHFMLVarHandler(int PIDopt);
        virtual ~AliHFMLVarHandler();

        //core methods --> implemented in each derived class
        virtual TTree* BuildTree(TString /*name*/, TString /*title*/) {return nullptr;}
        virtual bool SetVariables(AliAODRecoDecayHF* /*cand*/, float /*bfield*/, int /*masshypo*/, AliAODPidHF* /*pidrespo*/) {return false;} //to be called for each candidate
        //to be called for each candidate
        void SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected);
        void FillTree();
        
        //common methods
        void SetOptPID(int PIDopt) {fPidOpt=PIDopt;}
        void SetAddSingleTrackVars(bool add) {fAddSingleTrackVar=add;}
        void SetFillOnlySignal(bool fillopt=true) {fFillOnlySignal=fillopt;}

    protected:  
        //constant variables
        static const unsigned int knMaxProngs   = 4;
        static const unsigned int knMaxDet4Pid  = 2;
        static const unsigned int knMaxHypo4Pid = 3;

        const float kCSPEED = 2.99792457999999984e-02; // cm / ps

        //helper methods for derived clases (to be used in BuildTree and SetVariables functions)
        void AddCommonDmesonVarBranches();
        void AddSingleTrackBranches();
        void AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
        bool SetSingleTrackVars(AliAODTrack* prongtracks[]);
        bool SetPidVars(AliAODTrack* prongtracks[], AliAODPidHF* pidrespo, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF);
    
        //utils methods
        double CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF);
        float ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* cand, float bfield);
        float GetTOFmomentum(AliAODTrack* track, AliAODPidHF* pidrespo);

        TTree* fTreeVar; /// tree with variables
        unsigned int fNProngs; /// number of prongs
        int fCandType; ///flag for candidate type (bit map above)
        float fInvMass; ///candidate invariant mass
        float fPt; ///candidate pt
        float fDecayLength; ///candidate decay length
        float fDecayLengthXY; ///candidate decay length in the transverse plane
        float fNormDecayLengthXY; ///candidate normalised decay length in the transverse plane
        float fCosP; ///candidate cosine of pointing angle
        float fCosPXY; ///candidate cosine of pointing angle in the transcverse plane
        float fImpParXY; ///candidate impact parameter in the transverse plane
        float fDCA; ///DCA of candidates prongs
        float fPtProng[knMaxProngs]; ///prong pt
        float fTPCPProng[knMaxProngs]; ///prong TPC momentum
        int fNTPCclsPidProng[knMaxProngs]; ///prong track number of clusters in TPC used for PID
        float fTOFPProng[knMaxProngs]; ///prong TOF momentum
        float fTrackIntegratedLengthProng[knMaxProngs]; /// prong track integrated lengths
        float fStartTimeResProng[knMaxProngs]; /// prong track start time resolutions (for TOF)
        float fPIDNsigmaVector[knMaxProngs][knMaxDet4Pid+1][knMaxHypo4Pid]; ///PID nsigma variables
        float fPIDrawVector[knMaxProngs][knMaxDet4Pid]; ///raw PID variables
        int fPidOpt; ///option for PID variables
        bool fAddSingleTrackVar; //add single-track variables
        bool fFillOnlySignal; ///flag to enable only signal filling

    /// \cond CLASSIMP
    ClassDef(AliHFMLVarHandler,1); ///
    /// \endcond
};
#endif
