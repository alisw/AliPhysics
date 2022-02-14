#ifndef AliHFMLVARHANDLERNONPROMPTLC2V0BACHELOR_H
#define AliHFMLVARHANDLERNONPROMPTLC2V0BACHELOR_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandlerNonPromptLc2V0bachelor
// \brief helper class to handle a tree and variables for non-prompt Lc->V0bachelor ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"

class AliHFMLVarHandlerNonPromptLc2V0bachelor : public AliHFMLVarHandler
{
    public:

        enum massHypo {
            kpK0s,
            kpiL
        };

        enum channelopt {
            kLctopK0s = BIT(9),
            kLctopiL  = BIT(10)
        };

        enum massV0opt {
            kDeltaMassV0,
            kMassV0
        };

        AliHFMLVarHandlerNonPromptLc2V0bachelor();
        AliHFMLVarHandlerNonPromptLc2V0bachelor(int PIDopt, int massopt=kDeltaMassV0);
        virtual ~AliHFMLVarHandlerNonPromptLc2V0bachelor();

        AliHFMLVarHandlerNonPromptLc2V0bachelor(const AliHFMLVarHandlerNonPromptLc2V0bachelor &source) = delete;
        AliHFMLVarHandlerNonPromptLc2V0bachelor& operator=(const AliHFMLVarHandlerNonPromptLc2V0bachelor &source) = delete;

        virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
        virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int massHypo = 0, AliAODPidHF *pidrespo = nullptr);

        void SetMassV0Option(int opt) {fMassV0Opt = opt;}
        void SetUseKFParticleReco(bool useKF=true) {fUseKFParticle=useKF;}

    private:
        float fImpParProng[knMaxProngs] = {-999., -999., -999., -999.};  /// prong impact parameter
        float fNormd0MeasMinusExp = -999.;                               /// candidate topomatic variable
        float fMassK0s = -999.;                                          /// mass of the pipi in case of pK0s hypothesis
        float fMassL = -999.;                                            /// mass of the ppi in case of piL hypothesis
        float fImpParV0 = -999.;                                         /// impact parameter of V0
        float fDecayLengthV0 = -999.;                                    /// decay length of V0
        float fDCAV0 = -999.;                                            /// DCA V0 prongs
        float fPtV0 = -999.;                                             /// pT of V0
        float fcTauK0s = -999.;                                          /// cTau of K0s
        float fcTauL = -999.;                                            /// cTau of L
        float fCosPV0 = -999.;                                           /// cosine of pointing angle of V0
        float fsignd0 = -999.;                                           /// signed d0 proton (different from standard d0)
        float fArmqTOverAlphaV0 = -999.;                                 /// Armenteros qT/|alpha| of the V0
        int fMassV0Opt = kDeltaMassV0;                                   /// option for massV0 variables (mass or delta mass wrt K0s/L)
        bool fUseKFParticle = false;                                     /// flag to enable/disable recalculation of Lc decay vertex with KFParticle package
        float fKFTopoChi2 = -999.;                                       /// KF chi2 to primary vertex with topological constraint

        /// \cond CLASSIMP
        ClassDef(AliHFMLVarHandlerNonPromptLc2V0bachelor, 1); /// 
        /// \endcond
};
#endif
