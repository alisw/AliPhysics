#ifndef ALIANALYSISTASKSENONPROMPTLC_H
#define ALIANALYSISTASKSENONPROMPTLC_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSENonPromptLc
// \brief Analysis task to produce trees of Lc candidates for ML analyses of non-prompt Lc
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include <THnSparse.h>

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliHFMLVarHandler.h"
#include "AliNormalizationCounter.h"

class AliAnalysisTaskSENonPromptLc : public AliAnalysisTaskSE
{
public:

    enum
    {
        kLctopKpi = 0,
        kLctopK0s = 1,
        kLctopiL  = 2
    };

    AliAnalysisTaskSENonPromptLc();
    AliAnalysisTaskSENonPromptLc(const char *name, int fDecChannel, AliRDHFCuts *analysiscuts, bool createMLtree);
    virtual ~AliAnalysisTaskSENonPromptLc();

    void SetDecayChannel(bool dec = kLctopKpi)                                                    {fDecChannel = dec;}
    void SetReadMC(bool readMC = true)                                                            {fReadMC = readMC;}
    void SetAODMismatchProtection(int opt = 0)                                                    {fAODProtection = opt;}
    void SetAnalysisCuts(AliRDHFCuts *cuts)                                                       {fRDCuts = cuts;}
    void SetUseFinePtBinsForSparse(bool useFineBins = true)                                       {fUseFinPtBinsForSparse = useFineBins;}
    void SetFillNSparseAcceptanceLevel(bool fill = true)                                          {fFillAcceptanceLevel = fill;}
    /// methods for ML tree creation
    void SetCreateMLTree(bool flag = true)                                                        {fCreateMLtree = flag;}
    void SetMLTreePIDopt(int opt)                                                                 {fPIDopt = opt;} // default AliHFMLVarHandler::kNsigmaDetAndCombPID
    void SetMLTreeAddTrackVar(bool flag = true)                                                   {fAddSingleTrackVar = flag;}
    void SetKeepOnlyBkgFromHIJING(bool keeponlyhijing = true)                                     {fKeepOnlyBkgFromHIJING = keeponlyhijing;}
    void SetFillOnlySignalInMLtree(bool opt = true)
    {
        if (fReadMC)
            fFillOnlySignal = opt;
        else
        {
            if (opt)
                AliError("fReadMC has to be true");
        }
    }
    void EnableMLTreeEvtSampling(float fractokeep, unsigned long long seed, int option=0)
    {
        fEnableEvtSampling = true;
        fFracEvtToKeep = fractokeep;
        fSeedSampling = seed;
        fOptionSampling = option;
    }
    void EnableMLTreeCandSampling(float fractokeep, float maxptsampling)
    {
        fEnableCandSampling = true;
        fFracCandToKeep = fractokeep;
        fMaxCandPtSampling = maxptsampling;
    }
    void EnableKFRecoForV0bachelor(bool enable = true)                                            {fUseKFRecoForV0bachelor = enable;}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:
    enum
    {
        knVarForSparseAcc   = 3,
        knVarForSparseAccFD = 4,
        kMaxPtBins          = 100
    };

    AliAnalysisTaskSENonPromptLc(const AliAnalysisTaskSENonPromptLc &source);
    AliAnalysisTaskSENonPromptLc &operator=(const AliAnalysisTaskSENonPromptLc &source);

    int IsCandidateSelected(AliAODRecoDecayHF *&lc, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origownvtx);
    void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);
    bool CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau);
    void CreateEffSparses();

    AliAODEvent* fAOD = nullptr;                           /// AOD event

    TList *fOutput = nullptr;                              //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                          //!<! hist. for No. of events
    AliNormalizationCounter *fCounter = nullptr;           //!<! Counter for normalization
    THnSparseF *fnSparseMC[2] = {nullptr, nullptr};        //!<! THnSparse for MC
                                                           ///[0]: Acc step prompt Lc
                                                           ///[1]: Acc step FD Lc
    AliHFMLVarHandler *fMLhandler = nullptr;               //!<! object to handle ML tree creation and filling
    TTree *fMLtree = nullptr;                              //!<! tree with candidates for ML

    int fDecChannel = kLctopKpi;                           /// channel to analyse
    bool fReadMC = false;                                  /// flag for access to MC
    bool  fFillAcceptanceLevel = true;                     /// flag for filling true reconstructed Lc at acceptance level (see FillMCGenAccHistos)
    int fAODProtection = 0;                                /// flag to activate protection against AOD-dAOD mismatch.
                                                           /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    TList *fListCuts = nullptr;                            /// list of cuts
    AliRDHFCuts *fRDCuts = nullptr;                        /// Cuts for Analysis
    bool fUseFinPtBinsForSparse = true;                    /// flag to fill pt axis of sparse with 0.1 GeV/c wide bins
    bool fKeepOnlyBkgFromHIJING = false;                   /// flag to keep background from Hijing only

    bool fCreateMLtree = true;
    int fPIDopt = AliHFMLVarHandler::kNsigmaDetAndCombPID; /// option for PID variables
    bool fAddSingleTrackVar = true;                        /// option to store single track variables
    bool fFillOnlySignal = true;                           /// option to store only signal when using MC
    bool fEnableEvtSampling = true;                        /// flag to apply event sampling
    float fFracEvtToKeep = 1.1;                            /// fraction of events to be kept by event sampling
    unsigned long fSeedSampling = 0;                       /// seed for event sampling
    int fOptionSampling = 0;                               /// option for event sampling (0: keeps events with random < fracToKeep, 1: keep events with random > 1-fracToKeep)
    bool fEnableCandSampling = true;                       /// flag to apply candidate sampling
    float fFracCandToKeep = 1.1;                           /// fraction of candidates to be kept by sampling
    float fMaxCandPtSampling = 0.;                         /// maximun candidate pt to apply sampling
    bool fUseKFRecoForV0bachelor = true;                   /// flag to enable KFParticle reconstruction for Lc->V0bachelor in tree

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSENonPromptLc, 2); /// AliAnalysisTaskSE for non-prompt Lc tree
                                               /// \endcond
};

#endif
