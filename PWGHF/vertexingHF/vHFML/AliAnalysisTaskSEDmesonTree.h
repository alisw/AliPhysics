#ifndef ALIANALYSISTASKSEDMESONTREE_H
#define ALIANALYSISTASKSEDMESONTREE_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEDmesonTree
// \brief Analysis task to produce trees of D-meson candidates for ML analyses
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
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCHeader.h"
#include "AliNormalizationCounter.h"

class AliAnalysisTaskSEDmesonTree : public AliAnalysisTaskSE
{
public:

    enum
    {
        kD0toKpi = 0,
        kDplustoKpipi,
        kDstartoD0pi
    };

    AliAnalysisTaskSEDmesonTree();
    AliAnalysisTaskSEDmesonTree(const char *name, int fDecChannel, AliRDHFCuts *analysiscuts, bool createMLtree);
    virtual ~AliAnalysisTaskSEDmesonTree();

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

    void SetDecayChannel(bool dec = kD0toKpi)
    {
        fDecChannel = dec;
        switch(fDecChannel)
        {
            case kD0toKpi:
                fPdgD = 421;
                break;
            case kDplustoKpipi:
                fPdgD = 411;
                break;
            case kDstartoD0pi:
                fPdgD = 413;
                break;
        }
    }

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

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:
    enum
    {
        knVarForSparseAcc   = 2,
        knVarForSparseAccFD = 3,
        kMaxPtBins          = 100
    };

    AliAnalysisTaskSEDmesonTree(const AliAnalysisTaskSEDmesonTree &source);
    AliAnalysisTaskSEDmesonTree &operator=(const AliAnalysisTaskSEDmesonTree &source);

    int IsCandidateSelected(AliAODRecoDecayHF *&dMeson, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origownvtx);
    void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);
    bool CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau);
    void CreateEffSparses();

    AliAODEvent* fAOD = nullptr;                           /// AOD event

    TList *fOutput = nullptr;                              //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                          //!<! hist. for No. of events
    AliNormalizationCounter *fCounter = nullptr;           //!<! Counter for normalization
    THnSparseF *fnSparseMC[2] = {nullptr, nullptr};        //!<! THnSparse for MC
                                                           ///[0]: Acc step prompt D
                                                           ///[1]: Acc step FD D
    AliHFMLVarHandler *fMLhandler = nullptr;               //!<! object to handle ML tree creation and filling
    TTree *fMLtree = nullptr;                              //!<! tree with candidates for ML

    int fDecChannel = kD0toKpi;                            /// channel to analyse
    int fPdgD = 421;                                       /// pdg code of the D meson
    bool fReadMC = false;                                  /// flag for access to MC
    bool  fFillAcceptanceLevel = true;                     /// flag for filling true reconstructed D at acceptance level (see FillMCGenAccHistos)
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
    bool fAddNtrkl = false;                                /// flag to add number of tracklets in the tree
    bool fAddCentr = false;                                /// flag to add centrality percentile in the tree
    std::string fCentEstimator = "V0M";                    /// centrality estimator for tree

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEDmesonTree, 2); /// AliAnalysisTaskSE for production of D-meson trees
                                               /// \endcond
};

#endif
