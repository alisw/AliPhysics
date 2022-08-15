#ifndef ALIANALYSISTASKSEHFRESONANCEBUILDER_H
#define ALIANALYSISTASKSEHFRESONANCEBUILDER_H

/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEHFResonanceBuilder
// \brief Analysis task to produce trees of D-meson candidates combined with other particles
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include <TNtuple.h>

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
#include "AliHFMLResponse.h"

class AliAnalysisTaskSEHFResonanceBuilder : public AliAnalysisTaskSE
{
public:

    enum charmHadron
    {
        kD0toKpi = 0,
        kDplustoKpipi,
        kDstartoD0pi
    };

    enum bachelorID
    {
        kPion = 0,
        kKaon,
        kProton,
        kDeuteron,
        kNumBachIDs
    };

    AliAnalysisTaskSEHFResonanceBuilder();
    AliAnalysisTaskSEHFResonanceBuilder(const char *name, int decChannel, AliRDHFCuts *analysiscuts);
    virtual ~AliAnalysisTaskSEHFResonanceBuilder();

    void SetReadMC(bool readMC = true)                                                            {fReadMC = readMC;}
    void SetAODMismatchProtection(int opt = 0)                                                    {fAODProtection = opt;}
    void SetAnalysisCuts(AliRDHFCuts *cuts)                                                       {fRDCuts = cuts;}

    void SetDecayChannel(int dec = kDplustoKpipi)
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

    /// methods for ML application
    void SetDoMLApplication(bool flag = true, bool isMultiClass = false)                          {fApplyML = flag; fMultiClass = isMultiClass;}
    void SetMLConfigFile(std::string path = "")                                                   {fConfigPath = path;}
    void SetIsDependentOnMLSelector(bool flag=true, std::string name="MLSelector") {
        fDependOnMLSelector = flag;
        fMLSelectorName = name;
    }

    /// methods for bachelor selection
    void SetBachelorFB(int filterBit = 4)                                                         {fFilterBitBachelor = filterBit;}
    void SetPtBachelorSelection(double pt = 0.05)                                                 {fPtTrackMin = pt;}
    void SetNsigmaBachelorSelection(
        float nSigmaTPCPi = 3.,
        float nSigmaTPCKa = 0., 
        float nSigmaTPCPr = 0., 
        float nSigmaTPCDe = 0., 
        float nSigmaTOFPi = 3., 
        float nSigmaTOFKa = 0., 
        float nSigmaTOFPr = 0.,
        float nSigmaTOFDe = 0.
    ) {
        fNsigmaBachelorTPC = {nSigmaTPCPi, nSigmaTPCKa, nSigmaTPCPr, nSigmaTPCDe};
        fNsigmaBachelorTOF = {nSigmaTOFPi, nSigmaTOFKa, nSigmaTOFPr, nSigmaTOFDe};
    }

    /// methods for charm resonance selection
    void SetCharmResoMassWindows(
        std::vector<float> minMassPi,
        std::vector<float> maxMassPi,
        std::vector<float> minMassKa,
        std::vector<float> maxMassKa,
        std::vector<float> minMassPr,
        std::vector<float> maxMassPr,
        std::vector<float> minMassDe,
        std::vector<float> maxMassDe
    ) {
        fInvMassResoPiMin = minMassPi; fInvMassResoPiMax = maxMassPi;
        fInvMassResoKaMin = minMassKa; fInvMassResoKaMax = maxMassKa;
        fInvMassResoPrMin = minMassPr; fInvMassResoPrMax = maxMassPr;
        fInvMassResoDeMin = minMassDe; fInvMassResoDeMax = maxMassDe;
    }

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:

    AliAnalysisTaskSEHFResonanceBuilder(const AliAnalysisTaskSEHFResonanceBuilder &source);
    AliAnalysisTaskSEHFResonanceBuilder &operator=(const AliAnalysisTaskSEHFResonanceBuilder &source);

    int IsCandidateSelected(AliAODRecoDecayHF *&dMeson, AliAODRecoDecayHF *&dMesonWithVtx, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origownvtx, AliAODPidHF *&pidHF, std::size_t &iCand);
    int IsBachelorSelected(AliAODTrack *&track, AliAODPidHF *&pidHF);
    bool IsInvMassResoSelected(double &mass, int &bachId);

    std::array<int, kNumBachIDs> kPdgBachIDs = {211, 321, 2212, 1000010020};

    AliAODEvent* fAOD = nullptr;                                                          /// AOD event

    TList *fOutput = nullptr;                                                             //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                                                         //!<! hist. for No. of events
    std::array<TH2F*, kNumBachIDs> fHistNsigmaTPCSelBach{};                               //!<! array of histograms with NsigmaTPC vs. p for selected bachelor tracks
    std::array<TH2F*, kNumBachIDs> fHistNsigmaTOFSelBach{};                               //!<! array of histograms with NsigmaTOF vs. p for selected bachelor tracks
    std::array<TH1F*, 3> fHistBDTOutputScore{};                                           //!<! array of histograms with BDT output scores for D mesons
    TNtuple *fNtupleCharmReso = nullptr;                                                  //!<! ntuple for HF resonances
    AliNormalizationCounter *fCounter = nullptr;                                          //!<! Counter for normalization

    int fDecChannel = kDplustoKpipi;                                                      /// channel to analyse
    int fPdgD = 411;                                                                      /// pdg code of the D meson
    bool fReadMC = false;                                                                 /// flag for access to MC
    bool  fFillAcceptanceLevel = true;                                                    /// flag for filling true reconstructed D at acceptance level (see FillMCGenAccHistos)
    int fAODProtection = 0;                                                               /// flag to activate protection against AOD-dAOD mismatch.
                                                                                          /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    TList *fListCuts = nullptr;                                                           /// list of cuts
    AliRDHFCuts *fRDCuts = nullptr;                                                       /// Cuts for Analysis

    // ML application
    bool fApplyML = false;                                                                /// flag to enable ML application
    bool fMultiClass = false;                                                             /// flag to enable multi-class models (Bkg, Prompt, FD)
    std::string fConfigPath = "";                                                         /// path to ML config file
    AliHFMLResponse* fMLResponse = nullptr;                                               //!<! object to handle ML response
    bool fDependOnMLSelector = false;                                                     /// flag to read ML scores from a AliAnalysisTaskSECharmHadronMLSelector task
    std::vector<float> fPtLimsML{};                                                       /// pT bins in case application of ML model is done in MLSelector task
    std::vector<std::vector<double> > fMLScoreCuts{};                                     /// score cuts used in case application of ML model is done in MLSelector task
    std::vector<std::vector<std::string> > fMLOptScoreCuts{};                             /// score cut options (lower, upper) used in case application of ML model is done in MLSelector task
    std::string fMLSelectorName = "MLSelector";                                           /// name of MLSelector task
    std::vector<std::vector<double> > fScoresFromMLSelector{};                            /// scores from MLSelector task
    std::vector<std::vector<double> > fScoresFromMLSelectorSecond{};                      /// scores from MLSelector task for second mass hypothesis

    // bachelor selection
    int fFilterBitBachelor{4};                                                            /// filter bit for bachelor track
    std::array<float, kNumBachIDs> fNsigmaBachelorTPC = {0.f, 0.f, 0.f, 0.f};   /// Nsigma cuts for bachelor track in TPC
    std::array<float, kNumBachIDs> fNsigmaBachelorTOF = {0.f, 0.f, 0.f, 0.f};   /// Nsigma cuts for bachelor track in TOF
    float fPtTrackMin{0.05};                                                              /// minimum pT for bachelor track

    // resonance selection
    std::vector<float> fInvMassResoPiMin{2.0};                                            /// minimum invariant mass values for HF resonance (in case of pion combination)
    std::vector<float> fInvMassResoPiMax{2.5};                                            /// maximum invariant mass values for HF resonance (in case of pion combination)
    std::vector<float> fInvMassResoKaMin{2.2};                                            /// minimum invariant mass values for HF resonance (in case of kaon combination)
    std::vector<float> fInvMassResoKaMax{2.6};                                            /// maximum invariant mass values for HF resonance (in case of kaon combination)
    std::vector<float> fInvMassResoPrMin{2.7};                                            /// minimum invariant mass values for HF resonance (in case of proton combination)
    std::vector<float> fInvMassResoPrMax{3.3};                                            /// maximum invariant mass values for HF resonance (in case of proton combination)
    std::vector<float> fInvMassResoDeMin{3.6};                                            /// maximum invariant mass values for HF resonance (in case of deuteron combination)
    std::vector<float> fInvMassResoDeMax{5.0};                                            /// minimum invariant mass values for HF resonance (in case of deuteron combination)

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEHFResonanceBuilder, 3); /// AliAnalysisTaskSE for production of HF resonance trees
                                               /// \endcond
};

#endif
