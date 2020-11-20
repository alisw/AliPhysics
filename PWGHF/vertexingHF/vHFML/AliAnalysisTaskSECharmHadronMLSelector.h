#ifndef ALIANALYSISTASKSECHARMHADRONMLSELECTOR_H
#define ALIANALYSISTASKSECHARMHADRONMLSELECTOR_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSECharmHadronMLSelector
// \brief Analysis task to select charm-hadron candidates based on ML response for subsequent tasks
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include <vector>

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliHFMLResponse.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSECharmHadronMLSelector : public AliAnalysisTaskSE
{
public:

    enum
    {
        kDplustoKpipi = 0,
        kDstoKKpi     = 1
    };

    AliAnalysisTaskSECharmHadronMLSelector();
    AliAnalysisTaskSECharmHadronMLSelector(const char *name, int fDecChannel, AliRDHFCuts *analysiscuts);
    virtual ~AliAnalysisTaskSECharmHadronMLSelector();

    // Setters
    void SetDecayChannel(bool dec = kDplustoKpipi)                                                {fDecChannel = dec;}
    void SetAODMismatchProtection(int opt = 0)                                                    {fAODProtection = opt;}
    void SetTriggerInfo(TString trigClass, unsigned long long mask=0)                             {fTriggerClass = trigClass; fTriggerMask = mask;}
    void SetAnalysisCuts(AliRDHFCuts *cuts)                                                       {fRDCuts = cuts;}
    void SetMLConfigFile(TString path = "")                                                       {fConfigPath = path;}

    // Getters
    std::vector<int> GetSelectedCandidates() const                                                {return fChHadIdx;}
    std::vector<std::vector<double> > GetMLSCores() const                                         {return fMLScores;}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:

    AliAnalysisTaskSECharmHadronMLSelector(const AliAnalysisTaskSECharmHadronMLSelector &source);
    AliAnalysisTaskSECharmHadronMLSelector &operator=(const AliAnalysisTaskSECharmHadronMLSelector &source);

    int IsCandidateSelected(AliAODRecoDecayHF *&chHad, AliAnalysisVertexingHF *vHF, int absPdgMom, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> modelPred);

    AliAODEvent* fAOD = nullptr;                           /// AOD event

    TList *fOutput = nullptr;                              //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                          //!<! hist. for No. of events
    TH1F *fHistNallCand = nullptr;                         //!<! hist. for No. of all candidates
    TH1F *fHistNselCand = nullptr;                         //!<! hist. for No. of selected candidates
    TH1F *fHistBDTOutput[3] = {};                          //!<! hist. for distributions of BDT output scores (max 3)

    int fDecChannel = kDplustoKpipi;                       /// channel to analyse
    int fAODProtection = 0;                                /// flag to activate protection against AOD-dAOD mismatch.
                                                           /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    TString fTriggerClass = "";                            /// trigger class
    unsigned long long fTriggerMask = AliVEvent::kAny;     /// trigger mask
    TList *fListCuts = nullptr;                            /// list of cuts
    AliRDHFCuts *fRDCuts = nullptr;                        /// Cuts for Analysis
    TString fConfigPath = "";                              /// path to ML config file
    AliHFMLResponse* fMLResponse = nullptr;                //!<! object to handle ML response

    std::vector<int> fChHadIdx = {};                       /// vector with indexes of charm selected charm hadrons
    std::vector<std::vector<double> > fMLScores = {};      /// vector of vectors of ML output scores for each selected charm hadron

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSECharmHadronMLSelector, 1); /// AliAnalysisTaskSE for charm-hadron candidate selection with ML
                                                         /// \endcond
};

#endif
