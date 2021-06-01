#ifndef ALIANALYSISTASKSECHECKCHARMHADRONBKG_H
#define ALIANALYSISTASKSECHECKCHARMHADRONBKG_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSECheckCharmHadronBkg
// \brief Analysis task to study correlated backgrounds in D2H analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliHFMLResponse.h"

class AliAnalysisTaskSECheckCharmHadronBkg : public AliAnalysisTaskSE
{
public:
    enum
    {
        kUncorrelated,
        kDplusToKpipi,
        kDplusToKKpi,
        kDplusTopipipi,
        kDplusToKpie,
        kDplusToKpimu,
        kDzeroTopipipiX,
        kDzeroToKpipiX,
        kDstarToKpipi,
        kDsToKKpi,
        kDsToKminuspipluspiplus,
        kDsToKpluspipluspiminus,
        kBzeroToDzeroX,
        kBzeroToDplusX,
        kBzeroToDsX,
        kBzeroToLambdacX,
        kBzeroToJpsiX,
        kBplusToDzeroX,
        kBplusToDplusX,
        kBplusToDsX,
        kBplusToLambdacX,
        kBplusToJpsiX,
        kBsToDsX,
        kLambdacTopKpi,
        kLambdacToppipi
    };

    AliAnalysisTaskSECheckCharmHadronBkg();
    AliAnalysisTaskSECheckCharmHadronBkg(const char *name, AliRDHFCutsDplustoKpipi *analysiscuts);
    virtual ~AliAnalysisTaskSECheckCharmHadronBkg();

    void SetMassLimits(double lowlimit, double uplimit)
    {
        fLowMassLimit = lowlimit;
        fUpMassLimit = uplimit;
    }
    void SetMassBinWidth(double binwidth) { fMassBinWidth = binwidth; }
    void SetPtLimits(double lowlimit, double uplimit)
    {
        fPtMin = lowlimit;
        fPtMax = uplimit;
    }
    void SetPtBinWidth(double binwidth) { fPtBinWidth = binwidth; }
    void SetAODMismatchProtection(int opt = 1) { fAODProtection = opt; }

    void SetDoMLApplication(bool flag = true) {fApplyML = flag;}
    void SetMLConfigFile(TString path = "") { fConfigPath = path; }

    /// Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

private:
    AliAnalysisTaskSECheckCharmHadronBkg(const AliAnalysisTaskSECheckCharmHadronBkg &source);
    AliAnalysisTaskSECheckCharmHadronBkg &operator=(const AliAnalysisTaskSECheckCharmHadronBkg &source);
    void CreateSparse();
    int CheckMother(TClonesArray *mcArray, int daulab[3], int absdauPDG[3]);
    bool isPDGcodeAcceptable(int pdg);
    bool MatchToSpecie(int dauPDG[3], int pdg1, int pdg2, int pdg3, bool absPDGcode = true);
    bool IsInDecayChain(int abspdg, int motherlab, TClonesArray *mcArray);

    AliAODEvent *fAOD = nullptr;                               //!<! AOD event
    TList *fOutput = nullptr;                                  //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                              //!<! hist. for No. of events
    TH1F *fHistNCandidates = nullptr;                          //!<! hist. for No. of candidates
    THnSparseF *fMassVsPtVsOriginSparse = nullptr;             //!<! hist. for inv mass (topol+PID cuts)
    TH2F *fHistCentrality[3] = {nullptr, nullptr, nullptr};    //!<! hist. for cent distr (all,sel ev, )
    double fUpMassLimit = 1.6;                                 /// upper inv mass limit for histo (MeV/c^2)
    double fLowMassLimit = 2.2;                                /// lower inv mass limit for histo (MeV/c^2)
    double fPtMin = 0.;                                        /// upper pT limit for histo (GeV/c)
    double fPtMax = 50.;                                       /// lower pT limit for histo (GeV/c)
    int fNPtBins = 100;                                        /// Number of pT Bins
    int fNMassBins = 300;                                      /// Number of mass bins
    double fPtBinWidth = 0.5;                                  /// pT bin width (GeV/c)
    double fMassBinWidth = 0.002;                              /// mass bin width (MeV/c^2)
    TList *fListCuts = nullptr;                                /// list of cuts
    AliRDHFCutsDplustoKpipi *fRDCutsAnalysis = nullptr;        /// Cuts for Analysis
    int fAODProtection = 0;                                    /// flag to activate protection against AOD-dAOD mismatch.
                                                               /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    bool fUseBit = true;                                       /// flag to use bitmask
    bool fApplyML = false;                                     /// flag to apply ML models for charm-hadron selection
    TString fConfigPath = "";                                  /// path to ML config file
    AliHFMLResponse* fMLResponse = nullptr;                    //!<! object to handle ML response

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSECheckCharmHadronBkg, 1); /// AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
                                                       /// \endcond
};

#endif
