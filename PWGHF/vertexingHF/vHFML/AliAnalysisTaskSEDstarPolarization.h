#ifndef AliAnalysisTaskSEDstarPolarization_H
#define AliAnalysisTaskSEDstarPolarization_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEDstarPolarization
// \brief Analysis task to perform D*+ polarization analysis
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// S. Kundu, sourav.kundu@cern.ch
/////////////////////////////////////////////////////////////

#include <THnSparse.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

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
#include "AliHFMLResponseDstartoD0pi.h"
#include "AliHFMLResponseD0toKpi.h"

class AliAnalysisTaskSEDstarPolarization : public AliAnalysisTaskSE
{
public:

    enum
    {
        kDstartoD0pi = 0,
        kD0toKpi,
    };

    AliAnalysisTaskSEDstarPolarization();
    AliAnalysisTaskSEDstarPolarization(const char *name, AliRDHFCuts *analysiscuts);
    virtual ~AliAnalysisTaskSEDstarPolarization();

    void SetReadMC(bool readMC = true)                                                            {fReadMC = readMC;}
    void SetAODMismatchProtection(int opt = 0)                                                    {fAODProtection = opt;}
    void SetAnalysisCuts(AliRDHFCuts *cuts)                                                       {fRDCuts = cuts;}
    void SetUseFinePtBinsForSparse(bool useFineBins = true)                                       {fUseFinPtBinsForSparse = useFineBins;}
    void SetFillNSparseAcceptanceLevel(bool fill = true)                                          {fFillAcceptanceLevel = fill;}

    /// methods for ML application
    void SetDoMLApplication(bool flag = true)                                                     {fApplyML = flag;}
    void SetMLConfigFile(std::string path = "")                                                   {fConfigPath = path;}

    void SetIsDependentOnMLSelector(bool flag=true, std::string name="MLSelector") {
        fDependOnMLSelector = flag;
        fMLSelectorName = name;
    }

    void SetCheckWithD0()                                                                         {fDecChannel = kD0toKpi;}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:
    enum
    {
        knVarForSparseAcc    = 7,
        knVarForSparseReco   = 8,
    };

    AliAnalysisTaskSEDstarPolarization(const AliAnalysisTaskSEDstarPolarization &source);
    AliAnalysisTaskSEDstarPolarization &operator=(const AliAnalysisTaskSEDstarPolarization &source);

    int IsCandidateSelected(AliAODRecoDecayHF *&d, AliAODRecoDecayHF2Prong *&dZeroDau, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origownvtx,
                            std::vector<double> scoresFromMLSelector, std::vector<double> scoresFromMLSelectorSecond);
    void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, double centrality);
    bool CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau);
    void CreateEffSparses();
    void CreateRecoSparses();

    AliAODEvent* fAOD = nullptr;                                                    /// AOD event

    int fDecChannel = kDstartoD0pi;                                                 /// decay channel to analyse (main D*, D0 for analysis validation)

    TList *fOutput = nullptr;                                                       //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                                                   //!<! hist. for No. of events
    THnSparseF *fnSparseMC[2] = {nullptr, nullptr};                                 //!<! THnSparse for MC
    THnSparseF *fnSparseMCThetaPhiStar[2] = {nullptr, nullptr};                     //!<! THnSparse for MC
                                                                                    ///[0]: Acc step prompt D
                                                                                    ///[1]: Acc step FD D

    bool fReadMC = false;                                                           /// flag for access to MC
    bool  fFillAcceptanceLevel = true;                                              /// flag for filling true reconstructed D at acceptance level (see FillMCGenAccHistos)
    int fAODProtection = 0;                                                         /// flag to activate protection against AOD-dAOD mismatch.
                                                                                    /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    TList *fListCuts = nullptr;                                                     /// list of cuts
    AliRDHFCuts *fRDCuts = nullptr;                                                 /// Cuts for Analysis
    bool fUseFinPtBinsForSparse = true;                                             /// flag to fill pt axis of sparse with 0.1 GeV/c wide bins
                    
    // ML tree application
    THnSparseF* fnSparseReco[4] = {nullptr, nullptr, nullptr, nullptr};             //!<! THnSparse for reco candidates
    THnSparseF* fnSparseRecoThetaPhiStar[4] = {nullptr, nullptr, nullptr, nullptr}; //!<! THnSparse for reco candidates

    bool fApplyML = false;                                                          /// flag to enable ML application
    std::string fConfigPath = "";                                                   /// path to ML config file
    AliHFMLResponse* fMLResponse = nullptr;                                         //!<! object to handle ML response

    bool fDependOnMLSelector = false;                                               /// flag to read ML scores from a AliAnalysisTaskSECharmHadronMLSelector task
    std::vector<float> fPtLimsML{};                                                 /// pT bins in case application of ML model is done in MLSelector task   
    std::vector<std::vector<double> > fMLScoreCuts{};                               /// score cuts used in case application of ML model is done in MLSelector task   
    std::vector<std::vector<std::string> > fMLOptScoreCuts{};                       /// score cut options (lower, upper) used in case application of ML model is done in MLSelector task                                           
    std::string fMLSelectorName = "MLSelector";                                     /// name of MLSelector task

    ROOT::Math::PxPyPzMVector fourVecDstar{};                                       /// four vector for reconstructed D* in the lab
    ROOT::Math::PxPyPzMVector fourVecD0{};                                          /// four vector for reconstructed D0 in the lab
    ROOT::Math::PxPyPzMVector fourVecPi{};                                          /// four vector for reconstructed pion in the lab
    ROOT::Math::PxPyPzMVector fourVecPiCM{};                                        /// four vector for reconstructed pion in the D* RF

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEDstarPolarization, 5); /// AliAnalysisTaskSE for production of D-meson trees
                                               /// \endcond
};

#endif
