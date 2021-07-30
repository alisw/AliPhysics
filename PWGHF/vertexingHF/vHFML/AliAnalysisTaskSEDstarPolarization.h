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

class AliAnalysisTaskSEDstarPolarization : public AliAnalysisTaskSE
{
public:

    AliAnalysisTaskSEDstarPolarization();
    AliAnalysisTaskSEDstarPolarization(const char *name, AliRDHFCuts *analysiscuts);
    virtual ~AliAnalysisTaskSEDstarPolarization();

    void SetReadMC(bool readMC = true)                                                            {fReadMC = readMC;}
    void SetAODMismatchProtection(int opt = 0)                                                    {fAODProtection = opt;}
    void SetAnalysisCuts(AliRDHFCuts *cuts)                                                       {fRDCuts = cuts;}
    void SetUseFinePtBinsForSparse(bool useFineBins = true)                                       {fUseFinPtBinsForSparse = useFineBins;}
    void SetFillNSparseAcceptanceLevel(bool fill = true)                                          {fFillAcceptanceLevel = fill;}

    /// methods for ML application
    void SetDoMLApplication(bool flag = true, bool isMultiClass = false)                          {fApplyML = flag; fMultiClass = isMultiClass;}
    void SetMLConfigFile(std::string path = "")                                                   {fConfigPath = path;}
    void SetMLBinsForSparse(int nbins = 300, double min = 0.85, double max = 1.)                  {fNMLBins[0] = nbins; fMLOutputMin[0] = min; fMLOutputMax[0] = max;}
    void SetMultiClassMLBinsForSparse(int nbinsBkg = 100, int nbinsPrompt = 100, int nbinsFD = 100,
                                      double minBkg = 0., double maxBkg = 1.,
                                      double minPrompt = 0., double maxPrompt = 1.,
                                      double minFD = 0., double maxFD = 1.)
    {
        fNMLBins[0] = nbinsBkg; fNMLBins[1] = nbinsPrompt; fNMLBins[2] = nbinsFD;
        fMLOutputMin[0] = minBkg; fMLOutputMin[1] = minPrompt; fMLOutputMin[2] = minFD;
        fMLOutputMax[0] = maxBkg; fMLOutputMax[1] = maxPrompt; fMLOutputMax[2] = maxFD;
    }

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:
    enum
    {
        knVarForSparseAcc    = 4,
        knVarForSparseReco   = 4,
    };

    AliAnalysisTaskSEDstarPolarization(const AliAnalysisTaskSEDstarPolarization &source);
    AliAnalysisTaskSEDstarPolarization &operator=(const AliAnalysisTaskSEDstarPolarization &source);

    int IsCandidateSelected(AliAODRecoCascadeHF *&dStar, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origownvtx);
    void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, double centrality);
    bool CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau);
    void CreateEffSparses();
    void CreateRecoSparses();

    AliAODEvent* fAOD = nullptr;                                                /// AOD event

    TList *fOutput = nullptr;                                                   //!<! list send on output slot 0
    TH1F *fHistNEvents = nullptr;                                               //!<! hist. for No. of events
    THnSparseF *fnSparseMC[2] = {nullptr, nullptr};                             //!<! THnSparse for MC
                                                                                ///[0]: Acc step prompt D
                                                                                ///[1]: Acc step FD D

    bool fReadMC = false;                                                       /// flag for access to MC
    bool  fFillAcceptanceLevel = true;                                          /// flag for filling true reconstructed D at acceptance level (see FillMCGenAccHistos)
    int fAODProtection = 0;                                                     /// flag to activate protection against AOD-dAOD mismatch.
                                                                                /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    TList *fListCuts = nullptr;                                                 /// list of cuts
    AliRDHFCuts *fRDCuts = nullptr;                                             /// Cuts for Analysis
    bool fUseFinPtBinsForSparse = true;                                         /// flag to fill pt axis of sparse with 0.1 GeV/c wide bins
                    
    // ML tree application
    THnSparseF* fnSparseReco[4] = {nullptr, nullptr, nullptr, nullptr};         //!<! THnSparse for reco candidates
    bool fApplyML = false;                                                      /// flag to enable ML application
    bool fMultiClass = false;                                                   /// flag to enable multi-class models (Bkg, Prompt, FD)
    std::string fConfigPath = "";                                               /// path to ML config file
    AliHFMLResponseDstartoD0pi* fMLResponse = nullptr;                          //!<! object to handle ML response
    int fNMLBins[3] = {1000, 100, 100};                                         /// number of bins for ML output axis in THnSparse
    double fMLOutputMin[3] = {0., 0., 0.};                                      /// min for ML output axis in THnSparse
    double fMLOutputMax[3] = {1., 1., 1.};                                      /// max for ML output axis in THnSparse

    ROOT::Math::PxPyPzMVector fourVecDstar{};                                   /// four vector for reconstructed D* in the lab
    ROOT::Math::PxPyPzMVector fourVecD0{};                                      /// four vector for reconstructed D0 in the lab
    ROOT::Math::PxPyPzMVector fourVecPi{};                                      /// four vector for reconstructed pion in the lab
    ROOT::Math::PxPyPzMVector fourVecD0CM{};                                    /// four vector for reconstructed D0 in the D* RF
    ROOT::Math::PxPyPzMVector fourVecPiCM{};                                    /// four vector for reconstructed pion in the D* RF

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEDstarPolarization, 1); /// AliAnalysisTaskSE for production of D-meson trees
                                               /// \endcond
};

#endif
