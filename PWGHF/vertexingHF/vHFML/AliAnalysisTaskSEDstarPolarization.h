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
#include "AliVertexerTracks.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisFilter.h"

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

    void SetReadMC(bool readMC = true)                                                                          {fReadMC = readMC;}
    void SetAODMismatchProtection(int opt = 0)                                                                  {fAODProtection = opt;}
    void SetAnalysisCuts(AliRDHFCuts *cuts)                                                                     {fRDCuts = cuts;}
    void SetUseFinePtBinsForSparse(bool useFineBins = true)                                                     {fUseFinPtBinsForSparse = useFineBins;}
    void SetFillNSparseAcceptanceLevel(bool fill = true)                                                        {fFillAcceptanceLevel = fill;}
    void SetRecomputeDstarCombinatorial(bool recompute = true)                                                  {fRecomputeDstarCombinatorial = recompute;}

    /// methods for ML application
    void SetDoMLApplication(bool flag = true)                                                                   {fApplyML = flag;}
    void SetMLConfigFile(std::string path = "")                                                                 {fConfigPath = path;}

    void SetIsDependentOnMLSelector(bool flag=true, std::string name="MLSelector") {
        fDependOnMLSelector = flag;
        fMLSelectorName = name;
    }

    void SetCheckWithD0()                                                                                       {fDecChannel = kD0toKpi;}
    void SetFillBkgSparse()                                                                                     {fFillBkgSparse = true;}

    void SetBinsForMLAxes(std::array<int, 3> nBins, std::array<double, 3> mins, std::array<double, 3> maxs) {
        fNBinsML = nBins;
        fMLOutputMin = mins;
        fMLOutputMax = maxs;
    }

    void EnableTrackCutVariation(bool flag = true)                                                              {fApplyTrackCutVariations = flag;}

    // method to set Qn-vector calculation
    void SetComputeQnVector(bool flag = true)                                                                   {fComputeQnVectors = flag;}
    void SetQnVecTaskName(std::string name)                                                                     {fTenderTaskName = name;}
    void SetQnCalibFileName(std::string name)                                                                   {fQnCalibFileName = name;}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

private:
    enum
    {
        knVarForSparseAcc    = 11,
        knVarForSparseReco   = 19,
    };

    AliAnalysisTaskSEDstarPolarization(const AliAnalysisTaskSEDstarPolarization &source);
    AliAnalysisTaskSEDstarPolarization &operator=(const AliAnalysisTaskSEDstarPolarization &source);

    AliAODRecoCascadeHF *MakeCascade(AliAODRecoDecayHF2Prong* trackD0, AliAODTrack *track, AliESDtrack *esdTrackPi, AliESDVertex *fV1);
    bool RecomputeDstarCombinatorial(TClonesArray *array2Prongs, TClonesArray *arrayTracks, std::vector<AliAODRecoCascadeHF*> &arrayDstar);
    bool SelectInvMassAndPtDstarD0pi(double *px, double *py, double *pz);

    int IsCandidateSelected(AliAODRecoDecayHF *&d, AliAODRecoDecayHF2Prong *&dZeroDau, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origownvtx,
                            std::vector<double> scoresFromMLSelector, std::vector<double> scoresFromMLSelectorSecond, std::vector<double> &scores, std::vector<double> &scoresSecond, int trackCutFlags[4]);
    void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, double centrality);
    bool CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau);
    void CreateEffSparses();
    void CreateRecoSparses();
    double GetPhiInRange(double phi);
    double GetDeltaPsiSubInRange(double psi1, double psi2);

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

    TH1F* fHistEvPlane[3] = {nullptr, nullptr, nullptr};                            //!<! Histograms with event plane angle
    TH2F* fHistEvPlaneResol[3] = {nullptr, nullptr, nullptr};                        //!<! Histograms for event plane resolution estimation

    bool fApplyML = false;                                                          /// flag to enable ML application
    std::string fConfigPath = "";                                                   /// path to ML config file
    AliHFMLResponse* fMLResponse = nullptr;                                         //!<! object to handle ML response

    bool fDependOnMLSelector = false;                                               /// flag to read ML scores from a AliAnalysisTaskSECharmHadronMLSelector task
    std::vector<float> fPtLimsML{};                                                 /// pT bins in case application of ML model is done in MLSelector task   
    std::vector<std::vector<double> > fMLScoreCuts{};                               /// score cuts used in case application of ML model is done in MLSelector task   
    std::vector<std::vector<std::string> > fMLOptScoreCuts{};                       /// score cut options (lower, upper) used in case application of ML model is done in MLSelector task                                           
    std::string fMLSelectorName = "MLSelector";                                     /// name of MLSelector task

    std::array<int, 3> fNBinsML = {100, 100, 100};                                  /// number of bins for the ML output scores in the THnSparse
    std::array<double, 3> fMLOutputMin = {0., 0., 0.};                              /// minimum axis value for the ML output scores in the THnSparse
    std::array<double, 3> fMLOutputMax = {1., 1., 1.};                              /// maximum axis value for the ML output scores in the THnSparse

    ROOT::Math::PxPyPzMVector fourVecDstar{};                                       /// four vector for reconstructed D* in the lab
    ROOT::Math::PxPyPzMVector fourVecD0{};                                          /// four vector for reconstructed D0 in the lab
    ROOT::Math::PxPyPzMVector fourVecPi{};                                          /// four vector for reconstructed pion in the lab
    ROOT::Math::PxPyPzMVector fourVecPiCM{};                                        /// four vector for reconstructed pion in the D* RF

    bool fFillBkgSparse = false;                                                    /// flag to fill or not fill bkg sparses
    bool fRecomputeDstarCombinatorial = false;                                      /// flag to recompute D* combinatorial (D*->D0pi)
    double fBzkG = 0.;                                                              /// z componenent of field in kG

    bool fApplyTrackCutVariations = false;                                          /// flag to enable track-cut variation to be stored in ThnSparse
    AliRDHFCuts *fRDCutsTrackVariations[4] = {nullptr, nullptr, nullptr, nullptr};  /// Cuts for Analysis with alternative track selections

    AliESDtrackCuts *fEsdTrackCutsSoftPi = nullptr;                                 /// ESD track cuts for soft pion in case of combinatoric recomputation
    AliAnalysisFilter *fTrkFilterSoftPi = nullptr;                                  /// track filter for soft pion in case of combinatoric recomputation

    // event-plane calculation
    bool fComputeQnVectors = false;                                                 /// flag to enable computation of Qn-vectors
    std::string fTenderTaskName = "HFTenderQnVectors";                              /// name of tender task needed to get the calibrated Qn vectors
    std::string fQnCalibFileName = "";                                              /// AODB file name for calibrations (if Qn-framework not used)

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEDstarPolarization, 11); /// AliAnalysisTaskSE for production of D-meson trees
                                               /// \endcond
};

#endif
