#ifndef AliAnalysisV0Lam_cxx
#define AliAnalysisV0Lam_cxx
class TH1F;
class TH1D;
class TH2D;
class TH3D;
//class TProfile;
/* class AliESDEvent; */
class AliAODEvent;
class AliAODv0;
class AliAODMCParticle;
/* class AliESDtrackCuts; */
/* class AliESDpid; */
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisV0LamEventCollection.h"
#include "AliAODpidUtil.h"
/* #include "AliESDpid.h" */
#include "AliCentrality.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisV0LamCutProcessing.h"

// Author: Jai Salzwedel, jai.salzwedel@cern.ch

class AliAnalysisV0Lam : public AliAnalysisTaskSE {
  public:
    AliAnalysisV0Lam();
    AliAnalysisV0Lam(const char *name, Int_t varCutType, Int_t nomincalCutIndex, Bool_t flattenCent);
    virtual ~AliAnalysisV0Lam();
    virtual void UserCreateOutputObjects();
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *);
    enum 
    {
      nCentBins   = 20,
      zVertexBins = 10,
      nEventsToMix =  5,
    };

  
  private:
    AliAODEvent    *fAOD; //! AOD object
    TList          *fOutputList; //! Compact Output list
    AliAODpidUtil  *fpidAOD; //!
    AliAnalysisV0LamCutProcessing *fCutProcessor; //!
    
    // Variables that are constant, preset parameters (e.g. cut values)
    int    fMaxV0Mult;                    //maximum number of V0s, array size
    double fEtaDaughter;	          //maximum eta for daughter particles
    double fMassWindowLam;                //Accepted V0 mass allowance
    double fTOFLow;                       //boundary for TOF usage
    double fSigmaCutTOFPion;              //max NSigma allowed
    double fSigmaCutTPCPion;
    double fSigmaCutTOFProton;
    double fSigmaCutTPCProton; 
    double fPDGLambda, fPDGProton, fPDGPion; //true PDG masses
    int    fEventCount;
    int    fNumberOfVariableCutValues; //this is > 1 if checking systematics of a
    int    fNumberOfCfVariableCutValues; //Only different from above if doing variable avg sep cuts
    int    fDefaultCutType;             //DCA, CosP, Mass... which is being varied
    int    fNominalCutIndex;    //Index of nominal cut value
    int    fNumberVariableAvgSepCuts;
    bool   fIsUsingVariableAvgSepCut;
    bool   fIsMCEvent;
    bool   fFlattenCent;
    AliAnalysisV0LamEventCollection ***fEC; //!
    AliAnalysisV0LamEvent *fEvt; //!

    // Variables used for local debugging and output
    int fTotalLambda; //counts number of v0s found across all events
    int fTotalAntiLambda;
    int fV0Candidates;
    
    // Histograms created in the analysis
    TH2F *fTPCVsPPosLam; //!
    TH2F *fTPCVsPNegLam; //!
    TH2F *fTPCVsPPosALam; //!
    TH2F *fTPCVsPNegALam; //!
    TH1F *fMultDistRough; //!
    TH2F *fMultDistLambda; //!
    TH2F *fMultDistAntiLambda; //!
    TH3F *fMultCentLambda; //!
    TH3F *fMultCentAntiLambda; //!
    TH1F *fCentrality; //!
    TH2F *fRemainingFromBeginningToV0Finder; //!
    TH3F *fRemainingFromBeginningToRecon; //!
    TH3F *fRemainingFromV0FinderToRecon; //!
    TH1F *fMCTruthOfOriginalParticles; //!
    TH1F *fMCTruthOfV0FinderParticles; //!
    TH2F *fMCTruthOfReconstructedParticles; //!
    TH1F *fMCFakeParticleIdentity; //!
    TH1F *fMCOtherV0Identity; //!
    //V0 Shared daughter culling statistics (mostly for MC, not working)
    /* TH1F *fDataCompeted; //!  */
    /* TH1F *fDataCulled; //! */
    /* TH1F *fRemainingV0s; //! */
    /* TH1F *fRemainingFrac; //! */
    
    //Momentum resolution (pre-)correction analysis
    TH1F *fHistPsmearingKreconMinusKtruthLL; //!
    TH2F *fHistPsmearingKreconVsKtruthLL; //!
    TH1F *fHistPsmearingKreconMinusKtruthAA; //!
    TH2F *fHistPsmearingKreconVsKtruthAA; //!
    TH1F *fHistPsmearingKreconMinusKtruthLA; //!
    TH2F *fHistPsmearingKreconVsKtruthLA; //!

    // New momentum smearing analysis
    /* vector<TH2F *> fHistsSignalMomResTruthLL; */
    /* vector<TH2F *> fHistsSignalMomResRecLL; */
    /* vector<TH2F *> fHistsBkgMomResTruthLL; */
    /* vector<TH2F *> fHistsBkgMomResRecLL; */
    /* vector<TH2F *> fHistsSignalMomResTruthAA; */
    /* vector<TH2F *> fHistsSignalMomResRecAA; */
    /* vector<TH2F *> fHistsBkgMomResTruthAA; */
    /* vector<TH2F *> fHistsBkgMomResRecAA; */
    /* vector<TH2F *> fHistsSignalMomResTruthLA; */
    /* vector<TH2F *> fHistsSignalMomResRecLA; */
    /* vector<TH2F *> fHistsBkgMomResTruthLA; */
    /* vector<TH2F *> fHistsBkgMomResRecLA; */
    /* vector<TH1D *> fLednickyWeightsLA; */
    /* vector<TH1D *> fLednickyWeightsLL; */
    /* TH2F *fSignalMomResTruthUnweightedLL; //! */
    /* TH2F *fSignalMomResTruthUnweightedAA; //! */
    /* TH2F *fSignalMomResTruthUnweightedLA; //! */
    /* TH2F *fBkgMomResTruthUnweightedLL; //! */
    /* TH2F *fBkgMomResTruthUnweightedAA; //! */
    /* TH2F *fBkgMomResTruthUnweightedLA; //! */

    /* TH3F *fSignalMomResTruthVsPtLL; //! */
    /* TH3F *fSignalMomResTruthVsPtAA; //! */
    /* TH3F *fSignalMomResTruthVsPtLA; //! */
    /* TH3F *fBkgMomResTruthVsPtLL; //! */
    /* TH3F *fBkgMomResTruthVsPtAA; //! */
    /* TH3F *fBkgMomResTruthVsPtLA; //! */

    // MC info
    TH1F *fMCTruthPtLam; //!
    TH1F *fMCTruthPtALam; //!
    TH1F *fMCTruthPhiLam; //!
    TH1F *fMCTruthPhiALam; //!
    TH1F *fMCTruthEtaLam; //!
    TH1F *fMCTruthEtaALam; //!
    
    //Pair kT Tracking
    TH3F *fKtLamLamSig; //!
    TH3F *fKtALamALamSig; //!
    TH3F *fKtLamALamSig; //!
    TH3F *fKtLamLamBkg; //!
    TH3F *fKtALamALamBkg; //!
    TH3F *fKtLamALamBkg; //!
    //Basic correlation functions
    TH3F *fSignalLamLam; //!
    TH3F *fBkgLamLam; //!
    TH3F *fSignalALamALam; //!
    TH3F *fBkgALamALam; //!
    TH3F *fSignalLamALam; //!
    TH3F *fBkgLamALam; //!
    //Daughter separation distance
    TH3F *fSignalLamLamProtSep; //!
    TH3F *fSignalLamLamPiMinusSep; //!
    TH3F *fSignalALamALamAntiProtSep; //!
    TH3F *fSignalALamALamPiPlusSep; //!
    TH3F *fSignalLamALamProtPiPlusSep; //!
    TH3F *fSignalLamALamAntiProtPiMinusSep; //!
    TH3F *fBkgLamLamProtSep; //!
    TH3F *fBkgLamLamPiMinusSep; //!
    TH3F *fBkgALamALamAntiProtSep; //!
    TH3F *fBkgALamALamPiPlusSep; //!
    TH3F *fBkgLamALamProtPiPlusSep; //!
    TH3F *fBkgLamALamAntiProtPiMinusSep; //!


    // opposite charged pair separation 
    TH3F *fSignalLamLamPlusMinusSep; //!
    TH3F *fSignalALamALamPlusMinusSep; //!
    TH3F *fSignalLamALamProtSep; //!
    TH3F *fSignalLamALamPionSep; //!

    TH3F *fBkgLamLamPlusMinusSep; //!
    TH3F *fBkgALamALamPlusMinusSep; //!
    TH3F *fBkgLamALamProtSep; //!
    TH3F *fBkgLamALamPionSep; //!

    //Functions
    void MyInit();
    vector<TVector3> GetGlobalPositionAtGlobalRadiiThroughTPC(const AliAODTrack *track, const Float_t bfield);
    Double_t GetAverageSeparation(vector<TVector3> &globalPositions1st, vector<TVector3> &globalPositions2nd);
    TVector3 GetEmissionPoint(const AliAODMCParticle * const track, TVector3 primVertex);
    void DoV0JudgmentCuts(const AliAnalysisV0LamEvent * const event, const int totalV0s);
    int DetermineWhichV0IsWorse(const AliAnalysisV0LamEvent * const event, const int V01, const int V02, const int Criterion, const int cutIndex);
    AliReconstructedV0::MCV0Origin_t DetermineV0Origin(AliAODv0 *v0, TClonesArray *mcArray);
    /* Double_t GetLednickyWeight(UInt_t histIndex, Double_t kstar, Int_t pairType); */
    int GetV0MCParticleID(AliAODv0 *v0, TClonesArray *mcArray);
    AliReconstructedV0::MCV0Origin_t DeterminePdgCodeOfMcParticle(AliAODMCParticle *mcParticle, TClonesArray *mcArray);
    void GetMCParticleMomentumTruth(TVector3 &pTruth, AliAODv0 *v0, TClonesArray *mcArray);
    void BinOriginInformationForMCParticles(TH1F *mcOriginalV0Hist, TH1F *mcV0FinderHist, TH2F *mcV0PassedCutsHist);
    void SetBinsOnOriginHists(TH1 *mcHist);
    //    void SetBinsOnOriginHists(TH2 *mcHist);
    void SetBinsOnOriginHists(TH3 *mcHist);
    TH1F *CreateLambdaOriginHist(TClonesArray *mcArray, Int_t nLastHijingLabel);
    void FillReconstructedV0MCOrigin(const AliReconstructedV0 * v0, TH2F *histPassedCutsOrigin);
    bool IsInjectedParticle(AliAODv0 *v0, TClonesArray *mcArray, Int_t nLastHijingLabel);
    bool IsCorrectEventTrigger();
    void AddV0ToMultiplicityCounts(AliReconstructedV0 *v0, vector<int> & lambdaCount, vector<int> & antiLambdaCount);
    void HistogramEventMultiplicities(vector<int> & lambdaCount, vector<int> & antiLambdaCount, int centralityBin);
    void FillTPCSignalHists(const AliReconstructedV0 *v0, const double posDaughterP, const double posDaughterTPCSignal, const double negDaughterP, const double negDaughterTPCSignal);
    void CheckForFakeV0s(const AliReconstructedV0 *v0, TH1F *mcFakeParticleIdentity, TH1F *mcOtherV0Identity, const AliReconstructedV0::MCV0Origin_t mcV0Origin);
    double CalculateKstar(TVector3 p1, TVector3 p2, double mass1, double mass2);
    void BinMomentumSmearing(TVector3 v01MomentumRecon,TVector3 v01MomentumTruth, TVector3 v02MomentumRecon,TVector3 v02MomentumTruth, int pairType);
    /* void DoMomResCorrelationWeighting(const AliReconstructedV0 &v01, const AliReconstructedV0 &v02, Bool_t isMixedEvent, Int_t centBin, Int_t pairType); */
    void DoPairStudies(const AliAnalysisV0LamEvent *event, const int centralityBin);
    bool RejectEventCentFlat(float MagField, float CentPercent);
    

    AliAnalysisV0Lam(const AliAnalysisV0Lam&); // not implemented
    AliAnalysisV0Lam& operator=(const AliAnalysisV0Lam&); // not implemented
    ClassDef(AliAnalysisV0Lam, 1); 
};

#endif
