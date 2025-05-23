#ifndef AliAnalysisTaskCVEPIDCMEDiff_cxx
#define AliAnalysisTaskCVEPIDCMEDiff_cxx
#include <TH3.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "TList.h"
#include "TProfile.h"
#include "TProfile3D.h"

enum class PairType {
    kLambdaProton = 0,
    kLambdaAntiProton,
    kAntiLambdaProton,
    kAntiLambdaAntiProton,
    kNumPairs = 4
};

class AliAnalysisTaskCVEPIDCMEDiff : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCVEPIDCMEDiff();
  AliAnalysisTaskCVEPIDCMEDiff(const char* name);
  virtual ~AliAnalysisTaskCVEPIDCMEDiff();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  // Switch
  void IfCalculateLambdaProton(bool b) { this->isCalculateLambdaProton = b; }
  void IfCalculateLambdaHadron(bool b) { this->isCalculateLambdaHadron = b; }
  void IfCalculateLambdaLambda(bool b) { this->isCalculateLambdaLambda = b; }
  void IfCalculateLambdaPion(bool b) { this->isCalculateLambdaPion = b; }
  void IfDoNUE(bool bDoNUE) { this->isDoNUE = bDoNUE; }
  void IfDoLambdaNUE(bool bDoLambdaNUE) { this->isDoLambdaNUE = bDoLambdaNUE; }
  void IfNarrowDcaCuts768(bool bNarrowDcaCuts768) { this->isNarrowDcaCuts768 = bNarrowDcaCuts768; }
  void IfProtonCustomizedDCACut(bool bProtonCustomizedDCACut) { this->isProtonCustomizedDCACut = bProtonCustomizedDCACut; }
  void IfTightPileUp(bool bTightPileUp) { this->isTightPileUp = bTightPileUp; }

  // read in
  void SetListForNUE(TList* flist) { this->fListNUE = (TList*)flist->Clone(); }
  void SetListForNUA(TList* flist) { this->fListNUA = (TList*)flist->Clone(); }
  void SetListForVZEROCalib(TList* flist) { this->fListVZEROCalib = (TList*)flist->Clone(); }

  // Global
  void SetPeriod(TString period) { this->fPeriod = period; }
  // Event
  void SetVzCut(float vzCut) { this->fVzCut = vzCut; }
  // Plane
  void SetPlaneEstimator(TString planeEstimator) { this->fPlaneEstimator = planeEstimator; }
  // Track
  void SetFilterBit(int filterBit) { this->fFilterBit = filterBit; }
  void SetNclsCut(int nclsCut) { this->fNclsCut = nclsCut; }
  void SetChi2Max(float chi2Max) { this->fChi2Max = chi2Max; }
  void SetChi2Min(float chi2Min) { this->fChi2Min = chi2Min; }
  void SetSpecialHadronDCAzMax(float hadronDCAzMax) { this->fSpecialHadronDCAzMax = hadronDCAzMax; }
  // Proton
  void SetSpecialProtonDCAzMax(float protonDCAzMax) { this->fSpecialProtonDCAzMax = protonDCAzMax; }
  // PID
  void SetNSigmaTPC(float nSigmaTPC) { this->fNSigmaTPC = nSigmaTPC; }
  void SetNSigmaRMS(float nSigmaRMS) { this->fNSigmaRMS = nSigmaRMS; }
  // V0
  void SetV0CPAMin(float v0CPAMin) { this->fV0CPAMin = v0CPAMin; }
  void SetV0DecayLengthMax(float v0DecayLengthMax) { this->fV0DecayLengthMax = v0DecayLengthMax; }
  void SetV0DecayLengthMin(float v0DecayLengthMin) { this->fV0DecayLengthMin = v0DecayLengthMin; }
  void SetV0DcaBetweenDaughtersMax(float v0DcaBetweenDaughtersMax) { this->fV0DcaBetweenDaughtersMax = v0DcaBetweenDaughtersMax; }
  // V0 Daughter Cut
  void SetDaughtersTPCNclsMin(float daughtersTPCNclsMin) { this->fDaughtersTPCNclsMin = daughtersTPCNclsMin; }
  void SetDaughtersDCAToPrimVtxMin(float daughtersDCAToPrimVtxMin) { this->fDaughtersDCAToPrimVtxMin = daughtersDCAToPrimVtxMin; }
  void SetRatioCrossedRowsFindable(float ratioCrossedRowsFindable) { this->fRatioCrossedRowsFindable = ratioCrossedRowsFindable; }
  void SetDaughtersNSigmaTPC(float daughtersNSigmaTPC) { this -> fDaughtersNSigmaTPC = daughtersNSigmaTPC; }

 private:
  bool fDebug;
  ////////////////////////
  // Procedural function
  ////////////////////////
  float GetTPCPlane();
  float GetV0CPlane(float centSPD1);
  void ResetVectors();
  bool LoopTracks();
  bool LoopV0s();
  bool PairV0Trk();
  bool PairV0V0();

  ////////////////////////
  // Functional function
  ////////////////////////
  // Read in
  bool LoadCalibHistForThisRun(); //deal with all the readin
  // Pile-up
  bool RejectEvtTFFit(float centSPD0);
  // Track
  bool AcceptAODTrack(AliAODTrack* track);
  bool CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck);
  float GetNUECor(int charge, float pt);
  float GetPIDNUECor(int pdgcode, float pt);
  float GetNUACor(int charge, float phi, float eta, float vz);
  // V0
  bool IsGoodV0(AliAODv0* aodV0);
  bool IsGoodDaughterTrack(const AliAODTrack* track);
  int GetLambdaCode(const AliAODTrack* pTrack, const AliAODTrack* ntrack);
  // Plane
  float GetTPCPlaneNoAutoCorr(std::vector<int> vec_id);
  inline float GetEventPlane(float qx, float qy, float harmonic);
  // Get DCA
  bool GetDCA(float &dcaxy, float &dcaz, AliAODTrack* track);
  // Sum pT bin, Delta eta bin
  inline float GetSumPtBin(float sumPt);
  inline float GetDeltaEtaBin(float deltaEta);
  inline float GetDCABin(float dca);

  //////////////////////
  // Switch           //
  //////////////////////
  bool isCalculateLambdaHadron;
  bool isCalculateLambdaPion;
  bool isCalculateLambdaProton;
  bool isCalculateLambdaLambda;

  bool isDoNUE;
  bool isDoLambdaNUE;
  bool isNarrowDcaCuts768;
  bool isProtonCustomizedDCACut;
  bool isTightPileUp;

  //////////////////////
  // Cuts and options //
  //////////////////////
  // Global
  TString fTrigger; //
  TString fPeriod;  // period
  // Event
  float fVzCut;      // vz cut
  // Plane
  TString fPlaneEstimator; // TPC or V0C
  // Track
  int fFilterBit;          // AOD filter bit selection
  int fNclsCut;            // ncls cut for all tracks
  float fChi2Max;          // upper limmit for chi2
  float fChi2Min;          // lower limmit for chi2
  // Hadron
  float fSpecialHadronDCAzMax;     // upper limit for track DCAz
  // Proton
  float fSpecialProtonDCAzMax;    // upper limit for proton DCAz
  // PID
  float fNSigmaTPC;
  float fNSigmaRMS;
  // V0
  float fV0CPAMin;                 //
  float fV0DecayLengthMin;         //
  float fV0DecayLengthMax;         //
  float fV0DcaBetweenDaughtersMax; //
  // V0 Daughter
  float fDaughtersTPCNclsMin;      //
  float fDaughtersDCAToPrimVtxMin; //
  float fRatioCrossedRowsFindable; //
  float fDaughtersNSigmaTPC;       //

  ///////////////////The following files are from the data//////////////////////////////////
  /////////////
  // Handles //
  /////////////
  AliAODEvent* fAOD;            //!<! aod Event
  AliPIDResponse* fPIDResponse; //!<! PID Handler
  AliMultSelection* fMultSel;   //!<!

  ////////////////////////////////
  // Global Variables from data //
  ////////////////////////////////
  std::array<double, 3> fVertex{};
  int fRunNum;       // runnumber
  int fOldRunNum;    // latest runnumber
  int fRunNumBin;    // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
  int fVzBin;        // vertex z bin
  float fCent;       // centrality
  int fCentBin;      // centrality bin: 0-10
  // Variable to get TPC Plane
  float fSumQ2xTPC;
  float fSumQ2yTPC;
  float fWgtMultTPC;
  // Plane
  float fPsi2;

  // Plane tracks Map key:id value:(phi,weight)
  std::unordered_map<int, std::vector<float>> mapTPCTrksIDPhiWgt;

  // Vector for particles from Tracks [pt,eta,phi,id,pdgcode,weight,pidweight,dcaxy]
  std::vector<std::array<float,8>> vecParticle;
  // Vector for V0s [pt,eta,phi,id,pdgcode,weight,mass,id1,id2]
  std::vector<std::array<float,9>> vecParticleV0;

  ///////////////////The following files are read from external sources////////////////////
  ////////////////////////
  // Pile up Function
  ////////////////////////
  std::unique_ptr<TF1> fSPDCutPU;
  std::unique_ptr<TF1> fV0CutPU;
  std::unique_ptr<TF1> fCenCutLowPU;
  std::unique_ptr<TF1> fCenCutHighPU;
  std::unique_ptr<TF1> fMultCutPU;
  ////////////////////////
  // NUE
  ////////////////////////
  TList* fListNUE; // read list for NUE
  ////////////////////////
  // NUA
  ////////////////////////
  TList* fListNUA; // read lists for NUA
  TH3F* hCorrectNUAPos; //
  TH3F* hCorrectNUANeg; //
  ////////////////////////
  // VZERO
  ////////////////////////
  TList* fListVZEROCalib; // read list for V0 Calib
  //for recenter
  AliOADBContainer* contQxncm;
  AliOADBContainer* contQyncm;
  TH1D* hQx2mV0C;
  TH1D* hQy2mV0C;
  // for gain equalization
  TH2F* fHCorrectV0ChWeghts;


  ///////////////////The following files will be saved//////////////////////////////////
  //////////////
  // QA Plots //
  //////////////
  TList* fQAList;  //!<!
  // General QA
  // Event-wise
  TH1D* fEvtCount;  //!<!
  std::map<int, int>* runNumList;  //!<!
  TH1I* fHistRunNumBin;  //!<!
  std::array<TH1D*, 2> fHistCent{};  //!<!
  std::array<TH1D*, 2> fHistVz{}; //!<!
  std::array<TH2D*, 6> fHist2CentQA{}; //!<!
  std::array<TH2D*, 2> fHist2MultCentQA{}; //!<!
  std::array<TH2D*, 6> fHist2MultMultQA{}; //!<!
  // Track-wise
  TH1D* fHistPt;  //!<!
  TH1D* fHistEta;  //!<!
  TH1D* fHistNhits;  //!<!
  TH2D* fHist2PDedx;  //!<!
  TH1D* fHistDcaXY;  //!<!
  TH1D* fHistDcaZ;  //!<!
  std::array<TH1D*, 2> fHistPhi{}; //!<!

  //Proton QA
  TH1D* fHistProtonPt;  //!<!
  TH1D* fHistProtonEta;  //!<!
  TH1D* fHistProtonPhi;  //!<!
  TH2D* fHistProtonPtDcaXY;  //!<!
  TH2D* fHistProtonPtDcaZ;  //!<!
  TH1D* fHistAntiProtonPt;  //!<!
  TH1D* fHistAntiProtonEta;  //!<!
  TH1D* fHistAntiProtonPhi;  //!<!
  TH2D* fHistAntiProtonPtDcaXY;  //!<!
  TH2D* fHistAntiProtonPtDcaZ;  //!<!

  // V0s QA
  TH1D* fHistV0Pt;              //!<! Raw V0s' pT
  TH1D* fHistV0Eta;             //!<! Raw V0s' eta
  TH1D* fHistV0DcatoPrimVertex; //!<! Raw V0s' DcatoPV
  TH1D* fHistV0CPA;             //!<! Raw V0s' CPA(cosine pointing angle)
  TH1D* fHistV0DecayLength;     //!<! Raw V0s' DecayLength
  TH1D* fHistV0NegDaughterDca;  //!<! Raw V0s' NegDaughterDca
  TH1D* fHistV0PosDaughterDca;  //!<! Raw V0s' PosDaughterDca
  // Lambda QA
  //[0]:Before the Mass Cut [1]:After the Mass Cut
  std::array<TH1D*,2> fHistLambdaPt{};                  //!<!
  std::array<TH1D*,2> fHistLambdaEta{};                 //!<!
  std::array<TH1D*,2> fHistLambdaPhi{};                 //!<!
  std::array<TH1D*,2> fHistLambdaDcaToPrimVertex{};     //!<!
  std::array<TH1D*,2> fHistLambdaNegDaugtherDca{};      //!<!
  std::array<TH1D*,2> fHistLambdaPosDaugtherDca{};      //!<!
  std::array<TH1D*,2> fHistLambdaCPA{};                 //!<!
  std::array<TH1D*,2> fHistLambdaDecayLength{};         //!<!
  std::array<TH3D*,2> fHist3LambdaCentPtMass{};         //!<!
  std::array<TH2D*,2> fHist2LambdaMassPtY{};            //!<!
  std::array<TH1D*,2> fHistAntiLambdaPt{};              //!<!
  std::array<TH1D*,2> fHistAntiLambdaEta{};             //!<!
  std::array<TH1D*,2> fHistAntiLambdaPhi{};             //!<!
  std::array<TH1D*,2> fHistAntiLambdaDcaToPrimVertex{}; //!<!
  std::array<TH1D*,2> fHistAntiLambdaNegDaugtherDca{};  //!<!
  std::array<TH1D*,2> fHistAntiLambdaPosDaugtherDca{};  //!<!
  std::array<TH1D*,2> fHistAntiLambdaCPA{};             //!<!
  std::array<TH1D*,2> fHistAntiLambdaDecayLength{};     //!<!
  std::array<TH3D*,2> fHist3AntiLambdaCentPtMass{};     //!<!
  std::array<TH2D*,2> fHist2AntiLambdaMassPtY{};        //!<!


  /////////////
  // Results //
  /////////////

  TList* fResultsList;  //!<!
  // Plane
  TH2D* fHist2Psi2;

  // Lambda - Proton
  // Inv Mass
  std::array<TH3D*, 4> fHist3LambdaProtonMassIntg{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  std::array<TH3D*, 4> fHist3LambdaProtonMassSPt{};  //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  std::array<TH3D*, 4> fHist3LambdaProtonMassDEta{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  // Diff δ
  std::array<TProfile3D*,4> fProfile3DDeltaLambdaProtonMassIntg{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  std::array<TProfile3D*,4> fProfile3DDeltaLambdaProtonMassSPt{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  std::array<TProfile3D*,4> fProfile3DDeltaLambdaProtonMassDEta{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  // Diff γ
  std::array<TProfile3D*,4> fProfile3DGammaLambdaProtonMassIntg{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  std::array<TProfile3D*,4> fProfile3DGammaLambdaProtonMassSPt{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  std::array<TProfile3D*,4> fProfile3DGammaLambdaProtonMassDEta{}; //!<!  [0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar

  // Lambda - Hadron
  // Inv Mass
  std::array<TH3D*, 4> fHist3LambdaHadronMassIntg{}; //!<!
  // Diff δ
  std::array<TProfile3D*,4> fProfile3DDeltaLambdaHadronMassIntg{}; //!<!  [0]:Λ-h+ [1]:Λ-h-   [2]:Λbar-h+ [3]:Λbar-h-
  // Diff γ
  std::array<TProfile3D*,4> fProfile3DGammaLambdaHadronMassIntg{}; //!<!  [0]:Λ-h+ [1]:Λ-h-   [2]:Λbar-h+ [3]:Λbar-h-

  // Lambda - Pion
  std::array<TH3D*, 4> fHist3LambdaPionMassIntg{}; //!<!
  // Diff δ
  std::array<TProfile3D*,4> fProfile3DDeltaLambdaPionMassIntg{}; //!<!  [0]:Λ-h+ [1]:Λ-h-   [2]:Λbar-h+ [3]:Λbar-h-
  // Diff γ
  std::array<TProfile3D*,4> fProfile3DGammaLambdaPionMassIntg{}; //!<!  [0]:Λ-h+ [1]:Λ-h-   [2]:Λbar-h+ [3]:Λbar-h-

  // Lambda - Lambda
  std::array<TH3D*, 4> fHist3LambdaLambdaMassMass{}; //!<!
  // Diff δ
  std::array<TProfile3D*,4> fProfile3DDeltaLambdaLambdaMassMass{}; //!<!  [0]:Λ-Λ [1]:Λ-Λbar   [2]:Λbar-Λ [3]:Λbar-Λbar
  // Diff γ
  std::array<TProfile3D*,4> fProfile3DGammaLambdaLambdaMassMass{}; //!<!  [0]:Λ-Λ [1]:Λ-Λbar   [2]:Λbar-Λ [3]:Λbar-Λbar

  AliAnalysisTaskCVEPIDCMEDiff(const AliAnalysisTaskCVEPIDCMEDiff&);
  AliAnalysisTaskCVEPIDCMEDiff& operator=(const AliAnalysisTaskCVEPIDCMEDiff&);

  ClassDef(AliAnalysisTaskCVEPIDCMEDiff, 3);
};

#endif
