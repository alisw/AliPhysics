#ifndef AliAnalysisTaskLambdaProtonCVE_cxx
#define AliAnalysisTaskLambdaProtonCVE_cxx
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODv0.h"
#include "AliAODZDC.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TComplex.h"
#include "TList.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliEventCuts.h"

class AliAnalysisTaskLambdaProtonCVE : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskLambdaProtonCVE();
  AliAnalysisTaskLambdaProtonCVE(const char *name);
  virtual ~AliAnalysisTaskLambdaProtonCVE();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  //read in
  void SetListForNUE(TList *flist)         {this->fListNUE        = (TList*)flist->Clone();}
  void SetListForNUA(TList *flist)         {this->fListNUA        = (TList*)flist->Clone();}
  void SetListForVZEROCalib(TList *flist)  {this->fListVZEROCalib = (TList*)flist->Clone();}
  void SetListForZDCCalib(TList *flist)    {this->fListZDCCalib   = (TList*)flist->Clone();}
  //Global
  void SetDebug(int debug)                 {this->fDebug   =   debug;}
  void SetTrigger(TString trigger)         {this->fTrigger = trigger;}
  void SetPeriod(TString period)           {this->fPeriod  =  period;}
  //Event
  void SetVzCut(double vzCut)              {this->fVzCut       =       vzCut;}
  void SetCentCut(float centDiffCut)       {this->fCentDiffCut = centDiffCut;}
  //Calib
  void IfTPCCalibOn(bool bTPCCalibOn)      {this->IsTPCCalibOn   =   bTPCCalibOn;}
  void IfVZEROCalibOn(bool bVZEROCalibOn)  {this->IsVZEROCalibOn = bVZEROCalibOn;}
  void IfZDCCalibOn(bool bZDCCalibOn)      {this->IsZDCCalibOn   =   bZDCCalibOn;}
  //Plane
  void SetPlanePtMin(float planePtMin)     {this->fPlanePtMin = planePtMin;}
  void SetPlanePtMax(float planePtMax)     {this->fPlanePtMax = planePtMax;}
  void SetPlaneEtaGapPos(float etaGapPos)  {this->fEtaGapPos  =  etaGapPos;}
  void SetPlaneEtaGapNeg(float etaGapNeg)  {this->fEtaGapNeg  =  etaGapNeg;}
  //Track
  void SetFilterBit(int filterBit)         {this->fFilterBit   =   filterBit;}
  void SetNclsCut(int nclsCut)             {this->fNclsCut     =     nclsCut;}
  void SetChi2Max(float chi2Max)           {this->fChi2Max     =     chi2Max;}
  void SetChi2Min(float chi2Min)           {this->fChi2Min     =     chi2Min;}
  void SetDCAcutZ(float dcaCutz)           {this->fDcaCutz     =     dcaCutz;}
  void SetDCAcutXY(float dcaCutxy)         {this->fDcaCutxy    =    dcaCutxy;}
  void SetPtMin(float ptMin)               {this->fPtMin       =       ptMin;}
  void SetPtMax(float ptMax)               {this->fPtMax       =       ptMax;}
  void SetNUEOn(bool doNUE)                {this->IsDoNUE      =       doNUE;}
  void SetNUAOn(bool doNUA)                {this->IsDoNUA      =       doNUA;}
  void SetProtonPtMax(double protonPtMax)  {this->fProtonPtMax = protonPtMax;}
  //PID
  void SetNSigmaTPCCut(double nSigmaTPC) {this->fNSigmaTPCCut = nSigmaTPC;}
  void SetNSigmaTOFCut(double nSigmaTOF) {this->fNSigmaTOFCut = nSigmaTOF;}
  //V0
  void SetV0PtMin(double v0PtMin)                                   {this->fV0PtMin                  =                  v0PtMin;}
  void SetV0CPAMin(double v0CPAMin)                                 {this->fV0CPAMin                 =                 v0CPAMin;}
  void SetV0RapidityMax(double v0RapidityMax)                       {this->fV0RapidityMax            =            v0RapidityMax;}
  void SetV0DecayLengthMax(double v0DecayLengthMax)                 {this->fV0DecayLengthMax         =         v0DecayLengthMax;}
  void SetV0DecayLengthMin(double v0DecayLengthMin)                 {this->fV0DecayLengthMin         =         v0DecayLengthMin;}
  void SetV0DCAToPrimVtxMax(double v0DCAToPrimVtxMax)               {this->fV0DCAToPrimVtxMax        =        v0DCAToPrimVtxMax;}
  void SetV0DcaBetweenDaughtersMax(double v0DcaBetweenDaughtersMax) {this->fV0DcaBetweenDaughtersMax = v0DcaBetweenDaughtersMax;}
  //V0 Daughter Cut
  void SetDaughtersPtMax(double daughtersPtMax)                     {this->fDaughtersPtMax           =           daughtersPtMax;}
  void SetDaughtersEtaMax(double daughtersEtaMax)                   {this->fDaughtersEtaMax          =          daughtersEtaMax;}
  void SetDaughtersTPCNclsMin(double daughtersTPCNclsMin)           {this->fDaughtersTPCNclsMin      =      daughtersTPCNclsMin;}
  void SetDaughtersDCAToPrimVtxMin(double daughtersDCAToPrimVtxMin) {this->fDaughtersDCAToPrimVtxMin = daughtersDCAToPrimVtxMin;}
  void SetPosProtonTPCNsigma(float V0PosProtonTPCNsigma)            {this->fV0PosProtonTPCNsigma     =     V0PosProtonTPCNsigma;}
  void SetNegPionTPCNsigma(float   V0NegPionTPCNsigma)              {this->fV0NegPionTPCNsigma       =       V0NegPionTPCNsigma;}
  void SetNegProtonTPCNsigma(float V0NegProtonTPCNsigma)            {this->fV0NegProtonTPCNsigma     =     V0NegProtonTPCNsigma;}
  void SetPosPionTPCNsigma(float   V0PosPionTPCNsigma)              {this->fV0PosPionTPCNsigma       =       V0PosPionTPCNsigma;}
  void IfDaughtersPIDUseTOF(bool daughterPIDUseTOF)                 {this->IsV0DaughterUseTOF        =        daughterPIDUseTOF;}
  void SetPosProtonTOFNsigma(float V0PosProtonTOFNsigma)            {this->fV0PosProtonTOFNsigma     =     V0PosProtonTOFNsigma;}
  void SetNegPionTOFNsigma(float   V0NegPionTOFNsigma)              {this->fV0NegPionTOFNsigma       =       V0NegPionTOFNsigma;}
  void SetNegProtonTOFNsigma(float V0NegProtonTOFNsigma)            {this->fV0NegProtonTOFNsigma     =     V0NegProtonTOFNsigma;}
  void SetPosPionTOFNsigma(float   V0PosPionTOFNsigma)              {this->fV0PosPionTOFNsigma       =       V0PosPionTOFNsigma;}
  //Lambda Mass Cut
  void SetMassMean(double massMean)           {this->fMassMean      =      massMean;}
  void SetLambdaMassCut(double lambdaMassCut) {this->fLambdaMassCut = lambdaMassCut;}
  //QA switch
  void IfQATPC(bool bQATPC)                   {this->IsQATPC        =        bQATPC;}
  void IfQAVZERO(bool bQAVZERO)               {this->IsQAVZERO      =      bQAVZERO;}
  void IfQAZDC(bool bQAZDC)                   {this->IsQAZDC        =        bQAZDC;}
  void IfCheckPIDFlow(bool bCheckPIDFlow)     {this->IsCheckPIDFlow = bCheckPIDFlow;}

private:
  ////////////////////////
  // Procedural function
  ////////////////////////
  bool    GetVZEROPlane();
  bool      GetZDCPlane();
  bool GetZDCPlaneLsFit();
  bool       LoopTracks();
  bool      GetTPCPlane();
  bool          LoopV0s();
  bool       PairLambda();
  void     ResetVectors();

  ////////////////////////
  //Functional function
  ////////////////////////
  // Read in
  bool      LoadCalibHistForThisRun();
  // Pile-up
  bool           RejectEvtMultComp();
  bool              RejectEvtTFFit();
  bool      RejectEvtTPCITSfb32TOF();
  bool              AODPileupCheck();
  bool           PileUpMultiVertex();
  bool              RemovalForRun1();
  double                  GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  // Track
  bool              AcceptAODTrack(AliAODTrack *track);
  bool          CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck);
  double                 GetNUECor(int charge, double pt);
  double                 GetNUACor(int charge, double phi, double eta, double vz);
  // V0
  bool                    IsGoodV0(AliAODv0 *aodV0);
  bool         IsGoodDaughterTrack(const AliAODTrack *track);
  int                GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *ntrack);
  // Plane
  double             GetEventPlane(double qx, double qy, double harmonic);

  //////////////////////
  // Cuts and options //
  //////////////////////
  //Global
  int                       fDebug; // debug level controls amount of output statements
  TString                 fTrigger; //
  TString                  fPeriod; // period

  //Event
  double                    fVzCut; // vz cut
  float               fCentDiffCut; // centrality restriction for V0M and TRK

  //Calib
  bool              IsVZEROCalibOn; // switch for VZERO qn calib
  bool                IsZDCCalibOn;
  bool                IsTPCCalibOn; // switch for TPC qn calib
  bool                   IsQAVZERO; // flag for V0 qn QA
  bool                     IsQAZDC;
  bool                     IsQATPC;

  //Plane
  float                fPlanePtMin;
  float                fPlanePtMax;
  float                 fEtaGapPos; // value for the Eta Gap Pos
  float                 fEtaGapNeg;

  //Track
  int                   fFilterBit; // AOD filter bit selection
  int                     fNclsCut; // ncls cut for all tracks
  float                   fChi2Max; // upper limmit for chi2
  float                   fChi2Min; // lower limmit for chi2
  float                   fDcaCutz; // dcaz cut for all tracks
  float                  fDcaCutxy; // dcaxy cut for all tracks
  float                     fPtMin; // minimum pt for track
  float                     fPtMax; // maximum pt for track
  bool                     IsDoNUE; // switch for NUE
  bool                     IsDoNUA; // switch for NUA
  double              fProtonPtMax; // Max pt for proton
  //PID
  float             fNSigmaTPCCut;
  float             fNSigmaTOFCut;

  //V0
  double                  fV0PtMin; //
  double                 fV0CPAMin; //
  double            fV0RapidityMax; //
  double         fV0DecayLengthMin; //
  double         fV0DecayLengthMax; //
  double        fV0DCAToPrimVtxMax; //
  double fV0DcaBetweenDaughtersMax; //

  //V0 Daughter
  double           fDaughtersPtMax; //
  double          fDaughtersEtaMax; //
  double      fDaughtersTPCNclsMin; //
  double fDaughtersDCAToPrimVtxMin; //
  float      fV0PosProtonTPCNsigma; //
  float        fV0NegPionTPCNsigma; //
  float      fV0NegProtonTPCNsigma; //
  float        fV0PosPionTPCNsigma; //
  bool          IsV0DaughterUseTOF; //
  float      fV0PosProtonTOFNsigma; //
  float        fV0NegPionTOFNsigma; //
  float      fV0NegProtonTOFNsigma; //
  float        fV0PosPionTOFNsigma; //

  //Lambda Mass
  double            fLambdaMassCut; //

  //Check PID Flow
  bool              IsCheckPIDFlow;

  // Global Variables Unchanged in an Evt
  double                 fMassMean; //
  const float              fEtaCut; // eta cut
  const float             fDedxCut; // dedx cut

  ///////////////////The following files are from the data//////////////////////////////////
  /////////////
  // Handles //
  /////////////
  AliAODEvent*                fAOD; // aod Event
  AliPIDResponse*     fPIDResponse; // PID Handler
  AliAnalysisUtils*         fUtils; // Event Selection Options

  ////////////////////////////////
  // Global Variables from data //
  ////////////////////////////////
  std::map<int,int>*    runNumList;
  double                fVertex[3];
  int                      fRunNum; // runnumber
  int                   fOldRunNum;
  int                   fRunNumBin; // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
  int                       fVzBin; // vertex z bin
  int                     fCentBin; // centrality bin: 0-10
  double                     fCent; // value of centrality

  double                  fPsi1ZNC;
  double                  fPsi1ZNA;
  double                  fPsi2V0C;
  double                  fPsi2V0A;
  double               fPsi2TPCPos;
  double               fPsi2TPCNeg;

  double              SumQ2xTPCPos;
  double              SumQ2yTPCPos;
  double            fWgtMultTPCPos;

  double              SumQ2xTPCNeg;
  double              SumQ2yTPCNeg;
  double            fWgtMultTPCNeg;

  //Plane Vector
  std::vector<int>           vecPosEPTrkID;
  std::vector<int>           vecNegEPTrkID;
  std::vector<double>       vecPosEPTrkPhi;
  std::vector<double>       vecNegEPTrkPhi;
  std::vector<double>    vecPosEPTrkNUAWgt;
  std::vector<double>    vecNegEPTrkNUAWgt;

  //Track Vector
  std::vector<int>            vecPDGCode;
  std::vector<int>                 vecID;
  std::vector<double>             vecPhi;
  std::vector<double>             vecEta;
  std::vector<double>              vecPt;
  std::vector<double>       vecNUAWeight;
  std::vector<double>       vecNUEWeight;
  std::vector<double>    vecNUAWeightPID;
  std::vector<double>    vecNUEWeightPID;

  //Lambda Vector
  std::vector<int>          vecLambdaCode;
  std::vector<double>        vecLambdaPhi;
  std::vector<double>         vecLambdaPt;
  std::vector<int>       vecDaughterPosID; // Pos Daughter ID
  std::vector<int>       vecDaughterNegID; // Neg Daughter ID

  ///////////////////The following files are read from external sources////////////////////
  ////////////////////////
  // Pile up Function
  ////////////////////////
  TF1*                       fSPDCutPU;
  TF1*                        fV0CutPU;
  TF1*                    fCenCutLowPU;
  TF1*                   fCenCutHighPU;
  TF1*                      fMultCutPU;
  ////////////////////////
  // NUE
  ////////////////////////
  //10h/15o
  TList*                      fListNUE; // read list for NUE
  TH1D*                 hNUEweightPlus;
  TH1D*                hNUEweightMinus;
  ////////////////////////
  // NUA
  ////////////////////////
  TList*          fListNUA; // read lists for NUA
  //10h
  TH2D*     hNUAweightPlus;
  TH2D*    hNUAweightMinus;
  //15o
  TH3F*    hCorrectNUAPos; // Protty
  TH3F*    hCorrectNUANeg; // Protty
  ////////////////////////
  // VZERO
  ////////////////////////
  TList*         fListVZEROCalib; // read list for V0 Calib
  //10h
  TH2D*              hMultV0Read;
  TProfile3D*    pV0XMeanRead[3];
  TProfile3D*    pV0YMeanRead[3];
  //15o
  TH1D*                   hMultV0; //Dobrin
  AliOADBContainer*      contMult;
  AliOADBContainer*     contQxncm;
  AliOADBContainer*     contQyncm;
  AliOADBContainer*     contQxnam;
  AliOADBContainer*     contQynam;
  TH1D*                hQx2mV0[2];
  TH1D*                hQy2mV0[2];
  //18q
  TH2F*       fHCorrectV0ChWeghts;
  ////////////////////////
  // ZDC
  ////////////////////////
  TList*               fListZDCCalib;
  // 10h
  TTree*                        tree;
  float                 vtxQuant1[3];
  float                 vtxQuant2[3];
  TProfile*         fProfileForZNCGE;
  TProfile*         fProfileForZNAGE;
  THnSparseF*        fHn4DForZNCQxRC;
  THnSparseF*        fHn4DForZNCQyRC;
  THnSparseF*        fHn4DForZNCMtRC;
  THnSparseF*        fHn4DForZNAQxRC;
  THnSparseF*        fHn4DForZNAQyRC;
  THnSparseF*        fHn4DForZNAMtRC;
  THnSparseI*    fHn4DForZNCCountsRC;
  THnSparseI*    fHn4DForZNACountsRC;
  TProfile2D*      fProfile2DForCosC;
  TProfile2D*      fProfile2DForSinC;
  TProfile2D*      fProfile2DForCosA;
  TProfile2D*      fProfile2DForSinA;
  //15o //TODO
  //18q 18r
  TH1D*             fHZDCCparameters;
  TH1D*             fHZDCAparameters;

  ///////////////////The following files will be saved//////////////////////////////////
  ////////////
  //QA Plots//
  ////////////
  // Event-wise QA
  TList*             fOutputList;
  TH1D*                fEvtCount;
  TH1I*           fHistRunNumBin;
  TH1D*             fHistCent[2];
  TH1D*               fHistVz[2];
  TH2D*         fHist2DCentQA[8];
  TH2D*     fHist2DMultCentQA[2];
  TH2D*     fHist2DMultMultQA[6];
  // Track-wise QA
  TH1D*             fHistPt;
  TH1D*            fHistEta;
  TH1D*          fHistNhits;
  TH2D*        fHist2DPDedx;
  TH1D*          fHistDcaXY;
  TH1D*           fHistDcaZ;
  TH1D*         fHistPhi[2];
  TH2D*    fHist2DEtaPhi[2];
  // Psi QA
  //V0C [0]GE [1]RC
  TProfile*          fProfileV0CQxCent[2];
  TProfile*          fProfileV0CQyCent[2];
  TProfile*           fProfileV0CQxVtx[2];
  TProfile*           fProfileV0CQyVtx[2];
  TH2D*        fHist2DCalibPsi2V0CCent[2];
  //V0A [0]GE [1]RC
  TProfile*          fProfileV0AQxCent[2];
  TProfile*          fProfileV0AQyCent[2];
  TProfile*           fProfileV0AQxVtx[2];
  TProfile*           fProfileV0AQyVtx[2];
  TH2D*        fHist2DCalibPsi2V0ACent[2];
  //ZNC [0]Raw [1]GE
  TProfile* fProfileZNCTowerMeanEnegry[2];
  // [0]GE [1]RC
  TProfile*          fProfileZNCQxCent[2];
  TProfile*          fProfileZNCQyCent[2];
  // [0]GE [1]RC [2]SF
  TH2D*        fHist2DCalibPsi1ZNCCent[3];
  //ZNA
  TProfile* fProfileZNATowerMeanEnegry[2];
  TProfile*          fProfileZNAQxCent[2];
  TProfile*          fProfileZNAQyCent[2];
  TH2D*        fHist2DCalibPsi1ZNACent[3];
  //ZNC-ZNA
  TProfile*    fProfileZDCQxAQxCCent[2];
  TProfile*    fProfileZDCQxAQyCCent[2];
  TProfile*    fProfileZDCQyAQxCCent[2];
  TProfile*    fProfileZDCQyAQyCCent[2];
  //V0s QA
  TH1D*                 fHistV0Pt; // Raw V0s' pT
  TH1D*                fHistV0Eta; // Raw V0s' eta
  TH1D*    fHistV0DcatoPrimVertex; // Raw V0s' DcatoPV
  TH1D*                fHistV0CPA; // Raw V0s' CPA(cosine pointing angle)
  TH1D*        fHistV0DecayLength; // Raw V0s' DecayLength
  ///Lambda QA
  //[0]:Before the Mass Cut [1]:After the Mass Cut
  TH1D*                        fHistLambdaPt[2]; //
  TH1D*                       fHistLambdaEta[2]; //
  TH1D*                       fHistLambdaPhi[2]; //
  TH1D*           fHistLambdaDcaToPrimVertex[2]; //
  TH1D*                       fHistLambdaCPA[2]; //
  TH1D*               fHistLambdaDecayLength[2]; //
  TH1D*                      fHistLambdaMass[2]; //
  TH1D*                    fHistAntiLambdaPt[2]; //
  TH1D*                   fHistAntiLambdaEta[2]; //
  TH1D*                   fHistAntiLambdaPhi[2]; //
  TH1D*       fHistAntiLambdaDcaToPrimVertex[2]; //
  TH1D*                   fHistAntiLambdaCPA[2]; //
  TH1D*           fHistAntiLambdaDecayLength[2]; //
  TH1D*                  fHistAntiLambdaMass[2]; //
  TProfile*           fProfileLambdaMassVsPt[2]; //
  TProfile*       fProfileAntiLambdaMassVsPt[2]; //
  //Flow
  //[0]TPC [1]V0C [2]V0A [3]ZNC [4]ZNA
  TProfile2D*          fProfile2DRawFlowCentPthPos[5];
  TProfile2D*          fProfile2DRawFlowCentPthNeg[5];
  TProfile2D*        fProfile2DRawFlowCentPtProton[5];
  TProfile2D*    fProfile2DRawFlowCentPtAntiProton[5];
  TProfile2D*        fProfile2DRawFlowCentPtLambda[5];
  TProfile2D*    fProfile2DRawFlowCentPtAntiLambda[5];

  /////////////
  // Results //
  /////////////
  // Plane we used
  TH2D*    fHist2DPsi2TPCPosCent;
  TH2D*    fHist2DPsi2TPCNegCent;
  TH2D*       fHist2DPsi2V0CCent;
  TH2D*       fHist2DPsi2V0ACent;
  TH2D*       fHist2DPsi1ZNCCent;
  TH2D*       fHist2DPsi1ZNACent;
  //Res
  TProfile*          fProfileTPCPsi2Correlation; //
  TProfile*          fProfileV0MPsi2Correlation; //
  TProfile*          fProfileZDCPsi1Correlation; //
  TProfile*          fProfileZDCPsi2Correlation; //
  TProfile*    fProfileV0CTPCPosPsi2Correlation; //
  TProfile*    fProfileV0ATPCPosPsi2Correlation; //
  TProfile*    fProfileV0CTPCNegPsi2Correlation; //
  TProfile*    fProfileV0ATPCNegPsi2Correlation; //

  ///Lambda-X correlators
  TProfile*                 fProfileDelta_Lambda_hPos; //
  TProfile*                 fProfileDelta_Lambda_hNeg; //
  TProfile*               fProfileDelta_Lambda_Proton; //
  TProfile*           fProfileDelta_Lambda_AntiProton; //
  TProfile*             fProfileDelta_AntiLambda_hPos; //
  TProfile*             fProfileDelta_AntiLambda_hNeg; //
  TProfile*           fProfileDelta_AntiLambda_Proton; //
  TProfile*       fProfileDelta_AntiLambda_AntiProton; //

  TProfile*              fProfileGammaTPC_Lambda_hPos; //
  TProfile*              fProfileGammaTPC_Lambda_hNeg; //
  TProfile*            fProfileGammaTPC_Lambda_Proton; //
  TProfile*        fProfileGammaTPC_Lambda_AntiProton; //
  TProfile*          fProfileGammaTPC_AntiLambda_hPos; //
  TProfile*          fProfileGammaTPC_AntiLambda_hNeg; //
  TProfile*        fProfileGammaTPC_AntiLambda_Proton; //
  TProfile*    fProfileGammaTPC_AntiLambda_AntiProton; //

  TProfile*              fProfileGammaV0C_Lambda_hPos; //
  TProfile*              fProfileGammaV0C_Lambda_hNeg; //
  TProfile*            fProfileGammaV0C_Lambda_Proton; //
  TProfile*        fProfileGammaV0C_Lambda_AntiProton; //
  TProfile*          fProfileGammaV0C_AntiLambda_hPos; //
  TProfile*          fProfileGammaV0C_AntiLambda_hNeg; //
  TProfile*        fProfileGammaV0C_AntiLambda_Proton; //
  TProfile*    fProfileGammaV0C_AntiLambda_AntiProton; //

  TProfile*              fProfileGammaV0A_Lambda_hPos; //
  TProfile*              fProfileGammaV0A_Lambda_hNeg; //
  TProfile*            fProfileGammaV0A_Lambda_Proton; //
  TProfile*        fProfileGammaV0A_Lambda_AntiProton; //
  TProfile*          fProfileGammaV0A_AntiLambda_hPos; //
  TProfile*          fProfileGammaV0A_AntiLambda_hNeg; //
  TProfile*        fProfileGammaV0A_AntiLambda_Proton; //
  TProfile*    fProfileGammaV0A_AntiLambda_AntiProton; //

  TProfile*              fProfileGammaZNC_Lambda_hPos; //
  TProfile*              fProfileGammaZNC_Lambda_hNeg; //
  TProfile*            fProfileGammaZNC_Lambda_Proton; //
  TProfile*        fProfileGammaZNC_Lambda_AntiProton; //
  TProfile*          fProfileGammaZNC_AntiLambda_hPos; //
  TProfile*          fProfileGammaZNC_AntiLambda_hNeg; //
  TProfile*        fProfileGammaZNC_AntiLambda_Proton; //
  TProfile*    fProfileGammaZNC_AntiLambda_AntiProton; //

  TProfile*              fProfileGammaZNA_Lambda_hPos; //
  TProfile*              fProfileGammaZNA_Lambda_hNeg; //
  TProfile*            fProfileGammaZNA_Lambda_Proton; //
  TProfile*        fProfileGammaZNA_Lambda_AntiProton; //
  TProfile*          fProfileGammaZNA_AntiLambda_hPos; //
  TProfile*          fProfileGammaZNA_AntiLambda_hNeg; //
  TProfile*        fProfileGammaZNA_AntiLambda_Proton; //
  TProfile*    fProfileGammaZNA_AntiLambda_AntiProton; //

  AliAnalysisTaskLambdaProtonCVE(const AliAnalysisTaskLambdaProtonCVE&);
  AliAnalysisTaskLambdaProtonCVE& operator=(const AliAnalysisTaskLambdaProtonCVE&);

  ClassDef(AliAnalysisTaskLambdaProtonCVE, 1);
};

#endif
