#ifndef AliAnalysisTaskSigma1385PM_H
#define AliAnalysisTaskSigma1385PM_H

#include <THnSparse.h>
#include <TNtupleD.h>

#include <deque>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskTrackMixer.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

class AliAnalysisTaskSigma1385PM : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSigma1385PM();
  AliAnalysisTaskSigma1385PM(const char* name, Bool_t MCcase);
  virtual ~AliAnalysisTaskSigma1385PM();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void SetFillQAPlot(Bool_t input) { fFillQAPlot = input; }
  void SetMixerTask(AliAnalysisTaskTrackMixer* mixerTask) {
    fMixingPool = mixerTask;
  }
  void SetUseBuiltinMixer(Bool_t builtinmixer) {
    fUseBuiltinMixer = builtinmixer;
  }
  void SetMixing(Bool_t setmixing) { fSetMixing = setmixing; }
  void SetnMix(Int_t nMix) { fnMix = nMix; }
  void SetIsPrimaryMC(Bool_t isprimarymc) { fIsPrimaryMC = isprimarymc; }
  void SetINEL(Bool_t input) { fIsINEL = input; }
  void SetHighMult(Bool_t input) { fIsHM = input; }
  void SetAsymmCut(Bool_t input) { fUseAsymmCut = input; }
  void SetOnlyUseOnTheFlyV0(Bool_t input) { fOnlyUseOnTheFlyV0 = input; }
  void SetFillnTuple(Bool_t fillntuple) { fFillnTuple = fillntuple; }

  // Setter for cut variables
  void SetFilterbitSigmaStarPion(Double_t lParameter) {
    fFilterBit = lParameter;
  }
  void SetMaxNsigSigmaStarPion(Double_t lParameter) {
    fTPCNsigSigmaStarPionCut = lParameter;
  }
  void SetMaxEtaSigmaStarPion(Double_t lParameter) {
    fSigmaStarPionEtaCut = lParameter;
  }
  void SetMaxVertexZSigmaStarPion(Double_t lParameter) {
    fSigmaStarPionZVertexCut = lParameter;
  }
  void SetMaxVertexXYsigSigmaStarPion(Double_t lParameter) {
    fSigmaStarPionXYVertexSigmaCut = lParameter;
  }

  void SetMaxNsigV0Proton(Double_t lParameter) {
    fTPCNsigLambdaProtonCut = lParameter;
  }
  void SetMaxNsigV0Pion(Double_t lParameter) {
    fTPCNsigLambdaPionCut = lParameter;
  }
  void SetMaxDCAV0daughters(Double_t lParameter) {
    fDCADistLambdaDaughtersCut = lParameter;
  }
  void SetMaxDCAPVV0(Double_t lParameter) { fDCArDistLambdaPVCut = lParameter; }
  void SetMaxDCAPVV0PosDaughter(Double_t lParameter) {
    fDCAPositiveTrack = lParameter;
  }
  void SetMaxDCAPVV0NegDaughter(Double_t lParameter) {
    fDCANegativeTrack = lParameter;
  }
  void SetMinCPAV0(Double_t lParameter) {
    fV0CosineOfPointingAngleCut = lParameter;
  }
  void SetMaxRapidityV0(Double_t lParameter) {
    fMaxLambdaRapidity = lParameter;
  }
  void SetLowRadiusV0(Double_t lParameter) { fLambdaLowRadius = lParameter; }
  void SetHighRadiusV0(Double_t lParameter) { fLambdaHighRadius = lParameter; }
  void SetLifetimeV0(Double_t lParameter) { fLambdaLifetime = lParameter; }
  void SetMaxMassWindowV0(Double_t lParameter) {
    fV0MassWindowCut = lParameter;
  }

  void SetSigmaStarRapidityCutHigh(Double_t lParameter) {
    fSigmaStarYCutHigh = lParameter;
  }
  void SetSigmaStarRapidityCutLow(Double_t lParameter) {
    fSigmaStarYCutLow = lParameter;
  }
  void SetSigmaStarAsymmCutHigh(Double_t lParameter) {
    fSigmaStarAsymmCutHigh = lParameter;
  }
  void SetSigmaStarAsymmCutLow(Double_t lParameter) {
    fSigmaStarAsymmCutLow = lParameter;
  }

  Bool_t GoodTracksSelection();
  Bool_t GoodV0Selection();
  void FillTracks();
  void FillNtuples();
  void FillMCinput(AliMCEvent* fMCEvent, int Fillbin = 0);
  void FillTrackToEventPool();
  Bool_t IsTrueSigmaStar(UInt_t v0, UInt_t pion, UInt_t BkgCheck = 0);
  double GetTPCnSigma(AliVTrack* track, AliPID::EParticleType type);
  void GetImpactParam(AliVTrack* track, Float_t p[2], Float_t cov[3]);

  Bool_t IsSelectedTPCGeoCut(AliAODTrack* track);
  Bool_t IsSelectedTPCGeoCut(AliESDtrack* track);
  void SetCutOpen();

  // helper
  THnSparse* CreateTHnSparse(TString name,
                             TString title,
                             Int_t ndim,
                             std::vector<TAxis> bins,
                             Option_t* opt = "");
  THnSparse* CreateTHnSparse(TString name,
                             TString title,
                             TString templ,
                             Option_t* opt = "");
  Long64_t FillTHnSparse(TString name,
                         std::vector<Double_t> x,
                         Double_t w = 1.);
  Long64_t FillTHnSparse(THnSparse* h,
                         std::vector<Double_t> x,
                         Double_t w = 1.);
  TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
  TAxis AxisVar(TString name, std::vector<Double_t> bin);
  TAxis AxisStr(TString name, std::vector<TString> bin);

  AliEventCuts fEventCuts;  // Event cuts
  // TPC GeoCut
  Bool_t fCheckTPCGeo;                   //
  Double_t fTPCActiveLengthCutDeltaY;    //
  Double_t fTPCActiveLengthCutDeltaZ;    //
  Double_t fRequireCutGeoNcrNclLength;   //
  Double_t fRequireCutGeoNcrNclGeom1Pt;  //
  Double_t fCutGeoNcrNclFractionNcr;     //
  Double_t fCutGeoNcrNclFractionNcl;     //

 private:
  typedef std::vector<AliVTrack*> tracklist;
  typedef std::deque<tracklist> eventpool;
  typedef std::vector<std::vector<eventpool>> mixingpool;

  AliESDtrackCuts* fTrackCuts;             //!
  AliPIDResponse* fPIDResponse;            //!
  AliAnalysisTaskTrackMixer* fMixingPool;  //!

  AliVEvent* fEvt;             //!
  AliMCEvent* fMCEvent;        //!
  THistManager* fHistos;       //!
  AliAODVertex* fVertex;       //!
  TNtupleD* fNtupleSigma1385;  //!
  TClonesArray* fMCArray;      //!

  Bool_t fIsAOD;              //!
  Bool_t fIsNano;             //!
  Bool_t fSetMixing;          //
  Bool_t fUseBuiltinMixer;    //
  Bool_t fFillQAPlot;         //
  Bool_t fIsMC;               //
  Bool_t fIsPrimaryMC;        //
  Bool_t fFillnTuple;         //
  Bool_t fIsINEL;             //
  Bool_t fIsHM;               //
  Bool_t fUseAsymmCut;        //
  Bool_t fOnlyUseOnTheFlyV0;  //

  mixingpool fEMpool;  //!
  TAxis fBinCent;      //!
  TAxis fBinZ;         //!
  Double_t fPosPV[3];  //!
  Double_t fMagField;  //!

  Double_t fCent;  //!
  Int_t fnMix;     //!
  Int_t fCentBin;  //!
  Int_t fZbin;     //!

  // Pion cuts
  UInt_t fFilterBit;                        //
  Double_t fTPCNsigSigmaStarPionCut;        //
  Double_t fSigmaStarPionEtaCut;            //
  Double_t fSigmaStarPionZVertexCut;        //
  Double_t fSigmaStarPionXYVertexSigmaCut;  //

  // Lambda cuts
  Double_t fTPCNsigLambdaProtonCut;      //
  Double_t fTPCNsigLambdaPionCut;        //
  Double_t fDCADistLambdaDaughtersCut;   //
  Double_t fDCArDistLambdaPVCut;         //
  Double_t fDCAPositiveTrack;            //
  Double_t fDCANegativeTrack;            //
  Double_t fV0CosineOfPointingAngleCut;  //
  Double_t fMaxLambdaRapidity;           //
  Double_t fLambdaLowRadius;             //
  Double_t fLambdaHighRadius;            //
  Double_t fLambdaLifetime;              //
  Double_t fV0MassWindowCut;             //

  // Sigma Star cut
  Double_t fSigmaStarYCutHigh;      //
  Double_t fSigmaStarYCutLow;       //
  Double_t fSigmaStarAsymmCutHigh;  //
  Double_t fSigmaStarAsymmCutLow;   //

  // Good track/v0 vector array
  std::vector<UInt_t> fGoodTrackArray;
  std::vector<std::vector<UInt_t>> fGoodV0Array;

  ClassDef(AliAnalysisTaskSigma1385PM, 12);
  // Add rapidity/radius/Lifetime/Y cut of lambda
  // Add NanoOption
  // 4: Add GetImpactParm function for nano
  // 5: Seprate MC Sparse, INEL study capability
  // 6: Update some of deafult vaules
  // 7: Add skipping option for QA histos
  // 8: Rebuild MC part
  // 9: Update class format
  // 10: Fix streamer issue
  // 11: Add Asymm cut option
  // 12: Add OnTheFlyV0 option
};

#endif