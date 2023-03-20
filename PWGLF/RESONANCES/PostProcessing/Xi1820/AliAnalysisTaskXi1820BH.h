#ifndef AliAnalysisTaskXi1820BH_H
#define AliAnalysisTaskXi1820BH_H

#include <THnSparse.h>
#include <TNtupleD.h>

#include <deque>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

class AliAnalysisTaskXi1820BH : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskXi1820BH();
  AliAnalysisTaskXi1820BH(const char *name, Bool_t MCcase);
  virtual ~AliAnalysisTaskXi1820BH();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetFillQAPlot(Bool_t input) { fFillQAPlot = input; }
  void SetMixing(Bool_t setmixing) { fSetMixing = setmixing; }
  void SetnMix(Int_t nMix) { fnMix = nMix; }
  void SetIsPrimaryMC(Bool_t isprimarymc) { fIsPrimaryMC = isprimarymc; }
  void SetINEL(Bool_t input) { fIsINEL = input; }
  void SetRunXi1820Zero(Bool_t input) {fSetZero = input;}
  void SetHighMult(Bool_t input) { fIsHM = input; }
  void SetAsymmCut(Bool_t input) { fUseAsymmCut = input; }
  void SetOnlyUseOnTheFlyV0(Bool_t input) { fOnlyUseOnTheFlyV0 = input; }

  // Setter for cut variables
  void SetFilterbitXi1820ZeroKaon(Double_t lParameter) { fFilterBit = lParameter; }
  void SetMaxNsigXi1820ZeroKaon(Double_t lParameter) { fTPCNsigXi1820ZeroKaonCut = lParameter; }
  void SetMaxEtaXi1820ZeroKaon(Double_t lParameter) { fXi1820ZeroKaonEtaCut = lParameter; }
  void SetMaxVertexZXi1820ZeroKaon(Double_t lParameter) { fXi1820ZeroKaonZVertexCut = lParameter; }
  void SetMaxVertexXYsigXi1820ZeroKaon(Double_t lParameter) { fXi1820KaonZeroXYVertexSigmaCut = lParameter; }
  void SetMaxNsigLambdaProton(Double_t lParameter) { fTPCNsigLambdaProtonCut = lParameter; }
  void SetMaxNsigLambdaPion(Double_t lParameter) { fTPCNsigLambdaPionCut = lParameter; }
  void SetMaxDCALambdadaughters(Double_t lParameter) { fDCADistLambdaDaughtersCut = lParameter; }
  void SetMaxDCAPVK0s(Double_t lParameter) { fDCArDistK0sPVCut = lParameter; }
  void SetMaxDCAPVLambda(Double_t lParameter) { fDCArDistLambdaPVCut = lParameter; }
  void SetMaxDCAPVLambdaPosDaughter(Double_t lParameter) { fDCAPositiveTrackLambda = lParameter; }
  void SetMaxDCAPVLambdaNegDaughter(Double_t lParameter) { fDCANegativeTrackLambda = lParameter; }
  void SetMaxDCAPVK0sPosDaughter(Double_t lParameter) { fDCAPositiveTrackK0s = lParameter; }
  void SetMaxDCAPVK0sNegDaughter(Double_t lParameter) { fDCANegativeTrackK0s = lParameter; }
  void SetMinCPALambda(Double_t lParameter) { fLambdaCosineOfPointingAngleCut = lParameter; }
  void SetMaxRapidityLambda(Double_t lParameter) { fMaxLambdaRapidity = lParameter; }
  void SetLowRadiusLambda(Double_t lParameter) { fLambdaLowRadius = lParameter; }
  void SetHighRadiusLambda(Double_t lParameter) { fLambdaHighRadius = lParameter; }
  void SetLifetimeLambda(Double_t lParameter) { fLambdaLifetime = lParameter; }
  void SetMaxMassWindowLambda(Double_t lParameter) { fLambdaMassWindowCut = lParameter; }

  void SetXi1820RapidityCutHigh(Double_t lParameter) { fXi1820YCutHigh = lParameter; }
  void SetXi1820RapidityCutLow(Double_t lParameter) { fXi1820YCutLow = lParameter; }
  void SetXi1820AsymmCutHigh(Double_t lParameter) { fXi1820AsymmCutHigh = lParameter; }
  void SetXi1820AsymmCutLow(Double_t lParameter) { fXi1820AsymmCutLow = lParameter; }

  Bool_t GoodKaonSelection();
  Bool_t GoodK0sSelection();
  Bool_t GoodLambdaSelection();
  void FillKaonV0();
  void FillK0sV0();
  void FillMCEventProperties(AliMCEvent *fMCEvent);
  void FillMCinput(AliMCEvent *fMCEvent, int Fillbin = 0);
  void FillKaonToEventPool();
  void FillK0sToEventPool();
  Bool_t IsTrueXi1820PM(UInt_t v0, UInt_t kaon, UInt_t BkgCheck = 0);
  Bool_t IsTrueXi1820Zero(UInt_t v0, UInt_t k0s, UInt_t BkgCheck = 0);
  double GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type);
  double GetTOFnSigma(AliVTrack *track, AliPID::EParticleType type);
  void GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3]);

  Bool_t IsSelectedTPCGeoCut(AliAODTrack *track);
  Bool_t IsSelectedTPCGeoCut(AliESDtrack *track);

  // helper
  THnSparse *CreateTHnSparse(TString name,
                             TString title,
                             Int_t ndim,
                             std::vector<TAxis> bins,
                             Option_t *opt = "");
  THnSparse *CreateTHnSparse(TString name,
                             TString title,
                             TString templ,
                             Option_t *opt = "");
  Long64_t FillTHnSparse(TString name,
                         std::vector<Double_t> x,
                         Double_t w = 1.);
  Long64_t FillTHnSparse(THnSparse *h,
                         std::vector<Double_t> x,
                         Double_t w = 1.);
  TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
  TAxis AxisVar(TString name, std::vector<Double_t> bin);
  TAxis AxisStr(TString name, std::vector<TString> bin);

  AliEventCuts fEventCuts; // Event cuts
  // TPC GeoCut
  Bool_t fCheckTPCGeo;                  //
  Double_t fTPCActiveLengthCutDeltaY;   //
  Double_t fTPCActiveLengthCutDeltaZ;   //
  Double_t fRequireCutGeoNcrNclLength;  //
  Double_t fRequireCutGeoNcrNclGeom1Pt; //
  Double_t fCutGeoNcrNclFractionNcr;    //
  Double_t fCutGeoNcrNclFractionNcl;    //

private:
  typedef std::vector<AliVTrack *> tracklist;
  typedef std::vector<AliAODv0 *> v0list;
  typedef std::deque<tracklist> eventpool;
  typedef std::deque<v0list> v0pool;
  typedef std::vector<std::vector<eventpool>> mixingpool;
  typedef std::vector<std::vector<v0pool>> v0mixingpool;

  AliESDtrackCuts *fTrackCuts;  //!
  AliPIDResponse *fPIDResponse; //!

  AliVEvent *fEvt;        //!
  AliMCEvent *fMCEvent;   //!
  THistManager *fHistos;  //!
  AliAODVertex *fVertex;  //!
  TClonesArray *fMCArray; //!

  Bool_t fIsAOD;             //!
  Bool_t fIsNano;            //!
  Bool_t fSetMixing;         //
  Bool_t fSetZero;           //
  Bool_t fFillQAPlot;        //
  Bool_t fIsMC;              //
  Bool_t fIsPrimaryMC;       //
  Bool_t fIsINEL;            //
  Bool_t fIsHM;              //
  Bool_t fUseAsymmCut;       //
  Bool_t fOnlyUseOnTheFlyV0; //

  mixingpool fEMpoolKaon;  //!
  v0mixingpool fEMpoolK0s; //!
  TAxis fBinCent;          //!
  TAxis fBinZ;             //!
  Double_t fPosPV[3];      //!
  Double_t fMagField;      //!

  Double_t fCent; //!
  Int_t fnMix;    //!
  Int_t fCentBin; //!
  Int_t fZbin;    //!

  // Kaon cuts
  UInt_t fFilterBit;                        //
  Double_t fTPCNsigXi1820ZeroKaonCut;       //
  Double_t fTOFNsigXi1820ZeroKaonCut;       //
  Double_t fXi1820ZeroKaonEtaCut;           //
  Double_t fXi1820ZeroKaonZVertexCut;       //
  Double_t fXi1820KaonZeroXYVertexSigmaCut; //

  // K0s cuts
  Double_t fTPCNsigK0sPionCut;           //
  Double_t fDCADistK0sDaughtersCut;      //
  Double_t fDCArDistK0sPVCut;            //
  Double_t fDCAPositiveTrackK0s;         //
  Double_t fDCANegativeTrackK0s;         //
  Double_t fK0sCosineOfPointingAngleCut; //
  Double_t fMaxK0sRapidity;              //
  Double_t fK0sLowRadius;                //
  Double_t fK0sHighRadius;               //
  Double_t fK0sLifetime;                 //
  Double_t fK0sMassWindowCut;            //
  // Lambda cuts
  Double_t fTPCNsigLambdaProtonCut;         //
  Double_t fTPCNsigLambdaPionCut;           //
  Double_t fDCADistLambdaDaughtersCut;      //
  Double_t fDCArDistLambdaPVCut;            //
  Double_t fDCAPositiveTrackLambda;         //
  Double_t fDCANegativeTrackLambda;         //
  Double_t fLambdaCosineOfPointingAngleCut; //
  Double_t fMaxLambdaRapidity;              //
  Double_t fLambdaLowRadius;                //
  Double_t fLambdaHighRadius;               //
  Double_t fLambdaLifetime;                 //
  Double_t fLambdaMassWindowCut;            //

  // Xi1820 cuts
  Double_t fXi1820YCutHigh;     //
  Double_t fXi1820YCutLow;      //
  Double_t fXi1820AsymmCutHigh; //
  Double_t fXi1820AsymmCutLow;  //

  // Good kaon/k0s/lambda vector array
  std::vector<UInt_t> fGoodKaonArray;
  std::vector<UInt_t> fGoodK0sArray;
  std::vector<std::vector<UInt_t>> fGoodLambdaArray;

  ClassDef(AliAnalysisTaskXi1820BH, 1);
  // First class of Xi1820 analysis
};

#endif