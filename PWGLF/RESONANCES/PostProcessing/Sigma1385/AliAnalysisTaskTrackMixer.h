#ifndef AliAnalysisTaskTrackMixer_H
#define AliAnalysisTaskTrackMixer_H

#include <deque>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

class AliAnalysisTaskTrackMixer : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTrackMixer();
  AliAnalysisTaskTrackMixer(const char* name);
  virtual ~AliAnalysisTaskTrackMixer();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  void GoodTracksSelection();
  void FillTrackToEventPool();

  std::vector<std::vector<std::deque<std::vector<AliVTrack*>>>> GetMixingPool() { return fEMpool; }
  std::deque<std::vector<AliVTrack*>>* GetEventPool(int Centbin, int zVertexbin) {
    return &fEMpool[Centbin][zVertexbin];
  }

  void SetFillQAPlot(Bool_t input) { fFillQAPlot = input; }
  void SetnMix(Int_t nMix) { fnMix = nMix; }
  void SetHighMult(Bool_t input) { fIsHM = input; }

  // Setter for cut variables
  void SetTrackFilterbit(Double_t linput) { fTrackFilterBit = linput; }
  void SetTrackMaxNsig(Double_t linput) { fTrackTPCNsigCut = linput; }
  void SetTrackMaxEta(Double_t linput) { fTrackEtaCut = linput; }
  void SetTrackMaxVertexZ(Double_t linput) { fTrackZVertexCut = linput; }
  void SetTrackMaxVertexXYsig(Double_t linput) { fTrackXYVertsigCut = linput; }

  double GetTPCnSigma(AliVTrack* track, AliPID::EParticleType type);
  void GetImpactParam(AliVTrack* track, Float_t p[2], Float_t cov[3]);
  Bool_t IsSelectedTPCGeoCut(AliAODTrack* track);
  Bool_t IsSelectedTPCGeoCut(AliESDtrack* track);

  // helper
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

  AliESDtrackCuts* fTrackCuts;   //!
  AliPIDResponse* fPIDResponse;  //!

  AliVEvent* fEvt;        //!
  THistManager* fHistos;  //!
  AliAODVertex* fVertex;  //!

  Bool_t fIsAOD;       //!
  Bool_t fIsNano;      //!
  Bool_t fFillQAPlot;  //
  Bool_t fIsHM;        //

  mixingpool fEMpool;  //!
  TAxis fBinCent;      //!
  TAxis fBinZ;         //!
  Double_t fPosPV[3];  //!
  Double_t fMagField;  //!

  Double_t fCent;  //!
  Int_t fnMix;     //!
  Int_t fCentBin;  //!
  Int_t fZbin;     //!

  // Track cuts
  UInt_t fTrackFilterBit;             //
  Double_t fTrackTPCNsigCut;          //
  Double_t fTrackEtaCut;              //
  Double_t fTrackZVertexCut;          //
  Double_t fTrackXYVertsigCut;        //
  AliPID::EParticleType fTPCPIDType;  //

  // Good track/v0 std::vector array
  std::vector<UInt_t> fGoodTrackArray;

  ClassDef(AliAnalysisTaskTrackMixer, 1);
  // 1: First version.
};

#endif