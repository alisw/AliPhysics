#ifndef ALIANALYSISTASKHFSIMPLEVERTICES_H
#define ALIANALYSISTASKHFSIMPLEVERTICES_H

/* Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//
//*************************************************************************

#include "AliFJWrapper.h"
#include <map>
#include <string>
#include "DCAFitterN.h"
#include "AliVEvent.h"
#include "AliEventCuts.h"

class TList;
class AliESDEvent;

class AliAnalysisTaskHFSimpleVertices : public AliAnalysisTaskSE
{

 public:
  AliAnalysisTaskHFSimpleVertices();
  virtual ~AliAnalysisTaskHFSimpleVertices();

  virtual void UserExec(Option_t* option);
  virtual void UserCreateOutputObjects();
  virtual void Terminate(Option_t* option);
  void InitFromJson(TString filename);
  void SetUseCutOnSPDVsTrackVtx(Bool_t opt) { fCutOnSPDVsTrackVtx = opt; }
  void SetZVertexMaxRange(Double_t zmax) { fMaxZVert = zmax; }
  void SetUseVertexerTracks() { fSecVertexerAlgo = 0; }
  void SetUseO2Vertexer() { fSecVertexerAlgo = 1; }
  void SetUseKFParticleVertexer() { fSecVertexerAlgo = 2; }
  void SetFindVertexForCascades(Bool_t opt) { fFindVertexForCascades = opt; }
  void SetReadMC(Bool_t read) { fReadMC = read; }
  void SetUseCandidateAnalysisCuts() { fCandidateCutLevel = 2; }
  void SetUseCandidateSkimCuts() { fCandidateCutLevel = 1; }
  void SetUseNoCandidateCuts() { fCandidateCutLevel = 0; }
  void SetUsePtDependentFiducialAcceptance() { fMaxRapidityCand = -999.; }
  void SetMaxRapidityForFiducialAcceptance(Double_t ymax) { fMaxRapidityCand = ymax; }
  void SetUseAliEventCuts(Bool_t opt) { fUseAliEventCuts = opt; }
  void EnableCPUTimeCheck(Bool_t enable = kTRUE, Bool_t milliseconds = kFALSE)
  {
    fEnableCPUTimeCheck = enable;
    fCountTimeInMilliseconds = milliseconds;
  }
  void SetDoJetFinding(Bool_t set) {doJetFinding=set;}
  void SetJetMatching(Bool_t set) {matchedJetsOnly = set;}
  void SetDoJetSubstructure(Bool_t set) {doJetSubstructure=set;}


 private:
  AliAnalysisTaskHFSimpleVertices(const AliAnalysisTaskHFSimpleVertices& source);
  AliAnalysisTaskHFSimpleVertices& operator=(const AliAnalysisTaskHFSimpleVertices& source);

  std::string GetJsonString(const char* jsonFileName, const char* section, const char* key);
  int GetJsonInteger(const char* jsonFileName, const char* section, const char* key);
  int GetJsonBool(const char* jsonFileName, const char* section, const char* key);
  float GetJsonFloat(const char* jsonFileName, const char* section, const char* key);
  float* GetJsonArray(const char* jsonFileName, const char* section, const char* key, int& size);
  float** GetJsonMatrix(const char* jsonFileName, const char* section, const char* key, int& size1, int& size2);
  void InitDefault();
  Int_t GetPtBin(Double_t ptCand, Double_t* ptBinLims, Double_t nPtBins);
  Int_t GetPtBinSingleTrack(Double_t ptTrack);
  void ProcessTriplet(TObjArray* threeTrackArray, AliAODRecoDecay* rd4massCalc3, AliESDVertex* primVtxTrk, AliAODVertex* vertexAODp, float bzkG, double dist12, AliMCEvent* mcEvent);
  Bool_t GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG);
  Int_t SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG, Double_t d0track[2]);
  Bool_t SelectV0(AliESDv0* v0, AliESDVertex* primvtx);
  AliESDVertex* ReconstructSecondaryVertex(TObjArray* trkArray, AliESDVertex* primvtx);
  AliAODVertex* ConvertToAODVertex(AliESDVertex* trkv);
  Int_t SelectInvMassAndPt2prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc2);
  Int_t SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3);
  AliAODRecoDecayHF2Prong* Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG);
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG);
  AliAODRecoCascadeHF* MakeCascade(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG);
  Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y) const;
  Int_t DzeroSkimCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t JpsiSkimCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t DplusSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t DsSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t LcSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t DzeroSelectionCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t DplusSelectionCuts(AliAODRecoDecayHF3Prong* cand, Double_t bzkG);
  Int_t DsSelectionCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t JpsiSelectionCuts(AliAODRecoDecayHF2Prong* cand, AliESDtrack* trk_p, AliESDtrack* trk_n, AliESDVertex* primvtx, float bzkG);
  Int_t LcSelectionCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t MatchToMC(AliAODRecoDecay* rd, Int_t pdgabs, AliMCEvent* mcEvent, Int_t ndgCk, const TObjArray* trkArray, const Int_t* pdgDg) const;
  void FindJets(AliESDEvent *esd, Int_t iNegTrack_0, Int_t iPosTrack_0, AliAODRecoDecayHF2Prong* the2Prong, AliESDVertex* primVtx, AliMCEvent* mcEvent = nullptr, Int_t labD=-100);
  void FindGenJets(AliMCEvent* mcEvent, Int_t labD);
  enum ESelBits2prong { kbitDzero = 0,
                        kbitDzerobar,
                        kbitJpsi };
  enum ESelBits3prong { kbitDplus = 0,
                        kbitDs,
                        kbitLc };
  enum { kMaxNPtBinsDzero = 100,
         kNCutVarsDzero = 14 };
  enum { kMaxNPtBinsJpsi = 100,
         kNCutVarsJpsi = 5 };
  enum { kMaxNPtBinsLc = 100,
         kNCutVarsLc = 7 };
  enum { kMaxNPtBinsDplus = 100,
         kNCutVarsDplus = 8 };
  enum { kMaxNPtBinsDs = 100,
         kNCutVarsDs = 10 };
  enum { kMaxNPtBinsSingleTrack = 100,
         kNCutVarsSingleTrack = 2 };
  enum { kMaxNPtBins2ProngsSkims = 100,
         kNCutVars2Prong = 4 };
  enum { kMaxNPtBins3ProngsSkims = 100,
         kNCutVars3Prong = 4 };

  /// list of triggers avaliable in O2 that are also avaliable in AliPhysics
  std::map<std::string, Int_t> triggerMask = {{"kINT7", AliVEvent::kINT7},
                                              {"kEMC7", AliVEvent::kEMC7},
                                              {"kINT7inMUON", AliVEvent::kINT7inMUON},
                                              {"kMuonSingleLowPt7", AliVEvent::kMuonSingleLowPt7},
                                              {"kMuonSingleHighPt7", AliVEvent::kMuonSingleHighPt7},
                                              {"kMuonUnlikeLowPt7", AliVEvent::kMuonUnlikeLowPt7},
                                              {"kMuonLikeLowPt7", AliVEvent::kMuonLikeLowPt7}};

  TList* fOutput;                     //!<!  list of output histos
  TH1F* fHistNEvents;                 //!<!  histo with N of events
  TH1F* fHistTrackStatus;             //!<!  histo with counts of tracks passing cuts
  TH1F* fHistPtAllTracks;             //!<!  histo with pt all tracks
  TH1F* fHistPtSelTracks;             //!<!  histo with pt selected tracks
  TH1F* fHistTglAllTracks;            //!<!  histo with tgl all tracks
  TH1F* fHistTglSelTracks;            //!<!  histo with tgl selected tracks
  TH1F* fHistEtaAllTracks;            //!<!  histo with eta all tracks
  TH1F* fHistPtSelTracks2prong;       //!<!  histo with pt selected tracks
  TH1F* fHistPtSelTracks3prong;       //!<!  histo with pt selected tracks
  TH1F* fHistPtSelTracksbachelor;     //!<!  histo with pt selected tracks
  TH1F* fHistEtaSelTracks2prong;      //!<!  histo with eta selected tracks
  TH1F* fHistEtaSelTracks3prong;      //!<!  histo with eta selected tracks
  TH1F* fHistEtaSelTracksbachelor;    //!<!  histo with eta selected tracks
  TH1F* fHistImpParAllTracks;         //!<!  histo with d0 all tracks
  TH1F* fHistImpParSelTracks2prong;   //!<!  histo with d0 selected tracks
  TH1F* fHistImpParSelTracks3prong;   //!<!  histo with d0 selected tracks
  TH1F* fHistImpParSelTracksbachelor; //!<!  histo with d0 selected tracks
  TH1F* fHistITSmapAllTracks;         //!<!  histo with its map all tracks
  TH1F* fHistITSmapSelTracks;         //!<!  histo withits map selected tracks

  TH1F* fHistPrimVertX;     //!<!  histo of prim vertex x
  TH1F* fHistPrimVertY;     //!<!  histo of prim vertex y
  TH1F* fHistPrimVertZ;     //!<!  histo of prim vertex z
  TH1F* fHistPrimVertContr; //!<!  histo of prim vertex N contributors
  TH1F* fHist2ProngVertX;   //!<!  histo of 2-prong vertex x
  TH1F* fHist2ProngVertY;   //!<!  histo of 2-prong vertex y
  TH1F* fHist2ProngVertZ;   //!<!  histo of 2-prong vertex z
  TH1F* fHist3ProngVertX;   //!<!  histo of 3-prong vertex x
  TH1F* fHist3ProngVertY;   //!<!  histo of 3-prong vertex y
  TH1F* fHist3ProngVertZ;   //!<!  histo of 3-prong vertex z
  TH1F* fHistDist12LcpKpi;  //!<!  histo of LcpKpi+ distance between primary and secondary vertex reconstructed from the pair of tracks

  TH1F* fHistInvMassD0;           //!<!  histo with D0 inv mass
  TH1F* fHistPtD0;                //!<!  histo with D0 pt
  TH2F* fHistYPtD0;               //!<!  histo with D0 y vs pt
  TH1F* fHistPtD0Dau0;            //!<!  histo with D0 prong pt
  TH1F* fHistPtD0Dau1;            //!<!  histo with D0 prong pt
  TH1F* fHistImpParD0Dau0;        //!<!  histo with D0 prong d0
  TH1F* fHistImpParD0Dau1;        //!<!  histo with D0 prong d0
  TH1F* fHistImpParErrD0Dau0;     //!<!  histo with D0 prong d0 err
  TH1F* fHistImpParErrD0Dau1;     //!<!  histo with D0 prong d0 err
  TH1F* fHistd0Timesd0;           //!<!  histo with d0xd0
  TH1F* fHistCosPointD0;          //!<!  histo with D0 cosine of pointing angle
  TH1F* fHistDecLenD0;            //!<!  histo with D0 decay length
  TH1F* fHistDecLenXYD0;          //!<!  histo with D0 decay length XY
  TH1F* fHistDecLenErrD0;         //!<!  histo with D0 decay length err
  TH1F* fHistDecLenXYErrD0;       //!<!  histo with D0 decay length XY err
  TH1F* fHistNormDecLenD0;        //!<!  histo with D0 normalised decay length
  TH1F* fHistNormDecLenXYD0;      //!<!  histo with D0 normalised decay length XY
  TH1F* fHistCovMatPrimVXX2Prong; //!<!  histo with cov mat prim vert for the 2-prong candidate
  TH1F* fHistCovMatSecVXX2Prong;  //!<!  histo with cov mat sec vert for the 2-prong candidate
  TH1F* fHistCovMatPrimVYY2Prong; //!<!  histo with cov mat prim vert for the 2-prong candidate
  TH1F* fHistCovMatSecVYY2Prong;  //!<!  histo with cov mat sec vert for the 2-prong candidate
  TH1F* fHistCovMatPrimVXZ2Prong; //!<!  histo with cov mat prim vert for the 2-prong candidate
  TH1F* fHistCovMatSecVXZ2Prong;  //!<!  histo with cov mat sec vert for the 2-prong candidate
  TH1F* fHistCovMatPrimVZZ2Prong; //!<!  histo with cov mat prim vert for the 2-prong candidate
  TH1F* fHistCovMatSecVZZ2Prong;  //!<!  histo with cov mat sec vert for the 2-prong candidate
  TH1F* fHistD0SignalVertX;       //!<!  histo of D0 (MC truth) vertex x
  TH1F* fHistD0SignalVertY;       //!<!  histo of D0 (MC truth) vertex y
  TH1F* fHistD0SignalVertZ;       //!<!  histo of D0 (MC truth) vertex z
  TH1F* fHistInvMassD0Signal;     //!<!  histo with D0 (MC truth) inv mass
  TH1F* fHistInvMassD0Refl;       //!<!  histo with D0 (reflection) inv mass

  TH1F* fHistInvMassJpsi;       //!<!  histo with Jpsi inv mass
  TH1F* fHistPtJpsi;            //!<!  histo with Jpsi pt
  TH1F* fHistPtJpsiDau0;        //!<!  histo with Jpsi prong pt
  TH1F* fHistPtJpsiDau1;        //!<!  histo with Jpsi prong pt
  TH1F* fHistImpParJpsiDau0;    //!<!  histo with Jpsi prong d0
  TH1F* fHistImpParJpsiDau1;    //!<!  histo with Jpsi prong d0
  TH1F* fHistd0Timesd0Jpsi;     //!<!  histo with d0xd0
  TH1F* fHistCosPointJpsi;      //!<!  histo with Jpsi cosine of pointing angle
  TH1F* fHistDecLenJpsi;        //!<!  histo with Jpsi decay length
  TH1F* fHistDecLenXYJpsi;      //!<!  histo with Jpsi decay length XY
  TH1F* fHistDecLenErrJpsi;     //!<!  histo with Jpsi decay length err
  TH1F* fHistDecLenXYErrJpsi;   //!<!  histo with Jpsi decay length XY err
  TH1F* fHistJpsiSignalVertX;   //!<!  histo of Jpsi (MC truth) vertex x
  TH1F* fHistJpsiSignalVertY;   //!<!  histo of Jpsi (MC truth) vertex y
  TH1F* fHistJpsiSignalVertZ;   //!<!  histo of Jpsi (MC truth) vertex z
  TH1F* fHistInvMassJpsiSignal; //!<!  histo with Jpsi (MC truth) inv mass

  TH1F* fHistInvMassDplus;         //!<!  histo with D+ inv mass
  TH1F* fHistInvMassDplusSignal;   //!<!  histo with D+ inv mass (only signal)
  TH1F* fHistPtDplus;              //!<!  histo with D+ pt
  TH2F* fHistYPtDplus;             //!<!  histo with D+ y vs pt
  TH1F* fHistPtDplusDau0;          //!<!  histo with D+ prong pt
  TH1F* fHistPtDplusDau1;          //!<!  histo with D+ prong pt
  TH1F* fHistPtDplusDau2;          //!<!  histo with D+ prong pt
  TH1F* fHistImpParDplusDau0;      //!<!  histo with D+ prong d0
  TH1F* fHistImpParDplusDau1;      //!<!  histo with D+ prong d0
  TH1F* fHistImpParDplusDau2;      //!<!  histo with D+ prong d0
  TH1F* fHistDecLenDplus;          //!<!  histo with D+ decay length
  TH1F* fHistDecLenXYDplus;        //!<!  histo with D+ decay length XY
  TH1F* fHistNormDecLenXYDplus;    //!<!  histo with D+ normalized decay length XY
  TH1F* fHistImpParErrDplusDau;    //!<!  histo with D+ prong d0 err
  TH1F* fHistDecLenErrDplus;       //!<!  histo with D+ decay length err
  TH1F* fHistDecLenXYErrDplus;     //!<!  histo with D+ decay length XY err
  TH1F* fHistCosPointDplus;        //!<!  histo with D+ cosine of pointing angle
  TH1F* fHistCosPointXYDplus;      //!<!  histo with D+ cosine of pointing angle XY
  TH1F* fHistImpParXYDplus;        //!<!  histo with D+ impact parameter XY
  TH1F* fHistNormIPDplus;          //!<!  histo with max difference between prong observed and expeceted impact parameters
  TH1F* fHistoSumSqImpParDplusDau; //!<!  histo with squared sum of prong impact parameters
  TH1F* fHistCovMatPrimVXX3Prong;  //!<!  histo with cov mat prim vert for the 3-prong candidate
  TH1F* fHistCovMatSecVXX3Prong;   //!<!  histo with cov mat sec vert for the 3-prong candidate
  TH1F* fHistCovMatPrimVYY3Prong;  //!<!  histo with cov mat prim vert for the 3-prong candidate
  TH1F* fHistCovMatSecVYY3Prong;   //!<!  histo with cov mat sec vert for the 3-prong candidate
  TH1F* fHistCovMatPrimVXZ3Prong;  //!<!  histo with cov mat prim vert for the 3-prong candidate
  TH1F* fHistCovMatSecVXZ3Prong;   //!<!  histo with cov mat sec vert for the 3-prong candidate
  TH1F* fHistCovMatPrimVZZ3Prong;  //!<!  histo with cov mat prim vert for the 3-prong candidate
  TH1F* fHistCovMatSecVZZ3Prong;   //!<!  histo with cov mat sec vert for the 3-prong candidate

  TH1F* fHistInvMassDs;       //!<!  histo with Ds->KKpi inv mass
  TH1F* fHistInvMassDsSignal; //!<!  histo with Ds->KKpi inv mass (signal)
  TH1F* fHistInvMassDsRefl;   //!<!  histo with Ds->KKpi inv mass (reflection)
  TH1F* fHistPtDs;            //!<!  histo with Ds pt
  TH2F* fHistYPtDs;           //!<!  histo with Ds y vs pt
  TH1F* fHistPtDsDau0;        //!<!  histo with Ds prong pt
  TH1F* fHistPtDsDau1;        //!<!  histo with Ds prong pt
  TH1F* fHistPtDsDau2;        //!<!  histo with Ds prong pt
  TH1F* fHistImpParDsDau0;    //!<!  histo with Ds prong d0
  TH1F* fHistImpParDsDau1;    //!<!  histo with Ds prong d0
  TH1F* fHistImpParDsDau2;    //!<!  histo with Ds prong d0
  TH1F* fHistDecLenDs;        //!<!  histo with Ds decay length
  TH1F* fHistDecLenXYDs;      //!<!  histo with Ds decay length XY
  TH1F* fHistNormDecLenXYDs;  //!<!  histo with Ds normalized decay length XY
  TH1F* fHistCosPointDs;      //!<!  histo with Ds cosine of pointing angle
  TH1F* fHistCosPointXYDs;    //!<!  histo with Ds cosine of pointing angle XY
  TH1F* fHistImpParXYDs;      //!<!  histo with Ds impact parameter XY
  TH1F* fHistNormIPDs;        //!<!  histo with max difference between prong observed and expeceted impact parameters
  TH1F* fHistDeltaMassPhiDs;  //!<!  histo with abs. mass difference between KK pair and phi meson 
  TH1F* fHistCos3PiKDs;       //!<!  histo with cube of cosine of angle between pion and kaon
  TH1F* fHistAbsCos3PiKDs;    //!<!  histo with absolute value of cube of cosine of angle between pion and kaon

  TH1F* fHistInvMassLc;             //!<!  histo with LcpKpi+ inv mass
  TH1F* fHistPtLc;                  //!<!  histo with LcpKpi+ pt
  TH1F* fHistEtaLc;                 //!<!  histo with LcpKpi+ eta
  TH1F* fHistPhiLc;                 //!<!  histo with LcpKpi+ phi
  TH2F* fHistYPtLc;                 //!<!  histo with LcpKpi+ y vs pt
  TH1F* fHistPtLcDau0;              //!<!  histo with LcpKpi+ prong pt
  TH1F* fHistPtLcDau1;              //!<!  histo with LcpKpi+ prong pt
  TH1F* fHistPtLcDau2;              //!<!  histo with LcpKpi+ prong pt
  TH1F* fHistDecLenLc;              //!<!  histo with LcpKpi+ decay length
  TH1F* fHistDecLenXYLc;            //!<!  histo with LcpKpi+ decay length xy
  TH1F* fHistCosPointLc;            //!<!  histo with LcpKpi+ cosine of pointing angle
  TH1F* fHistCosPointXYLc;          //!<!  histo with LcpKpi+ cosine of pointing angle xy
  TH1F* fHistCtLc;                  //!<!  histo with LcpKpi+ proper decay length
  TH1F* fHistImpParLcDau0;          //!<!  histo with LcpKpi+ prong0 d0
  TH1F* fHistImpParLcDau1;          //!<!  histo with LcpKpi+ prong1 d0
  TH1F* fHistImpParLcDau2;          //!<!  histo with LcpKpi+ prong2 d0
  TH1F* fHistInvMassLcPrompt;       //!<!  histo with prompt LcpKpi+ inv mass
  TH1F* fHistEtaLcPrompt;           //!<!  histo with prompt LcpKpi+ eta
  TH1F* fHistPhiLcPrompt;           //!<!  histo with prompt LcpKpi+ phi
  TH2F* fHistYPtLcPrompt;           //!<!  histo with prompt LcpKpi+ y vs pt
  TH1F* fHistPtLcDau0Prompt;        //!<!  histo with prompt LcpKpi+ prong pt
  TH1F* fHistPtLcDau1Prompt;        //!<!  histo with prompt LcpKpi+ prong pt
  TH1F* fHistPtLcDau2Prompt;        //!<!  histo with prompt LcpKpi+ prong pt
  TH1F* fHistDecLenLcPrompt;        //!<!  histo with prompt LcpKpi+ decay length
  TH1F* fHistDecLenXYLcPrompt;      //!<!  histo with prompt LcpKpi+ decay length xy
  TH1F* fHistCtLcPrompt;            //!<!  histo with prompt LcpKpi+ proper decay length
  TH1F* fHistCosPointLcPrompt;      //!<!  histo with prompt LcpKpi+ cosine of pointing angle
  TH1F* fHistCosPointXYLcPrompt;    //!<!  histo with prompt LcpKpi+ cosine of pointing angle xy
  TH1F* fHistImpParLcDau0Prompt;    //!<!  histo with prompt LcpKpi+ prong0 d0
  TH1F* fHistImpParLcDau1Prompt;    //!<!  histo with prompt LcpKpi+ prong1 d0
  TH1F* fHistImpParLcDau2Prompt;    //!<!  histo with prompt LcpKpi+ prong2 d0
  TH1F* fHistInvMassLcNonPrompt;    //!<!  histo with nonprompt LcpKpi+ inv mass
  TH1F* fHistEtaLcNonPrompt;        //!<!  histo with nonprompt LcpKpi+ eta
  TH1F* fHistPhiLcNonPrompt;        //!<!  histo with nonprompt LcpKpi+ phi
  TH2F* fHistYPtLcNonPrompt;        //!<!  histo with nonprompt LcpKpi+ y vs pt
  TH1F* fHistPtLcDau0NonPrompt;     //!<!  histo with nonprompt LcpKpi+ prong pt
  TH1F* fHistPtLcDau1NonPrompt;     //!<!  histo with nonprompt LcpKpi+ prong pt
  TH1F* fHistPtLcDau2NonPrompt;     //!<!  histo with nonprompt LcpKpi+ prong pt
  TH1F* fHistDecLenLcNonPrompt;     //!<!  histo with nonprompt LcpKpi+ decay length
  TH1F* fHistDecLenXYLcNonPrompt;   //!<!  histo with nonprompt LcpKpi+ decay length xy
  TH1F* fHistCtLcNonPrompt;         //!<!  histo with nonprompt LcpKpi+ proper decay length
  TH1F* fHistCosPointLcNonPrompt;   //!<!  histo with nonprompt LcpKpi+ cosine of pointing angle
  TH1F* fHistCosPointXYLcNonPrompt; //!<!  histo with nonprompt LcpKpi+ cosine of pointing angle xy
  TH1F* fHistImpParLcDau0NonPrompt; //!<!  histo with nonprompt LcpKpi+ prong0 d0
  TH1F* fHistImpParLcDau1NonPrompt; //!<!  histo with nonprompt LcpKpi+ prong1 d0
  TH1F* fHistImpParLcDau2NonPrompt; //!<!  histo with nonprompt LcpKpi+ prong2 d0

  TH1F* fHistInvMassK0s;    //!<!  histo with K0s inv mass
  TH1F* fHistInvMassLcK0sp; //!<!  histo with LcpKpi+ inv mass
  TH1F* fHistPtLcK0sp;      //!<!  histo with LcpKpi+ pt

  TH1F* fHistPtGenPrompt[5];        //!<! histos for efficiency (prompt)
  TH1F* fHistPtGenFeeddw[5];        //!<! histos for efficiency (from B)
  TH1F* fHistPtGenLimAccPrompt[5];  //!<! histos for efficiency (prompt)
  TH1F* fHistPtGenLimAccFeeddw[5];  //!<! histos for efficiency (from B)
  TH1F* fHistEtaGenLimAccPrompt[5]; //!<! histos for efficiency (prompt)
  TH1F* fHistEtaGenLimAccFeeddw[5]; //!<! histos for efficiency (from B)
  TH1F* fHistPhiGenLimAccPrompt[5]; //!<! histos for efficiency (prompt)
  TH1F* fHistPhiGenLimAccFeeddw[5]; //!<! histos for efficiency (from B)
  TH1F* fHistPtRecoGenPtPrompt[5];  //!<! histos for efficiency (prompt)
  TH1F* fHistPtRecoGenPtFeeddw[5];  //!<! histos for efficiency (from B)
  TH1F* fHistPtRecoPrompt[5];       //!<! histos for efficiency (prompt)
  TH1F* fHistPtRecoFeeddw[5];       //!<! histos for efficiency (from B)

  TH2F* fHistCPUTimeTrackVsNTracks;  //!<! histo with CPU time for track selection vs number of selected tracks for candidates
  TH2F* fHistCPUTimeCandVsNTracks;   //!<! histo with CPU time for candidate selection vs number of selected tracks for candidates
  TH2F* fHistWallTimeTrackVsNTracks; //!<! histo with wall time for track selection vs number of selected tracks for candidates
  TH2F* fHistWallTimeCandVsNTracks;  //!<! histo with wall time for candidate selection vs number of selected tracks for candidates

  TH1F* fHistJetPt; //!<! histo with Jet transverse momentum
  TH1F* fHistJetE; //!<! histo with Jet energy
  TH1F* fHistJetPhi; //!<! histo with Jet phi
  TH1F* fHistJetEta; //!<! histo with Jet pseudorapidity
  TH1F* fHistJetNConstituents; //!<! histo with number of jet constituents
  TH1F* fHistJetCandPt; //!<! histo with HF candidate transverese momentum in Jet
  TH1F* fHistJetZg; //!<! histo with Zg
  TH1F* fHistJetRg; //!<! histo with Rg
  TH1F* fHistJetNsd; //!<! histo with Nsd

  TH1F* fHistJetPt_Part; //!<! histo with Jet transverse momentum at particle level
  TH1F* fHistJetE_Part; //!<! histo with Jet energy at particle level
  TH1F* fHistJetPhi_Part; //!<! histo with Jet phi at particle level
  TH1F* fHistJetEta_Part; //!<! histo with Jet pseudorapidity at particle level
  TH1F* fHistJetNConstituents_Part; //!<! histo with number of jet constituents at particle level
  TH1F* fHistJetCandPt_Part; //!<! histo with HF candidate transverese momentum in Jet at particle level
  TH1F* fHistJetZg_Part; //!<! histo with Zg at particle level
  TH1F* fHistJetRg_Part; //!<! histo with Rg at particle level
  TH1F* fHistJetNsd_Part; //!<! histo with Nsd at particle level

  Bool_t fReadMC;              // flag for access to MC
  Bool_t fUsePhysSel;          // flag use/not use phys sel
  Bool_t fUseAliEventCuts;     // flag to use/not use default AliEventCuts
  Int_t fTriggerMask;          // mask used in physics selection
  Bool_t fSelectOnCentrality;  // flag to activate cut on centrality
  Double_t fMinCentrality;     // centrality: lower limit
  Double_t fMaxCentrality;     // centrality: upper limit
  TString fCentrEstimator;     // centrality: estimator
  Bool_t fCutOnSPDVsTrackVtx;  // flag to activate cut on SPD-track vertex
  Double_t fMaxZVert;          // cut on z vertex position
  Bool_t fDo3Prong;            // flag yes/no for 3 prongs
  Double_t fMaxDecVertRadius2; // square of max radius of decay vertex

  Double_t fMassDzero;   // D0 mass from PDG
  Double_t fMassJpsi;    // Jpsi mass from PDG
  Double_t fMassDplus;   // D+ mass from PDG
  Double_t fMassDs;      // D_s mass from PDG
  Double_t fMassLambdaC; // Lc mass from PDG
  Double_t fMassK0s;     // K0s mass from PDG

  Int_t fSecVertexerAlgo;                      // Algorithm for secondary vertex
  AliVertexerTracks* fVertexerTracks;          // Run-2 vertexer
  o2::vertexing::DCAFitter2 fO2Vertexer2Prong; // o2 vertexer
  o2::vertexing::DCAFitter3 fO2Vertexer3Prong; // o2 vertexer
  Bool_t fVertexerPropagateToPCA;              // bool to propagate to PCA
  Double_t fVertexerMaxR;                      // max r vertexer
  Double_t fVertexerMaxDZIni;                  // max DZIni vertexer
  Double_t fVertexerMinParamChange;            // min ParamChange vertexer
  Double_t fVertexerMinRelChi2Change;          // min rel chi2 change vertexer
  Bool_t fVertexerUseAbsDCA;                   // use abs DCA vertexer
  AliEventCuts fEventCuts;                     // Standard AliEvent cuts
  AliESDtrackCuts* fTrackCuts2pr;              // Track cut object for 2 prongs
  AliESDtrackCuts* fTrackCuts3pr;              // Track cut object for 3 prongs
  AliESDtrackCuts* fTrackCutsBach;             // Track cut object for bachelor
  AliESDtrackCuts* fTrackCutsV0Dau;            // Track cut object for V0 daughters
  Int_t fMaxTracksToProcess;                   // Max n. of tracks, to limit test duration

  Int_t fNPtBinsSingleTrack;                                                     // Number of pt bins for single track cuts
  Double_t fPtBinLimsSingleTrack[kMaxNPtBinsDzero];                              // [fNPtBinsSingleTrack+1] limits of pt bins for single track cuts
  Double_t fSingleTrackCuts2Prong[kMaxNPtBinsSingleTrack][kNCutVarsSingleTrack]; // 2-prong single track cuts
  Double_t fSingleTrackCuts3Prong[kMaxNPtBinsSingleTrack][kNCutVarsSingleTrack]; // 3-prong single track cuts

  Double_t fPtBinLimsDzeroSkims[kMaxNPtBins2ProngsSkims]; // pt bins for D0 skim selections
  Double_t fPtBinLimsJpsiSkims[kMaxNPtBins2ProngsSkims];  // pt bins for Jpsi skim selections
  Double_t fPtBinLimsDplusSkims[kMaxNPtBins3ProngsSkims]; // pt bins for D+ skim selections
  Double_t fPtBinLimsDsSkims[kMaxNPtBins3ProngsSkims];    // pt bins for Ds+ skim selections
  Double_t fPtBinLimsLcSkims[kMaxNPtBins3ProngsSkims];    // pt bins for Lc+ skim selections
  Double_t fPtBinLimsXicSkims[kMaxNPtBins3ProngsSkims];   // pt bins for Xic+ skim selections
  Double_t fPtWithoutVtxToll;                             // pt tolerance for selections before vertex reconstruction

  Double_t fPtBinLims[kMaxNPtBinsDzero];                             // limits of pt bins
  Double_t fMinPtDzero;                                              // D0 min pt
  Double_t fMaxPtDzero;                                              // D0 max pt
  Double_t fMinPtDplus;                                              // D+ min pt
  Double_t fMaxPtDplus;                                              // D+ max pt
  Double_t fMinPtDs;                                                 // Ds min pt
  Double_t fMaxPtDs;                                                 // Ds max pt
  Double_t fMinPtJpsi;                                               // Jpsi min pt
  Double_t fMaxPtJpsi;                                               // Jpsi max pt
  Double_t fMinPtLc;                                                 // Lc min pt
  Double_t fMaxPtLc;                                                 // Lc max pt
  Int_t fCandidateCutLevel;                                          // Cuts: 0 = no, 1 = skim, 2 = analysis
  Double_t fDzeroSkimCuts[kMaxNPtBins2ProngsSkims][kNCutVars2Prong]; // D0 skimming cuts
  Double_t fJpsiSkimCuts[kMaxNPtBins2ProngsSkims][kNCutVars2Prong];  // Jpsi skimming cuts
  Double_t fDplusSkimCuts[kMaxNPtBins3ProngsSkims][kNCutVars3Prong]; // D+ skimming cuts
  Double_t fDsSkimCuts[kMaxNPtBins3ProngsSkims][kNCutVars3Prong];    // Ds skimming cuts
  Double_t fLcSkimCuts[kMaxNPtBins3ProngsSkims][kNCutVars3Prong];    // Lc skimming cuts
  Double_t fXicSkimCuts[kMaxNPtBins3ProngsSkims][kNCutVars3Prong];   // Xic skimming cuts
  Double_t fDzeroCuts[kMaxNPtBinsDzero][kNCutVarsDzero];             // D0 cuts
  Double_t fDplusCuts[kMaxNPtBinsDplus][kNCutVarsDplus];             // D+ cuts
  Double_t fDsCuts[kMaxNPtBinsDs][kNCutVarsDs];                      // Ds cuts
  Double_t fJpsiCuts[kMaxNPtBinsJpsi][kNCutVarsJpsi];                // Jpsi cuts
  Double_t fLcCuts[kMaxNPtBinsLc][kNCutVarsLc];                      // LcpKpi+ cuts
  Double_t fMinPt3Prong;                                             // Min pt for 3 prong candidate
  Double_t fMaxRapidityCand;                                         // Max rapidity cut (if -999 use pt dependent cut)
  static constexpr Int_t fNPtBinsDzero = 25;                         // Number of pt bins Dzero
  static constexpr Int_t fNPtBinsDplus = 12;                         // Number of pt bins Dplus
  static constexpr Int_t fNPtBinsDs = 8;                            // Number of pt bins Ds
  static constexpr Int_t fNPtBinsJpsi = 9;                           // Number of pt bins Jpsi
  static constexpr Int_t fNPtBinsLc = 10;                            // Number of pt bins Lc
  Double_t fPtBinLimsDzero[kMaxNPtBinsDzero];                        // [fNPtBinsDzero+1] limits of pt bins
  Double_t fPtBinLimsDplus[kMaxNPtBinsDplus];                        // [fNPtBinsDplus+1] limits of pt bins
  Double_t fPtBinLimsDs[kMaxNPtBinsDs];                              // [fNPtBinsDs+1] limits of pt bins
  Double_t fPtBinLimsJpsi[kMaxNPtBinsJpsi];                          // [fNPtBinsJpsi+1] limits of pt bins
  Double_t fPtBinLimsLc[kMaxNPtBinsLc];                              // [fNPtBinsLc+1] limits of pt bins
  Int_t fSelectD0;                                                   // flag to activate cuts for D0
  Int_t fSelectD0bar;                                                // flag to activate cuts for D0bar
  Int_t fSelectDplus;                                                // flag to activate cuts for Dplus
  Int_t fSelectDs;                                                   // flag to activate cuts for Ds
  Int_t fSelectJpsi;                                                 // flag to activate cuts for Jpsi
  Int_t fSelectLcpKpi;                                               // flag to activate cuts for LcpKpi
  Bool_t fFindVertexForCascades;                                     // flag to activate vertex reco for V0+bachelor
  Double_t fMinPtV0;                                                 // minimum V0 pt for Lc->pK0s
  Double_t fMinCosPointV0;                                           // K0s cuts for Lc->pK0s
  Double_t fCutOnK0sMass;                                            // K0s cuts for Lc->pK0s

  Bool_t fEnableCPUTimeCheck;      // flag to enable CPU time benchmark
  Bool_t fCountTimeInMilliseconds; // flag to switch from seconds (default) to milliseconds
  
  Bool_t doJetFinding;
  Bool_t matchedJetsOnly;
  Bool_t doJetSubstructure;
  AliFJWrapper* fFastJetWrapper;
  AliFJWrapper* fFastJetWrapper_Part;

  Double_t jetPtMin;                         // minimum jet pT
  Double_t jetPtMax;                         // maximum jet pT
  Double_t jetR;                             // jet resolution parameter
  Double_t jetTrackPtMin;                   // minimum track pT for jet finding
  Double_t jetTrackPtMax;                   // maximum track pT for jet finding
  Double_t jetTrackEtaMin;                  // maximum track eta for jet finding
  Double_t jetTrackEtaMax;                  // minimum track eta for jet finding
  Double_t jetCandPtMin;                   // maximum candidate pT for jet finding
  Double_t jetCandPtMax;                   // minimum candidate pT for jet find
  Double_t jetCandYMin;                    // maximum candidate rapidity for jet find
  Double_t jetCandYMax;                    // minimum candidate rapidity for jet find

  Double_t jetPtMin_Part;
  Double_t jetPtMax_Part;
  Double_t jetR_Part;
  Double_t jetTrackPtMin_Part;
  Double_t jetTrackPtMax_Part;
  Double_t jetTrackEtaMin_Part;
  Double_t jetTrackEtaMax_Part;
  Double_t jetCandPtMin_Part;
  Double_t jetCandPtMax_Part;
  Double_t jetCandYMin_Part;
  Double_t jetCandYMax_Part;

  Double_t zCut;
  Double_t beta;
  Double_t zCut_Part;
  Double_t beta_Part;

  ClassDef(AliAnalysisTaskHFSimpleVertices, 29);
};

#endif
