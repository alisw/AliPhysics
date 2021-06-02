#ifndef ALIANALYSISTASKHFSIMPLEVERTICES_H
#define ALIANALYSISTASKHFSIMPLEVERTICES_H

/* Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//          
//*************************************************************************

#include <map>
#include <string>
#include "DCAFitterN.h"
#include "AliVEvent.h"

class TList;
class AliESDEvent;

class AliAnalysisTaskHFSimpleVertices : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskHFSimpleVertices();
  virtual ~AliAnalysisTaskHFSimpleVertices();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);
  void InitFromJson(TString filename);
  void SetUseCutOnSPDVsTrackVtx(Bool_t opt){fCutOnSPDVsTrackVtx=opt;}
  void SetZVertexMaxRange(Double_t zmax){fMaxZVert=zmax;}
  void SetUseVertexerTracks(){fSecVertexerAlgo=0;}
  void SetUseO2Vertexer(){fSecVertexerAlgo=1;}
  void SetUseKFParticleVertexer(){fSecVertexerAlgo=2;}
  void SetReadMC(Bool_t read){fReadMC=read;}
  void SetUseCandidateAnalysisCuts(){fCandidateCutLevel=2;}
  void SetUseCandidateSkimCuts(){fCandidateCutLevel=1;}
  void SetUseNoCandidateCuts(){fCandidateCutLevel=0;}
  void SetUsePtDependentFiducialAcceptance(){fMaxRapidityCand=-999.;}
  void SetMaxRapidityForFiducialAcceptance(Double_t ymax){fMaxRapidityCand=ymax;}
  void EnableCPUTimeCheck(Bool_t enable=kTRUE, Bool_t milliseconds=kFALSE) {fEnableCPUTimeCheck=enable; fCountTimeInMilliseconds=milliseconds;}

 private:

  AliAnalysisTaskHFSimpleVertices(const AliAnalysisTaskHFSimpleVertices &source);
  AliAnalysisTaskHFSimpleVertices& operator=(const AliAnalysisTaskHFSimpleVertices &source);

  std::string GetJsonString(const char* jsonFileName, const char* key);
  int GetJsonInteger(const char* jsonFileName, const char* key);
  int GetJsonBool(const char* jsonFileName, const char* key);
  float GetJsonFloat(const char* jsonFileName, const char* key);
  float* GetJsonArray(const char* jsonFileName, const char* key, int &size);
  float** GetJsonMatrix(const char* jsonFileName, const char* key, int &size1, int &size2);
  void InitDefault();
  Int_t GetPtBin(Double_t ptCand);
  Int_t GetPtBinSingleTrack(Double_t ptTrack);
  void ProcessTriplet(TObjArray* threeTrackArray, AliAODRecoDecay* rd4massCalc3, AliESDVertex* primVtxTrk, AliAODVertex *vertexAODp, float bzkG, double dist12, AliMCEvent* mcEvent);
  Bool_t GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG);
  Int_t SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG, Double_t d0track[2]);
  AliESDVertex* ReconstructSecondaryVertex(TObjArray* trkArray, AliESDVertex* primvtx);
  AliAODVertex* ConvertToAODVertex(AliESDVertex* trkv);
  Int_t SelectInvMassAndPt2prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc2);
  Int_t SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3);
  AliAODRecoDecayHF2Prong* Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG);
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG);
  AliAODRecoCascadeHF* MakeCascade(TObjArray *twoTrackArray, AliAODVertex* secVert, Double_t bzkG);
  Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y) const;
  Int_t DzeroSkimCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t JpsiSkimCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t DplusSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t DsSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t LcSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t DzeroSelectionCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t DplusSelectionCuts(AliAODRecoDecayHF3Prong *cand, Double_t bzkG);
  Int_t JpsiSelectionCuts(AliAODRecoDecayHF2Prong* cand,AliESDtrack* trk_p,AliESDtrack* trk_n,AliESDVertex* primvtx,float bzkG);
  Int_t LcSelectionCuts(AliAODRecoDecayHF3Prong *cand);
  Int_t MatchToMC(AliAODRecoDecay* rd, Int_t pdgabs, AliMCEvent* mcEvent,Int_t ndgCk, const TObjArray *trkArray, const Int_t *pdgDg) const;
  
  enum ESelBits2prong {kbitDzero = 0,kbitDzerobar,kbitJpsi};
  enum ESelBits3prong {kbitDplus = 0,kbitDs,kbitLc};
  enum {kMaxNPtBins = 100, kNCutVarsDzero=11};
  enum {kMaxNPtBinsJpsi = 9, kNCutVarsJpsi=4};
  enum {kMaxNPtBinsLc = 10, kNCutVarsLc = 8 };
  enum {kMaxNPtBinsDplus = 50, kNCutVarsDplus = 8};
  enum {kMaxNPtBinsSingleTrack = 10, kNCutVarsSingleTrack = 5};

  /// list of triggers avaliable in O2 that are also avaliable in AliPhysics
  std::map<std::string, Int_t> triggerMask = {{"kINT7", AliVEvent::kINT7},
                                              {"kEMC7", AliVEvent::kEMC7},
                                              {"kINT7inMUON", AliVEvent::kINT7inMUON},
                                              {"kMuonSingleLowPt7", AliVEvent::kMuonSingleLowPt7},
                                              {"kMuonSingleHighPt7", AliVEvent::kMuonSingleHighPt7},
                                              {"kMuonUnlikeLowPt7", AliVEvent::kMuonUnlikeLowPt7},
                                              {"kMuonLikeLowPt7", AliVEvent::kMuonLikeLowPt7}};

  TList*  fOutput;                   //!<!  list of output histos
  TH1F* fHistNEvents;                //!<!  histo with N of events
  TH1F* fHistPtAllTracks;            //!<!  histo with pt all tracks
  TH1F* fHistPtSelTracks;            //!<!  histo with pt selected tracks
  TH1F* fHistTglAllTracks;           //!<!  histo with tgl all tracks
  TH1F* fHistTglSelTracks;           //!<!  histo with tgl selected tracks
  TH1F* fHistEtaAllTracks;           //!<!  histo with eta all tracks
  TH1F* fHistEtaSelTracks2prong;     //!<!  histo with eta selected tracks
  TH1F* fHistEtaSelTracks3prong;     //!<!  histo with eta selected tracks
  TH1F* fHistImpParAllTracks;        //!<!  histo with d0 all tracks
  TH1F* fHistImpParSelTracks2prong;  //!<!  histo with d0 selected tracks
  TH1F* fHistImpParSelTracks3prong;  //!<!  histo with d0 selected tracks
  TH1F* fHistITSmapAllTracks;        //!<!  histo with its map all tracks
  TH1F* fHistITSmapSelTracks;        //!<!  histo withits map selected tracks
  
  TH1F* fHistPrimVertX;              //!<!  histo of prim vertex x
  TH1F* fHistPrimVertY;              //!<!  histo of prim vertex y
  TH1F* fHistPrimVertZ;              //!<!  histo of prim vertex z
  TH1F* fHist2ProngVertX;            //!<!  histo of 2-prong vertex x
  TH1F* fHist2ProngVertY;            //!<!  histo of 2-prong vertex y
  TH1F* fHist2ProngVertZ;            //!<!  histo of 2-prong vertex z
  TH1F* fHist3ProngVertX;           //!<!  histo of 3-prong vertex x
  TH1F* fHist3ProngVertY;           //!<!  histo of 3-prong vertex y
  TH1F* fHist3ProngVertZ;           //!<!  histo of 3-prong vertex z
  TH1F *fHistDist12LcpKpi;           //!<!  histo of LcpKpi+ distance between primary and secondary vertex reconstructed from the pair of tracks

  TH1F* fHistInvMassD0;              //!<!  histo with D0 inv mass
  TH1F* fHistPtD0;                   //!<!  histo with D0 pt
  TH2F *fHistYPtD0;                  //!<!  histo with D0 y vs pt
  TH1F* fHistPtD0Dau0;               //!<!  histo with D0 prong pt
  TH1F* fHistPtD0Dau1;               //!<!  histo with D0 prong pt
  TH1F* fHistImpParD0Dau0;           //!<!  histo with D0 prong d0
  TH1F* fHistImpParD0Dau1;           //!<!  histo with D0 prong d0
  TH1F* fHistd0Timesd0;              //!<!  histo with d0xd0
  TH1F* fHistCosPointD0;             //!<!  histo with D0 cosine of pointing angle
  TH1F* fHistDecLenD0;               //!<!  histo with D0 decay length
  TH1F* fHistDecLenXYD0;             //!<!  histo with D0 decay length XY
  TH1F* fHistImpParErrD0Dau;         //!<!  histo with D0 prong d0 err
  TH1F* fHistDecLenErrD0;            //!<!  histo with D0 decay length err
  TH1F* fHistDecLenXYErrD0;          //!<!  histo with D0 decay length XY err
  TH1F* fHistCovMatPrimVXX2Prong;    //!<!  histo with cov mat prim vert for the 2-prong candidate
  TH1F* fHistCovMatSecVXX2Prong;     //!<!  histo with cov mat sec vert for the 2-prong candidate
  TH1F* fHistD0SignalVertX;          //!<!  histo of D0 (MC truth) vertex x
  TH1F* fHistD0SignalVertY;          //!<!  histo of D0 (MC truth) vertex y
  TH1F* fHistD0SignalVertZ;          //!<!  histo of D0 (MC truth) vertex z
  TH1F* fHistInvMassD0Signal;        //!<!  histo with D0 (MC truth) inv mass
  TH1F* fHistInvMassD0Refl;          //!<!  histo with D0 (reflection) inv mass
 
  TH1F* fHistInvMassJpsi;              //!<!  histo with Jpsi inv mass
  TH1F* fHistPtJpsi;                   //!<!  histo with Jpsi pt
  TH1F* fHistPtJpsiDau0;               //!<!  histo with Jpsi prong pt
  TH1F* fHistPtJpsiDau1;
  TH1F* fHistImpParJpsiDau0;           //!<!  histo with Jpsi prong d0
  TH1F* fHistImpParJpsiDau1;
  TH1F* fHistd0Timesd0Jpsi;              //!<!  histo with d0xd0
  TH1F* fHistCosPointJpsi;             //!<!  histo with Jpsi cosine of pointing angle
  TH1F* fHistDecLenJpsi;               //!<!  histo with Jpsi decay length
  TH1F* fHistDecLenXYJpsi;             //!<!  histo with Jpsi decay length XY
  TH1F* fHistDecLenErrJpsi;            //!<!  histo with Jpsi decay length err
  TH1F* fHistDecLenXYErrJpsi;          //!<!  histo with Jpsi decay length XY err
  TH1F* fHistJpsiSignalVertX;          //!<!  histo of Jpsi (MC truth) vertex x
  TH1F* fHistJpsiSignalVertY;          //!<!  histo of Jpsi (MC truth) vertex y
  TH1F* fHistJpsiSignalVertZ;          //!<!  histo of Jpsi (MC truth) vertex z
  TH1F* fHistInvMassJpsiSignal;        //!<!  histo with Jpsi (MC truth) inv mass

  TH1F* fHistInvMassDplus;           //!<!  histo with D+ inv mass
  TH1F* fHistInvMassDplusSignal;     //!<!  histo with D+ inv mass (only signal)
  TH1F* fHistPtDplus;                //!<!  histo with D+ pt
  TH2F *fHistYPtDplus;               //!<!  histo with D+ y vs pt
  TH1F* fHistPtDplusDau0;            //!<!  histo with D+ prong pt
  TH1F* fHistPtDplusDau1;            //!<!  histo with D+ prong pt
  TH1F* fHistPtDplusDau2;            //!<!  histo with D+ prong pt
  TH1F* fHistImpParDplusDau0;        //!<!  histo with D+ prong d0
  TH1F* fHistImpParDplusDau1;        //!<!  histo with D+ prong d0
  TH1F* fHistImpParDplusDau2;        //!<!  histo with D+ prong d0
  TH1F* fHistDecLenDplus;            //!<!  histo with D+ decay length
  TH1F* fHistDecLenXYDplus;          //!<!  histo with D+ decay length XY
  TH1F* fHistNormDecLenXYDplus;      //!<!  histo with D+ normalized decay length XY
  TH1F* fHistImpParErrDplusDau;      //!<!  histo with D+ prong d0 err
  TH1F* fHistDecLenErrDplus;         //!<!  histo with D+ decay length err
  TH1F* fHistDecLenXYErrDplus;       //!<!  histo with D+ decay length XY err
  TH1F* fHistCosPointDplus;          //!<!  histo with D+ cosine of pointing angle
  TH1F* fHistCosPointXYDplus;        //!<!  histo with D+ cosine of pointing angle XY
  TH1F* fHistImpParXYDplus;          //!<!  histo with D+ impact parameter XY
  TH1F* fHistNormIPDplus;            //!<!  histo with max difference between prong observed and expeceted impact parameters
  TH1F* fHistoSumSqImpParDplusDau;   //!<!  histo with squared sum of prong impact parameters
  TH1F* fHistCovMatPrimVXX3Prong;    //!<!  histo with cov mat prim vert for the 3-prong candidate
  TH1F* fHistCovMatSecVXX3Prong;     //!<!  histo with cov mat sec vert for the 3-prong candidate

  TH1F *fHistInvMassDs;              //!<!  histo with Ds->KKpi inv mass
  TH1F *fHistInvMassDsSignal;        //!<!  histo with Ds->KKpi inv mass (signal)
  TH1F *fHistInvMassDsRefl;          //!<!  histo with Ds->KKpi inv mass (reflection)
  TH1F *fHistPtDs;                   //!<!  histo with Ds pt
  TH2F *fHistYPtDs;                  //!<!  histo with Ds y vs pt
  TH1F *fHistDecLenDs;               //!<!  histo with Ds decay length
  TH1F *fHistCosPointDs;             //!<!  histo with Ds cosine of pointing angle

  TH1F *fHistInvMassLc;              //!<!  histo with LcpKpi+ inv mass
  TH1F *fHistPtLc;                   //!<!  histo with LcpKpi+ pt
  TH2F *fHistYPtLc;                  //!<!  histo with LcpKpi+ y vs pt
  TH1F *fHistPtLcDau0;               //!<!  histo with LcpKpi+ prong pt
  TH1F *fHistPtLcDau1;               //!<!  histo with LcpKpi+ prong pt
  TH1F *fHistPtLcDau2;               //!<!  histo with LcpKpi+ prong pt
  TH1F *fHistDecLenLc;               //!<!  histo with LcpKpi+ decay length
  TH1F *fHistCosPointLc;             //!<!  histo with LcpKpi+ cosine of pointing angle

  TH1F *fHistInvMassLcK0sp;         //!<!  histo with LcpKpi+ inv mass
  TH1F *fHistPtLcK0sp;              //!<!  histo with LcpKpi+ pt

  TH1F* fHistPtGenPrompt[5];        //!<! histos for efficiency (prompt)
  TH1F* fHistPtGenFeeddw[5];        //!<! histos for efficiency (from B)
  TH1F* fHistPtGenLimAccPrompt[5];  //!<! histos for efficiency (prompt)
  TH1F* fHistPtGenLimAccFeeddw[5];  //!<! histos for efficiency (from B)
  TH1F* fHistPtRecoGenPtPrompt[5];  //!<! histos for efficiency (prompt)
  TH1F* fHistPtRecoGenPtFeeddw[5];  //!<! histos for efficiency (from B)
  TH1F* fHistPtRecoPrompt[5];       //!<! histos for efficiency (prompt)
  TH1F* fHistPtRecoFeeddw[5];       //!<! histos for efficiency (from B)
  
  TH2F* fHistCPUTimeTrackVsNTracks;  //!<! histo with CPU time for track selection vs number of selected tracks for candidates
  TH2F* fHistCPUTimeCandVsNTracks;   //!<! histo with CPU time for candidate selection vs number of selected tracks for candidates
  TH2F* fHistWallTimeTrackVsNTracks; //!<! histo with wall time for track selection vs number of selected tracks for candidates
  TH2F* fHistWallTimeCandVsNTracks;  //!<! histo with wall time for candidate selection vs number of selected tracks for candidates

  Bool_t  fReadMC;             // flag for access to MC
  Bool_t  fUsePhysSel;         // flag use/not use phys sel
  Int_t   fTriggerMask;        // mask used in physics selection
  Bool_t  fSelectOnCentrality; // flag to activate cut on centrality
  Double_t fMinCentrality;     // centrality: lower limit
  Double_t fMaxCentrality;     // centrality: upper limit
  TString fCentrEstimator;     // centrality: estimator
  Bool_t fCutOnSPDVsTrackVtx;  // flag to activate cut on SPD-track vertex
  Double_t fMaxZVert;          // cut on z vertex position
  Bool_t  fDo3Prong;           // flag yes/no for 3 prongs
  Double_t fMaxDecVertRadius2; // square of max radius of decay vertex
  
  Double_t fMassDzero;         // D0 mass from PDG
  Double_t fMassJpsi;          // Jpsi mass from PDG
  Double_t fMassDplus;         // D+ mass from PDG
  Double_t fMassDs;            // D_s mass from PDG
  Double_t fMassLambdaC;       // Lc mass from PDG

  Int_t fSecVertexerAlgo;                  // Algorithm for secondary vertex
  AliVertexerTracks* fVertexerTracks;             // Run-2 vertexer
  o2::vertexing::DCAFitter2 fO2Vertexer2Prong;    // o2 vertexer
  o2::vertexing::DCAFitter3 fO2Vertexer3Prong;    // o2 vertexer
  Bool_t fVertexerPropagateToPCA;
  Double_t fVertexerMaxR;
  Double_t fVertexerMaxDZIni;
  Double_t fVertexerMinParamChange;
  Double_t fVertexerMinRelChi2Change;
  Bool_t fVertexerUseAbsDCA;
  AliESDtrackCuts* fTrackCuts2pr;  // Track cut object for 2 prongs
  AliESDtrackCuts* fTrackCuts3pr;  // Track cut object for 3 prongs
  Int_t fMaxTracksToProcess;       // Max n. of tracks, to limit test duration

  Int_t fNPtBinsSingleTrack;   // Number of pt bins for single track cuts
  Double_t fPtBinLimsSingleTrack[kMaxNPtBins];   // [fNPtBinsSingleTrack+1] limits of pt bins for single track cuts
  Double_t fSingleTrackCuts2Prong[kMaxNPtBinsSingleTrack][kNCutVarsSingleTrack]; // 2-prong single track cuts
  Double_t fSingleTrackCuts3Prong[kMaxNPtBinsSingleTrack][kNCutVarsSingleTrack]; // 3-prong single track cuts

  Int_t fNPtBins;                     // Number of pt bins
  Double_t fPtBinLims[kMaxNPtBins];   // [fNPtBins+1] limits of pt bins
  Double_t fMinPtDzero;               // D0 min pt
  Double_t fMaxPtDzero;               // D0 max pt
  Double_t fMinPtDplus;               // D+ min pt
  Double_t fMaxPtDplus;               // D+ max pt
  Double_t fMinPtJpsi;                // Jpsi min pt
  Double_t fMaxPtJpsi;                // Jpsi max pt
  Int_t    fCandidateCutLevel;        // Cuts: 0 = no, 1 = skim, 2 = analysis
  Double_t fDzeroSkimCuts[5];         // D0 skimming cuts
  Double_t fJpsiSkimCuts[5];          // Jpsi skimming cuts
  Double_t fDplusSkimCuts[5];         // D0 skimming cuts
  Double_t fDsSkimCuts[5];            // D0 skimming cuts
  Double_t fLcSkimCuts[5];            // D0 skimming cuts
  Double_t fDzeroCuts[kMaxNPtBins][kNCutVarsDzero]; // D0 cuts
  Double_t fDplusCuts[kMaxNPtBinsDplus][kNCutVarsDplus]; // D+ cuts
  Double_t fJpsiCuts[kMaxNPtBinsJpsi][kNCutVarsJpsi]; // Jpsi cuts
  Int_t fSelectD0;                    // flag to activate cuts for D0
  Int_t fSelectD0bar;                 // flag to activate cuts for D0bar
  Double_t fMinPt3Prong;              // Min pt for 3 prong candidate
  Double_t fMaxRapidityCand;          // Max rapidity cut (if -999 use pt dependent cut)
  Int_t fNPtBinsDplus;                          // Number of pt bins Dplus
  Int_t fNPtBinsJpsi;                           // Number of pt bins Jpsi
  Int_t fNPtBinsLc;                             // Number of pt bins Lc
  Int_t fSelectDplus;                           // flag to activate cuts for Dplus
  Int_t fSelectJpsi;                            // flag to activate cuts for Jpsi
  Double_t fPtBinLimsDplus[kMaxNPtBinsDplus];   // [fNPtBinsDplus+1] limits of pt bins
  Double_t fPtBinLimsLc[kMaxNPtBinsLc];         // [fNPtBinsLc+1] limits of pt bins
  Double_t fPtBinLimsJpsi[kMaxNPtBinsJpsi];     // [fNPtBinsJpsi+1] limits of pt bins
  Double_t fLcCuts[kMaxNPtBinsLc][kNCutVarsLc]; // LcpKpi+ cuts
  Int_t fSelectLcpKpi;                          // flag to activate cuts for LcpKpi
  Double_t fMinPtV0;                            // minimum V0 pt for Lc->pK0s

  Bool_t fEnableCPUTimeCheck;                   //flag to enable CPU time benchmark
  Bool_t fCountTimeInMilliseconds;              // flag to switch from seconds (default) to milliseconds
  
  ClassDef(AliAnalysisTaskHFSimpleVertices,21);
};


#endif
