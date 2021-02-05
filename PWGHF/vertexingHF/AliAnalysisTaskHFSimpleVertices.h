#ifndef ALIANALYSISTASKHFSIMPLEVERTICES_H
#define ALIANALYSISTASKHFSIMPLEVERTICES_H

/* Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//          
//*************************************************************************

#include "DCAFitterN.h"

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
  void SetReadMC(Bool_t read){fReadMC=read;}
  void SetUseCandidateAnalysisCuts(){fCandidateCutLevel=2;}
  void SetUseCandidateSkimCuts(){fCandidateCutLevel=1;}
  void SetUseNoCandidateCuts(){fCandidateCutLevel=0;}
  
 private:

  AliAnalysisTaskHFSimpleVertices(const AliAnalysisTaskHFSimpleVertices &source);
  AliAnalysisTaskHFSimpleVertices& operator=(const AliAnalysisTaskHFSimpleVertices &source);


  char* GetJsonString(const char* jsonFileName, const char* key);
  int GetJsonInteger(const char* jsonFileName, const char* key);
  bool GetJsonBool(const char* jsonFileName, const char* key);
  float GetJsonFloat(const char* jsonFileName, const char* key);
  
  void InitDefault();
  Int_t GetPtBin(Double_t ptCand);
  void ProcessTriplet(TObjArray* threeTrackArray, AliAODRecoDecay* rd4massCalc3, AliESDVertex* primVtxTrk, AliAODVertex *vertexAODp, float bzkG, double dist12);
  Bool_t GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG);
  Int_t SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG);
  AliESDVertex* ReconstructSecondaryVertex(TObjArray* trkArray, AliESDVertex* primvtx);
  AliAODVertex* ConvertToAODVertex(AliESDVertex* trkv);
  Int_t SelectInvMassAndPt2prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc2);
  Int_t SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3);
  AliAODRecoDecayHF2Prong* Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG);
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG);


  Int_t DzeroSkimCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t JpsiSkimCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t DplusSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t DsSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t LcSkimCuts(AliAODRecoDecayHF3Prong* cand);
  Int_t DzeroSelectionCuts(AliAODRecoDecayHF2Prong* cand);
  Int_t JpsiSelectionCuts(AliAODRecoDecayHF2Prong* cand,AliESDtrack* trk_p,AliESDtrack* trk_n,AliESDVertex* primvtx,float bzkG);
  Int_t LcSelectionCuts(AliAODRecoDecayHF3Prong *cand);
  Int_t MatchToMC(AliAODRecoDecay* rd, Int_t pdgabs, AliMCEvent* mcEvent,Int_t ndgCk, const TObjArray *trkArray, const Int_t *pdgDg) const;
  
  enum ESelBits2prong {kbitDzero = 0,kbitDzerobar,kbitJpsi};
  enum ESelBits3prong {kbitDplus = 0,kbitDs,kbitLc};
  enum {kMaxNPtBins = 100, kNCutVarsDzero=11};
  enum {kMaxNPtBinsJpsi = 9, kNCutVarsJpsi=4};
  enum { kMaxNPtBinsLc = 10, kNCutVarsLc = 8 };

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
  TH1F* fHist2ProngVertX;            //!<!  histo of D0 vertex x
  TH1F* fHist2ProngVertY;            //!<!  histo of D0 vertex y
  TH1F* fHist2ProngVertZ;            //!<!  histo of D0 vertex z
  TH1F* fHistDplusVertX;             //!<!  histo of D+ vertex x
  TH1F* fHistDplusVertY;             //!<!  histo of D+ vertex y
  TH1F* fHistDplusVertZ;             //!<!  histo of D+ vertex z
  TH1F *fHistLcpKpiVertX;            //!<!  histo of LcpKpi+ vertex x
  TH1F *fHistLcpKpiVertY;            //!<!  histo of LcpKpi+ vertex y
  TH1F *fHistLcpKpiVertZ;            //!<!  histo of LcpKpi+ vertex z
  TH1F *fHistDist12LcpKpi;           //!<!  histo of LcpKpi+ distance between primary and secondary vertex reconstructed from the pair of tracks

  TH1F* fHistInvMassD0;              //!<!  histo with D0 inv mass
  TH1F* fHistPtD0;                   //!<!  histo with D0 pt
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
  TH1F* fHistPtDPlus;                //!<!  histo with D+ pt
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
  TH1F *fHistPtDs;                   //!<!  histo with Ds pt
  TH1F *fHistDecLenDs;               //!<!  histo with Ds decay length
  TH1F *fHistCosPointDs;             //!<!  histo with Ds cosine of pointing angle

  TH1F *fHistInvMassLc;              //!<!  histo with LcpKpi+ inv mass
  TH1F *fHistPtLc;                   //!<!  histo with LcpKpi+ pt
  TH1F *fHistPtLcDau0;               //!<!  histo with LcpKpi+ prong pt
  TH1F *fHistPtLcDau1;               //!<!  histo with LcpKpi+ prong pt
  TH1F *fHistPtLcDau2;               //!<!  histo with LcpKpi+ prong pt
  TH1F *fHistDecLenLc;               //!<!  histo with LcpKpi+ decay length
  TH1F *fHistCosPointLc;             //!<!  histo with LcpKpi+ cosine of pointing angle

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

  AliESDtrackCuts* fTrackCuts2pr;  // Track cut object for 2 prongs
  AliESDtrackCuts* fTrackCuts3pr;  // Track cut object for 3 prongs
  Int_t fMaxTracksToProcess;       // Max n. of tracks, to limit test duration


  Int_t fNPtBins;                     // Number of pt bins
  Double_t fPtBinLims[kMaxNPtBins];   // [fNPtBins+1] limits of pt bins
  Double_t fMinPtDzero;               // D0 min pt
  Double_t fMaxPtDzero;               // D0 max pt
  Double_t fMinPtJpsi;                // Jpsi min pt
  Double_t fMaxPtJpsi;                // Jpsi max pt
  Int_t    fCandidateCutLevel;        // Cuts: 0 = no, 1 = skim, 2 = analysis
  Double_t fDzeroSkimCuts[5];         // D0 skimming cuts
  Double_t fJpsiSkimCuts[5];          // Jpsi skimming cuts
  Double_t fDplusSkimCuts[5];         // D0 skimming cuts
  Double_t fDsSkimCuts[5];            // D0 skimming cuts
  Double_t fLcSkimCuts[5];            // D0 skimming cuts
  Double_t fDzeroCuts[kMaxNPtBins][kNCutVarsDzero]; // D0 cuts
  Double_t fJpsiCuts[kMaxNPtBinsJpsi][kNCutVarsJpsi]; // Jpsi cuts
  Int_t fSelectD0;                    // flag to activate cuts for D0
  Int_t fSelectD0bar;                 // flag to activate cuts for D0bar
  Double_t fMinPt3Prong;              // Min pt for 3 prong candidate

  Int_t fNPtBinsJpsi;
  Int_t fNPtBinsLc;                             // Number of pt bins
  Int_t fSelectJpsi;
  Double_t fPtBinLimsLc[kMaxNPtBinsLc];         // [fNPtBins+1] limits of pt bins
  Double_t fPtBinLimsJpsi[kMaxNPtBinsJpsi];
  Double_t fLcCuts[kMaxNPtBinsLc][kNCutVarsLc]; // LcpKpi+ cuts
  Int_t fSelectLcpKpi;                          // flag to activate cuts for LcpKpi

  ClassDef(AliAnalysisTaskHFSimpleVertices,12);
};


#endif
