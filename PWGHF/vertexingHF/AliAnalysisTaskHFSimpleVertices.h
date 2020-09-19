#ifndef ALIANALYSISTASKHFSIMPLEVERTICES_H
#define ALIANALYSISTASKHFSIMPLEVERTICES_H

/* Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//          
//*************************************************************************

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

 private:

  AliAnalysisTaskHFSimpleVertices(const AliAnalysisTaskHFSimpleVertices &source);
  AliAnalysisTaskHFSimpleVertices& operator=(const AliAnalysisTaskHFSimpleVertices &source);


  char* GetJsonString(const char* jsonFileName, const char* key);
  int GetJsonInteger(const char* jsonFileName, const char* key);
  bool GetJsonBool(const char* jsonFileName, const char* key);
  float GetJsonFloat(const char* jsonFileName, const char* key);
  
  void InitDefault();
  Bool_t GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG);
  Bool_t SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG);
  AliESDVertex* ReconstructSecondaryVertex(AliVertexerTracks* vt, TObjArray* trkArray, AliESDVertex* primvtx);
  AliAODVertex* ConvertToAODVertex(AliESDVertex* trkv);
  Int_t SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3);
  AliAODRecoDecayHF2Prong* Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG);
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG);

  enum ESelBits3prong {kbitDplus = 0,kbitDs,kbitLc};
  
  TList*  fOutput;                   //!<!  list of output histos
  TH1F* fHistNEvents;                //!<!  histo with N of events
  TH1F* fHistPtAllTracks;            //!<!  histo with pt all tracks
  TH1F* fHistPtSelTracks;            //!<!  histo with pt selected tracks
  TH1F* fHistTglAllTracks;           //!<!  histo with tgl all tracks
  TH1F* fHistTglSelTracks;           //!<!  histo with tgl selected tracks
  TH1F* fHistImpParAllTracks;        //!<!  histo with d0 all tracks
  TH1F* fHistImpParSelTracks;        //!<!  histo with d0 selected tracks
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

  TH1F* fHistInvMassD0;              //!<!  histo with D0 inv mass
  TH1F* fHistPtD0;                   //!<!  histo with D0 pt
  TH1F* fHistPtD0Dau0;               //!<!  histo with D0 prong pt
  TH1F* fHistPtD0Dau1;               //!<!  histo with D0 prong pt
  TH1F* fHistImpParD0Dau0;           //!<!  histo with D0 prong d0
  TH1F* fHistImpParD0Dau1;           //!<!  histo with D0 prong d0
  TH1F* fHistd0Timesd0;              //!<!  histo with d0xd0
  TH1F* fHistDecLenD0;               //!<!  histo with D0 decay length
  TH1F* fHistDecLenXYD0;             //!<!  histo with D0 decay length XY
  TH1F* fHistInvMassDplus;           //!<!  histo with D+ inv mass

  Bool_t  fUsePhysSel;         // flag use/not use phys sel
  Int_t   fTriggerMask;        // mask used in physics selection
  Bool_t  fSelectOnCentrality; // flag to activate cut on centrality
  Double_t fMinCentrality;     // centrality: lower limit
  Double_t fMaxCentrality;     // centrality: upper limit
  TString fCentrEstimator;     // centrality: estimator
  Bool_t  fDo3Prong;           // flag yes/no for 3 prongs
  
  Double_t fMassDzero;         // D0 mass from PDG
  Double_t fMassDplus;         // D+ mass from PDG
  Double_t fMassDs;            // D_s mass from PDG
  Double_t fMassLambdaC;       // Lc mass from PDG

  AliESDtrackCuts* fTrackCuts;  // Track cut object
  Int_t fMaxTracksToProcess;    // Max n. of tracks, to limit test duration
  
  ClassDef(AliAnalysisTaskHFSimpleVertices,2);
};


#endif
