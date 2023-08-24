#ifndef ALIJHFTAGGING_H
#define ALIJHFTAGGING_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliAnalysisTaskEmcalJet.h"

class TRandom3;
class AliEmcalJet;
class AliHFJetsTagging;
class AliVertexerTracks;
class AliFJWrapper;
class AliJetContainer;
class AliHFJetsTaggingVertex;
class AliRDHFJetsCutsVertex;
//==============================================================

class AliJHFTagging : public AliAnalysisTaskEmcalJet
{
public:
  enum TTypeImpPar
  {
    kXY,
    kXYSig,
    kXYZ,
    kXYZSig,
    kXYZSigmaOnly
  };

  AliJHFTagging();
  AliJHFTagging(const char *name);
  AliJHFTagging(const AliJHFTagging &ap);
  AliJHFTagging &operator=(const AliJHFTagging &ap);
  virtual ~AliJHFTagging();

  // methods to fill from AliAnalysisTaskSE.
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() { Init(); }
  //   virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *opt = "");
  virtual Bool_t Run();

  TClonesArray *GetInputList() const { return fInputList; }
  TClonesArray *GetInputListb() const { return fInputListb; }
  TClonesArray *GetInputListc() const { return fInputListc; }
  TClonesArray *GetInputListlf() const { return fInputListlf; }

  double GetDeltaPtRandomCone();
  double GetDeltaPtPerpEmbedding(double signalEta, double signalPhi);
  void SetEmbeddPerpendicular(bool EmbeddPerpendicular = true) { fEmbeddPerpendicular = EmbeddPerpendicular; }; // EMB_clus

  void UsePartonDefinition(bool val = true) { fUsePartonDef = val; }

  void FillCandidateJet(int CutIndex, int JetFlavor);
  void FillCandidateHFJet(AliEmcalJet *jet, TClonesArray *inputList);
  void FillCandidateHFJetMC(AliEmcalJet *jet, short jetFlavor);

  void SetDoSVAnalysis(bool value) { fDoSVAnalysis = value; }
  void SetDoTCAnalysis(bool value) { fDoTrackCountingAnalysis = value; }

  void SetTaggerWorkingPoint(double value) { fThresholdIP = value; }
  void SetMaxDispersion(double value) { fMaxDespersion = value; }
  void SetThresholdLxy(double value) { fThresholdLxy = value; }

  // B jet tracks selection
  void SetTrackMinPt(double val) { fTCMinTrackPt = val; }
  void SetTPCClusterMin(int val) { fTCMinClusTPC = val; }
  void SetITSHitsMin(int val) { fTCMinHitsITS = val; }
  void SetMaxChi2perNDF(double val) { fTCMaxChi2pNDF = val; }
  void SetMaxIPxy(double val) { fTCMaxIPxy = val; }
  void SetMaxIPz(double val) { fTCMaxIPz = val; }
  void SetMaxbDecayLength(double val) { fTCMaxDecayLength = val; }
  void SetMaxDCATrackJet(double val) { fTCMaxDCATrackJet = val; }

  void SetFillControlHistos(bool fillHist) { fFillControlHists = fillHist; }

private:
  bool CalculateJetSignedTrackImpactParameter(AliAODTrack *track, AliEmcalJet *jet, double *impar, double *cov, double &sign, double &dcajetrack, double &lineardecaylength);
  double GetValImpactParameter(TTypeImpPar type, double *impar, double *cov);

  bool IsTrackAcceptedBJetCuts(AliAODTrack *track, int jetFlavour);

  bool MatchJetsGeometricDefault();                                                        // jet matching function 1/4
  void DoJetLoop();                                                                        // jet matching function 2/4
  void SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, int matching = 0);           // jet matching function 3/4
  void GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, double &d) const; // jet matching function 4/4

  void MakeControlHistograms();

  TClonesArray *fInputList;   //! HF Jet List
  TClonesArray *fInputListb;  //! b-jet list
  TClonesArray *fInputListc;  //! c-jet list
  TClonesArray *fInputListlf; //! lf-jet list

  std::unique_ptr<AliHFJetsTagging> fHFJetUtils; //!

  bool fFillControlHists;    //
  double fThresholdIP;       //
  double fMaxDespersion;     //
  double fThresholdLxy;      //
  bool fEmbeddPerpendicular; // EMB_clus use perpendicular track embedding
  unsigned int fEvtNum;      //

  AliESDVertex *fDiamond;       //!
  AliVertexerTracks *fVertexer; //!

  TRandom3 *fRandom;             //! Random cone input
  AliFJWrapper *fFastJetWrapper; ///< EMB_clus wrapper for fast jet finding
  TRandom *fTrackGenerator;      ///< EMB_clus generator for track perpendicular to signal jet
  TClonesArray *fMCArray;        //!

  // Bjet cuts
  double fTCMinTrackPt;     //
  int fTCMinClusTPC;        //
  int fTCMinHitsITS;        //
  double fTCMaxChi2pNDF;    //
  double fTCMaxIPxy;        //
  double fTCMaxIPz;         //
  double fTCMaxDecayLength; //
  double fTCMaxDCATrackJet; //

  TH1D *fhistInclusiveJetCuts; //!
  TH1D *fhistbJetCuts;         //!
  TH1D *fhistcJetCuts;         //!
  TH1D *fhistlfJetCuts;        //!

  TH1D *fh1dJetRecPtAcceptedunCorr; //!

  TH2D *f2histRhoVsDeltaPt;       //!
  TH2D *f2histRhoVsDeltaPtFirst;  //!
  TH2D *f2histRhoVsDeltaPtSecond; //!
  TH2D *f2histRhoVsDeltaPtThird;  //!

  TH2D *fh2DeltaPtEmbeddCorrelationPerpendicular; //!

  TH1D *fh1dJetGenPt;             //! Generator level jet momentum for unfolding
  TH1D *fh1dJetGenPtUnidentified; //!
  TH1D *fh1dJetGenPtudsg;         //!
  TH1D *fh1dJetGenPtc;            //!
  TH1D *fh1dJetGenPtb;            //!

  TH1D *fh1dJetRecPt;             //! Detector level jets
  TH1D *fh1dJetRecPtAccepted;     //! Detector level jets accepted
  TH2D *fh1dJetRecEtaPhiAccepted; //! Detector level jets accepted

  TH1D *fh1dJetRecPtUnidentified;         //!
  TH1D *fh1dJetRecPtudsg;                 //!
  TH1D *fh1dJetRecPtc;                    //!
  TH1D *fh1dJetRecPtb;                    //!
  TH1D *fh1dJetRecPtUnidentifiedAccepted; //!
  TH1D *fh1dJetRecPtudsgAccepted;         //!
  TH1D *fh1dJetRecPtcAccepted;            //!
  TH1D *fh1dJetRecPtbAccepted;            //!

  TH2D *fh2dJetGenPtVsJetRecPt;     //! raw momentum response matrix
  TH2D *fh2dJetGenPtVsJetRecPtb;    //! b momentum response matrix
  TH2D *fh2dJetGenPtVsJetRecPtc;    //! c momentum response matrix
  TH2D *fh2dJetGenPtVsJetRecPtudsg; //! udsg momentum response matrix

  TH2D *fh2dJetSignedImpParXY;             //!
  TH2D *fh2dJetSignedImpParXYUnidentified; //!
  TH2D *fh2dJetSignedImpParXYudsg;         //!
  TH2D *fh2dJetSignedImpParXYb;            //!
  TH2D *fh2dJetSignedImpParXYc;            //!

  TH2D *fh2dJetSignedImpParXYSignificance;             //!
  TH2D *fh2dJetSignedImpParXYSignificanceUnidentified; //!
  TH2D *fh2dJetSignedImpParXYSignificanceudsg;         //!
  TH2D *fh2dJetSignedImpParXYSignificanceb;            //!
  TH2D *fh2dJetSignedImpParXYSignificancec;            //!

  TH2D *fh2dJetSignedImpParXYZ;             //!
  TH2D *fh2dJetSignedImpParXYZUnidentified; //!
  TH2D *fh2dJetSignedImpParXYZudsg;         //!
  TH2D *fh2dJetSignedImpParXYZb;            //!
  TH2D *fh2dJetSignedImpParXYZc;            //!

  TH2D *fh2dJetSignedImpParXYZSignificance;             //!
  TH2D *fh2dJetSignedImpParXYZSignificanceUnidentified; //!
  TH2D *fh2dJetSignedImpParXYZSignificanceudsg;         //!
  TH2D *fh2dJetSignedImpParXYZSignificanceb;            //!
  TH2D *fh2dJetSignedImpParXYZSignificancec;            //!

  TH2D *fh2dJetSignedImpParXYSignificanceFirst;             //!
  TH2D *fh2dJetSignedImpParXYSignificanceUnidentifiedFirst; //!
  TH2D *fh2dJetSignedImpParXYSignificanceudsgFirst;         //!
  TH2D *fh2dJetSignedImpParXYSignificancebFirst;            //!
  TH2D *fh2dJetSignedImpParXYSignificancecFirst;            //!

  // Second
  TH2D *fh2dJetSignedImpParXYSignificanceSecond;             //!
  TH2D *fh2dJetSignedImpParXYSignificanceUnidentifiedSecond; //!
  TH2D *fh2dJetSignedImpParXYSignificanceudsgSecond;         //!
  TH2D *fh2dJetSignedImpParXYSignificancebSecond;            //!
  TH2D *fh2dJetSignedImpParXYSignificancecSecond;            //!

  // Third
  TH2D *fh2dJetSignedImpParXYSignificanceThird;             //!
  TH2D *fh2dJetSignedImpParXYSignificanceUnidentifiedThird; //!
  TH2D *fh2dJetSignedImpParXYSignificanceudsgThird;         //!
  TH2D *fh2dJetSignedImpParXYSignificancebThird;            //!
  TH2D *fh2dJetSignedImpParXYSignificancecThird;            //!

  AliJetContainer *fJetContainerMC;   //! Container with reconstructed jets
  AliJetContainer *fJetContainerData; //! Container with reconstructed jets
  AliAODEvent *fAODIn;                //! AOD Input Event
  AliAODVertex *fPrimaryVertex;       //! Event Primary Vertex

  // Secondary Vertex
  bool fDoSVAnalysis;            //
  bool fDoTrackCountingAnalysis; //
  bool fUsePartonDef;

  AliHFJetsTaggingVertex *fVtxTagger3Prong; //!
  AliHFJetsTaggingVertex *fVtxTagger2Prong; //!

  std::unique_ptr<AliRDHFJetsCutsVertex> fjetCuts3Prong; //! SV cuts
  std::unique_ptr<AliRDHFJetsCutsVertex> fjetCuts2Prong; //! SV cuts

  AliESDtrackCuts *fEsdTrackCuts; //! cuts used on the track selected for the SV reconstruction

  TH3D *fHistSV2Prong;             //!
  TH3D *fHistSV2ProngUnidentified; //!
  TH3D *fHistSV2Prongb;            //!
  TH3D *fHistSV2Prongc;            //!
  TH3D *fHistSV2Pronglf;           //!

  TH2D *fHistDispersion2Prong;             //!
  TH2D *fHistDispersion2ProngUnidentified; //!
  TH2D *fHistDispersion2Prongb;            //!
  TH2D *fHistDispersion2Prongc;            //!
  TH2D *fHistDispersion2Pronglf;           //!

  TH3D *fHistSV3Prong;             //!
  TH3D *fHistSV3ProngUnidentified; //!
  TH3D *fHistSV3Prongb;            //!
  TH3D *fHistSV3Prongc;            //!
  TH3D *fHistSV3Pronglf;           //!

  TH2D *fHistDispersion3Prong;             //!
  TH2D *fHistDispersion3ProngUnidentified; //!
  TH2D *fHistDispersion3Prongb;            //!
  TH2D *fHistDispersion3Prongc;            //!
  TH2D *fHistDispersion3Pronglf;           //!

  ClassDef(AliJHFTagging, 1);
};
#endif // AliJHFTagging_H
