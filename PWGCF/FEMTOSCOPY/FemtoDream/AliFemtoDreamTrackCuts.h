/*
 * AliFemtoDreamTrackCuts.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMTRACKCUTS_H_
#define ALIFEMTODREAMTRACKCUTS_H_
#include "Rtypes.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackMCHist.h"
#include "AliFemtoDreamTrackHist.h"

class AliFemtoDreamTrackCuts {
 public:
  AliFemtoDreamTrackCuts();
  AliFemtoDreamTrackCuts(const AliFemtoDreamTrackCuts& cuts);
  AliFemtoDreamTrackCuts &operator=(const AliFemtoDreamTrackCuts& cuts);
  virtual ~AliFemtoDreamTrackCuts();
  static AliFemtoDreamTrackCuts *PrimProtonCuts(bool isMC, bool DCAPlots,
                                                bool CombSigma,
                                                bool ContribSplitting);
  static AliFemtoDreamTrackCuts *PrimKaonCuts(bool isMC, bool DCAPlots = false,
                                              bool CombSigma = false,
                                              bool ContribSplitting = false);
  static AliFemtoDreamTrackCuts *PrimDeuteronCuts(bool isMC, bool DCAPlots,
                                                bool CombSigma,
                                                bool ContribSplitting);
  static AliFemtoDreamTrackCuts *DecayProtonCuts(bool isMC, bool PileUpRej,
                                                 bool ContribSplitting);
  static AliFemtoDreamTrackCuts *DecayPionCuts(bool isMC, bool PileUpRej,
                                               bool ContribSplitting);
  static AliFemtoDreamTrackCuts *Xiv0PionCuts(bool isMC, bool PileUpRej,
                                              bool ContribSplitting);
  static AliFemtoDreamTrackCuts *Xiv0ProtonCuts(bool isMC, bool PileUpRej,
                                                bool ContribSplitting);
  static AliFemtoDreamTrackCuts *XiBachPionCuts(bool isMC, bool PileUpRej,
                                                bool ContribSplitting);
  static AliFemtoDreamTrackCuts *OmegaBachKaonCuts(bool isMC, bool PileUpRej,
                                                   bool ContribSplitting);
  //  static AliFemtoDreamTrackCuts *OmegaKaonCuts(bool isMC,
  //                                               bool ContribSplitting);

  //Setters for Plots
  void SetMinimalBooking(bool doIt) {
    fMinimalBooking = doIt;
  }
  ;
  bool GetMinimalBooking() {
    return fMinimalBooking;
  }
  ;
  void SetPlotDCADist(bool plot) {
    fDCAPlots = plot;
  }
  ;
  void SetPlotTOFMass(bool mass) {
    fTOFM = mass;
  }
  ;
  void SetOriginMultiplicityHists(bool plot) {
    fDoMultBinning = plot;
  }
  void CheckParticleMothers(bool plot) {
    fCheckMother = plot;
  }
  void SetPlotCombSigma(bool plot) {
    fCombSigma = plot;
  }
  ;
  void SetPlotContrib(bool plot) {
    fContribSplitting = plot;
  }
  ;
  void SetIsMonteCarlo(bool isMC) {
    fMCData = isMC;
  }
  ;
  void SetFillQALater(bool FillLater) {
    fFillQALater = FillLater;
  }
  ;
  bool GetIsMonteCarlo() {
    return fMCData;
  }
  ;
  //Setters for the Track Cuts
  void SetCheckFilterBit(bool check) {
    fCheckFilterBit = check;
  }
  ;
  void SetFilterBit(UInt_t FilterBit) {
    fFilterBit = FilterBit;
    fCheckFilterBit = kTRUE;
  }
  ;
  void SetCheckESDFiltering(bool check) {
    fCheckESDFiltering = check;
  }
  ;
  void SetPtRange(float pmin, float pmax) {
    fpTmin = pmin;
    fpTmax = pmax;
    fcutPt = kTRUE;
  }
  ;
  void SetPtExclusion(float exmin, float exmax) {
    fpTexmin = exmin;
    fpTexmax = exmax;
    fexclPt = kTRUE;
  }
  ;
  void SetEtaRange(float etamin, float etamax) {
    fetamin = etamin;
    fetamax = etamax;
    fcutEta = kTRUE;
  }
  ;
  float GetEtaMin() {
    return fetamin;
  }
  float GetEtaMax() {
    return fetamax;
  }
  void SetCutCharge(int charge) {
    fcutCharge = kTRUE;
    fCharge = charge;
  }
  ;
  void SetCheckPileUpITS(bool check) {
    fCheckPileUpITS = check;
  }
  ;
  void SetCheckPileUpSPDTOF(bool check) {
    fCheckPileUpSPDTOF = check;
  }
  ;
  void SetCheckPileUpTOF(bool check) {
    fCheckPileUpTOF = check;
  }
  ;
  void SetCheckPileUp(bool check) {
    fCheckPileUp = check;
  }
  ;
  void SetNClsTPC(int nCls) {
    fnTPCCls = nCls;
    fcutnTPCCls = kTRUE;
  }
  ;
  void SetMaxSharedClsTPC(int nSCls) {
    fMaxSharedClsTPC = nSCls;
    fCutSharedClsTPC = true;
  }
  ;
  void SetChi2Cut(float MinChi2, float MaxChi2) {
    fMinCutChi2 = MinChi2;
    fMaxCutChi2 = MaxChi2;
    fCutChi2 = true;
  }
  ;
  void SetDCAReCalculation(bool which) {
    fDCAProp = which;
  }
  ;
  void SetDCAVtxXY(float dcaXY) {
    fDCAToVertexXY = dcaXY;
    fCutDCAToVtxXY = kTRUE;
  }
  ;
  void SetCutDCAVtxXY(bool cutdcaXY) {
    fCutDCAToVtxXY = cutdcaXY;
  }
  ;
  void SetDCAVtxZ(float dcaZ) {
    fDCAToVertexZ = dcaZ;
    fCutDCAToVtxZ = kTRUE;
  }
  ;
  void SetCutDCAVtxZ(bool cutdcaZ) {
    fCutDCAToVtxZ = cutdcaZ;
  }
  ;
  void SetCutSharedCls(bool cutit) {
    fCutSharedCls = cutit;
  }
  ;
  void SetCheckTPCRefit(bool cutit) {
    fCheckTPCRefit = cutit;
  }
  ;
  void SetCutTPCCrossedRows(bool cutit, int CrossedRows, float ratio) {
    fCutTPCCrossedRows = cutit;
    fCrossedRows = CrossedRows;
    fRatioCrossedRows = ratio;
  }
  ;
  void SetPID(AliPID::EParticleType pid, float pTPCThresh, float sigVal = 3.,
              bool AllowITSonly = false, float sigValITS = 3., bool CombITSTPC = false) {
    fParticleID = pid;
    fPIDPTPCThreshold = pTPCThresh;
    fNSigValue = sigVal;
    fAllowITSonly = AllowITSonly;
    fNSigValueITS = sigValITS;
    fCombITSTPC = CombITSTPC;
    fCutPID = kTRUE;
  }
  ;
  void SetRejLowPtPionsTOF(bool use) {
    fRejectPions = use;
  }
  ;
  void SetCutSmallestSig(bool cutit) {
    fCutHighPtSig = cutit;
  }
  ;
  int GetPDGCode();
  //selection Methods
  bool isSelected(AliFemtoDreamTrack *Track);
  void BookQA(AliFemtoDreamTrack *Track);
  void BookMC(AliFemtoDreamTrack *Track);
  void FillGenerated(float pT) {
    if (fMCHists)
      fMCHists->FillMCGen(pT);
  }
  ;
  //  void FillSharedClusterQA(AliFemtoDreamTrack *Track);
  //Histogram things
  void Init(TString name = TString("MinimalBooking"));
  TList *GetQAHists() {
    return fHists->GetHistList();
  }
  ;
  UInt_t GetFilterBit() const {
    return fFilterBit;
  }

  TList *GetMCQAHists() {
    return fMCHists ? fMCHists->GetHistList() : nullptr;
  }
  ;
  TString ClassName() {
    return "AliFemtoDreamTrackCuts";
  }
  ;
  void SetName(TString name) {
    if (fHists)
      fHists->SetName(name.Data());
  }
  ;
  void SetMCName(TString name) {
    if (fMCHists)
      fMCHists->SetName(name.Data());
  }
  ;
 private:
  bool TrackingCuts(AliFemtoDreamTrack *Track);
  bool PIDCuts(AliFemtoDreamTrack *Track);
  bool SmallestNSig(AliFemtoDreamTrack *Track);
  bool DCACuts(AliFemtoDreamTrack *Track);
  void BookTrackCuts();
  void FillMCContributions(AliFemtoDreamTrack *Track);
  AliFemtoDreamTrackMCHist *fMCHists;  //!
  AliFemtoDreamTrackHist *fHists;     //!
  bool fMinimalBooking;               //
  bool fMCData;                       //
  bool fDCAPlots;                     //
  bool fTOFM;                         //
  bool fDoMultBinning;                //
  bool fCheckMother;                  //
  bool fCombSigma;                    //
  bool fContribSplitting;             //
  bool fFillQALater;                  //
  bool fCheckFilterBit;               // This one is used for AODs
  bool fCheckESDFiltering;  // This one checks if the filtering of ESDs to AODs with FB128 passes
  bool fCheckPileUpITS;               //
  bool fCheckPileUpSPDTOF;               //
  bool fCheckPileUpTOF;               //
  bool fCheckPileUp;                //  Should only be used for Daughters of v0s
  UInt_t fFilterBit;                  //
  float fpTmin;                      //
  float fpTmax;                      //
  float fpTexmin;                   //
  float fpTexmax;                   //
  bool fcutPt;                        //
  bool fexclPt;                        //
  float fetamin;                     //
  float fetamax;                     //
  bool fcutEta;                       //
  bool fcutCharge;                    //
  int fCharge;                        //
  int fnTPCCls;                       //
  bool fcutnTPCCls;                   //
  int fMaxSharedClsTPC;               //
  bool fCutSharedClsTPC;              //
  bool fCutChi2;                      //
  float fMinCutChi2;                  //
  float fMaxCutChi2;                  //
  bool fDCAProp;  //  kTRUE means that the DCA gets recalculated by PropagateToDCA, kFALSE just uses the info stored in the AOD
  float fDCAToVertexXY;              //
  bool fCutDCAToVtxXY;                //
  float fDCAToVertexZ;               //
  bool fCutDCAToVtxZ;                 //
  bool fCutSharedCls;                 //
  bool fCheckTPCRefit;                //
  bool fCutTPCCrossedRows;            //
  int fCrossedRows;                 //
  float fRatioCrossedRows;            //
  bool fCutPID;                       //
  bool fAllowITSonly;                       //
  bool fCutHighPtSig;  // Reject tracks which have a lower Sigma for other particles (implemented for electrons, pion, kaons and protons)
  AliPID::EParticleType fParticleID;  //
  float fNSigValue;                  // defaults to 3
  float fNSigValueITS;                  // defaults to 3
  bool  fCombITSTPC;                    // defaults is false 
  float fPIDPTPCThreshold;           // defaults to 0
  bool fRejectPions;  // Supress Pions at low pT with the TOF, if information is available
ClassDef(AliFemtoDreamTrackCuts,9)
  ;
};

#endif /* ALIFEMTODREAMTRACKCUTS_H_ */
