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
  virtual ~AliFemtoDreamTrackCuts();
  static AliFemtoDreamTrackCuts *PrimProtonCuts(bool isMC,bool DCAPlots,
                                                bool CombSigma,
                                                bool ContribSplitting);
  static AliFemtoDreamTrackCuts *DecayProtonCuts(bool isMC,
                                                 bool ContribSplitting);
  static AliFemtoDreamTrackCuts *DecayPionCuts(bool isMC,
                                               bool ContribSplitting);
  static AliFemtoDreamTrackCuts *Xiv0PionCuts(bool isMC,
                                               bool ContribSplitting);
  static AliFemtoDreamTrackCuts *Xiv0ProtonCuts(bool isMC,
                                               bool ContribSplitting);
  static AliFemtoDreamTrackCuts *XiBachPionCuts(bool isMC,
                                               bool ContribSplitting);
//  static AliFemtoDreamTrackCuts *OmegaKaonCuts(bool isMC,
//                                               bool ContribSplitting);

  //Setters for Plots
  void SetPlotDCADist(bool plot) {fDCAPlots=plot;};
  void SetPlotCombSigma(bool plot) {fCombSigma=plot;};
  void SetPlotContrib(bool plot) {fContribSplitting=plot;};
  void SetIsMonteCarlo(bool isMC) {fMCData=isMC;};
  void SetFillQALater(bool FillLater){fFillQALater=FillLater;};
  bool GetIsMonteCarlo() {return fMCData;};
  //Setters for the Track Cuts
  void SetCheckFilterBit(bool check){fCheckFilterBit = check;};
  void SetFilterBit(UInt_t FilterBit){fFilterBit = FilterBit; fCheckFilterBit = kTRUE;};
  void SetPtRange(double pmin, double pmax){fpTmin = pmin; fpTmax = pmax; fcutPt = kTRUE;};
  void SetEtaRange(double etamin, double etamax){fetamin=etamin; fetamax=etamax; fcutEta = kTRUE;};
  double GetEtaMin() {return fetamin;}
  double GetEtaMax() {return fetamax;}
  void SetCutCharge(int charge){fcutCharge = kTRUE; fCharge = charge;};
  void SetCheckPileUpITS(bool check){fCheckPileUpITS=check;};
  void SetCheckPileUpTOF(bool check){fCheckPileUpTOF=check;};
  void SetCheckPileUp(bool check){fCheckPileUp=check;};
  void SetNClsTPC(int nCls){fnTPCCls = nCls; fcutnTPCCls = kTRUE;};
  void SetDCAReCalculation(bool which){fDCAProp = which;};
  void SetDCAVtxXY(double dcaXY){fDCAToVertexXY = dcaXY; fCutDCAToVtxXY = kTRUE;};
  void SetCutDCAVtxXY(bool cutdcaXY){fCutDCAToVtxXY = cutdcaXY;};
  void SetDCAVtxZ(double dcaZ){fDCAToVertexZ = dcaZ; fCutDCAToVtxZ = kTRUE;};
  void SetCutDCAVtxZ(bool cutdcaZ){fCutDCAToVtxZ = cutdcaZ;};
  void SetCutSharedCls(bool cutit){fCutSharedCls = cutit;};
  void SetCheckTPCRefit(bool cutit){fCheckTPCRefit=cutit;};
  void SetCutTPCCrossedRows(bool cutit){fCutTPCCrossedRows = cutit;};
  void SetPID(AliPID::EParticleType pid, double pTPChresh, double sigVal = 3.)
  {fParticleID = pid; fPIDPTPCThreshold = pTPChresh; fNSigValue = sigVal; fCutPID = kTRUE;};
  void SetRejLowPtPionsTOF(bool use){fRejectPions = use;};
  void SetCutSmallestSig(bool cutit){fCutHighPtSig = cutit;};
  //selection Methods
  bool isSelected(AliFemtoDreamTrack *Track);
  void BookQA(AliFemtoDreamTrack *Track);
  void BookMC(AliFemtoDreamTrack *Track);
//  void FillSharedClusterQA(AliFemtoDreamTrack *Track);
  //Histogram things
  void Init();
  TList *GetQAHists() {return fHists->GetHistList();};
  TList *GetMCQAHists() {return fMCHists->GetHistList();};
  TString ClassName() {return "AliFemtoDreamTrackCuts";};
  void SetName(TString name){fHists->SetName(name.Data());};
 private:
  bool TrackingCuts(AliFemtoDreamTrack *Track);
  bool PIDAODCuts(AliFemtoDreamTrack *Track);
  bool SmallestNSig(AliFemtoDreamTrack *Track);
  bool DCACuts(AliFemtoDreamTrack *Track);
  void BookTrackCuts();
  void FillMCContributions(AliFemtoDreamTrack *Track);
  AliFemtoDreamTrackMCHist *fMCHists; //!
  AliFemtoDreamTrackHist *fHists;     //!
  bool fMCData;                       //
  bool fDCAPlots;                     //
  bool fCombSigma;                    //
  bool fContribSplitting;             //
  bool fFillQALater;                  //
  bool fCheckFilterBit;               //
  bool fCheckPileUpITS;               //
  bool fCheckPileUpTOF;               //
  bool fCheckPileUp;                  //  Should only be used for Daughters of v0s
  UInt_t fFilterBit;                  //
  double fpTmin;                      //
  double fpTmax;                      //
  bool fcutPt;                        //
  double fetamin;                     //
  double fetamax;                     //
  bool fcutEta;                       //
  bool fcutCharge;                    //
  int fCharge;                        //
  int fnTPCCls;                       //
  bool fcutnTPCCls;                   //
  bool fDCAProp;                      //  kTRUE means that the DCA gets recalculated by PropagateToDCA, kFALSE just uses the info stored in the AOD
  double fDCAToVertexXY;              //
  bool fCutDCAToVtxXY;                //
  double fDCAToVertexZ;               //
  bool fCutDCAToVtxZ;                 //
  bool fCutSharedCls;                 //
  bool fCheckTPCRefit;                //
  bool fCutTPCCrossedRows;            //
  bool fCutPID;                       //
  bool fCutHighPtSig;                 // Reject tracks which have a lower Sigma for other particles (implemented for electrons, pion, kaons and protons)
  AliPID::EParticleType fParticleID;  //
  double fNSigValue;                  // defaults to 3
  double fPIDPTPCThreshold;           // defaults to 0
  bool fRejectPions;                  // Supress Pions at low pT with the TOF, if information is available
  ClassDef(AliFemtoDreamTrackCuts,1);
};

#endif /* ALIFEMTODREAMTRACKCUTS_H_ */
