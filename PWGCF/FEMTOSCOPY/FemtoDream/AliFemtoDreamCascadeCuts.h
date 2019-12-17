/*
 * AliFemtoDreamCascadeCuts.h
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMCASCADECUTS_H_
#define ALIFEMTODREAMCASCADECUTS_H_
#include "Rtypes.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeHist.h"
#include "AliFemtoDreamv0MCHist.h"
#include "AliFemtoDreamTrackCuts.h"
class AliFemtoDreamCascadeCuts {
 public:
  AliFemtoDreamCascadeCuts();
  AliFemtoDreamCascadeCuts(const AliFemtoDreamCascadeCuts& cuts);
  AliFemtoDreamCascadeCuts &operator=(const AliFemtoDreamCascadeCuts& cuts);
  virtual ~AliFemtoDreamCascadeCuts();
  static AliFemtoDreamCascadeCuts *XiCuts(bool isMC, bool contribSplitting);
  static AliFemtoDreamCascadeCuts *XiFor1530Cuts(bool isMC, bool contribSplitting);
  static AliFemtoDreamCascadeCuts *OmegaCuts(bool isMC, bool contribSplitting);
  void SetMinimalBooking(bool doIt) {
    fMinimalBooking = doIt;
  }
  ;
  bool GetMinimalBooking() {
    return fMinimalBooking;
  }
  ;
  void Setv0Negcuts(AliFemtoDreamTrackCuts *cuts) {
    fNegCuts = cuts;
  }
  ;
  void Setv0PosCuts(AliFemtoDreamTrackCuts *cuts) {
    fPosCuts = cuts;
  }
  ;
  void SetBachCuts(AliFemtoDreamTrackCuts *cuts) {
    fBachCuts = cuts;
  }
  ;
  void SetIsMonteCarlo(bool isMC) {
    fMCData = isMC;
  }
  ;
  bool GetIsMonteCarlo() {
    return fMCData;
  }
  ;
  void SetContributionSplitting(bool contrSplit) {
    fContribSplitting = contrSplit;
  }
  ;
  void SetRunNumberQA(int iMinRun, int iMaxRun) {
    fRunNumberQA = true;
    fMinRunNumber = iMinRun;
    fMaxRunNumber = iMaxRun;
  }
  void SetXiMassRange(float mass, float width, float widthExcl = -1.) {
    fcutXiMass = true;
    fXiMass = mass;
    fXiMassWidth = width;
    fXiMassWidthExcl = widthExcl;
  }
  ;
  void SetXiCharge(int charge) {
    fcutXiCharge = true;
    fXiCharge = charge;
  }
  void SetCutXiDaughterDCA(float maxDCA) {
    fcutDCAXiDaug = true;
    fMaxDCAXiDaug = maxDCA;
  }
  ;
  void SetCutXiMinDistBachToPrimVtx(float minDist) {
    fcutMinDistVtxBach = true;
    fMinDistVtxBach = minDist;
  }
  ;
  void SetCutXiCPA(float cpa) {
    fcutCPAXi = true;
    fCPAXi = cpa;
  }
  ;
  void SetCutXiTransverseRadius(float minRad, float maxRad) {
    fcutXiTransRadius = true;
    fMinXiTransRadius = minRad;
    fMaxXiTransRadius = maxRad;
  }
  void Setv0MassRange(float mass, float width) {
    fcutv0Mass = true;
    fv0Mass = mass;
    fv0Width = width;
  }
  void SetCutv0MaxDaughterDCA(float maxDCA) {
    fcutv0MaxDCADaug = true;
    fv0MaxDCADaug = maxDCA;
  }
  void SetCutv0CPA(float CPA) {
    fcutCPAv0 = true;
    fCPAv0 = CPA;
  }
  void SetCutv0TransverseRadius(float minRad, float maxRad) {
    fcutv0TransRadius = true;
    fMinv0TransRadius = minRad;
    fMaxv0TransRadius = maxRad;
  }
  void SetCutv0MinDistToPrimVtx(float minDist) {
    fcutv0MinDistVtx = true;
    fv0MinDistVtx = minDist;
  }
  void SetCutv0MinDaugDistToPrimVtx(float minDist) {
    fcutv0DaugMinDistVtx = true;
    fv0DaugMinDistVtx = minDist;
  }
  ;
  void SetRejectOmegas(float mass, float width) {
    fRejOmega = true;
    fRejOmegaMass = mass;
    fRejOmegaWidth = width;
  }
  void SetRejectXis(float mass, float width) {
    fRejXi = true;
    fRejXiMass = mass;
    fRejXiWidth = width;
  }
  void SetPtRangeXi(float PtMin, float PtMax) {
    fPtMin = PtMin;
    fPtMax = PtMax;
    fCutPt = true;
  }
  void SetPtRangev0(float PtMin, float PtMax) {
    fPtMinv0 = PtMin;
    fPtMaxv0 = PtMax;
    fCutPtv0 = true;
  }
  void SetPDGCodeCasc(int PDG) {
    fPDGCasc = PDG;
  }
  ;
  int GetPDGCodeCasc() {
    return fPDGCasc;
  }
  ;
  void SetPDGCodePosDaug(int PDG) {
    fPDGPosDaug = PDG;
  }
  ;
  int GetPDGCodePosDaug() {
    return fPDGPosDaug;
  }
  ;
  void SetPDGCodeNegDaug(int PDG) {
    fPDGNegDaug = PDG;
  }
  ;
  int GetPDGCodeNegDaug() {
    return fPDGNegDaug;
  }
  ;
  void SetPDGCodeBach(int PDG) {
    fPDGBachDaug = PDG;
  }
  ;
  int GetPDGCodeBach() {
    return fPDGBachDaug;
  }
  ;
  void SetPDGCodev0(int PDG) {
    fPDGv0 = PDG;
  }
  ;
  int GetPDGv0() {
    return fPDGv0;
  }
  ;
  TString ClassName() const {
    return "AliFemtoDreamCascadeCuts";
  }
  ;
  void Init();
  bool isSelected(AliFemtoDreamCascade *casc);
  void BookQA(AliFemtoDreamCascade *casc);
  void BookCuts();
  void BookMCQA(AliFemtoDreamCascade *casc);
  void FillMCContributions(AliFemtoDreamCascade *casc);
  void SetName(TString OutputName) {
    if (fHistList)
      fHistList->SetName(OutputName.Data());
  }
  ;
  void SetMCName(TString OutputName) {
    if (fMCHistList)
      fMCHistList->SetName(OutputName.Data());
  }
  ;
  TList *GetQAHists() {
    return fHistList;
  }
  ;
  TList *GetMCQAHists() {
    return fMCHistList;
  }
  ;
  void FillGenerated(float pT) {
    if (fMCHist)
      fMCHist->FillMCGen(pT);
  }
  ;
 private:
  AliFemtoDreamCascadeHist *fHist;            //!
  AliFemtoDreamv0MCHist *fMCHist;             //!
  AliFemtoDreamTrackCuts *fNegCuts;   //
  AliFemtoDreamTrackCuts *fPosCuts;   //
  AliFemtoDreamTrackCuts *fBachCuts;  //
  TList *fHistList;  //!
  TList *fMCHistList;                 //!
  bool fMinimalBooking;               //
  bool fMCData;                       //
  bool fContribSplitting;             //
  bool fRunNumberQA;                  //
  int fMinRunNumber;                  //
  int fMaxRunNumber;                  //
  bool fCutPt;                        //
  double fPtMin;                  //
  double fPtMax;          //
  bool fCutPtv0;          //
  double fPtMinv0;  //
  double fPtMaxv0;    //
  bool fcutXiMass;    //
  float fXiMass;    //
  float fXiMassWidth;  //
  float fXiMassWidthExcl;  //
  bool fcutXiCharge;  //
  int fXiCharge;  //
  bool fcutDCAXiDaug;  //
  float fMaxDCAXiDaug;  //
  bool fcutMinDistVtxBach;  //
  float fMinDistVtxBach;    //
  bool fcutCPAXi;     //
  float fCPAXi;   //
  bool fcutXiTransRadius;  //
  float fMinXiTransRadius;  //
  float fMaxXiTransRadius;  //
  bool fcutv0Mass;  //
  float fv0Mass;    //
  float fv0Width;   //
  bool fcutv0MaxDCADaug;  //
  float fv0MaxDCADaug;    //
  bool fcutCPAv0;     //
  float fCPAv0;       //
  bool fcutv0TransRadius;       //
  float fMinv0TransRadius;      //
  float fMaxv0TransRadius;      //
  bool fcutv0MinDistVtx;        //
  float fv0MinDistVtx;          //
  bool fcutv0DaugMinDistVtx;    //
  float fv0DaugMinDistVtx;      //
  bool fRejOmega;               //
  float fRejOmegaMass;          //
  float fRejOmegaWidth;         //
  bool fRejXi;               //
  float fRejXiMass;          //
  float fRejXiWidth;         //
  int fPDGCasc;                 //
  int fPDGv0;                   //
  int fPDGPosDaug;              //
  int fPDGNegDaug;              //
  int fPDGBachDaug;             //
ClassDef(AliFemtoDreamCascadeCuts,3)
};

#endif /* ALIFEMTODREAMCASCADECUTS_H_ */
