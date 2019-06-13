/*
 * AliFemtoDreamv0Cuts.h
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMV0CUTS_H_
#define ALIFEMTODREAMV0CUTS_H_
#include "Rtypes.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0MCHist.h"
#include "AliFemtoDreamv0Hist.h"
#include "AliFemtoDreamv0.h"
class AliFemtoDreamv0Cuts {
 public:
  enum PileUpRejectionMode {
    BothDaughtersCombined = 0,  ///< Both daughters, ITS SPD/SSD or TOF
    OneDaughterCombined = 1,    ///< One daughters, ITS SPD/SSD or TOF
    BothDaughtersITSonly = 2,   ///< Both daughters, ITS SPD/SSD only
    BothDaughtersTOFonly = 3,   ///< Both daughters, TOF only
    OneDaughterITSonly = 4,     ///< One daughter, ITS SPD/SSD only
    OneDaughterTOFonly = 5,     ///< One daughter, TOF only
    BothDaughersSPDonly = 6,    ///< Both daughters, SPD only
    OneDaugherSPDonly = 7,      ///< One daughter, SPD only
    BothDaughersSPDTOF = 8,     ///< Both daughters, SPD or TOF
    OneDaugherSPDTOF = 9        ///< One daughter, SPD or TOF
  };

  AliFemtoDreamv0Cuts();
  AliFemtoDreamv0Cuts(const AliFemtoDreamv0Cuts& cuts);
  AliFemtoDreamv0Cuts& operator=(const AliFemtoDreamv0Cuts& cuts);
  virtual ~AliFemtoDreamv0Cuts();
  static AliFemtoDreamv0Cuts* LambdaCuts(bool isMC, bool CPAPlots,
                                         bool SplitContrib);
  static AliFemtoDreamv0Cuts *LambdaSigma0Cuts(bool isMC, bool CPAPlots,
                                               bool SplitContrib);
  //Setters for plots
  void SetMinimalBooking(bool doIt) {
    fMinimalBooking = doIt;
  }
  ;
  bool GetMinimalBooking() {
    return fMinimalBooking;
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
  void SetPlotCPADist(bool plot) {
    fCPAPlots = plot;
  }
  ;
  void SetPlotContrib(bool plot) {
    fContribSplitting = plot;
  }
  ;
  void SetOriginMultiplicityHists(bool plot) {
    fDoMultBinning = plot;
  }
  void CheckParticleMothers(bool plot) {
    fCheckMother = plot;
  }
  void SetAxisInvMassPlots(int nBins, float minMass, float maxMass) {
    fNumberXBins = nBins;
    fAxisMinMass = minMass;
    fAxisMaxMass = maxMass;
  }
  void SetRunNumberQA(int iMinRun, int iMaxRun) {
    fRunNumberQA = true;
    fMinRunNumber = iMinRun;
    fMaxRunNumber = iMaxRun;
  }
  //Setters for the daughter track cuts
  void SetPosDaugterTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fPosCuts = cuts;
    fCutDaughters = true;
  }
  ;
  void SetNegDaugterTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fNegCuts = cuts;
    fCutDaughters = true;
  }
  ;
  void SetCutDaughters(bool cut) {
    fCutDaughters = cut;
  }
  bool GetCutDaughters() const { return fCutDaughters;};
  //Setters for PDG Codes of the daughters+v0
  void SetPDGCodev0(int pdgCode) {
    fPDGv0 = pdgCode;
  }
  ;
  int GetPDGv0() const {
    return fPDGv0;
  }
  ;
  void SetPDGCodePosDaug(int pdgCode) {
    fPDGDaugP = pdgCode;
  }
  ;
  int GetPDGPosDaug() const {
    return fPDGDaugP;
  }
  ;
  void SetPDGCodeNegDaug(int pdgCode) {
    fPDGDaugN = pdgCode;
  }
  ;
  int GetPDGNegDaug() const {
    return fPDGDaugN;
  }
  ;
  //Setters v0 cuts
  void SetCheckOnFlyStatus(bool val) {
    fOnFlyStatus = val;
    fCutOnFlyStatus = true;
  }
  ;
  void SetCutCharge(int charge) {
    fCutCharge = true;
    fCharge = charge;
  }
  ;
  void SetPtRange(float pmin, float pmax) {
    fpTmin = pmin;
    fpTmax = pmax;
    fCutPt = true;
  }
  void SetKaonRejection(float MassLow, float MassUp) {
    fKaonRejection = true;
    fKaonRejLow = MassLow;
    fKaonRejUp = MassUp;
  }
  ;
  void SetCutMaxDecayVtx(float maxDecayVtx) {
    fCutDecayVtxXYZ = true;
    fMaxDecayVtxXYZ = maxDecayVtx;
  }
  ;
  void SetCutTransverseRadius(float minRadius, float maxRadius) {
    fMinTransRadius = minRadius;
    fMaxTransRadius = maxRadius;
    fCutTransRadius = true;
  }
  void SetCutDCADaugToPrimVtx(float minDCA) {
    fMinDCADaugToPrimVtx = minDCA;
    fCutMinDCADaugPrimVtx = true;
  }
  ;
  void SetCutDCADaugTov0Vtx(float maxDCA) {
    fMaxDCADaugToDecayVtx = maxDCA;
    fCutMaxDCADaugToDecayVtx = true;
  }
  void SetCutCPA(float cpa) {
    fMinCPA = cpa;
    fCutCPA = true;
  }
  ;
  void SetCutInvMass(float width) {
    fInvMassCutWidth = width;
    fCutInvMass = true;
  }
  ;
  void SetCutWindow(float down, float up) {
    fCutInvMass = false;
    fCutInvMassSidebands = true;
    fInvMassCutSBdown = down;
    fInvMassCutSBup = up;
  }
  void Init();
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
  void FillGenerated(float pT) {
    if (fMCHist)
      fMCHist->FillMCGen(pT);
  }
  ;
  TList *GetQAHists() {
    return fHistList;
  }
  ;
//  TList *GetQAHistsPosDaug() {return fPosCuts->GetQAHists();};
//  TList *GetQAHistsNegDaug() {return fNegCuts->GetQAHists();};
  TList *GetMCQAHists() {
    return fMCHistList;
  }
  ;
//  TList *GetMCQAHistsPosDaug() {return fPosCuts->GetMCQAHists();};
//  TList *GetMCQAHistsNegDaug() {return fNegCuts->GetMCQAHists();};
  bool isSelected(AliFemtoDreamv0 *v0);
  TString ClassName() {
    return "v0Cuts";
  }
  ;
  void SetArmenterosCut(bool doIt) {
    fDoArmenterosCut = doIt;
  }
  void SetArmenterosCut(float qtLow, float qtUp, float alphaLow,
                        float alphaUp) {
    fDoArmenterosCut = true;
    fArmenterosQtLow = qtLow;
    fArmenterosQtUp = qtUp;
    fArmenterosAlphaLow = alphaLow;
    fArmenterosAlphaUp = alphaUp;
  }
  void SetDaughterTimingCut(PileUpRejectionMode mode) {
    fDoCombinedTimingCut = true;
    fCombinedTiming = mode;
  }

 private:
  bool TimingCut(AliFemtoDreamv0 *v0);
  bool RejectAsKaon(AliFemtoDreamv0 *v0);
  bool ArmenterosSelection(AliFemtoDreamv0 *v0);
  bool DaughtersPassCuts(AliFemtoDreamv0 *v0);
  bool MotherPassCuts(AliFemtoDreamv0 *v0);
  bool CPAandMassCuts(AliFemtoDreamv0 *v0);
  void BookQA(AliFemtoDreamv0 *v0);
  void BookMC(AliFemtoDreamv0 *v0);
  void BookTrackCuts();
  void FillMCContributions(AliFemtoDreamv0 *v0);
  float CalculateInvMass(AliFemtoDreamv0 *v0, int PDGPosDaug, int PDGNegDaug);
  TList *fHistList;                   //!
  TList *fMCHistList;                 //!
  AliFemtoDreamv0MCHist *fMCHist;     //!
  AliFemtoDreamv0Hist *fHist;         //!
  //Here the cut values for the Daughters are stored
  AliFemtoDreamTrackCuts *fPosCuts;   //
  AliFemtoDreamTrackCuts *fNegCuts;   //
  bool fCutDaughters;
  //These are all the cuts directly linked to the v0
  bool fMinimalBooking;               //
  bool fMCData;                       //
  bool fCPAPlots;                     //
  bool fContribSplitting;             //
  bool fDoMultBinning;               //
  bool fCheckMother;                 //
  bool fRunNumberQA;                  //
  int fMinRunNumber;                  //
  int fMaxRunNumber;                  //
  bool fCutOnFlyStatus;               //
  bool fOnFlyStatus;                  //
  bool fCutCharge;                    //
  int fCharge;                        //
  bool fDoCombinedTimingCut;          //
  PileUpRejectionMode fCombinedTiming;//
  bool fDoArmenterosCut;              //
  float fArmenterosQtLow;             //
  float fArmenterosQtUp;              //
  float fArmenterosAlphaLow;          //
  float fArmenterosAlphaUp;           //
  bool fCutPt;                        //
  float fpTmin;                      //
  float fpTmax;                      //
  bool fKaonRejection;                //
  float fKaonRejLow;                 //
  float fKaonRejUp;                  //
  bool fCutDecayVtxXYZ;               //
  float fMaxDecayVtxXYZ;             //
  bool fCutTransRadius;               //
  float fMinTransRadius;             //
  float fMaxTransRadius;             //
  bool fCutMinDCADaugPrimVtx;         //
  float fMinDCADaugToPrimVtx;        //
  bool fCutMaxDCADaugToDecayVtx;      //
  float fMaxDCADaugToDecayVtx;       //
  bool fCutCPA;                       //
  float fMinCPA;                     //
  bool fCutInvMass;                   //
  float fInvMassCutWidth;            //
  bool fCutInvMassSidebands;         //
  float fInvMassCutSBdown;           //
  float fInvMassCutSBup;             //

  //Range for the axis of the hists
  float fAxisMinMass;                //
  float fAxisMaxMass;                //
  int fNumberXBins;                   //
  //PDG Codes of the Mother and the Daughter needed for Inv Mass Calc. and
  //matching in the MC Sample
  int fPDGv0;                         //
  int fPDGDaugP;                      //
  int fPDGDaugN;                      //
ClassDef(AliFemtoDreamv0Cuts,6)
};

#endif /* ALIFEMTODREAMV0CUTS_H_ */
