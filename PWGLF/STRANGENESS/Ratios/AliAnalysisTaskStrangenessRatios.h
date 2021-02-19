/// \class AliAnalysisTaskStrangenessRatios

#ifndef __AliAnalysisTaskStrangenessRatios__
#define __AliAnalysisTaskStrangenessRatios__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include "AliEventCuts.h"

class AliPIDResponse;
class TH2F;
class TList;
class TTree;

struct MiniCascade {
  Double32_t pt;
  Double32_t eta;
  Double32_t mass;
  Double32_t ct;
  Double32_t radius; //[0,25.4,8]
  Double32_t radiusV0; //[0,25.4,8]
  Double32_t dcaBachPV; //[0,2.54,8]
  Double32_t dcaV0PV; //[0,2.54,8]
  Double32_t dcaV0piPV; //[0,2.54,8]
  Double32_t dcaV0prPV; //[0,2.54,8]
  Double32_t dcaV0tracks; //[0,2.54,8]
  Double32_t dcaBachV0; //[0,2.54,8]
  Double32_t cosPA; //[0.95,1,16]
  Double32_t cosPAV0; //[0.95,1,16]
  Double32_t V0invMassDelta; //[-0.01,0.01,8]
  Double32_t tpcNsigmaBach; //[-5,5,8]
  Double32_t tpcNsigmaV0Pr; //[-5,5,8]
  Double32_t tpcNsigmaV0Pi; //[-5,5,8]
  Double32_t competingMass; //[0,0.254,8]
  Double32_t bachBarCosPA;  //[0.9999, 1., 8]
  unsigned char tpcClBach;
  unsigned char tpcClV0Pr;
  unsigned char tpcClV0Pi;
  unsigned char centrality;
  bool matter;
  bool hasTOFhit;
  bool hasITSrefit;
};

struct MiniCascadeMC : public MiniCascade {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  int pdg;
  bool isReconstructed;
};

class AliAnalysisTaskStrangenessRatios : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskStrangenessRatios(bool isMC = false, TString taskname = "StrangenessRatios");
  static AliAnalysisTaskStrangenessRatios* AddTask(bool isMC, TString tskname, TString suffix);
  virtual ~AliAnalysisTaskStrangenessRatios();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *) {}

  AliEventCuts  fEventCut; ///<

  //Setters for topological cuts
  void SetRadiusCut(float xi = 1.2, float omega = 1.) {fCutRadius[0]=xi; fCutRadius[1]=omega;}
  void SetRadiusV0Cut(float cut = 3.0) {fCutRadiusV0=cut;}
  void SetDCABachToPVCut(float cut = 0.1) {fCutDCABachToPV=cut;}
  void SetDCAV0toPVCut(float cut = 0.1) {fCutDCAV0toPV=cut;}
  void SetDCAV0piToPVCut(float cut = 0.2) {fCutDCAV0piToPV=cut;}
  void SetDCAV0prToPVCut(float cut = 0.2) {fCutDCAV0prToPV=cut;}
  void SetDCAV0tracksCut(float cut = 1.0) {fCutDCAV0tracks=cut;}
  void SetCutDCABachToV0XiCut(float xi = 1.0, float omega = 0.6) {fCutDCABachToV0[0]=xi; fCutDCABachToV0[0]=omega;}
  void SetCosPACut(float cut = 0.95) {fCutCosPA=cut;}
  void SetCosPAV0Cut(float cut = 0.95) {fCutCosPAV0=cut;}
  void SetV0MassWindowCut(float cut = 0.005) {fCutV0MassWindow=cut;}
  void SetYCut(float cut = 0.5) {fCutY=cut;}
  void SetYDaughtCut(float cut = 0.8) {fCutYDaught=cut;}
  void SetNsigmaTPCCut(float cut = 4.0) {fCutNsigmaTPC=cut;}
  void SetCtCut(float cut = 4 /*cm*/) {fCutCt=cut;}
  void SetCtV0Cut(float cut = 30 /*cm*/) {fCutCtV0=cut;}
  void SetCompetingMassCut(float cut = 0.008 /*GeV/c^2*/) {fCutCompetingMass=cut;}
  void SetTPCcluCut(int cut = 70) {fCutTPCclu=cut;}
  void SetSaveOnlyTrueCandidates(bool cut = true) { fOnlyTrueCandidates = cut; }

private:
  AliAnalysisTaskStrangenessRatios (const AliAnalysisTaskStrangenessRatios &source);
  AliAnalysisTaskStrangenessRatios &operator=(const AliAnalysisTaskStrangenessRatios &source);

  TList*          fList = nullptr;             //!<! List of the output histograms
  TTree*          fTreeXi = nullptr;           //!<! Tree for Xis
  TTree*          fTreeOmega = nullptr;        //!<! Tree for Omegas

  MiniCascade* fRecCascade = nullptr;          //!<! Transient fRecCascade
  MiniCascadeMC fGenCascade;
  AliPIDResponse* fPID = nullptr;              //!<! ALICE PID framework
  bool fMC;
  bool fOnlyTrueCandidates = false;  ///< Save only true Xi and Omegas in MC

  //configurable cuts
  float fCutRadius[2] = {1.2, 1.0};
  float fCutRadiusV0 = 3.0;
  float fCutDCABachToPV = 0.1;
  float fCutDCAV0toPV = 0.1;
  float fCutDCAV0piToPV = 0.2;
  float fCutDCAV0prToPV = 0.2;
  float fCutDCAV0tracks = 1.2;
  float fCutDCABachToV0[2]{1.0, 1.0};
  float fCutCosPA = 0.95;
  float fCutCosPAV0 = 0.95;
  float fCutV0MassWindow = 0.005;
  float fCutY = 0.5;
  float fCutYDaught = 0.8;
  float fCutNsigmaTPC = 4.0;
  float fCutCt = 4;
  float fCutCtV0 = 30;
  float fCutCompetingMass = 0.;
  int fCutTPCclu = 70;

  bool IsTopolSelected(bool isXi = true);
  float Eta2y(float pt, float m, float eta) const;


  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskStrangenessRatios, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskStrangenessRatios__) */
