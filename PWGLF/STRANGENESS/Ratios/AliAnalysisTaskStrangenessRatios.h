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
  AliAnalysisTaskStrangenessRatios(bool isMC, TString taskname = "StrangenessRatios");
  static AliAnalysisTaskStrangenessRatios* AddTask(bool isMC, TString tskname, TString suffix);
  virtual ~AliAnalysisTaskStrangenessRatios();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *) {}


  AliEventCuts  fEventCut; ///<

private:
  AliAnalysisTaskStrangenessRatios (const AliAnalysisTaskStrangenessRatios &source);
  AliAnalysisTaskStrangenessRatios &operator=(const AliAnalysisTaskStrangenessRatios &source);

  TList*          fList;             //!<! List of the output histograms
  TTree*          fTreeXi;           //!<! List of the output histograms
  TTree*          fTreeOmega;         //!<! List of the output histograms

  MiniCascade* fRecCascade; 
  MiniCascadeMC fGenCascade; 
  AliPIDResponse* fPID;              //!<! ALICE PID framework
  bool fMC;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskStrangenessRatios, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskStrangenessRatios__) */
