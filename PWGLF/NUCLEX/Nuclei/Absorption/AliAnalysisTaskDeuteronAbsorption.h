/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDeuteronAbsorption_H
#define AliAnalysisTaskDeuteronAbsorption_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"
#include "AliPID.h"
#include <string>

class AliPIDResponse;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliESDtrackCuts;

#define kNabsSpecies 5

class AliAnalysisTaskDeuteronAbsorption : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskDeuteronAbsorption(const char *name = "AliAnalysisTaskDeuteronAbsorption");
  virtual ~AliAnalysisTaskDeuteronAbsorption();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option) {}

  double GetMindEdx() const { return fMindEdx; }
  void SetMindEdx(double opt) { fMindEdx = opt; }
  void SetTreeFlag(Bool_t tmode) {fTreemode = tmode;};

  double GetMinTPCsignalN() const { return fMinTPCsignalN; }
  void SetMinTPCsignalN(double signalN = 50) { fMinTPCsignalN = signalN; }

  void SetESDtrackCuts(const AliESDtrackCuts& cuts) { fESDtrackCuts = cuts; }

  static const AliPID::EParticleType fgkSpecies[kNabsSpecies];
  static const std::string fgkParticleNames[kNabsSpecies];
  static const double fgkPhiParamPos[4][4];
  static const double fgkPhiParamNeg[4][4];

  bool fUseTRDboundariesCut;
  float fNtpcSigmas;
  AliEventCuts fEventCuts;

private:
  double fMindEdx = 100.0; /// Cut on the minimum dE/dx in TPC
  int    fMinTPCsignalN = 50; /// Minimum number of PID clusters in the TPC
  Bool_t fTreemode = kFALSE;    // Flag for filling the tree mode

  AliPIDResponse *fPIDResponse;   //! pid response
  AliESDtrackCuts fESDtrackCuts;  // input track cuts
                                  //
  TList *fOutputList;             //! output list

  TTree *fTreeTrack;  //! tree for some track parameters
  
  // Variables for the tree
  //Double_t tP;
  Double_t tPt;
  Double_t tEta;
  Double_t tPhi;
  Double_t tnsigTPC;
  Double_t tnsigTOF;
  Double_t tmass2;
  Int_t tnPIDclsTPC;
  Double_t tTOFsigDx;
  Double_t tTOFsigDz;
  Int_t tTOFclsN;
  Int_t tID;
  
  //
  TH1F *fHistZv;      //! Primary vertex z distribution
  TH3F *fHist3TPCpid[kNabsSpecies];  //! QA TPC dE/dx per species
  TH3F *fHist3TPCpidAll;  //! QA TPC dE/dx no species selection
  TH3F *fHist3TOFpid[kNabsSpecies];  //! QA TOF beta per species
  TH3F *fHist3TOFpidAll;  //! QA TOF beta no species selection
  TH3F *fHist3TOFmass[kNabsSpecies]; //! QA TOF mass per species
  TH3F *fHist3TOFmassAll; //! QA TOF mass no species selection

  TH1F *fHist1AcceptanceAll[2][2][2]; //! Number of tracks vs p, negative (0) and positive (1), without(0) and with (1) TRD, without (0) and with (1) TOF matching
  TH2F *fHist2Matching[kNabsSpecies][2][2]; //! TOF mass per species vs p, negative (0) and positive (1), without(0) and with (1) TRD 
  TH2F *fHist2Phi[2][2]; //! phi vs pt, negative (0) and positive (1), without(0) and with (1) TRD
  TH2F *fHist2TPCnSigma[kNabsSpecies][2][2]; //! Number of tracks per species vs p, negative (0) and positive (1), without(0) and with (1) TRD

  TH2F *fHist2MatchingMC[kNabsSpecies][2][2]; //! TOF mass per species vs p, negative (0) and positive (1), without(0) and with (1) TRD 

  TF1 *fTRDboundariesPos[4]; //! Function with the phi limits of TRD boundaries as a function of pt
  TF1 *fTRDboundariesNeg[4]; //! Function with the phi limits of TRD boundaries as a function of pt

  //
  AliAnalysisTaskDeuteronAbsorption(const AliAnalysisTaskDeuteronAbsorption &);            // not implemented
  AliAnalysisTaskDeuteronAbsorption &operator=(const AliAnalysisTaskDeuteronAbsorption &); // not implemented

  ClassDef(AliAnalysisTaskDeuteronAbsorption, 1);
};

#endif
