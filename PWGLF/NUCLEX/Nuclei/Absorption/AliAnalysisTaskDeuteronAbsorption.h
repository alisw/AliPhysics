/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDeuteronAbsorption_H
#define AliAnalysisTaskDeuteronAbsorption_H

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
#include <string>

class AliPIDResponse;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliESDtrackCuts;

class AliAnalysisTaskDeuteronAbsorption : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskDeuteronAbsorption(const char *name = "AliAnalysisTaskDeuteronAbsorption");
  virtual ~AliAnalysisTaskDeuteronAbsorption();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option) {}

  double GetMindEdx() const { return fMindEdx; }
  void SetMindEdx(double dedx = 100.) { fMindEdx = dedx; }

  double GetMinTPCsignalN() const { return fMinTPCsignalN; }
  void SetMinTPCsignalN(double signalN = 50) { fMinTPCsignalN = signalN; }

  void SetESDtrackCuts(AliESDtrackCuts * cuts) { fESDtrackCuts = cuts; }

  static const AliPID::EParticleType fgkSpecies[4];
  static const std::string fgkParticleNames[4];
  static const double fgkPhiParamPos[4][4];
  static const double fgkPhiParamNeg[4][4];

private:
  double fMindEdx; /// Cut on the minimum dE/dx in TPC
  int    fMinTPCsignalN; /// Minimum number of PID clusters in the TPC

  AliPIDResponse *fPIDResponse;   //! pid response
  AliESDtrackCuts *fESDtrackCuts; //-> input track cuts
                                  //
  TList *fOutputList;             //! output list
  //
  TH1F *fHistZv;      //! Primary vertex z distribution
  TH3F *fHist3TPCpid[4];  //! QA TPC dE/dx per species
  TH3F *fHist3TPCpidAll;  //! QA TPC dE/dx no species selection
  TH3F *fHist3TOFpid[4];  //! QA TOF beta per species
  TH3F *fHist3TOFpidAll;  //! QA TOF beta no species selection
  TH3F *fHist3TOFmass[4]; //! QA TOF mass per species
  TH3F *fHist3TOFmassAll; //! QA TOF mass no species selection

  TH1F *fHist1AcceptanceAll[2][2][2]; //! Number of tracks vs p, negative (0) and positive (1), without(0) and with (1) TRD, without (0) and with (1) TOF matching
  TH2F *fHist2Matching[4][2][2]; //! TOF mass per species vs p, negative (0) and positive (1), without(0) and with (1) TRD 
  TH2F *fHist2Phi[2][2]; //! phi vs pt, negative (0) and positive (1), without(0) and with (1) TRD
  TH2F *fHist2TPCnSigma[4][2][2]; //! Number of tracks per species vs p, negative (0) and positive (1), without(0) and with (1) TRD

  TH2F *fHist2MatchingMC[4][2][2]; //! TOF mass per species vs p, negative (0) and positive (1), without(0) and with (1) TRD 

  TF1 *fTRDboundariesPos[4]; //! Function with the phi limits of TRD boundaries as a function of pt
  TF1 *fTRDboundariesNeg[4]; //! Function with the phi limits of TRD boundaries as a function of pt

  //
  AliAnalysisTaskDeuteronAbsorption(const AliAnalysisTaskDeuteronAbsorption &);            // not implemented
  AliAnalysisTaskDeuteronAbsorption &operator=(const AliAnalysisTaskDeuteronAbsorption &); // not implemented

  ClassDef(AliAnalysisTaskDeuteronAbsorption, 1);
};

#endif
