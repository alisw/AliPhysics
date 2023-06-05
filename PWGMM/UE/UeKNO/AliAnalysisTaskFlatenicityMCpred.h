/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskFlatenicityMCpred_H
#define AliAnalysisTaskFlatenicityMCpred_H

class AliESDtrackCuts;
class AliESDEvent;
class AliESDAD;
class TList;
class TH1D;
class TH2D;
class TH1I;
class TProfile;

#include "AliAnalysisTaskSE.h"

#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"

class AliAnalysisTaskFlatenicityMCpred : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskFlatenicityMCpred();
  AliAnalysisTaskFlatenicityMCpred(const char *name);
  virtual ~AliAnalysisTaskFlatenicityMCpred();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  double GetFlatenicityV0();
  double GetFlatenicityMC();
  void ExtractMultiplicities();
  void ExtractMultiplicitiesMC();
  void SetNbinsNch(int nbins = 600) { fNchBin = nbins; }
  void SetMinFlat(double rhomin = 0.) { fminrho = rhomin; }
  void SetMaxFlat(double rhomax = 0.) { fmaxrho = rhomax; }

  bool HasRecVertex();
  void FillHistos1(Int_t multGen, const std::vector<double> &ptGen,
                   const std::vector<double> &yGen,
                   const std::vector<double> &etaGen,
                   const std::vector<int> &idGen);

  int FillMCarray(std::vector<double> &pt, std::vector<double> &yGen,
                  std::vector<double> &etaGen, std::vector<int> &id);

protected:
private:
  AliESDEvent *fESD; //! input ESD event
  AliEventCuts fEventCuts;
  AliStack *fMCStack; //! MC stack
  AliMCEvent *fMC;    //! MC Event
  double fVtxz;
  double fminrho;
  double fmaxrho;
  int fNchBin;
  int fV0Mindex;
  double fv0multalice;
  double flatenicity_m;
  double fmultV0A_m;
  double fmultV0C_m;
  double flatenicity_t;
  int fmultV0A_t;
  int fmultV0C_t;
  int fnchAll;

  TList *fOutputList; //! output list in the root file
  double fEtaCut;
  double fv0mpercentile;
  AliMultSelection *fMultSelection;
  TH1D *hCounter;
  TH1D *hV0MBadruns;

  TH2D *hMultVsFlat1;
  TH2D *hflatVsNchV01;

  TH3D *hMultVsFlatVsPt1[4];
  TH2D *hFlatVsPt1[4];

  TH1D *hMultV0MPerc[9];
  TH1D *hFlatV0MPerc[9];

  TH2D *hFlatVsPt1V0MPerc[4][9];

  AliAnalysisTaskFlatenicityMCpred(
      const AliAnalysisTaskFlatenicityMCpred &); // not implemented
  AliAnalysisTaskFlatenicityMCpred &
  operator=(const AliAnalysisTaskFlatenicityMCpred &); // not implemented

  ClassDef(AliAnalysisTaskFlatenicityMCpred, 3);
};

#endif
