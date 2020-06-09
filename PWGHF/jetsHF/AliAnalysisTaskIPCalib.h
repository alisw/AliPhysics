#ifndef ALIANALYSISTASKIPCALIB_H
#define ALIANALYSISTASKIPCALIB_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

class AliAnalysisTaskIPCalib : public AliAnalysisTaskEmcalJet
{
public:
  AliAnalysisTaskIPCalib();
  AliAnalysisTaskIPCalib(const char *name);
  virtual ~AliAnalysisTaskIPCalib();

  void UserCreateOutputObjects();

  //   virtual Bool_t              UserNotify()                                      ;

  void Terminate(Option_t *option);

  void SetReadMC(Bool_t readMC) { fReadMC = readMC; }

  void SetCorrectResolutionPscat(Bool_t correct) { fCorrectResPscat = correct; }
  void SetCorrectResolutionNvtxContrib(Bool_t correct) { fCorrectResNvtxContrib = correct; }

  void SetCorrectionFunctionPscat(TF1 *corrFunc, Int_t nITS) { fCorrectionFactorsPscat[nITS] = corrFunc; }

  void SetCorrectionFunctionNvtxContrib(TF1 *corrFunc, Int_t nITS) { fCorrectionFactorsNvtxContrib[nITS] = corrFunc; }

  Double_t CorrectPscatIPs(Double_t sIP, Double_t pScat, Int_t nITS);
  Double_t CorrectNvtxContribIPs(Double_t sIP, Int_t nVtxContrib, Int_t nITS);

  static AliAnalysisTaskIPCalib *AddTaskIPCalib(const char *ntracks = "usedefault",
                                                  TString pathToCorrFuncPscat = "",
                                                  TString pathToCorrFuncNvtxContrib = "",
                                                  const char *suffix = "");

protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();
  //   Bool_t                      RetrieveEventObjects()                            ;

  void AllocateTrackHistograms();

  void DoTrackLoop();

  THistManager fHistManager;     ///< Histogram manager
  Bool_t fReadMC;                // Flag whether to analyze MC or Data
  Bool_t fCorrectResPscat;       // Flag whether to correct the IP resolution interms of the tracks' pScat
  Bool_t fCorrectResNvtxContrib; // Flag whether to correct the IP resolution interms of the primary vertex contributors

  TF1 *fCorrectionFactorsPscat[5];       //
  TF1 *fCorrectionFactorsNvtxContrib[5]; //

  Float_t fAvgTrials; //! Average number of trials

private:
  AliAnalysisTaskIPCalib(const AliAnalysisTaskIPCalib &);            // not implemented
  AliAnalysisTaskIPCalib &operator=(const AliAnalysisTaskIPCalib &); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskIPCalib, 4);
  /// \endcond
};
#endif
