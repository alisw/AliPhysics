#ifndef AliAnalysisTaskIPCalib_H
#define AliAnalysisTaskIPCalib_H

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

  void SetCorrectResolution(Bool_t correct) { fCorrectRes = correct; }

  void SetCorrectionFunction(TF1 *corrFunc, Int_t nITS) { fCorrectionFactors[nITS] = corrFunc; }

  Double_t CorrectIPs(Double_t sIP, Double_t pScat, Int_t nITS);

  static AliAnalysisTaskIPCalib *AddTaskIPCalib(const char *ntracks = "usedefault",
                                                  TString pathToCorrFunc = "",
                                                  const char *suffix = "");

protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();
  //   Bool_t                      RetrieveEventObjects()                            ;

  void AllocateTrackHistograms();

  void DoTrackLoop();

  THistManager fHistManager; ///< Histogram manager
  Bool_t fReadMC;            // Flag whether to analyze MC or Data
  Bool_t fCorrectRes;        // Flag whether to correct the IP resolution

  TF1 *fCorrectionFactors[5]; //

  Float_t fAvgTrials; //! Average number of trials

private:
  AliAnalysisTaskIPCalib(const AliAnalysisTaskIPCalib &);            // not implemented
  AliAnalysisTaskIPCalib &operator=(const AliAnalysisTaskIPCalib &); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskIPCalib, 3);
  /// \endcond
};
#endif
