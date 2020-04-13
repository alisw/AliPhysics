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
  AliAnalysisTaskIPCalib(const char* name);
  virtual ~AliAnalysisTaskIPCalib();

  void UserCreateOutputObjects();

  //   virtual Bool_t              UserNotify()                                      ;

  void Terminate(Option_t* option);

  void SetReadMC(Bool_t readMC) {fReadMC = readMC;  }

  static AliAnalysisTaskIPCalib* AddTaskIPCalib(const char* ntracks = "usedefault",
                                                const char* suffix = "");

 protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();
  //   Bool_t                      RetrieveEventObjects()                            ;

  void AllocateTrackHistograms();

  void DoTrackLoop();

  THistManager fHistManager; ///< Histogram manager
  Bool_t fReadMC; // Flag whether to analyze MC or Data

  Float_t fAvgTrials; //! Average number of trials

 private:
  AliAnalysisTaskIPCalib(const AliAnalysisTaskIPCalib&);            // not implemented
  AliAnalysisTaskIPCalib& operator=(const AliAnalysisTaskIPCalib&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskIPCalib, 1);
  /// \endcond
};
#endif
