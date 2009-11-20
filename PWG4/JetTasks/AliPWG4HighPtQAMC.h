/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// This class compares the global reconstruction with the MC information
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HIGHPTQAMC_H
#define ALIPWG4HIGHPTQAMC_H

#include "AliAnalysisTask.h"

class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDEvent;
class AliESDtrackCuts;

class AliPWG4HighPtQAMC: public AliAnalysisTask {

 public:
  AliPWG4HighPtQAMC();
  AliPWG4HighPtQAMC(const char *name);
  ~AliPWG4HighPtQAMC() {;}

  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}
  void SetCutsITS(AliESDtrackCuts* trackCutsITS) {fTrackCutsITS = trackCutsITS;}

 protected:

 private:

  void InitHistPointers();
  AliPWG4HighPtQAMC(const AliPWG4HighPtQAMC&);
  AliPWG4HighPtQAMC& operator=(const AliPWG4HighPtQAMC&);

  AliESDEvent *fESD;              //! ESD object
  AliESDtrackCuts *fTrackCuts;    // TrackCuts for global reconstructed vs MC comparison
  AliESDtrackCuts *fTrackCutsITS; // TrackCuts including ITSrefit

  
  TH1F *fNEvent;                               //! Event counter
  TH1F *fPtAll;                                //! Pt spectrum all charged particles
  TH1F *fPtSel;                                //! Pt spectrum all selected charged particles by fTrackCuts
  TH2F *fPtAllminPtMCvsPtAll;                  //! Momentum resolution (global vs MC)
  TH3F *fPtAllminPtMCvsPtAllNPointTPC;         //! Momentum resolution vs NPointTPC
  TH3F *fPtAllminPtMCvsPtAllDCAR;              //! Momentum resolution vs DCAR
  TH3F *fPtAllminPtMCvsPtAllDCAZ;              //! Momentum resolution vs DCAZ
  TH3F *fPtAllminPtMCvsPtAllPhi;               //! Momentum resolution vs Phi
  TH3F *fPtAllminPtMCvsPtAllNPointITS;         //! Momentum resolution vs NPointITS
  TH3F *fPtAllminPtMCvsPtAllNSigmaToVertex;    //! Momentum resolution vs NSigmaToVertes
  TH3F *fPtAllminPtMCvsPtAllChi2C;             //! Momentum resolution vs Chi2Constrained
  TH3F *fPtAllminPtMCvsPtAllRel1PtUncertainty; //! Momentum resolution vs relUncertainty1Pt

  TH1F *fPtAllMC;     //! Pt spectrum all charged particles
  TH1F *fPtSelMC;     //! Pt spectrum all selected charged particles by fTrackCuts
  TH1F *fPtSelMCITS;  //! Pt spectrum all selected charged particles by fTrackCutsITS

  TList *fHistList; //! List of Histograms
  
  TH1F *fPtSelITS;                              //! Pt spectrum all selected charged particles by fTrackCutsITS
  TH2F *fPtITSminPtMCvsPtITS;                  //! Momentum resolution (global with ITSrefit vs MC)
  TH3F *fPtITSminPtMCvsPtITSNPointTPC;         //! Momentum resolution vs NPointTPC 
  TH3F *fPtITSminPtMCvsPtITSDCAR;              //! Momentum resolution vs DCAR
  TH3F *fPtITSminPtMCvsPtITSDCAZ;              //! Momentum resolution vs DCAZ
  TH3F *fPtITSminPtMCvsPtITSPhi;               //! Momentum resolution vs Phi
  TH3F *fPtITSminPtMCvsPtITSNPointITS;         //! Momentum resolution vs NPointITS
  TH3F *fPtITSminPtMCvsPtITSNSigmaToVertex;    //! Momentum resolution vs NSigmaToVertex
  TH3F *fPtITSminPtMCvsPtITSChi2C;             //! Momentum resolution vs Chi2Constrained
  TH3F *fPtITSminPtMCvsPtITSRel1PtUncertainty; //! Momentum resolution vs relUncertainty1Pt

  TList *fHistListITS; //! List of Histograms

 
  ClassDef(AliPWG4HighPtQAMC,1) 
  
};
#endif
