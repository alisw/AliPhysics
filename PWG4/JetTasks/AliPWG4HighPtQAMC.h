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
class TProfile;
class TList;
class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class AliStack;
class AliESDVertex;
class AliGenPythiaEventHeader;
//class AliAnalysisHelperJetTasks;

class AliPWG4HighPtQAMC: public AliAnalysisTask {

 public:
  AliPWG4HighPtQAMC();
  AliPWG4HighPtQAMC(const char *name);
  ~AliPWG4HighPtQAMC() {;}
 
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify(); //Copied from AliAnalysisTaskJetSpectrum2

  Bool_t SelectEvent();    //decides if event is used for analysis

  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}

  void SetTrackType(Int_t trackType) {fTrackType = trackType;}
  void SetSigmaConstrainedMax(Double_t sigma) {fSigmaConstrainedMax=sigma;}
  void SetPtMax(Float_t ptmax) {fPtMax = ptmax;}
  Float_t GetPtMax()           {return fPtMax;}
  void SetPtBinEdges(Int_t region, Double_t ptmax, Double_t ptBinWidth);

  static AliGenPythiaEventHeader*  GetPythiaEventHeader(AliMCEvent *mcEvent);
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root

 protected:

 private:
  AliPWG4HighPtQAMC(const AliPWG4HighPtQAMC&);
  AliPWG4HighPtQAMC& operator=(const AliPWG4HighPtQAMC&);

  AliESDEvent *fESD;              //! ESD object
  AliMCEvent  *fMC;               //! MC event object
  AliStack    *fStack;            //! stack object

  const AliESDVertex   *fVtx;     //! vertex object

  AliESDtrackCuts *fTrackCuts;    // TrackCuts for global reconstructed vs MC comparison

  Int_t   fTrackType;             // 0: global track; 1:TPConly track 2: TPConly constrained track 3: global ITSrefit

  Double_t fSigmaConstrainedMax;  // max sigma on constrained fit
  Float_t fPtMax;                 // Maximum pT for histograms
  Float_t fPtBinEdges[3][2];      // 3 regions total with different binning for pT axis of histos

  Float_t fAvgTrials;             // Average number of trials
  
  TH1F *fNEventAll;                            //! Event counter
  TH1F *fNEventSel;                            //! Event counter
  TH1F *fNEventReject;                         //! Book keeping of reason of rejecting events
 
  TProfile*     fh1Xsec;                       //! pythia cross section and trials
  TH1F*         fh1Trials;                     //! trials which are added
  TH1F*         fh1PtHard;                     //! pt hard of the event
  TH1F*         fh1PtHardTrials;               //! pt hard of the event

  TH1F *fPtAll;                                //! Pt spectrum all charged particles
  TH1F *fPtSel;                                //! Pt spectrum all selected charged particles by fTrackCuts
  TH1F *fPtSelFakes;                           //! Pt distributions for tracks with negative label (=fake tracks)
  TH1F *fNPointTPCFakes;                       //! NTPCCluster of fake tracks
  TH1F *fPtSelLargeLabel;                      //! Filled if label is larger than nMCtracks
  TH1F *fMultRec;                              //! Bookkeeping of multiple times reconstructed tracks
  TH1F *fNPointTPCMultRec;                     //! NTPCClusters of multiple reconstructed tracks
  TH2F *fDeltaPtMultRec;                       //! Delta pT versus pT of first track for multiple reconstructed tracks

  TH2F *fPtAllvsPtMC;                          //! Reconstructed momentum vs generated momentum
  TH2F *fPtAllminPtMCvsPtAll;                  //! Momentum resolution (global vs MC)
  TH3F *fPtAllvsPtMCvsMult;                    //! Reconstructed momentum vs generated momentum vs multiplicity
  TH3F *fPtAllminPtMCvsPtAllvsMult;            //! Momentum resolution (global vs MC) vs multiplicity
  TH3F *fPtAllminPtMCvsPtAllNPointTPC;         //! Momentum resolution vs NPointTPC
  TH3F *fPtAllminPtMCvsPtAllNPointTPCIter1;    //! Momentum resolution vs NPointTPC Iter1
  TH3F *fPtAllminPtMCvsPtAllChi2TPC;           //! Momentum resolution vs Chi2TPC
  TH3F *fPtAllminPtMCvsPtAllChi2TPCIter1;      //! Momentum resolution vs Chi2TPC Iter1
  TH3F *fPtAllminPtMCvsPtAllDCAR;              //! Momentum resolution vs DCAR
  TH3F *fPtAllminPtMCvsPtAllDCAZ;              //! Momentum resolution vs DCAZ
  TH3F *fPtAllminPtMCvsPtAllPhi;               //! Momentum resolution vs Phi
  TH3F *fPtAllminPtMCvsPtAllNPointITS;         //! Momentum resolution vs NPointITS
  TH3F *fPtAllminPtMCvsPtAllNSigmaToVertex;    //! Momentum resolution vs NSigmaToVertes
  TH3F *fPtAllminPtMCvsPtAllChi2C;             //! Momentum resolution vs Chi2Constrained
  TH3F *fPtAllminPtMCvsPtAllRel1PtUncertainty; //! Momentum resolution vs relUncertainty1Pt

  TH1F *fPtAllMC;     //! Pt spectrum all charged particles
  TH1F *fPtSelMC;     //! Pt spectrum all selected charged particles by fTrackCuts

  TList *fHistList; //! List of Histograms
  
  ClassDef(AliPWG4HighPtQAMC,1) 
  
};
#endif
