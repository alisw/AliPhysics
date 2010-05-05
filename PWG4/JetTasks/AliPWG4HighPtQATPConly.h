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
// This class compares the global reconstruction with the TPConly reconstruction
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HIGHPTQATPCONLY_H
#define ALIPWG4HIGHPTQATPCONLY_H

#include "AliAnalysisTask.h"

class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDEvent;
class AliESDfriend;
class AliESDfriendTrack;
class AliMCEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliESDtrack;

class AliPWG4HighPtQATPConly: public AliAnalysisTask {

 public:
  AliPWG4HighPtQATPConly();
  AliPWG4HighPtQATPConly(const char *name);
  ~AliPWG4HighPtQATPConly() {;}

  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  Bool_t IsCosmic(const AliESDtrack* track1 = 0x0, Int_t trackNumber = 0, Double_t ptMin = 6.);
  virtual void   Terminate(Option_t *);

  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}
  void SetCutsITS(AliESDtrackCuts* trackCutsITS) {fTrackCutsITS = trackCutsITS;}
  void SetMaxCosmicAngle(Double_t angle) {fMaxCosmicAngle = angle;}

 protected:

 private:

  void InitHistPointers();
  AliPWG4HighPtQATPConly(const AliPWG4HighPtQATPConly&);
  AliPWG4HighPtQATPConly& operator=(const AliPWG4HighPtQATPConly&);

  AliESDEvent *fESD;              //! ESD object
  AliESDfriend *fESDfriend;       //! ESD friend object
  AliMCEvent *fMC;                //! MC event object
  AliESDtrackCuts *fTrackCuts;    // TrackCuts for global vs TPConly comparison
  AliESDtrackCuts *fTrackCutsITS; // TrackCuts including ITSrefit
  
  Double_t fMaxCosmicAngle;       // Max deviation from pi (angle between two tracks) in case of cosmic candidate

  TH1F *fNEventAll;                             //! Event counter
  TH1F *fNEventSel;                             //! Event counter: Selected events for analysis
  TH1F *fPtAll;                                 //! Pt spectrum all charged particles
  TH1F *fPtSel;                                 //! Pt spectrum all selected charged particles by fTrackCuts
  TH2F *fPtAllminPtTPCvsPtAll;                  //! Momentum resolution (global vs TPConly)
  TH3F *fPtAllminPtTPCvsPtAllEtaPos;            //! Momentum resolution (global vs TPConly) vs Eta for positive particles
  TH3F *fPtAllminPtTPCvsPtAllEtaNeg;            //! Momentum resolution (global vs TPConly) vs Eta for negative particles
  TH3F *fPtAllminPtTPCvsPtAllNPointTPC;         //! Momentum resolution vs NPointTPC
  TH3F *fPtAllminPtTPCvsPtAllNPointTPCS;        //! Momentum resolution vs NPointTPCShared/NPointTPC
  TH3F *fPtAllminPtTPCvsPtAllDCAR;              //! Momentum resolution vs DCAR
  TH3F *fPtAllminPtTPCvsPtAllDCAZ;              //! Momentum resolution vs DCAZ
  TH3F *fPtAllminPtTPCvsPtAllPhi;               //! Momentum resolution vs Phi
  TH3F *fPtAllminPtTPCvsPtAllNPointITS;         //! Momentum resolution vs NPointITS
  TH3F *fPtAllminPtTPCvsPtAllNSigmaToVertex;    //! Momentum resolution vs NSigmaToVertes
  TH3F *fPtAllminPtTPCvsPtAllChi2C;             //! Momentum resolution vs Chi2Constrained
  TH3F *fPtAllminPtTPCvsPtAllRel1PtUncertainty; //! Momentum resolution vs relUncertainty1Pt
  TH3F *fPtAllminPtTPCvsPtAllChi2PerNClusTPC;   //! Momentum resolution vs Chi2PerNClusTPC
  TH3F *fPtAllminPtTPCvsPtAllChi2PerNClusITS;   //! Momentum resolution vs Chi2PerNClusITS

  TH3F *fPtAllminPtTPCvsNPointTPCPhi;           //! Momentum resolution vs NPointTPC vs Phi
  TH3F *fPtAllminPtTPCvsNPointITSPhi;           //! Momentum resolution vs NPointITS vs Phi
  TH3F *fPtAllminPtTPCvsRel1PtUncertaintyPhi;   //! Momentum resolution vs Rel1PtUncertainty vs Phi

  TH2F *fEtaPhiOutliers;                        //! Eta Phi for outliers in momentum resolution
 
  TH1F *fPtSelITSouter;                         //! Pt at ITS outer wall for all selected charged particles by fTrackCuts
  TH2F *fPtITSouterminPtTPCvsPtAll;                  //! Momentum resolution (global vs ITSouter-TPCinner)
  TH3F *fPtITSouterminPtTPCvsPtAllEtaPos;            //! Momentum resolution (global vs ITSouter-TPCinner) vs Eta for positive particles
  TH3F *fPtITSouterminPtTPCvsPtAllEtaNeg;            //! Momentum resolution (global vs ITSouter-TPCinner) vs Eta for negative particles
  TH3F *fPtITSouterminPtTPCvsPtAllNPointTPC;         //! Momentum resolution vs NPointTPC
  TH3F *fPtITSouterminPtTPCvsPtAllNPointTPCS;        //! Momentum resolution vs NPointTPCS
  TH3F *fPtITSouterminPtTPCvsPtAllDCAR;              //! Momentum resolution vs DCAR
  TH3F *fPtITSouterminPtTPCvsPtAllDCAZ;              //! Momentum resolution vs DCAZ
  TH3F *fPtITSouterminPtTPCvsPtAllPhi;               //! Momentum resolution vs Phi
  TH3F *fPtITSouterminPtTPCvsPtAllNPointITS;         //! Momentum resolution vs NPointITS
  TH3F *fPtITSouterminPtTPCvsPtAllNSigmaToVertex;    //! Momentum resolution vs NSigmaToVertes
  TH3F *fPtITSouterminPtTPCvsPtAllChi2C;             //! Momentum resolution vs Chi2Constrained
  TH3F *fPtITSouterminPtTPCvsPtAllRel1PtUncertainty; //! Momentum resolution vs relUncertainty1Pt
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC;   //! Momentum resolution vs Chi2PerNClusTPC
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITS;   //! Momentum resolution vs Chi2PerNClusITS
 
  TH2F *fPtITSouterminPtTPCvsPtAllITSLayer0;                  //! Track has at least 1st SPD layer
  TH2F *fPtITSouterminPtTPCvsPtAllITSLayer1;                  //! Track has at least 2nd SPD layer and not 1st SPD
  TH2F *fPtITSouterminPtTPCvsPtAllITSLayer2;                  //! Track has at least 1st SDD layer and not SPD layers
  TH2F *fPtITSouterminPtTPCvsPtAllITSLayer3;                  //! Track has at least 1st SDD layer and not SPD layers and not 1st SDD
  TH2F *fPtITSouterminPtTPCvsPtAllITSLayer4;                  //! Track has at least 1st SSD layer and not SPD or SDD layers
  TH2F *fPtITSouterminPtTPCvsPtAllITSLayer5;                  //! Track has at least 1st SSD layer and not SPD or SDD layers or 1st SSD

  TH2F *fPtITSouterminPtTPCvsPtAllNoSPD;                  //! Track has no signal in SPD layers
  TH2F *fPtITSouterminPtTPCvsPtAllNoSDD;                  //! Track has no signal in SDD layers
  TH2F *fPtITSouterminPtTPCvsPtAllNoSSD;                  //! Track has no signal in SSD layers

  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0;                  //! Track has at least 1st SPD layer
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1;                  //! Track has at least 2nd SPD layer and not 1st SPD
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2;                  //! Track has at least 1st SDD layer and not SPD layers
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3;                  //! Track has at least 1st SDD layer and not SPD layers and not 1st SDD
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4;                  //! Track has at least 1st SSD layer and not SPD or SDD layers
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5;                  //! Track has at least 1st SSD layer and not SPD or SDD layers or 1st SSD

  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD;                  //! Track has no signal in SPD layers
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD;                  //! Track has no signal in SDD layers
  TH3F *fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD;                  //! Track has no signal in SSD layers

  TList *fHistList; //! List of Histograms
  
  TH1F *fPtAllTPC;     //! Pt spectrum all charged particles
  TH1F *fPtSelTPC;     //! Pt spectrum all selected charged particles by fTrackCuts
  TH1F *fPtSelTPCITS;  //! Pt spectrum all selected charged particles by fTrackCutsITS
  TList *fHistListTPC; //! List of Histograms
  
  TH1F *fPtSelITS;                              //! Pt spectrum all selected charged particles by fTrackCutsITS
  TH2F *fPtITSminPtTPCvsPtITS;                  //! Momentum resolution (global with ITSrefit vs TPConly)
  TH3F *fPtITSminPtTPCvsPtITSEtaPos;            //! Momentum resolution (global with ITSrefit vs TPConly) vs Eta for positive particles
  TH3F *fPtITSminPtTPCvsPtITSEtaNeg;            //! Momentum resolution (global with ITSrefit vs TPConly) vs Eta for negative particles
  TH3F *fPtITSminPtTPCvsPtITSNPointTPC;         //! Momentum resolution vs NPointTPC 
  TH3F *fPtITSminPtTPCvsPtITSNPointTPCS;        //! Momentum resolution vs NPointTPCS 
  TH3F *fPtITSminPtTPCvsPtITSDCAR;              //! Momentum resolution vs DCAR
  TH3F *fPtITSminPtTPCvsPtITSDCAZ;              //! Momentum resolution vs DCAZ
  TH3F *fPtITSminPtTPCvsPtITSPhi;               //! Momentum resolution vs Phi
  TH3F *fPtITSminPtTPCvsPtITSNPointITS;         //! Momentum resolution vs NPointITS
  TH3F *fPtITSminPtTPCvsPtITSNSigmaToVertex;    //! Momentum resolution vs NSigmaToVertex
  TH3F *fPtITSminPtTPCvsPtITSChi2C;             //! Momentum resolution vs Chi2Constrained
  TH3F *fPtITSminPtTPCvsPtITSRel1PtUncertainty; //! Momentum resolution vs relUncertainty1Pt
  TH3F *fPtITSminPtTPCvsPtITSChi2PerNClusTPC;   //! Momentum resolution vs Chi2PerNClusTPC
  TH3F *fPtITSminPtTPCvsPtITSChi2PerNClusITS;   //! Momentum resolution vs Chi2PerNClusITS

  TH3F *fPtITSminPtTPCvsNPointTPCPhi;           //! Momentum resolution vs NPointTPC vs Phi
  TH3F *fPtITSminPtTPCvsNPointITSPhi;           //! Momentum resolution vs NPointITS vs Phi
  TH3F *fPtITSminPtTPCvsRel1PtUncertaintyPhi;   //! Momentum resolution vs Rel1PtUncertainty vs Phi

  TH3F *fPtRel1PtUncertaintyChi2PerClusTPC;     //! Global Pt vs relUncertainty1Pt vs Chi2PerClusTPC
  TH3F *fPtNPointTPCSChi2PerClusTPC;            //! Global Pt vs NPointTPCShared/NPointTPC vs Chi2PerClusTPC
  TH3F *fPtNPointTPCSRel1PtUncertainty;         //! Global Pt vs NPointTPCShared/NPointTPC vs relUncertainty1Pt
  TH3F *fPtRel1PtUncertaintyChi2PerClusTPCSharedSel;     //! Global Pt vs relUncertainty1Pt vs Chi2PerClusTPC for NPointTPCShared/NPointTPC>0.05

  TList *fHistListITS; //! List of Histograms

  TH1F *fPtCosmicCandidates;                    //! Cosmic Candidates
  TH1F *fDeltaPtCosmicCandidates;               //! Cosmic Candidates Delta Pt
  TH1F *fDeltaPhi;                              //! Cosmic Candidates Delta Phi
  TH1F *fDeltaEta;                              //! Cosmic Candidates Delta Eta

  TList *fHistListCosmics;                      //! List of Histograms for cosmic candidates

  ClassDef(AliPWG4HighPtQATPConly,1) 
  
};
#endif
