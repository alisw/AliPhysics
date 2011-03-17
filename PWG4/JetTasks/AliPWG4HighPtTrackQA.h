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
// This class stores QA variables as function of pT for different type
// of tracks and track selection criteria
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HIGHPTTRACKQA_H
#define ALIPWG4HIGHPTTRACKQA_H

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TList;
class TArrayF;

class AliVEvent;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDVertex;
class AliAODTrack;

class AliGenPythiaEventHeader;
class AliMCEvent;
//class AliAnalysisHelperJetTasks;

class AliPWG4HighPtTrackQA: public AliAnalysisTaskSE {

 public:
  AliPWG4HighPtTrackQA();
  AliPWG4HighPtTrackQA(const char *name);
  virtual ~AliPWG4HighPtTrackQA() {;}
 
  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify(); //Copied from AliAnalysisTaskJetSpectrum2

  enum DataType {kESD,kAOD};

  Bool_t IsPbPb() {return fIsPbPb;}  //is PbPb data?
  Bool_t SelectEvent();              //decides if event is used for analysis
  Int_t CalculateCentrality(AliVEvent *ev);
  Int_t CalculateCentrality(AliESDEvent *esd);
  Int_t CalculateCentrality(AliAODEvent *aod);
  void DoAnalysisESD();
  void DoAnalysisAOD();
  void FillHistograms();


  //Setters
  void SetDataType(DataType d)             {fDataType = d;}
  void SetIsPbPb(Bool_t cs)                {fIsPbPb = cs;}
  void SetCentralityClass(int cent)        {fCentClass=cent;}
  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}
  void SetTrackType(Int_t trackType) {fTrackType = trackType;}
  void SetFilterMask(UInt_t filterMask)    {fFilterMask = filterMask;}

  void SetSigmaConstrainedMax(Double_t sigma) {fSigmaConstrainedMax=sigma;}
  void SetPtMax(Float_t ptmax) {fPtMax = ptmax;}
  void SetNVariables(Int_t nv) {fNVariables = nv;}

  Float_t GetPtMax()           {return fPtMax;}
  Float_t GetTPCClusterInfo(AliAODTrack *tr,Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159) const;

  static AliGenPythiaEventHeader*  GetPythiaEventHeader(AliMCEvent *mcEvent);
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root

 protected:

 private:
  AliPWG4HighPtTrackQA(const AliPWG4HighPtTrackQA&);
  AliPWG4HighPtTrackQA& operator=(const AliPWG4HighPtTrackQA&);

  DataType fDataType;             //! kESD or kAOD

  AliVEvent   *fEvent;
  AliESDEvent *fESD;      //! ESD object
  const AliESDVertex   *fVtx;     //! vertex object

  AliESDtrackCuts *fTrackCuts;    // TrackCuts
  Int_t   fTrackType;             // 0: global track; 1:TPConly track 2: TPConly constrained track 3: global ITSrefit
  UInt_t fFilterMask;             //! Select tracks from specific track cuts belonging to certain filter mask for AOD analysis

  Double_t fSigmaConstrainedMax;  // max sigma on constrained fit
  Float_t fPtMax;                 // Maximum pT for histograms

  Bool_t   fIsPbPb;               //  kTRUE if PbPb
  Int_t fCentClass;               // Select only events from predefined centrality class

  /*
  0: pt
  1: phi
  2: eta
  3: dca2D
  4: dcaZ 
  5: nClustersTPC
  6: nPointITS   
  7: chi2C       
  8: nSigmaToVertex
  9: relUncertainty1Pt
  10: chi2PerClusterTPC
  11: #crossed rows
  12: (#crossed rows)/(#findable clusters)
  */
  Int_t fNVariables;             // Number of variables
  TArrayF *fVariables;           // QA variables

  Float_t fAvgTrials;             // Average number of trials
  
  TH1F *fNEventAll;                            //! Event counter
  TH1F *fNEventSel;                            //! Event counter
  TH1F *fNEventReject;                         //! Book keeping of reason of rejecting events
 
  TH1F *fh1Centrality;                         //! Centrality

  TProfile*     fh1Xsec;                       //! pythia cross section and trials
  TH1F*         fh1Trials;                     //! trials which are added
  TH1F*         fh1PtHard;                     //! pt hard of the event
  TH1F*         fh1PtHardTrials;               //! pt hard of the event

  TH1F *fh1NTracksAll;                         //! All tracks
  TH1F *fh1NTracksReject;                      //! Reason why track was rejected
  TH1F *fh1NTracksSel;                         //! Number of accepted tracks

  TH1F *fPtAll;                                //! Pt spectrum all charged particles
  TH1F *fPtSel;                                //! Pt spectrum all selected charged particles by fTrackCuts
  TH2F *fPtPhi;                                //! Pt vs Phi
  TH2F *fPtEta;                                //! Pt vs Eta
  TH2F *fPtDCA2D;                              //! Pt vs DCA2D
  TH2F *fPtDCAZ;                               //! Pt vs DCAZ
  TH2F *fPtNClustersTPC;                       //! Pt vs nClustersTPC
  TH2F *fPtNPointITS;                          //! Pt vs nPointITS
  TH2F *fPtChi2C;                              //! Pt vs Chi2C
  TH2F *fPtNSigmaToVertex;                     //! Pt vs nSigmaToVertex
  TH2F *fPtRelUncertainty1Pt;                  //! Pt vs relUncertainty1Pt
  TH2F *fPtUncertainty1Pt;                   //! Pt vs Uncertainty1Pt
  TH2F *fPtChi2PerClusterTPC;                  //! Pt vs Chi2PerClusterTPC
  TH2F *fPtNCrossedRows;                       //! Pt vs NCrossedRows
  TH2F *fPtNCrossedRowsNClusF;                 //! Pt vs NCrossedRows/NClusF
  TH3F *fPtNCrRNCrRNClusF;                     //! Pt vs NCrossedRows vs NCrossedRows/NClusF 

  //histos for covariance matrix elements
  TH2F *fPtSigmaY2;                            //! Pt vs sigma(y)^2 extCov[0]
  TH2F *fPtSigmaZ2;                            //! Pt vs sigma(z)^2 extCov[2]
  TH2F *fPtSigmaSnp2;                          //! Pt vs sigma(Snp)^2 extCov[5]
  TH2F *fPtSigmaTgl2;                          //! Pt vs sigma(Tgl)^2 extCov[9]
  TH2F *fPtSigma1Pt2;                          //! Pt vs sigma(1/pT)^2 extCov[14]

  TList *fHistList; //! List of Histograms
  
 
  ClassDef(AliPWG4HighPtTrackQA,1) 
  
};
#endif
