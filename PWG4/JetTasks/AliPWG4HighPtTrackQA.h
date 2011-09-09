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
class AliESDtrack;

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
  Int_t GetCentralityClass(Float_t cent=-1.);
  void DoAnalysisESD();
  void DoAnalysisAOD();
  void FillHistograms();

  void FillSystematicCutHist(AliESDtrack *track);

  //Setters
  void SetDataType(DataType d)             {fDataType = d;}
  void SetIsPbPb(Bool_t cs)                {fIsPbPb = cs;}
  void SetCentralityClass(int cent)        {fCentClass=cent;}
  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}
  void SetCutsITSLoose(AliESDtrackCuts* trackCuts) {fTrackCutsITSLoose = trackCuts;}
  void SetCutsTPConly(AliESDtrackCuts* trackCuts) {fTrackCutsTPConly = trackCuts;}
  void SetTrackType(Int_t trackType) {fTrackType = trackType;}
  void SetFilterMask(UInt_t filterMask)    {fFilterMask = filterMask;}

  void SetSigmaConstrainedMax(Double_t sigma) {fSigmaConstrainedMax=sigma;}
  void SetPtMax(Float_t ptmax) {fPtMax = ptmax;}
  void SetPtBinEdges(Int_t region, Double_t ptmax, Double_t ptBinWidth);
  void SetNVariables(Int_t nv) {fNVariables = nv;}

  Float_t GetPtMax()           {return fPtMax;}
  Float_t GetTPCClusterInfo(AliAODTrack *tr,Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159) const;
  Int_t   GetTrackLengthTPC(AliESDtrack *track);
  Int_t   GetTrackLengthTPC(AliAODTrack *track);

  static AliGenPythiaEventHeader*  GetPythiaEventHeader(AliMCEvent *mcEvent);
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root

 protected:

 private:
  AliPWG4HighPtTrackQA(const AliPWG4HighPtTrackQA&);
  AliPWG4HighPtTrackQA& operator=(const AliPWG4HighPtTrackQA&);

  DataType fDataType;             //! kESD or kAOD

  AliVEvent   *fEvent;            //! AliVEvent object
  AliESDEvent *fESD;              //! ESD object
  const AliESDVertex   *fVtx;     //! vertex object

  AliESDtrackCuts *fTrackCuts;         // TrackCuts
  AliESDtrackCuts *fTrackCutsITSLoose; // Loose ITS track cuts
  AliESDtrackCuts *fTrackCutsTPConly;  // TPC track cuts
  Int_t   fTrackType;                  // 0: global track; 1:TPConly track 2: TPConly constrained track 3: global ITSrefit 4: TPConly constrained track with QA selection based on global track
  UInt_t fFilterMask;                  //! Select tracks from specific track cuts belonging to certain filter mask for AOD analysis

  Double_t fSigmaConstrainedMax;  // max sigma on constrained fit
  Float_t fPtMax;                 // Maximum pT for histograms
  Float_t fPtBinEdges[3][2];      // 3 regions total with different binning for pT axis of histos

  Bool_t   fIsPbPb;               // kTRUE if PbPb
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
  TH2F *fPtNClustersTPCIter1;                  //! Pt vs nClustersTPCIter1
  TH2F *fPtNClustersTPCShared;                 //! Pt vs nClustersTPCShared
  TH2F *fPtNClustersTPCSharedFrac;             //! Pt vs nClustersTPCSharedFrac
  TH2F *fPtNPointITS;                          //! Pt vs nPointITS
  TH2F *fPtChi2C;                              //! Pt vs Chi2C
  TH2F *fPtNSigmaToVertex;                     //! Pt vs nSigmaToVertex
  TH2F *fPtRelUncertainty1Pt;                  //! Pt vs relUncertainty1Pt
  TH3F *fPtRelUncertainty1PtNClus;             //! Pt vs relUncertainty1Pt vs NClustersTPC
  TH3F *fPtRelUncertainty1PtNClusIter1;        //! Pt vs relUncertainty1Pt vs NClustersTPCIter1
  TH3F *fPtRelUncertainty1PtChi2;              //! Pt vs relUncertainty1Pt vs Chi2TPC/NClus
  TH3F *fPtRelUncertainty1PtChi2Iter1;         //! Pt vs relUncertainty1Pt vs Chi2TPC/NClusIter1
  TH3F *fPtRelUncertainty1PtPhi;               //! Pt vs relUncertainty1PtPhi
  TH3F *fPtRelUncertainty1PtTrkLength;         //! Pt vs relUncertainty1Pt vs track length in TPC
  TH2F *fPtUncertainty1Pt;                     //! Pt vs Uncertainty1Pt
  TH2F *fPtChi2PerClusterTPC;                  //! Pt vs Chi2PerClusterTPC
  TH2F *fPtChi2PerClusterTPCIter1;                  //! Pt vs Chi2PerClusterTPCIter1
  TH2F *fPtNCrossedRows;                       //! Pt vs NCrossedRows
  TH2F *fPtNCrossedRowsNClusF;                 //! Pt vs NCrossedRows/NClusF
  TH3F *fPtNCrRNCrRNClusF;                     //! Pt vs NCrossedRows vs NCrossedRows/NClusF 

  //histos for covariance matrix elements
  TH2F *fPtSigmaY2;                            //! 1/Pt vs sigma(y) extCov[0]
  TH2F *fPtSigmaZ2;                            //! 1/Pt vs sigma(z) extCov[2]
  TH2F *fPtSigmaSnp2;                          //! 1/Pt vs sigma(Snp) extCov[5]
  TH2F *fPtSigmaTgl2;                          //! 1/Pt vs sigma(Tgl) extCov[9]
  TH2F *fPtSigma1Pt2;                          //! 1/Pt vs sigma(1/pT) extCov[14]

  //profiles for covariance matrix elements
  TProfile *fProfPtSigmaY2;                    //! 1/Pt vs sigma(y) extCov[0]
  TProfile *fProfPtSigmaZ2;                    //! 1/Pt vs sigma(z) extCov[2]
  TProfile *fProfPtSigmaSnp2;                  //! 1/Pt vs sigma(Snp) extCov[5]
  TProfile *fProfPtSigmaTgl2;                  //! 1/Pt vs sigma(Tgl) extCov[9]
  TProfile *fProfPtSigma1Pt2;                  //! 1/Pt vs sigma(1/pT) extCov[14]

  TProfile *fProfPtSigma1Pt;                   //! pT vs sigma(1/Pt)
  TProfile *fProfPtPtSigma1Pt;                 //! pT vs pT*sigma(1/Pt)

  TH1F *fSystTrackCuts;                        //! Bookkeeping for track cuts

  TList *fHistList; //! List of Histograms
  
 
  ClassDef(AliPWG4HighPtTrackQA,1) 
  
};
#endif
