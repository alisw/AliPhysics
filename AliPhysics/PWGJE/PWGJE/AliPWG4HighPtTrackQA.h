/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
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
class AliAODVertex;
class AliAODTrack;
class AliESDtrack;
class AliAnalysisUtils;

class AliGenPythiaEventHeader;
class AliMCEvent;
//class AliAnalysisHelperJetTasks;

/**
 * @class AliPWG4HighPtTrackQAs
 * @brief Track-based QA for the jet working group
 * @author Marta Verweij - UU
 * 
 * This class stores QA variables as function of pT for different type
 * of tracks and track selection criteria
 * Output: Histograms for different set of cuts
 */
class AliPWG4HighPtTrackQA: public AliAnalysisTaskSE {

 public:
  AliPWG4HighPtTrackQA();
  AliPWG4HighPtTrackQA(const char *name);
  virtual ~AliPWG4HighPtTrackQA() {;}
 
  enum DataType {kESD,kAOD};

  Bool_t IsPbPb() {return fIsPbPb;}  //is PbPb data?
  Bool_t SelectEvent();              //decides if event is used for analysis
  Int_t CalculateCentrality(AliVEvent *ev);
  Int_t CalculateCentrality(AliESDEvent *esd);
  Int_t CalculateCentrality(const AliAODEvent *aod);
  Int_t GetCentralityClass(Float_t cent=-1.) const;
  void DoAnalysisESD();
  void DoAnalysisAOD();
  void FillHistograms();

  //Setters
  void SetDataType(DataType d)                     {fDataType=d;}
  void SetIsPbPb(Bool_t cs)                        {fIsPbPb=cs;}
  void SetEMCALTrigger(const char *emcaltrigger)   {fEmcalTrigger=emcaltrigger;}
  void SetCentralityClass(int cent)                {fCentClass=cent;}
  void SetCuts(AliESDtrackCuts* trackCuts)         {fTrackCuts         = trackCuts;}
  void SetCutsITSLoose(AliESDtrackCuts* trackCuts) {fTrackCutsITSLoose = trackCuts;}
  void SetCutsTPConly(AliESDtrackCuts* trackCuts)  {fTrackCutsTPConly  = trackCuts;}
  void SetTrackType(Int_t trackType)               {fTrackType=trackType;}
  void SetFilterMask(UInt_t filterMask)            {fFilterMask=filterMask;}
  void SetIncludeNoITS(Bool_t f)                   {fIncludeNoITS=f;}
  void SetUseSPDTrackletVsClusterBG(Bool_t b)      {fTklVsClusSPDCut=b;}
  void SetZvertexDiffValue(Double_t cut)           {fZvertexDiff=cut;}

  void SetSigmaConstrainedMax(Double_t sigma)      {fSigmaConstrainedMax=sigma;}
  void SetPtMax(Float_t ptmax)                     {fPtMax = ptmax;}
  void SetPtBinEdges(Int_t region, Double_t ptmax, Double_t ptBinWidth);
  void SetNVariables(Int_t nv)                     {fNVariables = nv;}

  Float_t GetPtMax()           {return fPtMax;}
  Float_t GetTPCClusterInfo(const AliAODTrack *tr,Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159, Bool_t useFitMap=kFALSE) const;
  Float_t GetTPCClusterInfoFitMap(const AliESDtrack *tr,Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159) const;
  Int_t   GetTrackLengthTPC(const AliESDtrack *track) const;
  Int_t   GetTrackLengthTPC(const AliAODTrack *track) const;
  Float_t GetGoldenChi2(AliESDtrack *origtrack);
  Float_t GetGGCChi2(AliESDtrack *origtrack);

  static AliGenPythiaEventHeader*  GetPythiaEventHeader(const AliMCEvent *mcEvent);
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root

 protected:
  virtual void   Init();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify(); //Copied from AliAnalysisTaskJetSpectrum2

 private:
  /**
   * @enum Variable_t  
   * @brief Index of observables in the array
   */
  enum Variable_t {
    kVarPT = 0,                   ///< pt
    kVarPhi = 1,                  ///< phi
    kVarEta = 2,                  ///< eta
    kVarDCAXY = 3,                ///< dcaXY
    kVarDCAZ = 4,                 ///< dcaZ
    kVarClustersTPC = 5,          ///< nClustersTPC
    kVarPointsITS = 6,            ///< nPointsITS
    kVarChi2C = 7,                ///< chi2C
    kVarNsigVtx = 8,              ///< nSigmaToVertex
    kVarRelUncertaintyInvPt = 9,  ///< relUncertainty1Pt
    kVarChi2ClusterTPC = 10,      ///< chi2PerClusterTPC
    kVarCrossedRowsTPC = 11,      ///< #crossed rows
    kVarCrossedRowsFindTPC = 12,  ///< (#crossed rows)/(#findable clusters)
    kVarSigmaY2 = 13,             ///< SigmaY2
    kVarSigmaZ2 = 14,             ///< SigmaZ2
    kVarSigmaSnp2 = 15,           ///< SigmaSnp2
    kVarSigmaTgl2 = 16,           ///< SigmaTgl2
    kVarSigmaInvPt2 = 17,         ///< Sigma1Pt2
    kVarClustersTPCIter1 = 18,    ///< NClustersTPCIter1
    kVarChi2TPCIter1 = 19,        ///< TPCChi2Iter1
    kVarClustersTPCShared = 20,   ///< NClustersTPCShared
    kVarChi2Gold = 21,            ///< Chi2Gold (TPC constrained vs global)
    kVarChi2ConstrainedGlob = 22, ///< Chi2GGC (global constrained vs global)
    kVarNCrossedRowsFit = 23,     ///< NCrossed rows from fit map
    kVarCrossedRowsFindFit = 24,  ///< (#crossed rows)/(#findable clusters) from fit map
    kVarChi2ITS = 25,             ///< ITS chi2
    kVarTPCclustersdEdx = 26      ///< TPC number of clusters for dE/dx
  };

  AliPWG4HighPtTrackQA(const AliPWG4HighPtTrackQA&);
  AliPWG4HighPtTrackQA& operator=(const AliPWG4HighPtTrackQA&);

  DataType                        fDataType;                              ///< kESD or kAOD

  AliVEvent                       *fEvent;                                //!<! AliVEvent object
  AliESDEvent                     *fESD;                                  //!<! ESD object
  const AliESDVertex              *fVtx;                                  //!<! vertex object
  const AliAODVertex              *fVtxAOD;                               //!<! vertex object AOD

  AliESDtrackCuts                 *fTrackCuts;                            ///< TrackCuts
  AliESDtrackCuts                 *fTrackCutsITSLoose;                    ///< Loose ITS track cuts
  AliESDtrackCuts                 *fTrackCutsTPConly;                     ///< TPC track cuts
  Int_t                           fTrackType;                             ///< 0: global track; 1:TPConly track 2: TPConly constrained track 3: global ITSrefit 4: TPConly constrained track with QA selection based on global track
  UInt_t                          fFilterMask;                            ///< Select tracks from specific track cuts belonging to certain filter mask for AOD analysis
  Bool_t                          fIncludeNoITS;                          ///< includes tracks with failed ITS refit

  Double_t                        fSigmaConstrainedMax;                   ///< max sigma on constrained fit
  Float_t                         fPtMax;                                 ///< Maximum pT for histograms
  Float_t                         fPtBinEdges[4][2];                      ///< 4 regions total with different binning for pT axis of histos

  Bool_t                          fIsPbPb;                                ///< kTRUE if PbPb
  Int_t                           fCentClass;                             ///< Select only events from predefined centrality class
  TString                         fEmcalTrigger;                          ///< EMCAL trigger class (data only)
  Bool_t                          fInit;                                  ///< true after initialization (relevant for ESD analysis)
  AliAnalysisUtils                *fAliAnalysisUtils;                     //!<!vertex selection (optional)
  Double_t                        fVertex[3];                             //!<!event vertex
  Double_t                        fVertexSPD[3];                          //!<!event Svertex
  Int_t                           fNVertCont;                             //!<!event vertex number of contributors
  Int_t                           fNVertSPDCont;                          //!<!event SPD vertex number of contributors

  Bool_t                          fTklVsClusSPDCut;                       ///< Apply tracklet-vs-cluster SPD cut to reject background events in pp
  Double_t                        fZvertexDiff;                           ///< upper limit for distance between primary and SPD vertex

  Int_t                           fNVariables;                            ///< Number of variables
  TArrayF                         *fVariables;                            ///< QA variables

  UChar_t                         fITSClusterMap;                         ///< map of clusters, one bit per a layer

  Float_t                         fAvgTrials;                             ///< Average number of trials
  
  TH1F                            *fNEventAll;                            //!<! Event counter
  TH1F                            *fNEventSel;                            //!<! Event counter
  TH1F                            *fNEventReject;                         //!<! Book keeping of reason of rejecting events
  TH1F                            *fhEvMult;                              //!<! histo for event multplicity
  TH1F                            *fhTrackletsMult;                       //!<! histo for tracklets multplicity
  
  TH1F                            *fh1Centrality;                         //!<! Centrality

  TProfile                        *fh1Xsec;                               //!<! pythia cross section and trials
  TH1F                            *fh1Trials;                             //!<! trials which are added
  TH1F                            *fh1PtHard;                             //!<! pt hard of the event
  TH1F                            *fh1PtHardTrials;                       //!<! pt hard of the event

  TH1F                            *fh1NTracksAll;                         //!<! All tracks
  TH1F                            *fh1NTracksReject;                      //!<! Reason why track was rejected
  TH1F                            *fh1NTracksSel;                         //!<! Number of accepted tracks

  TH1F                            *fPtAll;                                //!<! Pt spectrum all charged particles
  TH1F                            *fPtSel;                                //!<! Pt spectrum all selected charged particles by fTrackCuts
  TH2F                            *fPtPhi;                                //!<! Pt vs Phi
  TH2F                            *fPtEta;                                //!<! Pt vs Eta
  TH3F                            *fPtEtaPhi;                             //!<! Pt vs Eta vs Phi
  TH2F                            *fPtDCA2D;                              //!<! Pt vs DCA2D
  TH2F                            *fPtDCAZ;                               //!<! Pt vs DCAZ
  TH2F                            *fPtNClustersTPC;                       //!<! Pt vs nClustersTPC
  TH2F                            *fPtNClustersTPCPhi;                    //!<! Phi vs nClustersTPC
  TH2F                            *fPtNClustersTPCIter1;                  //!<! Pt vs nClustersTPCIter1
  TH3F                            *fPtNClustersTPCIter1Phi;               //!<! Pt vs nClustersTPCIter1 vs Phi
  TH2F                            *fPtNClustersTPCShared;                 //!<! Pt vs nClustersTPCShared
  TH2F                            *fPtNClustersTPCSharedFrac;             //!<! Pt vs nClustersTPCSharedFrac
  TH2F                            *fPtNPointITS;                          //!<! Pt vs nPointITS
  TH3F                            *fPtNPointITSPhi;                       //!<! Pt vs nPointITS vs phi
  TH2F                            *fPtChi2C;                              //!<! Pt vs Chi2C
  TH2F                            *fPtNSigmaToVertex;                     //!<! Pt vs nSigmaToVertex

  TH2F                            *fPtRelUncertainty1Pt;                  //!<! Pt vs relUncertainty1Pt
  TH3F                            *fPtRelUncertainty1PtNClus;             //!<! Pt vs relUncertainty1Pt vs NClustersTPC
  TH3F                            *fPtRelUncertainty1PtNClusIter1;        //!<! Pt vs relUncertainty1Pt vs NClustersTPCIter1
  TH3F                            *fPtRelUncertainty1PtNPointITS;         //!<! Pt vs relUncertainty1Pt vs NPointITS
  TH3F                            *fPtRelUncertainty1PtITSClusterMap;     //!<! Pt vs relUncertainty1Pt vs byte of ITS clustermap

  TH3F                            *fPtRelUncertainty1PtChi2;              //!<! Pt vs relUncertainty1Pt vs Chi2TPC/NClus
  TH3F                            *fPtRelUncertainty1PtChi2Iter1;         //!<! Pt vs relUncertainty1Pt vs Chi2TPC/NClusIter1
  TH3F                            *fPtRelUncertainty1PtPhi;               //!<! Pt vs relUncertainty1PtPhi

  TH2F                            *fPtChi2PerClusterTPC;                  //!<! Pt vs Chi2PerClusterTPC
  TH2F                            *fPtChi2PerClusterTPCIter1;             //!<! Pt vs Chi2PerClusterTPCIter1
  TH2F                            *fPtNCrossedRows;                       //!<! Pt vs NCrossedRows
  TH3F                            *fPtNCrossedRowsPhi;                    //!<! Pt vs NCrossedRows vs Phi
  TH3F                            *fPtNCrossedRowsNClusFPhi;              //!<! Pt vs NCrossedRows/NClusF vs Phi
  TH3F                            *fPtNCrRNCrRNClusF;                     //!<! Pt vs NCrossedRows vs NCrossedRows/NClusF 
  TH2F                            *fPtNCrossedRowsFit;                    //!<! Pt vs NCrossedRows from NClusterFitMap
  TH3F                            *fPtNCrossedRowsFitPhi;                 //!<! Pt vs NCrossedRows from FitMap vs Phi
  TH3F                            *fPtNCrossedRowsNClusFFitPhi;           //!<! Pt vs NCrossedRows/NClusF from NClusterFitMap
  TH2F                            *fNCrossedRowsNCrossedRowsFit;          //!<! NCrossedRows from ClusterMap vs NCrossedRows from NClusterFitMap
  TH2F                            *fNClustersNCrossedRows;                //!<! NClusters vs NCrossedRows
  TH2F                            *fNClustersNCrossedRowsFit;             //!<! NClusters vs NCrossedRows from fit map
  TH3F                            *fPtNClustersNClustersFitMap;           //!<! pT vs ncls vs ncls from fit map
  TH3F                            *fPtRelUncertainty1PtNCrossedRows;      //!<! Pt vs relUncertainty1Pt vs NCrossedRows
  TH3F                            *fPtRelUncertainty1PtNCrossedRowsFit;   //!<! Pt vs relUncertainty1Pt vs NCrossedRowsFitMap
  
  TH2F                            *fPtTPCSignalN;                         //!<! Pt vs TPCsignalN
  
  TH2F                            *fPtChi2Gold;                           //!<! Pt vs Chi2 between global and TPC constrained track
  TH2F                            *fPtChi2GGC;                            //!<! Pt vs Chi2 between global and global constrained track
  TH3F                            *fPtChi2GoldPhi;                        //!<! Pt vs Chi2 between global and TPC constrained track vs phi
  TH3F                            *fPtChi2GGCPhi;                         //!<! Pt vs Chi2 between global and global constrained track vs phi
  TH2F                            *fChi2GoldChi2GGC;                      //!<! Correlations between gold chi2 and GGC chi2

  TH3F                            *fPtChi2ITSPhi;                         //!<! Pt vs ITS chi2 vs phi

  //histos for covariance matrix elements
  TH2F                            *fPtSigmaY2;                            //!<! 1/Pt vs sigma(y) extCov[0]
  TH2F                            *fPtSigmaZ2;                            //!<! 1/Pt vs sigma(z) extCov[2]
  TH2F                            *fPtSigmaSnp2;                          //!<! 1/Pt vs sigma(Snp) extCov[5]
  TH2F                            *fPtSigmaTgl2;                          //!<! 1/Pt vs sigma(Tgl) extCov[9]
  TH2F                            *fPtSigma1Pt2;                          //!<! 1/Pt vs sigma(1/pT) extCov[14]

  //profiles for covariance matrix elements
  TProfile                        *fProfPtSigmaY2;                        //!<! 1/Pt vs sigma(y) extCov[0]
  TProfile                        *fProfPtSigmaZ2;                        //!<! 1/Pt vs sigma(z) extCov[2]
  TProfile                        *fProfPtSigmaSnp2;                      //!<! 1/Pt vs sigma(Snp) extCov[5]
  TProfile                        *fProfPtSigmaTgl2;                      //!<! 1/Pt vs sigma(Tgl) extCov[9]
  TProfile                        *fProfPtSigma1Pt2;                      //!<! 1/Pt vs sigma(1/pT) extCov[14]

  TProfile                        *fProfPtSigma1Pt;                       //!<! pT vs sigma(1/Pt)
  TProfile                        *fProfPtPtSigma1Pt;                     //!<! pT vs pT*sigma(1/Pt)

  TList                           *fHistList;                             //!<! List of Histograms
 
  ClassDef(AliPWG4HighPtTrackQA,10)
};
#endif
