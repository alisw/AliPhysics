#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$
#include "AliAnalysisTaskEmcal.h"
#include "AliEventPoolManager.h"
#include "AliEventCuts.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"
#include <THn.h>
#include <THnSparse.h>

class TH1;
class TH2;
class TH3;
class AliVVZERO;
class AliEvtPoolManager;

using std::vector;

class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcal {
 public:
	AliAnalysisTaskGammaHadron();
	AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputSameEventAnalysis);
virtual ~AliAnalysisTaskGammaHadron();

  //setters for the analysis
  void                        SetEffHistGamma(THnF *h)                              { fHistEffGamma    = h      ; }
  void                        SetEffHistHadron(THnF *h)                             { fHistEffHadron   = h      ; }
  void                        SetSavePool(Bool_t input)                             { fSavePool        = input  ; }
  void                        SetPlotMore(Bool_t input)                             { fPlotQA          = input  ; }
  void                        SetEvtTriggerType(UInt_t input)                       { fTriggerType     = input  ; }
  void                        SetEvtMixType(UInt_t input)                           { fMixingEventType = input  ; }
//  void                        SetCutsId(Id, cent, ptCl,Ecl,,.... UInt_t input)      { fMixingEventType = input  ; }
  void                        SetNLM(Int_t input)                                   { fMaxNLM = input;}
  void                        SetM02(Double_t inputMin,Double_t inputMax)           { fClShapeMin = inputMin; fClShapeMax = inputMax;}
  void                        SetRmvMatchedTrack(Bool_t input, Double_t dEta=-1, Double_t dPhi=-1) { fRmvMTrack  = input; fTrackMatchEta=dEta; fTrackMatchPhi=dPhi;}
  void                        SetUseManualEvtCuts(Bool_t input)                     { fUseManualEventCuts = input;}
  void                        SetDoRotBkg(Bool_t input)                             { fDoRotBkg = input;}
  void                        SetNRotBkgSamples(Int_t input)                        { fNRotBkgSamples = input;}


  //Functions for mixed event purposes
  void                        SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager*        GetEventPoolManager()                                 {return fPoolMgr;}
  // Set which pools will be saved
  void                        AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt);
 private:
  AliEventCuts                fEventCuts;                   ///< event selection utility
  AliFiducialCut*             fFiducialCuts;                ///< fiducial cuts for the EMCal and DCal in terms of eta and phi
  AliEMCALRecoUtils*          fFiducialCellCut;             ///< fiducial cut for EMCal+DCal in terms of rows and collumns

 protected:

  void                        InitArrays()                                                 ;
  // overwritten EMCal base class functions
  Bool_t                      Run()                             	                          ;
  void                        ExecOnce()         									      ;
  Bool_t                      IsEventSelected()											  ;
  void                        UserCreateOutputObjects()        		                      ;

  //Functions for mixed event purposes
  void                        InitEventMixer()											  ;
  TObjArray*                  CloneToCreateTObjArray(AliParticleContainer* tracks)          ;

  //..Correlate and fill
  Bool_t                      FillHistograms()                                              ;
  Int_t                       CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  Int_t                       CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  void                        FillPi0CandsHist(AliTLorentzVector CaloClusterVec,AliTLorentzVector CaloClusterVec2,AliTLorentzVector CaloClusterVecPi0,Double_t fMaxClusM02,Double_t Weight);
  void                        FillGhHisograms(Int_t identifier,AliTLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t Weight);
  void                        FillQAHisograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec);
  Bool_t                      AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster);
  Bool_t                      DetermineMatchedTrack(AliVCluster* caloCluster,Double_t &etadiff,Double_t & phidiff);

  //..Delta phi does also exist in AliAnalysisTaskEmcal. It is overwritten here (ask Raymond)
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec) ;
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,Double_t phi_EVP)       ;
  Double_t                    GetEff(AliTLorentzVector ParticleVec)                         ;
  void                        GetDistanceToSMBorder(AliVCluster* caloCluster,Int_t &etaCellDist,Int_t &phiCellDist);
  AliVCluster*                GetLeadingCluster(const char* opt, AliClusterContainer* clusters);

  Bool_t                      fGammaOrPi0;               ///< This tells me whether the correltation and the filling of histograms is done for gamma or pi0
  Bool_t                      fSEvMEv;                   ///< This option performs the analysis either for same event or for mixed event analysis
  Bool_t                      fDebug;			        ///< Can be set for debugging
  Bool_t                      fSavePool;                 ///< Defines whether to save output pools in a root file
  Bool_t                      fPlotQA;                   ///< plot additional QA histograms
  Bool_t                      fUseManualEventCuts;       ///< Use manual cuts if automatic setup is not available for the period

  //..Input histograms
  THnF                       *fHistEffGamma;             ///< ??input efficiency for trigger particles
  THnF                       *fHistEffHadron;            ///< ??input efficiency for associate particles

  //..Constants
  Double_t                    fRtoD;                     ///< conversion of rad to degree
  static const Int_t          kNIdentifier=3;            ///< number of different versions of the same histogram type, can later be used for centrality or mixed event eg.
  static const Int_t          kNvertBins=20;             ///< vertex bins in which the ME are mixed
  static const Int_t          kNcentBins=8;              ///< centrality bins in which the ME are mixed
  static const Int_t          kNoGammaBins=9;            ///< Bins in gamma pT
  static const Int_t          kNoZtBins=7;               ///< Bins in Zt
  static const Int_t          kNoXiBins=8;               ///< Bins in Xi
  Double_t                    fArray_G_Bins[10];         ///< 10=kNoGammaBins+1
  Double_t                    fArray_ZT_Bins[8];         ///< 8=kNoZtBins+1
  Double_t                    fArray_XI_Bins[9];         ///< 9=kNoXiBins+1
  Double_t                    fArrayNVertBins[21];       ///< 21=kNvertBins+1

  //..cuts
  Double_t                    fClShapeMin;               ///< Minimum cluster shape
  Double_t                    fClShapeMax;               ///< Maximum cluster shape
  Int_t                       fMaxNLM;                   ///< Maximum number of local maxima
  Bool_t                      fRmvMTrack;                ///< Switch to enable removing clusters with a matched track
  Double_t                    fTrackMatchEta;            ///< eta range in which a track is called a match to a cluster
  Double_t                    fTrackMatchPhi;            ///< phi range in which a track is called a match to a cluster
  //..Event pool variables
  TAxis                      *fMixBCent;                 ///< Number of centrality bins for the mixed event
  TAxis                      *fMixBZvtx;                 ///< Number of vertex bins for the mixed event
  AliEventPoolManager        *fPoolMgr;                  ///< event pool manager
  Int_t                       fTrackDepth;               ///<  #tracks to fill pool
  Int_t                       fPoolSize;                 ///<  Maximum number of events
  vector<vector<Double_t> >   fEventPoolOutputList;      //!<! ???vector representing a list of pools (given by value range) that will be saved
  //..Event selection types
  UInt_t                      fTriggerType;              ///<  Event types that are used for the trigger (gamma or pi0)
  UInt_t                      fMixingEventType;          ///<  Event types that are used for the tracks in the mixed event
  UInt_t                      fCurrentEventTrigger;      //!<! Trigger of the current event

  // MC stuff
  Bool_t                      fParticleLevel;            ///< Set particle level analysis
  Bool_t                      fIsMC;                     ///< Trigger, MC analysis
  UInt_t                      fAODfilterBits[2];         ///< AOD track filter bit map

  // Other stuff
  TList                      *fEventCutList;           //!<! Output list for event cut histograms
  TList                      *fOutputListQA;           //!<! Output list

  // Histograms -

  TH2                       *fHistClusPairInvarMasspT; //!<! Tyler's histogram
  TH1 					    *fHistPi0;                 //!<! Tyler's histogram
  TH2                       *fMAngle;                  //!<! Tyler's histogram
  TH2                       *fPtAngle;                 //!<! Tyler's histogram


  TRandom3    *fRand;//!<! Random number generator.  Initialzed by rot background
  TH1         *fClusEnergy; //!<! Energy of clusters accepted for pi0 analysis
  Bool_t       fDoRotBkg; ///< Whether or not to calculate the rotational background
  Int_t        fNRotBkgSamples; ///< How many samples to use in the rotational background
  THnSparseF  *fPi0Cands; //!<! Michael's THnSparse for pi0 Candidates


  TH1 					    *fHistEvsPt;               //!<! E vs pT
  TH1 					   **fHistBinCheckPt;          //!<! plot Pt distribution for ideal binning
  TH1 					   **fHistBinCheckZt;          //!<! plot Zt distribution for ideal binning
  TH1 					   **fHistBinCheckXi;          //!<! plot Xi distribution for ideal binning
  TH2                      **fHistDEtaDPhiGammaQA;     //!<! Distribution of gammas in delta phi delta eta
  TH2                      **fHistDEtaDPhiTrackQA;     //!<! Distribution of tracks in delta phi delta eta
  TH2                      **fHistClusterTime;         //!<! Cluster time vs energy
  THnSparseF                *fCorrVsManyThings;        //!<! Thn sparse filled with delta phi, delta eta,Eg,zt,xi,vertex Z,centrality...
  THnSparseF                *fClusterProp;             //!<! Thn sparse filled with cluster properties
  TH2                	    *fHPoolReady;              //!<! Check how many Jobs start mixing
  //

 private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 11) // Class to analyse gamma hadron correlations
};
#endif
