#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$
#include "AliAnalysisTaskEmcal.h"
#include "AliEventPoolManager.h"
#include <THn.h>
#include "AliStaObjects.h"//<<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
#include "AliEventCuts.h"
#include <THnSparse.h>

class TH1;
class TH2;
class TH3;
class AliVVZERO;
class AliEvtPoolManager;

using std::vector;

// class to store EMC hits    <<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
class EmcHitPi0 {

  TLorentzVector thishit;
  Short_t hittype; // 100 for bg event/real event, 1 for added pi0, 2 for added eta, -1 for not needed; from primary pi0: 101, from secondary pi0: 102, from K0: 103, from material: 104
  Byte_t smno;
  Int_t imo; // index of original mother in monte carlo stack
  Int_t pid; // particle ID
  Float_t weight; // weight from mother particle
  Bool_t bclean; // clean if only one contributor

  Int_t NCells;
  std::vector<int> CellRay;

   // UShort_t *CellIDArray;
  // Int_t GetNCells() {return NCells;}

  //  void SetCellIDArray(const UShort_t *cellid);
  // UShort_t GetCellIDArray() {return *CellIDArray;}

public:
  //virtual ~EmcHit();
  EmcHitPi0();
  friend class EmcEventPi0;
  friend class AliAnalysisTaskGammaHadron;

  void Print(){Printf("E=%.2f, type=%d, MoID=%d, PID=%d, w=%.3f",thishit.E(),hittype,imo,pid,weight);   }
};
// class to store old events
class EmcEventPi0 {

  //    Int_t fCenPercent;
  //    Int_t fVtx;
	Float_t TrigPhi; // phi of highest pT hit on EMCal
  Float_t TrigTheta; // eta of highest pT hit ...

  const static int nMaxHit = 1000;

  int nHits;
  EmcHitPi0 hit[nMaxHit];

public:
  EmcEventPi0();
  EmcEventPi0(const EmcEventPi0 &obj);
  //virtual ~EmcEvent();
  //    void SetGlobalInfo(const Int_t&, const Int_t&, const Int_t&, const Int_t&, const Double_t&, const Double_t&);
  void SetGlobalInfo(const Int_t&, const Float_t&, const Float_t&);
  int evsize() {return nHits;}
  void Reset();
  void Print();

  friend class AliAnalysisTaskGammaHadron;
};
//                              <<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>


class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcal {
 public:
	AliAnalysisTaskGammaHadron();
	AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputSameEventAnalysis, Bool_t InitMCorData);
virtual ~AliAnalysisTaskGammaHadron();

  //setters for the analysis, currently not many implemented
  void                        SetEffHistGamma(THnF *h)                              { fHistEffGamma    = h      ; }
  void                        SetEffHistHadron(THnF *h)                             { fHistEffHadron   = h      ; }
  void                        SetSavePool(Bool_t input)                             { fSavePool        = input  ; }
  void                        SetEvtTriggerType(UInt_t input)                       { fTriggerType     = input  ; }
  void                        SetEvtMixType(UInt_t input)                           { fMixingEventType = input  ; }
//  void                        SetCutsId(Id, cent, ptCl,Ecl,,.... UInt_t input)      { fMixingEventType = input  ; }
  void                        SetNLM(Int_t input)                                   { fMaxNLM = input;}
  void                        SetM02(Double_t inputMin,Double_t inputMax)           { fClShapeMin = inputMin; fClShapeMax = inputMax;}
  void                        SetRmvMatchedTrack(Bool_t input, Double_t dEta=-1, Double_t dPhi=-1) { fRmvMTrack  = input; fTrackMatchEta=dEta; fTrackMatchPhi=dPhi;}
  void                        SetUseManualEvtCuts(Bool_t input)                     { fUseManualEventCuts = input;}

  //Functions for mixed event purposes
  void                        SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager*        GetEventPoolManager()                                 {return fPoolMgr;}
  // Set which pools will be saved
  void                        AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt);
 private:
  AliEventCuts                fEventCuts;                   ///< event selection utility

 protected:

  virtual void                ProcessMC()                                                   ;
  void                        InitArrays()                                                 ;
  // overwritten EMCal framework functions
  Bool_t                      Run()                             	                          ;
  void                        ExecOnce()         									      ;
  Bool_t                      IsEventSelected()											  ;
  void                        UserCreateOutputObjects()        		                      ;

  //Functions for mixed event purposes
  void                        InitEventMixer()											  ;
  TObjArray*                  CloneToCreateTObjArray(AliParticleContainer* tracks)          ;

  Bool_t                      FillHistograms()                                              ;
  Int_t                       CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  Int_t                       CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  void                        FillGhHisograms(Int_t identifier,AliTLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Weight);
  void                        FillQAHisograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec);
  Bool_t                      AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster);
  Bool_t                      DetermineMatchedTrack(AliVCluster* caloCluster);
  //<<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
    TObjArray*                  CloneClustersTObjArray(AliClusterContainer* clusters)          ;
    void GetMulClassPi0(Int_t&);
    void AddMixEventPi0(const Int_t, const Int_t, const Int_t, Int_t&, const Float_t&, const Float_t&);
  //<<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>

  //..Delta phi does also exist in AliAnalysisTaskEmcal. It is overwritten here (ask Raymond)
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec) ;
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,Double_t phi_EVP)       ;
  Double_t                    GetEff(AliTLorentzVector ParticleVec)                         ;

  Bool_t                      fGammaOrPi0;               ///< This tells me whether the correltation and the filling of histograms is done for gamma or pi0
  Bool_t                      fDoMixing;                 ///< This option enables mixed events being used in the analysi
  Bool_t                      fMCorData;                 //<Are we looking at simulations or at the real thing
  Bool_t                      fDebug;			        ///< Can be set for debugging
  Bool_t                      fSavePool;                 ///< Defines whether to save output pools in a root file
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
  TList                      *fOutputList1;            //!<! Output list
  TList                      *fOutputListTrAs;         //!<! Output list
  TList                      *fOutputListGamma;        //!<! Output list
  TList                      *fOutputListXi;           //!<! Output list
  TList                      *fOutputListZeta;         //!<! Output list
  TList                      *fOutputListQA;           //!<! Output list

  // Histograms -

  TH1  					    *fHistNoClusPt;            //!<! ?No of calorimeter Clusters as a function of p_T
  TH1					   **fHistptAssHadronG[3];     //!<! pt distr. of the associated hadron as a function of Eg
  TH1					   **fHistptAssHadronZt[3];    //!<! pt distr. of the associated hadron as a function of Zt
  TH1					   **fHistptAssHadronXi[3];    //!<! pt distr. of the associated hadron as a function of Xi
  TH1					   **fHistptTriggG[3];         //!<! pt distr. of the trigger as a function of Eg
  TH1					   **fHistptTriggZt[3];        //!<! pt distr. of the trigger as a function of Zt
  TH1					   **fHistptTriggXi[3];        //!<! pt distr. of the trigger as a function of Xi

  TH2                       *fHistClusPairInvarMasspT; //!<! Tyler's histogram

  TH1 					    *fHistPi0;                 //!<!
  TH2                                       *fHistClusPairInvarMasspTMIX; //!<!
  TH2                                       *fHistClusPairInvarMasspTMIXolap; //!<!
  TH2                                       *fMAngle; //!<!
  TH2                                       *fMAngleMIX; //!<!
  TH2                                       *fPtAngle; //!<!
  TH2                                       *fPtAngleMIX; //!<!
  TH1 					    *fHistEvsPt;               //!<! E vs pT
  TH1 					   **fHistBinCheckPt;          //!<! plot Pt distribution for ideal binning
  TH1 					   **fHistBinCheckZt;          //!<! plot Zt distribution for ideal binning
  TH1 					   **fHistBinCheckXi;          //!<! plot Xi distribution for ideal binning
  TH2					   **fHistDEtaDPhiG[3][10];    //!<! No of g-h pairs in the deta eta delta phi plane for certain gamma energies
  TH2					   **fHistDEtaDPhiZT[3][8];    //!<! No of g-h pairs in the deta eta delta phi plane for certain zT values
  TH2					   **fHistDEtaDPhiXI[3][9];    //!<! No of g-h pairs in the deta eta delta phi plane for certain Xi values
  TH2                      **fHistDEtaDPhiGammaQA;     //!<! Distribution of gammas in delta phi delta eta
  TH2                      **fHistDEtaDPhiTrackQA;     //!<! Distribution of tracks in delta phi delta eta
  TH2                      **fHistCellsCluster;        //!<! Number of cells in cluster as function of energy
  TH2                      **fHistClusterShape;        //!<! Cluster shape vs energy
  TH2                      **fHistClusterShape0;       //!<! Cluster shape vs energy
  TH2                      **fHistClusterShape1;       //!<! Cluster shape vs energy
  TH2                      **fHistClusterShape2;       //!<! Cluster shape vs energy
  TH2                      **fHistClusterShape3;       //!<! Cluster shape vs energy
  TH2                      **fHistClusterShape4;       //!<! Cluster shape vs energy
  TH2                       *fHistMatchEtaPhiAllCl2;   //!<! matched track distance for 2 cell clusters
  TH2                       *fHistMatchEtaPhiAllCl3;   //!<! matched track distance for 3 cell clusters
  TH2                       *fHistMatchEtaPhiAllCl4;   //!<! matched track distance for 4 cell clusters
  TH2                      **fHistClusterTime;         //!<! Cluster time vs energy
  THnSparseF                *fCorrVsManyThings;        //!<! Thn sparse filled with delta phi, delta eta,Eg,zt,xi,vertex Z,centrality...
  THnSparseF                *fCorrVsManyThingsME;      //!<! Thn sparse filled with delta phi, delta eta,Eg,zt,xi,vertex Z,centrality...
  THnSparseF                *fClusterProp;             //!<! Thn sparse filled with cluster properties
  TH2                	    *fHPoolReady;              //!<! Check how many Jobs start mixing
  //
  //


  const static int nMulClass =   5;  // <<<<<><<<<<<<<<><<<<<<<<<<><<<<<<<<<<<<><<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<>
  const static int nZClass   =   3;
  const static int nPtClass = 1;
  int iEvt[nMulClass][nZClass][nPtClass];
  const static int nEvt      =   3;//30; // mixing "depth"

  EmcEventPi0 evt;
  EmcEventPi0 EmcEventList[nMulClass][nZClass][nPtClass][nEvt];

  EmcEventPi0 thisEvent;

 private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 10) // Class to analyse gamma hadron correlations
};
#endif
