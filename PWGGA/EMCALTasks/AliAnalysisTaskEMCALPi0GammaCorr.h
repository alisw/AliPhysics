#ifndef AliAnalysisTaskEMCALPi0GammaCorr_H
#define AliAnalysisTaskEMCALPi0GammaCorr_H

// $Id$
#include "AliAnalysisTaskEmcal.h"
#include "AliEventPoolManager.h"
#include <THn.h>
#include "AliStaObjects.h"
#include "AliEventCuts.h"

class TH1;
class TH2;
class TH3;
class THnSparse;
class THnSparse;
class AliVVZERO;
class AliEvtPoolManager;

using std::vector;

class AliAnalysisTaskEMCALPi0GammaCorr : public AliAnalysisTaskEmcal {
 public:
	AliAnalysisTaskEMCALPi0GammaCorr();
	AliAnalysisTaskEMCALPi0GammaCorr(Bool_t InputSameEventAnalysis);
virtual ~AliAnalysisTaskEMCALPi0GammaCorr();

  void SetEffHistGamma(THnF *h)                              { fHistEffGamma    = h      ; }
  void SetEffHistHadron(THnF *h)                             { fHistEffHadron   = h      ; }
  void SetSavePool(Bool_t input)                             { fSavePool        = input  ; }
  void SetEvtTriggerType(UInt_t input)                       { fTriggerType     = input  ; }
  void SetEvtMixType(UInt_t input)                           { fMixingEventType = input  ; }
  void SetUseManualEvtCuts(Bool_t input)                     { fUseManualEventCuts = input;}

  void SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager*        GetEventPoolManager()                                 {return fPoolMgr;}
  void AddEventPoolsToOutput(double minCent, double maxCent,  double minZvtx, double maxZvtx, double minPt, double maxPt);
  
  private:
  AliEventCuts fEventCuts;                  
  
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
  
  void                        FillClusterHisto(AliVCluster* cluster, THnSparse* histo);
  void                        FillPionHisto(AliVCluster* cluster1, AliVCluster* cluster2, THnSparse* histo);
  void                        FillPionCorrelation(AliVCluster* cluster1, AliVCluster* cluster2, AliVParticle* track, THnSparse* histo, double weight);
  void                        FillPhotonCorrelation(AliVCluster* cluster, AliVParticle* track, THnSparse* histo, double weight);
  int                         CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, double Weight);
  Bool_t                      PassedCuts(AliVCluster* caloCluster);
  double                      GetIsolation_Track(AliVCluster* cluster);
  
  TObjArray*                  CloneClustersTObjArray(AliClusterContainer* clusters)          ;
  double                    GetEff(AliTLorentzVector ParticleVec)                         ;

  Bool_t                      fSavePool;                 ///< Defines whether to save output pools in a root file
  Bool_t                      fUseManualEventCuts;       ///< Use manual cuts if automatic setup is not available for the period
  
  //..Input histograms
  THnF                       *fHistEffGamma;             ///< ??input efficiency for trigger particles
  THnF                       *fHistEffHadron;            ///< ??input efficiency for associate particles

  //..Constant
  static const int          kNvertBins=20;             ///< vertex bins in which the ME are mixed
  static const int          kNcentBins=8;              ///< centrality bins in which the ME are mixed
  double                    fArrayNVertBins[21];       ///< 21=kNvertBins+1

  //..Event pool variables
  TAxis                      *fMixBCent;                 ///< Number of centrality bins for the mixed event
  TAxis                      *fMixBZvtx;                 ///< Number of vertex bins for the mixed event
  AliEventPoolManager        *fPoolMgr;                  ///< event pool manager
  int                       fTrackDepth;               ///<  #tracks to fill pool
  int                       fPoolSize;                 ///<  Maximum number of events
  vector<vector<double> >   fEventPoolOutputList;      //!<! ???vector representing a list of pools (given by value range) that will be saved

  UInt_t                      fTriggerType;              ///<  Event types that are used for the trigger (gamma or pi0)
  UInt_t                      fMixingEventType;          ///<  Event types that are used for the tracks in the mixed event
  UInt_t                      fCurrentEventTrigger;      //!<! Trigger of the current event

  THnSparse                 *h_Track;                   //!<!
  THnSparse                 *h_Cluster;                 //!<!
  THnSparse                 *h_ClusterTrack;                 //!<! THnSparse with info on cluster and track.
  THnSparse                 *h_ClusterTrack_Mixed;                 //!<!
  THnSparse                 *h_Pi0;                 //!<!
  THnSparse                 *h_Pi0Track;                 //!<! THnSparse with info on pi0 and track.
  THnSparse                 *h_Pi0Track_Mixed;                 //!<!

  TList                      *fEventCutList;           //!<! Output list for event cut histograms
 // TList                      *OutputList;            //!<! Output list
  
  const static int nEvt      =   10;//30; // mixing "depth"

  private:
  AliAnalysisTaskEMCALPi0GammaCorr(const AliAnalysisTaskEMCALPi0GammaCorr&);            // not implemented
  AliAnalysisTaskEMCALPi0GammaCorr &operator=(const AliAnalysisTaskEMCALPi0GammaCorr&); // not implemented
  ClassDef(AliAnalysisTaskEMCALPi0GammaCorr, 10) // Class to analyse gamma hadron correlations
};
#endif



