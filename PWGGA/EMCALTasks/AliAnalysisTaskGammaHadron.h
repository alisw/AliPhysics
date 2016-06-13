#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$
#include "AliAnalysisTaskEmcal.h"
#include "AliEventPoolManager.h"
#include <THn.h>
//#include "AliPool.h"

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;
class AliEvtPoolManager;

using std::vector;

class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcal {
 public:
	AliAnalysisTaskGammaHadron();
	AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputSameEventAnalysis);
virtual ~AliAnalysisTaskGammaHadron();

  //setters for the analysis, currently not many implemented
  void                        SetEffHistGamma(THnF *h)                              { fHistEffGamma       = h        ; }
  void                        SetEffHistHadron(THnF *h)                             { fHistEffHadron      = h        ; }
  void                        SetSavePool(Bool_t input)                             { fSavePool           = input    ; }


  //Functions for mixed event purposes
  void                        SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager*        GetEventPoolManager()                                 {return fPoolMgr;}
  // Set which pools will be saved
  void                        AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt);

 protected:

  void                        InitArrays()                                                  ;
  // EMCal framework functions
  Bool_t                      Run()                                                         ;
  void                        ExecOnce()         										  ;
  void                        UserCreateOutputObjects()                                     ;

  //Functions for mixed event purposes
  void                        InitEventMixer()											  ;
  TObjArray*                  CloneToCreateTObjArray(AliParticleContainer* tracks)          ;

  Bool_t                      FillHistograms()                                              ;
  Int_t                       CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  Int_t                       CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  void                        FillGhHisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Anglecut, Double_t Weight);
  void                        FillQAHisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut);
  Bool_t                      AccClusterForAna(AliVCluster* cluster)                        ;
  //Delta phi does also exist in AliAnalysisTaskEmcal. It is overwritten here (ask Raymond)
  Double_t                    DeltaPhi(TLorentzVector ClusterVec,AliVParticle* TrackVec)    ;
  Double_t                    GetEff(TLorentzVector ParticleVec)                            ;
  Bool_t                      fGammaOrPi0;               // This tells me whether the correltation and the filling of histograms is done for gamma or pi0
  Bool_t                      fDoMixing;                 // This option enables mixed events being used in the analysi
  Bool_t                      fDebug;			        // Can be set for debugging
  Bool_t                      fUsePerTrigWeight;		    // Sets whether you want to look at absolute yields or per trigger yields
  Bool_t                      fSavePool;                 // Defines whether to save output pools in a root file

  // Input histograms
  THnF                       *fHistEffGamma;             // input efficiency for trigger particles
  THnF                       *fHistEffHadron;            // input efficiency for associate particles

  // Constants
  Double_t                    fRtoD;                     // conversion of rad to degree
  static const Int_t          kNIdentifier=3;            // number of different versions of the same histogram type, can later be used for centrality or mixed event eg.
  static const Int_t          kNDPhistos=31;             // =  nbins[0];
  vector<Int_t>               fVector_G_Bins;            // vector that contains the bins of the G historgram
  vector<Double_t>            fVector_ZT_Bins;           // vector that contains the bins of the Zt historgram
  vector<Double_t>            fVector_XI_Bins;           // vector that contains the bins of the Xi historgram
  Double_t                    fZtStep;                   // Bin width for the zT histograms
  Double_t                    fXiStep;                   // Bin width for the Xi histograms

  // Event pool variables
  TAxis                      *fMixBCent;                 //! Number of centrality bins for the mixed event
  TAxis                      *fMixBZvtx;                 //! Number of vertex bins for the mixed event
  AliEventPoolManager        *fPoolMgr;                  //! event pool manager
  Int_t                       fTrackDepth;               //  #tracks to fill pool
  Int_t                       fPoolSize;                 //  Maximum number of events
  vector<vector<Double_t> >   fEventPoolOutputList;      //  vector representing a list of pools (given by value range) that will be saved
  // Event selection types
  UInt_t                      fTriggerType;              ///<  Event types that are used for the trigger (gamma or pi0)
  UInt_t                      fMixingEventType;          ///<  Event types that are used for the tracks in the mixed event
  UInt_t                      fCurrentEventTrigger;      ///<  Trigger of the current event

  // MC stuff
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  UInt_t                      fAODfilterBits[2];         // AOD track filter bit map

  // Other stuff
  TList                      *fOutputList1;            //! Output list
  TList                      *fOutputList2;            //! Output list
  TList                      *fOutputList3;            //! Output list
  TList                      *fOutputListGamma;        //! Output list
  TList                      *fOutputListXi;           //! Output list
  TList                      *fOutputListZeta;         //! Output list
  TList                      *fOutputListQA;           //! Output list

  // Histograms -
  TH1  					    *fHistNoClusPt;            //! No of calorimeter Clusters as a function of p_T
  TH1					   **fHistNoClusPtH;           //! No of calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH1 					    *fHistPi0;                 //!
  TH2					   **fHistDEtaDPhiG[3];        //! No of g-h pairs in the deta eta delta phi plane for certain gamma energies
  TH2					   **fHistDEtaDPhiZT[3];       //! No of g-h pairs in the deta eta delta phi plane for certain zT values
  TH2					   **fHistDEtaDPhiXI[3];       //! No of g-h pairs in the deta eta delta phi plane for certain Xi values
  TH2                      **fHistDEtaDPhiGammaQA;     //! Distribution of gammas in delta phi delta eta
  TH2                      **fHistDEtaDPhiTrackQA;     //! Distribution of tracks in delta phi delta eta

  TH1					   **fHistptAssHadron[3];      //! pt distributions of the associated hadron in a certain p_t bin of the gamma
  TH1					   **fHistDpGh[3];             //! delta phi g-h distribution fro a given p_t gamma bin
  TH2                	    *fHPoolReady;              //! Check how many Jobs start mixing
  //
  //

 private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 8) // Class to analyse gamma hadron correlations
};
#endif
