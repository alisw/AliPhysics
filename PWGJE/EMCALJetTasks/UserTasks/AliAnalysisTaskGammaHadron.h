#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEventPoolManager.h"
#include <THn.h>
//#include "AliPool.h"

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;
class AliEvtPoolManager;


class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcalJet {
 public:
	AliAnalysisTaskGammaHadron();
  virtual ~AliAnalysisTaskGammaHadron();

  void                        UserCreateOutputObjects();

  //Set things for the analyis
  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut      ; }
  void                        SetMaxCellsInCluster(Int_t b)                        { fMaxCellsInCluster  = b        ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s        ; }
  void                        SetMC(Bool_t m)                                      { fIsMC               = m        ; }
  void                        SetAdditionalCentEst(const char* meth2)              { fCentMethod_alt     = meth2    ; }
  void                        SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)         { fAODfilterBits[0]   = b0  ; fAODfilterBits[1] = b1  ; }
  void                        SetEffHistGamma(THnF *h)                             { fHistEff_Gamma      = h        ; }
  void                        SetEffHistHadron(THnF *h)                            { fHistEff_Hadron     = h        ; }

 protected:

  void                        ExecOnce()         										  ;
  void                        InitEventMixer()											  ;
  Bool_t                      RetrieveEventObjects()                                        ;
  TObjArray*                  CloneToCreateTObjArray(AliParticleContainer* tracks)          ;
  Bool_t                      Run()                                                         ;
  Bool_t                      FillHistograms()                                              ;
  Int_t                       CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  Bool_t                      AccClusterForAna(AliVCluster* cluster)                        ;
  Double_t                    DeltaPhi(TLorentzVector ClusterVec,AliVParticle* TrackVec)    ;
  void                        Fill_GH_Hisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Anglecut, Double_t Weight);

  // Input histograms
  THnF                      *fHistEff_Gamma;            // input efficiency for trigger particles
  THnF                      *fHistEff_Hadron;           // input efficiency for associate particles

  // Constants
  Double_t                    fRtoD;                     // conversion of rad to degree

  TAxis*                      fMixBCent;                 // Number of centrality bins for the mixed event
  TAxis*                      fMixBZvtx;                 // Number of vertex bins for the mixed event
  TString                     fCentMethod_alt;           // alternative centrality selection method
  Double_t                    fCent_alt;                 // alternative centrality
  AliEventPoolManager        *fPoolMgr;                 //! event mixer
  Int_t                       fTrackDepth;               //  #tracks to fill pool
  Int_t                       fPoolSize;                 //  Maximum number of events

  // Cuts
  Float_t                     fCellEnergyCut;            // Energy cell cut
  Int_t                       fMaxCellsInCluster;        // Maximum number (approx) of cells in a cluster

  // MC stuff
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  UInt_t                      fAODfilterBits[2];         // AOD track filter bit map

  // Other stuff
  TList                      *fOutputList1;            //! Output list
  TList                      *fOutputList2;            //! Output list
  TList                      *fOutputList3;            //! Output list

  // Histograms -
  TH1  					    *fHistNoClus_pt_Trigger;   //! No of calorimeter Clusters as a function of p_T
  TH1  					    *fHistNoClus_pt;           //! No of calorimeter Clusters as a function of p_T
  TH1					   **fHistNoClus_ptH;          //! No of calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH2					   **fHist_dEta_dPhi;          //! No of g-h pairs in the deta eta delta phi plane
  TH2					   **fHist_dEta_dPhi_low;      //! No of g-h pairs in the deta eta delta phi plane
  TH2					   **fHist_dEta_dPhi_med;      //! No of g-h pairs in the deta eta delta phi plane
  TH2					   **fHist_dEta_dPhi_high;     //! No of g-h pairs in the deta eta delta phi plane
  TH1 					    *fHistpi0;                 //!

  TH1					  **fHistpt_assHadron[3];      //! pt distributions of the associated hadron in a certain p_t bin of the gamma
  TH1					  **fHist_DP_gh[3];            //! delta phi g-h distribution fro a given p_t gamma bin
  TH2                	   *fHPoolReady;               //! Check how many Jobs start mixing
  //
  //

 private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 6) // Class to analyse gamma hadron correlations
};
#endif
