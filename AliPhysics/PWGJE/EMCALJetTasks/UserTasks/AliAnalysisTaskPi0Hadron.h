#ifndef ALIANALYSISTASKPI0HADRON_H
#define ALIANALYSISTASKPI0HADRON_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskPi0Hadron : public AliAnalysisTaskEmcalJet {
 public:
	AliAnalysisTaskPi0Hadron();
  virtual ~AliAnalysisTaskPi0Hadron();

  void                        UserCreateOutputObjects();

  //Set things for the analyis
  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetMaxCellsInCluster(Int_t b)                        { fMaxCellsInCluster  = b          ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s          ; }
  void                        SetMC(Bool_t m)                                      { fIsMC               = m          ; }
  void                        SetAdditionalCentEst(const char* meth2, const char* meth3="") { fCentMethod2 = meth2; fCentMethod3 = meth3; }
  void                        SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)         { fAODfilterBits[0]   = b0  ; fAODfilterBits[1] = b1  ; }

 protected:

  void                        ExecOnce()         										  ;
  Bool_t                      RetrieveEventObjects()                                        ;
  Bool_t                      FillHistograms()                                              ;
  Int_t                       DoCellLoop(Float_t &sum)                                      ;
  Int_t                       DoClusterLoop(Float_t &sum, AliVCluster* &leading)            ;
  Int_t                       DoTrackLoop(Float_t &sum, AliVParticle* &leading)             ;
  Bool_t                      AccClusterForAna(AliVCluster* cluster)                        ;
  Double_t                    DeltaPhi(TLorentzVector ClusterVec,AliVParticle* TrackVec)    ;
  void                        Fill_GH_Hisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Anglecut);

  //Constants
  Double_t                    fRtoD;                     // conversion of rad to degree

  // Cuts
  Float_t                     fCellEnergyCut;            // Energy cell cut
  Int_t                       fMaxCellsInCluster;        // Maximum number (approx) of cells in a cluster

  // MC stuff
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  TString                     fCentMethod2;              // Centrality method 2
  TString                     fCentMethod3;              // Centrality method 3
  UInt_t                      fAODfilterBits[2];         // AOD track filter bit map

  //Other stuff
  TList                      *fOutputList1;            //! Output list
  TList                      *fOutputList2;            //! Output list
  TList                      *fOutputList3;            //! Output list

  // Histograms
  TH1  						*fHistNoClus_pt;           //! No of calorimeter Clusters as a function of p_T
  TH1  						*fHistNoClus_pt_tr;        //! No of calorimeter trigger Clusters as a function of p_T
  TH1					    *fHistNoClus_ptH;          //! No of calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH1					    *fHistNoClus_ptH_tr;       //! No of calorimeter trigger Clusters as a function of p_T with a hadron in the second hemisphere
  TH2					    *fHist_DetaDphi;           //! No of g-h pairs in the deta eta delta phi plane
  TH1 						*fHistpi0;                 //!

  TH1                                       *fHistClusPairInvarMass;   //!    -<()>-
  TH2                                       *fHistClusPairInvarMasspT;   //!
  TH2                                       *fHistClusPairInvarMassE;   //!
  TH2                                       *fHistClusPairInvarMassPlay;   //!
  TH2                                       *fHistClusPairInvarMasspTSlice;   //!
  TH1                                       *fHistReadout;   //!

  TH1					   **fHistpt_assHadron;        //! pt distributions of the associated hadron in a certain p_t bin of the gamma
  TH1					   **fHistpt_assHadron_tr;     //! pt distributions of the associated hadron in a certain p_t bin of the gamma that triggered the event
  TH1					   **fHist_DP_gh;              //! delta phi g-h distribution fro a given p_t gamma bin
  //
  //

 private:
  AliAnalysisTaskPi0Hadron(const AliAnalysisTaskPi0Hadron&);            // not implemented
  AliAnalysisTaskPi0Hadron &operator=(const AliAnalysisTaskPi0Hadron&); // not implemented

  ClassDef(AliAnalysisTaskPi0Hadron, 4) // Class to analyse gamma hadron correlations
};
#endif
