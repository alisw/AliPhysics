#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcalJet {
 public:
	AliAnalysisTaskGammaHadron();
  virtual ~AliAnalysisTaskGammaHadron();

  void                        UserCreateOutputObjects();

  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s          ; }
  void                        SetMC(Bool_t m)                                      { fIsMC               = m          ; }
  void                        SetAdditionalCentEst(const char* meth2, const char* meth3="") { fCentMethod2 = meth2; fCentMethod3 = meth3; }
  void                        SetDoV0QA(Int_t b)                                   { fDoV0QA             = b          ; }
  void                        SetDoEPQA(Int_t b)                                   { fDoEPQA             = b          ; }
  void                        SetMaxCellsInCluster(Int_t b)                        { fMaxCellsInCluster  = b          ; }
  void                        SetDoLeadingObjectPosition(Int_t b)                  { fDoLeadingObjectPosition = b     ; }
  void                        SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)         { fAODfilterBits[0]   = b0  ; fAODfilterBits[1] = b1  ; }

 protected:

  void                        ExecOnce()        											  ;
  Bool_t                      RetrieveEventObjects()                                        ;
  Bool_t                      FillHistograms()                                              ;

  Int_t                       DoCellLoop(Float_t &sum)                                      ;
  Int_t                       DoClusterLoop(Float_t &sum, AliVCluster* &leading)            ;
  Int_t                       DoTrackLoop(Float_t &sum, AliVParticle* &leading)             ;

  Float_t                     fCellEnergyCut;            // Energy cell cut
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  TString                     fCentMethod2;              // Centrality method 2
  TString                     fCentMethod3;              // Centrality method 3
  Int_t                       fDoV0QA;                   // Add V0 QA histograms
  Int_t                       fDoEPQA;                   // Add event plane QA histograms
  Int_t                       fDoLeadingObjectPosition;  // Add axis for leading object position (eta-phi)
  Int_t                       fMaxCellsInCluster;        // Maximum number (approx) of cells in a cluster
  UInt_t                      fAODfilterBits[2];         // AOD track filter bit map
  Double_t                    fCent2;                    //!Event centrality with method 2
  Double_t                    fCent3;                    //!Event centrality with method 3
  AliVVZERO                  *fVZERO;                    //!Event V0 object
  Double_t                    fV0ATotMult;               //!Event V0A total multiplicity
  Double_t                    fV0CTotMult;               //!Event V0C total multiplicity
 

  // ELIANE
  TH1  						*fHistNoClus_pt;           //! No of calorimeter Clusters as a function of p_T
  TH1  						*fHistNoClus_pt_tr;        //! No of calorimeter trigger Clusters as a function of p_T
  TH1					    *fHistNoClus_ptH;          //! No of calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH1					    *fHistNoClus_ptH_tr;       //! No of calorimeter trigger Clusters as a function of p_T with a hadron in the second hemisphere
  TH1					    *fHistNoClus_ptLeadH;      //! No of calorimeter Clusters as a function of p_T with a leading hadron in the second hemisphere
  TH1					    *fHistNoClus_Leadpt;       //! No of leading calorimeter Clusters as a function of p_T
  TH1					    *fHistNoClus_LeadptH;      //! No of leading calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH1					    *fHistNoClus_LeadptLeadH;  //! No of leading calorimeter Clusters as a function of p_T with a leading hadron in the second hemisphere

  TH1					    *fHistNoClus_xEH;          //! No of calorimeter Clusters as a function of x_E with a hadron in the second hemisphere
  TH1					    *fHistNoClus_LeadxEH;      //! No of calorimeter Clusters as a function of x_E with a leading hadron in the second hemisphere
  TH1					    *fHistNoClus_xELeadH;      //! No of leading calorimeter Clusters as a function of x_E with a hadron in the second hemisphere
  TH1					    *fHistNoClus_LeadxELeadH;  //! No of leading calorimeter Clusters as a function of x_E with a leading hadron in the second hemisphere
  TH1					   **fHistpt_assHadron;        //! pt distributions of the associated hadron in a certain p_t bin of the gamma
  TH1					   **fHistpt_assHadron_tr;     //! pt distributions of the associated hadron in a certain p_t bin of the gamma that triggered the event
  //
  //

 private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 3) // Quality task for Emcal analysis
};
#endif
