#ifndef ALIANALYSISTASKEMCALJETPROPERTIES_H
#define ALIANALYSISTASKEMCALJETPROPERTIES_H
/**
 * \file AliAnalysisTaskEmcalJetProperties.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetProperties
 *
 * In this header file the class AliAnalysisTaskEmcalJetProperties is declared.
 * This is a task used for the measurement of charged jet properties.
 *
 * \author Sidharth Kumar Prasad <sidharth.kumar.prasad@cern.ch>, Bose Institute
 * \author Debjani Banerjee <banerjee.debjani@cern.ch>, Bose Institute
 * \author Prottoy Das <prottoy.das@cern.ch>, Bose Institute
 * 
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// This code is the modified version of AliAnalysisTaskJetProperties.h


#include <TF1.h>                      
class TString;                        
class TF1;                            
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TProfile;
class AliRhoParameter;                

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEmcalJetProperties
 * \brief Implementation of a sample jet analysis task.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 */

class AliAnalysisTaskEmcalJetProperties : public AliAnalysisTaskEmcalJet {
 public:
  
  AliAnalysisTaskEmcalJetProperties();
  AliAnalysisTaskEmcalJetProperties(const char *name);
  virtual ~AliAnalysisTaskEmcalJetProperties();
  
  void UserCreateOutputObjects();
  void Terminate(Option_t *option);
  

 protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();

  Bool_t IsRec = kTRUE;//True for matched rec jet list, false for matched gen jet list

  TClonesArray *fJets = new TClonesArray();     //!<! jets
  UInt_t fNExclLeadJets = 2; // number of leading jets to be excluded from the median calculation

  Double_t BkgSubtracted(Int_t KtJetCont, Int_t AktJetCont);
  Double_t GetNch(Bool_t Isrec);
  Double_t GetNchUEPC(Bool_t Isrec);
  Double_t GetFF(Bool_t Isrec, AliEmcalJet *jet, Int_t iTrack);
  Double_t GetFFUEPC(Bool_t Isrec, AliEmcalJet *jet, TList *UETrkList, Int_t iTrack);
  Double_t* GetPtSum(Bool_t Isrec, AliEmcalJet *jet, Double_t* PtSum);
  Double_t* GetPtSumUE(Bool_t Isrec, TList *UETrkList, Double_t etaTilted, Double_t phiTilted, Double_t* PtSum);
  Int_t GetMLeadingJet();
  
  void SetJetContainerBase(Int_t c){ fContainerBase = c;}
  void SetJetContainerTag(Int_t c){ fContainerTag  = c;}
  Int_t fContainerBase = -1;
  Int_t fContainerTag = -1;
  
  Int_t fContainerKtRec = -1;
  Int_t fContainerKtGen = -1;
  Int_t fContainerAktRec = -1;
  Int_t fContainerAktGen = -1;
  void SetJetContainerKtRec(Int_t c){ fContainerKtRec = c;}
  void SetJetContainerKtGen(Int_t c){ fContainerKtGen = c;}
  void SetJetContainerAktRec(Int_t c){ fContainerAktRec = c;}
  void SetJetContainerAktGen(Int_t c){ fContainerAktGen = c;}

  const Double_t fMaxDist = 0.24;
  void MatchJets(Int_t c1, Int_t c2, Double_t maxDist, TList* matchedListJet1, TList* matchedListJet2);
  // Int_t MyMatchJets(Int_t c1, Int_t c2, Double_t maxDist, TList* matchedListParton, TList* matchedListJet);

  void AllocateJetHistograms();
  void  FillJetProperties();
  void  FillJetShape();
  void  FillJetShapeUE2(const Double_t);
  void  GetTracksTiltedwrpJetAxis(Float_t alpha, 
				  TList* outputlist, 
				  const AliEmcalJet* jet, 
				  Double_t radius,Double_t& sumPt,Int_t& Ntracks_in_UE);
  void  FillJetShapeUE3(TList* outputlist, Double_t EtaTilt, Double_t PhiTilt);
  void RandomConeCalculation(Bool_t CheckOverlappingWithRC);
  //void RandomConeCalculation();
  Bool_t      IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);
  Double_t fCfactor = 1.0;
  void FillMy2DResponseMatrices();
  //  void FillMy2DResponseMatricesTest();


  THistManager fHistManager;///< Histogram manager

 private:
  AliAnalysisTaskEmcalJetProperties(const AliAnalysisTaskEmcalJetProperties&)           ; 
  AliAnalysisTaskEmcalJetProperties &operator=(const AliAnalysisTaskEmcalJetProperties&); 
  TList*    fTrackListUE; //! List of tracks in jet cone UE

  TList *matchedGenJetList; //! List of matched Gen level Jets
  TList *matchedRecJetList; //! List of matched Rec level Jets

  //FillJetProperties
  TH1F* fh1PtJet1D;          //!jet pt

  //FillJetShape, FillJetShapeUE2 (PC), FillJetShapeUE3 (RC)
  TH1F* fh1PtLeadingJet;          //!highest jet pt

  TH1F* fh1PtLeadingJet_PtHardBinNch;
  TH1F* fh1PtLeadingJet_PtHardBinFF;


  TH1F* fh1PtSumInJetConeUE2;
  TH1F* fh1PtSumInJetConeUE3;
  TH2F* fh2NtracksLeadingJet;
  TH2F* fh2NtracksLeadingJetUE2;
  TH2F* fh2NtracksLeadingJetUE3;
  TProfile* fProNtracksLeadingJet;
  TProfile* fProNtracksLeadingJetUE2;
  TProfile* fProNtracksLeadingJetUE3;
  TH2F* fh2PtTrack;
  TH2F* fh2PtTrackUE2;
  TH2F* fh2PtTrackUE3;
  TH2F* fh2FF;
  TH2F* fh2FFUE2;
  TH2F* fh2FFUE3;
  TH2F* fh2Ksi;
  TH2F* fh2KsiUE2;
  TH2F* fh2KsiUE3;
  TH2F* fh2DelR80pcPt;
  TH2F* fh2DelR80pcPtUE2;
  TH2F* fh2DelR80pcPtUE3;
  TProfile* fProDelR80pcPt;
  TProfile* fProDelR80pcPtUE2;
  TProfile* fProDelR80pcPtUE3;
  TH3F* fh3PtDelRPtSum;
  TH3F* fh3PtDelRPtSumUE2;
  TH3F* fh3PtDelRPtSumUE3;
  TProfile* fProDelRPtDenUEN2;
  TProfile* fProDelRPtDenUEN3;
  TProfile* fProDelRPtSum[21];         //!Pt sum vs R
  TProfile* fProDelRPtSumUE2[21];         //!Pt sum vs R
  TProfile* fProDelRPtSumUE3[21];         //!Pt sum vs R
  TProfile* fProDelRPtDenUE2[21];         //!Pt sum vs R
  TProfile* fProDelRPtDenUE3[21];         //!Pt sum vs R

  TH1F *fh1DelRRCJet1;
  TH1F *fh1DelRRCJet2;
  TH1F *fh1nRC;

  //BkgSubtracted
  TProfile* hProDelRPtDenUEN;
  TProfile* hhProDelRPtSumUE[21];       //!Pt sum vs R
  TProfile* hhProDelRPtDenUE[21];       //!Pt density vs R
  TH1F*     fh1rho;                   //!rho distribution
  TH1F*     fh1rhoN;                   //!rho for Ntracks distribution
  TH1F*     h1LeadingJetPtWUE;        //!highest jet pt
  TH1F*     h1LeadingJetPtWoUEWoCAll;       //!highest jet pt after UE subtraction
  TProfile* fProNtracksLeadingJetWUE; //!number of tracks in leading jet
  TProfile* fProNtracksLeadingJetUEWoCAll; //!number of tracks in leading jet (UE)
  TProfile* fProNtracksLeadingJetWoUEWoCAll;//!number of tracks in leading jet after UE subtraction
  TProfile* fProNtracksLeadingJetUEwrtSubtPtWoC; //!number of tracks in leading jet (UE) wrt subt jet pt
  TProfile* fProNtracksLeadingJetWoUEwrtSubtPtWoC;//!number of tracks in leading jet after UE subtraction wrt subt jet pt

  TH1F      *h1Cfactor; // C factor
  TH1F*     h1LeadingJetPtWoUEWCAll;       //!highest jet pt after UE subtraction
  TProfile* fProNtracksLeadingJetUEWCAll; //!number of tracks in leading jet (UE)
  TProfile* fProNtracksLeadingJetWoUEWCAll;//!number of tracks in leading jet after UE subtraction
  TProfile* fProNtracksLeadingJetUEwrtSubtPtWC; //!number of tracks in leading jet (UE) wrt subt jet pt
  TProfile* fProNtracksLeadingJetWoUEwrtSubtPtWC;//!number of tracks in leading jet after UE subtraction wrt subt jet pt

  THnSparseD* h4ResNch;
  TH2D* h2NchGen;
  TH2D* h2NchRec;

  THnSparseD* h4ResNchUE;
  TH2D* h2NchGenUE;
  TH2D* h2NchRecUE;

  THnSparseD* h4ResFF;
  TH2D* h2FFGen;
  TH2D* h2FFRec;
  THnSparseD* h4ResFFM;
  TH2D* h2FFMFake;
  TH2D* h2FFMMiss;
  TH2D* h2FFMGen;
  TH2D* h2FFMRec;

  THnSparseD* h4ResFFMUE;
  TH2D* h2FFMUEFake;
  TH2D* h2FFMUEMiss;
  TH2D* h2FFMUEGen;
  TH2D* h2FFMUERec;

  THnSparseD* h5ResPtSum;
  THnSparseD* h5ResPtSumUE;

  TH2D* h2ResJetPt;
  TH2D* h2ResJetPtNch;

   //______________________HM_________________________________//                                                                                              
  Double_t  fMultV0A;                                  //!  mult. V0A                                                                                        
  Double_t  fMultV0C;                                  //!  mult. V0C                                                                                        
  Double_t  fMultV0M;                                  //!  mult. V0A+V0C                                                                                    

  TH1D* fhistMultV0A;
  TH1D* fhistMultV0C;
  TH1D* fhistMultV0M;

  //______________________HM_________________________________//                                                                                               


  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetProperties, 1);
  /// \endcond
};
#endif

