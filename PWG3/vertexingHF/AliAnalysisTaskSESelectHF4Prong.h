#ifndef ALIANALYSISTASKSESELECTHF4PRONG_H
#define ALIANALYSISTASKSESELECTHF4PRONG_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSESelectHF4Prong
// AliAnalysisTaskSE for the selection of heavy-flavour decay candidates
// and creation of a stand-alone AOD for 4prong D0 decay
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//         F.Colamaria, fabio.colamaria@ba.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TChain.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliAODMCParticle.h"
class AliAODEvent;
class AliAnalysisManager;
class AliAnalysisVertexingHF;
class AliRDHFCuts;

class AliAnalysisTaskSESelectHF4Prong : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSESelectHF4Prong();
  AliAnalysisTaskSESelectHF4Prong(const char *name,AliRDHFCutsD0toKpipipi* cuts);
  virtual ~AliAnalysisTaskSESelectHF4Prong();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetPtBinH(Double_t* ptlimits);
  void PrintPtBinHandMCFlag();
  void SetMCTruth(int flag) {fMCTruth = flag;};
  Int_t GetMCTruth() const {return fMCTruth;}
  
  void AnalysisReflection(AliAODEvent* aodIn, AliAODRecoDecayHF4Prong* d, Int_t hyp01, Int_t hyp02, Int_t hyp03, Int_t hyp04);
  Int_t StudyMCTruth(TClonesArray* mcArray, AliAODRecoDecayHF4Prong* d);
  void FillReflHistos(AliAODRecoDecayHF4Prong* d, Int_t hyp01, Int_t hyp02, Int_t hyp03, Int_t hyp04, Int_t flagLayer1, Int_t flagLayer2, Int_t flagLayer3, Int_t flagLayer4);
  void FillReflD0Histos(AliAODRecoDecayHF4Prong* d, Int_t hyp01, Int_t hyp02, Int_t hyp03, Int_t hyp04, Int_t flagLayforD0_1, Int_t flagLayforD0_2, Int_t flagLayforD0_3, Int_t flagLayforD0_4);
  void FillReflD0barHistos(AliAODRecoDecayHF4Prong* d, Int_t hyp01, Int_t hyp02, Int_t hyp03, Int_t hyp04, Int_t flagLayforD0bar_1, Int_t flagLayforD0bar_2, Int_t flagLayforD0bar_3, Int_t flagLayforD0bar_4);
 
 private:

  AliAnalysisTaskSESelectHF4Prong(const AliAnalysisTaskSESelectHF4Prong &source);
  AliAnalysisTaskSESelectHF4Prong& operator=(const AliAnalysisTaskSESelectHF4Prong& source); 
  TClonesArray *fVerticesHFTClArr;     //! Array of heavy-flavour vertices
  TClonesArray *fCharm4ProngTClArr;        //! Array of D0->K3pi

  Double_t fPtBinH[6];			      //bin i has pt between values i and i+1
  Double_t fmassD0[2];			      //!To store invariant mass for D0 hypotheses
  Double_t fmassD0bar[2];		      //!To store invariant mass for D0bar hypotheses
  Int_t fSelected;			      //!Flag for selection of candidate
  Int_t fMCTruth;			      //flag for MC truth analysis

  TList *fOutput;                             //! list send on output slot 1
  TList *fOutput2;                            //! list send on output slot 2
  TList *fOutput3;                            //! list send on output slot 3
  TList *fOutput4;                            //! list send on output slot 4
  TList *fOutput5;                            //! list send on output slot 5
  TList *fOutputC;                            //! list send on output slot 6

  //output histograms
  TH1F *fhInvMassD0Sum10MevBin1;                        //! Invariant mass D01+D02 (good hyp) 10Mev BIN1
  TH1F *fhInvMassD0barSum10MevBin1;                     //! Invariant mass D0bar1+D0bar2 (good hyp) 10Mev
  TH1F *fhInvMassSumAll10MevBin1;                       //! Invariant mass superimpose (good hyp only) 10Mev
  TH1F *fhInvMassD0Sum5MevBin1;                         //! Invariant mass D01+D02 (good hyp) 5Mev
  TH1F *fhInvMassD0barSum5MevBin1;                      //! Invariant mass D0bar1+D0bar2 (good hyp) 5Mev
  TH1F *fhInvMassSumAll5MevBin1;                        //! Invariant mass superimpose (good hyp only) 5Mev

  TH1F *fhInvMassD0Sum10MevBin2;                        //! Invariant mass D01+D02 (good hyp) 10Mev BIN2
  TH1F *fhInvMassD0barSum10MevBin2;                     //! Invariant mass D0bar1+D0bar2 (good hyp) 10Mev
  TH1F *fhInvMassSumAll10MevBin2;                       //! Invariant mass superimpose (good hyp only) 10Mev
  TH1F *fhInvMassD0Sum5MevBin2;                         //! Invariant mass D01+D02 (good hyp) 5Mev
  TH1F *fhInvMassD0barSum5MevBin2;                      //! Invariant mass D0bar1+D0bar2 (good hyp) 5Mev
  TH1F *fhInvMassSumAll5MevBin2;                        //! Invariant mass superimpose (good hyp only) 5Mev

  TH1F *fhInvMassD0Sum10MevBin3;                        //! Invariant mass D01+D02 (good hyp) 10Mev BIN3
  TH1F *fhInvMassD0barSum10MevBin3;                     //! Invariant mass D0bar1+D0bar2 (good hyp) 10Mev
  TH1F *fhInvMassSumAll10MevBin3;                       //! Invariant mass superimpose (good hyp only) 10Mev
  TH1F *fhInvMassD0Sum5MevBin3;                         //! Invariant mass D01+D02 (good hyp) 5Mev
  TH1F *fhInvMassD0barSum5MevBin3;                      //! Invariant mass D0bar1+D0bar2 (good hyp) 5Mev
  TH1F *fhInvMassSumAll5MevBin3;                        //! Invariant mass superimpose (good hyp only) 5Mev

  TH1F *fhInvMassD0Sum10MevBin4;                        //! Invariant mass D01+D02 (good hyp) 10Mev BIN4
  TH1F *fhInvMassD0barSum10MevBin4;                     //! Invariant mass D0bar1+D0bar2 (good hyp) 10Mev
  TH1F *fhInvMassSumAll10MevBin4;                       //! Invariant mass superimpose (good hyp only) 10Mev
  TH1F *fhInvMassD0Sum5MevBin4;                         //! Invariant mass D01+D02 (good hyp) 5Mev
  TH1F *fhInvMassD0barSum5MevBin4;                      //! Invariant mass D0bar1+D0bar2 (good hyp) 5Mev
  TH1F *fhInvMassSumAll5MevBin4;                        //! Invariant mass superimpose (good hyp only) 5Mev

  TH1F *fhInvMassD0Sum10MevBin5;                        //! Invariant mass D01+D02 (good hyp) 10Mev BIN5
  TH1F *fhInvMassD0barSum10MevBin5;                     //! Invariant mass D0bar1+D0bar2 (good hyp) 10Mev
  TH1F *fhInvMassSumAll10MevBin5;                       //! Invariant mass superimpose (good hyp only) 10Mev
  TH1F *fhInvMassD0Sum5MevBin5;                         //! Invariant mass D01+D02 (good hyp) 5Mev
  TH1F *fhInvMassD0barSum5MevBin5;                      //! Invariant mass D0bar1+D0bar2 (good hyp) 5Mev
  TH1F *fhInvMassSumAll5MevBin5;                        //! Invariant mass superimpose (good hyp only) 5Mev

  TH2F *fhReflBin1;			 		//! Analysis for reflection - Bin1
  TH2F *fhReflBin2;			  		//! Analysis for reflection - Bin2
  TH2F *fhReflBin3;			  		//! Analysis for reflection - Bin3
  TH2F *fhReflBin4;			  		//! Analysis for reflection - Bin4
  TH2F *fhReflBin5;			  		//! Analysis for reflection - Bin5
  TH2F *fhReflD0Bin1;			 		//! Analysis for reflection - D0 study - Bin1
  TH2F *fhReflD0Bin2;			  		//! Analysis for reflection - D0 study - Bin2
  TH2F *fhReflD0Bin3;			  		//! Analysis for reflection - D0 study - Bin3
  TH2F *fhReflD0Bin4;			  		//! Analysis for reflection - D0 study - Bin4
  TH2F *fhReflD0Bin5;			  		//! Analysis for reflection - D0 study - Bin5
  TH2F *fhReflD0barBin1;				//! Analysis for reflection - D0bar study - Bin1
  TH2F *fhReflD0barBin2;				//! Analysis for reflection - D0bar study - Bin2
  TH2F *fhReflD0barBin3;				//! Analysis for reflection - D0bar study - Bin3
  TH2F *fhReflD0barBin4;				//! Analysis for reflection - D0bar study - Bin4
  TH2F *fhReflD0barBin5;				//! Analysis for reflection - D0bar study - Bin5

  TH2F *fScatterP4PID;					  //! K momentum vs like sign Pi momentum after PID
  TH2F *fPtVsY;						  //! Pt vs Y of selected candidates (by PPR cuts)
  TH2F *fPtVsYAll;					  //! Pt vs Y of all candidates

  TH1F *fEventCounter;				       //! Event Counter
  TH1F *fCutDCA;				       //! DCA histogram doubl.
  TH1F *fCutDCA3;				       //! DCA histogram trips
  TH1F *fCutDCA2;				       //! DCA histogram quads1
  TH1F *fCutDCA5;				       //! DCA histogram quads2
  TH1F *fCutVertexDist2;			       //! Vertex doubl. to primary distance
  TH1F *fCutVertexDist3;			       //! Vertex trips to primary distance
  TH1F *fCutVertexDist4;			       //! Vertex quads to primary distance
  TH1F *fCutCosinePoint;			       //! Cosine of pointing angle
  TH1F *fCutPt;				       	       //! Candidate D0 Pt
  TH1F *fCutY;				       	       //! Candidate D0 Y
  TH1F *fPIDSel;				       //! PID Selected
  TH1F *fPIDSelBin1;			               //! PID Selected Bin1
  TH1F *fPIDSelBin2;			               //! PID Selected Bin2
  TH1F *fPIDSelBin3;			               //! PID Selected Bin3
  TH1F *fPIDSelBin4;			               //! PID Selected Bin4
  TH1F *fPIDSelBin5;			               //! PID Selected Bin5
  TH2F *fMultipleHyps;				       //! Multiple hypotesis accepted counter
  TH2F *fMultipleHypsType;			       //! Multiple hypotesis accepted counter

  TH1F *fPtSel;				       	       //! Pt of selected candidates
 
  AliRDHFCutsD0toKpipipi *fCuts;			// Cuts container

  ClassDef(AliAnalysisTaskSESelectHF4Prong,2); // AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
};

#endif

