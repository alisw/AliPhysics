/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>             *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// The task:
// stores TRD PID quantities in a Tree
//
//  Author:
//  Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//

#ifndef ALITRDPIDTREE_H
#define ALITRDPIDTREE_H
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

class TTreeStream;
class AliInputEventHandler;
class TArrayF;
template <class X>
class THnSparseT;
typedef class THnSparseT<TArrayF> THnSparseF;
class TFile;
class AliESDEvent;
class AliMCEvent;
class AliESDtrackCuts;
class AliPIDResponse;
class AliESD;
class AliESDtrack;
class AliAnalysisTask;
class AliESDInputHandler;
class AliESDv0KineCuts;
class AliAnalysisManager;
//class AliCentrality;
class AliMultSelection;
class AliTRDgeometry;
class TTree;
class TSystem;
class TStyle;
class TROOT;
class Riostream;
class TChain;
class TH2;
class TF1;
class TH1;
class TObjArray;


class AliTRDPIDTree : public AliAnalysisTaskSE {

    public:
  typedef enum{
      kpp = 0,
      kpPb = 1,
      kPbPb = 2
  } ECollisionSystem_t;

  typedef enum{
      kNoPileUpCut=0,
      kLHC15o = 1,
      kLHC18q = 2
  } Period;

  AliTRDPIDTree(const char *name = "trd_pid_tree");
  virtual ~AliTRDPIDTree();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0);
  virtual void   Terminate(const Option_t*);
  Int_t          CompareFloat(Float_t f1=1, Float_t f2=0) const;
  Int_t          GetV0tag(Int_t trackIndex) const;
  Int_t          *fV0tags;  //! Pointer to array with tags for identified particles from V0 decays

  Bool_t Ispp() const { return fCollisionSystem.TestBitNumber(kpp); }
  Bool_t IspPb() const { return fCollisionSystem.TestBitNumber(kpPb); }
  Bool_t IsPbPb() const { return fCollisionSystem.TestBitNumber(kPbPb); }

  void SetCollisionSystem(ECollisionSystem_t system){
      fCollisionSystem.Clear();
      fCollisionSystem.SetBitNumber(system, kTRUE);
  }
  void SetppAnalysis(){
      fCollisionSystem.SetBitNumber(kpPb, kFALSE);
      fCollisionSystem.SetBitNumber(kPbPb, kFALSE);
      fCollisionSystem.SetBitNumber(kpp, kTRUE);
  }
  void SetpPbAnalysis() {
      fCollisionSystem.SetBitNumber(kpp, kFALSE);
      fCollisionSystem.SetBitNumber(kPbPb, kFALSE);
      fCollisionSystem.SetBitNumber(kpPb, kTRUE);
  }
  void SetPbPbAnalysis() {
      fCollisionSystem.SetBitNumber(kpp, kFALSE);
      fCollisionSystem.SetBitNumber(kpPb, kFALSE);
      fCollisionSystem.SetBitNumber(kPbPb, kTRUE);
  };

  void SetUseExtraPileupCut(Int_t UseExtraPileupCut=0){
      fUseExtraPileupCut=UseExtraPileupCut;
  };

  protected:
  static Double_t fgMinLayer;  //! Cut variable for min number of layers
  AliESDv0KineCuts *fV0cuts;           //! ESD V0 cuts
  TObjArray *fV0electrons;             //! array with pointer to identified particles from V0 decays (electrons)
  TObjArray *fV0pions;                 //! array with pointer to identified particles from V0 decays (pions)
  TObjArray *fV0protons;               //! array with pointer to identified particles from V0 decays (ptotons)

  void FillTree(AliESDtrack *track, Int_t pdg, Int_t runnumber, Int_t centralityvalue);
  void SetupV0qa();
  void FillV0PIDlist();
  void ClearV0PIDlist();
  Double_t GetPhi(AliESDtrack *const fTrack,Int_t iPl, Double_t& ar, Double_t& thetaar);
  Bool_t PassTrackCuts(AliESDtrack *fESDTrack=0);
  Bool_t HasMissingLayer(const AliVTrack *fESDTrack=0);
  Int_t  GetNTrackletsPID(const AliVTrack *fESDTrack=0) const;

  private:
  //
  //
  AliESDEvent *fESDEvent;              //! ESD object
  AliMCEvent  *fMCEvent;               //! MC object
  AliStack    *fMCStack;               //! MC stack
  TTree *fTreeTRDPID;                  //! Tree with V0
  AliPIDResponse *fPIDResponse;        //! PID handling
  TObjArray *fOutputContainer;         //! output data container
  AliESDtrackCuts *fESDtrackCuts;      //! basic cut variables for all non-V0 tracks
  AliESDtrackCuts *fESDtrackCutsV0;    //! basic cut variables for all V0 tracks
  //
  TList                 *fListQATRD;   //! List with QATRD histograms
  TList                 *fListQATRDV0; //! List with V0 kine cuts QATRD histograms

  Int_t fNumTagsStored;                //! Number of entries of fV0tags

  TBits fCollisionSystem;              //! Collision System;

  Float_t fpdg;                        //! particle type (pdg value)
  Int_t frun;                          //! run number

  Int_t fUseExtraPileupCut;           // cut on correlation of VZERO multiplicity & TPCout tracks 
  
  // TTree stuff for PID References
  Int_t frunnumber;                  //! Tree: Run number
  Double_t fcentrality;              //! Tree: Centrality
  Double_t fTRDslices[48];           //! Tree: dEdx slices (8 per layer)
  Double_t fTRDMomentum[6];          //! Tree: local track momentum at anode wire
  Int_t fTRDNcls;                    //! Tree: number of clusters
  Int_t fTRDntracklets;              //! Tree: number of tracking tracklets
  Int_t fTRDntrackletsPID;           //! Tree: number of pid tracklets
  Double_t fTRDphi[6];               //! Tree: local track inclination phi (from track extrapol)
  Double_t fTRDY[6];                 //! Tree: local track y-position (from track extrapol)
  Double_t fTRDtheta;                //! Tree: local track theta in layer 0 (from track extrapol)
  Double_t fTRDthetalayer[6];        //! Tree: local track theta (from track extrapol)
  Double_t fTRDTPCtgl;               //! Tree: TPC dip angle
  Double_t fTRDsignal;               //! Tree: Truncated mean signal
  Int_t fTRDnclsdEdx;                //! Tree: number of clusters dedx (truncated mean method)
  Int_t fTRDnch;                     //! Tree: number of chambers dedx (truncated mean method)
  Float_t fNSigmaTPC[3];             //! Tree: TPC nsigma ele, pion, proton
  Float_t fNSigmaTOF[3];             //! Tree: TOF nsigma ele, pion, proton
  Int_t fPDG;                        //! Tree: PDG value from V0 identification
  Int_t fPDGTRUE;                    //! Tree: true PDG value (only available in MC)
  Int_t fTrackCharge;                //! Tree: charge of track
  Float_t fDCA[2];                   //! Tree: DCA
  Float_t fChi2;                     //! Tree: Chi at TRD
  Double_t fsigmaTRD[5];             //! Tree: truncated mean sigma
  Double_t fdeltaTRD[5];             //! Tree: truncated delta
  Double_t fratioTRD[5];             //! Tree: truncated ratio

  // Histograms
  TH1F *fhtrackCuts;                 //! Track and Event Cuts - QA
  TH1F *fhEventCount;                //! count number of events analysed 
  TH2F *fhArmenteros;                //! 2D V0 QA Hist
  TH2F *fHistV0MvsTPCoutBeforePileUpCuts; //! histos to monitor pile up cuts
  TH2F *fHistV0MvsTPCoutAfterPileUpCuts;  //!

  AliTRDPIDTree(const AliTRDPIDTree&); // not implemented
  AliTRDPIDTree& operator=(const AliTRDPIDTree&); // not implemented
  
  ClassDef(AliTRDPIDTree, 5);
};
#endif
