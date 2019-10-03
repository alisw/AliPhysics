#ifndef ALIANALYSISTASKSEHFCJQA_H
#define ALIANALYSISTASKSEHFCJQA_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEHFCJqa
// AliAnalysisTask with QA plots for HFCJ analyses
//
// Authors: Andrea Rossi, andrea.rossi@cern.ch
//          Elena Bruna,  elena.bruna@to.infn.it
//*************************************************************************

class TH1F;
class TH2F;
class TH3F;
class AliAODDEvent;
class AliAODMCHeader;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF;
class AliAODMCParticle;
class AliAnalysisVertexingHF;
class AliRDHFCutsD0toKpi;
class AliNormalizationCounter;
class AliPIDResponse;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEHFCJqa : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSEHFCJqa();
  AliAnalysisTaskSEHFCJqa(const char* name);
  virtual ~AliAnalysisTaskSEHFCJqa();

 // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);  
  AliAODMCParticle* IsMCJet(TClonesArray *arrayMC,const AliAODJet *jet, Double_t &contribution);
  AliAODMCParticle* GetMCPartonOrigin(TClonesArray *arrayMC,AliAODMCParticle *p, Int_t &idx);
  void SetCutObject(AliRDHFCuts *cuts){fCuts=cuts;}
  void SetFilterBit(Int_t bit){ffilterbit=bit;}
  void SetLoadJet(Int_t ljet,TString strJetArray=""){fLoadJet=ljet;fJetArrayString=strJetArray;}
 
 private:
  AliAnalysisTaskSEHFCJqa(const AliAnalysisTaskSEHFCJqa&); // copy constructo not implemented yet
  AliAnalysisTaskSEHFCJqa& operator=(const AliAnalysisTaskSEHFCJqa&); // assignment operator not implemented yet
  Bool_t FillTrackHistosAndSelectTrack(AliAODTrack *aodtr,const AliESDVertex *primary,Double_t magfield); // method to filter tracks and fill related histograms
  void SetupPIDresponse();
  void FillJetRecoHisto(const AliAODJet *jet,Int_t partonnat,Double_t contribution,Double_t ptpart);
  void FillTrackHistosPID(AliAODTrack *aodtr);

  Bool_t fReadMC;                         // flag to read MC data
  Int_t ffilterbit;                       // selected filter bit
  Bool_t fKeepTrackNegID;                //  flag for rejecting track with neg ID
  AliPIDResponse *fpidResp;               // !pid response object
  AliRDHFCuts *fCuts;                     // cut object (temporary D2H cut object calss used)
  TH1F *fhEventCounter;                   //! histo with counter of event selected
  TH3F *fhImpParResolITSsel;              //! histo with imp par distribution as a function of pt and ITS clusters
  TH3F *fhImpParResolITSselGoodTracks;    //! histo with imp par distribution as a function of pt and ITS clusters for selected tracks
  THnSparseF *fhSparseFilterMask;          //! sparse histo with track information
  THnSparseF *fhSparseFilterMaskTrackAcc;    //! sparse with filter bits and track kine/geometrical properties
  THnSparseF *fhSparseFilterMaskImpPar;    //! sparse with kine/geometrical prop and imp par xy
  THnSparseF *fhSparseEoverPeleTPC;       //! sparse histo with TPC-EMCal PID electron information
  THnSparseF *fhSparseShowShapeEleTPC;    //! sparse histo with TPC-EMCAL PID electron info, including shower shape & Ncells
  TH3F *fhnSigmaTPCTOFEle;                //! sparse with TPC-TOF nsigma informations, centered on ele hypo
  TH3F *fhnSigmaTPCTOFPion;              //! sparse with TPC-TOF nsigma informations, centered on pion hypo
  TH3F *fhnSigmaTPCTOFKaon;               //! sparse with TPC-TOF nsigma informations, centered on kaon hypo
  TH3F *fhnSigmaTPCTOFProton;             //! sparse with TPC-TOF nsigma informations, centered on proton hypo
  THnSparseF *fhTrackEMCal;              //! sparse with EMCal cluster properties related to clusters matched to tracks 
  THnSparseF *fSparseRecoJets;           //! sparse histo with jet properties
  Int_t fLoadJet;                         // flag for reading jet array (0=no, 1=online reco jets, 2=from friend file)
  TString fJetArrayString;                // jet array name
  TList *fListTrackAndPID;                // list with single track and PID properties
  TList *fListJets;                       // list with jet properties
ClassDef(AliAnalysisTaskSEHFCJqa,2); // analysis task for MC study
};

#endif
