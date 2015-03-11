#ifndef ALIANALYSISTASKSEDMESONSFILTERCJ_H
#define ALIANALYSISTASKSEDMESONSFILTERCJ_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

//-----------------------------------------------------------------------
// Author : A. Grelli, Utrecht University
//          C. Bianchin, Utrecht University
//          X. Zhang, LBNL
//          S. Aiola, Yale University
//-----------------------------------------------------------------------


#include "AliAnalysisTaskSE.h"

class TH2;
class TString;
class TClonesArray;
class AliRDHFCuts;
class AliAODRecoCascadeHF;

class AliAnalysisTaskSEDmesonsFilterCJ : public AliAnalysisTaskSE 
{

 public:

  enum ECandidateType{ kD0toKpi, kDstartoKpipi };
  
  AliAnalysisTaskSEDmesonsFilterCJ();
  AliAnalysisTaskSEDmesonsFilterCJ(const Char_t* name,AliRDHFCuts* cuts,ECandidateType candtype);
  virtual ~AliAnalysisTaskSEDmesonsFilterCJ();

  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  virtual void     Init();
  virtual void     LocalInit() { Init(); }

  // inizializations
  Bool_t DefineHistoForAnalysis();

  // set MC usage
  void   SetMC(Bool_t theMCon) { fUseMCInfo = theMCon ; }
  Bool_t GetMC() const         { return fUseMCInfo    ; }
  
  // set usage of generated or reconstucted quantities (relevant for MC)
  void SetUseReco(Bool_t useReco=kTRUE) { fUseReco = useReco ; }
  Bool_t GetUseReco() const             { return fUseReco    ; }
 
  void SetMassLimits(Double_t range, Int_t pdg);
  void SetMassLimits(Double_t lowlimit, Double_t uplimit);

  // Array of D0 width for the Dstar
  Bool_t SetD0WidthForDStar(Int_t nptbins, Float_t *width);
  
  Float_t DeltaR(AliVParticle *p1, AliVParticle *p2) const;

 protected:
  void ExecOnce();
  void ProcessD0(AliAODRecoDecayHF* charmCand, Int_t isSelected, Int_t mcLabel);
  void ProcessDstar(AliAODRecoCascadeHF* dstar, Int_t mcLabel);
  
 private:
  
  AliAnalysisTaskSEDmesonsFilterCJ(const AliAnalysisTaskSEDmesonsFilterCJ &source);
  AliAnalysisTaskSEDmesonsFilterCJ& operator=(const AliAnalysisTaskSEDmesonsFilterCJ& source); 

  Bool_t          fUseMCInfo;              //  Use MC info
  Bool_t          fUseReco;                // use reconstructed tracks when running on MC
  UInt_t          fCandidateType;          //  Dstar or D0
  TString         fCandidateName;          //  Dstar or D0
  Int_t           fPDGmother;              //  PDG code of D meson
  Int_t           fNProngs;                //  number of prong of the decay channel  
  Int_t           fPDGdaughters[4];        //  PDG codes of daughters
  Float_t         fSigmaD0[30];            //  D0 sigma for Dstar
  TString         fBranchName;             //  AOD branch name
  AliRDHFCuts    *fCuts;                   //  cuts 
  Double_t        fMinMass;                //  mass lower limit histogram
  Double_t        fMaxMass;                //  mass upper limit histogram
  Bool_t          fInhibitTask;            //
  Bool_t          fInitOk;                 //!
  AliAODEvent    *fAodEvent;               //!
  TClonesArray   *fArrayDStartoD0pi;       //!
  TClonesArray   *fMCarray;                //!
  TClonesArray   *fCandidateArray;         //! contains candidates selected by AliRDHFCuts
  TClonesArray   *fSideBandArray;          //! contains candidates selected by AliRDHFCuts::IsSelected(kTracks), to be used for side bands (DStar case only!!)
  TList          *fOutput;                 //! user output
  TH1            *fHistStat;               //!
  TH1            *fHistNSBCandEv;          //!
  TH1            *fHistNCandEv;            //!
  TH2            *fHistImpParS;            //!
  TH2            *fHistImpParB;            //!
  TH1            *fHistPtPion;             //!
  TH2            *fHistInvMassPtD;         //!
  TH1            *fHistInvMassS;           //!
  TH1            *fHistInvMassB;           //!
  TH2            *fHistAlphaDDS;           //!
  TH2            *fHistAlphaDpisS;         //!
  TH2            *fHistAlphaDpiS;          //!
  TH2            *fHistAlphaDKS;           //!
  TH2            *fHistAlphaDDB;           //!
  TH2            *fHistAlphaDpisB;         //!
  TH2            *fHistAlphaDpiB;          //!
  TH2            *fHistAlphaDKB;           //!
  TH2            *fHistDeltaRDDS;          //!
  TH2            *fHistDeltaRDpisS;        //!
  TH2            *fHistDeltaRDpiS;         //!
  TH2            *fHistDeltaRDKS;          //!
  TH2            *fHistDeltaRDDB;          //!
  TH2            *fHistDeltaRDpisB;        //!
  TH2            *fHistDeltaRDpiB;         //!
  TH2            *fHistDeltaRDKB;          //!
  TH2            *fHistAlphaDpiR;          //!
  TH2            *fHistAlphaDKR;           //!
  TH2            *fHistDeltaRDpiR;         //!
  TH2            *fHistDeltaRDKR;          //!

  ClassDef(AliAnalysisTaskSEDmesonsFilterCJ, 4); // task for selecting D mesons to be used as an input for D-Jet correlations
};

#endif
