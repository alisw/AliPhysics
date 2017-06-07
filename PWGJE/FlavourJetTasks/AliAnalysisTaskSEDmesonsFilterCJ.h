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
//	    B. Trzeciak, Utrecht University
//-----------------------------------------------------------------------


#include "AliAnalysisTaskEmcal.h"

class TH2;
class TString;
class TClonesArray;
class AliAODMCHeader;
class AliRDHFCuts;
class AliAODRecoCascadeHF;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecay;
class AliStack;

class AliAnalysisTaskSEDmesonsFilterCJ : public AliAnalysisTaskEmcal 
{

 public:

  enum ECandidateType{ kD0toKpi, kDstartoKpipi };
  enum EParticleOrigin { kQuarkNotFound, kFromCharm, kFromBottom };
  enum EDecayChannel { kDecayOther, kDecayD0toKpi, kDecayDStartoKpipi };
  
  AliAnalysisTaskSEDmesonsFilterCJ();
  AliAnalysisTaskSEDmesonsFilterCJ(const Char_t* name,AliRDHFCuts* cuts,ECandidateType candtype);
  virtual ~AliAnalysisTaskSEDmesonsFilterCJ();

  void     UserCreateOutputObjects();
  Bool_t   Run();
  void     Init();
  void     LocalInit() { Init(); }

  // inizializations
  Bool_t DefineHistoForAnalysis();

  // set MC usage
  void   SetMC(Bool_t theMCon) { fUseMCInfo = theMCon ; }
  Bool_t GetMC() const         { return fUseMCInfo    ; }

  // set MC RM or eff
  void   SetBuildRMEff(Bool_t theRM)   { fBuildRMEff = theRM ; }
  Bool_t GetBuildRMEff() const         { return fBuildRMEff    ; }
  
  void   SetUsePythia(Bool_t theUsePythia) { fUsePythia = theUsePythia; }
  Bool_t GetUsePythia() const { return fUsePythia; }

  // set usage of generated or reconstucted quantities (relevant for MC)
  void SetUseReco(Bool_t useReco=kTRUE)    { fUseReco = useReco           ; }
  Bool_t GetUseReco() const                { return fUseReco              ; }

  void   SetCombineDmesons(Bool_t c)       { fCombineDmesons = c          ; }
  Bool_t GetCombineDmesons() const         { return fCombineDmesons       ; }

  void   SetMultipleCandidates(Bool_t c)   { fMultCand = c                ; }
  Bool_t GetMultipleCandidates() const     { return fMultCand             ; }

  void   SetAnalysedCandidate(Int_t c)     { fAnalyseCand = c             ; }
  Int_t  GetAnalysedCandidate() const      { return fAnalyseCand          ; }

  void   SetRejectQuarkNotFound(Bool_t c)  { fRejectQuarkNotFound = c     ; }
  Bool_t GetRejectQuarkNotFound() const    { return fRejectQuarkNotFound  ; }

  void   SetRejectDfromB(Bool_t c)         { fRejectDfromB = c            ; }
  Bool_t GetRejectDfromB() const           { return fRejectDfromB         ; }

  void   SetKeepOnlyDfromB(Bool_t c)       { fKeepOnlyDfromB = c          ; }
  Bool_t GetKeepOnlyDfromB() const         { return fKeepOnlyDfromB       ; }
 
  void SetMassLimits(Double_t range, Int_t pdg);
  void SetMassLimits(Double_t lowlimit, Double_t uplimit);

  // Array of D0 width for the Dstar
  Bool_t SetD0WidthForDStar(Int_t nptbins, Float_t *width);
  
  Float_t DeltaR(AliVParticle *p1, AliVParticle *p2) const;

  static Double_t AddDaughters(AliAODRecoDecay* cand, TObjArray& daughters);
  static Double_t AddMCDaughters(AliAODMCParticle* mcDmeson, TObjArray& mcdaughters,  TClonesArray* mcArray);

  static Int_t CheckOrigin(AliAODRecoDecay* cand, TClonesArray* mcArray); // AOD
  static Int_t CheckOrigin(AliAODMCParticle* part, TClonesArray* mcArray); // AOD
  static Int_t CheckOrigin(AliAODRecoDecay* cand, AliStack* stack); // ESD
  static Int_t CheckOrigin(Int_t ipart, AliStack* stack); // ESD
  
  static Int_t CheckDecayChannel(AliAODMCParticle* part, TClonesArray* mcArray); // AOD
  static Int_t CheckDecayChannel(Int_t ipart, AliStack* stack); // ESD
  
  void GetTrackPrimaryGenerator(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen);
  void GetMCTrackPrimaryGenerator(AliAODMCParticle *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen);
  Bool_t IsTrackInjected(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC);
  Bool_t IsMCTrackInjected(AliAODMCParticle *track,AliAODMCHeader *header,TClonesArray *arrayMC);


 protected:
  void ExecOnce();
  void ProcessD0(AliAODRecoDecayHF2Prong* charmCand, Int_t isSelected);
  void ProcessDstar(AliAODRecoCascadeHF* dstar, Int_t isSelected);
  void FillD0MCTruthKinHistos(AliAODRecoDecayHF2Prong* charmCand, Int_t isSelected, Int_t isD0);
  void FillDStarMCTruthKinHistos(AliAODRecoCascadeHF* dstar, Int_t isSelected, Int_t isDstar);
  void FillDstarSideBands(AliAODRecoCascadeHF* dstar);
  void AddEventTracks(TClonesArray* coll, AliParticleContainer* tracks);
  void AddMCEventTracks(TClonesArray* coll, AliParticleContainer* mctracks);
  

  Bool_t          fUseMCInfo;              //  Use MC info
  Bool_t 	  fBuildRMEff;		   //  MC RM or efficiency studies
  Bool_t 	  fUsePythia;		   //  Use Pythia info only for MC
  Bool_t          fUseReco;                //  use reconstructed tracks when running on MC
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
  Bool_t          fCombineDmesons;         //  create an additional collection with D meson candidates and the rest of the tracks (for jet finding)
  Bool_t          fMultCand;               //  In case of multiple candidates per event
  Int_t           fAnalyseCand;            //  Number of the candidate to be analysed
  Bool_t          fRejectQuarkNotFound;    //  reject D mesons for which the original charm or bottom quark could not be found (MC)
  Bool_t          fRejectDfromB;           //  reject D mesons coming from a B meson decay (MC)
  Bool_t          fKeepOnlyDfromB;         //  only accept D mesons coming from a B meson decay (MC)
  AliAODEvent    *fAodEvent;               //!
  AliAODMCHeader *fMCHeader;		   //!
  TClonesArray   *fArrayDStartoD0pi;       //!
  TClonesArray   *fMCarray;                //!
  TClonesArray   *fCandidateArray;         //! contains candidates selected by AliRDHFCuts
  TClonesArray   *fSideBandArray;          //! contains candidates selected by AliRDHFCuts::IsSelected(kTracks), to be used for side bands (DStar case only!!)
  TClonesArray   *fCombinedDmesons;        //! contains candidates selected by AliRDHFCuts and the rest of the event tracks
  TClonesArray   *fCombinedDmesonsBkg;     //! contains bkg candidates selected by AliRDHFCuts and the rest of the event tracks
  TClonesArray   *fMCCombinedDmesons;      //! contains MC D0 and MC event particles
  Int_t           fNCand;                  //! number of selected D candidates already added to fCandidateArray
  Int_t           fNSBCand;                //! number of selected side-band D candidates already added to fSideBandArray
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


 private:
  
  AliAnalysisTaskSEDmesonsFilterCJ(const AliAnalysisTaskSEDmesonsFilterCJ &source);
  AliAnalysisTaskSEDmesonsFilterCJ& operator=(const AliAnalysisTaskSEDmesonsFilterCJ& source); 

  ClassDef(AliAnalysisTaskSEDmesonsFilterCJ, 7); // task for selecting D mesons to be used as an input for D-Jet correlations
};

#endif
