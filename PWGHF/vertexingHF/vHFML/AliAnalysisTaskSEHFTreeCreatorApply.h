#ifndef AliAnalysisTaskSEHFTreeCreatorApply_H
#define AliAnalysisTaskSEHFTreeCreatorApply_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///*************************************************************************
/// \class AliAnalysisTaskSEHFTreeCreatorApply
///
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
///*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include "TProfile.h"

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliNormalizationCounter.h"
#include "AliPIDResponse.h"
#include "AliHFTreeHandlerApply.h"
#include "AliHFTreeHandlerApplyLc2V0bachelor.h"
#include "AliHFTreeHandlerApplyDstoKKpi.h"
#include "AliHFMLResponse.h"

class AliAODMCHeader;
class AliAODEvent;
class TClonesArray;
class AliCDBEntry;

class AliAnalysisTaskSEHFTreeCreatorApply : public AliAnalysisTaskSE
{
public:
  
  AliAnalysisTaskSEHFTreeCreatorApply();
  AliAnalysisTaskSEHFTreeCreatorApply(const char *name, TList *cutsList);
  virtual ~AliAnalysisTaskSEHFTreeCreatorApply();
  
  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetReadMC(Bool_t opt=kFALSE){fReadMC=opt;}
  Bool_t GetReadMC() const {return fReadMC;}
  void SetUseSelectionBit(Bool_t opt=kFALSE){fUseSelectionBit=opt;}
  Bool_t GetUseSelectionBit() const {return fUseSelectionBit;}
  void SetSystem(Int_t opt){fSys=opt;}
  Int_t GetSystem() const {return fSys;}
  void SetAODMismatchProtection(Int_t opt=0) {fAODProtection=opt;}
  Int_t GetAODMismatchProtection() const {return fAODProtection;}

  void SetFillDsTree(Int_t opt){fWriteVariableTreeDs=opt;}
  Int_t GetFillDsTree() const {return fWriteVariableTreeDs;}
  void SetFillLc2V0bachelorTree(Int_t opt){fWriteVariableTreeLc2V0bachelor=opt;}
  Int_t GetFillLc2V0bachelorTree() const {return fWriteVariableTreeLc2V0bachelor;}
  void SetPIDoptDsTree(Int_t opt){fPIDoptDs=opt;}
  Int_t GetPIDoptDsTree() const {return fPIDoptDs;}
  void SetPIDoptLc2V0bachelorTree(Int_t opt){fPIDoptLc2V0bachelor=opt;}
  Int_t GetPIDoptLc2V0bachelorTree() const {return fPIDoptLc2V0bachelor;}
  void SetWriteOnlySignalTree(Bool_t opt){fWriteOnlySignal=opt;}
  Bool_t GetWriteOnlySignalTree() const {return fWriteOnlySignal;}
  void SetFillMCGenTrees(Bool_t fillMCgen) {fFillMCGenTrees=fillMCgen;}
  Bool_t GetFillMCGenTrees() const {return fFillMCGenTrees;}

  std::string GetPeriod(const AliVEvent *ev);
  void SetMultiplVsZProfile(std::string period, TProfile *hprof)
  {
    delete fMultEstimatorAvg[period];
    fMultEstimatorAvg[period] = new TProfile(*hprof);
  }
  void SetMultiplVsZProfileSHM(std::string period, TProfile *hprof)
  {
    delete fMultEstimatorAvgSHM[period];
    fMultEstimatorAvgSHM[period] = new TProfile(*hprof);
  }

  void SetRefMult(Double_t refMult) { fRefMult = refMult; }
  Double_t GetRefMult() { return fRefMult; }
  void SetRefMultSHM(Double_t refMult) { fRefMultSHM = refMult; }
  Double_t GetRefMultSHM() { return fRefMultSHM; }
  void SetCorrNtrVtx(Bool_t corr = true) { fCorrNtrVtx = corr; }
  Bool_t GetCorrNtrVtx() const { return fCorrNtrVtx; }
  void SetCorrV0MVtx(Bool_t corr = true) { fCorrV0MVtx = corr; }
  Bool_t GetCorrV0MVtx() const { return fCorrV0MVtx; }
  
  void SetDsMassKKOption(AliHFTreeHandlerApplyDstoKKpi::massKKopt opt) {fDsMassKKOpt=opt;}
  Int_t GetDsMassKKOption() const {return fDsMassKKOpt;}
  void SetLc2V0bachelorCalcSecoVtx(Int_t opt=1) {fLc2V0bachelorCalcSecoVtx=opt;}
  Int_t GetLc2V0bachelorCalcSecoVtx() const {return fLc2V0bachelorCalcSecoVtx;}
  void SetLc2V0type(Int_t opt=1) {fV0typeForLc2V0bachelor=opt;}
  Int_t GetLc2V0type() const {return fV0typeForLc2V0bachelor;}
  void SetTreeSingleTrackVarsOpt(Int_t opt) {fTreeSingleTrackVarsOpt=opt;}
  Int_t GetTreeSingleTrackVarsOpt() const {return fTreeSingleTrackVarsOpt;}
  void SetReducePbPbBranches(Bool_t b) { fReducePbPbBranches = b; }
  Bool_t GetReducePbPbBranches() const { return fReducePbPbBranches; }
  void SetSaveSTDSelection(Bool_t b) { fSaveSTDSelection = b; }
  Bool_t GetSaveSTDSelection() const { return fSaveSTDSelection; }

  void SetGoodTrackFilterBit(Int_t i) { fGoodTrackFilterBit = i; }
  Int_t GetGoodTrackFilterBit() const { return fGoodTrackFilterBit; }
  void SetGoodTrackEtaRange(Double_t d) { fGoodTrackEtaRange = d; }
  Double_t GetGoodTrackEtaRange() const { return fGoodTrackEtaRange; }
  void SetGoodTrackMinPt(Double_t d) { fGoodTrackMinPt = d; }
  Double_t GetGoodTrackMinPt() const { return fGoodTrackMinPt; }
  void SetITSUpgradeProduction(Bool_t b) { fITSUpgradeProduction = b; }
  Bool_t GetITSUpgradeProduction() const { return fITSUpgradeProduction; }
  void SetITSUpgradePreSelect(Bool_t b) { fITSUpgradePreSelect = b; }
  Bool_t GetITSUpgradePreselect() const { return fITSUpgradePreSelect; }

  void ApplyPhysicsSelectionOnline(bool apply=true) { fApplyPhysicsSelOnline = apply; }
  Bool_t GetApplyPhysicsSelOnline() const { return fApplyPhysicsSelOnline; }
  void ApplyEventSelectionOnline(bool apply=true) { fApplyEventSelOnline = apply; }
  Bool_t GetApplyEventSelOnline() const { return fApplyEventSelOnline; }
  void EnableEventDownsampling(float fractokeep, unsigned long seed) {
    fEnableEventDownsampling = true;
    fFracToKeepEventDownsampling = fractokeep;
    fSeedEventDownsampling = seed;
  }
  Bool_t GetEnableEventDownsampling() const { return fEnableEventDownsampling; }
  Float_t GetFracToKeepEventDownsampling() const { return fFracToKeepEventDownsampling; }
  unsigned long GetSeedEventDownsampling() const { return fSeedEventDownsampling; }
  void EnableCandDownsampling(float fractokeep, float maxptsampling) {
    fEnableCandDownsampling = true;
    fFracToKeepCandDownsampling = fractokeep;
    fMaxPtCandDownsampling = maxptsampling;
  }
  Bool_t GetEnableCandDownsampling() const { return fEnableCandDownsampling; }
  Float_t GetFracToKeepCandDownsampling() const { return fFracToKeepCandDownsampling; }
  Float_t GetMaxPtCandDownsampling() const { return fMaxPtCandDownsampling; }

  void SetMLConfigFile(TString path = ""){fConfigPath = path;}
  TString GetMLConfigFile() const { return fConfigPath; }

  void SetNsigmaTPCDataDrivenCorrection(Int_t syst) {
    fEnableNsigmaTPCDataCorr=true;
    fSystemForNsigmaTPCDataCorr=syst;
  }
  Bool_t GetEnableNsigmaTPCDataCorr() const { return fEnableNsigmaTPCDataCorr; }
  Int_t GetSystemForNsigmaTPCDataCorr() const { return fSystemForNsigmaTPCDataCorr; }
  
  void Process3Prong(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
  void ProcessCasc(TClonesArray *arrayCasc, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
  void ProcessMCGen(TClonesArray *mcarray, AliAODMCHeader *mcHeader);

  Bool_t CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau, Bool_t ITSUpgradeStudy);
  Bool_t IsCandidateFromHijing(AliAODRecoDecayHF *cand, AliAODMCHeader *mcHeader, TClonesArray* arrMC, AliAODTrack *tr = 0x0);
  void SelectGoodTrackForReconstruction(AliAODEvent *aod, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags);
  AliAODVertex* ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion);
  unsigned int GetEvID();

private:

  AliAnalysisTaskSEHFTreeCreatorApply(const AliAnalysisTaskSEHFTreeCreatorApply&);
  AliAnalysisTaskSEHFTreeCreatorApply& operator=(const AliAnalysisTaskSEHFTreeCreatorApply&);

  unsigned int            fEventNumber;
  TH1F                    *fNentries;                            //!<!   histogram with number of events on output slot 1
  TH2F                    *fHistoNormCounter;                    //!<!   histogram with number of events on output slot 1
  TList                   *fListCounter;                         //!<!   list for normalization counter on output slot 3
  AliNormalizationCounter *fCounter;                             //!<!   AliNormalizationCounter

  TList                   *fListCuts;                            /// list of cuts sent to output slot 2
  AliRDHFCutsDstoKKpi     *fFiltCutsDstoKKpi;                    /// DstoKKpi filtering (or loose) cuts
  AliRDHFCutsLctoV0       *fFiltCutsLc2V0bachelor;               /// Lc2V0bachelor filtering (or loose) cuts
  AliRDHFCutsDstoKKpi     *fCutsDstoKKpi;                        /// DstoKKpi analysis cuts
  AliRDHFCutsLctoV0       *fCutsLc2V0bachelor;                   /// Lc2V0bachelor analysis cuts
  AliRDHFCuts             *fEvSelectionCuts;                     /// Event selection cuts

  Bool_t                  fReadMC;                               /// flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t                  fUseSelectionBit;                      /// flag to use selection bit when looping over TClonesArray
  Int_t                   fSys;                                  /// fSys=0 -> p-p; fSys=1 ->PbPb
  Int_t                   fAODProtection;                        /// flag to activate protection against AOD-dAOD mismatch.
  // -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

  Bool_t                  fWriteOnlySignal;                       /// flag to decide wether to write only signal or signal+background
  Bool_t                  fFillMCGenTrees;                        /// flag to enable fill of the generated trees
  Int_t                   fWriteVariableTreeDs;                   /// flag to decide whether to write the Ds variables on a tree variables
  Int_t                   fWriteVariableTreeLc2V0bachelor;        /// flag to decide whether to write the Lc->pK0s variables on a tree variables
  // 0 don't fill
  // 1 fill standard tree

  TTree                   *fVariablesTreeDs;                     //!<! tree of the candidate variables
  TTree                   *fVariablesTreeLc2V0bachelor;          //!<! tree of the candidate variables
  TTree                   *fGenTreeDs;                           //!<! tree of the gen Ds variables
  TTree                   *fGenTreeLc2V0bachelor;                //!<! tree of the gen Lc2V0bachelor variables
  TTree                   *fTreeEvChar;                          //!<! tree of event variables

  AliHFTreeHandlerApplyDstoKKpi       *fTreeHandlerDs;                //!<! handler object for the tree with topological variables
  AliHFTreeHandlerApplyLc2V0bachelor  *fTreeHandlerLc2V0bachelor;     //!<! handler object for the tree with topological variables
  AliHFTreeHandlerApplyDstoKKpi       *fTreeHandlerGenDs;             //!<! handler object for the tree with topological variables
  AliHFTreeHandlerApplyLc2V0bachelor  *fTreeHandlerGenLc2V0bachelor;  //!<! handler object for the tree with topological variables

  AliPIDResponse          *fPIDresp;                             /// PID response
  Int_t                   fPIDoptDs;                             /// PID option for Ds tree
  Int_t                   fPIDoptLc2V0bachelor;                  /// PID option for Lc2V0bachelor tree

  UShort_t                fBC;                                   /// bunch crossing number
  Int_t                   fOrbit;                                /// orbit
  Int_t                   fPeriod;                               /// period
  Int_t                   fEventID;                              /// event ID (unique when combined with run number)
  Int_t                   fEventIDExt;                           /// upper 32-bit of event ID
  Long64_t                fEventIDLong;                          /// single unique event id (long64)
  TString                 fFileName;                             /// Store filename for an unique event ID
  unsigned int            fDirNumber;                            /// Store directory number for an unique event ID

  Float_t                 fCentrality;                           /// event centrality
  Float_t                 fzVtxReco;                             /// reconstructed Zvtx
  Float_t                 fzVtxGen;                              /// generated Zvtx
  Int_t                   fNcontributors;                        /// number of contributors
  Int_t                   fNtracks;                              /// number of tracks
  Int_t                   fIsEvRej;                              /// flag with information about rejection of the event
  Int_t                   fIsEvRej_INT7;                         /// flag with information about rejection of the event
  Int_t                   fIsEvRej_HighMultSPD;                  /// flag with information about rejection of the event
  Int_t                   fIsEvRej_HighMultV0;                   /// flag with information about rejection of the event
  Int_t                   fRunNumber;                            /// run number
  Int_t                   fRunNumberCDB;                         /// run number (for OCDB)
  Int_t                   fnTracklets;                           /// number of tracklets
  Int_t                   fnTrackletsCorr;                       /// number of tracklets (corrected)
  Int_t                   fnTrackletsCorrSHM;                    /// number of tracklets (corrected)
  Int_t                   fnV0A;                                 /// V0A multiplicity
  Int_t                   fnTPCCls;                              /// TPC multiplicity
  Int_t                   fMultGen;                              /// generated multiplicity around mid-rapidity [-1,1]
  Int_t                   fMultGenV0A;                           /// generated multiplicity in V0A range
  Int_t                   fMultGenV0C;                           /// generated multiplicity in V0C range
  ULong64_t               fTriggerMask;                          /// Trigger mask bitmap
  Bool_t                  fTriggerOnlineINT7;                    /// Flag explicitly whether bitmap contains INT7
  Bool_t                  fTriggerOnlineHighMultSPD;             /// Flag explicitly whether bitmap contains HighMultSPD
  Bool_t                  fTriggerOnlineHighMultV0;              /// Flag explicitly whether bitmap kHighMultV0
  Bool_t                  fTriggerBitINT7;                       /// Flag explicitly whether bitmap contains INT7
  Bool_t                  fTriggerBitHighMultSPD;                /// Flag explicitly whether bitmap contains HighMultSPD
  Bool_t                  fTriggerBitHighMultV0;                 /// Flag explicitly whether bitmap kHighMultV0
  Bool_t                  fTriggerBitCentral;                    /// Flag explicitly whether bitmap contains kCentral
  Bool_t                  fTriggerBitSemiCentral;                /// Flag explicitly whether bitmap contains kSemiCentral
  TString                 fTriggerClasses;                       /// Collect all trigger classes
  Bool_t                  fTriggerClassINT7;                     /// Flag explicitly whether classes contain INT7
  Bool_t                  fTriggerClassHighMultSPD;              /// Flag explicitly whether classes contain HighMultSPD
  Bool_t                  fTriggerClassHighMultV0m;              /// Flag explicitly whether classes contain HighMultV0
  Int_t                   fnV0M;                                 /// V0M multiplicity
  Int_t                   fnV0MEq;                               /// V0M multiplicity (equalized)
  Int_t                   fnV0MCorr;                             /// V0M multiplicity (corrected)
  Int_t                   fnV0MEqCorr;                           /// V0M multiplicity (equalized + corrected)
  Float_t                 fPercV0M;                              /// V0M multiplicity percentile
  Float_t                 fMultV0M;                              /// V0M multiplicity from mult selection task

  Double_t                fRefMult;                              /// reference multiplicity
  Double_t                fRefMultSHM;                           /// reference multiplicity
  std::map<std::string, TProfile*> fMultEstimatorAvg;
  std::map<std::string, TProfile*> fMultEstimatorAvgSHM;
  Bool_t fCorrNtrVtx;
  Bool_t fCorrV0MVtx;
  AliCDBEntry *fCdbEntry;

  Int_t                   fDsMassKKOpt;                          /// option for Ds massKK (mass or delta mass)
  Int_t                   fLc2V0bachelorCalcSecoVtx;             /// option to calculate the secondary vertex for Lc2V0bachelor. False by default, has to be added to AddTask in case we want to start using it.
  Int_t                   fV0typeForLc2V0bachelor;               /// option to select Offline+OnTheFly (0), only Offline (1=default), only OnTheFly (2) V0's for the Lc->V0bachelor decay
  Int_t                   fTreeSingleTrackVarsOpt;               /// option for single-track variables to be filled in the trees
  
  Int_t                   fGoodTrackFilterBit;                   /// Setting filter bit for bachelor on-the-fly reconstruction candidate
  Double_t                fGoodTrackEtaRange;                    /// Setting eta-range for bachelor on-the-fly reconstruction candidate
  Double_t                fGoodTrackMinPt;                       /// Setting min pT for bachelor on-the-fly reconstruction candidate
  Bool_t                  fITSUpgradeProduction;                 /// Setting for analysing an ITS Upgrade production
  Bool_t                  fITSUpgradePreSelect;                  /// Setting to enable ITSUpgrade Preselect function

  Bool_t fEnableNsigmaTPCDataCorr;                               /// flag to enable data-driven NsigmaTPC correction
  Int_t fSystemForNsigmaTPCDataCorr;                             /// system for data-driven NsigmaTPC correction

  Bool_t fApplyPhysicsSelOnline;                                 /// flag to apply physics selection in the task
  Bool_t fApplyEventSelOnline;                                   /// flag to fill reco TTrees only with selected events
  Bool_t fEnableEventDownsampling;                               /// flag to apply event downsampling
  Float_t fFracToKeepEventDownsampling;                          /// fraction of events to be kept by event downsampling
  unsigned long fSeedEventDownsampling;                          /// seed for event downsampling
  Bool_t fEnableCandDownsampling;                                /// flag to apply cand downsampling
  Float_t fFracToKeepCandDownsampling;                           /// fraction of cands to be kept by cand downsampling
  Float_t fMaxPtCandDownsampling;                                /// max pT used for cand downsampling

  TString fConfigPath;                                           /// path to ML config file
  AliHFMLResponse* fMLResponse;                                  //!<! object to handle ML response
  Bool_t fReducePbPbBranches;                                    /// variable to disable unnecessary branches in PbPb
  Bool_t fSaveSTDSelection;                                      /// variable to store candidates that pass std cuts as well, even when ML < MLcut

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEHFTreeCreatorApply,7);
  /// \endcond
};

#endif
