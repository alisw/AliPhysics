#ifndef AliAnalysisTaskSEHFTreeCreator_H
#define AliAnalysisTaskSEHFTreeCreator_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///*************************************************************************
/// \class AliAnalysisTaskSEHFTreeCreator
/// 
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// J. Mulligan, james.mulligan@berkeley.edu
// N. Zardoshti, nima.zardoshti@cern.ch
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
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliNormalizationCounter.h"
#include "AliPIDResponse.h"
#include "AliHFTreeHandler.h"
#include "AliHFTreeHandlerD0toKpi.h"
#include "AliHFTreeHandlerDplustoKpipi.h"
#include "AliHFTreeHandlerDstoKKpi.h"
#include "AliHFTreeHandlerLctopKpi.h"
#include "AliHFTreeHandlerBplustoD0pi.h"
#include "AliHFTreeHandlerBstoDspi.h"
#include "AliHFTreeHandlerDstartoKpipi.h"
#include "AliHFTreeHandlerLc2V0bachelor.h"
#include "AliHFTreeHandlerLbtoLcpi.h"
#include "AliHFTreeHandlerInclusiveJet.h"
#include "AliJetTreeHandler.h"
#include "AliParticleTreeHandler.h"
#include "AliTrackletTreeHandler.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliMCParticleContainer.h"
#include "AliJetContainer.h"

class AliAODMCHeader;
class AliAODEvent;
class TClonesArray;
class AliEmcalJet;
class AliRhoParameter;
class AliCDBEntry;

class AliAnalysisTaskSEHFTreeCreator : public AliAnalysisTaskSE
{
public:
  
  enum JetAlgorithm{
    antikt,
    kt,
    ca
  };
   
    AliAnalysisTaskSEHFTreeCreator();
    AliAnalysisTaskSEHFTreeCreator(const char *name,TList *cutsList, int fillNJetTrees, bool fillJetConstituentTrees);
    virtual ~AliAnalysisTaskSEHFTreeCreator();
    
    
    /// Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void ExecOnce();
    virtual Bool_t RetrieveEventObjects();
    virtual void Terminate(Option_t *option);
    
    void SetRefMult(Double_t refMult) { fRefMult = refMult; }
    Double_t GetRefMult() { return fRefMult; }
    void SetRefMultSHM(Double_t refMult) { fRefMultSHM = refMult; }
    Double_t GetRefMultSHM() { return fRefMultSHM; }
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
    std::string GetPeriod(const AliVEvent *ev);
    void SetCorrNtrVtx(bool corr = true) { fCorrNtrVtx = corr; }
    bool GetCorrNtrVtx() const { return fCorrNtrVtx; }
    void SetCorrV0MVtx(bool corr = true) { fCorrV0MVtx = corr; }
    bool GetCorrV0MVtx() const { return fCorrV0MVtx; }
    
    void SetReadMC(Bool_t opt=kFALSE){fReadMC=opt;}
    void SetSystem(Int_t opt){fSys=opt;}
    void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
    void SetWriteOnlySignalTree(Bool_t opt){fWriteOnlySignal=opt;}
    void SetFillD0Tree(Int_t opt){fWriteVariableTreeD0=opt;}
    void SetFillDsTree(Int_t opt){fWriteVariableTreeDs=opt;}
    void SetFillDplusTree(Int_t opt){fWriteVariableTreeDplus=opt;}
    void SetFillLctopKpiTree(Int_t opt){fWriteVariableTreeLctopKpi=opt;}
    void SetFillBplusTree(Int_t opt){fWriteVariableTreeBplus=opt;}
    void SetFillBsTree(Int_t opt){fWriteVariableTreeBs=opt;}
    void SetFillDstarTree(Int_t opt){fWriteVariableTreeDstar=opt;}
    void SetFillLc2V0bachelorTree(Int_t opt){fWriteVariableTreeLc2V0bachelor=opt;}
    void SetFillLbTree(Int_t opt){fWriteVariableTreeLb=opt;}
    void SetFillInclusiveJetTree(Int_t opt){fWriteVariableTreeInclusiveJet=opt;}
    void SetPIDoptD0Tree(Int_t opt){fPIDoptD0=opt;}
    void SetPIDoptDsTree(Int_t opt){fPIDoptDs=opt;}
    void SetPIDoptDplusTree(Int_t opt){fPIDoptDplus=opt;}
    void SetPIDoptLctopKpiTree(Int_t opt){fPIDoptLctopKpi=opt;}
    void SetPIDoptBplusTree(Int_t opt){fPIDoptBplus=opt;}
    void SetPIDoptBsTree(Int_t opt){fPIDoptBs=opt;}
    void SetPIDoptDstarTree(Int_t opt){fPIDoptDstar=opt;}
    void SetPIDoptLc2V0bachelorTree(Int_t opt){fPIDoptLc2V0bachelor=opt;}
    void SetPIDoptLbTree(Int_t opt){fPIDoptLb=opt;}
    void SetFillMCGenTrees(Bool_t fillMCgen) {fFillMCGenTrees=fillMCgen;}
  
    void SetMinJetPtCorr(double pt) { fMinJetPtCorr = pt; }
    void SetFillJetEtaPhi(bool b) { fFillJetEtaPhi = b; }
    void SetFillPtCorr(bool b) { fFillPtCorr = b; }
    void SetFillPtUncorr(bool b) { fFillPtUncorr = b; }
    void SetFillArea(bool b) { fFillArea = b; }
    void SetFillNConstituents(bool b) { fFillNConstituents = b; }
    void SetFillZLeading(bool b) { fFillZLeading = b; }
    void SetFillRadialMoment(bool b) { fFillRadialMoment = b; }
    void SetFillpTD(bool b) { fFillpTD = b; }
    void SetFillMass(bool b) { fFillMass = b; }
    void SetFillMatchingJetID(bool b) { fFillMatchingJetID = b; }

    void SetFillJets(bool b) {fFillJets = b; }
    void SetDoJetSubstructure(bool b) {fDoJetSubstructure = b; }
    void SetJetRadius(Double_t d) {fJetRadius = d; }
    void SetJetSubRadius(Double_t d) {fSubJetRadius = d; }
    void SetJetAlgorithm(Int_t i) {fJetAlgorithm = i; }
    void SetSubJetAlgorithm(Int_t i) {fSubJetAlgorithm = i; }
    void SetMinJetPt(Double_t d) {fMinJetPt = d; }
    void SetSoftDropZCut(Double_t d) {fSoftDropZCut = d; }
    void SetSoftDropBeta(Double_t d) {fSoftDropBeta = d; }
    void SetTrackingEfficiency(Double_t d) {fTrackingEfficiency = d;}
  
    void SetGoodTrackFilterBit(Int_t i) { fGoodTrackFilterBit = i; }
    void SetGoodTrackEtaRange(Double_t d) { fGoodTrackEtaRange = d; }
    void SetGoodTrackMinPt(Double_t d) { fGoodTrackMinPt = d; }
    void SetITSUpgradeProduction(Bool_t b) { fITSUpgradeProduction = b; }
    void SetITSUpgradePreSelect(Bool_t b) { fITSUpgradePreSelect = b; }
    void SetStoreOnlyHIJINGBackground(Bool_t b) {fStoreOnlyHIJINGBackground = b; }
    void SetLctopKpiPreselection(Bool_t d){ fPreSelectLctopKpi = d; }
    void SetFillInjCandHijingTrackCombi(Bool_t b){ fFillInjCandHijingTrackCombi = b; }
    
    void SetDsMassKKOption(AliHFTreeHandlerDstoKKpi::massKKopt opt) {fDsMassKKOpt=opt;}
    void SetLc2V0bachelorCalcSecoVtx(Int_t opt=1) {fLc2V0bachelorCalcSecoVtx=opt;}
    void SetLc2V0type(Int_t opt=1) {fV0typeForLc2V0bachelor=opt;}
    void SetOnFlySelectionValues(float invmass, float pt, float impparprod, float cosp, float cospxy){
      fInvMassOnFlyCut = invmass;
      fPtOnFlyCut = pt;
      fImpParProdOnFlyCut = impparprod;
      fCosPOnFlyCut = cosp;
      fCosPXYOnFlyCut = cospxy;
    }

    void SetTreeSingleTrackVarsOpt(Int_t opt) {fTreeSingleTrackVarsOpt=opt;}
  
    Int_t  GetSystem() const {return fSys;}
    Bool_t GetWriteOnlySignalTree() const {return fWriteOnlySignal;}
    
    void Process2Prong(TClonesArray *array2prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
    void Process3Prong(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
    void ProcessDstar(TClonesArray *arrayDstar, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield);
    void ProcessCasc(TClonesArray *arrayCasc, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield);
    void ProcessBplus(TClonesArray *array2prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
    void ProcessBs(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
    void ProcessLb(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader);
    void ProcessInclusiveJet(AliAODEvent *aod, TClonesArray *arrMC);
    void ProcessMCGen(TClonesArray *mcarray);
    void ProcessMCGenInclusiveJet(TClonesArray *mcarray);
  
    Bool_t CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau, Bool_t ITSUpgradeStudy);
    Bool_t IsCandidateFromHijing(AliAODRecoDecayHF *cand, AliAODMCHeader *mcHeader, TClonesArray* arrMC, AliAODTrack *tr = 0x0);
    
    void SelectGoodTrackForReconstruction(AliAODEvent *aod, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags);
    AliAODVertex* ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion);
  
    void SetNsigmaTPCDataDrivenCorrection(Int_t syst) {
        fEnableNsigmaTPCDataCorr=true; 
        fSystemForNsigmaTPCDataCorr=syst; 
    }

    void ApplyPhysicsSelectionOnline(bool apply=true) { fApplyPhysicsSelOnline = apply; }

    void EnableEventDownsampling(float fractokeep, unsigned long seed) {
        fEnableEventDownsampling = true;
        fFracToKeepEventDownsampling = fractokeep;
        fSeedEventDownsampling = seed;
    }

    // Particles (tracks or MC particles)
    //-----------------------------------------------------------------------------------------------
    void                        SetFillParticleTree(Bool_t b) {fFillParticleTree = b;}
    void                        SetFillTrackletTree(Bool_t b) {fFillTrackletTree = b;}
    AliParticleContainer*       AddParticleContainer(const char *n);
    AliTrackContainer*          AddTrackContainer(const char *n);
    AliMCParticleContainer*     AddMCParticleContainer(const char *n);
    AliParticleContainer*       GetParticleContainer(Int_t i=0) const;
    AliParticleContainer*       GetParticleContainer(const char* name) const;
    void                        FillParticleTree();

    // Jets
    //-----------------------------------------------------------------------------------------------
    void SetFillNJetTrees(Int_t n){fWriteNJetTrees=n;}
    AliJetContainer* AddJetContainer(AliJetContainer::EJetType_t jetType, AliJetContainer::EJetAlgo_t jetAlgo, AliJetContainer::ERecoScheme_t recoScheme, Double_t radius, UInt_t accType, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag = "Jet");
    AliJetContainer* AddJetContainer(const char *n, UInt_t accType, Float_t jetRadius);
    AliJetContainer* GetJetContainer(Int_t i=0) const;
    void FillJetTree();
  
    
    unsigned long GetEvID();
    
private:
    
    AliAnalysisTaskSEHFTreeCreator(const AliAnalysisTaskSEHFTreeCreator&);
    AliAnalysisTaskSEHFTreeCreator& operator=(const AliAnalysisTaskSEHFTreeCreator&);
    
    unsigned int            fEventNumber;
    TH1F                    *fNentries;                            //!<!   histogram with number of events on output slot 1
    TH2F                    *fHistoNormCounter;                    //!<!   histogram with number of events on output slot 1
    TList                   *fListCuts;                            //      list of cuts sent to output slot 2
    AliRDHFCutsD0toKpi      *fFiltCutsD0toKpi;                     //      D0toKpi filtering (or loose) cuts
    AliRDHFCutsDstoKKpi     *fFiltCutsDstoKKpi;                    //      DstoKKpi filtering (or loose) cuts
    AliRDHFCutsDplustoKpipi *fFiltCutsDplustoKpipi;                //      DplustoKpipi filtering (or loose) cuts
    AliRDHFCutsLctopKpi     *fFiltCutsLctopKpi    ;                //      LctopKpi filtering (or loose) cuts
    AliRDHFCutsD0toKpi      *fFiltCutsBplustoD0pi;                 //      BplustoD0pi filtering (or loose) cuts
    AliRDHFCutsDstoKKpi     *fFiltCutsBstoDspi;                    //      BstoDspi filtering (or loose) cuts
    AliRDHFCutsDStartoKpipi *fFiltCutsDstartoKpipi;                //      DstartoKpipi filtering (or loose) cuts
    AliRDHFCutsLctoV0       *fFiltCutsLc2V0bachelor;               //      Lc2V0bachelor filtering (or loose) cuts
    AliRDHFCutsLctopKpi     *fFiltCutsLbtoLcpi;                    //      LbtoLcpi filtering (or loose) cuts
    AliRDHFCutsD0toKpi      *fCutsD0toKpi;                         //      D0toKpi analysis cuts
    AliRDHFCutsDstoKKpi     *fCutsDstoKKpi;                        //      DstoKKpi analysis cuts
    AliRDHFCutsDplustoKpipi *fCutsDplustoKpipi;                    //      DplustoKpipi analysis cuts
    AliRDHFCutsLctopKpi     *fCutsLctopKpi;                        //      LctopKpi analysis cuts
    AliRDHFCutsD0toKpi      *fCutsBplustoD0pi;                     //      BplustoD0pi analysis cuts
    AliRDHFCutsDstoKKpi     *fCutsBstoDspi;                        //      BstoDspi analysis cuts
    AliRDHFCutsDStartoKpipi *fCutsDstartoKpipi;                    //      DstartoKpipi analysis cuts
    AliRDHFCutsLctoV0       *fCutsLc2V0bachelor;                   //      Lc2V0bachelor analysis cuts
    AliRDHFCutsLctopKpi     *fCutsLbtoLcpi;                        //      LbtoLcpi analysis cuts
    AliRDHFCuts             *fEvSelectionCuts;                     //      Event selection cuts
    Bool_t                  fReadMC;                               //     flag for MC array: kTRUE = read it, kFALSE = do not read it
    TList                   *fListCounter;                         //!<!   list for normalization counter on output slot 3
    AliNormalizationCounter *fCounter;                             //!<!   AliNormalizationCounter
    Bool_t                  fUseSelectionBit;
    Int_t                   fSys;                                  // fSys=0 -> p-p; fSys=1 ->PbPb
    Int_t                   fAODProtection;                        // flag to activate protection against AOD-dAOD mismatch.
                                                                   // -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
    Int_t                   fWriteVariableTreeD0;                  // flag to decide whether to write the candidate variables on a tree variables
                                                                   // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeDs;                  // flag to decide whether to write the candidate variables on a tree variables
 														                                       // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeDplus;               // flag to decide whether to write the candidate variables on a tree variables
    													                                     // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeLctopKpi;            // flag to decide whether to write the candidate variables on a tree variables
    													                                     // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeBplus;              // flag to decide whether to write the candidate variables on a tree variables
                                                                   // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeBs;                   // flag to decide whether to write the candidate variables on a tree variables
                                                                   // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeDstar;                // flag to decide whether to write the candidate variables on a tree variables
    													                                     // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeLc2V0bachelor;        // flag to decide whether to write the candidate variables on a tree variables
    													                                     // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                   fWriteVariableTreeLb;                   // flag to decide whether to write the candidate variables on a tree variables
                                                                   // 0 don't fill
                                                                   // 1 fill standard tree

    Int_t                   fWriteVariableTreeInclusiveJet;        // flag to decide whether to write the candidate variables on a tree variables
                                                                   // 0 don't fill
                                                                   // 1 fill standard tree


    TTree                   *fVariablesTreeD0;                     //!<! tree of the candidate variables
    TTree                   *fVariablesTreeDs;                     //!<! tree of the candidate variables
    TTree                   *fVariablesTreeDplus;                  //!<! tree of the candidate variables
    TTree                   *fVariablesTreeLctopKpi;               //!<! tree of the candidate variables
    TTree                   *fVariablesTreeBplus;                  //!<! tree of the candidate variables
    TTree                   *fVariablesTreeBs;                     //!<! tree of the candidate variables
    TTree                   *fVariablesTreeDstar;                  //!<! tree of the candidate variables
    TTree                   *fVariablesTreeLc2V0bachelor;          //!<! tree of the candidate variables
    TTree                   *fVariablesTreeLb;                     //!<! tree of the candidate variables
    TTree                   *fVariablesTreeInclusiveJet;           //!<! tree of the candidate variables
    TTree                   *fGenTreeD0;                           //!<! tree of the gen D0 variables
    TTree                   *fGenTreeDs;                           //!<! tree of the gen Ds variables
    TTree                   *fGenTreeDplus;                        //!<! tree of the gen D+ variables
    TTree                   *fGenTreeLctopKpi;                     //!<! tree of the gen LctopKpi variables
    TTree                   *fGenTreeBplus;                        //!<! tree of the gen B+ variables
    TTree                   *fGenTreeBs;                           //!<! tree of the gen Bs variables
    TTree                   *fGenTreeDstar;                        //!<! tree of the gen Dstar variables
    TTree                   *fGenTreeLc2V0bachelor;                //!<! tree of the gen Lc2V0bachelor variables
    TTree                   *fGenTreeLb;                           //!<! tree of the gen Lb variables
    TTree                   *fGenTreeInclusiveJet;                 //!<! tree of the gen Inclusive variables
    TTree                   *fTreeEvChar;                          //!<! tree of event variables
    bool                    fWriteOnlySignal;
    AliHFTreeHandlerD0toKpi        *fTreeHandlerD0;                //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstoKKpi       *fTreeHandlerDs;                //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDplustoKpipi   *fTreeHandlerDplus;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLctopKpi       *fTreeHandlerLctopKpi;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerBplustoD0pi    *fTreeHandlerBplus;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerBstoDspi       *fTreeHandlerBs;                //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstartoKpipi   *fTreeHandlerDstar;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLc2V0bachelor  *fTreeHandlerLc2V0bachelor;     //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLbtoLcpi       *fTreeHandlerLb;                //!<! handler object for the tree with topological variables
    AliHFTreeHandlerInclusiveJet   *fTreeHandlerInclusiveJet;      //!<! handler object for the tree 
    AliHFTreeHandlerD0toKpi        *fTreeHandlerGenD0;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstoKKpi       *fTreeHandlerGenDs;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDplustoKpipi   *fTreeHandlerGenDplus;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLctopKpi       *fTreeHandlerGenLctopKpi;       //!<! handler object for the tree with topological variables
    AliHFTreeHandlerBplustoD0pi    *fTreeHandlerGenBplus;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerBstoDspi       *fTreeHandlerGenBs;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstartoKpipi   *fTreeHandlerGenDstar;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLc2V0bachelor  *fTreeHandlerGenLc2V0bachelor;  //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLbtoLcpi       *fTreeHandlerGenLb;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerInclusiveJet   *fTreeHandlerGenInclusiveJet;   //!<! handler object for the tree 
    AliPIDResponse          *fPIDresp;                             /// PID response
    int                     fPIDoptD0;                             /// PID option for D0 tree
    int                     fPIDoptDs;                             /// PID option for Ds tree
    int                     fPIDoptDplus;                          /// PID option for D+ tree
    int                     fPIDoptLctopKpi;                       /// PID option for Lc2pKpi tree
    int                     fPIDoptBplus;                          /// PID option for B+ tree
    int                     fPIDoptBs;                             /// PID option for Bs tree
    int                     fPIDoptDstar;                          /// PID option for D* tree
    int                     fPIDoptLc2V0bachelor;                  /// PID option for Lc2V0bachelor tree
    int                     fPIDoptLb;                             /// PID option for Lb tree
    Float_t                 fCentrality;                           /// event centrality
    Float_t                 fzVtxReco;                             /// reconstructed Zvtx
    Float_t                 fzVtxGen;                              /// generated Zvtx
    Int_t                   fNcontributors;                        /// number of contributors
    Int_t                   fNtracks;                              /// number of tracks
    Int_t                   fIsEvRej;                              /// flag with information about rejection of the event
    Int_t                   fIsEvRej_INT7;                         /// flag with information about rejection of the event
    Int_t                   fIsEvRej_HighMultSPD;                  /// flag with information about rejection of the event
    Int_t                   fIsEvRej_HighMultV0;                   /// flag with information about rejection of the event
    Int_t                   fIsEvRej_EMCEJE;                       /// flag with information about rejection of the event
    Bool_t                  fIsEvSel_INT7;                         /// boolean whether event accept for INT7
    Bool_t                  fIsEvSel_HighMultSPD;                  /// boolean whether event accept for SHM
    Bool_t                  fIsEvSel_HighMultV0;                   /// boolean whether event accept for VHM
    Bool_t                  fIsEvSel_EMCEJE;                       /// boolean whether event accept for EMCEJE
    Double_t                fCross_Section;                        /// cross section of pt hard bin event
    Int_t                   fTrials;                               /// Trials of pt hard bin event
    Double_t                fpthard;                               /// Pt hard bin of event
    Int_t                   fRunNumber;                            /// run number
    Int_t                   fRunNumberCDB;                         /// run number (for OCDB)
    UShort_t                fBC;                                   /// bunch crossing number
    Int_t                   fOrbit;                                /// orbit
    Int_t                   fPeriod;                               /// period
    Int_t                   fEventID;                              /// event ID (for guaranteed uniqueness combine with ext ID)
    Int_t                   fEventIDExt;                           /// upper 32-bit of event ID
    Long64_t                fEventIDLong;                          /// single unique event id (long64)
    TString                 fFileName;
    unsigned int            fDirNumber;
    Int_t                   fnTracklets;                           /// number of tracklets
    Int_t                   fnTrackletsCorr;                       /// number of tracklets (corrected)
    Int_t                   fnTrackletsCorrSHM;                    /// number of tracklets (corrected)
    Double_t                fRefMult;                              /// reference multiplicity
    Double_t                fRefMultSHM;                           /// reference multiplicity
    Int_t                   fnV0A;                                 /// V0A multiplicity
    Int_t                   fMultGen;                              /// generated multiplicity around mid-rapidity [-1,1]
    Int_t                   fMultGenV0A;                           /// generated multiplicity in V0A range
    Int_t                   fMultGenV0C;                           /// generated multiplicity in V0C range
    ULong64_t               fTriggerMask;                          /// Trigger mask bitmap
    Bool_t                  fTriggerOnlineINT7;                    /// Flag explicitly whether bitmap contains INT7
    Bool_t                  fTriggerOnlineHighMultSPD;             /// Flag explicitly whether bitmap contains HighMultSPD
    Bool_t                  fTriggerOnlineHighMultV0;              /// Flag explicitly whether bitmap kHighMultV0
    Bool_t                  fTriggerOnlineEMCALEJ1;                /// Flag explicitly whether bitmap kEMCEJ1
    Bool_t                  fTriggerOnlineEMCALEJ2;                /// Flag explicitly whether bitmap kEMCEJ2
    Bool_t                  fTriggerOnlineDCALDJ1;                 /// Flag explicitly whether bitmap DCALDJ1
    Bool_t                  fTriggerOnlineDCALDJ2;                 /// Flag explicitly whether bitmap DCALDJ2
    Bool_t                  fTriggerBitINT7;                       /// Flag explicitly whether bitmap contains INT7
    Bool_t                  fTriggerBitHighMultSPD;                /// Flag explicitly whether bitmap contains HighMultSPD
    Bool_t                  fTriggerBitHighMultV0;                 /// Flag explicitly whether bitmap kHighMultV0
    Bool_t                  fTriggerBitCentral;                    /// Flag explicitly whether bitmap contains kCentral
    Bool_t                  fTriggerBitSemiCentral;                /// Flag explicitly whether bitmap contains kSemiCentral
    Bool_t                  fTriggerBitEMCEJE;                     /// Flag explicitly whether bitmap contains kEMCEJE
    TString                 fTriggerClasses;                       /// Collect all trigger classes
    Bool_t                  fTriggerClassINT7;                     /// Flag explicitly whether classes contain INT7
    Bool_t                  fTriggerClassHighMultSPD;              /// Flag explicitly whether classes contain HighMultSPD
    Bool_t                  fTriggerClassHighMultV0m;              /// Flag explicitly whether classes contain HighMultV0
    Bool_t                  fTriggerClassEMCALEJ1;                 /// Flag explicitly whether classes contain EMCALEJ1 (contained in AliVEvent::kEMCEJE)
    Bool_t                  fTriggerClassEMCALEJ2;                 /// Flag explicitly whether classes contain EMCALEJ2 (conatined in AliVEvent::kEMCEJE)
    Bool_t                  fTriggerClassDCALDJ1;                  /// Flag explicitly whether classes contain DCALDJ1 (contained in AliVEvent::kEMCEJE)
    Bool_t                  fTriggerClassDCALDJ2;                  /// Flag explicitly whether classes contain DCALDJ2 (conatined in AliVEvent::kEMCEJE)
    Int_t                   fnV0M;                                 /// V0M multiplicity
    Int_t                   fnV0MEq;                               /// V0M multiplicity (equalized)
    Int_t                   fnV0MCorr;                             /// V0M multiplicity (corrected)
    Int_t                   fnV0MEqCorr;                           /// V0M multiplicity (equalized + corrected)
    Float_t                 fPercV0M;                              /// V0M multiplicity percentile
    Float_t                 fMultV0M;                              /// V0M multiplicity from mult selection task

    Bool_t                  fFillMCGenTrees;                       /// flag to enable fill of the generated trees
  
    Int_t                   fDsMassKKOpt;                          /// option for Ds massKK (mass or delta mass)
    Int_t                   fLc2V0bachelorCalcSecoVtx;             /// option to calculate the secondary vertex for Lc2V0bachelor. False by default, has to be added to AddTask in case we want to start using it.
    Int_t                   fV0typeForLc2V0bachelor;               /// option to select Offline+OnTheFly (0), only Offline (1=default), only OnTheFly (2) V0's for the Lc->V0bachelor decay
    Float_t                 fInvMassOnFlyCut;                      ///Cut on invariant mass for on fly hadron selection
    Float_t                 fPtOnFlyCut;                           ///Cut on pT for on fly hadron  selection
    Float_t                 fImpParProdOnFlyCut;                   ///Cut on d0xd0 for on fly hadron  selection
    Float_t                 fCosPOnFlyCut;                         ///Cut on cos pointing angle for on fly hadron  selection
    Float_t                 fCosPXYOnFlyCut;                       ///Cut on cos pointing angle xy for on fly hadron selection
  
    Int_t                   fTreeSingleTrackVarsOpt;               /// option for single-track variables to be filled in the trees

    Double_t                fJetRadius;                            /// Setting the radius for jet finding
    Double_t                fSubJetRadius;                         /// Setting the radius for subjet finding
    Int_t                   fJetAlgorithm;                         /// Setting the jet finding algorithm
    Int_t                   fSubJetAlgorithm;                      /// Setting the jet finding algorithm
    Double_t                fMinJetPt;                             /// Setting the jet finding min pT
    Double_t                fSoftDropZCut;                         /// setting the soft drop z parameter
    Double_t                fSoftDropBeta;                         /// setting the soft drop beta parameter
    Double_t                fTrackingEfficiency;                   /// Setting the jet finding tracking efficiency
  
    Int_t                   fGoodTrackFilterBit;                   /// Setting filter bit for bachelor on-the-fly reconstruction candidate
    Double_t                fGoodTrackEtaRange;                    /// Setting eta-range for bachelor on-the-fly reconstruction candidate
    Double_t                fGoodTrackMinPt;                       /// Setting min pT for bachelor on-the-fly reconstruction candidate
    Bool_t                  fITSUpgradeProduction;                 /// Setting for analysing an ITS Upgrade production
    Bool_t                  fITSUpgradePreSelect;                  /// Setting to enable ITSUpgrade Preselect function
    Bool_t                  fStoreOnlyHIJINGBackground;            /// Setting to store only HIJING background candidates
    Bool_t                  fPreSelectLctopKpi;                    /// Setting for reduce the filling of reconstucted candidates  in an ITS Upgrade production
    Bool_t                  fFillInjCandHijingTrackCombi;          /// Setting to store injected candidate + HIJING track for background shape studies

    // Particles (tracks / MC particles)
    // Add a single AliTrackContainer and/or AliMCParticleContainer to select particles
    // A separate (identical) AliParticleTreeHandler will be used to fill each tree.
    //-----------------------------------------------------------------------------------------------
    bool                    fFillParticleTree;                     ///< Store tree of all tracks inside the jet
    bool                    fFillTrackletTree;                     ///< Store tree of all tracklets
  
    TTree*                  fVariablesTreeParticle;                //!<! Particle tree
    TTree*                  fVariablesTreeTracklet;                //!<! Tracklet tree
    TTree*                  fVariablesTreeGenParticle;             //!<! MC particle tree
  
    AliParticleTreeHandler* fTreeHandlerParticle;                  //!<! handler object for particle tree
    AliTrackletTreeHandler* fTreeHandlerTracklet;                  //!<! handler object for tracklet tree
    AliParticleTreeHandler* fTreeHandlerGenParticle;               //!<! handler object for MC particle tree
  
    TObjArray               fParticleCollArray;                    ///< particle/track collection array
  
    // Jets
    //-----------------------------------------------------------------------------------------------
  
    // Write N jet trees, according to constructor argument. Should match number of jet containers added.
    // If fFillJetConstituentTrees is true, then also fill N separate trees of jet constituent info.
    Int_t                   fWriteNJetTrees;                       ///< number of jet trees to write
    bool                    fFillJetConstituentTrees;              ///< Store tree of all tracks inside the jet
  
    std::vector<TTree*>     fVariablesTreeJet;                     //!<! vector of jet trees
    std::vector<TTree*>     fVariablesTreeJetConstituent;          //!<! vector of jet constituent trees

    std::vector<AliJetTreeHandler*> fTreeHandlerJet;               //!<! vector of handler objects for jet tree
  
    // Jet container and array
    Bool_t                  fLocalInitialized;                     ///< whether or not the task has been already initialized
    TObjArray               fJetCollArray;                         ///< array of jet containers
    double                  fMinJetPtCorr;                         ///< Min jet Pt (background subtracted) to fill jet into tree
  
    // Jet background subtraction
    TString                 fRhoName;                              ///<  rho name
    AliRhoParameter        *fRho;                                  //!<! event rho
    Double_t                fRhoVal;                               //!<! event rho value
  
    // Fill jet tree according to the below flags. By default, it only contains: event id, jet id
    bool                    fFillJetEtaPhi;                        ///< Jet eta/phi
    bool                    fFillPtCorr;                           ///< Pt of the jet (GeV/c) (background subtracted)
    bool                    fFillPtUncorr;                         ///< Pt of the jet (GeV/c) (not background subtracted)
    bool                    fFillArea;                             ///< Area
    bool                    fFillNConstituents;                    ///< N constituents
    bool                    fFillZLeading;                         ///< ZLeading
    bool                    fFillRadialMoment;                     ///< Radial moment
    bool                    fFillpTD;                              ///< pT,D
    bool                    fFillMass;                             ///< Mass
    bool                    fFillMatchingJetID;                    ///< jet matching


    
    bool                    fFillJets;                             /// FillJetInfo
    bool                    fDoJetSubstructure;                    /// FillJetSubstructure
    
  
    bool fEnableNsigmaTPCDataCorr; /// flag to enable data-driven NsigmaTPC correction
    int fSystemForNsigmaTPCDataCorr; /// system for data-driven NsigmaTPC correction

    std::map<std::string, TProfile*> fMultEstimatorAvg;
    std::map<std::string, TProfile*> fMultEstimatorAvgSHM;
    bool fCorrNtrVtx;
    bool fCorrV0MVtx;

    bool fApplyPhysicsSelOnline;                                   /// flag to apply physics selection in the task
    bool fEnableEventDownsampling;                                 /// flag to apply event downsampling
    float fFracToKeepEventDownsampling;                            /// fraction of events to be kept by event downsampling
    unsigned long fSeedEventDownsampling;                          /// seed for event downsampling

    AliCDBEntry *fCdbEntry;

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEHFTreeCreator,29);
    /// \endcond
};

#endif

