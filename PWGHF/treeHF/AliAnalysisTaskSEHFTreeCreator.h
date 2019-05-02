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
///*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsBPlustoD0Pi.h"
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
#include "AliHFTreeHandlerDstartoKpipi.h"
#include "AliHFTreeHandlerLc2V0bachelor.h"
#include "AliJetTreeHandler.h"
#include "AliJetContainer.h"

class AliAODEvent;
class TClonesArray;
class AliEmcalJet;
class AliRhoParameter;

class AliAnalysisTaskSEHFTreeCreator : public AliAnalysisTaskSE
{
public:


    
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
    
    
    void SetReadMC(Bool_t opt=kFALSE){fReadMC=opt;}
    void SetSystem(Int_t opt){fSys=opt;}
    void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
    void SetWriteOnlySignalTree(Bool_t opt){fWriteOnlySignal=opt;}
    void SetFillD0Tree(Int_t opt){fWriteVariableTreeD0=opt;}
    void SetFillDsTree(Int_t opt){fWriteVariableTreeDs=opt;}
    void SetFillDplusTree(Int_t opt){fWriteVariableTreeDplus=opt;}
    void SetFillLctopKpiTree(Int_t opt){fWriteVariableTreeLctopKpi=opt;}
    void SetFillBplusTree(Int_t opt){fWriteVariableTreeBplus=opt;}
    void SetFillDstarTree(Int_t opt){fWriteVariableTreeDstar=opt;}
    void SetFillLc2V0bachelorTree(Int_t opt){fWriteVariableTreeLc2V0bachelor=opt;}
    void SetPIDoptD0Tree(Int_t opt){fPIDoptD0=opt;}
    void SetPIDoptDsTree(Int_t opt){fPIDoptDs=opt;}
    void SetPIDoptDplusTree(Int_t opt){fPIDoptDplus=opt;}
    void SetPIDoptLctopKpiTree(Int_t opt){fPIDoptLctopKpi=opt;}
    void SetPIDoptBplusTree(Int_t opt){fPIDoptBplus=opt;}
    void SetPIDoptDstarTree(Int_t opt){fPIDoptDstar=opt;}
    void SetPIDoptLc2V0bachelorTree(Int_t opt){fPIDoptLc2V0bachelor=opt;}
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
  
    void SetDsMassKKOption(AliHFTreeHandlerDstoKKpi::massKKopt opt) {fDsMassKKOpt=opt;}
    void SetLc2V0bachelorCalcSecoVtx(Int_t opt=1) {fLc2V0bachelorCalcSecoVtx=opt;}
  
    void SetTreeSingleTrackVarsOpt(Int_t opt) {fTreeSingleTrackVarsOpt=opt;}
  
    Int_t  GetSystem() const {return fSys;}
    Bool_t GetWriteOnlySignalTree() const {return fWriteOnlySignal;}
    
    void Process2Prong(TClonesArray *array2prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield);
    void Process3Prong(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield);
    void ProcessDstar(TClonesArray *arrayDstar, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield);
    void ProcessCasc(TClonesArray *arrayCasc, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield);
    void ProcessMCGen(TClonesArray *mcarray);
  
    Bool_t CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
    AliAODVertex* ReconstructBplusVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion);
  
    void SetNsigmaTPCDataDrivenCorrection(Int_t syst) {
        fEnableNsigmaTPCDataCorr=true; 
        fSystemForNsigmaTPCDataCorr=syst; 
    }

    // Jets
    //-----------------------------------------------------------------------------------------------
    void SetFillNJetTrees(Int_t n){fWriteNJetTrees=n;}
  
    AliJetContainer* AddJetContainer(AliJetContainer::EJetType_t jetType, AliJetContainer::EJetAlgo_t jetAlgo, AliJetContainer::ERecoScheme_t recoScheme, Double_t radius, UInt_t accType, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag = "Jet");
    AliJetContainer* AddJetContainer(const char *n, UInt_t accType, Float_t jetRadius);
    AliJetContainer* GetJetContainer(Int_t i=0) const;
    void FillJetTree();
  
    
    unsigned int GetEvID();
    
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
    AliRDHFCutsBPlustoD0Pi  *fFiltCutsBplustoD0pi;                 //      BplustoD0pi filtering (or loose) cuts
    AliRDHFCutsDStartoKpipi *fFiltCutsDstartoKpipi;                //      DstartoKpipi filtering (or loose) cuts
    AliRDHFCutsLctoV0       *fFiltCutsLc2V0bachelor;               //      Lc2V0bachelor filtering (or loose) cuts
    AliRDHFCutsD0toKpi      *fCutsD0toKpi;                         //      D0toKpi analysis cuts
    AliRDHFCutsDstoKKpi     *fCutsDstoKKpi;                        //      DstoKKpi analysis cuts
    AliRDHFCutsDplustoKpipi *fCutsDplustoKpipi;                    //      DplustoKpipi analysis cuts
    AliRDHFCutsLctopKpi     *fCutsLctopKpi;                        //      LctopKpi analysis cuts
    AliRDHFCutsBPlustoD0Pi  *fCutsBplustoD0pi;                     //      BplustoD0pi analysis cuts
    AliRDHFCutsDStartoKpipi *fCutsDstartoKpipi;                    //      DstartoKpipi analysis cuts
    AliRDHFCutsLctoV0       *fCutsLc2V0bachelor;                   //      Lc2V0bachelor analysis cuts
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
    Int_t                    fWriteVariableTreeBplus;              // flag to decide whether to write the candidate variables on a tree variables
                                                                   // 0 don't fill
                                                                   // 1 fill standard tree

    Int_t                  fWriteVariableTreeDstar;                // flag to decide whether to write the candidate variables on a tree variables
    													                                     // 0 don't fill
                                                                   // 1 fill standard tree
    Int_t                  fWriteVariableTreeLc2V0bachelor;        // flag to decide whether to write the candidate variables on a tree variables
    													                                     // 0 don't fill
                                                                   // 1 fill standard tree

    TTree                   *fVariablesTreeD0;                     //!<! tree of the candidate variables
    TTree                   *fVariablesTreeDs;                     //!<! tree of the candidate variables
    TTree                   *fVariablesTreeDplus;                  //!<! tree of the candidate variables
    TTree                   *fVariablesTreeLctopKpi;               //!<! tree of the candidate variables
    TTree                   *fVariablesTreeBplus;                  //!<! tree of the candidate variables
    TTree                   *fVariablesTreeDstar;                  //!<! tree of the candidate variables
    TTree                   *fVariablesTreeLc2V0bachelor;          //!<! tree of the candidate variables
    TTree                   *fGenTreeD0;                           //!<! tree of the gen D0 variables
    TTree                   *fGenTreeDs;                           //!<! tree of the gen Ds variables
    TTree                   *fGenTreeDplus;                        //!<! tree of the gen D+ variables
    TTree                   *fGenTreeLctopKpi;                     //!<! tree of the gen LctopKpi variables
    TTree                   *fGenTreeBplus;                        //!<! tree of the gen B+ variables
    TTree                   *fGenTreeDstar;                        //!<! tree of the gen Dstar variables
    TTree                   *fGenTreeLc2V0bachelor;                //!<! tree of the gen Lc2V0bachelor variables
    TTree                   *fTreeEvChar;                          //!<! tree of event variables
    bool                    fWriteOnlySignal;
    AliHFTreeHandlerD0toKpi        *fTreeHandlerD0;                //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstoKKpi       *fTreeHandlerDs;                //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDplustoKpipi   *fTreeHandlerDplus;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLctopKpi       *fTreeHandlerLctopKpi;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerBplustoD0pi    *fTreeHandlerBplus;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstartoKpipi   *fTreeHandlerDstar;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLc2V0bachelor  *fTreeHandlerLc2V0bachelor;     //!<! handler object for the tree with topological variables
    AliHFTreeHandlerD0toKpi        *fTreeHandlerGenD0;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstoKKpi       *fTreeHandlerGenDs;             //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDplustoKpipi   *fTreeHandlerGenDplus;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLctopKpi       *fTreeHandlerGenLctopKpi;       //!<! handler object for the tree with topological variables
    AliHFTreeHandlerBplustoD0pi    *fTreeHandlerGenBplus;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerDstartoKpipi   *fTreeHandlerGenDstar;          //!<! handler object for the tree with topological variables
    AliHFTreeHandlerLc2V0bachelor  *fTreeHandlerGenLc2V0bachelor;  //!<! handler object for the tree with topological variables
    AliPIDResponse          *fPIDresp;                             /// PID response
    int                     fPIDoptD0;                             /// PID option for D0 tree
    int                     fPIDoptDs;                             /// PID option for Ds tree
    int                     fPIDoptDplus;                          /// PID option for D+ tree
    int                     fPIDoptLctopKpi;                       /// PID option for Lc2pKpi tree
    int                     fPIDoptBplus;                          /// PID option for B+ tree
    int                     fPIDoptDstar;                          /// PID option for D* tree
    int                     fPIDoptLc2V0bachelor;                  /// PID option for Lc2V0bachelor tree
    Float_t                 fCentrality;                           /// event centrality
    Float_t                 fzVtxReco;                             /// reconstructed Zvtx
    Float_t                 fzVtxGen;                              /// generated Zvtx
    Int_t                   fNcontributors;                        /// number of contributors
    Int_t                   fNtracks;                              /// number of tracks
    Int_t                   fIsEvRej;                              /// flag with information about rejection of the event
    Int_t                   fRunNumber;                            /// run number
    UInt_t                  fEventID;                              /// event ID (unique when combined with run number)
    TString                 fFileName;
    unsigned int            fDirNumber;
    Int_t                   fnTracklets;                           /// number of tracklets
    Int_t                   fnV0A;                                 /// V0A multiplicity 

    Bool_t                  fFillMCGenTrees;                       /// flag to enable fill of the generated trees
  
    Int_t                   fDsMassKKOpt;                          /// option for Ds massKK (mass or delta mass)
    Int_t                   fLc2V0bachelorCalcSecoVtx;             /// option to calculate the secondary vertex for Lc2V0bachelor. False by default, has to be added to AddTask in case we want to start using it.
  
    Int_t                   fTreeSingleTrackVarsOpt;               /// option for single-track variables to be filled in the trees
  
  
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
  
    bool fEnableNsigmaTPCDataCorr; /// flag to enable data-driven NsigmaTPC correction
    int fSystemForNsigmaTPCDataCorr; /// system for data-driven NsigmaTPC correction

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSEHFTreeCreator,13);
    /// \endcond
};

#endif

