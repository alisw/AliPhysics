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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
//
//
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
////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TProcessID.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <THnSparse.h>
#include <TRandom3.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TChain.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliExternalTrackParam.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliNormalizationCounter.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskSE.h"
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
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliAnalysisTaskSEHFTreeCreator.h"
#include "AliAODPidHF.h"
#include "AliESDUtils.h"
#include "AliMultSelection.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerInput.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEHFTreeCreator);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreator::AliAnalysisTaskSEHFTreeCreator():
AliAnalysisTaskSEHFTreeCreator("", nullptr, 0, false)
{
  
  /// Default constructor
  
  fParticleCollArray.SetOwner(kTRUE);
  fJetCollArray.SetOwner(kTRUE);
  
}
//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreator::AliAnalysisTaskSEHFTreeCreator(const char *name, TList *cutsList, int fillNJetTrees, bool fillJetConstituentTrees):
AliAnalysisTaskSE(name),
fEventNumber(0),
fNentries(0x0),
fHistoNormCounter(0x0),
fListCuts(cutsList),
fFiltCutsD0toKpi(0x0),
fFiltCutsDstoKKpi(0x0),
fFiltCutsDplustoKpipi(0x0),
fFiltCutsLctopKpi(0x0),
fFiltCutsBplustoD0pi(0x0),
fFiltCutsBstoDspi(0x0),
fFiltCutsDstartoKpipi(0x0),
fFiltCutsLc2V0bachelor(0x0),
fFiltCutsLbtoLcpi(0x0),
fCutsD0toKpi(0x0),
fCutsDstoKKpi(0x0),
fCutsDplustoKpipi(0x0),
fCutsLctopKpi(0x0),
fCutsBplustoD0pi(0x0),
fCutsBstoDspi(0x0),
fCutsDstartoKpipi(0x0),
fCutsLc2V0bachelor(0x0),
fCutsLbtoLcpi(0x0),
fEvSelectionCuts(0x0),
fReadMC(0),
fListCounter(0x0),
fCounter(0x0),
fUseSelectionBit(kTRUE),
fSys(0),
fAODProtection(1),
fWriteVariableTreeD0(0),
fWriteVariableTreeDs(0),
fWriteVariableTreeDplus(0),
fWriteVariableTreeLctopKpi(0),
fWriteVariableTreeBplus(0),
fWriteVariableTreeBs(0),
fWriteVariableTreeDstar(0),
fWriteVariableTreeLc2V0bachelor(0),
fWriteVariableTreeLb(0),
fWriteVariableTreeInclusiveJet(0),
fVariablesTreeD0(0x0),
fVariablesTreeDs(0x0),
fVariablesTreeDplus(0x0),
fVariablesTreeLctopKpi(0x0),
fVariablesTreeBplus(0x0),
fVariablesTreeBs(0x0),
fVariablesTreeDstar(0x0),
fVariablesTreeLc2V0bachelor(0x0),
fVariablesTreeLb(0x0),
fVariablesTreeInclusiveJet(0x0),
fGenTreeD0(0x0),
fGenTreeDs(0x0),
fGenTreeDplus(0x0),
fGenTreeLctopKpi(0x0),
fGenTreeBplus(0x0),
fGenTreeBs(0x0),
fGenTreeDstar(0x0),
fGenTreeLc2V0bachelor(0x0),
fGenTreeLb(0x0),
fGenTreeInclusiveJet(0x0),
fTreeEvChar(0x0),
fWriteOnlySignal(kFALSE),
fTreeHandlerD0(0x0),
fTreeHandlerDs(0x0),
fTreeHandlerDplus(0x0),
fTreeHandlerLctopKpi(0x0),
fTreeHandlerBplus(0x0),
fTreeHandlerBs(0x0),
fTreeHandlerDstar(0x0),
fTreeHandlerLc2V0bachelor(0x0),
fTreeHandlerLb(0x0),
fTreeHandlerInclusiveJet(0x0),
fTreeHandlerGenD0(0x0),
fTreeHandlerGenDs(0x0),
fTreeHandlerGenDplus(0x0),
fTreeHandlerGenLctopKpi(0x0),
fTreeHandlerGenBplus(0x0),
fTreeHandlerGenBs(0x0),
fTreeHandlerGenDstar(0x0),
fTreeHandlerGenLc2V0bachelor(0x0),
fTreeHandlerGenLb(0x0),
fTreeHandlerGenInclusiveJet(0x0),
fPIDresp(0x0),
fPIDoptD0(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDs(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDplus(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptLctopKpi(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptBplus(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptBs(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDstar(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptLc2V0bachelor(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptLb(AliHFTreeHandler::kRawAndNsigmaPID),
fCentrality(-999.),
fzVtxReco(0.),
fzVtxGen(0.),
fNcontributors(0),
fNtracks(0),
fIsEvRej(0),
fIsEvRej_INT7(0),
fIsEvRej_HighMultSPD(0),
fIsEvRej_HighMultV0(0),
fIsEvRej_EMCEJE(0),
fIsEvSel_INT7(false),
fIsEvSel_HighMultSPD(false),
fIsEvSel_HighMultV0(false),
fIsEvSel_EMCEJE(false),
fCross_Section(-1.),
fTrials(-1),
fpthard(-1.),
fRunNumber(0),
fRunNumberCDB(0),
fBC(0),
fOrbit(0),
fPeriod(0),
fEventID(0),
fEventIDExt(0),
fEventIDLong(0),
fFileName(""),
fDirNumber(0),
fnTracklets(0),
fnTrackletsCorr(0),
fnTrackletsCorrSHM(0),
fRefMult(9.26),
fRefMultSHM(9.26),
fnV0A(0),
fMultGen(0),
fMultGenV0A(0),
fMultGenV0C(0),
fTriggerMask(0),
fTriggerOnlineINT7(false),
fTriggerOnlineHighMultSPD(false),
fTriggerOnlineHighMultV0(false),
fTriggerOnlineEMCALEJ1(false),
fTriggerOnlineEMCALEJ2(false),
fTriggerOnlineDCALDJ1(false),
fTriggerOnlineDCALDJ2(false),
fTriggerBitINT7(false),
fTriggerBitHighMultSPD(false),
fTriggerBitHighMultV0(false),
fTriggerBitCentral(false),
fTriggerBitSemiCentral(false),
fTriggerBitEMCEJE(false),
fTriggerClasses(""),
fTriggerClassINT7(false),
fTriggerClassHighMultSPD(false),
fTriggerClassHighMultV0m(false),
fTriggerClassEMCALEJ1(false),
fTriggerClassEMCALEJ2(false),
fTriggerClassDCALDJ1(false),
fTriggerClassDCALDJ2(false),
fnV0M(0),
fnV0MEq(0),
fnV0MCorr(0),
fnV0MEqCorr(0),
fPercV0M(0.),
fMultV0M(0.),
fFillMCGenTrees(kTRUE),
fDsMassKKOpt(1),
fLc2V0bachelorCalcSecoVtx(0),
fV0typeForLc2V0bachelor(1),
fInvMassOnFlyCut(0.3),
fPtOnFlyCut(-1),
fImpParProdOnFlyCut(9999.),
fCosPOnFlyCut(-9999.),
fCosPXYOnFlyCut(-9999.),
fTreeSingleTrackVarsOpt(AliHFTreeHandler::kRedSingleTrackVars),
fJetRadius(0.4),
fSubJetRadius(0.0),
fJetAlgorithm(JetAlgorithm::antikt),
fSubJetAlgorithm(JetAlgorithm::ca),
fMinJetPt(0.0),
fSoftDropZCut(0.1),
fSoftDropBeta(0.0),
fTrackingEfficiency(1.0),
fGoodTrackFilterBit(-1),
fGoodTrackEtaRange(999.),
fGoodTrackMinPt(0.),
fITSUpgradeProduction(0),
fITSUpgradePreSelect(0),
fStoreOnlyHIJINGBackground(0),
fPreSelectLctopKpi(false),
fFillInjCandHijingTrackCombi(false),
fFillParticleTree(false),
fFillTrackletTree(false),
fVariablesTreeParticle(0),
fVariablesTreeTracklet(0),
fVariablesTreeGenParticle(0),
fTreeHandlerParticle(nullptr),
fTreeHandlerTracklet(nullptr),
fTreeHandlerGenParticle(nullptr),
fParticleCollArray(),
fWriteNJetTrees(fillNJetTrees),
fFillJetConstituentTrees(fillJetConstituentTrees),
fVariablesTreeJet(0),
fVariablesTreeJetConstituent(0),
fTreeHandlerJet(0),
fLocalInitialized(kFALSE),
fJetCollArray(),
fMinJetPtCorr(0.),
fRhoName(),
fRho(0),
fRhoVal(0),
fFillJetEtaPhi(false),
fFillPtCorr(false),
fFillPtUncorr(false),
fFillArea(false),
fFillNConstituents(false),
fFillZLeading(false),
fFillRadialMoment(false),
fFillpTD(false),
fFillMass(false),
fFillMatchingJetID(false),
fFillJets(false),
fDoJetSubstructure(false),
fDoPtHard(false),
fEnableNsigmaTPCDataCorr(false),
fSystemForNsigmaTPCDataCorr(AliAODPidHF::kNone),
fCorrNtrVtx(false),
fCorrV0MVtx(false),
fMultEstimatorAvg(),
fMultEstimatorAvgSHM(),
fApplyPhysicsSelOnline(false),
fEnableEventDownsampling(false),
fFracToKeepEventDownsampling(1.1),
fSeedEventDownsampling(0),
fCdbEntry(nullptr)
{
  fParticleCollArray.SetOwner(kTRUE);
  fJetCollArray.SetOwner(kTRUE);
  
  if (fListCuts) {
    fFiltCutsD0toKpi      =(AliRDHFCutsD0toKpi*)fListCuts->FindObject("D0toKpiFilteringCuts");
    fFiltCutsDstoKKpi     =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("DstoKKpiFilteringCuts");
    fFiltCutsDplustoKpipi =(AliRDHFCutsDplustoKpipi*)fListCuts->FindObject("DplustoKpipiFilteringCuts");
    fFiltCutsLctopKpi     =(AliRDHFCutsLctopKpi*)fListCuts->FindObject("LctopKpiFilteringCuts");
    fFiltCutsBplustoD0pi  =(AliRDHFCutsD0toKpi*)fListCuts->FindObject("BplustoD0piFilteringCuts");
    fFiltCutsBstoDspi     =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("BstoDspiFilteringCuts");
    fFiltCutsDstartoKpipi =(AliRDHFCutsDStartoKpipi*)fListCuts->FindObject("DstartoKpipiFilteringCuts");
    fFiltCutsLc2V0bachelor=(AliRDHFCutsLctoV0*)fListCuts->FindObject("Lc2V0bachelorFilteringCuts");
    fFiltCutsLbtoLcpi     =(AliRDHFCutsLctopKpi*)fListCuts->FindObject("LbtoLcpiFilteringCuts");
    fCutsD0toKpi          =(AliRDHFCutsD0toKpi*)fListCuts->FindObject("D0toKpiAnalysisCuts");
    fCutsDstoKKpi         =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("DstoKKpiAnalysisCuts");
    fCutsDplustoKpipi     =(AliRDHFCutsDplustoKpipi*)fListCuts->FindObject("DplustoKpipiAnalysisCuts");
    fCutsLctopKpi         =(AliRDHFCutsLctopKpi*)fListCuts->FindObject("LctopKpiAnalysisCuts");
    fCutsBplustoD0pi      =(AliRDHFCutsD0toKpi*)fListCuts->FindObject("BplustoD0piAnalysisCuts");
    fCutsBstoDspi         =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("BstoDspiAnalysisCuts");
    fCutsDstartoKpipi     =(AliRDHFCutsDStartoKpipi*)fListCuts->FindObject("DstartoKpipiAnalysisCuts");
    fCutsLc2V0bachelor    =(AliRDHFCutsLctoV0*)fListCuts->FindObject("Lc2V0bachelorAnalysisCuts");
    fCutsLbtoLcpi         =(AliRDHFCutsLctopKpi*)fListCuts->FindObject("LbtoLcpiAnalysisCuts");
    
    if(fWriteVariableTreeD0) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsD0toKpi->Clone();
    else if(fWriteVariableTreeDplus) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDplustoKpipi->Clone();
    else if(fWriteVariableTreeDs) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDstoKKpi->Clone();
    else if(fWriteVariableTreeDstar) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDstartoKpipi->Clone();
    else if(fWriteVariableTreeBplus) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsBplustoD0pi->Clone();
    else if(fWriteVariableTreeBs) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsBstoDspi->Clone();
    else if(fWriteVariableTreeLctopKpi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLctopKpi->Clone();
    else if(fWriteVariableTreeLc2V0bachelor) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLc2V0bachelor->Clone();
    else if(fWriteVariableTreeLb) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLbtoLcpi->Clone();
    else {
      if(fFiltCutsD0toKpi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsD0toKpi->Clone();
      else if(fFiltCutsDstoKKpi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDstoKKpi->Clone();
      else if(fFiltCutsDplustoKpipi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDplustoKpipi->Clone();
      else if(fFiltCutsDstartoKpipi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDstartoKpipi->Clone();
      else if(fFiltCutsBplustoD0pi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsBplustoD0pi->Clone();
      else if(fFiltCutsBstoDspi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsBstoDspi->Clone();
      else if(fFiltCutsLctopKpi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLctopKpi->Clone();
      else if(fFiltCutsLc2V0bachelor) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLc2V0bachelor->Clone();
      else if(fFiltCutsLbtoLcpi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLbtoLcpi->Clone();
      else printf("AliAnalysisTaskSEHFTreeCreator:: Constructor: No event selection cuts could be stored, code will crash!\n");
    }
  }
  
  DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TH1F container (number of events)
  DefineOutput(1,TH1F::Class());
  // Output slot #2 writes into a TH2F container (number of events)
  DefineOutput(2,TH2F::Class());
  // Output slot #3 writes into a TList container (cuts)
  DefineOutput(3,TList::Class());
  // Output slot #4 writes Normalization Counter
  DefineOutput(4,TList::Class());
  // Output slot #5 stores the tree of the event-characterisation variables
  DefineOutput(5,TTree::Class());
  // Output slot #6 stores the tree of the D0 candidate variables after track selection
  DefineOutput(6,TTree::Class());
  // Output slot #7 stores the tree of the gen D0 variables
  DefineOutput(7,TTree::Class());
  // Output slot #8 stores the tree of the D+ candidate variables after track selection
  DefineOutput(8,TTree::Class());
  // Output slot #9 stores the tree of the gen D0 variables
  DefineOutput(9,TTree::Class());
  // Output slot #10 stores the tree of the Ds+ candidate variables after track selection
  DefineOutput(10,TTree::Class());
  // Output slot #11 stores the tree of the gen Ds+ variables
  DefineOutput(11,TTree::Class());
  // Output slot #12 stores the tree of the LctopKpi candidate variables after track selection
  DefineOutput(12,TTree::Class());
  // Output slot #13 stores the tree of the gen LctopKpi variables
  DefineOutput(13,TTree::Class());
  // Output slot #14 stores the tree of the B+ candidate variables after track selection
  DefineOutput(14,TTree::Class());
  // Output slot #15 stores the tree of the gen B+ variables
  DefineOutput(15,TTree::Class());
  // Output slot #16 stores the tree of the Dstar candidate variables after track selection
  DefineOutput(16,TTree::Class());
  // Output slot #17 stores the tree of the gen Dstar variables
  DefineOutput(17,TTree::Class());
  // Output slot #18 stores the tree of the Lc2V0bachelor candidate variables after track selection
  DefineOutput(18,TTree::Class());
  // Output slot #19 stores the tree of the gen Lc2V0bachelor variables
  DefineOutput(19,TTree::Class());
  // Output slot #20 stores the tree of the Bs candidate variables after track selection
  DefineOutput(20,TTree::Class());
  // Output slot #21 stores the tree of the gen Bs variables
  DefineOutput(21,TTree::Class());
  // Output slot #22 stores the tree of the Lb candidate variables after track selection
  DefineOutput(22,TTree::Class());
  // Output slot #23 stores the tree of the gen Lb variables
  DefineOutput(23,TTree::Class());
  // Output slot #24 stores the tree of the inclusive jet variables 
  DefineOutput(24,TTree::Class());
  // Output slot #25 stores the tree of the gen inclusive jet variables
  DefineOutput(25,TTree::Class()); 
  // Output slot #26 stores the tree of the track variables after track selection
  DefineOutput(26,TTree::Class());
  // Output slot #27 stores the tree of the MC particle variables
  DefineOutput(27,TTree::Class());
  // Output slot #28 stores tracklets
  DefineOutput(28,TTree::Class());
  
  // Set up separate output slot for each jet tree
  // (for simplicity, keep the jet trees in the last slots)
  for (int i=0; i<fillNJetTrees; i++) {
    // Output slot #29 stores the tree of the jet variables
    DefineOutput(29+i,TTree::Class());
  }
  
  // Set up separate output slot for each jet constituent tree (if enabled)
  if (fillJetConstituentTrees) {
    for (int i=0; i<fillNJetTrees; i++) {
      // Output slot #29 stores the tree of the jet variables
      DefineOutput(29+fillNJetTrees+i,TTree::Class());
    }
  }
  
}

//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreator::~AliAnalysisTaskSEHFTreeCreator()
{
  delete fListCuts;
  delete fFiltCutsD0toKpi;
  delete fFiltCutsDstoKKpi;
  delete fFiltCutsDplustoKpipi;
  delete fFiltCutsLctopKpi;
  delete fFiltCutsBplustoD0pi;
  delete fFiltCutsBstoDspi;
  delete fFiltCutsDstartoKpipi;
  delete fFiltCutsLc2V0bachelor;
  delete fFiltCutsLbtoLcpi;
  delete fCutsD0toKpi;
  delete fCutsDstoKKpi;
  delete fCutsDplustoKpipi;
  delete fCutsLctopKpi;
  delete fCutsBplustoD0pi;
  delete fCutsBstoDspi;
  delete fCutsDstartoKpipi;
  delete fCutsLc2V0bachelor;
  delete fCutsLbtoLcpi;
  delete fEvSelectionCuts;
  delete fNentries;
  delete fHistoNormCounter;
  delete fListCounter;
  delete fCounter;
  delete fTreeHandlerD0;
  delete fTreeHandlerDs;
  delete fTreeHandlerDplus;
  delete fTreeHandlerLctopKpi;
  delete fTreeHandlerBplus;
  delete fTreeHandlerBs;
  delete fTreeHandlerDstar;
  delete fTreeHandlerLc2V0bachelor;
  delete fTreeHandlerLb;
  delete fTreeHandlerParticle;
  delete fTreeHandlerTracklet;
  
  for(auto& thj : fTreeHandlerJet)
    delete thj;
  delete fTreeHandlerGenD0;
  delete fTreeHandlerGenDs;
  delete fTreeHandlerGenDplus;
  delete fTreeHandlerGenLctopKpi;
  delete fTreeHandlerGenBplus;
  delete fTreeHandlerGenBs;
  delete fTreeHandlerGenDstar;
  delete fTreeHandlerGenLc2V0bachelor;
  delete fTreeHandlerGenLb;
  delete fTreeHandlerGenParticle;
  delete fTreeEvChar;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::Init()
{
  /// Initialization
  
  if(fDebug > 1) printf("AliAnalysisTaskSEHFTreeCreator::Init() \n");
  
  PostData(3,fListCuts);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::UserCreateOutputObjects()
{
  
  /// Create the output container
  //
  if(fDebug > 1) printf("AliAnalysisTaskSEHFTreeCreator::UserCreateOutputObjects() \n");
  
  const char* nameoutput=GetOutputSlot(1)->GetContainer()->GetName();
  fNentries=new TH1F(nameoutput, "Number of events", 44,-0.5,43.5);
  fNentries->GetXaxis()->SetBinLabel(1,"n. evt. read");
  fNentries->GetXaxis()->SetBinLabel(2,"n. evt. matched dAOD");
  fNentries->GetXaxis()->SetBinLabel(3,"n. evt. mismatched dAOD");
  fNentries->GetXaxis()->SetBinLabel(4,"n. evt. analised");
  fNentries->GetXaxis()->SetBinLabel(5,"n. evt. passing IsEvSelected (including pileup)");
  fNentries->GetXaxis()->SetBinLabel(6,"n. evt. rejected due to trigger");
  fNentries->GetXaxis()->SetBinLabel(7,"n. evt. rejected due to not reco vertex");
  fNentries->GetXaxis()->SetBinLabel(8,"n. evt. rejected for contr vertex");
  fNentries->GetXaxis()->SetBinLabel(9,"n. evt. rejected for vertex out of accept");
  fNentries->GetXaxis()->SetBinLabel(10,"n. evt. rejected for pileup events");
  fNentries->GetXaxis()->SetBinLabel(11,"n. evt. of out centrality events");
  fNentries->GetXaxis()->SetBinLabel(12,"n. of 2 prong candidates");
  fNentries->GetXaxis()->SetBinLabel(13,"n. D0 after filtering");
  fNentries->GetXaxis()->SetBinLabel(14,"n. D0 after selection");
  fNentries->GetXaxis()->SetBinLabel(15,"n. of not on-the-fly rec D0");
  fNentries->GetXaxis()->SetBinLabel(16,"n. of 3 prong candidates");
  fNentries->GetXaxis()->SetBinLabel(17,"n. Ds after filtering");
  fNentries->GetXaxis()->SetBinLabel(18,"n. Ds after selection");
  fNentries->GetXaxis()->SetBinLabel(19,"n. of not on-the-fly rec Ds");
  fNentries->GetXaxis()->SetBinLabel(20,"n. Dplus after filtering");
  fNentries->GetXaxis()->SetBinLabel(21,"n. Dplus after selection");
  fNentries->GetXaxis()->SetBinLabel(22,"n. of not on-the-fly rec Dplus");
  fNentries->GetXaxis()->SetBinLabel(23,"n. LctopKpi after filtering");
  fNentries->GetXaxis()->SetBinLabel(24,"n. LctopKpi after selection");
  fNentries->GetXaxis()->SetBinLabel(25,"n. of not on-the-fly rec LctopKpi");
  fNentries->GetXaxis()->SetBinLabel(26,"n. Bplus D0's after filtering");
  fNentries->GetXaxis()->SetBinLabel(27,"n. Bplus after selection");
  fNentries->GetXaxis()->SetBinLabel(28,"n. of not on-the-fly rec Bplus");
  fNentries->GetXaxis()->SetBinLabel(29,"n. of Dstar candidates");
  fNentries->GetXaxis()->SetBinLabel(30,"n. Dstar after filtering");
  fNentries->GetXaxis()->SetBinLabel(31,"n. Dstar after selection");
  fNentries->GetXaxis()->SetBinLabel(32,"n. of not on-the-fly rec Dstar");
  fNentries->GetXaxis()->SetBinLabel(33,"n. of cascade candidates");
  fNentries->GetXaxis()->SetBinLabel(34,"n. Lc2V0bachelor after filtering");
  fNentries->GetXaxis()->SetBinLabel(35,"n. Lc2V0bachelor after selection");
  fNentries->GetXaxis()->SetBinLabel(36,"n. of not on-the-fly rec Lc2V0bachelor");
  fNentries->GetXaxis()->SetBinLabel(37,"n. Bs Ds's after filtering");
  fNentries->GetXaxis()->SetBinLabel(38,"n. Bs after selection");
  fNentries->GetXaxis()->SetBinLabel(39,"n. of not on-the-fly rec Bs");
  fNentries->GetXaxis()->SetBinLabel(40,"n. Lb Lc's after filtering");
  fNentries->GetXaxis()->SetBinLabel(41,"n. Lb after selection");
  fNentries->GetXaxis()->SetBinLabel(42,"n. of not on-the-fly rec Lb");
  
  nameoutput=GetOutputSlot(2)->GetContainer()->GetName();
  fHistoNormCounter=new TH2F(nameoutput, "Number of events for norm;;centrality", 5,-0.5,4.5,102,-1.,101.);
  fHistoNormCounter->GetXaxis()->SetBinLabel(1,"n. evt. w primary V");
  fHistoNormCounter->GetXaxis()->SetBinLabel(2,"n. evt. w/o primary V");
  fHistoNormCounter->GetXaxis()->SetBinLabel(3,"n. evt. w primary V z>10cm");
  fHistoNormCounter->GetXaxis()->SetBinLabel(4,"n. evt. for norm.");
  fHistoNormCounter->GetXaxis()->SetBinLabel(5,"n. evt. pileup");
  
  
  fListCounter=new TList();
  fListCounter->SetOwner(kTRUE);
  fListCounter->SetName("NormCounter");
  fCounter = new AliNormalizationCounter("norm_counter");
  fCounter->SetStudyMultiplicity(kTRUE,1.);
  fCounter->Init();
  fListCounter->Add(fCounter);
  
  //count number of enabled trees
  Int_t nEnabledTrees = 1; // event tree always enabled
  if(fWriteVariableTreeD0) nEnabledTrees++;
  if(fWriteVariableTreeDs) nEnabledTrees++;
  if(fWriteVariableTreeDplus) nEnabledTrees++;
  if(fWriteVariableTreeLctopKpi) nEnabledTrees++;
  if(fWriteVariableTreeBplus) nEnabledTrees++;
  if(fWriteVariableTreeBs) nEnabledTrees++;
  if(fWriteVariableTreeDstar) nEnabledTrees++;
  if(fWriteVariableTreeLc2V0bachelor) nEnabledTrees++;
  if(fWriteVariableTreeLb) nEnabledTrees++;
  if (fFillParticleTree) nEnabledTrees++;
  if (fFillTrackletTree) nEnabledTrees++;
  if(fReadMC && fFillMCGenTrees) {
    nEnabledTrees = (nEnabledTrees-1)*2+1;
  }
  nEnabledTrees += fWriteNJetTrees;
  if (fFillJetConstituentTrees) {
    nEnabledTrees += fWriteNJetTrees;
  }
  
  
  //
  // Output slot 4-25 : trees of the candidate and event-characterization variables
  //
  OpenFile(5);
  fTreeEvChar = new TTree("tree_event_char","tree_event_char");
  //set variables
  fTreeEvChar->Branch("centrality", &fCentrality);
  fTreeEvChar->Branch("z_vtx_reco", &fzVtxReco);
  fTreeEvChar->Branch("n_vtx_contributors", &fNcontributors);
  fTreeEvChar->Branch("n_tracks", &fNtracks);
  fTreeEvChar->Branch("is_ev_rej", &fIsEvRej);
  fTreeEvChar->Branch("is_ev_rej_INT7", &fIsEvRej_INT7);
  fTreeEvChar->Branch("is_ev_rej_HighMultSPD", &fIsEvRej_HighMultSPD);
  fTreeEvChar->Branch("is_ev_rej_HighMultV0", &fIsEvRej_HighMultV0);
  fTreeEvChar->Branch("is_ev_rej_EMCEJE", &fIsEvRej_EMCEJE);
  fTreeEvChar->Branch("run_number", &fRunNumber);
  fTreeEvChar->Branch("ev_id", &fEventID);
  fTreeEvChar->Branch("ev_id_ext", &fEventIDExt);
  fTreeEvChar->Branch("ev_id_long", &fEventIDLong);
  fTreeEvChar->Branch("n_tracklets", &fnTracklets);
  fTreeEvChar->Branch("V0Amult", &fnV0A);
  fTreeEvChar->Branch("trigger_bitmap", &fTriggerMask);
  fTreeEvChar->Branch("trigger_online_INT7", &fTriggerOnlineINT7);
  fTreeEvChar->Branch("trigger_online_HighMultSPD", &fTriggerOnlineHighMultSPD);
  fTreeEvChar->Branch("trigger_online_HighMultV0", &fTriggerOnlineHighMultV0);
  fTreeEvChar->Branch("trigger_online_EMCALEJ1", &fTriggerOnlineEMCALEJ1);
  fTreeEvChar->Branch("trigger_online_EMCALEJ2", &fTriggerOnlineEMCALEJ2);
  fTreeEvChar->Branch("trigger_online_DCALDJ1", &fTriggerOnlineDCALDJ1);
  fTreeEvChar->Branch("trigger_online_DCALDJ2", &fTriggerOnlineDCALDJ2);
  fTreeEvChar->Branch("trigger_hasbit_INT7", &fTriggerBitINT7);
  fTreeEvChar->Branch("trigger_hasbit_HighMultSPD", &fTriggerBitHighMultSPD);
  fTreeEvChar->Branch("trigger_hasbit_HighMultV0", &fTriggerBitHighMultV0);
  fTreeEvChar->Branch("trigger_hasbit_Central", &fTriggerBitCentral);
  fTreeEvChar->Branch("trigger_hasbit_SemiCentral", &fTriggerBitSemiCentral);
  fTreeEvChar->Branch("trigger_hasbit_EMCEJE", &fTriggerBitEMCEJE);
  fTreeEvChar->Branch("trigger_classes", &fTriggerClasses);
  fTreeEvChar->Branch("trigger_hasclass_INT7", &fTriggerClassINT7);
  fTreeEvChar->Branch("trigger_hasclass_HighMultSPD", &fTriggerClassHighMultSPD);
  fTreeEvChar->Branch("trigger_hasclass_HighMultV0", &fTriggerClassHighMultV0m);
  fTreeEvChar->Branch("trigger_hasclass_EMCALEJ1", &fTriggerClassEMCALEJ1);
  fTreeEvChar->Branch("trigger_hasclass_EMCALEJ2", &fTriggerClassEMCALEJ2);
  fTreeEvChar->Branch("trigger_hasclass_DCALDJ1", &fTriggerClassDCALDJ1);
  fTreeEvChar->Branch("trigger_hasclass_DCALDJ2", &fTriggerClassDCALDJ2);
  fTreeEvChar->Branch("z_vtx_gen", &fzVtxGen);
  fTreeEvChar->Branch("n_tracklets_corr", &fnTrackletsCorr);
  fTreeEvChar->Branch("n_tracklets_corr_shm", &fnTrackletsCorrSHM);
  fTreeEvChar->Branch("v0m", &fnV0M);
  fTreeEvChar->Branch("v0m_eq", &fnV0MEq);
  fTreeEvChar->Branch("v0m_corr", &fnV0MCorr);
  fTreeEvChar->Branch("v0m_eq_corr", &fnV0MEqCorr);
  fTreeEvChar->Branch("mult_gen", &fMultGen);
  fTreeEvChar->Branch("mult_gen_v0a", &fMultGenV0A);
  fTreeEvChar->Branch("mult_gen_v0c", &fMultGenV0C);
  fTreeEvChar->Branch("perc_v0m", &fPercV0M);
  fTreeEvChar->Branch("mult_v0m", &fMultV0M);
  fTreeEvChar->Branch("is_ev_sel_int7", &fIsEvSel_INT7);
  fTreeEvChar->Branch("is_ev_sel_shm", &fIsEvSel_HighMultSPD);
  fTreeEvChar->Branch("is_ev_sel_vhm", &fIsEvSel_HighMultV0);
  fTreeEvChar->Branch("is_ev_sel_EMCEJE", &fIsEvSel_EMCEJE);
  fTreeEvChar->Branch("cross_section", &fCross_Section);
  fTreeEvChar->Branch("trials", &fTrials);
  fTreeEvChar->Branch("pthard", &fpthard);
  fTreeEvChar->SetMaxVirtualSize(1.e+8/nEnabledTrees);
  
  if(fWriteVariableTreeD0){
    OpenFile(6);
    TString nameoutput = "tree_D0";
    fTreeHandlerD0 = new AliHFTreeHandlerD0toKpi(fPIDoptD0);
    fTreeHandlerD0->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerD0->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerD0->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerD0->SetFillJets(fFillJets);
    fTreeHandlerD0->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerD0->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerD0->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerD0->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeD0 = (TTree*)fTreeHandlerD0->BuildTree(nameoutput,nameoutput);
    fVariablesTreeD0->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeD0);
    
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(7);
      TString nameoutput = "tree_D0_gen";
      fTreeHandlerGenD0 = new AliHFTreeHandlerD0toKpi(0);
      fTreeHandlerGenD0->SetFillJets(fFillJets);
      fTreeHandlerGenD0->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenD0->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenD0->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeD0 = (TTree*)fTreeHandlerGenD0->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeD0->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeD0);
    }
  }
  if(fWriteVariableTreeDs){
    OpenFile(8);
    TString nameoutput = "tree_Ds";
    fTreeHandlerDs = new AliHFTreeHandlerDstoKKpi(fPIDoptDs);
    fTreeHandlerDs->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerDs->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerDs->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerDs->SetMassKKOption(fDsMassKKOpt);
    fTreeHandlerDs->SetFillJets(fFillJets);
    fTreeHandlerDs->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerDs->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerDs->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerDs->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeDs = (TTree*)fTreeHandlerDs->BuildTree(nameoutput,nameoutput);
    fVariablesTreeDs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeDs);
    
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(9);
      TString nameoutput = "tree_Ds_gen";
      fTreeHandlerGenDs = new AliHFTreeHandlerDstoKKpi(0);
      fTreeHandlerGenDs->SetFillJets(fFillJets);
      fTreeHandlerGenDs->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenDs->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenDs->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeDs = (TTree*)fTreeHandlerGenDs->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeDs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeDs);
    }
  }
  if(fWriteVariableTreeDplus){
    OpenFile(10);
    TString nameoutput = "tree_Dplus";
    fTreeHandlerDplus = new AliHFTreeHandlerDplustoKpipi(fPIDoptDplus);
    fTreeHandlerDplus->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerDplus->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerDplus->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerDplus->SetFillJets(fFillJets);
    fTreeHandlerDplus->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerDplus->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerDplus->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerDplus->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeDplus = (TTree*)fTreeHandlerDplus->BuildTree(nameoutput,nameoutput);
    fVariablesTreeDplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeDplus);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(11);
      TString nameoutput = "tree_Dplus_gen";
      fTreeHandlerGenDplus = new AliHFTreeHandlerDplustoKpipi(0);
      fTreeHandlerGenDplus->SetFillJets(fFillJets);
      fTreeHandlerGenDplus->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenDplus->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenDplus->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeDplus = (TTree*)fTreeHandlerGenDplus->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeDplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeDplus);
    }
  }
  if(fWriteVariableTreeLctopKpi){
    OpenFile(12);
    TString nameoutput = "tree_LctopKpi";
    fTreeHandlerLctopKpi = new AliHFTreeHandlerLctopKpi(fPIDoptLctopKpi);
    fTreeHandlerLctopKpi->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerLctopKpi->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerLctopKpi->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerLctopKpi->SetFillJets(fFillJets);
    fTreeHandlerLctopKpi->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerLctopKpi->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerLctopKpi->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerLctopKpi->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeLctopKpi = (TTree*)fTreeHandlerLctopKpi->BuildTree(nameoutput,nameoutput);
    fVariablesTreeLctopKpi->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeLctopKpi);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(13);
      TString nameoutput = "tree_LctopKpi_gen";
      fTreeHandlerGenLctopKpi = new AliHFTreeHandlerLctopKpi(0);
      fTreeHandlerGenLctopKpi->SetFillJets(fFillJets);
      fTreeHandlerGenLctopKpi->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenLctopKpi->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenLctopKpi->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeLctopKpi = (TTree*)fTreeHandlerGenLctopKpi->BuildTreeMCGen(nameoutput,nameoutput);
      fTreeHandlerGenLctopKpi->AddBranchResonantDecay(fGenTreeLctopKpi);
      fGenTreeLctopKpi->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeLctopKpi);
    }
  }
  if(fWriteVariableTreeBplus){
    OpenFile(14);
    TString nameoutput = "tree_Bplus";
    fTreeHandlerBplus = new AliHFTreeHandlerBplustoD0pi(fPIDoptBplus);
    fTreeHandlerBplus->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerBplus->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerBplus->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerBplus->SetFillJets(fFillJets);
    fTreeHandlerBplus->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerBplus->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerBplus->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerBplus->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeBplus = (TTree*)fTreeHandlerBplus->BuildTree(nameoutput,nameoutput);
    fVariablesTreeBplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeBplus);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(15);
      TString nameoutput = "tree_Bplus_gen";
      fTreeHandlerGenBplus = new AliHFTreeHandlerBplustoD0pi(0);
      fTreeHandlerGenBplus->SetFillJets(fFillJets);
      fTreeHandlerGenBplus->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenBplus->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenBplus->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeBplus = (TTree*)fTreeHandlerGenBplus->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeBplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeBplus);
    }
  }
  if(fWriteVariableTreeDstar){
    OpenFile(16);
    TString nameoutput = "tree_Dstar";
    fTreeHandlerDstar = new AliHFTreeHandlerDstartoKpipi(fPIDoptDstar);
    fTreeHandlerDstar->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerDstar->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerDstar->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerDstar->SetFillJets(fFillJets);
    fTreeHandlerDstar->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerDstar->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerDstar->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerDstar->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeDstar = (TTree*)fTreeHandlerDstar->BuildTree(nameoutput,nameoutput);
    fVariablesTreeDstar->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeDstar);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(17);
      TString nameoutput = "tree_Dstar_gen";
      fTreeHandlerGenDstar = new AliHFTreeHandlerDstartoKpipi(0);
      fTreeHandlerGenDstar->SetFillJets(fFillJets);
      fTreeHandlerGenDstar->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenDstar->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenDstar->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeDstar = (TTree*)fTreeHandlerGenDstar->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeDstar->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeDstar);
    }
  }
  if(fWriteVariableTreeLc2V0bachelor){
    OpenFile(18);
    TString nameoutput = "tree_Lc2V0bachelor";
    fTreeHandlerLc2V0bachelor = new AliHFTreeHandlerLc2V0bachelor(fPIDoptLc2V0bachelor);
    fTreeHandlerLc2V0bachelor->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerLc2V0bachelor->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerLc2V0bachelor->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerLc2V0bachelor->SetCalcSecoVtx(fLc2V0bachelorCalcSecoVtx);
    fTreeHandlerLc2V0bachelor->SetFillJets(fFillJets);
    fTreeHandlerLc2V0bachelor->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerLc2V0bachelor->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerLc2V0bachelor->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerLc2V0bachelor->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeLc2V0bachelor = (TTree*)fTreeHandlerLc2V0bachelor->BuildTree(nameoutput,nameoutput);
    fVariablesTreeLc2V0bachelor->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeLc2V0bachelor);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(19);
      TString nameoutput = "tree_Lc2V0bachelor_gen";
      fTreeHandlerGenLc2V0bachelor = new AliHFTreeHandlerLc2V0bachelor(0);
      fTreeHandlerGenLc2V0bachelor->SetFillJets(fFillJets);
      fTreeHandlerGenLc2V0bachelor->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenLc2V0bachelor->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenLc2V0bachelor->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeLc2V0bachelor = (TTree*)fTreeHandlerGenLc2V0bachelor->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeLc2V0bachelor->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeLc2V0bachelor);
    }
  }
  if(fWriteVariableTreeBs){
    OpenFile(20);
    TString nameoutput = "tree_Bs";
    fTreeHandlerBs = new AliHFTreeHandlerBstoDspi(fPIDoptBs);
    fTreeHandlerBs->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerBs->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerBs->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerBs->SetBsSelectionValues(fInvMassOnFlyCut,fPtOnFlyCut,fImpParProdOnFlyCut,fCosPOnFlyCut,fCosPXYOnFlyCut);
    fTreeHandlerBs->SetFillJets(fFillJets);
    fTreeHandlerBs->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerBs->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerBs->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerBs->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeBs = (TTree*)fTreeHandlerBs->BuildTree(nameoutput,nameoutput);
    fVariablesTreeBs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeBs);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(21);
      TString nameoutput = "tree_Bs_gen";
      fTreeHandlerGenBs = new AliHFTreeHandlerBstoDspi(0);
      fTreeHandlerGenBs->SetFillJets(fFillJets);
      fTreeHandlerGenBs->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenBs->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenBs->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeBs = (TTree*)fTreeHandlerGenBs->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeBs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeBs);
    }
  }
  if(fWriteVariableTreeLb){
    OpenFile(22);
    TString nameoutput = "tree_Lb";
    fTreeHandlerLb = new AliHFTreeHandlerLbtoLcpi(fPIDoptLb);
    fTreeHandlerLb->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerLb->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerLb->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerLb->SetFillJets(fFillJets);
    fTreeHandlerLb->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerLb->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerLb->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerLb->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeLb = (TTree*)fTreeHandlerLb->BuildTree(nameoutput,nameoutput);
    fVariablesTreeLb->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeLb);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(23);
      TString nameoutput = "tree_Lb_gen";
      fTreeHandlerGenLb = new AliHFTreeHandlerLbtoLcpi(0);
      fTreeHandlerGenLb->SetFillJets(fFillJets);
      fTreeHandlerGenLb->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenLb->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenLb->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeLb = (TTree*)fTreeHandlerGenLb->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeLb->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeLb);
    }
  }


  if(fWriteVariableTreeInclusiveJet){
    OpenFile(24);
    TString nameoutput = "tree_InclusiveJet";
    fTreeHandlerInclusiveJet = new AliHFTreeHandlerInclusiveJet();
    fTreeHandlerInclusiveJet->SetFillJets(fFillJets);
    fTreeHandlerInclusiveJet->SetDoJetSubstructure(fDoJetSubstructure);
    fTreeHandlerInclusiveJet->SetTrackingEfficiency(fTrackingEfficiency);
    fTreeHandlerInclusiveJet->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
    fTreeHandlerInclusiveJet->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
    fVariablesTreeInclusiveJet = (TTree*)fTreeHandlerInclusiveJet->BuildTree(nameoutput,nameoutput);
    fVariablesTreeInclusiveJet->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeInclusiveJet);
    
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(25);
      TString nameoutput = "tree_InclusiveJet_gen";
      fTreeHandlerGenInclusiveJet = new AliHFTreeHandlerInclusiveJet();
      fTreeHandlerGenInclusiveJet->SetFillJets(fFillJets);
      fTreeHandlerGenInclusiveJet->SetDoJetSubstructure(fDoJetSubstructure);
      fTreeHandlerGenInclusiveJet->SetJetProperties(fJetRadius,fJetAlgorithm,fMinJetPt);
      fTreeHandlerGenInclusiveJet->SetSubJetProperties(fSubJetRadius,fSubJetAlgorithm,fSoftDropZCut,fSoftDropBeta);
      fGenTreeInclusiveJet = (TTree*)fTreeHandlerGenInclusiveJet->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeInclusiveJet->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeInclusiveJet);
    }
  }

  
  if(fFillParticleTree){
    OpenFile(26);
    TString nameoutput = "tree_Particle";
    fTreeHandlerParticle = new AliParticleTreeHandler();
    fTreeHandlerParticle->SetParticleContainer(GetParticleContainer(0));
    fVariablesTreeParticle = (TTree*)fTreeHandlerParticle->BuildTree(nameoutput,nameoutput);
    fVariablesTreeParticle->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeParticle);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(27);
      TString nameoutput = "tree_Particle_gen";
      fTreeHandlerGenParticle = new AliParticleTreeHandler();
      fTreeHandlerGenParticle->SetParticleContainer(GetParticleContainer(1));
      fVariablesTreeGenParticle = (TTree*)fTreeHandlerGenParticle->BuildTree(nameoutput,nameoutput);
      fVariablesTreeGenParticle->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fVariablesTreeGenParticle);
    }
  }
  if(fFillTrackletTree){
    OpenFile(28);
    TString nameoutput = "tree_Tracklet";
    fTreeHandlerTracklet = new AliTrackletTreeHandler();
    fVariablesTreeTracklet = (TTree*)fTreeHandlerTracklet->BuildTree(nameoutput,nameoutput);
    fVariablesTreeTracklet->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeTracklet); 
  }
  if(fWriteNJetTrees > 0){
    for (int i=0; i<fJetCollArray.GetEntriesFast(); i++) {
      OpenFile(29 + i);
      
      // Create jet tree handlers and configure them
      fTreeHandlerJet.push_back(new AliJetTreeHandler());
      fTreeHandlerJet.at(i)->SetFillJetConstituentTree(fFillJetConstituentTrees);
      fTreeHandlerJet.at(i)->SetJetContainer(GetJetContainer(i));
      fTreeHandlerJet.at(i)->SetMinJetPtCorr(fMinJetPtCorr);
      
      fTreeHandlerJet.at(i)->SetFillJetEtaPhi(fFillJetEtaPhi);
      fTreeHandlerJet.at(i)->SetFillPtCorr(fFillPtCorr);
      fTreeHandlerJet.at(i)->SetFillPtUncorr(fFillPtUncorr);
      fTreeHandlerJet.at(i)->SetFillArea(fFillArea);
      fTreeHandlerJet.at(i)->SetFillNConstituents(fFillNConstituents);
      fTreeHandlerJet.at(i)->SetFillZLeading(fFillZLeading);
      fTreeHandlerJet.at(i)->SetFillRadialMoment(fFillRadialMoment);
      fTreeHandlerJet.at(i)->SetFillpTD(fFillpTD);
      fTreeHandlerJet.at(i)->SetFillMass(fFillMass);
      fTreeHandlerJet.at(i)->SetFillMatchingJetID(fFillMatchingJetID);
      
      // Build jet trees
      TString nameoutput = GetJetContainer(i)->GetName();
      fVariablesTreeJet.push_back((TTree*)fTreeHandlerJet.at(i)->BuildJetTree(nameoutput,nameoutput));
      fVariablesTreeJet.at(i)->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fVariablesTreeJet.at(i));
      
      // Build jet constituent trees (if enabled)
      if (fFillJetConstituentTrees) {
        OpenFile(29 + fWriteNJetTrees + i);
        TString nameoutput = Form("Constituents_%s", GetJetContainer(i)->GetName());
        fVariablesTreeJetConstituent.push_back((TTree*)fTreeHandlerJet.at(i)->BuildJetConstituentTree(nameoutput,nameoutput));
        fVariablesTreeJetConstituent.at(i)->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeJetConstituent.at(i));
      }
    }
  }
  
  // Post the data
  PostData(1,fNentries);
  PostData(2,fHistoNormCounter);
  PostData(4,fListCounter);
  PostData(5,fTreeEvChar);
  if(fWriteVariableTreeD0){
    PostData(6,fVariablesTreeD0);
    if(fFillMCGenTrees && fReadMC) PostData(7,fGenTreeD0);
  }
  if(fWriteVariableTreeDs){
    PostData(8,fVariablesTreeDs);
    if(fFillMCGenTrees && fReadMC) PostData(9,fGenTreeDs);
  }
  if(fWriteVariableTreeDplus){
    PostData(10,fVariablesTreeDplus);
    if(fFillMCGenTrees && fReadMC) PostData(11,fGenTreeDplus);
  }
  if(fWriteVariableTreeLctopKpi){
    PostData(12,fVariablesTreeLctopKpi);
    if(fFillMCGenTrees && fReadMC) PostData(13,fGenTreeLctopKpi);
  }
  if(fWriteVariableTreeBplus){
    PostData(14,fVariablesTreeBplus);
    if(fFillMCGenTrees && fReadMC) PostData(15,fGenTreeBplus);
  }
  if(fWriteVariableTreeDstar){
    PostData(16,fVariablesTreeDstar);
    if(fFillMCGenTrees && fReadMC) PostData(17,fGenTreeDstar);
  }
  if(fWriteVariableTreeLc2V0bachelor){
    PostData(18,fVariablesTreeLc2V0bachelor);
    if(fFillMCGenTrees && fReadMC) PostData(19,fGenTreeLc2V0bachelor);
  }
  if(fWriteVariableTreeBs){
    PostData(20,fVariablesTreeBs);
    if(fFillMCGenTrees && fReadMC) PostData(21,fGenTreeBs);
  }
  if(fWriteVariableTreeLb){
    PostData(22,fVariablesTreeLb);
    if(fFillMCGenTrees && fReadMC) PostData(23,fGenTreeLb);
  }
  if(fWriteVariableTreeInclusiveJet){
    PostData(24,fVariablesTreeInclusiveJet);
    if(fFillMCGenTrees && fReadMC) PostData(25,fGenTreeInclusiveJet);
  }
  if(fFillParticleTree){
    PostData(26,fVariablesTreeParticle);
    if(fFillMCGenTrees && fReadMC) PostData(27,fVariablesTreeGenParticle);
  }
  if(fFillTrackletTree){
    PostData(28,fVariablesTreeTracklet);
  }
  if(fWriteNJetTrees > 0){
    // Post each jet tree to a separate output slot (for simplicity, keep the jet tree in the last slots)
    const int nJetCollections = fJetCollArray.GetEntriesFast();
    for (int i=0; i<nJetCollections; i++) {
      PostData(29+i,fVariablesTreeJet.at(i));
    }
    // Post jet constituent trees (if enabled)
    if (fFillJetConstituentTrees) {
      for (int i=0; i<nJetCollections; i++) {
        PostData(29+nJetCollections+i,fVariablesTreeJetConstituent.at(i));
      }
    }
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::FillJetTree() {
  
  // If it is the first event, then execute ExecOnce()
  if (!fLocalInitialized){
    ExecOnce();
  }
  
  // Retrieve particles/jets corresponding to each particle/jet container
  if (!RetrieveEventObjects()) {
    return;
  }
  
  // Fill particle trees
  if (fFillParticleTree) {
    fTreeHandlerParticle->FillTree(fRunNumber, fEventID, fEventIDExt, fEventIDLong);
    
    if (fTreeHandlerGenParticle) {
      fTreeHandlerGenParticle->FillTree(fRunNumber, fEventID, fEventIDExt, fEventIDLong);
    }
    
  }
  
  // Set Jet ID for all jets, as index of accepted jets in each event: 0, 1, 2, ..., N
  // Note: We assign the ID before filling the trees, to ensure jet matches can be set in MC
  for(Int_t i =0; i<fJetCollArray.GetEntriesFast(); i++) {
    fTreeHandlerJet.at(i)->SetJetLabels();
  }
  
  // Loop through jet containers, set jet variables for each, and fill each tree
  for(Int_t i =0; i<fJetCollArray.GetEntriesFast(); i++) {
    fTreeHandlerJet.at(i)->FillTree(fRunNumber, fEventID, fEventIDExt, fEventIDLong);
  }
  
}

std::string AliAnalysisTaskSEHFTreeCreator::GetPeriod(const AliVEvent* event){
  Int_t runNo  = event->GetRunNumber();
  // pp: 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  // pPb 2013: 0-LHC13b, 1-LHC13c
  // pPb 2016: 0-LHC16q: 265499->265525; 265309->265387, 1-LHC16q:265435, 2-LHC16q:265388->265427, LHC16t: 267163->267166
  
  if (runNo>195343 && runNo<195484) return "LHC13b";
  if (runNo>195528 && runNo<195678) return "LHC13c";
  if ((runNo>=265499 && runNo<=265525) || (runNo>=265309 && runNo<=265387)) return "LHC16q_0";
  if (runNo == 265435) return "LHC16q_1";
  if (runNo>=265388 && runNo<=265427) return "LHC16q_2";
  if (runNo>=267163 && runNo<=267166) return "LHC16t";
  
  // if(runNo>114930 && runNo<117223) period = 0;
  // if(runNo>119158 && runNo<120830) period = 1;
  // if(runNo>122373 && runNo<126438) period = 2;
  // if(runNo>127711 && runNo<130851) period = 3;
  
  if (runNo>=252235 && runNo<=252375) return "LHC16d";
  if (runNo>=252603 && runNo<=253591) return "LHC16e";
  if (runNo>=254124 && runNo<=254332) return "LHC16g";
  if (runNo>=254378 && runNo<=255469) return "LHC16h_1";
  if (runNo>=254418 && runNo<=254422) return "LHC16h_2";
  if (runNo>=256146 && runNo<=256420) return "LHC16j";
  if (runNo>=256504 && runNo<=258537) return "LHC16k";
  if (runNo>=258883 && runNo<=260187) return "LHC16l";
  if (runNo>=262395 && runNo<=264035) return "LHC16o";
  if (runNo>=264076 && runNo<=264347) return "LHC16p";
  
  if (runNo>=270822 && runNo<=270830) return "LHC17e";
  if (runNo>=270854 && runNo<=270865) return "LHC17f";
  if (runNo>=271868 && runNo<=273103) return "LHC17h";
  if (runNo>=273591 && runNo<=274442) return "LHC17i";
  if (runNo>=274593 && runNo<=274671) return "LHC17j";
  if (runNo>=274690 && runNo<=276508) return "LHC17k";
  if (runNo>=276551 && runNo<=278216) return "LHC17l";
  if (runNo>=278914 && runNo<=280140) return "LHC17m";
  if (runNo>=280282 && runNo<=281961) return "LHC17o";
  if (runNo>=282504 && runNo<=282704) return "LHC17r";
  
  if(runNo>=285008 && runNo<=285447) return "LHC18b";
  if(runNo>=285978 && runNo<=286350) return "LHC18d";
  if(runNo>=286380 && runNo<=286937) return "LHC18e";
  if(runNo>=287000 && runNo<=287977) return "LHC18f";
  if(runNo>=288619 && runNo<=288750) return "LHC18g";
  if(runNo>=288804 && runNo<=288806) return "LHC18h";
  if(runNo>=288861 && runNo<=288909) return "LHC18i";
  if(runNo==288943) return "LHC18j";
  if(runNo>=289165 && runNo<=289201) return "LHC18k";
  if(runNo>=289240 && runNo<=289971) return "LHC18l";
  if(runNo>=290222 && runNo<=292839) return "LHC18m";
  if(runNo>=293357 && runNo<=293359) return "LHC18n";
  if(runNo>=293368 && runNo<=293898) return "LHC18o";
  if(runNo>=294009 && runNo<=294925) return "LHC18p";
  
  return "undefined";
}
//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::UserExec(Option_t */*option*/)
{ 
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  
  //setbuf(stdout, NULL);
  //fflush(stdout);
  /// Execute analysis for current event:
  
  if(fEnableEventDownsampling) {
        gRandom->SetSeed(fSeedEventDownsampling);
        if(gRandom->Rndm() > fFracToKeepEventDownsampling)
            return;
  }

  fBC = aod->GetBunchCrossNumber();
  fOrbit = (Int_t)(aod->GetOrbitNumber());
  fPeriod = (Int_t)(aod->GetPeriodNumber());
  fEventIDLong = (Long64_t)(aod->GetHeader()->GetEventIdAsLong());
  if (!fEventIDLong)
    fEventIDLong = (Long64_t(aod->GetTimeStamp()) << 32) + Long64_t((aod->GetNumberOfTPCClusters()<<5) | (aod->GetNumberOfTPCTracks()));
  fEventIDExt = Int_t(fEventIDLong >> 32);
  fEventID    = Int_t(fEventIDLong & 0xffffffff);

  fNentries->Fill(0); // all events
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fNentries->Fill(2);
      return;
    }
    fNentries->Fill(1);
  }
  
  
  TString candidates2prongArrayName="D0toKpi";
  TString candidates3prongArrayName="Charm3Prong";
  TString candidatesCascArrayDstarName="Dstar";
  TString candidatesCascArrayName="CascadesHF";
  TClonesArray *array2prong=0;
  TClonesArray *array3Prong=0;
  TClonesArray *arrayDstar=0;
  TClonesArray *arrayCasc=0;
  
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      array2prong=(TClonesArray*)aodFromExt->GetList()->FindObject(candidates2prongArrayName.Data());
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject(candidates3prongArrayName.Data());
      arrayDstar=(TClonesArray*)aodFromExt->GetList()->FindObject(candidatesCascArrayDstarName.Data());
      arrayCasc=(TClonesArray*)aodFromExt->GetList()->FindObject(candidatesCascArrayName.Data());
    }
  } else if(aod) {
    array2prong=(TClonesArray*)aod->GetList()->FindObject(candidates2prongArrayName.Data());
    array3Prong=(TClonesArray*)aod->GetList()->FindObject(candidates3prongArrayName.Data());
    arrayDstar=(TClonesArray*)aod->GetList()->FindObject(candidatesCascArrayDstarName.Data());
    arrayCasc=(TClonesArray*)aod->GetList()->FindObject(candidatesCascArrayName.Data());
  }
  
  if(!array2prong || !array3Prong || !arrayDstar || !arrayCasc || !aod) {
    printf("AliAnalysisTaskSEHFTreeCreator::UserExec: input branches not found!\n");
    return;
  }
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;
  fNentries->Fill(3); // count events
  
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;
  
  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC particles branch not found!\n");
      return;
    }
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
      return;
    }
    fzVtxGen = mcHeader->GetVtxZ();     
  }
  
  Bool_t isSameEvSelD0=kTRUE;
  Bool_t isSameEvSelDs=kTRUE;
  Bool_t isSameEvSelDplus=kTRUE;
  Bool_t isSameEvSelLctopKpi=kTRUE;
  Bool_t isSameEvSelBplus=kTRUE;
  Bool_t isSameEvSelBs=kTRUE;
  Bool_t isSameEvSelDstar=kTRUE;
  Bool_t isSameEvSelLc2V0bachelor=kTRUE;
  Bool_t isSameEvSelLb=kTRUE;
  if(fWriteVariableTreeD0)
    isSameEvSelD0=!((fFiltCutsD0toKpi->IsEventSelected(aod) && !fCutsD0toKpi->IsEventSelected(aod))||(!fFiltCutsD0toKpi->IsEventSelected(aod) && fCutsD0toKpi->IsEventSelected(aod)));
  if(fWriteVariableTreeDs)
    isSameEvSelDs=!((fFiltCutsDstoKKpi->IsEventSelected(aod) && !fCutsDstoKKpi->IsEventSelected(aod))||(!fFiltCutsDstoKKpi->IsEventSelected(aod) && fCutsDstoKKpi->IsEventSelected(aod)));
  if(fWriteVariableTreeDplus)
    isSameEvSelDplus=!((fFiltCutsDplustoKpipi->IsEventSelected(aod) && !fCutsDplustoKpipi->IsEventSelected(aod))||(!fFiltCutsDplustoKpipi->IsEventSelected(aod) && fCutsDplustoKpipi->IsEventSelected(aod)));
  if(fWriteVariableTreeLctopKpi)
    isSameEvSelLctopKpi=!((fFiltCutsLctopKpi->IsEventSelected(aod) && !fCutsLctopKpi->IsEventSelected(aod))||(!fFiltCutsLctopKpi->IsEventSelected(aod) && fCutsLctopKpi->IsEventSelected(aod)));
  if(fWriteVariableTreeBplus)
    isSameEvSelBplus = !((fFiltCutsBplustoD0pi->IsEventSelected(aod) && !fCutsBplustoD0pi->IsEventSelected(aod)) || (!fFiltCutsBplustoD0pi->IsEventSelected(aod) && fCutsBplustoD0pi->IsEventSelected(aod)));
  if(fWriteVariableTreeBs)
    isSameEvSelBs = !((fFiltCutsBstoDspi->IsEventSelected(aod) && !fCutsBstoDspi->IsEventSelected(aod)) || (!fFiltCutsBstoDspi->IsEventSelected(aod) && fCutsBstoDspi->IsEventSelected(aod)));
  if(fWriteVariableTreeDstar)
    isSameEvSelDstar=!((fFiltCutsDstartoKpipi->IsEventSelected(aod) && !fCutsDstartoKpipi->IsEventSelected(aod))||(!fFiltCutsDstartoKpipi->IsEventSelected(aod) && fCutsDstartoKpipi->IsEventSelected(aod)));
  if(fWriteVariableTreeLc2V0bachelor)
    isSameEvSelLc2V0bachelor=!((fFiltCutsLc2V0bachelor->IsEventSelected(aod) && !fCutsLc2V0bachelor->IsEventSelected(aod))||(!fFiltCutsLc2V0bachelor->IsEventSelected(aod) && fCutsLc2V0bachelor->IsEventSelected(aod)));
  if(fWriteVariableTreeLb)
    isSameEvSelLb = !((fFiltCutsLbtoLcpi->IsEventSelected(aod) && !fCutsLbtoLcpi->IsEventSelected(aod)) || (!fFiltCutsLbtoLcpi->IsEventSelected(aod) && fCutsLbtoLcpi->IsEventSelected(aod)));
  
  Bool_t isSameEvSel = isSameEvSelD0 && isSameEvSelDs && isSameEvSelDplus && isSameEvSelLctopKpi && isSameEvSelBplus && isSameEvSelBs && isSameEvSelDstar && isSameEvSelLc2V0bachelor && isSameEvSelLb;
  if(!isSameEvSel) {
    Printf("AliAnalysisTaskSEHFTreeCreator::UserExec: differences in the event selection cuts same meson");
    return;
  }
  if((fWriteVariableTreeD0 && fWriteVariableTreeDs && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsDstoKKpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeDplus && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsDplustoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeDplus && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsDplustoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeLctopKpi && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsLctopKpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeLctopKpi && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsLctopKpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDplus && fWriteVariableTreeLctopKpi && (fFiltCutsDplustoKpipi->IsEventSelected(aod)!=fFiltCutsLctopKpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeBplus && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsBplustoD0pi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeBplus && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsBplustoD0pi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDplus && fWriteVariableTreeBplus && (fFiltCutsDplustoKpipi->IsEventSelected(aod)!=fFiltCutsBplustoD0pi->IsEventSelected(aod))) ||
     (fWriteVariableTreeLctopKpi && fWriteVariableTreeBplus && (fFiltCutsLctopKpi->IsEventSelected(aod)!=fFiltCutsBplustoD0pi->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeDstar && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsDstartoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDplus && fWriteVariableTreeDstar && (fFiltCutsDplustoKpipi->IsEventSelected(aod)!=fFiltCutsDstartoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeDstar && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsDstartoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeLctopKpi && fWriteVariableTreeDstar && (fFiltCutsLctopKpi->IsEventSelected(aod)!=fFiltCutsDstartoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeBplus && fWriteVariableTreeDstar && (fFiltCutsBplustoD0pi->IsEventSelected(aod)!=fFiltCutsDstartoKpipi->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeLc2V0bachelor && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod))) ||
     (fWriteVariableTreeDplus && fWriteVariableTreeLc2V0bachelor && (fFiltCutsDplustoKpipi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeLc2V0bachelor && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod))) ||
     (fWriteVariableTreeLctopKpi && fWriteVariableTreeLc2V0bachelor && (fFiltCutsLctopKpi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod))) ||
     (fWriteVariableTreeBplus && fWriteVariableTreeLc2V0bachelor && (fFiltCutsBplustoD0pi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod))) ||
     (fWriteVariableTreeDstar && fWriteVariableTreeLc2V0bachelor && (fFiltCutsDstartoKpipi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeBs && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDplus && fWriteVariableTreeBs && (fFiltCutsDplustoKpipi->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeBs && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeLctopKpi && fWriteVariableTreeBs && (fFiltCutsLctopKpi->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeBplus && fWriteVariableTreeBs && (fFiltCutsBplustoD0pi->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDstar && fWriteVariableTreeBs && (fFiltCutsDstartoKpipi->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeLc2V0bachelor && fWriteVariableTreeBs && (fFiltCutsLc2V0bachelor->IsEventSelected(aod)!=fFiltCutsBstoDspi->IsEventSelected(aod))) ||
     (fWriteVariableTreeD0 && fWriteVariableTreeLb && (fFiltCutsD0toKpi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDplus && fWriteVariableTreeLb && (fFiltCutsDplustoKpipi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDs && fWriteVariableTreeLb && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeLctopKpi && fWriteVariableTreeLb && (fFiltCutsLctopKpi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeBplus && fWriteVariableTreeLb && (fFiltCutsBplustoD0pi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeDstar && fWriteVariableTreeLb && (fFiltCutsDstartoKpipi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeLc2V0bachelor && fWriteVariableTreeLb && (fFiltCutsLc2V0bachelor->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod))) ||
     (fWriteVariableTreeBs && fWriteVariableTreeLb && (fFiltCutsBstoDspi->IsEventSelected(aod)!=fFiltCutsLbtoLcpi->IsEventSelected(aod)))
     ){
    Printf("AliAnalysisTaskSEHFTreeCreator::UserExec: differences in the event selection cuts different meson");
    return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx = (AliAODVertex*)aod->GetPrimaryVertex();
  fNcontributors = vtx->GetNContributors();
  fzVtxReco = vtx->GetZ();
  fNtracks = aod->GetNumberOfTracks();
  fRunNumber=aod->GetRunNumber();

  //n tracklets
  AliAODTracklets* tracklets=aod->GetTracklets();
  Int_t nTr=tracklets->GetNumberOfTracklets();
  Int_t countTreta1=0;
  for(Int_t iTr=0; iTr<nTr; iTr++){
    Double_t theta=tracklets->GetTheta(iTr);
    Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
    if(eta>-1.0 && eta<1.0) countTreta1++;//count at central rapidity
  }
  fnTracklets=countTreta1;
  fnTrackletsCorr = -1.;

  TProfile *estimatorAvg = fMultEstimatorAvg[GetPeriod(aod)];
  if (fCorrNtrVtx && estimatorAvg)
    fnTrackletsCorr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg, countTreta1, vtx->GetZ(), fRefMult));
  TProfile *estimatorAvgSHM = fMultEstimatorAvgSHM[GetPeriod(aod)];
  if (fCorrNtrVtx && estimatorAvgSHM)
    fnTrackletsCorrSHM = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvgSHM, countTreta1, vtx->GetZ(), fRefMultSHM));

  fCounter->StoreEvent(aod,fEvSelectionCuts,fReadMC,fnTrackletsCorr);

  Bool_t isEvSel=fEvSelectionCuts->IsEventSelected(aod);
      
  if(fEvSelectionCuts->IsEventRejectedDueToTrigger())fNentries->Fill(5);
  if(fEvSelectionCuts->IsEventRejectedDueToNotRecoVertex())fNentries->Fill(6);
  if(fEvSelectionCuts->IsEventRejectedDueToVertexContributors())fNentries->Fill(7);
  if(fEvSelectionCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fNentries->Fill(8);
  if(fEvSelectionCuts->IsEventRejectedDueToPileup())fNentries->Fill(9);
  if(fEvSelectionCuts->IsEventRejectedDueToCentrality())fNentries->Fill(10);
  
  fCentrality = fEvSelectionCuts->GetCentrality(aod);
  if(fCentrality<0) fCentrality=-1.;
  bool isEvRejPhysSel = fEvSelectionCuts->IsEventRejectedDuePhysicsSelection();
  //normalisation counter
  if(!isEvRejPhysSel){
    if(isEvSel){
      //selected events with primary vertex
      fHistoNormCounter->Fill(0.,fCentrality);
    }
    else{
      if(fEvSelectionCuts->GetWhyRejection()==0){
        //rejected events bc no primary vertex
        fHistoNormCounter->Fill(1.,fCentrality);
      }
      //rejected events bc good primary vertex but >10cm
      if(fEvSelectionCuts->GetWhyRejection()==6){
        //nPrimaryV++;
        fHistoNormCounter->Fill(0.,fCentrality);
        //nzVtxGT10++;
        fHistoNormCounter->Fill(2.,fCentrality);
      }
      if(fEvSelectionCuts->GetWhyRejection()==1){
        //nPileup++;
        fHistoNormCounter->Fill(4.,fCentrality);
      }
    }
    if(fEvSelectionCuts->CountEventForNormalization()){
      //nCountForNorm++;
      fHistoNormCounter->Fill(3.,fCentrality);
    }
  }
  Bool_t isEvRejCent  = fEvSelectionCuts->IsEventRejectedDueToCentrality();
  
  if(!isEvSel && (isEvRejCent || (fApplyPhysicsSelOnline && isEvRejPhysSel))){
    return; //cut only centrality and physics selection if enabled, else tag only
  }
  if(isEvSel) fNentries->Fill(4);
  fIsEvRej = fEvSelectionCuts->GetEventRejectionBitMap();
    
  auto trig_mask_cuts = fEvSelectionCuts->GetTriggerMask();
    
  fEvSelectionCuts->SetTriggerMask(AliVEvent::kINT7);
  fIsEvSel_INT7 = fEvSelectionCuts->IsEventSelected(aod);
  fIsEvRej_INT7 = fEvSelectionCuts->GetEventRejectionBitMap();

  fEvSelectionCuts->SetTriggerMask(AliVEvent::kHighMultSPD);
  fIsEvSel_HighMultSPD = fEvSelectionCuts->IsEventSelected(aod);
  fIsEvRej_HighMultSPD = fEvSelectionCuts->GetEventRejectionBitMap();

  fEvSelectionCuts->SetTriggerMask(AliVEvent::kHighMultV0);
  fIsEvSel_HighMultV0 = fEvSelectionCuts->IsEventSelected(aod);
  fIsEvRej_HighMultV0 = fEvSelectionCuts->GetEventRejectionBitMap();
    
  fEvSelectionCuts->SetTriggerMask(AliVEvent::kEMCEJE);
  fIsEvSel_EMCEJE = fEvSelectionCuts->IsEventSelected(aod);
  fIsEvRej_EMCEJE = fEvSelectionCuts->GetEventRejectionBitMap();

  fEvSelectionCuts->SetTriggerMask(trig_mask_cuts);

  if(fDoPtHard){
    TString currentfilepath_pyxsec = ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()))->GetName();
    for (Int_t s_end=0; s_end<11; s_end++) currentfilepath_pyxsec.Remove(currentfilepath_pyxsec.Length()-1);
    TString pyxsec_name = "pyxsec_hists.root";
    TFile *F_xsection = new TFile(currentfilepath_pyxsec+pyxsec_name);
    TList *L_xsection = (TList *) F_xsection->Get("cFilterList");
    TProfile *fh_xsection = (TProfile *) L_xsection->FindObject("h1Xsec");
    TH1D *fh_Trials = (TH1D *) L_xsection->FindObject("h1Trials");

    //fCross_Section = mcHeader->GetCrossSection();
    //fTrials = mcHeader->GetTrials();
    fCross_Section = fh_xsection->GetBinContent(1);
    fTrials = fh_Trials->GetBinContent(1);
    fpthard = mcHeader->GetPtHard();

    delete F_xsection;
    delete L_xsection;
    delete fh_xsection;
    delete fh_Trials;

  }

  //V0 multiplicities
  AliAODVZERO *vzeroAOD = (AliAODVZERO*)aod->GetVZEROData();
  Double_t vzeroA = vzeroAOD ? vzeroAOD->GetMTotV0A() : 0.;
  Double_t vzeroC = vzeroAOD ? vzeroAOD->GetMTotV0C() : 0.;
  Double_t vzeroAEq = AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(aod);
  Double_t vzeroCEq = AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(aod);
  fnV0A = static_cast<Int_t>(vzeroA);
  fnV0M = static_cast<Int_t>(vzeroA + vzeroC);
  fnV0MEq = static_cast<Int_t>(vzeroAEq + vzeroCEq);
  fnV0MCorr = -1;
  fnV0MEqCorr = -1;
  if (fCorrV0MVtx) {
    fnV0MCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroA, vtx->GetZ()) + AliESDUtils::GetCorrV0C(vzeroC, vtx->GetZ()));
    fnV0MEqCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroAEq, vtx->GetZ()) + AliESDUtils::GetCorrV0C(vzeroCEq, vtx->GetZ()));
  }

  // multiplicity percentiles
  const auto multSel = static_cast<AliMultSelection*>(aod->FindListObject("MultSelection"));
  fPercV0M = multSel ? multSel->GetMultiplicityPercentile("V0M") : -1.;
  // multiplicity from mult selection task
  const auto multEst = multSel ? multSel->GetEstimator("V0M") : nullptr;
  fMultV0M = multEst ? multEst->GetValue() : -1.;

  // generated multiplicity
  fMultGen = -1;
  fMultGenV0A = -1;
  fMultGenV0C = -1;
  if (fReadMC) {
    TClonesArray *arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    fMultGen = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-1.0,1.0);
    fMultGenV0A = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,2.8,5.1);
    fMultGenV0C = AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(arrayMC,-3.7,-1.7);
  }


  // Extract fired triggers
  fTriggerMask = static_cast<AliVAODHeader*>(aod->GetHeader())->GetOfflineTrigger();
  fTriggerBitINT7 = static_cast<bool>(fTriggerMask & AliVEvent::kINT7);
  fTriggerBitHighMultSPD = static_cast<bool>(fTriggerMask & AliVEvent::kHighMultSPD);
  fTriggerBitHighMultV0 = static_cast<bool>(fTriggerMask & AliVEvent::kHighMultV0);
  fTriggerBitCentral = static_cast<bool>(fTriggerMask & AliVEvent::kCentral);
  fTriggerBitSemiCentral = static_cast<bool>(fTriggerMask & AliVEvent::kSemiCentral);
  fTriggerBitEMCEJE = static_cast<bool>(fTriggerMask & AliVEvent::kEMCEJE);
  
  fTriggerClasses = aod->GetFiredTriggerClasses();
  fTriggerClassINT7 = fTriggerClasses.Contains("CINT7-B");
  fTriggerClassHighMultSPD = fTriggerClasses.Contains("CVHMSH2-B");
  fTriggerClassHighMultV0m = fTriggerClasses.Contains("CVHMV0M-B");
  fTriggerClassEMCALEJ1 = fTriggerClasses.Contains("EJ1");
  fTriggerClassEMCALEJ2 = fTriggerClasses.Contains("EJ2");
  fTriggerClassDCALDJ1 = fTriggerClasses.Contains("DJ1");
  fTriggerClassDCALDJ2 = fTriggerClasses.Contains("DJ2");
  
  // bits for CTP inputs
  if (fRunNumberCDB != fRunNumber) {
    fCdbEntry = AliCDBManager::Instance()->Get("GRP/CTP/Config", fRunNumber);
    fRunNumberCDB = fRunNumber;
  }

  AliTriggerConfiguration *trgCfg = fCdbEntry ? static_cast<AliTriggerConfiguration*>(fCdbEntry->GetObject()) : nullptr;
  TObjArray inputs;
  if (trgCfg)
    inputs = trgCfg->GetInputs();
  const auto inputSHM = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("0SHM")) : nullptr;
  const auto inputV0M = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("0VHM")) : nullptr;
  const auto inputV0A = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("0V0A")) : nullptr;
  const auto inputV0C = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("0V0C")) : nullptr;
  const auto inputEMCALEJ1 = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("1EJ1")) : nullptr;
  const auto inputEMCALEJ2 = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("1EJ2")) : nullptr;
  const auto inputDCALDJ1 = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("1DJ1")) : nullptr;
  const auto inputDCALDJ2 = trgCfg ? static_cast<AliTriggerInput*>(inputs.FindObject("1DJ2")) : nullptr;
  const auto triggerBits = aod->GetHeader()->GetL0TriggerInputs();
  fTriggerOnlineHighMultSPD = inputSHM ? TESTBIT(triggerBits, inputSHM->GetIndexCTP() - 1) : -1;
  fTriggerOnlineHighMultV0 = inputV0M ? TESTBIT(triggerBits, inputV0M->GetIndexCTP() - 1) : -1;
  fTriggerOnlineINT7 = (inputV0C && inputV0A) ?
                       (TESTBIT(triggerBits, inputV0C->GetIndexCTP() - 1) &&
                        TESTBIT(triggerBits, inputV0A->GetIndexCTP() - 1)) : -1;

  fTriggerOnlineEMCALEJ1 = inputEMCALEJ1 ? TESTBIT(triggerBits, inputEMCALEJ1->GetIndexCTP() - 1) : -1;
  fTriggerOnlineEMCALEJ2 = inputEMCALEJ2 ? TESTBIT(triggerBits, inputEMCALEJ2->GetIndexCTP() - 1) : -1;
  fTriggerOnlineDCALDJ1 = inputDCALDJ1 ? TESTBIT(triggerBits, inputDCALDJ1->GetIndexCTP() - 1) : -1;
  fTriggerOnlineDCALDJ2 = inputDCALDJ2 ? TESTBIT(triggerBits, inputDCALDJ2->GetIndexCTP() - 1) : -1;
  
  fTreeEvChar->Fill(); 
  //get PID response
  if(!fPIDresp) fPIDresp = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  
  if(fWriteVariableTreeD0) Process2Prong(array2prong,aod,mcArray,aod->GetMagneticField(),mcHeader);
  if(fWriteVariableTreeDs || fWriteVariableTreeDplus || fWriteVariableTreeLctopKpi) Process3Prong(array3Prong,aod,mcArray,aod->GetMagneticField(),mcHeader);
  if(fWriteVariableTreeDstar) ProcessDstar(arrayDstar,aod,mcArray,aod->GetMagneticField());
  if(fWriteVariableTreeLc2V0bachelor) ProcessCasc(arrayCasc,aod,mcArray,aod->GetMagneticField());
  if(fWriteVariableTreeBplus) ProcessBplus(array2prong,aod,mcArray,aod->GetMagneticField(),mcHeader);
  if(fWriteVariableTreeBs) ProcessBs(array3Prong,aod,mcArray,aod->GetMagneticField(),mcHeader);
  if(fWriteVariableTreeLb) ProcessLb(array3Prong,aod,mcArray,aod->GetMagneticField(),mcHeader);
  if(fWriteVariableTreeInclusiveJet){
    ProcessInclusiveJet(aod,mcArray);
    if(fFillMCGenTrees && fReadMC) ProcessMCGenInclusiveJet(mcArray);
  }
  if(fFillMCGenTrees && fReadMC) ProcessMCGen(mcArray);
  
  
  // Fill the jet tree
  if (fWriteNJetTrees > 0 || fFillParticleTree) {
    FillJetTree();
  }
  if (fFillTrackletTree){
    fTreeHandlerTracklet->SetTrackletContainer(aod->GetTracklets());
    fTreeHandlerTracklet->FillTree(fRunNumber, fEventID, fEventIDExt, fEventIDLong);
  }
  
  // Post the data
  PostData(1,fNentries);
  PostData(2,fHistoNormCounter);
  PostData(4,fListCounter);
  PostData(5,fTreeEvChar);
  if(fWriteVariableTreeD0){
    PostData(6,fVariablesTreeD0);
    if(fFillMCGenTrees && fReadMC) PostData(7,fGenTreeD0);
  }
  if(fWriteVariableTreeDs){
    PostData(8,fVariablesTreeDs);
    if(fFillMCGenTrees && fReadMC) PostData(9,fGenTreeDs);
  }
  if(fWriteVariableTreeDplus){
    PostData(10,fVariablesTreeDplus);
    if(fFillMCGenTrees && fReadMC) PostData(11,fGenTreeDplus);
  }
  if(fWriteVariableTreeLctopKpi){
    PostData(12,fVariablesTreeLctopKpi);
    if(fFillMCGenTrees && fReadMC) PostData(13,fGenTreeLctopKpi);
  }
  if(fWriteVariableTreeBplus){
    PostData(14,fVariablesTreeBplus);
    if(fFillMCGenTrees && fReadMC) PostData(15,fGenTreeBplus);
  }
  if(fWriteVariableTreeDstar){
    PostData(16,fVariablesTreeDstar);
    if(fFillMCGenTrees && fReadMC) PostData(17,fGenTreeDstar);
  }
  if(fWriteVariableTreeLc2V0bachelor){
    PostData(18,fVariablesTreeLc2V0bachelor);
    if(fFillMCGenTrees && fReadMC) PostData(19,fGenTreeLc2V0bachelor);
  }
  if(fWriteVariableTreeBs){
    PostData(20,fVariablesTreeBs);
    if(fFillMCGenTrees && fReadMC) PostData(21,fGenTreeBs);
  }
  if(fWriteVariableTreeLb){
    PostData(22,fVariablesTreeLb);
    if(fFillMCGenTrees && fReadMC) PostData(23,fGenTreeLb);
  }
  if(fWriteVariableTreeInclusiveJet){
    PostData(24,fVariablesTreeInclusiveJet);
    if(fFillMCGenTrees && fReadMC) PostData(25,fGenTreeInclusiveJet);
  }
  if(fFillParticleTree){
    PostData(26,fVariablesTreeParticle);
    if(fFillMCGenTrees && fReadMC) PostData(27,fVariablesTreeGenParticle);
  } 
  if(fFillTrackletTree){
    PostData(28,fVariablesTreeTracklet);
  }
  if(fWriteNJetTrees > 0){
    // Post each jet tree to a separate output slot (for simplicity, keep the jet tree in the last slots)
    const int nJetCollections = fJetCollArray.GetEntriesFast();
    for (int i=0; i<nJetCollections; i++) {
      PostData(29+i,fVariablesTreeJet.at(i));
    }
    // Post jet constituent trees (if enabled)
    if (fFillJetConstituentTrees) {
      for (int i=0; i<nJetCollections; i++) {
        PostData(29+nJetCollections+i,fVariablesTreeJetConstituent.at(i));
      }
    }
  }
  
  return;
}

/**
 * Perform steps needed to initialize the analysis.
 * This function relies on the presence of an input
 * event (ESD or AOD event). Consequently it is called
 * internally by UserExec for the first event.
 *
 * This function connects all containers attached to
 * this task to the corresponding arrays in the
 * input event.
 *
 * (Copied from AliAnalysisTaskEmcal / AliAnalysisTaskEmcalJet).
 */
void AliAnalysisTaskSEHFTreeCreator::ExecOnce()
{
  
  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }
  
  //Load all requested track branches - each container knows name already
  for (Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetArray(InputEvent());
  }
  
  if (fParticleCollArray.GetEntriesFast()>0) {
    
    AliParticleContainer *cont = GetParticleContainer(0);
    if (!cont) {
      AliError(Form("%s: Particle container %d not found",GetName(),0));
      return;
    }
    TString contName = cont->GetArrayName();
    TClonesArray* tracks = cont->GetArray();
    if (!tracks) {
      AliError(Form("%s: Could not retrieve first track branch!", GetName()));
      return;
    }
  }
  
  //Load all requested jet branches - each container knows name already
  if(fJetCollArray.GetEntriesFast()==0) {
    AliWarning("There are no jet collections");
  }
  else {
    
    // get rho from the event
    if (!fRhoName.IsNull() && !fRho) {
      fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
      if (!fRho) {
        AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
        fLocalInitialized = kFALSE;
        return;
      }
    }
    
    for(Int_t i =0; i<fJetCollArray.GetEntriesFast(); i++) {
      AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
      cont->SetRunNumber(InputEvent()->GetRunNumber());
      cont->SetArray(InputEvent());
      cont->LoadRho(InputEvent());
    }
    
    //Get Jets, cuts and rho for first jet container
    AliJetContainer *cont = GetJetContainer(0);
    
    if (!cont->GetArrayName().IsNull()) {
      TClonesArray *jets = cont->GetArray();
      if(!jets && fJetCollArray.GetEntriesFast()>0) {
        AliErrorStream() << GetName() << ": Could not retrieve first jet branch!\n";
        std::stringstream foundbranches;
        bool first(true);
        for(auto e : *(InputEvent()->GetList())){
          if(first){
            // Skip printing a comma on the first time through
            first = false;
          }
          else {
            foundbranches << ", ";
          }
          foundbranches << e->GetName();
        }
        std::string fbstring = foundbranches.str();
        AliErrorStream() << "Found branches: " << fbstring << std::endl;
        fLocalInitialized = kFALSE;
        return;
      }
    }
    
    if (!fRho) { // if rho name is not provided, tries to use the rho object of the first jet branch
      fRhoName = cont->GetRhoName();
      fRho = cont->GetRhoParameter();
    }
  }
  
  fLocalInitialized = kTRUE;
}

/**
 * Retrieve objects from event. This operation needs to be performed for every event.
 * @return kTRUE if successful, kFALSE otherwise
 */
Bool_t AliAnalysisTaskSEHFTreeCreator::RetrieveEventObjects()
{
  
  AliEmcalContainer* cont = nullptr;
  TIter nextPartColl(&fParticleCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextPartColl()))){
    cont->NextEvent(InputEvent());
  }
  
  if (fRho) fRhoVal = fRho->GetVal();
  
  AliEmcalContainer* contJet = nullptr;
  TIter nextJetColl(&fJetCollArray);
  while ((contJet = static_cast<AliEmcalContainer*>(nextJetColl()))) contJet->NextEvent(InputEvent());
  
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AliAnalysisTaskSEHFTreeCreator: Terminate() \n");
  
  
  fNentries = dynamic_cast<TH1F*>(GetOutputData(1));
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }
  fHistoNormCounter = dynamic_cast<TH2F*>(GetOutputData(2));
  if(!fHistoNormCounter){
    printf("ERROR: fHistoNormCounter not available\n");
    return;
  }
  fListCuts = dynamic_cast<TList*>(GetOutputData(3));
  if(!fListCuts){
    printf("ERROR: fListCuts not available\n");
    return;
  }
  fListCounter = dynamic_cast<TList*>(GetOutputData(4));
  if(!fListCounter){
    printf("ERROR: fListCounter not available\n");
    return;
  }
  return;
}
//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::Process2Prong(TClonesArray *array2prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t n2prong = array2prong->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",n2prong);
  
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t nSelectedD0=0;
  Int_t nFilteredD0=0;
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  
  for (Int_t i2prong = 0; i2prong < n2prong; i2prong++) {
    fNentries->Fill(11);
    
    //D0
    Bool_t isD0tagged=kTRUE;
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)array2prong->UncheckedAt(i2prong);
    if(fUseSelectionBit && d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
      isD0tagged=kFALSE;
    }
    
    if(isD0tagged && fWriteVariableTreeD0){

      Int_t preSelectedD0 = -1;
      if(fITSUpgradePreSelect){
        TObjArray arrTracks(2);
        for(Int_t ipr=0;ipr<2;ipr++){
          AliAODTrack *tr=vHF->GetProng(aod,d,ipr);
          arrTracks.AddAt(tr,ipr);
        }
        preSelectedD0 = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 2, 421, pdgDgD0toKpi);
        if(preSelectedD0 == 0) continue; //Mixture hijing + injected
        if(preSelectedD0 == 2) continue; //Only MatchedToMC injected signal
        if(fWriteOnlySignal == 1 && preSelectedD0 != 1) continue; //Only matched signal when only signal is requested
      }

      fNentries->Fill(12);
      nFilteredD0++;
      if((vHF->FillRecoCand(aod,d))) {//Fill the data members of the candidate only if they are empty.
        
        //filtering cuts
        Int_t isSelectedFilt     = fFiltCutsD0toKpi->IsSelected(d,AliRDHFCuts::kAll,aod); //selected
        Int_t isSelectedPidFilt  = fFiltCutsD0toKpi->IsSelectedPID(d);
        bool isUsePidFilt = fFiltCutsD0toKpi->GetIsUsePID();
        if(isUsePidFilt) fFiltCutsD0toKpi->SetUsePID(kFALSE);
        Int_t isSelectedTopoFilt = fFiltCutsD0toKpi->IsSelected(d,AliRDHFCuts::kAll,aod);
        fFiltCutsD0toKpi->SetUsePID(isUsePidFilt);
        
        if(isSelectedFilt){
          fNentries->Fill(13);
          nSelectedD0++;
          
          //analysis cuts
          Int_t isSelectedAnalysis = fCutsD0toKpi->IsSelected(d,AliRDHFCuts::kAll,aod); //selected
          Bool_t isSelAnCutsD0=kFALSE;
          Bool_t isSelAnCutsD0bar=kFALSE;
          if(isSelectedAnalysis==1 || isSelectedAnalysis==3) isSelAnCutsD0=kTRUE;
          if(isSelectedAnalysis==2 || isSelectedAnalysis==3) isSelAnCutsD0bar=kTRUE;
          Int_t isSelectedPidAnalysis = fCutsD0toKpi->IsSelectedPID(d); //selected
          Bool_t isSelPidAnCutsD0=kFALSE;
          Bool_t isSelPidAnCutsD0bar=kFALSE;
          if(isSelectedPidAnalysis==1 || isSelectedPidAnalysis==3) isSelPidAnCutsD0=kTRUE;
          if(isSelectedPidAnalysis==2 || isSelectedPidAnalysis==3) isSelPidAnCutsD0bar=kTRUE;
          bool isUsePidAn = fCutsD0toKpi->GetIsUsePID();
          if(isUsePidAn) fCutsD0toKpi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsD0toKpi->IsSelected(d,AliRDHFCuts::kAll,aod);
          Bool_t isSelTopoAnCutsD0=kFALSE;
          Bool_t isSelTopoAnCutsD0bar=kFALSE;
          if(isSelectedTopoAnalysis==1 || isSelectedTopoAnalysis==3) isSelTopoAnCutsD0=kTRUE;
          if(isSelectedTopoAnalysis==2 || isSelectedTopoAnalysis==3) isSelTopoAnCutsD0bar=kTRUE;
          fCutsD0toKpi->SetUsePID(isUsePidAn);
          Bool_t isSelTracksAnCuts=kFALSE;
          Int_t isSelectedTrackAnalysis = fCutsD0toKpi->IsSelected(d,AliRDHFCuts::kTracks,aod);
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          fTreeHandlerD0->SetIsDzeroDzeroBar(isSelectedAnalysis, isSelectedTopoAnalysis, isSelectedPidAnalysis, isSelectedFilt, isSelectedTopoFilt, isSelectedPidFilt);
          
          Bool_t unsetvtx=kFALSE;
          if(!d->GetOwnPrimaryVtx()){
            d->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsD0toKpi->GetIsPrimaryWithoutDaughters()){
            if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
            if(fFiltCutsD0toKpi->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
            else fFiltCutsD0toKpi->CleanOwnPrimaryVtx(d,aod,origownvtx);
          }
          
          Int_t labD0 = -1;
          Int_t pdgD0 = -99;
          Int_t origin= -1;
          Float_t ptGenD0 = -99.;
          
          AliAODMCParticle *partD0=0x0;
          if(fReadMC) {
            labD0 = d->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
            if(labD0>=0){
              partD0 = (AliAODMCParticle*)arrMC->At(labD0);
              pdgD0 = partD0->GetPdgCode();
              ptGenD0 = partD0->Pt();
              origin = AliVertexingHFUtils::CheckOrigin(arrMC,partD0,!fITSUpgradeProduction);
            }
          }
          
          bool issignal = kFALSE;
          bool isbkg =    kFALSE;
          bool isFD =     kFALSE;
          bool isprompt = kFALSE;
          bool isrefl =   kFALSE;
          Int_t masshypo = 0;
          
          if (isSelectedFilt==1 || isSelectedFilt==3) { //D0
            masshypo=0;
            if(fReadMC){
              if(labD0>=0){
                if(origin==4 || origin==5) {
                  if(origin==4) isprompt=kTRUE;
                  else if(origin==5) isFD=kTRUE;
                  if(pdgD0==421){
                    issignal=kTRUE;
                  }
                  else {
                    isrefl=kTRUE;
                  }
                }
              }//end labD0check
              else{//background
                if(fStoreOnlyHIJINGBackground){
                  Bool_t isHijing = kFALSE;
                  if (mcHeader) isHijing = IsCandidateFromHijing(d,mcHeader,arrMC);
                  if (isHijing) isbkg = kTRUE;
                } else {
                  isbkg=kTRUE;
                }
              }
              if(issignal || isbkg || isrefl) fTreeHandlerD0->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
            }//end read MC
            if(!fReadMC || (issignal || isbkg || isrefl)) {
              fTreeHandlerD0->SetIsSelectedStd(isSelAnCutsD0, isSelTopoAnCutsD0, isSelPidAnCutsD0, isSelTracksAnCuts);
              fTreeHandlerD0->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenD0,d,bfield,masshypo,fPIDresp);
              if (fFillJets) fTreeHandlerD0->SetJetVars(aod->GetTracks(),d,d->InvMassD0(),arrMC,partD0);
              fTreeHandlerD0->FillTree();
            }
          }//end D0
          if (isSelectedFilt>1){//D0bar
            issignal = kFALSE;
            isbkg =    kFALSE;
            isFD =     kFALSE;
            isprompt = kFALSE;
            isrefl =   kFALSE;
            masshypo = 1;
            if(fReadMC){
              if(labD0>=0){
                if(origin==4 || origin==5) {
                  if(origin==4) isprompt=kTRUE;
                  else if(origin==5) isFD=kTRUE;
                  if(pdgD0==-421){
                    issignal=kTRUE;
                  }
                  else {
                    isrefl=kTRUE;
                  }
                }
              } //end label check
              else{ //background MC
                if(fStoreOnlyHIJINGBackground){
                  Bool_t isHijing = kFALSE;
                  if (mcHeader) isHijing = IsCandidateFromHijing(d,mcHeader,arrMC);
                  if (isHijing) isbkg = kTRUE;
                } else {
                  isbkg=kTRUE;
                }
              }
              if(issignal || isbkg || isrefl) fTreeHandlerD0->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
            }//end readMC
            if(!fReadMC || (issignal || isbkg || isrefl)) {
              fTreeHandlerD0->SetIsSelectedStd(isSelAnCutsD0bar, isSelTopoAnCutsD0bar, isSelPidAnCutsD0bar, isSelTracksAnCuts);
              fTreeHandlerD0->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenD0,d,bfield,masshypo,fPIDresp);
              if (fFillJets) fTreeHandlerD0->SetJetVars(aod->GetTracks(),d,d->InvMassD0bar(),arrMC,partD0);
              fTreeHandlerD0->FillTree();
            }
          }//end D0bar
          if(recVtx)fFiltCutsD0toKpi->CleanOwnPrimaryVtx(d,aod,origownvtx);
          if(unsetvtx) d->UnsetOwnPrimaryVtx();
        }//end is selected filt
      }
      else {
        fNentries->Fill(14); //monitor how often this fails
      }
    }//end D0
  }//end loop on candidates
  
  delete vHF;
  return;
}

//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::Process3Prong(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t n3prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of 3prongs: %d\n",n3prong);
  
  Int_t pdgDstoKKpi[3]={321,321,211};
  Int_t nSelectedDs=0;
  Int_t nFilteredDs=0;
  
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  Int_t nSelectedDplus=0;
  Int_t nFilteredDplus=0;
  
  Int_t pdgLctopKpi[3]={2212,321,211};
  Int_t nSelectedLctopKpi=0;
  Int_t nFilteredLctopKpi=0;
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  
  for (Int_t i3prong = 0; i3prong < n3prong; i3prong++) {
    fNentries->Fill(15);
    
    //Ds
    Bool_t isDstagged=kTRUE;
    AliAODRecoDecayHF3Prong *ds    = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3prong);
    if(fUseSelectionBit && !(ds->HasSelectionBit(AliRDHFCuts::kDsCuts))){
      isDstagged=kFALSE;
    }
    
    if(isDstagged && fWriteVariableTreeDs){

      Int_t preSelectedDs = -1;
      if(fITSUpgradePreSelect){
        TObjArray arrTracks(3);
        for(Int_t ipr=0;ipr<3;ipr++){
          AliAODTrack *tr=vHF->GetProng(aod,ds,ipr);
          arrTracks.AddAt(tr,ipr);
        }
        preSelectedDs = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 3, 431, pdgDstoKKpi);
        if(preSelectedDs == 0) continue; //Mixture hijing + injected
        if(preSelectedDs == 2) continue; //Only MatchedToMC injected signal
        if(fWriteOnlySignal == 1 && preSelectedDs != 1) continue; //Only matched signal when only signal is requested
      }

      fNentries->Fill(16);
      nFilteredDs++;
      if((vHF->FillRecoCand(aod,ds))) {////Fill the data members of the candidate only if they are empty.
        
        Int_t isSelectedFilt=fFiltCutsDstoKKpi->IsSelected(ds,AliRDHFCuts::kAll,aod);
        Int_t isKKpi=isSelectedFilt&1;
        Int_t ispiKK=isSelectedFilt&2;
        Int_t isPhiKKpi=isSelectedFilt&4;
        Int_t isPhipiKK=isSelectedFilt&8;
        Int_t isK0starKKpi=isSelectedFilt&16;
        Int_t isK0starpiKK=isSelectedFilt&32;
        
        if(isSelectedFilt>0){
          fNentries->Fill(17);
          nSelectedDs++;
          
          //test analysis cuts
          Bool_t isSelAnCutsKKpi=kFALSE;
          Bool_t isSelAnCutspiKK=kFALSE;
          Bool_t isSelAnPidCutsKKpi=kFALSE;
          Bool_t isSelAnPidCutspiKK=kFALSE;
          Bool_t isSelAnTopoCutsKKpi=kFALSE;
          Bool_t isSelAnTopoCutspiKK=kFALSE;
          Int_t isSelectedAnalysis=fCutsDstoKKpi->IsSelected(ds,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis=fCutsDstoKKpi->IsSelectedPID(ds);
          bool isUsePidAn = fCutsDstoKKpi->GetIsUsePID();
          if(isUsePidAn) fCutsDstoKKpi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsDstoKKpi->IsSelected(ds,AliRDHFCuts::kAll,aod);
          fCutsDstoKKpi->SetUsePID(isUsePidAn);
          Bool_t isSelTracksAnCuts=kFALSE;
          Int_t isSelectedTrackAnalysis = fCutsDstoKKpi->IsSelected(ds,AliRDHFCuts::kTracks,aod);
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          if(fWriteVariableTreeDs==1) {
            if(isSelectedAnalysis&4) isSelAnCutsKKpi=kTRUE;
            if(isSelectedAnalysis&8) isSelAnCutspiKK=kTRUE;
            if(isSelectedPidAnalysis==1 || isSelectedPidAnalysis==3) isSelAnPidCutsKKpi=kTRUE;
            if(isSelectedPidAnalysis==2 || isSelectedPidAnalysis==3) isSelAnPidCutspiKK=kTRUE;
            if(isSelectedTopoAnalysis&4) isSelAnTopoCutsKKpi=kTRUE;
            if(isSelectedTopoAnalysis&8) isSelAnTopoCutspiKK=kTRUE;
          }
          else if(fWriteVariableTreeDs==2) {
            if(isSelectedAnalysis&16) isSelAnCutsKKpi=kTRUE;
            if(isSelectedAnalysis&32) isSelAnCutspiKK=kTRUE;
            if(isSelectedPidAnalysis==1 || isSelectedPidAnalysis==3) isSelAnPidCutsKKpi=kTRUE;
            if(isSelectedPidAnalysis==2 || isSelectedPidAnalysis==3) isSelAnPidCutspiKK=kTRUE;
            if(isSelectedTopoAnalysis&16) isSelAnTopoCutsKKpi=kTRUE;
            if(isSelectedTopoAnalysis&32) isSelAnTopoCutspiKK=kTRUE;
          }
          else if(fWriteVariableTreeDs==3) {
            if(isSelectedAnalysis&1) isSelAnCutsKKpi=kTRUE;
            if(isSelectedAnalysis&2) isSelAnCutspiKK=kTRUE;
            if(isSelectedPidAnalysis==1 || isSelectedPidAnalysis==3) isSelAnPidCutsKKpi=kTRUE;
            if(isSelectedPidAnalysis==2 || isSelectedPidAnalysis==3) isSelAnPidCutspiKK=kTRUE;
            if(isSelectedTopoAnalysis&1) isSelAnTopoCutsKKpi=kTRUE;
            if(isSelectedTopoAnalysis&2) isSelAnTopoCutspiKK=kTRUE;
          }
          
          Bool_t unsetvtx=kFALSE;
          if(!ds->GetOwnPrimaryVtx()){
            ds->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsDstoKKpi->GetIsPrimaryWithoutDaughters()){
            if(ds->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*ds->GetOwnPrimaryVtx());
            if(fFiltCutsDstoKKpi->RecalcOwnPrimaryVtx(ds,aod))recVtx=kTRUE;
            else fFiltCutsDstoKKpi->CleanOwnPrimaryVtx(ds,aod,origownvtx);
          }
          
          Int_t labDs=-1;
          Int_t labDplus=-1;
          Int_t pdgCode0=-999;
          Int_t orig=0;
          Float_t ptGenDs = -99.;
          //checking origin
          AliAODMCParticle *partDs = 0x0;
          if(fReadMC){
            labDs = ds->MatchToMC(431,arrMC,3,pdgDstoKKpi);
            labDplus = ds->MatchToMC(411,arrMC,3,pdgDstoKKpi);
            
            if(labDs>=0){
              Int_t labDau0=((AliAODTrack*)ds->GetDaughter(0))->GetLabel();
              AliAODMCParticle* p=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDau0));
              pdgCode0=TMath::Abs(p->GetPdgCode());
              partDs = (AliAODMCParticle*)arrMC->At(labDs);
              ptGenDs = partDs->Pt();
            }
            else{
              if(labDplus>=0) {
                Int_t labDau0=((AliAODTrack*)ds->GetDaughter(0))->GetLabel();
                AliAODMCParticle* p=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDau0));
                pdgCode0=TMath::Abs(p->GetPdgCode());
                partDs = (AliAODMCParticle*)arrMC->At(labDplus);
                ptGenDs = partDs->Pt();
              }
            }
            if(partDs) orig = AliVertexingHFUtils::CheckOrigin(arrMC,partDs,!fITSUpgradeProduction);
          }
          
          //filling the Ds tree
          if ((fWriteVariableTreeDs==1 && (isPhiKKpi || isPhipiKK)) || (fWriteVariableTreeDs==2 && (isK0starKKpi || isK0starpiKK)) || (fWriteVariableTreeDs==3 && (isKKpi || ispiKK))){
            
            bool issignal = kFALSE;
            bool isbkg = kFALSE;
            bool isprompt = kFALSE;
            bool isFD = kFALSE;
            bool isrefl = kFALSE;
            
            if((fWriteVariableTreeDs==3 && isKKpi) || (fWriteVariableTreeDs==1 && isPhiKKpi) || (fWriteVariableTreeDs==2 && isK0starKKpi)) {
              if(fReadMC) {
                if(labDs>=0) {
                  if(orig==4 || orig==5) {
                    if(pdgCode0==321) issignal = kTRUE;
                    else if(pdgCode0==211) isrefl = kTRUE;
                    if(orig==4) isprompt = kTRUE;
                    else if(orig==5) isFD = kTRUE;
                  }
                }
                else {
                  if(fStoreOnlyHIJINGBackground){
                    Bool_t isHijing = kFALSE;
                    if (mcHeader) isHijing = IsCandidateFromHijing(ds,mcHeader,arrMC);
                    if (isHijing) isbkg = kTRUE;
                  } else {
                    isbkg=kTRUE;
                  }
                  if(labDplus>=0) fTreeHandlerDs->SetIsDplustoKKpi(kTRUE);//put also D+ -->KKpi in bkg
                }
                //do not apply cuts, but enable flag if is selected
                if(issignal || isbkg || isrefl) fTreeHandlerDs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
              }
              if(!fReadMC || (issignal || isbkg || isrefl)) {
                fTreeHandlerDs->SetIsSelectedStd(isSelAnCutsKKpi,isSelAnTopoCutsKKpi,isSelAnPidCutsKKpi,isSelTracksAnCuts);
                fTreeHandlerDs->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenDs,ds,bfield,0,fPIDresp);
                if (fFillJets) fTreeHandlerDs->SetJetVars(aod->GetTracks(),ds,ds->InvMassDsKKpi(),arrMC,partDs);
                fTreeHandlerDs->FillTree();
              }
            }
            issignal = kFALSE;
            isbkg = kFALSE;
            isprompt = kFALSE;
            isFD = kFALSE;
            isrefl = kFALSE;
            if((fWriteVariableTreeDs==3 && ispiKK) || (fWriteVariableTreeDs==1 && isPhipiKK) || (fWriteVariableTreeDs==2 && isK0starpiKK)) {
              if(fReadMC) {
                if(labDs>=0) {
                  if(orig==4 || orig==5) {
                    if(pdgCode0==211) issignal = kTRUE;
                    else if(pdgCode0==321) isrefl = kTRUE;
                    if(orig==4) isprompt = kTRUE;
                    else if(orig==5) isFD = kTRUE;
                  }
                }
                else {
                  if(fStoreOnlyHIJINGBackground){
                    Bool_t isHijing = kFALSE;
                    if (mcHeader) isHijing = IsCandidateFromHijing(ds,mcHeader,arrMC);
                    if (isHijing) isbkg = kTRUE;
                  } else {
                    isbkg=kTRUE;
                  }
                  if(labDplus>=0) fTreeHandlerDs->SetIsDplustoKKpi(kTRUE);//put also D+ -->KKpi in bkg
                }
                //do not apply cuts, but enable flag if is selected
                if(issignal || isbkg || isrefl) fTreeHandlerDs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
              }
              if(!fReadMC || (issignal || isbkg || isrefl)) {
                fTreeHandlerDs->SetIsSelectedStd(isSelAnCutspiKK,isSelAnTopoCutspiKK,isSelAnPidCutspiKK,isSelTracksAnCuts);
                fTreeHandlerDs->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenDs,ds,bfield,1,fPIDresp);
                if (fFillJets) fTreeHandlerDs->SetJetVars(aod->GetTracks(),ds,ds->InvMassDspiKK(),arrMC,partDs);
                fTreeHandlerDs->FillTree();
              }
            }
          }//end fill tree
          if(recVtx)fFiltCutsDstoKKpi->CleanOwnPrimaryVtx(ds,aod,origownvtx);
          if(unsetvtx) ds->UnsetOwnPrimaryVtx();
        }//end is selected
      }
      else{
        fNentries->Fill(18); //monitor how often this fails
      }
    }//end Ds
    
    
    //*************************************
    //Dplus
    Bool_t isDplustagged=kTRUE;
    AliAODRecoDecayHF3Prong *dplus = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3prong);
    if(fUseSelectionBit && !(dplus->HasSelectionBit(AliRDHFCuts::kDplusCuts))){
      isDplustagged=kFALSE;
    }
    if(isDplustagged && fWriteVariableTreeDplus){

      Int_t preSelectedDplus = -1;
      if(fITSUpgradePreSelect){
        TObjArray arrTracks(3);
        for(Int_t ipr=0;ipr<3;ipr++){
          AliAODTrack *tr=vHF->GetProng(aod,dplus,ipr);
          arrTracks.AddAt(tr,ipr);
        }
        preSelectedDplus = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 3, 411, pdgDgDplustoKpipi);
        if(preSelectedDplus == 0) continue; //Mixture hijing + injected
        if(preSelectedDplus == 2) continue; //Only MatchedToMC injected signal
        if(fWriteOnlySignal == 1 && preSelectedDplus != 1) continue; //Only matched signal when only signal is requested
      }

      nFilteredDplus++;
      fNentries->Fill(19);
      if((vHF->FillRecoCand(aod,dplus))) {////Fill the data members of the candidate only if they are empty.
        
        Int_t isSelectedFilt = fFiltCutsDplustoKpipi->IsSelected(dplus,AliRDHFCuts::kAll,aod);
        
        if(isSelectedFilt){
          fNentries->Fill(20);
          nSelectedDplus++;
          
          //test analysis cuts
          Bool_t isSelAnCuts = kFALSE;
          Bool_t isSelAnPidCuts = kFALSE;
          Bool_t isSelAnTopolCuts = kFALSE;
          Int_t isSelectedAnalysis = fCutsDplustoKpipi->IsSelected(dplus,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis = fCutsDplustoKpipi->IsSelectedPID(dplus);
          bool isUsePidAn = fCutsDplustoKpipi->GetIsUsePID();
          if(isUsePidAn) fCutsDplustoKpipi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsDplustoKpipi->IsSelected(dplus,AliRDHFCuts::kAll,aod);
          if(isSelectedAnalysis) isSelAnCuts = kTRUE;
          if(isSelectedPidAnalysis) isSelAnPidCuts = kTRUE;
          if(isSelectedTopoAnalysis) isSelAnTopolCuts = kTRUE;
          fCutsDplustoKpipi->SetUsePID(isUsePidAn);
          Bool_t isSelTracksAnCuts=kFALSE;
          Int_t isSelectedTrackAnalysis = fCutsDplustoKpipi->IsSelected(dplus,AliRDHFCuts::kTracks,aod);
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          Bool_t unsetvtx=kFALSE;
          if(!dplus->GetOwnPrimaryVtx()){
            dplus->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsDplustoKpipi->GetIsPrimaryWithoutDaughters()){
            if(dplus->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*dplus->GetOwnPrimaryVtx());
            if(fFiltCutsDplustoKpipi->RecalcOwnPrimaryVtx(dplus,aod))recVtx=kTRUE;
            else fFiltCutsDplustoKpipi->CleanOwnPrimaryVtx(dplus,aod,origownvtx);
          }
          
          Int_t labDp=-1;
          bool isPrimary=kFALSE;
          bool isFeeddown=kFALSE;
          bool issignal=kFALSE;
          bool isbkg=kFALSE;
          Int_t pdgCode=-2;
          Float_t ptGenDplus = -99.;
          //read MC
          AliAODMCParticle *partDp=0x0;
          if(fReadMC){
            labDp = dplus->MatchToMC(411,arrMC,3,pdgDgDplustoKpipi);
            if(labDp>=0){
              partDp = (AliAODMCParticle*)arrMC->At(labDp);
              ptGenDplus = partDp->Pt();
              Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partDp,!fITSUpgradeProduction);//Prompt = 4, FeedDown = 5
              if(orig==4 || orig==5) {
                issignal=kTRUE;
                pdgCode=TMath::Abs(partDp->GetPdgCode());
                if(orig==4){
                  isPrimary=kTRUE;
                  isFeeddown=kFALSE;
                }
                else if(orig==5){
                  isPrimary=kFALSE;
                  isFeeddown=kTRUE;
                }
              }
            }
            else {
              if(fStoreOnlyHIJINGBackground){
                Bool_t isHijing = kFALSE;
                if (mcHeader) isHijing = IsCandidateFromHijing(dplus,mcHeader,arrMC);
                if (isHijing) isbkg = kTRUE;
              } else {
                isbkg=kTRUE;
              }
            }
            if(issignal || isbkg) fTreeHandlerDplus->SetCandidateType(issignal,isbkg,isPrimary,isFeeddown,kFALSE);
          } //end read MC
          
          // fill tree
          if(!fReadMC || (issignal || isbkg)) {
            fTreeHandlerDplus->SetIsSelectedStd(isSelAnCuts,isSelAnTopolCuts,isSelAnPidCuts,isSelTracksAnCuts);
            fTreeHandlerDplus->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenDplus,dplus,bfield,0,fPIDresp);
            if (fFillJets) fTreeHandlerDplus->SetJetVars(aod->GetTracks(),dplus,dplus->InvMassDplus(),arrMC,partDp);
            fTreeHandlerDplus->FillTree();
          }
          //end fill tree
          if(recVtx)fFiltCutsDplustoKpipi->CleanOwnPrimaryVtx(dplus,aod,origownvtx);
          if(unsetvtx) dplus->UnsetOwnPrimaryVtx();
        } //end topol and PID cuts
      }//end ok fill reco cand
      else{
        fNentries->Fill(21); //monitor how often this fails
      }
    }//end Dplus
    
    
    //*************************************
    //LctopKpi
    Bool_t isLctopKpitagged=kTRUE;
    AliAODRecoDecayHF3Prong *lctopkpi = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3prong);
    if(fUseSelectionBit && !(lctopkpi->HasSelectionBit(AliRDHFCuts::kLcCuts))){
      isLctopKpitagged=kFALSE;
    }
    if(isLctopKpitagged && fWriteVariableTreeLctopKpi){

      TObjArray arrTracks(3);
      if(fITSUpgradePreSelect || fPreSelectLctopKpi){
        for(Int_t ipr=0;ipr<3;ipr++){
          AliAODTrack *tr=vHF->GetProng(aod,lctopkpi,ipr);
          arrTracks.AddAt(tr,ipr);
        }
      }

      Int_t preSelectedLc = -1;
      if(fITSUpgradePreSelect){
        preSelectedLc = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 3, 4122, pdgLctopKpi);
        if(preSelectedLc == 0) continue; //Mixture hijing + injected
        if(preSelectedLc == 2) continue; //Only MatchedToMC injected signal
        if(fWriteOnlySignal == 1 && preSelectedLc != 1) continue; //Only matched signal when only signal is requested
      }

      if(fPreSelectLctopKpi){
        Bool_t preSelectedLcMass=kTRUE;
        preSelectedLcMass = fFiltCutsLctopKpi->PreSelectMass(arrTracks);
        if (!preSelectedLcMass) continue;
      }

      nFilteredLctopKpi++;
      fNentries->Fill(22);
      if((vHF->FillRecoCand(aod,lctopkpi))) {////Fill the data members of the candidate only if they are empty.
        
        Int_t isSelectedFilt    = fFiltCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);

        if(isSelectedFilt){
          fNentries->Fill(23);
          nSelectedLctopKpi++;
          
          // check analysis cuts
          Bool_t isSelAnCutspKpi=kFALSE;
          Bool_t isSelAnCutspiKp=kFALSE;
          Bool_t isSelPIDpKpi=kFALSE;
          Bool_t isSelPIDpiKp=kFALSE;
          Bool_t isSelTopopKpi=kFALSE;
          Bool_t isSelTopopiKp=kFALSE;
          Bool_t ispKpi=kFALSE;
          Bool_t ispiKp=kFALSE;
          Int_t isSelectedAnalysis= fCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis = fCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kPID,aod);
          Bool_t isUsePidAn = fCutsLctopKpi->GetIsUsePID();
          if(isUsePidAn) fCutsLctopKpi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
          fCutsLctopKpi->SetUsePID(isUsePidAn);
          Bool_t isSelTracksAnCuts=kFALSE;
          Int_t isSelectedTrackAnalysis = fCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kTracks,aod);
          if(isSelectedTrackAnalysis > 0)                            isSelTracksAnCuts=kTRUE;
          if(isSelectedAnalysis==1 || isSelectedAnalysis==3)         isSelAnCutspKpi=kTRUE;
          if(isSelectedAnalysis>=2)                                  isSelAnCutspiKp=kTRUE;
          if(isSelectedTopoAnalysis==1 || isSelectedTopoAnalysis==3) isSelTopopKpi=kTRUE;
          if(isSelectedTopoAnalysis>=2)                              isSelTopopiKp=kTRUE;
          if(isSelectedPidAnalysis==1 || isSelectedPidAnalysis==3)   isSelPIDpKpi=kTRUE;
          if(isSelectedPidAnalysis>=2)                               isSelPIDpiKp=kTRUE;
          if(isSelectedFilt==1 || isSelectedFilt==3)                 ispKpi=kTRUE;
          if(isSelectedFilt>=2)                                      ispiKp=kTRUE;
          
          Bool_t unsetvtx=kFALSE;
          if(!lctopkpi->GetOwnPrimaryVtx()){
            lctopkpi->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsLctopKpi->GetIsPrimaryWithoutDaughters()){
            if(lctopkpi->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*lctopkpi->GetOwnPrimaryVtx());
            if(fFiltCutsLctopKpi->RecalcOwnPrimaryVtx(lctopkpi,aod))recVtx=kTRUE;
            else fFiltCutsLctopKpi->CleanOwnPrimaryVtx(lctopkpi,aod,origownvtx);
          }

          bool isPrimary=kFALSE;
          bool isFeeddown=kFALSE;
          bool issignal=kFALSE;
          bool isbkg=kFALSE;
          bool isrefl=kFALSE;
          Int_t labDp=-1;
          Float_t ptGenLcpKpi = -99.;
          int restype = -1;
          if(ispKpi) {
            //read MC
            AliAODMCParticle *partDp=0x0;
            if(fReadMC){
              labDp = lctopkpi->MatchToMC(4122,arrMC,3,pdgLctopKpi);
              if(labDp>=0){
                partDp = (AliAODMCParticle*)arrMC->At(labDp);
                ptGenLcpKpi = partDp->Pt();
                Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partDp,!fITSUpgradeProduction);//Prompt = 4, FeedDown = 5
                if(orig==4 || orig==5) {
                  issignal=kTRUE;
                  if(orig==4){
                    isPrimary=kTRUE;
                    isFeeddown=kFALSE;
                  }
                  else if(orig==5){
                    isPrimary=kFALSE;
                    isFeeddown=kTRUE;
                  }
                }
                //check daughters
                Int_t labDauLc0=((AliAODTrack*)lctopkpi->GetDaughter(0))->GetLabel();
                Int_t labDauLc1=((AliAODTrack*)lctopkpi->GetDaughter(1))->GetLabel();
                Int_t labDauLc2=((AliAODTrack*)lctopkpi->GetDaughter(2))->GetLabel();
                AliAODMCParticle* pDauLc0=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDauLc0));
                AliAODMCParticle* pDauLc1=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDauLc1));
                AliAODMCParticle* pDauLc2=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDauLc2));
                Int_t pdgDauLc0=TMath::Abs(pDauLc0->GetPdgCode());
                Int_t pdgDauLc1=TMath::Abs(pDauLc1->GetPdgCode());
                Int_t pdgDauLc2=TMath::Abs(pDauLc2->GetPdgCode());
                if(pdgDauLc0==211 && pdgDauLc1==321 && pdgDauLc2==2212) isrefl=kTRUE;
                restype = fTreeHandlerLctopKpi->GetLcResonantDecay(arrMC,partDp);
              }
              else {
                if(fStoreOnlyHIJINGBackground){
                  Bool_t isHijing = kFALSE;
                  if (mcHeader) isHijing = IsCandidateFromHijing(lctopkpi,mcHeader,arrMC);
                  if (isHijing) isbkg = kTRUE;
                } else {
                  isbkg=kTRUE;
                }
              }
              if(issignal || isbkg || isrefl) fTreeHandlerLctopKpi->SetCandidateType(issignal,isbkg,isPrimary,isFeeddown,isrefl);
            } //end read MC
            
            // fill tree
            if(!fReadMC || (issignal || isbkg || isrefl)) {
              fTreeHandlerLctopKpi->SetIsSelectedStd(isSelAnCutspKpi,isSelTopopKpi,isSelPIDpKpi,isSelTracksAnCuts);
              fTreeHandlerLctopKpi->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenLcpKpi,lctopkpi,bfield,1,fPIDresp);
              fTreeHandlerLctopKpi->SetVariableResonantDecay(restype);
              if (fFillJets) fTreeHandlerLctopKpi->SetJetVars(aod->GetTracks(),lctopkpi,lctopkpi->InvMassLcpKpi(),arrMC,partDp);
              fTreeHandlerLctopKpi->FillTree();
            }
          } // end pKpi
          isPrimary=kFALSE;
          isFeeddown=kFALSE;
          issignal=kFALSE;
          isbkg=kFALSE;
          isrefl=kFALSE;
          labDp=-1;
          restype = -1;
          if(ispiKp) {
            //read MC
            AliAODMCParticle *partDp=0x0;
            if(fReadMC){
              labDp = lctopkpi->MatchToMC(4122,arrMC,3,pdgLctopKpi);
              if(labDp>=0){
                partDp = (AliAODMCParticle*)arrMC->At(labDp);
                ptGenLcpKpi = partDp->Pt();
                Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partDp,!fITSUpgradeProduction);//Prompt = 4, FeedDown = 5
                if(orig==4 || orig==5) {
                  issignal=kTRUE;
                  if(orig==4){
                    isPrimary=kTRUE;
                    isFeeddown=kFALSE;
                  }
                  else if(orig==5){
                    isPrimary=kFALSE;
                    isFeeddown=kTRUE;
                  }
                }
                //check daughters
                Int_t labDauLc0=((AliAODTrack*)lctopkpi->GetDaughter(0))->GetLabel();
                Int_t labDauLc1=((AliAODTrack*)lctopkpi->GetDaughter(1))->GetLabel();
                Int_t labDauLc2=((AliAODTrack*)lctopkpi->GetDaughter(2))->GetLabel();
                AliAODMCParticle* pDauLc0=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDauLc0));
                AliAODMCParticle* pDauLc1=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDauLc1));
                AliAODMCParticle* pDauLc2=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDauLc2));
                Int_t pdgDauLc0=TMath::Abs(pDauLc0->GetPdgCode());
                Int_t pdgDauLc1=TMath::Abs(pDauLc1->GetPdgCode());
                Int_t pdgDauLc2=TMath::Abs(pDauLc2->GetPdgCode());
                if(pdgDauLc0==2212 && pdgDauLc1==321 && pdgDauLc2==211) isrefl=kTRUE;
                restype = fTreeHandlerLctopKpi->GetLcResonantDecay(arrMC,partDp);
              }
              else {
                if(fStoreOnlyHIJINGBackground){
                  Bool_t isHijing = kFALSE;
                  if (mcHeader) isHijing = IsCandidateFromHijing(lctopkpi,mcHeader,arrMC);
                  if (isHijing) isbkg = kTRUE;
                } else {
                  isbkg=kTRUE;
                }
              }
              if(issignal || isbkg || isrefl) fTreeHandlerLctopKpi->SetCandidateType(issignal,isbkg,isPrimary,isFeeddown,isrefl);
            } //end read MC
            
            // fill tree
            if(!fReadMC || (issignal || isbkg || isrefl)) {
              fTreeHandlerLctopKpi->SetIsSelectedStd(isSelAnCutspiKp,isSelTopopiKp,isSelPIDpiKp,isSelTracksAnCuts);
              fTreeHandlerLctopKpi->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenLcpKpi,lctopkpi,bfield,2,fPIDresp);
              fTreeHandlerLctopKpi->SetVariableResonantDecay(restype);
              if (fFillJets) fTreeHandlerLctopKpi->SetJetVars(aod->GetTracks(),lctopkpi,lctopkpi->InvMassLcpiKp(),arrMC,partDp);
              fTreeHandlerLctopKpi->FillTree();
            }
          } // end fill piKpi
          if(recVtx)fFiltCutsLctopKpi->CleanOwnPrimaryVtx(lctopkpi,aod,origownvtx);
          if(unsetvtx) lctopkpi->UnsetOwnPrimaryVtx();
        } //end topol and PID cuts
      }//end ok fill reco cand
      else{
        fNentries->Fill(24); //monitor how often this fails
      }
    }//end LctopKpi
    
    
    
    
  }//end loop on cadidates
  
  delete vHF;
  return;
}

//________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::ProcessDstar(TClonesArray *arrayDstar, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield){
  
  if(fStoreOnlyHIJINGBackground) AliInfo("Storing HIJING background only not implemented for Dstar");

  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t nCasc = arrayDstar->GetEntriesFast();
  if(fDebug>2) printf("Number of D*-> D0 pi: %d\n",nCasc);
  
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t nSelectedDstar=0;
  Int_t nFilteredDstar=0;    
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  
  for (Int_t iCasc = 0; iCasc < nCasc; iCasc++) {
    fNentries->Fill(28);
    
    //Dstar
    Bool_t isDstartagged=kTRUE;
    AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*)arrayDstar->UncheckedAt(iCasc);
    if(fUseSelectionBit && d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kDstarCuts)){
      isDstartagged=kFALSE;
    }
    
    if (isDstartagged && fWriteVariableTreeDstar){
      
      fNentries->Fill(29);
      nFilteredDstar++;
      if((vHF->FillRecoCasc(aod,d,kTRUE))) {//Fill the data members of the candidate only if they are empty.
        
        Int_t isSelectedFilt = fFiltCutsDstartoKpipi->IsSelected(d,AliRDHFCuts::kAll,aod); //selected
        if(isSelectedFilt > 0){
          fNentries->Fill(30);
          nSelectedDstar++;
          
          //test analysis cuts
          Bool_t isSelAnCuts = kFALSE;
          Bool_t isSelAnPidCuts = kFALSE;
          Bool_t isSelAnTopolCuts = kFALSE;
          Int_t isSelectedAnalysis = fCutsDstartoKpipi->IsSelected(d,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis = fCutsDstartoKpipi->IsSelectedPID(d);
          Bool_t isUsePidAn = fCutsDstartoKpipi->GetIsUsePID();
          if(isUsePidAn) fCutsDstartoKpipi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsDstartoKpipi->IsSelected(d,AliRDHFCuts::kAll,aod);
          Bool_t isSelTracksAnCuts=kFALSE;
          Int_t isSelectedTrackAnalysis = fCutsDstartoKpipi->IsSelected(d,AliRDHFCuts::kTracks,aod);
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          if(isSelectedAnalysis) isSelAnCuts = kTRUE;
          if(isSelectedPidAnalysis) isSelAnPidCuts = kTRUE;
          if(isSelectedTopoAnalysis) isSelAnTopolCuts = kTRUE;
          fCutsDstartoKpipi->SetUsePID(isUsePidAn);
          
          //For Dstar the removal of daughters is done only for D0 prongs.
          AliAODRecoDecayHF2Prong* dd = (AliAODRecoDecayHF2Prong*)d->Get2Prong();
          Bool_t unsetvtx=kFALSE;
          if(!dd->GetOwnPrimaryVtx()){
            dd->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsDstartoKpipi->GetIsPrimaryWithoutDaughters()){
            if(dd->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*dd->GetOwnPrimaryVtx());
            if(fFiltCutsDstartoKpipi->RecalcOwnPrimaryVtx(dd,aod))recVtx=kTRUE;
            else fFiltCutsDstartoKpipi->CleanOwnPrimaryVtx(dd,aod,origownvtx);
          }

          Int_t labDstar = -1;
          Int_t pdgDstar = -99;
          Int_t origin= -1;
          Float_t ptGenDstar = -99;
          
          AliAODMCParticle *partDstar=0x0;
          if(fReadMC) {
            labDstar = d->MatchToMC(413,421,pdgDgDStartoD0pi, pdgDgD0toKpi, arrMC);
            if(labDstar>=0){
              partDstar = (AliAODMCParticle*)arrMC->At(labDstar);
              pdgDstar = TMath::Abs(partDstar->GetPdgCode());
              ptGenDstar = partDstar->Pt();
              origin = AliVertexingHFUtils::CheckOrigin(arrMC,partDstar,!fITSUpgradeProduction);
            }
          }
          
          bool issignal = kFALSE;
          bool isbkg =    kFALSE;
          bool isFD =     kFALSE;
          bool isprompt = kFALSE;
          bool isrefl =   kFALSE;
          Int_t masshypo = 0;
          
          if(fReadMC){
            if(labDstar>=0){
              if(origin==4 || origin==5) {
                if(origin==4) isprompt=kTRUE;
                else if(origin==5) isFD=kTRUE;
                if(pdgDstar==413){
                  issignal=kTRUE;
                }
              }
            }//end labDstar check
            else{//background
              isbkg=kTRUE;
            }
            if(issignal || isbkg || isrefl) fTreeHandlerDstar->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
          }//end read MC
          if(!fReadMC || (issignal || isbkg || isrefl)) {
            fTreeHandlerDstar->SetIsSelectedStd(isSelAnCuts,isSelAnTopolCuts,isSelAnPidCuts,isSelTracksAnCuts);
            fTreeHandlerDstar->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenDstar,d,bfield,masshypo,fPIDresp);
            if (fFillJets) fTreeHandlerDstar->SetJetVars(aod->GetTracks(),d,d->InvMassDstarKpipi(),arrMC,partDstar);	
            fTreeHandlerDstar->FillTree();
          }
          if(recVtx)fFiltCutsDstartoKpipi->CleanOwnPrimaryVtx(dd,aod,origownvtx);
          if(unsetvtx) dd->UnsetOwnPrimaryVtx();
        }//end is selected filt
      }
      else {
        fNentries->Fill(31); //monitor how often this fails
      }
    }//end Dstar
  }//end loop on candidates
  
  delete vHF;
  return;
}
//________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::ProcessCasc(TClonesArray *arrayCasc, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield){
  
  if(fStoreOnlyHIJINGBackground) AliInfo("Storing HIJING background only not implemented for Lc->pK0s");

  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t nCasc = arrayCasc->GetEntriesFast();
  if(fDebug>2) printf("Number of Cascade: %d\n",nCasc);
  
  Int_t pdgDgLc2K0Spr[2]={2212,310};
  Int_t pdgDgK0stoDaughters[2]={211,211};
  Int_t nSelectedLc2V0bachelor=0;
  Int_t nFilteredLc2V0bachelor=0;
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  
  for (Int_t iCasc = 0; iCasc < nCasc; iCasc++) {
    fNentries->Fill(32);
    
    //Lc2V0bachelor
    Bool_t isLc2V0bachelortagged=kTRUE;
    AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*)arrayCasc->UncheckedAt(iCasc);
    if(fUseSelectionBit && d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kLctoV0Cuts)){
      isLc2V0bachelortagged=kFALSE;
    }
    
    if (isLc2V0bachelortagged && fWriteVariableTreeLc2V0bachelor){
      
      AliAODv0 * v0part;
      if(d->GetIsFilled() == 0) v0part = (AliAODv0*)(aod->GetV0(d->GetProngID(1)));
      else                      v0part = d->Getv0();
      Bool_t isOnFlyV0 = v0part->GetOnFlyStatus();

      if(fV0typeForLc2V0bachelor==1 && isOnFlyV0 == kTRUE) continue;
      if(fV0typeForLc2V0bachelor==2 && isOnFlyV0 == kFALSE) continue;

      if(fFiltCutsLc2V0bachelor->GetUsePreselect() && d->GetIsFilled() == 0){
        TObjArray arrTracks(2);
        for(Int_t ipr=0;ipr<2;ipr++){
          AliAODTrack *tr;
          if(ipr==0) tr=vHF->GetProng(aod,d,ipr);
          else tr = (AliAODTrack*)(aod->GetV0(d->GetProngID(1)));
          arrTracks.AddAt(tr,ipr);
        }
        Int_t PreSelectLc = fFiltCutsLc2V0bachelor->PreSelect(arrTracks);
        if(PreSelectLc==0) continue;
      }

      fNentries->Fill(33);
      nFilteredLc2V0bachelor++;
      if((vHF->FillRecoCasc(aod,d,kFALSE,fLc2V0bachelorCalcSecoVtx))) {//Fill the data members of the candidate only if they are empty.
        
        //To calculate secondary vertex for pp/pPb if requested
        //Remember to run also CleanUpTask!
        if(d->GetIsFilled()==1 && fLc2V0bachelorCalcSecoVtx) vHF->RecoSecondaryVertexForCascades(aod, d);
        
        Int_t isSelectedFilt = fFiltCutsLc2V0bachelor->IsSelected(d,AliRDHFCuts::kAll,aod); //selected
        if(isSelectedFilt > 0){
          fNentries->Fill(34);
          nSelectedLc2V0bachelor++;
          
          //test analysis cuts
          Bool_t isSelAnCutstopK0s = kFALSE;
          Bool_t isSelAnCutstoLpi = kFALSE;
          Bool_t isSelAnPidCutstopK0s = kFALSE;
          Bool_t isSelAnPidCutstoLpi = kFALSE;
          Bool_t isSelAnTopolCutstopK0s = kFALSE;
          Bool_t isSelAnTopolCutstoLpi = kFALSE;
          Int_t isSelectedAnalysis = fCutsLc2V0bachelor->IsSelected(d,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis = fCutsLc2V0bachelor->IsSelectedPID(d);
          Bool_t isUsePidAn = fCutsLc2V0bachelor->GetIsUsePID();
          if(isUsePidAn) fCutsLc2V0bachelor->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsLc2V0bachelor->IsSelected(d,AliRDHFCuts::kAll,aod);
          
          //Standard selection Lc->pK0s, but keep also Lc->Lpi (different bit)
          if( (isSelectedAnalysis&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr))                                                                                                       isSelAnCutstopK0s = kTRUE;
          if( ((isSelectedAnalysis&(AliRDHFCutsLctoV0::kLcToLpi)) == (AliRDHFCutsLctoV0::kLcToLpi)) || ((isSelectedAnalysis&(AliRDHFCutsLctoV0::kLcToLBarpi)) == (AliRDHFCutsLctoV0::kLcToLBarpi)) )         isSelAnCutstoLpi = kTRUE;
          if( (isSelectedPidAnalysis&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr))                                                                                                    isSelAnPidCutstopK0s = kTRUE;
          if( ((isSelectedPidAnalysis&(AliRDHFCutsLctoV0::kLcToLpi)) == (AliRDHFCutsLctoV0::kLcToLpi)) || ((isSelectedPidAnalysis&(AliRDHFCutsLctoV0::kLcToLBarpi)) == (AliRDHFCutsLctoV0::kLcToLBarpi)))    isSelAnPidCutstoLpi = kTRUE;
          if( (isSelectedTopoAnalysis&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr))                                                                                                   isSelAnTopolCutstopK0s = kTRUE;
          if( ((isSelectedTopoAnalysis&(AliRDHFCutsLctoV0::kLcToLpi)) == (AliRDHFCutsLctoV0::kLcToLpi)) || ((isSelectedTopoAnalysis&(AliRDHFCutsLctoV0::kLcToLBarpi)) == (AliRDHFCutsLctoV0::kLcToLBarpi)) ) isSelAnTopolCutstoLpi = kTRUE;
          
          fCutsLc2V0bachelor->SetUsePID(isUsePidAn);
          Bool_t isSelTracksAnCuts=kFALSE;
          Int_t isSelectedTrackAnalysis = fCutsLc2V0bachelor->IsSelected(d,AliRDHFCuts::kTracks,aod);
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          Bool_t unsetvtx=kFALSE;
          if(!d->GetOwnPrimaryVtx()){
            d->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsLc2V0bachelor->GetIsPrimaryWithoutDaughters()){
            if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
            if(fFiltCutsLc2V0bachelor->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
            else fFiltCutsLc2V0bachelor->CleanOwnPrimaryVtx(d,aod,origownvtx);
          }

          Int_t labLc2V0bachelor = -1;
          Int_t pdgLc2V0bachelor = -99;
          Int_t origin= -1;
          Float_t ptGenLc2V0bachelor = -99;
          
          AliAODMCParticle *partLc2V0bachelor=0x0;
          if(fReadMC) {
            labLc2V0bachelor = d->MatchToMC(4122,310,pdgDgLc2K0Spr, pdgDgK0stoDaughters, arrMC, kTRUE);
            if(labLc2V0bachelor>=0){
              partLc2V0bachelor = (AliAODMCParticle*)arrMC->At(labLc2V0bachelor);
              pdgLc2V0bachelor = TMath::Abs(partLc2V0bachelor->GetPdgCode());
              ptGenLc2V0bachelor = partLc2V0bachelor->Pt();
              origin = AliVertexingHFUtils::CheckOrigin(arrMC,partLc2V0bachelor,!fITSUpgradeProduction);
            }
          }
          
          bool issignal = kFALSE;
          bool isbkg =    kFALSE;
          bool isFD =     kFALSE;
          bool isprompt = kFALSE;
          bool isrefl =   kFALSE;
          Int_t masshypo = 0;
          
          if(fReadMC){
            if(labLc2V0bachelor>=0){
              if(origin==4 || origin==5) {
                if(origin==4) isprompt=kTRUE;
                else if(origin==5) isFD=kTRUE;
                if(pdgLc2V0bachelor==4122){
                  issignal=kTRUE;
                }
              }
            }//end labLc2V0bachelor check
            else{//background
              isbkg=kTRUE;
            }
            if(issignal || isbkg || isrefl) fTreeHandlerLc2V0bachelor->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
          }//end read MC
          if(!fReadMC || (issignal || isbkg || isrefl)) {
            fTreeHandlerLc2V0bachelor->SetIsSelectedStd(isSelAnCutstopK0s,isSelAnTopolCutstopK0s,isSelAnPidCutstopK0s,isSelTracksAnCuts);
            fTreeHandlerLc2V0bachelor->SetIsLctoLpi(isSelAnCutstoLpi, isSelAnTopolCutstoLpi, isSelAnPidCutstoLpi);
            fTreeHandlerLc2V0bachelor->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenLc2V0bachelor,d,bfield,masshypo,fPIDresp);
            if (fFillJets) fTreeHandlerLc2V0bachelor->SetJetVars(aod->GetTracks(),d,d->InvMassLctoK0sP(),arrMC,partLc2V0bachelor);
            fTreeHandlerLc2V0bachelor->FillTree();
          }
          if(recVtx)fFiltCutsLc2V0bachelor->CleanOwnPrimaryVtx(d,aod,origownvtx);
          if(unsetvtx) d->UnsetOwnPrimaryVtx();
        }//end is selected filt
      }
      else {
        fNentries->Fill(35); //monitor how often this fails
      }
    }//end Lc2V0bachelor
  }//end loop on candidates
  
  delete vHF;
  return;
}
//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::ProcessBplus(TClonesArray *array2prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t n2prong = array2prong->GetEntriesFast();
  if(fDebug>2) printf("Number of 2prongs: %d\n",n2prong);
  
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t pdgDgD0topiK[2]={211,321};
  Int_t nSelectedD0=0;
  Int_t nFilteredD0=0;
  
  Int_t pdgDgBplustoD0piInt[2] = {211,421};
  UInt_t pdgDgBplustoD0piUInt[2] = {211,421};
  Int_t nSelectedBplus = 0;
  Int_t nFilteredBplus = 0;
  
  // Select good track before hand to save time
  Int_t nTracks= aod->GetNumberOfTracks();
  if (nTracks==0) return;
  Bool_t seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectGoodTrackForReconstruction(aod, nTracks, nSeleTrks, seleTrkFlags);
  if(nSeleTrks<1) return;
  
  // Get AliAODPidHF Bplus object for IsBsPionSelected function in BplusTreeHandler
  AliAODPidHF* fPidHFD0 = fFiltCutsBplustoD0pi->GetPidHF();

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  
  for (Int_t i2prong = 0; i2prong < n2prong; i2prong++) {
    
    //D0 from Bplus
    Bool_t isD0fromBplustagged=kTRUE;
    AliAODRecoDecayHF2Prong *dfromB = (AliAODRecoDecayHF2Prong*)array2prong->UncheckedAt(i2prong);
    if(fUseSelectionBit && dfromB->GetSelectionMap()) if(!dfromB->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
      isD0fromBplustagged=kFALSE;
    }

    Int_t preSelectedBplus = -1;
    if(fITSUpgradePreSelect){
      TObjArray arrTracks(2);
      for(Int_t ipr=0;ipr<2;ipr++){
        AliAODTrack *tr=vHF->GetProng(aod,dfromB,ipr);
        arrTracks.AddAt(tr,ipr);
      }
      preSelectedBplus = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 2, 421, pdgDgD0toKpi);
      if(preSelectedBplus == 0) continue; //Mixture hijing + injected
      if(preSelectedBplus == 2) continue; //Only MatchedToMC injected signal
      if(fWriteOnlySignal == 1 && preSelectedBplus != 1) continue; //Only matched signal when only signal is requested
    }

    if(isD0fromBplustagged && fWriteVariableTreeBplus){
      nFilteredD0++;
      if ((vHF->FillRecoCand(aod, dfromB))) { //Fill the data members of the candidate only if they are empty.

        //To significantly speed up the task when only signal was requested
        if(fWriteOnlySignal && !fITSUpgradePreSelect){
          if(dfromB->MatchToMC(421,arrMC,2,pdgDgD0toKpi) < 0) continue;
        }
        
        Int_t isSelectedFilt=fFiltCutsBplustoD0pi->IsSelected(dfromB,AliRDHFCuts::kAll,aod);
        if(isSelectedFilt > 0){
          fNentries->Fill(25);
          nSelectedD0++;

          //test analysis cuts
          Bool_t isSelAnCuts=kFALSE;
          Bool_t isSelAnPidCuts=kFALSE;
          Bool_t isSelAnTopoCuts=kFALSE;
          Bool_t isSelTracksAnCuts=kFALSE;

          //analysis cuts
          Int_t isSelectedAnalysis = fCutsBplustoD0pi->IsSelected(dfromB,AliRDHFCuts::kAll,aod); //selected
          Int_t isSelectedPidAnalysis = fCutsBplustoD0pi->IsSelectedPID(dfromB); //selected
          bool isUsePidAn = fCutsBplustoD0pi->GetIsUsePID();
          if(isUsePidAn) fCutsBplustoD0pi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsBplustoD0pi->IsSelected(dfromB,AliRDHFCuts::kAll,aod);
          fCutsBplustoD0pi->SetUsePID(isUsePidAn);
          Int_t isSelectedTrackAnalysis = fCutsBplustoD0pi->IsSelected(dfromB,AliRDHFCuts::kTracks,aod);
          
          if(isSelectedAnalysis > 0) isSelAnCuts=kTRUE;
          if(isSelectedPidAnalysis > 0) isSelAnPidCuts=kTRUE;
          if(isSelectedTopoAnalysis > 0) isSelAnTopoCuts=kTRUE;
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          Bool_t unsetvtx=kFALSE;
          if(!dfromB->GetOwnPrimaryVtx()){
            dfromB->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsBplustoD0pi->GetIsPrimaryWithoutDaughters()){
            if(dfromB->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*dfromB->GetOwnPrimaryVtx());
            if(fFiltCutsBplustoD0pi->RecalcOwnPrimaryVtx(dfromB,aod))recVtx=kTRUE;
            else fFiltCutsBplustoD0pi->CleanOwnPrimaryVtx(dfromB,aod,origownvtx);
          }
          
          //loop over all selected tracks for pion from Bplus
          for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
            if(!seleTrkFlags[iTrack]) continue;

            AliAODTrack* pionTrack = dynamic_cast<AliAODTrack*>(aod->GetTrack(iTrack));
            if (!pionTrack) AliFatal("Not a standard AOD");

            Bool_t trackInjected = kFALSE;
            if(fReadMC){
              Bool_t trackInjected = AliVertexingHFUtils::IsTrackInjected(pionTrack,mcHeader,arrMC);
              if(preSelectedBplus==3 && trackInjected) continue;
            }

            nFilteredBplus++;
            if(fTreeHandlerBplus->IsBplusPionSelected(pionTrack, fFiltCutsBplustoD0pi, fPidHFD0, aod, vtx1)){

              //we check if the IDs of the tracks are different
              AliAODTrack* twoProngdaughter0 = (AliAODTrack*)dfromB->GetDaughter(0);
              AliAODTrack* twoProngdaughter1 = (AliAODTrack*)dfromB->GetDaughter(1);
              UShort_t idProng0 = twoProngdaughter0->GetID();
              UShort_t idProng1 = twoProngdaughter1->GetID();
              if (pionTrack->GetID() != idProng0 && pionTrack->GetID() != idProng1){
                
                //we use the Bplus pion and D0 tracks to reconstruct the vertex for the Bplus
                //(order swapped wrt to Bs and Lb because of already existing MatchToMC function)
                AliExternalTrackParam firstTrack;
                firstTrack.CopyFromVTrack(pionTrack);
                AliExternalTrackParam secondTrack;
                secondTrack.CopyFromVTrack(dfromB);
                
                // we calculate the vertex of the mother candidate
                TObjArray daughterTracks;
                daughterTracks.Add(&firstTrack);
                daughterTracks.Add(&secondTrack);
                Double_t dispersion = 0;
                AliAODVertex *vertexMother = ReconstructDisplVertex(vtx1, &daughterTracks, bfield, dispersion);
                
                if (vertexMother){ //check if calculation vertex Bplus succeeded
                  
                  //use the new vertex to create the Bplus candidate
                  Double_t xdummy = 0., ydummy = 0.;
                  Double_t d0z0[2], covd0z0[3], d0[2], d0err[2];
                  
                  firstTrack.PropagateToDCA(vertexMother, bfield, 100., d0z0, covd0z0);
                  secondTrack.PropagateToDCA(vertexMother, bfield, 100., d0z0, covd0z0);
                  
                  //we reconstruct the mother decay prong
                  Double_t px[2], py[2], pz[2];
                  px[0] = firstTrack.Px();
                  py[0] = firstTrack.Py();
                  pz[0] = firstTrack.Pz();
                  px[1] = secondTrack.Px();
                  py[1] = secondTrack.Py();
                  pz[1] = secondTrack.Pz();
                  
                  UShort_t id[2];
                  id[0] = firstTrack.GetID();
                  id[1] = secondTrack.GetID();
                  
                  firstTrack.PropagateToDCA(vtx1, bfield, 100., d0z0, covd0z0);
                  d0[0] = d0z0[0];
                  d0err[0] = TMath::Sqrt(covd0z0[0]);
                  secondTrack.PropagateToDCA(vtx1, bfield, 100., d0z0, covd0z0);
                  d0[1] = d0z0[0];
                  d0err[1] = TMath::Sqrt(covd0z0[0]);
                  
                  Double_t dca = secondTrack.GetDCA(&firstTrack, bfield, xdummy, ydummy);
                  Short_t chargeMother = dfromB->Charge() + pionTrack->Charge();
                  if(chargeMother == 0) AliWarning("Bplus does not have charge, please check!");
                  
                  //Using filtering cuts, too many AliAODRecoDecay objects are built per event
                  //The maximum number of TRef for a given TProcesssID is 2^24=16777216 can be reached.
                  //Save current Object count, and set it back after filling the TTree to not get AliFatal
                  Int_t ObjectNumber = TProcessID::GetObjectCount();
                  AliAODRecoDecayHF2Prong trackBplus(vertexMother, px, py, pz, d0, d0err, dca);
                  
                  trackBplus.SetCharge(chargeMother);
                  trackBplus.GetSecondaryVtx()->AddDaughter(pionTrack);
                  trackBplus.GetSecondaryVtx()->AddDaughter(dfromB);
                  trackBplus.SetPrimaryVtxRef((AliAODVertex*)aod->GetPrimaryVertex());
                  trackBplus.SetProngIDs(2, id);
                  
                  /*Add some hardcoded cuts Bplus object for the moment. To be moved to TreeHandler when there are more*/
                  Double_t invmassBplus = trackBplus.InvMass(2, pdgDgBplustoD0piUInt);
                  Double_t massBplusPDG = TDatabasePDG::Instance()->GetParticle(521)->Mass();
                  if(TMath::Abs(invmassBplus-massBplusPDG) < 0.3){
                    
                    fNentries->Fill(26);
                    nSelectedBplus++;
                  
                    Int_t labBplus = -1;
                    bool issignal = kFALSE;
                    bool isbkg =    kFALSE;
                    bool isFD =     kFALSE;
                    bool isprompt = kTRUE; //beauty, so "always" prompt
                    bool isrefl =   kFALSE; //To check, do we have reflections?
                    Float_t ptGenBplus = -99.;
                    //read mc
                    AliAODMCParticle *partBplus=0x0;
                    if (fReadMC){
                      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
                      //Effect is per mille level, but big enough to be rejected by MatchToMC() function.
                      //Fix implemented, loosening cut from 0.001% to 2.0%. If > 2.0%, negative label is returned.
                      labBplus = trackBplus.MatchToMCB2Prong(521,421,pdgDgBplustoD0piInt,pdgDgD0topiK,arrMC);
                      
                      //Use feed-down bit to flag these candidates for offline study.
                      //Can be good for analysis or coming from an "incomplete simulated" decay
                      if(labBplus < -1){ isFD=kTRUE; labBplus=TMath::Abs(labBplus); }

                      if(labBplus >= 0) {
                        partBplus = (AliAODMCParticle*)arrMC->At(labBplus);
                        ptGenBplus = partBplus->Pt();
                        issignal = kTRUE;
                      } else {
                        if(fStoreOnlyHIJINGBackground){
                          Bool_t isHijing = kFALSE;
                          if (mcHeader) isHijing = IsCandidateFromHijing(dfromB,mcHeader,arrMC,pionTrack);
                          if (isHijing) isbkg = kTRUE;
                        } else {
                          isbkg=kTRUE;
                        }
                      }
                      
                      if(issignal || isbkg || isrefl) fTreeHandlerBplus->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                    } //end read MC

                    // fill tree
                    if(!fReadMC || (issignal || isbkg || isrefl)) {
                      fTreeHandlerBplus->SetIsSelectedStd(isSelAnCuts,isSelAnTopoCuts,isSelAnPidCuts,isSelTracksAnCuts);
                      fTreeHandlerBplus->SetVariables(fRunNumber, fEventID, fEventIDExt, fEventIDLong, ptGenBplus, &trackBplus, bfield, 0, fPIDresp);
                      if (fFillJets) fTreeHandlerBplus->SetJetVars(aod->GetTracks(),&trackBplus,trackBplus.InvMass(2,pdgDgBplustoD0piUInt),arrMC,partBplus);
                      fTreeHandlerBplus->FillTree();
                    }
                  }//end hardcoded cuts Bplus object
                  //Restore Object count, to save space in the table keeping track of all referenced objects
                  TProcessID::SetObjectCount(ObjectNumber);
                }//end vertex Bplus check
                delete vertexMother; vertexMother = nullptr;
              }//end check pion ID with D0 daughters
            }//end Bplus pion filtering cuts
          }//end track-loop
          if(recVtx)fFiltCutsBplustoD0pi->CleanOwnPrimaryVtx(dfromB,aod,origownvtx);
          if(unsetvtx) dfromB->UnsetOwnPrimaryVtx();
        }//end D0 filtering cuts
      } else {
        fNentries->Fill(27); //monitor how often this fails
      }
    }//end Bplus
  }//end loop on candidates
  
  delete vHF;
  return;
}
//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::ProcessBs(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t n3prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of 3prongs: %d\n",n3prong);
  
  Int_t pdgDstoKKpi[3]={321,321,211};
  Int_t nSelectedDs=0;
  Int_t nFilteredDs=0;
  
  Int_t pdgDgBstoDspi[2] = {431,211};
  UInt_t pdgDgBstoDspiUInt[2] = {431,211};
  Int_t nSelectedBs = 0;
  Int_t nFilteredBs = 0;
  
  // Select good track before hand to save time
  Int_t nTracks= aod->GetNumberOfTracks();
  if (nTracks==0) return;
  Bool_t seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectGoodTrackForReconstruction(aod, nTracks, nSeleTrks, seleTrkFlags);
  if(nSeleTrks<1) return;
  
  // Get AliAODPidHF Bs object for IsBsPionSelected function in BsTreeHandler
  AliAODPidHF* fPidHFDs = fFiltCutsBstoDspi->GetPidHF();

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  
  for (Int_t i3prong = 0; i3prong < n3prong; i3prong++) {
    
    //Ds
    Bool_t isDstagged=kTRUE;
    AliAODRecoDecayHF3Prong *ds    = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3prong);
    if(fUseSelectionBit && !(ds->HasSelectionBit(AliRDHFCuts::kDsCuts))){
      isDstagged=kFALSE;
    }

    TObjArray arrTracks(3);
    if(fITSUpgradePreSelect || fFiltCutsBstoDspi->GetUsePreselect()){
      for(Int_t ipr=0;ipr<3;ipr++){
        AliAODTrack *tr=vHF->GetProng(aod,ds,ipr);
        arrTracks.AddAt(tr,ipr);
      }
    }
    Int_t preSelectedBs = -1;
    if(fITSUpgradePreSelect){
      Int_t preSelectedBs = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 3, 431, pdgDstoKKpi);
      if(preSelectedBs == 0) continue; //Mixture hijing + injected
      if(preSelectedBs == 2) continue; //Only MatchedToMC injected signal
      if(fWriteOnlySignal == 1 && preSelectedBs != 1) continue; //Only matched signal when only signal is requested
    }
    Int_t preSelectedDsCuts = -1;
    if(fFiltCutsBstoDspi->GetUsePreselect() && ds->GetIsFilled() == 0){
      preSelectedDsCuts = fFiltCutsBstoDspi->PreSelect(arrTracks);
      if(preSelectedDsCuts==0) continue;
    }

    if(isDstagged && fWriteVariableTreeBs){
      nFilteredDs++;
      if((vHF->FillRecoCand(aod,ds))) {////Fill the data members of the candidate only if they are empty.

        //To significantly speed up the task when only signal was requested
        if(fWriteOnlySignal && !fITSUpgradePreSelect){
          if(ds->MatchToMC(431,arrMC,3,pdgDstoKKpi) < 0) continue;
        }

        //Only looking at Ds -> phi pi -> KK pi.
        Int_t isSelectedFilt=fFiltCutsBstoDspi->IsSelected(ds,AliRDHFCuts::kAll,aod);
        Int_t isPhiKKpi=isSelectedFilt&4;
        Int_t isPhipiKK=isSelectedFilt&8;

        if(isPhiKKpi || isPhipiKK){
          fNentries->Fill(36);
          nSelectedDs++;
          
          //Set mass hypothesis for Ds, for filling TTree in correct way
          //Assuming here that candidate never satisfies isPhiKKpi and isPhipiKK (due to tight Ds mass cut)
          Int_t masshypoDs = -1;
          if(isPhiKKpi) masshypoDs=0;
          else          masshypoDs=1;

          //test analysis cuts
          Bool_t isSelAnCuts=kFALSE;
          Bool_t isSelAnPidCuts=kFALSE;
          Bool_t isSelAnTopoCuts=kFALSE;
          Bool_t isSelTracksAnCuts=kFALSE;
          
          Int_t isSelectedAnalysis=fCutsBstoDspi->IsSelected(ds,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis=fCutsBstoDspi->IsSelectedPID(ds);
          bool isUsePidAn = fCutsBstoDspi->GetIsUsePID();
          if(isUsePidAn) fCutsBstoDspi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsBstoDspi->IsSelected(ds,AliRDHFCuts::kAll,aod);
          fCutsBstoDspi->SetUsePID(isUsePidAn);
          Int_t isSelectedTrackAnalysis = fCutsBstoDspi->IsSelected(ds,AliRDHFCuts::kTracks,aod);
          
          //See Ds (Process3Prong) for all different options, for now just accept all of them
          if(isSelectedAnalysis > 0) isSelAnCuts=kTRUE;
          if(isSelectedPidAnalysis > 0) isSelAnPidCuts=kTRUE;
          if(isSelectedTopoAnalysis > 0) isSelAnTopoCuts=kTRUE;
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          Bool_t unsetvtx=kFALSE;
          if(!ds->GetOwnPrimaryVtx()){
            ds->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsBstoDspi->GetIsPrimaryWithoutDaughters()){
            if(ds->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*ds->GetOwnPrimaryVtx());
            if(fFiltCutsBstoDspi->RecalcOwnPrimaryVtx(ds,aod))recVtx=kTRUE;
            else fFiltCutsBstoDspi->CleanOwnPrimaryVtx(ds,aod,origownvtx);
          }
          
          //loop over all selected tracks for pion from Bs
          for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
            if(!seleTrkFlags[iTrack]) continue;
            
            AliAODTrack* pionTrack = dynamic_cast<AliAODTrack*>(aod->GetTrack(iTrack));
            if (!pionTrack) AliFatal("Not a standard AOD");
            if (ds->Charge() != -1.*pionTrack->Charge()) continue;

            Bool_t trackInjected = kFALSE;
            if(fReadMC){
              trackInjected = AliVertexingHFUtils::IsTrackInjected(pionTrack,mcHeader,arrMC);
              if(preSelectedBs==3 && trackInjected) continue;
            }

            nFilteredBs++;
            if(fTreeHandlerBs->IsBsPionSelected(pionTrack, fFiltCutsBstoDspi, fPidHFDs, aod, vtx1)){

              //we check if the IDs of the tracks are different
              AliAODTrack* threeProngdaughter0 = (AliAODTrack*)ds->GetDaughter(0);
              AliAODTrack* threeProngdaughter1 = (AliAODTrack*)ds->GetDaughter(1);
              AliAODTrack* threeProngdaughter2 = (AliAODTrack*)ds->GetDaughter(2);
              UShort_t idProng0 = threeProngdaughter0->GetID();
              UShort_t idProng1 = threeProngdaughter1->GetID();
              UShort_t idProng2 = threeProngdaughter2->GetID();
              if (pionTrack->GetID() != idProng0 && pionTrack->GetID() != idProng1 && pionTrack->GetID() != idProng2){
                
                //we use the Bs pion and Ds tracks to reconstruct the vertex for the Bs
                AliExternalTrackParam firstTrack;
                firstTrack.CopyFromVTrack(ds);
                AliExternalTrackParam secondTrack;
                secondTrack.CopyFromVTrack(pionTrack);
                
                // we calculate the vertex of the mother candidate
                TObjArray daughterTracks;
                daughterTracks.Add(&firstTrack);
                daughterTracks.Add(&secondTrack);
                Double_t dispersion = 0;
                AliAODVertex *vertexMother = ReconstructDisplVertex(vtx1, &daughterTracks, bfield, dispersion);
                
                //check if calculation vertex Bs succeeded
                if (vertexMother){
                  
                  //use the new vertex to create the Bs candidate
                  Double_t xdummy = 0., ydummy = 0.;
                  Double_t d0z0[2], covd0z0[3], d0[2], d0err[2];
                  
                  firstTrack.PropagateToDCA(vertexMother, bfield, 100., d0z0, covd0z0);
                  secondTrack.PropagateToDCA(vertexMother, bfield, 100., d0z0, covd0z0);
                  
                  //we reconstruct the mother decay prong
                  Double_t px[2], py[2], pz[2];
                  px[0] = firstTrack.Px();
                  py[0] = firstTrack.Py();
                  pz[0] = firstTrack.Pz();
                  px[1] = secondTrack.Px();
                  py[1] = secondTrack.Py();
                  pz[1] = secondTrack.Pz();
                  
                  UShort_t id[2];
                  id[0] = firstTrack.GetID();
                  id[1] = secondTrack.GetID();
                  
                  firstTrack.PropagateToDCA(vtx1, bfield, 100., d0z0, covd0z0);
                  d0[0] = d0z0[0];
                  d0err[0] = TMath::Sqrt(covd0z0[0]);
                  secondTrack.PropagateToDCA(vtx1, bfield, 100., d0z0, covd0z0);
                  d0[1] = d0z0[0];
                  d0err[1] = TMath::Sqrt(covd0z0[0]);
                  
                  Double_t dca = secondTrack.GetDCA(&firstTrack, bfield, xdummy, ydummy);
                  Short_t chargeMother = ds->Charge() + pionTrack->Charge();
                  if(chargeMother != 0) AliWarning("Bs0 got charge, please check!");
                  
                  //Using filtering cuts, too many AliAODRecoDecay objects are built per event
                  //The maximum number of TRef for a given TProcesssID is 2^24=16777216 can be reached.
                  //Save current Object count, and set it back after filling the TTree to not get AliFatal
                  Int_t ObjectNumber = TProcessID::GetObjectCount();
                  AliAODRecoDecayHF2Prong trackBs(vertexMother, px, py, pz, d0, d0err, dca);
                  
                  trackBs.SetCharge(chargeMother);
                  trackBs.GetSecondaryVtx()->AddDaughter(ds);
                  trackBs.GetSecondaryVtx()->AddDaughter(pionTrack);
                  trackBs.SetPrimaryVtxRef((AliAODVertex*)aod->GetPrimaryVertex());
                  trackBs.SetProngIDs(2, id);
                  
                  if(fTreeHandlerBs->IsBsSelected(&trackBs)){
                    fNentries->Fill(37);
                    nSelectedBs++;
                    
                    Int_t labBs = -1;
                    bool issignal = kFALSE;
                    bool isbkg =    kFALSE;
                    bool isFD =     kFALSE;
                    bool isprompt = kTRUE; //beauty, so "always" prompt
                    bool isrefl =   kFALSE;
                    Float_t ptGenBs = -99.;

                    //for storing injected candidate + HIJING track for background shape studies
                    //NB: using reflection bit to get these candidates saved for offline study.
                    bool isDsPrompt  = kFALSE;
                    bool isDsFDBplus = kFALSE;
                    bool isDsFDB0    = kFALSE;
                    bool isDsFDLb0   = kFALSE;
                    bool isDsFDBs0   = kFALSE;

                    //read mc
                    AliAODMCParticle *partBs=0x0;
                    if (fReadMC){
                      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
                      //Effect is per mille level, but big enough to be rejected by MatchToMC() function.
                      //Fix implemented, loosening cut from 0.001% to 2.0%. If > 2.0%, negative label is returned.
                      labBs = trackBs.MatchToMCB3Prong(531,431,pdgDgBstoDspi,pdgDstoKKpi,arrMC);

                      //Use feed-down bit to flag these candidates for offline study.
                      //Can be good for analysis or coming from an "incomplete simulated" decay
                      if(labBs < -1){ isFD=kTRUE; labBs=TMath::Abs(labBs); }

                      if(labBs >= 0) {
                        partBs = (AliAODMCParticle*)arrMC->At(labBs);
                        ptGenBs = partBs->Pt();
                        Bool_t isHijing = IsCandidateFromHijing(ds,mcHeader,arrMC,pionTrack);;
                        if (!isHijing) issignal = kTRUE;
                      } else {
                        if(fStoreOnlyHIJINGBackground || fFillInjCandHijingTrackCombi){
                          Bool_t isHijing=kFALSE;
                          if (mcHeader){
                            isHijing= IsCandidateFromHijing(ds,mcHeader,arrMC,pionTrack);
                            if(!isHijing && !trackInjected && fFillInjCandHijingTrackCombi){
                              // Store injected Ds and HIJING pion track for background shape studies
                              Int_t labDs = ds->MatchToMC(431,arrMC,3,pdgDstoKKpi);
                              if(labDs >= 0){
                                AliAODMCParticle *partDs = (AliAODMCParticle*)arrMC->At(labDs);
                                Int_t orig = AliVertexingHFUtils::CheckOrigin(arrMC,partDs,!fITSUpgradeProduction);
                                if(orig == 4) isDsPrompt = kTRUE;
                                if(orig == 5){
                                  Int_t labMother = partDs->GetMother();
                                  if(labMother >= 0){
                                    AliAODMCParticle *partMother = (AliAODMCParticle*)arrMC->At(labMother);
                                    Int_t pdgMother = TMath::Abs(partMother->GetPdgCode());
                                    if(pdgMother == 511) isDsFDB0 = kTRUE;
                                    if(pdgMother == 521) isDsFDBplus = kTRUE;
                                    if(pdgMother == 531) isDsFDBs0 = kTRUE;
                                    if(pdgMother == 5122) isDsFDLb0 = kTRUE;
                                  }
                                }
                              }
                            }
                          }
                          if (isHijing) isbkg = kTRUE;
                          if (isDsPrompt || isDsFDB0 || isDsFDBplus || isDsFDBs0 || isDsFDLb0) isrefl = kTRUE;
                        } else {
                          isbkg=kTRUE;
                        }
                      }
                      
                      if(issignal || isbkg || isrefl) fTreeHandlerBs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                    } //end read MC

                    // fill tree
                    if(!fReadMC || (issignal || isbkg || isrefl)) {
                      fTreeHandlerBs->SetIsSelectedStd(isSelAnCuts,isSelAnTopoCuts,isSelAnPidCuts,isSelTracksAnCuts);
                      fTreeHandlerBs->SetDsBackgroundShapeType(isDsPrompt, isDsFDBplus, isDsFDB0, isDsFDLb0, isDsFDBs0);
                      fTreeHandlerBs->SetVariables(fRunNumber, fEventID, fEventIDExt, fEventIDLong, ptGenBs, &trackBs, bfield, masshypoDs, fPIDresp);
                      if (fFillJets) fTreeHandlerBs->SetJetVars(aod->GetTracks(),&trackBs,trackBs.InvMass(2,pdgDgBstoDspiUInt),arrMC,partBs);
                      fTreeHandlerBs->FillTree();
                    }
                  }//end hardcoded cuts Bs object
                  //Restore Object count, to save space in the table keeping track of all referenced objects
                  TProcessID::SetObjectCount(ObjectNumber);
                }//end vertex Bs check
                delete vertexMother; vertexMother = nullptr;
              }//end check pion ID with Ds daughters
            }//end Bs pion filtering cuts
          }//end track-loop
          if(recVtx)fFiltCutsBstoDspi->CleanOwnPrimaryVtx(ds,aod,origownvtx);
          if(unsetvtx) ds->UnsetOwnPrimaryVtx();
        }//end Ds filtering cuts
      } else{
        fNentries->Fill(38); //monitor how often this fails
      }
    }//end Ds-loop
  }//end loop on cadidates
  
  delete vHF;
  return;
}
//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::ProcessLb(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t n3prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of 3prongs: %d\n",n3prong);
  
  Int_t pdgLctopKpi[3]={2212,321,211};
  Int_t nSelectedLctopKpi=0;
  Int_t nFilteredLctopKpi=0;
  
  Int_t pdgDgLbtoLcpi[2] = {4122,211};
  UInt_t pdgDgLbtoLcpiUInt[2] = {4122,211};
  Int_t nSelectedLb = 0;
  Int_t nFilteredLb = 0;
  
  // Select good track before hand to save time
  Int_t nTracks= aod->GetNumberOfTracks();
  if (nTracks==0) return;
  Bool_t seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectGoodTrackForReconstruction(aod, nTracks, nSeleTrks, seleTrkFlags);
  if(nSeleTrks<1) return;
  
  // Get AliAODPidHF Lb object for IsLbPionSelected function in LbTreeHandler
  AliAODPidHF* fPidHFLc = fFiltCutsLbtoLcpi->GetPidHF();
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  
  for (Int_t i3prong = 0; i3prong < n3prong; i3prong++) {
    
    //Lc->pKpi
    Bool_t isLctopKpitagged=kTRUE;
    AliAODRecoDecayHF3Prong *lctopkpi = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3prong);
    if(fUseSelectionBit && !(lctopkpi->HasSelectionBit(AliRDHFCuts::kLcCuts))){
      isLctopKpitagged=kFALSE;
    }

    TObjArray arrTracks(3);
    if(fITSUpgradePreSelect || fPreSelectLctopKpi){
      for(Int_t ipr=0;ipr<3;ipr++){
        AliAODTrack *tr=vHF->GetProng(aod,lctopkpi,ipr);
        arrTracks.AddAt(tr,ipr);
      }
    }

    Int_t preSelectedLb = -1;
    if(fITSUpgradePreSelect){
      Int_t preSelectedLb = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 3, 4122, pdgLctopKpi);
      if(preSelectedLb == 0) continue; //Mixture hijing + injected
      if(preSelectedLb == 2) continue; //Only MatchedToMC injected signal
      if(fWriteOnlySignal == 1 && preSelectedLb != 1) continue; //Only matched signal when only signal is requested
    }

    if(fPreSelectLctopKpi){
      Bool_t preSelectedLbMass=kTRUE;
      preSelectedLbMass=fFiltCutsLbtoLcpi->PreSelectMass(arrTracks);
      if (!preSelectedLbMass) continue;
    }

    if(isLctopKpitagged && fWriteVariableTreeLb){
      nFilteredLctopKpi++;

      if((vHF->FillRecoCand(aod,lctopkpi))) {////Fill the data members of the candidate only if they are empty.

        //To significantly speed up the task when only signal was requested
        if(fWriteOnlySignal && !fITSUpgradePreSelect){
          if(lctopkpi->MatchToMC(4122,arrMC,3,pdgLctopKpi) < 0) continue;
        }
        
        //Only looking at Lc -> p K pi
        Int_t isSelectedFilt = fFiltCutsLbtoLcpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
        if(isSelectedFilt){
          fNentries->Fill(39);
          nSelectedLctopKpi++;
          
          //Set mass hypothesis for Lc, for filling TTree in correct way
          //Assuming here that candidate never satisfies pKpi and piKp (due to tight Lc mass cut)
          Int_t masshypoLc = -1;
          if(isSelectedFilt==1) masshypoLc=1;
          else                  masshypoLc=2;
          
          //test analysis cuts
          Bool_t isSelAnCuts=kFALSE;
          Bool_t isSelAnPidCuts=kFALSE;
          Bool_t isSelAnTopoCuts=kFALSE;
          Bool_t isSelTracksAnCuts=kFALSE;
          
          Int_t isSelectedAnalysis=fCutsLbtoLcpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
          Int_t isSelectedPidAnalysis=fCutsLbtoLcpi->IsSelectedPID(lctopkpi);
          bool isUsePidAn = fCutsLbtoLcpi->GetIsUsePID();
          if(isUsePidAn) fCutsLbtoLcpi->SetUsePID(kFALSE);
          Int_t isSelectedTopoAnalysis = fCutsLbtoLcpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
          fCutsLbtoLcpi->SetUsePID(isUsePidAn);
          Int_t isSelectedTrackAnalysis = fCutsLbtoLcpi->IsSelected(lctopkpi,AliRDHFCuts::kTracks,aod);
          
          //See Lc (Process3Prong) for the different options, for now just accept all of them
          if(isSelectedAnalysis > 0) isSelAnCuts=kTRUE;
          if(isSelectedPidAnalysis > 0) isSelAnPidCuts=kTRUE;
          if(isSelectedTopoAnalysis > 0) isSelAnTopoCuts=kTRUE;
          if(isSelectedTrackAnalysis > 0) isSelTracksAnCuts=kTRUE;
          
          Bool_t unsetvtx=kFALSE;
          if(!lctopkpi->GetOwnPrimaryVtx()){
            lctopkpi->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
            // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
            // Pay attention if you use continue inside this loop!!!
          }
          Bool_t recVtx=kFALSE;
          AliAODVertex *origownvtx=0x0;
          if(fFiltCutsLbtoLcpi->GetIsPrimaryWithoutDaughters()){
            if(lctopkpi->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*lctopkpi->GetOwnPrimaryVtx());
            if(fFiltCutsLbtoLcpi->RecalcOwnPrimaryVtx(lctopkpi,aod))recVtx=kTRUE;
            else fFiltCutsLbtoLcpi->CleanOwnPrimaryVtx(lctopkpi,aod,origownvtx);
          }

          //loop over all selected tracks for pion from Lb
          for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
            if(!seleTrkFlags[iTrack]) continue;

            AliAODTrack* pionTrack = dynamic_cast<AliAODTrack*>(aod->GetTrack(iTrack));
            if (!pionTrack) AliFatal("Not a standard AOD");
            if (lctopkpi->Charge() != -1.*pionTrack->Charge()) continue;

            Bool_t trackInjected = kFALSE;
            if(fReadMC){
              trackInjected = AliVertexingHFUtils::IsTrackInjected(pionTrack,mcHeader,arrMC);
              if(preSelectedLb==3 && trackInjected) continue;
            }

            nFilteredLb++;
            if(fTreeHandlerLb->IsLbPionSelected(pionTrack, fFiltCutsLbtoLcpi, fPidHFLc, aod, vtx1)){

              //we check if the IDs of the tracks are different
              AliAODTrack* threeProngdaughter0 = (AliAODTrack*)lctopkpi->GetDaughter(0);
              AliAODTrack* threeProngdaughter1 = (AliAODTrack*)lctopkpi->GetDaughter(1);
              AliAODTrack* threeProngdaughter2 = (AliAODTrack*)lctopkpi->GetDaughter(2);
              UShort_t idProng0 = threeProngdaughter0->GetID();
              UShort_t idProng1 = threeProngdaughter1->GetID();
              UShort_t idProng2 = threeProngdaughter2->GetID();
              if (pionTrack->GetID() != idProng0 && pionTrack->GetID() != idProng1 && pionTrack->GetID() != idProng2){
                
                //we use the Lb pion and Lc tracks to reconstruct the vertex for the Lb
                AliExternalTrackParam firstTrack;
                firstTrack.CopyFromVTrack(lctopkpi);
                AliExternalTrackParam secondTrack;
                secondTrack.CopyFromVTrack(pionTrack);
                
                // we calculate the vertex of the mother candidate
                TObjArray daughterTracks;
                daughterTracks.Add(&firstTrack);
                daughterTracks.Add(&secondTrack);
                Double_t dispersion = 0;
                AliAODVertex *vertexMother = ReconstructDisplVertex(vtx1, &daughterTracks, bfield, dispersion);
                
                //check if calculation vertex Lb succeeded
                if (vertexMother){
                  
                  //use the new vertex to create the Lb candidate
                  Double_t xdummy = 0., ydummy = 0.;
                  Double_t d0z0[2], covd0z0[3], d0[2], d0err[2];
                  
                  firstTrack.PropagateToDCA(vertexMother, bfield, 100., d0z0, covd0z0);
                  secondTrack.PropagateToDCA(vertexMother, bfield, 100., d0z0, covd0z0);
                  
                  //we reconstruct the mother decay prong
                  Double_t px[2], py[2], pz[2];
                  px[0] = firstTrack.Px();
                  py[0] = firstTrack.Py();
                  pz[0] = firstTrack.Pz();
                  px[1] = secondTrack.Px();
                  py[1] = secondTrack.Py();
                  pz[1] = secondTrack.Pz();
                  
                  UShort_t id[2];
                  id[0] = firstTrack.GetID();
                  id[1] = secondTrack.GetID();
                  
                  firstTrack.PropagateToDCA(vtx1, bfield, 100., d0z0, covd0z0);
                  d0[0] = d0z0[0];
                  d0err[0] = TMath::Sqrt(covd0z0[0]);
                  secondTrack.PropagateToDCA(vtx1, bfield, 100., d0z0, covd0z0);
                  d0[1] = d0z0[0];
                  d0err[1] = TMath::Sqrt(covd0z0[0]);
                  
                  Double_t dca = secondTrack.GetDCA(&firstTrack, bfield, xdummy, ydummy);
                  Short_t chargeMother = lctopkpi->Charge() + pionTrack->Charge();
                  if(chargeMother != 0) AliWarning("Lb0 got charge, please check!");
                  
                  //Using filtering cuts, too many AliAODRecoDecay objects are built per event
                  //The maximum number of TRef for a given TProcesssID is 2^24=16777216 can be reached.
                  //Save current Object count, and set it back after filling the TTree to not get AliFatal
                  Int_t ObjectNumber = TProcessID::GetObjectCount();
                  AliAODRecoDecayHF2Prong trackLb(vertexMother, px, py, pz, d0, d0err, dca);
                  
                  trackLb.SetCharge(chargeMother);
                  trackLb.GetSecondaryVtx()->AddDaughter(lctopkpi);
                  trackLb.GetSecondaryVtx()->AddDaughter(pionTrack);
                  trackLb.SetPrimaryVtxRef((AliAODVertex*)aod->GetPrimaryVertex());
                  trackLb.SetProngIDs(2, id);
                  
                  /*Add some hardcoded cuts Lb object for the moment. To be moved to TreeHandler when there are more*/
                  Double_t invmassLb = trackLb.InvMass(2, pdgDgLbtoLcpiUInt);
                  Double_t massLbPDG = TDatabasePDG::Instance()->GetParticle(5122)->Mass();
                  if(TMath::Abs(invmassLb-massLbPDG) < 0.3){
                    
                    fNentries->Fill(40);
                    nSelectedLb++;
                    
                    Int_t labLb = -1;
                    bool issignal = kFALSE;
                    bool isbkg =    kFALSE;
                    bool isFD =     kFALSE;
                    bool isprompt = kTRUE; //beauty, so "always" prompt
                    bool isrefl =   kFALSE; //To check, do we have reflections?
                    Float_t ptGenLb = -99.;
                    //read mc
                    AliAODMCParticle *partLb=0x0;
                    if (fReadMC){
                      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
                      //Effect is per mille level, but big enough to be rejected by MatchToMC() function.
                      //Fix implemented, loosening cut from 0.001% to 2.0%. If > 2.0%, negative label is returned.
                      labLb = trackLb.MatchToMCB3Prong(5122,4122,pdgDgLbtoLcpi,pdgLctopKpi,arrMC);
                      
                      //Use feed-down bit to flag these candidates for offline study.
                      //Can be good for analysis or coming from an "incomplete simulated" decay
                      if(labLb < -1){ isFD=kTRUE; labLb=TMath::Abs(labLb); }

                      if(labLb >= 0) {
                        partLb = (AliAODMCParticle*)arrMC->At(labLb);
                        ptGenLb = partLb->Pt();
                        issignal = kTRUE;
                      } else{
                        if(fStoreOnlyHIJINGBackground){
                          Bool_t isHijing = kFALSE;
                          if (mcHeader) isHijing = IsCandidateFromHijing(lctopkpi,mcHeader,arrMC,pionTrack);
                          if (isHijing) isbkg = kTRUE;
                        } else {
                          isbkg=kTRUE;
                        }
                      }
                      
                      if(issignal || isbkg || isrefl) fTreeHandlerLb->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                    } //end read MC
                    
                    // fill tree
                    if(!fReadMC || (issignal || isbkg || isrefl)) {
                      fTreeHandlerLb->SetIsSelectedStd(isSelAnCuts,isSelAnTopoCuts,isSelAnPidCuts,isSelTracksAnCuts);
                      fTreeHandlerLb->SetVariables(fRunNumber, fEventID, fEventIDExt, fEventIDLong, ptGenLb, &trackLb, bfield, masshypoLc, fPIDresp);
                      if (fFillJets) fTreeHandlerLb->SetJetVars(aod->GetTracks(),&trackLb,trackLb.InvMass(2,pdgDgLbtoLcpiUInt),arrMC,partLb);
                      fTreeHandlerLb->FillTree();
                    }
                  }//end hardcoded cuts Lb object
                  //Restore Object count, to save space in the table keeping track of all referenced objects
                  TProcessID::SetObjectCount(ObjectNumber);
                }//end vertex Lb check
                delete vertexMother; vertexMother = nullptr;
              }//end check pion ID with Lc daughters
            }//end Lb pion filtering cuts
          }//end track-loop
          if(recVtx)fFiltCutsLbtoLcpi->CleanOwnPrimaryVtx(lctopkpi,aod,origownvtx);
          if(unsetvtx) lctopkpi->UnsetOwnPrimaryVtx();
        }//end Lc filtering cuts
      } else{
        fNentries->Fill(41); //monitor how often this fails
      }
    }//end Lc-loop
  }//end loop on cadidates
  
  delete vHF;
  return;
}




//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::ProcessInclusiveJet(AliAODEvent *aod, TClonesArray *arrMC){
  

  fTreeHandlerInclusiveJet->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, 0.0, NULL, 0.0, 0, NULL);
  fTreeHandlerInclusiveJet->SetAndFillInclusiveJetVars(aod->GetTracks(),arrMC);

  return;
}

void AliAnalysisTaskSEHFTreeCreator::ProcessMCGenInclusiveJet(TClonesArray *arrayMC){

  fTreeHandlerGenInclusiveJet->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong);
  fTreeHandlerGenInclusiveJet->SetAndFillInclusiveGenJetVars(arrayMC);

}

//_________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::ProcessMCGen(TClonesArray *arrayMC){
  /// Fill MC gen trees
  
  for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
    
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
    Int_t absPDG = TMath::Abs(mcPart->GetPdgCode());
    
    if(absPDG == 411 || absPDG == 421 || absPDG == 431 || absPDG == 4122 || absPDG == 521 || absPDG == 413 || absPDG == 531 || absPDG == 5122) {
      Bool_t isPrimary = kFALSE;
      Bool_t isFeeddown = kFALSE;
      //Bplus&Bs&Lb are "always" primary
      if(absPDG != 521 && absPDG != 531 && absPDG != 5122) {
        Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,!fITSUpgradeProduction);//Prompt = 4, FeedDown = 5
        if(orig!=4 && orig!=5) continue; //keep only prompt or feed-down
        
        if(orig==4){
          isPrimary = kTRUE;
          isFeeddown = kFALSE;
        }
        else if(orig==5){
          isPrimary = kFALSE;
          isFeeddown = kTRUE;
        }
      } else {
        isPrimary = kTRUE;
        isFeeddown = kFALSE;
      }
      
      Int_t  deca = 0;
      Int_t  deca2 = 0;
      Int_t  labDau[3] = {-1,-1,-1};
      Int_t  labDau2[3] = {-1,-1,-1}; //Needed for 2nd decay same particle
      Int_t  labDau4pr[4] = {-1,-1,-1,-1}; //Needed for 4 prong decays (like Bs and Lb)
      Bool_t isDaugInAcc = kFALSE;
      
      if(absPDG == 411 && fWriteVariableTreeDplus) {
        deca = AliVertexingHFUtils::CheckDplusDecay(arrayMC,mcPart,labDau);
        if(deca<1 || labDau[0]<0 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
        fTreeHandlerGenDplus->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenDplus->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenDplus->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenDplus->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenDplus->FillTree();
      }
      else if(absPDG == 421 && fWriteVariableTreeD0) {
        deca = AliVertexingHFUtils::CheckD0Decay(arrayMC,mcPart,labDau);
        if(deca!=1 || labDau[0]<0 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,2,labDau,fITSUpgradeProduction);
        fTreeHandlerGenD0->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenD0->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenD0->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenD0->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenD0->FillTree();
      }
      else if(absPDG == 431 && fWriteVariableTreeDs) {
        deca = AliVertexingHFUtils::CheckDsDecay(arrayMC,mcPart,labDau);
        if(deca!=fWriteVariableTreeDs || labDau[0]<0 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
        fTreeHandlerGenDs->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenDs->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenDs->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenDs->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenDs->FillTree();
      }
      else if(absPDG == 521 && fWriteVariableTreeBplus) {
        deca = AliVertexingHFUtils::CheckBplusDecay(arrayMC,mcPart,labDau);
        
        //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
        //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*deca - 1) is returned.
        //Use feed-down bit to flag these candidates for offline study.
        //Can be good for analysis or coming from an "incomplete simulated" decay
        if(deca == -2){ isFeeddown=kTRUE; deca=1; }

        if(deca!=1 || labDau[0]==-1 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
        fTreeHandlerGenBplus->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenBplus->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenBplus->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenBplus->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenBplus->FillTree();
      }
      else if(absPDG == 531 && fWriteVariableTreeBs) {
        deca = AliVertexingHFUtils::CheckBsDecay(arrayMC,mcPart,labDau4pr,kTRUE);
        //Only accept Bs-> pi Ds(->phipi->KKpi) decays

        //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
        //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*deca - 1) is returned.
        //Use feed-down bit to flag these candidates for offline study.
        //Can be good for analysis or coming from an "incomplete simulated" decay
        if(deca == -2){ isFeeddown=kTRUE; deca=1; }
        
        if(deca!=1 || labDau4pr[0]==-1 || labDau4pr[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,4,labDau4pr,fITSUpgradeProduction);
        fTreeHandlerGenBs->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenBs->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenBs->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenBs->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenBs->FillTree();
      }
      else if(absPDG == 413 && fWriteVariableTreeDstar) {
        deca = AliVertexingHFUtils::CheckDstarDecay(arrayMC,mcPart,labDau);
        if(deca!=1 || labDau[0]<0 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
        fTreeHandlerGenDstar->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenDstar->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenDstar->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenDstar->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenDstar->FillTree();
      }
      else if(absPDG == 4122 && (fWriteVariableTreeLctopKpi || fWriteVariableTreeLc2V0bachelor)) {
        deca = AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC,mcPart,labDau);
        deca2 = AliVertexingHFUtils::CheckLcV0bachelorDecay(arrayMC,mcPart,labDau2);
        if(deca<1 || labDau[0]==-1 || labDau[1]<0){
          if(deca2!=1 || labDau2[0]<0 || labDau2[1]<0 || !fWriteVariableTreeLc2V0bachelor) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau2,fITSUpgradeProduction);
          fTreeHandlerGenLc2V0bachelor->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenLc2V0bachelor->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenLc2V0bachelor->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
          if(fFillJets) fTreeHandlerGenLc2V0bachelor->SetGenJetVars(arrayMC,mcPart);
          fTreeHandlerGenLc2V0bachelor->FillTree();
        } else if(fWriteVariableTreeLctopKpi) {
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
          fTreeHandlerGenLctopKpi->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenLctopKpi->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenLctopKpi->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
          int resdecay = fTreeHandlerGenLctopKpi->GetLcResonantDecay(arrayMC,mcPart);
          fTreeHandlerGenLctopKpi->SetMCGenVariableResonantDecay(resdecay);
          if(fFillJets) fTreeHandlerGenLctopKpi->SetGenJetVars(arrayMC,mcPart);
          fTreeHandlerGenLctopKpi->FillTree();
        }
      } else if(absPDG == 5122 && fWriteVariableTreeLb) {
        deca = AliVertexingHFUtils::CheckLbDecay(arrayMC,mcPart,labDau4pr);
        
        //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
        //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*deca - 1) is returned.
        //Use feed-down bit to flag these candidates for offline study.
        //Can be good for analysis or coming from an "incomplete simulated" decay
        if(deca < -1){ isFeeddown=kTRUE; deca=-1*deca-1; }
        
        if(deca<1 || labDau4pr[0]==-1 || labDau4pr[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,4,labDau4pr,fITSUpgradeProduction);
        fTreeHandlerGenLb->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenLb->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenLb->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        if(fFillJets) fTreeHandlerGenLb->SetGenJetVars(arrayMC,mcPart);
        fTreeHandlerGenLb->FillTree();
      }
    }
  }  
}

//--------------------------------------------------------
Bool_t AliAnalysisTaskSEHFTreeCreator::CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau, Bool_t ITSUpgradeStudy){
  /// check if the decay products are in the good eta and pt range
    
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) {
      return kFALSE;
    }
    Double_t eta = mcPartDaughter->Eta();
    if(ITSUpgradeStudy){
      if (TMath::Abs(eta) >= 1) return kFALSE;
    } else {
      Double_t pt = mcPartDaughter->Pt();
      if (TMath::Abs(eta) > 0.9 || pt < 0.1) return kFALSE;
    }
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSEHFTreeCreator::IsCandidateFromHijing(AliAODRecoDecayHF *cand, AliAODMCHeader *mcHeader, TClonesArray* arrMC, AliAODTrack *tr){
  //
  // Returns true if candidate (or cand + track for on-the-fly reconstruction) are from Hijing, so not injected.
  //
  
  Bool_t candNotHijing = AliVertexingHFUtils::IsCandidateInjected(cand, mcHeader,arrMC);
  if(tr){
    Bool_t trackNotHijing = AliVertexingHFUtils::IsTrackInjected(tr,mcHeader,arrMC);
    if(!candNotHijing && !trackNotHijing) return kTRUE;
  } else {
    if(!candNotHijing) return kTRUE;
  }
  return kFALSE;
  
}
//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::SelectGoodTrackForReconstruction(AliAODEvent *aod, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags)
{
  //
  // Select good tracks for on the fly reconstruction (Bplus, Bs, Lb, ..) and return the array of their ids
  //
  
  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;
    
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(i));
    
    // skip tracks with negative ID
    // (these are duplicated TPC-only AOD tracks, for jet analysis...)
    if(track->GetID()<0) continue;
    
    // TEMPORARY: check that the cov matrix is there
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;

    // skip pure ITS SA tracks
    if(track->GetStatus()&AliESDtrack::kITSpureSA) continue;
    // skip tracks without ITS
    if(!(track->GetStatus()&AliESDtrack::kITSin)) continue;
    // avoid ghost TPC tracks
    if(!(track->GetStatus()&AliESDtrack::kTPCrefit)) continue;
    
    if(fGoodTrackFilterBit>=0 && !track->TestFilterMask(BIT(fGoodTrackFilterBit))) continue;
    if(fabs(track->Eta())>fGoodTrackEtaRange) continue;
    if(track->Pt()<fGoodTrackMinPt) continue;
    
    seleFlags[i]=kTRUE;
    nSeleTrks++;

  } // end loop on tracks
}
//________________________________________________________________
AliAODVertex* AliAnalysisTaskSEHFTreeCreator::ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion) {
  //
  // Helper function to recalculate a vertex.
  //
  
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  AliVertexerTracks vertexer;
  vertexer.SetFieldkG(bField);
  
  vertexer.SetVtxStart((AliESDVertex*)primary); //primary vertex
  vertexESD = (AliESDVertex*)vertexer.VertexForSelectedESDTracks(tracks);
  
  if (!vertexESD) return vertexAOD;
  
  if (vertexESD->GetNContributors() != tracks->GetEntriesFast()){
    delete vertexESD; vertexESD = nullptr;
    return vertexAOD;
  }
  
  // convert to AliAODVertex
  Double_t pos[3], cov[6], chi2perNDF;
  for (Int_t a = 0; a < 3; a++)pos[a] = 0.;
  for (Int_t b = 0; b < 6; b++)cov[b] = 0.;
  chi2perNDF = 0;
  
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  
  Double_t vertRadius2 = pos[0] * pos[0] + pos[1] * pos[1];
  if (vertRadius2 > 8.){ //(2.82)^2 radius beam pipe
    delete vertexESD; vertexESD = nullptr;
    return vertexAOD;
  }
  
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD = nullptr;
  Int_t nprongs = tracks->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, nprongs);
  
  return vertexAOD;
}

//________________________________________________________________
unsigned long AliAnalysisTaskSEHFTreeCreator::GetEvID() {
  TString currentfilename = ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()))->GetName();
  if(!fFileName.EqualTo(currentfilename)) {
    fEventNumber = 0;
    fFileName = currentfilename;
    TObjArray *path = fFileName.Tokenize("/");
    TString s = ((TObjString*)path->At( ((path->GetLast())-1) ))->GetString();
    fDirNumber = (unsigned int)s.Atoi();
    delete path;
  }
  Long64_t ev_number = Entry();
  if(fReadMC){
    ev_number = fEventNumber;
  }
  fEventNumber++;

  unsigned long evID = fPeriod & 0xfffffff;
  evID = (evID << 24) | (fOrbit & 0xffffff);
  evID = (evID << 12) | (fBC & 0xfff);
  return evID;
}

/**
 * @brief Create new container for MC particles and attach it to the task.
 * (copied from AliAnalysisTaskEmcal)
 *
 * The name provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new container for MC particles
 */
AliMCParticleContainer* AliAnalysisTaskSEHFTreeCreator::AddMCParticleContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;
  
  AliMCParticleContainer* cont = new AliMCParticleContainer(n);
  fParticleCollArray.Add(cont);
  
  return cont;
}

/**
 * @brief Create new track container and attach it to the task.
 * (copied from AliAnalysisTaskEmcal)
 *
 * The name provided to this function must match the name of the array
 * attached to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new track container
 */
AliTrackContainer* AliAnalysisTaskSEHFTreeCreator::AddTrackContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;
  
  AliTrackContainer* cont = new AliTrackContainer(n);
  fParticleCollArray.Add(cont);
  
  return cont;
}

/**
 * @brief Create new particle container and attach it to the task.
 * (copied from AliAnalysisTaskEmcal)
 *
 * The name provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new particle container
 */
AliParticleContainer* AliAnalysisTaskSEHFTreeCreator::AddParticleContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;
  
  AliParticleContainer* cont = new AliParticleContainer(n);
  fParticleCollArray.Add(cont);
  
  return cont;
}

/**
 * @brief Get \f$ i^{th} \f$ particle container attached to this task
 * (copied from AliAnalysisTaskEmcal)
 * @param[in] i Index of the particle container
 * @return Particle container found for the given index (NULL if no particle container exists for that index)
 */
AliParticleContainer* AliAnalysisTaskSEHFTreeCreator::GetParticleContainer(Int_t i) const
{
  if (i<0 || i>fParticleCollArray.GetEntriesFast()) return 0;
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  return cont;
}

/**
 * @brief Find particle container attached to this task according to its name
 * (copied from AliAnalysisTaskEmcal)
 * @param[in] name Name of the particle container
 * @return Particle container found under the given name
 */
AliParticleContainer* AliAnalysisTaskSEHFTreeCreator::GetParticleContainer(const char *name) const
{
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.FindObject(name));
  return cont;
}


/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * (copied from AliAnalysisTaskEmcalJet)
 * @param[in] jetType One of the AliJetContainer::EJetType_t enumeration values (charged, full, neutral)
 * @param[in] jetAlgo One of the AliJetContainer::EJetAlgo_t enumeration values (anti-kt, kt, ...)
 * @param[in] recoScheme One of the AliJetContainer::ERecoScheme_t enumeration values (pt-scheme, ...)
 * @param[in] radius Resolution parameter (0.2, 0.4, ...)
 * @param[in] accType One of the AliEmcalJet::JetAcceptanceType enumeration values (kTPC, kEMCAL, kDCAL, ...),
 * or a combination using bitwise OR: For example, (kEMCAL | kDCAL) will select all jets in either EMCal or DCal.
 * @param[in] partCont Particle container of the objects used to generate the jets
 * @param[in] clusCont Cluster container of the objects used to generate the jets
 * @param[in] tag Label to distinguish different jet branches (defaul is 'Jet')
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskSEHFTreeCreator::AddJetContainer(AliJetContainer::EJetType_t jetType, AliJetContainer::EJetAlgo_t jetAlgo, AliJetContainer::ERecoScheme_t recoScheme, Double_t radius, UInt_t accType, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag)
{
  AliJetContainer *cont = new AliJetContainer(jetType, jetAlgo, recoScheme, radius, partCont, clusCont, tag);
  cont->SetJetAcceptanceType(accType);
  fJetCollArray.Add(cont);
  
  return cont;
}

/**
 * Create new jet container and attach it to the task. This method is usually called in the add task macro.
 * (copied from AliAnalysisTaskEmcalJet)
 * @param[in] n Name of the jet branch
 * @param[in] accType One of the AliEmcalJet::JetAcceptanceType enumeration values (kTPC, kEMCAL, kDCAL, ...),
 * or a combination using bitwise OR: For example, (kEMCAL | kDCAL) will select all jets in either EMCal or DCal.
 * @param[in] jetRadius Resolution parameter (0.2, 0.4, ...)
 * @return Pointer to the new jet container
 */
AliJetContainer* AliAnalysisTaskSEHFTreeCreator::AddJetContainer(const char *n, UInt_t accType, Float_t jetRadius)
{
  if (TString(n).IsNull()) return 0;
  
  AliJetContainer *cont = new AliJetContainer(n);
  cont->SetJetRadius(jetRadius);
  cont->SetJetAcceptanceType(accType);
  fJetCollArray.Add(cont);
  
  return cont;
}

/**
 * Get \f$ i^{th} \f$ jet container attached to this task
 * (copied from AliAnalysisTaskEmcalJet)
 * @param[in] i Index of the jet container
 * @return Jet container found for the given index (NULL if no jet container exists for that index)
 */
AliJetContainer* AliAnalysisTaskSEHFTreeCreator::GetJetContainer(Int_t i) const
{
  if (i < 0 || i >= fJetCollArray.GetEntriesFast()) return 0;
  AliJetContainer *cont = static_cast<AliJetContainer*>(fJetCollArray.At(i));
  return cont;
}
