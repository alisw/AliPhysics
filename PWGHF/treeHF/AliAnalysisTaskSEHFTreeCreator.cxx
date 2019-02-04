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
////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <THnSparse.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TChain.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
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
#include "AliHFTreeHandlerDstartoKpipi.h"
#include "AliAnalysisTaskSEHFTreeCreator.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEHFTreeCreator);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreator::AliAnalysisTaskSEHFTreeCreator():
AliAnalysisTaskSE(),
fNentries(0x0),
fHistoNormCounter(0x0),
fListCuts(0x0),
fFiltCutsD0toKpi(0x0),
fFiltCutsDstoKKpi(0x0),
fFiltCutsDplustoKpipi(0x0),
fFiltCutsLctopKpi(0x0),
fFiltCutsBplustoD0pi(0x0),
fFiltCutsDstartoKpipi(0x0),
fCutsD0toKpi(0x0),
fCutsDstoKKpi(0x0),
fCutsDplustoKpipi(0x0),
fCutsLctopKpi(0x0),
fCutsBplustoD0pi(0x0),
fCutsDstartoKpipi(0x0),
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
fWriteVariableTreeDstar(0),
fVariablesTreeD0(0x0),
fVariablesTreeDs(0x0),
fVariablesTreeDplus(0x0),
fVariablesTreeLctopKpi(0x0),
fVariablesTreeBplus(0x0),
fVariablesTreeDstar(0x0),
fGenTreeD0(0x0),
fGenTreeDs(0x0),
fGenTreeDplus(0x0),
fGenTreeLctopKpi(0x0),
fGenTreeBplus(0x0),
fGenTreeDstar(0x0),
fTreeEvChar(0x0),
fWriteOnlySignal(kFALSE),
fTreeHandlerD0(0x0),
fTreeHandlerDs(0x0),
fTreeHandlerDplus(0x0),
fTreeHandlerLctopKpi(0x0),
fTreeHandlerBplus(0x0),
fTreeHandlerDstar(0x0),
fTreeHandlerGenD0(0x0),
fTreeHandlerGenDs(0x0),
fTreeHandlerGenDplus(0x0),
fTreeHandlerGenLctopKpi(0x0),
fTreeHandlerGenBplus(0x0),
fTreeHandlerGenDstar(0x0),
fPIDresp(0x0),
fPIDoptD0(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDs(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDplus(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptLctopKpi(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptBplus(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDstar(AliHFTreeHandler::kRawAndNsigmaPID),
fCentrality(-999.),
fzVtxReco(0.),
fzVtxGen(0.),
fNcontributors(0),
fNtracks(0),
fIsEvRej(0),
fRunNumber(0),
fFillMCGenTrees(kTRUE),
fDsMassKKOpt(1),
fTreeSingleTrackVarsOpt(AliHFTreeHandler::kRedSingleTrackVars)
{

/// Default constructor

}
//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreator::AliAnalysisTaskSEHFTreeCreator(const char *name, TList *cutsList):
AliAnalysisTaskSE(name),
fNentries(0x0),
fHistoNormCounter(0x0),
fListCuts(0x0),
fFiltCutsD0toKpi(0x0),
fFiltCutsDstoKKpi(0x0),
fFiltCutsDplustoKpipi(0x0),
fFiltCutsLctopKpi(0x0),
fFiltCutsBplustoD0pi(0x0),
fFiltCutsDstartoKpipi(0x0),
fCutsD0toKpi(0x0),
fCutsDstoKKpi(0x0),
fCutsDplustoKpipi(0x0),
fCutsLctopKpi(0x0),
fCutsBplustoD0pi(0x0),
fCutsDstartoKpipi(0x0),
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
fWriteVariableTreeDstar(0),
fVariablesTreeD0(0x0),
fVariablesTreeDs(0x0),
fVariablesTreeDplus(0x0),
fVariablesTreeLctopKpi(0x0),
fVariablesTreeBplus(0x0),
fVariablesTreeDstar(0x0),
fGenTreeD0(0x0),
fGenTreeDs(0x0),
fGenTreeDplus(0x0),
fGenTreeLctopKpi(0x0),
fGenTreeBplus(0x0),
fGenTreeDstar(0x0),
fTreeEvChar(0x0),
fWriteOnlySignal(kFALSE),
fTreeHandlerD0(0x0),
fTreeHandlerDs(0x0),
fTreeHandlerDplus(0x0),
fTreeHandlerLctopKpi(0x0),
fTreeHandlerBplus(0x0),
fTreeHandlerDstar(0x0),
fTreeHandlerGenD0(0x0),
fTreeHandlerGenDs(0x0),
fTreeHandlerGenDplus(0x0),
fTreeHandlerGenLctopKpi(0x0),
fTreeHandlerGenBplus(0x0),
fTreeHandlerGenDstar(0x0),
fPIDresp(0x0),
fPIDoptD0(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDs(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDplus(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptLctopKpi(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptBplus(AliHFTreeHandler::kRawAndNsigmaPID),
fPIDoptDstar(AliHFTreeHandler::kRawAndNsigmaPID),
fCentrality(-999.),
fzVtxReco(0.),
fzVtxGen(0.),
fNcontributors(0),
fNtracks(0),
fIsEvRej(0),
fRunNumber(0),
fFillMCGenTrees(kTRUE),
fDsMassKKOpt(1),
fTreeSingleTrackVarsOpt(AliHFTreeHandler::kRedSingleTrackVars)
{
    /// Standard constructor
    
    
    if(fFiltCutsD0toKpi){
    delete fFiltCutsD0toKpi;fFiltCutsD0toKpi=NULL;
    }
    if(fFiltCutsDstoKKpi){
    delete fFiltCutsDstoKKpi;fFiltCutsDstoKKpi=NULL;
    }
    if(fFiltCutsDplustoKpipi){
    delete fFiltCutsDplustoKpipi;fFiltCutsDplustoKpipi=NULL;
    }
    if(fFiltCutsLctopKpi){
    delete fFiltCutsLctopKpi;fFiltCutsLctopKpi=NULL;
    }
    if(fFiltCutsBplustoD0pi){
    delete fFiltCutsBplustoD0pi;fFiltCutsBplustoD0pi=NULL;
    }
    if(fFiltCutsDstartoKpipi){
      delete fFiltCutsDstartoKpipi;fFiltCutsDstartoKpipi=NULL;
    }
    if(fCutsD0toKpi){
    delete fCutsD0toKpi;fCutsD0toKpi=NULL;
    }
    if(fCutsDstoKKpi){
    delete fCutsDstoKKpi;fCutsDstoKKpi=NULL;
    }
    if(fCutsDplustoKpipi){
    delete fCutsDplustoKpipi;fCutsDplustoKpipi=NULL;
    }
    if(fCutsLctopKpi){
    delete fCutsLctopKpi;fCutsLctopKpi=NULL;
    }
    if(fCutsBplustoD0pi){
    delete fCutsBplustoD0pi;fCutsBplustoD0pi=NULL;
    }
    if(fCutsDstartoKpipi){
      delete fCutsDstartoKpipi;fCutsDstartoKpipi=NULL;
    }
    fListCuts=cutsList;
    
    fFiltCutsD0toKpi     =(AliRDHFCutsD0toKpi*)fListCuts->FindObject("D0toKpiFilteringCuts");
    fFiltCutsDstoKKpi    =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("DstoKKpiFilteringCuts");
    fFiltCutsDplustoKpipi=(AliRDHFCutsDplustoKpipi*)fListCuts->FindObject("DplustoKpipiFilteringCuts");
    fFiltCutsLctopKpi    =(AliRDHFCutsLctopKpi*)fListCuts->FindObject("LctopKpiFilteringCuts");
    fFiltCutsBplustoD0pi =(AliRDHFCutsBPlustoD0Pi*)fListCuts->FindObject("BplustoD0piFilteringCuts");
    fFiltCutsDstartoKpipi=(AliRDHFCutsDStartoKpipi*)fListCuts->FindObject("DstartoKpipiFilteringCuts");
    fCutsD0toKpi         =(AliRDHFCutsD0toKpi*)fListCuts->FindObject("D0toKpiAnalysisCuts");
    fCutsDstoKKpi        =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("DstoKKpiAnalysisCuts");
    fCutsDplustoKpipi    =(AliRDHFCutsDplustoKpipi*)fListCuts->FindObject("DplustoKpipiAnalysisCuts");
    fCutsLctopKpi        =(AliRDHFCutsLctopKpi*)fListCuts->FindObject("LctopKpiAnalysisCuts");
    fCutsBplustoD0pi     =(AliRDHFCutsBPlustoD0Pi*)fListCuts->FindObject("BplustoD0piAnalysisCuts");
    fCutsDstartoKpipi    =(AliRDHFCutsDStartoKpipi*)fListCuts->FindObject("DstartoKpipiAnalysisCuts");

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
}

//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreator::~AliAnalysisTaskSEHFTreeCreator()
{
    if (fListCuts) {
        delete fListCuts;
        fListCuts = 0x0;
    }
    if (fFiltCutsD0toKpi) {
        delete fFiltCutsD0toKpi;
        fFiltCutsD0toKpi = 0x0;
    }
    if (fFiltCutsDstoKKpi) {
        delete fFiltCutsDstoKKpi;
        fFiltCutsDstoKKpi = 0x0;
    }
    if (fFiltCutsDplustoKpipi) {
        delete fFiltCutsDplustoKpipi;
        fFiltCutsDplustoKpipi = 0x0;
    }
    if (fFiltCutsLctopKpi) {
        delete fFiltCutsLctopKpi;
        fFiltCutsLctopKpi = 0x0;
    }
    if (fFiltCutsBplustoD0pi) {
        delete fFiltCutsBplustoD0pi;
        fFiltCutsBplustoD0pi = 0x0;
    }
    if (fFiltCutsDstartoKpipi) {
      delete fFiltCutsDstartoKpipi;
      fFiltCutsDstartoKpipi = 0x0;
    }
    if (fCutsD0toKpi) {
        delete fCutsD0toKpi;
        fCutsD0toKpi = 0x0;
    }
    if (fCutsDstoKKpi) {
        delete fCutsDstoKKpi;
        fCutsDstoKKpi = 0x0;
    }
    if (fCutsDplustoKpipi) {
        delete fCutsDplustoKpipi;
        fCutsDplustoKpipi = 0x0;
    }
    if (fCutsLctopKpi) {
        delete fCutsLctopKpi;
        fCutsLctopKpi = 0x0;
    }
    if (fCutsBplustoD0pi) {
        delete fCutsBplustoD0pi;
        fCutsBplustoD0pi = 0x0;
    }
    if (fCutsDstartoKpipi) {
      delete fCutsDstartoKpipi;
      fCutsDstartoKpipi = 0x0;
    }
    if (fNentries){
        delete fNentries;
        fNentries = 0x0;
    }
    if (fHistoNormCounter){
        delete fHistoNormCounter;
        fHistoNormCounter = 0x0;
    }
    if (fListCounter) {
        delete fListCounter;
        fListCounter = 0x0;
    }
    if(fCounter){
        delete fCounter;
        fCounter=0x0;
    }
    if(fTreeHandlerD0) {
      delete fTreeHandlerD0;
      fTreeHandlerD0 = 0x0;
    }
    if(fTreeHandlerDs) {
      delete fTreeHandlerDs;
      fTreeHandlerDs = 0x0;
    }
    if(fTreeHandlerDplus) {
      delete fTreeHandlerDplus;
      fTreeHandlerDplus = 0x0;
    }
    if(fTreeHandlerLctopKpi) {
      delete fTreeHandlerLctopKpi;
      fTreeHandlerLctopKpi = 0x0;
    }
    if(fTreeHandlerBplus) {
        delete fTreeHandlerBplus;
        fTreeHandlerBplus = 0x0;
    }
    if(fTreeHandlerDstar) {
      delete fTreeHandlerDstar;
      fTreeHandlerDstar = 0x0;
    }
    if(fTreeHandlerGenD0) {
      delete fTreeHandlerGenD0;
      fTreeHandlerGenD0 = 0x0;
    }
    if(fTreeHandlerGenDs) {
      delete fTreeHandlerGenDs;
      fTreeHandlerGenDs = 0x0;
    }
    if(fTreeHandlerGenDplus) {
      delete fTreeHandlerGenDplus;
      fTreeHandlerGenDplus = 0x0;
    }
    if(fTreeHandlerGenLctopKpi) {
      delete fTreeHandlerGenLctopKpi;
      fTreeHandlerGenLctopKpi = 0x0;
    }
    if(fTreeHandlerGenBplus) {
        delete fTreeHandlerGenBplus;
        fTreeHandlerGenBplus = 0x0;
    }
    if(fTreeHandlerGenDstar) {
        delete fTreeHandlerGenDstar;
        fTreeHandlerGenDstar = 0x0;
    }
    if(fTreeEvChar) {
        delete fTreeEvChar;
        fTreeEvChar = 0x0;
    }
    
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
    fNentries=new TH1F(nameoutput, "Number of events", 34,-0.5,33.5);
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
    fNentries->GetXaxis()->SetBinLabel(26,"n. Bplus after filtering");
    fNentries->GetXaxis()->SetBinLabel(27,"n. Bplus after selection");
    fNentries->GetXaxis()->SetBinLabel(28,"n. of not on-the-fly rec Bplus");
    fNentries->GetXaxis()->SetBinLabel(29,"n. of cascade candidates");
    fNentries->GetXaxis()->SetBinLabel(30,"n. Dstar after filtering");
    fNentries->GetXaxis()->SetBinLabel(31,"n. Dstar after selection");
    fNentries->GetXaxis()->SetBinLabel(32,"n. of not on-the-fly rec Dstar");
    
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
    fCounter->Init();
    fListCounter->Add(fCounter);
    
    //count number of enabled trees
    Int_t nEnabledTrees = 1; // event tree always enabled
    if(fWriteVariableTreeD0) nEnabledTrees++;
    if(fWriteVariableTreeDs) nEnabledTrees++;
    if(fWriteVariableTreeDplus) nEnabledTrees++;
    if(fWriteVariableTreeLctopKpi) nEnabledTrees++;
    if(fWriteVariableTreeBplus) nEnabledTrees++;
    if(fWriteVariableTreeDstar) nEnabledTrees++;
    if(fReadMC && fFillMCGenTrees) {
        nEnabledTrees = (nEnabledTrees-1)*2+1;
    }

    //
    // Output slot 4-17 : trees of the candidate and event-characterization variables
    //
    OpenFile(5);
    fTreeEvChar = new TTree("tree_event_char","tree_event_char");
    //set variables
    TString varnames[7] = {"centrality", "z_vtx_reco", "n_vtx_contributors", "n_tracks", "is_ev_rej", "run_number", "z_vtx_gen"};
    fTreeEvChar->Branch(varnames[0].Data(),&fCentrality,Form("%s/F",varnames[0].Data()));
    fTreeEvChar->Branch(varnames[1].Data(),&fzVtxReco,Form("%s/F",varnames[1].Data()));
    fTreeEvChar->Branch(varnames[2].Data(),&fNcontributors,Form("%s/I",varnames[2].Data()));
    fTreeEvChar->Branch(varnames[3].Data(),&fNtracks,Form("%s/I",varnames[3].Data()));
    fTreeEvChar->Branch(varnames[4].Data(),&fIsEvRej,Form("%s/I",varnames[4].Data()));
    fTreeEvChar->Branch(varnames[5].Data(),&fRunNumber,Form("%s/I",varnames[5].Data()));
    if(fReadMC) fTreeEvChar->Branch(varnames[6].Data(),&fzVtxGen,Form("%s/F",varnames[6].Data()));
    //fTreeEvChar->SetMaxVirtualSize(1.e+8/nEnabledTrees);

    if(fWriteVariableTreeD0){
        OpenFile(6);
        TString nameoutput = "tree_D0";
        fTreeHandlerD0 = new AliHFTreeHandlerD0toKpi(fPIDoptD0);
        fTreeHandlerD0->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
        if(fReadMC && fWriteOnlySignal) fTreeHandlerD0->SetFillOnlySignal(fWriteOnlySignal);
        fVariablesTreeD0 = (TTree*)fTreeHandlerD0->BuildTree(nameoutput,nameoutput);
        //fVariablesTreeD0->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeD0);
      
        if(fFillMCGenTrees && fReadMC) {
          OpenFile(7);
          TString nameoutput = "tree_D0_gen";
          fTreeHandlerGenD0 = new AliHFTreeHandlerD0toKpi(0);
          fGenTreeD0 = (TTree*)fTreeHandlerGenD0->BuildTreeMCGen(nameoutput,nameoutput);
          //fGenTreeD0->SetMaxVirtualSize(1.e+8/nEnabledTrees);
          fTreeEvChar->AddFriend(fGenTreeD0);
        }
    }
    if(fWriteVariableTreeDs){
        OpenFile(8);
        TString nameoutput = "tree_Ds";
        fTreeHandlerDs = new AliHFTreeHandlerDstoKKpi(fPIDoptDs);
        fTreeHandlerDs->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
        if(fReadMC && fWriteOnlySignal) fTreeHandlerDs->SetFillOnlySignal(fWriteOnlySignal);
        fTreeHandlerDs->SetMassKKOption(fDsMassKKOpt);
        fVariablesTreeDs = (TTree*)fTreeHandlerDs->BuildTree(nameoutput,nameoutput);
        //fVariablesTreeDs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeDs);
      
      if(fFillMCGenTrees && fReadMC) {
        OpenFile(9);
        TString nameoutput = "tree_Ds_gen";
        fTreeHandlerGenDs = new AliHFTreeHandlerDstoKKpi(0);
        fGenTreeDs = (TTree*)fTreeHandlerGenDs->BuildTreeMCGen(nameoutput,nameoutput);
        //fGenTreeDs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fGenTreeDs);
      }
    }
    if(fWriteVariableTreeDplus){
        OpenFile(10);
        TString nameoutput = "tree_Dplus";
        fTreeHandlerDplus = new AliHFTreeHandlerDplustoKpipi(fPIDoptDplus);
        fTreeHandlerDplus->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
        if(fReadMC && fWriteOnlySignal) fTreeHandlerDplus->SetFillOnlySignal(fWriteOnlySignal);
        fVariablesTreeDplus = (TTree*)fTreeHandlerDplus->BuildTree(nameoutput,nameoutput);
        //fVariablesTreeDplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeDplus);
      if(fFillMCGenTrees && fReadMC) {
        OpenFile(11);
        TString nameoutput = "tree_Dplus_gen";
        fTreeHandlerGenDplus = new AliHFTreeHandlerDplustoKpipi(0);
        fGenTreeDplus = (TTree*)fTreeHandlerGenDplus->BuildTreeMCGen(nameoutput,nameoutput);
        //fGenTreeDplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fGenTreeDplus);
      }
    }
    if(fWriteVariableTreeLctopKpi){
        OpenFile(12);
        TString nameoutput = "tree_LctopKpi";
        fTreeHandlerLctopKpi = new AliHFTreeHandlerLctopKpi(fPIDoptLctopKpi);
        fTreeHandlerLctopKpi->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
        if(fReadMC && fWriteOnlySignal) fTreeHandlerLctopKpi->SetFillOnlySignal(fWriteOnlySignal);
        fVariablesTreeLctopKpi = (TTree*)fTreeHandlerLctopKpi->BuildTree(nameoutput,nameoutput);
        //fVariablesTreeLctopKpi->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeLctopKpi);
      if(fFillMCGenTrees && fReadMC) {
        OpenFile(13);
        TString nameoutput = "tree_LctopKpi_gen";
        fTreeHandlerGenLctopKpi = new AliHFTreeHandlerLctopKpi(0);
        fGenTreeLctopKpi = (TTree*)fTreeHandlerGenLctopKpi->BuildTreeMCGen(nameoutput,nameoutput);
        //fGenTreeLctopKpi->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fGenTreeLctopKpi);
      }
    }
    if(fWriteVariableTreeBplus){
        OpenFile(14);
        TString nameoutput = "tree_Bplus";
        fTreeHandlerBplus = new AliHFTreeHandlerBplustoD0pi(fPIDoptBplus);
        fTreeHandlerBplus->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
        if(fReadMC && fWriteOnlySignal) fTreeHandlerBplus->SetFillOnlySignal(fWriteOnlySignal);
        fVariablesTreeBplus = (TTree*)fTreeHandlerBplus->BuildTree(nameoutput,nameoutput);
        //fVariablesTreeBplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeBplus);
        if(fFillMCGenTrees && fReadMC) {
            OpenFile(15);
            TString nameoutput = "tree_Bplus_gen";
            fTreeHandlerGenBplus = new AliHFTreeHandlerBplustoD0pi(0);
            fGenTreeBplus = (TTree*)fTreeHandlerGenBplus->BuildTreeMCGen(nameoutput,nameoutput);
            //fGenTreeBplus->SetMaxVirtualSize(1.e+8/nEnabledTrees);
            fTreeEvChar->AddFriend(fGenTreeBplus);
        }
    }
    if(fWriteVariableTreeDstar){
        OpenFile(16);
        TString nameoutput = "tree_Dstar";
        fTreeHandlerDstar = new AliHFTreeHandlerDstartoKpipi(fPIDoptDstar);
        fTreeHandlerDstar->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
        if(fReadMC && fWriteOnlySignal) fTreeHandlerDstar->SetFillOnlySignal(fWriteOnlySignal);
        fVariablesTreeDstar = (TTree*)fTreeHandlerDstar->BuildTree(nameoutput,nameoutput);
        //fVariablesTreeDstar->SetMaxVirtualSize(1.e+8/nEnabledTrees);
        fTreeEvChar->AddFriend(fVariablesTreeDstar);
        if(fFillMCGenTrees && fReadMC) {
            OpenFile(17);
            TString nameoutput = "tree_Dstar_gen";
            fTreeHandlerGenDstar = new AliHFTreeHandlerDstartoKpipi(0);
            fGenTreeDstar = (TTree*)fTreeHandlerGenDstar->BuildTreeMCGen(nameoutput,nameoutput);
            //fGenTreeDstar->SetMaxVirtualSize(1.e+8/nEnabledTrees);
            fTreeEvChar->AddFriend(fGenTreeDstar);
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
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::UserExec(Option_t */*option*/)

{
    /// Execute analysis for current event:
    
    AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
    
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
    TString candidatesCascArrayName="Dstar";
    TClonesArray *array2prong=0;
    TClonesArray *array3Prong=0;
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
            arrayCasc=(TClonesArray*)aodFromExt->GetList()->FindObject(candidatesCascArrayName.Data());
        }
    } else if(aod) {
        array2prong=(TClonesArray*)aod->GetList()->FindObject(candidates2prongArrayName.Data());
        array3Prong=(TClonesArray*)aod->GetList()->FindObject(candidates3prongArrayName.Data());
        arrayCasc=(TClonesArray*)aod->GetList()->FindObject(candidatesCascArrayName.Data());
    }
    
    if(!array2prong || !array3Prong || !aod || !arrayCasc) {
        printf("AliAnalysisTaskSEHFTreeCreator::UserExec: input branches not found!\n");
        return;
    }
    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;
    
    fNentries->Fill(3); // count events
    
    AliRDHFCuts* evCuts = 0x0;
    if(fWriteVariableTreeD0) evCuts = (AliRDHFCuts*)fFiltCutsD0toKpi->Clone();
    else if(fWriteVariableTreeDplus) evCuts = (AliRDHFCuts*)fFiltCutsDplustoKpipi->Clone();
    else if(fWriteVariableTreeDs) evCuts = (AliRDHFCuts*)fFiltCutsDstoKKpi->Clone();
    else if(fWriteVariableTreeDstar) evCuts = (AliRDHFCuts*)fFiltCutsDstartoKpipi->Clone();
    else if(fWriteVariableTreeBplus) evCuts = (AliRDHFCuts*)fFiltCutsBplustoD0pi->Clone();
    else if(fWriteVariableTreeLctopKpi) evCuts = (AliRDHFCuts*)fFiltCutsLctopKpi->Clone();

    if(!evCuts){
      //Setting one of the filtering cuts as event selection cuts if evCuts is still empty
      //(which is the case if all meson TTrees are disabled).
      //In this way, the task doesn't crash and fTreeEvChar is still filled.
      if(fFiltCutsD0toKpi) evCuts = (AliRDHFCuts*)fFiltCutsD0toKpi->Clone();
      else if(fFiltCutsDplustoKpipi) evCuts = (AliRDHFCuts*)fFiltCutsDplustoKpipi->Clone();
      else if(fFiltCutsDstartoKpipi) evCuts = (AliRDHFCuts*)fFiltCutsDstartoKpipi->Clone();
      else if(fFiltCutsBplustoD0pi) evCuts = (AliRDHFCuts*)fFiltCutsBplustoD0pi->Clone();
      else if(fFiltCutsLctopKpi) evCuts = (AliRDHFCuts*)fFiltCutsLctopKpi->Clone();
      else{
        printf("AliAnalysisTaskSEHFTreeCreator::UserExec: No event selection cuts found!\n");
        return;
      }
    }

    if(evCuts->IsEventRejectedDueToTrigger())fNentries->Fill(5);
    if(evCuts->IsEventRejectedDueToNotRecoVertex())fNentries->Fill(6);
    if(evCuts->IsEventRejectedDueToVertexContributors())fNentries->Fill(7);
    if(evCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fNentries->Fill(8);
    if(evCuts->IsEventRejectedDueToPileup())fNentries->Fill(9);
    if(evCuts->IsEventRejectedDueToCentrality())fNentries->Fill(10);
    
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
    Bool_t isSameEvSelDstar=kTRUE;

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
    if(fWriteVariableTreeDstar)
      isSameEvSelDstar=!((fFiltCutsDstartoKpipi->IsEventSelected(aod) && !fCutsDstartoKpipi->IsEventSelected(aod))||(!fFiltCutsDstartoKpipi->IsEventSelected(aod) && fCutsDstartoKpipi->IsEventSelected(aod)));

    Bool_t isSameEvSel = isSameEvSelD0 && isSameEvSelDs && isSameEvSelDplus && isSameEvSelLctopKpi && isSameEvSelBplus && isSameEvSelDstar;
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
       (fWriteVariableTreeBplus && fWriteVariableTreeDstar && (fFiltCutsBplustoD0pi->IsEventSelected(aod)!=fFiltCutsDstartoKpipi->IsEventSelected(aod)))
        ){
      Printf("AliAnalysisTaskSEHFTreeCreator::UserExec: differences in the event selection cuts different meson");
      return;
    }
    
    fCounter->StoreEvent(aod,evCuts,fReadMC);
    Bool_t isEvSel=evCuts->IsEventSelected(aod);
    fCentrality = evCuts->GetCentrality(aod);
    if(fCentrality<0) fCentrality=-1.;
    //normalisation counter
    if(!evCuts->IsEventRejectedDuePhysicsSelection()){
      if(isEvSel){
      //selected events with primary vertex
      fHistoNormCounter->Fill(0.,fCentrality);
      }
      else{
      	if(evCuts->GetWhyRejection()==0){
      	//rejected events bc no primary vertex
      	fHistoNormCounter->Fill(1.,fCentrality);
      	}
      	//rejected events bc good primary vertex but >10cm
        if(evCuts->GetWhyRejection()==6){
        //nPrimaryV++;
        fHistoNormCounter->Fill(0.,fCentrality);
        //nzVtxGT10++;
        fHistoNormCounter->Fill(2.,fCentrality);
        }
        if(evCuts->GetWhyRejection()==1){
        //nPileup++;
        fHistoNormCounter->Fill(4.,fCentrality);
        }
      }
      if(evCuts->CountEventForNormalization()){
        //nCountForNorm++;
        fHistoNormCounter->Fill(3.,fCentrality);
      }
    }
    
    Bool_t isEvRejCent  = evCuts->IsEventRejectedDueToCentrality();

    if(!isEvSel && isEvRejCent) return; //cut only centrality, else tag only
    if(isEvSel) fNentries->Fill(4);
    // AOD primary vertex
    AliAODVertex *vtx = (AliAODVertex*)aod->GetPrimaryVertex();
    fNcontributors = vtx->GetNContributors();
    fzVtxReco = vtx->GetZ();
    fNtracks = aod->GetNumberOfTracks();
    fIsEvRej = evCuts->GetEventRejectionBitMap();
    fRunNumber=aod->GetRunNumber();
    fTreeEvChar->Fill();
    
    //get PID response
    if(!fPIDresp) fPIDresp = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();

    if(fWriteVariableTreeD0 || fWriteVariableTreeBplus) Process2Prong(array2prong,aod,mcArray,aod->GetMagneticField());
    if(fWriteVariableTreeDs || fWriteVariableTreeDplus || fWriteVariableTreeLctopKpi) Process3Prong(array3Prong,aod,mcArray,aod->GetMagneticField());
    if(fWriteVariableTreeDstar) ProcessCasc(arrayCasc,aod,mcArray,aod->GetMagneticField());
    if(fFillMCGenTrees && fReadMC) ProcessMCGen(mcArray);
  
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

    delete evCuts;
    evCuts = 0x0;

    return;
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
void AliAnalysisTaskSEHFTreeCreator::Process2Prong(TClonesArray *array2prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield){
    
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    //Needed separate one for Bplus, because D0 deletes its daughters from PV
    AliAODVertex *vtx1_Bplus = (AliAODVertex*)aod->GetPrimaryVertex();

    Int_t n2prong = array2prong->GetEntriesFast();
    if(fDebug>2) printf("Number of D0->Kpi: %d\n",n2prong);

    Int_t pdgDgD0toKpi[2]={321,211};
    Int_t pdgDgD0topiK[2]={211,321};
    Int_t nSelectedD0=0;
    Int_t nFilteredD0=0;

    Int_t pdgDgBplustoD0piInt[2] = {211,421};
    Int_t nSelectedBplus = 0;
    Int_t nFilteredBplus = 0;

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

                    if(fReadMC) {
                        AliAODMCParticle *partD0=0x0;
                        labD0 = d->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
                        if(labD0>=0){
                            partD0 = (AliAODMCParticle*)arrMC->At(labD0);
                            pdgD0 = partD0->GetPdgCode();
                            origin = AliVertexingHFUtils::CheckOrigin(arrMC,partD0,kTRUE);
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
                                if(origin==4) isprompt=kTRUE;
                                else if(origin==5) isFD=kTRUE;
                                if(pdgD0==421){
                                    issignal=kTRUE;
                                }
                                else {
                                    isrefl=kTRUE;
                                }
                            }//end labD0check
                            else{//background
                                isbkg=kTRUE;
                            }
                            fTreeHandlerD0->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                        }//end read MC
                        fTreeHandlerD0->SetIsSelectedStd(isSelAnCutsD0, isSelTopoAnCutsD0, isSelPidAnCutsD0);
                        fTreeHandlerD0->SetVariables(d,bfield,masshypo,fPIDresp);
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
                                if(origin==4) isprompt=kTRUE;
                                else if(origin==5) isFD=kTRUE;
                                if(pdgD0==-421){
                                    issignal=kTRUE;
                                }
                                else {
                                    isrefl=kTRUE;
                                }
                            } //end label check
                            else{ //background MC
                                isbkg=kTRUE;
                            }
                            fTreeHandlerD0->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                        }//end readMC
                        fTreeHandlerD0->SetIsSelectedStd(isSelAnCutsD0bar, isSelTopoAnCutsD0bar, isSelPidAnCutsD0bar);
                        fTreeHandlerD0->SetVariables(d,bfield,masshypo,fPIDresp);
                    }//end D0bar
                    if(recVtx)fFiltCutsD0toKpi->CleanOwnPrimaryVtx(d,aod,origownvtx);
                    if(unsetvtx) d->UnsetOwnPrimaryVtx();
                }//end is selected filt
            }
            else {
                fNentries->Fill(14); //monitor how often this fails
            }
        }//end D0
        
        //*************************************
        //Bplus (need to be reconstructed on the fly using a D0 candidate)
        
        //D0 from Bplus
        Bool_t isD0fromBplustagged=kTRUE;
        AliAODRecoDecayHF2Prong *dfromB = (AliAODRecoDecayHF2Prong*)array2prong->UncheckedAt(i2prong);
        if(fUseSelectionBit && dfromB->GetSelectionMap()) if(!dfromB->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
            isD0fromBplustagged=kFALSE;
        }
        
        if(isD0fromBplustagged && fWriteVariableTreeBplus){
            fNentries->Fill(25);
            nFilteredBplus++;
            if ((vHF->FillRecoCand(aod, dfromB))) { //Fill the data members of the candidate only if they are empty.

                Int_t isSelectedFilt=fFiltCutsBplustoD0pi->IsD0forD0ptbinSelectedMVA(dfromB, AliRDHFCuts::kAll, aod, vtx1_Bplus, bfield);
                if(isSelectedFilt > 0){

                  //loop over all tracks for pion from Bplus
                  for (Int_t iTrack = 0; iTrack < aod->GetNumberOfTracks(); ++iTrack){
                      
                    AliAODTrack* pionTrack = dynamic_cast<AliAODTrack*>(aod->GetTrack(iTrack));
                    if (!pionTrack) AliFatal("Not a standard AOD");

                    if(fFiltCutsBplustoD0pi->IsBplusPionSelectedMVA(pionTrack, AliRDHFCuts::kAll, aod, vtx1_Bplus, bfield) > 0){

                      //we check if the IDs of the tracks are different
                      AliAODTrack* twoProngdaughter0 = (AliAODTrack*)dfromB->GetDaughter(0);
                      AliAODTrack* twoProngdaughter1 = (AliAODTrack*)dfromB->GetDaughter(1);
                      UShort_t idProng0 = twoProngdaughter0->GetID();
                      UShort_t idProng1 = twoProngdaughter1->GetID();
                      if (pionTrack->GetID() != idProng0 && pionTrack->GetID() != idProng1){

                        //Pre reconstructing vertex cuts (to speed things up)
                        if(fFiltCutsBplustoD0pi->IsD0SelectedPreRecVtxMVA(dfromB,pionTrack,vtx1_Bplus,bfield,0) > 0){
                        
                          //we use the BPlus pion and D0 tracks to reconstruct the vertex for the BPlus
                          AliExternalTrackParam firstTrack;
                          firstTrack.CopyFromVTrack(pionTrack);
                          AliExternalTrackParam secondTrack;
                          secondTrack.CopyFromVTrack(dfromB);
                          
                          // we calculate the vertex of the mother candidate
                          TObjArray daughterTracks;
                          daughterTracks.Add(&firstTrack);
                          daughterTracks.Add(&secondTrack);
                          Double_t dispersion = 0;
                          AliAODVertex *vertexMother = ReconstructBplusVertex(vtx1_Bplus, &daughterTracks, bfield, dispersion);

                          if (vertexMother){ //check if calculation vertex Bplus succeeded

                            //use the new vertex to create the BPlus candidate
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
                            
                            firstTrack.PropagateToDCA(vtx1_Bplus, bfield, 100., d0z0, covd0z0);
                            d0[0] = d0z0[0];
                            d0err[0] = TMath::Sqrt(covd0z0[0]);
                            secondTrack.PropagateToDCA(vtx1_Bplus, bfield, 100., d0z0, covd0z0);
                            d0[1] = d0z0[0];
                            d0err[1] = TMath::Sqrt(covd0z0[0]);

                            Double_t dca = secondTrack.GetDCA(&firstTrack, bfield, xdummy, ydummy);
                            Short_t chargeMother = dfromB->Charge() + pionTrack->Charge();

                            AliAODRecoDecayHF2Prong trackBPlus(vertexMother, px, py, pz, d0, d0err, dca);
                            
                            trackBPlus.SetCharge(chargeMother);
                            trackBPlus.GetSecondaryVtx()->AddDaughter(pionTrack);
                            trackBPlus.GetSecondaryVtx()->AddDaughter(dfromB);
                            trackBPlus.SetPrimaryVtxRef((AliAODVertex*)aod->GetPrimaryVertex());
                            trackBPlus.SetProngIDs(2, id);

                            if(fFiltCutsBplustoD0pi->IsSelected(&trackBPlus, 0, aod) > 0){

                              Bool_t isSelAnCuts=kFALSE;
                              Int_t isSelectedAnalysis=fCutsBplustoD0pi->IsD0forD0ptbinSelectedMVA(dfromB, AliRDHFCuts::kAll, aod, vtx1_Bplus, bfield);
                              if(isSelectedAnalysis > 0) isSelAnCuts=kTRUE;

                              Bool_t isSelAnCutsBplus = kTRUE;
                              if(fCutsBplustoD0pi->IsD0forD0ptbinSelectedMVA(dfromB, AliRDHFCuts::kAll, aod, vtx1_Bplus, bfield) < 1) isSelAnCutsBplus = kFALSE;
                              if(fCutsBplustoD0pi->IsBplusPionSelectedMVA(pionTrack, AliRDHFCuts::kAll, aod, vtx1_Bplus, bfield) < 1) isSelAnCutsBplus = kFALSE;
                              if(fCutsBplustoD0pi->IsD0SelectedPreRecVtxMVA(dfromB,pionTrack,vtx1_Bplus,bfield,0) > 0) isSelAnCutsBplus = kFALSE;
                              if(fCutsBplustoD0pi->IsSelected(&trackBPlus, 0, aod) < 1) isSelAnCutsBplus = kFALSE;

                              fNentries->Fill(26);
                              nSelectedBplus++;
                                
                              Int_t labBplus = -1;
                              if (fReadMC) labBplus = trackBPlus.MatchToMCB2Prong(521,421,pdgDgBplustoD0piInt,pdgDgD0topiK,arrMC);

                              bool issignal = kFALSE;
                              bool isbkg =    kFALSE;
                              bool isFD =     kFALSE;
                              bool isprompt = kTRUE;
                              bool isrefl =   kFALSE;
                              Int_t masshypo = 0;
                                
                              if(labBplus >= 0) {issignal = kTRUE;}
                              else isbkg = kTRUE;
                                
                              if (fReadMC) {
                                fTreeHandlerBplus->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                              }
                              fTreeHandlerBplus->SetIsSelectedStd(isSelAnCutsBplus,isSelAnCutsBplus,isSelAnCutsBplus);
                              fTreeHandlerBplus->SetVariables(&trackBPlus, bfield, masshypo, fPIDresp);

                            } // end Bplus is selected filt
                          } // end calculation vertex Bplus
                          
                          delete vertexMother; vertexMother = nullptr;
                          
                        } // end pre RecVtx selection
                      } // end track ID check
                    } //end Bplus pion pre selection
                  } //end loop over pion track
                } //end D0 filt pre selection
            } else {
              fNentries->Fill(27); //monitor how often this fails
            }
        }//end Bplus

    }//end loop on candidates

    if(fWriteVariableTreeD0) fTreeHandlerD0->FillTree();
    if(fWriteVariableTreeBplus) fTreeHandlerBplus->FillTree();

    delete vHF;
    return;
}

//--------------------------------------------------------
void AliAnalysisTaskSEHFTreeCreator::Process3Prong(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield){
    
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

                        if(fWriteVariableTreeDs==1) {
                          if(isSelectedAnalysis&4) isSelAnCutsKKpi=kTRUE;
                          if(isSelectedAnalysis&8) isSelAnCutspiKK=kTRUE;
                          if(isSelectedPidAnalysis&4) isSelAnPidCutsKKpi=kTRUE;
                          if(isSelectedPidAnalysis&8) isSelAnPidCutspiKK=kTRUE;
                          if(isSelectedTopoAnalysis&4) isSelAnTopoCutsKKpi=kTRUE;
                          if(isSelectedTopoAnalysis&8) isSelAnTopoCutspiKK=kTRUE;
                        }
                        else if(fWriteVariableTreeDs==2) {
                          if(isSelectedAnalysis&16) isSelAnCutsKKpi=kTRUE;
                          if(isSelectedAnalysis&32) isSelAnCutspiKK=kTRUE;
                          if(isSelectedPidAnalysis&16) isSelAnPidCutsKKpi=kTRUE;
                          if(isSelectedPidAnalysis&32) isSelAnPidCutspiKK=kTRUE;
                          if(isSelectedTopoAnalysis&16) isSelAnTopoCutsKKpi=kTRUE;
                          if(isSelectedTopoAnalysis&32) isSelAnTopoCutspiKK=kTRUE;
                        }
                        else if(fWriteVariableTreeDs==3) {
                          if(isSelectedAnalysis&1) isSelAnCutsKKpi=kTRUE;
                          if(isSelectedAnalysis&2) isSelAnCutspiKK=kTRUE;
                          if(isSelectedPidAnalysis&1) isSelAnPidCutsKKpi=kTRUE;
                          if(isSelectedPidAnalysis&2) isSelAnPidCutspiKK=kTRUE;
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
                        
                        //checking origin
                        if(fReadMC){
                            AliAODMCParticle *partDs = 0x0;
                            labDs = ds->MatchToMC(431,arrMC,3,pdgDstoKKpi);
                            labDplus = ds->MatchToMC(411,arrMC,3,pdgDstoKKpi);
                            
                            if(labDs>=0){
                                Int_t labDau0=((AliAODTrack*)ds->GetDaughter(0))->GetLabel();
                                AliAODMCParticle* p=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDau0));
                                pdgCode0=TMath::Abs(p->GetPdgCode());
                                partDs = (AliAODMCParticle*)arrMC->At(labDs);
                                
                            }
                            else{
                                if(labDplus>=0) {
                                    Int_t labDau0=((AliAODTrack*)ds->GetDaughter(0))->GetLabel();
                                    AliAODMCParticle* p=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDau0));
                                    pdgCode0=TMath::Abs(p->GetPdgCode());
                                    partDs = (AliAODMCParticle*)arrMC->At(labDplus);
                                }
                            }
                            if(partDs) orig = AliVertexingHFUtils::CheckOrigin(arrMC,partDs,kTRUE);
                        }
                        
                        //filling the Ds tree
                        if ((fWriteVariableTreeDs==1 && (isPhiKKpi || isPhipiKK)) || (fWriteVariableTreeDs==2 && (isK0starKKpi || isK0starpiKK)) || (fWriteVariableTreeDs==3 && (isKKpi || ispiKK))){
                            
                            bool issignal = kFALSE;
                            bool isbkg = kFALSE;
                            bool isprompt = kFALSE;
                            bool isFD = kFALSE;
                            bool isrefl = kFALSE;
                            
                            if(isKKpi || isPhiKKpi || isK0starKKpi) {
                              if(fReadMC) {
                                if(labDs>=0) {
                                  if(pdgCode0==321) issignal = kTRUE;
                                  else if(pdgCode0==211) isrefl = kTRUE;
                                  if(orig==4) isprompt = kTRUE;
                                  else if(orig==5) isFD = kTRUE;
                                }
                                else {
                                  isbkg = kTRUE;
                                  if(labDplus>=0) fTreeHandlerDs->SetIsDplustoKKpi(kTRUE);//put also D+ -->KKpi in bkg
                                }
                                //do not apply cuts, but enable flag if is selected
                                fTreeHandlerDs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                              }
                              fTreeHandlerDs->SetIsSelectedStd(isSelAnCutsKKpi,isSelAnTopoCutsKKpi,isSelAnPidCutsKKpi);
                              fTreeHandlerDs->SetVariables(ds,bfield,0,fPIDresp);
                            }
                          issignal = kFALSE;
                          isbkg = kFALSE;
                          isprompt = kFALSE;
                          isFD = kFALSE;
                          isrefl = kFALSE;
                          if(ispiKK || isPhipiKK || isK0starpiKK) {
                            if(fReadMC) {
                              if(labDs>=0) {
                                if(pdgCode0==211) issignal = kTRUE;
                                else if(pdgCode0==321) isrefl = kTRUE;
                                if(orig==4) isprompt = kTRUE;
                                else if(orig==5) isFD = kTRUE;
                              }
                              else {
                                isbkg = kTRUE;
                                if(labDplus>=0) fTreeHandlerDs->SetIsDplustoKKpi(kTRUE);//put also D+ -->KKpi in bkg
                              }
                              //do not apply cuts, but enable flag if is selected
                                fTreeHandlerDs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                              }
                            fTreeHandlerDs->SetIsSelectedStd(isSelAnCutspiKK,isSelAnTopoCutspiKK,isSelAnPidCutspiKK);
                            fTreeHandlerDs->SetVariables(ds,bfield,1,fPIDresp);
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

                  Bool_t unsetvtx=kFALSE;
                  if(!dplus->GetOwnPrimaryVtx()){
                  dplus->SetOwnPrimaryVtx(vtx1);
                  unsetvtx=kTRUE;
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
                  //read MC
                  if(fReadMC){
                  labDp = dplus->MatchToMC(411,arrMC,3,pdgDgDplustoKpipi);
                  if(labDp>=0){
                    issignal=kTRUE;
                    AliAODMCParticle *partDp = (AliAODMCParticle*)arrMC->At(labDp);
                    Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partDp,kTRUE);//Prompt = 4, FeedDown = 5
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
                    else isbkg=kTRUE;
                    fTreeHandlerDplus->SetCandidateType(issignal,isbkg,isPrimary,isFeeddown,kFALSE);
                   } //end read MC
                   
                   // fill tree
                    fTreeHandlerDplus->SetIsSelectedStd(isSelAnCuts,isSelectedTopoAnalysis,isSelectedPidAnalysis);
                    fTreeHandlerDplus->SetVariables(dplus,bfield,0,fPIDresp);
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
            nFilteredLctopKpi++;
            fNentries->Fill(22);
            if((vHF->FillRecoCand(aod,lctopkpi))) {////Fill the data members of the candidate only if they are empty.
              
                Int_t isSelectedFilt    = fFiltCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
                //Printf("isSelectedFilt = %i isSelectedAnalysis = %i",isSelectedFilt,isSelectedAnalysis);
                if(isSelectedFilt){
                  fNentries->Fill(23);
                  nSelectedLctopKpi++;

                  // check analysis cuts
                  Bool_t isSelAnCuts=kFALSE;
                  Bool_t isSelPID=kFALSE;
                  Bool_t isSelTopo=kFALSE;
                  Bool_t ispKpi=kFALSE;
                  Bool_t ispiKp=kFALSE;
                  Int_t isSelectedAnalysis= fCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
                  Int_t isSelectedTopoAnalysis = fCutsLctopKpi->IsSelected(lctopkpi,AliRDHFCuts::kAll,aod);
                  Int_t isSelectedPidAnalysis = fCutsLctopKpi->IsSelectedPID(lctopkpi);
                  if(isSelectedAnalysis) isSelAnCuts=kTRUE;
                  if(isSelectedTopoAnalysis) isSelTopo=kTRUE;
                  if(isSelectedPidAnalysis) isSelPID=kTRUE;
                  if(isSelectedFilt==1 || isSelectedFilt==3)     ispKpi=kTRUE;
                  if(isSelectedFilt>2)                           ispiKp=kTRUE;

                  Bool_t unsetvtx=kFALSE;
                  if(!lctopkpi->GetOwnPrimaryVtx()){
                  lctopkpi->SetOwnPrimaryVtx(vtx1);
                  unsetvtx=kTRUE;
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
                  if(ispKpi) {
                    //read MC
                    if(fReadMC){
                      labDp = lctopkpi->MatchToMC(4122,arrMC,3,pdgLctopKpi);
                      if(labDp>=0){
                        issignal=kTRUE;
                        AliAODMCParticle *partDp = (AliAODMCParticle*)arrMC->At(labDp);
                        Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partDp,kTRUE);//Prompt = 4, FeedDown = 5
                        if(orig==4){
                          isPrimary=kTRUE;
                          isFeeddown=kFALSE;
                        }
                        else if(orig==5){
                          isPrimary=kFALSE;
                          isFeeddown=kTRUE;
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
                      }
                      else isbkg=kTRUE;
                      fTreeHandlerLctopKpi->SetCandidateType(issignal,isbkg,isPrimary,isFeeddown,isrefl);
                      //Printf("labLc = %i, issignal = %i, isPrimary = %i, isFeeddown = %i, isBkg = %i",labDp,issignal,isPrimary,isFeeddown,isbkg);
                    } //end read MC

                    // fill tree
                    fTreeHandlerLctopKpi->SetIsSelectedStd(isSelAnCuts,isSelTopo,isSelPID);
                    fTreeHandlerLctopKpi->SetVariables(lctopkpi,bfield,1,fPIDresp);
                  } // end pKpi
                  isPrimary=kFALSE;
                  isFeeddown=kFALSE;
                  issignal=kFALSE;
                  isbkg=kFALSE;
                  isrefl=kFALSE;
                  labDp=-1;
                  if(ispiKp) {
                    //read MC
                    if(fReadMC){
                      labDp = lctopkpi->MatchToMC(4122,arrMC,3,pdgLctopKpi);
                      if(labDp>=0){
                        issignal=kTRUE;
                        AliAODMCParticle *partDp = (AliAODMCParticle*)arrMC->At(labDp);
                        Int_t orig=AliVertexingHFUtils::CheckOrigin(arrMC,partDp,kTRUE);//Prompt = 4, FeedDown = 5
                        if(orig==4){
                          isPrimary=kTRUE;
                          isFeeddown=kFALSE;
                        }
                        else if(orig==5){
                          isPrimary=kFALSE;
                          isFeeddown=kTRUE;
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
                      }
                      else isbkg=kTRUE;
                      fTreeHandlerLctopKpi->SetCandidateType(issignal,isbkg,isPrimary,isFeeddown,isrefl);
                      //Printf("labLc = %i, issignal = %i, isPrimary = %i, isFeeddown = %i, isBkg = %i",labDp,issignal,isPrimary,isFeeddown,isbkg);
                    } //end read MC

                    // fill tree
                    fTreeHandlerLctopKpi->SetIsSelectedStd(isSelAnCuts,isSelTopo,isSelPID);
                    fTreeHandlerLctopKpi->SetVariables(lctopkpi,bfield,2,fPIDresp);
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
    
    if(fWriteVariableTreeDs) fTreeHandlerDs->FillTree();
    if(fWriteVariableTreeDplus) fTreeHandlerDplus->FillTree();
    if(fWriteVariableTreeLctopKpi) fTreeHandlerLctopKpi->FillTree();

    delete vHF;
    return;
}

//________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::ProcessCasc(TClonesArray *arrayCasc, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield){
    
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    
    Int_t nCasc = arrayCasc->GetEntriesFast();
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
        AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*)arrayCasc->UncheckedAt(iCasc);
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
                    if(isSelectedAnalysis) isSelAnCuts = kTRUE;
                    if(isSelectedPidAnalysis) isSelAnPidCuts = kTRUE;
                    if(isSelectedTopoAnalysis) isSelAnTopolCuts = kTRUE;
                    fCutsDstartoKpipi->SetUsePID(isUsePidAn);

                    Bool_t unsetvtx=kFALSE;
                    if(!d->GetOwnPrimaryVtx()){
                        d->SetOwnPrimaryVtx(vtx1);
                        unsetvtx=kTRUE;
                        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
                        // Pay attention if you use continue inside this loop!!!
                    }
                    Bool_t recVtx=kFALSE;
                    
                    AliAODVertex *origownvtx=0x0;
                    if(fFiltCutsDstartoKpipi->GetIsPrimaryWithoutDaughters()){
                        if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
                        if(fFiltCutsDstartoKpipi->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
                        else fFiltCutsDstartoKpipi->CleanOwnPrimaryVtx(d,aod,origownvtx);
                    }
            
                    Int_t labDstar = -1;
                    Int_t pdgDstar = -99;
                    Int_t origin= -1;
                    AliAODMCParticle *partDstar=0x0;

                    if(fReadMC) {
                        labDstar = d->MatchToMC(413,421,pdgDgDStartoD0pi, pdgDgD0toKpi, arrMC);
                        if(labDstar>=0){
                            partDstar = (AliAODMCParticle*)arrMC->At(labDstar);
                            pdgDstar = partDstar->GetPdgCode();
                            origin = AliVertexingHFUtils::CheckOrigin(arrMC,partDstar,kTRUE);
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
                            if(origin==4) isprompt=kTRUE;
                            else if(origin==5) isFD=kTRUE;
                            if(pdgDstar==413){
                                issignal=kTRUE;
                            }
                        }//end labDstar check
                        else{//background
                            isbkg=kTRUE;
                        }
                        fTreeHandlerDstar->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                    }//end read MC
                    fTreeHandlerDstar->SetIsSelectedStd(isSelAnCuts,isSelectedTopoAnalysis,isSelectedPidAnalysis);
                    fTreeHandlerDstar->SetVariables(d,bfield,masshypo,fPIDresp);
	    
                    if(recVtx)fFiltCutsDstartoKpipi->CleanOwnPrimaryVtx(d,aod,origownvtx);
                    if(unsetvtx) d->UnsetOwnPrimaryVtx();
                }//end is selected filt
            }
            else {
                fNentries->Fill(31); //monitor how often this fails
            }
        }//end Dstar
	}//end loop on candidates
  
    if(fWriteVariableTreeDstar) fTreeHandlerDstar->FillTree();

    delete vHF;
    return;
}

//_________________________________________________________________
void AliAnalysisTaskSEHFTreeCreator::ProcessMCGen(TClonesArray *arrayMC){
  /// Fill MC gen trees
  
  for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
      
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
    Int_t absPDG = TMath::Abs(mcPart->GetPdgCode());
    
      if(absPDG == 411 || absPDG == 421 || absPDG == 431 || absPDG == 4122 || absPDG == 521 || absPDG == 413) {
        Bool_t isPrimary = kFALSE;
        Bool_t isFeeddown = kFALSE;
        //Bplus will always end up with orig=4, so primary
        Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,kTRUE);//Prompt = 4, FeedDown = 5
        if(orig==4){
          isPrimary = kTRUE;
          isFeeddown = kFALSE;
        }
        else if(orig==5){
          isPrimary = kFALSE;
          isFeeddown = kTRUE;
        }

        Int_t  deca = 0;
        Int_t  labDau[3] = {-1,-1,-1};
        Bool_t isDaugInAcc = kFALSE;

        if(absPDG == 411 && fWriteVariableTreeDplus) {
          deca = AliVertexingHFUtils::CheckDplusDecay(arrayMC,mcPart,labDau);
          if(deca<1 || labDau[0]<0 || labDau[1]<0) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau);
          fTreeHandlerGenDplus->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenDplus->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenDplus->SetMCGenVariables(mcPart);
        }
        else if(absPDG == 421 && fWriteVariableTreeD0) {
          deca = AliVertexingHFUtils::CheckD0Decay(arrayMC,mcPart,labDau);
          if(deca!=1 || labDau[0]<0 || labDau[1]<0) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,2,labDau);
          fTreeHandlerGenD0->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenD0->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenD0->SetMCGenVariables(mcPart);
        }
        else if(absPDG == 431 && fWriteVariableTreeDs) {
          deca = AliVertexingHFUtils::CheckDsDecay(arrayMC,mcPart,labDau);
          if(deca!=1 || labDau[0]<0 || labDau[1]<0) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau);
          fTreeHandlerGenDs->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenDs->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenDs->SetMCGenVariables(mcPart);
        }
        if(absPDG == 4122 && fWriteVariableTreeLctopKpi) {
          deca = AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC,mcPart,labDau);
          if(deca<1 || labDau[0]==-1 || labDau[1]<0) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau);
          fTreeHandlerGenLctopKpi->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenLctopKpi->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenLctopKpi->SetMCGenVariables(mcPart);
        }
        else if(absPDG == 521 && fWriteVariableTreeBplus) {
          deca = AliVertexingHFUtils::CheckBplusDecay(arrayMC,mcPart,labDau);
          if(deca!=1 || labDau[0]==-1 || labDau[1]<0) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau);
          fTreeHandlerGenBplus->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenBplus->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenBplus->SetMCGenVariables(mcPart);
        }
        else if(absPDG == 413 && fWriteVariableTreeDstar) {
          deca = AliVertexingHFUtils::CheckDstarDecay(arrayMC,mcPart,labDau);
          if(deca!=1 || labDau[0]<0 || labDau[1]<0) continue;
          isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau);
          fTreeHandlerGenDstar->SetDauInAcceptance(isDaugInAcc);
          fTreeHandlerGenDstar->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
          fTreeHandlerGenDstar->SetMCGenVariables(mcPart);
        }
      }
    }
  
  if(fWriteVariableTreeD0) fTreeHandlerGenD0->FillTree();
  if(fWriteVariableTreeDs) fTreeHandlerGenDs->FillTree();
  if(fWriteVariableTreeDplus) fTreeHandlerGenDplus->FillTree();
  if(fWriteVariableTreeLctopKpi) fTreeHandlerGenLctopKpi->FillTree();
  if(fWriteVariableTreeBplus) fTreeHandlerGenBplus->FillTree();
  if(fWriteVariableTreeDstar) fTreeHandlerGenDstar->FillTree();
}

//--------------------------------------------------------
Bool_t AliAnalysisTaskSEHFTreeCreator::CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau){
  /// check if the decay products are in the good eta and pt range
  
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) {
      return kFALSE;
    }
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > 0.9 || pt < 0.1) {
      return kFALSE;
    }
  }
  return kTRUE;
}

//________________________________________________________________
AliAODVertex* AliAnalysisTaskSEHFTreeCreator::ReconstructBplusVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion) {
    //
    // Helper function to recalculate a vertex.
    //
    
    AliESDVertex *vertexESD = 0;
    AliAODVertex *vertexAOD = 0;
    
    AliVertexerTracks vertexer;
    vertexer.SetFieldkG(bField);
    
    vertexer.SetVtxStart((AliESDVertex*)primary); //primary vertex
    vertexESD = (AliESDVertex*)vertexer.VertexForSelectedESDTracks(tracks);
    
    // delete vertexer; vertexer=NULL;
    
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
    Int_t nprongs = 2; //tracks->GetEntriesFast();
    vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, nprongs);
    
    return vertexAOD;
}
