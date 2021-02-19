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
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TString.h>
#include <TObjString.h>
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
#include "AliHFTreeHandlerApply.h"
#include "AliHFTreeHandlerApplyDstoKKpi.h"
#include "AliHFTreeHandlerApplyLc2V0bachelor.h"
#include "AliHFMLResponseDstoKKpi.h"
#include "AliHFMLResponseLctoV0bachelor.h"
#include "AliAnalysisTaskSEHFTreeCreatorApply.h"
#include "AliAODPidHF.h"
#include "AliESDUtils.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerInput.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEHFTreeCreatorApply);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreatorApply::AliAnalysisTaskSEHFTreeCreatorApply():
AliAnalysisTaskSEHFTreeCreatorApply("", nullptr)
{

  /// Default constructor

}
//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreatorApply::AliAnalysisTaskSEHFTreeCreatorApply(const char *name, TList *cutsList):
AliAnalysisTaskSE(name),
fEventNumber(0),
fNentries(0x0),
fHistoNormCounter(0x0),
fListCounter(0x0),
fCounter(0x0),
fListCuts(cutsList),
fFiltCutsDstoKKpi(0x0),
fFiltCutsLc2V0bachelor(0x0),
fCutsDstoKKpi(0x0),
fCutsLc2V0bachelor(0x0),
fEvSelectionCuts(0x0),
fReadMC(0),
fUseSelectionBit(kTRUE),
fSys(0),
fAODProtection(0),
fWriteOnlySignal(kFALSE),
fFillMCGenTrees(kTRUE),
fWriteVariableTreeDs(0),
fWriteVariableTreeLc2V0bachelor(0),
fVariablesTreeDs(0x0),
fVariablesTreeLc2V0bachelor(0x0),
fGenTreeDs(0x0),
fGenTreeLc2V0bachelor(0x0),
fTreeEvChar(0x0),
fTreeHandlerDs(0x0),
fTreeHandlerLc2V0bachelor(0x0),
fTreeHandlerGenDs(0x0),
fTreeHandlerGenLc2V0bachelor(0x0),
fPIDresp(0x0),
fPIDoptDs(AliHFTreeHandlerApply::kRawAndNsigmaPID),
fPIDoptLc2V0bachelor(AliHFTreeHandlerApply::kRawAndNsigmaPID),
fBC(0),
fOrbit(0),
fPeriod(0),
fEventID(0),
fEventIDExt(0),
fEventIDLong(0),
fFileName(""),
fDirNumber(0),
fCentrality(-999.),
fzVtxReco(0.),
fzVtxGen(0.),
fNcontributors(0),
fNtracks(0),
fIsEvRej(0),
fIsEvRej_INT7(0),
fIsEvRej_HighMultSPD(0),
fIsEvRej_HighMultV0(0),
fRunNumber(0),
fRunNumberCDB(0),
fnTracklets(0),
fnTrackletsCorr(0),
fnTrackletsCorrSHM(0),
fnV0A(0),
fnTPCCls(0),
fMultGen(0),
fMultGenV0A(0),
fMultGenV0C(0),
fTriggerMask(0),
fTriggerOnlineINT7(false),
fTriggerOnlineHighMultSPD(false),
fTriggerOnlineHighMultV0(false),
fTriggerBitINT7(false),
fTriggerBitHighMultSPD(false),
fTriggerBitHighMultV0(false),
fTriggerBitCentral(false),
fTriggerBitSemiCentral(false),
fTriggerClasses(""),
fTriggerClassINT7(false),
fTriggerClassHighMultSPD(false),
fTriggerClassHighMultV0m(false),
fnV0M(0),
fnV0MEq(0),
fnV0MCorr(0),
fnV0MEqCorr(0),
fPercV0M(0.),
fMultV0M(0.),
fRefMult(9.26),
fRefMultSHM(9.26),
fMultEstimatorAvg(),
fMultEstimatorAvgSHM(),
fCorrNtrVtx(false),
fCorrV0MVtx(false),
fCdbEntry(nullptr),
fDsMassKKOpt(1),
fLc2V0bachelorCalcSecoVtx(0),
fV0typeForLc2V0bachelor(1),
fTreeSingleTrackVarsOpt(AliHFTreeHandlerApply::kRedSingleTrackVars),
fGoodTrackFilterBit(-1),
fGoodTrackEtaRange(999.),
fGoodTrackMinPt(0.),
fITSUpgradeProduction(0),
fITSUpgradePreSelect(0),
fEnableNsigmaTPCDataCorr(false),
fSystemForNsigmaTPCDataCorr(AliAODPidHF::kNone),
fApplyPhysicsSelOnline(false),
fApplyEventSelOnline(false),
fEnableEventDownsampling(false),
fFracToKeepEventDownsampling(1.1),
fSeedEventDownsampling(0),
fEnableCandDownsampling(false),
fFracToKeepCandDownsampling(1.1),
fMaxPtCandDownsampling(999.),
fConfigPath(""),
fMLResponse(0x0),
fReducePbPbBranches(false),
fSaveSTDSelection(false)
{

  if (fListCuts) {
    fFiltCutsDstoKKpi     =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("DstoKKpiFilteringCuts");
    fFiltCutsLc2V0bachelor=(AliRDHFCutsLctoV0*)fListCuts->FindObject("Lc2V0bachelorFilteringCuts");
    fCutsDstoKKpi         =(AliRDHFCutsDstoKKpi*)fListCuts->FindObject("DstoKKpiAnalysisCuts");
    fCutsLc2V0bachelor    =(AliRDHFCutsLctoV0*)fListCuts->FindObject("Lc2V0bachelorAnalysisCuts");

    if(fWriteVariableTreeDs) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDstoKKpi->Clone();
    else if(fWriteVariableTreeLc2V0bachelor) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLc2V0bachelor->Clone();
    else {
      if(fFiltCutsDstoKKpi) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsDstoKKpi->Clone();
      else if(fFiltCutsLc2V0bachelor) fEvSelectionCuts = (AliRDHFCuts*)fFiltCutsLc2V0bachelor->Clone();
      else printf("AliAnalysisTaskSEHFTreeCreatorApply:: Constructor: No event selection cuts could be stored, code will crash!\n");
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
  // Output slot #6 stores the tree of the Ds+ candidate variables after track selection
  DefineOutput(6,TTree::Class());
  // Output slot #7 stores the tree of the gen Ds+ variables
  DefineOutput(7,TTree::Class());
  // Output slot #8 stores the tree of the Lc2V0bachelor candidate variables after track selection
  DefineOutput(8,TTree::Class());
  // Output slot #9 stores the tree of the gen Lc2V0bachelor variables
  DefineOutput(9,TTree::Class());

}

//________________________________________________________________________
AliAnalysisTaskSEHFTreeCreatorApply::~AliAnalysisTaskSEHFTreeCreatorApply()
{
  delete fListCuts;
  delete fFiltCutsDstoKKpi;
  delete fFiltCutsLc2V0bachelor;
  delete fCutsDstoKKpi;
  delete fCutsLc2V0bachelor;
  delete fEvSelectionCuts;
  delete fNentries;
  delete fHistoNormCounter;
  delete fListCounter;
  delete fCounter;
  delete fTreeHandlerDs;
  delete fTreeHandlerLc2V0bachelor;

  delete fTreeHandlerGenDs;
  delete fTreeHandlerGenLc2V0bachelor;
  delete fTreeEvChar;
  
  delete fMLResponse;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreatorApply::Init()
{
  /// Initialization

  if(fDebug > 1) printf("AliAnalysisTaskSEHFTreeCreatorApply::Init() \n");

  AliInfo(Form(" Resetting NsigmaDataDrivenCorrection of cutobjects to %d for system %d \n", fEnableNsigmaTPCDataCorr,fSystemForNsigmaTPCDataCorr));
  if(fFiltCutsDstoKKpi) fFiltCutsDstoKKpi->EnableNsigmaDataDrivenCorrection(fEnableNsigmaTPCDataCorr,fSystemForNsigmaTPCDataCorr);
  if(fFiltCutsLc2V0bachelor) fFiltCutsLc2V0bachelor->EnableNsigmaDataDrivenCorrection(fEnableNsigmaTPCDataCorr,fSystemForNsigmaTPCDataCorr);
  if(fCutsDstoKKpi) fCutsDstoKKpi->EnableNsigmaDataDrivenCorrection(fEnableNsigmaTPCDataCorr,fSystemForNsigmaTPCDataCorr);
  if(fCutsLc2V0bachelor) fCutsLc2V0bachelor->EnableNsigmaDataDrivenCorrection(fEnableNsigmaTPCDataCorr,fSystemForNsigmaTPCDataCorr);
  
  PostData(3,fListCuts);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreatorApply::UserCreateOutputObjects()
{

  /// Create the output container

  if(fDebug > 1) printf("AliAnalysisTaskSEHFTreeCreatorApply::UserCreateOutputObjects() \n");

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
  fNentries->GetXaxis()->SetBinLabel(43,"n. evt. rejected by downsampling");
  
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
  if(fWriteVariableTreeDs) nEnabledTrees++;
  if(fWriteVariableTreeLc2V0bachelor) nEnabledTrees++;
  if(fReadMC && fFillMCGenTrees) {
    nEnabledTrees = (nEnabledTrees-1)*2+1;
  }
  
  //
  // Output slot 4-9 : trees of the candidate and event-characterization variables
  //
  OpenFile(5);
  fTreeEvChar = new TTree("tree_event_char","tree_event_char");
  //set variables
  fTreeEvChar->Branch("centrality", &fCentrality);
  fTreeEvChar->Branch("z_vtx_reco", &fzVtxReco);
  fTreeEvChar->Branch("n_vtx_contributors", &fNcontributors);
  fTreeEvChar->Branch("n_tracks", &fNtracks);
  fTreeEvChar->Branch("is_ev_rej", &fIsEvRej);
  if(!fReducePbPbBranches){
    fTreeEvChar->Branch("is_ev_rej_INT7", &fIsEvRej_INT7);
    fTreeEvChar->Branch("is_ev_rej_HighMultSPD", &fIsEvRej_HighMultSPD);
    fTreeEvChar->Branch("is_ev_rej_HighMultV0", &fIsEvRej_HighMultV0);
  }
  fTreeEvChar->Branch("run_number", &fRunNumber);
  fTreeEvChar->Branch("ev_id", &fEventID);
  fTreeEvChar->Branch("ev_id_ext", &fEventIDExt);
  fTreeEvChar->Branch("ev_id_long", &fEventIDLong);
  fTreeEvChar->Branch("n_tracklets", &fnTracklets);
  fTreeEvChar->Branch("V0Amult", &fnV0A);
  fTreeEvChar->Branch("n_tpc_cls", &fnTPCCls);
  if(!fReducePbPbBranches){
    fTreeEvChar->Branch("trigger_bitmap", &fTriggerMask);
    fTreeEvChar->Branch("trigger_online_INT7", &fTriggerOnlineINT7);
    fTreeEvChar->Branch("trigger_online_HighMultSPD", &fTriggerOnlineHighMultSPD);
    fTreeEvChar->Branch("trigger_online_HighMultV0", &fTriggerOnlineHighMultV0);
  }
  fTreeEvChar->Branch("trigger_hasbit_INT7", &fTriggerBitINT7);
  if(!fReducePbPbBranches){
    fTreeEvChar->Branch("trigger_hasbit_HighMultSPD", &fTriggerBitHighMultSPD);
    fTreeEvChar->Branch("trigger_hasbit_HighMultV0", &fTriggerBitHighMultV0);
  }
  fTreeEvChar->Branch("trigger_hasbit_Central", &fTriggerBitCentral);
  fTreeEvChar->Branch("trigger_hasbit_SemiCentral", &fTriggerBitSemiCentral);
  if(!fReducePbPbBranches){
    fTreeEvChar->Branch("trigger_classes", &fTriggerClasses);
    fTreeEvChar->Branch("trigger_hasclass_INT7", &fTriggerClassINT7);
    fTreeEvChar->Branch("trigger_hasclass_HighMultSPD", &fTriggerClassHighMultSPD);
    fTreeEvChar->Branch("trigger_hasclass_HighMultV0", &fTriggerClassHighMultV0m);
  }
  fTreeEvChar->Branch("z_vtx_gen", &fzVtxGen);
  if(!fReducePbPbBranches){
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
  }
  fTreeEvChar->SetMaxVirtualSize(1.e+8/nEnabledTrees);
  
  if(fWriteVariableTreeDs){
    if(fConfigPath != ""){
      fMLResponse = new AliHFMLResponseDstoKKpi("DstoKKpiMLResponse", "DstoKKpiMLResponse", fConfigPath.Data());
      fMLResponse->MLResponseInit();
    }

    OpenFile(6);
    TString nameoutput = "tree_Ds";
    fTreeHandlerDs = new AliHFTreeHandlerApplyDstoKKpi(fPIDoptDs);
    fTreeHandlerDs->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerDs->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerDs->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerDs->SetMassKKOption(fDsMassKKOpt);
    fVariablesTreeDs = (TTree*)fTreeHandlerDs->BuildTree(nameoutput,nameoutput);
    fVariablesTreeDs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeDs);
    
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(7);
      TString nameoutput = "tree_Ds_gen";
      fTreeHandlerGenDs = new AliHFTreeHandlerApplyDstoKKpi(0);
      fGenTreeDs = (TTree*)fTreeHandlerGenDs->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeDs->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeDs);
    }
  }
  if(fWriteVariableTreeLc2V0bachelor){
    if(fConfigPath != ""){
      fMLResponse = new AliHFMLResponseLctoV0bachelor("LctoV0bachelorMLResponse", "LctoV0bachelorMLResponse", fConfigPath.Data());
      fMLResponse->MLResponseInit();
    }

    OpenFile(8);
    TString nameoutput = "tree_Lc2V0bachelor";
    fTreeHandlerLc2V0bachelor = new AliHFTreeHandlerApplyLc2V0bachelor(fPIDoptLc2V0bachelor);
    fTreeHandlerLc2V0bachelor->SetOptSingleTrackVars(fTreeSingleTrackVarsOpt);
    if(fReadMC && fWriteOnlySignal) fTreeHandlerLc2V0bachelor->SetFillOnlySignal(fWriteOnlySignal);
    if(fEnableNsigmaTPCDataCorr) fTreeHandlerLc2V0bachelor->EnableNsigmaTPCDataDrivenCorrection(fSystemForNsigmaTPCDataCorr);
    fTreeHandlerLc2V0bachelor->SetCalcSecoVtx(fLc2V0bachelorCalcSecoVtx);
    fTreeHandlerLc2V0bachelor->SetReducePbPbBranches(fReducePbPbBranches);
    fVariablesTreeLc2V0bachelor = (TTree*)fTreeHandlerLc2V0bachelor->BuildTree(nameoutput,nameoutput);
    fVariablesTreeLc2V0bachelor->SetMaxVirtualSize(1.e+8/nEnabledTrees);
    fTreeEvChar->AddFriend(fVariablesTreeLc2V0bachelor);
    if(fFillMCGenTrees && fReadMC) {
      OpenFile(9);
      TString nameoutput = "tree_Lc2V0bachelor_gen";
      fTreeHandlerGenLc2V0bachelor = new AliHFTreeHandlerApplyLc2V0bachelor(0);
      fGenTreeLc2V0bachelor = (TTree*)fTreeHandlerGenLc2V0bachelor->BuildTreeMCGen(nameoutput,nameoutput);
      fGenTreeLc2V0bachelor->SetMaxVirtualSize(1.e+8/nEnabledTrees);
      fTreeEvChar->AddFriend(fGenTreeLc2V0bachelor);
    }
  }

  //Set seed of gRandom
  if(fEnableEventDownsampling) gRandom->SetSeed(fSeedEventDownsampling);

  // Post the data
  PostData(1,fNentries);
  PostData(2,fHistoNormCounter);
  PostData(4,fListCounter);
  PostData(5,fTreeEvChar);
  if(fWriteVariableTreeDs){
    PostData(6,fVariablesTreeDs);
    if(fFillMCGenTrees && fReadMC) PostData(7,fGenTreeDs);
  }
  if(fWriteVariableTreeLc2V0bachelor){
    PostData(8,fVariablesTreeLc2V0bachelor);
    if(fFillMCGenTrees && fReadMC) PostData(9,fGenTreeLc2V0bachelor);
  }

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreatorApply::UserExec(Option_t */*option*/){

  /// Execute analysis for current event:

  if(fEnableEventDownsampling && gRandom->Rndm() > fFracToKeepEventDownsampling) {
    fNentries->Fill(42);
    return;
  }
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  fBC = aod->GetBunchCrossNumber();
  fOrbit = (Int_t)(aod->GetOrbitNumber());
  fPeriod = (Int_t)(aod->GetPeriodNumber());
  fEventIDLong = (Long64_t)(aod->GetHeader()->GetEventIdAsLong());
  if (!fEventIDLong)
    fEventIDLong = (Long64_t(aod->GetTimeStamp()) << 32) + Long64_t((aod->GetNumberOfTPCClusters()<<5) | (aod->GetNumberOfTPCTracks()));
  fEventIDExt = Int_t(fEventIDLong >> 32);
  if(!fReadMC) {
    fEventID    = Int_t(fEventIDLong & 0xffffffff);
  } else {
    fEventID    = Int_t(GetEvID());
  }
  
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
    printf("AliAnalysisTaskSEHFTreeCreatorApply::UserExec: input branches not found!\n");
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
      printf("AliAnalysisTaskSEHFTreeCreatorApply::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEHFTreeCreatorApply::UserExec: MC header branch not found!\n");
      return;
    }
    fzVtxGen = mcHeader->GetVtxZ();
  }
  
  Bool_t isSameEvSelDs=kTRUE;
  Bool_t isSameEvSelLc2V0bachelor=kTRUE;

  if(fWriteVariableTreeDs)
    isSameEvSelDs=!((fFiltCutsDstoKKpi->IsEventSelected(aod) && !fCutsDstoKKpi->IsEventSelected(aod))||(!fFiltCutsDstoKKpi->IsEventSelected(aod) && fCutsDstoKKpi->IsEventSelected(aod)));
  if(fWriteVariableTreeLc2V0bachelor)
    isSameEvSelLc2V0bachelor=!((fFiltCutsLc2V0bachelor->IsEventSelected(aod) && !fCutsLc2V0bachelor->IsEventSelected(aod))||(!fFiltCutsLc2V0bachelor->IsEventSelected(aod) && fCutsLc2V0bachelor->IsEventSelected(aod)));
  
  Bool_t isSameEvSel = isSameEvSelDs && isSameEvSelLc2V0bachelor;
  if(!isSameEvSel) {
    Printf("AliAnalysisTaskSEHFTreeCreatorApply::UserExec: differences in the event selection cuts same meson");
    return;
  }
  
  if((fWriteVariableTreeDs && fWriteVariableTreeLc2V0bachelor && (fFiltCutsDstoKKpi->IsEventSelected(aod)!=fFiltCutsLc2V0bachelor->IsEventSelected(aod)))){
    Printf("AliAnalysisTaskSEHFTreeCreatorApply::UserExec: differences in the event selection cuts different meson");
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
  fIsEvRej_INT7 = fEvSelectionCuts->GetEventRejectionBitMap();
  
  fEvSelectionCuts->SetTriggerMask(AliVEvent::kHighMultSPD);
  fIsEvRej_HighMultSPD = fEvSelectionCuts->GetEventRejectionBitMap();
  
  fEvSelectionCuts->SetTriggerMask(AliVEvent::kHighMultV0);
  fIsEvRej_HighMultV0 = fEvSelectionCuts->GetEventRejectionBitMap();
  
  fEvSelectionCuts->SetTriggerMask(trig_mask_cuts);
    
  //TPC multiplicities
  fnTPCCls = aod->GetNumberOfTPCClusters();

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
  
  fTriggerClasses = aod->GetFiredTriggerClasses();
  fTriggerClassINT7 = fTriggerClasses.Contains("CINT7-B");
  fTriggerClassHighMultSPD = fTriggerClasses.Contains("CVHMSH2-B");
  fTriggerClassHighMultV0m = fTriggerClasses.Contains("CVHMV0M-B");
  
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
  const auto triggerBits = aod->GetHeader()->GetL0TriggerInputs();
  fTriggerOnlineHighMultSPD = inputSHM ? TESTBIT(triggerBits, inputSHM->GetIndexCTP() - 1) : -1;
  fTriggerOnlineHighMultV0 = inputV0M ? TESTBIT(triggerBits, inputV0M->GetIndexCTP() - 1) : -1;
  fTriggerOnlineINT7 = (inputV0C && inputV0A) ?
                       (TESTBIT(triggerBits, inputV0C->GetIndexCTP() - 1) &&
                        TESTBIT(triggerBits, inputV0A->GetIndexCTP() - 1)) : -1;
  
  fTreeEvChar->Fill();
  if(fFillMCGenTrees && fReadMC) ProcessMCGen(mcArray,mcHeader);
  
  //reduce output size by continuing only with selected events to reco TTrees
  if(!isEvSel && fApplyEventSelOnline) return;

  //get PID response
  if(!fPIDresp) fPIDresp = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  
  if(fWriteVariableTreeDs) Process3Prong(array3Prong,aod,mcArray,aod->GetMagneticField(),mcHeader);
  if(fWriteVariableTreeLc2V0bachelor) ProcessCasc(arrayCasc,aod,mcArray,aod->GetMagneticField(),mcHeader);
  
  // Post the data
  PostData(1,fNentries);
  PostData(2,fHistoNormCounter);
  PostData(4,fListCounter);
  PostData(5,fTreeEvChar);
  if(fWriteVariableTreeDs){
    PostData(6,fVariablesTreeDs);
    if(fFillMCGenTrees && fReadMC) PostData(7,fGenTreeDs);
  }
  if(fWriteVariableTreeLc2V0bachelor){
    PostData(8,fVariablesTreeLc2V0bachelor);
    if(fFillMCGenTrees && fReadMC) PostData(9,fGenTreeLc2V0bachelor);
  }
  
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFTreeCreatorApply::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AliAnalysisTaskSEHFTreeCreatorApply: Terminate() \n");
  
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
void AliAnalysisTaskSEHFTreeCreatorApply::Process3Prong(TClonesArray *array3Prong, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  Int_t n3prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of 3prongs: %d\n",n3prong);
  
  Int_t pdgDstoKKpi[3]={321,321,211};
  Int_t nSelectedDs=0;
  Int_t nFilteredDs=0;
  
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

      //PreSelection to speed up task significantly when needed (PbPb `18 or ITSUpgrade)
      TObjArray arrTracks(3);
      if(fITSUpgradePreSelect || fFiltCutsDstoKKpi->GetUsePreselect()){
        for(Int_t ipr=0;ipr<3;ipr++){
          AliAODTrack *tr=vHF->GetProng(aod,ds,ipr);
          arrTracks.AddAt(tr,ipr);
        }
      }
      Int_t preSelectedDs = -1;
      if(fITSUpgradePreSelect){
        preSelectedDs = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 3, 431, pdgDstoKKpi);
        if(preSelectedDs == 0) continue; //Mixture hijing + injected
        if(preSelectedDs == 2) continue; //Only MatchedToMC injected signal
        if(fWriteOnlySignal == 1 && preSelectedDs != 1) continue; //Only matched signal when only signal is requested
      }
      Int_t preSelectedDsCuts = -1;
      if(fFiltCutsDstoKKpi->GetUsePreselect() && ds->GetIsFilled() == 0){
        preSelectedDsCuts = fFiltCutsDstoKKpi->PreSelect(arrTracks);
        if(preSelectedDsCuts==0) continue;
      }

      fNentries->Fill(16);
      nFilteredDs++;
      if((vHF->FillRecoCand(aod,ds))) {////Fill the data members of the candidate only if they are empty.
        
        if(fEnableCandDownsampling){
          Double_t ptCand = ds->Pt();
          Double_t pseudoRand = ptCand * 1000. - (long)(ptCand * 1000);
          if (pseudoRand > fFracToKeepCandDownsampling && ptCand < fMaxPtCandDownsampling) continue;
        }

        Int_t isSelectedFilt=fFiltCutsDstoKKpi->IsSelected(ds,AliRDHFCuts::kAll,aod);
        Int_t isKKpi=isSelectedFilt&1;
        Int_t ispiKK=isSelectedFilt&2;
        Int_t isPhiKKpi=isSelectedFilt&4;
        Int_t isPhipiKK=isSelectedFilt&8;
        Int_t isK0starKKpi=isSelectedFilt&16;
        Int_t isK0starpiKK=isSelectedFilt&32;
        
        if(isSelectedFilt>0){
          
          Double_t modelPred0 = -1.;
          Double_t modelPred1 = -1.;
          AliAODPidHF *Pid_HF = fFiltCutsDstoKKpi->GetPidHF();
          Bool_t isSelectedMLFilt0 = kTRUE;
          Bool_t isSelectedMLFilt1 = kTRUE;
          if(fMLResponse){
            isSelectedMLFilt0 = fMLResponse->IsSelected(modelPred0, ds, aod->GetMagneticField(), Pid_HF, 0);
            isSelectedMLFilt1 = fMLResponse->IsSelected(modelPred1, ds, aod->GetMagneticField(), Pid_HF, 1);
          }
          if(isSelectedMLFilt0 || isSelectedMLFilt1){
            
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
            Bool_t isParticleFromOutOfBunchPileUpEvent = kFALSE;
            //checking origin
            AliAODMCParticle *partDs = 0x0;
            if(fReadMC){
              labDs = ds->MatchToMC(431,arrMC,3,pdgDstoKKpi);
              labDplus = ds->MatchToMC(411,arrMC,3,pdgDstoKKpi);
              
              if(labDs>=0){
                // PILEUP protection for PbPb2018: remove particles from pileup events in efficiency computation
                isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(labDs, mcHeader, arrMC);
                Int_t labDau0=((AliAODTrack*)ds->GetDaughter(0))->GetLabel();
                AliAODMCParticle* p=(AliAODMCParticle*)arrMC->UncheckedAt(TMath::Abs(labDau0));
                pdgCode0=TMath::Abs(p->GetPdgCode());
                partDs = (AliAODMCParticle*)arrMC->At(labDs);
                ptGenDs = partDs->Pt();
              }
              else{
                if(labDplus>=0) {
                  // PILEUP protection for PbPb2018: remove particles from pileup events in efficiency computation
                  isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(labDplus, mcHeader, arrMC);
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
                if(fReadMC && !isParticleFromOutOfBunchPileUpEvent) {
                  if(labDs>=0) {
                    if(orig==4 || orig==5) {
                      if(pdgCode0==321) issignal = kTRUE;
                      else if(pdgCode0==211) isrefl = kTRUE;
                      if(orig==4) isprompt = kTRUE;
                      else if(orig==5) isFD = kTRUE;
                    }
                  }
                  else {
                    isbkg = kTRUE;
                    if(labDplus>=0) fTreeHandlerDs->SetIsDplustoKKpi(kTRUE);//put also D+ -->KKpi in bkg
                  }
                  //do not apply cuts, but enable flag if is selected
                  if(issignal || isbkg || isrefl) fTreeHandlerDs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                }
                if(!fReadMC || (issignal || isbkg || isrefl)) {
                  fTreeHandlerDs->SetIsSelectedStd(isSelAnCutsKKpi,isSelAnTopoCutsKKpi,isSelAnPidCutsKKpi,isSelTracksAnCuts);
                  fTreeHandlerDs->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenDs,modelPred0,ds,bfield,0,fPIDresp,Pid_HF);
                  if(isSelectedMLFilt0) fTreeHandlerDs->FillTree();
                }
              }
              issignal = kFALSE;
              isbkg = kFALSE;
              isprompt = kFALSE;
              isFD = kFALSE;
              isrefl = kFALSE;
              if((fWriteVariableTreeDs==3 && ispiKK) || (fWriteVariableTreeDs==1 && isPhipiKK) || (fWriteVariableTreeDs==2 && isK0starpiKK)) {
                if(fReadMC && !isParticleFromOutOfBunchPileUpEvent) {
                  if(labDs>=0) {
                    if(orig==4 || orig==5) {
                      if(pdgCode0==211) issignal = kTRUE;
                      else if(pdgCode0==321) isrefl = kTRUE;
                      if(orig==4) isprompt = kTRUE;
                      else if(orig==5) isFD = kTRUE;
                    }
                  }
                  else {
                    isbkg = kTRUE;
                    if(labDplus>=0) fTreeHandlerDs->SetIsDplustoKKpi(kTRUE);//put also D+ -->KKpi in bkg
                  }
                  //do not apply cuts, but enable flag if is selected
                  if(issignal || isbkg || isrefl) fTreeHandlerDs->SetCandidateType(issignal,isbkg,isprompt,isFD,isrefl);
                }
                if(!fReadMC || (issignal || isbkg || isrefl)) {
                  fTreeHandlerDs->SetIsSelectedStd(isSelAnCutspiKK,isSelAnTopoCutspiKK,isSelAnPidCutspiKK,isSelTracksAnCuts);
                  fTreeHandlerDs->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenDs,modelPred1,ds,bfield,1,fPIDresp,Pid_HF);
                  if(isSelectedMLFilt1) fTreeHandlerDs->FillTree();
                }
              }
            }//end fill tree
            if(recVtx)fFiltCutsDstoKKpi->CleanOwnPrimaryVtx(ds,aod,origownvtx);
            if(unsetvtx) ds->UnsetOwnPrimaryVtx();
          }//end ML selection
        }//end is selected
      }
      else{
        fNentries->Fill(18); //monitor how often this fails
      }
    }//end Ds
  }//end loop on cadidates
  
  delete vHF;
  return;
}
//________________________________________________________________
void AliAnalysisTaskSEHFTreeCreatorApply::ProcessCasc(TClonesArray *arrayCasc, AliAODEvent *aod, TClonesArray *arrMC, Float_t bfield, AliAODMCHeader *mcHeader){
  
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
      
      //PreSelection to speed up task significantly when needed (PbPb `18 or ITSUpgrade)
      TObjArray arrTracks(2);
      if(fITSUpgradePreSelect || fFiltCutsLc2V0bachelor->GetUsePreselect()){
        for(Int_t ipr=0;ipr<2;ipr++){
          AliAODTrack *tr;
          if(ipr==0) tr = vHF->GetProng(aod,d,ipr);
          else       tr = (AliAODTrack*)(v0part);
          arrTracks.AddAt(tr,ipr);
        }
        if(fFiltCutsLc2V0bachelor->GetUsePreselect() > 1){
          arrTracks.Expand(3);
          arrTracks.AddAt(v0part,2);
        }
      }
      Int_t preSelectedLc = -1;
      if(fITSUpgradePreSelect){
        preSelectedLc = AliVertexingHFUtils::PreSelectITSUpgrade(arrMC, mcHeader, arrTracks, 2, 4122, pdgDgLc2K0Spr);
        if(preSelectedLc == 0) continue; //Mixture hijing + injected
        if(preSelectedLc == 2) continue; //Only MatchedToMC injected signal
        if(fWriteOnlySignal == 1 && preSelectedLc != 1) continue; //Only matched signal when only signal is requested
      }
      Int_t preSelectedLcCuts = -1;
      if(fFiltCutsLc2V0bachelor->GetUsePreselect() && d->GetIsFilled() == 0){
        preSelectedLcCuts = fFiltCutsLc2V0bachelor->PreSelect(arrTracks);
        if(preSelectedLcCuts==0) continue;
      }
      
      fNentries->Fill(33);
      nFilteredLc2V0bachelor++;
      if((vHF->FillRecoCasc(aod,d,kFALSE,fLc2V0bachelorCalcSecoVtx))) {//Fill the data members of the candidate only if they are empty.
        
        //To calculate secondary vertex for pp/pPb if requested
        //Remember to run also CleanUpTask!
        if(d->GetIsFilled()==1 && fLc2V0bachelorCalcSecoVtx) vHF->RecoSecondaryVertexForCascades(aod, d);
        
        if(fEnableCandDownsampling){
          Double_t ptCand = d->Pt();
          Double_t pseudoRand = ptCand * 1000. - (long)(ptCand * 1000);
          if (pseudoRand > fFracToKeepCandDownsampling && ptCand < fMaxPtCandDownsampling) continue;
        }

        Int_t isSelectedFilt = fFiltCutsLc2V0bachelor->IsSelected(d,AliRDHFCuts::kAll,aod); //selected
        
        if(isSelectedFilt > 0){
          
          Double_t modelPred = -1.;
          AliAODPidHF *Pid_HF = fFiltCutsLc2V0bachelor->GetPidHF();
          Bool_t isSelectedMLFilt = kTRUE;
          if(fMLResponse) isSelectedMLFilt = fMLResponse->IsSelected(modelPred, d, aod->GetMagneticField(), Pid_HF, 0);

          Int_t isSelectedCutAn = kFALSE;
          if(fSaveSTDSelection) isSelectedCutAn = fCutsLc2V0bachelor->IsSelected(d,AliRDHFCuts::kAll,aod);

          if(isSelectedMLFilt || (fSaveSTDSelection && isSelectedCutAn > 0)){
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
            Bool_t isParticleFromOutOfBunchPileUpEvent = kFALSE;
            AliAODMCParticle *partLc2V0bachelor=0x0;
            if(fReadMC) {
              labLc2V0bachelor = d->MatchToMC(4122,310,pdgDgLc2K0Spr, pdgDgK0stoDaughters, arrMC, kTRUE);
              if(labLc2V0bachelor>=0){
                // PILEUP protection for PbPb2018: remove particles from pileup events in efficiency computation
                isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(labLc2V0bachelor, mcHeader, arrMC);
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
            
            if(fReadMC && !isParticleFromOutOfBunchPileUpEvent){
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
              fTreeHandlerLc2V0bachelor->SetVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong,ptGenLc2V0bachelor,modelPred,d,bfield,masshypo,fPIDresp,Pid_HF);
              fTreeHandlerLc2V0bachelor->FillTree();
            }
            if(recVtx)fFiltCutsLc2V0bachelor->CleanOwnPrimaryVtx(d,aod,origownvtx);
            if(unsetvtx) d->UnsetOwnPrimaryVtx();
          }//end is selected filt
        } //end ML selection
      } // end Lc IsSelected
      else {
        fNentries->Fill(35); //monitor how often this fails
      }
    }//end Lc2V0bachelor
  }//end loop on candidates
  
  delete vHF;
  return;
}
//_________________________________________________________________
void AliAnalysisTaskSEHFTreeCreatorApply::ProcessMCGen(TClonesArray *arrayMC, AliAODMCHeader *mcHeader){
  /// Fill MC gen trees
  
  for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
    
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
    Int_t absPDG = TMath::Abs(mcPart->GetPdgCode());
    
    if(absPDG == 431 || absPDG == 4122) {
      Bool_t isPrimary = kFALSE;
      Bool_t isFeeddown = kFALSE;

      Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,!fITSUpgradeProduction);//Prompt = 4, FeedDown = 5
      if(orig!=4 && orig!=5) continue; //keep only prompt or feed-down
        
      // PILEUP protection for PbPb2018: remove particles from pileup events in efficiency computation
      Bool_t isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(mcPart->GetLabel(), mcHeader, arrayMC);
      if(isParticleFromOutOfBunchPileUpEvent) continue;

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
      
      if(absPDG == 431 && fWriteVariableTreeDs) {
        deca = AliVertexingHFUtils::CheckDsDecay(arrayMC,mcPart,labDau);
        if(deca!=fWriteVariableTreeDs || labDau[0]<0 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
        fTreeHandlerGenDs->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenDs->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenDs->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        fTreeHandlerGenDs->FillTree();
      } else if(absPDG == 4122 && fWriteVariableTreeLc2V0bachelor) {
        deca = AliVertexingHFUtils::CheckLcV0bachelorDecay(arrayMC,mcPart,labDau);
        if(deca!=1 || labDau[0]<0 || labDau[1]<0) continue;
        isDaugInAcc = CheckDaugAcc(arrayMC,3,labDau,fITSUpgradeProduction);
        fTreeHandlerGenLc2V0bachelor->SetDauInAcceptance(isDaugInAcc);
        fTreeHandlerGenLc2V0bachelor->SetCandidateType(kTRUE,kFALSE,isPrimary,isFeeddown,kFALSE);
        fTreeHandlerGenLc2V0bachelor->SetMCGenVariables(fRunNumber,fEventID,fEventIDExt,fEventIDLong, mcPart);
        fTreeHandlerGenLc2V0bachelor->FillTree();
      }
    }
  }
}
//________________________________________________________________________
std::string AliAnalysisTaskSEHFTreeCreatorApply::GetPeriod(const AliVEvent* event){

  /// Get the period name corresponding to the runnumber that is being analysed

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
//--------------------------------------------------------
Bool_t AliAnalysisTaskSEHFTreeCreatorApply::CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau, Bool_t ITSUpgradeStudy){
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
Bool_t AliAnalysisTaskSEHFTreeCreatorApply::IsCandidateFromHijing(AliAODRecoDecayHF *cand, AliAODMCHeader *mcHeader, TClonesArray* arrMC, AliAODTrack *tr){
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
void AliAnalysisTaskSEHFTreeCreatorApply::SelectGoodTrackForReconstruction(AliAODEvent *aod, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags)
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
AliAODVertex* AliAnalysisTaskSEHFTreeCreatorApply::ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, Double_t bField, Double_t dispersion) {
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
unsigned int AliAnalysisTaskSEHFTreeCreatorApply::GetEvID() {
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
  unsigned int evID = (unsigned int)ev_number + (unsigned int)(fDirNumber<<17);
  fEventNumber++;
  return evID;
}
