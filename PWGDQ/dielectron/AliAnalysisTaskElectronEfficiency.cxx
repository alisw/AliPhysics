/*************************************************************************
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//          Single Electron and Pair-Prefilter Efficiency Task           //
//                                     (description in .h file)          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TChain.h"
#include <math.h>
#include "TObjArray.h"
#include "TList.h"
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliKFParticle.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliVertex.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliDielectron.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronPID.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronHelper.h"

// todo: clean up includes...

#include "AliAnalysisTaskElectronEfficiency.h"
#include <bitset>

ClassImp(AliAnalysisTaskElectronEfficiency)

using std::cout;
using std::endl;


//________________________________________________________________________
AliAnalysisTaskElectronEfficiency::AliAnalysisTaskElectronEfficiency() : 
AliAnalysisTaskSE(),
fESD(0x0),
mcEvent(0x0),
fPIDResponse(0x0),
fPostPIDCntrdCorrTPC(0x0),
fPostPIDWdthCorrTPC(0x0),
fPostPIDCntrdCorrITS(0x0),
fPostPIDWdthCorrITS(0x0),
fUsedVars(0x0),
fSelectPhysics(kFALSE),
fTriggerMask(AliVEvent::kAny),
fEventFilter(0x0),
fRequireVtx(kFALSE),
fCheckV0daughterElectron(kFALSE),
fCutInjectedSignal(kFALSE),
fNminEleInEventForRej(2),
fSupportedCutInstance(0),
//fEventcount(0),
fRandomizeDaughters(kFALSE),
fRandom3(4357),//default seed 4357
fMaxVtxZ(0.),
fCentMin( -1.),
fCentMax(101.),
fEtaMinGEN(-10.),
fEtaMaxGEN( 10.),
fPtMinGEN( -1.),
fPtMaxGEN(100.),
fNptBins(0),
fNetaBins(0),
fNphiBins(0),
fPtBins(0x0),
fEtaBins(0x0),
fPhiBins(0x0),
fsRunBins(""),
fvTrackCuts(),
fvExtraTrackCuts(),
fvDoPrefilterEff(),
fvRejCutMee(),
fvRejCutTheta(),
fvRejCutPhiV(),
fNgen(0x0),
fvReco_Ele(),
fvReco_Ele_poslabel(),
fvAllPionsForRej(),
fvPionsRejByAllSigns(),
fvPionsRejByUnlike(),
fOutputList(0x0),
fOutputListSupportHistos(0x0),
fEventStat(0x0),
tracksT(0x0),
fWriteTree(kFALSE),
pxESD(-1.), // tree variables
pyESD(-1.),
pzESD(-1.),
pTPC(-1.),
chargeT(-999),
signalITS(-1.),
signalTPC(-1.),
beta(-1.),
kchi2ITS(-1.),
kNclsITS(-1),
kITSchi2Cl(-1.),
kNclsTPC(-1),
kTPCchi2Cl(-1.),
kNclsTPCdEdx(-1),
kNFclsTPCr(-1.),
kNFclsTPCfCross(-1.),
kNtrkltsTRD(-1),
kNtrkltsTRDPID(-1),
sigmaEleITS(-999.),
sigmaEleTPC(-999.),
sigmaEleTOF(-999.),
probEleTRD(-999.),
sigmaPioITS(-999.),
sigmaPioTPC(-999.),
sigmaPioTOF(-999.),
sigmaKaoITS(-999.),
sigmaKaoTPC(-999.),
sigmaProITS(-999.),
sigmaProTPC(-999.),
isGlobalT(kFALSE),
isGlobalSDD(kFALSE),
labelT(-1),
pdgT(-1),
labelmotherT(-1),
pdgmotherT(-1),
labelgrandmotherT(-1),
pxMC(-1.),
pyMC(-1.),
pzMC(-1.),
selectedByCut(0),
selectedByExtraCut(0)
{
  /// Default Constructor
}


//________________________________________________________________________
AliAnalysisTaskElectronEfficiency::AliAnalysisTaskElectronEfficiency(const char *name) : 
AliAnalysisTaskSE(name),
fESD(0x0),
mcEvent(0x0),
fPIDResponse(0x0),
fPostPIDCntrdCorrTPC(0x0),
fPostPIDWdthCorrTPC(0x0),
fPostPIDCntrdCorrITS(0x0),
fPostPIDWdthCorrITS(0x0),
fUsedVars(0x0),
fSelectPhysics(kFALSE),
fTriggerMask(AliVEvent::kAny),
fEventFilter(0x0),
fRequireVtx(kFALSE),
fCheckV0daughterElectron(kFALSE),
fCutInjectedSignal(kFALSE),
fNminEleInEventForRej(2),
fSupportedCutInstance(0),
//fEventcount(0),
fRandomizeDaughters(kFALSE),
fRandom3(4357),//default seed 4357
fMaxVtxZ(0.),
fCentMin( -1.),
fCentMax(101.),
fEtaMinGEN(-10.),
fEtaMaxGEN( 10.),
fPtMinGEN( -1.),
fPtMaxGEN(100.),
fNptBins(0),
fNetaBins(0),
fNphiBins(0),
fPtBins(0x0),
fEtaBins(0x0),
fPhiBins(0x0),
fsRunBins(""),
fvTrackCuts(),
fvExtraTrackCuts(),
fvDoPrefilterEff(),
fvRejCutMee(),
fvRejCutTheta(),
fvRejCutPhiV(),
fNgen(0x0),
fvReco_Ele(),
fvReco_Ele_poslabel(),
fvAllPionsForRej(),
fvPionsRejByAllSigns(),
fvPionsRejByUnlike(),
fOutputList(0x0),
fOutputListSupportHistos(0x0),
fEventStat(0x0),
tracksT(0x0),
fWriteTree(kFALSE),
pxESD(-1.), // tree variables
pyESD(-1.),
pzESD(-1.),
pTPC(-1.),
chargeT(-999),
signalITS(-1.),
signalTPC(-1.),
beta(-1.),
kchi2ITS(-1.),
kNclsITS(-1),
kITSchi2Cl(-1.),
kNclsTPC(-1),
kTPCchi2Cl(-1.),
kNclsTPCdEdx(-1),
kNFclsTPCr(-1.),
kNFclsTPCfCross(-1.),
kNtrkltsTRD(-1),
kNtrkltsTRDPID(-1),
sigmaEleITS(-999.),
sigmaEleTPC(-999.),
sigmaEleTOF(-999.),
probEleTRD(-999.),
sigmaPioITS(-999.),
sigmaPioTPC(-999.),
sigmaPioTOF(-999.),
sigmaKaoITS(-999.),
sigmaKaoTPC(-999.),
sigmaProITS(-999.),
sigmaProTPC(-999.),
isGlobalT(kFALSE),
isGlobalSDD(kFALSE),
labelT(-1),
pdgT(-1),
labelmotherT(-1),
pdgmotherT(-1),
labelgrandmotherT(-1),
pxMC(-1.),
pyMC(-1.),
pzMC(-1.),
selectedByCut(0),
selectedByExtraCut(0)
{
  /// Constructor
  
  fvTrackCuts.clear();
  fvExtraTrackCuts.clear();
  fvDoPrefilterEff.clear();
  fvRejCutMee.clear();
  fvRejCutTheta.clear();
  fvRejCutPhiV.clear();
  fvReco_Ele.clear();
  fvReco_Ele_poslabel.clear();
  fvAllPionsForRej.clear();
  fvPionsRejByAllSigns.clear();
  fvPionsRejByUnlike.clear();
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TH1D::Class());
  
  fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kP, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kPIn, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  // for ESD tracks, all variables are filled anyhow...
}


//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::UserCreateOutputObjects()
{
  /// Create histograms
  /// Called once
	AliInfo("Create the output objects");
	Printf("Now running: CreateOutputObjects()");
  
  Printf("  __________________________________________________");
  Printf("  - centrality range:  %f  to  %f", fCentMin, fCentMax);
  Printf("  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
  
  OpenFile(1, "RECREATE");
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner();
  fOutputList->Add(fNgen);
  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) {
    fOutputList->Add(fvReco_Ele.at(iCut));
    fOutputList->Add(fvReco_Ele_poslabel.at(iCut));
    if(fvDoPrefilterEff.at(iCut)){
      fOutputList->Add(fvAllPionsForRej.at(iCut));
      fOutputList->Add(fvPionsRejByAllSigns.at(iCut));
      fOutputList->Add(fvPionsRejByUnlike.at(iCut));
    }
    // be really careful if you need to implement this (see comments in UserExec):
    //    fOutputList->Add(fvReco_Pio.at(iCut));
    //    fOutputList->Add(fvReco_Kao.at(iCut));
    //    fOutputList->Add(fvReco_Pro.at(iCut));
  }
  PostData(1, fOutputList);
  
  OpenFile(2, "RECREATE");
  CreateSupportHistos();
  PostData(2, fOutputListSupportHistos);  
  
  OpenFile(3, "RECREATE");
  //Define Tree for Track Values
  tracksT = new TTree("tracksT","tracksT");
  // branches for both data and MC:
  tracksT->Branch("pxESD",              &pxESD);
  tracksT->Branch("pyESD",              &pyESD);
  tracksT->Branch("pzESD",              &pzESD);
  tracksT->Branch("pTPC",               &pTPC);
  tracksT->Branch("chargeT",            &chargeT);
  tracksT->Branch("signalITS",          &signalITS);
  tracksT->Branch("signalTPC",          &signalTPC);
  tracksT->Branch("beta",               &beta);
  tracksT->Branch("kchi2ITS",           &kchi2ITS);
  tracksT->Branch("kNclsITS",           &kNclsITS);
  tracksT->Branch("kITSchi2Cl",         &kITSchi2Cl); // redundant, for convenience
  tracksT->Branch("kNclsTPC",           &kNclsTPC);
  tracksT->Branch("kTPCchi2Cl",         &kTPCchi2Cl);
  tracksT->Branch("kNclsTPCdEdx",       &kNclsTPCdEdx);
  tracksT->Branch("kNFclsTPCr",         &kNFclsTPCr);
  tracksT->Branch("kNFclsTPCfCross",    &kNFclsTPCfCross);
  tracksT->Branch("kNtrkltsTRD",        &kNtrkltsTRD);
  tracksT->Branch("kNtrkltsTRDPID",     &kNtrkltsTRDPID);
  tracksT->Branch("sigmaEleITS",        &sigmaEleITS);
  tracksT->Branch("sigmaEleTPC",        &sigmaEleTPC);
  tracksT->Branch("sigmaEleTOF",        &sigmaEleTOF);
  //tracksT->Branch("sigmaEleTRD",        &sigmaEleTRD);
  tracksT->Branch("probEleTRD",         &probEleTRD);
  tracksT->Branch("sigmaPioITS",        &sigmaPioITS);
  tracksT->Branch("sigmaPioTPC",        &sigmaPioTPC);
  tracksT->Branch("sigmaPioTOF",        &sigmaPioTOF);
  tracksT->Branch("sigmaKaoITS",        &sigmaKaoITS);
  tracksT->Branch("sigmaKaoTPC",        &sigmaKaoTPC);
  tracksT->Branch("sigmaProITS",        &sigmaProITS);
  tracksT->Branch("sigmaProTPC",        &sigmaProTPC);
  tracksT->Branch("isGlobalT",          &isGlobalT);
  tracksT->Branch("isGlobalSDD",        &isGlobalSDD);
  // additional branches for MC: (not possible to set them inside the "if (fIsMC) {...}")
  tracksT->Branch("labelT",             &labelT);
  tracksT->Branch("pdgT",               &pdgT);
  tracksT->Branch("labelmotherT",       &labelmotherT);
  tracksT->Branch("pdgmotherT",         &pdgmotherT);
  tracksT->Branch("labelgrandmotherT",  &labelgrandmotherT);
  //tracksT->Branch("pMC",                &pMC);
  tracksT->Branch("pxMC",               &pxMC);
  tracksT->Branch("pyMC",               &pyMC);
  tracksT->Branch("pzMC",               &pzMC);
  tracksT->Branch("selectedByCut",      &selectedByCut);
  tracksT->Branch("selectedByExtraCut", &selectedByExtraCut);
  PostData(3, tracksT);
  
  OpenFile(4, "RECREATE");
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",kEventStatBins,0,kEventStatBins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel. (data only)");
    fEventStat->GetXaxis()->SetBinLabel(3,"After event filter");
  }
  PostData(4, fEventStat);
}


//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::UserExec(Option_t *)
{
  /// Strategy:  Process one cutInstance (or multiple ones) for analysis tracking&PID efficiency (as usual).
  ///            Process optional, separate cutInstance for prefilter efficiencies: it also produces the usual tracking&PID efficiency
  ///            (but of course for the specified prefilter track sample, so mainly for convenience and curiosity),
  ///            and then computes the pair rejection efficiency, using random rejection of "testparticles" with the selected electrons. (further info in 'CalcPrefilterEff()')
  
  AliESDInputHandler *inputHandlerESD = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandlerESD) { Printf("ERROR: Could not get ESDInputHandler\n"); }
  else fESD = inputHandlerESD->GetEvent();
  if (!fESD) { Printf("ERROR: fESD not available"); return; }

  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if (mcH) mcEvent = mcH->MCEvent();
  if (!mcEvent) { Printf("ERROR: mcEvent not available"); return; }
  if (!fPIDResponse) SetPIDResponse( ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse() );
  AliDielectronVarManager::SetPIDResponse(fPIDResponse);
  AliDielectronMC::Instance()->ConnectMCEvent();
  // set pid correction function to var manager
  if(fPostPIDCntrdCorrTPC) AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorrTPC);
  if(fPostPIDWdthCorrTPC)  AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorrTPC);
  if(fPostPIDCntrdCorrITS) AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorrITS);
  if(fPostPIDWdthCorrITS)  AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDWdthCorrITS);
  
  AliStack *fStack = mcEvent->Stack();
  
  Bool_t isESD=kTRUE;
  Bool_t isAOD=kFALSE; //currently not supported!
  // taken from AliAnalysisTaskMultiDielectron::UserExec():
  // -----
  // Was event selected ?
  ULong64_t isSelected = AliVEvent::kAny;
  //Bool_t isRejected = kFALSE;
  if( fSelectPhysics && inputHandlerESD){
    if((isESD && inputHandlerESD->GetEventSelection()) || isAOD){
      isSelected = inputHandlerESD->IsEventSelected();
      //if (fExcludeTriggerMask && (isSelected&fExcludeTriggerMask)) isRejected=kTRUE;
      //if (fTriggerLogic==kAny) isSelected&=fTriggerMask;
      //else if (fTriggerLogic==kExact) isSelected=((isSelected&fTriggerMask)==fTriggerMask);
      //cout << " isSelected = " << (bitset<32>)isSelected << " \t fTriggerMask = " << (bitset<32>)fTriggerMask << endl;
      isSelected&=fTriggerMask;
      //TString firedTriggerClasses=InputEvent()->GetFiredTriggerClasses();
      //if(!fFiredTrigger.IsNull()) isSelected=(firedTriggerClasses.Contains(fFiredTrigger))^fFiredExclude;
    }
  }
  // -----
  // Before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected==0) {
    PostData(4,fEventStat);
    return;
  }
  // After physics selection
  fEventStat->Fill(kPhysicsSelectionEvents);
  
  // Event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(fESD)) return; // fESD instead of 'InputEvent()'
  }
  // centrality
  Double_t centralityF=-1;
  AliCentrality *esdCentrality = fESD->GetCentrality();
  if (esdCentrality) centralityF = esdCentrality->GetCentralityPercentile("V0M");
  if (centralityF<fCentMin || centralityF>=fCentMax) return;
  // vertex, just for monitoring
  const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
  Double_t vtxZGlobal = -99.;
  Int_t nCtrb = -1;
  if (vtxESD) {
    vtxZGlobal = vtxESD->GetZ(); // was GetZv(); until Jan 2015
    nCtrb = vtxESD->GetNContributors();
  }
  
  // After event filter
  fEventStat->Fill(kFilteredEvents);
  PostData(4,fEventStat);
  
  
//  ++fEventcount;
//  Printf("__________ next Event selected ( %i ) __________", fEventcount);
  Int_t kRunNumber = fESD->GetRunNumber();
  
  // for selected events, fill some histos:
  (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(0)))->Fill(0.);
  (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(1)))->Fill(centralityF);
  (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(4)))->Fill(vtxZGlobal);//hVertexZ
  (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(5)))->Fill(nCtrb);//hNvertexCtrb
  //hNTrksEvent_cent->Fill(fESD->GetNumberOfTracks(),centralityF);
  
  Int_t Nacc = AliDielectronHelper::GetNacc(fESD);
  Double_t sigmaEleITS_PIDresp;
  Double_t sigmaEleTPC_PIDresp;
  Double_t sigmaEleITS_Raw;
  Double_t sigmaEleTPC_Raw;
  
  //
  // store if there is at least one cutset to be used for Prefilter efficiency determination.
  // store for each cutset if there is at least one electron selected by 'fvExtraTrackCuts'.
  //
  Bool_t              atleastOnePrefilterSetting=kFALSE;
  Bool_t              atleastOneEleExtraSelected=kFALSE;
  std::vector<Bool_t> fvAtleastOneEleExtra_perCut;
  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) {
    if (fvDoPrefilterEff.at(iCut) == kTRUE) atleastOnePrefilterSetting=kTRUE;
    fvAtleastOneEleExtra_perCut.push_back(kFALSE); // initialize vector contents
  }
  
  //
  // store all accepted track IDs in vectors to use them for the Prefilter efficiency determination.
  //
  std::vector<Int_t> vecTrkID;
  vecTrkID.push_back(0); // first element is neccessary to hook track vector into cut vector.
  std::vector< std::vector<Int_t> > vecEleCand_perCut; // one vector per cut setting
  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) { vecEleCand_perCut.push_back(vecTrkID); }
  
  Int_t NEleSelected = 0;
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) 
  {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) { Printf("ERROR: Could not receive track %d", iTracks); continue; }
    
    Double_t trackPt(-1.);
    trackPt = track->Pt();
    
    if (fCheckV0daughterElectron) {
      //Bool_t isV0daughterElectron = kFALSE;
      //isV0daughterElectron = track->TestBit(BIT(14));
      if (track->TestBit(BIT(14))) continue;
    }
    Int_t label = track->GetLabel();
    Int_t abslabel = TMath::Abs( track->GetLabel() );
    //if(label < 0) continue; // two sets of histograms will be filled to check difference.
    AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack( abslabel ));
    if (!mctrack) continue;
    //
    Int_t trackpdg = (Int_t) TMath::Abs(mctrack->PdgCode());
    if (!( trackpdg==Int_t(11) ) ) continue; // do not try to include hadrons, see below.
    //Printf(" found electron (iTracks=%i)", iTracks);
    // _______________
    // the following inclusion of some hadrons should NOT be used!
    // for some reason a few of the true electrons will get lost at some point and not filled into histograms!
    // tested at GSI with ./run user -o -l MC_PbPb/2.76ATeV/LHC12a17h_fix/169099.ana.txt -n 20 -f -c MC_PbPb
    //if (!( trackpdg==Int_t(11) || trackpdg==Int_t(211) || trackpdg==Int_t(321) || trackpdg==Int_t(2212) ) ) {
    //  continue; //Printf(" skip this track (trackpdg=%i) (iTracks=%i)", trackpdg, iTracks);
    //}
    
    Int_t motherlabel = mctrack->Particle()->GetFirstMother();
    if (motherlabel<0) continue; // motherlabel 0 already gives the first valid particle.
    AliMCParticle *mother = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherlabel));
    if (!mother) continue;
    if (!mother->Particle()->IsPrimary()) continue;
    // from Patrick:
    // for removal of injected signals.
    if (fCutInjectedSignal && IsInjectedSignal(mcEvent, abslabel)) {
      //Printf(" kicked out injected signal (trackpdg=%i) (iTracks=%i)", trackpdg, iTracks);
      continue;
    }
    
    Double_t mcPt(-1.),mcEta(-9.),mcPhi(-9.);
    mcPt = mctrack->Pt(); mcEta = mctrack->Eta(); mcPhi = mctrack->Phi();
    
    selectedByCut = 0;
    selectedByExtraCut = 0;
    
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) // loop over all specified cutInstances
    {
      //cutting logic taken from AliDielectron::FillTrackArrays()
      UInt_t selectedMask=(1<<fvTrackCuts.at(iCut)->GetCuts()->GetEntries())-1;
      //apply track cuts
      UInt_t cutmask=fvTrackCuts.at(iCut)->IsSelected(track);
      //cout << "   cutmask = " << cutmask << "   selectedMask = " << selectedMask << endl;
      if (cutmask!=selectedMask) continue;
      //cout << "   PDG code " << trackpdg << endl;
      //cout << "   cutmask = " << cutmask << "   selectedMask = " << selectedMask << endl;
      
      // for later pairing to get random rejection probability
      //cout << "adding ele for cut " << iCut << ": " << endl;
      vecEleCand_perCut.at(iCut).push_back(iTracks);
      //cout << "new ele in array: " << vecEleCand_perCut.at(iCut).at( vecEleCand_perCut.at(iCut).size()-1 ) << endl;
      
//      if (trackpdg==Int_t(11)) { //Electron
        fvReco_Ele.at(iCut)->Fill(mcPt,mcEta,mcPhi);
        if(label > 0) {
          fvReco_Ele_poslabel.at(iCut)->Fill(mcPt,mcEta,mcPhi);
        }
//      }
      // be really careful if you need to implement this (see comments in UserExec):
//      else if (trackpdg==Int_t(211)) { //PiPlus
//        fvReco_Pio.at(i)->Fill(mcPt,mcEta,mcPhi);
//      }
//      else if (trackpdg==Int_t(321)) { //KPlus
//        fvReco_Kao.at(i)->Fill(mcPt,mcEta,mcPhi);
//      }
//      else if (trackpdg==Int_t(2212)) { //Proton
//        fvReco_Pro.at(i)->Fill(mcPt,mcEta,mcPhi);
//      }
//      else continue;
      
      selectedByCut|=1<<(iCut); // store bitwise which cut settings the track survived.
      
      // store infos related to prefilter efficiency
      if (fvDoPrefilterEff.at(iCut) == kTRUE) {
        // check if one of the prefilter electrons also survives the analysis cuts of this cutset.
        //cout << "fvExtraTrackCuts.at("<<iCut<<") --- .size() = "<<fvExtraTrackCuts.size() << endl;
        UInt_t selectedMaskExtra = (1<<fvExtraTrackCuts.at(iCut)->GetCuts()->GetEntries())-1;
        UInt_t cutmaskExtra      = fvExtraTrackCuts.at(iCut)->IsSelected(track);
        if (cutmaskExtra==selectedMaskExtra) {
          selectedByExtraCut|=1<<(iCut); // this is just for the tree.
          atleastOneEleExtraSelected = kTRUE;
          fvAtleastOneEleExtra_perCut.at(iCut) = kTRUE;
        }
      }
      
    } //loop over cutInstances
    //cout << "selectedByCut = " << selectedByCut << endl;
    
    //________________________________________________________________
    if (selectedByCut==0) continue;
    // only go on if the track survived at least 1 of the cut settings!
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    // get track information
    // feel free to add more information to the tree, some more variables are already defined as branches...
    pxESD = track->Px();
    pyESD = track->Py();
    pzESD = track->Pz();
    chargeT   = track->Charge();
    signalITS = track->GetITSsignal();
    signalTPC = track->GetTPCsignal();
    //
    kNclsTPC  = track->GetTPCNcls();
    kTPCchi2Cl  = track->GetTPCchi2() / kNclsTPC; // only for ESD! // AOD like this: particle->Chi2perNDF()*(tpcNcls-5)/tpcNcls
    kNclsTPCdEdx = track->GetTPCsignalN(); // ("fTPCsignalN") number of points used for dEdx - maybe not in AOD...?
    kNFclsTPCr = track->GetTPCClusterInfo(2,1); // number of findable clusters(crossed rows) in the TPC with more robust definition //ESD & AOD!
    Float_t kNFclsTPC = track->GetTPCNclsF(); //tpcClusFindable // number of findable clusters in the TPC //ESD & AOD!
    kNFclsTPCfCross = (kNFclsTPC>0)?(kNFclsTPCr/kNFclsTPC):0; // fraction crossed rows/findable clusters in the TPC, as done in AliESDtrackCuts  //ESD & AOD!
    //kNtrkltsTRD = 0; // only exists for ESD, see below
    kNtrkltsTRDPID = track->GetTRDntrackletsPID();
    //
    // fill the AliDielectronVarManager to get some specific variables, e.g. which were used for track selection.
    Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
    AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.
    AliDielectronVarManager::Fill(track, values);
    //
    sigmaEleITS_Raw = values[AliDielectronVarManager::kITSnSigmaEleRaw];
    sigmaEleTPC_Raw = values[AliDielectronVarManager::kTPCnSigmaEleRaw];
    sigmaEleITS = values[AliDielectronVarManager::kITSnSigmaEle];
    sigmaEleTPC = values[AliDielectronVarManager::kTPCnSigmaEle];
    sigmaEleTOF = values[AliDielectronVarManager::kTOFnSigmaEle];
    //
    // temporary, as cross-check:
    sigmaEleITS_PIDresp = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
    sigmaEleTPC_PIDresp = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
    //sigmaEleTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
    
    // TODO: which momentum is better?
    //pTPC        = values[AliDielectronVarManager::kPIn]; // uses GetInnerParam() //Track parameters estimated at the inner wall of TPC
    pTPC        = ((AliESDtrack*)track)->GetTPCInnerParam()->P(); //Track parameters estimated at the inner wall of TPC using the TPC stand-alone
    kchi2ITS    = ((AliESDtrack*)track)->GetITSchi2();
    kNclsITS    = ((AliESDtrack*)track)->GetNcls(0);
    kITSchi2Cl  = kchi2ITS / kNclsITS;
    
    // get MC track information
    pdgT          = mctrack->PdgCode();
    labelmotherT  = motherlabel;
    pdgmotherT    = mother->PdgCode();
    
    
    if (fWriteTree) {
      // variable 'selectedByCut' is stored in the Tree to distinguish the cut settings...
      tracksT->Fill(); //Fill Track Tree
    } // if fWriteTree
    
    
    if (selectedByCut & 1<<fSupportedCutInstance) // only go on if the track survived in the cut setting for which you want to fill the support histograms.
    { 
      NEleSelected++;
      Float_t pESD = track->P();
      
      // for selected tracks, fill some histos:
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(6)))->Fill(values[AliDielectronVarManager::kPIn], pTPC);//hPInVarMgr_PInStandAlone
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(7)))->Fill(track->Pt());//hPt (reco)
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(8)))->Fill(pESD, pTPC);//hP_PIn
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(9)))->Fill(pESD, mctrack->P());//hP_Pgen
      // ITS
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(10)))->Fill(pESD, signalITS);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(11)))->Fill(pESD, sigmaEleITS);
      //(dynamic_cast<TH2F *>(fOutputListSupportHistos->At(12)))->Fill(pESD, );
      //(dynamic_cast<TH2F *>(fOutputListSupportHistos->At(13)))->Fill(pESD, );
      // TPC
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(14)))->Fill(pTPC, signalTPC);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(15)))->Fill(pTPC, sigmaEleTPC);
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(16)))->Fill(pTPC, sigmaEleTPC, signalTPC);
      // old way, directly from fPIDResponse object.
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(20)))->Fill(pESD, sigmaEleITS_PIDresp);//hITSnSigmaElePIDresp_P
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(21)))->Fill(pTPC, sigmaEleTPC_PIDresp);//hTPCnSigmaElePIDresp_P
      // run dependency
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(22)))->Fill(pTPC, signalTPC, kRunNumber);//hTPC_dEdx_P_run
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(23)))->Fill(pTPC, sigmaEleTPC, kRunNumber);//hTPCnSigmaEle_P_run
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(24)))->Fill(pESD, sigmaEleITS, kRunNumber);
      // TOF
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(25)))->Fill(pTPC, sigmaEleTOF);
      //(dynamic_cast<TH2F *>(fOutputListSupportHistos->At(26)))->Fill(pTPC, );
      // Eta and Phi
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(27)))->Fill(mcEta);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(28)))->Fill(mcPhi);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(29)))->Fill(mcEta, mcPhi);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(30)))->Fill(mcEta, signalTPC);
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(31)))->Fill(mcEta, signalTPC, pTPC);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(32)))->Fill(mcEta, sigmaEleTPC);
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(33)))->Fill(mcEta, sigmaEleTPC, pTPC);
      //(dynamic_cast<TH3F *>(fOutputListSupportHistos->At(34)))->Fill(mcEta, sigmaEleTPC, );//hTPCnSigmaEle_Eta_RefMultTPConly
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(35)))->Fill(mcEta, sigmaEleTPC, Nacc);//hTPCnSigmaEle_Eta_Nacc
      // DCA
      //(dynamic_cast<TH1F *>(fOutputListSupportHistos->At(36)))->Fill();
      //(dynamic_cast<TH1F *>(fOutputListSupportHistos->At(37)))->Fill();
      //(dynamic_cast<TH2F *>(fOutputListSupportHistos->At(38)))->Fill();
      // Quality
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(40)))->Fill(kNFclsTPCfCross);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(41)))->Fill(kNFclsTPCr);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(42)))->Fill(kNclsTPC);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(43)))->Fill(kNclsITS);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(44)))->Fill(kTPCchi2Cl);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(45)))->Fill(kITSchi2Cl);
      //(dynamic_cast<TH1F *>(fOutputListSupportHistos->At(46)))->Fill(kNclsSFracTPC);
      //(dynamic_cast<TH1F *>(fOutputListSupportHistos->At(47)))->Fill(kTPCclsDiff);
      (dynamic_cast<TH1F *>(fOutputListSupportHistos->At(48)))->Fill(kNclsTPCdEdx);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(49)))->Fill(kNclsTPC, kNFclsTPCr);
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(50)))->Fill(track->Pt(), kNFclsTPCr);
      
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(60)))->Fill(mcEta, sigmaEleTPC);
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(61)))->Fill(mcEta, sigmaEleTPC, pTPC);
      //(dynamic_cast<TH3F *>(fOutputListSupportHistos->At(62)))->Fill(mcEta, sigmaEleTPC, );//hITSnSigmaEle_Eta_RefMultTPConly
      (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(63)))->Fill(mcEta, sigmaEleTPC, Nacc);//hITSnSigmaEle_Eta_Nacc
      
      // check of post PID correction. raw values from AliDielectronVarManager
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(64)))->Fill(pESD, sigmaEleITS_Raw);//hITSnSigmaEleRaw_P
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(65)))->Fill(pTPC, sigmaEleTPC_Raw);//hTPCnSigmaEleRaw_P
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(66)))->Fill(mcEta, sigmaEleITS_Raw);//hITSnSigmaEleRaw_Eta
      (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(67)))->Fill(mcEta, sigmaEleTPC_Raw);//hTPCnSigmaEleRaw_Eta
      
    } //fSupportedCutInstance
    
  } //track loop
  
  //(dynamic_cast<TH2F *>(fOutputListSupportHistos->At(2)))->Fill(centralityF,0);//hNTrksEvent_cent
  (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(3)))->Fill(centralityF,NEleSelected);//hNEleEvent_cent
  
//  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) {
//    cout << " iCut=" << iCut << ", ele candidate IDs: " << flush;
//    for (UInt_t iele=1; iele<vecEleCand_perCut.at(iCut).size(); iele++) { // iele=0 is the default entry which is filled evertime, but not used (see init of 'vecEleCand_perCut').
//      cout << "  " << vecEleCand_perCut.at(iCut).at(iele) << flush;
//    }
//  } //cout << " " << endl;
  
  //
  // call the Prefilter efficiency calculation
  // please see comments above the function declaration for further info.
  //
  if (atleastOnePrefilterSetting && atleastOneEleExtraSelected) CalcPrefilterEff(mcEvent, vecEleCand_perCut, fvAtleastOneEleExtra_perCut);
  
  
  if(mcEvent){
    Int_t nMCtracks = mcEvent->GetNumberOfTracks();
    //AliStack *fStack = mcEvent->Stack();
    
    for(Int_t iMCtrack = 0; iMCtrack < nMCtracks; iMCtrack++)
    {
      //if(!fStack->IsPhysicalPrimary(iMCtrack)) continue;
      AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iMCtrack));
      if( TMath::Abs(mctrack->PdgCode()) !=  11 ) continue;
      
      // from Theo:
      Double_t mcPt(-1.),mcEta(-9.),mcPhi(-9.);
      mcPt = mctrack->Pt(); mcEta = mctrack->Eta(); mcPhi = mctrack->Phi();
      if(mcPt  < fPtMinGEN  || mcPt  > fPtMaxGEN)  continue;
      if(mcEta < fEtaMinGEN || mcEta > fEtaMaxGEN) continue;
      
      // check which kinematically ok tracks are skipped by "IsPhysicalPrimary".
      // checking here will already be clean of low-pt e- from scattering and bremsstrahlung etc.
      if(!fStack->IsPhysicalPrimary(iMCtrack)) { //Printf("track not ->IsPhysicalPrimary()"); 
        //        TParticle* mothertemp      = fStack->Particle(mctrack->GetMother());
        //        if (mothertemp->GetPdgCode() != 22) { // just to reduce the amount of printouts
        //          Printf("  track: %s \t mother: %s \t mother->IsPrimary(): %d", GetParticleName(mctrack->PdgCode()), GetParticleName(mothertemp->GetPdgCode()), mothertemp->IsPrimary());
        //        }
        continue; // excludes huge number of e from gamma. also from: K+- [some of them are "IsPrimary()"], K_L0, pi0 (not prim). (nothing else seen in LHC12a17h_fix/169099.ana.txt -n 5)
      }
      
      // from Theo:
      Int_t motherlabel = mctrack->Particle()->GetFirstMother();
      if(!motherlabel) { Printf("bad motherlabel"); continue; }
      AliMCParticle *mother = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherlabel));
      if(!mother) { Printf("bad mother"); continue; }
      if(!mother->Particle()->IsPrimary()) { //Printf("mother not ->IsPrimary()"); 
        //Printf("  mother: %s \t fromBGEvent: %d \t isPhysPrim: %d", GetParticleName(mother->PdgCode()), mcEvent->IsFromBGEvent(motherlabel), mcEvent->IsPhysicalPrimary(motherlabel));
        continue; // kicks out e from gamma, pi0, K_L0, K+, D-, D+, (...?)
        // when checking "fStack->IsPhysicalPrimary(i)" before, then this kicks out e from: 
        // pi0, D+-, D0, D0_bar, B-, B0_bar, eta, tau-, D_s+ (all "fromBGEvent: 1 & isPhysPrim: 0");
        // J/psi, Xi_c0, eta, D0 (all "fromBGEvent: 0 & isPhysPrim: 0")
      }
      
      // from Patrick:
      // for removal of injected signals.
      if (fCutInjectedSignal && IsInjectedSignal(mcEvent, iMCtrack)) continue;
      // this still has a strong effect after "fStack->IsPhysicalPrimary(i)" and "mother->Particle()->IsPrimary()". skips e from J/psi
      
      fNgen->Fill(mcPt,mcEta,mcPhi);
    } //track loop
  } //MC loop
  //Printf("__________ end of Event ( %i ) __________", fEventcount);
}


//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CalcPrefilterEff(AliMCEvent* mcEventLocal, const std::vector< std::vector<Int_t> > & vvEleCand, const std::vector<Bool_t> & vbEleExtra) 
{
  /// Determination of random electron rejection efficiency due to pair-prefiltering (used for photon conversion + Dalitz rejection).
  /// It is estimated by pairing a sample of "testparticles" (using primary, non-injected, charged pions) with the electrons of the prefilter sample (stored in 'vvEleCand').
  /// With this method the pair has no real correlation, just like in the case of random rejection of a signal electron.
  /// The pairing is only done for cut settings for which in the current event at least one final analysis electron was found (stored in 'vbEleExtra').
  /// Probably a second electron of any type should be in the event, otherwise the chance of rejecting an electron-electron pair is zero. Controlled by variable 'fNminEleInEventForRej' (=2 by default).
  //  @TODO: question: how to treat ULS and LS pairing in this context?
  /// The prefilter pair cuts are applied to these random pairs. All and rejected testparticles are stored in 3D histograms.
  /// The advantage of this method is a huge number of random pairs, the disadvantage may be a systematic difference in the implicitly available number of pairs...
  /// (by Patrick)
  
  Int_t nMCtracks  = mcEventLocal->GetNumberOfTracks();
  AliStack *fStack = mcEventLocal->Stack();
  Double_t eleMass = AliPID::ParticleMass(AliPID::kElectron);
  
  for(Int_t iMCtrack = 0; iMCtrack < nMCtracks; iMCtrack++)
  {
    // select only primary, non-injected, charged pions:
    if(!fStack->IsPhysicalPrimary(iMCtrack)) continue;
    AliMCParticle *mcPion = dynamic_cast<AliMCParticle *>(mcEventLocal->GetTrack(iMCtrack));
    if( TMath::Abs(mcPion->PdgCode()) !=  211 ) continue;
    if (fCutInjectedSignal && IsInjectedSignal(mcEventLocal, iMCtrack)) continue;
    
    Double_t mcPt(-1.),mcEta(-9.),mcPhi(-9.);
    mcPt = mcPion->Pt(); mcEta = mcPion->Eta(); mcPhi = mcPion->Phi();
    if(mcPt  < fPtMinGEN  || mcPt  > fPtMaxGEN)  continue;
    if(mcEta < fEtaMinGEN || mcEta > fEtaMaxGEN) continue;
    
    // need TLorentzVectors for the pairing
		Double_t vMomPi[3]; mcPion->PxPyPz(vMomPi);
    Double_t momPi    = mcPion->P();
		Double_t energyPi = TMath::Sqrt(momPi*momPi + eleMass*eleMass);	// need to use hardcoded electron mass.
    TLorentzVector dau1;
    dau1.SetPxPyPzE(vMomPi[0],vMomPi[1],vMomPi[2],energyPi);
    Int_t chargePion = mcPion->Charge();
    
    // loop over all cut settings
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
    {
      if (vbEleExtra.at(iCut) == kFALSE) continue; // if there is no global electron in the event, no prefilter pairing should be done.
      if (vvEleCand.at(iCut).size()-1 < fNminEleInEventForRej) continue; // size-1 because vec[0] is not used (always filled with '0').
      if (fvDoPrefilterEff.at(iCut) == kFALSE) continue; // may not even be needed anymore...
      // flags for the pion
      Bool_t rejByAllSigns = kFALSE;
      Bool_t rejByUnlike   = kFALSE;
      
      // pair the testparticle with all selected electron candidates and check for close pairs. 
      // pairing starts at vec[1] because vec[0] is always set to 0, since this is needed to initialize the 2D vector (see UserExec).
      for (UInt_t iEle=1; iEle<vvEleCand.at(iCut).size(); iEle++) 
      {
        AliMCParticle *mcEle = dynamic_cast<AliMCParticle *>(mcEventLocal->GetTrack(vvEleCand.at(iCut).at(iEle)));
        Double_t vMomEle[3]; mcEle->PxPyPz(vMomEle);
        Double_t momEle    = mcEle->P();
        Double_t energyEle = TMath::Sqrt(momEle*momEle + eleMass*eleMass);	// need to use hardcoded electron mass.
        TLorentzVector dau2;
        dau2.SetPxPyPzE(vMomEle[0],vMomEle[1],vMomEle[2],energyEle);
        Int_t chargeEle = mcEle->Charge();
        
        // Invariant Mass:
        Double_t invmass = (dau1+dau2).M();
        // Pair Pt:
        Double_t pairpt  = (dau1+dau2).Pt();
        // Opening Angle:
        Double_t theta   = dau1.Angle(dau2.Vect());
        // PhiV:
        Double_t phiv    = PhivPair(fESD->GetMagneticField(), chargePion, chargeEle, dau1.Vect(), dau2.Vect());
        //cout << " invmass = " << invmass << "  theta = " << theta << "  phiv = " << phiv << endl;
        
        (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(51)))->Fill(invmass,pairpt);//hMeePtee
        (dynamic_cast<TH3F *>(fOutputListSupportHistos->At(52)))->Fill(invmass,phiv,theta);//hMeePhiVOpen
        (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(53)))->Fill(invmass,theta);//hMeeOpen
        (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(54)))->Fill(invmass,phiv);//hMeePhiV
        (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(55)))->Fill(pairpt, theta);//hPteeOpen
        (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(56)))->Fill(pairpt, phiv);//hPteePhiV
        (dynamic_cast<TH2F *>(fOutputListSupportHistos->At(57)))->Fill(theta,  phiv);//hOpenPhiV
        
        // tag the pion as rejected, if the pair falls within following cut:
        if ( (invmass < fvRejCutMee.at(iCut)   || fvRejCutMee.at(iCut)   < 0.) // -> within cut or cut disabled
            && (theta < fvRejCutTheta.at(iCut) || fvRejCutTheta.at(iCut) < 0.)
            && (phiv  > fvRejCutPhiV.at(iCut)  || fvRejCutPhiV.at(iCut)  > 3.14) )
        {
          //cout << " pion rejected! invmass = " << invmass << "  theta = " << theta << "  phiv = " << phiv << endl;
          rejByAllSigns = kTRUE;
          if (chargePion!=chargeEle) rejByUnlike = kTRUE;
        }
      } //electron loop
      
      //fill histograms per cut setting to determine rejection efficiency
      fvAllPionsForRej.at(iCut)->Fill(mcPt,mcEta,mcPhi); //reference histogram (denominator)
      if (rejByAllSigns) fvPionsRejByAllSigns.at(iCut)->Fill(mcPt,mcEta,mcPhi);
      if (rejByUnlike)   fvPionsRejByUnlike.at(iCut)->Fill(mcPt,mcEta,mcPhi);
    } //cut loop
    
  } //pion loop
}


//______________________________________________
Double_t AliAnalysisTaskElectronEfficiency::PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2) //const
{
  /// Following the idea to use opening of collinear pairs in magnetic field from e.g. PHENIX
  /// to identify conversions. Angle between ee plane and magnetic field is calculated (0 to pi).
  /// Due to tracking to the primary vertex, conversions with no intrinsic opening angle 
  /// always end up as pair in "cowboy" configuration. The function as defined here then 
  /// returns values close to pi.
  /// Correlated Like Sign pairs (from double conversion / dalitz + conversion) may show up 
  /// at pi or at 0 depending on which leg has the higher momentum. (not checked yet)
  /// This expected ambiguity is not seen due to sorting of track arrays in this framework. 
  /// To reach the same result as for ULS (~pi), the legs are flipped for LS.
  /// from PWGDQ/dielectron/AliDielectronPair.cxx

  //Define local buffer variables for leg properties
  Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
  Double_t px2=-9999.,py2=-9999.,pz2=-9999.;

  TVector3 fD1=dau1;
  TVector3 fD2=dau2;
  Int_t    d1Q=charge1;
  //Int_t    d2Q=charge2;
  if (fRandomizeDaughters) { // randomize daughters if requested
    if (fRandom3.Rndm()>0.5) {
      fD1=dau2;
      fD2=dau1;
      d1Q=charge2;
      //d2Q=charge1;
    }
  }
  else { // sort particles according to pt, as done by default in AliDielectronPair
    if (dau1.Pt() < dau2.Pt()) {
      fD1=dau2;
      fD2=dau1;
      d1Q=charge2;
      //d2Q=charge1;
    }
  }
  
  if (charge1*charge2 > 0.) { // Like Sign
    if(MagField<0){ // inverted behaviour
      if(d1Q>0){
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }
    }
  }
  else { // Unlike Sign
    if(MagField>0){ // regular behaviour
      if(d1Q>0){
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();
        
        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();
        
        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();
        
        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();
        
        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }
    }
  }

  Double_t px = px1+px2;
  Double_t py = py1+py2;
  Double_t pz = pz1+pz2;
  Double_t dppair = TMath::Sqrt(px*px+py*py+pz*pz);
  
  //unit vector of (pep+pem) 
  Double_t pl = dppair;
  Double_t ux = px/pl;
  Double_t uy = py/pl;
  Double_t uz = pz/pl;
  Double_t ax = uy/TMath::Sqrt(ux*ux+uy*uy);
  Double_t ay = -ux/TMath::Sqrt(ux*ux+uy*uy); 
  
  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  //Double_t ptep = iep->Px()*ax + iep->Py()*ay; 
  //Double_t ptem = iem->Px()*ax + iem->Py()*ay; 
  
  Double_t pxep = px1;
  Double_t pyep = py1;
  Double_t pzep = pz1;
  Double_t pxem = px2;
  Double_t pyem = py2;
  Double_t pzem = pz2;
  
  //vector product of pep X pem 
  Double_t vpx = pyep*pzem - pzep*pyem; 
  Double_t vpy = pzep*pxem - pxep*pzem; 
  Double_t vpz = pxep*pyem - pyep*pxem; 
  Double_t vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz); 
  //Double_t thev = acos(vpz/vp); 
  
  //unit vector of pep X pem 
  Double_t vx = vpx/vp; 
  Double_t vy = vpy/vp; 
  Double_t vz = vpz/vp; 
  
  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz) 
  Double_t wx = uy*vz - uz*vy; 
  Double_t wy = uz*vx - ux*vz; 
  //Double_t wz = ux*vy - uy*vx; 
  //Double_t wl = sqrt(wx*wx+wy*wy+wz*wz); 
  // by construction, (wx,wy,wz) must be a unit vector. 
  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them 
  // should be small if the pair is conversion.
  // this function then returns values close to pi!
  Double_t cosPhiV = wx*ax + wy*ay; 
  Double_t phiv = TMath::ACos(cosPhiV); 
  
  return phiv;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskElectronEfficiency::IsInjectedSignal(AliMCEvent* mcEventLocal, Int_t tracklabel) 
{
  /// Function to identify injected signals. (by Patrick)
  
  Bool_t isInjected  = kFALSE;
  Bool_t fromBGEvent = 0;
  Bool_t isPhysPrim  = 0;
  TParticle* mother  = NULL;
  AliStack *fStack   = mcEventLocal->Stack();
  AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEventLocal->GetTrack(tracklabel));
  
  fromBGEvent = mcEventLocal->IsFromBGEvent(tracklabel);
  isPhysPrim  = mcEventLocal->IsPhysicalPrimary(tracklabel);
  // fromBGEvent:  particle comes from the generated MC event (?). true for most pi, K, p, some e.
  // !fromBGEvent: particle could be injected or produced in Geant (from weak decay / photon conv).
  //               weak decays and photon conv are handled by Geant -> will NOT be fromBGEvent!
  // isPhysPrim: in a chain of decaying particles, denotes the first which is stable under strong and EM force.
  // - injected: the physical primary will NOT be fromBGEvent! (<= that is the relevant criterion!)
  // - Geant:    the physical primary will be fromBGEvent!
  // 
  // example: for injected J/psi -> ee: the e+- will be phys prim (because of EM decay), but not fromBGEvent!
  // goal: skip injected signals (e.g. J/psi -> ee)
  // but keep signals from Geant (e.g. photon conversions, but again not from injected J/psi)
  //Printf(" pdg: %s \t fromBGEvent: %d \t isPhysPrim: %d", GetParticleName(pdgT), fromBGEvent, isPhysPrim);
  //
  // if the particle isn't already the physical primary, then iteratively go back to it:
  labelmotherT = mctrack->GetMother();
  while (!isPhysPrim) 
  {
    fromBGEvent  = mcEventLocal->IsFromBGEvent(labelmotherT);     // save property of mother
    isPhysPrim   = mcEventLocal->IsPhysicalPrimary(labelmotherT); // save property of mother
    mother      = fStack->Particle(labelmotherT);
    labelmotherT = mother->GetMother(0); // prepare to get properties of grandmother
    //
    //Printf("  mother: %s \t fromBGEvent: %d \t isPhysPrim: %d", GetParticleName(mother->GetPdgCode()), fromBGEvent, isPhysPrim);
  }
  // now check if the physical primary is from background event. if not, it was injected.
  if (!fromBGEvent) { //Printf("  physical primary is not from background event => injected! skipping. ( %s )", GetParticleName(mother->PdgCode()));
    isInjected = kTRUE;
  }
  return isInjected;
}


//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf(" Now running: Terminate()");

  //get output data and draw 'fHistPt'
//  if (!GetOutputData(0)) return;
//  TH1F *hist=(TH1F*)(((TObjArray*)GetOutputData(0))->FindObject("fHistPt"));
//  if (hist) hist->Draw();
}

//________________________________________________________________________
AliAnalysisTaskElectronEfficiency::~AliAnalysisTaskElectronEfficiency()
{
  /// Destructor
  Printf(" Now running: ~Destructor");
  
  Printf("deleting TList");
  if (fOutputList) {
    fOutputList->Clear();
    delete fOutputList;
  }
  Printf("setting pointer to 0x0");
  fOutputList=0x0;
  
  if (fOutputListSupportHistos) {
    fOutputListSupportHistos->Clear();
    delete fOutputListSupportHistos;
  }
  fOutputListSupportHistos=0x0;
  
  // other objects (pointers) may be deleted like this:
  //  if (fTrackCutsITSSA)  delete fTrackCutsITSSA;
  //  fTrackCutsITSSA=0x0;
}

//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CreateHistograms(TString names, Int_t cutInstance)
{
	Printf(" Now running: CreateHistograms()");

  TObjArray *arrNames=names.Tokenize(";");
  TString name=Form("%02d",cutInstance);
  if (cutInstance<arrNames->GetEntriesFast()){
    name=arrNames->At(cutInstance)->GetName();
  }

  //Printf("%i\t %f\t %i\t %f\t %i\t %f\t ",fNptBins,fPtBins[1],fNetaBins,fEtaBins[1],fNphiBins,fPhiBins[1]);
  TH3F *hNreco_Ele = new TH3F(Form("Nreco_Ele_%s",name.Data()),Form("Nreco_Ele_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  TH3F *hNreco_Ele_poslabel = new TH3F(Form("Nreco_Ele_poslabel_%s",name.Data()),Form("Nreco_Ele_poslabel_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  fvReco_Ele.push_back(hNreco_Ele);
  fvReco_Ele_poslabel.push_back(hNreco_Ele_poslabel);
  
  // one needs the histogram 'hAllPionsForRej' for each cutInstance independently, because empty events may differ between the cutsets which run together.
  TH3F *hAllPionsForRej = new TH3F(Form("AllPionsForRej_%s",name.Data()),Form("AllPionsForRej_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  TH3F *hNPionsRejByAllSigns = new TH3F(Form("NPionsRejByAllSigns_%s",name.Data()),Form("NPionsRejByAllSigns_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  TH3F *hNPionsRejByUnlike = new TH3F(Form("NPionsRejByUnlike_%s",name.Data()),Form("NPionsRejByUnlike_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  fvAllPionsForRej.push_back(hAllPionsForRej);
  fvPionsRejByAllSigns.push_back(hNPionsRejByAllSigns);
  fvPionsRejByUnlike.push_back(hNPionsRejByUnlike);
  
  // be really careful if you need to implement this (see comments in UserExec):
  //  TH3F *hNreco_Pio = new TH3F(Form("Nreco_Pio_%s",name.Data()),Form("Nreco_Pio_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  //  TH3F *hNreco_Kao = new TH3F(Form("Nreco_Kao_%s",name.Data()),Form("Nreco_Kao_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  //  TH3F *hNreco_Pro = new TH3F(Form("Nreco_Pro_%s",name.Data()),Form("Nreco_Pro_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  //  fvReco_Pio.push_back(hNreco_Pio);
  //  fvReco_Kao.push_back(hNreco_Kao);
  //  fvReco_Pro.push_back(hNreco_Pro);
}

//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CreateHistoGen()
{
	Printf(" Now running: CreateHistoGen()");
  
  Printf("fNptBins=%i\t fPtBins[1]=%f\t fNetaBins=%i\t fEtaBins[1]=%f\t fNphiBins=%i\t fPhiBins[1]=%f\t ",fNptBins,fPtBins[1],fNetaBins,fEtaBins[1],fNphiBins,fPhiBins[1]);
  fNgen               = new TH3F("fNgen","fNgen",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
}

//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CreateSupportHistos()
{
	Printf(" Now running: CreateSupportHistos()");
  
  fOutputListSupportHistos = new TList();
  fOutputListSupportHistos->SetName(GetName());
  fOutputListSupportHistos->SetOwner();
  
  // Event variables
  TH1F* hnEvents         = new TH1F("nEvents","Number of processed events after cuts;;N_{events}",1,0.,1.);//,AliDielectronVarManager::kNevents);
  TH1F* hCent            = new TH1F("centrality","N. events vs. centrality;centrality [%];N_{events}",100,0,100);//,AliDielectronVarManager::kCentrality);
  TH2F* hNTrksEvent_cent = new TH2F("hNTrksEvent_cent", "N. tracks vs. centrality (not used);centrality [%];N_{Tracks}", 50,0.,100., 50,-0.5,49.5);
  TH2F* hNEleEvent_cent  = new TH2F("hNEleEvent_cent", "N. selected electrons per event vs. centrality;centrality [%];N_{Tracks}", 50,0.,100., 50,-0.5,49.5);
  TH1F* hVertexZ         = new TH1F("hVertexZ","hVertexZ",300,-15.,15.);
  TH1F* hNvertexCtrb     = new TH1F("hNvertexCtrb","hNvertexCtrb",5000,-0.5,4999.5);
  fOutputListSupportHistos->AddAt(hnEvents, 0);
  fOutputListSupportHistos->AddAt(hCent, 1);
  fOutputListSupportHistos->AddAt(hNTrksEvent_cent, 2);
  fOutputListSupportHistos->AddAt(hNEleEvent_cent, 3);
  fOutputListSupportHistos->AddAt(hVertexZ, 4);
  fOutputListSupportHistos->AddAt(hNvertexCtrb, 5);
  
//  TH1F* hUseless06 = new TH1F("hUseless06","hUseless06", 1,0.,1.);
//  fOutputListSupportHistos->AddAt(hUseless06, 6);
  TH2F* hPInVarMgr_PInStandAlone = new TH2F("PInVarMgr_PInStandAlone","GetTPCInnerParam()->P() vs GetInnerParam->P() (VarMgr method); PIn global tracking [GeV]; PIn TPC standalone [GeV]",
                                            160,0.,8.,160,0.,8.);
  fOutputListSupportHistos->AddAt(hPInVarMgr_PInStandAlone, 6);
  
  // Track variables
  TH1F* hPt      = new TH1F("Pt","Pt;Pt [GeV];#tracks",200,0,10.);//,AliDielectronVarManager::kPt);
  TH2F* hP_PIn   = new TH2F("P_PIn","TPC inner P vs P; P [GeV]; TPC inner P [GeV]",
                        160,0.,8.,160,0.,8.);//,AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  TH2F* hP_Pgen  = new TH2F("P_Pgen","P gen vs P; P [GeV]; P gen [GeV]",
                        160,0.,8.,160,0.,8.);//,AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  fOutputListSupportHistos->AddAt(hPt,     7);
  fOutputListSupportHistos->AddAt(hP_PIn,  8);
  fOutputListSupportHistos->AddAt(hP_Pgen, 9);
  
  
  // ITS
  TH2F* hITS_dEdx_P = new TH2F("ITS_dEdx_P","ITS dEdx;P [GeV];ITS signal (arb units)",
                        160,0.,8.,700,0.,700.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal,makeLogx);
  TH2F* hITSnSigmaEle_P = new TH2F("ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV];ITS number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  TH2F* hITSnSigmaPio_P = new TH2F("ITSnSigmaPio_P","ITS number of sigmas Pions;P [GeV];ITS number of sigmas Pions",
                        160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio,makeLogx);
  TH2F* hITSnSigmaKao_P = new TH2F("ITSnSigmaKao_P","ITS number of sigmas Kaons;P [GeV];ITS number of sigmas Kaons",
                        160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao,makeLogx);
  fOutputListSupportHistos->AddAt(hITS_dEdx_P, 10);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_P, 11);
  fOutputListSupportHistos->AddAt(hITSnSigmaPio_P, 12);
  fOutputListSupportHistos->AddAt(hITSnSigmaKao_P, 13);
  
  // TPC
  TH2F* hTPC_dEdx_P = new TH2F("TPC_dEdx_P","TPC dEdx;P [GeV];TPC signal (arb units)",
                        160,0.,8.,120,0.,120.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,makeLogx);
  TH2F* hTPCnSigmaEle_P = new TH2F("TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
  TH3F* hTPCnSigmaEle_P_dEdx = new TH3F("TPCnSigmaEle_P_dEdx","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;TPC signal (arb units)",
                        80,0.,4.,80,-4.,4.,50,50.,100.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal,makeLogx);
  TH2F* hTPCnSigmaPio_P = new TH2F("TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions",
                        160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,makeLogx);
  TH2F* hTPCnSigmaKao_P = new TH2F("TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons",
                        160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,makeLogx);
  TH2F* hTPCnSigmaPro_P = new TH2F("TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons",
                        160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,makeLogx);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_P, 14);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P, 15);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P_dEdx, 16);
  fOutputListSupportHistos->AddAt(hTPCnSigmaPio_P, 17);
  fOutputListSupportHistos->AddAt(hTPCnSigmaKao_P, 18);
  fOutputListSupportHistos->AddAt(hTPCnSigmaPro_P, 19);
  
  TH2F* hITSnSigmaElePIDresp_P = new TH2F("ITSnSigmaElePIDresp_P","ITS nSigmaEle from VarManager;P [GeV];ITS nSigmaEle from VarManager",
                                         160,0.,8.,100,-5.,5.);
  TH2F* hTPCnSigmaElePIDresp_P = new TH2F("TPCnSigmaElePIDresp_P","TPC nSigmaEle from VarManager;P [GeV];TPC nSigmaEle from VarManager",
                                         160,0.,8.,100,-5.,5.);
  fOutputListSupportHistos->AddAt(hITSnSigmaElePIDresp_P, 20);
  fOutputListSupportHistos->AddAt(hTPCnSigmaElePIDresp_P, 21);
  
  // run dependency
  // run string "fsRunBins" must be sorted in increasing order!
  TObjArray *objaRuns=fsRunBins.Tokenize(", ");
  const Int_t nRuns=objaRuns->GetEntries();
  if (nRuns<1) { fsRunBins = "-99"; cout << "warning: bins for run dependence not specified!" << endl; }
  fsRunBins.Append(Form(", %i", (Int_t) (atoi(objaRuns->At(nRuns-1)->GetName()) + 10))); // create upper limit for bin of last run!
  
  Double_t* dRunBinning = (AliDielectronHelper::MakeArbitraryBinning(fsRunBins))->GetMatrixArray();
  //for (int i=0; i<nRuns+1; i++) { cout << i << " \t " << dRunBinning[i] << " \t "; }
  
  Int_t nbinsPIn=80, nbinsTPCsig=50;
  TH3F* hTPC_dEdx_P_run = new TH3F("TPC_dEdx_P_run","TPC dEdx;P [GeV];TPC signal (arb units);run number",
                                   nbinsPIn   , (AliDielectronHelper::MakeLinBinning(nbinsPIn   ,0.,4.))->GetMatrixArray(), 
                                   nbinsTPCsig, (AliDielectronHelper::MakeLinBinning(nbinsTPCsig,50.,100.))->GetMatrixArray(), 
                                   nRuns      , dRunBinning
                                   );//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kRunNumber);
  TH3F* hTPCnSigmaEle_P_run = new TH3F("TPCnSigmaEle_P_run","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;run number",
                                       nbinsPIn   , (AliDielectronHelper::MakeLinBinning(nbinsPIn   ,0.,4.))->GetMatrixArray(), 
                                       50         , (AliDielectronHelper::MakeLinBinning(50         ,-5.,5.))->GetMatrixArray(), 
                                       nRuns      , dRunBinning
                                       );//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
  TH3F* hITSnSigmaEle_P_run = new TH3F("ITSnSigmaEle_P_run","ITS number of sigmas Electrons;P [GeV];ITS number of sigmas Electrons;run number",
                                       nbinsPIn   , (AliDielectronHelper::MakeLinBinning(nbinsPIn   ,0.,4.))->GetMatrixArray(), 
                                       50         , (AliDielectronHelper::MakeLinBinning(50         ,-5.,5.))->GetMatrixArray(), 
                                       nRuns      , dRunBinning
                                       );//,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kRunNumber);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_P_run, 22);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P_run, 23);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_P_run, 24);
  
  // TOF
  TH2F* hTOFnSigmaEle_P = new TH2F("TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
  TH2F* hTOFbeta = new TH2F("TOFbeta","TOF beta;P [GeV];TOF beta",
                        160,0.,8.,120,0.,1.2);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,makeLogx);
  fOutputListSupportHistos->AddAt(hTOFnSigmaEle_P, 25);
  fOutputListSupportHistos->AddAt(hTOFbeta, 26);
  
  // Eta and Phi
  TH1F* hEta = new TH1F("Eta","Eta; Eta;#tracks",
                        200,-2,2);//,AliDielectronVarManager::kEta);
  TH1F* hPhi = new TH1F("Phi","Phi; Phi;#tracks",
                        320,0.,6.4);//,AliDielectronVarManager::kPhi);
  TH2F* hEta_Phi = new TH2F("Eta_Phi","Eta Phi Map; Eta; Phi",
                        100,-1,1,320,0,6.4);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  fOutputListSupportHistos->AddAt(hEta, 27);
  fOutputListSupportHistos->AddAt(hPhi, 28);
  fOutputListSupportHistos->AddAt(hEta_Phi, 29);
  
  TH2F* hTPC_dEdx_Eta = new TH2F("TPC_dEdx_Eta","TPC dEdx;Eta;TPC signal (arb units)",
                        100,-1,1,120,0.,120.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  TH3F* hTPC_dEdx_Eta_P = new TH3F("TPC_dEdx_Eta_P","TPC dEdx;Eta;TPC signal (arb units); TPC inner P [GeV]",
                        100,-1,1,60,0.,120.,80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
  TH2F* hTPCnSigmaEle_Eta = new TH2F("TPCnSigmaEle_Eta","TPC number of sigmas Electrons; Eta; TPC number of sigmas Electrons",
                        100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  // 3D, may be used for dEdx eta correction
  TH3F* hTPCnSigmaEle_Eta_P = new TH3F("TPCnSigmaEle_Eta_P","TPC number of sigmas Electrons; Eta; TPC number of sigmas Electrons; TPC inner P [GeV]",
                        50,-1,1, 50,-5.,5., 80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
  TH3F* hTPCnSigmaEle_Eta_RefMultTPConly = new TH3F("TPCnSigmaEle_Eta_RefMultTPConly","TPC Ref Mult from AOD header';Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                        50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  TH3F* hTPCnSigmaEle_Eta_Nacc = new TH3F("TPCnSigmaEle_Eta_Nacc","Nacc from 'AliDielectronHelper::GetNacc()';Eta;n#sigma_{ele}^{TPC};N_{acc}",
                        50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_Eta, 30);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_Eta_P, 31);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta, 32);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta_P, 33);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta_RefMultTPConly, 34);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta_Nacc, 35);
//  TH2F* hTPCnSigmaKao_Eta = new TH2F("TPCnSigmaKao_Eta","TPC number of sigmas Kaons; Eta; TPC number of sigmas Kaons",
//                        100,-1,1,200,-10.,10.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
//  TH3F* hTPCnSigmaKao_Eta_P = new TH3F("TPCnSigmaKao_Eta_P","TPC number of sigmas Kaons; Eta; TPC number of sigmas Kaons; TPC inner P [GeV]",
//                        50,-1,1,100,-10.,10.,80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao,AliDielectronVarManager::kPIn);
//  fOutputListSupportHistos->AddAt(hTPCnSigmaKao_Eta, 34);
//  fOutputListSupportHistos->AddAt(hTPCnSigmaKao_Eta_P, 35);
  
  // DCA
  TH1F* hdXY = new TH1F("dXY","dXY;dXY [cm];#tracks",
                        200,-2.,2.);//.,AliDielectronVarManager::kImpactParXY);
  TH1F* hdZ = new TH1F("dZ","dZ;dZ [cm];#tracks",
                        400,-4.,4.);//.,AliDielectronVarManager::kImpactParZ);
  TH2F* hdXY_dZ = new TH2F("dXY_dZ","dXY dZ Map; dXY; dZ",
                        100,-1.,1.,300,-3.,3.);//.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  fOutputListSupportHistos->AddAt(hdXY, 36);
  fOutputListSupportHistos->AddAt(hdZ, 37);
  fOutputListSupportHistos->AddAt(hdXY_dZ, 38);
  
  TH1F* hUseless39 = new TH1F("hUseless","hUseless",1,0,1);
  fOutputListSupportHistos->AddAt(hUseless39, 39);

  // Quality
  TH1F* hTPCcrossedRowsOverFindable = new TH1F("TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",120,0.,1.2);//,AliDielectronVarManager::kNFclsTPCfCross);
  TH1F* hTPCcrossedRows = new TH1F("TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNFclsTPCr);
  TH1F* hTPCnCls = new TH1F("TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC);
  TH1F* hITSnCls = new TH1F("ITSnCls","Number of Clusters ITS;ITS number clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsITS);
  TH1F* hTPCchi2 = new TH1F("TPCchi2","TPC chi2 per Cluster;TPC chi2/Cl;#tracks",100,0.,10.);//.,AliDielectronVarManager::kTPCchi2Cl);
  TH1F* hITSchi2 = new TH1F("ITSchi2","ITS chi2 per Cluster;ITS chi2/Cl;#tracks",100,0.,10.);//.,AliDielectronVarManager::kITSchi2Cl);
  TH1F* hNclsSFracTPC = new TH1F("NclsSFracTPC","Fraction of shared clusters assigned in the TPC;TPC fraction of shared clusters;#tracks",200,0,10.);//.,AliDielectronVarManager::kNclsSFracTPC);
  TH1F* hTPCclsDiff = new TH1F("TPCclsDiff","TPC cluster difference;TPC cluster difference;#tracks",200,0,20.);//.,AliDielectronVarManager::kTPCclsDiff);
  TH1F* hTPCsignalN = new TH1F("TPCsignalN","Number of PID Clusters TPC;TPC number PID clusters;#tracks",160,-0.5,159.5);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
  fOutputListSupportHistos->AddAt(hTPCcrossedRowsOverFindable, 40);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows, 41);
  fOutputListSupportHistos->AddAt(hTPCnCls, 42);
  fOutputListSupportHistos->AddAt(hITSnCls, 43);
  fOutputListSupportHistos->AddAt(hTPCchi2, 44);
  fOutputListSupportHistos->AddAt(hITSchi2, 45);
  fOutputListSupportHistos->AddAt(hNclsSFracTPC, 46);
  fOutputListSupportHistos->AddAt(hTPCclsDiff, 47);
  fOutputListSupportHistos->AddAt(hTPCsignalN, 48);
  
  TH2F* hTPCcrossedRows_TPCnCls = new TH2F("TPCcrossedRows_TPCnCls","TPC crossed rows vs TPC number clusters;TPC number clusters;TPC crossed rows",
                        160,-0.5,159.5,160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
  TH2F* hTPCcrossedRows_Pt = new TH2F("TPCcrossedRows_Pt","TPC crossed rows vs Pt;Pt [GeV];TPC crossed rows",
                        160,0.,8.,160,-0.5,159.5);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_TPCnCls, 49);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_Pt, 50);
  
  
  
  // pair-prefilter histograms
  //2D and 3D histograms
  TH2F* hMeePtee    = new TH2F("InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                               100,0.,1., 100,0.,2.);//AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  TH3F* hMeePhiVOpen= new TH3F("InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                               100,0.,0.5, 100,0.,TMath::Pi(), 50,0.,TMath::Pi()/2.);//AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
  //opening angle and PhiV
  TH2F* hMeeOpen    = new TH2F("InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                               100,0.,1., 100,0.,TMath::Pi());//AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
  TH2F* hMeePhiV    = new TH2F("InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                               100,0.,1., 100,0.,TMath::Pi());//AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
  TH2F* hPteeOpen   = new TH2F("PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                               100,0.,2., 100,0.,TMath::Pi());//AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  TH2F* hPteePhiV   = new TH2F("PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                               100,0.,2., 100,0.,TMath::Pi());//AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
  TH2F* hOpenPhiV   = new TH2F("OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                               100,0.,TMath::Pi(), 100,0.,TMath::Pi());//AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
  fOutputListSupportHistos->AddAt(hMeePtee,     51);
  fOutputListSupportHistos->AddAt(hMeePhiVOpen, 52);
  fOutputListSupportHistos->AddAt(hMeeOpen,     53);
  fOutputListSupportHistos->AddAt(hMeePhiV,     54);
  fOutputListSupportHistos->AddAt(hPteeOpen,    55);
  fOutputListSupportHistos->AddAt(hPteePhiV,    56);
  fOutputListSupportHistos->AddAt(hOpenPhiV,    57);
  
  TH1F* hUseless58 = new TH1F("hUseless","hUseless",1,0,1);
  fOutputListSupportHistos->AddAt(hUseless58, 58);
  TH1F* hUseless59 = new TH1F("hUseless","hUseless",1,0,1);
  fOutputListSupportHistos->AddAt(hUseless59, 59);
  
  
  // histograms for ITS eta and dEdx correction
  TH2F* hITSnSigmaEle_Eta = new TH2F("ITSnSigmaEle_Eta","ITS number of sigmas Electrons; Eta; ITS number of sigmas Electrons",
                                     100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  // 3D, may be used for dEdx eta correction
  TH3F* hITSnSigmaEle_Eta_P = new TH3F("ITSnSigmaEle_Eta_P","ITS number of sigmas Electrons; Eta; ITS number of sigmas Electrons; ITS inner P [GeV]",
                                       50,-1,1, 50,-5.,5., 80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kPIn);
  TH3F* hITSnSigmaEle_Eta_RefMultTPConly = new TH3F("ITSnSigmaEle_Eta_RefMultTPConly","TPC Ref Mult from AOD header';Eta;n#sigma_{ele}^{ITS};N_{TPC ref}",
                                                    50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  TH3F* hITSnSigmaEle_Eta_Nacc = new TH3F("ITSnSigmaEle_Eta_Nacc","Nacc from 'AliDielectronHelper::GetNacc()';Eta;n#sigma_{ele}^{ITS};N_{acc}",
                                          50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kNacc);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta,      60);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta_P,    61);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta_RefMultTPConly, 62);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta_Nacc, 63);
  
  // raw sigma Ele in ITS and TPC
  TH2F* hITSnSigmaEleRaw_P = new TH2F("ITSnSigmaEleRaw_P","ITS number of sigmas Electrons (raw);P [GeV];ITS number of sigmas Electrons (raw)",
                                      160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEleRaw);
  TH2F* hTPCnSigmaEleRaw_P = new TH2F("TPCnSigmaEleRaw_P","TPC number of sigmas Electrons (raw);P [GeV];TPC number of sigmas Electrons (raw)",
                                      160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEleRaw);
  TH2F* hITSnSigmaEleRaw_Eta = new TH2F("ITSnSigmaEleRaw_Eta","ITS number of sigmas Electrons (raw); Eta; ITS number of sigmas Electrons (raw)",
                                        100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEleRaw);
  TH2F* hTPCnSigmaEleRaw_Eta = new TH2F("TPCnSigmaEleRaw_Eta","TPC number of sigmas Electrons (raw); Eta; TPC number of sigmas Electrons (raw)",
                                        100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw);
  fOutputListSupportHistos->AddAt(hITSnSigmaEleRaw_P,     64);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEleRaw_P,     65);
  fOutputListSupportHistos->AddAt(hITSnSigmaEleRaw_Eta,   66);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEleRaw_Eta,   67);
  
}


// taken from AliDielectron.cxx
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetCentroidCorrFunction(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
//  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorrTPC = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDCntrdCorrTPC)  {
    fPostPIDCntrdCorrTPC->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST TPC PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTPC->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrTPC->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTPC->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrTPC->GetXaxis()->GetName());
    }
    printf("\n");
//    fUsedVars->SetBitNumber(varx, kTRUE); // this is also done in AliDielectronPID::IsSelected()
//    fUsedVars->SetBitNumber(vary, kTRUE);
//    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetWidthCorrFunction(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
//  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorrTPC = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDWdthCorrTPC)  {
    fPostPIDWdthCorrTPC->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST TPC PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTPC->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrTPC->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTPC->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrTPC->GetXaxis()->GetName());
    }
    printf("\n");
//    fUsedVars->SetBitNumber(varx, kTRUE); // this is also done in AliDielectronPID::IsSelected()
//    fUsedVars->SetBitNumber(vary, kTRUE);
//    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetCentroidCorrFunctionITS(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  //  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrdITS%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorrITS = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDCntrdCorrITS)  {
    fPostPIDCntrdCorrITS->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST ITS PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    //    fUsedVars->SetBitNumber(varx, kTRUE); // this is also done in AliDielectronPID::IsSelected()
    //    fUsedVars->SetBitNumber(vary, kTRUE);
    //    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetWidthCorrFunctionITS(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  //  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdthITS%d%d%d",varx,vary,varz);
  fPostPIDWdthCorrITS = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDWdthCorrITS)  {
    fPostPIDWdthCorrITS->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST ITS PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    //    fUsedVars->SetBitNumber(varx, kTRUE); // this is also done in AliDielectronPID::IsSelected()
    //    fUsedVars->SetBitNumber(vary, kTRUE);
    //    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
