 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

/*  AliAnalysisTaskPidPidCorrelations.cxx
 *  
 *  Analysis task for PID-PID two-particle correlations using PID form  ITS, TPC, TOF, HMPID
 * 
 *  Authors:	Gyula BENCEDI (WIGNER RCP) gyula.bencedi@cern.ch
 * 		Levente MOLNAR (CNRS-IPHC) levente.molnar@cern.ch
 * 
 *  NOTE: running only on AODs, in PbPb
 * 
 * 
 *  Last modified: Wed Jan 29 15:08:56 CET 2014
 * 
 * 
*/

//_____ ROOT fMyAODHeaders
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "THashList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TList.h"
#include "TFormula.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TPDGCode.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include <vector>
#include <algorithm>
#include <functional>

//_____ ALIROOT headers
#include "AliAnalysisFilter.h"

//_____  Manager, Handler
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//_____ Event pool includes
#include "AliEventPoolManager.h"
// #include "AliPool.h"

//_____ AOD includes
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHMPIDrings.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODInputHandler.h"
#include "AliAODpidUtil.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

//_____ PID includes
#include "AliAODPid.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliPIDCombined.h"

//_____ Additional includes
#include "AliLog.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliGenEventHeader.h"
#include "AliTHn.h" 
#include "AliAnalysisUtils.h"

//_____ CORRFW includes
#include "AliCFContainer.h"
// #include "AliCFGridSparse.h"
// #include "AliCFEffGrid.h"
// #include "AliCFManager.h"
// #include "AliCFCutBase.h"

//_____ AnalysisTask headers
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPidPidCorrelations.h"

using namespace myAliPID;
using namespace std;

ClassImp( AliAnalysisTaskPidPidCorrelations )
ClassImp( AliPidPidCorrelationReducedTrack )

#define sign(x) = (x>0) ? 1 : -1

#define length(a) ( sizeof(a) / sizeof(*a) )

#define pi TMath::Pi()

//________________________________________________________________________
AliAnalysisTaskPidPidCorrelations::AliAnalysisTaskPidPidCorrelations() // All data members should be initialised here
  : AliAnalysisTaskSE()

  , fUseMC ( kFALSE )

  , fMyAODEvent ( 0 )
  , fMyAODHeader ( 0 )
  , fMyAODTrack ( 0 )
  , fPIDResponse ( 0 )
  , fMyPrimVtx ( 0 )
  , fMyMcArray ( 0 )
  , fMyMcHeader ( 0 )
  , fMyMcHandler ( 0 )
  , fPoolMgr ( 0x0 )
  , fMyCFCont ( 0x0 )
  , fMyprimRecoTracksPID ( 0x0 )
  , fMyprimMCParticlesPID ( 0x0 )
  , fMyprimRecoTracksMatchedPID ( 0x0 )
  
//   , fWriteCorrTree ( kTRUE )
//   , fVariablesTreeCorr ( 0x0 )
//   , fCorrVariables ( 0 )

  , fTriggerType ( 1 )
  , fMyMcType ( -1 )
//   , fFilterBit ( 128 )
  , fFilterType ( 2 )
  , fCentrality ( -1 )
  , fCentralityPercentileMin ( 0. )
  , fCentralityPercentileMax ( 10. )
  , fNbinsCent ( 1 )
  , fCentAxis ( 0x0 )
  , fNbinsZvtx ( 1 )
  , fZvtxAxis ( 0x0 )
  , fNbinsPt ( 1 )
  , fPtAxis ( 0x0 )
  , fNbinsEta ( 1 )
  , fEtaAxis ( 0x0 )
  , fCentralityEstimator ( "V0M" )

  , fTrackEtaMin ( -0.9 )
  , fTrackEtaMax ( 0.9 )
  , fTrackPtMin ( 0.2 )
  , fTrackPtMax ( 10.0 )
  , fTrackStatus ( 0 )
  , fnTracksVertex ( 1 )
  , fRejectZeroTrackEvents ( kTRUE )
  , fEtaCut ( 0.9 )
  , fVxMax ( 0.3 )
  , fVyMax ( 0.3 )
  , fVzMax ( 10. )
  , fRemoveWeakDecays ( kFALSE )
  , fRemoveDuplicates ( kFALSE )
  , fDeltaEtaMax ( 2. )
	
  , fSelectCharge ( 0 )
  , fTriggerSelectCharge ( 0 )
  , fAssociatedSelectCharge ( 0 )
  , fTriggerRestrictEta ( -1 )
  , fCutConversions ( kFALSE )
  , fCutResonances ( kFALSE )
  , fRejectResonanceDaughters ( -1 )
  , fOnlyOneEtaSide ( 0 )
  , fWeightPerEvent ( kFALSE )
  , fPtOrder ( kTRUE )
  , fQCut ( kFALSE )
  , fDeltaPtMin ( 0. )
  
  , fPIDtrigType ( 999 )
  , fPIDassocType ( 999 )
  
  , fCustomBinning ( "" )
  , fBinningString ( "" )

  , fMinNumTrack ( 2000 )
  , fPoolSize ( 1000 )
  , fMinNEventsToMix ( 5 )

  , fFillpT ( kFALSE )
  , fTwoTrackEfficiencyCut ( 1 )
  , twoTrackEfficiencyCutValue ( 0.02 )
  , fTwoTrackCutMinRadius ( 0.8 )
  , fEtaOrdering ( kFALSE )
  , fFillMixed ( kTRUE )
  
  , fList ( 0x0 )
  , fHistNev ( 0x0 )
  , fHistTriggerStats ( 0x0 )
  , fHistTriggerRun ( 0x0 )
  , fHistEventStat ( 0x0 )  
  , fHistRefTracks ( 0x0 )
  , fHistCentStats ( 0x0 )
  , fHistCentralityPercentile ( 0x0 )
  , fHistCentralityClass10 ( 0x0 )
  , fHistCentralityClass5 ( 0x0 )
  , fHistV0M ( 0x0 )
  , fHistMcGenerator ( 0x0 )
  , fHist_HijingBg ( 0x0 )
  , fHistNumOfPartPerEvt ( 0x0 )
  , fHistMcStats ( 0x0 )
  , fHistMcAllPt ( 0x0 )
  , fHistMcAllPt_Hijing ( 0x0 )
  , fHistMcAllPt_Dec ( 0x0 )
  , fHistMcAllPt_Inj ( 0x0 )
  , fHistMcAllEta_NotHijing ( 0x0 )
  , fHistMcAllEta_Hijing ( 0x0 )
  , fHistMcAllEta ( 0x0 )
//   , fHistCorrSingle ( 0x0 )
//   , fHistPoolMgrQA ( 0x0 )
//   , fHistCorrPairSame ( 0x0 )
//   , fHistCorrPairMixed ( 0x0 )
  , fHistControlConvResoncances ( 0x0 )
  
  , fHistTriggerWeighting ( 0x0 )
  , fTriggerWeighting ( 0x0 )
  
  , fHistHBTbefore ( 0x0 )
  , fHistHBTafter ( 0x0 )
//   , fHistConversionbefore ( 0x0 )
//   , fHistConversionafter ( 0x0 )
//   , fHistPsiMinusPhi ( 0x0 )
//   , fHistResonancesBefore ( 0x0 )
//   , fHistResonancesRho ( 0x0 )
//   , fHistResonancesK0 ( 0x0 )
//   , fHistResonancesLambda ( 0x0 )
//   , fHistQbefore ( 0x0 )
//   , fHistQafter ( 0x0 )

{
  // Dummy constructor ALWAYS needed for I/O.
}
  
//________________________________________________________________________
AliAnalysisTaskPidPidCorrelations::AliAnalysisTaskPidPidCorrelations(const char *name) // All data members should be initialised here
  : AliAnalysisTaskSE(name)

  , fUseMC ( kFALSE )

  , fMyAODEvent ( 0 )
  , fMyAODHeader ( 0 )
  , fMyAODTrack ( 0 )
  , fPIDResponse ( 0 )
  , fMyPrimVtx ( 0 )
  , fMyMcArray ( 0 )
  , fMyMcHeader ( 0 )
  , fMyMcHandler ( 0 )
  , fPoolMgr ( 0x0 )
  , fMyCFCont ( 0x0 )
  , fMyprimRecoTracksPID ( 0x0 )
  , fMyprimMCParticlesPID ( 0x0 )
  , fMyprimRecoTracksMatchedPID ( 0x0 )

//   , fWriteCorrTree ( kTRUE )
//   , fVariablesTreeCorr ( 0x0 )
//   , fCorrVariables ( 0 )
  
  , fTriggerType ( 1 )
  , fMyMcType ( -1 )
//   , fFilterBit ( 128 )
  , fFilterType ( 2 )
  , fCentrality ( -1 )
  , fCentralityPercentileMin ( 0. )
  , fCentralityPercentileMax ( 10. )
  , fNbinsCent ( 1 )
  , fCentAxis ( 0x0 )
  , fNbinsZvtx ( 1 )
  , fZvtxAxis ( 0x0 )
  , fNbinsPt ( 1 )
  , fPtAxis ( 0x0 )
  , fNbinsEta ( 1 )
  , fEtaAxis ( 0x0 )
  , fCentralityEstimator ( "V0M" )

  , fTrackEtaMin ( -0.9 )
  , fTrackEtaMax ( 0.9 )
  , fTrackPtMin ( 0.2 )
  , fTrackPtMax ( 10.0 )
  , fTrackStatus ( 0 )
  , fnTracksVertex ( 1 )
  , fRejectZeroTrackEvents ( kTRUE )
  , fEtaCut ( 0.9 )
  , fVxMax ( 0.3 )
  , fVyMax ( 0.3 )
  , fVzMax ( 10. )
  , fRemoveWeakDecays ( kFALSE )
  , fRemoveDuplicates ( kFALSE )
  , fDeltaEtaMax ( 2. )
	
  , fSelectCharge ( 0 )
  , fTriggerSelectCharge ( 0 )
  , fAssociatedSelectCharge ( 0 )
  , fTriggerRestrictEta ( -1 )
  , fCutConversions ( kFALSE )
  , fCutResonances ( kFALSE )
  , fRejectResonanceDaughters ( -1 )
  , fOnlyOneEtaSide ( 0 )
  , fWeightPerEvent ( kFALSE )
  , fPtOrder ( kTRUE )
  , fQCut ( kFALSE )
  , fDeltaPtMin ( 0. )
  
  , fPIDtrigType ( 999 )
  , fPIDassocType ( 999 )
  
  , fCustomBinning ( "" )
  , fBinningString ( "" )

  , fMinNumTrack ( 2000 )
  , fPoolSize ( 1000 )
  , fMinNEventsToMix ( 5 )

  , fFillpT ( kFALSE )
  , fTwoTrackEfficiencyCut ( 1 )
  , twoTrackEfficiencyCutValue ( 0.02 )
  , fTwoTrackCutMinRadius ( 0.8 )
  , fEtaOrdering ( kFALSE )
  , fFillMixed ( kTRUE )
  
  , fList ( 0x0 )
  , fHistNev ( 0x0 )
  , fHistTriggerStats ( 0x0 )
  , fHistTriggerRun ( 0x0 )
  , fHistEventStat ( 0x0 )  
  , fHistRefTracks ( 0x0 )
  , fHistCentStats ( 0x0 )
  , fHistCentralityPercentile ( 0x0 )
  , fHistCentralityClass10 ( 0x0 )
  , fHistCentralityClass5 ( 0x0 )
  , fHistV0M ( 0x0 )
  , fHistMcGenerator ( 0x0 )
  , fHist_HijingBg ( 0x0 )  
  , fHistNumOfPartPerEvt ( 0x0 )
  , fHistMcStats ( 0x0 )
  , fHistMcAllPt ( 0x0 )
  , fHistMcAllPt_Hijing ( 0x0 )
  , fHistMcAllPt_Dec ( 0x0 )
  , fHistMcAllPt_Inj ( 0x0 )
  , fHistMcAllEta_NotHijing ( 0x0 )
  , fHistMcAllEta_Hijing ( 0x0 )
  , fHistMcAllEta ( 0x0 )
//   , fHistCorrSingle ( 0x0 )
//   , fHistPoolMgrQA ( 0x0 )
//   , fHistCorrPairSame ( 0x0 )
//   , fHistCorrPairMixed ( 0x0 )
  , fHistControlConvResoncances ( 0x0 )
  
  , fHistTriggerWeighting ( 0x0 )
  , fTriggerWeighting ( 0x0 )
  
  , fHistHBTbefore ( 0x0 )
  , fHistHBTafter ( 0x0 )
//   , fHistConversionbefore ( 0x0 )
//   , fHistConversionafter ( 0x0 )
//   , fHistPsiMinusPhi ( 0x0 )
//   , fHistResonancesBefore ( 0x0 )
//   , fHistResonancesRho ( 0x0 )
//   , fHistResonancesK0 ( 0x0 )
//   , fHistResonancesLambda ( 0x0 )
//   , fHistQbefore ( 0x0 )
//   , fHistQafter ( 0x0 )

{
  // Constructor
  
  //____ THnSparse
  for(Int_t i=0; i<2; i++) { fHistCorrPair[i] = 0x0; }

  for (Int_t i=0; i<fNMaxBinsPt; i++)
    for(Int_t j=0; j<2; j++)
      fHistTwoTrackDistancePt[i][j] = 0x0;
  
  for (Int_t iSpecies=0; iSpecies<AliPID::kSPECIES; iSpecies++)
  {
    nsigmaITS[iSpecies] = -999.;
    nsigmaTPC[iSpecies] = -999.;
    nsigmaTOF[iSpecies] = -999.;
    nsigmaHMPID[iSpecies] = -999.;

    for (Int_t i=0; i<4; i++)
      fnsigmas[iSpecies][i] = 0x0;
  }
  
  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++)
  {
    //_____ QA
    for (Int_t i=0; i < 14; i++) fHistQA[iCent][i] = 0x0;
    for (Int_t i=0; i<6; i++) fHistRefTracksCent[iCent][i] = 0x0;

    //_____ PID
    for (Int_t iSpecies=0; iSpecies<AliPID::kSPECIES; iSpecies++)
    {
      fHistNSigmaTPC[iCent][iSpecies] = 0x0;
      fHistNSigmaTOF[iCent][iSpecies] = 0x0;
    }
    
    fHistTPCdEdx[iCent]		= 0x0;
    fHistTOFbeta[iCent]		= 0x0;
    fHistTPCdEdx_selected[iCent] 	= 0x0;
    fHistTOFbeta_selected[iCent] 	= 0x0;
   
    fHistoNSigma[iCent] = 0x0;
    
//     for (Int_t iPtBinTrig=0; iPtBinTrig<fNMaxBinsPt; iPtBinTrig++)
//     {
//       for (Int_t iPtBinAssoc=0; iPtBinAssoc<fNMaxBinsPt; iPtBinAssoc++)
//       {
// 	fHistDeltaPhi[iCent][iPtBinTrig][iPtBinAssoc]	= 0x0;
// 	fHistDeltaPhiMix[iCent][iPtBinTrig][iPtBinAssoc]	= 0x0;
// 	fHistDphiDeta[iCent][iPtBinTrig][iPtBinAssoc]	= 0x0;
// 	fHistDphiDetaMix[iCent][iPtBinTrig][iPtBinAssoc]	= 0x0;
//       }
//     }
    fHistSingleHadronsPt[iCent]		= 0x0;
//     fHistSingleHadronsPt_Mixed[iCent]	= 0x0;
    fHistSingleHadronsEtaPt[iCent]	= 0x0;
//     fHistSingleHadronsEtaPt_Mixed[iCent]	= 0x0;
  }  

  DefineOutput(1,TList::Class());
  DefineOutput(2,AliCFContainer::Class());
  DefineOutput(3,TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskPidPidCorrelations::~AliAnalysisTaskPidPidCorrelations()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.

  if (fCentAxis)	{ delete fCentAxis; fCentAxis = 0x0; }
  if (fZvtxAxis)	{ delete fZvtxAxis; fZvtxAxis = 0x0; }
  if (fPtAxis) 	{ delete fPtAxis; fPtAxis = 0x0; }
  if (fEtaAxis) 	{ delete fEtaAxis; fEtaAxis = 0x0; }  

  if (fList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fList; fList = 0x0; }
  if (fMyCFCont && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fMyCFCont; fMyCFCont = 0x0;  }
  if (fPoolMgr) { delete fPoolMgr; fPoolMgr = 0x0; }
  if (fPIDResponse) { delete fPIDResponse; fPIDResponse = 0x0; }
}

//_______________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::UserExec(Option_t *)
{
  // _____ get input event
  AliVEvent* event = InputEvent();
  if (!event) { AliError("UserExec(). Could not retrieve event"); return; }

  static int iEvent = 0; iEvent++; 
  fHistNev -> SetBinContent(1,iEvent);

  // _____ get AOD event
  fMyAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fMyAODEvent) { AliError("Cannot get the AOD event"); return; }  
  // _____ get AOD fMyAODHeader
  fMyAODHeader = dynamic_cast<AliAODHeader*>(fMyAODEvent->GetHeader());
  if(!fMyAODHeader) { AliError("Cannot get the AOD fMyAODHeader"); return; }

  if (fUseMC)
  {
    fMyMcArray = dynamic_cast<TClonesArray*>(fMyAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fMyMcArray) { AliError("AliAnalysisTaskPidPidCorrelations::UserExec: Could not find Monte-Carlo in AOD"); return; }
    fMyMcHeader = (AliAODMCHeader*)fMyAODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!fMyMcHeader) { AliError("AliAnalysisTaskPidPidCorrelations::UserExec: MC header branch not found!\n"); return; }
    // MC handler
    fMyMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  }

  //____________ get correlation and correction
  Analyse();
  
  PostData(1, fList);
  PostData(2, fMyCFCont);
//   if (fWriteCorrTree) PostData(3,fVariablesTreeCorr);
}

//____________________________________________________________________
void  AliAnalysisTaskPidPidCorrelations::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("PIDAnalysisSettings","Analysis Settings in PID correlation");
  settingsTree->Branch("fnTracksVertex", &fnTracksVertex,"nTracksVertex/I");
//   settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fTrackStatus", &fTrackStatus,"TrackStatus/I");
  settingsTree->Branch("fRejectZeroTrackEvents", &fRejectZeroTrackEvents,"RejectZeroTrackEvents/O");
  settingsTree->Branch("fRemoveWeakDecays", &fRemoveWeakDecays,"RemoveWeakDecays/O");
  settingsTree->Branch("fRemoveDuplicates", &fRemoveDuplicates,"RemoveDuplicates/O");
  settingsTree->Branch("fVzMax", &fVzMax,"ZVertex/D");
  settingsTree->Branch("fTrackEtaMin", &fTrackEtaMin, "fTrackEtaMin/D");
  settingsTree->Branch("fOnlyOneEtaSide", &fOnlyOneEtaSide,"OnlyOneEtaSide/I");
  settingsTree->Branch("fTrackPtMin", &fTrackPtMin, "fTrackPtMin/D");
  settingsTree->Branch("fTriggerType", &fTriggerType,"fTriggerType/I");
  settingsTree->Branch("fSelectCharge", &fSelectCharge,"SelectCharge/I");
  settingsTree->Branch("fTriggerSelectCharge", &fTriggerSelectCharge,"TriggerSelectCharge/I");
  settingsTree->Branch("fAssociatedSelectCharge", &fAssociatedSelectCharge,"fAssociatedSelectCharge/I");  
  settingsTree->Branch("fTriggerRestrictEta", &fTriggerRestrictEta,"TriggerRestrictEta/D");
  settingsTree->Branch("fEtaOrdering", &fEtaOrdering,"EtaOrdering/O");
  settingsTree->Branch("fCutConversions", &fCutConversions,"CutConversions/O");
  settingsTree->Branch("fCutResonances", &fCutResonances,"CutResonances/O");
  settingsTree->Branch("fRejectResonanceDaughters", &fRejectResonanceDaughters,"RejectResonanceDaughters/I");
  settingsTree->Branch("fFillpT", &fFillpT,"FillpT/O");
  settingsTree->Branch("fMinNumTrack", &fMinNumTrack,"fMinNumTrack/I");
  settingsTree->Branch("fWeightPerEvent", &fWeightPerEvent,"WeightPerEvent/O");
  settingsTree->Branch("fPtOrder", &fPtOrder,"PtOrder/O");
  settingsTree->Branch("fTwoTrackEfficiencyCut", &fTwoTrackEfficiencyCut,"TwoTrackEfficiencyCut/D");
  settingsTree->Branch("twoTrackEfficiencyCutValue", &twoTrackEfficiencyCutValue,"twoTrackEfficiencyCutValue/D");
  settingsTree -> Fill();
  fList -> Add(settingsTree);
}

//_______________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::Analyse() 
{
  fHistEventStat -> Fill("Total",1); // all events

  // store offline trigger bits
  fHistTriggerStats -> Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());  
  
//_________Trigger selection
  Bool_t isSelected = kTRUE; //Bool_t lTrigger = kTRUE;
  // instead of task->SelectCollisionCandidate(mask) in AddTask macro
  isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerType);
  /*  
  Printf("Trigger classes: %s:", fAOD -> GetFiredTriggerClasses().Data());
  Bool_t lTrigger1 = (fAOD -> GetFiredTriggerClasses().Contains("CVLN_R1-B-NOPF-ALLNOTRD")) ? 1 : 0;
  Bool_t lTrigger2 = (fAOD -> GetFiredTriggerClasses().Contains("CINT7WU-B-NOPF-ALL")) ? 1 : 0;
  lTrigger = lTrigger1 && lTrigger2;
  */
  if(!isSelected/* && !lTrigger*/) { /*AliInfo(" === Event Rejected!");*/ fHistEventStat -> Fill("No trigger",1); PostData(1,fList); return; }

  TString trgClasses = fMyAODEvent -> GetFiredTriggerClasses();
//   Printf("GetFiredTriggerClasses: %s",trgClasses.Data());
  
  if (fMyAODEvent->GetFiredTriggerClasses().Contains("ALL")) fHistTriggerRun -> Fill(0);
  if (fMyAODEvent->GetFiredTriggerClasses().Contains("ALLNOTRD")) fHistTriggerRun -> Fill(1);
  if (fMyAODEvent->GetFiredTriggerClasses().Contains("MUON")) fHistTriggerRun -> Fill(2);
  if (fMyAODEvent->GetFiredTriggerClasses().Contains("CENT")) fHistTriggerRun -> Fill(3);
  if (fMyAODEvent->GetFiredTriggerClasses().Contains("CENTNOTRD")) fHistTriggerRun -> Fill(4);
  if (fMyAODEvent->GetFiredTriggerClasses().Contains("TPC")) fHistTriggerRun -> Fill(5);

    
  //_________Centrality
  if (fCentralityEstimator.Length() > 0)
  {
    AliCentrality* centralityObj = 0x0;

    Int_t lCentralityClass10 	= 0;
    Int_t lCentralityClass5 	= 0;  
  
    fMyAODHeader = dynamic_cast<AliAODHeader*>(fMyAODEvent->GetHeader());
      
    // QA for centrality estimators
    fHistCentStats->Fill(0.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("V0M"));
    fHistCentStats->Fill(1.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("FMD"));
    fHistCentStats->Fill(2.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("TRK"));
    fHistCentStats->Fill(3.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("TKL"));
    fHistCentStats->Fill(4.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("CL0"));
    fHistCentStats->Fill(5.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("CL1"));
    fHistCentStats->Fill(6.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
    fHistCentStats->Fill(7.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
    fHistCentStats->Fill(8.,fMyAODHeader->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));

    // centrality QA (V0M)
    fHistV0M->Fill(fMyAODEvent->GetVZEROData()->GetMTotV0A(), fMyAODEvent->GetVZEROData()->GetMTotV0C());
	
    // centrality QA (reference tracks)
    fHistRefTracks->Fill(0.,fMyAODHeader->GetRefMultiplicity());
    fHistRefTracks->Fill(1.,fMyAODHeader->GetRefMultiplicityPos());
    fHistRefTracks->Fill(2.,fMyAODHeader->GetRefMultiplicityNeg());
    fHistRefTracks->Fill(3.,fMyAODHeader->GetTPConlyRefMultiplicity());
    fHistRefTracks->Fill(4.,fMyAODHeader->GetNumberOfITSClusters(0));
    fHistRefTracks->Fill(5.,fMyAODHeader->GetNumberOfITSClusters(1));
  
//     if (fMyAODEvent)
//       centralityObj = fMyAODEvent -> GetCentrality();  	
    centralityObj = fMyAODHeader -> GetCentralityP();

    Bool_t lCentralityInClass = kFALSE;
    if (centralityObj)
    {
//       Printf("======= GetQuality = %f\n",centralityObj->GetQuality());      
//       fCentrality = centralityObj -> GetCentralityPercentileUnchecked(fCentralityEstimator.Data());
      fCentrality = centralityObj -> GetCentralityPercentile(fCentralityEstimator.Data()); 	// returns the centrality percentile, a float Form 0 to 100 (or to the trigger efficiency)
      lCentralityClass10	= centralityObj -> GetCentralityClass10(fCentralityEstimator.Data());//"FMD"); // returns centrality class for 10% width (a integer Form 0 to 10)
      lCentralityClass5	= centralityObj -> GetCentralityClass5(fCentralityEstimator.Data());//"TKL"); // returns centrality class for 5% width (a integer Form 0 to 20)
      lCentralityInClass	= centralityObj -> IsEventInCentralityClass(fCentralityPercentileMin,fCentralityPercentileMax,fCentralityEstimator.Data()); // returns kTRUE if the centrality of the event is between a and b, otherwise kFALSE
      if ( !lCentralityInClass ) { fHistEventStat -> Fill("Wrong centrality", 1); /*Printf("fCentrality = %f\n",fCentrality);*/ PostData(1, fList); return; }
    
      //AliInfo(Form("Centrality is %f",fCentrality));
      fHistCentralityPercentile -> Fill( fCentrality );
      fHistCentralityClass10 -> Fill( lCentralityClass10 );
      fHistCentralityClass5 -> Fill( lCentralityClass5 );
    }
    else
    {
      Printf("WARNING: Centrality object is 0");
      fCentrality = -1.;
    }
  }

  if (fCentrality < 0) return;
  Int_t centBin = GetCentBin(fCentrality);
  if (centBin<0) return;
    
  //___ QA
  fHistRefTracksCent[centBin][0]->Fill(fCentrality,fMyAODHeader->GetRefMultiplicity());
  fHistRefTracksCent[centBin][1]->Fill(fCentrality,fMyAODHeader->GetRefMultiplicityPos());
  fHistRefTracksCent[centBin][2]->Fill(fCentrality,fMyAODHeader->GetRefMultiplicityNeg());
  fHistRefTracksCent[centBin][3]->Fill(fCentrality,fMyAODHeader->GetTPConlyRefMultiplicity());
  fHistRefTracksCent[centBin][4]->Fill(fCentrality,fMyAODHeader->GetNumberOfITSClusters(0));
  fHistRefTracksCent[centBin][5]->Fill(fCentrality,fMyAODHeader->GetNumberOfITSClusters(1));

  //_________Vertex selection
  if(!VertexSelection(fMyAODEvent, fnTracksVertex, centBin, fVxMax, fVyMax, fVzMax)) return;

  //_________check the PIDResponse handler
  if (!fPIDResponse) return;
    
  fHistEventStat -> Fill("Analyzed",1);


  if (!fUseMC)
  {
    //_________get magnetic field
    Double_t bSign		= (fMyAODEvent -> GetMagneticField() > 0) ? 1 : -1;  // for two track efficiency cut
    //Double_t bSignAux	= fMyAODHeader -> GetMagneticField() ; // for dca cut  

    //____ Get AOD RECO tracks
    fMyprimRecoTracksPID = AcceptTracks(centBin,0x0,/*kFALSE,*/kTRUE); // !FIXME onlyprimaries part
    //Printf("Accepted %d fMyprimRecoTracksPID", tracks->GetEntries());
    if(fMyprimRecoTracksPID==0x0) { AliInfo(" ==== fMyprimRecoTracksPID: Zero track pointer"); return; }
    //AliDebug(1, Form("Single Event analysis : nTracks = %4d", fMyprimRecoTracksPID -> GetEntries()));

    if (fRejectZeroTrackEvents && fMyprimRecoTracksPID->GetEntriesFast() == 0)
    {
      AliInfo(Form("========== Rejecting event because it has no tracks: %f %d", fCentrality, fMyprimRecoTracksPID->GetEntriesFast()));
      fHistEventStat -> Fill("NO tracks per event",1);
      //delete fMyprimRecoTracksPID; fMyprimRecoTracksPID=0x0;
      return;
    }
    //delete fMyprimRecoTracksPID;

    const AliVVertex* vertex = fMyAODEvent->GetPrimaryVertex();
    Double_t zVtx = vertex->GetZ();
    
    //_________ Fill Correction AliCFContainer
    FillCFcontainers(0x0, fMyprimRecoTracksPID, 0x0, fCentrality/*, zVtx*/);
    
    //_________ Fill Correlations  
    Double_t weight = 1.;
    if (fFillpT)
      weight = -1.;
  
    //____ step 0
    //FillCorrelations(fMyprimRecoTracksPID, 0x0, fCentrality, zVtx, 0.0, kFALSE, 0.02, kFALSE, /*0,*/ weight);
   
    //____ step 1
    if (fTwoTrackEfficiencyCut > 0)
      FillCorrelations(fMyprimRecoTracksPID, 0x0, fCentrality, zVtx, bSign, kTRUE, twoTrackEfficiencyCutValue, /*1,*/ weight);

    if (fFillMixed)
    {
      // event mixing
    
      // 1. First get an event pool corresponding in mult (cent) and
      //    zvertex to the current event. Once initialized, the pool
      //    should contain nMix (reduced) events. This routine does not
      //    pre-scan the chain. The first several events of every chain
      //    will be skipped until the needed pools are filled to the
      //    specified depth. If the pool categories are not too rare, this
      //    should not be a problem. If they are rare, you could lose
      //    statistics.

      // 2. Collect the whole pool's content of tracks into one TObjArray
      //    (bgTracks), which is effectively a single background super-event.

      // 3. The reduced and bgTracks arrays must both be passed into
      //    FillCorrelations(). Also nMix should be passed in, so a weight
      //    of 1./nMix can be applied.

      AliEventPool* pool = fPoolMgr->GetEventPool(fCentrality, zVtx);
    
      if (!pool)
	AliInfo(Form("No pool found for centrality = %f, zVtx = %f", fCentrality, zVtx));
      else
      {
//       pool->SetDebug(1);
//       pool->PrintInfo();
//       PrintPoolManagerContents();
  
	//if (pool->IsReady())
	if (pool->IsReady() || pool->NTracksInPool() > fMinNumTrack / 10 || pool->GetCurrentNEvents() >= fMinNEventsToMix)
	{
	  Int_t nMix = pool->GetCurrentNEvents();
//       cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
      
	  // ! FIXME later
//       ((TH1F*) fList->FindObject("eventStat"))->Fill(2);
//       ((TH2F*) fList->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
//       if (pool->IsReady())
// 	((TH1F*) fList->FindObject("eventStat"))->Fill(3);
    
	  // Fill mixed-event histos here  
	  for (Int_t jMix=0; jMix<nMix; jMix++)
	  {
	    TObjArray* bgTracks = pool->GetEvent(jMix);
	    if (!bgTracks) continue;
	  
	    if (fTwoTrackEfficiencyCut > 0)
	      FillCorrelations(fMyprimRecoTracksPID, bgTracks, fCentrality, zVtx, bSign, kTRUE, twoTrackEfficiencyCutValue, /*2,*/ 1./nMix);
	  }
	}
	pool->UpdatePool(fMyprimRecoTracksPID);
	//pool->PrintInfo();
      }//_____ pool NULL check
    }//_______ run mixing
    else
      delete fMyprimRecoTracksPID;
  }
  else
  {
    //______ check MC
    Bool_t isMCok = kFALSE;
    isMCok = CheckMcDistributions(fMyMcArray,fMyMcHeader);
    if (isMCok == kFALSE) return;

    //____ Get MC primaries
    fMyprimMCParticlesPID = AcceptMcTracks(centBin,kTRUE,kTRUE);
    if (fMyprimMCParticlesPID==0x0) { AliInfo(" ==== fMyprimMCParticlesPID: Zero track pointer"); return; }
    //CleanUp(primMCParticles, fMyMcArray, skipParticlesAbove);
//     delete fMyprimMCParticlesPID; fMyprimMCParticlesPID=0x0;
    //____ Get MC RECO-matched
    fMyprimRecoTracksMatchedPID = AcceptTracks(centBin,fMyMcArray,/*kTRUE,*/kTRUE);
    if (fMyprimRecoTracksMatchedPID==0x0) { AliInfo(" ==== Zero track pointer"); return; }
    //TObjArray* fMyprimRecoTracksMatchedPID	= AcceptMcRecoMachedTracks(centBin,kTRUE,kTRUE);
    //CleanUp(fMyprimRecoTracksMatchedPID, fMyMcArray, skipParticlesAbove);
//     delete fMyprimRecoTracksMatchedPID; fMyprimRecoTracksMatchedPID=0x0;

    //___ check pdg
    for (Int_t i=0; i<fMyMcArray->GetEntriesFast(); i++)
      ((TH1F*) fList->FindObject("pids"))->Fill(((AliAODMCParticle*)fMyMcArray->At(i))->GetPdgCode());

    //______ Fill Corrections
    FillCFcontainers(fMyprimMCParticlesPID, 0x0, fMyprimRecoTracksMatchedPID, fCentrality/*, zVtx*/);
//     Printf("%d %d %d" fMyprimMCParticlesPID->GetEntries(), fMyprimRecoTracksPID->GetEntries(), fMyprimRecoTracksMatchedPID->GetEntries());
//     delete fMyprimMCParticlesPID;
//     delete fMyprimRecoTracksMatchedPID;
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::FillCorrelations(TObjArray* particles, TObjArray* particlesMixed, Double_t centrality, Double_t zVtx, Double_t bSign, Bool_t twoTrackEfficiencyCut, Double_t fTwoTrackEfficiencyCutValue, /*Int_t step,*/ Double_t weight)
{
  if (!particles){ AliInfo(" ==============  particles TObjArray is NULL pointer, return"); return; }

//   Double_t* trackVariablesSingle = new Double_t[kTrackVariablesSingle];
  Double_t* trackVariablesPair = new Double_t[kTrackVariablesPair];
  
  Bool_t fillpT = kFALSE;
  if (weight < 0) fillpT = kTRUE;

  // define end of particle loops
  Int_t iMax = particles -> GetEntriesFast();
  Int_t jMax = iMax;
  if (particlesMixed)
    jMax = particlesMixed -> GetEntriesFast();

  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* input = (particlesMixed) ? particlesMixed : particles;

  TArrayF secondEta(jMax);
  TArrayF secondPhi(jMax);
  TArrayF secondPt(jMax);
//   TArrayS secondCharge(jMax);
  TArrayF secondPID(jMax);

  for (Int_t i=0;i < jMax;i++)
  {
    secondEta[i]		= ((AliPidPidCorrelationReducedTrack*)input->At(i))->Eta();
    secondPhi[i]		= ((AliPidPidCorrelationReducedTrack*)input->At(i))->Phi();
    secondPt[i]		= ((AliPidPidCorrelationReducedTrack*)input->At(i))->Pt();
//     secondCharge[i]	= (Short_t)((AliPidPidCorrelationReducedTrack*)input->At(i))->Charge();
    secondPID[i]		= ((AliPidPidCorrelationReducedTrack*)input->At(i))->GetMyPartID();
  }


  //_______ Reject ResonanceDaughters 
	
  // identify K, Lambda candidates and flag those particles
  // a TObject bit is used for this
  const UInt_t kResonanceDaughterFlag = 1 << 14;
  if (fRejectResonanceDaughters > 0)
  {
      Double_t resonanceMass = -1;
      Double_t massDaughter1 = -1;
      Double_t massDaughter2 = -1;
      const Double_t interval = 0.02;
      
      switch (fRejectResonanceDaughters)
      {
	case 1: resonanceMass = 1.2; 	massDaughter1 = 0.1396; massDaughter2 = 0.9383; 		break; // method test
	case 2: resonanceMass = 0.4976;	massDaughter1 = 0.1396; massDaughter2 = massDaughter1; 	break; // k0
	case 3: resonanceMass = 1.115; 	massDaughter1 = 0.1396; massDaughter2 = 0.9383; 		break; // lambda
	default: AliFatal(Form("Invalid setting %d", fRejectResonanceDaughters));
      }

      for (Int_t i=0; i<particles->GetEntriesFast(); i++)
	particles->At(i)->ResetBit(kResonanceDaughterFlag);
      
      if (particlesMixed)
	for (Int_t i=0; i<jMax; i++)
	  particlesMixed->At(i)->ResetBit(kResonanceDaughterFlag);
      
      for (Int_t i=0; i<particles->GetEntriesFast(); i++)
      {
	//AliVParticle* triggerParticle = (AliVParticle*) particles->UncheckedAt(i);
	AliPidPidCorrelationReducedTrack* triggerParticle = (AliPidPidCorrelationReducedTrack*) particles->At(i);

	for (Int_t j=0; j<jMax; j++)
	{
	  if (!particlesMixed && i == j)	continue;
	
	  AliVParticle* particle = 0;
	  if (!particlesMixed)
	  {
	    particle = (AliPidPidCorrelationReducedTrack*) particles->At(j);
	    //particle = (AliVParticle*) particles->UncheckedAt(j);
	  }
	  else
	  {
	    particle = (AliPidPidCorrelationReducedTrack*) particlesMixed->At(j);
	    //particle = (AliVParticle*) particlesMixed->UncheckedAt(j);
	  }

	  // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event)
// 	  if (particlesMixed && triggerParticle->IsEqual(particle)) continue;
	 
	  if (triggerParticle->Charge() * particle->Charge() > 0) continue;
      
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerParticle->Eta(), triggerParticle->Phi(), particle->Pt(), particle->Eta(), particle->Phi(), massDaughter1, massDaughter2);
	      
	  if (TMath::Abs(mass - resonanceMass*resonanceMass) < interval*5)
	  {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerParticle->Eta(), triggerParticle->Phi(), particle->Pt(), particle->Eta(), particle->Phi(), massDaughter1, massDaughter2);

	    if (mass > (resonanceMass-interval)*(resonanceMass-interval) && mass < (resonanceMass+interval)*(resonanceMass+interval))
	    {
	      triggerParticle->SetBit(kResonanceDaughterFlag);
	      particle->SetBit(kResonanceDaughterFlag);	      
	      //Printf("Flagged %d %d %f", i, j, TMath::Sqrt(mass));
	    }
	  }
	}
      }
    }//________ fRejectResonanceDaughters

  //_______ fHistTriggerWeighting
  if (fWeightPerEvent)
  {
    for (Int_t i=0; i<particles->GetEntriesFast(); i++)
    {
      AliVParticle* triggerParticle = (AliVParticle*) particles->At(i);
	
      // some optimization
      Double_t triggerEta = triggerParticle->Eta();

      if (fTriggerRestrictEta > 0 && TMath::Abs(triggerEta) > fTriggerRestrictEta)
	continue;
      
      if (fOnlyOneEtaSide != 0)
	if (fOnlyOneEtaSide * triggerEta < 0)
	  continue;
      if (fTriggerSelectCharge != 0)
	if (triggerParticle->Charge() * fTriggerSelectCharge < 0)
	  continue;
	
      fHistTriggerWeighting->Fill(triggerParticle->Pt());
     }
  }

  //_______ trigger particle loop
  for (Int_t i=0; i < iMax; i++)
  {
//     AliVParticle* triggerParticle = (AliVParticle*) particles -> UncheckedAt(i);
    AliPidPidCorrelationReducedTrack* triggerParticle = dynamic_cast<AliPidPidCorrelationReducedTrack*>(particles -> At(i));
    if (!triggerParticle) { AliInfo(" ===== triggerParticle: zero pointer. Continue...\n"); continue; }

    Double_t triggerEta	= triggerParticle -> Eta();
    Double_t triggerPhi	= triggerParticle -> Phi();
    Double_t triggerPt	= triggerParticle -> Pt();
    Int_t triggerPID	= triggerParticle -> GetMyPartID();
    
    if (fTriggerRestrictEta > 0 && TMath::Abs(triggerEta) > fTriggerRestrictEta)
      continue;
	
    if (fOnlyOneEtaSide != 0)
      if (fOnlyOneEtaSide * triggerEta < 0)
	continue;
      
    if (fTriggerSelectCharge != 0)
      if (triggerParticle->Charge() * fTriggerSelectCharge < 0)
	continue;

    if (fRejectResonanceDaughters > 0)
      if (triggerParticle->TestBit(kResonanceDaughterFlag))
      {
	//Printf("Skipped i=%d", i);
	continue;
      }
		
      //Int_t CentBin = fHistCorrPairSame->GetGrid(0)->GetGrid()->GetAxis(7)->FindBin(centrality);//GetCentBin(centrality);
      //Int_t ZvtxBin = fHistCorrPairSame->GetGrid(0)->GetGrid()->GetAxis(6)->FindBin(zVtx);//GetZvtxBin(zVtx);

    //fEtaSpectraTrig[CentBin][ZvtxBin] -> Fill(triggerEta);

    // ___ trigger eta PID
    //if ( triggerPID == fPartPionPlus || triggerPID == fPartPionMinus ) 		fEtaSpectraTrigPion[CentBin][ZvtxBin] -> Fill(triggerEta);
    //if ( triggerPID == fPartKaonPlus || triggerPID == fPartKaonMinus ) 		fEtaSpectraTrigKaon[CentBin][ZvtxBin] -> Fill(triggerEta);
    //if ( triggerPID == fPartProtonPlus || triggerPID == fPartProtonMinus )	fEtaSpectraTrigProton[CentBin][ZvtxBin] -> Fill(triggerEta);

    //___ fill single particle histograms
/*
    trackVariablesSingle[0]    = triggerPID;
    trackVariablesSingle[1]    = triggerPt;
    trackVariablesSingle[2]    = zVtx;
    trackVariablesSingle[3]    = centrality;
    trackVariablesSingle[4]    = triggerParticle->Charge()>0 ? 7 : 8;
    */
//     if( fWriteCorrTree == kTRUE && !particlesMixed)
//     {
//       fCorrVariables[0] = triggerPID;
//       fCorrVariables[1] = triggerPt;
//       fCorrVariables[2] = zVtx;
//       fCorrVariables[3] = centrality;
//       fCorrVariables[4] = triggerParticle->Charge()>0 ? 7 : 8;
//       
//       fVariablesTreeCorr->Fill();
//     }
    
    if (fillpT) weight = triggerPt;	
    Double_t useWeight = weight;

    if (!particlesMixed)
    {
//       fHistCorrSingle -> Fill(trackVariablesSingle,step,useWeight);
//       fHistCorrSingle -> Fill(trackVariablesSingle,useWeight);
      
      fHistSingleHadronsPt[GetCentBin(centrality)]->Fill(triggerPt);
      fHistSingleHadronsEtaPt[GetCentBin(centrality)]->Fill(triggerPt,triggerEta);      
    }
    
    //_______ assoc particle loop
    for (Int_t j=0; j < jMax; j++)
    {
      // no auto correlations (only for non mixing)
      if (!particlesMixed && i == j) continue;

      //AliVParticle* assocParticle = 0x0;
      AliPidPidCorrelationReducedTrack* assocParticle = 0x0;
      if (!particlesMixed)
      {
	//assocParticle = (AliVParticle*) particles->At(j);
	assocParticle = (AliPidPidCorrelationReducedTrack*) particles->At(j);
      }
      else
      {
	//assocParticle = (AliVParticle*) mixed->At(j);
      	assocParticle = (AliPidPidCorrelationReducedTrack*) particlesMixed->At(j);
      }

      if (!assocParticle) continue;

      // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event)
      //! FIXME
      //       if (particlesMixed && triggerParticle -> IsEqual(assocParticle))
      // 	continue;
        
      if (fPtOrder)
	if (assocParticle->Pt() >= triggerParticle->Pt())
	  continue;
	
      if (fAssociatedSelectCharge != 0)
	if (assocParticle->Charge() * fAssociatedSelectCharge < 0)
	  continue;

      if (fSelectCharge > 0)
      {
      	// skip like sign
        if (fSelectCharge == 1 && assocParticle->Charge() * triggerParticle->Charge() > 0)	continue;
	// skip unlike sign
        if (fSelectCharge == 2 && assocParticle->Charge() * triggerParticle->Charge() < 0)	continue;
      }
        
      if (fEtaOrdering)
      {
	if (triggerEta < 0 && secondEta[j] < triggerEta)	continue;
	if (triggerEta > 0 && secondEta[j] > triggerEta)	continue;
      }

      if (fRejectResonanceDaughters > 0)
	if (triggerParticle->TestBit(kResonanceDaughterFlag))
	{
	  //Printf("Skipped i=%d", i);
	  continue;
	}

      // conversions
      if (fCutConversions && assocParticle->Charge() * triggerParticle->Charge() < 0)
      {
	Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.510e-3, 0.510e-3);
	  
	if (mass < 0.1)
	{
	  mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.510e-3, 0.510e-3);
	  fHistControlConvResoncances->Fill(0.0, mass);

	  if (mass < 0.04*0.04) continue;
	}
      }

      // K0s
      if (fCutResonances && assocParticle->Charge() * triggerParticle->Charge() < 0)
      {
	Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.1396, 0.1396);
	  
	const Float_t kK0smass = 0.4976;
	  
	if (TMath::Abs(mass - kK0smass*kK0smass) < 0.1)
	{
	  mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.1396, 0.1396);
	  fHistControlConvResoncances->Fill(1, mass - kK0smass*kK0smass);

	  if (mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))
	    continue;
	}
      }

      // Lambda
      if (fCutResonances && assocParticle->Charge() * triggerParticle->Charge() < 0)
      {
	Float_t mass1 = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.1396, 0.9383);
	Float_t mass2 = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.9383, 0.1396);
	  
	const Float_t kLambdaMass = 1.115;

	if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < 0.1)
	{
	  mass1 = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.1396, 0.9383);
	  fHistControlConvResoncances->Fill(2, mass1 - kLambdaMass*kLambdaMass);
	    
	  if (mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
	    continue;
	}
	if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < 0.1)
	{
	  mass2 = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), assocParticle->Pt(), secondEta[j], assocParticle->Phi(), 0.9383, 0.1396);
	  fHistControlConvResoncances->Fill(2, mass2 - kLambdaMass*kLambdaMass);

	  if (mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
	    continue;
	}
      }

      //_______ two-Track Efficiency Cut , HBT-like cut
      if (twoTrackEfficiencyCut)
      {
	// the variables & cuthave been developed by the HBT group 
	// see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
	
	Double_t phi1	= triggerParticle -> Phi();
	//Double_t phi1	= triggerParticle -> Phi()*TMath::DegToRad();
	Double_t pt1	= triggerParticle -> Pt();
	Double_t charge1	= triggerParticle -> Charge();
		    
	//Double_t phi2	= assocParticle -> Phi()*TMath::DegToRad();
	Double_t phi2	= assocParticle -> Phi();
	Double_t pt2	= assocParticle -> Pt();
	Double_t charge2	= assocParticle -> Charge();

	Double_t deta = triggerEta - secondEta[j];
	Double_t dphi = triggerPhi - secondPhi[j];
		      
	fHistHBTbefore -> Fill(deta,dphi);

	// optimization
	if (TMath::Abs(deta) < fTwoTrackEfficiencyCutValue * 2.5 * 3) // twoTrackEfficiencyCutValue = 0.02 [default for dphicorrelations]
	{
	  // check first boundaries to see if is worth to loop and find the minimum
	  Double_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
	  Double_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
		  	  
	  const Double_t kLimit = fTwoTrackEfficiencyCutValue * 3;
		
	  Double_t dphistarminabs = 1e5;
	  Double_t dphistarmin = 1e5;
	  if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	  {
	    for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
	    {
	      Double_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);
	      Double_t dphistarabs = TMath::Abs(dphistar);
	      
	      if (dphistarabs < dphistarminabs)
	      {
		dphistarmin = dphistar;
		dphistarminabs = dphistarabs;
	      }
	    }

	    Int_t PtBin = GetPtBin(pt2);
// 	    Printf("================ pt = %f \t <--> \t PtBin = %d\n",pt2,PtBin);
	    //fHistTwoTrackDistancePt[0] -> Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2)); // !FIXME
	    fHistTwoTrackDistancePt[PtBin][0] -> Fill(deta, dphistarmin);
 
	    if (dphistarminabs < fTwoTrackEfficiencyCutValue && TMath::Abs(deta) < fTwoTrackEfficiencyCutValue)
	    {
	      // Printf("Removed track pair %d %d with %f %f %f %f %f %f %f %f %f", i, j, deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
	      // AliInfo(Form("HBT: Removed track pair %d %d with [[%f %f]] %f %f %f | %f %f %d %f %f %d %f", i, j, deta, dphi, dphistarminabs, dphistar1, dphistar2, phi1rad, pt1, charge1, phi2rad, pt2, charge2, bSign));
	      continue;
	    }
		
	    //fHistTwoTrackDistancePt[1] -> Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2)); // !FIXME
	    fHistTwoTrackDistancePt[PtBin][1] -> Fill(deta, dphistarmin);
	  }
	}
	fHistHBTafter -> Fill(deta,dphi);
      }//____ HBT cut

      trackVariablesPair[0] =  triggerPt;				// pt trigger
      trackVariablesPair[1] =  secondPt[j];			// pt assoc
      trackVariablesPair[2] =  DeltaPhi(triggerPhi-secondPhi[j]);	// delta phi
      trackVariablesPair[3] =  triggerEta - secondEta[j];		// delta eta
      trackVariablesPair[4] =  zVtx;				// z of the primary vertex
      
//       trackVariablesPair[0] =  triggerPID; 			// trigger PID
// //       trackVariablesPair[1] =  secondPID[j];			// assoc PID
//       trackVariablesPair[1] =  triggerPt;				// pt trigger
//       trackVariablesPair[2] =  secondPt[j];			// pt assoc
//       trackVariablesPair[3] =  DeltaPhi(triggerPhi-secondPhi[j]);	// delta phi
//       trackVariablesPair[4] =  triggerEta - secondEta[j];		// delta eta
//       trackVariablesPair[5] =  zVtx;				// z of the primary vertex
// //       trackVariablesPair[7] =  triggerParticle->Charge()>0 ? 7 : 8;
//       trackVariablesPair[6] =  assocParticle->Charge()>0 ? 7 : 8;
      
      //_______ momentum difference cut - suppress femtoscopic effects
      if (fQCut)
      {
	//Double_t ptMin        = 0.1; //const for the time being (should be changeable later on)
	Double_t ptDifference = TMath::Abs( triggerPt - secondPt[j]);
// 	fHistQbefore -> Fill(trackVariablesPair[1],trackVariablesPair[2],ptDifference);
	if (ptDifference < fDeltaPtMin) continue;
// 	fHistQafter -> Fill(trackVariablesPair[1],trackVariablesPair[2],ptDifference);
      }

      if (fWeightPerEvent)
      {
	Int_t weightBin = fHistTriggerWeighting->GetXaxis()->FindBin(trackVariablesPair[0]);
	//Printf("Using weight %f", fHistTriggerWeighting->GetBinContent(weightBin));
	useWeight /= fHistTriggerWeighting->GetBinContent(weightBin);
      }
    
      if (!particlesMixed)
      {
	if (triggerPID == fPIDtrigType && secondPID[j] == fPIDassocType) fHistCorrPair[0] -> Fill(trackVariablesPair,useWeight);
// 	fHistCorrPairSame -> Fill(trackVariablesPair,step,useWeight);
/*
	Int_t ptBinTrackTrig  = fPtAxis -> FindBin( triggerParticle->Pt() );
	Int_t ptBinTrackAssoc = fPtAxis -> FindBin( assocParticle->Pt() );

	if (ptBinTrackTrig<1 || ptBinTrackTrig>fNbinsPt || ptBinTrackAssoc<1 || ptBinTrackAssoc>fNbinsPt) continue;

	fHistDeltaPhi[GetCentBin(centrality)][ptBinTrackTrig-1][ptBinTrackAssoc-1]	-> Fill(DeltaPhi(triggerParticle->Phi()-assocParticle->Phi()));
	fHistDphiDeta[GetCentBin(centrality)][ptBinTrackTrig-1][ptBinTrackAssoc-1]	-> Fill(DeltaPhi(triggerParticle->Phi()-assocParticle->Phi()),triggerParticle->Eta()-assocParticle->Eta());
	fHistTracksEtaTrigVsEtaAssoc[GetCentBin(centrality)]			-> Fill(triggerParticle -> Eta(),assocParticle -> Eta());*/
      }
      else if (particlesMixed)
      {
	if (triggerPID == fPIDtrigType && secondPID[j] == fPIDassocType) fHistCorrPair[1] -> Fill(trackVariablesPair,useWeight);
// 	fHistCorrPairMixed -> Fill(trackVariablesPair,step,useWeight);
	
// 	Int_t ptBinTrackTrigM  = fPtAxis -> FindBin( triggerParticle->Pt() );
// 	Int_t ptBinTrackAssocM = fPtAxis -> FindBin( assocParticle->Pt() );
// 
// 	if (ptBinTrackTrigM<1 || ptBinTrackTrigM>fNbinsPt || ptBinTrackAssocM<1 || ptBinTrackAssocM>fNbinsPt) continue;
// 
// 	fHistDeltaPhiMix[GetCentBin(centrality)][ptBinTrackTrigM-1][ptBinTrackAssocM-1]	-> Fill(DeltaPhi(triggerParticle->Phi()-assocParticle->Phi()));
// 	fHistDphiDetaMix[GetCentBin(centrality)][ptBinTrackTrigM-1][ptBinTrackAssocM-1]	-> Fill(DeltaPhi(triggerParticle->Phi()-assocParticle->Phi()),triggerParticle->Eta()-assocParticle->Eta());
// 	fHistTracksEtaTrigVsEtaAssocMixed[GetCentBin(centrality)]				-> Fill(triggerParticle->Eta(),assocParticle->Eta());	
      }
    }//_____ assoc loop
  }//_____ trigger loop
//   delete[] trackVariablesSingle;
  delete[] trackVariablesPair;
}

//_____________________________________________________
Bool_t AliAnalysisTaskPidPidCorrelations::CheckMcDistributions(TClonesArray* arrayMC, AliAODMCHeader* mcHeader)
{
  Int_t myPDGproton = 2212;
  Int_t myPDGkaon = 321;
  Int_t myPDGpion = 211;
//   Int_t tmpPDG = -1;

  //____ counters
  Int_t kCntr =0;
  Int_t kCntrInjected = 0;
  Int_t kCntrFromDecay = 0;
   
  Int_t entrMC = arrayMC->GetEntriesFast();
  if ( entrMC < 1 ) return kFALSE;
    //
  AliAODMCParticle* mcPart = 0x0;
        
  TString kGenName="";
//     Bool_t isFromHijing = kFALSE;
    
  for(Int_t ie=0 ;ie < entrMC;ie++)
  {
//     isFromHijing = kFALSE;

    mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ie));
    if(!ie) continue;
        
    if ( (TMath::Abs(mcPart->GetPdgCode()) != myPDGproton) && (TMath::Abs(mcPart->GetPdgCode()) != myPDGkaon) && (TMath::Abs(mcPart->GetPdgCode()) != myPDGpion) ) continue;

    if( mcPart->IsPrimary() == kFALSE ) continue;

    kCntr++;

    ((TH1F*)(fList->FindObject("fHistMcStats")))->Fill(0);
    //
    kGenName = GetGenerator(ie,mcHeader);
    //AliInfo(Form("====> IsPrimary: %d IsPhysicsPrimary %d Gen %s",mcPart->IsPrimary(),mcPart->IsPhysicalPrimary(),kGenName.Data()));
		
    FillMcGeneratorHistos(kGenName);

    if ( kGenName.Contains("Hijing")) fMyMcType = 0;
    else if ( kGenName.Contains("Pythia")) fMyMcType = 1;
    else fMyMcType = -1;
                
    Int_t idxmother = mcPart->GetMother();
    if (idxmother == -1) kCntrInjected++;
    else kCntrFromDecay++;
        
    ((TH1F*)(fList->FindObject("fHistMcAllPt")))->Fill(mcPart->Pt());
    if (idxmother == -1 )	((TH1F*)(fList->FindObject("fHistMcAllPt_Inj")))->Fill(mcPart->Pt());
    if (idxmother > -1 )	((TH1F*)(fList->FindObject("fHistMcAllPt_Dec")))->Fill(mcPart->Pt());

    ((TH1F*)(fList->FindObject("fHistMcAllEta")))->Fill(mcPart->Eta());
        
    if (kGenName.Contains("ijing")) //_____ from HJING
    {
//     isFromHijing = kTRUE;
      ((TH1F*)(fList->FindObject("fHistMcAllPt_Hijing")))->Fill(mcPart->Pt());
      ((TH1F*)(fList->FindObject("fHistMcAllEta_Hijing")))->Fill(mcPart->Eta());            
    }
    else
      ((TH1F*)(fList->FindObject("fHistMcAllEta_NotHijing")))->Fill(mcPart->Eta());        
  } // mc particle loop

  ((TH1F*)(fList->FindObject("fHistNumOfPartPerEvt")))->Fill(kCntr);
  //AliInfo(Form(" === MC: \t injected \t %d, \t Decay: \t %d\n",kCntrInjected,kCntrFromDecay));
    
  return kTRUE;
}

//________________________________________________________________________
TString AliAnalysisTaskPidPidCorrelations::GetGenerator(Int_t label, AliAODMCHeader* MCheader)
{
  // taken from .../PWGHF/vertexingHF/AliVertexingHFUtils.cxx

  // get the name of the generator that produced a given particle
  Int_t nsumpart = 0;
  TList*lh = MCheader -> GetCocktailHeaders();
  Int_t nh = lh -> GetEntries();
  for(Int_t i=0;i < nh;i++)
  {
    AliGenEventHeader* gh = (AliGenEventHeader*)lh -> At(i);
    TString genname = gh -> GetName();
    Int_t npart = gh -> NProduced();
    if ( label >= nsumpart && label < (nsumpart+npart) )	return genname;
    nsumpart += npart;
  }
  TString empty="";
  return empty;
}

//________________________________________________
Bool_t AliAnalysisTaskPidPidCorrelations::IsFromHijingBg(Int_t mcTrackLabel)
{
  Bool_t isFromHijingBg = kFALSE;
  
  if( mcTrackLabel < 0 ) return isFromHijingBg;
  
//   AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fMyMcArray->At(mcTrackLabel));
  
  if ( !mcTrackLabel ) return isFromHijingBg;
  if ( mcTrackLabel < 0 ) return isFromHijingBg;
  
  TString nameGen = GetGenerator(mcTrackLabel,fMyMcHeader);
  
  if ( nameGen.Contains("ijing") ) isFromHijingBg = kTRUE;
  if ( isFromHijingBg == kTRUE )
  {
    AliInfo(" === isFromHijingBg = kTRUE === \n\n");
    ((TH1F*)(fList->FindObject("fHist_HijingBg"))) -> Fill(0);
  }
  else
    ((TH1F*)(fList->FindObject("fHist_HijingBg"))) -> Fill(1);
    
  return isFromHijingBg;
}

//________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::FillMcGeneratorHistos(TString genLabel)
{
  Int_t fillVal = 0;
    
  if ( genLabel.Contains("hijing_0"))  	fillVal = 0;
  else if ( genLabel.Contains("pythia_1")) fillVal = 1;
  else
  {
    for(Int_t i=3;i<=20;i++)
    {
      if( ((TH1F*)(fList->FindObject("fHistMcGenerator")))->GetBinContent(i) < 1 )
      {
	fillVal = i--;
        ((TH1F*)(fList->FindObject("fHistMcGenerator")))->GetXaxis()->SetBinLabel(i,Form("%s",genLabel.Data()));
        break;
      }
    }
  }

  ((TH1F*)(fList->FindObject("fHistMcGenerator")))->Fill(fillVal);
  
  return;
}

//________________________________________________________________________
Bool_t  AliAnalysisTaskPidPidCorrelations::VertexSelection(TObject* obj, Int_t ntracks, Int_t centBin, Double_t gVxMax, Double_t gVyMax, Double_t gVzMax)
{
	if (obj->InheritsFrom("AliAODEvent"))
	{
  	Int_t nVertex = ((AliAODEvent*)obj) -> GetNumberOfVertices();
  
		if ( nVertex > 0 )
  	{
	fMyPrimVtx = (AliAODVertex*)((AliAODEvent*)obj)->GetPrimaryVertex();
    	if(!fMyPrimVtx) return kFALSE;

    	Bool_t lReconstructedVertexPresent = kFALSE;
    	Bool_t lVertexAcceptable = kFALSE;  
    	Double_t fCov[6] = {0.};
   	 	fMyPrimVtx -> GetCovarianceMatrix(fCov);

    	lReconstructedVertexPresent = ( ( fMyPrimVtx->GetNContributors() > 0 ) && ( fCov[5] != 0 ) );  
    	if(!lReconstructedVertexPresent) { Printf("Bad vertex params"); fHistEventStat -> Fill("Bad vertex params",1); PostData(1,fList); return kFALSE; }

    	// before vertex cut |zvtx|<10 cm, for fMyPrimVtx only
    	fHistQA[centBin][0] -> Fill((fMyPrimVtx->GetX()));
    	fHistQA[centBin][1] -> Fill((fMyPrimVtx->GetY()));
    	fHistQA[centBin][2] -> Fill((fMyPrimVtx->GetZ()));

    	// count events having a proper vertex before vertex cut
    	fHistEventStat -> Fill("Proper vertex",1);
    	lVertexAcceptable = ( fMyPrimVtx->GetNContributors() < ntracks || (TMath::Abs(fMyPrimVtx->GetX()) < gVxMax && TMath::Abs(fMyPrimVtx->GetY()) < gVyMax && TMath::Abs(fMyPrimVtx->GetZ()) < gVzMax) );
    	if (!lVertexAcceptable) { Printf("Vertex out of range"); fHistEventStat -> Fill("Bad vertex position",1); PostData(1,fList); return kFALSE; }

	// do zvtx binning
  		Int_t vtx = GetZvtxBin(fMyPrimVtx->GetZ());
  		if (vtx < 0) return kFALSE;
        
    	// count events after vertex cut
    	fHistEventStat -> Fill("Proper vertex within 10cm",1);

    	// after vertex cut, for fMyPrimVtx only
    	fHistQA[centBin][3] -> Fill((fMyPrimVtx->GetX()));
    	fHistQA[centBin][4] -> Fill((fMyPrimVtx->GetY()));
    	fHistQA[centBin][5] -> Fill((fMyPrimVtx->GetZ()));
    	//___ end vertex
  	}
  	else
  	{
    	AliInfo(" Primary-vertex Selection: event REJECTED ...");
    	return kFALSE;
  	}
	}

	//_____ AOD MC , NOT needed...
    /*
  if (obj->InheritsFrom("AliAODMCHeader"))
  {
    if (TMath::Abs(((AliAODMCHeader*) obj)->GetVtxZ()) >= fVzMax)
    {
      AliInfo(" Primary-vertex Selection: event (based on MC) REJECTED ...");
      return kFALSE;
    }
  }
*/
	return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::CleanUp(TObjArray* tracks, TObject* mcObj)
{
  if (!tracks) return;
  
  if (fRemoveWeakDecays)	RemoveWeakDecays(tracks, mcObj);
  if (fRemoveDuplicates)	RemoveDuplicates(tracks);
}

//________________________________________________________________________
 //________________________________________________________________________
TObjArray* AliAnalysisTaskPidPidCorrelations::AcceptTracks(Int_t centBin, TObject* arrayMC, /*Bool_t onlyprimaries,*/ Bool_t useCuts)
{
  //_____ number of RECO AOD tracks
  Int_t nTracks = fMyAODEvent->GetNumberOfTracks();
//   AliDebug(5, Form("#tracks= %d %d", nTracks, centBin));
//   AliInfo(Form("There are %d tracks in this event",nTracks));
  
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted -> SetOwner(kTRUE);

  //_____ Track loop
  for (Int_t itrk=0; itrk < nTracks; itrk++)
  {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fMyAODEvent->GetTrack(itrk));
    if (!track) { AliInfo("AcceptTracks: Could not receive track"); continue; }

    /* FIXME
    // select only primary tracks
    if (onlyprimaries) { if (fMyAODEvent->GetType() != AliAODTrack::kPrimary) continue; }

    AliVParticle* vtrack = dynamic_cast<AliVParticle*>(fMyAODEvent->GetTrack(itrk));
    //____use only tracks from primaries
    if(onlyprimaries)
    {
      if (vtrack->GetLabel()<=0)	continue;
      if (!(static_cast<AliAODMCParticle*>(fMyMcArray->At(vtrack->GetLabel()))->IsPhysicalPrimary()))		continue;
    }
    delete vtrack; // vtrack needed just for correction for underestimation of strangeness in Monte Carlos, see below, FIXME later...
    */
    //______track selection cuts, hybrid track cuts

    Bool_t bGood = kFALSE;
    if		(fFilterType == 0)	bGood = kTRUE;
    else if	(fFilterType == 1)	bGood = track -> IsHybridTPCConstrainedGlobal();
    else if	(fFilterType == 2)	bGood = track -> IsHybridGlobalConstrainedGlobal();
//     if ( !(track -> TestFilterBit(fFilterBit)) ) { AliInfo(" === TestFilterBit: not passed..."); return 0x0; } 
//     if ( (fFilterBit>0) && ( (!track->TestFilterBit(fFilterBit) || (!bGood)) )) return 0;
    if ( !bGood )
    {
//       AliInfo(" ==== IsHybridGlobalConstrainedGlobal = kFALSE\n");
//       Printf("%s:%d Not matching filter %d/%d %d/%d",(Char_t*)__FILE__,__LINE__,itrk,fMyAODEvent->GetNumberOfTracks(),fFilterBit,track->GetFilterMap());	
      continue;
    }
    //if (fTrackStatus != 0 && !CheckTrack(track)) return 0; //clm
   
    Double_t dxy	  = 0.;
    Double_t dz	  = 0.;
    Double_t chi2ndf = 0.;
    Double_t nCrossedRowsTPC = 0.;
    Double_t ratioCrossedRowsOverFindableClustersTPC = 0.;

    dxy = track -> DCA();
    dz = track -> ZAtDCA();
    fHistQA[centBin][6] -> Fill(dxy);
    fHistQA[centBin][7] -> Fill(dz);
    chi2ndf = track -> Chi2perNDF();
    fHistQA[centBin][11] -> Fill(chi2ndf);  
    nCrossedRowsTPC = track -> GetTPCClusterInfo(2,1);
    fHistQA[centBin][12] -> Fill(nCrossedRowsTPC); 
    //Double_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  
    if (track -> GetTPCNclsF() > 0)
    {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC /1.0/track->GetTPCNclsF();
      fHistQA[centBin][13] -> Fill(ratioCrossedRowsOverFindableClustersTPC);
    }

    //_____ only charged
    if (track->Charge()==0) continue;

    if (useCuts)
    {
      if ( track->Eta() < fTrackEtaMin || track->Eta() > fTrackEtaMax || track->Pt() < fTrackPtMin || track->Pt() > fTrackPtMax )
	continue;
    }

    //______ QA
    fHistQA[centBin][8] -> Fill(track->Pt());
    fHistQA[centBin][9] -> Fill(track->Eta());
    fHistQA[centBin][10] -> Fill(track->Phi());	

    //!FIXME -- correction for underestimation of strangeness in Monte Carlos...
    //.....
	
    //_____ MC-RECO matched
    if (arrayMC)
    {
      Int_t label = track->GetLabel();
      if (label<0) continue;
      
      //___ check HIJING Bg
      Bool_t kIsFromHijingBg = kFALSE;
      kIsFromHijingBg = IsFromHijingBg(label);
      if (!kIsFromHijingBg) continue;
	//_____ redefine track as a matched MC-particle
	//track = dynamic_cast<AliAODTrack*>(arrayMC->At(TMath::Abs(label)));
	//if (!track) continue;
    }

    //______ PID 
    Int_t myPartID = -1;
    myPartID = GetParticleID((AliVTrack*)track,centBin,kTRUE);
//     Printf("=========================== myPartID = %d",myPartID);
    
    tracksAccepted -> Add(new AliPidPidCorrelationReducedTrack(myPartID,track->Eta(),track->Phi(),track->Pt(),track->Charge()));
    //tracksAccepted -> SetUniqueID(eventno * 100000 + TMath::Abs(part->GetLabel())); //______ !FIXME
    //tracksAccepted -> SetUniqueID(part->GetUniqueID()); //______ !FIXME
    //tracksAccepted -> Add(part);
  }//_____ track loop
  
  return tracksAccepted;
}

//________________________________________________________________________
//________________________________________________________________________
TObjArray* AliAnalysisTaskPidPidCorrelations::AcceptMcTracks(Int_t centBin, Bool_t onlyprimaries, Bool_t useCuts)
{
  //_____ number of MC AOD particles
  Int_t nMCPart = fMyMcArray->GetEntriesFast();
  //if(fVerbose>0) Printf("AOD tracks: %d", nTracks); //! FIXME

  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted -> SetOwner(kTRUE); //! FIXME
  
  for (Int_t iMC=0; iMC < nMCPart; iMC++)
  {
    AliAODMCParticle* partMC = dynamic_cast<AliAODMCParticle*>(fMyMcArray->At(iMC));
    if (!partMC) { AliWarning("AcceptMcTracks: Could not receive particle"); continue; }
    
    //_____ eventually only primaries
    if (onlyprimaries)
    {
      if (!partMC->IsPhysicalPrimary()) continue;
    }
    //_____ eventually only hadrons
    Int_t pdgCode = partMC->GetPdgCode();
    Bool_t isHadron = TMath::Abs(pdgCode)==211 || TMath::Abs(pdgCode)==2212 || TMath::Abs(pdgCode)==321;
    if (!isHadron) continue;
    //_____ eventually only charged
    if (!partMC->Charge()) continue;

    if (useCuts)
    {
      if ( partMC->Eta() < fTrackEtaMin || partMC->Eta() > fTrackEtaMax || partMC->Pt() < fTrackPtMin || partMC->Pt() > fTrackPtMax )
        continue;
    }

    //______ QA
    fHistQA[centBin][8] -> Fill(partMC->Pt());
    fHistQA[centBin][9] -> Fill(partMC->Eta());
    fHistQA[centBin][10] -> Fill(partMC->Phi());

    //______ PID 
    Int_t myPartID = fPartUndefined;
    myPartID = GetParticleIDMC((AliVTrack*)partMC,centBin,kTRUE);
//     Printf("myPartID = %d",myPartID);

    tracksAccepted -> Add(new AliPidPidCorrelationReducedTrack(myPartID,partMC->Eta(),partMC->Phi(),partMC->Pt(),partMC->Charge()));
    //tracksAccepted -> SetUniqueID(eventno * 100000 + TMath::Abs(part->GetLabel())); //______ !FIXME
    //tracksAccepted -> SetUniqueID(part->GetUniqueID()); //______ !FIXME
    //tracksAccepted -> Add(part);

  }//_____ track loop

  return tracksAccepted;
}

//________________________________________________________________________
//________________________________________________________________________
TObjArray* AliAnalysisTaskPidPidCorrelations::AcceptMcRecoMachedTracks(Int_t centBin, Bool_t onlyprimaries, Bool_t useCuts)
{
  //_____ number of RECO AOD tracks
  Int_t nTracks = fMyAODEvent->GetNumberOfTracks();
  //if(fVerbose>0) Printf("AOD tracks: %d", nTracks); //! FIXME

  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted -> SetOwner(kFALSE); //! FIXME

  //fMyPrimVtx = dynamic_cast<AliAODVertex*>(fMyAODEvent->GetPrimaryVertexSPD());//GetPrimaryVertex()
  //Double_t vzAOD = fMyPrimVtx->GetZ();

  //______ track loop
  for (Int_t itrk=0; itrk < nTracks; itrk++)
  {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fMyAODEvent->GetTrack(itrk));
    if (!track) { AliWarning("AcceptMcRecoMachedTracks: Could not receive track"); continue; }

    if (track->GetLabel()<0) continue;

    AliVParticle* vtrack = dynamic_cast<AliVParticle*>(fMyMcArray->At(track->GetLabel()));
    if (!vtrack) { AliWarning("AcceptMcRecoMachedTracks: Could not receive vtrack"); continue; }

    //use only tracks from primaries
    if (onlyprimaries && !((static_cast<AliAODMCParticle*>(vtrack))->IsPhysicalPrimary())) continue; // !FIXME gives seg fault...

    //______track selection cuts, hybrid track cuts
    Bool_t bGood = kFALSE;
    if	(fFilterType == 0)	bGood = kTRUE;
    else if	(fFilterType == 1)	bGood = track -> IsHybridTPCConstrainedGlobal();
    else if	(fFilterType == 2)	bGood = track -> IsHybridGlobalConstrainedGlobal();
    //if ( !(track -> TestFilterBit(fFilterBit)) ) return 0; 
    //if ( (fFilterBit>0) && ((!track) -> TestFilterBit(fFilterBit) || (!bGood))) )
    if ( !bGood )
    {
      //AliInfo(" ==== IsHybridGlobalConstrainedGlobal = kFALSE\n");
      //Printf("%s:%d Not matching filter %d/%d %d/%d",(Char_t*)__FILE__,__LINE__,it,fMyAODEvent->GetNumberOfTracks(),fFilterBit,track->GetFilterMap());	
      continue;
    }
    
//     if (fTrackStatus != 0 && !CheckTrack(track)) return 0; //clm
/*
    if(vtrack->GetLabel()<0)
	continue;
    else
    {
      AliAODMCParticle* partOfTrack = dynamic_cast<AliAODMCParticle*>(fMyMcArray->At(vtrack->GetLabel()));
      if(!partOfTrack) Printf("label=%d", vtrack->GetLabel());
      delete partOfTrack;delete vtrack;

      //correction for underestimation of strangeness in Monte Carlos !FIXME
      //...
    }
*/
    //_____ only charged
    if (track->Charge()==0) continue;

    if (useCuts)
    {
      if ( track->Eta() < fTrackEtaMin || track->Eta() > fTrackEtaMax || track->Pt() < fTrackPtMin || track->Pt() > fTrackPtMax )
	continue;
    }

    //______ QA
    fHistQA[centBin][8] -> Fill(track->Pt());
    fHistQA[centBin][9] -> Fill(track->Eta());		
    fHistQA[centBin][10] -> Fill(track->Phi());

    //______ PID 
    Int_t myPartID = fPartUndefined;
    myPartID = GetParticleID(dynamic_cast<AliVTrack*>(track),centBin,kTRUE);
    //Printf("myPartID = %d",myPartID);

    tracksAccepted -> Add(new AliPidPidCorrelationReducedTrack(myPartID,track->Eta(),track->Phi(),track->Pt(),track->Charge()));
    //tracksAccepted -> SetUniqueID(eventno * 100000 + TMath::Abs(part->GetLabel())); //______ !FIXME
    //tracksAccepted -> SetUniqueID(part->GetUniqueID()); //______ !FIXME
    //tracksAccepted -> Add(part);
  }
  
  //Double_t vzMC = fMyMcHeader->GetVtxZ();
  return tracksAccepted;
}

//________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::UserCreateOutputObjects()
{
//   fList = new THashList;
  fList = new TList();
  fList -> SetOwner(kTRUE);
  fList -> SetName(GetOutputListName());

  fPIDResponse = (dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager() -> GetInputEventHandler())) -> GetPIDResponse(); 
  if( !fPIDResponse ) { AliError("Cannot get the PID response"); return; }
  
  AddSettingsTree();

  const Char_t* kMyPIDTypeName[] = {"ITS","TPC","TOF","HMPID"} ;
  const Char_t* kMyParticleSpeciesName[] = { "Electrons","Muons","Pions","Kaons","Protons","Undefined" } ;


  fHistNev = new TH1I("fHistNev","Nev.",1,0,1); fHistNev -> Sumw2();
  fList -> Add(fHistNev); 
  fHistTriggerStats = new TH1F("fHistTriggerStats","Trigger statistics;TriggerBit;N_{events}",130,0,130);
  fList -> Add(fHistTriggerStats);
  fHistTriggerRun = new TH1F("fHistTriggerRun","Trigger;trigger",10,-0.5,9.5);
  fList -> Add(fHistTriggerRun);

  const Int_t nEventStatBins = 6;
  TString gCutName[nEventStatBins] = {"Total","No trigger","Wrong centrality","Bad vertex params","Bad vertex position","Analyzed"};
  fHistEventStat = new TH1F("fHistEventStat","Event statistics;;N_{events}",nEventStatBins,-0.5,nEventStatBins-0.5);
  for(Int_t i = 1; i <= nEventStatBins; i++) fHistEventStat -> GetXaxis() -> SetBinLabel(i,gCutName[i-1].Data());
  fList -> Add(fHistEventStat);
  
  TString gCentName[9] = {"V0M","FMD","TRK","TKL","CL0","CL1","V0MvsFMD","TKLvsV0M","ZEMvsZDC"};
  fHistCentStats = new TH2F("fHistCentStats","Centrality statistics;;Cent percentile",9,-0.5,8.5,110,-5,105);
  for(Int_t i = 1; i <= 9; i++) fHistCentStats->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  fList -> Add(fHistCentStats);
  
  fHistCentralityPercentile = new TH1F("fHistCentralityPercentile","Centrality Percentile;Centrality;Entries",101,-1,100);
  fList -> Add(fHistCentralityPercentile);
  fHistCentralityClass10 = new TH1F("fHistCentralityClass10","Centrality Class 10;Centrality;Entries",10,-0.5,9.5);
  fList -> Add(fHistCentralityClass10);
  fHistCentralityClass5 = new TH1F("fHistCentralityClass5","Centrality Class 5;Centrality;Entries",20,-0.5,19.5);
  fList -> Add(fHistCentralityClass5);  

  fHistV0M = new TH2F("fHistV0M","V0 Multiplicity C vs. A",5, 0, 20000, 5, 0, 20000);
  fList -> Add(fHistV0M);
  
  TString gRefTrackName[6] = {"tracks","tracksPos","tracksNeg","tracksTPConly","clusITS0","clusITS1"};
  fHistRefTracks  = new TH2F("fHistRefTracks","Nr of Ref tracks/event vs. ref track estimator;;Nr of tracks",6, 0, 6, 20, 0, 20000);
  for(Int_t i = 1; i <= 6; i++) fHistRefTracks->GetXaxis()->SetBinLabel(i,gRefTrackName[i-1].Data());
  fList -> Add(fHistRefTracks);
  
  //______ MC QA

  if (fUseMC)
  {
  
  fHistMcGenerator = new TH1F("fHistMcGenerator","All MC generators;;Entries",20,0,20);
  fHistMcGenerator->GetXaxis()->SetBinLabel(1,"Hijing");
  fHistMcGenerator->GetXaxis()->SetBinLabel(2,"Pythia");
  fHistMcGenerator->GetXaxis()->SetBinLabel(3,"Unknown");
  fHistMcGenerator->Sumw2();
  fList->Add(fHistMcGenerator);

  fHist_HijingBg = new TH1F("fHist_HijingBg","MC generator;;Entries",2,0,2);
  fHist_HijingBg->GetXaxis()->SetBinLabel(1,"Hijing BG");
  fHist_HijingBg->GetXaxis()->SetBinLabel(2,"No Hijing BG");
  fHist_HijingBg->Sumw2();
  fList->Add(fHist_HijingBg);  
  
  fHistNumOfPartPerEvt = new TH1F("fHistNumOfPartPerEvt","Number of part per event;;Entries",1,0,1);
  fHistNumOfPartPerEvt->GetXaxis()->SetBinLabel(1,"# particles per event");
  fHistNumOfPartPerEvt->Sumw2();
  fList->Add(fHistNumOfPartPerEvt);

  fHistMcStats = new TH1F("fHistMcStats","MC stats;;Entries",6,0,6);
  fHistMcStats->Sumw2();
  fHistMcStats->GetXaxis() -> SetBinLabel(1,"MC all");
  fList->Add(fHistMcStats);
  
  fHistMcAllPt = new TH1F("fHistMcAllPt","MC All;pT (GeV/c);Entries",1000,0,100);
  fHistMcAllPt->Sumw2();
  fList->Add(fHistMcAllPt);
  fHistMcAllPt_Hijing= new TH1F("fHistMcAllPt_Hijing","Hijing;pT (GeV/c);Entries",1000,0,100);
  fHistMcAllPt_Hijing->Sumw2();
  fList->Add(fHistMcAllPt_Hijing);

  fHistMcAllPt_Dec= new TH1F("fHistMcAllPt_Dec","Decay;pT (GeV/c);Entries",1000,0,100);
  fHistMcAllPt_Dec->Sumw2();
  fList->Add(fHistMcAllPt_Dec);
  fHistMcAllPt_Inj= new TH1F("fHistMcAllPt_Inj","Injected;pT (GeV/c);Entries",1000,0,100);
  fHistMcAllPt_Inj->Sumw2();
  fList->Add(fHistMcAllPt_Inj);

  fHistMcAllEta_NotHijing= new TH1F("fHistMcAllEta_NotHijing","not Hijing;#eta;Entries",2000,-10,10);
  fHistMcAllEta_NotHijing->Sumw2();
  fList->Add(fHistMcAllEta_NotHijing);
  fHistMcAllEta_Hijing= new TH1F("fHistMcAllEta_Hijing","Hijing;#eta;Entries",2000,-10,10);
  fHistMcAllEta_Hijing->Sumw2();
  fList->Add(fHistMcAllEta_Hijing);
  fHistMcAllEta= new TH1F("fHistMcAllEta","All;#eta;Entries",2000,-10,10);
  fHistMcAllEta->Sumw2();
  fList->Add(fHistMcAllEta);

  fList->Add(new TH1F("pids", ";pdg;tracks", 6001, -3000.5, 3000.5));

  }
  
  for (Int_t iCent=0; iCent < fNbinsCent; iCent++)
  {
    fHistQA[iCent][0] = new TH1F(Form("fHistQAvx_Cent%02d",iCent), Form("Histo Vx All, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -5., 5.);
    fHistQA[iCent][1] = new TH1F(Form("fHistQAvy_Cent%02d",iCent), Form("Histo Vy All, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -5., 5.);
    fHistQA[iCent][2] = new TH1F(Form("fHistQAvz_Cent%02d",iCent), Form("Histo Vz All, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -25., 25.);  
    fHistQA[iCent][3] = new TH1F(Form("fHistQAvxA_Cent%02d",iCent), Form("Histo Vx  After Cut, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -5., 5.);
    fHistQA[iCent][4] = new TH1F(Form("fHistQAvyA_Cent%02d",iCent), Form("Histo Vy After Cut, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -5., 5.);
    fHistQA[iCent][5] = new TH1F(Form("fHistQAvzA_Cent%02d",iCent), Form("Histo Vz After Cut, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -25., 25.);
    fHistQA[iCent][6] = new TH1F(Form("fHistQADcaXyC_Cent%02d",iCent), Form("Histo DCAxy after cut, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -5., 5.);
    fHistQA[iCent][7] = new TH1F(Form("fHistQADcaZC_Cent%02d",iCent), Form("Histo DCAz after cut, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))), 50, -5., 5.);   
    fHistQA[iCent][8] = new TH1F(Form("fHistQAPt_Cent%02d",iCent),Form("p_{T} distribution, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))),900,0.,9.);
    fHistQA[iCent][9] = new TH1F(Form("fHistQAEta_Cent%02d",iCent),Form("#eta distribution, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))),360,-1.8,1.8);
    fHistQA[iCent][10] = new TH1F(Form("fHistQAPhi_Cent%02d",iCent),Form("#phi distribution, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))),340,0,6.8);
    fHistQA[iCent][11] = new TH1F(Form("fHistQAChi2_Cent%02d",iCent),Form("Chi2 per NDF, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))),100,0,10);
    fHistQA[iCent][12] = new TH1F(Form("fHistQAnCrossedRowsTPC_Cent%02d",iCent),Form("Number of TPC ccrossed rows, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))),200,0,200);
    fHistQA[iCent][13] = new TH1F(Form("fHistQAratioCrossedRowsOverFindableClustersTPC_Cent%02d",iCent),Form("Number of TPC ccrossed rows find clusters, CentBin %d-%d",Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),Int_t(fCentAxis -> GetBinUpEdge(iCent+1))),200,0,2);
    for(Int_t i=0; i<14; i++) fList -> Add(fHistQA[iCent][i]);
    
    for (Int_t i=0;i<6;i++)
    {
      fHistRefTracksCent[iCent][i] = new TH2F(Form("fHistRefTracksCent_%02d_%d",iCent,i),
					      Form("Nr of tracksNr of Ref tracks/event in ref track estimator %s for CentBin %d-%d;centrality [%%];Nr of Ref tracks/event",
						   gRefTrackName[i-1].Data(),
						   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
						   Int_t(fCentAxis->GetBinUpEdge(iCent+1))),10.*(fCentralityPercentileMax-fCentralityPercentileMin),fCentralityPercentileMin,fCentralityPercentileMax,20,0,20000);
      fList->Add(fHistRefTracksCent[iCent][i]);
    }
  
/*    
    for (Int_t iPtBinTrig=0; iPtBinTrig < fNbinsPt; iPtBinTrig++)
    {
      for (Int_t iPtBinAssoc=0; iPtBinAssoc < fNbinsPt; iPtBinAssoc++)
      {
	fHistDeltaPhi[iCent][iPtBinTrig][iPtBinAssoc] = new TH1F(Form("fHistDeltaPhi_Cent%02d_PtBin%02d_%02d",iCent,iPtBinTrig,iPtBinAssoc), 
							       Form("%d-%d %%, %3.1f<p_{T,trig}<%3.1f, %3.1f<p_{T,assoc}<%3.1f",
								    Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),
								    Int_t(fCentAxis -> GetBinUpEdge(iCent+1)),
								    fPtAxis -> GetBinLowEdge(iPtBinTrig+1),
								    fPtAxis -> GetBinUpEdge(iPtBinTrig+1),
								    fPtAxis -> GetBinLowEdge(iPtBinAssoc+1),
								    fPtAxis -> GetBinUpEdge(iPtBinAssoc+1)),
								    72, -0.5*TMath::Pi(), 1.5*TMath::Pi());
	
	fHistDeltaPhiMix[iCent][iPtBinTrig][iPtBinAssoc] = new TH1F(Form("fHistDeltaPhiMix_Cent%02d_PtBin%02d_%02d",iCent,iPtBinTrig,iPtBinAssoc), 
								  Form("%d-%d %%, %3.1f<p_{T,trig}<%3.1f, %3.1f<p_{T,assoc}<%3.1f MIXED",
								       Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis -> GetBinUpEdge(iCent+1)),
								       fPtAxis -> GetBinLowEdge(iPtBinTrig+1),
								       fPtAxis -> GetBinUpEdge(iPtBinTrig+1),
								       fPtAxis -> GetBinLowEdge(iPtBinAssoc+1),
								       fPtAxis -> GetBinUpEdge(iPtBinAssoc+1)),
									72, -0.5*TMath::Pi(), 1.5*TMath::Pi());
	
	fHistDeltaPhi[iCent][iPtBinTrig][iPtBinAssoc]    -> SetXTitle("#Delta#varphi [rad]");
	fHistDeltaPhiMix[iCent][iPtBinTrig][iPtBinAssoc] -> SetXTitle("#Delta#varphi [rad]");

	fHistDeltaPhi[iCent][iPtBinTrig][iPtBinAssoc]    -> Sumw2();
	fHistDeltaPhiMix[iCent][iPtBinTrig][iPtBinAssoc] -> Sumw2();
	
	fList -> Add(fHistDeltaPhi[iCent][iPtBinTrig][iPtBinAssoc]);
	fList -> Add(fHistDeltaPhiMix[iCent][iPtBinTrig][iPtBinAssoc]);

	fHistDphiDeta[iCent][iPtBinTrig][iPtBinAssoc] = new TH2F(Form("fHistDphiDeta_Cent%02d_PtBin%02d_%02d",iCent,iPtBinTrig,iPtBinAssoc), 
								  Form("%d-%d %%, %3.1f<p_{T,trig}<%3.1f, %3.1f<p_{T,assoc}<%3.1f",
								       Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis -> GetBinUpEdge(iCent+1)),
								       fPtAxis -> GetBinLowEdge(iPtBinTrig+1),
								       fPtAxis -> GetBinUpEdge(iPtBinTrig+1),
								       fPtAxis -> GetBinLowEdge(iPtBinAssoc+1),
								       fPtAxis -> GetBinUpEdge(iPtBinAssoc+1)),
									72, -0.5*TMath::Pi(), 1.5*TMath::Pi(),
									36, -1.8, 1.8);
	
	fHistDphiDetaMix[iCent][iPtBinTrig][iPtBinAssoc] = new TH2F(Form("fHistDphiDetaMix_Cent%02d_PtBin%02d_%02d",iCent,iPtBinTrig,iPtBinAssoc), 
								     Form("%d-%d %%, %3.1f<p_{T,trig}<%3.1f, %3.1f<p_{T,assoc}<%3.1f MIXED",
									  Int_t(fCentAxis -> GetBinLowEdge(iCent+1)),
									  Int_t(fCentAxis -> GetBinUpEdge(iCent+1)),
									  fPtAxis -> GetBinLowEdge(iPtBinTrig+1),
									  fPtAxis -> GetBinUpEdge(iPtBinTrig+1),
									  fPtAxis -> GetBinLowEdge(iPtBinAssoc+1),
									  fPtAxis -> GetBinUpEdge(iPtBinAssoc+1)),
									  72, -0.5*TMath::Pi(), 1.5*TMath::Pi(),
									  36, -1.8, 1.8);
	
	fHistDphiDeta[iCent][iPtBinTrig][iPtBinAssoc]    -> SetXTitle("#Delta#varphi [rad]");
	fHistDphiDetaMix[iCent][iPtBinTrig][iPtBinAssoc] -> SetXTitle("#Delta#varphi [rad]");
	fHistDphiDeta[iCent][iPtBinTrig][iPtBinAssoc]    -> SetYTitle("#Delta#eta");
	fHistDphiDetaMix[iCent][iPtBinTrig][iPtBinAssoc] -> SetYTitle("#Delta#eta");

	fHistDphiDeta[iCent][iPtBinTrig][iPtBinAssoc]    -> Sumw2();
	fHistDphiDetaMix[iCent][iPtBinTrig][iPtBinAssoc] -> Sumw2();
	
	fList -> Add(fHistDphiDeta[iCent][iPtBinTrig][iPtBinAssoc]);
	fList -> Add(fHistDphiDetaMix[iCent][iPtBinTrig][iPtBinAssoc]);
      }
    }*/

    fHistTracksEtaTrigVsEtaAssoc[iCent]		= new TH2F(Form("fHistTracksEtaTrigVsEtaAssoc_%02d",iCent),"#eta_{trig} vs #eta_{assoc}",20,-1.,1.,20,-1.,1.);
    fList -> Add(fHistTracksEtaTrigVsEtaAssoc[iCent]);
    fHistTracksEtaTrigVsEtaAssocMixed[iCent]	= new TH2F(Form("fHistTracksEtaTrigVsEtaAssocMixed_%02d",iCent),"#eta_{trig} vs #eta_{assoc}",20,-1.,1.,20,-1.,1.);
    fList -> Add(fHistTracksEtaTrigVsEtaAssocMixed[iCent]);
    
    fHistSingleHadronsPt[iCent]		= new TH1F(Form("fHistSingleHadronsPt_Cent%02d",iCent),"p_{T}", fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray());
    fList -> Add(fHistSingleHadronsPt[iCent]);
//     fHistSingleHadronsPt_Mixed[iCent]	= new TH1F(Form("fHistSingleHadronsPt_Mixed_Cent%02d",iCent), "p_{T}", fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray());
//     fList -> Add(fHistSingleHadronsPt_Mixed[iCent]);

    fHistSingleHadronsEtaPt[iCent]	= new TH2F(Form("fHistSingleHadronsEtaPt_Cent%02d",iCent), "#eta vs p_{T}",fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray(),fEtaAxis->GetNbins(),(Double_t*)fEtaAxis->GetXbins()->GetArray());
    fList -> Add(fHistSingleHadronsEtaPt[iCent]);
//     fHistSingleHadronsEtaPt_Mixed[iCent]	= new TH2F(Form("fHistSingleHadronsEtaPt_Mixed_Cent%02d",iCent), "#eta vs p_{T}",fNbinsPt, (Double_t*)fPtAxis->GetXbins()->GetArray(),fEtaAxis->GetNbins(),(Double_t*)fEtaAxis->GetXbins()->GetArray());
//     fList -> Add(fHistSingleHadronsEtaPt_Mixed[iCent]);
    
    
    fHistTPCdEdx[iCent] = new TH2F(Form("fHistTPCdEdx_%d",iCent), ";p_{TPCin} (GeV/c);dE/dx (au.)", 5, 0., 5., 50, 0., 500.);
    fList -> Add(fHistTPCdEdx[iCent]);
    fHistTOFbeta[iCent] = new TH2F(Form("fHistTOFbeta_%d",iCent), ";p (GeV/c);v/c", 5, 0., 5., 50, 0.1, 1.1);
    fList -> Add(fHistTOFbeta[iCent]);

    fHistTPCdEdx_selected[iCent] = new TH2F(Form("fHistTPCdEdx_selected_%d",iCent), ";p_{TPCin} (GeV/c);dE/dx (au.)",  5, 0., 5., 50, 0., 500.);
    fList -> Add(fHistTPCdEdx_selected[iCent]);
    fHistTOFbeta_selected[iCent] = new TH2F(Form("fHistTOFbeta_selected_%d",iCent), ";p (GeV/c);v/c",  5, 0., 5., 50, 0.1, 1.1);
    fList -> Add(fHistTOFbeta_selected[iCent]);

    for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
    {
      fHistNSigmaTPC[iCent][iSpecies] = new TH2F(Form("fHistNSigmaTPC_%d_%s",iCent,AliPID::ParticleName(iSpecies)), Form(";p_{T} (GeV/c); N_{#sigma-TPC}^{%s}", AliPID::ParticleLatexName(iSpecies)), 5, 0., 5., 20, -10., 10.);
      fList -> Add(fHistNSigmaTPC[iCent][iSpecies]);
      fHistNSigmaTOF[iCent][iSpecies] = new TH2F(Form("fHistNSigmaTOF_%d_%s",iCent,AliPID::ParticleName(iSpecies)), Form(";p_{T} (GeV/c); N_{#sigma-TOF}^{%s}", AliPID::ParticleLatexName(iSpecies)), 5, 0., 5., 20, -10., 10.);
      fList -> Add(fHistNSigmaTOF[iCent][iSpecies]);
    }

    
    Double_t ptmin[] = { 0.,0.,0.,0. };
    Double_t ptmax[] = { 2.,3.,5.,6. };
    
    //_____ nsigma plot
    for(Int_t ipart=0;ipart < AliPID::kSPECIES;ipart++)
    {
      for(Int_t ipid=0;ipid < kMyNSigmaPIDType;ipid++)
      {
	Double_t miny = -10.;
	Double_t maxy =  10.;
	
        //if (ipid == kMyNSigmaTPCTOF) { miny=0; maxy=50; }
        fHistoNSigma[iCent] = new TH2F(Form("NSigma_Cent%02d_%d_%d",iCent,ipart,ipid),
				       Form("n#sigma %d-%d %% %s %s",Int_t(fCentAxis->GetBinLowEdge(iCent+1)),Int_t(fCentAxis->GetBinUpEdge(iCent+1)),kMyParticleSpeciesName[ipart],kMyPIDTypeName[ipid]),
				       10.*(ptmax[ipid]-ptmin[ipid]),ptmin[ipid],ptmax[ipid],2000,miny,maxy);
        fHistoNSigma[iCent] -> GetXaxis() -> SetTitle("p_{T} (GeV / c)");
        fHistoNSigma[iCent] -> GetYaxis() -> SetTitle(Form("n#sigma %s %s",kMyParticleSpeciesName[ipart],kMyPIDTypeName[ipid]));
        fList -> Add(fHistoNSigma[iCent]);
      }
    }
    
    if (fUseMC)
    {
      //_____ nsigmaRec plot
      for(Int_t ipart=0;ipart < AliPID::kSPECIES;ipart++)
      {
	for(Int_t ipid=0;ipid < kMyNSigmaPIDType;ipid++)
	{
	  Double_t miny	=-10.;
          Double_t maxy	= 10.;
          //if(ipid == kMyNSigmaTPCTOF) { miny=0; maxy=20; }
        fHistoNSigma[iCent] = new TH2F(Form("NSigmaRec_Cent%02d_%d_%d",iCent,ipart,ipid),
				       Form("n#sigma reconstructed %d-%d %% %s %s",Int_t(fCentAxis->GetBinLowEdge(iCent+1)),Int_t(fCentAxis->GetBinUpEdge(iCent+1)),kMyParticleSpeciesName[ipart],kMyPIDTypeName[ipid]),
				       10.*(ptmax[ipid]-ptmin[ipid]),ptmin[ipid],ptmax[ipid],2000,miny,maxy);
        fHistoNSigma[iCent] -> GetXaxis() -> SetTitle("p_{T} (GeV / c)");
        fHistoNSigma[iCent] -> GetYaxis() -> SetTitle(Form("n#sigma %s %s",kMyParticleSpeciesName[ipart],kMyPIDTypeName[ipid]));
        fList -> Add(fHistoNSigma[iCent]);
	}
      }
      
      //_____ nsigmaMC plot
      for(Int_t ipart=0;ipart<fPartNSpeciesQA;ipart++)
      {
	for(Int_t ipid=0;ipid < kMyNSigmaPIDType;ipid++)
	{
	  Double_t miny	= -10.;
	  Double_t maxy	=  10.;
	  //if(ipid == kMyNSigmaTPCTOF) { miny=0; maxy=50; }
        fHistoNSigma[iCent] = new TH2F(Form("NSigmaMC_Cent%02d_%d_%d",iCent,ipart,ipid),
				       Form("n#sigma %d-%d %% %s %s",Int_t(fCentAxis->GetBinLowEdge(iCent+1)),Int_t(fCentAxis->GetBinUpEdge(iCent+1)),kMyParticleSpeciesName[ipart],kMyPIDTypeName[ipid]),
				       10.*(ptmax[ipid]-ptmin[ipid]),ptmin[ipid],ptmax[ipid],2000,miny,maxy);
        fHistoNSigma[iCent] -> GetXaxis() -> SetTitle("p_{T} (GeV / c)");
        fHistoNSigma[iCent] -> GetYaxis() -> SetTitle(Form("n#sigma %s %s",kMyParticleSpeciesName[ipart],kMyPIDTypeName[ipid]));
        fList -> Add(fHistoNSigma[iCent]);
	}
      }    
    }// MC
  }// iCent
         

  fHistControlConvResoncances = new TH2F("fHistControlConvResoncances", ";id;delta mass", 3, -0.5, 2.5, 20, -0.1, 0.1);
  fList -> Add(fHistControlConvResoncances);

/*	
  Int_t binsH[4]   = { fNbinsZvtx		,  fNbinsCent			, fPoolSize	,  fMinNumTrack	};
  Double_t minH[4] = { -10		,  fCentralityPercentileMin	, 0	 	,  0		};
  Double_t maxH[4] = { 10  	 	,  fCentralityPercentileMax  	, fPoolSize	,  fMinNumTrack	};
  
  fHistPoolMgrQA = new THnSparseD("fHistPoolMgrQA","vtxz:cent:Nevents:Ntracks",4,binsH,minH,maxH);
  fHistPoolMgrQA->Sumw2();
  fList -> Add(fHistPoolMgrQA);*/

	//_______ FillCorrelation _______

//   Int_t anaSteps = 3;       // analysis steps

	// single particle histograms
//   Int_t iBinSingle[kTrackVariablesSingle];        // binning for track variables
//   Double_t* dBinsSingle[kTrackVariablesSingle];   // bins for track variables  
//   TString axisTitleSingle[kTrackVariablesSingle]; // axis titles for track variables
  
  // two particle histograms
  Int_t iBinPair[kTrackVariablesPair];         // binning for track variables
  Double_t* dBinsPair[kTrackVariablesPair];    // bins for track variables  
//   TString axisTitlePair[kTrackVariablesPair];  // axis titles for track variables

  //_____ dim 1,2 -- trig, assoc PID bins
//   const Int_t kNpidBins = 6;
//   Double_t pidBins[kNpidBins+1] = {1.,2.,3.,4.,5.,6.,7.};
//   iBinSingle[0]    		= kNpidBins;
//   dBinsSingle[0]   		= pidBins;
//   axisTitleSingle[0]  = "PID_{trig.}"; 
//   iBinPair[0]      		= kNpidBins;
//   dBinsPair[0]     		= pidBins;
//   axisTitlePair[0]    = "PID_{trig.}"; 
//   iBinPair[1]      		= kNpidBins;
//   dBinsPair[1]     		= pidBins;
//   axisTitlePair[1]    = "PID_{assoc.}"; 
  //_____ dim 3,4 -- pT,trig; pT,assoc bins
  //const Int_t kNPtBins = 16;
  //Double_t ptBins[kNPtBins+1] = {0.2,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.,12.,15.,20.};
//   iBinSingle[1]     = fNbinsPt;//kNPtBins;
//   dBinsSingle[1]    = (Double_t*)fPtAxis->GetXbins()->GetArray();//ptBins;    
//   axisTitleSingle[1]  = "p_{T,trig.} (GeV/c)";
  iBinPair[0]       = fNbinsPt;//kNPtBins;
  dBinsPair[0]      = (Double_t*)fPtAxis->GetXbins()->GetArray();//ptBins;
//   axisTitlePair[2]    = "p_{T,trig.} (GeV/c)"; 
  iBinPair[1]       = fNbinsPt;//kNPtBins;
  dBinsPair[1]      = (Double_t*)fPtAxis->GetXbins()->GetArray();//ptBins;
//   axisTitlePair[3]    = "p_{T,assoc.} (GeV/c)"; 
  //_____ dim 5 -- dphi bins
  const Int_t kNDeltaPhiBins = 72;
  Double_t deltaPhiBins[kNDeltaPhiBins+1];
  for(Int_t i = 0; i < kNDeltaPhiBins+1; i++) { deltaPhiBins[i] = -TMath::Pi()/2. + i * 5.*TMath::Pi()/180.; } 
  iBinPair[2]       = kNDeltaPhiBins;
  dBinsPair[2]      = deltaPhiBins;
//   axisTitlePair[4]  = "#Delta#varphi (rad)"; 
  //_____ dim 6 -- deta bins
  const Int_t kNDeltaEtaBins = 80;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for(Int_t i=0; i < kNDeltaEtaBins+1; i++) { deltaEtaBins[i] = - fDeltaEtaMax + i * 2 * fDeltaEtaMax / (Double_t)kNDeltaEtaBins; } 
  iBinPair[3]       = kNDeltaEtaBins;
  dBinsPair[3]      = deltaEtaBins;
//   axisTitlePair[5]  	= "#Delta#eta";
  //_____ dim 7 -- vertex z bins !! same as in AliEventPoolManager()
  //const Int_t kNVertexZBins  = 10;
  //Double_t vertexZBins[kNZvtxBins+1] = { -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10 };
  //const Int_t kNVertexZBins = 9;
  //Double_t vertexZBins[kNVertexZBins+1] = {-10., -7., -5., -3., -1., 1., 3., 5., 7., 10.};
//   iBinSingle[2]    = fNbinsZvtx;//kNVertexZBins;
//   dBinsSingle[2]   = (Double_t*)fZvtxAxis->GetXbins()->GetArray();//vertexZBins;
  iBinPair[4]      = fNbinsZvtx;//kNVertexZBins;
  dBinsPair[4]     = (Double_t*)fZvtxAxis->GetXbins()->GetArray();//vertexZBins;
//   axisTitleSingle[2]  = "v_{z} (cm)"; 
//   axisTitlePair[6]    = "v_{z} (cm)";
  //_____ dim 8 -- centrality bins !! same as in AliEventPoolManager()
  //const Int_t kNCentralityBins = 9;
  //Double_t centralityBins[kNCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
//   iBinSingle[3]  = fNbinsCent;//kNCentralityBins;
//   dBinsSingle[3] = (Double_t*)fCentAxis->GetXbins()->GetArray();//centralityBins;	  
//   axisTitleSingle[3]  = "centrality [#%]";
//   iBinPair[7]    = fNbinsCent;//kNCentralityBins;
//   dBinsPair[7]   = (Double_t*)fCentAxis->GetXbins()->GetArray();//centralityBins;
//   axisTitlePair[7]  	= "centrality [#%]"; 

/*
	TString histName;
  histName = "fHistCorrSingle"; 
  if (fCentralityEstimator) histName += fCentralityEstimator.Data();
  const TString title = histName+";#part;#p_{T,trig};vertex Z (bin);centr (bin);";
  fHistCorrSingle = new AliTHn(histName.Data(),title.Data(),anaSteps,kTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<kTrackVariablesSingle; j++)
	{
    fHistCorrSingle -> SetBinLimits(j, dBinsSingle[j]);
    fHistCorrSingle -> SetVarTitle(j, axisTitleSingle[j]);
  }
  fList->Add(fHistCorrSingle);
  //______ THn for triggers and assoc, pair histo SAME evt
  histName = "fHistCorrPairSame";
  if (fCentralityEstimator) histName += fCentralityEstimator.Data();
  fHistCorrPairSame = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++)
	{
    fHistCorrPairSame -> SetBinLimits(j, dBinsPair[j]);
    fHistCorrPairSame -> SetVarTitle(j, axisTitlePair[j]);
  }
  fList->Add(fHistCorrPairSame);
  //______ THn for triggers and assoc, pair histo MIXED evt
  histName = "fHistCorrPairMixed";
  if (fCentralityEstimator) histName += fCentralityEstimator.Data();
  fHistCorrPairMixed = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++)
	{
    fHistCorrPairMixed -> SetBinLimits(j, dBinsPair[j]);
    fHistCorrPairMixed -> SetVarTitle(j, axisTitlePair[j]);
  }
  fList->Add(fHistCorrPairMixed);

  AliInfo("Finished setting up the AliTHn");*/

  //_______ THnSparse histos

  Int_t ptMin 	= fPtAxis->GetXmin();
  Int_t ptMax 	= fPtAxis->GetXmax();
  Int_t zvtxMin	= fZvtxAxis->GetXmin();
  Int_t zvtxMax	= fZvtxAxis->GetXmax();
//   Int_t centMin	= fCentAxis->GetXmin();
//   Int_t centMax	= fCentAxis->GetXmax();

  //_____ single histo
/*
  Int_t binsHisto4[kTrackVariablesSingle]   = { iBinSingle[0]	,  iBinSingle[1]	, iBinSingle[2]	 ,  iBinSingle[3]	, 2};
  Double_t minHisto4[kTrackVariablesSingle] = { 1		 	,  ptMin	  	, zvtxMin	 ,  centMin	  	, 7};
  Double_t maxHisto4[kTrackVariablesSingle] = { 8  	 	,  ptMax  	, zvtxMax	 ,  centMax		, 9};
*/
  //_____ pair histo
//   Int_t binsHistoPair[kTrackVariablesPair]   = { iBinPair[0]	, iBinPair[1]	, iBinPair[2]	, iBinPair[3]	,	iBinPair[4]	,	iBinPair[5]	,	iBinPair[6]	, iBinPair[7]	, 2,2};
//   Double_t minHistoPair[kTrackVariablesPair] = { 1	 	, 1		, ptMin		, ptMin	 	,	-0.5*TMath::Pi()	,	-fDeltaEtaMax	,	zvtxMin		, centMin	, 7,7};
//   Double_t maxHistoPair[kTrackVariablesPair] = { 8 		, 8		, ptMax		, ptMax	 	,	1.5*TMath::Pi()	,	fDeltaEtaMax	,	zvtxMax		, centMax	, 9,9};

//   Int_t binsHistoPair[kTrackVariablesPair]   = { iBinPair[0]	, iBinPair[1]	, iBinPair[2]	, iBinPair[3]	,	iBinPair[4]	,	iBinPair[5]	,	iBinPair[6]	, 2,2};
//   Double_t minHistoPair[kTrackVariablesPair] = { 1	 	, 1		, ptMin		, ptMin	 	,	-0.5*TMath::Pi()	,	-fDeltaEtaMax	,	zvtxMin		, 7,7};
//   Double_t maxHistoPair[kTrackVariablesPair] = { 7 		, 7		, ptMax		, ptMax	 	,	1.5*TMath::Pi()	,	fDeltaEtaMax	,	zvtxMax		, 9,9};

//   Int_t binsHistoPair[kTrackVariablesPair]   = { iBinPair[0]	, iBinPair[1]	, iBinPair[2]	,	iBinPair[3]	,	iBinPair[4]	,	iBinPair[5]	, 2};
//   Double_t minHistoPair[kTrackVariablesPair] = { 1	 	, ptMin		, ptMin	 	,	-0.5*TMath::Pi()	,	-fDeltaEtaMax	,	zvtxMin		, 7};
//   Double_t maxHistoPair[kTrackVariablesPair] = { 7 		, ptMax		, ptMax	 	,	1.5*TMath::Pi()	,	fDeltaEtaMax	,	zvtxMax		, 9};
  
  

  Int_t binsHistoPair[kTrackVariablesPair]   = { iBinPair[0]	, iBinPair[1]	,	iBinPair[2]	,	iBinPair[3]	,	iBinPair[4]	};
  Double_t minHistoPair[kTrackVariablesPair] = { static_cast<Double_t>(ptMin)		, static_cast<Double_t>(ptMin)	 	,	-0.5*TMath::Pi()	,	static_cast<Double_t>(-fDeltaEtaMax)	,	static_cast<Double_t>(zvtxMin)		};
  Double_t maxHistoPair[kTrackVariablesPair] = { static_cast<Double_t>(ptMax)		, static_cast<Double_t>(ptMax)	 	,	1.5*TMath::Pi()	,	static_cast<Double_t>(fDeltaEtaMax)	,	static_cast<Double_t>(zvtxMax)		};


/*
  Double_t ptMin = 0.0, ptMax = 6.0;
//   Double_t binsPt[]  = {0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0};
  Double_t binsPt[]  = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const Int_t nPtBins = length(binsPt)-1;

  Double_t zvtxMin = -10., zvtxMax = 10.;
  Double_t binsZvtx[] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
  const Int_t nZvtxBins = length(binsZvtx)-1;

  Double_t centMin = 0., centMax = 10.;
//   Double_t binsCent[] = {0.,1.,2.,4.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.1}; // 90
  Double_t binsCent[] = {0.,1.,2.,5.,10.}; // 90
  const Int_t nCentBins = length(binsCent)-1;


  //_____ single histo, 4 dim matrix
  Int_t binsHisto4[4+1]   = { 6 ,	nPtBins	, nZvtxBins	, nCentBins	,2};
  Double_t minHisto4[4+1] = { 1 ,	ptMin	, zvtxMin	, centMin     	,7};
  Double_t maxHisto4[4+1] = { 7 ,	ptMax  	, zvtxMax     	, centMax	,9};

  //_____ pair histo, 8 dim matrix
  Int_t binsHistoPair[8+2]   = { 6 , 6	 , nPtBins	, nPtBins ,   		     72 	,  Int_t(fDeltaEtaMax*10) 	,  nZvtxBins	, nCentBins	,2,2};
  Double_t minHistoPair[8+2] = { 1 , 1	 , ptMin		, ptMin	  ,  	-0.5*TMath::Pi()	, -fDeltaEtaMax   	 	,  zvtxMin	, centMin	,7,7};
  Double_t maxHistoPair[8+2] = { 7 , 7	 , ptMax		, ptMax	  ,	 1.5*TMath::Pi()	,  fDeltaEtaMax   	 	,  zvtxMax	, centMax	,9,9};
*/
/*
  fHistCorrSingle = new THnSparseD("fHistCorrSingle","pdg_trig:pt_trig:Zvtx:centr:chhdr_trig",4+1,binsHisto4,minHisto4,maxHisto4);
  fHistCorrSingle->GetAxis(0)->SetTitle("partType_{trig}");
  fHistCorrSingle->SetBinEdges(1,dBinsSingle[1]);
//   fHistCorrSingle->SetBinEdges(1,binsPt);
  fHistCorrSingle->GetAxis(1)->SetTitle("p_{T, trig} (GeV/c)");
  fHistCorrSingle->SetBinEdges(2,dBinsSingle[2]);
//   fHistCorrSingle->SetBinEdges(2,binsZvtx);
  fHistCorrSingle->GetAxis(2)->SetTitle("v_{z} (cm)");
  fHistCorrSingle->SetBinEdges(3,dBinsSingle[3]);
//   fHistCorrSingle->SetBinEdges(3,binsCent);
  fHistCorrSingle->GetAxis(3)->SetTitle("Centrality [#%]");
  fHistCorrSingle->GetAxis(4)->SetTitle("hadron_trig");

  fHistCorrSingle->Sumw2();
  fList -> Add(fHistCorrSingle);
*/

  for(Int_t i=0;i<2;i++)
  {
    fHistCorrPair[i] = new THnSparseD(Form("fHistCorrPair_%d",i),"pt_trig:pt_assoc:DeltaPhi:DeltaEta:Zvtx",5,binsHistoPair,minHistoPair,maxHistoPair);
    fHistCorrPair[i]->SetBinEdges(0,dBinsPair[0]);
    fHistCorrPair[i]->SetBinEdges(1,dBinsPair[1]);
    fHistCorrPair[i]->GetAxis(0)->SetTitle("p_{T, trig} (GeV/c)");
    fHistCorrPair[i]->GetAxis(1)->SetTitle("p_{T, assoc} (GeV/c)");
    fHistCorrPair[i]->SetBinEdges(2,dBinsPair[2]);
    fHistCorrPair[i]->SetBinEdges(3,dBinsPair[3]);
    fHistCorrPair[i]->GetAxis(2)->SetTitle("#Delta #varphi");
    fHistCorrPair[i]->GetAxis(3)->SetTitle("#Delta #eta");
    fHistCorrPair[i]->SetBinEdges(4,dBinsPair[4]);
    fHistCorrPair[i]->GetAxis(4)->SetTitle("v_{z} (cm)");
    fHistCorrPair[i]->Sumw2();
    fList -> Add(fHistCorrPair[i]);
  }

  //     TAxis* axis = fHistCorrPairSame->GetGrid(0)->GetGrid()->GetAxis(2);
  fTriggerWeighting = fHistCorrPair[0]->GetAxis(0);
  fHistTriggerWeighting = new TH1D("fHistTriggerWeighting","", fTriggerWeighting->GetNbins(),fTriggerWeighting->GetXbins()->GetArray());
    


/*
  for(Int_t i=0;i<2;i++)
  {
    fHistCorrPair[i] = new THnSparseD(Form("fHistCorrPair_%d",i),"pdg_trig:pt_trig:pt_assoc:DeltaPhi:DeltaEta:Zvtx:chhdr_assoc",7,binsHistoPair,minHistoPair,maxHistoPair);
    fHistCorrPair[i]->GetAxis(0)->SetTitle("partType_{trig}");
//     fHistCorrPair[i]->GetAxis(1)->SetTitle("partType_{assoc}");
    fHistCorrPair[i]->SetBinEdges(0,dBinsPair[0]);
//     fHistCorrPair[i]->SetBinEdges(1,dBinsPair[1]);
//     fHistCorrPair[i]->SetBinEdges(2,binsPt);
//     fHistCorrPair[i]->SetBinEdges(3,binsPt);
    fHistCorrPair[i]->GetAxis(1)->SetTitle("p_{T, trig} (GeV/c)");
    fHistCorrPair[i]->GetAxis(2)->SetTitle("p_{T, assoc} (GeV/c)");
    fHistCorrPair[i]->SetBinEdges(1,dBinsPair[1]);
    fHistCorrPair[i]->SetBinEdges(2,dBinsPair[2]);
    fHistCorrPair[i]->GetAxis(3)->SetTitle("#Delta #varphi");
    fHistCorrPair[i]->GetAxis(4)->SetTitle("#Delta #eta");
    fHistCorrPair[i]->SetBinEdges(3,dBinsPair[3]);
    fHistCorrPair[i]->SetBinEdges(4,dBinsPair[4]);
//     fHistCorrPair[i]->SetBinEdges(6,binsZvtx);
//     fHistCorrPair[i]->SetBinEdges(7,binsCent);
    fHistCorrPair[i]->GetAxis(5)->SetTitle("v_{z} (cm)");
    fHistCorrPair[i]->SetBinEdges(5,dBinsPair[5]);
//     fHistCorrPair[i]->GetAxis(7)->SetTitle("hadron_trig");
    fHistCorrPair[i]->GetAxis(6)->SetTitle("hadron_assoc");

    fHistCorrPair[i]->Sumw2();
    fList -> Add(fHistCorrPair[i]);
  }
*/
/*
  if( fWriteCorrTree == kTRUE )
  {
    fVariablesTreeCorr = new TTree("CorrTree","2PartCorrPID tree");
    Int_t nVar = 5;
    fCorrVariables = new Float_t [nVar];
    TString * fCorrVariableNames = new TString[nVar];
    fCorrVariableNames[0] = "fMyPIDtrig";
    fCorrVariableNames[1] = "fMyPTtrig";
    fCorrVariableNames[2] = "fMyZvtx";
    fCorrVariableNames[3] = "fMyCent";
    fCorrVariableNames[4] = "fMyChardHdr";
    
    for(Int_t ivar=0; ivar<nVar; ivar++)
      fVariablesTreeCorr->Branch(fCorrVariableNames[ivar].Data(),&fCorrVariables[ivar],Form("%s/f",fCorrVariableNames[ivar].Data()));
    PostData(3,fVariablesTreeCorr);
  }*/
  
  for (Int_t i=0; i<fNMaxBinsPt; i++)
  {
    fHistTwoTrackDistancePt[i][0] = new TH2F(Form("fHistTwoTrackDistancePt[%d][0]",i),Form("fHistTwoTrackDistancePt[%d][0] -- PtBin %.2f-%.2f;#Delta#eta;#Delta#varphi^{*}_{min}",i,fPtAxis->GetBinLowEdge(i+1),fPtAxis->GetBinUpEdge(i+1)), 100, -0.05, 0.05, 400, -0.2, 0.2);
    fList->Add(fHistTwoTrackDistancePt[i][0]);
    fHistTwoTrackDistancePt[i][1] = (TH2F*)fHistTwoTrackDistancePt[i][0]->Clone(Form("fHistTwoTrackDistancePt[%d][1]",i));
    fList->Add(fHistTwoTrackDistancePt[i][1]);
  }
  
  // QA
  fHistHBTbefore        = new TH2F("fHistHBTbefore","before HBT cut",20,0,2,20,0,2.*TMath::Pi()); fList->Add(fHistHBTbefore);
  fHistHBTafter         = new TH2F("fHistHBTafter","after HBT cut",20,0,2,20,0,2.*TMath::Pi()); fList->Add(fHistHBTafter);
/*  
  fHistConversionbefore = new TH2F("fHistConversionbefore","before Conversion cut",200,0,2,200,0,2.*TMath::Pi());
  fHistConversionafter  = new TH2F("fHistConversionafter","after Conversion cut",200,0,2,200,0,2.*TMath::Pi());
  fHistPsiMinusPhi      = new TH2F("fHistPsiMinusPhi","",4,-0.5,3.5,100,0,2.*TMath::Pi());
  fHistResonancesBefore = new TH3D("fHistResonancesBefore","before resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesRho    = new TH3D("fHistResonancesRho","after #rho resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesK0     = new TH3D("fHistResonancesK0","after #rho, K0 resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesLambda = new TH3D("fHistResonancesLambda","after #rho, K0, Lambda resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistQbefore          = new TH3D("fHistQbefore","before momentum difference cut;#Delta#eta;#Delta#phi;|#Delta p_{T}| (GeV/c)",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistQafter           = new TH3D("fHistQafter","after momentum difference cut;#Delta#eta;#Delta#phi;|#Delta p_{T}| (GeV/c)",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
*/

  //_____ PoolManager setup for mixing
  SetupForMixing();

  //_____ AliCFContainer settings
  Int_t nTrackBin[kNvars];
  Double_t* trackBins[kNvars];
  TString trackAxisTitle[kNvars];

  TString defaultBinningStr;
  defaultBinningStr = 	"eta: -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\n"
// 			"phi: 0., 0.62831, 1.25662, 1.88493, 2.51324, 3.14156, 3.76987, 4.39818, 5.02649, 5.65480, 6.28312\n"
			"p_T: 0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0\n"
			"PID: 0., 1., 2., 3., 4., 5., 6., 7., 8., 9.\n"
			"centrality: 0.0, 5.0\n"
			"vertex: -10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.\n";
 // =========================================================
  // Customization (adopted from AliUEHistograms)
  // =========================================================
  TObjArray* lines = defaultBinningStr.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    TString tag = line(0, line.Index(":")+1);
    if (!fCustomBinning.BeginsWith(tag) && !fCustomBinning.Contains(TString("\n") + tag))
      fBinningString += line + "\n";
    else
      AliInfo(Form("Using custom binning for %s", tag.Data()));
  }
  delete lines;
  fBinningString += fCustomBinning;
  
  AliInfo(Form("Used THn Binning:\n%s",fBinningString.Data()));

  
  // eta
  trackBins[0] = GetBinning(fBinningString, "eta", nTrackBin[0]);
  trackAxisTitle[0] = "#eta";
//   Int_t nEtaBins = 20;
//   Double_t etaMin = -2., etaMax = -2.;
  // phi
//   trackBins[1] = GetBinning(fBinningString, "phi", nTrackBin[1]);  
//   trackAxisTitle[1] = "#phi (rad)";
//   Int_t nPhiBins = 36;
//   Double_t phiMin = 0.; Double_t phiMax = 2.*TMath::Pi();
  // pT
  trackBins[1] = GetBinning(fBinningString, "p_T", nTrackBin[1]);
  trackAxisTitle[1] = "p_{T} (GeV/c)";
//   Int_t nPtBins = 20;
//   Double_t ptMin = 0., ptMax = 6.  
  // particle species
  trackBins[2] = GetBinning(fBinningString, "PID", nTrackBin[2]);
  trackAxisTitle[2] = "particle species";
//   Int_t nPidBins = 9;
//   Double_t pidMin = -0.5, pidMax = 8.5;
  // centrality
  trackBins[3] = GetBinning(fBinningString, "centrality", nTrackBin[3]);
  trackAxisTitle[3] = "cent (%)";
//   Double_t centMin = 0., centMax = 100.1;
//   Double_t binsCent[] = {0.,1.,2.,4.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.1};
//   Int_t nCentBins = length(binsCent)-1;  
  // vtx-z axis
  trackBins[4] = GetBinning(fBinningString, "vertex", nTrackBin[4]);
  trackAxisTitle[4] = "z-vtx (cm)";  
//   Int_t nVzBins = 40;
//   Double_t vzMin = -10., vzMax = 10.;  
  
//   Int_t nbins[kNvars] = {nEtaBins, nPhiBins, nPtBins, nPidBins, nCentBins, nVzBins};
//   Double_t xmin[kNvars] = {etaMin, phiMin, ptMin, pidMin, centMin, vzMin};
//   Double_t xmax[kNvars] = {etaMax, phiMax, ptMax, pidMax, centMax, vzMax};  
  
  fMyCFCont = new AliCFContainer("PidPidCorrCFCont","CF for PidPid Corr",kNsteps,kNvars,nTrackBin);
//   fMyCFCont = new AliCFContainer("PidPidCorrCFCont","CF for PidPid Corr",kNsteps,kNvars,nbins);
  
  for ( Int_t idim=0; idim < kNvars; idim++)
  {
    fMyCFCont->SetVarTitle(idim, Form("%s",trackAxisTitle[idim].Data()));
    fMyCFCont->SetBinLimits(idim, trackBins[idim]);
//     fMyCFCont->SetBinLimits(idim, xmin[idim], xmax[idim]);
  }
  
  TString stepTitle[kNsteps] = {"generated","reco","recomatch"};
  for (Int_t istep=0; istep<kNsteps; istep++)
    fMyCFCont->SetStepTitle(istep, stepTitle[istep].Data());
  
  PostData(1, fList);
  PostData(2, fMyCFCont);  
}

//________________________________________________________________________
Double_t* AliAnalysisTaskPidPidCorrelations::GetBinning(const Char_t* configuration, const Char_t* tag, Int_t& nBins)
{
  // takes the binning from <configuration> identified by <tag>
  // configuration syntax example:
  // eta: 2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
  // phi: .....
  //
  // returns bin edges which have to be deleted by the caller
  
  TString config(configuration);
  TObjArray* lines = config.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    if (line.BeginsWith(TString(tag) + ":"))
    {
      line.Remove(0, strlen(tag) + 1);
      line.ReplaceAll(" ", "");
      TObjArray* binning = line.Tokenize(",");
      Double_t* bins = new Double_t[binning->GetEntriesFast()];
      for (Int_t j=0; j<binning->GetEntriesFast(); j++)
	bins[j] = TString(binning->At(j)->GetName()).Atof();
      
      nBins = binning->GetEntriesFast() - 1;

      delete binning;
      delete lines;
      return bins;
    }
  }
  
  delete lines;
  AliFatal(Form("Tag %s not found in %s", tag, configuration));
  return 0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPidPidCorrelations::CheckTrack(AliAODTrack* track)
{
  // check if the track status flags are set  
  UInt_t status = track -> GetStatus();
  if ((status & fTrackStatus) == fTrackStatus)
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::CalculateNSigmas(AliAODTrack* track, Int_t centBin, Bool_t* pidFlag, Bool_t fillQA)
{
  if(!fPIDResponse)
    fPIDResponse = (dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager() -> GetInputEventHandler())) -> GetPIDResponse(); 
  if(!fPIDResponse)
    AliFatal("AliAnalysisTaskPidPidCorrelations::CalculateNSigmas(): Cannot get pid response");

  // ____________ set cuts on pT, nsigma for ITS, TPC, TOF, HMPID
  Double_t fgTPCPIDmomcut[AliPID::kSPECIES]	= { 0., 0., 0.5, 0.5, 0.7 };	// _____ TPC pT
  Double_t fgTPCPIDsigmacut[AliPID::kSPECIES]	= { 0., 0., 3. , 3. , 3.  };	// _____ TPC nsigma	
  Double_t fgTOFPIDmomcut[AliPID::kSPECIES]	= { 0., 0., 1.5, 1.5, 2.0 };	// _____ TOF pT
  Double_t fgTOFPIDsigmacut[AliPID::kSPECIES]	= { 0., 0., 2. , 2. , 2.  };	// _____ TOF nsgima
// 
  Double_t p	= track -> P();
  Double_t pt	= track -> Pt();
  Double_t ptpc	= track -> GetTPCmomentum();

  //___ nsigmas
//   Double_t nsigmaITS[AliPID::kSPECIES]	= {-999,-999,-999,-999,-999};	// 	ITS
//   Double_t nsigmaTPC[AliPID::kSPECIES]	= {-999,-999,-999,-999,-999};	// 	TPC
//   Double_t nsigmaTOF[AliPID::kSPECIES]	= {-999,-999,-999,-999,-999};	//	TOF
//   Double_t nsigmaHMPID[AliPID::kSPECIES]	= {-999,-999,-999,-999,-999};	//	HMPID
  
  //___ auxiliary nsigma for QA
  Double_t auxNsigmaTPC[AliPID::kSPECIES]	= {-999,-999,-999,-999,-999};	//	TPC
  Double_t auxNsigmaTOF[AliPID::kSPECIES]	= {-999,-999,-999,-999,-999};	//	TOF

  Double_t dEdx = MakeTPCPID(track, auxNsigmaTPC);
  Double_t beta = MakeTOFPID(track, auxNsigmaTOF);

  Bool_t isOnTPCPID = dEdx > 0.;
  Bool_t isOnTOFPID = beta > 0.;

  Bool_t ITSnSigmaIsOk[AliPID::kSPECIES] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
  Bool_t TPCnSigmaIsOk[AliPID::kSPECIES] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
  Bool_t TOFnSigmaIsOk[AliPID::kSPECIES] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
  Bool_t HMPnSigmaIsOk[AliPID::kSPECIES] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};

  for(Int_t ispecies=0; ispecies < AliPID::kSPECIES; ispecies++)
  {
    //_____ ITS
    nsigmaITS[ispecies] = fPIDResponse -> NumberOfSigmasITS((AliVTrack*)track,(AliPID::EParticleType)ispecies);
    if (nsigmaITS[ispecies] != -999.) ITSnSigmaIsOk[ispecies] = kTRUE;
    //_____ TPC
    nsigmaTPC[ispecies] = fPIDResponse -> NumberOfSigmasTPC((AliVTrack*)track,(AliPID::EParticleType)ispecies);
    if (nsigmaTPC[ispecies] != -999.) TPCnSigmaIsOk[ispecies] = kTRUE;
    //_____ TOF
    nsigmaTOF[ispecies] = fPIDResponse -> NumberOfSigmasTOF((AliVTrack*)track,(AliPID::EParticleType)ispecies);
    if (nsigmaTOF[ispecies] != -999.) TOFnSigmaIsOk[ispecies] = kTRUE;
    //_____ HMPID
    nsigmaHMPID[ispecies] = fPIDResponse -> NumberOfSigmasHMPID((AliVTrack*)track,(AliPID::EParticleType)ispecies);
    if (nsigmaHMPID[ispecies] != -999.) HMPnSigmaIsOk[ispecies] = kTRUE;
  }

  for(Int_t ispecies=0; ispecies < AliPID::kSPECIES; ispecies++)
  {
    if (ITSnSigmaIsOk[ispecies]) fnsigmas[ispecies][0] = nsigmaITS[ispecies];
    if (TPCnSigmaIsOk[ispecies]) fnsigmas[ispecies][1] = nsigmaTPC[ispecies];
    if (TOFnSigmaIsOk[ispecies]) fnsigmas[ispecies][2] = nsigmaTOF[ispecies];
    if (HMPnSigmaIsOk[ispecies]) fnsigmas[ispecies][3] = nsigmaHMPID[ispecies];
  }

  //______ fill QA
  if (fillQA)
  {
    TH2F* h[fNMaxBinsCentrality];

    for (Int_t iSpecies=0; iSpecies < AliPID::kSPECIES; iSpecies++)
    {
      pidFlag[iSpecies] = kFALSE;
      
      //_______________ PID signal plots for TPC, TOF
      if (isOnTPCPID)
      {
     	fHistNSigmaTPC[centBin][iSpecies] -> Fill(pt, nsigmaTPC[iSpecies]);
	fHistTPCdEdx[centBin] -> Fill(ptpc, dEdx);
      }
      if (isOnTOFPID)
      {
      	fHistNSigmaTOF[centBin][iSpecies] -> Fill(pt, nsigmaTOF[iSpecies]);
    	fHistTOFbeta[centBin] -> Fill(p, beta);
      }

      //_______________ PID nsigma plots for TPCTOF, TPC, TOF
      if (isOnTPCPID && isOnTOFPID)
      {
	// get TPC & TOF nsigma -- combined
	if (pt > fgTOFPIDmomcut[iSpecies]) // cut on TOF pT
// 	nsigmaTPCTOF[iSpecies] = TMath::Sqrt(nsigmaTPC[iSpecies]*nsigmaTPC[iSpecies]+nsigmaTOF[iSpecies]*nsigmaTOF[iSpecies]);

      	if (pt < fgTOFPIDmomcut[iSpecies] && TMath::Abs(nsigmaTOF[iSpecies]) < fgTOFPIDsigmacut[iSpecies] && TMath::Abs(nsigmaTPC[iSpecies]) < fgTPCPIDmomcut[iSpecies])
      	{
	  fHistTOFbeta_selected[centBin] -> Fill(p, beta);
	  pidFlag[iSpecies] = kTRUE;
	}
      }
      //_____ TPC-only PID cuts
      else if (isOnTPCPID && !isOnTOFPID)
      {
	if (pt < fgTPCPIDmomcut[iSpecies] && TMath::Abs(nsigmaTPC[iSpecies]) < fgTPCPIDsigmacut[iSpecies])
	{
	  fHistTPCdEdx_selected[centBin] -> Fill(ptpc, dEdx);
	  pidFlag[iSpecies] = kTRUE;
	}
      }
      //_____ TOF-only PID cuts
      else if (!isOnTPCPID && isOnTOFPID)
      {
      	if (pt < fgTOFPIDmomcut[iSpecies] && TMath::Abs(nsigmaTOF[iSpecies]) < fgTOFPIDsigmacut[iSpecies])
      	{
	  fHistTOFbeta_selected[centBin] -> Fill(p, beta);
	  pidFlag[iSpecies] = kTRUE;
	}
      }
  
      //__________ fill nsigma vs pT for all the detectors
      for (int iPid=0; iPid < kMyNSigmaPIDType; iPid++)
      {
	if (fUseMC)
	  h[centBin] = GetHisto2D(Form("NSigmaRec_Cent%02d_%d_%d",centBin,iSpecies,iPid));
	else
	  h[centBin] = GetHisto2D(Form("NSigma_Cent%02d_%d_%d",centBin,iSpecies,iPid));
	
	h[centBin] -> Fill(pt,fnsigmas[iSpecies][iPid]);
      }
    }//____ iSpecies   
  }//____ fill QA
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::FindNSigma(AliAODTrack* track) 
{
  //_____ ITS
  Double_t fgITSPIDmomcutMin[AliPID::kSPECIES]	= { 0., 0., 0.2 , 0.6 , 0.4 };	// pT_min
  Double_t fgITSPIDsigmacutMin[AliPID::kSPECIES]	= { 0., 0.,-1.5 ,-2.0 ,-6.0 };	// nsigma_min
  Double_t fgITSPIDmomcutMax[AliPID::kSPECIES]	= { 0., 0., 1.2 , 1.4 , 2.0 };	// pT_max
  Double_t fgITSPIDsigmacutMax[AliPID::kSPECIES]	= { 0., 0., 1.5 , 1.0 , 1.0 };	// nsigma_max
  //_____ TPC
  Double_t fgTPCPIDmomcutMin[AliPID::kSPECIES]	= { 0., 0., 0.2 , 0.6 , 0.6 };	// pT_min
  Double_t fgTPCPIDsigmacutMin[AliPID::kSPECIES]	= { 0., 0.,-2.0 ,-2.0 ,-4.0 };	// nsigma_min
  Double_t fgTPCPIDmomcutMax[AliPID::kSPECIES]	= { 0., 0., 1.0 , 2.0 , 4.0 };	// pT_max
  Double_t fgTPCPIDsigmacutMax[AliPID::kSPECIES]	= { 0., 0., 2.0 , 3.0 , 3.0 };	// nsigma_max
  //_____ TOF
  Double_t fgTOFPIDmomcutMin[AliPID::kSPECIES]	= { 0., 0., 0.4 , 0.6 , 0.6 };	// pT_min
  Double_t fgTOFPIDsigmacutMin[AliPID::kSPECIES]	= { 0., 0.,-1.0 ,-2.0 ,-2.0 };	// nsgima_min
  Double_t fgTOFPIDmomcutMax[AliPID::kSPECIES]	= { 0., 0., 1.2 , 2.6 , 4.0 };	// pT_max
  Double_t fgTOFPIDsigmacutMax[AliPID::kSPECIES]	= { 0., 0., 1.0 , 2.0 , 3.0 };	// nsgima_max
  //_____ HMPID
  Double_t fgHMPIDPIDmomcutMin[AliPID::kSPECIES]		= { 0., 0., 0.5 , 0.8 , 1.4 };	// pT_min
  Double_t fgHMPIDPIDsigmacutMin[AliPID::kSPECIES]	= { 0., 0.,-2.0 ,-2.0 ,-2.0 };	// nsigma_min
  Double_t fgHMPIDPIDmomcutMax[AliPID::kSPECIES]		= { 0., 0., 1.5 , 1.5 , 4.0 };	// pT_max
  Double_t fgHMPIDPIDsigmacutMax[AliPID::kSPECIES]	= { 0., 0., 3.0 , 6.0 , 4.0 };	// nsigma_max

  Bool_t isITSok[AliPID::kSPECIES]	=  {kFALSE};
  Bool_t isTPCok[AliPID::kSPECIES] 	=  {kFALSE};
  Bool_t isTOFok[AliPID::kSPECIES] 	=  {kFALSE};
  Bool_t isHMPIDok[AliPID::kSPECIES]   	=  {kFALSE};

  //Int_t trackCharge = (track->Charge()>0) ? 1 : -1;
//Int_t trackCharge = 0;
//trackCharge = sign(track->Charge());

//if (trackCharge>0) 		return fPartHadronPlus;
//else if (trackCharge<0) 	return fPartHadronMinus;

  //Printf(" --->>>  %f <--> %f <--> %f <--> %f <--> %f <--> %f",track->Pt(),fgITSPIDmomcutMin[2], fgITSPIDmomcutMax[2], fgITSPIDsigmacutMin[2], fgITSPIDsigmacutMax[2],fnsigmas[2][0]);
  
  for(Int_t iSpecies=0; iSpecies < AliPID::kSPECIES; iSpecies++)
  {
    if ( (track->Pt() > fgITSPIDmomcutMin[iSpecies] && track->Pt() < fgITSPIDmomcutMax[iSpecies]) && (fnsigmas[iSpecies][0] > fgITSPIDsigmacutMin[iSpecies] && fnsigmas[iSpecies][0] < fgITSPIDsigmacutMax[iSpecies]) )
      isITSok[iSpecies] = kTRUE;
    if ( (track->Pt() > fgTPCPIDmomcutMin[iSpecies] && track->Pt() < fgTPCPIDmomcutMax[iSpecies]) && (fnsigmas[iSpecies][1] > fgTPCPIDsigmacutMin[iSpecies] && fnsigmas[iSpecies][1] < fgTPCPIDsigmacutMax[iSpecies]) )
      isTPCok[iSpecies] = kTRUE;
    if ( (track->Pt() > fgTOFPIDmomcutMin[iSpecies] && track->Pt() < fgTOFPIDmomcutMax[iSpecies]) && (fnsigmas[iSpecies][2] > fgTOFPIDsigmacutMin[iSpecies] && fnsigmas[iSpecies][2] < fgTOFPIDsigmacutMax[iSpecies]) )
      isTOFok[iSpecies] = kTRUE;
    if ( (track->Pt() > fgHMPIDPIDmomcutMin[iSpecies] && track->Pt() < fgHMPIDPIDmomcutMax[iSpecies]) && (fnsigmas[iSpecies][3] > fgHMPIDPIDsigmacutMin[iSpecies] && fnsigmas[iSpecies][3] < fgHMPIDPIDsigmacutMax[iSpecies]) )
      isHMPIDok[iSpecies] = kTRUE;
    
    switch(iSpecies)
    {
      case 2:
	if ( (isITSok[iSpecies] == kTRUE || isTPCok[iSpecies] == kTRUE || isTOFok[iSpecies] == kTRUE || isHMPIDok[iSpecies] == kTRUE) && track->Charge()>0 )
	  return fPartPionPlus;
	else if ( (isITSok[iSpecies] == kTRUE || isTPCok[iSpecies] == kTRUE || isTOFok[iSpecies] == kTRUE || isHMPIDok[iSpecies] == kTRUE) && track->Charge()<0 )
	  return fPartPionMinus;
      break;
      case 3:
	if ( (isITSok[iSpecies] == kTRUE || isTPCok[iSpecies] == kTRUE || isTOFok[iSpecies] == kTRUE || isHMPIDok[iSpecies] == kTRUE) && track->Charge()>0 )
	  return fPartKaonPlus;
	else if ( (isITSok[iSpecies] == kTRUE || isTPCok[iSpecies] == kTRUE || isTOFok[iSpecies] == kTRUE || isHMPIDok[iSpecies] == kTRUE) && track->Charge()<0 )
	  return fPartKaonMinus;
      break;
      case 4:
	if ( (isITSok[iSpecies] == kTRUE || isTPCok[iSpecies] == kTRUE || isTOFok[iSpecies] == kTRUE || isHMPIDok[iSpecies] == kTRUE) && track->Charge()>0 )
	  return fPartProtonPlus;
	else if ( (isITSok[iSpecies] == kTRUE || isTPCok[iSpecies] == kTRUE || isTOFok[iSpecies] == kTRUE || isHMPIDok[iSpecies] == kTRUE) && track->Charge()<0 )
	  return fPartProtonMinus;
      break;
    }
  }//iSpecies
  
  return fPartUndefined;
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::GetParticleID(AliVTrack* trk, Int_t centbin, Bool_t fillQA)
{
  if (!trk) { AliInfo(" ==== zero track pointer."); return fPartUndefined; }

  Bool_t* pidFlag; pidFlag = new Bool_t[AliPID::kSPECIES];
  CalculateNSigmas((AliAODTrack*)trk, centbin, pidFlag, fillQA);
  delete[] pidFlag;
  Int_t mypid = -1;
  mypid = FindNSigma((AliAODTrack*)trk);
//     Printf(" >>>>>>>>>>>>>>>>> mypid = %d",mypid);
  return mypid;
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::GetParticleIDMC(AliVTrack* trk, Int_t centbin, Bool_t fillQA)
{
  if (!trk) { AliInfo(" ==== zero track pointer."); return fPartUndefined; }

  if (fUseMC)//______ MC
  {
    if(!fMyMcArray) AliFatal("Error: AliAnalysisTaskPidPidCorrelations::GetParticleID called on data\n");
    
    Int_t pdgCode = 999;
//     AliAODMCParticle* partMC = dynamic_cast<AliAODMCParticle*>(fMyMcArray->At(TMath::Abs(trk->GetLabel())));
    AliAODMCParticle* partMC = (AliAODMCParticle*)fMyMcArray->At(TMath::Abs(trk->GetLabel()));
    if (!partMC){ AliError("Cannot get MC particle"); return fPartUndefined; }
    pdgCode = partMC->GetPdgCode();
//     delete partMC; partMC=0x0;

    //Int_t qPart = 0;
    //qPart  = (Int_t)partMC->GetPDG()->Charge();
    //if (qPart<0) return fPartHadronMinus;
    //else if (qPart<0) return fPartHadronPlus;
		
    switch(pdgCode)
    {
      case 2212:
	if (fillQA) {for(Int_t ipid=0;ipid<kMyNSigmaPIDType;ipid++){TH2F* h = GetHisto2D(Form("NSigmaMC_Cent%02d_%d_%d",centbin,fPartProtonQA,ipid));h -> Fill(trk->Pt(),fnsigmas[fPartProtonQA][ipid]);}}
	  return fPartProtonPlus;
      break;
      case -2212:
	if (fillQA) {for(Int_t ipid=0;ipid<kMyNSigmaPIDType;ipid++){TH2F* h = GetHisto2D(Form("NSigmaMC_Cent%02d_%d_%d",centbin,fPartProtonQA,ipid));h -> Fill(trk->Pt(),fnsigmas[fPartProtonQA][ipid]);}}
	  return fPartProtonPlus;
      break;
      case 321:
	if (fillQA) {for(Int_t ipid=0;ipid<kMyNSigmaPIDType;ipid++){TH2F* h = GetHisto2D(Form("NSigmaMC_Cent%02d_%d_%d",centbin,fPartKaonQA,ipid));h -> Fill(trk->Pt(),fnsigmas[fPartKaonQA][ipid]);}}
	  return fPartKaonPlus;
      break;
      case -321:
	if (fillQA) {for(Int_t ipid=0;ipid<kMyNSigmaPIDType;ipid++){TH2F* h = GetHisto2D(Form("NSigmaMC_Cent%02d_%d_%d",centbin,fPartKaonQA,ipid));h -> Fill(trk->Pt(),fnsigmas[fPartKaonQA][ipid]);}}
	  return fPartKaonMinus;
      break;
      case 211:
	if (fillQA) {for(Int_t ipid=0;ipid<kMyNSigmaPIDType;ipid++){TH2F* h = GetHisto2D(Form("NSigmaMC_Cent%02d_%d_%d",centbin,fPartPionQA,ipid));h -> Fill(trk->Pt(),fnsigmas[fPartPionQA][ipid]);}}
	  return fPartPionPlus;
      break;
      case -211:
	if (fillQA) {for(Int_t ipid=0;ipid<kMyNSigmaPIDType;ipid++){TH2F* h = GetHisto2D(Form("NSigmaMC_Cent%02d_%d_%d",centbin,fPartPionQA,ipid));h -> Fill(trk->Pt(),fnsigmas[fPartPionQA][ipid]);}}
	  return fPartPionMinus;
      break;
      default:
	return fPartUndefined;
    }
  } // MC
  return fPartUndefined;
}

//___________________________________________________________
Double_t AliAnalysisTaskPidPidCorrelations::MakeTPCPID(AliAODTrack* track, Double_t* nSigma)
{
  	//___ fills nSigma array with TPC nsigmas for e, mu, pi, K, p
  
  	/* check TPC PID */
  	if (!HasTPCPID(track)) return -1.;

  	/* get TPC info */
  	Double_t ptpc	= track -> GetTPCmomentum();
  	Double_t dEdx	= track -> GetTPCsignal();
  	Int_t dEdxN	= track -> GetTPCsignalN();
    
  	/* loop over particles */
  	for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
  	{
    		Double_t bethe = fPIDResponse -> GetTPCResponse().GetExpectedSignal(ptpc, (AliPID::EParticleType)iSpecies);
    		Double_t diff = dEdx - bethe;
    		Double_t sigma = fPIDResponse -> GetTPCResponse().GetExpectedSigma(ptpc, dEdxN, (AliPID::EParticleType)iSpecies);
    		nSigma[iSpecies] = diff/sigma;
  	}
  	return dEdx;
}

//___________________________________________________________
Double_t AliAnalysisTaskPidPidCorrelations::MakeTOFPID(AliAODTrack* track, Double_t* nSigma)
{
  	//___ fills nSigma array with TOF nsigmas for e, mu, pi, K, p

  	/* check TOF PID */
  	if (!HasTOFPID(track)) return -1.;

  	/* get TOF info */
  	Double_t p	= track -> P();
  	Double_t time	= track -> GetTOFsignal() - fPIDResponse -> GetTOFResponse().GetStartTime(p);
  	Double_t timei[5];
  	track -> GetIntegratedTimes(timei,AliPID::kSPECIES);
  
 	/* loop over particles */
  	for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
  	{
    		Double_t timez = time - timei[iSpecies];
    		Double_t sigma = fPIDResponse -> GetTOFResponse().GetExpectedSigma(p, timei[iSpecies], AliPID::ParticleMass(iSpecies));
    		nSigma[iSpecies] = timez/sigma;
  	}
  	//return timei[0]/time;
	return GetBeta(track);
}

//___________________________________________________________
Bool_t AliAnalysisTaskPidPidCorrelations::HasTPCPID(AliAODTrack* track) const
{
  // check PID signal 
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse -> CheckPIDStatus(AliPIDResponse::kTPC,track);
  if (statusTPC != AliPIDResponse::kDetPidOk)
    return kFALSE;
  return kTRUE;
}

//___________________________________________________________
Bool_t AliAnalysisTaskPidPidCorrelations::HasTOFPID(AliAODTrack* track) const
{
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse -> CheckPIDStatus(AliPIDResponse::kTOF,track);
  if(statusTOF != AliPIDResponse::kDetPidOk)
    return kFALSE;
  if(0)
  {
    Int_t startTimeMask = fPIDResponse -> GetTOFResponse().GetStartTimeMask(track->P());
    if (startTimeMask < 0) return kFALSE; 
  }
  return kTRUE;
}

//___________________________________________________________
Double_t AliAnalysisTaskPidPidCorrelations ::GetBeta(AliAODTrack* track)
{
	//it is called only when TOF PID is available
	Double_t p	= track -> P();
	Double_t time	= track -> GetTOFsignal() - fPIDResponse -> GetTOFResponse().GetStartTime(p);
	Double_t timei[5];
	track -> GetIntegratedTimes(timei,AliPID::kSPECIES);
	return timei[0]/time;
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::RemoveDuplicates(TObjArray* tracks)
{
	// remove particles with the same label
  	Int_t before = tracks -> GetEntriesFast();

  	for (Int_t i=0; i<before; ++i) 
  	{
    		AliVParticle* part = (AliVParticle*) tracks -> At(i);
    
    		for (Int_t j=i+1; j<before; ++j) 
    		{
      	AliVParticle* part2 = (AliVParticle*) tracks -> At(j);
      
      	if (part -> GetLabel() == part2 -> GetLabel())
      	{
        		Printf("Removing %d with label %d (duplicated in %d)", i, part -> GetLabel(), j); part -> Dump(); part2 -> Dump();
        		TObject* object = tracks -> RemoveAt(i);
        		if (tracks -> IsOwner())
          		delete object;
        	break;
      	}
    		}
  	}
  	tracks -> Compress();
  
  	if (before > tracks -> GetEntriesFast())
    		AliInfo(Form("Reduced Form %d to %d", before, tracks -> GetEntriesFast())); 
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::RemoveWeakDecays(TObjArray* tracks, TObject* mcObj)
{
  // remove particles Form weak decays
  // <tracks> can be the following cases:
  // a. tracks: in this case the label is taken and then case b.
  // b. particles: it is checked if IsSecondaryFromWeakDecay is true
  // <mcObj> can be AOD (TClonesArray) or ESD (AliMCEvent)
  
  TClonesArray* arrayMC = 0;
  AliMCEvent* mcEvent = 0;
  if (mcObj -> InheritsFrom("AliMCEvent"))
    mcEvent = static_cast<AliMCEvent*>(mcObj);
  else if (mcObj -> InheritsFrom("TClonesArray"))
    arrayMC = static_cast<TClonesArray*>(mcObj);
  else
  {
    mcObj -> Dump();
    AliFatal("Invalid object passed");
  }
  
  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) 
  {
    AliVParticle* part = (AliVParticle*) tracks->At(i);
    
    if (part->InheritsFrom("AliESDtrack") || part->InheritsFrom("AliAODTrack"))
      part = ((mcEvent) ? mcEvent->GetTrack(TMath::Abs(part->GetLabel())) : (AliVParticle*)arrayMC->At(TMath::Abs(part->GetLabel())));
    
    if (part->InheritsFrom("AliAODMCParticle"))
    {
      if (!((AliAODMCParticle*) part)->IsSecondaryFromWeakDecay())
		continue;
    }
    else if (part->InheritsFrom("AliMCParticle") && mcEvent)
    {
      if (!(mcEvent->Stack()->IsSecondaryFromWeakDecay(part->GetLabel())))
		continue;
    }
    else
    {
      part -> Dump();
      AliFatal("Unknown particle");
    }
    
//     Printf("Removing %d with label %d", i, part->GetLabel()); part->Dump();
    TObject* object = tracks -> RemoveAt(i);
    if (tracks -> IsOwner())
      delete object;
  }
 
  tracks -> Compress();
  
  if (before > tracks->GetEntriesFast())
    AliInfo(Form("Reduced Form %d to %d", before, tracks->GetEntriesFast())); 
}

//________________________________________________________________________
Double_t AliAnalysisTaskPidPidCorrelations::DeltaPhi(Double_t Dphi) const
{
  if (Dphi < -0.5*TMath::Pi())	Dphi += TMath::TwoPi();
  if (Dphi > 1.5*TMath::Pi())	Dphi -= TMath::TwoPi();
  return Dphi;  
}

//___________________________________________________________
TH2F* AliAnalysisTaskPidPidCorrelations::GetHisto2D(const Char_t* name)
{
  return (TH2F*) fList -> FindObject(name);
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::GetCentBin(Double_t percentile)
{
//   Double_t percentile = fMyAODEvent -> GetCentrality() -> GetCentralityPercentile(fCentralityEstimator.Data());
  Int_t bin = fCentAxis -> FindBin(percentile) - 1;
  if (bin >= fNbinsCent) bin = -1;
  return bin;
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::SetCentBinning(Int_t nBins, Double_t* limits)
{
  if (nBins > fNMaxBinsCentrality)
  {
    AliInfo(Form("WARNING : only %d centrality bins (out of the %d proposed) will be considered",fNMaxBinsCentrality,nBins));
    nBins = fNMaxBinsCentrality;
  }
  if (nBins <= 0)
  {
    AliInfo("WARNING : at least one centrality bin must be considered");
    nBins = 1;
  }
  fNbinsCent = nBins;
  fCentAxis  = new TAxis(fNbinsCent, limits);
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::GetZvtxBin(Double_t zvtx)
{
  Int_t bin = fZvtxAxis -> FindBin(zvtx) - 1;
  if (bin >= fNbinsZvtx) bin = -1;
  return bin;
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::SetZvtxBinning(Int_t nBins, Double_t* limits)
{
  if (nBins > fNMaxBinsZvtx)
  {
    AliInfo(Form("WARNING : only %d z-vertex bins (out of the %d proposed) will be considered",fNMaxBinsZvtx,nBins));
    nBins = fNMaxBinsZvtx;
  }
  if (nBins <= 0)
  {
    AliInfo("WARNING : at least one centrality bin must be considered");
    nBins = 1;
  }
  fNbinsZvtx = nBins;
  fZvtxAxis  = new TAxis(fNbinsZvtx, limits);
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::GetPtBin(Double_t valPt)
{
  Int_t bin = fPtAxis -> FindBin(valPt) - 1;
  if (bin >= fNbinsPt) bin = -1;
  return bin;
}

//__________________________________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::SetPtBinning(Int_t nBins, Double_t* limits)
{
  if (nBins > fNMaxBinsPt)
  {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsPt,nBins));
    nBins = fNMaxBinsPt;
  }
  if (nBins <= 0)
  {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  fNbinsPt = nBins;
  fPtAxis  = new TAxis(fNbinsPt, limits);
}

//___________________________________________________________
Int_t AliAnalysisTaskPidPidCorrelations::GetEtaBin(Double_t valEta)
{
  Int_t bin = fEtaAxis -> FindBin(valEta) - 1;
  if (bin >= fNbinsEta) bin = -1;
  return bin;
}

//__________________________________________________________________________________________________
void AliAnalysisTaskPidPidCorrelations::SetEtaBinning(Int_t nBins, Double_t* limits)
{
  if (nBins > fNMaxBinsEta)
  {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsEta,nBins));
    nBins = fNMaxBinsEta;
  }
  if (nBins <= 0)
  {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  fNbinsEta = nBins;
  fEtaAxis  = new TAxis(nBins, limits);
}

//________________________________________________________________________
Double_t AliAnalysisTaskPidPidCorrelations::GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius, Double_t bSign)
{ 
  //
  // calculates dphistar
  //
  
  Double_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}

//________________________________________________________________________ copy&paste from $ALICE_ROOT/PWGCF/Correlations/Base/AliUEHistograms.h
Float_t AliAnalysisTaskPidPidCorrelations::GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared
  // same can be achieved, but with more computing time with
  /*TLorentzVector photon, p1, p2;
  p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
  p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
  photon = p1+p2;
  photon.M()*/
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = TMath::Exp(-eta1);
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = TMath::Exp(-eta2);
    tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2 ) ) );
  
//   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
  
  return mass2;
}

//________________________________________________________________________ copy&paste from $ALICE_ROOT/PWGCF/Correlations/Base/AliUEHistograms.h
Float_t AliAnalysisTaskPidPidCorrelations::GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared approximately
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = 1.0-eta2+eta2*eta2/2-eta2*eta2*eta2/6+eta2*eta2*eta2*eta2/24;
    tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  // fold onto 0...pi
  Float_t deltaPhi = TMath::Abs(phi1 - phi2);
  while (deltaPhi > TMath::TwoPi())
    deltaPhi -= TMath::TwoPi();
  if (deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  
  Float_t cosDeltaPhi = 0;
  if (deltaPhi < TMath::Pi()/3)
    cosDeltaPhi = 1.0 - deltaPhi*deltaPhi/2 + deltaPhi*deltaPhi*deltaPhi*deltaPhi/24;
  else if (deltaPhi < 2*TMath::Pi()/3)
    cosDeltaPhi = -(deltaPhi - TMath::Pi()/2) + 1.0/6 * TMath::Power((deltaPhi - TMath::Pi()/2), 3);
  else
    cosDeltaPhi = -1.0 + 1.0/2.0*(deltaPhi - TMath::Pi())*(deltaPhi - TMath::Pi()) - 1.0/24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( cosDeltaPhi + 1.0 / tantheta1 / tantheta2 ) ) );
  
//   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
  
  return mass2;
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::PrintPoolManagerContents()
{
  // Determine the number of pools in the pool manager.
  AliEventPool* poolin = fPoolMgr -> GetEventPool(0,0);
  Int_t NPoolsCentrality = 0;
  while (poolin)
  {
    NPoolsCentrality++;
    poolin = fPoolMgr -> GetEventPool(NPoolsCentrality,0);
  } 

  poolin = fPoolMgr -> GetEventPool(0,0);
  Int_t NPoolsVtxZ = 0;
  while (poolin)
  {
    NPoolsVtxZ++;
    poolin = fPoolMgr -> GetEventPool(0,NPoolsVtxZ);
  } 

  // Loop over all Pools in the matrix of the pool manager.
  std::cout<<" Pool manager contents: (Nevt,NTrack)"<<std::endl;
  for (Int_t iCentrality = 0; iCentrality < NPoolsCentrality; iCentrality++)
  {
    std::cout<<Form("Centrality Bin: %2i --> ", iCentrality);

    for (Int_t iVtxZ = 0; iVtxZ < NPoolsVtxZ; iVtxZ++)
    {
      poolin = fPoolMgr->GetEventPool(iCentrality, iVtxZ);
      std::cout<<Form("(%2i,%4i) ",poolin->GetCurrentNEvents(), poolin->NTracksInPool());
      
//       Double_t vars[] = { iCentrality, iVtxZ, poolin->GetCurrentNEvents(), poolin->NTracksInPool() };
//       fHistPoolMgrQA -> Fill(vars);
    }
    std::cout<<std::endl;
  }
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::SetupForMixing()
{
  const Int_t trackDepth 	= fMinNumTrack;
  const Int_t poolsize 		= fPoolSize;
  Double_t centralityBins[] = { fCentralityPercentileMin, fCentralityPercentileMax };
  const Int_t nCentralityBins = length(centralityBins)-1;
  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins/*fCentAxis->GetNbins()+1*/, centralityBins/*const_cast<Double_t*>(fCentAxis->GetXbins()->GetArray())*/, fZvtxAxis->GetNbins()+1, const_cast<Double_t*>(fZvtxAxis->GetXbins()->GetArray()));
//   fPoolMgr -> SetTargetValues(fMinNumTrack, 0.1, 5);
}

//___________________________________________________________
TObjArray* AliAnalysisTaskPidPidCorrelations::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliPidPidCorrelationReducedTrack which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone -> SetOwner(kTRUE);
  
  for (Int_t i=0; i < tracks->GetEntriesFast(); i++)
  {
    AliPidPidCorrelationReducedTrack* particle = (AliPidPidCorrelationReducedTrack*) tracks->UncheckedAt(i);
    AliPidPidCorrelationReducedTrack* copy = new AliPidPidCorrelationReducedTrack(particle->GetMyPartID(),particle->Eta(),particle->Phi(),particle->Pt(),particle->Charge());
    copy -> SetUniqueID(particle->GetUniqueID());
    tracksClone -> Add(copy);
  }
  return tracksClone;
}
//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::FillCFcontainers(TObjArray* mc, TObjArray* reco, TObjArray* recomatch, Double_t cent/*, Double_t zvtx*//*const AliVEvent* event*/)
{
  if (!fMyCFCont) return;
/*
  Double_t zvtx=0., zvtxMC=0.;
  if ( event->IsA() == AliAODEvent::Class() )
  {
    const AliVVertex* vertex = (AliVVertex*)static_cast<const AliAODEvent*>(event)->GetPrimaryVertexSPD();
    zvtx = vertex->GetZ();
  }
  if ( event->IsA() == AliAODEvent::Class() )
  {
    AliAODMCHeader* aodMCHeader = static_cast<AliAODMCHeader*> (static_cast<const AliAODEvent*>(event)->FindListObject(AliAODMCHeader::StdBranchName()));
    zvtxMC = aodMCHeader->GetVtxZ();
  }
  */
  for (Int_t istep=0; istep < kNsteps; istep++)
  {
    TObjArray* list = mc;
    if (istep == 1)
      list = reco;
    else if (istep == 2)
      list = recomatch;
    
    if (!list)
      continue;

    Double_t cfValue[kNvars];
    for(Int_t ie = 0 ; ie < list->GetEntriesFast(); ie++)
    {
      AliPidPidCorrelationReducedTrack* mytrack = dynamic_cast<AliPidPidCorrelationReducedTrack*>(list->At(ie));
  
      cfValue[kVarEta]	= mytrack->Eta();
//       cfValue[kVarPhi]	= mytrack->Phi();
      cfValue[kVarPt]	= mytrack->Pt();
      cfValue[kVarPID]	= mytrack->GetMyPartID();
      cfValue[kVarCent]	= cent;
      cfValue[kVarZvtx]	= ( istep == kStepRec ) ? (fMyAODEvent->GetPrimaryVertex()->GetZ()) : (fMyMcHeader->GetVtxZ()); //! FIXME
//       cfValue[kVarZvtx]	= zvtx;
      
      fMyCFCont->Fill(cfValue,istep);
    }
  }
}

//___________________________________________________________
TString AliAnalysisTaskPidPidCorrelations::GetOutputListName() const
{
  TString listName("listPIDCorr");
  listName += TString::Format("_%smix",         fFillMixed ? "" : "no");
  listName += TString::Format("_cent%.0f%.0f", fCentralityPercentileMin, fCentralityPercentileMax);
  listName += TString::Format("_ptMin%.0fMeV",  1e3*fTrackPtMin);
  return listName;
}

//___________________________________________________________
void AliAnalysisTaskPidPidCorrelations::Terminate(Option_t *) 
{
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();

  fList   = dynamic_cast<TList*>(GetOutputData(1));
  if(!fList) { Printf("ERROR: could not retrieve TList fList"); return; }

  fMyCFCont = dynamic_cast<AliCFContainer*> (GetOutputData(2));
  if ( !fMyCFCont ) { AliError("Cannot find container in file"); return; }  
  
  if(fPoolMgr) delete fPoolMgr;

  return;
}

