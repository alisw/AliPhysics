/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, proviyaded that the above copyright notice appears in all *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is         *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//                        Analysis for Jet Hadrochemistry                                 //
//                                                                                        //
//    This analysis extracts pT-spectra of charged kaons, protons, and pions              //
//                      for the inclusive event and in jets.                              //
//   It is based on particles identification via the dE/dx signal of the TPC              //
//                   and the time of flight nsigma from the TOF.                          //
//                              This is the ESD version.                                  //
//                                                                                        //
// Author: Sierra Cantway (Weyhmiller) <sierra.lisa.weyhmiller@cern.ch>, Yale University  //
//      Author: Mesut Arslandok <mesut.arslandok@cern.ch>, Yale University                //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TVectorF.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliPDG.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliLumiTools.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliVParticle.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliAnalysisJetHadro.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliFJWrapper.h"
#include "AliEmcalJet.h"
#include "AliEmcalJetTask.h"
#include "AliRhoParameter.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
#include <TObjString.h>
using namespace std;
using std::cout;
using std::setw;

ClassImp(AliAnalysisJetHadro)

#define USE_STREAMER 1


// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisJetHadro::AliAnalysisJetHadro()
: AliAnalysisTaskEmcalJet("JetHadro"), fEventCuts(0), fPIDResponse(0), fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_2015(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCuts_Bit768_v(0),
fPIDCombined(0x0),
fMCStack(0x0),
fVertex(0x0),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreeMC(0x0),
fTreeCuts(0x0),
fTreejetsEMCBGconst(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJBGconst(0x0),
fTreejetEvents(0x0),
fTreejetsEMC(0x0),
fTreejetsEMCBG(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPercentageOfEvents(0),
fSmallOut(kFALSE),
fMCtrue(kFALSE),
fIncludeITS(kTRUE),
fFilljetsFJBGTree(kTRUE),
fFillJetsFJBGConst(kTRUE),
fDoIncTracks(kTRUE),
fDoFastJet(kTRUE),
fDoEMCJet(kTRUE),
fFillFastJet(kTRUE),
fFillJetsFJConst(kTRUE),
fFillEMCJet(kTRUE),
fFillIncTracks(kTRUE),
fFill_TPC(kTRUE),
fFillpTPC_pT(kTRUE),
fFillp_pT(kTRUE),
fFillpTPC_p(kTRUE),
fFill_TOF(kTRUE),
fFill_TOF_expecs(kTRUE),
fFill_TPC_expecs(kTRUE),
fFillJetsEMCConst(kTRUE),
fFillJetsEMCBGConst(kTRUE),
fcent_min(0.0),
fcent_max(100.0),
fjetMinPtSub(-1000.0),
fjetMinArea(-1000.0),
fFillTreeMC(kFALSE),
fUseCouts(kFALSE),
fSetTPCmom(0),
fSetTOFmom(0),
fSetBetamom(0),
fSetEta(0),
fNTPCMom_Bins(2000),
fTPCMom_Bins(0),
fTOFMom_NBins(1000),
fTOFMom_Bins(0),
fNEta_Bins(9),
fEta_Bins(0),
fNdEdxBins(2000),
fdEdxBins(0),
fNBetaBins(1000),
fBetaBins(0),
fNTOFNSigmaBins(1000),
fTOFNSigmaBins(0),
fNSigmasElTOF(0),
fNSigmasMuTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventCountInFile(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fPhi(0),
fSign(0),
fJetContainer(0),
fbgJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
frhoFJ(0),
fisGoodIncEvent(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
fNumRealFJJets(0),
fNumRealEMCJets(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fHistCentrality(0),
fHistImpParam(0),
fHistVertex(0),
fHistIncTracks_dEdx(0),
fHistIncTracks_moms(0),
fHistIncTracks_moms_p(0),
fHistIncTracks_moms_pTPC_p(0),
fHistIncTracks_kin(0),
fHistIncTracks_beta(0),
fHistIncTracks_t0(0),
fHistIncTracks_TOFpi_nsigma(0),
fHistIncTracks_TOFka_nsigma(0),
fHistIncTracks_TOFpr_nsigma(0),
fHistIncTracks_TOFpi_nsigma_1cls(0),
fHistIncTracks_TOFka_nsigma_1cls(0),
fHistIncTracks_TOFpr_nsigma_1cls(0),
fHistJetTracks_dEdx(0),
fHistJetTracks_moms(0),
fHistJetTracks_moms_p(0),
fHistJetTracks_moms_pTPC_p(0),
fHistJetTracks_kin(0),
fHistJetTracks_beta(0),
fHistJetTracks_TOFpi_nsigma(0),
fHistJetTracks_TOFka_nsigma(0),
fHistJetTracks_TOFpr_nsigma(0),
fHistJetTracks_TOFpi_nsigma_1cls(0),
fHistJetTracks_TOFka_nsigma_1cls(0),
fHistJetTracks_TOFpr_nsigma_1cls(0),
fHistIncTracks_mpi(0),
fHistBetaExpec_pi(0),
fHistBetaExpec_ka(0),
fHistBetaExpec_pr(0),
fHistjet_BetaExpec_pi(0),
fHistjet_BetaExpec_ka(0),
fHistjet_BetaExpec_pr(0),
fHist_pi_mismatch(0),
fHist_ka_mismatch(0),
fHist_pr_mismatch(0),
fHist_jet_pi_mismatch(0),
fHist_jet_ka_mismatch(0),
fHist_jet_pr_mismatch(0),
fHist_elExpec_pihyp(0),
fHist_muExpec_pihyp(0),
fHist_kaExpec_pihyp(0),
fHist_prExpec_pihyp(0),
fHist_piExpec_kahyp(0),
fHist_prExpec_kahyp(0),
fHist_piExpec_prhyp(0),
fHist_kaExpec_prhyp(0),
fHist_deExpec_prhyp(0),
fHist_jet_elExpec_pihyp(0),
fHist_jet_muExpec_pihyp(0),
fHist_jet_kaExpec_pihyp(0),
fHist_jet_prExpec_pihyp(0),
fHist_jet_piExpec_kahyp(0),
fHist_jet_prExpec_kahyp(0),
fHist_jet_piExpec_prhyp(0),
fHist_jet_kaExpec_prhyp(0),
fHist_jet_deExpec_prhyp(0),
fHistTOFSigmaExpec_pi(0),
fHistTOFSigmaExpec_ka(0),
fHistTOFSigmaExpec_pr(0),
fHistjet_TOFSigmaExpec_pi(0),
fHistjet_TOFSigmaExpec_ka(0),
fHistjet_TOFSigmaExpec_pr(0),
fHistIncTracks_spi(0),
fHistIncTracks_mel(0),
fHistIncTracks_sel(0),
fHistIncTracks_mka(0),
fHistIncTracks_ska(0),
fHistIncTracks_mpr(0),
fHistIncTracks_spr(0),
fHistJet_ptsub_v_area(0),
fHistJet_kin(0),
fHistJet_moms(0)
{
  // default Constructor
}

//________________________________________________________________________
AliAnalysisJetHadro::AliAnalysisJetHadro(const char *name)
: AliAnalysisTaskEmcalJet(name), fEventCuts(0), fPIDResponse(0), fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_2015(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCuts_Bit768_v(0),
fPIDCombined(0x0),
fMCStack(0x0),
fVertex(0x0),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreejetsEMCBGconst(0x0),
fTreeMC(0x0),
fTreeCuts(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJBGconst(0x0),
fTreejetEvents(0x0),
fTreejetsEMC(0x0),
fTreejetsEMCBG(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPercentageOfEvents(0),
fSmallOut(kFALSE),
fMCtrue(kFALSE),
fIncludeITS(kTRUE),
fFilljetsFJBGTree(kTRUE),
fFillJetsFJBGConst(kTRUE),
fDoIncTracks(kTRUE),
fDoFastJet(kTRUE),
fDoEMCJet(kTRUE),
fFillFastJet(kTRUE),
fFillJetsFJConst(kTRUE),
fFillEMCJet(kTRUE),
fFillJetsEMCConst(kTRUE),
fFillJetsEMCBGConst(kTRUE),
fFillIncTracks(kTRUE),
fFill_TPC(kTRUE),
fFillpTPC_pT(kTRUE),
fFillp_pT(kTRUE),
fFillpTPC_p(kTRUE),
fFill_TOF(kTRUE),
fFill_TOF_expecs(kTRUE),
fFill_TPC_expecs(kTRUE),
fcent_min(0.0),
fcent_max(0.0),
fjetMinPtSub(-1000.0),
fjetMinArea(-1000.0),
fFillTreeMC(kFALSE),
fUseCouts(kFALSE),
fSetTPCmom(0),
fSetTOFmom(0),
fSetBetamom(0),
fSetEta(0),
fNTPCMom_Bins(2000),
fTPCMom_Bins(0),
fTOFMom_NBins(1000),
fTOFMom_Bins(0),
fNEta_Bins(9),
fEta_Bins(0),
fNdEdxBins(2000),
fdEdxBins(0),
fNBetaBins(1000),
fBetaBins(0),
fNTOFNSigmaBins(1000),
fTOFNSigmaBins(0),
fNSigmasElTOF(0),
fNSigmasMuTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventCountInFile(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fPhi(0),
fSign(0),
fJetContainer(0),
fbgJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
frhoFJ(0),
fisGoodIncEvent(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
fNumRealFJJets(0),
fNumRealEMCJets(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fHistCentrality(0),
fHistImpParam(0),
fHistVertex(0),
fHistIncTracks_dEdx(0),
fHistIncTracks_moms(0),
fHistIncTracks_moms_p(0),
fHistIncTracks_moms_pTPC_p(0),
fHistIncTracks_kin(0),
fHistIncTracks_beta(0),
fHistIncTracks_t0(0),
fHistIncTracks_TOFpi_nsigma(0),
fHistIncTracks_TOFka_nsigma(0),
fHistIncTracks_TOFpr_nsigma(0),
fHistIncTracks_TOFpi_nsigma_1cls(0),
fHistIncTracks_TOFka_nsigma_1cls(0),
fHistIncTracks_TOFpr_nsigma_1cls(0),
fHistJetTracks_dEdx(0),
fHistJetTracks_moms(0),
fHistJetTracks_moms_p(0),
fHistJetTracks_moms_pTPC_p(0),
fHistJetTracks_kin(0),
fHistJetTracks_beta(0),
fHistJetTracks_TOFpi_nsigma(0),
fHistJetTracks_TOFka_nsigma(0),
fHistJetTracks_TOFpr_nsigma(0),
fHistJetTracks_TOFpi_nsigma_1cls(0),
fHistJetTracks_TOFka_nsigma_1cls(0),
fHistJetTracks_TOFpr_nsigma_1cls(0),
fHistIncTracks_mpi(0),
fHistBetaExpec_pi(0),
fHistBetaExpec_ka(0),
fHistBetaExpec_pr(0),
fHistjet_BetaExpec_pi(0),
fHistjet_BetaExpec_ka(0),
fHistjet_BetaExpec_pr(0),
fHist_pi_mismatch(0),
fHist_ka_mismatch(0),
fHist_pr_mismatch(0),
fHist_jet_pi_mismatch(0),
fHist_jet_ka_mismatch(0),
fHist_jet_pr_mismatch(0),
fHist_elExpec_pihyp(0),
fHist_muExpec_pihyp(0),
fHist_kaExpec_pihyp(0),
fHist_prExpec_pihyp(0),
fHist_piExpec_kahyp(0),
fHist_prExpec_kahyp(0),
fHist_piExpec_prhyp(0),
fHist_kaExpec_prhyp(0),
fHist_deExpec_prhyp(0),
fHist_jet_elExpec_pihyp(0),
fHist_jet_muExpec_pihyp(0),
fHist_jet_kaExpec_pihyp(0),
fHist_jet_prExpec_pihyp(0),
fHist_jet_piExpec_kahyp(0),
fHist_jet_prExpec_kahyp(0),
fHist_jet_piExpec_prhyp(0),
fHist_jet_kaExpec_prhyp(0),
fHist_jet_deExpec_prhyp(0),
fHistTOFSigmaExpec_pi(0),
fHistTOFSigmaExpec_ka(0),
fHistTOFSigmaExpec_pr(0),
fHistjet_TOFSigmaExpec_pi(0),
fHistjet_TOFSigmaExpec_ka(0),
fHistjet_TOFSigmaExpec_pr(0),
fHistIncTracks_spi(0),
fHistIncTracks_mel(0),
fHistIncTracks_sel(0),
fHistIncTracks_mka(0),
fHistIncTracks_ska(0),
fHistIncTracks_mpr(0),
fHistIncTracks_spr(0),
fHistJet_ptsub_v_area(0),
fHistJet_kin(0),
fHistJet_moms(0)
{
  //
  //         standard constructur which should be used
  //
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:***************** CONSTRUCTOR CALLED: AliAnalysisJetHadro  *****************"<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  // ==========================================
  //
  // ==========================================
  //
  // Define outputs
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());
  DefineOutput(10, TTree::Class());
  DefineOutput(11, TTree::Class());
  DefineOutput(12, TTree::Class());
  // ==========================================

}
//________________________________________________________________________
AliAnalysisJetHadro::~AliAnalysisJetHadro()
{

  //
  // Destructor
  //
  std::cout << " Info::siweyhmi: ===== In the Destructor ===== " << std::endl;
  if (fHistCentrality)      delete fHistCentrality;
  if (fHistImpParam)        delete fHistImpParam;
  if (fHistVertex)          delete fHistVertex;
  if (fHistIncTracks_dEdx)          delete fHistIncTracks_dEdx;
  if (fHistIncTracks_moms)          delete fHistIncTracks_moms;
  if (fHistIncTracks_moms_p)          delete fHistIncTracks_moms_p;
  if (fHistIncTracks_moms_pTPC_p)          delete fHistIncTracks_moms_pTPC_p;
  if (fHistIncTracks_kin)          delete fHistIncTracks_kin;
  if (fHistIncTracks_beta)          delete fHistIncTracks_beta;
  if (fHistJetTracks_beta)          delete fHistJetTracks_beta;
  if (fHistIncTracks_TOFpi_nsigma)          delete fHistIncTracks_TOFpi_nsigma;
  if (fHistIncTracks_TOFka_nsigma)          delete fHistIncTracks_TOFka_nsigma;
  if (fHistIncTracks_TOFpr_nsigma)          delete fHistIncTracks_TOFpr_nsigma;
  if (fHistIncTracks_TOFpi_nsigma_1cls)          delete fHistIncTracks_TOFpi_nsigma_1cls;
  if (fHistIncTracks_TOFka_nsigma_1cls)          delete fHistIncTracks_TOFka_nsigma_1cls;
  if (fHistIncTracks_TOFpr_nsigma_1cls)          delete fHistIncTracks_TOFpr_nsigma_1cls;
  if (fHistJetTracks_dEdx)          delete fHistJetTracks_dEdx;
  if (fHistJetTracks_moms)          delete fHistJetTracks_moms;
  if (fHistJetTracks_moms_p)          delete fHistJetTracks_moms_p;
  if (fHistJetTracks_moms_pTPC_p)          delete fHistJetTracks_moms_pTPC_p;
  if (fHistJetTracks_kin)          delete fHistJetTracks_kin;
  if (fHistIncTracks_t0)          delete fHistIncTracks_t0;
  if (fHistJetTracks_TOFpi_nsigma)          delete fHistJetTracks_TOFpi_nsigma;
  if (fHistJetTracks_TOFka_nsigma)          delete fHistJetTracks_TOFka_nsigma;
  if (fHistJetTracks_TOFpr_nsigma)          delete fHistJetTracks_TOFpr_nsigma;
  if (fHistJetTracks_TOFpi_nsigma_1cls)          delete fHistJetTracks_TOFpi_nsigma_1cls;
  if (fHistJetTracks_TOFka_nsigma_1cls)          delete fHistJetTracks_TOFka_nsigma_1cls;
  if (fHistJetTracks_TOFpr_nsigma_1cls)          delete fHistJetTracks_TOFpr_nsigma_1cls;
  if (fHistBetaExpec_pi) delete fHistBetaExpec_pi;
  if (fHistBetaExpec_ka) delete fHistBetaExpec_ka;
  if (fHistBetaExpec_pr) delete fHistBetaExpec_pr;
  if (fHistjet_BetaExpec_pi) delete fHistjet_BetaExpec_pi;
  if (fHistjet_BetaExpec_ka) delete fHistjet_BetaExpec_ka;
  if (fHistjet_BetaExpec_pr) delete fHistjet_BetaExpec_pr;
  if (fHist_pi_mismatch) delete fHist_pi_mismatch;
  if (fHist_ka_mismatch) delete fHist_ka_mismatch;
  if (fHist_pr_mismatch) delete fHist_pr_mismatch;
  if (fHist_jet_pi_mismatch) delete fHist_jet_pi_mismatch;
  if (fHist_jet_ka_mismatch) delete fHist_jet_ka_mismatch;
  if (fHist_jet_pr_mismatch) delete fHist_jet_pr_mismatch;
  if (fHist_elExpec_pihyp) delete fHist_elExpec_pihyp;
  if (fHist_muExpec_pihyp) delete fHist_muExpec_pihyp;
  if (fHist_kaExpec_pihyp) delete fHist_kaExpec_pihyp;
  if (fHist_prExpec_pihyp) delete fHist_prExpec_pihyp;
  if (fHist_piExpec_kahyp) delete fHist_piExpec_kahyp;
  if (fHist_prExpec_kahyp) delete fHist_prExpec_kahyp;
  if (fHist_piExpec_prhyp) delete fHist_piExpec_prhyp;
  if (fHist_kaExpec_prhyp) delete fHist_kaExpec_prhyp;
  if (fHist_deExpec_prhyp) delete fHist_deExpec_prhyp;
  if (fHist_jet_elExpec_pihyp) delete fHist_jet_elExpec_pihyp;
  if (fHist_jet_muExpec_pihyp) delete fHist_jet_muExpec_pihyp;
  if (fHist_jet_kaExpec_pihyp) delete fHist_jet_kaExpec_pihyp;
  if (fHist_jet_prExpec_pihyp) delete fHist_jet_prExpec_pihyp;
  if (fHist_jet_piExpec_kahyp) delete fHist_jet_piExpec_kahyp;
  if (fHist_jet_prExpec_kahyp) delete fHist_jet_prExpec_kahyp;
  if (fHist_jet_piExpec_prhyp) delete fHist_jet_piExpec_prhyp;
  if (fHist_jet_kaExpec_prhyp) delete fHist_jet_kaExpec_prhyp;
  if (fHist_jet_deExpec_prhyp) delete fHist_jet_deExpec_prhyp;
  if (fHistTOFSigmaExpec_pi) delete fHistTOFSigmaExpec_pi;
  if (fHistTOFSigmaExpec_ka) delete fHistTOFSigmaExpec_ka;
  if (fHistTOFSigmaExpec_pr) delete fHistTOFSigmaExpec_pr;
  if (fHistjet_TOFSigmaExpec_pi) delete fHistjet_TOFSigmaExpec_pi;
  if (fHistjet_TOFSigmaExpec_ka) delete fHistjet_TOFSigmaExpec_ka;
  if (fHistjet_TOFSigmaExpec_pr) delete fHistjet_TOFSigmaExpec_pr;
  if (fHistIncTracks_mpi)          delete fHistIncTracks_mpi;
  if (fHistIncTracks_spi)          delete fHistIncTracks_spi;
  if (fHistIncTracks_mel)          delete fHistIncTracks_mel;
  if (fHistIncTracks_sel)          delete fHistIncTracks_sel;
  if (fHistIncTracks_mka)          delete fHistIncTracks_mka;
  if (fHistIncTracks_ska)          delete fHistIncTracks_ska;
  if (fHistIncTracks_mpr)          delete fHistIncTracks_mpr;
  if (fHistIncTracks_spr)          delete fHistIncTracks_spr;
  if (fHistJet_ptsub_v_area)          delete fHistJet_ptsub_v_area;
  if (fHistJet_kin)          delete fHistJet_kin;
  if (fHistJet_moms)          delete fHistJet_moms;
  if (fPIDCombined) delete fPIDCombined;
  if (fESDtrackCuts)          delete fESDtrackCuts;
  if (fESDtrackCuts_Bit128)   delete fESDtrackCuts_Bit128;
  if (fESDtrackCuts_Bit768)   delete fESDtrackCuts_Bit768;
  if (fESDtrackCuts_Bit768_v)   delete fESDtrackCuts_Bit768_v;
  if (fESDtrackCuts_2015)   delete fESDtrackCuts_2015;
  if (fTreeSRedirector)       delete fTreeSRedirector;

}
//
// ---------------------------------------------------------------------------------
//                                     Functions
// ---------------------------------------------------------------------------------
//
void AliAnalysisJetHadro::Initialize()
{
  //
  std::cout << " Info::siweyhmi: ===== In the Initialize ===== " << std::endl;
  //
  // ------------------------------------------------
  //
  //
  // TPC only tracks
  fESDtrackCuts_Bit128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fESDtrackCuts_Bit128->SetName("Bit128");
  fESDtrackCuts_Bit128->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit128->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit128->SetMinNClustersTPC(70);
  fESDtrackCuts_Bit128->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit128->SetDCAToVertex2D(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(4);
  }

  //Hybrid track cuts
  fESDtrackCuts_Bit768 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCuts_Bit768->SetName("Bit768");
  fESDtrackCuts_Bit768->SetEtaRange(-0.9,0.9);
  fESDtrackCuts_Bit768->SetPtRange(0.15,1000.);

  fESDtrackCuts_Bit768->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit768->SetDCAToVertex2D(kTRUE);

  fESDtrackCuts_Bit768->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit768->SetMaxFractionSharedTPCClusters(0.4);

  fESDtrackCuts_Bit768->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

  fESDtrackCuts_Bit768_v = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCuts_Bit768_v->SetName("Bit768_v");
  fESDtrackCuts_Bit768_v->SetEtaRange(-0.9,0.9);
  fESDtrackCuts_Bit768_v->SetPtRange(0.15,1000.);

  fESDtrackCuts_Bit768_v->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit768_v->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit768_v->SetDCAToVertex2D(kTRUE);

  fESDtrackCuts_Bit768_v->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit768_v->SetMaxFractionSharedTPCClusters(0.4);

  fESDtrackCuts_Bit768_v->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  // Default track cuts
  fESDtrackCuts = new AliESDtrackCuts("esdTrackCuts","");
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  fESDtrackCuts->SetPtRange(0.15,1000.);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNCrossedRowsTPC(80);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");
  fESDtrackCuts->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetDCAToVertex2D(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if (fIncludeITS) {
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  }

  //TOF inc pi,k,p 2015 spectrum track cuts
  fESDtrackCuts_2015 = new AliESDtrackCuts("fESDtrackCuts_2015","");
  fESDtrackCuts_2015 = (AliESDtrackCuts*) AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(kFALSE, 1);
  fESDtrackCuts_2015->SetEtaRange(-0.9,0.9);
  fESDtrackCuts_2015->SetPtRange(0.15,1000.);

  //
  //
  std::cout << " Info::siweyhmi: ===================================================== " << std::endl;
  std::cout << " Info::siweyhmi: =============== Summary of Track Cuts =============== " << std::endl;
  std::cout << " Info::siweyhmi: ===================================================== " << std::endl;
}
//________________________________________________________________________
void AliAnalysisJetHadro::UserCreateOutputObjects()
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  Initialize();
  std::cout << " Info::siweyhmi: ===== In the UserCreateOutputObjects ===== " << std::endl;
  // ------------  setup PIDCombined  ---------------
  fPIDCombined = new AliPIDCombined("pidCombined","");
  //
  // **********************   Input handler to get the PID object *********************
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler)
  AliFatal("Input handler needed");
  else {
    fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
    if (!fPIDResponse) std::cout << " Info::siweyhmi: ======= PIDResponse object was not created ====== " << std::endl;
  }
  //
  // ************************************************************************
  //   OpenFile output --> one can open several files
  // ************************************************************************
  //
  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  fEventCuts.AddQAplotsToList(fListHist);

  Float_t tpc_mom_bins[58] = {};
  for (int i=0; i < (57+1); i++){
    if (i < 20) tpc_mom_bins[i] = 0.0 + i*0.05;
    else if (i < 30) tpc_mom_bins[i] = 1.0 + (i-20)*0.1;
    else if (i < 40) tpc_mom_bins[i] = 2.0 + (i-30)*0.2;
    else if (i < 46) tpc_mom_bins[i] = 4.0 + (i-40)*0.5;
    else if (i < 55) tpc_mom_bins[i] = 7.0 + (i-46)*1.0;
    else tpc_mom_bins[i] = 16.0 + (i-55)*2.0;
  }
  Int_t n_tpc_mom_bins = sizeof(tpc_mom_bins)/sizeof(Float_t) - 1;

  Float_t TOF_mom_bins[48] = {};
  for (int i=0; i < (47+1); i++){
    if (i < 20) TOF_mom_bins[i] = 0.0 + i*0.05;
    else if (i < 30) TOF_mom_bins[i] = 1.0 + (i-20)*0.1;
    else if (i < 40) TOF_mom_bins[i] = 2.0 + (i-30)*0.2;
    else if (i < 46) TOF_mom_bins[i] = 4.0 + (i-40)*0.5;
    else TOF_mom_bins[i] = 7.0 + (i-46)*1.0;
  }
  Int_t n_TOF_mom_bins = sizeof(TOF_mom_bins)/sizeof(Float_t) - 1;

  Float_t dedx_bins[2001] = {};
  for (int i=0; i < (2000+1); i++){
    dedx_bins[i] = 0.0 + i*1.0;
  }
  Int_t n_dedx_bins = sizeof(dedx_bins)/sizeof(Float_t) - 1;

  Float_t beta_bins[1001] = {};
  for (int i=0; i < (1000+1); i++){
    beta_bins[i] = 0.0 + i*0.002;
  }
  Int_t n_beta_bins = sizeof(beta_bins)/sizeof(Float_t) - 1;

  Float_t tof_n_sigma_bins[2001] = {};
  for (int i=0; i < (2000+1); i++){
    tof_n_sigma_bins[i] = -100.0 + i*0.1;
  }
  Int_t n_tof_n_sigma_bins = sizeof(tof_n_sigma_bins)/sizeof(Float_t) - 1;

  Float_t eta_bins[5] = {};
  for (int i=0; i < (4+1); i++){
    if (i < 3) eta_bins[i] = 0.0 + i*0.2;
    else eta_bins[i] = 0.6 + (i-3)*0.3;
  }
  Int_t n_eta_bins = sizeof(eta_bins)/sizeof(Float_t) - 1;

  //
  // ************************************************************************
  //   Event histograms
  // ************************************************************************
  //
  fHistCentrality        = new TH1F("hCentrality",           "control histogram for centrality"           , 10,  0., 100.);
  fHistVertex            = new TH1F("hVertex",               "control histogram for vertex Z position"    , 200, -20., 20.);

  if (fFill_TPC) {
    fHistIncTracks_dEdx    = new TH3F("fHistIncTracks_dEdx",   "dEdx histogram for inclusive tracks"        , n_tpc_mom_bins,tpc_mom_bins,n_dedx_bins,dedx_bins, n_eta_bins,eta_bins);
    fHistJetTracks_dEdx    = new TH3F("fHistJetTracks_dEdx",   "dEdx histogram for Jet tracks"            , n_tpc_mom_bins,tpc_mom_bins,n_dedx_bins,dedx_bins, n_eta_bins,eta_bins);
  }

  if (fFillpTPC_pT){
    fHistIncTracks_moms    = new TH2F("fHistIncTracks_moms",   "All mom types for inclusive tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
    fHistJetTracks_moms    = new TH2F("fHistJetTracks_moms",   "All mom types for Jet tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
  }
  if (fFillp_pT) {
    fHistIncTracks_moms_p    = new TH2F("fHistIncTracks_moms_p",   "All mom types for inclusive tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
    fHistJetTracks_moms_p    = new TH2F("fHistJetTracks_moms_p",   "All mom types for Jet tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
  }
  if (fFillpTPC_p) {
    fHistIncTracks_moms_pTPC_p    = new TH2F("fHistIncTracks_moms_pTPC_p",   "All mom types for inclusive tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
    fHistJetTracks_moms_pTPC_p    = new TH2F("fHistJetTracks_moms_pTPC_p",   "All mom types for Jet tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
  }

  fHistIncTracks_kin    = new TH3F("fHistIncTracks_kin",     "Kinematics histogram for inclusive tracks"  , 200, 0., 20., 9, 0.0, 0.9, 64, -3.2, 3.2);

  if (fFill_TOF){
    fHistIncTracks_beta    = new TH2F("fHistIncTracks_beta",   "Beta histogram for inclusive tracks"        , n_TOF_mom_bins,TOF_mom_bins,   n_beta_bins,beta_bins);
    fHistIncTracks_TOFpi_nsigma    = new TH3F("fHistIncTracks_TOFpi_nsigma",   "TOF Nsigma histogram for inclusive tracks under the pion hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistIncTracks_TOFka_nsigma    = new TH3F("fHistIncTracks_TOFka_nsigma",   "TOF Nsigma histogram for inclusive tracks under the kaon hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistIncTracks_TOFpr_nsigma    = new TH3F("fHistIncTracks_TOFpr_nsigma",   "TOF Nsigma histogram for inclusive tracks under the proton hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    fHistIncTracks_TOFpi_nsigma_1cls    = new TH3F("fHistIncTracks_TOFpi_nsigma_1cls",   "TOF Nsigma histogram for inclusive tracks under the pion hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistIncTracks_TOFka_nsigma_1cls    = new TH3F("fHistIncTracks_TOFka_nsigma_1cls",   "TOF Nsigma histogram for inclusive tracks under the kaon hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistIncTracks_TOFpr_nsigma_1cls    = new TH3F("fHistIncTracks_TOFpr_nsigma_1cls",   "TOF Nsigma histogram for inclusive tracks under the proton hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    fHistJetTracks_beta    = new TH2F("fHistJetTracks_beta",   "Beta histogram for jet tracks"        , n_TOF_mom_bins,TOF_mom_bins,   n_beta_bins,beta_bins);
    fHistJetTracks_TOFpi_nsigma    = new TH3F("fHistJetTracks_TOFpi_nsigma",   "TOF Nsigma histogram for jet tracks under the pion hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistJetTracks_TOFka_nsigma    = new TH3F("fHistJetTracks_TOFka_nsigma",   "TOF Nsigma histogram for jet tracks under the kaon hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistJetTracks_TOFpr_nsigma    = new TH3F("fHistJetTracks_TOFpr_nsigma",   "TOF Nsigma histogram for jet tracks under the proton hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    fHistJetTracks_TOFpi_nsigma_1cls    = new TH3F("fHistJetTracks_TOFpi_nsigma_1cls",   "TOF Nsigma histogram for jet tracks under the pion hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistJetTracks_TOFka_nsigma_1cls    = new TH3F("fHistJetTracks_TOFka_nsigma_1cls",   "TOF Nsigma histogram for jet tracks under the kaon hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHistJetTracks_TOFpr_nsigma_1cls    = new TH3F("fHistJetTracks_TOFpr_nsigma_1cls",   "TOF Nsigma histogram for jet tracks under the proton hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    Float_t t0_bins[141];
    for (int i=0; i < (140+1); i++){
      t0_bins[i] = -700.0 + i*10.0;
    }
    fHistIncTracks_t0   = new TH2F("fHistIncTracks_t0",   "TO histogram for inc tracks", n_TOF_mom_bins,TOF_mom_bins,   140,t0_bins);
  }

  fHistJetTracks_kin    = new TH3F("fHistJetTracks_kin",     "Kinematics histogram for Jet tracks"  , 200, 0., 20., 9, 0., 0.9, 64, -3.2, 3.2);

  if (fFill_TOF_expecs){
    fHistBetaExpec_pi    = new TH2F("fHistBetaExpec_pi",   "Expected Beta histogram for pions", n_TOF_mom_bins,TOF_mom_bins, n_beta_bins,beta_bins);
    fHistBetaExpec_ka    = new TH2F("fHistBetaExpec_ka",   "Expected Beta histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,  n_beta_bins,beta_bins);
    fHistBetaExpec_pr    = new TH2F("fHistBetaExpec_pr",   "Expected Beta histogram for protons", n_TOF_mom_bins,TOF_mom_bins,  n_beta_bins,beta_bins);

    fHistjet_BetaExpec_pi    = new TH2F("fHistjet_BetaExpec_pi",   "Expected Beta histogram for pions", n_TOF_mom_bins,TOF_mom_bins, n_beta_bins,beta_bins);
    fHistjet_BetaExpec_ka    = new TH2F("fHistjet_BetaExpec_ka",   "Expected Beta histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,  n_beta_bins,beta_bins);
    fHistjet_BetaExpec_pr    = new TH2F("fHistjet_BetaExpec_pr",   "Expected Beta histogram for protons", n_TOF_mom_bins,TOF_mom_bins, n_beta_bins,beta_bins);

    fHist_pi_mismatch    = new TH3F("fHist_pi_mismatch",   "Pion mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_ka_mismatch    = new TH3F("fHist_ka_mismatch",   "Kaon mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_pr_mismatch    = new TH3F("fHist_pr_mismatch",   "Proton mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    fHist_jet_pi_mismatch    = new TH3F("fHist_jet_pi_mismatch",   "Jet Pion mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_ka_mismatch    = new TH3F("fHist_jet_ka_mismatch",   "Jet Kaon mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_pr_mismatch    = new TH3F("fHist_jet_pr_mismatch",   "Jet Proton mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    fHist_elExpec_pihyp    = new TH3F("fHist_elExpec_pihyp",   "Expected Nsigma histogram for electrons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_muExpec_pihyp    = new TH3F("fHist_muExpec_pihyp",   "Expected Nsigma histogram for muons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_kaExpec_pihyp    = new TH3F("fHist_kaExpec_pihyp",   "Expected Nsigma histogram for kaons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_prExpec_pihyp    = new TH3F("fHist_prExpec_pihyp",   "Expected Nsigma histogram for protons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_piExpec_kahyp    = new TH3F("fHist_piExpec_kahyp",   "Expected Nsigma histogram for pions under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_prExpec_kahyp    = new TH3F("fHist_prExpec_kahyp",   "Expected Nsigma histogram for protons under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins,n_eta_bins,eta_bins);
    fHist_piExpec_prhyp    = new TH3F("fHist_piExpec_prhyp",   "Expected Nsigma histogram for pions under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_kaExpec_prhyp    = new TH3F("fHist_kaExpec_prhyp",   "Expected Nsigma histogram for kaons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_deExpec_prhyp    = new TH3F("fHist_deExpec_prhyp",   "Expected Nsigma histogram for deuterons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    fHist_jet_elExpec_pihyp    = new TH3F("fHist_jet_elExpec_pihyp",   "Expected Nsigma histogram for electrons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_muExpec_pihyp    = new TH3F("fHist_jet_muExpec_pihyp",   "Expected Nsigma histogram for muons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_kaExpec_pihyp    = new TH3F("fHist_jet_kaExpec_pihyp",   "Expected Nsigma histogram for kaons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_prExpec_pihyp    = new TH3F("fHist_jet_prExpec_pihyp",   "Expected Nsigma histogram for protons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_piExpec_kahyp    = new TH3F("fHist_jet_piExpec_kahyp",   "Expected Nsigma histogram for pions under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_prExpec_kahyp    = new TH3F("fHist_jet_prExpec_kahyp",   "Expected Nsigma histogram for protons under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_piExpec_prhyp    = new TH3F("fHist_jet_piExpec_prhyp",   "Expected Nsigma histogram for pions under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_kaExpec_prhyp    = new TH3F("fHist_jet_kaExpec_prhyp",   "Expected Nsigma histogram for kaons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
    fHist_jet_deExpec_prhyp    = new TH3F("fHist_jet_deExpec_prhyp",   "Expected Nsigma histogram for deuterons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

    Float_t pi_expec_bins[401];
    for (int i=0; i < (400+1); i++){
      pi_expec_bins[i] = 60.0 + i*0.1;
    }
    Float_t ka_expec_bins[5401];
    for (int i=0; i < (5400+1); i++){
      ka_expec_bins[i] = 60.0 + i*0.1;
    }
    Float_t pr_expec_bins[5401];
    for (int i=0; i < (5400+1); i++){
      pr_expec_bins[i] = 60.0 + i*0.1;
    }

    fHistTOFSigmaExpec_pi    = new TH3F("fHistTOFSigmaExpec_pi",   "Expected TOF Sigma histogram for pions", n_TOF_mom_bins,TOF_mom_bins,   400,pi_expec_bins, n_eta_bins,eta_bins);
    fHistTOFSigmaExpec_ka    = new TH3F("fHistTOFSigmaExpec_ka",   "Expected TOF Sigma histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,   5400,ka_expec_bins, n_eta_bins,eta_bins);
    fHistTOFSigmaExpec_pr    = new TH3F("fHistTOFSigmaExpec_pr",   "Expected TOF Sigma histogram for protons", n_TOF_mom_bins,TOF_mom_bins,   5400,pr_expec_bins, n_eta_bins,eta_bins);

    fHistjet_TOFSigmaExpec_pi    = new TH3F("fHistjet_TOFSigmaExpec_pi",   "Expected TOF Sigma histogram for pions", n_TOF_mom_bins,TOF_mom_bins,  400,pi_expec_bins, n_eta_bins,eta_bins);
    fHistjet_TOFSigmaExpec_ka    = new TH3F("fHistjet_TOFSigmaExpec_ka",   "Expected TOF Sigma histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,   5400,ka_expec_bins, n_eta_bins,eta_bins);
    fHistjet_TOFSigmaExpec_pr    = new TH3F("fHistjet_TOFSigmaExpec_pr",   "Expected TOF Sigma histogram for protons", n_TOF_mom_bins,TOF_mom_bins,   5400,pr_expec_bins, n_eta_bins,eta_bins);
  }

  if (fFill_TPC_expecs){
    Float_t mean_expec_bins[361];
    for (int i=0; i < (360+1); i++){
      mean_expec_bins[i] = 40.0 + i*1.0;
    }
    Float_t sigma_expec_bins[26];
    for (int i=0; i < (25+1); i++){
      sigma_expec_bins[i] = 0.0 + i*1.0;
    }

    fHistIncTracks_mpi  = new TH3F("fHistIncTracks_mpi",     "Expected mean pion histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
    fHistIncTracks_spi  = new TH3F("fHistIncTracks_spi",     "Expected sigma pion histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);

    fHistIncTracks_mel  = new TH3F("fHistIncTracks_mel",     "Expected mean electron histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
    fHistIncTracks_sel  = new TH3F("fHistIncTracks_sel",     "Expected sigma electron histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);

    fHistIncTracks_mka  = new TH3F("fHistIncTracks_mka",     "Expected mean kaon histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
    fHistIncTracks_ska  = new TH3F("fHistIncTracks_ska",     "Expected sigma kaon histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);

    fHistIncTracks_mpr  = new TH3F("fHistIncTracks_mpr",     "Expected mean proton histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
    fHistIncTracks_spr  = new TH3F("fHistIncTracks_spr",     "Expected sigma proton histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);
  }

  fHistJet_ptsub_v_area  = new TH2F("fHistJet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);
  fHistJet_kin  = new TH3F("fHistJet_kin", "Kinematics histogram for Jets"  , 300, 0., 300., 48, -0.6, 0.6, 130, 0., 6.5);
  fHistJet_moms  = new TH2F("fHistJet_moms", "All mom types for jets"  , 300, 0., 300., 600, 0., 600.);

  if (fFill_TPC) {
    fListHist->Add(fHistIncTracks_dEdx);
    fListHist->Add(fHistJetTracks_dEdx);
  }
  if (fFillpTPC_pT){
    fListHist->Add(fHistIncTracks_moms);
    fListHist->Add(fHistJetTracks_moms);
  }

  if (fFillp_pT){
    fListHist->Add(fHistIncTracks_moms_p);
    fListHist->Add(fHistJetTracks_moms_p);
  }

  if (fFillpTPC_p){
    fListHist->Add(fHistIncTracks_moms_pTPC_p);
    fListHist->Add(fHistJetTracks_moms_pTPC_p);
  }

    fListHist->Add(fHistIncTracks_kin);

  if (fFill_TOF){
    fListHist->Add(fHistIncTracks_beta);
    fListHist->Add(fHistIncTracks_TOFpi_nsigma);
    fListHist->Add(fHistIncTracks_TOFka_nsigma);
    fListHist->Add(fHistIncTracks_TOFpr_nsigma);

    fListHist->Add(fHistIncTracks_TOFpi_nsigma_1cls);
    fListHist->Add(fHistIncTracks_TOFka_nsigma_1cls);
    fListHist->Add(fHistIncTracks_TOFpr_nsigma_1cls);

    fListHist->Add(fHistJetTracks_beta);
    fListHist->Add(fHistJetTracks_TOFpi_nsigma);
    fListHist->Add(fHistJetTracks_TOFka_nsigma);
    fListHist->Add(fHistJetTracks_TOFpr_nsigma);

    fListHist->Add(fHistJetTracks_TOFpi_nsigma_1cls);
    fListHist->Add(fHistJetTracks_TOFka_nsigma_1cls);
    fListHist->Add(fHistJetTracks_TOFpr_nsigma_1cls);

    fListHist->Add(fHistIncTracks_t0);
  }

  fListHist->Add(fHistJetTracks_kin);

  if (fFill_TOF_expecs){
    fListHist->Add(fHistBetaExpec_pi);
    fListHist->Add(fHistBetaExpec_ka);
    fListHist->Add(fHistBetaExpec_pr);
    fListHist->Add(fHistjet_BetaExpec_pi);
    fListHist->Add(fHistjet_BetaExpec_ka);
    fListHist->Add(fHistjet_BetaExpec_pr);
    fListHist->Add(fHist_pi_mismatch);
    fListHist->Add(fHist_ka_mismatch);
    fListHist->Add(fHist_pr_mismatch);
    fListHist->Add(fHist_jet_pi_mismatch);
    fListHist->Add(fHist_jet_ka_mismatch);
    fListHist->Add(fHist_jet_pr_mismatch);
    fListHist->Add(fHist_elExpec_pihyp);
    fListHist->Add(fHist_muExpec_pihyp);
    fListHist->Add(fHist_kaExpec_pihyp);
    fListHist->Add(fHist_prExpec_pihyp);
    fListHist->Add(fHist_piExpec_kahyp);
    fListHist->Add(fHist_prExpec_kahyp);
    fListHist->Add(fHist_piExpec_prhyp);
    fListHist->Add(fHist_kaExpec_prhyp);
    fListHist->Add(fHist_deExpec_prhyp);

    fListHist->Add(fHist_jet_elExpec_pihyp);
    fListHist->Add(fHist_jet_muExpec_pihyp);
    fListHist->Add(fHist_jet_kaExpec_pihyp);
    fListHist->Add(fHist_jet_prExpec_pihyp);
    fListHist->Add(fHist_jet_piExpec_kahyp);
    fListHist->Add(fHist_jet_prExpec_kahyp);
    fListHist->Add(fHist_jet_piExpec_prhyp);
    fListHist->Add(fHist_jet_kaExpec_prhyp);
    fListHist->Add(fHist_jet_deExpec_prhyp);

    fListHist->Add(fHistTOFSigmaExpec_pi);
    fListHist->Add(fHistTOFSigmaExpec_ka);
    fListHist->Add(fHistTOFSigmaExpec_pr);
    fListHist->Add(fHistjet_TOFSigmaExpec_pi);
    fListHist->Add(fHistjet_TOFSigmaExpec_ka);
    fListHist->Add(fHistjet_TOFSigmaExpec_pr);
  }

  fListHist->Add(fHistJet_ptsub_v_area);
  fListHist->Add(fHistJet_kin);
  fListHist->Add(fHistJet_moms);

  if (fFill_TPC_expecs){
    fListHist->Add(fHistIncTracks_mpi);
    fListHist->Add(fHistIncTracks_spi);

    fListHist->Add(fHistIncTracks_mel);
    fListHist->Add(fHistIncTracks_sel);

    fListHist->Add(fHistIncTracks_mka);
    fListHist->Add(fHistIncTracks_ska);

    fListHist->Add(fHistIncTracks_mpr);
    fListHist->Add(fHistIncTracks_spr);
  }

  //
  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fTreejetsEMCconst  = ((*fTreeSRedirector)<<"jetsEMCconst").GetTree();
  fTreejetsEMCBGconst  = ((*fTreeSRedirector)<<"jetsEMCBGconst").GetTree();
  fTreejetsFJ        = ((*fTreeSRedirector)<<"jetsFJ").GetTree();
  fTreejetsFJBG      = ((*fTreeSRedirector)<<"jetsFJBG").GetTree();
  fTreejetsEMC       = ((*fTreeSRedirector)<<"jetsEMC").GetTree();
  fTreejetsEMCBG       = ((*fTreeSRedirector)<<"jetsEMCBG").GetTree();
  fTreejetsFJconst   = ((*fTreeSRedirector)<<"jetsFJconst").GetTree();
  fTreejetsFJBGconst   = ((*fTreeSRedirector)<<"jetsFJBGconst").GetTree();
  fTreejetEvents     = ((*fTreeSRedirector)<<"jeteventInfo").GetTree();
  fTreeMC        = ((*fTreeSRedirector)<<"fTreeMC").GetTree();
  fTreeCuts      = ((*fTreeSRedirector)<<"tracks").GetTree();
  //
  // ************************************************************************
  //   Send output objects to container
  // ************************************************************************
  //
  PostData(1, fListHist);
  PostData(2, fTreejetsEMCconst);
  PostData(3, fTreejetsEMCBGconst);
  PostData(4, fTreejetsFJ);
  PostData(5, fTreejetsFJBG);
  PostData(6, fTreejetsFJconst);
  PostData(7, fTreejetsFJBGconst);
  PostData(8, fTreejetEvents);
  PostData(9, fTreejetsEMC);
  PostData(10, fTreejetsEMCBG);
  PostData(11, fTreeMC);
  PostData(12, fTreeCuts);

  std::cout << " Info::siweyhmi: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::Run()
{
  //
  // main event loop
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the UserExec ===== " << std::endl;
  //
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) fMCtrue = kFALSE;
  fEventCountInFile++;
  //
  //  get the filename
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf(" Error::siweyhmi: Could not receive input chain"); return kFALSE;}
  TObjString fileName(chain->GetCurrentFile()->GetName());
  //
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  //
  if (eventHandler) fMCEvent = eventHandler->MCEvent();
  AliGenEventHeader* genHeader = 0x0;
  TString genheaderName;
  if (fMCEvent){
    genHeader = fMCEvent->GenEventHeader();
    genheaderName = genHeader->GetName();
    if(!genHeader){ Printf(" Error::siweyhmi: Event generator header not available!!!\n"); return kFALSE; }
  }
  //
  // If ESDs exist get some event variables
  //
  fCentrality = -5;
  fCentImpBin =-10.;
  AliCentrality    *esdCentrality = 0x0;
  AliMultSelection *MultSelection = 0x0;
  AliMultSelectionTask *MultSelectionTask = 0x0;
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (fESD)
  {
    //
    // Init magnetic filed for golden chi2 cut
    fESD->InitMagneticField();
    //
    // event selection
      if ( (fPassIndex==3 || fPassIndex==2) && fYear>2013){
        //
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1); // standard
        if (!fEventCuts.AcceptEvent(fESD)) {return kFALSE;}
      }
    //
    //
    esdCentrality = fESD->GetCentrality();
    MultSelection = (AliMultSelection*) fESD-> FindListObject("MultSelection");
  }

  //
  // Get rid of "E-AliESDpid::GetTPCsignalTunedOnData: Tune On Data requested, but MC event not set. Call SetCurrentMCEvent before!" errors
  //if (!fPIDResponse) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (fMCEvent && fPIDResponse) fPIDResponse->SetCurrentMCEvent(fMCEvent);
  //
  if(fMCtrue)
  {
    //
    // Add different generator particles to PDG Data Base to avoid problems when reading MC generator particles
    AliPDG::AddParticlesToPdgDataBase();
    //
    // ========================== MC =========================
    //
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) { Printf(" Error::siweyhmi: No MC stack available !!!\n"); return kFALSE;}
    //
    if (MultSelection) {
      fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::siweyhmi: Centralitity is taken from MultSelection " << fCentrality << std::endl;
    } else if (esdCentrality) {
      fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::siweyhmi: Centralitity is taken from esdCentrality " << fCentrality << std::endl;
    }
    //
    // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
    // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
    Float_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
    AliGenHijingEventHeader *lHIJINGHeader = 0x0;  // event header for HIJING

    if (!TMath::IsNaN(((AliGenHijingEventHeader*) genHeader)->ImpactParameter()) && fMCEvent){
      //
      // If HIJING based MC
      lHIJINGHeader = (AliGenHijingEventHeader*) genHeader;
      fNHardScatters = lHIJINGHeader->HardScatters();
      fNProjectileParticipants = lHIJINGHeader->ProjectileParticipants();
      fNTargetParticipants     = lHIJINGHeader->TargetParticipants();
      fNNColl   = lHIJINGHeader->NN();
      fNNwColl  = lHIJINGHeader->NNw();
      fNwNColl  = lHIJINGHeader->NwN();
      fNwNwColl = lHIJINGHeader->NwNw();
      fMCImpactParameter = lHIJINGHeader->ImpactParameter();
      fHistImpParam->Fill(fMCImpactParameter);
      if (fUseCouts)  std::cout << " Info::siweyhmi: HIJING: Centralitity is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenHijingEventHeader*) genHeader)->GetName() << std::endl;
      if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
      if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
      if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
      if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
      if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
      if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
      if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
      if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
      if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
      if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
    }
    if (fCentrality<0) fCentrality=fCentImpBin;
    //
    if (fUseCouts) {
      std::cout << " Info::siweyhmi: ========================================================================================== " << std::endl;
      std::cout << " Info::siweyhmi: " << fEventCountInFile << std::endl;
      std::cout << " Info::siweyhmi: Centrality = " << fCentrality << " ----- Impact Param = " << fMCImpactParameter << " fCentralityImp = " << fCentImpBin << std::endl;
      std::cout << " Info::siweyhmi: ========================================================================================== " << std::endl;
    }
  }
    //
    // ========================== Real =========================
    //
    if (!fESD)          { Printf(" Error::siweyhmi: fESD not available"); return kFALSE; }
    if (!fESDtrackCuts) { Printf(" Error::siweyhmi: fESDtrackCuts not available"); return kFALSE; }
    //
    // ------------------------------------------------
    // ------- monitor vertex position =---------------
    // ------------------------------------------------
    //
    Bool_t isVertexOk = kTRUE;
    fVertex = fESD->GetPrimaryVertexTracks();
    if( fVertex->GetNContributors()<1) isVertexOk = kFALSE;
    if( fVertex->GetNContributors()>1) {
      fVz    = fVertex->GetZ();
    }
    fNContributors   = fVertex->GetNContributors();
    //
    // ------------------------------------------------
    // ------- event vertex cut along Z ---------------
    // ------------------------------------------------
    //
    // if (fMCtrue && TMath::Abs(fVz) > 10) return;   // For MC put fixed cut
    if (TMath::Abs(fVz)>10) return kFALSE;
    //
    if (fVertex && isVertexOk) fHistVertex->Fill(fVz);
    else return kFALSE;
    //
    // ------------------------------------------------
    // ---------- Centrality definition ---------------
    // ------------------------------------------------
    //
      if (MultSelection) {
        fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      } else if (esdCentrality) {
        fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      } else {
        std::cout << " Info::siweyhmi: Error: There is no cent info " << std::endl;
      }
    //
    if (fUseCouts) {
      std::cout << " Info::siweyhmi: =============================================================================================== " << std::endl;
      std::cout << " Info::siweyhmi: Event counter = " << fEventCountInFile << " - cent =  " << fCentrality << std::endl;
      std::cout << " Info::siweyhmi: =============================================================================================== " << std::endl;
    }
  fHistCentrality->Fill(fCentrality);  // count events after physics and vertex selection
  //
  // in case small stat is enough
  if (fPercentageOfEvents>0 && (fEventCountInFile%fPercentageOfEvents)!=0) return kFALSE;

  // if centrality estimation failed
  if (fCentrality > 100.) return kFALSE;
  //
  // ======================================================================
  //   Filling part
  // ======================================================================
  //
  // Real Data Analysis
  //
  if (!fMCtrue && fESD && fCentrality>=fcent_min && fCentrality<fcent_max){
    fisGoodIncEvent = 0;
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = 0.0;
    fjetRhoVal = 0.0;
    fNumRealFJJets = 0;
    fNumRealEMCJets = 0;

    if (fDoIncTracks){
      FillIncTracksReal();
    }
    if (fDoEMCJet){
      FindJetsEMC();
    }
    if (fDoFastJet){
      FindJetsFJ();
    }
    FillEventTree();
    if (fUseCouts)  std::cout << " Info::siweyhmi: (Real Data Analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return kTRUE;
  }
  //
  // full MC analysis
  //
  if (fMCtrue && fESD  && fCentrality>=fcent_min && fCentrality<fcent_max){
    fisGoodIncEvent = 0;
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = 0.0;
    fjetRhoVal = 0.0;
    fNumRealFJJets = 0;
    fNumRealEMCJets = 0;

    FillTreeMC();
    if (fDoEMCJet){
      FindJetsEMC();
    }
    if (fDoFastJet){
      FindJetsFJ();
    }
    FillEventTree();
    if (fUseCouts)  std::cout << " Info::siweyhmi: (full MC analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return kTRUE;
  }
  return kTRUE;

}
//________________________________________________________________________
void AliAnalysisJetHadro::FindJetsEMC()
{

  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsEMC ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();
  //
  Float_t pT_sub_min = fjetMinPtSub;
  if (fUseCouts) cout << "Minimum jet pt after subtraction is " << fjetMinPtSub << endl;
  // Get the jet container
  fJetContainer = this->GetJetContainer("detJets");
  TString fRhoName = fJetContainer->GetRhoName();
  if (fUseCouts) cout << "Rho Name is " << fRhoName << endl;

  if (fJetContainer->GetRhoParameter()) fjetRhoVal = fJetContainer->GetRhoVal();
  if (fUseCouts) cout << "In FindJetsEMC Rho value is " << fjetRhoVal << endl;

  //
  // Get some jet container properties
  Int_t Njets         = fJetContainer->GetNJets();
  Int_t NAcceptedjets = fJetContainer->GetNAcceptedJets();
  Float_t jetRadius  = fJetContainer->GetJetRadius();
  Float_t jetEtaMin = fJetContainer->GetJetEtaMin();
  Float_t jetEtaMax = fJetContainer->GetJetEtaMax();
  UInt_t jetAcceptanceType = fJetContainer->GetAcceptanceType();
  //
  // loop over jets
  for(auto jet : fJetContainer->accepted())
  {
    //if (TMath::Abs(jet->Eta()) >= 0.9-jetRadius) continue; //this is not needed when using jetcontainer->accepted
    fhasAcceptedEMCjet = 1;
    fJetPt = jet->Pt();
    fJetEta = jet->Eta();
    fJetPhi = jet->Phi();
    Float_t JetM = jet->M();
    Int_t JetLabel = jet->GetLabel();
    Float_t JetAreaPt = jet->AreaPt(); //jet transverse area
    Int_t JetNumberOfConstituents = jet->GetNumberOfParticleConstituents(); //tracks + clusters = particles

    Float_t particleMaxChargedPt = jet->MaxChargedPt(); //max pT of charged particle in jet
    Float_t jetptsub = jet->PtSub(fjetRhoVal, kFALSE);

    if (jetptsub > pT_sub_min && JetAreaPt > fjetMinArea)
    {
      fhasRealEMCjet = 1;
      fNumRealEMCJets +=1;
    }

    if(!fTreeSRedirector) return;
    if (fFillEMCJet){
      (*fTreeSRedirector)<<"jetsEMC"<<
      "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
      "jetRhoVal=" << fjetRhoVal <<
      "jetpt="     << fJetPt    << // jetPt
      "jeteta="    << fJetEta   << // jetEta
      "jetphi="    << fJetPhi   << // jetPhi
      "jetArea="   << JetAreaPt  <<
      "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
      "particleMaxChargedPt="  << particleMaxChargedPt <<
      "jetptsub="  << jetptsub <<
      "cent="      << fCentrality  <<  //  centrality
      "\n";
    }

    fHistJet_ptsub_v_area->Fill(JetAreaPt, jetptsub);

    if (JetAreaPt <= fjetMinArea) continue;
    fHistJet_kin->Fill(jetptsub, fJetEta, fJetPhi);
    fHistJet_moms->Fill(jetptsub, fJetPt);

    if (jetptsub <= pT_sub_min) continue;

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->Track(i);
      AliESDtrack* esdtrack = (AliESDtrack*)(particle);
      if (!esdtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      //
      if (!esdtrack->GetInnerParam()) continue;               // Ask if track is in the TPC
      //if (!fESDtrackCuts->AcceptTrack(esdtrack))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      if (!(esdtrack->GetTPCsignalN()>0)) continue;
      //
      // Get the track variables
      GetExpecteds(esdtrack);
      SetCutBitsAndSomeTrackVariables(esdtrack);
      Float_t length    = esdtrack->GetIntegratedLength();
      Float_t tofSignal = esdtrack->GetTOFsignal();
      Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
      Float_t beta = -.05;
      int nTOFClusters = esdtrack->GetNTOFclusters(); //All matchable clusters
      if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

      Float_t fTPCmom_choice = fPtot;
      if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
      if (fSetTPCmom == 2) fTPCmom_choice = fPt;

      Float_t fTOFmom_choice = fPt;
      if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
      if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

      Float_t fBetamom_choice = fPVertex;
      if (fSetBetamom == 1) fTPCmom_choice = fPt;
      if (fSetBetamom == 2) fTPCmom_choice = fPtot;

      Float_t fEta_choice = fEta;
      if (fSetEta == 1) fEta_choice = fY;

      if (TMath::Abs(fEta) < 0.9){
        if (fFill_TPC) fHistJetTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
        if (fFillpTPC_pT) fHistJetTracks_moms->Fill(fPt,fPtot);
        if (fFillp_pT) fHistJetTracks_moms_p->Fill(fPt,fPVertex);
        if (fFillpTPC_p) fHistJetTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
        fHistJetTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
      }

      //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
      Bool_t fTOFout = kFALSE;
      if ((esdtrack->GetStatus() & AliESDtrack::kTOFout) != 0) {
      fTOFout = kTRUE;
      }

      //
      //kTIME flag
      Bool_t fTime = kFALSE;
      if ((esdtrack->GetStatus() & AliESDtrack::kTIME) != 0) {
      fTime = kTRUE;
      }

      Double_t inttime[6]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, 5 = deuteron
      esdtrack->GetIntegratedTimes(inttime, 6); // Returns the array with integrated times for each particle hypothesis

      Float_t fTOFMismatchTime = -20000000.;
      fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(esdtrack->Eta());
      if (fTOFMismatchTime <= 0){
        fTOFMismatchTime = -20000000.;
      }

      for (Int_t i = 0; i < 3; i++) {
        const Float_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
        if (i==0){
          fHistjet_BetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
          fHist_jet_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHistjet_TOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasElTOF) < 2.0){
            fHist_jet_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasMuTOF) < 2.0){
            fHist_jet_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_jet_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_jet_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
        }
        if (i==1){
          fHistjet_BetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
          fHist_jet_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHistjet_TOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_jet_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_jet_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
        }
        if (i==2){
          fHistjet_BetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
          fHist_jet_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHistjet_TOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_jet_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_jet_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasDeTOF) < 2.0){
            fHist_jet_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
        }
      }

      if (TMath::Abs(fEta) < 0.9  && fTOFout && fTime && fFill_TOF){
        fHistJetTracks_beta->Fill(fBetamom_choice,beta);
        fHistJetTracks_TOFpi_nsigma->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
        fHistJetTracks_TOFka_nsigma->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
        fHistJetTracks_TOFpr_nsigma->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
        if (nTOFClusters < 2) { //ESD only
          fHistJetTracks_TOFpi_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
          fHistJetTracks_TOFka_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
          fHistJetTracks_TOFpr_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
        }
      }

      if (fFillJetsEMCConst && fTOFout && fTime){

        (*fTreeSRedirector)<<"jetsEMCconst"<<
        //
        "jetRadius=" << jetRadius << // jet Radius
        "Njets="     << Njets << // Number of jets in event
        "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
        "jetRhoVal=" << fjetRhoVal <<
        "jetpt="     << fJetPt    << // jetPt
        "jeteta="    << fJetEta   << // jetEta
        "jetphi="    << fJetPhi   << // jetEta
        "jetArea="   << JetAreaPt  <<
        "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
        "particleMaxChargedPt="  << particleMaxChargedPt <<
        "jetptsub="  << jetptsub <<
        "dEdx="      << fTPCSignal            <<  //  dEdx of the track
        "sign="      << fSign                 <<  //  charge
        "ptot="      << fPtot                 <<  //  TPC momentum
        "p="         << fPVertex              <<  //  momentum at vertex
        "pT="        << fPt                   <<  // transverse momentum
        "eta="       << fEta                  <<  //  eta
        "phi="       << fPhi                  <<  //  phi
        "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
        "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
        "dEdxMeanPi="  << fDEdxPi              <<
        "dEdxSigmaPi=" << fSigmaPi            <<
        "dEdxMeanKa="  << fDEdxKa              <<
        "dEdxSigmaKa=" << fSigmaKa            <<
        "dEdxMeanPr="  << fDEdxPr              <<
        "dEdxSigmaPr=" << fSigmaPr            <<
        "cent="      << fCentrality;
        if (!fSmallOut){
          (*fTreeSRedirector)<<"jetsEMCconst"<<
          "primMult="  << fNContributors <<          //  #prim tracks
          "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
          "pitpcpid="  << fNSigmasPiTPC         <<
          "katpcpid="  << fNSigmasKaTPC         <<
          "prtpcpid="  << fNSigmasPrTPC         <<
          "tofSignal=" << tofSignal         <<
          "prtofpid="  << fNSigmasPrTOF<<
          "dEdxMeanDe="  << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<
          "beta="      << beta;
        }
        (*fTreeSRedirector)<<"jetsEMCconst"<<"\n";
      }
    }

  }

  if (fFillJetsEMCBG || fFillJetsEMCBGConst){
  fbgJetContainer = this->GetJetContainer("bgJets");
  //
  // Get some jet container properties
  Int_t Njetsbg         = fbgJetContainer->GetNJets();
  Int_t NAcceptedjetsbg = fbgJetContainer->GetNAcceptedJets();
  Float_t jetRadiusbg  = fbgJetContainer->GetJetRadius();
  Float_t jetEtaMinbg = fbgJetContainer->GetJetEtaMin();
  Float_t jetEtaMaxbg = fbgJetContainer->GetJetEtaMax();
  //
  // loop over jets
  for(auto jet : fbgJetContainer->accepted())
  {
    fJetPt = jet->Pt();
    fJetEta = jet->Eta();
    fJetPhi = jet->Phi();
    Float_t JetM = jet->M();
    Int_t JetLabel = jet->GetLabel();
    Float_t JetAreaPt = jet->AreaPt(); //jet transverse area
    Int_t JetNumberOfConstituents = jet->GetNumberOfParticleConstituents(); //tracks + clusters
    Bool_t IsJetMc = jet->IsMC(); //using >0.0 of jet area inside emcal
    Float_t particleMaxChargedPt = jet->MaxChargedPt(); //max pT of charged particle in jet
    Float_t jetMCPt = jet->MCPt();
    Float_t jetptsub = jet->PtSub(fjetRhoVal, kFALSE); //bg sub pt

    if (fFillJetsEMCBG){
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"jetsEMCBG"<<
    "jetRadius=" << jetRadiusbg << // jet Radius
    "Njets="     << Njetsbg << // Number of jets in event
    "NAcceptedjets=" << NAcceptedjetsbg << // Number of accepted jets in event
    "jetEtaMinbg=" << jetEtaMinbg <<
    "jetEtaMaxbg=" << jetEtaMaxbg <<
    "jetRhoVal=" << fjetRhoVal <<
    "jetpt="     << fJetPt    << // jetPt
    "jeteta="    << fJetEta   << // jetEta
    "jetphi="    << fJetPhi   << // jetPhi
    "jetArea="   << JetAreaPt  <<
    "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
    "particleMaxChargedPt="  << particleMaxChargedPt <<
    "jetptsub="  << jetptsub <<
    "cent="      << fCentrality  <<  //  centrality
    "\n";
    }

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->Track(i);
      AliESDtrack* esdtrack = (AliESDtrack*)(particle);
      if (!esdtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      //
      if (!esdtrack->GetInnerParam()) continue;               // Ask if track is in the TPC
      //if (!fESDtrackCuts->AcceptTrack(esdtrack))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      if (!(esdtrack->GetTPCsignalN()>0)) continue;
      //
      // Get the track variables
      Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(esdtrack);
      GetExpecteds(esdtrack);
      SetCutBitsAndSomeTrackVariables(esdtrack);
      Float_t length    = esdtrack->GetIntegratedLength();
      Float_t tofSignal = esdtrack->GetTOFsignal();
      Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
      Float_t beta = -.05;
      int nTOFClusters = esdtrack->GetNTOFclusters(); //All matchable clusters
      if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

      // Fill track constituent information
      if (fFillJetsEMCBGConst){

        (*fTreeSRedirector)<<"jetsEMCBGconst"<<
        //
        "jetRadius=" << jetRadiusbg << // jet Radius
        "Njets="     << Njetsbg << // Number of jets in event
        "NAcceptedjets=" << NAcceptedjetsbg << // Number of accepted jets in event
        "jetRhoVal=" << fjetRhoVal <<
        "jetpt="     << fJetPt    << // jetPt
        "jeteta="    << fJetEta   << // jetEta
        "jetphi="    << fJetPhi   << // jetPhi
        "jetArea="   << JetAreaPt  <<
        "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
        "particleMaxChargedPt="  << particleMaxChargedPt <<
        "jetptsub="  << jetptsub <<
        "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
        "dEdx="      << fTPCSignal            <<  //  dEdx of the track
        "sign="      << fSign                 <<  //  charge
        "ptot="      << fPtot                 <<  //  TPC momentum
        "p="         << fPVertex              <<  //  momentum at vertex
        "pT="        << fPt                   <<  // transverse momentum
        "eta="       << fEta                  <<  //  eta
        "phi="       << fPhi                  <<  //  phi
        "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
        "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
        "dEdxMeanPi="  << fDEdxPi              <<
        "dEdxSigmaPi=" << fSigmaPi            <<
        "dEdxMeanKa="  << fDEdxKa              <<
        "dEdxSigmaKa=" << fSigmaKa            <<
        "dEdxMeanPr="  << fDEdxPr              <<
        "dEdxSigmaPr=" << fSigmaPr            <<
        "cent="      << fCentrality;
        if (!fSmallOut){
          (*fTreeSRedirector)<<"jetsEMCBGconst"<<
          "primMult="  << fNContributors <<          //  #prim tracks
          "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
          "pitpcpid="  << fNSigmasPiTPC         <<
          "katpcpid="  << fNSigmasKaTPC         <<
          "prtpcpid="  << fNSigmasPrTPC         <<
          "tofSignal=" << tofSignal         <<
          "prtofpid="  << fNSigmasPrTOF<<
          "dEdxMeanDe="  << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<
          "beta="      << beta;
        }
        (*fTreeSRedirector)<<"jetsEMCBGconst"<<"\n";
      }
    }

  }
  }

}
//________________________________________________________________________
void AliAnalysisJetHadro::FindJetsFJ()
{
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsFJ ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();
  //
  // Create jetwrapper with the same settings used in FindJetsEMC
  int nJetRadiusBins = 1;
  float fTrackPt = 0.15;
  float fTrackEta = 0.9;
  float fGhostArea = 0.005;
  float bgJetAbsEtaCut = 0.7;           // fixed
  float bgJetRadius = 0.2;              // fixed
  //
  float pT_sub_min = fjetMinPtSub;    // can be 40, 60, 80
  std::vector<float>  fJetRadius;
  std::vector<float>  fPtSubMin;
  fJetRadius.resize(nJetRadiusBins);
  fJetRadius[0] = 0.4;
  fJetRadius[1] = 0.2;
  fJetRadius[2] = 0.6;
  fPtSubMin.resize(nJetRadiusBins);

    for (int iJetRadius=0; iJetRadius<nJetRadiusBins; iJetRadius++){
      Float_t jetAbsEtaCut = fTrackEta-fJetRadius[iJetRadius];   // fixed

      //
      //SOME CODE FROM NIMA
      AliFJWrapper *fFastJetWrapper;
      fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
      fFastJetWrapper->Clear();
      fFastJetWrapper->SetR(fJetRadius[iJetRadius]);
      fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
      fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
      //fFastJetWrapper->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
      fFastJetWrapper->SetGhostArea(fGhostArea);
      fFastJetWrapper->SetAreaType(fastjet::AreaType::active_area);
      fFastJetWrapper->SetMaxRap(1.0);

      AliFJWrapper *fFastJetWrapperBG;
      fFastJetWrapperBG = new AliFJWrapper("fFastJetWrapperBG","fFastJetWrapperBG");
      fFastJetWrapperBG->Clear();
      fFastJetWrapperBG->SetR(bgJetRadius);
      fFastJetWrapperBG->SetAlgorithm(fastjet::JetAlgorithm::kt_algorithm);
      fFastJetWrapperBG->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
      //fFastJetWrapperBG->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
      fFastJetWrapperBG->SetGhostArea(fGhostArea);
      fFastJetWrapperBG->SetAreaType(fastjet::AreaType::active_area_explicit_ghosts);
      fFastJetWrapperBG->SetMaxRap(1.0);

      //std::vector<fastjet::PseudoJet> particlesEmbeddedSubtracted; //will be filled with your subtracted event
      float particleEtaCut = 0.9;
      //
      // loop over esd tracks and add their four vector to wrapper --> identical to track container in EMC jet
      for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
        AliESDtrack* track = fESD->GetTrack(iTrack);
        if (!(fESDtrackCuts->AcceptTrack(track)) ) continue; // default cuts which should match EMC jets
        if (track->Pt() < fTrackPt || TMath::Abs(track->Eta()) >= particleEtaCut) continue;
        fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
        fFastJetWrapperBG->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
      }
      //
      // background jet definitions
      fastjet::JetMedianBackgroundEstimator bgE;
      fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
      bgE.set_selector(selectorBG);

      fastjet::Selector selectorBGjets = !fastjet::SelectorIsPureGhost() * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackPt, 1000.0);
      fFastJetWrapperBG->Run();
      std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBGjets(fFastJetWrapperBG->GetInclusiveJets()));

      if (jetsBG.size() > 0) {
        bgE.set_jets(jetsBG);  // give the kT jets to the background estimator
        frhoFJ = bgE.rho();
      }
      if (fUseCouts) std::cout << "frhoFJ is " << frhoFJ << std::endl;
      //
      // start of background jet loop
      if (fFilljetsFJBGTree || fFillJetsFJBGConst){
        for (Int_t ijet=0; ijet<Int_t(jetsBG.size()); ijet++) {
          fastjet::PseudoJet jet = jetsBG[ijet];
          //if (jet.pt() < fTrackPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= bgJetAbsEtaCut) continue; //redundant because of the selector
          Float_t jetpt = jet.pt();
          Float_t jetphi = jet.phi();
          Float_t jeteta = jet.eta();
          Float_t jetArea = jet.area();
          Float_t jetptsub = jetpt - frhoFJ*jetArea;
          Int_t jetNum = jetsBG.size();
          Int_t Njets = jetsBG.size();
          std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(jet.constituents());
          Int_t nConstituents = constituents.size();
          fastjet::PseudoJet highestpt_const = constituents[0];
          if (highestpt_const.pt() > 100.0) continue;

          if (fFilljetsFJBGTree)
          {
            (*fTreeSRedirector)<<"jetsFJBG"<<
            "bjJetRadius="    << bgJetRadius << // jet Radius
            "bgJetAbsEtaCut=" << bgJetAbsEtaCut << //abs eta cut for jet
            "jetNum="         << jetNum <<    //  number of jets
            "jetpt="          << jetpt <<
            "jetphi="         << jetphi <<
            "jeteta="         << jeteta <<
            "jetptsub="       << jetptsub << //bg sub jet pt (pt - rho*Area)
            "nConstituents="  << nConstituents << //num consts
            "rhoFJ="          << frhoFJ << //event rho
            "jetArea="        << jetArea << //jet area
            "cent="           << fCentrality  <<  //  centrality
            "\n";
          }

          for(Int_t i = 0; i < nConstituents; i++)
          {
            fastjet::PseudoJet &constituent = constituents[i];
            Float_t pt = constituent.pt();
            if (pt<1.e-10) continue;
            Float_t phi = constituent.phi();
            Float_t eta = constituent.eta();

            if (fFillJetsFJBGConst){
              (*fTreeSRedirector)<<"jetsFJBGconst"<<
              "pT="        << pt                   <<  // transverse momentum
              "eta="       << eta                  <<  //  eta
              "phi="       << phi                  <<  //  phi
              "\n";
            }
          }
        } // end of background jet loop
      }
      //
      // background subtraction on the constituent level TODO
      // fastjet::contrib::ConstituentSubtractor subtractorConstituent(&bgE); //add the background estimator to the correct subtractor
      // subtractorConstituent.set_common_bge_for_rho_and_rhom(true); //CHECK : should not be the case since particles have mass
      // subtractorConstituent.set_max_standardDeltaR(0.25); // set the max event wise subtraction distance
      // particlesEmbeddedSubtracted = subtractorConstituent.subtract_event(particlesEmbedded, particleEtaCut); //perform subtraction and fill the subtracted event container
      //

      // run jet finder using wrapper
      fFastJetWrapper->Run();
      std::vector<fastjet::PseudoJet> jets = sorted_by_pt(fFastJetWrapper->GetInclusiveJets());
      // start of jet loop
      for (Int_t ijet=0; ijet<Int_t(jets.size()); ijet++){
        fastjet::PseudoJet jet = jets[ijet];
        if (jet.pt() < fTrackPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
        Float_t jetpt = jet.pt();
        Float_t jetphi = jet.phi();
        Float_t jeteta = jet.eta();
        Float_t jetArea = jet.area();
        Float_t jetptsub = jetpt - frhoFJ*jetArea;
        Int_t jetNum = jets.size();
        std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(jets[ijet].constituents());

        Int_t nConstituents = constituents.size();
        if (nConstituents<1) continue;
        fastjet::PseudoJet highestpt_const = constituents[0];
        if (highestpt_const.pt() > 100.0) continue;
        fhasAcceptedFJjet = 1;

        if (jetptsub > pT_sub_min && jetArea > fjetMinArea)
        {
          fhasRealFJjet = 1;
          fNumRealFJJets +=1;
        }

        if (fFillFastJet){
          (*fTreeSRedirector)<<"jetsFJ"<<
          "jetNum="       << jetNum <<    //  number of jets
          "jetpt="        << jetpt <<
          "jetphi="       << jetphi <<
          "jeteta="       << jeteta <<
          "nConst="       << nConstituents <<    //  global event ID //this includes ghosts! Be careful
          "cent="         << fCentrality           <<  //  centrality
          "jetptsub="     << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="        << frhoFJ << //event rho
          "jetArea="      << jetArea << //jet area
          "\n";
        }

        fHistJet_ptsub_v_area->Fill(jetArea, jetptsub);

        if (jetArea <= fjetMinArea) continue;
        fHistJet_kin->Fill(jetptsub, jeteta, jetphi);
        fHistJet_moms->Fill(jetptsub, jetpt);

        if (jetptsub <= pT_sub_min) continue;

        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituent = constituents[i];
          Float_t pt = constituent.pt();
          if (pt<1.e-10) continue;

          Int_t trackIndex = constituent.user_index();
          AliESDtrack* trackConst = fESD->GetTrack(trackIndex);

          //Track cuts start
          fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
          fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
          //
          // --------------------------------------------------------------
          //      Get relevant track info and set cut bits
          // --------------------------------------------------------------
          //
          Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackConst);
          Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackConst);
          Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackConst) || fESDtrackCuts_Bit768_v->AcceptTrack(trackConst) ;
          if (!trackConst->GetInnerParam()) continue;               // Ask if track is in the TPC
          if (!(trackConst->GetTPCsignalN()>0)) continue;
          //
          // Get the track variables
          GetExpecteds(trackConst);
          SetCutBitsAndSomeTrackVariables(trackConst);

          Float_t length    = trackConst->GetIntegratedLength();
          Float_t tofSignal = trackConst->GetTOFsignal();
          Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
          Float_t beta = -.05;
          int nTOFClusters = trackConst->GetNTOFclusters(); //All matchable clusters
          if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

          Float_t fTPCmom_choice = fPtot;
          if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
          if (fSetTPCmom == 2) fTPCmom_choice = fPt;

          Float_t fTOFmom_choice = fPt;
          if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
          if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

          Float_t fBetamom_choice = fPVertex;
          if (fSetBetamom == 1) fTPCmom_choice = fPt;
          if (fSetBetamom == 2) fTPCmom_choice = fPtot;

          Float_t fEta_choice = fEta;
          if (fSetEta == 1) fEta_choice = fY;

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9){
            if (fFill_TPC) fHistJetTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
            if (fFillpTPC_pT) fHistJetTracks_moms->Fill(fPt,fPtot);
            if (fFillp_pT) fHistJetTracks_moms_p->Fill(fPt,fPVertex);
            if (fFillpTPC_p) fHistJetTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
            fHistJetTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
          }

          //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
          Bool_t fTOFout = kFALSE;
          if ((trackConst->GetStatus() & AliESDtrack::kTOFout) != 0) { //Track has the kTOFout flag
          fTOFout = kTRUE;
          }

          //
          //kTIME flag
          Bool_t fTime = kFALSE;
          if ((trackConst->GetStatus() & AliESDtrack::kTIME) != 0) { //Track has the kTIME flag
          fTime = kTRUE;
          }

          Double_t inttime[6]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, 5 = deuteron
          trackConst->GetIntegratedTimes(inttime, 6); // Returns the array with integrated times for each particle hypothesis

          Float_t fTOFMismatchTime = -20000000.;
          fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(trackConst->Eta());
          if (fTOFMismatchTime <= 0){
            fTOFMismatchTime = -20000000.;
          }

          for (Int_t i = 0; i < 3; i++) {
            const Float_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
            if (i==0){
              fHistjet_BetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
              fHist_jet_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHistjet_TOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasElTOF) < 2.0){
                fHist_jet_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasMuTOF) < 2.0){
                fHist_jet_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHist_jet_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHist_jet_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
            }
            if (i==1){
              fHistjet_BetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
              fHist_jet_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHistjet_TOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHist_jet_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHist_jet_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              }
            }
            if (i==2){
              fHistjet_BetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
              fHist_jet_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHistjet_TOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHist_jet_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHist_jet_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasDeTOF) < 2.0){
                fHist_jet_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
            }
          }

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9  && fTOFout && fTime && fFill_TOF){
            fHistJetTracks_beta->Fill(fBetamom_choice,beta);
            fHistJetTracks_TOFpi_nsigma->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
            fHistJetTracks_TOFka_nsigma->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
            fHistJetTracks_TOFpr_nsigma->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
            if (nTOFClusters < 2) { //ESD only
              fHistJetTracks_TOFpi_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
              fHistJetTracks_TOFka_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
              fHistJetTracks_TOFpr_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
            }
          }

          if (fFillJetsFJConst && fTOFout && fTime && ifDefaultCuts)
          {

          (*fTreeSRedirector)<<"jetsFJconst"<<
          "jetRadius=" << fJetRadius[iJetRadius] << // jet Radius
          "jetNum="    << jetNum <<    //  number of jets
          "jetpt="     << jetpt <<     //  global event ID
          "jetphi="    << jetphi <<    //  global event ID
          "jeteta="    << jeteta <<    //  global event ID
          "nConst="    << nConstituents <<    //  global event ID

          "jetptsub="  << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="     << frhoFJ << //event rho
          "jetArea="   << jetArea << //jet area

          "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand

          "dEdx="      << fTPCSignal            <<  //  dEdx of the track
          "sign="      << fSign                 <<  //  charge
          "ptot="      << fPtot                 <<  //  TPC momentum
          "p="         << fPVertex              <<  //  momentum at vertex
          "pT="        << fPt                   <<  // transverse momentum
          "eta="       << fEta                  <<  //  eta
          "cent="      << fCentrality           <<  //  centrality
          "phi="       << fPhi                  <<  //  phi

          "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
          "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
          "dEdxMeanPi=" << fDEdxPi              <<
          "dEdxSigmaPi=" << fSigmaPi            <<
          "dEdxMeanKa=" << fDEdxKa              <<
          "dEdxSigmaKa=" << fSigmaKa            <<
          "dEdxMeanPr=" << fDEdxPr              <<
          "dEdxSigmaPr=" << fSigmaPr            <<
          "dEdxMeanDe=" << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<

          "beta=" << beta         <<

          "eltofpid="  << fNSigmasElTOF         <<  // nsigma TPC for electrons
          "pitofpid="  << fNSigmasPiTOF         <<
          "katofpid="  << fNSigmasKaTOF         <<
          "prtofpid="  << fNSigmasPrTOF         <<
          "detofpid="  << fNSigmasDeTOF;

          if (!fSmallOut){
            (*fTreeSRedirector)<<"jetsFJconst"<<
            "bit128="    << fBit128 <<        // TPC only tracks cuts
            "bit768="    << fBit768 <<        // Hybrid track cuts
            "primMult="  << fNContributors <<          //  #prim tracks
            "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
            "pitpcpid="  << fNSigmasPiTPC         <<
            "katpcpid="  << fNSigmasKaTPC         <<
            "prtpcpid="  << fNSigmasPrTPC         <<
            "detpcpid="  << fNSigmasDeTPC         <<

            "tofSignal=" << tofSignal;
          }
        (*fTreeSRedirector)<<"jetsFJconst"<<"\n";
        }
        }
      } // end of jet loop

      delete fFastJetWrapper;
      delete fFastJetWrapperBG;
    }
}
//________________________________________________________________________
void AliAnalysisJetHadro::FillIncTracksReal()
{
  //
  // Fill dEdx information for the TPC and also clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillIncTracksReal ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();
  // --------------------------------------------------------------
  // Get the event
  AliVEvent *event=InputEvent();
  if (CountEmptyEvents()) return;
  fisGoodIncEvent = 1;
  //
  // --------------------------------------------------------------
  //  Main track loop
  // --------------------------------------------------------------
  //
  //
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop

    fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
    fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
    //
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track) continue;
    //
    // --------------------------------------------------------------
    //      Get relevant track info and set cut bits
    // --------------------------------------------------------------
    //
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(track);
    Bool_t cuts_2015 = fESDtrackCuts_2015->AcceptTrack(track);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(track);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(track) || fESDtrackCuts_Bit768_v->AcceptTrack(track) ;
    if (!track->GetInnerParam()) continue;               // Ask if track is in the TPC
    if (!(track->GetTPCsignalN()>0)) continue;
    //
    // Get the track variables
    GetExpecteds(track);
    SetCutBitsAndSomeTrackVariables(track);
    //
    // --------------------------------------------------------------
    //   Fill the trees
    // --------------------------------------------------------------
    //

    Float_t length = track->GetIntegratedLength();
    Float_t tofSignal = track->GetTOFsignal();
    Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
    Float_t beta = -.05;
    int nTOFClusters = track->GetNTOFclusters(); //All matchable clusters
    if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

    if (fPt>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

    Float_t fTPCmom_choice = fPtot;
    if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
    if (fSetTPCmom == 2) fTPCmom_choice = fPt;

    Float_t fTOFmom_choice = fPt;
    if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
    if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

    Float_t fBetamom_choice = fPVertex;
    if (fSetBetamom == 1) fTPCmom_choice = fPt;
    if (fSetBetamom == 2) fTPCmom_choice = fPtot;

    Float_t fEta_choice = fEta;
    if (fSetEta == 1) fEta_choice = fY;

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9){
      if (fFill_TPC) fHistIncTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
      if (fFillpTPC_pT) fHistIncTracks_moms->Fill(fPt,fPtot);
      if (fFillp_pT) fHistIncTracks_moms_p->Fill(fPt,fPVertex);
      if (fFillpTPC_p) fHistIncTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
      fHistIncTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);

      if (fFill_TPC_expecs){
        fHistIncTracks_mpi->Fill(fTPCmom_choice,fDEdxPi,TMath::Abs(fEta_choice));
        fHistIncTracks_spi->Fill(fTPCmom_choice,fSigmaPi,TMath::Abs(fEta_choice));

        fHistIncTracks_mel->Fill(fTPCmom_choice,fDEdxEl,TMath::Abs(fEta_choice));
        fHistIncTracks_sel->Fill(fTPCmom_choice,fSigmaEl,TMath::Abs(fEta_choice));

        fHistIncTracks_mka->Fill(fTPCmom_choice,fDEdxKa,TMath::Abs(fEta_choice));
        fHistIncTracks_ska->Fill(fTPCmom_choice,fSigmaKa,TMath::Abs(fEta_choice));

        fHistIncTracks_mpr->Fill(fTPCmom_choice,fDEdxPr,TMath::Abs(fEta_choice));
        fHistIncTracks_spr->Fill(fTPCmom_choice,fSigmaPr,TMath::Abs(fEta_choice));
      }
    }
    //

    //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
    Bool_t fTOFout = kFALSE;
    if ((track->GetStatus() & AliESDtrack::kTOFout) != 0) { //Track has the kTOFout flag
    fTOFout = kTRUE;
    }

    //kTIME flag
    Bool_t fTime = kFALSE;
    if ((track->GetStatus() & AliESDtrack::kTIME) != 0) { //Track has the kTIME flag
    fTime = kTRUE;
    }

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9 && fTOFout && fTime && fFill_TOF_expecs){
      Double_t inttime[6]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, 5 = deuteron
      track->GetIntegratedTimes(inttime, 6); // Returns the array with integrated times for each particle hypothesis

      Float_t fTOFMismatchTime = -20000000.;
      fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());
      if (fTOFMismatchTime <= 0){
        fTOFMismatchTime = -20000000.;
      }

      for (Int_t i = 0; i < 3; i++) {
        const Double_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
        if (i==0){
          fHistBetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
          fHist_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHistTOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasElTOF) < 2.0){
            fHist_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasMuTOF) < 2.0){
            fHist_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
        }
        if (i==1){
          fHistBetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
          fHist_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHistTOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
        }
        if (i==2){
          fHistBetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
          fHist_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHistTOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasDeTOF) < 2.0){
            fHist_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
        }
      }
    }

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9 && fTOFout && fTime && fFill_TOF){
      fHistIncTracks_beta->Fill(fBetamom_choice,beta);
      fHistIncTracks_t0->Fill(fBetamom_choice,t0);
      fHistIncTracks_TOFpi_nsigma->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
      fHistIncTracks_TOFka_nsigma->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
      fHistIncTracks_TOFpr_nsigma->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
      if (nTOFClusters < 2) { //ESD only
        fHistIncTracks_TOFpi_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
        fHistIncTracks_TOFka_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
        fHistIncTracks_TOFpr_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
      }
    }

	  if (fFillIncTracks && fTOFout && fTime && ifDefaultCuts)
    {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"tracks"<<
      //
      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "bit768="    << fBit768 <<        // Hybrid track cuts
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "phi="       << fPhi                  <<  //  phi
      "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
      "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
      "dEdxMeanPi=" << fDEdxPi              <<
      "dEdxSigmaPi=" << fSigmaPi            <<
      "dEdxMeanKa=" << fDEdxKa              <<
      "dEdxSigmaKa=" << fSigmaKa            <<
      "dEdxMeanPr=" << fDEdxPr              <<
      "dEdxSigmaPr=" << fSigmaPr            <<
      "cent="      << fCentrality           <<
      "eltofpid="  << fNSigmasElTOF         <<  // nsigma TPC for electrons
      "pitofpid="  << fNSigmasPiTOF         <<
      "katofpid="  << fNSigmasKaTOF         <<
      "prtofpid="  << fNSigmasPrTOF         <<
      "detofpid="  << fNSigmasDeTOF         <<
      "beta="      << beta;

      if (!fSmallOut){
        (*fTreeSRedirector)<<"tracks"<<
        "bit128="    << fBit128 <<        // TPC only tracks cuts
        "primMult="  << fNContributors <<          //  #prim tracks
        "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
        "pitpcpid="  << fNSigmasPiTPC         <<
        "katpcpid="  << fNSigmasKaTPC         <<
        "prtpcpid="  << fNSigmasPrTPC         <<
        "tofSignal=" << tofSignal         <<
        "prtofpid="  << fNSigmasPrTOF<<
        "dEdxMeanDe=" << fDEdxDe              <<
        "dEdxSigmaDe=" << fSigmaDe;

      }
      (*fTreeSRedirector)<<"tracks"<<"\n";
    }

  }// end of track loop

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillTreeMC()
{
  //
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();

  Int_t trackOrigin = -10;
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillTreeMC ===== " << std::endl;
  //
  // ======================================================================
  // ------   reconstructed MC particles with dEdx information-------------
  // ======================================================================
  //
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  {
    //
    // Esd track
    //
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    //
    // Get generated track info
    Int_t lab = TMath::Abs(trackReal->GetLabel());
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
    Int_t pdg = trackMCgen->Particle()->GetPdgCode();
    //
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(lab,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    //
    // check the origin of the track
    trackOrigin = -10;
    if (fMCStack->IsPhysicalPrimary(lab))        trackOrigin = 0;
    if (fMCStack->IsSecondaryFromMaterial(lab))  trackOrigin = 1;
    if (fMCStack->IsSecondaryFromWeakDecay(lab)) trackOrigin = 2;
    if (trackOrigin<-1) continue; // TODO
    //
    // Track cuts from dtector
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackReal);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackReal) || fESDtrackCuts_Bit768_v->AcceptTrack(trackReal) ;
    //
    if (!trackReal->GetInnerParam()) continue;     // TODO        // Ask if track is in the TPC
    if (!fESDtrackCuts->AcceptTrack(trackReal))  continue;    // TODO
    if (!(trackReal->GetTPCsignalN()>0)) continue; // TODO
    //
    // match the track with mc track
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
    if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
    if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
    if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
    if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
    if (iPart == -10) continue; // TODO
    //
    GetExpecteds(trackReal);
    SetCutBitsAndSomeTrackVariables(trackReal);
    //
    if (trackReal-> GetInnerParam()){
      fPtotMC       = trackReal->GetInnerParam()->GetP();
      fTPCSignalMC  = trackReal->GetTPCsignal();
    }
    fEtaMC        = trackReal->Eta();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();
    //
    //
    //
    Float_t dca[2], covar[3];
    trackReal->GetImpactParameters(dca, covar);
    Float_t length = trackReal->GetIntegratedLength();
    Float_t tofSignal = trackReal->GetTOFsignal();
    Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
    Float_t beta = -.05;
    int nTOFClusters = trackReal->GetNTOFclusters(); //All matchable clusters
    if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"fTreeMC"<<
    "orig="     << trackOrigin <<   // origin of the track
    "part="      << iPart <<
    "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
    "sign="      << fSignMC <<         // sign
    "ptot="      << fPtotMC <<         // tpc momentum
    "p="         << pMC <<             // vertex momentum
    "pT="        << fPtMC <<           // transverse momentum
    "eta="       << fEtaMC <<          // mc eta
    "phi="       << fPhiMC <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "centimp="   << fCentImpBin <<
    "vZ="        << fVz <<
    "nsigmatofka="  << fNSigmasKaTOF         <<  // interaction rate
    "nsigmatofpr="  << fNSigmasPrTOF         <<  // interaction rate
    //
    "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
    "tofSignal=" << tofSignal         <<
    "beta=" << beta         <<
    "defCut="    << ifDefaultCuts <<  // default cut
    "bit128="    << fBit128 <<  // run Number
    "bit768="    << fBit768 <<  // run Number
    "primmult="  << fNContributors <<  //  #prim tracks
    "fCdd="      << covar[0] <<
    "fCdz="      << covar[1] <<
    "fCzz="      << covar[2] <<
    "tpcpileup=" << isTPCPileup <<
    "itspileup=" << isITSPileup <<
    "\n";

  } // ======= end of track loop for MC dEdx filling =======

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillEventTree()
{
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEventTree ===== " << std::endl;
  (*fTreeSRedirector)<<"jeteventInfo"<<
  "rhoFJ="      << frhoFJ << //event rho
  "rhoEMC="  << fjetRhoVal <<
  "cent="      << fCentrality  <<  //  centrality
  "isGoodIncEvent="   << fisGoodIncEvent <<
  "hasAcceptedFJjet="   << fhasAcceptedFJjet <<
  "hasRealFJjet="   << fhasRealFJjet <<
  "hasAcceptedEMCjet="   << fhasAcceptedEMCjet <<
  "hasRealEMCjet="   << fhasRealEMCjet <<
  "NumRealFJJets="   << fNumRealFJJets <<
  "NumRealEMCJets="   << fNumRealEMCJets <<
  "\n";
}
//________________________________________________________________________
void AliAnalysisJetHadro::GetExpecteds(AliESDtrack *track)
{

  //
  // bettaGamma is not well deifned below bg=0.01 --> below 200MeV protons and deuterons
  Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
  Float_t ptotForBetaGammaThr = 0.2;
  // if (ptotForBetaGamma<ptotForBetaGammaThr) return;
  //
  // --------------------------------------------------------------
  //  Calculates expected sigma and dEdx for a given track and returns colesest expected particle and its index
  // --------------------------------------------------------------
  //
  fNSigmasElTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fNSigmasPiTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fNSigmasKaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //
  fNSigmasElTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasMuTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion,     fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,     fPIDResponse->GetTOFResponse().GetTimeZero());
  //
  if (ptotForBetaGamma>ptotForBetaGammaThr){
    fNSigmasPrTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fNSigmasDeTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    fNSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,   fPIDResponse->GetTOFResponse().GetTimeZero());
    fNSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron, fPIDResponse->GetTOFResponse().GetTimeZero());
  }

  //
  //
  //
  // Electron Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasElTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxEl  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kElectron);
    fSigmaEl = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kElectron);
  }
  //
  // Pion Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasPiTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxPi  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kPion);
    fSigmaPi = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kPion);
  }
  //
  // Kaon Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasKaTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxKa  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kKaon);
    fSigmaKa = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kKaon);
  }
  //
  // Proton Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasPrTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxPr  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kProton);
    fSigmaPr = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kProton);
  }
  //
  // Deuteron Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasDeTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxDe  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
    fSigmaDe = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
  }

}
//________________________________________________________________________
void AliAnalysisJetHadro::SetCutBitsAndSomeTrackVariables(AliESDtrack *track)
{
  //
  // Set some track variables
  //
  fPVertex = track->P();
  fSign= track->GetSign();
  fPt  =track->Pt();
  fY   =track->Y();
  fPhi =track->Phi()-TMath::Pi();
  fEta = track->Eta();

  //
  // TPC related quantities
  if (track->GetInnerParam()){
    fPtot      = track->GetInnerParam()->GetP();
    fTPCSignal = track->GetTPCsignal();
    //
    Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    //
    // --------------------------------------------------------------
    //      Bayesian PID part
    // --------------------------------------------------------------
    //
    if (ptotForBetaGamma>ptotForBetaGammaThr) {
      fPIDCombined->SetDefaultTPCPriors();
      Double_t probTPC[AliPID::kSPECIES]={0.};
      Double_t probTOF[AliPID::kSPECIES]={0.};
      // Get TPC probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      fTrackProbPiTPC = probTPC[AliPID::kPion];
      fTrackProbKaTPC = probTPC[AliPID::kKaon];
      fTrackProbPrTPC = probTPC[AliPID::kProton];
      //fTrackProbDeTPC = probTPC[AliPID::kDeuteron];
      // Get TOF probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
      fTrackProbPiTOF = probTOF[AliPID::kPion];
      fTrackProbKaTOF = probTOF[AliPID::kKaon];
      fTrackProbPrTOF = probTOF[AliPID::kProton];
    }
  }
  //

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::CountEmptyEvents()
{
  //
  // count Empty Events
  //
  Bool_t emptyEvent= kTRUE;
  for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack) {   // Track loop
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (track->GetInnerParam() && fESDtrackCuts->AcceptTrack(track) && track->GetTPCsignalN()>0 && track->Pt()<100.0) {
      emptyEvent= kFALSE;
      break;
    }
  }
  //
  // check if the event is empty
  if (emptyEvent) {
    std::cout << " Info::siweyhmi: Empty event " << std::endl;
  }
  else {
    if (fUseCouts) std::cout << " Info::siweyhmi: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  }
  return emptyEvent;

}
//________________________________________________________________________
void AliAnalysisJetHadro::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  std::cout << " Info::siweyhmi: ===== In the Terminate ===== " << std::endl;

  /*
  AliTPCPIDResponse fTPCresponse = fPIDResponse->GetTPCResponse();
  Double_t temp_mip = fTPCresponse.GetMIP();
  std::cout << "temp_mip is " << temp_mip << std::endl;
  for (int j=0; j<8; j++){
    AliPID::EParticleType species = AliPID::kPion;
    if (j==1) {
      species = AliPID::kElectron;
      std::cout << "el" << std::endl;
    }
    if (j==2) {
      species = AliPID::kKaon;
      std::cout << "ka" << std::endl;
    }
    if (j==3) {
      species = AliPID::kProton;
      std::cout << "pr" << std::endl;
    }
    if (j==4) {
      species = AliPID::kDeuteron;
      std::cout << "d" << std::endl;
    }
    if (j==5) {
      species = AliPID::kTriton;
      std::cout << "t" << std::endl;
    }
    if (j==6) {
      species = AliPID::kHe3;
      std::cout << "He3" << std::endl;
    }
    if (j==7) {
      species = AliPID::kAlpha;
      std::cout << "He4" << std::endl;
    }
    TSpline3* resp_fxn = fTPCresponse.GetResponseFunction(species, AliTPCPIDResponse::ETPCgainScenario::kDefault);
    Double_t mass=AliPID::ParticleMassZ(species);
    std::cout << "The mass is " << AliPID::ParticleMassZ(species) << std::endl;
    //charge factor. BB goes with z^2, however in reality it is slightly larger (calibration, threshold effects, ...)
    // !!! Splines for light nuclei need to be normalised to this factor !!!
    const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(species),2.3);
    std::cout << "The particle charge is " << AliPID::ParticleCharge(species) << std::endl;

    stringstream name;
    name << "resp_fxn_" << j;
    //gSystem->ChangeDirectory("/Users/sierraweyhmiller/Sierra/Caines_Research_Yr2/JetHadro/most_recent_code");
    //resp_fxn->Write();
    resp_fxn->SaveAs(name.str().c_str());
  }
  */

}
