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
// AliAnalysisTaskSE for HF quality assurance
//
// Author: Chiara Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliCounterCollection.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliVertexingHFUtils.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"

#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowVector.h"

#include "AliTRDTriggerAnalysis.h"

#include "AliAnalysisTaskSEHFQA.h"

using std::cout;
using std::endl;


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEHFQA);
/// \endcond

//____________________________________________________________________________

AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA():AliAnalysisTaskSE()
  , fOutputEntries(0x0)
  , fOutputPID(0x0)
  , fOutputTrack(0x0)
  , fOutputCounters(0x0)
  , fOutputCheckCentrality(0x0)
  , fOutputEvSelection(0x0)
  , fOutputFlowObs(0x0)
  , fDecayChannel(AliAnalysisTaskSEHFQA::kD0toKpi)
  , fCuts(0x0)
  , fFlowEvent(0x0)
  , fRFPcuts(0x0)
  , fEstimator(AliRDHFCuts::kCentTRK)
  , fReadMC(kFALSE)
  , fSimpleMode(kFALSE)
  , fUseSelectionBit(kTRUE)
  , fOnOff()
  , fFillDistrTrackEffChecks(kFALSE)
  , fAODProtection(1)
  , fHisNentries(0)
  , fHisNentriesSelBit(0)
  , fHisTOFflags(0)
  , fHisTOFsig(0)
  , fHisTOFstartTimeMask(0)
  , fHisTOFstartTimeRes(0)
  , fHisTOFstartTimeDistrib(0)
  , fHisTOFtime(0)
  , fHisTOFtimeKaonHyptime(0)
  , fHisTOFtimeKaonHyptimeAC(0)
  , fHisTOFsigmaKSigPid(0)
  , fHisTOFsigmaPionSigPid(0)
  , fHisTOFsigmaProtonSigPid(0)
  , fHisTOFsigPid3sigPion(0)
  , fHisTOFsigPid3sigKaon(0)
  , fHisTOFsigPid3sigProton(0)
  , fHisTPCsig(0)
  , fHisTPCsigvsp(0)
  , fHisTPCsigvspAC(0)
  , fHisTPCsigmaK(0)
  , fHisTPCsigmaPion(0)
  , fHisTPCsigmaProton(0)
  , fHisTPCsigNvsPtAllTracks(0)
  , fHisTPCsigNvsPhiAllTracks(0)
  , fHisTPCsigNvsEtaAllTracks(0)
  , fHisTPCsigNvsPtDaughters(0)
  , fHisTPCsigNvsPhiDaughters(0)
  , fHisTPCsigNvsEtaDaughters(0)
  , fHisTOFsigmaMCKSigPid(0)
  , fHisTOFsigmaMCPionSigPid(0)
  , fHisTOFsigmaMCProtonSigPid(0)
  , fHisTPCsigmaMCK(0)
  , fHisTPCsigmaMCPion(0)
  , fHisTPCsigmaMCProton(0)
  , fHisnClsITS(0)
  , fHisnClsITSselTr(0)
  , fHisnClsITSSA(0)
  , fHisnClsITSSAspdAny(0)
  , fHisnClsITSSAspdIn(0)
  , fHisnClsITSSAspdOut(0)
  , fHisnLayerITS(0)
  , fHisnLayerITSselTr(0)
  , fHisnLayerITSsa(0)
  , fHisnClsSPD(0)
  , fHisptGoodTr(0)
  , fHisptGoodTrFromDaugh(0)
  , fHisptGoodTrFromDaugh_filt(0)
  , fHisdistrGoodTr(0)
  , fHisdistrSelTr(0)
  , fHisd0dau(0)
  , fHisd0dau_filt(0)
  , fHisd0dauphi(0)
  , fHisd0dauphi_filt(0)
  , fHisd0zdau(0)
  , fHisd0zdau_filt(0)
  , fHisd0zdauphi(0)
  , fHisd0zdauphi_filt(0)
  , fHisd0TracksSPDin(0)
  , fHisd0TracksSPDany(0)
  , fHisd0TracksFilterBit4(0)
  , fHisd0TracksTPCITSSPDany(0)
  , fHisPtDaughters(0)
  , fHisPhiDaughters(0)
  , fHisEtaDaughters(0)
  , fHisEtavsPhiDaughters(0)
  , fHisNTPCclsvsPtDaughters(0)
  , fHisNTPCclsvsPhiDaughters(0)
  , fHisNTPCclsvsEtaDaughters(0)
  , fHisNTPCCrossedRowsvsPtDaughters(0)
  , fHisNTPCCrossedRowsvsPhiDaughters(0)
  , fHisNTPCCrossedRowsvsEtaDaughters(0)
  , fHisRatioCRowsOverFclsvsPtDaughters(0)
  , fHisRatioCRowsOverFclsvsPhiDaughters(0)
  , fHisRatioCRowsOverFclsvsEtaDaughters(0)
  , fHisNITSclsvsPtDaughters(0)
  , fHisNITSclsvsPhiDaughters(0)
  , fHisNITSclsvsEtaDaughters(0)
  , fHisSPDclsDaughters(0)
  , fHisPtAllTracks(0)
  , fHisPhiAllTracks(0)
  , fHisEtaAllTracks(0)
  , fHisEtavsPhiAllTracks(0)
  , fHisNTPCclsvsPtAllTracks(0)
  , fHisNTPCclsvsPhiAllTracks(0)
  , fHisNTPCclsvsEtaAllTracks(0)
  , fHisNTPCCrossedRowsvsPtAllTracks(0)
  , fHisNTPCCrossedRowsvsPhiAllTracks(0)
  , fHisNTPCCrossedRowsvsEtaAllTracks(0)
  , fHisRatioCRowsOverFclsvsPtAllTracks(0)
  , fHisRatioCRowsOverFclsvsPhiAllTracks(0)
  , fHisRatioCRowsOverFclsvsEtaAllTracks(0)
  , fHisNITSclsvsPtAllTracks(0)
  , fHisNITSclsvsPhiAllTracks(0)
  , fHisNITSclsvsEtaAllTracks(0)
  , fHisSPDclsAllTracks(0)
  , fHisdistrFakeTr(0)
  , fHisd0f(0)
  , fHisd0f_filt(0)
  , fHisptFakeTr(0)
  , fHisptFakeTrFromDaugh(0)
  , fHisptFakeTrFromDaughFilt(0)
  , fHisNtracklets(0)
  , fHisNtracklets01(0)
  , fHisNtracklets01AllEv(0)
  , fHisMult(0)
  , fHisMultFBit4(0)
  , fHisMultComb05(0)
  , fHisMultComb08(0)
  , fHisNtrackletsIn(0)
  , fHisMultIn(0)
  , fHisNtrackletsOut(0)
  , fHisMultOut(0)
  , fHisMultvsPercentile(0)
  , fHisntrklvsPercentile(0)
  , fHisntrklvsPercentile01(0)
  , fHisntrklvsPercentile01AllEv(0)
  , fHisnTPCTracksvsPercentile(0)
  , fHisnTPCITSTracksvsPercentile(0)
  , fHisnTPCITS1SPDTracksvsPercentile(0)
  , fHisStdEstimSignalPercentile(0)
  , fHisStdEstimSignalNtrackletsIn(0)
  , fHisStdEstimSignal(0)
  , fHisStdPercentileSecondPercentile(0)
  , fHisStdSignalSecondSignal(0)
  , fHisStdPercentileOldFrwPercentile(0)
  , fHisStdPercentileOldFrwPercentileDev(0)
  , fHisxvtx(0)
  , fHisyvtx(0)
  , fHiszvtx(0)
  , fHisxvtxSelEv(0)
  , fHisyvtxSelEv(0)
  , fHiszvtxSelEv(0)
  , fHisWhichVert(0)
  , fHisWhichVertSelEv(0)
  , fHisnClsITSvsNtrackletsSel(0)
  , fHiszvtxvsSPDzvtx(0)
  , fHiszvtxvsSPDzvtxSel(0)
  , fHisTrigCent(0)
  , fHisTrigMul(0)
  , fHisTrigCentSel(0)
  , fHisTrigMulSel(0)
  , fHisWhyEvRejected(0)
  , fHisFEvents(0)
  , fHisTPCVZE_AngleQ(0)
  , fHisCentVsMultRPS(0)
{
  /// default constructor
  fOnOff[0]=kTRUE;
  fOnOff[1]=kTRUE;
  fOnOff[2]=kTRUE;
  fOnOff[3]=kTRUE;
  fOnOff[4]=kTRUE;

  for (Int_t iii=0; iii<3; iii++) {
    fHisQ[iii]=0;
    fHisAngleQ[iii]=0;
    fHisPhiEta[iii]=0;
  }
}

//____________________________________________________________________________
AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA(const char *name, AliAnalysisTaskSEHFQA::DecChannel ch,AliRDHFCuts* cuts):
  AliAnalysisTaskSE(name)
  , fOutputEntries(0x0)
  , fOutputPID(0x0)
  , fOutputTrack(0x0)
  , fOutputCounters(0x0)
  , fOutputCheckCentrality(0x0)
  , fOutputEvSelection(0x0)
  , fOutputFlowObs(0x0)
  , fDecayChannel(ch)
  , fCuts(0x0)
  , fFlowEvent(0x0)
  , fRFPcuts(0x0)
  , fEstimator(AliRDHFCuts::kCentTRK)
  , fReadMC(kFALSE)
  , fSimpleMode(kFALSE)
  , fUseSelectionBit(kTRUE)
  , fOnOff()
  , fFillDistrTrackEffChecks(kFALSE)
  , fAODProtection(1)
  , fHisNentries(0)
  , fHisNentriesSelBit(0)
  , fHisTOFflags(0)
  , fHisTOFsig(0)
  , fHisTOFstartTimeMask(0)
  , fHisTOFstartTimeRes(0)
  , fHisTOFstartTimeDistrib(0)
  , fHisTOFtime(0)
  , fHisTOFtimeKaonHyptime(0)
  , fHisTOFtimeKaonHyptimeAC(0)
  , fHisTOFsigmaKSigPid(0)
  , fHisTOFsigmaPionSigPid(0)
  , fHisTOFsigmaProtonSigPid(0)
  , fHisTOFsigPid3sigPion(0)
  , fHisTOFsigPid3sigKaon(0)
  , fHisTOFsigPid3sigProton(0)
  , fHisTPCsig(0)
  , fHisTPCsigvsp(0)
  , fHisTPCsigvspAC(0)
  , fHisTPCsigmaK(0)
  , fHisTPCsigmaPion(0)
  , fHisTPCsigmaProton(0)
  , fHisTPCsigNvsPtAllTracks(0)
  , fHisTPCsigNvsPhiAllTracks(0)
  , fHisTPCsigNvsEtaAllTracks(0)
  , fHisTPCsigNvsPtDaughters(0)
  , fHisTPCsigNvsPhiDaughters(0)
  , fHisTPCsigNvsEtaDaughters(0)
  , fHisTOFsigmaMCKSigPid(0)
  , fHisTOFsigmaMCPionSigPid(0)
  , fHisTOFsigmaMCProtonSigPid(0)
  , fHisTPCsigmaMCK(0)
  , fHisTPCsigmaMCPion(0)
  , fHisTPCsigmaMCProton(0)
  , fHisnClsITS(0)
  , fHisnClsITSselTr(0)
  , fHisnClsITSSA(0)
  , fHisnClsITSSAspdAny(0)
  , fHisnClsITSSAspdIn(0)
  , fHisnClsITSSAspdOut(0)
  , fHisnLayerITS(0)
  , fHisnLayerITSselTr(0)
  , fHisnLayerITSsa(0)
  , fHisnClsSPD(0)
  , fHisptGoodTr(0)
  , fHisptGoodTrFromDaugh(0)
  , fHisptGoodTrFromDaugh_filt(0)
  , fHisdistrGoodTr(0)
  , fHisdistrSelTr(0)
  , fHisd0dau(0)
  , fHisd0dau_filt(0)
  , fHisd0dauphi(0)
  , fHisd0dauphi_filt(0)
  , fHisd0zdau(0)
  , fHisd0zdau_filt(0)
  , fHisd0zdauphi(0)
  , fHisd0zdauphi_filt(0)
  , fHisd0TracksSPDin(0)
  , fHisd0TracksSPDany(0)
  , fHisd0TracksFilterBit4(0)
  , fHisd0TracksTPCITSSPDany(0)
  , fHisPtDaughters(0)
  , fHisPhiDaughters(0)
  , fHisEtaDaughters(0)
  , fHisEtavsPhiDaughters(0)
  , fHisNTPCclsvsPtDaughters(0)
  , fHisNTPCclsvsPhiDaughters(0)
  , fHisNTPCclsvsEtaDaughters(0)
  , fHisNTPCCrossedRowsvsPtDaughters(0)
  , fHisNTPCCrossedRowsvsPhiDaughters(0)
  , fHisNTPCCrossedRowsvsEtaDaughters(0)
  , fHisRatioCRowsOverFclsvsPtDaughters(0)
  , fHisRatioCRowsOverFclsvsPhiDaughters(0)
  , fHisRatioCRowsOverFclsvsEtaDaughters(0)
  , fHisNITSclsvsPtDaughters(0)
  , fHisNITSclsvsPhiDaughters(0)
  , fHisNITSclsvsEtaDaughters(0)
  , fHisSPDclsDaughters(0)
  , fHisPtAllTracks(0)
  , fHisPhiAllTracks(0)
  , fHisEtaAllTracks(0)
  , fHisEtavsPhiAllTracks(0)
  , fHisNTPCclsvsPtAllTracks(0)
  , fHisNTPCclsvsPhiAllTracks(0)
  , fHisNTPCclsvsEtaAllTracks(0)
  , fHisNTPCCrossedRowsvsPtAllTracks(0)
  , fHisNTPCCrossedRowsvsPhiAllTracks(0)
  , fHisNTPCCrossedRowsvsEtaAllTracks(0)
  , fHisRatioCRowsOverFclsvsPtAllTracks(0)
  , fHisRatioCRowsOverFclsvsPhiAllTracks(0)
  , fHisRatioCRowsOverFclsvsEtaAllTracks(0)
  , fHisNITSclsvsPtAllTracks(0)
  , fHisNITSclsvsPhiAllTracks(0)
  , fHisNITSclsvsEtaAllTracks(0)
  , fHisSPDclsAllTracks(0)
  , fHisdistrFakeTr(0)
  , fHisd0f(0)
  , fHisd0f_filt(0)
  , fHisptFakeTr(0)
  , fHisptFakeTrFromDaugh(0)
  , fHisptFakeTrFromDaughFilt(0)
  , fHisNtracklets(0)
  , fHisNtracklets01(0)
  , fHisNtracklets01AllEv(0)
  , fHisMult(0)
  , fHisMultFBit4(0)
  , fHisMultComb05(0)
  , fHisMultComb08(0)
  , fHisNtrackletsIn(0)
  , fHisMultIn(0)
  , fHisNtrackletsOut(0)
  , fHisMultOut(0)
  , fHisMultvsPercentile(0)
  , fHisntrklvsPercentile(0)
  , fHisntrklvsPercentile01(0)
  , fHisntrklvsPercentile01AllEv(0)
  , fHisnTPCTracksvsPercentile(0)
  , fHisnTPCITSTracksvsPercentile(0)
  , fHisnTPCITS1SPDTracksvsPercentile(0)
  , fHisStdEstimSignalPercentile(0)
  , fHisStdEstimSignalNtrackletsIn(0)
  , fHisStdEstimSignal(0)
  , fHisStdPercentileSecondPercentile(0)
  , fHisStdSignalSecondSignal(0)
  , fHisStdPercentileOldFrwPercentile(0)
  , fHisStdPercentileOldFrwPercentileDev(0)
  , fHisxvtx(0)
  , fHisyvtx(0)
  , fHiszvtx(0)
  , fHisxvtxSelEv(0)
  , fHisyvtxSelEv(0)
  , fHiszvtxSelEv(0)
  , fHisWhichVert(0)
  , fHisWhichVertSelEv(0)
  , fHisnClsITSvsNtrackletsSel(0)
  , fHiszvtxvsSPDzvtx(0)
  , fHiszvtxvsSPDzvtxSel(0)
  , fHisTrigCent(0)
  , fHisTrigMul(0)
  , fHisTrigCentSel(0)
  , fHisTrigMulSel(0)
  , fHisWhyEvRejected(0)
  , fHisFEvents(0)
  , fHisTPCVZE_AngleQ(0)
  , fHisCentVsMultRPS(0)
{
  /// constructor

  //SetCutObject(cuts);
  fCuts=cuts;

  fOnOff[0]=kTRUE;
  fOnOff[1]=kTRUE;
  fOnOff[2]=kTRUE;
  fOnOff[3]=kTRUE;
  fOnOff[4]=kTRUE;

  for (Int_t iii=0; iii<3; iii++) {
    fHisQ[iii]=0;
    fHisAngleQ[iii]=0;
    fHisPhiEta[iii]=0;
  }

  // Output slot #1 writes into a TList container (number of events)
  DefineOutput(1,TList::Class());
  // Output slot #2 writes into a TList container (PID)
  if (fOnOff[1]) DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TList container (Tracks)
  if (fOnOff[0]) DefineOutput(3,TList::Class());  //My private output
  // Output slot #4 writes into a AliRDHFCuts container (cuts)
  switch(fDecayChannel){
  case 0:
    DefineOutput(4,AliRDHFCutsDplustoKpipi::Class());  //My private output
    break;
  case 1:
    DefineOutput(4,AliRDHFCutsD0toKpi::Class());  //My private output
    break;
  case 2:
    DefineOutput(4,AliRDHFCutsDStartoKpipi::Class());  //My private output
    break;
  case 3:
    DefineOutput(4,AliRDHFCutsDstoKKpi::Class());  //My private output
    break;
  case 4:
    DefineOutput(4,AliRDHFCutsD0toKpipipi::Class());  //My private output
    break;
  case 5:
    DefineOutput(4,AliRDHFCutsLctopKpi::Class());  //My private output
    break;
  case kLambdactoV0:
    DefineOutput(4,AliRDHFCutsLctoV0::Class());  //My private output
    break;
  }
  if (fOnOff[2]) {
    // Output slot #5 writes into a TList container (AliCounterCollection)
    DefineOutput(5,TList::Class());  //My private output
    // Output slot #6 writes into a TList container (TH1F)
    DefineOutput(6,TList::Class());  //My private output
  }

  if(fOnOff[3]) DefineOutput(7,TList::Class());  //My private output
  if(fOnOff[4]) DefineOutput(8,TList::Class());  //My private output

}

//___________________________________________________________________________
AliAnalysisTaskSEHFQA::~AliAnalysisTaskSEHFQA()
{
  /// destructor

  delete fOutputEntries;

  delete fOutputPID;

  delete fOutputTrack;

  delete fOutputCounters;

  delete fOutputCheckCentrality;

  delete fOutputEvSelection;

  if(fOnOff[4]) {
    delete fOutputFlowObs;
    delete fFlowEvent;
  }

}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::Init(){

  /// initialization
  if(fDebug > 1) printf("AnalysisTaskSEHFQA::Init() \n");
  AliRDHFCuts *copycut = 0x0;

  switch(fDecayChannel){
  case 0:
    {
      copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fCuts)));
    }
    break;
  case 1:
    {
      copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
    }
    break;
  case 2:
    {
      copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
    }
    break;
  case 3:
    {
      copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fCuts)));
    }
    break;
  case 4:
    {
      copycut=new AliRDHFCutsD0toKpipipi(*(static_cast<AliRDHFCutsD0toKpipipi*>(fCuts)));
    }
    break;
  case 5:
    {
      copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fCuts)));
    }
    break;
  case kLambdactoV0:
    {
      copycut=new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0*>(fCuts)));
    }
    break;
  default:
    AliFatal("Bad initialization for the decay channe - Exiting...");
    break;
  }

  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
  if (copycut){
    copycut->SetName(nameoutput);

    // Post the data
    PostData(4,copycut);
  }else{
    AliFatal("Failing initializing AliRDHFCuts object - Exiting...");
  }

  return;

}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::UserCreateOutputObjects()
{

  /// create the output container
  if(fDebug > 1) printf("AnalysisTaskSEHFQA::UserCreateOutputObjects() \n");

  //count events
  fOutputEntries=new TList();
  fOutputEntries->SetOwner();
  fOutputEntries->SetName(GetOutputSlot(1)->GetContainer()->GetName());


  TString hnameEntries="hNentries";
  fHisNentries=new TH1F(hnameEntries.Data(), "Counts the number of events", 15,-0.5,14.5);
  fHisNentries->GetXaxis()->SetBinLabel(1,"nEventsRead");
  fHisNentries->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
  fHisNentries->GetXaxis()->SetBinLabel(3,"Mismatched dAOD (Event numbers)");
  fHisNentries->GetXaxis()->SetBinLabel(4,"Mismatched dAOD (TProcessID)");
  fHisNentries->GetXaxis()->SetBinLabel(5,"Mismatched Old New Centrality");
  fHisNentries->GetXaxis()->SetBinLabel(6,"nEventsAnal");
  fHisNentries->GetXaxis()->SetBinLabel(7,"Pile-up Rej");
  fHisNentries->GetXaxis()->SetBinLabel(8,"No VertexingHF");
  fHisNentries->GetXaxis()->SetBinLabel(9,"nCandidates(AnCuts)");
  fHisNentries->GetXaxis()->SetBinLabel(10,"EventsWithGoodVtx");
  fHisNentries->GetXaxis()->SetBinLabel(11,"N candidates");
  if(fReadMC){
    fHisNentries->GetXaxis()->SetBinLabel(12,"MC Cand from c");
    fHisNentries->GetXaxis()->SetBinLabel(13,"MC Cand from b");
    fHisNentries->GetXaxis()->SetBinLabel(14,"N fake Trks");
    fHisNentries->GetXaxis()->SetBinLabel(15,"N true Trks");
  }

  fHisNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  hnameEntries="HasSelBit";
  fHisNentriesSelBit=new TH2F(hnameEntries.Data(), "Counts the number of events with SelectionBit", 5,0.,5.,100,0.,30.);
  fHisNentriesSelBit->GetXaxis()->SetBinLabel(1,"Dplus");
  fHisNentriesSelBit->GetXaxis()->SetBinLabel(2,"Ds");
  fHisNentriesSelBit->GetXaxis()->SetBinLabel(3,"LcKpi");
  fHisNentriesSelBit->GetXaxis()->SetBinLabel(4,"D0toKpi");
  fHisNentriesSelBit->GetXaxis()->SetBinLabel(5,"Dstar");

  fOutputEntries->Add(fHisNentries);
  fOutputEntries->Add(fHisNentriesSelBit);


  //PID
  if(fOnOff[1]){
    fOutputPID=new TList();
    fOutputPID->SetOwner();
    fOutputPID->SetName(GetOutputSlot(2)->GetContainer()->GetName());

    //TOF pid
    fHisTOFflags=new TH1F("hTOFflags","TOF flags",7,-0.5,6.5);
    fHisTOFflags->SetMinimum(0.);
    fHisTOFflags->GetXaxis()->SetBinLabel(1,"All Tracks");
    fHisTOFflags->GetXaxis()->SetBinLabel(2,"kTPCout");
    fHisTOFflags->GetXaxis()->SetBinLabel(3,"kTOFout");
    fHisTOFflags->GetXaxis()->SetBinLabel(4,"kTIME");
    fHisTOFflags->GetXaxis()->SetBinLabel(5,"kTOFpid");
    fHisTOFflags->GetXaxis()->SetBinLabel(6,"kTOFmismatch");
    fHisTOFflags->GetXaxis()->SetBinLabel(7,"kDetPidOK");

    TString hname="hTOFsig";
    fHisTOFsig=new TH1F(hname.Data(),"Distribution of TOF signal;TOF time [ps];Entries", 100, -2.e3,40.e3);

    hname="hTOFstartTimeMask";
    fHisTOFstartTimeMask=new TH1F(hname.Data(),"TOF start time mask; Mask ;Entries", 8, -0.5,7.5);
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(1,"FILL");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(2,"TOF");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(3,"T0A");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(4,"TOF.and.T0A");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(5,"T0C");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(6,"TOF.and.T0C");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(7,"T0AC");
    fHisTOFstartTimeMask->GetXaxis()->SetBinLabel(8,"TOF.and.T0AC");

    hname="hTOFstartTimeRes";
    fHisTOFstartTimeRes=new TH1F(hname.Data(),"TOF start time resolution; Resolution (ps) ;Entries", 100, 0.,300.);

    hname="hTOFstartTimeDistrib";
    fHisTOFstartTimeDistrib=new TH1F(hname.Data(),"TOF start time distribution; Start time ;Entries", 400, -1000.,1000.);

    hname="hTOFtime";
    fHisTOFtime=new TH1F(hname.Data(),"Distribution of TOF time Kaon;TOF time(Kaon) [ps];Entries", 1000, 0.,50000.);

    hname="hTOFtimeKaonHyptime";
    fHisTOFtimeKaonHyptime=new TH2F(hname.Data(),"TOFtime - timeHypothesisForKaon;p[GeV/c];TOFtime - timeHypothesisForKaon [ps]",500,0.,10.,1000,-20000.,20000.);

    hname="hTOFtimeKaonHyptimeAC";
    fHisTOFtimeKaonHyptimeAC=new TH2F(hname.Data(),"TOFtime - timeHypothesisForKaon;p[GeV/c];TOFtime - timeHypothesisForKaon [ps]",500,0.,10.,1000,-20000.,20000.);

    hname="hTOFsigmaKSigPid";
    fHisTOFsigmaKSigPid=new TH2F(hname.Data(),"(TOFsignal-timeK)/tofSigPid;p[GeV/c];(TOFsignal-timeK)/tofSigPid",500,0.,10.,400,-20,20);

    hname="hTOFsigmaPionSigPid";
    fHisTOFsigmaPionSigPid=new TH2F(hname.Data(),"(TOFsignal-time#pi)/tofSigPid;p[GeV/c];(TOFsignal-time#pi)/tofSigPid",500,0.,10.,400,-20,20);

    hname="hTOFsigmaProtonSigPid";
    fHisTOFsigmaProtonSigPid=new TH2F(hname.Data(),"(TOFsignal-timep)/tofSigPid;p[GeV/c];(TOFsignal-time p)/tofSigPid",500,0.,10.,400,-20,20);

    hname="hTOFsigPid3sigPion";
    fHisTOFsigPid3sigPion=new TH1F(hname.Data(),"TOF PID resolution (#pi) [ps]",500,0.,1000.);

    hname="hTOFsigPid3sigKaon";
    fHisTOFsigPid3sigKaon=new TH1F(hname.Data(),"TOF PID resolution (K) [ps]",500,0.,1000.);

    hname="hTOFsigPid3sigProton";
    fHisTOFsigPid3sigProton=new TH1F(hname.Data(),"TOF PID resolution (p) [ps]",500,0.,1000.);


    //TPC pid
    hname="hTPCsig";
    fHisTPCsig=new TH1F(hname.Data(),"Distribution of TPC signal;TPC sig;Entries", 100, 35.,100.);

    hname="hTPCsigvsp";
    fHisTPCsigvsp=new TH2F(hname.Data(),"TPCsig vs p;TPC p[GeV/c];TPCsig",500,0.,10.,1000,35.,100.);

    hname="hTPCsigvspAC";
    fHisTPCsigvspAC=new TH2F(hname.Data(),"TPCsig vs p;TPCp[GeV/c];TPCsig",500,0.,10.,1000,35.,100.);

    hname="hTPCsigmaK";
    fHisTPCsigmaK=new TH2F(hname.Data(),"TPC Sigma for K as a function of momentum;p[GeV/c];Sigma Kaon",500,0.,10.,400,-20,20);

    hname="hTPCsigmaPion";
    fHisTPCsigmaPion=new TH2F(hname.Data(),"TPC Sigma for #pi as a function of momentum;p[GeV/c];Sigma #pi",500,0.,10.,400,-20,20);

    hname="hTPCsigmaProton";
    fHisTPCsigmaProton=new TH2F(hname.Data(),"TPC Sigma for proton as a function of momentum;p[GeV/c];Sigma Proton",500,0.,10.,400,-20,20);


    fOutputPID->Add(fHisTOFflags);
    fOutputPID->Add(fHisTOFsig);
    fOutputPID->Add(fHisTPCsig);
    fOutputPID->Add(fHisTOFstartTimeMask);
    fOutputPID->Add(fHisTOFstartTimeRes);
    fOutputPID->Add(fHisTOFstartTimeDistrib);
    fOutputPID->Add(fHisTOFtime);
    fOutputPID->Add(fHisTOFtimeKaonHyptime);
    fOutputPID->Add(fHisTOFtimeKaonHyptimeAC);
    fOutputPID->Add(fHisTOFsigmaKSigPid);
    fOutputPID->Add(fHisTOFsigmaPionSigPid);
    fOutputPID->Add(fHisTOFsigmaProtonSigPid);
    fOutputPID->Add(fHisTOFsigPid3sigPion);
    fOutputPID->Add(fHisTOFsigPid3sigKaon);
    fOutputPID->Add(fHisTOFsigPid3sigProton);
    fOutputPID->Add(fHisTPCsigvsp);
    fOutputPID->Add(fHisTPCsigvspAC);
    fOutputPID->Add(fHisTPCsigmaK);
    fOutputPID->Add(fHisTPCsigmaPion);
    fOutputPID->Add(fHisTPCsigmaProton);

    if(fFillDistrTrackEffChecks){

      hname="hTPCsigNvsPtAllTracks";
      fHisTPCsigNvsPtAllTracks=new TH2F(hname.Data(),"Distribution of n. points used for TPC dE/dx vs. p_{T};p_{T} [GeV/c]; n. points", 200, 0.,20.,161,-0.5,160.5);

      hname="hTPCsigNvsPhiAllTracks";
      fHisTPCsigNvsPhiAllTracks=new TH2F(hname.Data(),"Distribution of n. points used for TPC dE/dx vs. #phi;#phi [rad]; n. points", 100, 0.,2*TMath::Pi(),161,-0.5,160.5);

      hname="hTPCsigNvsEtaAllTracks";
      fHisTPCsigNvsEtaAllTracks=new TH2F(hname.Data(),"Distribution of n. points used for TPC dE/dx vs. #eta;eta; n. points", 80,-2.,2.,161,-0.5,160.5);

      hname="hTPCsigNvsPtDaughters";
      fHisTPCsigNvsPtDaughters=new TH2F(hname.Data(),"Distribution of n. points used for TPC dE/dx vs. p_{T};p_{T} [GeV/c]; n. points", 200, 0.,20.,161,-0.5,160.5);

      hname="hTPCsigNvsPhiDaughters";
      fHisTPCsigNvsPhiDaughters=new TH2F(hname.Data(),"Distribution of n. points used for TPC dE/dx vs. #phi;#phi [rad]; n. points", 100, 0.,2*TMath::Pi(),161,-0.5,160.5);

      hname="hTPCsigNvsEtaDaughters";
      fHisTPCsigNvsEtaDaughters=new TH2F(hname.Data(),"Distribution of n. points used for TPC dE/dx vs. #eta;eta; n. points", 80,-2.,2.,161,-0.5,160.5);

      fOutputPID->Add(fHisTPCsigNvsPtAllTracks);
      fOutputPID->Add(fHisTPCsigNvsPhiAllTracks);
      fOutputPID->Add(fHisTPCsigNvsEtaAllTracks);
      fOutputPID->Add(fHisTPCsigNvsPtDaughters);
      fOutputPID->Add(fHisTPCsigNvsPhiDaughters);
      fOutputPID->Add(fHisTPCsigNvsEtaDaughters);
    }


    if(fReadMC){
      //TOF
      hname="hTOFsigmaMCKSigPid";
      fHisTOFsigmaMCKSigPid=new TH2F(hname.Data(),"(TOFsignal-timeK)/tofSigPid;p[GeV/c];(TOFsignal-timeK)/tofSigPid",500,0.,10.,400,-20,20);

      hname="hTOFsigmaMCPionSigPid";
      fHisTOFsigmaMCPionSigPid=new TH2F(hname.Data(),"(TOFsignal-time#pi)/tofSigPid;p[GeV/c];(TOFsignal-time#pi)/tofSigPid",500,0.,10.,400,-20,20);

      hname="hTOFsigmaMCProtonSigPid";
      fHisTOFsigmaMCProtonSigPid=new TH2F(hname.Data(),"(TOFsignal-timep)/tofSigPid;p[GeV/c];(TOFsignal-time p)/tofSigPid",500,0.,10.,400,-20,20);

      //TPC
      hname="hTPCsigmaMCK";
      fHisTPCsigmaMCK=new TH2F(hname.Data(),"TPC Sigma for K as a function of momentum;p[GeV/c];Sigma Kaon",500,0.,10.,400,-20,20);

      hname="hTPCsigmaMCPion";
      fHisTPCsigmaMCPion=new TH2F(hname.Data(),"TPC Sigma for #pi as a function of momentum;p[GeV/c];Sigma #pi",500,0.,10.,400,-20,20);

      hname="hTPCsigmaMCProton";
      fHisTPCsigmaMCProton=new TH2F(hname.Data(),"TPC Sigma for proton as a function of momentum;p[GeV/c];Sigma Proton",500,0.,10.,400,-20,20);

      fOutputPID->Add(fHisTOFsigmaMCKSigPid);
      fOutputPID->Add(fHisTOFsigmaMCPionSigPid);
      fOutputPID->Add(fHisTOFsigmaMCProtonSigPid);
      fOutputPID->Add(fHisTPCsigmaMCK);
      fOutputPID->Add(fHisTPCsigmaMCPion);
      fOutputPID->Add(fHisTPCsigmaMCProton);

    }
  }

  //quality of the tracks
  if(fOnOff[0]){
    fOutputTrack=new TList();
    fOutputTrack->SetOwner();
    fOutputTrack->SetName(GetOutputSlot(3)->GetContainer()->GetName());

    TString hname="hnClsITS";
    fHisnClsITS=new TH1F(hname.Data(),"Distribution of number of ITS clusters;nITScls;Entries",7,-0.5,6.5);

    hname="hnClsITSselTr";
    fHisnClsITSselTr=new TH1F(hname.Data(),"Distribution of number of ITS clusters selected tracks;nITScls;Entries",7,-0.5,6.5);

    hname="hnClsITS-SA";
    fHisnClsITSSA=new TH1F(hname.Data(),"Distribution of number of ITS clusters(ITS-SA);nITScls;Entries",7,-0.5,6.5);
    hname="hnClsITS-SA-SPDAny";
    fHisnClsITSSAspdAny=new TH1F(hname.Data(),"Distribution of number of ITS clusters(ITS-SA) - SPD kAny;nITScls;Entries",7,-0.5,6.5);
    hname="hnClsITS-SA-SPDIn";
    fHisnClsITSSAspdIn=new TH1F(hname.Data(),"Distribution of number of ITS clusters(ITS-SA) - SPDin;nITScls;Entries",7,-0.5,6.5);
    hname="hnClsITS-SA-SPDOut";
    fHisnClsITSSAspdOut=new TH1F(hname.Data(),"Distribution of number of ITS clusters(ITS-SA) - SPDout;nITScls;Entries",7,-0.5,6.5);


    hname="hnLayerITS";
    fHisnLayerITS=new TH1F(hname.Data(),"Number of tracks with point in layer;ITS layer;",7,-1.5,5.5);
    fHisnLayerITS->GetXaxis()->SetBinLabel(1,"n tracks");
    fHisnLayerITS->GetXaxis()->SetBinLabel(2,"SPDin");
    fHisnLayerITS->GetXaxis()->SetBinLabel(3,"SPDout");
    fHisnLayerITS->GetXaxis()->SetBinLabel(4,"SDDin");
    fHisnLayerITS->GetXaxis()->SetBinLabel(5,"SDDout");
    fHisnLayerITS->GetXaxis()->SetBinLabel(6,"SSDin");
    fHisnLayerITS->GetXaxis()->SetBinLabel(7,"SSDout");

    hname="hnLayerITSselTr";
    fHisnLayerITSselTr=new TH1F(hname.Data(),"Number of selected tracks with point in layer;ITS layer;",7,-1.5,5.5);
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(1,"n tracks");
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(2,"SPDin");
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(3,"SPDout");
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(4,"SDDin");
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(5,"SDDout");
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(6,"SSDin");
    fHisnLayerITSselTr->GetXaxis()->SetBinLabel(7,"SSDout");

    hname="hnLayerITSsa";
    fHisnLayerITSsa=new TH1F(hname.Data(),"Number of ITSsa tracks with point in layer;ITS layer;",7,-1.5,5.5);
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(1,"n tracks");
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(2,"SPDin");
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(3,"SPDout");
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(4,"SDDin");
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(5,"SDDout");
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(6,"SSDin");
    fHisnLayerITSsa->GetXaxis()->SetBinLabel(7,"SSDout");

    hname="hnClsSPD";
    fHisnClsSPD=new TH1F(hname.Data(),"Distribution of number of SPD clusters;nSPDcls;Entries",3,-0.5,2.5);

    hname="hptGoodTr";
    fHisptGoodTr=new TH1F(hname.Data(),"Pt distribution of 'good' tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
    fHisptGoodTr->SetTitleOffset(1.3,"Y");

    if(!fSimpleMode){
      hname="hptGoodTrFromDaugh";
      fHisptGoodTrFromDaugh=new TH1F(hname.Data(),"Pt distribution of 'good' candidate's daughters;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
      fHisptGoodTrFromDaugh->SetTitleOffset(1.3,"Y");
      fOutputTrack->Add(fHisptGoodTrFromDaugh);
      hname="hptGoodTrFromDaugh_filt";
      fHisptGoodTrFromDaugh_filt=new TH1F(hname.Data(),"Pt distribution of 'good' candidate's daughters, cuts level;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
      fHisptGoodTrFromDaugh_filt->SetTitleOffset(1.3,"Y");
      fOutputTrack->Add(fHisptGoodTrFromDaugh_filt);
    }

    hname="hdistrGoodTr";
    fHisdistrGoodTr=new TH1F(hname.Data(),"Distribution of number of 'good' candidate's daughters per event;no.good-tracks/ev;Entries",4000,-0.5,3999.5);
    fHisdistrGoodTr->SetTitleOffset(1.3,"Y");

    hname="hdistrSelTr";
    fHisdistrSelTr=new TH1F(hname.Data(),"Distribution of number of Selected tracks per event;no.good-tracks/ev;Entries",4000,-0.5,3999.5);
    fHisdistrSelTr->SetTitleOffset(1.3,"Y");

    hname="hd0dau";
    fHisd0dau=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of D daughter tracks;d_{0rphi}[cm];Entries/10^{3} cm",200,-0.1,0.1);

    hname="hd0dau_filt";
    fHisd0dau_filt=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of D daughter tracks, cut level;d_{0rphi}[cm];Entries/10^{3} cm",200,-0.1,0.1);

    hname="hd0dauphi";
    fHisd0dauphi=new TH2F(hname.Data(), "Impact parameter (rphi) distribution of D daughter tracks versus #phi; #phi [rad]; d_{0rphi} [cm]",400,0,6.3,200,-0.1,0.1);

    hname="hd0dauphi_filt";
    fHisd0dauphi_filt=new TH2F(hname.Data(), "Impact parameter (rphi) distribution of D daughter tracks versus #phi, cut level; #phi [rad]; d_{0rphi} [cm]",400,0,6.3,200,-0.1,0.1);

    hname="hd0zdau";
    fHisd0zdau=new TH1F(hname.Data(),"Impact parameter (z) distribution of D daughter tracks;d_{0z}[cm];Entries/10^{3} cm",200,-0.1,0.1);

    hname="hd0zdau_filt";
    fHisd0zdau_filt=new TH1F(hname.Data(),"Impact parameter (z) distribution of D daughter tracks, cut level;d_{0z}[cm];Entries/10^{3} cm",200,-0.1,0.1);


    hname="hd0zdauphi";
    fHisd0zdauphi=new TH2F(hname.Data(), "Impact parameter (z) distribution of D daughter tracks versus #phi; #phi [rad]; d_{0z} [cm]",400,0,6.3,200,-0.1,0.1);

    hname="hd0zdauphi_filt";
    fHisd0zdauphi_filt=new TH2F(hname.Data(), "Impact parameter (z) distribution of D daughter tracks versus #phi, filtering level; #phi [rad]; d_{0z} [cm]",400,0,6.3,200,-0.1,0.1);

    hname="hd0TracksSPDin";
    fHisd0TracksSPDin=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of AOD tracks kITSrefit, SPDinner; d_{0rphi}[cm];Entries",200,-0.5,0.5);

    hname="hd0TracksSPDany";
    fHisd0TracksSPDany=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of AOD tracks kITSrefit, SPDany; d_{0rphi}[cm];Entries",200,-0.5,0.5);

    hname="hd0TracksFilterBit4";
    fHisd0TracksFilterBit4=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of AOD tracks FilterBit4; d_{0rphi}[cm];Entries",200,-0.5,0.5);

    hname="hd0TracksTPCITSSPDany";
    fHisd0TracksTPCITSSPDany=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of AOD tracks TPC+ITScuts+SPDany; d_{0rphi}[cm];Entries",200,-0.5,0.5);


    if(fFillDistrTrackEffChecks){
      hname="hPtDaughters";
      fHisPtDaughters=new TH1F(hname.Data(),"p_{T} distributions of the daughter tracks;p_{T} [GeV/c];Entries",200,0.,20.);

      hname="hPhiDaughters";
      fHisPhiDaughters=new TH1F(hname.Data(),"#phi distribution of the daughter tracks;#phi [rad];Entries",100,0.,2*(TMath::Pi()));

      hname="hEtaDaughters";
      fHisEtaDaughters=new TH1F(hname.Data(),"#eta distribution of the daughter tracks;#eta;Entries",80,-2.,2.);

      hname="hEtavsPhiDaughters";
      fHisEtavsPhiDaughters=new TH2F(hname.Data(),"#eta vs #phi distribution of the daughter tracks;#phi;#eta",100,0.,2*(TMath::Pi()),80,-2.,2.);

      hname="hNTPCclsvsPtDaughters";
      fHisNTPCclsvsPtDaughters=new TH2F(hname.Data(),"N TPC clusters vs p_{T} distribution of the daughter tracks;p_{T} [GeV/c];N TPC cls",200,0.,20.,85,-0.5,169.5);

      hname="hNTPCclsvsPhiDaughters";
      fHisNTPCclsvsPhiDaughters=new TH2F(hname.Data(),"N TPC clusters vs #phi distribution of the daughter tracks;#phi [rad];N TPC cls",100,0.,2*(TMath::Pi()),85,-0.5,169.5);

      hname="hNTPCclsvsEtaDaughters";
      fHisNTPCclsvsEtaDaughters=new TH2F(hname.Data(),"N TPC clusters vs #eta distribution of the daughter tracks;#eta;N TPC cls",80,-2.,2.,85,-0.5,169.5);

      hname="hNTPCCrossedRowsvsPtDaughters";
      fHisNTPCCrossedRowsvsPtDaughters=new TH2F(hname.Data(),"N TPC crossed rows vs p_{T} distribution of the daughter tracks;p_{T} [GeV/c];N TPC cros. rows",200,0.,20.,100,-0.5,199.5);

      hname="hNTPCCrossedRowsvsPhiDaughters";
      fHisNTPCCrossedRowsvsPhiDaughters=new TH2F(hname.Data(),"N TPC crossed rows vs #phi distribution of the daughter tracks;#phi [rad];N TPC cros. rows",100,0.,2*(TMath::Pi()),100,-0.5,199.5);

      hname="hNTPCCrossedRowsvsEtaDaughters";
      fHisNTPCCrossedRowsvsEtaDaughters=new TH2F(hname.Data(),"N TPC crossed rows vs #eta distribution of the daughter tracks;#eta;N TPC cros. rows",80,-2.,2.,100,-0.5,199.5);

      hname="hRatioCRowsOverFclsvsPtDaughters";
      fHisRatioCRowsOverFclsvsPtDaughters=new TH2F(hname.Data(),"CrossedRows/FindableClusters vs p_{T} distribution of the daughter tracks;p_{T} [GeV/c];CRows/FCls",200,0.,20,100,0.,1.);

      hname="hRatioCRowsOverFclsvsPhiDaughters";
      fHisRatioCRowsOverFclsvsPhiDaughters=new TH2F(hname.Data(),"CrossedRows/FindableClusters vs #phi distribution of the daughter tracks;#phi [rad];CRows/FCls",100,0.,2*(TMath::Pi()),100,0.,1.);

      hname="hRatioCRowsOverFclsvsEtaDaughters";
      fHisRatioCRowsOverFclsvsEtaDaughters=new TH2F(hname.Data(),"CrossedRows/FindableClusters vs #eta distribution of the daughter tracks;#eta;CRows/FCls",80,-2.,2.,100,0.,1.);

      hname="hNITSclsvsPtDaughters";
      fHisNITSclsvsPtDaughters=new TH2F(hname.Data(),"N ITS clusters vs p_{T} distribution of the daughter tracks;p_{T} [GeV/c];N ITS cls",200,0.,20,7,-0.5,6.5);

      hname="hNITSclsvsPhiDaughters";
      fHisNITSclsvsPhiDaughters=new TH2F(hname.Data(),"N ITS clusters vs #phi distribution of the daughter tracks;#phi [rad];N ITS cls",100,0.,2*(TMath::Pi()),7,-0.5,6.5);

      hname="hNITSclsvsEtaDaughters";
      fHisNITSclsvsEtaDaughters=new TH2F(hname.Data(),"N ITS clusters vs #eta distribution of the daughter tracks;#eta;N ITS cls",80,-2.,2.,7,-0.5,6.5);

      hname="hSPDclsDaughters";
      fHisSPDclsDaughters = new TH1F(hname.Data(),"N SPD points distribution;;Entries",4,-0.5,3.5);
      fHisSPDclsDaughters->GetXaxis()->SetBinLabel(1, "no SPD");
      fHisSPDclsDaughters->GetXaxis()->SetBinLabel(2, "kOnlyFirst");
      fHisSPDclsDaughters->GetXaxis()->SetBinLabel(3, "kOnlySecond");
      fHisSPDclsDaughters->GetXaxis()->SetBinLabel(4, "kBoth");

      hname="hPtAllTracks";
      fHisPtAllTracks=new TH1F(hname.Data(),"p_{T} distributions of the AOD tracks (ID>0);p_{T} [GeV/c];Entries",200,0.,20.);

      hname="hPhiAllTracks";
      fHisPhiAllTracks=new TH1F(hname.Data(),"#phi distribution of the AOD tracks (ID>0);#phi [rad];Entries",100,0.,2*(TMath::Pi()));

      hname="hEtaAllTracks";
      fHisEtaAllTracks=new TH1F(hname.Data(),"#eta distribution of the AOD tracks (ID>0);#eta;Entries",80,-2.,2.);

      hname="hEtavsPhiAllTracks";
      fHisEtavsPhiAllTracks=new TH2F(hname.Data(),"#eta vs #phi distribution of the AOD tracks (ID>0);#phi;#eta",100,0.,2*(TMath::Pi()),80,-2.,2.);

      hname="hNTPCclsvsPtAllTracks";
      fHisNTPCclsvsPtAllTracks=new TH2F(hname.Data(),"N TPC clusters vs p_{T} distribution of the AOD tracks (ID>0);p_{T} [GeV/c];N TPC cls",200,0.,20,85,-0.5,169.5);

      hname="hNTPCclsvsPhiAllTracks";
      fHisNTPCclsvsPhiAllTracks=new TH2F(hname.Data(),"N TPC clusters vs #phi distribution of the AOD tracks (ID>0);#phi [rad];N TPC cls",100,0.,2*(TMath::Pi()),85,-0.5,169.5);

      hname="hNTPCclsvsEtaAllTracks";
      fHisNTPCclsvsEtaAllTracks=new TH2F(hname.Data(),"N TPC clusters vs #eta distribution of the AOD tracks (ID>0);#eta;N TPC cls",80,-2.,2.,85,-0.5,169.5);

      hname="hNTPCCrossedRowsvsPtAllTracks";
      fHisNTPCCrossedRowsvsPtAllTracks=new TH2F(hname.Data(),"N TPC crossed rows vs p_{T} distribution of the AOD tracks;p_{T} [GeV/c];N TPC cros. rows",200,0.,20.,100,-0.5,199.5);

      hname="hNTPCCrossedRowsvsPhiAllTracks";
      fHisNTPCCrossedRowsvsPhiAllTracks=new TH2F(hname.Data(),"N TPC crossed rows vs #phi distribution of the AOD tracks;#phi [rad];N TPC cros. rows",100,0.,2*(TMath::Pi()),100,-0.5,199.5);

      hname="hNTPCCrossedRowsvsEtaAllTracks";
      fHisNTPCCrossedRowsvsEtaAllTracks=new TH2F(hname.Data(),"N TPC crossed rows vs #eta distribution of the AOD tracks;#eta;N TPC cros. rows",80,-2.,2.,100,-0.5,199.5);

      hname="hRatioCRowsOverFclsvsPtAllTracks";
      fHisRatioCRowsOverFclsvsPtAllTracks=new TH2F(hname.Data(),"CrossedRows/FindableClusters vs p_{T} distribution of the AOD tracks (ID>0);p_{T} [GeV/c];CRows/FCls",200,0.,20,100,0.,1.);

      hname="hRatioCRowsOverFclsvsPhiAllTracks";
      fHisRatioCRowsOverFclsvsPhiAllTracks=new TH2F(hname.Data(),"CrossedRows/FindableClusters vs #phi distribution of the AOD tracks (ID>0);#phi [rad];CRows/FCls",100,0.,2*(TMath::Pi()),100,0.,1.);

      hname="hRatioCRowsOverFclsvsEtaAllTracks";
      fHisRatioCRowsOverFclsvsEtaAllTracks=new TH2F(hname.Data(),"CrossedRows/FindableClusters vs #eta distribution of the AOD tracks (ID>0);#eta;CRows/FCls",80,-2.,2.,100,0.,1.);

      hname="hNITSclsvsPtAllTracks";
      fHisNITSclsvsPtAllTracks=new TH2F(hname.Data(),"N ITS clusters vs p_{T} distribution of the AOD tracks (ID>0);p_{T} [GeV/c];N ITS cls",200,0.,20,7,-0.5,6.5);

      hname="hNITSclsvsPhiAllTracks";
      fHisNITSclsvsPhiAllTracks=new TH2F(hname.Data(),"N ITS clusters vs #phi distribution of the AOD tracks (ID>0);#phi [rad];N ITS cls",100,0.,2*(TMath::Pi()),7,-0.5,6.5);

      hname="hNITSclsvsEtaAllTracks";
      fHisNITSclsvsEtaAllTracks=new TH2F(hname.Data(),"N ITS clusters vs #eta distribution of the AOD tracks (ID>0);#eta;N ITS cls",80,-2.,2.,7,-0.5,6.5);

      hname="hSPDclsAllTracks";
      fHisSPDclsAllTracks = new TH1F(hname.Data(),"N SPD points distribution AOD tracks (ID>0);;Entries",4,-0.5,3.5);
      fHisSPDclsAllTracks->GetXaxis()->SetBinLabel(1, "no SPD");
      fHisSPDclsAllTracks->GetXaxis()->SetBinLabel(2, "kOnlyFirst");
      fHisSPDclsAllTracks->GetXaxis()->SetBinLabel(3, "kOnlySecond");
      fHisSPDclsAllTracks->GetXaxis()->SetBinLabel(4, "kBoth");


      fOutputTrack->Add(fHisPtDaughters);
      fOutputTrack->Add(fHisPhiDaughters);
      fOutputTrack->Add(fHisEtaDaughters);
      fOutputTrack->Add(fHisEtavsPhiDaughters);
      fOutputTrack->Add(fHisNTPCclsvsPtDaughters);
      fOutputTrack->Add(fHisNTPCclsvsPhiDaughters);
      fOutputTrack->Add(fHisNTPCclsvsEtaDaughters);
      fOutputTrack->Add(fHisNTPCCrossedRowsvsPtDaughters);
      fOutputTrack->Add(fHisNTPCCrossedRowsvsPhiDaughters);
      fOutputTrack->Add(fHisNTPCCrossedRowsvsEtaDaughters);
      fOutputTrack->Add(fHisRatioCRowsOverFclsvsPtDaughters);
      fOutputTrack->Add(fHisRatioCRowsOverFclsvsPhiDaughters);
      fOutputTrack->Add(fHisRatioCRowsOverFclsvsEtaDaughters);
      fOutputTrack->Add(fHisNITSclsvsPtDaughters);
      fOutputTrack->Add(fHisNITSclsvsPhiDaughters);
      fOutputTrack->Add(fHisNITSclsvsEtaDaughters);
      fOutputTrack->Add(fHisSPDclsDaughters);
      fOutputTrack->Add(fHisPtAllTracks);
      fOutputTrack->Add(fHisPhiAllTracks);
      fOutputTrack->Add(fHisEtaAllTracks);
      fOutputTrack->Add(fHisEtavsPhiAllTracks);
      fOutputTrack->Add(fHisNTPCclsvsPtAllTracks);
      fOutputTrack->Add(fHisNTPCclsvsPhiAllTracks);
      fOutputTrack->Add(fHisNTPCclsvsEtaAllTracks);
      fOutputTrack->Add(fHisNTPCCrossedRowsvsPtAllTracks);
      fOutputTrack->Add(fHisNTPCCrossedRowsvsPhiAllTracks);
      fOutputTrack->Add(fHisNTPCCrossedRowsvsEtaAllTracks);
      fOutputTrack->Add(fHisRatioCRowsOverFclsvsPtAllTracks);
      fOutputTrack->Add(fHisRatioCRowsOverFclsvsPhiAllTracks);
      fOutputTrack->Add(fHisRatioCRowsOverFclsvsEtaAllTracks);
      fOutputTrack->Add(fHisNITSclsvsPtAllTracks);
      fOutputTrack->Add(fHisNITSclsvsPhiAllTracks);
      fOutputTrack->Add(fHisNITSclsvsEtaAllTracks);
      fOutputTrack->Add(fHisSPDclsAllTracks);

    }

    fOutputTrack->Add(fHisnClsITS);
    fOutputTrack->Add(fHisnClsITSselTr);
    fOutputTrack->Add(fHisnClsITSSA);
    fOutputTrack->Add(fHisnClsITSSAspdAny);
    fOutputTrack->Add(fHisnClsITSSAspdIn);
    fOutputTrack->Add(fHisnClsITSSAspdOut);
    fOutputTrack->Add(fHisnLayerITS);
    fOutputTrack->Add(fHisnLayerITSselTr);
    fOutputTrack->Add(fHisnLayerITSsa);
    fOutputTrack->Add(fHisnClsSPD);
    fOutputTrack->Add(fHisptGoodTr);
    fOutputTrack->Add(fHisdistrGoodTr);
    fOutputTrack->Add(fHisdistrSelTr);
    fOutputTrack->Add(fHisd0TracksSPDin);
    fOutputTrack->Add(fHisd0TracksSPDany);
    fOutputTrack->Add(fHisd0TracksFilterBit4);
    fOutputTrack->Add(fHisd0TracksTPCITSSPDany);
    fOutputTrack->Add(fHisd0dau);
    fOutputTrack->Add(fHisd0dauphi);
    fOutputTrack->Add(fHisd0zdau);
    fOutputTrack->Add(fHisd0zdauphi);
    fOutputTrack->Add(fHisd0dau_filt);
    fOutputTrack->Add(fHisd0dauphi_filt);
    fOutputTrack->Add(fHisd0zdau_filt);
    fOutputTrack->Add(fHisd0zdauphi_filt);


    if(fReadMC){
      hname="hdistrFakeTr";
      fHisdistrFakeTr=new TH1F(hname.Data(),"Distribution of number of fake tracks per event;no.fake-tracks/ev;Entries",4000,-0.5,3999.5);
      fHisdistrFakeTr->SetTitleOffset(1.3,"Y");

      hname="hd0f";
      fHisd0f=new TH1F(hname.Data(),"Impact parameter distribution of fake tracks;d_{0}[cm];Entries/10^{3} cm",200,-0.1,0.1);

      hname="hd0f_filt";
      fHisd0f_filt=new TH1F(hname.Data(),"Impact parameter distribution of fake tracks, cut level;d_{0}[cm];Entries/10^{3} cm",200,-0.1,0.1);


      hname="hptFakeTr";
      fHisptFakeTr=new TH1F(hname.Data(),"Pt distribution of fake tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
      fHisptFakeTr->SetTitleOffset(1.3,"Y");
      if(!fSimpleMode){
	hname="hptFakeTrFromDaugh";
	fHisptFakeTrFromDaugh=new TH1F(hname.Data(),"Pt distribution of fake tracks from daughters;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
	fHisptFakeTrFromDaugh->SetTitleOffset(1.3,"Y");
	fOutputTrack->Add(fHisptFakeTrFromDaugh);

	hname="hptFakeTrFromDaugh_filt";
	fHisptFakeTrFromDaughFilt=new TH1F(hname.Data(),"Pt distribution of fake tracks from daughters, cut level;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
	fHisptFakeTrFromDaughFilt->SetTitleOffset(1.3,"Y");
	fOutputTrack->Add(fHisptFakeTrFromDaughFilt);
      }

      fOutputTrack->Add(fHisptFakeTr);
      fOutputTrack->Add(fHisdistrFakeTr);
      fOutputTrack->Add(fHisd0f);
      fOutputTrack->Add(fHisd0f_filt);
    }
  }


  if(fOnOff[2] && fCuts->GetUseCentrality()){

    //Centrality (Counters)
    fOutputCounters=new TList();
    fOutputCounters->SetOwner();
    fOutputCounters->SetName(GetOutputSlot(5)->GetContainer()->GetName());

    AliCounterCollection *stdEstimator=new AliCounterCollection("stdEstimator");
    stdEstimator->AddRubric("run",500000);
    stdEstimator->AddRubric("centralityclass","-10_0/0_10/10_20/20_30/30_40/40_50/50_60/60_70/70_80/80_90/90_100/-990_-980");
    stdEstimator->Init();
    AliCounterCollection *secondEstimator=new AliCounterCollection("secondEstimator");
    secondEstimator->AddRubric("run",500000);
    secondEstimator->AddRubric("centralityclass","-10_0/0_10/10_20/20_30/30_40/40_50/50_60/60_70/70_80/80_90/90_100/-990_-980");
    secondEstimator->Init();

    fOutputCounters->Add(stdEstimator);
    fOutputCounters->Add(secondEstimator);

    //Centrality (Checks)
    fOutputCheckCentrality=new TList();
    fOutputCheckCentrality->SetOwner();
    fOutputCheckCentrality->SetName(GetOutputSlot(6)->GetContainer()->GetName());

    TString hname="hNtrackletsIn";
    fHisNtrackletsIn=new TH1F(hname.Data(),"Number of tracklets in Centrality range;ntracklets;Entries",5000,-0.5,4999.5);

    hname="hMultIn";
    fHisMultIn=new TH1F(hname.Data(),"Multiplicity;multiplicity in Centrality range;Entries",10000,-0.5,9999.5);

    hname="hNtrackletsOut";
    fHisNtrackletsOut=new TH1F(hname.Data(),"Number of tracklets out of Centrality range;ntracklets;Entries",5000,-0.5,4999.5);

    hname="hMultOut";
    fHisMultOut=new TH1F(hname.Data(),"Multiplicity out of Centrality range;multiplicity;Entries",10000,-0.5,9999.5);

    hname="hMultvsPercentile";
    fHisMultvsPercentile=new TH2F(hname.Data(),"Multiplicity vs Percentile;multiplicity;percentile",10000,-0.5,9999.5,240,-10.,110);

    hname="hntrklvsPercentile";
    fHisntrklvsPercentile=new TH2F(hname.Data(),"N tracklets vs Percentile;ntracklets;percentile",5000,-0.5,4999.5,240,-10.,110);

    hname="hntrklvsPercentile01";
    fHisntrklvsPercentile01=new TH2F(hname.Data(),"N tracklets vs Percentile |#eta|<1;ntracklets;percentile",5000,-0.5,4999.5,240,-10.,110);

    hname="hntrklvsPercentile01AllEv";
    fHisntrklvsPercentile01AllEv=new TH2F(hname.Data(),"N tracklets vs Percentile |#eta|<1 - All Events;ntracklets;percentile",5000,-0.5,4999.5,240,-10.,110);

    hname="hnTPCTracksvsPercentile";
    fHisnTPCTracksvsPercentile=new TH2F(hname.Data(),"N TPC tracks vs Percentile;nTPCTracks;percentile",5000,-0.5,9999.5,240,-10.,110);

    hname="hnTPCITSTracksvsPercentile";
    fHisnTPCITSTracksvsPercentile=new TH2F(hname.Data(),"N TPC+ITS tracks vs Percentile;nTPCITSTracks;percentile",5000,-0.5,9999.5,240,-10.,110);

    hname="hnTPCITS1SPDTracksvsPercentile";
    fHisnTPCITS1SPDTracksvsPercentile=new TH2F(hname.Data(),"N TPC+ITS+1SPD tracks vs Percentile;nTPCITS1SPDTracks;percentile",5000,-0.5,9999.5,240,-10.,110);

    hname="hStdEstimSignalPercentile";
    fHisStdEstimSignalPercentile = new TH2F(hname.Data(),"Std estimator signal vs Percentile;Std estimator signal;percentile",1000,-0.5,9999.5,120,-10.,110);

    hname="hStdEstimSignalNtrackletsIn";
    fHisStdEstimSignalNtrackletsIn = new TH2F(hname.Data(),"Std estimator signal vs Number of tracklets in the CC;Std estimator signal;number of tracklets",1000,-0.5,9999.5,5000,-0.5,4999.5);

    hname="hStdEstimSignal";
    fHisStdEstimSignal = new TH1F(hname.Data(),"Std estimator signal",700,0,1400);

    hname="hStdPercentileSecondPercentile";
    fHisStdPercentileSecondPercentile = new TH2F(hname.Data(),"Std estimator Percentile Vs Second Estimator Percentile;Std estimator percentile;Second estimator percentile",120,-10.,110,120,-10.,110);

    hname="hStdSignalSecondSignal";
    fHisStdSignalSecondSignal = new TH2F(hname.Data(),"Std estimator signal Vs Second Estimator signal;Std estimator;Second estimator",1000,-0.5,9999.5,1000,-0.5,9999.5);

    hname="hStdPercentileOldFrwPercentile";
    fHisStdPercentileOldFrwPercentile = new TH2F(hname.Data(),"Std estimator Percentile Vs Old Framework estimator Percentile;Std estimator percentile;Old framework estimator percentile",120,-10.,110,120,-10.,110);

    hname="hStdPercentileOldFrwPercentileDev";
    fHisStdPercentileOldFrwPercentileDev = new TH1F(hname.Data(),"Std estimator Percentile Vs Old Framework estimator Percentile deviation",5,-0.5,4.5);
    fHisStdPercentileOldFrwPercentileDev->GetXaxis()->SetBinLabel(1, "nEvts with >1% difference");
    fHisStdPercentileOldFrwPercentileDev->GetXaxis()->SetBinLabel(2, "nEvts with >0.5% difference");
    fHisStdPercentileOldFrwPercentileDev->GetXaxis()->SetBinLabel(3, "nEvts default <20, new >20");
    fHisStdPercentileOldFrwPercentileDev->GetXaxis()->SetBinLabel(4, "nEvts default >20, new <20");
    fHisStdPercentileOldFrwPercentileDev->GetXaxis()->SetBinLabel(5, "nEvts rejected due to centrality mismatch");


    fOutputCheckCentrality->Add(fHisNtrackletsIn);
    fOutputCheckCentrality->Add(fHisNtrackletsOut);
    fOutputCheckCentrality->Add(fHisMultIn);
    fOutputCheckCentrality->Add(fHisMultOut);
    fOutputCheckCentrality->Add(fHisMultvsPercentile);
    fOutputCheckCentrality->Add(fHisntrklvsPercentile);
    fOutputCheckCentrality->Add(fHisntrklvsPercentile01);
    fOutputCheckCentrality->Add(fHisntrklvsPercentile01AllEv);
    fOutputCheckCentrality->Add(fHisnTPCTracksvsPercentile);
    fOutputCheckCentrality->Add(fHisnTPCITSTracksvsPercentile);
    fOutputCheckCentrality->Add(fHisnTPCITS1SPDTracksvsPercentile);
    fOutputCheckCentrality->Add(fHisStdEstimSignalPercentile);
    fOutputCheckCentrality->Add(fHisStdEstimSignal);
    fOutputCheckCentrality->Add(fHisStdEstimSignalNtrackletsIn);
    fOutputCheckCentrality->Add(fHisStdPercentileSecondPercentile);
    fOutputCheckCentrality->Add(fHisStdSignalSecondSignal);
    fOutputCheckCentrality->Add(fHisStdPercentileOldFrwPercentile);
    fOutputCheckCentrality->Add(fHisStdPercentileOldFrwPercentileDev);

    PostData(6,fOutputCheckCentrality);

  } else{
    if(fOnOff[0]){
      TString hname="hNtracklets";
      fHisNtracklets=new TH1F(hname.Data(),"Number of tracklets;ntracklets;Entries",5000,-0.5,4999.5);
      hname="hNtracklets01";
      fHisNtracklets01=new TH1F(hname.Data(),"Number of tracklets |#eta|<1;ntracklets;Entries",5000,-0.5,4999.5);
      hname="hNtracklets01AllEv";
      fHisNtracklets01AllEv=new TH1F(hname.Data(),"Number of tracklets |#eta|<1 - All events;ntracklets;Entries",5000,-0.5,4999.5);
      hname="hMult";
      fHisMult=new TH1F(hname.Data(),"Multiplicity;multiplicity;Entries",10000,-0.5,9999.5);
      hname="hMultFBit4";
      fHisMultFBit4=new TH1F(hname.Data(),"Multiplicity (global+tracklet) with filter bit 4;multiplicity;Entries",10000,-0.5,9999.5);
      hname="hMultComb05";
      fHisMultComb05=new TH1F(hname.Data(),"Multiplicity (global+tracklet) in |#eta|<0.5;multiplicity;Entries",10000,-0.5,9999.5);
      hname="hMultComb08";
      fHisMultComb08=new TH1F(hname.Data(),"Multiplicity (global+tracklet) in |#eta|<0.8;multiplicity;Entries",10000,-0.5,9999.5);

      fOutputTrack->Add(fHisNtracklets);
      fOutputTrack->Add(fHisNtracklets01);
      fOutputTrack->Add(fHisNtracklets01AllEv);
      fOutputTrack->Add(fHisMult);
      fOutputTrack->Add(fHisMultFBit4);
      fOutputTrack->Add(fHisMultComb05);
      fOutputTrack->Add(fHisMultComb08);
    }
  }

  //event selection (z vertex for the moment)
  if(fOnOff[3]){
    fOutputEvSelection=new TList();
    fOutputEvSelection->SetOwner();
    fOutputEvSelection->SetName(GetOutputSlot(7)->GetContainer()->GetName());
    AliCounterCollection *evselection=new AliCounterCollection("evselection");
    evselection->AddRubric("run",500000);
    evselection->AddRubric("evnonsel","zvtx");
    evselection->Init();

    fHisxvtx=new TH1F("hxvtx", "Distribution of x_{VTX};x_{VTX} [cm];Entries",800,-1,1);
    fHisyvtx=new TH1F("hyvtx", "Distribution of y_{VTX};y_{VTX} [cm];Entries",800,-1,1);
    fHiszvtx=new TH1F("hzvtx", "Distribution of z_{VTX};z_{VTX} [cm];Entries",800,-30,30);
    fHisxvtxSelEv=new TH1F("hxvtxSelEv", "Distribution of x_{VTX} Selected Ev;x_{VTX} [cm];Entries",800,-1,1);
    fHisyvtxSelEv=new TH1F("hyvtxSelEv", "Distribution of y_{VTX} Selected Ev;y_{VTX} [cm];Entries",800,-1,1);
    fHiszvtxSelEv=new TH1F("hzvtxSelEv", "Distribution of z_{VTX} Selected Ev;z_{VTX} [cm];Entries",800,-30,30);
    fHisWhichVert=new TH1F("hWhichVert","Vertex Type",4,-1.5,2.5);
    fHisWhichVert->GetXaxis()->SetBinLabel(1,"Not found");
    fHisWhichVert->GetXaxis()->SetBinLabel(2,"Track");
    fHisWhichVert->GetXaxis()->SetBinLabel(3,"SPD-3D");
    fHisWhichVert->GetXaxis()->SetBinLabel(4,"SPD-z");
    fHisWhichVertSelEv=new TH1F("hWhichVertSelEv","Vertex Type",4,-1.5,2.5);
    fHisWhichVertSelEv->GetXaxis()->SetBinLabel(1,"Not found");
    fHisWhichVertSelEv->GetXaxis()->SetBinLabel(2,"Track");
    fHisWhichVertSelEv->GetXaxis()->SetBinLabel(3,"SPD-3D");
    fHisWhichVertSelEv->GetXaxis()->SetBinLabel(4,"SPD-z");

    fHisnClsITSvsNtrackletsSel=new TH2F("hnClsITSvsNtrackletsSel","number of SPD clusters vs number of SPD tracklets; n. SPD clusters; Ntracklets",200,0,6000,500,0,20000); // max values should be changed for pp data to about 200 and 1000 respectively
    fHiszvtxvsSPDzvtx=new TH2F("hzvtxvsSPDzvtx","event primary z-vertex vs SPD z-vertex - before event selection; PV z-vertex [cm]; SPD z-vertex [cm]",800,-30,30,800,-30,30);
    fHiszvtxvsSPDzvtxSel=new TH2F("hzvtxvsSPDzvtxSel","event primary z-vertex vs SPD z-vertex - after event selection; PV z-vertex [cm]; SPD z-vertex [cm]",800,-30,30,800,-30,30);

    fHisTrigCent=new TH2F("hTrigCent","Centrality vs. Trigger types",24,-1.5,22.5,12,-10,110);
    fHisTrigCent->GetXaxis()->SetBinLabel(1,"All");
    fHisTrigCent->GetXaxis()->SetBinLabel(2,"kAny");
    fHisTrigCent->GetXaxis()->SetBinLabel(3,"kMB");
    fHisTrigCent->GetXaxis()->SetBinLabel(4,"kINT7");
    fHisTrigCent->GetXaxis()->SetBinLabel(5,"kINT8");
    fHisTrigCent->GetXaxis()->SetBinLabel(6,"kCINT5");
    fHisTrigCent->GetXaxis()->SetBinLabel(7,"kCent");
    fHisTrigCent->GetXaxis()->SetBinLabel(8,"kSemiCent");
    fHisTrigCent->GetXaxis()->SetBinLabel(9,"kEMC1");
    fHisTrigCent->GetXaxis()->SetBinLabel(10,"kEMC7");
    fHisTrigCent->GetXaxis()->SetBinLabel(11,"kEMC8");
    fHisTrigCent->GetXaxis()->SetBinLabel(12,"kEMCJET7");
    fHisTrigCent->GetXaxis()->SetBinLabel(13,"kEMCGAMMA7");
    fHisTrigCent->GetXaxis()->SetBinLabel(14,"kEMCJET8");
    fHisTrigCent->GetXaxis()->SetBinLabel(15,"kEMCGAMMA8");
    fHisTrigCent->GetXaxis()->SetBinLabel(16,"Muons");
    fHisTrigCent->GetXaxis()->SetBinLabel(17,"PHOS");
    fHisTrigCent->GetXaxis()->SetBinLabel(18,"TRD");
    fHisTrigCent->GetXaxis()->SetBinLabel(19,"TRDHJT");
    fHisTrigCent->GetXaxis()->SetBinLabel(20,"TRDHSE");
    fHisTrigCent->GetXaxis()->SetBinLabel(21,"HighMult");
    fHisTrigCent->GetXaxis()->SetBinLabel(22,"SPI7");
    fHisTrigCent->GetXaxis()->SetBinLabel(23,"SPI8");
    fHisTrigCent->GetXaxis()->SetBinLabel(24,"Others");

    fHisTrigMul=new TH2F("hTrigMul","Multiplicity vs. Trigger types",24,-1.5,22.5,1000,0.,10000.);
    fHisTrigMul->GetXaxis()->SetBinLabel(1,"All");
    fHisTrigMul->GetXaxis()->SetBinLabel(2,"kAny");
    fHisTrigMul->GetXaxis()->SetBinLabel(3,"kMB");
    fHisTrigMul->GetXaxis()->SetBinLabel(4,"kINT7");
    fHisTrigMul->GetXaxis()->SetBinLabel(5,"kINT8");
    fHisTrigMul->GetXaxis()->SetBinLabel(6,"kCINT5");
    fHisTrigMul->GetXaxis()->SetBinLabel(7,"kCent");
    fHisTrigMul->GetXaxis()->SetBinLabel(8,"kSemiCent");
    fHisTrigMul->GetXaxis()->SetBinLabel(9,"kEMC1");
    fHisTrigMul->GetXaxis()->SetBinLabel(10,"kEMC7");
    fHisTrigMul->GetXaxis()->SetBinLabel(11,"kEMC8");
    fHisTrigMul->GetXaxis()->SetBinLabel(12,"kEMCJET7");
    fHisTrigMul->GetXaxis()->SetBinLabel(13,"kEMCGAMMA7");
    fHisTrigMul->GetXaxis()->SetBinLabel(14,"kEMCJET8");
    fHisTrigMul->GetXaxis()->SetBinLabel(15,"kEMCGAMMA8");
    fHisTrigMul->GetXaxis()->SetBinLabel(16,"Muons");
    fHisTrigMul->GetXaxis()->SetBinLabel(17,"PHOS");
    fHisTrigMul->GetXaxis()->SetBinLabel(18,"TRD");
    fHisTrigMul->GetXaxis()->SetBinLabel(19,"TRDHJT");
    fHisTrigMul->GetXaxis()->SetBinLabel(20,"TRDHSE");
    fHisTrigMul->GetXaxis()->SetBinLabel(21,"HighMult");
    fHisTrigMul->GetXaxis()->SetBinLabel(22,"SPI7");
    fHisTrigMul->GetXaxis()->SetBinLabel(23,"SPI8");
    fHisTrigMul->GetXaxis()->SetBinLabel(24,"Others");

    fHisTrigCentSel=new TH2F("hTrigCentSel","Trigger types",24,-1.5,22.5,12,-10,110);
    fHisTrigCentSel->GetXaxis()->SetBinLabel(1,"All");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(2,"kAny");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(3,"kMB");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(4,"kINT7");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(5,"kINT8");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(6,"kCINT5");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(7,"kCent");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(8,"kSemiCent");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(9,"kEMC1");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(10,"kEMC7");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(11,"kEMC8");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(12,"kEMCJET7");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(13,"kEMCGAMMA7");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(14,"kEMCJET8");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(15,"kEMCGAMMA8");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(16,"Muons");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(17,"PHOS");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(18,"TRD");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(19,"TRDHJT");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(20,"TRDHSE");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(21,"HighMult");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(22,"SPI7");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(23,"SPI8");
    fHisTrigCentSel->GetXaxis()->SetBinLabel(24,"Others");

    fHisTrigMulSel=new TH2F("hTrigMulSel","Multiplicity after selection vs. Trigger types",24,-1.5,22.5,1000,0.,10000.);
    fHisTrigMulSel->GetXaxis()->SetBinLabel(1,"All");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(2,"kAny");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(3,"kMB");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(4,"kINT7");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(5,"kINT8");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(6,"kCINT5");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(7,"kCent");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(8,"kSemiCent");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(9,"kEMC1");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(10,"kEMC7");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(11,"kEMC8");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(12,"kEMCJET7");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(13,"kEMCGAMMA7");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(14,"kEMCJET8");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(15,"kEMCGAMMA8");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(16,"Muons");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(17,"PHOS");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(18,"TRD");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(19,"TRDHJT");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(20,"TRDHSE");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(21,"HighMult");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(22,"SPI7");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(23,"SPI8");
    fHisTrigMulSel->GetXaxis()->SetBinLabel(24,"Others");

    AliCounterCollection *trigCounter=new AliCounterCollection("trigCounter");
    trigCounter->AddRubric("run",500000);
    trigCounter->AddRubric("triggerType","All/Any/MB/Cent/SemiCent/EMCAL/MUON/NoPhysSelMUON/NoPhysSelEvNot7/NoPhysSelCMUP1/NoPhysSelMB/NoPhysSelCent/NoPhysSelSemiCent/CINT7/INT8");
    trigCounter->Init();

    AliCounterCollection *trigCounter2=new AliCounterCollection("trigCounter2");
    trigCounter2->AddRubric("run",500000);
    trigCounter2->AddRubric("triggerType","All/Any/MB/CINT7/INT8/NoPhysSelEvNot7/NoPhysSelMB/HighMult/SPI7/SPI8/EMC1/EMC7/EMC8/EMCJET7/EMCJET8/EMCGAMMA/TRD/TRDHJT/TRDHSE");
    trigCounter2->Init();

    fHisWhyEvRejected=new TH1F("hWhyEvRejected", "Why Event rejected",7,-1.5,5.5);

    fHisWhyEvRejected->GetXaxis()->SetBinLabel(1,"N events");
    fHisWhyEvRejected->GetXaxis()->SetBinLabel(2,"pileup");
    fHisWhyEvRejected->GetXaxis()->SetBinLabel(3,"centrality");
    fHisWhyEvRejected->GetXaxis()->SetBinLabel(4,"Vertex not found");
    fHisWhyEvRejected->GetXaxis()->SetBinLabel(5,"trigger");
    fHisWhyEvRejected->GetXaxis()->SetBinLabel(6,"z vertex out of 10 cm");
    fHisWhyEvRejected->GetXaxis()->SetBinLabel(7,"physics sel");


    fOutputEvSelection->Add(evselection);
    fOutputEvSelection->Add(fHisxvtx);
    fOutputEvSelection->Add(fHisyvtx);
    fOutputEvSelection->Add(fHiszvtx);
    fOutputEvSelection->Add(fHisxvtxSelEv);
    fOutputEvSelection->Add(fHisyvtxSelEv);
    fOutputEvSelection->Add(fHiszvtxSelEv);
    fOutputEvSelection->Add(fHisWhichVert);
    fOutputEvSelection->Add(fHisWhichVertSelEv);
    fOutputEvSelection->Add(fHisTrigCent);
    fOutputEvSelection->Add(fHisTrigMul);
    fOutputEvSelection->Add(fHisTrigMulSel);
    fOutputEvSelection->Add(fHisTrigCentSel);
    fOutputEvSelection->Add(trigCounter);
    fOutputEvSelection->Add(trigCounter2);
    fOutputEvSelection->Add(fHisWhyEvRejected);
    fOutputEvSelection->Add(fHisnClsITSvsNtrackletsSel);
    fOutputEvSelection->Add(fHiszvtxvsSPDzvtx);
    fOutputEvSelection->Add(fHiszvtxvsSPDzvtxSel);

  }
  if(fOnOff[4]){ // FLOW OBSERVABLES
    fOutputFlowObs=new TList();
    fOutputFlowObs->SetOwner();
    fOutputFlowObs->SetName(GetOutputSlot(8)->GetContainer()->GetName());

    fFlowEvent = new AliFlowEvent(3000);
    fRFPcuts = new AliFlowTrackCuts("rfpCuts");

    fHisFEvents = new TH2F("hFlowEvents","FlowEvent Selection",7,0,7,7,-10,60);
    fHisFEvents->GetXaxis()->SetBinLabel(1,"REACHED");
    fHisFEvents->GetXaxis()->SetBinLabel(2,"TRIGGERED");
    fHisFEvents->GetXaxis()->SetBinLabel(3,"kMB");
    fHisFEvents->GetXaxis()->SetBinLabel(4,"kCent");
    fHisFEvents->GetXaxis()->SetBinLabel(5,"kSemiC");
    fHisFEvents->GetXaxis()->SetBinLabel(6,"Triggered + vtx cut");
    fHisFEvents->GetXaxis()->SetBinLabel(7,"UnexpectedBehaviour");
    fOutputFlowObs->Add(fHisFEvents);

    TString ref[3] = {"FB1","FB128","VZE"};
    Int_t etabin[3] = {40,40,20};
    Int_t etamax[3] = { 1, 1, 5};
    for(Int_t i=0; i<3; ++i) {
      fHisQ[i]= new TProfile2D( Form("h%s_Q",ref[i].Data()),
			     Form("Q_{2} components for %s",ref[i].Data()),
			     4,0,4,12,0,60,"s");
      fHisQ[i]->GetXaxis()->SetBinLabel(1,"Qx^{-}");
      fHisQ[i]->GetXaxis()->SetBinLabel(2,"Qy^{-}");
      fHisQ[i]->GetXaxis()->SetBinLabel(3,"Qx^{+}");
      fHisQ[i]->GetXaxis()->SetBinLabel(4,"Qy^{+}");
      fHisQ[i]->GetYaxis()->SetTitle("Centrality");
      fOutputFlowObs->Add(fHisQ[i]);

      fHisAngleQ[i] = new TH2F( Form("h%s_AngleQ",ref[i].Data()),
			     Form("#Psi_{2} for %s",ref[i].Data()),
			     72,0,TMath::Pi(),12,0,60);
      fHisAngleQ[i]->GetXaxis()->SetTitle( Form("#Psi_{2}^{%s}",ref[i].Data()) );
      fHisAngleQ[i]->GetYaxis()->SetTitle("Centrality");
      fOutputFlowObs->Add(fHisAngleQ[i]);

      fHisPhiEta[i] = new TH3F( Form("h%s_PhiEta",ref[i].Data()),
			     Form("Eta vs Phi for %s",ref[i].Data()),
			     144,0,TMath::TwoPi(),etabin[i],-1.0*etamax[i],+1.0*etamax[i],12,0,60);
      fHisPhiEta[i]->GetXaxis()->SetTitle("Phi");
      fHisPhiEta[i]->GetYaxis()->SetTitle("Eta");
      fHisPhiEta[i]->GetZaxis()->SetTitle("Centrality");
      fOutputFlowObs->Add(fHisPhiEta[i]);

    }
    fHisTPCVZE_AngleQ = new TH3F("hTPCVZE_AngleQ","#Psi_{2}^{VZERO} vs #Psi_{2}^{TPC}",   72,0,TMath::Pi(),72,0,TMath::Pi(),12,0,60);
    fHisTPCVZE_AngleQ->GetXaxis()->SetTitle("#Psi_{2}^{TPC}");
    fHisTPCVZE_AngleQ->GetYaxis()->SetTitle("#Psi_{2}^{VZE}");
    fHisTPCVZE_AngleQ->GetZaxis()->SetTitle("Centrality");
    fOutputFlowObs->Add(fHisTPCVZE_AngleQ);

    fHisCentVsMultRPS = new TH2F("hCentVsMultRPS", " Centrality Vs. Multiplicity RPs",5000, 0, 5000.,12,0,60 );
    fHisCentVsMultRPS->GetXaxis()->SetTitle("Multiplicity RPs");
    fHisCentVsMultRPS->GetYaxis()->SetTitle("Centrality");
    fOutputFlowObs->Add(fHisCentVsMultRPS);
  }

  // Post the data
  PostData(1,fOutputEntries);

  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) PostData(5,fOutputCounters);
  if(fOnOff[3]) PostData(7,fOutputEvSelection);
  if(fOnOff[4]) PostData(8,fOutputFlowObs);

  if(!fOnOff[0] && !fOnOff[1] && !fOnOff[2]) AliError("Nothing will be filled!");
}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(fDebug>2) printf("Analysing decay %d\n",fDecayChannel);
  // Post the data already here
  PostData(1,fOutputEntries);
  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) {
    PostData(5,fOutputCounters);
    if(fCuts->GetUseCentrality()) PostData(6,fOutputCheckCentrality);
  }


  fHisNentries->Fill(0);
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel==-1) fHisNentries->Fill(2);
    if (matchingAODdeltaAODlevel==0)  fHisNentries->Fill(3);
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return;
    }
    fHisNentries->Fill(1);
  }

  if( fCuts->IsEventRejectedDueToMismatchOldNewCentrality() ) fHisNentries->Fill(4);


  TClonesArray *arrayProng =0;

  // Load all the branches of the DeltaAOD - needed for SelectionBit counting
  TClonesArray *arrayProngD0toKpi  =0;
  TClonesArray *arrayProng3Prong   =0;
  TClonesArray *arrayProng4Prong   =0;
  TClonesArray *arrayProngDstar    =0;
  TClonesArray *arrayProngCascades =0;

  Int_t pdg=0;
  Int_t *pdgdaughters=0x0;
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

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
      AliAODEvent *aodFromExt = ext->GetAOD();

      arrayProng3Prong  =(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      arrayProng4Prong  =(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
      arrayProngD0toKpi =(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      arrayProngDstar   =(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
      arrayProngCascades=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");

    switch(fDecayChannel){

    case 0:
      arrayProng=arrayProng3Prong;
      pdg=411;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=211;//pi
	  pdgdaughters[1]=321;//K
	  pdgdaughters[2]=211;//pi
	}
	break;
    case 1:
      arrayProng=arrayProngD0toKpi;
      pdg=421;
	if(fReadMC){
	  pdgdaughters =new Int_t[2];
	  pdgdaughters[0]=211;//pi
	  pdgdaughters[1]=321;//K
	}
	break;
    case 2:
      arrayProng=arrayProngDstar;
      pdg=413;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[1]=211;//pi
	  pdgdaughters[0]=321;//K
	  pdgdaughters[2]=211;//pi (soft?)
	}
	break;
    case 3:
      arrayProng=arrayProng3Prong;
      pdg=431;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=321;//K
	  pdgdaughters[1]=321;//K
	  pdgdaughters[2]=211;//pi
	}
	break;
    case 4:
      arrayProng=arrayProng4Prong;
      pdg=421;
	if(fReadMC){
	  pdgdaughters =new Int_t[4];
	  pdgdaughters[0]=321;
	  pdgdaughters[1]=211;
	  pdgdaughters[2]=211;
	  pdgdaughters[3]=211;
	}
	break;
    case 5:
      arrayProng=arrayProng3Prong;
      pdg=4122;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=2212;//p
	  pdgdaughters[1]=321;//K
	  pdgdaughters[2]=211;//pi
	}
	break;
    case kLambdactoV0:
      arrayProng=arrayProngCascades;
	pdg=4122;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=2212;//p
	  pdgdaughters[1]=211;//pi
	  pdgdaughters[2]=211;//pi
	}
	break;
      }
    }
  } else if(aod) {

    arrayProng3Prong  =(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayProng4Prong  =(TClonesArray*)aod->GetList()->FindObject("Charm4Prong");
    arrayProngD0toKpi =(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    arrayProngDstar   =(TClonesArray*)aod->GetList()->FindObject("Dstar");
    arrayProngCascades=(TClonesArray*)aod->GetList()->FindObject("CascadesHF");

    switch(fDecayChannel){

    case 0:
      arrayProng=arrayProng3Prong;
      pdg=411;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=211;//pi
	pdgdaughters[1]=321;//K
	pdgdaughters[2]=211;//pi
      }
      break;
    case 1:
      arrayProng=arrayProngD0toKpi;
      pdg=421;
      if(fReadMC){
	pdgdaughters =new Int_t[2];
	pdgdaughters[0]=211;//pi
	pdgdaughters[1]=321;//K
      }
      break;
    case 2:
      arrayProng=arrayProngDstar;
      pdg=413;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[1]=211;//pi
	pdgdaughters[0]=321;//K
	pdgdaughters[2]=211;//pi (soft?)
      }
      break;
    case 3:
      arrayProng=arrayProng3Prong;
      pdg=431;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=321;//K
	pdgdaughters[1]=321;//K
	pdgdaughters[2]=211;//pi
      }
      break;
    case 4:
      arrayProng=arrayProng4Prong;
      pdg=421;
      if(fReadMC){
	pdgdaughters =new Int_t[4];
	pdgdaughters[0]=321;
	pdgdaughters[1]=211;
	pdgdaughters[2]=211;
	pdgdaughters[3]=211;
      }
      break;
    case 5:
      arrayProng=arrayProng3Prong;
      pdg=4122;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=2212;//p
	pdgdaughters[1]=321;//K
	pdgdaughters[2]=211;//pi
      }
      break;
    case kLambdactoV0:
      arrayProng=arrayProngCascades;
      pdg=4122;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=2212;//p
	pdgdaughters[1]=211;//pi
	pdgdaughters[2]=211;//pi
      }
      break;
    }
  }
  Bool_t isSimpleMode=fSimpleMode;
  if(!arrayProng) {
    AliInfo("Branch not found! The output will contain only track related histograms\n");
    isSimpleMode=kTRUE;
    fHisNentries->Fill(7);
  }

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(!aod) {
    delete [] pdgdaughters;
    return;
  }

  //check if MC
  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSEHFQA::UserExec: MC particles branch not found!\n");
      delete [] pdgdaughters;
      return;
    }

    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEHFQA::UserExec: MC header branch not found!\n");
      delete [] pdgdaughters;
      return;
    }
  }


  UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Double_t centrality=fCuts->GetCentrality(aod);
  AliAODHeader * header = dynamic_cast<AliAODHeader*>(aod->GetHeader());
  if(!header) AliFatal("Not a standard AOD");

  Double_t multiplicity=header->GetRefMultiplicity();
  Int_t runNumber = aod->GetRunNumber();
  TString trigClass=aod->GetFiredTriggerClasses();
  Int_t nAODtracks=aod->GetNumberOfTracks();
  Int_t nSelTracksTPCOnly=0;
  Int_t nSelTracksTPCITS=0;
  Int_t nSelTracksTPCITS1SPD=0;
  Int_t ntracksFBit4=0;

  AliTRDTriggerAnalysis trdSelection;
  trdSelection.CalcTriggers(aod);

  if(fReadMC) {
    if(aod->GetTriggerMask()==0 &&
       (runNumber>=195344 && runNumber<=195677)){
      AliDebug(3,"Event rejected because of null trigger mask");
      delete [] pdgdaughters;
      return;
    }
  }

  for (Int_t k=0;k<nAODtracks;k++){
    AliAODTrack* track=dynamic_cast<AliAODTrack*>(aod->GetTrack(k));
    if(!track) AliFatal("Not a standard AOD");
    if(track->GetID()<0) continue;
    Int_t nclsTot=0,nclsSPD=0;
    for(Int_t l=0;l<6;l++) {
      if(TESTBIT(track->GetITSClusterMap(),l)) {
	nclsTot++; if(l<2) nclsSPD++;
      }
    }
    UShort_t nTPCClus=track->GetTPCClusterMap().CountBits();
    if(TMath::Abs(track->Eta())<0.8 && nTPCClus>=70 && track->GetStatus()&AliESDtrack::kTPCrefit){
      if(track->TestFilterBit(1))  nSelTracksTPCOnly++;
      if(track->GetStatus()&AliESDtrack::kITSrefit){
	nSelTracksTPCITS++;
	if(nclsSPD>0) nSelTracksTPCITS1SPD++;
      }
    }
    if(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){
      ntracksFBit4++;
    }
  }

  if(fOnOff[4]) {
    FillFlowObs(aod);
    PostData(8,fOutputFlowObs);
  }
  if(fOnOff[3]){
    AliCounterCollection* trigCount=(AliCounterCollection*)fOutputEvSelection->FindObject("trigCounter");
    AliCounterCollection* trigCount2=(AliCounterCollection*)fOutputEvSelection->FindObject("trigCounter2");

    fHisTrigCent->Fill(-1.,centrality);
    fHisTrigMul->Fill(-1.,multiplicity);
    trigCount->Count(Form("triggerType:All/Run:%d",runNumber));
    trigCount2->Count(Form("triggerType:All/Run:%d",runNumber));
    if(evSelMask==0){
      if(aod->GetEventType()!=7 || trigClass.Contains("BEAMB")){
	trigCount->Count(Form("triggerType:NoPhysSelEvNot7/Run:%d",runNumber));
	trigCount2->Count(Form("triggerType:NoPhysSelEvNot7/Run:%d",runNumber));
      }else if(trigClass.Contains("CMUP1")){
	trigCount->Count(Form("triggerType:NoPhysSelCMUP1/Run:%d",runNumber));
      }else if(trigClass.Contains("MUON")){
	trigCount->Count(Form("triggerType:NoPhysSelMUON/Run:%d",runNumber));
      }else if(trigClass.Contains("CPBI2_B1-B") || trigClass.Contains(" CPBI2WU_B1-B")){
	trigCount->Count(Form("triggerType:NoPhysSelMB/Run:%d",runNumber));
	trigCount2->Count(Form("triggerType:NoPhysSelMB/Run:%d",runNumber));
      }else if(trigClass.Contains("CCENT") || trigClass.Contains("CVHN")){
	trigCount->Count(Form("triggerType:NoPhysSelCent/Run:%d",runNumber));
      }else if(trigClass.Contains("CSEMI") || trigClass.Contains("CVLN")){
	trigCount->Count(Form("triggerType:NoPhysSelSemiCent/Run:%d",runNumber));
      }
    }
    if(evSelMask & AliVEvent::kAny){
      fHisTrigCent->Fill(0.,centrality);
      fHisTrigMul->Fill(0.,multiplicity);
      trigCount->Count(Form("triggerType:Any/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:Any/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kMB){
      fHisTrigCent->Fill(1.,centrality);
      fHisTrigMul->Fill(1.,multiplicity);
      trigCount->Count(Form("triggerType:MB/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:MB/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kINT7){
      fHisTrigCent->Fill(2.,centrality);
      fHisTrigMul->Fill(2.,multiplicity);
      trigCount->Count(Form("triggerType:CINT7/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:CINT7/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kINT8){
      fHisTrigCent->Fill(3.,centrality);
      fHisTrigMul->Fill(3.,multiplicity);
      trigCount->Count(Form("triggerType:INT8/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:INT8/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kCINT5){
      fHisTrigCent->Fill(4.,centrality);
      fHisTrigMul->Fill(4.,multiplicity);
    }
    if(evSelMask & AliVEvent::kCentral){
      fHisTrigCent->Fill(5.,centrality);
      fHisTrigMul->Fill(5.,multiplicity);
      trigCount->Count(Form("triggerType:Cent/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kSemiCentral){
      fHisTrigCent->Fill(6.,centrality);
      fHisTrigMul->Fill(6.,multiplicity);
      trigCount->Count(Form("triggerType:SemiCent/Run:%d",runNumber));
    }

    if(evSelMask & AliVEvent::kEMC1){
      fHisTrigCent->Fill(7.,centrality);
      fHisTrigMul->Fill(7.,multiplicity);
      trigCount->Count(Form("triggerType:EMCAL/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:EMC1/Run:%d",runNumber));
    }
    if((evSelMask & AliVEvent::kEMC7) && trigClass.Contains("CEMC7")){
      fHisTrigCent->Fill(8.,centrality);
      fHisTrigMul->Fill(8.,multiplicity);
      trigCount->Count(Form("triggerType:EMCAL/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:EMC7/Run:%d",runNumber));
    }
    if((evSelMask & AliVEvent::kEMC8) && trigClass.Contains("CEMC8")){
      fHisTrigCent->Fill(9.,centrality);
      fHisTrigMul->Fill(9.,multiplicity);
      trigCount->Count(Form("triggerType:EMCAL/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:EMC8/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kEMCEJE){
       trigCount->Count(Form("triggerType:EMCAL/Run:%d",runNumber));
       if(trigClass.Contains("CEMC7EJE")) {
	 trigCount2->Count(Form("triggerType:EMCJET7/Run:%d",runNumber));
	 fHisTrigCent->Fill(10.,centrality);
	 fHisTrigMul->Fill(10.,multiplicity);

       }
       else if(trigClass.Contains("CEMC8EJE")) {
	 trigCount2->Count(Form("triggerType:EMCJET8/Run:%d",runNumber));
	 fHisTrigCent->Fill(12.,centrality);
	 fHisTrigMul->Fill(12.,multiplicity);
       }
    }
    if(evSelMask & AliVEvent::kEMCEGA){
      if(trigClass.Contains("CEMC7EGA")) {
	fHisTrigCent->Fill(11.,centrality);
	fHisTrigMul->Fill(11.,multiplicity);
      } else if (trigClass.Contains("CEMC8EGA")){
	fHisTrigCent->Fill(13.,centrality);
	fHisTrigMul->Fill(13.,multiplicity);

      }
      trigCount->Count(Form("triggerType:EMCAL/Run:%d",runNumber));
      trigCount2->Count(Form("triggerType:EMCGAMMA/Run:%d",runNumber));
    }
    if(evSelMask & (((AliVEvent::kCMUS5 | AliVEvent::kMUSH7) | (AliVEvent::kMUL7 | AliVEvent::kMUU7)) |  (AliVEvent::kMUS7 | AliVEvent::kMUON))){
      fHisTrigCent->Fill(14.,centrality);
      fHisTrigMul->Fill(14.,multiplicity);
      trigCount->Count(Form("triggerType:MUON/Run:%d",runNumber));
    }
    if(evSelMask & (AliVEvent::kPHI1 | AliVEvent::kPHI7)){
      fHisTrigCent->Fill(15.,centrality);
      fHisTrigMul->Fill(15.,multiplicity);
    }
    if(evSelMask & (AliVEvent::kTRD)){
      fHisTrigCent->Fill(16.,centrality);
      fHisTrigMul->Fill(16.,multiplicity);
      trigCount2->Count(Form("triggerType:TRD/Run:%d",runNumber));
    }
    if((evSelMask & AliVEvent::kTRD) && trdSelection.IsFired(AliTRDTriggerAnalysis::kHJT)){
      fHisTrigCent->Fill(17.,centrality);
      fHisTrigMul->Fill(17.,multiplicity);
      trigCount2->Count(Form("triggerType:TRDHJT/Run:%d",runNumber));
    }
    if((evSelMask & AliVEvent::kTRD) && trdSelection.IsFired(AliTRDTriggerAnalysis::kHSE)){
      fHisTrigCent->Fill(18.,centrality);
      fHisTrigMul->Fill(18.,multiplicity);
      trigCount2->Count(Form("triggerType:TRDHSE/Run:%d",runNumber));
    }
    if(evSelMask & (AliVEvent::kHighMult)){
      fHisTrigCent->Fill(19.,centrality);
      fHisTrigMul->Fill(19.,multiplicity);
      trigCount2->Count(Form("triggerType:HighMult/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kSPI7){
      if(trigClass.Contains("CSPI7")) {
	fHisTrigCent->Fill(20.,centrality);
	fHisTrigMul->Fill(20.,multiplicity);
	trigCount2->Count(Form("triggerType:SPI7/Run:%d",runNumber));
      }
    }
    if(evSelMask & AliVEvent::kSPI){
      if(trigClass.Contains("CSPI8")) {
	fHisTrigCent->Fill(21.,centrality);
	fHisTrigMul->Fill(21.,multiplicity);
        trigCount2->Count(Form("triggerType:SPI8/Run:%d",runNumber));
      }
    }
    if(evSelMask & (AliVEvent::kDG5 | AliVEvent::kZED)){
      fHisTrigCent->Fill(22.,centrality);
      fHisTrigMul->Fill(22.,multiplicity);
    }
  }


  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  Double_t magField  = aod->GetMagneticField();
  if(!aod->GetPrimaryVertex() || TMath::Abs(magField)<0.001) {
    delete [] pdgdaughters;
    return;
  }

  // count event
  fHisNentries->Fill(5);
  //count events with good vertex
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  Double_t pos[3],cov[6];
  vtx1->GetXYZ(pos);
  vtx1->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) fHisNentries->Fill(9);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  //TString trigclass=aod->GetFiredTriggerClasses();
  //if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fHisNentries->Fill(9); //tmp




  Bool_t evSelbyCentrality=kTRUE,evSelected=kTRUE,evSelByVertex=kTRUE,evselByPileup=kTRUE,evSelByPS=kTRUE;

  if(fOnOff[3]){
     fHisWhyEvRejected->Fill(-1);
  }

  //select event
  if(!fCuts->IsEventSelected(aod)) {
    evSelected=kFALSE;
    if(fCuts->IsEventRejectedDueToPileup()) {
      fHisWhyEvRejected->Fill(0);
      evselByPileup=kFALSE;
    }// rejected for pileup
    if(fCuts->IsEventRejectedDueToCentrality()) {
      fHisWhyEvRejected->Fill(1);
      evSelbyCentrality=kFALSE; //rejected by centrality
    }
    if(fCuts->IsEventRejectedDueToNotRecoVertex() ||
       fCuts->IsEventRejectedDueToVertexContributors()){
      evSelByVertex=kFALSE;
      fHisWhyEvRejected->Fill(2);
    }
    if(fCuts->IsEventRejectedDueToTrigger()){
      fHisWhyEvRejected->Fill(3);
    }
    if(fCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) {
      evSelByVertex=kFALSE;
      if(fOnOff[3]) ((AliCounterCollection*)fOutputEvSelection->FindObject("evselection"))->Count(Form("evnonsel:zvtx/Run:%d",runNumber));
      fHisWhyEvRejected->Fill(4);
    }
    if(fCuts->IsEventRejectedDuePhysicsSelection()) {
      evSelByPS=kFALSE;
      fHisWhyEvRejected->Fill(5);
    }
  }
  if(evSelected && fOnOff[3]){
    fHisTrigCentSel->Fill(-1.,centrality);
    fHisTrigMulSel->Fill(-1.,multiplicity);
    if(evSelMask & AliVEvent::kAny) {
      fHisTrigCentSel->Fill(0.,centrality);
      fHisTrigMulSel->Fill(0.,multiplicity);}
    if(evSelMask & AliVEvent::kMB) {
      fHisTrigCentSel->Fill(1.,centrality);
      fHisTrigMulSel->Fill(1.,multiplicity);}
    if(evSelMask & AliVEvent::kINT7){
      fHisTrigCentSel->Fill(2.,centrality);
      fHisTrigMulSel->Fill(2.,multiplicity);}
    if(evSelMask & AliVEvent::kINT8){
      fHisTrigCentSel->Fill(3.,centrality);
      fHisTrigMulSel->Fill(3.,multiplicity);}
    if(evSelMask & AliVEvent::kCINT5){
      fHisTrigCentSel->Fill(4.,centrality);
      fHisTrigMulSel->Fill(4.,multiplicity);}
    if(evSelMask & AliVEvent::kCentral){
      fHisTrigCentSel->Fill(5.,centrality);
      fHisTrigMulSel->Fill(5.,multiplicity);}
    if(evSelMask & AliVEvent::kSemiCentral){
      fHisTrigCentSel->Fill(6.,centrality);
      fHisTrigMulSel->Fill(6.,multiplicity);}
    if(evSelMask & AliVEvent::kEMC1){
      fHisTrigCentSel->Fill(7.,centrality);
      fHisTrigMulSel->Fill(7.,multiplicity);
    }
    if((evSelMask & AliVEvent::kEMC7) && trigClass.Contains("CEMC7")){
      fHisTrigCentSel->Fill(8.,centrality);
      fHisTrigMulSel->Fill(8.,multiplicity);
    }
    if((evSelMask & AliVEvent::kEMC8) && trigClass.Contains("CEMC8")){
      fHisTrigCentSel->Fill(9.,centrality);
      fHisTrigMulSel->Fill(9.,multiplicity);
    }
      if((evSelMask & AliVEvent::kEMCEJE) && trigClass.Contains("CEMC7EJE")){
      fHisTrigCentSel->Fill(10.,centrality);
      fHisTrigMulSel->Fill(10.,multiplicity);
    }
    if((evSelMask & AliVEvent::kEMCEGA) && trigClass.Contains("CEMC7EGA")){
      fHisTrigCentSel->Fill(11.,centrality);
      fHisTrigMulSel->Fill(11.,multiplicity);
    }
    if((evSelMask & AliVEvent::kEMCEJE) && trigClass.Contains("CEMC8EJE")){
      fHisTrigCentSel->Fill(12.,centrality);
      fHisTrigMulSel->Fill(12.,multiplicity);
    }
    if((evSelMask & AliVEvent::kEMCEGA) && trigClass.Contains("CEMC8EGA")){
      fHisTrigCentSel->Fill(13.,centrality);
      fHisTrigMulSel->Fill(13.,multiplicity);
    }
    if(evSelMask & (((AliVEvent::kCMUS5 | AliVEvent::kMUSH7) | (AliVEvent::kMUL7 | AliVEvent::kMUU7)) |  (AliVEvent::kMUS7 | AliVEvent::kMUON))){
      fHisTrigCentSel->Fill(14.,centrality);
      fHisTrigMulSel->Fill(14.,multiplicity);}
    if(evSelMask & (AliVEvent::kPHI1 | AliVEvent::kPHI7)){
      fHisTrigCentSel->Fill(15.,centrality);
      fHisTrigMulSel->Fill(15.,multiplicity);}
    if(evSelMask & (AliVEvent::kTRD)){
      fHisTrigCentSel->Fill(16.,centrality);
      fHisTrigMulSel->Fill(16.,multiplicity);
    }
    if((evSelMask & AliVEvent::kTRD) && trdSelection.IsFired(AliTRDTriggerAnalysis::kHJT)){
      fHisTrigCentSel->Fill(17.,centrality);
      fHisTrigMulSel->Fill(17.,multiplicity);
    }
    if((evSelMask & AliVEvent::kTRD) && trdSelection.IsFired(AliTRDTriggerAnalysis::kHSE)){
      fHisTrigCentSel->Fill(18.,centrality);
      fHisTrigMulSel->Fill(18.,multiplicity);
    }
    if(evSelMask & (AliVEvent::kHighMult)){
      fHisTrigCentSel->Fill(19.,centrality);
      fHisTrigMulSel->Fill(19.,multiplicity);}
    if(evSelMask & AliVEvent::kSPI7){
      if(trigClass.Contains("CSPI7")) {
	fHisTrigCentSel->Fill(20.,centrality);
	fHisTrigMulSel->Fill(20.,multiplicity);
      }
    }
    if(evSelMask & AliVEvent::kSPI){
      if(trigClass.Contains("CSPI8")) {
	fHisTrigCentSel->Fill(21.,centrality);
	fHisTrigMulSel->Fill(21.,multiplicity);
      }
    }
    if(evSelMask & (AliVEvent::kDG5 | AliVEvent::kZED)){
      fHisTrigCentSel->Fill(22.,centrality);
      fHisTrigMulSel->Fill(22.,multiplicity);}
  }

  if(evSelected || (!evSelbyCentrality && evSelByVertex && evselByPileup && evSelByPS)){ //events selected or not selected because of centrality
    if(fOnOff[2] && fCuts->GetUseCentrality()){

      Float_t stdCentf=fCuts->GetCentrality(aod);
      Int_t stdCent = (Int_t)(stdCentf+0.5);
      Float_t secondCentf =fCuts->GetCentrality(aod,fEstimator);
      Int_t secondCent = (Int_t)(secondCentf+0.5);
      Int_t mincent=stdCent-stdCent%10;
      Float_t stdSignal = 0.;
      Float_t secondSignal = 0.;
      AliAODVZERO *vzeroAOD = (AliAODVZERO*)aod->GetVZEROData();
      AliAODZDC *zdcAOD = (AliAODZDC*)aod->GetZDCData();
      const Double_t *towerZNASignal = zdcAOD->GetZNATowerEnergy();
      switch(fCuts->GetUseCentrality())
      {
         case AliRDHFCuts::kCentV0M:
            stdSignal = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();
            break;
         case AliRDHFCuts::kCentV0A:
            stdSignal = vzeroAOD->GetMTotV0A();
            break;
         case AliRDHFCuts::kCentZNA:
            stdSignal = towerZNASignal[0];
            break;
         default:
            stdSignal = 0.;
            break;
      }
      switch(fEstimator)
      {
         case AliRDHFCuts::kCentV0M:
            secondSignal = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();
            break;
         case AliRDHFCuts::kCentV0A:
            secondSignal = vzeroAOD->GetMTotV0A();
            break;
         case AliRDHFCuts::kCentZNA:
            secondSignal = towerZNASignal[0];
            break;
         default:
            secondSignal = 0.;
            break;
      }
      //AliCentrality *aodcent = aod->GetCentrality();
      // Float_t spdCentf = aodcent->GetCentralityPercentile("CL1");
      if(stdCentf==-1) {
         mincent=-10;
         stdCent=-1;
      }
      if(mincent==100)mincent--;
      ((AliCounterCollection*)fOutputCounters->FindObject("stdEstimator"))->Count(Form("centralityclass:%d_%d/Run:%d",mincent,mincent+10,runNumber));

      mincent=secondCent-secondCent%10;
      if(secondCentf==-1) {
	mincent=-10;
	secondCent=-1;
      }
      if(mincent==100)mincent--;
      ((AliCounterCollection*)fOutputCounters->FindObject("secondEstimator"))->Count(Form("centralityclass:%d_%d/Run:%d",mincent,mincent+10,runNumber));

      if(stdCent<fCuts->GetMinCentrality() || stdCent>fCuts->GetMaxCentrality()){
	fHisNtrackletsOut->Fill(aod->GetTracklets()->GetNumberOfTracklets());
	fHisMultOut->Fill(header->GetRefMultiplicity());
      }else{
	fHisNtrackletsIn->Fill(aod->GetTracklets()->GetNumberOfTracklets());
	fHisMultIn->Fill(header->GetRefMultiplicity());
      }
      fHisMultvsPercentile->Fill(header->GetRefMultiplicity(),stdCentf);
      fHisntrklvsPercentile->Fill(aod->GetTracklets()->GetNumberOfTracklets(),stdCentf);
      fHisntrklvsPercentile01->Fill(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.),stdCentf);
      fHisnTPCTracksvsPercentile->Fill(nSelTracksTPCOnly,stdCentf);
      fHisnTPCITSTracksvsPercentile->Fill(nSelTracksTPCITS,stdCentf);
      fHisnTPCITS1SPDTracksvsPercentile->Fill(nSelTracksTPCITS1SPD,stdCentf);
      fHisStdEstimSignalPercentile->Fill(stdSignal,stdCentf);
      fHisStdEstimSignal->Fill(stdSignal);
      fHisStdEstimSignalNtrackletsIn->Fill(stdSignal,aod->GetTracklets()->GetNumberOfTracklets());
      fHisStdPercentileSecondPercentile->Fill(stdCentf,secondCentf);
      fHisStdSignalSecondSignal->Fill(stdSignal,secondSignal);

      if( (fCuts->GetMultSelectionObjectName()).CompareTo("MultSelection")!=0 ){
          AliMultSelection *multSelectionOld = (AliMultSelection*) aod -> FindListObject("MultSelection");
          Float_t oldFrmwCent = multSelectionOld ? multSelectionOld->GetMultiplicityPercentile("V0M") : -1;
	  //          Printf(" HFQA task: new centrality %3.3f vs old centrality %3.2f \n",stdCentf,oldFrmwCent);
          fHisStdPercentileOldFrwPercentile->Fill(stdCentf,oldFrmwCent);
          //
          Double_t differenceCent = TMath::Abs( stdCentf - oldFrmwCent );
          if(differenceCent>1) fHisStdPercentileOldFrwPercentileDev->Fill(0.);
          if(differenceCent>0.5) fHisStdPercentileOldFrwPercentileDev->Fill(1.);
          if( (stdCentf<20) && (oldFrmwCent>20) ) fHisStdPercentileOldFrwPercentileDev->Fill(2.);
          if( (stdCentf>20) && (oldFrmwCent<20) ) fHisStdPercentileOldFrwPercentileDev->Fill(3.);
          if( fCuts->IsEventRejectedDueToMismatchOldNewCentrality() ) fHisStdPercentileOldFrwPercentileDev->Fill(4.);
      }

      PostData(6,fOutputCheckCentrality);

    } else{
      if(fOnOff[0]){
	fHisNtracklets->Fill(aod->GetTracklets()->GetNumberOfTracklets());
	fHisNtracklets01->Fill(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.));
	fHisMult->Fill(header->GetRefMultiplicity());
	fHisMultFBit4->Fill(ntracksFBit4);
	fHisMultComb05->Fill(header->GetRefMultiplicityComb05());
	fHisMultComb08->Fill(header->GetRefMultiplicityComb08());
      }
    }
  }

  if(evSelected || (!evSelbyCentrality && evSelByVertex && evselByPileup && evSelByPS) || (!evSelByVertex && evselByPileup && evSelByPS)){ //events selected or not selected because of centrality
    if(fOnOff[2] && fCuts->GetUseCentrality()){
      fHisntrklvsPercentile01AllEv->Fill(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.),fCuts->GetCentrality(aod));
    }else{
      if(fOnOff[0]){
  	fHisNtracklets01AllEv->Fill(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.));
      }
    }
  }

  if(fOnOff[3]){
    const AliVVertex *vertex = aod->GetPrimaryVertex();
    Double_t xvtx=vertex->GetX();
    Double_t yvtx=vertex->GetY();
    Double_t zvtx=vertex->GetZ();
    Int_t vtxTyp=0;
    if(vertex->GetNContributors()<=0) vtxTyp=-1;
    TString title=vertex->GetTitle();
    if(title.Contains("Z")) vtxTyp=2;
    if(title.Contains("3D")) vtxTyp=1;
    fHisxvtx->Fill(xvtx);
    fHisyvtx->Fill(yvtx);
    fHiszvtx->Fill(zvtx);
    fHisWhichVert->Fill(vtxTyp);
    const AliVVertex *vSPD = aod->GetPrimaryVertexSPD();
    fHiszvtxvsSPDzvtx->Fill(vSPD->GetZ(),zvtx);
    if(evSelected){
      fHisxvtxSelEv->Fill(xvtx);
      fHisyvtxSelEv->Fill(yvtx);
      fHiszvtxSelEv->Fill(zvtx);
      fHisWhichVertSelEv->Fill(vtxTyp);
      fHiszvtxvsSPDzvtxSel->Fill(vSPD->GetZ(),zvtx);
    }

    Int_t nCls = aod->GetNumberOfITSClusters(0) + aod->GetNumberOfITSClusters(1);
    fHisnClsITSvsNtrackletsSel->Fill(nCls,aod->GetTracklets()->GetNumberOfTracklets());

  }

  if(!evSelected) {
    delete [] pdgdaughters;
    return; //discard all events not selected (vtx and/or centrality)
  }


  AliAODPidHF* pidHF=fCuts->GetPidHF();
  if(!pidHF) {
    delete [] pdgdaughters;
    return;
  }
  //load all the  branches, re-fill the candidates and fill the SelectionBit histo

  Int_t nCand3Prong = arrayProng3Prong->GetEntriesFast();
  Int_t nCandD0toKpi = arrayProngD0toKpi->GetEntriesFast();
  Int_t nCandDstar = arrayProngDstar->GetEntriesFast();
  Int_t nCandCasc=0;
  Int_t n4Prong=0;
  if(arrayProngCascades)nCandCasc = arrayProngCascades->GetEntriesFast();
  if(arrayProng4Prong)n4Prong = arrayProng4Prong->GetEntriesFast();
  // D+, Ds and Lc
  AliAODRecoDecayHF *d;
  for (Int_t iCand = 0; iCand < nCand3Prong; iCand++) {
    d = (AliAODRecoDecayHF*)arrayProng3Prong->UncheckedAt(iCand);
    if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF3Prong*)d))continue;

    if(fUseSelectionBit && !isSimpleMode){
      Double_t ptCand_selBit = d->Pt();
      if(fUseSelectionBit && d->GetSelectionMap()) {
        if(d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) fHisNentriesSelBit->Fill(0.0,ptCand_selBit);
        if(d->HasSelectionBit(AliRDHFCuts::kDsCuts)) fHisNentriesSelBit->Fill(1.0,ptCand_selBit);
        if(d->HasSelectionBit(AliRDHFCuts::kLcCuts)) fHisNentriesSelBit->Fill(2.0,ptCand_selBit);
      }
    }
  }
  // D0kpi
  for (Int_t iCand = 0; iCand < nCandD0toKpi; iCand++) {
    d = (AliAODRecoDecayHF*)arrayProngD0toKpi->UncheckedAt(iCand);
    if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF2Prong*)d))continue;
    if(fUseSelectionBit && !isSimpleMode){
      Double_t ptCand_selBit = d->Pt();
      if(fUseSelectionBit && d->GetSelectionMap()) {
        if(d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) fHisNentriesSelBit->Fill(3.0,ptCand_selBit);
      }
    }
  }
  // Dstar
  for (Int_t iCand = 0; iCand < nCandDstar; iCand++) {
    d = (AliAODRecoDecayHF*)arrayProngDstar->UncheckedAt(iCand);
    if(!vHF->FillRecoCasc(aod,((AliAODRecoCascadeHF*)d),kTRUE))continue;
    if(fUseSelectionBit && !isSimpleMode){
      Double_t ptCand_selBit = d->Pt();
      if(fUseSelectionBit && d->GetSelectionMap()) {
        if(d->HasSelectionBit(AliRDHFCuts::kDstarCuts)) fHisNentriesSelBit->Fill(4.0,ptCand_selBit);
      }
    }
  }
//   //Cascade
  if(arrayProngCascades){
    for (Int_t iCand = 0; iCand < nCandCasc; iCand++) {
      d=(AliAODRecoDecayHF*)arrayProngCascades->UncheckedAt(iCand);
      if(!vHF->FillRecoCasc(aod,((AliAODRecoCascadeHF*)d),kFALSE))continue;
    }
  }
  //end refill

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();

  //AliPIDResponse* respF=pidHF->GetPidResponse();
  AliTPCPIDResponse* tpcres=new AliTPCPIDResponse();
  Bool_t oldPID=pidHF->GetOldPid();
  if(oldPID){
    Double_t alephParameters[5];
    pidHF->GetTPCBetheBlochParams(alephParameters);
    tpcres->SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);
  }


  Int_t ntracks=0;
  Int_t isGoodTrack=0, isFakeTrack=0, isSelTrack=0;

  if(aod) ntracks=aod->GetNumberOfTracks();

  if(fOnOff[0] || fOnOff[1]){
    //loop on tracks in the event
    for (Int_t k=0;k<ntracks;k++){
      AliAODTrack* track=dynamic_cast<AliAODTrack*>(aod->GetTrack(k));
      if(!track) AliFatal("Not a standard AOD");
      // Track selection cuts
      if(track->GetID()<0) continue;
      Double_t d0z0[2],covd0z0[3];
      if(!track->PropagateToDCA(vtx1,magField,99999.,d0z0,covd0z0)) continue;
      if(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){
	fHisd0TracksFilterBit4->Fill(d0z0[0]);
      }
      ULong_t trStatus=track->GetStatus();
      if(trStatus&AliESDtrack::kITSrefit){
	if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)){
	  fHisd0TracksSPDany->Fill(d0z0[0]);
	  if(track->HasPointOnITSLayer(0)){
	    fHisd0TracksSPDin->Fill(d0z0[0]);
	  }
	}
      }

      Bool_t selTrack=kTRUE;
      if (!((trStatus & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
	  !((trStatus & AliVTrack::kITSrefit) == AliVTrack::kITSrefit)){
	selTrack=kFALSE;
      }
      if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){ // BIT(4) standard cuts with very loose DCA
      	 selTrack=kFALSE;
      }
      if(TMath::Abs(track->Eta())>0.9){
      	 selTrack=kFALSE;
      }
      Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
      Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
      if (track->GetTPCNclsF()>0) {
	ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
      }

      Bool_t selTrackNoSPD=selTrack;
      if(selTrack){
	if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)){
	  fHisd0TracksTPCITSSPDany->Fill(d0z0[0]);
	}
      }

      AliAODPid *pid = track->GetDetPid();
      if(!pid && fDebug>1) cout<<"No AliAODPid found"<<endl;

      if(pid && fOnOff[1]){
	Double_t times[AliPID::kSPECIES];
	pid->GetIntegratedTimes(times);

	Double_t tofRes[AliPID::kSPECIES];
	pid->GetTOFpidResolution(tofRes);

	//check TOF
	fHisTOFflags->Fill(0.);
	if (trStatus&AliESDtrack::kTPCout) fHisTOFflags->Fill(1.);
	if (trStatus&AliESDtrack::kTOFout) fHisTOFflags->Fill(2.);
	if (trStatus&AliESDtrack::kTIME) fHisTOFflags->Fill(3.);
	if (trStatus&AliESDtrack::kTOFpid) fHisTOFflags->Fill(4.);
	if (trStatus&AliESDtrack::kTOFmismatch) fHisTOFflags->Fill(5.);

	Bool_t isTOFok=kFALSE;
	if(pidResp){
	  Double_t prob[AliPID::kSPECIES];
	  if(pidResp->ComputeTOFProbability(track,AliPID::kSPECIES,prob)==AliPIDResponse::kDetPidOk){
	    isTOFok=kTRUE;
	    fHisTOFflags->Fill(6.);
	  }
	}

	if(selTrack && isTOFok){
	  Double_t tofTime=pid->GetTOFsignal();
	  AliTOFHeader* tofH=(AliTOFHeader*)aod->GetTOFHeader();
	  if (tofH && (TMath::Abs(tofRes[0]) <= 1.E-16) ) { // new AOD
            // with new AOD we need to retrieve startTime, subtract it and retrieve correctly TOF PID resolutions  *PA*
	    AliTOFPIDResponse tofResp=pidResp->GetTOFResponse();
	    Double_t startTime = tofResp.GetStartTime(track->P());
	    Float_t startTimeRes = tofResp.GetStartTimeRes(track->P());
	    Int_t startTimeMask = tofResp.GetStartTimeMask(track->P());
	    fHisTOFstartTimeDistrib->Fill(startTime);
	    fHisTOFstartTimeMask->Fill(startTimeMask);
	    fHisTOFstartTimeRes->Fill(startTimeRes);
	    tofTime-=startTime;
	    for (Int_t type=0;type<AliPID::kSPECIES;type++) tofRes[type]=tofResp.GetExpectedSigma(track->P(),times[type],AliPID::ParticleMassZ(type));
	  }
	  fHisTOFtime->Fill(times[AliPID::kProton]);
	  fHisTOFtimeKaonHyptime->Fill(track->P(),tofTime-times[3]); //3 is kaon
	  fHisTOFsig->Fill(tofTime);
	  if (pid->GetTOFsignal()< 0) fHisTOFsig->Fill(-1);

	  Double_t nsigma[3]={-10,-10,-10};
	  nsigma[0]=pidResp->NumberOfSigmasTOF(track,AliPID::kPion);
	  nsigma[1]=pidResp->NumberOfSigmasTOF(track,AliPID::kKaon);
	  nsigma[2]=pidResp->NumberOfSigmasTOF(track,AliPID::kProton);

	  fHisTOFsigmaKSigPid->Fill(track->P(),nsigma[1]);
	  fHisTOFsigmaPionSigPid->Fill(track->P(),nsigma[0]);
	  fHisTOFsigmaProtonSigPid->Fill(track->P(),nsigma[2]);
	  if(fReadMC){
	    Int_t label=track->GetLabel();
	    if(label<=0) continue;
	    AliMCParticle* mcpart=(AliMCParticle*)mcArray->At(label);
	    if(mcpart){
	      Int_t abspdgcode=TMath::Abs(mcpart->PdgCode());
	      if(abspdgcode==211) fHisTOFsigmaMCPionSigPid->Fill(track->P(),nsigma[0]);
	      if(abspdgcode==321) fHisTOFsigmaMCKSigPid->Fill(track->P(),nsigma[1]);
	      if(abspdgcode==2212) fHisTOFsigmaMCProtonSigPid->Fill(track->P(),nsigma[2]);

	    }
	  }

	  for (Int_t iS=2; iS<5; iS++){ //we plot TOF Pid resolution for 3-sigma identified particles
	    if ( TMath::Abs(nsigma[iS-2])<3.){
	      switch (iS) {
	      case AliPID::kPion:
		fHisTOFsigPid3sigPion->Fill(tofRes[iS]);
		break;
	      case AliPID::kKaon:
		fHisTOFsigPid3sigKaon->Fill(tofRes[iS]);
		break;
	      case AliPID::kProton:
		fHisTOFsigPid3sigProton->Fill(tofRes[iS]);
		break;
	      default:
		break;
	      }
	    }
	  }
	}//if TOF status
	//}

	if(pidHF && pidHF->CheckStatus(track,"TPC") && selTrack){

	  Double_t TPCp=pid->GetTPCmomentum();
	  Double_t TPCsignal=pid->GetTPCsignal();
	  UShort_t TPCsignalN=pid->GetTPCsignalN();
	  fHisTPCsig->Fill(TPCsignal);
	  fHisTPCsigvsp->Fill(TPCp,TPCsignal);
	  //if (pidHF->IsKaonRaw(track, "TOF"))
	  Double_t nsigma[3]={-10,-10,-10};
	  pidHF->GetnSigmaTPC(track,(Int_t)AliPID::kPion,nsigma[0]);
	  pidHF->GetnSigmaTPC(track,(Int_t)AliPID::kKaon,nsigma[1]);
	  pidHF->GetnSigmaTPC(track,(Int_t)AliPID::kProton,nsigma[2]);

	  fHisTPCsigmaK->Fill(TPCp,nsigma[1]);

	  fHisTPCsigmaPion->Fill(TPCp,nsigma[0]);
	  fHisTPCsigmaProton->Fill(TPCp,nsigma[2]);

	  if(fReadMC){
	    Int_t label=track->GetLabel();
	    if(label<=0) continue;
	    AliMCParticle* mcpart=(AliMCParticle*)mcArray->At(label);
	    if(mcpart){
	      Int_t abspdgcode=TMath::Abs(mcpart->PdgCode());
	      if(abspdgcode==211) fHisTPCsigmaMCPion->Fill(track->P(),nsigma[0]);
	      if(abspdgcode==321) fHisTPCsigmaMCK->Fill(track->P(),nsigma[1]);
	      if(abspdgcode==2212) fHisTPCsigmaMCProton->Fill(track->P(),nsigma[2]);

	    }

	  }
	  if(fFillDistrTrackEffChecks && track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kTPCrefit){
	    fHisTPCsigNvsPtAllTracks->Fill(track->Pt(),(Float_t)TPCsignalN);
	    fHisTPCsigNvsPhiAllTracks->Fill(track->Phi(),(Float_t)TPCsignalN);
	    fHisTPCsigNvsEtaAllTracks->Fill(track->Eta(),(Float_t)TPCsignalN);
	  }

	}//if TPC status
      } //end PID histograms

      Int_t nclsTot=0,nclsSPD=0;

      //check clusters of the tracks
      if(fOnOff[0]){

	fHisnLayerITS->Fill(-1);
	if(selTrackNoSPD) fHisnLayerITSselTr->Fill(-1);
	for(Int_t l=0;l<6;l++) {
	  if(TESTBIT(track->GetITSClusterMap(),l)) {
	    fHisnLayerITS->Fill(l);
	    if(selTrackNoSPD) fHisnLayerITSselTr->Fill(l);
	    nclsTot++; if(l<2) nclsSPD++;
	  }
	}
	fHisnClsITS->Fill(nclsTot);
	fHisnClsSPD->Fill(nclsSPD);

	if(fFillDistrTrackEffChecks && track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kTPCrefit){

	  fHisPtAllTracks->Fill(track->Pt());
	  fHisPhiAllTracks->Fill(track->Phi());
	  fHisEtaAllTracks->Fill(track->Eta());
	  fHisEtavsPhiAllTracks->Fill(track->Phi(),track->Eta());
	  fHisNTPCclsvsPtAllTracks->Fill(track->Pt(),track->GetTPCNcls());
	  fHisNTPCclsvsPhiAllTracks->Fill(track->Phi(),track->GetTPCNcls());
	  fHisNTPCclsvsEtaAllTracks->Fill(track->Eta(),track->GetTPCNcls());

	  fHisNTPCCrossedRowsvsPtAllTracks->Fill(track->Pt(),nCrossedRowsTPC);
	  fHisNTPCCrossedRowsvsPhiAllTracks->Fill(track->Phi(),nCrossedRowsTPC);
	  fHisNTPCCrossedRowsvsEtaAllTracks->Fill(track->Eta(),nCrossedRowsTPC);

	  fHisRatioCRowsOverFclsvsPtAllTracks->Fill(track->Pt(),ratioCrossedRowsOverFindableClustersTPC);
	  fHisRatioCRowsOverFclsvsPhiAllTracks->Fill(track->Phi(),ratioCrossedRowsOverFindableClustersTPC);
	  fHisRatioCRowsOverFclsvsEtaAllTracks->Fill(track->Eta(),ratioCrossedRowsOverFindableClustersTPC);

	  if(!(track->HasPointOnITSLayer(0)) && !(track->HasPointOnITSLayer(1))){ //no SPD points
	    fHisSPDclsAllTracks->Fill(0);
	  }
	  if(track->HasPointOnITSLayer(0) && !(track->HasPointOnITSLayer(1))){ //kOnlyFirst
	    fHisSPDclsAllTracks->Fill(1);
	  }
	  if(!(track->HasPointOnITSLayer(0)) && track->HasPointOnITSLayer(1)){ //kOnlySecond
	    fHisSPDclsAllTracks->Fill(2);
	  }
	  if(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1)){ //kBoth
	    fHisSPDclsAllTracks->Fill(3);
	  }
	  fHisNITSclsvsPtAllTracks->Fill(track->Pt(), nclsTot);
	  fHisNITSclsvsPhiAllTracks->Fill(track->Phi(), nclsTot);
	  fHisNITSclsvsEtaAllTracks->Fill(track->Eta(), nclsTot);

	}

	if(track->Pt()>0.3 &&
	   TMath::Abs(track->Eta())<0.8 &&
	   track->GetStatus()&AliESDtrack::kITSrefit &&
	   track->GetStatus()&AliESDtrack::kTPCrefit &&
	   nclsSPD>0){
	  fHisnClsITSselTr->Fill(nclsTot);
	}
	if(!(track->GetStatus()&AliESDtrack::kTPCin) && track->GetStatus()&AliESDtrack::kITSrefit && !(track->GetStatus()&AliESDtrack::kITSpureSA)){//tracks retrieved in the ITS and not reconstructed in the TPC
	  fHisnClsITSSA->Fill(nclsTot);
     if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))
       fHisnClsITSSAspdAny->Fill(nclsTot);
     if(track->HasPointOnITSLayer(0))
       fHisnClsITSSAspdIn->Fill(nclsTot);
     if(track->HasPointOnITSLayer(1))
       fHisnClsITSSAspdOut->Fill(nclsTot);
	  fHisnLayerITSsa->Fill(-1);
	  for(Int_t l=0;l<6;l++) {
	    if(TESTBIT(track->GetITSClusterMap(),l)) {
	      fHisnLayerITSsa->Fill(l);
	    }
	  }
	}
	Int_t label=0;
	if(fReadMC){
	  label=track->GetLabel();
	  if (label<0) fHisNentries->Fill(13);
	  else fHisNentries->Fill(14);
	}


	if (track->Pt()>0.3 &&
	    track->GetStatus()&AliESDtrack::kTPCrefit &&
	    track->GetStatus()&AliESDtrack::kITSrefit &&
	    /*nclsTot>3 &&*/
	    nclsSPD>0) {//count good tracks


	  if(fReadMC && label<0) {
	    fHisptFakeTr->Fill(track->Pt());
	    isFakeTrack++;
	  } else {
	    fHisptGoodTr->Fill(track->Pt());
	    isGoodTrack++;
	  }

	  if(fCuts->IsDaughterSelected(track,&vESD,fCuts->GetTrackCuts(),aod)){
	    isSelTrack++;
	  }//select tracks for our analyses

	}
      } //fill track histos
    } //end loop on tracks

    //fill once per event
    if(fOnOff[0]){
      if (fReadMC) fHisdistrFakeTr->Fill(isFakeTrack);
      fHisdistrGoodTr->Fill(isGoodTrack);
      fHisdistrSelTr->Fill(isSelTrack);
    }

    if(!isSimpleMode){
      // loop over candidates
      Int_t nCand = arrayProng->GetEntriesFast();
      Int_t ndaugh=3;
      if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi) ndaugh=2;
      if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpipipi) ndaugh=4;

      for (Int_t iCand = 0; iCand < nCand; iCand++) {
	AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
         if(d->GetIsFilled()<1) continue;//0 = not refilled. or refill failed. 1 = standard dAODs. 2 = succesfully refilled
       	if(fUseSelectionBit && d->GetSelectionMap()) {
	  if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi && !d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) continue; //skip the D0 from Dstar
	  if(fDecayChannel==AliAnalysisTaskSEHFQA::kDplustoKpipi && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) continue; //skip the 3 prong !D+
	}
  if(fDecayChannel==AliAnalysisTaskSEHFQA::kLambdactoV0 && !((dynamic_cast<AliAODRecoCascadeHF*>(d))->CheckCascadeFlags())) continue;

	if(fReadMC){

	  Int_t labD = -1;
	  if (fDecayChannel==AliAnalysisTaskSEHFQA::kLambdactoV0 && (dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0()) {

	    Int_t pdgDgLctoV0bachelor[2]={2212,310};
	    Int_t pdgDgV0toDaughters[2]={211,211};
	    Int_t mcLabelK0S = (dynamic_cast<AliAODRecoCascadeHF*>(d))->MatchToMC(pdg,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE); // Lc->K0S+p and cc
	    pdgDgLctoV0bachelor[1]=3122, pdgDgLctoV0bachelor[0]=211;
	    pdgDgV0toDaughters[0]=2212,  pdgDgV0toDaughters[1]=211;
	    Int_t mcLabelLambda = (dynamic_cast<AliAODRecoCascadeHF*>(d))->MatchToMC(pdg,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE); // Lc->Lambda+pi and cc
	    if (mcLabelK0S!=-1 || mcLabelLambda!=-1) AliInfo(Form("mcLabelK0S=%d - mcLabelLambda=%d",mcLabelK0S,mcLabelLambda));

	    if (mcLabelK0S!=-1 && mcLabelLambda!=-1)
	      AliInfo("Strange: current Lc->V0+bachelor candidate has two MC different labels!");
	    else if (mcLabelK0S>-1 && mcLabelLambda==-1)
	      labD = mcLabelK0S;
	    else if (mcLabelLambda>-1 && mcLabelK0S==-1)
	      labD = mcLabelLambda;
	  }
	  else
	    labD = d->MatchToMC(pdg,mcArray,ndaugh,pdgdaughters);

	  if(labD>=0){
	    AliAODMCParticle *partD = (AliAODMCParticle*)mcArray->At(labD);
	    Int_t label=partD->GetMother();
	    AliAODMCParticle *mot = (AliAODMCParticle*)mcArray->At(label);
	    while(label>=0){//get first mother
	      mot = (AliAODMCParticle*)mcArray->At(label);
	      label=mot->GetMother();
	    }
	    if(mot){
	      Int_t pdgMotCode = mot->GetPdgCode();

	      if(TMath::Abs(pdgMotCode)==4) fHisNentries->Fill(11); //from primary charm
	      if(TMath::Abs(pdgMotCode)==5) fHisNentries->Fill(12); //from beauty
	    }
	  }
	}//end MC
	fHisNentries->Fill(10);//count the candidates (data and MC)

	for(Int_t id=0;id<ndaugh;id++){
	  //other histograms to be filled when the cut object is given
	  AliAODTrack* track=0;

	  if (fDecayChannel==AliAnalysisTaskSEHFQA::kLambdactoV0 && (dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0()) {
	    if (id==0)
	      track=(AliAODTrack*)(dynamic_cast<AliAODRecoCascadeHF*>(d))->GetBachelor();
	    else if (id==1)
	      track=(AliAODTrack*)(dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0PositiveTrack();
	    else if (id==2)
	      track=(AliAODTrack*)(dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0NegativeTrack();
	  }
	  else
	    track=(AliAODTrack*)d->GetDaughter(id);


	  // filtering cut level

	  if (fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(pdg)) && fCuts->IsSelected(d,AliRDHFCuts::kAll,aod)) {

	    Int_t label=0;
	    if(fReadMC)label=track->GetLabel();
	    if(fOnOff[0]){

	      if(fReadMC && label<0) {
		isFakeTrack++;
		fHisptFakeTrFromDaughFilt->Fill(track->Pt());

		fHisd0f_filt->Fill(d->Getd0Prong(id));
	      } else {
		fHisptGoodTrFromDaugh_filt->Fill(track->Pt());
		fHisd0dau_filt->Fill(d->Getd0Prong(id));
                Double_t phidaughter = d->PhiProng(id);
		if(phidaughter<0) phidaughter=2.0*TMath::Pi()+phidaughter;
		fHisd0dauphi_filt->Fill(phidaughter, d->Getd0Prong(id));
		Double_t d0rphiz[2],covd0[3];
		Bool_t isDCA=track->PropagateToDCA(vtx1,magField,9999.,d0rphiz,covd0);
		if(isDCA){
		  fHisd0zdau_filt->Fill(d0rphiz[1]);
		  fHisd0zdauphi_filt->Fill(phidaughter,d0rphiz[1]);
		}
	      }
	    }
	  }



	  //track quality

	  if (fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(pdg)) && fCuts->IsSelected(d,AliRDHFCuts::kTracks,aod)) {

	    Int_t label=0;
	    if(fReadMC)label=track->GetLabel();
	    if(fOnOff[0]){

	      if(fReadMC && label<0) {
		isFakeTrack++;
		fHisptFakeTrFromDaugh->Fill(track->Pt());

		fHisd0f->Fill(d->Getd0Prong(id));
	      } else {
		fHisptGoodTrFromDaugh->Fill(track->Pt());
		fHisd0dau->Fill(d->Getd0Prong(id));
                Double_t phidaughter = d->PhiProng(id);
		if(phidaughter<0) phidaughter=2.0*TMath::Pi()+phidaughter;
		fHisd0dauphi->Fill(phidaughter, d->Getd0Prong(id));
		Double_t d0rphiz[2],covd0[3];
		Bool_t isDCA=track->PropagateToDCA(vtx1,magField,9999.,d0rphiz,covd0);
		if(isDCA){
		  fHisd0zdau->Fill(d0rphiz[1]);
		  fHisd0zdauphi->Fill(phidaughter,d0rphiz[1]);
		}
	      }
	    }


	    if(fFillDistrTrackEffChecks){
	      Int_t nITScls = 0;
	      Double_t nTPCCrossedRows = track->GetTPCClusterInfo(2,1);
	      Double_t ratioCrossedRowsOverFcls = 1.0;
	      if(track->GetTPCNclsF()>0){
		ratioCrossedRowsOverFcls = (nTPCCrossedRows)/(track->GetTPCNclsF());
	      }
	      for(Int_t l=0;l<6;l++) {
		if(TESTBIT(track->GetITSClusterMap(),l)) {
		  nITScls++;
		}
	      }

	      fHisPtDaughters->Fill(track->Pt());
	      fHisPhiDaughters->Fill(track->Phi());
	      fHisEtaDaughters->Fill(track->Eta());
	      fHisEtavsPhiDaughters->Fill(track->Phi(),track->Eta());

	      fHisNTPCclsvsPtDaughters->Fill(track->Pt(),track->GetTPCNcls());
	      fHisNTPCclsvsPhiDaughters->Fill(track->Phi(),track->GetTPCNcls());
	      fHisNTPCclsvsEtaDaughters->Fill(track->Eta(),track->GetTPCNcls());

	      fHisNTPCCrossedRowsvsPtDaughters->Fill(track->Pt(),nTPCCrossedRows);
	      fHisNTPCCrossedRowsvsPhiDaughters->Fill(track->Phi(),nTPCCrossedRows);
	      fHisNTPCCrossedRowsvsEtaDaughters->Fill(track->Eta(),nTPCCrossedRows);

	      fHisRatioCRowsOverFclsvsPtDaughters->Fill(track->Pt(),ratioCrossedRowsOverFcls);
	      fHisRatioCRowsOverFclsvsPhiDaughters->Fill(track->Phi(),ratioCrossedRowsOverFcls);
	      fHisRatioCRowsOverFclsvsEtaDaughters->Fill(track->Eta(),ratioCrossedRowsOverFcls);

	      fHisNITSclsvsPtDaughters->Fill(track->Pt(), nITScls);
	      fHisNITSclsvsPhiDaughters->Fill(track->Phi(), nITScls);
	      fHisNITSclsvsEtaDaughters->Fill(track->Eta(), nITScls);
	      if(!(track->HasPointOnITSLayer(0)) && !(track->HasPointOnITSLayer(1))){ //no SPD points
		fHisSPDclsDaughters->Fill(0);
	      }
	      if(track->HasPointOnITSLayer(0) && !(track->HasPointOnITSLayer(1))){ //kOnlyFirst
		fHisSPDclsDaughters->Fill(1);
	      }
	      if(!(track->HasPointOnITSLayer(0)) && track->HasPointOnITSLayer(1)){ //kOnlySecond
		fHisSPDclsDaughters->Fill(2);
	      }
	      if(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1)){ //kBoth
		fHisSPDclsDaughters->Fill(3);
	      }


	      if(fOnOff[1]){
		AliAODPid *pid = track->GetDetPid();
		if(pid){
		  if(pidHF && pidHF->CheckStatus(track,"TPC")){
		    fHisTPCsigNvsPtDaughters->Fill(track->Pt(),pid->GetTPCsignalN());
		    fHisTPCsigNvsPhiDaughters->Fill(track->Phi(),pid->GetTPCsignalN());
		    fHisTPCsigNvsEtaDaughters->Fill(track->Eta(),pid->GetTPCsignalN());
		  }
		}
	      }
	    }


	    if (fCuts->IsSelected(d,AliRDHFCuts::kAll,aod) && fOnOff[1]){
	       fHisNentries->Fill(8); //candidates passing analysis cuts

	    AliAODPid *pid = track->GetDetPid();
	      if(pid){
		Double_t times[5];
		pid->GetIntegratedTimes(times);
		if(pidHF && pidHF->CheckStatus(track,"TOF")){
		  Double_t tofTime=pid->GetTOFsignal();
		  AliTOFHeader* tofH=(AliTOFHeader*)aod->GetTOFHeader();
		  Double_t tofRes[AliPID::kSPECIES];
		  pid->GetTOFpidResolution(tofRes);
		  if (tofH && (TMath::Abs(tofRes[0]) <= 1.E-16) ) { // new AOD
		    AliTOFPIDResponse tofResp=pidHF->GetPidResponse()->GetTOFResponse();
		    Double_t startTime=tofResp.GetStartTime(track->P());
		    tofTime-=startTime;
		  }
		  fHisTOFtimeKaonHyptimeAC->Fill(track->P(),tofTime-times[AliPID::kKaon]);
		}
		if(pidHF && pidHF->CheckStatus(track,"TPC")) fHisTPCsigvspAC->Fill(pid->GetTPCmomentum(),pid->GetTPCsignal());
	      }

	    } //end analysis cuts
	  } //end acceptance and track cuts
	} //end loop on tracks in the candidate
      } //end loop on candidates

    }
  } //end if on pid or track histograms
  delete vHF;
  delete tpcres;
  delete [] pdgdaughters;
  PostData(1,fOutputEntries);
  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) PostData(5,fOutputCounters);
  //Post data 6 done in case of centrality on
}

//____________________________________________________________________________
void AliAnalysisTaskSEHFQA::FillFlowObs(AliAODEvent *aod){
  /// fills the flow observables
  Double_t cc;
  cc = fCuts->GetCentrality(aod);
  fHisFEvents->Fill(0., cc);

  UInt_t mask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  UInt_t trigger=AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  if(mask & trigger) {
    fHisFEvents->Fill(1.,cc); // fired
    if (mask & AliVEvent::kMB) fHisFEvents->Fill(2.,cc);
    if (mask & AliVEvent::kCentral) fHisFEvents->Fill(3.,cc);
    if (mask & AliVEvent::kSemiCentral) fHisFEvents->Fill(4.,cc);
    Bool_t rejected=false;
    if(cc<0 || cc>60) rejected=true;
    const AliVVertex *vertex = aod->GetPrimaryVertex();
    Double_t zvtx=vertex->GetZ();
    if(TMath::Abs(zvtx)>fCuts->GetMaxVtxZ()) rejected=true;
    if(rejected) return; //not interesting for flow QA
  } else {
    return;
  }

  // event accepted
  fHisFEvents->Fill(5.,cc);
  fRFPcuts->SetParamType(AliFlowTrackCuts::kGlobal);
  fRFPcuts->SetPtRange(0.2,5.);
  fRFPcuts->SetEtaRange(-0.8,0.8);
  fRFPcuts->SetMinNClustersTPC(70);
  fRFPcuts->SetMinChi2PerClusterTPC(0.2);
  fRFPcuts->SetMaxChi2PerClusterTPC(4.0);
  fRFPcuts->SetAcceptKinkDaughters(kFALSE);
  fRFPcuts->SetEvent(aod);

  // "FB1" (i=0), "FB128" (i=1), "VZE" (i=2)
  Double_t psi[3];
  for(Int_t i=0; i!=3; ++i) {
    if(i==0) { // switching to bit 1
      fRFPcuts->SetMinimalTPCdedx(10.);
      fRFPcuts->SetAODfilterBit(1);
    } else { // switching to bit 128
      fRFPcuts->SetMinimalTPCdedx(-1);
      fRFPcuts->SetAODfilterBit(128);
    }
    if(i>1) {
      fRFPcuts->SetParamType(AliFlowTrackCuts::kVZERO);
      fRFPcuts->SetEtaRange(-5,+5);
      fRFPcuts->SetPhiMin(0);
      fRFPcuts->SetPhiMax(TMath::TwoPi());
    }
    fFlowEvent->Fill(fRFPcuts,fRFPcuts);
    fFlowEvent->TagSubeventsInEta(-5,0,0,+5);
    // getting informationt
    AliFlowVector vQ, vQaQb[2];
    fFlowEvent->Get2Qsub(vQaQb,2);
    vQ = vQaQb[0]+vQaQb[1];
    Double_t dMa=vQaQb[0].GetMult();
    Double_t dMb=vQaQb[1].GetMult();
    if( dMa<2 || dMb<2 ) {
      fHisFEvents->Fill(6.,cc); //???
      continue;
    }
    psi[i] = vQ.Phi()/2;
    // publishing
    fHisQ[i]->Fill(0,cc,vQaQb[0].X()/dMa,dMa); // Qx-
    fHisQ[i]->Fill(1,cc,vQaQb[0].Y()/dMa,dMa); // Qy-
    fHisQ[i]->Fill(2,cc,vQaQb[1].X()/dMb,dMb); // Qx+
    fHisQ[i]->Fill(3,cc,vQaQb[1].Y()/dMb,dMb); // Qy+
    fHisAngleQ[i]->Fill(psi[i],cc); // Psi
    AliFlowTrackSimple *track;
    for(Int_t t=0; t!=fFlowEvent->NumberOfTracks(); ++t) {
      track = (AliFlowTrackSimple*) fFlowEvent->GetTrack(t);
      if(!track) continue;
      if(!track->InRPSelection()) continue;
      fHisPhiEta[i]->Fill(track->Phi(),track->Eta(),cc,track->Weight()); //PhiEta
    }

  //histo filled only for TPCFB1
  if (i==0) {
    fHisCentVsMultRPS->Fill(fFlowEvent->GetNumberOfRPs(),cc);
  }
  }
  // TPC vs VZERO
  fHisTPCVZE_AngleQ->Fill(psi[0],psi[2],cc);
}

//____________________________________________________________________________
void AliAnalysisTaskSEHFQA::Terminate(Option_t */*option*/){
  /// terminate analysis

 fOutputEntries = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputEntries && fOnOff[1]) {
    printf("ERROR: %s not available\n",GetOutputSlot(1)->GetContainer()->GetName());
    return;
  }

  fOutputPID = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputPID && fOnOff[1]) {
    printf("ERROR: %s not available\n",GetOutputSlot(2)->GetContainer()->GetName());
    return;
  }

  fOutputTrack = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputTrack && fOnOff[0]) {
    printf("ERROR: %s not available\n",GetOutputSlot(3)->GetContainer()->GetName());
    return;
  }

}



