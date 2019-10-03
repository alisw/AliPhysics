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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTask to extract from ESD tracks the information
// on ITS tracking efficiency and resolutions.
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TStyle.h>
#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>  
#include <TCanvas.h>
#include <TParticle.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliTriggerConfiguration.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliVertexerTracks.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliESDtrackCuts.h"
#include "AliTrackPointArray.h"
#include "AliITSRecPoint.h"

#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliLog.h"

#include "AliGenEventHeader.h" 
#include "AliAnalysisTaskITSTrackingCheck.h"

ClassImp(AliAnalysisTaskITSTrackingCheck)
AliAnalysisTaskITSTrackingCheck::AliAnalysisTaskITSTrackingCheck() : 
AliAnalysisTaskSE(), 
fNITSLayers(6),
fReadMC(kFALSE),
fReadRPLabels(kFALSE),
fFillNtuples(kFALSE),
fUseITSSAforNtuples(kFALSE),
fUsePhysSel(kFALSE),
fRequireSPDvtx(kTRUE),
fRejPileupSPD(kFALSE),
fCheckSDDIsIn(0),
fESD(0), 
fMinMult(0),
fMaxMult(1000000),
fTriggerClass(""),
fTrigConfig(0),
fUsePtBinsForBox(kFALSE),
fOCDBPath(""),
fOutput(0), 
fHistNEvents(0),
fHistNEventsFrac(0),
fHistNtracks(0),
fHistNclsITSMI(0),
fHistNclsITSSA(0),
fHistNclsITSSAInAcc(0),
fHistClusterMapITSMI(0),
fHistClusterMapITSMIok(0),
fHistClusterMapITSMIbad(0),
fHistClusterMapITSMIskipped(0),
fHistClusterMapITSMIoutinz(0),
fHistClusterMapITSMInorefit(0),
fHistClusterMapITSMInocls(0),
fHistClusterMapITSMIokoutinzbad(0),
fHistClusterMapITSSA(0),
fHistClusterMapITSSAok(0),
fHistClusterMapITSSAbad(0),
fHistClusterMapITSSAskipped(0),
fHistClusterMapITSSAoutinz(0),
fHistClusterMapITSSAnorefit(0),
fHistClusterMapITSSAnocls(0),
fHistClusterMapITSSAokoutinzbad(0),
fHistClusterMapITSSAInAcc(0),
fHistClusterMapITSSAokInAcc(0),
fHistClusterMapITSSAbadInAcc(0),
fHistClusterMapITSSAskippedInAcc(0),
fHistClusterMapITSSAoutinzInAcc(0),
fHistClusterMapITSSAnorefitInAcc(0),
fHistClusterMapITSSAnoclsInAcc(0),
fHistClusterMapITSSAokoutinzbadInAcc(0),
fHistClusterMapModuleITSSAokInAcc(0),
fHistClusterMapModuleITSSAbadInAcc(0),
fHistClusterMapModuleITSSAnoclsInAcc(0),
fHistClusterMapModuleITSMIokInAcc(0),
fHistClusterMapModuleITSMIbadInAcc(0),
fHistClusterMapModuleITSMInoclsInAcc(0),
fHistNClustersMapModule(0),
fHistZatSPDouter0ok(0),
fHistZatSPDouter1ok(0),
fHistZatSPDouter2ok(0),
fHistZatSPDouter3ok(0),
fHistZatSPDouter0notok(0),
fHistZatSPDouter1notok(0),
fHistZatSPDouter2notok(0),
fHistZatSPDouter3notok(0),
fHistxlocSDDok(0),
fHistzlocSDDok(0),
fHistxlocVSmodSDDok(0),
fHistxlocSDDall(0),
fHistzlocSDDall(0),
fHistPhiTPCInAcc(0),
fHistEtaTPCInAcc(0),
fHistPtTPC(0),
fHistPtTPCInAcc(0),
fHistPtTPCInAccTOFbc0(0),
fHistPtTPCInAccwSDD(0),
fHistPtTPCInAccTOFbc0wSDD(0),
fHistPtTPCInAccwoSDD(0),
fHistPtTPCInAccTOFbc0woSDD(0),
fHistPtTPCInAccMCtwoSPD(0),
fHistPtTPCInAccMConeSPD(0),
fHistdEdxVSPtTPCInAcc(0),
fHistdEdxVSPtITSTPCsel(0),
fHistTPCclsVSPtTPCInAcc(0),
fHistTPCclsVSPtITSMISPDInAcc(0),
fHistPtVSphiTPCInAcc(0),
fHistPtTPCInAccNoTRDout(0),
fHistPtTPCInAccNoTOFout(0),
fHistPtTPCInAccWithPtTPCAtInnerWall(0),
fHistPtTPCInAccWithPtTPCAtVtx(0),
fHistDeltaPtTPC(0),
fHistPtTPCInAccP(0),
fHistPtTPCInAccS(0),
fHistPtTPCInAccPfromStrange(0),
fHistPtTPCInAccSfromStrange(0),
fHistPtTPCInAccSfromMat(0),
fHistPtITSMI2(0),
fHistPtITSMI3(0),
fHistPtITSMI4(0),
fHistPtITSMI5(0),
fHistPtITSMI6(0),
fHistPtITSMI7(0),
fHistPtITSMISPD(0),
fHistPtITSMIoneSPD(0),
fHistPtITSMItwoSPD(0),
fHistPtITSMI2InAcc(0),
fHistPtITSMI3InAcc(0),
fHistPtITSMI4InAcc(0),
fHistPtITSMI5InAcc(0),
fHistPtITSMI6InAcc(0),
fHistPtITSMI7InAcc(0),
fHistPtITSMISPDInAcc(0),
fHistPtITSMIoneSPDInAcc(0),
fHistPtITSMItwoSPDInAcc(0),
fHistPtITSMI2InAccTOFbc0(0),
fHistPtITSMI3InAccTOFbc0(0),
fHistPtITSMI4InAccTOFbc0(0),
fHistPtITSMI5InAccTOFbc0(0),
fHistPtITSMI6InAccTOFbc0(0),
fHistPtITSMISPDInAccTOFbc0(0),
fHistPtITSMIoneSPDInAccTOFbc0(0),
fHistPtITSMI2InAccwSDD(0),
fHistPtITSMI3InAccwSDD(0),
fHistPtITSMI4InAccwSDD(0),
fHistPtITSMI5InAccwSDD(0),
fHistPtITSMI6InAccwSDD(0),
fHistPtITSMISPDInAccwSDD(0),
fHistPtITSMIoneSPDInAccwSDD(0),
fHistPtITSMI2InAccTOFbc0wSDD(0),
fHistPtITSMI3InAccTOFbc0wSDD(0),
fHistPtITSMI4InAccTOFbc0wSDD(0),
fHistPtITSMI5InAccTOFbc0wSDD(0),
fHistPtITSMI6InAccTOFbc0wSDD(0),
fHistPtITSMISPDInAccTOFbc0wSDD(0),
fHistPtITSMIoneSPDInAccTOFbc0wSDD(0),
fHistPtITSMI2InAccwoSDD(0),
fHistPtITSMI3InAccwoSDD(0),
fHistPtITSMI4InAccwoSDD(0),
fHistPtITSMI5InAccwoSDD(0),
fHistPtITSMI6InAccwoSDD(0),
fHistPtITSMISPDInAccwoSDD(0),
fHistPtITSMIoneSPDInAccwoSDD(0),
fHistPtITSMI2InAccTOFbc0woSDD(0),
fHistPtITSMI3InAccTOFbc0woSDD(0),
fHistPtITSMI4InAccTOFbc0woSDD(0),
fHistPtITSMI5InAccTOFbc0woSDD(0),
fHistPtITSMI6InAccTOFbc0woSDD(0),
fHistPtITSMISPDInAccTOFbc0woSDD(0),
fHistPtITSMIoneSPDInAccTOFbc0woSDD(0),
fHistPtITSMIoneSPDInAccShared(0),
fHistPtITSMIoneSPDInAccSharedSPD(0),
fHistPtITSMISPD1InAccShared(0),
fHistPtITSMISPD2InAccShared(0),
fHistPtITSMIoneSPDInAccSharedFake(0),
fHistPtITSMIoneSPDInAccSharedSPDFake(0),
fHistPhiITSMI2InAcc(0),
fHistPhiITSMI3InAcc(0),
fHistPhiITSMI4InAcc(0),
fHistPhiITSMI5InAcc(0),
fHistPhiITSMI6InAcc(0),
fHistPhiITSMI7InAcc(0),
fHistPhiITSMISPDInAcc(0),
fHistPhiITSMIoneSPDInAcc(0),
fHistPhiITSMItwoSPDInAcc(0),
fHistEtaITSMI2InAcc(0),
fHistEtaITSMI3InAcc(0),
fHistEtaITSMI4InAcc(0),
fHistEtaITSMI5InAcc(0),
fHistEtaITSMI6InAcc(0),
fHistEtaITSMI7InAcc(0),
fHistEtaITSMISPDInAcc(0),
fHistEtaITSMIoneSPDInAcc(0),
fHistEtaITSMItwoSPDInAcc(0),
fHistPtITSMI2InAccFake(0),
fHistPtITSMI3InAccFake(0),
fHistPtITSMI4InAccFake(0),
fHistPtITSMI5InAccFake(0),
fHistPtITSMI6InAccFake(0),
fHistPtITSMI7InAccFake(0),
fHistPtITSMISPDInAccFake(0),
fHistPtITSMIoneSPDInAccFake(0),
fHistPtITSMItwoSPDInAccFake(0),
fHistPtITSMIoneSPDthreeSDDSSDInAcc(0),
fHistPtITSTPCsel(0),
fHistPtITSTPCselTOFbc0(0),
fHistPtITSTPCselwSDD(0),
fHistPtITSTPCselTOFbc0wSDD(0),
fHistPtITSTPCselwoSDD(0),
fHistPtITSTPCselTOFbc0woSDD(0),
fHistPtITSTPCselP(0),
fHistPtITSTPCselS(0),
fHistPtITSTPCselFake(0),
fHistPtITSTPCselPfromStrange(0),
fHistPtITSTPCselSfromStrange(0),
fHistPtITSTPCselSfromMat(0),
fHistPtITSMI2InAccP(0),
fHistPtITSMI3InAccP(0),
fHistPtITSMI4InAccP(0),
fHistPtITSMI5InAccP(0),
fHistPtITSMI6InAccP(0),
fHistPtITSMI7InAccP(0),
fHistPtITSMISPDInAccP(0),
fHistPtITSMIoneSPDInAccP(0),
fHistPtITSMItwoSPDInAccP(0),
fHistPtITSMI2InAccS(0),
fHistPtITSMI3InAccS(0),
fHistPtITSMI4InAccS(0),
fHistPtITSMI5InAccS(0),
fHistPtITSMI6InAccS(0),
fHistPtITSMI7InAccS(0),
fHistPtITSMISPDInAccS(0),
fHistPtITSMIoneSPDInAccS(0),
fHistPtITSMItwoSPDInAccS(0),
fHistPtITSMIokbadoutinz6(0),
fHistPtITSMIokbadoutinz4InAcc(0),
fHistPtITSMIokbadoutinz5InAcc(0),
fHistPtITSMIokbadoutinz6InAcc(0),
fHistPhiITSMIokbadoutinz6InAcc(0),
fHistRProdVtxInAccP(0),
fHistRProdVtxInAccS(0),
fHistd0rphiTPCInAccP150200(0),
fHistd0rphiTPCInAccP500700(0),
fHistd0rphiTPCInAccP10001500(0),
fHistd0rphiTPCInAccS150200(0),
fHistd0rphiTPCInAccS500700(0),
fHistd0rphiTPCInAccS10001500(0),
fHistd0rphiITSMISPDInAccP150200(0),
fHistd0rphiITSMISPDInAccP500700(0),
fHistd0rphiITSMISPDInAccP10001500(0),
fHistd0rphiITSMISPDInAccS150200(0),
fHistd0rphiITSMISPDInAccS500700(0),
fHistd0rphiITSMISPDInAccS10001500(0),
fHistd0rphiITSMIoneSPDInAccP150200(0),
fHistd0rphiITSMIoneSPDInAccP350450(0),
fHistd0rphiITSMIoneSPDInAccP500700(0),
fHistd0rphiITSMIoneSPDInAccP10001500(0),
fHistd0rphiITSMIoneSPDInAccP25004000(0),
fHistd0rphiITSMIoneSPDInAccP40008000(0),
fHistd0zITSMIoneSPDInAccP150200(0),
fHistd0zITSMIoneSPDInAccP500700(0),
fHistd0zITSMIoneSPDInAccP10001500(0),
fHistd0zVSetaTPCInAccP10001500(0),
fHistd0rphiVSphiITSMIoneSPDInAccP10001500(0),
fHistd0rphiVSetaITSMIoneSPDInAccP10001500(0),
fHistd0rphiITSMIoneSPDInAccS150200(0),
fHistd0rphiITSMIoneSPDInAccS350450(0),
fHistd0rphiITSMIoneSPDInAccS500700(0),
fHistd0rphiITSMIoneSPDInAccS500700from22(0),
fHistd0rphiITSMIoneSPDInAccS500700from211(0),
fHistd0rphiITSMIoneSPDInAccS500700from310(0),
fHistd0rphiITSMIoneSPDInAccS500700from321(0),
fHistd0rphiITSMIoneSPDInAccS500700from3122(0),
fHistd0rphiITSMIoneSPDInAccS10001500(0),
fHistd0rphiITSMIoneSPDInAccS25004000(0),
fHistd0rphiITSMIoneSPDInAccS40008000(0),
fHistd0rphiITSMIoneSPDInAccS150200fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS150200fromMat(0),
fHistd0rphiITSMIoneSPDInAccS350450fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS350450fromMat(0),
fHistd0rphiITSMIoneSPDInAccS500700fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS500700fromMat(0),
fHistd0rphiITSMIoneSPDInAccS10001500fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS10001500fromMat(0),
fHistd0rphiITSMIoneSPDInAccS25004000fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS25004000fromMat(0),
fHistd0rphiITSMIoneSPDInAccS40008000fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS40008000fromMat(0),
fHistd0zITSMIoneSPDInAccS150200(0),
fHistd0zITSMIoneSPDInAccS500700(0),
fHistd0zITSMIoneSPDInAccS10001500(0),
fHistPDGMoth(0),
fHistPDGMoth150200(0),
fHistPDGMoth500700(0),
fHistPDGMoth10001500(0),
fHistPDGTrk(0),
fHistITSRedChi2NonFakePt02(0),
fHistITSRedChi2FakePt02(0),
fHistITSRedChi2NonFakePt05(0),
fHistITSRedChi2FakePt05(0),
fHistITSRedChi2NonFakePt1(0),
fHistITSRedChi2FakePt1(0),
fNtupleESDTracks(0),
fNtupleITSAlignExtra(0),
fNtupleITSAlignSPDTracklets(0),
fESDtrackCutsTPC(0),
fESDtrackCutsITSTPC(0)
{
  // Constructor
  for(Int_t i=0; i<11; i++) fCountsPerPtBin[i]=0;

}

//________________________________________________________________________
AliAnalysisTaskITSTrackingCheck::AliAnalysisTaskITSTrackingCheck(const char *name) : 
AliAnalysisTaskSE(name), 
fNITSLayers(6),
fReadMC(kFALSE),
fReadRPLabels(kFALSE),
fFillNtuples(kFALSE),
fUseITSSAforNtuples(kFALSE),
fUsePhysSel(kFALSE),
fRequireSPDvtx(kTRUE),
fRejPileupSPD(kFALSE),
fCheckSDDIsIn(0),
fESD(0), 
fMinMult(0),
fMaxMult(1000000),
fTriggerClass(""),
fTrigConfig(0),
fUsePtBinsForBox(kFALSE),
fOCDBPath(""),
fOutput(0), 
fHistNEvents(0),
fHistNEventsFrac(0),
fHistNtracks(0),
fHistNclsITSMI(0),
fHistNclsITSSA(0),
fHistNclsITSSAInAcc(0),
fHistClusterMapITSMI(0),
fHistClusterMapITSMIok(0),
fHistClusterMapITSMIbad(0),
fHistClusterMapITSMIskipped(0),
fHistClusterMapITSMIoutinz(0),
fHistClusterMapITSMInorefit(0),
fHistClusterMapITSMInocls(0),
fHistClusterMapITSMIokoutinzbad(0),
fHistClusterMapITSSA(0),
fHistClusterMapITSSAok(0),
fHistClusterMapITSSAbad(0),
fHistClusterMapITSSAskipped(0),
fHistClusterMapITSSAoutinz(0),
fHistClusterMapITSSAnorefit(0),
fHistClusterMapITSSAnocls(0),
fHistClusterMapITSSAokoutinzbad(0),
fHistClusterMapITSSAInAcc(0),
fHistClusterMapITSSAokInAcc(0),
fHistClusterMapITSSAbadInAcc(0),
fHistClusterMapITSSAskippedInAcc(0),
fHistClusterMapITSSAoutinzInAcc(0),
fHistClusterMapITSSAnorefitInAcc(0),
fHistClusterMapITSSAnoclsInAcc(0),
fHistClusterMapITSSAokoutinzbadInAcc(0),
fHistClusterMapModuleITSSAokInAcc(0),
fHistClusterMapModuleITSSAbadInAcc(0),
fHistClusterMapModuleITSSAnoclsInAcc(0),
fHistClusterMapModuleITSMIokInAcc(0),
fHistClusterMapModuleITSMIbadInAcc(0),
fHistClusterMapModuleITSMInoclsInAcc(0),
fHistNClustersMapModule(0),
fHistZatSPDouter0ok(0),
fHistZatSPDouter1ok(0),
fHistZatSPDouter2ok(0),
fHistZatSPDouter3ok(0),
fHistZatSPDouter0notok(0),
fHistZatSPDouter1notok(0),
fHistZatSPDouter2notok(0),
fHistZatSPDouter3notok(0),
fHistxlocSDDok(0),
fHistzlocSDDok(0),
fHistxlocVSmodSDDok(0),
fHistxlocSDDall(0),
fHistzlocSDDall(0),
fHistPhiTPCInAcc(0),
fHistEtaTPCInAcc(0),
fHistPtTPC(0),
fHistPtTPCInAcc(0),
fHistPtTPCInAccTOFbc0(0),
fHistPtTPCInAccwSDD(0),
fHistPtTPCInAccTOFbc0wSDD(0),
fHistPtTPCInAccwoSDD(0),
fHistPtTPCInAccTOFbc0woSDD(0),
fHistPtTPCInAccMCtwoSPD(0),
fHistPtTPCInAccMConeSPD(0),
fHistdEdxVSPtTPCInAcc(0),
fHistdEdxVSPtITSTPCsel(0),
fHistTPCclsVSPtTPCInAcc(0),
fHistTPCclsVSPtITSMISPDInAcc(0),
fHistPtVSphiTPCInAcc(0),
fHistPtTPCInAccNoTRDout(0),
fHistPtTPCInAccNoTOFout(0),
fHistPtTPCInAccWithPtTPCAtInnerWall(0),
fHistPtTPCInAccWithPtTPCAtVtx(0),
fHistDeltaPtTPC(0),
fHistPtTPCInAccP(0),
fHistPtTPCInAccS(0),
fHistPtTPCInAccPfromStrange(0),
fHistPtTPCInAccSfromStrange(0),
fHistPtTPCInAccSfromMat(0),
fHistPtITSMI2(0),
fHistPtITSMI3(0),
fHistPtITSMI4(0),
fHistPtITSMI5(0),
fHistPtITSMI6(0),
fHistPtITSMI7(0),
fHistPtITSMISPD(0),
fHistPtITSMIoneSPD(0),
fHistPtITSMItwoSPD(0),
fHistPtITSMI2InAcc(0),
fHistPtITSMI3InAcc(0),
fHistPtITSMI4InAcc(0),
fHistPtITSMI5InAcc(0),
fHistPtITSMI6InAcc(0),
fHistPtITSMI7InAcc(0),
fHistPtITSMISPDInAcc(0),
fHistPtITSMIoneSPDInAcc(0),
fHistPtITSMItwoSPDInAcc(0),
fHistPtITSMI2InAccTOFbc0(0),
fHistPtITSMI3InAccTOFbc0(0),
fHistPtITSMI4InAccTOFbc0(0),
fHistPtITSMI5InAccTOFbc0(0),
fHistPtITSMI6InAccTOFbc0(0),
fHistPtITSMISPDInAccTOFbc0(0),
fHistPtITSMIoneSPDInAccTOFbc0(0),
fHistPtITSMI2InAccwSDD(0),
fHistPtITSMI3InAccwSDD(0),
fHistPtITSMI4InAccwSDD(0),
fHistPtITSMI5InAccwSDD(0),
fHistPtITSMI6InAccwSDD(0),
fHistPtITSMISPDInAccwSDD(0),
fHistPtITSMIoneSPDInAccwSDD(0),
fHistPtITSMI2InAccTOFbc0wSDD(0),
fHistPtITSMI3InAccTOFbc0wSDD(0),
fHistPtITSMI4InAccTOFbc0wSDD(0),
fHistPtITSMI5InAccTOFbc0wSDD(0),
fHistPtITSMI6InAccTOFbc0wSDD(0),
fHistPtITSMISPDInAccTOFbc0wSDD(0),
fHistPtITSMIoneSPDInAccTOFbc0wSDD(0),
fHistPtITSMI2InAccwoSDD(0),
fHistPtITSMI3InAccwoSDD(0),
fHistPtITSMI4InAccwoSDD(0),
fHistPtITSMI5InAccwoSDD(0),
fHistPtITSMI6InAccwoSDD(0),
fHistPtITSMISPDInAccwoSDD(0),
fHistPtITSMIoneSPDInAccwoSDD(0),
fHistPtITSMI2InAccTOFbc0woSDD(0),
fHistPtITSMI3InAccTOFbc0woSDD(0),
fHistPtITSMI4InAccTOFbc0woSDD(0),
fHistPtITSMI5InAccTOFbc0woSDD(0),
fHistPtITSMI6InAccTOFbc0woSDD(0),
fHistPtITSMISPDInAccTOFbc0woSDD(0),
fHistPtITSMIoneSPDInAccTOFbc0woSDD(0),
fHistPtITSMIoneSPDInAccShared(0),
fHistPtITSMIoneSPDInAccSharedSPD(0),
fHistPtITSMISPD1InAccShared(0),
fHistPtITSMISPD2InAccShared(0),
fHistPtITSMIoneSPDInAccSharedFake(0),
fHistPtITSMIoneSPDInAccSharedSPDFake(0),
fHistPhiITSMI2InAcc(0),
fHistPhiITSMI3InAcc(0),
fHistPhiITSMI4InAcc(0),
fHistPhiITSMI5InAcc(0),
fHistPhiITSMI6InAcc(0),
fHistPhiITSMI7InAcc(0),
fHistPhiITSMISPDInAcc(0),
fHistPhiITSMIoneSPDInAcc(0),
fHistPhiITSMItwoSPDInAcc(0),
fHistEtaITSMI2InAcc(0),
fHistEtaITSMI3InAcc(0),
fHistEtaITSMI4InAcc(0),
fHistEtaITSMI5InAcc(0),
fHistEtaITSMI6InAcc(0),
fHistEtaITSMI7InAcc(0),
fHistEtaITSMISPDInAcc(0),
fHistEtaITSMIoneSPDInAcc(0),
fHistEtaITSMItwoSPDInAcc(0),
fHistPtITSMI2InAccFake(0),
fHistPtITSMI3InAccFake(0),
fHistPtITSMI4InAccFake(0),
fHistPtITSMI5InAccFake(0),
fHistPtITSMI6InAccFake(0),
fHistPtITSMI7InAccFake(0),
fHistPtITSMISPDInAccFake(0),
fHistPtITSMIoneSPDInAccFake(0),
fHistPtITSMItwoSPDInAccFake(0),
fHistPtITSMIoneSPDthreeSDDSSDInAcc(0),
fHistPtITSTPCsel(0),
fHistPtITSTPCselTOFbc0(0),
fHistPtITSTPCselwSDD(0),
fHistPtITSTPCselTOFbc0wSDD(0),
fHistPtITSTPCselwoSDD(0),
fHistPtITSTPCselTOFbc0woSDD(0),
fHistPtITSTPCselP(0),
fHistPtITSTPCselS(0),
fHistPtITSTPCselFake(0),
fHistPtITSTPCselPfromStrange(0),
fHistPtITSTPCselSfromStrange(0),
fHistPtITSTPCselSfromMat(0),
fHistPtITSMI2InAccP(0),
fHistPtITSMI3InAccP(0),
fHistPtITSMI4InAccP(0),
fHistPtITSMI5InAccP(0),
fHistPtITSMI6InAccP(0),
fHistPtITSMI7InAccP(0),
fHistPtITSMISPDInAccP(0),
fHistPtITSMIoneSPDInAccP(0),
fHistPtITSMItwoSPDInAccP(0),
fHistPtITSMI2InAccS(0),
fHistPtITSMI3InAccS(0),
fHistPtITSMI4InAccS(0),
fHistPtITSMI5InAccS(0),
fHistPtITSMI6InAccS(0),
fHistPtITSMI7InAccS(0),
fHistPtITSMISPDInAccS(0),
fHistPtITSMIoneSPDInAccS(0),
fHistPtITSMItwoSPDInAccS(0),
fHistPtITSMIokbadoutinz6(0),
fHistPtITSMIokbadoutinz4InAcc(0),
fHistPtITSMIokbadoutinz5InAcc(0),
fHistPtITSMIokbadoutinz6InAcc(0),
fHistPhiITSMIokbadoutinz6InAcc(0),
fHistRProdVtxInAccP(0),
fHistRProdVtxInAccS(0),
fHistd0rphiTPCInAccP150200(0),
fHistd0rphiTPCInAccP500700(0),
fHistd0rphiTPCInAccP10001500(0),
fHistd0rphiTPCInAccS150200(0),
fHistd0rphiTPCInAccS500700(0),
fHistd0rphiTPCInAccS10001500(0),
fHistd0rphiITSMISPDInAccP150200(0),
fHistd0rphiITSMISPDInAccP500700(0),
fHistd0rphiITSMISPDInAccP10001500(0),
fHistd0rphiITSMISPDInAccS150200(0),
fHistd0rphiITSMISPDInAccS500700(0),
fHistd0rphiITSMISPDInAccS10001500(0),
fHistd0rphiITSMIoneSPDInAccP150200(0),
fHistd0rphiITSMIoneSPDInAccP350450(0),
fHistd0rphiITSMIoneSPDInAccP500700(0),
fHistd0rphiITSMIoneSPDInAccP10001500(0),
fHistd0rphiITSMIoneSPDInAccP25004000(0),
fHistd0rphiITSMIoneSPDInAccP40008000(0),
fHistd0zITSMIoneSPDInAccP150200(0),
fHistd0zITSMIoneSPDInAccP500700(0),
fHistd0zITSMIoneSPDInAccP10001500(0),
fHistd0zVSetaTPCInAccP10001500(0),
fHistd0rphiVSphiITSMIoneSPDInAccP10001500(0),
fHistd0rphiVSetaITSMIoneSPDInAccP10001500(0),
fHistd0rphiITSMIoneSPDInAccS150200(0),
fHistd0rphiITSMIoneSPDInAccS350450(0),
fHistd0rphiITSMIoneSPDInAccS500700(0),
fHistd0rphiITSMIoneSPDInAccS500700from22(0),
fHistd0rphiITSMIoneSPDInAccS500700from211(0),
fHistd0rphiITSMIoneSPDInAccS500700from310(0),
fHistd0rphiITSMIoneSPDInAccS500700from321(0),
fHistd0rphiITSMIoneSPDInAccS500700from3122(0),
fHistd0rphiITSMIoneSPDInAccS10001500(0),
fHistd0rphiITSMIoneSPDInAccS25004000(0),
fHistd0rphiITSMIoneSPDInAccS40008000(0),
fHistd0rphiITSMIoneSPDInAccS150200fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS150200fromMat(0),
fHistd0rphiITSMIoneSPDInAccS350450fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS350450fromMat(0),
fHistd0rphiITSMIoneSPDInAccS500700fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS500700fromMat(0),
fHistd0rphiITSMIoneSPDInAccS10001500fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS10001500fromMat(0),
fHistd0rphiITSMIoneSPDInAccS25004000fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS25004000fromMat(0),
fHistd0rphiITSMIoneSPDInAccS40008000fromStrange(0),
fHistd0rphiITSMIoneSPDInAccS40008000fromMat(0),
fHistd0zITSMIoneSPDInAccS150200(0),
fHistd0zITSMIoneSPDInAccS500700(0),
fHistd0zITSMIoneSPDInAccS10001500(0),
fHistPDGMoth(0),
fHistPDGMoth150200(0),
fHistPDGMoth500700(0),
fHistPDGMoth10001500(0),
fHistPDGTrk(0),
fHistITSRedChi2NonFakePt02(0),
fHistITSRedChi2FakePt02(0),
fHistITSRedChi2NonFakePt05(0),
fHistITSRedChi2FakePt05(0),
fHistITSRedChi2NonFakePt1(0),
fHistITSRedChi2FakePt1(0),
fNtupleESDTracks(0),
fNtupleITSAlignExtra(0),
fNtupleITSAlignSPDTracklets(0),
fESDtrackCutsTPC(0),
fESDtrackCutsITSTPC(0)
{
  // Constructor

  for(Int_t i=0; i<11; i++) fCountsPerPtBin[i]=0;

  // Define input and output slots here
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());  //My private output
}
//________________________________________________________________________
AliAnalysisTaskITSTrackingCheck::~AliAnalysisTaskITSTrackingCheck()
{
  // Destructor

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  gStyle->SetHistLineWidth(2);

  Int_t nPtBins=34;
  Float_t xPtBinsDefault[35]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.5,0.6,0.7,0.8,1.0,1.5,2.,2.5,3,4,5,6,8,10,15,20,25,30};
  Float_t xPtBinsForBox[35]={0,0.025,0.05,0.075,0.08,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.280,0.325,0.35,0.375,0.4,0.45,0.55,0.7,0.75,0.9,1.2,1.7,2.5,3.5,4.5,5,6,8,13,15,23,28,35};
  Float_t xPtBins[35];
  for(Int_t ib=0; ib<=nPtBins; ib++) {
    if(fUsePtBinsForBox) {
      xPtBins[ib]=xPtBinsForBox[ib];
    } else {
      xPtBins[ib]=xPtBinsDefault[ib];
    }
  }


  for(Int_t i=0; i<11; i++) fCountsPerPtBin[i]=0;

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();

  fHistPDGMoth = new TH1F("fHistPDGMoth","",4000,0,4000);
  fOutput->Add(fHistPDGMoth);
  fHistPDGMoth150200 = new TH1F("fHistPDGMoth150200","",4000,0,4000);
  fOutput->Add(fHistPDGMoth150200);
  fHistPDGMoth500700 = new TH1F("fHistPDGMoth500700","",4000,0,4000);
  fOutput->Add(fHistPDGMoth500700);
  fHistPDGMoth10001500 = new TH1F("fHistPDGMoth10001500","",4000,0,4000);
  fOutput->Add(fHistPDGMoth10001500);
  fHistPDGTrk = new TH1F("fHistPDGTrk","",4000,0,4000);
  fOutput->Add(fHistPDGTrk);

  fHistNEvents = new TH1F("fHistNEvents", "Events: -1 tot, 0 sel, 1 vSPD3D, 2 vSPDZ, 3 vSPD, 4 vTRK; Type; N Events",12, -1.5, 10.5);
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistNEventsFrac = (TH1F*)fHistNEvents->Clone("fHistNEventsFrac");
  fOutput->Add(fHistNEventsFrac);

  fHistNtracks = new TH1F("fHistNtracks", "N ESD tracks; N tracks; Events",5000, -0.5, 4999.5);
  fHistNtracks->Sumw2();
  fHistNtracks->SetMinimum(0);
  fOutput->Add(fHistNtracks);

  fHistNclsITSMI = new TH1F("fHistNclsITSMI", "N ITS clusters per track (MI); N clusters; Counts",8, -0.5, 7.5);
  fHistNclsITSMI->Sumw2();
  fHistNclsITSMI->SetMinimum(0);
  fOutput->Add(fHistNclsITSMI);

  fHistNclsITSSAInAcc = new TH1F("fHistNclsITSSAInAcc", "N ITS clusters per track (SA); N clusters; Counts",7, -0.5, 6.5);
  fHistNclsITSSAInAcc->Sumw2();
  fHistNclsITSSAInAcc->SetMinimum(0);
  fOutput->Add(fHistNclsITSSAInAcc);  

  fHistNclsITSSA = new TH1F("fHistNclsITSSA", "N ITS clusters per track (SA); N clusters; Counts",8, -0.5, 7.5);
  fHistNclsITSSA->Sumw2();
  fHistNclsITSSA->SetMinimum(0);
  fOutput->Add(fHistNclsITSSA);  

  fHistClusterMapITSMI = new TH1F("fHistClusterMapITSMI", "N tracks with point on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMI->Sumw2();
  fHistClusterMapITSMI->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMI);
  
  fHistClusterMapITSSA = new TH1F("fHistClusterMapITSSA", "N tracks with point on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSA->Sumw2();
  fHistClusterMapITSSA->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSA);

  fHistClusterMapITSSAInAcc = new TH1F("fHistClusterMapITSSAInAcc", "N tracks with point on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAInAcc->Sumw2();
  fHistClusterMapITSSAInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAInAcc);

  fHistClusterMapITSMIok = new TH1F("fHistClusterMapITSMIok", "N tracks with ok on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMIok->Sumw2();
  fHistClusterMapITSMIok->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIok);
  
  

  fHistClusterMapITSSAokInAcc = new TH1F("fHistClusterMapITSSAokInAcc", "N tracks with ok on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAokInAcc->Sumw2();
  fHistClusterMapITSSAokInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAokInAcc);

  fHistClusterMapModuleITSSAokInAcc = new TH1F("fHistClusterMapModuleITSSAokInAcc", "N tracks with ok on Module (SA); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSSAokInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSSAokInAcc);

  fHistClusterMapModuleITSMIokInAcc = new TH1F("fHistClusterMapModuleITSMIokInAcc", "N tracks with ok on Module (MI); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSMIokInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSMIokInAcc);

  fHistNClustersMapModule = new TH1F("fHistNClustersMapModule", "N clusters; Module; N tracks",2198, -0.5, 2197.5);
  fHistNClustersMapModule->SetMinimum(0);
  fOutput->Add(fHistNClustersMapModule);

  fHistClusterMapITSSAok = new TH1F("fHistClusterMapITSSAok", "N tracks with ok on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAok->Sumw2();
  fHistClusterMapITSSAok->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAok);

  fHistClusterMapITSMIbad = new TH1F("fHistClusterMapITSMIbad", "N tracks with bad on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMIbad->Sumw2();
  fHistClusterMapITSMIbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIbad);
  
  fHistClusterMapITSSAbadInAcc = new TH1F("fHistClusterMapITSSAbadInAcc", "N tracks with bad on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAbadInAcc->Sumw2();
  fHistClusterMapITSSAbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAbadInAcc);

  fHistClusterMapModuleITSSAbadInAcc = new TH1F("fHistClusterMapModuleITSSAbadInAcc", "N tracks with bad on Module (SA); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSSAbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSSAbadInAcc);

  fHistClusterMapModuleITSMIbadInAcc = new TH1F("fHistClusterMapModuleITSMIbadInAcc", "N tracks with bad on Module (MI); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSMIbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSMIbadInAcc);

  fHistClusterMapITSSAbad = new TH1F("fHistClusterMapITSSAbad", "N tracks with bad on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAbad->Sumw2();
  fHistClusterMapITSSAbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAbad);

  fHistClusterMapITSMIskipped = new TH1F("fHistClusterMapITSMIskipped", "N tracks with skip on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMIskipped->Sumw2();
  fHistClusterMapITSMIskipped->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIskipped);
  
  fHistClusterMapITSSAskippedInAcc = new TH1F("fHistClusterMapITSSAskippedInAcc", "N tracks with skip on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAskippedInAcc->Sumw2();
  fHistClusterMapITSSAskippedInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAskippedInAcc);

  fHistClusterMapITSSAskipped = new TH1F("fHistClusterMapITSSAskipped", "N tracks with skip on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAskipped->Sumw2();
  fHistClusterMapITSSAskipped->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAskipped);

  fHistClusterMapITSMIoutinz = new TH1F("fHistClusterMapITSMIoutinz", "N tracks out in z on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMIoutinz->Sumw2();
  fHistClusterMapITSMIoutinz->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIoutinz);
  
  fHistClusterMapITSSAoutinzInAcc = new TH1F("fHistClusterMapITSSAoutinzInAcc", "N tracks with out in z on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAoutinzInAcc->Sumw2();
  fHistClusterMapITSSAoutinzInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAoutinzInAcc);

  fHistClusterMapITSSAoutinz = new TH1F("fHistClusterMapITSSAoutinz", "N tracks with out in z on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAoutinz->Sumw2();
  fHistClusterMapITSSAoutinz->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAoutinz);

  fHistClusterMapITSSAokoutinzbad = new TH1F("fHistClusterMapITSSAokoutinzbad", "N tracks with cluster or bad zone or out in z (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAokoutinzbad->Sumw2();
  fHistClusterMapITSSAokoutinzbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAokoutinzbad);

  fHistClusterMapITSMIokoutinzbad = new TH1F("fHistClusterMapITSMIokoutinzbad", "N tracks with cluster or bad zone or out in z (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMIokoutinzbad->Sumw2();
  fHistClusterMapITSMIokoutinzbad->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMIokoutinzbad);

  fHistClusterMapITSSAokoutinzbadInAcc = new TH1F("fHistClusterMapITSSAokoutinzbadInAcc", "N tracks with cluster or bad zone or out in z (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAokoutinzbadInAcc->Sumw2();
  fHistClusterMapITSSAokoutinzbadInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAokoutinzbadInAcc);

  fHistClusterMapITSMInorefit = new TH1F("fHistClusterMapITSMInorefit", "N tracks with norefit on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMInorefit->Sumw2();
  fHistClusterMapITSMInorefit->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMInorefit);
  
  fHistClusterMapITSSAnorefitInAcc = new TH1F("fHistClusterMapITSSAnorefitInAcc", "N tracks with norefit on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAnorefitInAcc->Sumw2();
  fHistClusterMapITSSAnorefitInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnorefitInAcc);

  fHistClusterMapITSSAnorefit = new TH1F("fHistClusterMapITSSAnorefit", "N tracks with norefit on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAnorefit->Sumw2();
  fHistClusterMapITSSAnorefit->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnorefit);

  fHistClusterMapITSMInocls = new TH1F("fHistClusterMapITSMInocls", "N tracks with nocls on Layer (MI); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSMInocls->Sumw2();
  fHistClusterMapITSMInocls->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSMInocls);
  
  fHistClusterMapITSSAnoclsInAcc = new TH1F("fHistClusterMapITSSAnoclsInAcc", "N tracks with nocls on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAnoclsInAcc->Sumw2();
  fHistClusterMapITSSAnoclsInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnoclsInAcc);
  
  fHistClusterMapModuleITSSAnoclsInAcc = new TH1F("fHistClusterMapModuleITSSAnoclsInAcc", "N tracks with nocls on Module (SA); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSSAnoclsInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSSAnoclsInAcc);

  fHistClusterMapModuleITSMInoclsInAcc = new TH1F("fHistClusterMapModuleITSMInoclsInAcc", "N tracks with nocls on Module (MI); Module; N tracks",2198, -0.5, 2197.5);
  fHistClusterMapModuleITSMInoclsInAcc->SetMinimum(0);
  fOutput->Add(fHistClusterMapModuleITSMInoclsInAcc);

  fHistClusterMapITSSAnocls = new TH1F("fHistClusterMapITSSAnocls", "N tracks with nocls on Layer (SA); Layer; N tracks",7, -0.5, 6.5);
  fHistClusterMapITSSAnocls->Sumw2();
  fHistClusterMapITSSAnocls->SetMinimum(0);
  fOutput->Add(fHistClusterMapITSSAnocls);
  
  
  fHistZatSPDouter0ok = new TH1F("fHistZatSPDouter0ok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter0ok->Sumw2();
  fHistZatSPDouter0ok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter0ok);
  fHistZatSPDouter1ok = new TH1F("fHistZatSPDouter1ok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter1ok->Sumw2();
  fHistZatSPDouter1ok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter1ok);
  fHistZatSPDouter2ok = new TH1F("fHistZatSPDouter2ok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter2ok->Sumw2();
  fHistZatSPDouter2ok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter2ok);
  fHistZatSPDouter3ok = new TH1F("fHistZatSPDouter3ok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter3ok->Sumw2();
  fHistZatSPDouter3ok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter3ok);
  fHistZatSPDouter0notok = new TH1F("fHistZatSPDouter0notok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter0notok->Sumw2();
  fHistZatSPDouter0notok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter0notok);
  fHistZatSPDouter1notok = new TH1F("fHistZatSPDouter1notok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter1notok->Sumw2();
  fHistZatSPDouter1notok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter1notok);
  fHistZatSPDouter2notok = new TH1F("fHistZatSPDouter2notok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter2notok->Sumw2();
  fHistZatSPDouter2notok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter2notok);
  fHistZatSPDouter3notok = new TH1F("fHistZatSPDouter3notok", "Track z at SPDouter radius; z [cm]; N tracks",3200,-16,16);
  fHistZatSPDouter3notok->Sumw2();
  fHistZatSPDouter3notok->SetMinimum(0);
  fOutput->Add(fHistZatSPDouter3notok);
  
  fHistxlocSDDok = new TH1F("fHistxlocSDDok", "SDD points; xloc [cm]; N tracks",75, -3.75, 3.75);
  fHistxlocSDDok->Sumw2();
  fHistxlocSDDok->SetMinimum(0);
  fOutput->Add(fHistxlocSDDok);

  fHistxlocVSmodSDDok = new TH2F("fHistxlocVSmodSDDok", "SDD points; module; xloc [cm]",260,239.5,499.5,25, -3.75, 3.75);
  fOutput->Add(fHistxlocVSmodSDDok);

  fHistzlocSDDok = new TH1F("fHistzlocSDDok", "SDD points; zloc [cm]; N tracks",77, -3.85, 3.85);
  fHistzlocSDDok->Sumw2();
  fHistzlocSDDok->SetMinimum(0);
  fOutput->Add(fHistzlocSDDok);
  
  fHistxlocSDDall = new TH1F("fHistxlocSDDall", "SDD points; xloc [cm]; N tracks",75, -3.75, 3.75);
  fHistxlocSDDall->Sumw2();
  fHistxlocSDDall->SetMinimum(0);
  fOutput->Add(fHistxlocSDDall);

  fHistzlocSDDall = new TH1F("fHistzlocSDDall", "SDD points; zloc [cm]; N tracks",77, -3.85, 3.85);
  fHistzlocSDDall->Sumw2();
  fHistzlocSDDall->SetMinimum(0);
  fOutput->Add(fHistzlocSDDall);


  fHistPhiTPCInAcc = new TH1F("fHistPhiTPCInAcc","Azimuthal distribution of TPC tracks; #phi; N tracks",100, 0, 2.*3.1415);
  fHistPhiTPCInAcc->Sumw2();
  fHistPhiTPCInAcc->SetMinimum(0);
  fOutput->Add(fHistPhiTPCInAcc);

  fHistEtaTPCInAcc = new TH1F("fHistEtaTPCInAcc","Eta distribution of TPC tracks; #eta; N tracks",100, -1.5, +1.5);
  fHistEtaTPCInAcc->Sumw2();
  fHistEtaTPCInAcc->SetMinimum(0);
  fOutput->Add(fHistEtaTPCInAcc);


  fHistPhiITSMIokbadoutinz6InAcc = new TH1F("fHistPhiITSMIokbadoutinz6InAcc","Azimuthal distribution of ITSMI tracks with 6 layers OK; #phi; N tracks",100,0,2.*3.1415);
  fHistPhiITSMIokbadoutinz6InAcc->Sumw2();
  fHistPhiITSMIokbadoutinz6InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMIokbadoutinz6InAcc);
  
  fHistPtTPC = new TH1F("fHistPtTPC","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPC->Sumw2();
  fHistPtTPC->SetMinimum(0);
  fOutput->Add(fHistPtTPC);
  
  fHistPtITSMI7 = new TH1F("fHistPtITSMI7","pt distribution of ITSMI7 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI7->Sumw2();
  fHistPtITSMI7->SetMinimum(0);
  fOutput->Add(fHistPtITSMI7);
  
  fHistPtITSMI6 = new TH1F("fHistPtITSMI6","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6->Sumw2();
  fHistPtITSMI6->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6);
  
  fHistPtITSMI5 = new TH1F("fHistPtITSMI5","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5->Sumw2();
  fHistPtITSMI5->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5);
  
  fHistPtITSMI4 = new TH1F("fHistPtITSMI4","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4->Sumw2();
  fHistPtITSMI4->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4);
  
  fHistPtITSMI3 = new TH1F("fHistPtITSMI3","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3->Sumw2();
  fHistPtITSMI3->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3);
  
  fHistPtITSMI2 = new TH1F("fHistPtITSMI2","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2->Sumw2();
  fHistPtITSMI2->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2);
  
  fHistPtITSMISPD = new TH1F("fHistPtITSMISPD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPD->Sumw2();
  fHistPtITSMISPD->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPD);

  fHistPtITSMIoneSPD = new TH1F("fHistPtITSMIoneSPD","pt distribution of ITSMIoneSPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPD->Sumw2();
  fHistPtITSMIoneSPD->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPD);

  fHistPtITSMItwoSPD = new TH1F("fHistPtITSMItwoSPD","pt distribution of ITSMItwoSPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMItwoSPD->Sumw2();
  fHistPtITSMItwoSPD->SetMinimum(0);
  fOutput->Add(fHistPtITSMItwoSPD);

  fHistPtTPCInAcc = new TH1F("fHistPtTPCInAcc","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAcc->Sumw2();
  fHistPtTPCInAcc->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAcc);

  fHistPtTPCInAccTOFbc0 = new TH1F("fHistPtTPCInAccTOFbc0","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccTOFbc0->Sumw2();
  fHistPtTPCInAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccTOFbc0);

  fHistPtTPCInAccwSDD = new TH1F("fHistPtTPCInAccwSDD","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccwSDD->Sumw2();
  fHistPtTPCInAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccwSDD);

  fHistPtTPCInAccTOFbc0wSDD = new TH1F("fHistPtTPCInAccTOFbc0wSDD","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccTOFbc0wSDD->Sumw2();
  fHistPtTPCInAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccTOFbc0wSDD);
  
  fHistPtTPCInAccwoSDD = new TH1F("fHistPtTPCInAccwoSDD","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccwoSDD->Sumw2();
  fHistPtTPCInAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccwoSDD);

  fHistPtTPCInAccTOFbc0woSDD = new TH1F("fHistPtTPCInAccTOFbc0woSDD","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccTOFbc0woSDD->Sumw2();
  fHistPtTPCInAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccTOFbc0woSDD);
  
  fHistPtTPCInAccMCtwoSPD = new TH1F("fHistPtTPCInAccMCtwoSPD","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccMCtwoSPD->Sumw2();
  fHistPtTPCInAccMCtwoSPD->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccMCtwoSPD);
  
  fHistPtTPCInAccMConeSPD = new TH1F("fHistPtTPCInAccMConeSPD","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccMConeSPD->Sumw2();
  fHistPtTPCInAccMConeSPD->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccMConeSPD);
  
  fHistdEdxVSPtTPCInAcc = new TH2F("fHistdEdxVSPtTPCInAcc","dE/dx distribution of TPC tracks; p_{t} [GeV/c]; TPC dE/dx",100,0,1,50,0,500);
  fOutput->Add(fHistdEdxVSPtTPCInAcc);
  
  fHistdEdxVSPtITSTPCsel = new TH2F("fHistdEdxVSPtITSTPCsel","dE/dx distribution of TPC tracks; p_{t} [GeV/c]; ITS dE/dx",100,0,1,50,0,500);
  fOutput->Add(fHistdEdxVSPtITSTPCsel);
  
  fHistTPCclsVSPtTPCInAcc = new TH2F("fHistTPCclsVSPtTPCInAcc","TPC ncls; p_{t} [GeV/c]; nclsTPC",50,0,10,80,0,160);
  fOutput->Add(fHistTPCclsVSPtTPCInAcc);
  fHistTPCclsVSPtITSMISPDInAcc = new TH2F("fHistTPCclsVSPtITSMISPDInAcc","TPC ncls; p_{t} [GeV/c]; nclsTPC",50,0,10,80,0,160);
  fOutput->Add(fHistTPCclsVSPtITSMISPDInAcc);

  fHistPtVSphiTPCInAcc = new TH2F("fHistPtVSphiTPCInAcc","pt distribution of TPC tracks; phi; p_{t} [GeV/c]",18,0,6.28,20,0,0.5);
  fHistPtVSphiTPCInAcc->SetMinimum(0);
  fOutput->Add(fHistPtVSphiTPCInAcc);
  
  fHistPtTPCInAccNoTRDout = new TH1F("fHistPtTPCInAccNoTRDout","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccNoTRDout->Sumw2();
  fHistPtTPCInAccNoTRDout->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccNoTRDout);
  
  fHistPtTPCInAccNoTOFout = new TH1F("fHistPtTPCInAccNoTOFout","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccNoTOFout->Sumw2();
  fHistPtTPCInAccNoTOFout->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccNoTOFout);
  
  fHistPtTPCInAccWithPtTPCAtVtx = new TH1F("fHistPtTPCInAccWithPtTPCAtVtx","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccWithPtTPCAtVtx->Sumw2();
  fHistPtTPCInAccWithPtTPCAtVtx->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccWithPtTPCAtVtx);
  
  fHistPtTPCInAccWithPtTPCAtInnerWall = new TH1F("fHistPtTPCInAccWithPtTPCAtInnerWall","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccWithPtTPCAtInnerWall->Sumw2();
  fHistPtTPCInAccWithPtTPCAtInnerWall->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccWithPtTPCAtInnerWall);
  
  fHistDeltaPtTPC = new TH2F("fHistDeltaPtTPC","pt distribution of TPC tracks; p_{t} [GeV/c]; p_{t} TPC at vtx - p_{t} at inner wall [GeV/c]",100,0,1,50,-0.2,0.2);
  fHistDeltaPtTPC->SetMinimum(0);
  fOutput->Add(fHistDeltaPtTPC);
  
  fHistPtITSMI7InAcc = new TH1F("fHistPtITSMI7InAcc","pt distribution of ITSMI7 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI7InAcc->Sumw2();
  fHistPtITSMI7InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMI7InAcc);
  
  fHistPtITSMI6InAcc = new TH1F("fHistPtITSMI6InAcc","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAcc->Sumw2();
  fHistPtITSMI6InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAcc);
  
  fHistPtITSMI5InAcc = new TH1F("fHistPtITSMI5InAcc","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAcc->Sumw2();
  fHistPtITSMI5InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAcc);
  
  fHistPtITSMI4InAcc = new TH1F("fHistPtITSMI4InAcc","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAcc->Sumw2();
  fHistPtITSMI4InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAcc);
  
  fHistPtITSMI3InAcc = new TH1F("fHistPtITSMI3InAcc","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAcc->Sumw2();
  fHistPtITSMI3InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAcc);
  
  fHistPtITSMI2InAcc = new TH1F("fHistPtITSMI2InAcc","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAcc->Sumw2();
  fHistPtITSMI2InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAcc);
  
  fHistPtITSMISPDInAcc = new TH1F("fHistPtITSMISPDInAcc","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAcc->Sumw2();
  fHistPtITSMISPDInAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAcc);
  
  fHistPtITSMIoneSPDInAcc = new TH1F("fHistPtITSMIoneSPDInAcc","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAcc->Sumw2();
  fHistPtITSMIoneSPDInAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAcc);

  fHistPtITSMItwoSPDInAcc = new TH1F("fHistPtITSMItwoSPDInAcc","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMItwoSPDInAcc->Sumw2();
  fHistPtITSMItwoSPDInAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMItwoSPDInAcc);

  fHistPtITSMI6InAccTOFbc0 = new TH1F("fHistPtITSMI6InAccTOFbc0","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccTOFbc0->Sumw2();
  fHistPtITSMI6InAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccTOFbc0);
  
  fHistPtITSMI5InAccTOFbc0 = new TH1F("fHistPtITSMI5InAccTOFbc0","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccTOFbc0->Sumw2();
  fHistPtITSMI5InAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccTOFbc0);
  
  fHistPtITSMI4InAccTOFbc0 = new TH1F("fHistPtITSMI4InAccTOFbc0","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccTOFbc0->Sumw2();
  fHistPtITSMI4InAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccTOFbc0);
  
  fHistPtITSMI3InAccTOFbc0 = new TH1F("fHistPtITSMI3InAccTOFbc0","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccTOFbc0->Sumw2();
  fHistPtITSMI3InAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccTOFbc0);
  
  fHistPtITSMI2InAccTOFbc0 = new TH1F("fHistPtITSMI2InAccTOFbc0","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccTOFbc0->Sumw2();
  fHistPtITSMI2InAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccTOFbc0);
  
  fHistPtITSMISPDInAccTOFbc0 = new TH1F("fHistPtITSMISPDInAccTOFbc0","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccTOFbc0->Sumw2();
  fHistPtITSMISPDInAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccTOFbc0);
  
  fHistPtITSMIoneSPDInAccTOFbc0 = new TH1F("fHistPtITSMIoneSPDInAccTOFbc0","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccTOFbc0->Sumw2();
  fHistPtITSMIoneSPDInAccTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccTOFbc0);

  fHistPtITSMI6InAccwSDD = new TH1F("fHistPtITSMI6InAccwSDD","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccwSDD->Sumw2();
  fHistPtITSMI6InAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccwSDD);
  
  fHistPtITSMI5InAccwSDD = new TH1F("fHistPtITSMI5InAccwSDD","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccwSDD->Sumw2();
  fHistPtITSMI5InAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccwSDD);
  
  fHistPtITSMI4InAccwSDD = new TH1F("fHistPtITSMI4InAccwSDD","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccwSDD->Sumw2();
  fHistPtITSMI4InAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccwSDD);
  
  fHistPtITSMI3InAccwSDD = new TH1F("fHistPtITSMI3InAccwSDD","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccwSDD->Sumw2();
  fHistPtITSMI3InAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccwSDD);
  
  fHistPtITSMI2InAccwSDD = new TH1F("fHistPtITSMI2InAccwSDD","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccwSDD->Sumw2();
  fHistPtITSMI2InAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccwSDD);
  
  fHistPtITSMISPDInAccwSDD = new TH1F("fHistPtITSMISPDInAccwSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccwSDD->Sumw2();
  fHistPtITSMISPDInAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccwSDD);
  
  fHistPtITSMIoneSPDInAccwSDD = new TH1F("fHistPtITSMIoneSPDInAccwSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccwSDD->Sumw2();
  fHistPtITSMIoneSPDInAccwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccwSDD);

  fHistPtITSMI6InAccTOFbc0wSDD = new TH1F("fHistPtITSMI6InAccTOFbc0wSDD","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccTOFbc0wSDD->Sumw2();
  fHistPtITSMI6InAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccTOFbc0wSDD);
  
  fHistPtITSMI5InAccTOFbc0wSDD = new TH1F("fHistPtITSMI5InAccTOFbc0wSDD","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccTOFbc0wSDD->Sumw2();
  fHistPtITSMI5InAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccTOFbc0wSDD);
  
  fHistPtITSMI4InAccTOFbc0wSDD = new TH1F("fHistPtITSMI4InAccTOFbc0wSDD","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccTOFbc0wSDD->Sumw2();
  fHistPtITSMI4InAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccTOFbc0wSDD);
  
  fHistPtITSMI3InAccTOFbc0wSDD = new TH1F("fHistPtITSMI3InAccTOFbc0wSDD","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccTOFbc0wSDD->Sumw2();
  fHistPtITSMI3InAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccTOFbc0wSDD);
  
  fHistPtITSMI2InAccTOFbc0wSDD = new TH1F("fHistPtITSMI2InAccTOFbc0wSDD","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccTOFbc0wSDD->Sumw2();
  fHistPtITSMI2InAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccTOFbc0wSDD);
  
  fHistPtITSMISPDInAccTOFbc0wSDD = new TH1F("fHistPtITSMISPDInAccTOFbc0wSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccTOFbc0wSDD->Sumw2();
  fHistPtITSMISPDInAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccTOFbc0wSDD);
  
  fHistPtITSMIoneSPDInAccTOFbc0wSDD = new TH1F("fHistPtITSMIoneSPDInAccTOFbc0wSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccTOFbc0wSDD->Sumw2();
  fHistPtITSMIoneSPDInAccTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccTOFbc0wSDD);

  fHistPtITSMI6InAccwoSDD = new TH1F("fHistPtITSMI6InAccwoSDD","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccwoSDD->Sumw2();
  fHistPtITSMI6InAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccwoSDD);
  
  fHistPtITSMI5InAccwoSDD = new TH1F("fHistPtITSMI5InAccwoSDD","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccwoSDD->Sumw2();
  fHistPtITSMI5InAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccwoSDD);
  
  fHistPtITSMI4InAccwoSDD = new TH1F("fHistPtITSMI4InAccwoSDD","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccwoSDD->Sumw2();
  fHistPtITSMI4InAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccwoSDD);
  
  fHistPtITSMI3InAccwoSDD = new TH1F("fHistPtITSMI3InAccwoSDD","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccwoSDD->Sumw2();
  fHistPtITSMI3InAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccwoSDD);
  
  fHistPtITSMI2InAccwoSDD = new TH1F("fHistPtITSMI2InAccwoSDD","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccwoSDD->Sumw2();
  fHistPtITSMI2InAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccwoSDD);
  
  fHistPtITSMISPDInAccwoSDD = new TH1F("fHistPtITSMISPDInAccwoSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccwoSDD->Sumw2();
  fHistPtITSMISPDInAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccwoSDD);
  
  fHistPtITSMIoneSPDInAccwoSDD = new TH1F("fHistPtITSMIoneSPDInAccwoSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccwoSDD->Sumw2();
  fHistPtITSMIoneSPDInAccwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccwoSDD);

  fHistPtITSMI6InAccTOFbc0woSDD = new TH1F("fHistPtITSMI6InAccTOFbc0woSDD","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccTOFbc0woSDD->Sumw2();
  fHistPtITSMI6InAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccTOFbc0woSDD);
  
  fHistPtITSMI5InAccTOFbc0woSDD = new TH1F("fHistPtITSMI5InAccTOFbc0woSDD","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccTOFbc0woSDD->Sumw2();
  fHistPtITSMI5InAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccTOFbc0woSDD);
  
  fHistPtITSMI4InAccTOFbc0woSDD = new TH1F("fHistPtITSMI4InAccTOFbc0woSDD","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccTOFbc0woSDD->Sumw2();
  fHistPtITSMI4InAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccTOFbc0woSDD);
  
  fHistPtITSMI3InAccTOFbc0woSDD = new TH1F("fHistPtITSMI3InAccTOFbc0woSDD","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccTOFbc0woSDD->Sumw2();
  fHistPtITSMI3InAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccTOFbc0woSDD);
  
  fHistPtITSMI2InAccTOFbc0woSDD = new TH1F("fHistPtITSMI2InAccTOFbc0woSDD","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccTOFbc0woSDD->Sumw2();
  fHistPtITSMI2InAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccTOFbc0woSDD);
  
  fHistPtITSMISPDInAccTOFbc0woSDD = new TH1F("fHistPtITSMISPDInAccTOFbc0woSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccTOFbc0woSDD->Sumw2();
  fHistPtITSMISPDInAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccTOFbc0woSDD);
  
  fHistPtITSMIoneSPDInAccTOFbc0woSDD = new TH1F("fHistPtITSMIoneSPDInAccTOFbc0woSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccTOFbc0woSDD->Sumw2();
  fHistPtITSMIoneSPDInAccTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccTOFbc0woSDD);

  fHistPtITSMIoneSPDInAccShared = new TH1F("fHistPtITSMIoneSPDInAccShared","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccShared->Sumw2();
  fHistPtITSMIoneSPDInAccShared->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccShared);

  fHistPtITSMIoneSPDInAccSharedSPD = new TH1F("fHistPtITSMIoneSPDInAccSharedSPD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccSharedSPD->Sumw2();
  fHistPtITSMIoneSPDInAccSharedSPD->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccSharedSPD);

  fHistPtITSMISPD1InAccShared = new TH1F("fHistPtITSMISPD1InAccShared","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPD1InAccShared->Sumw2();
  fHistPtITSMISPD1InAccShared->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPD1InAccShared);

  fHistPtITSMISPD2InAccShared = new TH1F("fHistPtITSMISPD2InAccShared","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPD2InAccShared->Sumw2();
  fHistPtITSMISPD2InAccShared->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPD2InAccShared);

  fHistPtITSMIoneSPDInAccSharedFake = new TH1F("fHistPtITSMIoneSPDInAccSharedFake","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccSharedFake->Sumw2();
  fHistPtITSMIoneSPDInAccSharedFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccSharedFake);

  fHistPtITSMIoneSPDInAccSharedSPDFake = new TH1F("fHistPtITSMIoneSPDInAccSharedSPDFake","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccSharedSPDFake->Sumw2();
  fHistPtITSMIoneSPDInAccSharedSPDFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccSharedSPDFake);

  fHistPhiITSMI7InAcc = new TH1F("fHistPhiITSMI7InAcc","phi distribution of ITSMI7 tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMI7InAcc->Sumw2();
  fHistPhiITSMI7InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMI7InAcc);

  fHistPhiITSMI6InAcc = new TH1F("fHistPhiITSMI6InAcc","phi distribution of ITSMI6 tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMI6InAcc->Sumw2();
  fHistPhiITSMI6InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMI6InAcc);
  
  fHistPhiITSMI5InAcc = new TH1F("fHistPhiITSMI5InAcc","phi distribution of ITSMI5 tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMI5InAcc->Sumw2();
  fHistPhiITSMI5InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMI5InAcc);
  
  fHistPhiITSMI4InAcc = new TH1F("fHistPhiITSMI4InAcc","phi distribution of ITSMI4 tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMI4InAcc->Sumw2();
  fHistPhiITSMI4InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMI4InAcc);
  
  fHistPhiITSMI3InAcc = new TH1F("fHistPhiITSMI3InAcc","phi distribution of ITSMI3 tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMI3InAcc->Sumw2();
  fHistPhiITSMI3InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMI3InAcc);
  
  fHistPhiITSMI2InAcc = new TH1F("fHistPhiITSMI2InAcc","phi distribution of ITSMI2 tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMI2InAcc->Sumw2();
  fHistPhiITSMI2InAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMI2InAcc);
  
  fHistPhiITSMISPDInAcc = new TH1F("fHistPhiITSMISPDInAcc","phi distribution of ITSMISPD tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMISPDInAcc->Sumw2();
  fHistPhiITSMISPDInAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMISPDInAcc);
  
  fHistPhiITSMIoneSPDInAcc = new TH1F("fHistPhiITSMIoneSPDInAcc","phi distribution of ITSMISPD tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMIoneSPDInAcc->Sumw2();
  fHistPhiITSMIoneSPDInAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMIoneSPDInAcc);

  fHistPhiITSMItwoSPDInAcc = new TH1F("fHistPhiITSMItwoSPDInAcc","phi distribution of ITSMISPD tracks; #phi; N tracks",100,0,2*3.1415);
  fHistPhiITSMItwoSPDInAcc->Sumw2();
  fHistPhiITSMItwoSPDInAcc->SetMinimum(0);
  fOutput->Add(fHistPhiITSMItwoSPDInAcc);

  fHistEtaITSMI6InAcc = new TH1F("fHistEtaITSMI6InAcc","eta distribution of ITSMI6 tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMI6InAcc->Sumw2();
  fHistEtaITSMI6InAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMI6InAcc);

  fHistEtaITSMI7InAcc = new TH1F("fHistEtaITSMI7InAcc","eta distribution of ITSMI7 tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMI7InAcc->Sumw2();
  fHistEtaITSMI7InAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMI7InAcc);
  
  fHistEtaITSMI5InAcc = new TH1F("fHistEtaITSMI5InAcc","eta distribution of ITSMI5 tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMI5InAcc->Sumw2();
  fHistEtaITSMI5InAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMI5InAcc);
  
  fHistEtaITSMI4InAcc = new TH1F("fHistEtaITSMI4InAcc","eta distribution of ITSMI4 tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMI4InAcc->Sumw2();
  fHistEtaITSMI4InAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMI4InAcc);
  
  fHistEtaITSMI3InAcc = new TH1F("fHistEtaITSMI3InAcc","eta distribution of ITSMI3 tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMI3InAcc->Sumw2();
  fHistEtaITSMI3InAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMI3InAcc);
  
  fHistEtaITSMI2InAcc = new TH1F("fHistEtaITSMI2InAcc","eta distribution of ITSMI2 tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMI2InAcc->Sumw2();
  fHistEtaITSMI2InAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMI2InAcc);
  
  fHistEtaITSMISPDInAcc = new TH1F("fHistEtaITSMISPDInAcc","eta distribution of ITSMISPD tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMISPDInAcc->Sumw2();
  fHistEtaITSMISPDInAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMISPDInAcc);
  
  fHistEtaITSMIoneSPDInAcc = new TH1F("fHistEtaITSMIoneSPDInAcc","eta distribution of ITSMISPD tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMIoneSPDInAcc->Sumw2();
  fHistEtaITSMIoneSPDInAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMIoneSPDInAcc);

  fHistEtaITSMItwoSPDInAcc = new TH1F("fHistEtaITSMItwoSPDInAcc","eta distribution of ITSMISPD tracks; #eta; N tracks",100,-1.5,+1.5);
  fHistEtaITSMItwoSPDInAcc->Sumw2();
  fHistEtaITSMItwoSPDInAcc->SetMinimum(0);
  fOutput->Add(fHistEtaITSMItwoSPDInAcc);


  fHistPtITSMI7InAccFake = new TH1F("fHistPtITSMI7InAccFake","pt distribution of ITSMI7 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI7InAccFake->Sumw2();
  fHistPtITSMI7InAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMI7InAccFake);
  
  fHistPtITSMI6InAccFake = new TH1F("fHistPtITSMI6InAccFake","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccFake->Sumw2();
  fHistPtITSMI6InAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccFake);
  
  fHistPtITSMI5InAccFake = new TH1F("fHistPtITSMI5InAccFake","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccFake->Sumw2();
  fHistPtITSMI5InAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccFake);
  
  fHistPtITSMI4InAccFake = new TH1F("fHistPtITSMI4InAccFake","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccFake->Sumw2();
  fHistPtITSMI4InAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccFake);
  
  fHistPtITSMI3InAccFake = new TH1F("fHistPtITSMI3InAccFake","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccFake->Sumw2();
  fHistPtITSMI3InAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccFake);
  
  fHistPtITSMI2InAccFake = new TH1F("fHistPtITSMI2InAccFake","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccFake->Sumw2();
  fHistPtITSMI2InAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccFake);
  
  fHistPtITSMISPDInAccFake = new TH1F("fHistPtITSMISPDInAccFake","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccFake->Sumw2();
  fHistPtITSMISPDInAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccFake);
  
  fHistPtITSMIoneSPDInAccFake = new TH1F("fHistPtITSMIoneSPDInAccFake","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccFake->Sumw2();
  fHistPtITSMIoneSPDInAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccFake);

  fHistPtITSMItwoSPDInAccFake = new TH1F("fHistPtITSMItwoSPDInAccFake","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMItwoSPDInAccFake->Sumw2();
  fHistPtITSMItwoSPDInAccFake->SetMinimum(0);
  fOutput->Add(fHistPtITSMItwoSPDInAccFake);

  fHistPtITSMIoneSPDthreeSDDSSDInAcc = new TH1F("fHistPtITSMIoneSPDthreeSDDSSDInAcc","pt distribution of ITSMI tracks (>0 in SPD, >2 in SDD+SSD); p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDthreeSDDSSDInAcc->Sumw2();
  fHistPtITSMIoneSPDthreeSDDSSDInAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDthreeSDDSSDInAcc);

  fHistPtITSTPCsel = new TH1F("fHistPtITSTPCsel","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCsel->Sumw2();
  fHistPtITSTPCsel->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCsel);

  fHistPtITSTPCselTOFbc0 = new TH1F("fHistPtITSTPCselTOFbc0","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselTOFbc0->Sumw2();
  fHistPtITSTPCselTOFbc0->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselTOFbc0);

  fHistPtITSTPCselwSDD = new TH1F("fHistPtITSTPCselwSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselwSDD->Sumw2();
  fHistPtITSTPCselwSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselwSDD);

  fHistPtITSTPCselTOFbc0wSDD = new TH1F("fHistPtITSTPCselTOFbc0wSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselTOFbc0wSDD->Sumw2();
  fHistPtITSTPCselTOFbc0wSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselTOFbc0wSDD);

  fHistPtITSTPCselwoSDD = new TH1F("fHistPtITSTPCselwoSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselwoSDD->Sumw2();
  fHistPtITSTPCselwoSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselwoSDD);

  fHistPtITSTPCselTOFbc0woSDD = new TH1F("fHistPtITSTPCselTOFbc0woSDD","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselTOFbc0woSDD->Sumw2();
  fHistPtITSTPCselTOFbc0woSDD->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselTOFbc0woSDD);

  fHistPtITSTPCselP = new TH1F("fHistPtITSTPCselP","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselP->Sumw2();
  fHistPtITSTPCselP->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselP);

  fHistPtITSTPCselS = new TH1F("fHistPtITSTPCselS","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselS->Sumw2();
  fHistPtITSTPCselS->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselS);

  fHistPtITSTPCselFake = new TH1F("fHistPtITSTPCselFake","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselFake->Sumw2();
  fHistPtITSTPCselFake->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselFake);

  fHistPtITSTPCselSfromStrange = new TH1F("fHistPtITSTPCselSfromStrange","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselSfromStrange->Sumw2();
  fHistPtITSTPCselSfromStrange->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselSfromStrange);

  fHistPtITSTPCselPfromStrange = new TH1F("fHistPtITSTPCselPfromStrange","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselPfromStrange->Sumw2();
  fHistPtITSTPCselPfromStrange->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselPfromStrange);

  fHistPtITSTPCselSfromMat = new TH1F("fHistPtITSTPCselSfromMat","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSTPCselSfromMat->Sumw2();
  fHistPtITSTPCselSfromMat->SetMinimum(0);
  fOutput->Add(fHistPtITSTPCselSfromMat);

  fHistPtTPCInAccP = new TH1F("fHistPtTPCInAccP","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccP->Sumw2();
  fHistPtTPCInAccP->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccP);

  fHistPtTPCInAccPfromStrange = new TH1F("fHistPtTPCInAccPfromStrange","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccPfromStrange->Sumw2();
  fHistPtTPCInAccPfromStrange->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccPfromStrange);
  
  fHistPtITSMI7InAccP = new TH1F("fHistPtITSMI7InAccP","pt distribution of ITSMI7 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI7InAccP->Sumw2();
  fHistPtITSMI7InAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMI7InAccP);
  
  fHistPtITSMI6InAccP = new TH1F("fHistPtITSMI6InAccP","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccP->Sumw2();
  fHistPtITSMI6InAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccP);
  
  fHistPtITSMI5InAccP = new TH1F("fHistPtITSMI5InAccP","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccP->Sumw2();
  fHistPtITSMI5InAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccP);
  
  fHistPtITSMI4InAccP = new TH1F("fHistPtITSMI4InAccP","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccP->Sumw2();
  fHistPtITSMI4InAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccP);
  
  fHistPtITSMI3InAccP = new TH1F("fHistPtITSMI3InAccP","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccP->Sumw2();
  fHistPtITSMI3InAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccP);
  
  fHistPtITSMI2InAccP = new TH1F("fHistPtITSMI2InAccP","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccP->Sumw2();
  fHistPtITSMI2InAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccP);
  
  fHistPtITSMISPDInAccP = new TH1F("fHistPtITSMISPDInAccP","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccP->Sumw2();
  fHistPtITSMISPDInAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccP);
  
  fHistPtITSMIoneSPDInAccP = new TH1F("fHistPtITSMIoneSPDInAccP","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccP->Sumw2();
  fHistPtITSMIoneSPDInAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccP);

  fHistPtITSMItwoSPDInAccP = new TH1F("fHistPtITSMItwoSPDInAccP","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMItwoSPDInAccP->Sumw2();
  fHistPtITSMItwoSPDInAccP->SetMinimum(0);
  fOutput->Add(fHistPtITSMItwoSPDInAccP);

  fHistPtTPCInAccS = new TH1F("fHistPtTPCInAccS","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccS->Sumw2();
  fHistPtTPCInAccS->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccS);
  
  fHistPtTPCInAccSfromStrange = new TH1F("fHistPtTPCInAccSfromStrange","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccSfromStrange->Sumw2();
  fHistPtTPCInAccSfromStrange->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccSfromStrange);
  
  fHistPtTPCInAccSfromMat = new TH1F("fHistPtTPCInAccSfromMat","pt distribution of TPC tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtTPCInAccSfromMat->Sumw2();
  fHistPtTPCInAccSfromMat->SetMinimum(0);
  fOutput->Add(fHistPtTPCInAccSfromMat);
  
  fHistPtITSMI7InAccS = new TH1F("fHistPtITSMI7InAccS","pt distribution of ITSMI7 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI7InAccS->Sumw2();
  fHistPtITSMI7InAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMI7InAccS);
  
  fHistPtITSMI6InAccS = new TH1F("fHistPtITSMI6InAccS","pt distribution of ITSMI6 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI6InAccS->Sumw2();
  fHistPtITSMI6InAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMI6InAccS);
  
  fHistPtITSMI5InAccS = new TH1F("fHistPtITSMI5InAccS","pt distribution of ITSMI5 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI5InAccS->Sumw2();
  fHistPtITSMI5InAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMI5InAccS);
  
  fHistPtITSMI4InAccS = new TH1F("fHistPtITSMI4InAccS","pt distribution of ITSMI4 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI4InAccS->Sumw2();
  fHistPtITSMI4InAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMI4InAccS);
  
  fHistPtITSMI3InAccS = new TH1F("fHistPtITSMI3InAccS","pt distribution of ITSMI3 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI3InAccS->Sumw2();
  fHistPtITSMI3InAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMI3InAccS);
  
  fHistPtITSMI2InAccS = new TH1F("fHistPtITSMI2InAccS","pt distribution of ITSMI2 tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMI2InAccS->Sumw2();
  fHistPtITSMI2InAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMI2InAccS);
  
  fHistPtITSMISPDInAccS = new TH1F("fHistPtITSMISPDInAccS","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMISPDInAccS->Sumw2();
  fHistPtITSMISPDInAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMISPDInAccS);
  
  fHistPtITSMIoneSPDInAccS = new TH1F("fHistPtITSMIoneSPDInAccS","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIoneSPDInAccS->Sumw2();
  fHistPtITSMIoneSPDInAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMIoneSPDInAccS);

  fHistPtITSMItwoSPDInAccS = new TH1F("fHistPtITSMItwoSPDInAccS","pt distribution of ITSMISPD tracks; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMItwoSPDInAccS->Sumw2();
  fHistPtITSMItwoSPDInAccS->SetMinimum(0);
  fOutput->Add(fHistPtITSMItwoSPDInAccS);
  
  fHistPtITSMIokbadoutinz6 = new TH1F("fHistPtITSMIokbadoutinz6","pt distribution of ITSMI tracks with 6 layers OK; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIokbadoutinz6->Sumw2();
  fHistPtITSMIokbadoutinz6->SetMinimum(0);
  fOutput->Add(fHistPtITSMIokbadoutinz6);

  fHistPtITSMIokbadoutinz4InAcc = new TH1F("fHistPtITSMIokbadoutinz4InAcc","pt distribution of ITSMI tracks with 4 layers OK; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIokbadoutinz4InAcc->Sumw2();
  fHistPtITSMIokbadoutinz4InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMIokbadoutinz4InAcc);

  fHistPtITSMIokbadoutinz5InAcc = new TH1F("fHistPtITSMIokbadoutinz5InAcc","pt distribution of ITSMI tracks with 5 layers OK; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIokbadoutinz5InAcc->Sumw2();
  fHistPtITSMIokbadoutinz5InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMIokbadoutinz5InAcc);

  fHistPtITSMIokbadoutinz6InAcc = new TH1F("fHistPtITSMIokbadoutinz6InAcc","pt distribution of ITSMI tracks with 6 layers OK; p_{t} [GeV/c]; N tracks",nPtBins,xPtBins);
  fHistPtITSMIokbadoutinz6InAcc->Sumw2();
  fHistPtITSMIokbadoutinz6InAcc->SetMinimum(0);
  fOutput->Add(fHistPtITSMIokbadoutinz6InAcc);

  fHistRProdVtxInAccP = new TH1F("fHistRProdVtxInAccP","Radius of production vertex for primaries; r [cm]; N tracks",100,0,10);
  fHistRProdVtxInAccP->Sumw2();
  fHistRProdVtxInAccP->SetMinimum(0);
  fOutput->Add(fHistRProdVtxInAccP);

  fHistRProdVtxInAccS = new TH1F("fHistRProdVtxInAccS","Radius of production vertex for secondaries; r [cm]; N tracks",100,0,10);
  fHistRProdVtxInAccS->Sumw2();
  fHistRProdVtxInAccS->SetMinimum(0);
  fOutput->Add(fHistRProdVtxInAccS);
  

  fHistd0rphiTPCInAccP150200 = new TH1F("fHistd0rphiTPCInAccP150200","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-5,5);
  fHistd0rphiTPCInAccP150200->Sumw2();
  fHistd0rphiTPCInAccP150200->SetMinimum(0);
  fOutput->Add(fHistd0rphiTPCInAccP150200);

  fHistd0rphiTPCInAccP500700 = new TH1F("fHistd0rphiTPCInAccP500700","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-5,5);
  fHistd0rphiTPCInAccP500700->Sumw2();
  fHistd0rphiTPCInAccP500700->SetMinimum(0);
  fOutput->Add(fHistd0rphiTPCInAccP500700);

  fHistd0rphiTPCInAccP10001500 = new TH1F("fHistd0rphiTPCInAccP10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-5,5);
  fHistd0rphiTPCInAccP10001500->Sumw2();
  fHistd0rphiTPCInAccP10001500->SetMinimum(0);
  fOutput->Add(fHistd0rphiTPCInAccP10001500);

  fHistd0rphiTPCInAccS150200 = new TH1F("fHistd0rphiTPCInAccS150200","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-5,5);
  fHistd0rphiTPCInAccS150200->Sumw2();
  fHistd0rphiTPCInAccS150200->SetMinimum(0);
  fOutput->Add(fHistd0rphiTPCInAccS150200);

  fHistd0rphiTPCInAccS500700 = new TH1F("fHistd0rphiTPCInAccS500700","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-5,5);
  fHistd0rphiTPCInAccS500700->Sumw2();
  fHistd0rphiTPCInAccS500700->SetMinimum(0);
  fOutput->Add(fHistd0rphiTPCInAccS500700);

  fHistd0rphiTPCInAccS10001500 = new TH1F("fHistd0rphiTPCInAccS10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-5,5);
  fHistd0rphiTPCInAccS10001500->Sumw2();
  fHistd0rphiTPCInAccS10001500->SetMinimum(0);
  fOutput->Add(fHistd0rphiTPCInAccS10001500);

  fHistd0rphiITSMISPDInAccP150200 = new TH1F("fHistd0rphiITSMISPDInAccP150200","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMISPDInAccP150200->Sumw2();
  fHistd0rphiITSMISPDInAccP150200->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMISPDInAccP150200);

  fHistd0rphiITSMISPDInAccP500700 = new TH1F("fHistd0rphiITSMISPDInAccP500700","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMISPDInAccP500700->Sumw2();
  fHistd0rphiITSMISPDInAccP500700->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMISPDInAccP500700);

  fHistd0rphiITSMISPDInAccP10001500 = new TH1F("fHistd0rphiITSMISPDInAccP10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMISPDInAccP10001500->Sumw2();
  fHistd0rphiITSMISPDInAccP10001500->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMISPDInAccP10001500);

  fHistd0rphiITSMISPDInAccS150200 = new TH1F("fHistd0rphiITSMISPDInAccS150200","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMISPDInAccS150200->Sumw2();
  fHistd0rphiITSMISPDInAccS150200->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMISPDInAccS150200);

  fHistd0rphiITSMISPDInAccS500700 = new TH1F("fHistd0rphiITSMISPDInAccS500700","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMISPDInAccS500700->Sumw2();
  fHistd0rphiITSMISPDInAccS500700->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMISPDInAccS500700);

  fHistd0rphiITSMISPDInAccS10001500 = new TH1F("fHistd0rphiITSMISPDInAccS10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMISPDInAccS10001500->Sumw2();
  fHistd0rphiITSMISPDInAccS10001500->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMISPDInAccS10001500);

  fHistd0rphiITSMIoneSPDInAccP150200 = new TH1F("fHistd0rphiITSMIoneSPDInAccP150200","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccP150200->Sumw2();
  fHistd0rphiITSMIoneSPDInAccP150200->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccP150200);

  fHistd0rphiITSMIoneSPDInAccP350450 = new TH1F("fHistd0rphiITSMIoneSPDInAccP350450","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccP350450->Sumw2();
  fHistd0rphiITSMIoneSPDInAccP350450->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccP350450);

  fHistd0rphiITSMIoneSPDInAccP500700 = new TH1F("fHistd0rphiITSMIoneSPDInAccP500700","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccP500700->Sumw2();
  fHistd0rphiITSMIoneSPDInAccP500700->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccP500700);

  fHistd0rphiITSMIoneSPDInAccP10001500 = new TH1F("fHistd0rphiITSMIoneSPDInAccP10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccP10001500->Sumw2();
  fHistd0rphiITSMIoneSPDInAccP10001500->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccP10001500);

  fHistd0rphiITSMIoneSPDInAccP25004000 = new TH1F("fHistd0rphiITSMIoneSPDInAccP25004000","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccP25004000->Sumw2();
  fHistd0rphiITSMIoneSPDInAccP25004000->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccP25004000);

  fHistd0rphiITSMIoneSPDInAccP40008000 = new TH1F("fHistd0rphiITSMIoneSPDInAccP40008000","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccP40008000->Sumw2();
  fHistd0rphiITSMIoneSPDInAccP40008000->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccP40008000);

  fHistd0zITSMIoneSPDInAccP150200 = new TH1F("fHistd0zITSMIoneSPDInAccP150200","Longitudinal imp. par. to VertexTracks for primaries; d_{0} z [cm]; N tracks",300,-1.5,1.5);
  fHistd0zITSMIoneSPDInAccP150200->Sumw2();
  fHistd0zITSMIoneSPDInAccP150200->SetMinimum(0);
  fOutput->Add(fHistd0zITSMIoneSPDInAccP150200);

  fHistd0zITSMIoneSPDInAccP500700 = new TH1F("fHistd0zITSMIoneSPDInAccP500700","Longitudinal imp. par. to VertexTracks for primaries; d_{0} z [cm]; N tracks",300,-1.5,1.5);
  fHistd0zITSMIoneSPDInAccP500700->Sumw2();
  fHistd0zITSMIoneSPDInAccP500700->SetMinimum(0);
  fOutput->Add(fHistd0zITSMIoneSPDInAccP500700);

  fHistd0zITSMIoneSPDInAccP10001500 = new TH1F("fHistd0zITSMIoneSPDInAccP10001500","Longitudinal imp. par. to VertexTracks for primaries; d_{0} z [cm]; N tracks",300,-1.5,1.5);
  fHistd0zITSMIoneSPDInAccP10001500->Sumw2();
  fHistd0zITSMIoneSPDInAccP10001500->SetMinimum(0);
  fOutput->Add(fHistd0zITSMIoneSPDInAccP10001500);

  fHistd0zVSetaTPCInAccP10001500 = new TH2F("fHistd0zVSetaTPCInAccP10001500","Long. imp. par. to VertexSPD for TPC tracks; d_{0} z [cm]; eta",100,-4,4,10,-1,1);
  fOutput->Add(fHistd0zVSetaTPCInAccP10001500);

  fHistd0rphiVSphiITSMIoneSPDInAccP10001500 = new TH2F("fHistd0rphiVSphiITSMIoneSPDInAccP10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; phi",30,-0.3,0.3,40,0,2*3.1415);
  fOutput->Add(fHistd0rphiVSphiITSMIoneSPDInAccP10001500);

  fHistd0rphiVSetaITSMIoneSPDInAccP10001500 = new TH2F("fHistd0rphiVSetaITSMIoneSPDInAccP10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; eta",30,-0.3,0.3,10,-1,1);
  fOutput->Add(fHistd0rphiVSetaITSMIoneSPDInAccP10001500);

  fHistd0rphiITSMIoneSPDInAccS150200 = new TH1F("fHistd0rphiITSMIoneSPDInAccS150200","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS150200->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS150200->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS150200);

  fHistd0rphiITSMIoneSPDInAccS350450 = new TH1F("fHistd0rphiITSMIoneSPDInAccS350450","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS350450->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS350450->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS350450);

  fHistd0rphiITSMIoneSPDInAccS500700 = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700);

  fHistd0rphiITSMIoneSPDInAccS500700from22 = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700from22","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700from22->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700from22->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700from22);

  fHistd0rphiITSMIoneSPDInAccS500700from211 = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700from211","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700from211->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700from211->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700from211);

  fHistd0rphiITSMIoneSPDInAccS500700from310 = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700from310","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700from310->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700from310->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700from310);

  fHistd0rphiITSMIoneSPDInAccS500700from321 = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700from321","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700from321->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700from321->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700from321);

  fHistd0rphiITSMIoneSPDInAccS500700from3122 = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700from3122","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700from3122->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700from3122->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700from3122);

  fHistd0rphiITSMIoneSPDInAccS10001500 = new TH1F("fHistd0rphiITSMIoneSPDInAccS10001500","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS10001500->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS10001500->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS10001500);

  fHistd0rphiITSMIoneSPDInAccS25004000 = new TH1F("fHistd0rphiITSMIoneSPDInAccS25004000","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS25004000->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS25004000->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS25004000);

  fHistd0rphiITSMIoneSPDInAccS40008000 = new TH1F("fHistd0rphiITSMIoneSPDInAccS40008000","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS40008000->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS40008000->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS40008000);

  fHistd0rphiITSMIoneSPDInAccS150200fromStrange = new TH1F("fHistd0rphiITSMIoneSPDInAccS150200fromStrange","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS150200fromStrange->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS150200fromStrange->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS150200fromStrange);

  fHistd0rphiITSMIoneSPDInAccS150200fromMat = new TH1F("fHistd0rphiITSMIoneSPDInAccS150200fromMat","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS150200fromMat->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS150200fromMat->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS150200fromMat);

  fHistd0rphiITSMIoneSPDInAccS350450fromStrange = new TH1F("fHistd0rphiITSMIoneSPDInAccS350450fromStrange","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS350450fromStrange->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS350450fromStrange->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS350450fromStrange);

  fHistd0rphiITSMIoneSPDInAccS350450fromMat = new TH1F("fHistd0rphiITSMIoneSPDInAccS350450fromMat","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS350450fromMat->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS350450fromMat->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS350450fromMat);

  fHistd0rphiITSMIoneSPDInAccS500700fromStrange = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700fromStrange","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700fromStrange->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700fromStrange->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700fromStrange);

  fHistd0rphiITSMIoneSPDInAccS500700fromMat = new TH1F("fHistd0rphiITSMIoneSPDInAccS500700fromMat","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS500700fromMat->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS500700fromMat->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS500700fromMat);

  fHistd0rphiITSMIoneSPDInAccS10001500fromStrange = new TH1F("fHistd0rphiITSMIoneSPDInAccS10001500fromStrange","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS10001500fromStrange->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS10001500fromStrange->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS10001500fromStrange);

  fHistd0rphiITSMIoneSPDInAccS10001500fromMat = new TH1F("fHistd0rphiITSMIoneSPDInAccS10001500fromMat","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS10001500fromMat->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS10001500fromMat->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS10001500fromMat);

  fHistd0rphiITSMIoneSPDInAccS25004000fromStrange = new TH1F("fHistd0rphiITSMIoneSPDInAccS25004000fromStrange","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS25004000fromStrange->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS25004000fromStrange->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS25004000fromStrange);

  fHistd0rphiITSMIoneSPDInAccS25004000fromMat = new TH1F("fHistd0rphiITSMIoneSPDInAccS25004000fromMat","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS25004000fromMat->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS25004000fromMat->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS25004000fromMat);

  fHistd0rphiITSMIoneSPDInAccS40008000fromStrange = new TH1F("fHistd0rphiITSMIoneSPDInAccS40008000fromStrange","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS40008000fromStrange->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS40008000fromStrange->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS40008000fromStrange);

  fHistd0rphiITSMIoneSPDInAccS40008000fromMat = new TH1F("fHistd0rphiITSMIoneSPDInAccS40008000fromMat","Transverse imp. par. to VertexTracks for primaries; d_{0} rphi [cm]; N tracks",300,-1.5,1.5);
  fHistd0rphiITSMIoneSPDInAccS40008000fromMat->Sumw2();
  fHistd0rphiITSMIoneSPDInAccS40008000fromMat->SetMinimum(0);
  fOutput->Add(fHistd0rphiITSMIoneSPDInAccS40008000fromMat);

  
  fHistd0zITSMIoneSPDInAccS150200 = new TH1F("fHistd0zITSMIoneSPDInAccS150200","Longitudinal imp. par. to VertexTracks for secondaries; d_{0} z [cm]; N tracks",300,-1.5,1.5);
  fHistd0zITSMIoneSPDInAccS150200->Sumw2();
  fHistd0zITSMIoneSPDInAccS150200->SetMinimum(0);
  fOutput->Add(fHistd0zITSMIoneSPDInAccS150200);

  fHistd0zITSMIoneSPDInAccS500700 = new TH1F("fHistd0zITSMIoneSPDInAccS500700","Longitudinal imp. par. to VertexTracks for secondaries; d_{0} z [cm]; N tracks",300,-1.5,1.5);
  fHistd0zITSMIoneSPDInAccS500700->Sumw2();
  fHistd0zITSMIoneSPDInAccS500700->SetMinimum(0);
  fOutput->Add(fHistd0zITSMIoneSPDInAccS500700);

  fHistd0zITSMIoneSPDInAccS10001500 = new TH1F("fHistd0zITSMIoneSPDInAccS10001500","Longitudinal imp. par. to VertexTracks for secondaries; d_{0} z [cm]; N tracks",300,-1.5,1.5);
  fHistd0zITSMIoneSPDInAccS10001500->Sumw2();
  fHistd0zITSMIoneSPDInAccS10001500->SetMinimum(0);
  fOutput->Add(fHistd0zITSMIoneSPDInAccS10001500);

  fHistITSRedChi2NonFakePt02 = new TH1F("fHistITSRedChi2NonFakePt02","ITS chi2/cluster; ITS #chi^{2}/nclusters; tracks",500,0,100);
  fOutput->Add(fHistITSRedChi2NonFakePt02);
  fHistITSRedChi2FakePt02 = new TH1F("fHistITSRedChi2FakePt02","ITS chi2/cluster; ITS #chi^{2}/nclusters; tracks",500,0,100);
  fOutput->Add(fHistITSRedChi2FakePt02);
  fHistITSRedChi2NonFakePt05 = new TH1F("fHistITSRedChi2NonFakePt05","ITS chi2/cluster; ITS #chi^{2}/nclusters; tracks",500,0,100);
  fOutput->Add(fHistITSRedChi2NonFakePt05);
  fHistITSRedChi2FakePt05 = new TH1F("fHistITSRedChi2FakePt05","ITS chi2/cluster; ITS #chi^{2}/nclusters; tracks",500,0,100);
  fOutput->Add(fHistITSRedChi2FakePt05);
  fHistITSRedChi2NonFakePt1 = new TH1F("fHistITSRedChi2NonFakePt1","ITS chi2/cluster; ITS #chi^{2}/nclusters; tracks",500,0,100);
  fOutput->Add(fHistITSRedChi2NonFakePt1);
  fHistITSRedChi2FakePt1 = new TH1F("fHistITSRedChi2FakePt1","ITS chi2/cluster; ITS #chi^{2}/nclusters; tracks",500,0,100);
  fOutput->Add(fHistITSRedChi2FakePt1);


  // ntuples
  //
  fNtupleESDTracks = new TNtuple("fNtupleESDTracks","tracks","pt:eta:phi:d0:z0:sigmad0:sigmaz0:ptMC:pdgMC:pdgMCMoth:d0MC:d0MCv:z0MCv:sigmad0MCv:sigmaz0MCv:ITSflag:isPrimary:isTPCSel");  
  fOutput->Add(fNtupleESDTracks);

  fNtupleITSAlignExtra = new TNtuple("fNtupleITSAlignExtra","ITS alignment checks: extra clusters","layer:x:y:z:dxy:dz:xloc:zloc:npoints:pt");  
  fOutput->Add(fNtupleITSAlignExtra);

  fNtupleITSAlignSPDTracklets = new TNtuple("fNtupleITSAlignSPDTracklets","ITS alignment checks: SPD tracklets wrt SPD vertex","phi:theta:z:dxy:dz:pt");  
  fOutput->Add(fNtupleITSAlignSPDTracklets);

  PostData(1, fOutput);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event


  fESD = dynamic_cast<AliESDEvent*>(InputEvent());

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  //if(fESD->GetEventNumberInFile()<7) return;

  // only events in the requested multiplicity range
  if(!IsSelectedCentrality()) return;

  fHistNEvents->Fill(-1);


  Bool_t isSelected = kTRUE;
  if(fUsePhysSel) {
    isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  }
  if(!isSelected) return;

  TString firedTriggerClasses=fESD->GetFiredTriggerClasses();  
  if(!firedTriggerClasses.Contains(fTriggerClass.Data())) return;
  

  Bool_t sddIsIn=kTRUE;

      AliCDBManager* man = AliCDBManager::Instance();
      if(fOCDBPath.Contains("OCDB")) { // when running in the QAtrain this is not called (OCBD is already set)
	man->SetDefaultStorage(fOCDBPath.Data());
	man->SetRun(fESD->GetRunNumber());
      }
      if(!man) {       
	AliFatal("CDB not set but needed by AliAnalysisTaskITSTrackingCheck");
	return;    
      }   
      AliCDBEntry* eT=(AliCDBEntry*)man->Get("GRP/CTP/Config");    
      if(eT) {      
	fTrigConfig=(AliTriggerConfiguration*)eT->GetObject();    
      }    
      if(!eT || !fTrigConfig) {      
	AliError("Cannot retrieve CDB entry for GRP/CTP/Config");      
	return;     
      }
    

    sddIsIn=fESD->IsDetectorInTriggerCluster("ITSSDD",fTrigConfig);
  


  //if(fESD->GetEventType()!=7) return;

  // ***********  MC info ***************
  TArrayF mcVertex(3);
  mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
  Float_t dNchdy=-999.;

  TParticle *part=0;
  AliESDVertex *vertexMC=0;
  AliMCEvent* mcEvent=0x0;
  if (fReadMC) {
    AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    
    AliHeader* header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    AliGenEventHeader* genHeader = header->GenEventHeader();
    genHeader->PrimaryVertex(mcVertex);
    //printf("MC vertex: %f %f %f\n",mcVertex[0],mcVertex[1],mcVertex[2]);

    Int_t ngenpart = (Int_t)mcEvent->GetNumberOfTracks();
    //printf("# generated particles = %d\n",ngenpart);
    dNchdy=0;
    for(Int_t ip=0; ip<ngenpart; ip++) {
      part = mcEvent->Particle(ip);
      if (!part) continue;
      // keep only electrons, muons, pions, kaons and protons
      Int_t apdg = TMath::Abs(part->GetPdgCode());
      if(apdg!=11 && apdg!=13 && apdg!=211 && apdg!=321 && apdg!=2212) continue;      
      // reject secondaries
      //if(TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1]))<.01) printf("%f\n",TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1])));
      if(TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1]))>0.0010) continue;
      // reject incoming protons
      Double_t energy  = part->Energy();
      if(energy>900.) continue;
      Double_t pz = part->Pz();
      Double_t y = 100.;
      if(energy-pz+1.e-13 > 0.) y = 0.5*TMath::Log((energy+pz+1.e-13)/(energy-pz+1.e-13));
      if(TMath::Abs(y)<1.0) dNchdy += 0.5; // count 1/2 of particles in |y|<1
    }
    //printf("# primary particles = %7.1f\n",dNchdy);
  } 
  // ***********  MC info ***************
  Double_t mcVtxPos[3]={mcVertex[0],mcVertex[1],mcVertex[2]},mcVtxSigma[3]={0,0,0};
  vertexMC = new AliESDVertex(mcVtxPos,mcVtxSigma);



  //------- event selection --------
  Int_t   mincontrSPDvtx=1;
  Double_t maxzSPDvtx=5.;
  Double_t maxrSPDvtx=0.5;
  Bool_t useTRKvtx=kTRUE;
  Bool_t useSPDvtxifNotTRK=kTRUE;
  Int_t   mincontrTRKvtx=1;
  //Double_t maxzTRKvtx=20.;
  //Double_t maxrTRKvtx=1.;
  //------- TPC track selection --------
  Int_t    minclsTPC=70;
  Double_t maxchi2perTPCcl=4.;
  Double_t minEtaInAcc=-0.8; // -0.8
  Double_t maxEtaInAcc=0.8; // 0.8
  Double_t maxdcaxy=2.4;//2.4;
  Double_t maxdcaz=3.2;//3.2;
  AliESDtrackCuts* esdtrackCutsTPC = new AliESDtrackCuts("esdtrackCutsTPC");
  esdtrackCutsTPC->SetMaxDCAToVertexXY(maxdcaxy);
  esdtrackCutsTPC->SetMaxDCAToVertexZ(maxdcaz);
  esdtrackCutsTPC->SetDCAToVertex2D(kTRUE);
  esdtrackCutsTPC->SetRequireSigmaToVertex(kFALSE);
  esdtrackCutsTPC->SetRequireTPCRefit(kFALSE);// cannot do it because status not copied in AliESDtrack::FillTPCOnlyTrack
  esdtrackCutsTPC->SetAcceptKinkDaughters(kFALSE);
  esdtrackCutsTPC->SetMinNClustersTPC(minclsTPC);
  esdtrackCutsTPC->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
  esdtrackCutsTPC->SetEtaRange(minEtaInAcc,maxEtaInAcc);
  SetESDtrackCutsTPC(esdtrackCutsTPC);
  //------- ITS+TPC track selection --------
  Double_t maxdcaxyITSTPC=0.2;
  Double_t maxdcazITSTPC=1.e6;
  AliESDtrackCuts* esdtrackCutsITSTPC = new AliESDtrackCuts("esdtrackCutsITSTPC");
  esdtrackCutsITSTPC->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
  esdtrackCutsITSTPC->SetMaxDCAToVertexZ(maxdcazITSTPC);
  esdtrackCutsITSTPC->SetDCAToVertex2D(kFALSE);
  esdtrackCutsITSTPC->SetRequireSigmaToVertex(kFALSE);
  esdtrackCutsITSTPC->SetRequireITSRefit(kTRUE);
  esdtrackCutsITSTPC->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					       AliESDtrackCuts::kAny);
  esdtrackCutsITSTPC->SetAcceptKinkDaughters(kFALSE);
  esdtrackCutsITSTPC->SetMinNClustersTPC(minclsTPC);
  esdtrackCutsITSTPC->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
  esdtrackCutsITSTPC->SetEtaRange(minEtaInAcc,maxEtaInAcc);
  SetESDtrackCutsITSTPC(esdtrackCutsITSTPC);
  //---------------------------------------
  

  //
  
  /*  
  // **********  Trigger *****************
  ULong64_t triggerMask;
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);
  
  triggerMask=fESD->GetTriggerMask();
  // MB1: SPDFO || V0L || V0R
  Bool_t eventTriggered = (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right))); 
  //MB2: GFO && V0R
  //triggerMask & spdFO && ((triggerMask&v0left) || (triggerMask&v0right))
  // ************ Trigger ******************
  if(!eventTriggered) return;
  */


  fHistNEvents->Fill(0);

  // SPD vertex
  const AliESDVertex *spdv=fESD->GetPrimaryVertexSPD();
  // Select good SPD vertices
  TString spdvtitle=spdv->GetTitle();
  //if(!spdvtitle.Contains("3D")) return;
  if(fRequireSPDvtx) {
    if(spdv->GetNContributors()<mincontrSPDvtx ||
       TMath::Abs(spdv->GetZ())>maxzSPDvtx ||  
       spdv->GetX()*spdv->GetX()+spdv->GetY()*spdv->GetY()>maxrSPDvtx) {
      delete esdtrackCutsTPC; esdtrackCutsTPC=0;
      delete esdtrackCutsITSTPC; esdtrackCutsITSTPC=0;
      return;
    }
  }
  //spdv->Print();

  if(fRejPileupSPD) if(fESD->IsPileupFromSPD()) return;

  //
  // Tracks vertex
  const AliESDVertex *vertexESD = fESD->GetPrimaryVertexTracks();
  //vertexESD->Print();

  if(spdvtitle.Contains("3D")) {
    fHistNEvents->Fill(1);
    fHistNEvents->Fill(3);
  } else {
    fHistNEvents->Fill(2);
    fHistNEvents->Fill(3);
  }
  if(vertexESD) {
    if(vertexESD->GetStatus()) fHistNEvents->Fill(4);
  }
    


  Int_t ntracks = fESD->GetNumberOfTracks();
  //printf("Tracks # = %d\n",fESD->GetNumberOfTracks());

  fHistNtracks->Fill(ntracks);
  // Post the data already here
  PostData(1, fOutput);

  FillNClustersModuleMap();

  Int_t idet,status; Float_t xloc,zloc;
  Double_t rSPDouter=7.6,rSDDouter=23.9,rSSDouter=43.1;  
  Double_t zSPDouter=14.1,zSDDouter=29.7,zSSDouter=48.9;  


  // loop on tracks to check if there are 7 ITS layers (upgrade)
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack *track = fESD->GetTrack(itr);
    if(track->HasPointOnITSLayer(6)) {
      fNITSLayers=7;
      break;
    }
  }

    
    const AliESDVertex* spdMeanV = (const AliESDVertex*)man->Get("GRP/Calib/MeanVertexSPD")->GetObject();
  Double_t  sigZLR = spdMeanV->GetZRes();
 Double_t   zIP = spdMeanV->GetZ();
  // loop on tracks
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack *track = fESD->GetTrack(itr);
    //if(track->Charge()>0) continue;
    // remove kink daughters
    if(track->GetKinkIndex(0)>0) continue;


    // remove tracks not reco in ITS or TPC
    if (!(track->GetStatus() & AliESDtrack::kITSin) &&
	!(track->GetStatus() & AliESDtrack::kTPCin)) continue;

    // check if tracks if flagged as from bc0 by TOF
      Float_t dca[2],dcaCov[3];
      track->GetImpactParametersTPC(dca,dcaCov);
    
      
    Bool_t isTOFbc0 = kFALSE;
    if((track->IsOn(AliESDtrack::kTPCrefit))&&(TMath::Abs(track->GetTPCInnerParam()->GetX())<3)&&(track->GetNcls(1)>80)&&(track->GetTPCchi2()/track->GetNcls(1)<4)&&(TMath::Abs(track->GetTPCInnerParam()->GetZ()-zIP)<3*sigZLR)&&(dcaCov[0]>1e-9)&&(dca[0]*dca[0]/(dcaCov[0]+1e-2)<50)&&(track->GetTOFBunchCrossing()==0))
          
          isTOFbc0=kTRUE;
      

    //
    Bool_t isPrimary=kTRUE,isFromMat=kFALSE,isFromStrange=kFALSE;
    Double_t rProdVtx=0;
    //    Double_t zProdVtx=0;
    Int_t pdgTrk=0,pdgMoth=0;
    Int_t nClsMCSPD=0;
  
    Int_t trkLabel = TMath::Abs(track->GetLabel());
    Bool_t isFake=(track->GetLabel()>=0 ? kFALSE : kTRUE);
    // Bool_t hasShared=kFALSE,hasSharedSPD=kFALSE;
    // for(Int_t i=0;i<fNITSLayers;i++) {
    //   if(track->HasSharedPointOnITSLayer(i)) {
    // 	hasShared=kTRUE;
    // 	if(i<(fNITSLayers-4)) hasSharedSPD=kTRUE;
    //   }
    // }

    // check if it is primary
    if(fReadMC && mcEvent) {
      part = mcEvent->Particle(trkLabel);
      if (!part) continue;
      isPrimary = mcEvent->IsPhysicalPrimary(trkLabel);
      rProdVtx = TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1]));
      //      zProdVtx = TMath::Abs(part->Vz()-mcVertex[2]);
      //if(rProdVtx<2.8) isPrimary=kTRUE; // this could be tried
      pdgTrk = TMath::Abs(part->GetPdgCode());
      AliMCParticle* mcPart=(AliMCParticle*)mcEvent->GetTrack(trkLabel);
      Int_t idMother=mcPart->GetMother();
      if(idMother>=0) {
	TParticle* mm = mcEvent->Particle(idMother);
	if(mm) pdgMoth = TMath::Abs(mm->GetPdgCode());
      }
      if(pdgMoth==310 || pdgMoth==321 || pdgMoth==3122 || pdgMoth==3312) isFromStrange=kTRUE;
      if(pdgMoth==211 || pdgMoth==22 || pdgMoth==2112 || pdgMoth==2212) isFromMat=kTRUE;
      nClsMCSPD = NumberOfITSClustersMC(trkLabel,240);
    }

    Bool_t itsrefit=kFALSE,tpcrefit=kFALSE,itsfindable=kFALSE,itsfindableAcc=kFALSE;
    if ((track->GetStatus() & AliESDtrack::kITSrefit)) itsrefit=kTRUE;
    if ((track->GetStatus() & AliESDtrack::kTPCrefit)) tpcrefit=kTRUE;
    //if ((track->GetStatus() & AliESDtrack::kTPCin)) tpcrefit=kTRUE;



    // remove tracks with kTPCin and not kTPCrefit 
    if ((track->GetStatus() & AliESDtrack::kTPCin) && !tpcrefit) continue;

    AliESDtrack *trackTPC = 0;
    if(tpcrefit) trackTPC = AliESDtrackCuts::GetTPCOnlyTrack(fESD,itr);
    //    if(trackTPC) trackTPC->RelateToVertex(spdv,fESD->GetMagneticField(),100.); // relate it to the SPD vertex
   

    Int_t nclsITS = track->GetNcls(0);
    Int_t nclsSDDSSD = 0;
    Int_t nclsokbadoutinzITS = 0;
    Bool_t outInZ=kFALSE;

    for(Int_t layer=0; layer<fNITSLayers; layer++) {

      if(TESTBIT(track->GetITSClusterMap(),layer)) {
	if(tpcrefit) {
	  fHistClusterMapITSMI->Fill(layer);
	} else {
	  fHistClusterMapITSSA->Fill(layer);
	  if(!outInZ) fHistClusterMapITSSAInAcc->Fill(layer);
	}
      }

      if(layer>=(fNITSLayers-4) && track->HasPointOnITSLayer(layer)) nclsSDDSSD++;
      if(layer==0 && !track->HasPointOnITSLayer(1)) continue;
      if(layer==1 && !track->HasPointOnITSLayer(0)) continue;
      Bool_t goodIndex=track->GetITSModuleIndexInfo(layer,idet,status,xloc,zloc);
      if(!goodIndex) continue;
      if(status<=0) continue; // status cannot be zero
      if(idet<0) continue; // wrong module retrieved
      if(layer==0 && idet>=80) continue; // Inner SPD has modules from 0 to 79
      else if(layer==1 && (idet<80 || idet>=240)) continue; // Outer SPD has modules from 80 to 239
      else if(layer==2 && idet>=84) continue; // Inner SDD from 0 to 83
      else if(layer==3 && (idet<84 || idet>=260)) continue; // Outer SDD from 84 to 259
      else if(layer==4 && idet>=748) continue; // Inner SSD from 0 to 747
      else if(layer==5 && (idet<748 || idet>=1698)) continue; // Outer SSD from 748 to 1697
   
      if(layer>=2) idet+=240; // add n SPD modules
      if(layer>=4) idet+=260; // add n SDD modules
      if(status==4) outInZ=kTRUE;
      
      if(layer==1 && tpcrefit && itsrefit) {
	Double_t xyz[3];
	track->GetXYZAt(rSPDouter,fESD->GetMagneticField(),xyz);
	if(status==1) {
	  if(idet%4==0) fHistZatSPDouter0ok->Fill(xyz[2]);
	  if(idet%4==1) fHistZatSPDouter1ok->Fill(xyz[2]);
	  if(idet%4==2) fHistZatSPDouter2ok->Fill(xyz[2]);
	  if(idet%4==3) fHistZatSPDouter3ok->Fill(xyz[2]);
	}
	if(status==2 || status==5) {
	  if(idet%4==0) fHistZatSPDouter0notok->Fill(xyz[2]);
	  if(idet%4==1) fHistZatSPDouter1notok->Fill(xyz[2]);
	  if(idet%4==2) fHistZatSPDouter2notok->Fill(xyz[2]);
	  if(idet%4==3) fHistZatSPDouter3notok->Fill(xyz[2]);
	}
      }
      // kill SDD clusters in MC
      if(fReadMC && (layer==2 || layer==3)) {
	Double_t qITS[4]; track->GetITSdEdxSamples(qITS);
	Double_t conv=3.34/5.43; 
	if(layer==2 && status==1 && qITS[0]*conv < 30.) status=5;
	if(layer==3 && status==1 && qITS[1]*conv < 30.) status=5;
      }
      if(tpcrefit) {
	if(trackTPC) if(!fESDtrackCutsTPC->AcceptTrack(trackTPC)) continue;
	if(track->Pt()>0.7 && track->Pt()<2.0) {
	  if(status==1) {
	    fHistClusterMapITSMIok->Fill(layer);
	  }
	  if(status==2) fHistClusterMapITSMIbad->Fill(layer);
	  if(status==3) fHistClusterMapITSMIskipped->Fill(layer);
	  if(status==4) fHistClusterMapITSMIoutinz->Fill(layer);
	  if(status==5) fHistClusterMapITSMInocls->Fill(layer);
	  if(status==6) fHistClusterMapITSMInorefit->Fill(layer);
	  if(status==1 && !outInZ) fHistClusterMapModuleITSMIokInAcc->Fill(idet);
	  if((status==2 || status==7) && !outInZ) fHistClusterMapModuleITSMIbadInAcc->Fill(idet);
	  if(status==5 && !outInZ) fHistClusterMapModuleITSMInoclsInAcc->Fill(idet);


	}
	if(status==1 || status==2 || status==4) {
	  fHistClusterMapITSMIokoutinzbad->Fill(layer);
	  nclsokbadoutinzITS++;
	}
	/*if(idet>=80 && idet<120) {
	  Double_t xyz[3];
	  track->GetXYZAt(6.8,fESD->GetMagneticField(),xyz);
	  printf("det %d status %d phi %f z %f zloc %f xloc %f\n",idet,status,TMath::ATan2(xyz[1],xyz[0]),xyz[2],zloc,xloc);
	  }*/
	if((layer==2 || layer==3) && status!=2 && status!=4) {
	  //Float_t xlocCls;
	  //Int_t nClsInMod = NumberOfITSClusters(idet,xlocCls);
	  fHistxlocSDDall->Fill(xloc);
	  fHistzlocSDDall->Fill(zloc);
	  if(status==1  /*nClsInMod>0*/) {
	    fHistxlocSDDok->Fill(xloc);
	    fHistxlocVSmodSDDok->Fill(idet,xloc);
	    fHistzlocSDDok->Fill(zloc);

	  } 	
	}
      }
    }  

    // TPC track findable in ITS
    if(tpcrefit && trackTPC) { 
      if(fESDtrackCutsTPC->AcceptTrack(trackTPC)) {
	itsfindable=kTRUE;
	Double_t zAtSSDouter=100,zAtSDDouter=100,zAtSPDouter=100;
	track->GetZAt(rSSDouter,fESD->GetMagneticField(),zAtSSDouter);
	track->GetZAt(rSDDouter,fESD->GetMagneticField(),zAtSDDouter);
	track->GetZAt(rSPDouter,fESD->GetMagneticField(),zAtSPDouter);
	fHistPtTPC->Fill(track->Pt());  
	if(TMath::Abs(zAtSSDouter)<1.*zSSDouter &&
	   TMath::Abs(zAtSDDouter)<1.*zSDDouter &&
	   TMath::Abs(zAtSPDouter)<1.*zSPDouter) {
	  itsfindableAcc=kTRUE; 
	  fHistdEdxVSPtTPCInAcc->Fill(track->Pt(),track->GetTPCsignal());
	  fHistPtTPCInAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtTPCInAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtTPCInAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtTPCInAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtTPCInAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtTPCInAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  fHistTPCclsVSPtTPCInAcc->Fill(trackTPC->Pt(),trackTPC->GetNcls(1));
	  fHistPtVSphiTPCInAcc->Fill(track->Phi(),track->Pt());
	  if(!(track->GetStatus()&AliESDtrack::kTRDout)) fHistPtTPCInAccNoTRDout->Fill(track->Pt()); 
	  if(!(track->GetStatus()&AliESDtrack::kTOFout)) fHistPtTPCInAccNoTOFout->Fill(track->Pt()); 
	  fHistPtTPCInAccWithPtTPCAtVtx->Fill(trackTPC->Pt());
	  Double_t pTPCinnerwall[3];
	  track->GetInnerPxPyPz(pTPCinnerwall);
	  Double_t ptTPCinnerwall=TMath::Sqrt(pTPCinnerwall[0]*pTPCinnerwall[0]+pTPCinnerwall[1]*pTPCinnerwall[1]);
	  fHistPtTPCInAccWithPtTPCAtInnerWall->Fill(ptTPCinnerwall);
	  if(track->GetStatus()&AliESDtrack::kITSrefit && !(track->GetStatus()&AliESDtrack::kTOFout)) fHistDeltaPtTPC->Fill(trackTPC->P(),track->P()-trackTPC->P());
	  if(track->Pt()>1) {
	    fHistPhiTPCInAcc->Fill(track->Phi());  
	    fHistEtaTPCInAcc->Fill(track->Eta());
	  }  
	  if(isPrimary) {
	    fHistPtTPCInAccP->Fill(track->Pt());
	    if(pdgTrk==321) fHistPtTPCInAccPfromStrange->Fill(track->Pt());
	  } else {
	    fHistPtTPCInAccS->Fill(track->Pt());
	    if(isFromMat) fHistPtTPCInAccSfromMat->Fill(track->Pt());
	    if(isFromStrange) fHistPtTPCInAccSfromStrange->Fill(track->Pt());
	  }  
	  if(nClsMCSPD==2) fHistPtTPCInAccMCtwoSPD->Fill(track->Pt());
	  if(nClsMCSPD>=1) fHistPtTPCInAccMConeSPD->Fill(track->Pt());

	  //if(isPrimary) {fHistRProdVtxInAccP->Fill(rProdVtx);} else {fHistRProdVtxInAccS->Fill(rProdVtx);}
	}
      }
    }


    if(useTRKvtx) {
      // we need the vertex to compute the impact parameters
      if(!vertexESD) {delete trackTPC; trackTPC=0; continue;}
      if(!(vertexESD->GetStatus()) || vertexESD->GetNContributors()<mincontrTRKvtx) {
	if(useSPDvtxifNotTRK) {
	  vertexESD = fESD->GetPrimaryVertexSPD();
	} else {
	  delete trackTPC; trackTPC=0; continue;
	}
      }
    }

    //if(hasSharedSPD) continue;
    //if(track->GetITSFakeFlag()) continue;

    // track prolonged in ITS with different conditions
    if(itsrefit) {
      if(itsfindable) {
	if(nclsITS==7) fHistPtITSMI7->Fill(track->Pt());
	if(nclsITS==6) fHistPtITSMI6->Fill(track->Pt());
	if(nclsITS==5) fHistPtITSMI5->Fill(track->Pt());
	if(nclsITS==4) fHistPtITSMI4->Fill(track->Pt());
	if(nclsITS==3) fHistPtITSMI3->Fill(track->Pt());
	if(nclsITS==2) fHistPtITSMI2->Fill(track->Pt());
	if(ConditionSPD(track)) {
	  fHistPtITSMISPD->Fill(track->Pt());
	  fHistTPCclsVSPtITSMISPDInAcc->Fill(track->Pt(),track->GetNcls(1));
	}
	if(ConditionSPDone(track)) {
	  fHistPtITSMIoneSPD->Fill(track->Pt());
	}
	if(ConditionSPDtwo(track)) {
	  fHistPtITSMItwoSPD->Fill(track->Pt());
	}
	if(nclsokbadoutinzITS==6) fHistPtITSMIokbadoutinz6->Fill(track->Pt());
      }
      if(itsfindableAcc) {
	if(nclsITS==7) {
	  fHistPtITSMI7InAcc->Fill(track->Pt());
	  if(track->Pt()>1) {
	    fHistPhiITSMI7InAcc->Fill(track->Phi());
	    fHistEtaITSMI7InAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMI7InAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMI7InAccP->Fill(track->Pt());} else {fHistPtITSMI7InAccS->Fill(track->Pt());}  
	}
	if(nclsITS==6) {
	  fHistPtITSMI6InAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI6InAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI6InAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMI6InAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI6InAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI6InAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMI6InAcc->Fill(track->Phi());
	    fHistEtaITSMI6InAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMI6InAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMI6InAccP->Fill(track->Pt());} else {fHistPtITSMI6InAccS->Fill(track->Pt());}  
	}
	if(nclsITS==5) {
	  fHistPtITSMI5InAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI5InAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI5InAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMI5InAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI5InAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI5InAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMI5InAcc->Fill(track->Phi());
	    fHistEtaITSMI5InAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMI5InAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMI5InAccP->Fill(track->Pt());} else {fHistPtITSMI5InAccS->Fill(track->Pt());}  
	}
	if(nclsITS==4) {
	  fHistPtITSMI4InAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI4InAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI4InAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMI4InAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI4InAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI4InAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMI4InAcc->Fill(track->Phi());
	    fHistEtaITSMI4InAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMI4InAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMI4InAccP->Fill(track->Pt());} else {fHistPtITSMI4InAccS->Fill(track->Pt());}  
	}
	if(nclsITS==3) {
	  fHistPtITSMI3InAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI3InAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI3InAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMI3InAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI3InAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI3InAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMI3InAcc->Fill(track->Phi());
	    fHistEtaITSMI3InAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMI3InAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMI3InAccP->Fill(track->Pt());} else {fHistPtITSMI3InAccS->Fill(track->Pt());}  
	}
	if(nclsITS==2) {
	  fHistPtITSMI2InAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI2InAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI2InAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMI2InAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMI2InAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMI2InAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMI2InAcc->Fill(track->Phi());
	    fHistEtaITSMI2InAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMI2InAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMI2InAccP->Fill(track->Pt());} else {fHistPtITSMI2InAccS->Fill(track->Pt());}  
	}
	if(ConditionSPD(track)) {
	  fHistPtITSMISPDInAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMISPDInAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMISPDInAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMISPDInAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMISPDInAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMISPDInAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMISPDInAcc->Fill(track->Phi());
	    fHistEtaITSMISPDInAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMISPDInAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMISPDInAccP->Fill(track->Pt());} else {fHistPtITSMISPDInAccS->Fill(track->Pt());}  
	}
	if(ConditionSPDtwo(track)) {
	  fHistPtITSMItwoSPDInAcc->Fill(track->Pt());
	  if(track->Pt()>1) {
	    fHistPhiITSMItwoSPDInAcc->Fill(track->Phi());
	    fHistEtaITSMItwoSPDInAcc->Fill(track->Eta());
	  }
	  if(isFake) fHistPtITSMItwoSPDInAccFake->Fill(track->Pt());
	  if(isPrimary) {fHistPtITSMItwoSPDInAccP->Fill(track->Pt());} else {fHistPtITSMItwoSPDInAccS->Fill(track->Pt());}  
	}
	if(ConditionSPDone(track)) {
	  fHistPtITSMIoneSPDInAcc->Fill(track->Pt());
	  if(fCheckSDDIsIn && sddIsIn) fHistPtITSMIoneSPDInAccwSDD->Fill(track->Pt());
	  if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMIoneSPDInAccwoSDD->Fill(track->Pt());
	  if(isTOFbc0) {
	    fHistPtITSMIoneSPDInAccTOFbc0->Fill(track->Pt());
	    if(fCheckSDDIsIn && sddIsIn) fHistPtITSMIoneSPDInAccTOFbc0wSDD->Fill(track->Pt());
	    if(fCheckSDDIsIn && !sddIsIn) fHistPtITSMIoneSPDInAccTOFbc0woSDD->Fill(track->Pt());
	  }
	  if(track->Pt()>1) {
	    fHistPhiITSMIoneSPDInAcc->Fill(track->Phi());
	    fHistEtaITSMIoneSPDInAcc->Fill(track->Eta());
	  }
	  if(track->HasSharedPointOnITSLayer(0) || 
	     track->HasSharedPointOnITSLayer(1) || 
	     track->HasSharedPointOnITSLayer(2) || 
	     track->HasSharedPointOnITSLayer(3) || 
	     track->HasSharedPointOnITSLayer(4) || 
	     track->HasSharedPointOnITSLayer(5)) {
	    fHistPtITSMIoneSPDInAccShared->Fill(track->Pt());
	    if(isFake) fHistPtITSMIoneSPDInAccSharedFake->Fill(track->Pt());
	  }
	  if(track->HasSharedPointOnITSLayer(0) || track->HasSharedPointOnITSLayer(1)) {
	    fHistPtITSMIoneSPDInAccSharedSPD->Fill(track->Pt());
	    if(isFake) fHistPtITSMIoneSPDInAccSharedSPDFake->Fill(track->Pt());
	  }
	  if(track->HasSharedPointOnITSLayer(0)) fHistPtITSMISPD1InAccShared->Fill(track->Pt());
	  if(track->HasSharedPointOnITSLayer(1)) fHistPtITSMISPD2InAccShared->Fill(track->Pt());


	  Float_t chi2redITS = 0.;
	  if(track->GetNcls(0)>0) chi2redITS = track->GetITSchi2()/(Float_t)track->GetNcls(0);

	  if(track->Pt()>0.18 && track->Pt()<0.25)  fHistITSRedChi2NonFakePt02->Fill(chi2redITS);
	  if(track->Pt()>0.4 && track->Pt()<0.6)  fHistITSRedChi2NonFakePt05->Fill(chi2redITS);
	  if(track->Pt()>0.9 && track->Pt()<1.5)  fHistITSRedChi2NonFakePt1->Fill(chi2redITS);
	  if(isFake) {
	    fHistPtITSMIoneSPDInAccFake->Fill(track->Pt());
	    if(track->Pt()>0.18 && track->Pt()<0.25)  fHistITSRedChi2FakePt02->Fill(chi2redITS);
	    if(track->Pt()>0.4 && track->Pt()<0.6)  fHistITSRedChi2FakePt05->Fill(chi2redITS);
	    if(track->Pt()>0.9 && track->Pt()<1.5)  fHistITSRedChi2FakePt1->Fill(chi2redITS);
	  }
	  if(nclsSDDSSD>=3) fHistPtITSMIoneSPDthreeSDDSSDInAcc->Fill(track->Pt());

	  if(isPrimary) {fHistPtITSMIoneSPDInAccP->Fill(track->Pt());} else {fHistPtITSMIoneSPDInAccS->Fill(track->Pt());}  
	}
	if(nclsokbadoutinzITS==6) fHistPtITSMIokbadoutinz6InAcc->Fill(track->Pt());
	if(nclsokbadoutinzITS==5) fHistPtITSMIokbadoutinz5InAcc->Fill(track->Pt());
	if(nclsokbadoutinzITS==4) fHistPtITSMIokbadoutinz4InAcc->Fill(track->Pt());
	if(nclsokbadoutinzITS==6 && track->Pt()>1) fHistPhiITSMIokbadoutinz6InAcc->Fill(track->Phi());  
      }
    }

    if(tpcrefit) {
      fHistNclsITSMI->Fill(nclsITS);
    } else {
      fHistNclsITSSA->Fill(nclsITS);
      if(!outInZ) fHistNclsITSSAInAcc->Fill(nclsITS);
    }

  
    if((tpcrefit && fUseITSSAforNtuples) || // only ITS-SA for ntuples
       (!tpcrefit && !fUseITSSAforNtuples)) // only ITS-TPC for ntuples
      { delete trackTPC; trackTPC=0; continue; }

    // impact parameter to VertexTracks
    Float_t d0z0[2]={0,0},covd0z0[3]={0,0,0};
    Double_t d0z0TPC[2]={0,0},covd0z0TPC[3]={0,0,0};
    if(useTRKvtx) {
      if(!track->RelateToVertex(vertexESD,fESD->GetMagneticField(),kVeryBig)) { delete trackTPC; trackTPC=0; continue; }
    } else { 
      //if(!track->RelateToVertex(vertexMC,fESD->GetMagneticField(),kVeryBig)) continue;
      if(!track->RelateToVertex(spdv,fESD->GetMagneticField(),kVeryBig)) { delete trackTPC; trackTPC=0; continue; }
    } 
    track->GetImpactParameters(d0z0,covd0z0);
    if(trackTPC) trackTPC->PropagateToDCA(spdv,fESD->GetMagneticField(),kVeryBig,d0z0TPC,covd0z0TPC);
    if(covd0z0[0]<0. || covd0z0[2]<0. || covd0z0TPC[0]<0. || covd0z0TPC[2]<0.) { delete trackTPC; trackTPC=0; continue; }


    // track that passes final ITS+TPC cuts
    if(itsfindableAcc && fESDtrackCutsITSTPC->AcceptTrack(track)) {
      fHistPtITSTPCsel->Fill(track->Pt());
      if(fCheckSDDIsIn && sddIsIn) fHistPtITSTPCselwSDD->Fill(track->Pt());
      if(fCheckSDDIsIn && !sddIsIn) fHistPtITSTPCselwoSDD->Fill(track->Pt());
      if(isTOFbc0) {
	fHistPtITSTPCselTOFbc0->Fill(track->Pt());
	if(fCheckSDDIsIn && sddIsIn) fHistPtITSTPCselTOFbc0wSDD->Fill(track->Pt());
	if(fCheckSDDIsIn && !sddIsIn) fHistPtITSTPCselTOFbc0woSDD->Fill(track->Pt());
      }
      if(isFake) fHistPtITSTPCselFake->Fill(track->Pt());
      fHistdEdxVSPtITSTPCsel->Fill(track->Pt(),track->GetITSsignal());
      if(isPrimary) {
	fHistPtITSTPCselP->Fill(track->Pt());
	if(pdgTrk==321) fHistPtITSTPCselPfromStrange->Fill(track->Pt());
      } else {
	fHistPtITSTPCselS->Fill(track->Pt());
	if(isFromMat) fHistPtITSTPCselSfromMat->Fill(track->Pt());
	if(isFromStrange) fHistPtITSTPCselSfromStrange->Fill(track->Pt());
      }  
    }
    

    // fill d0 histos
    if(((!fUseITSSAforNtuples&&itsfindableAcc) || fUseITSSAforNtuples) /*&& track->Charge()>0*/) {
      if(track->Pt()>0.150 && track->Pt()<0.200) {
	if(isPrimary) {
	  fHistd0rphiTPCInAccP150200->Fill(d0z0TPC[0]);
	  if(ConditionSPD(track) && itsrefit)
	    fHistd0rphiITSMISPDInAccP150200->Fill(d0z0[0]);
	  if((ConditionSPDone(track)) && itsrefit) {
	    //if(((track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) || (!track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) fHistd0rphiITSMIoneSPDInAccP150200->Fill(d0z0[0]);
	    if(TMath::Abs(d0z0[0])<maxdcaxyITSTPC) fHistd0zITSMIoneSPDInAccP150200->Fill(d0z0[1]);
	  }
	} else {
	  fHistd0rphiTPCInAccS150200->Fill(d0z0TPC[0]);
	  if(ConditionSPD(track) && itsrefit)
	    fHistd0rphiITSMISPDInAccS150200->Fill(d0z0[0]);
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) {
	      fHistd0rphiITSMIoneSPDInAccS150200->Fill(d0z0[0]);
	      if(isFromStrange) fHistd0rphiITSMIoneSPDInAccS150200fromStrange->Fill(d0z0[0]);
	      if(isFromMat) fHistd0rphiITSMIoneSPDInAccS150200fromMat->Fill(d0z0[0]);
	    }
	    if(TMath::Abs(d0z0[0])<maxdcaxyITSTPC) fHistd0zITSMIoneSPDInAccS150200->Fill(d0z0[1]);
	    if(TMath::Abs(d0z0[0])<0.1) fHistPDGMoth150200->Fill(pdgMoth);
	  }
	}
      }
      if(track->Pt()>0.350 && track->Pt()<0.450) {
	if(isPrimary) {
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) fHistd0rphiITSMIoneSPDInAccP350450->Fill(d0z0[0]);
	  }
	} else {
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) {
	      fHistd0rphiITSMIoneSPDInAccS350450->Fill(d0z0[0]);
	      if(isFromStrange) fHistd0rphiITSMIoneSPDInAccS350450fromStrange->Fill(d0z0[0]);
	      if(isFromMat) fHistd0rphiITSMIoneSPDInAccS350450fromMat->Fill(d0z0[0]);
	    }
	  }
	}
      }
      if(track->Pt()>0.500 && track->Pt()<0.700) {
	if(isPrimary) {
	  fHistd0rphiTPCInAccP500700->Fill(d0z0TPC[0]);
	  if(ConditionSPD(track) && itsrefit)
	    fHistd0rphiITSMISPDInAccP500700->Fill(d0z0[0]);
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) fHistd0rphiITSMIoneSPDInAccP500700->Fill(d0z0[0]);
	    if(TMath::Abs(d0z0[0])<maxdcaxyITSTPC) fHistd0zITSMIoneSPDInAccP500700->Fill(d0z0[1]);
	  }
	} else {
	  fHistd0rphiTPCInAccS500700->Fill(d0z0TPC[0]);
	  if(ConditionSPD(track) && itsrefit)
	    fHistd0rphiITSMISPDInAccS500700->Fill(d0z0[0]);
	  if((ConditionSPDone(track)) && itsrefit) {
	    fHistPDGTrk->Fill(pdgTrk);
	    fHistPDGMoth->Fill(pdgMoth);
	    if(TMath::Abs(d0z0[0])<0.1) fHistPDGMoth500700->Fill(pdgMoth);
 	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) {
	      fHistd0rphiITSMIoneSPDInAccS500700->Fill(d0z0[0]);
	      if(isFromStrange) fHistd0rphiITSMIoneSPDInAccS500700fromStrange->Fill(d0z0[0]);
	      if(isFromMat) fHistd0rphiITSMIoneSPDInAccS500700fromMat->Fill(d0z0[0]);
	    }
 	    if(TMath::Abs(d0z0[0])<maxdcaxyITSTPC) fHistd0zITSMIoneSPDInAccS500700->Fill(d0z0[1]);
	    if(pdgMoth==310) fHistd0rphiITSMIoneSPDInAccS500700from310->Fill(d0z0[0]);
	    if(pdgMoth==321) fHistd0rphiITSMIoneSPDInAccS500700from321->Fill(d0z0[0]);
	    if(pdgMoth==3122) fHistd0rphiITSMIoneSPDInAccS500700from3122->Fill(d0z0[0]);
	    if(pdgMoth==211) fHistd0rphiITSMIoneSPDInAccS500700from211->Fill(d0z0[0]);
	    if(pdgMoth==22) fHistd0rphiITSMIoneSPDInAccS500700from22->Fill(d0z0[0]);
	    if(pdgMoth==211 && TMath::Abs(d0z0[0])>0.5) fHistRProdVtxInAccS->Fill(rProdVtx);
	  }
	}
      }
      if(track->Pt()>1.000 && track->Pt()<1.500) {
	if(isPrimary) {
	  fHistd0rphiTPCInAccP10001500->Fill(d0z0TPC[0]);
	  fHistd0zVSetaTPCInAccP10001500->Fill(d0z0TPC[1],trackTPC->Eta());

	  if(ConditionSPD(track) && itsrefit)
	    fHistd0rphiITSMISPDInAccP10001500->Fill(d0z0[0]);
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) fHistd0rphiITSMIoneSPDInAccP10001500->Fill(d0z0[0]);
	    if(TMath::Abs(d0z0[0])<maxdcaxyITSTPC) fHistd0zITSMIoneSPDInAccP10001500->Fill(d0z0[1]);
	    fHistd0rphiVSphiITSMIoneSPDInAccP10001500->Fill(d0z0[0],track->Phi());
	    fHistd0rphiVSetaITSMIoneSPDInAccP10001500->Fill(d0z0[0],track->Eta());
	  }
	} else {
	  fHistd0rphiTPCInAccS10001500->Fill(d0z0TPC[0]);
	  if(ConditionSPD(track) && itsrefit)
	    fHistd0rphiITSMISPDInAccS10001500->Fill(d0z0[0]);
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) {
	      fHistd0rphiITSMIoneSPDInAccS10001500->Fill(d0z0[0]);
	      if(isFromStrange) fHistd0rphiITSMIoneSPDInAccS10001500fromStrange->Fill(d0z0[0]);
	      if(isFromMat) fHistd0rphiITSMIoneSPDInAccS10001500fromMat->Fill(d0z0[0]);
	    }
	    if(TMath::Abs(d0z0[0])<maxdcaxyITSTPC) fHistd0zITSMIoneSPDInAccS10001500->Fill(d0z0[1]);
	    if(TMath::Abs(d0z0[0])<0.1) fHistPDGMoth10001500->Fill(pdgMoth);
	  }
	}
      }
      if(track->Pt()>2.500 && track->Pt()<4.000) {
	if(isPrimary) {
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) fHistd0rphiITSMIoneSPDInAccP25004000->Fill(d0z0[0]);
	  }
	} else {
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) {
	      fHistd0rphiITSMIoneSPDInAccS25004000->Fill(d0z0[0]);
	      if(isFromStrange) fHistd0rphiITSMIoneSPDInAccS25004000fromStrange->Fill(d0z0[0]);
	      if(isFromMat) fHistd0rphiITSMIoneSPDInAccS25004000fromMat->Fill(d0z0[0]);
	    }
	  }
	}
      }
      if(track->Pt()>4.000 && track->Pt()<8.000) {
	if(isPrimary) {
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) fHistd0rphiITSMIoneSPDInAccP40008000->Fill(d0z0[0]);
	  }
	} else {
	  if((ConditionSPDone(track)) && itsrefit) {
	    if(TMath::Abs(d0z0[1])<maxdcazITSTPC) {
	      fHistd0rphiITSMIoneSPDInAccS40008000->Fill(d0z0[0]);
	      if(isFromStrange) fHistd0rphiITSMIoneSPDInAccS40008000fromStrange->Fill(d0z0[0]);
	      if(isFromMat) fHistd0rphiITSMIoneSPDInAccS40008000fromMat->Fill(d0z0[0]);
	    }
	  }
	}
      }

    }

    if(trackTPC) { delete trackTPC; trackTPC=0; }

    // encode ITS cluster map, including MC info
    Int_t iITSflag=MakeITSflag(track);

    // if MC info is available: get particle properties
    Float_t ptMC=-999.,pdgMC=-999.,d0MC=-999.;
    Double_t d0z0MCv[2]={-999.,-999.},covd0z0MCv[3]={1.,1.,1.};
    if(fReadMC) {
      part = mcEvent->Particle(trkLabel);
      if (!part) continue;
      ptMC=part->Pt();
      pdgMC=part->GetPdgCode();
      d0MC=ParticleImpParMC(part,vertexMC,0.1*fESD->GetMagneticField());
      track->PropagateToDCA(vertexMC,fESD->GetMagneticField(),100.,d0z0MCv,covd0z0MCv);
      if(covd0z0MCv[0]<0. || covd0z0MCv[2]<0.) continue;
    }

    Double_t sigmad0MCv=TMath::Sqrt(covd0z0MCv[0]);
    if(!itsrefit) sigmad0MCv *= -1.;
    Float_t isPrimaryFl = (isPrimary ? 1. : 0.);
    Float_t isTPCSelFl  = (itsfindableAcc ? 1. : 0.);

    // fill ntuple with track properties
    if(fFillNtuples && SelectPt(track->Pt())) {
      Float_t fillArray[22]={static_cast<Float_t>(track->Pt()),static_cast<Float_t>(track->Eta()),static_cast<Float_t>(track->Phi()),static_cast<Float_t>(d0z0[0]),static_cast<Float_t>(d0z0[1]),static_cast<Float_t>(TMath::Sqrt(covd0z0[0])),static_cast<Float_t>(TMath::Sqrt(covd0z0[2])),ptMC,static_cast<Float_t>(pdgMC),static_cast<Float_t>(pdgMoth),static_cast<Float_t>(d0MC),static_cast<Float_t>(d0z0MCv[0]),static_cast<Float_t>(d0z0MCv[1]),static_cast<Float_t>(sigmad0MCv),static_cast<Float_t>(TMath::Sqrt(covd0z0MCv[2])),static_cast<Float_t>(iITSflag),static_cast<Float_t>(isPrimaryFl),static_cast<Float_t>(isTPCSelFl)};
      fNtupleESDTracks->Fill(fillArray);
    }

    //---------------------------------------------    
    // AliTrackPoints: alignment checks
    // 
    if(!fFillNtuples) continue;

    const AliTrackPointArray *array = track->GetTrackPointArray();
    if(!array) continue;
    AliTrackPoint point;
    Int_t pointOnLayer[6]={0,0,0,0,0,0};
    Int_t indexAssociated[6]={-1,-1,-1,-1,-1,-1},indexExtra=-1;
    Bool_t extra=kFALSE;
    Int_t layerId,layerExtra=-1;
    for(Int_t ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      Float_t r = TMath::Sqrt(point.GetX()*point.GetX()+point.GetY()*point.GetY());

      if(r>3 && r<6) {
	layerId = 0;
      } else if(r>6 && r<8) {
	layerId = 1;
      } else if(r>8 && r<18) {
	layerId = 2;
      } else if(r>18 && r<30) {
	layerId = 3;
      } else if(r>30 && r<40) {
	layerId = 4;
      } else if(r>40 && r<50) {
	layerId = 5;
      } else {
	layerId=100;
      }

      // only ITS points
      if(layerId>5) continue;

      if(!point.IsExtra()) {
	pointOnLayer[layerId]++;
	indexAssociated[layerId]=ipt;
      } else {
	// this is an extra cluster
	extra=kTRUE;
	layerExtra=layerId;
	indexExtra=ipt;
      }
    } // end loop on AliTrackPoints

    TString vtitle = spdv->GetTitle();
    if(!vtitle.Contains("3D")) continue; 

    // SPD tracklet
    if(indexAssociated[0]>=0 && indexAssociated[1]>=0) {
      AliTrackPoint pointSPD1,pointSPD2;
      array->GetPoint(pointSPD1,indexAssociated[0]);
      array->GetPoint(pointSPD2,indexAssociated[1]);
      Float_t phi=TMath::ATan2(pointSPD2.GetY()-pointSPD1.GetY(),pointSPD2.GetX()-pointSPD1.GetX());
      Float_t lambda=TMath::ATan((pointSPD2.GetZ()-pointSPD1.GetZ())/TMath::Sqrt((pointSPD2.GetX()-pointSPD1.GetX())*(pointSPD2.GetX()-pointSPD1.GetX())+(pointSPD2.GetY()-pointSPD1.GetY())*(pointSPD2.GetY()-pointSPD1.GetY())));
      Float_t theta=0.5*TMath::Pi()-lambda;
      TParticle particle(211,0,0,0,0,0,TMath::Cos(phi),TMath::Sin(phi),TMath::Tan(lambda),10.,pointSPD1.GetX(),pointSPD1.GetY(),pointSPD1.GetZ(),0);
      AliESDtrack tracklet(&particle);
      Float_t dz[2];
      // distance to primary SPD (only if 3D and high multiplicity)
      if(spdv->GetNContributors()>10) { 
	tracklet.GetDZ(spdv->GetX(),spdv->GetY(),spdv->GetZ(),0,dz);
	//tracklet.GetDZ(-0.07,0.25,spdv->GetZ(),0,dz);
	fNtupleITSAlignSPDTracklets->Fill(phi,theta,0.5*(pointSPD1.GetZ()+pointSPD2.GetZ()),dz[0],dz[1],track->Pt());
      }
    }

    // distance to extra
    if(extra && spdv->GetNContributors()>4 && indexAssociated[layerExtra]>-1) {
      AliTrackPoint pointExtra,pointAssociated;
      array->GetPoint(pointAssociated,indexAssociated[layerExtra]);
      array->GetPoint(pointExtra,indexExtra);
      Float_t phiExtra = TMath::ATan2(pointExtra.GetY()-spdv->GetY(),pointExtra.GetX()-spdv->GetX());
      Float_t phiAssociated = TMath::ATan2(pointAssociated.GetY()-spdv->GetY(),pointAssociated.GetX()-spdv->GetX());
      Float_t rExtra = TMath::Sqrt((pointExtra.GetX()-spdv->GetX())*(pointExtra.GetX()-spdv->GetX())+(pointExtra.GetY()-spdv->GetY())*(pointExtra.GetY()-spdv->GetY()));
      Float_t rAssociated = TMath::Sqrt((pointAssociated.GetX()-spdv->GetX())*(pointAssociated.GetX()-spdv->GetX())+(pointAssociated.GetY()-spdv->GetY())*(pointAssociated.GetY()-spdv->GetY()));
      Float_t dzExtra[2];
      dzExtra[0] = (phiExtra-phiAssociated)*0.5*(rExtra+rAssociated);
      dzExtra[1] = pointExtra.GetZ()-pointAssociated.GetZ()-(rExtra-rAssociated)*(pointAssociated.GetZ()-spdv->GetZ())/rAssociated;
      Float_t xlocExtra=-100.,zlocExtra=-100.;
      fNtupleITSAlignExtra->Fill(layerExtra,pointExtra.GetX(),pointExtra.GetY(),pointExtra.GetZ(),dzExtra[0],dzExtra[1],xlocExtra,zlocExtra,nclsITS,track->Pt());  
    }
    
    if(trackTPC) { delete trackTPC; trackTPC=0; }
  } // end loop on tracks
  
  if(vertexMC) { delete vertexMC; vertexMC=0; }

  delete esdtrackCutsTPC;    esdtrackCutsTPC=0;
  delete esdtrackCutsITSTPC; esdtrackCutsITSTPC=0;

  PostData(1, fOutput);

  return;
}      

//________________________________________________________________________
void AliAnalysisTaskITSTrackingCheck::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  fHistNEventsFrac = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEventsFrac"));
  if(fHistNEvents && fHistNEventsFrac) {
    for(Int_t ibin=2; ibin<=fHistNEvents->GetNbinsX(); ibin++) {
      if(fHistNEvents->GetBinContent(2)) fHistNEventsFrac->SetBinContent(ibin,fHistNEvents->GetBinContent(ibin)/fHistNEvents->GetBinContent(2)); 
    }
  }

  return;
}
//---------------------------------------------------------------------------
Int_t AliAnalysisTaskITSTrackingCheck::NumberOfITSClustersMC(Int_t label,Int_t nModules) const
{
  //
  // Return number of ITS clusters produced by MC particle with given label
  //
  
  AliESDInputHandlerRP *esdHRP = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdHRP) return -1;
  TTree *cTree = (TTree*)esdHRP->GetTreeR("ITS");
  if(!cTree) return -1;
  TClonesArray *clusters=0;
  cTree->SetBranchAddress("ITSRecPoints",&clusters);
  if(!clusters) return -1;

  AliITSRecPoint *c=0;
  Int_t i,n,icl,lay,ilab;
  Int_t ncls[7]={0,0,0,0,0,0,0};
  Int_t nclstot=0;

  for(i=0; i<nModules; i++) {
    cTree->GetEvent(i);
    n=clusters->GetEntriesFast();
    for (icl=0; icl<n; icl++) {
      c=(AliITSRecPoint*)clusters->UncheckedAt(icl);
      lay=c->GetLayer();
      //printf("mod %d lay %d xloc %f zloc %f\n",i,lay,c->GetDetLocalX(),c->GetDetLocalZ());
      for(ilab=0;ilab<3;ilab++) {
        if(c->GetLabel(ilab)==label) ncls[lay]++;
      }
    }
  }
  for(i=0;i<fNITSLayers;i++) { 
    //printf("%d\n",ncls[i]); 
    if(ncls[i]) nclstot++; 
  }

  return nclstot;
}
//---------------------------------------------------------------------------
Int_t AliAnalysisTaskITSTrackingCheck::NumberOfITSClusters(Int_t idet,Float_t &xloc) const
{
  //
  // Return number of ITS clusters produced by MC particle with given label
  //
  
  AliESDInputHandlerRP *esdHRP = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdHRP) return -1;
  TTree *cTree = (TTree*)esdHRP->GetTreeR("ITS");
  if(!cTree) return -1;
  TClonesArray *clusters=0;   // new TClonesArray("AliITSRecPoint",10000);
  cTree->SetBranchAddress("ITSRecPoints",&clusters);
  if(!clusters) return -1;

  AliITSRecPoint *c=0;
  Int_t n,icl;

  cTree->GetEvent(idet);
  n=clusters->GetEntriesFast();

  if(n==1) {
    for (icl=0; icl<n; icl++) {
      c=(AliITSRecPoint*)clusters->UncheckedAt(icl);
      xloc = c->GetDetLocalX();
    }
  }
  return n;
}
//---------------------------------------------------------------------------
void AliAnalysisTaskITSTrackingCheck::FillNClustersModuleMap()
{
  //
  // Number of clusters per module
  //
  
  AliESDInputHandlerRP *esdHRP = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdHRP) return;
  TTree *cTree = (TTree*)esdHRP->GetTreeR("ITS");
  if(!cTree) return;
  TClonesArray *clusters=0;   // new TClonesArray("AliITSRecPoint",10000);
  cTree->SetBranchAddress("ITSRecPoints",&clusters);

  for(Int_t idet=0;idet<2198;idet++) {
    cTree->GetEvent(idet);
    if(!clusters) continue;
    fHistNClustersMapModule->SetBinContent(idet,clusters->GetEntriesFast());
  }

  return;
}
//---------------------------------------------------------------------------
Double_t AliAnalysisTaskITSTrackingCheck::ParticleImpParMC(TParticle *part,
							   AliESDVertex *vert,
							   Double_t bzT) const
{
  //
  // Return the MC value of the impact parameter
  //

  bzT = TMath::Abs(bzT);

 
  Double_t vx=part->Vx()-vert->GetX();
  Double_t vy=part->Vy()-vert->GetY();
      
  Double_t pt=part->Pt();     
  Double_t px=part->Px();     
  Double_t py=part->Py();     
  Double_t charge = (part->GetPdgCode()>0. ? 1. : -1.);
  if(TMath::Abs(part->GetPdgCode())<100) charge*=-1.;


  if(TMath::Abs(px)<0.000001) px=0.000001;     
  Double_t rAnd=((10./2.99792458)*pt/bzT)*100.;
  Double_t center[3],d0;
  center[0]=vx-(1./charge)*rAnd*(py/pt);
  center[1]=vy+(1./charge)*rAnd*(px/pt);
  center[2]=TMath::Sqrt(center[0]*center[0]+center[1]*center[1]);
  d0 = -center[2]+rAnd;

  //if(pt>29) printf("%f %f %f %d %f %f %f %f %f\n",px,py,pt,part->GetPdgCode(),center[0],center[1],center[2],rAnd,d0);

  return d0;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskITSTrackingCheck::SelectPt(Double_t pt)
{
  //
  // Keep only tracks in given pt bins
  // 
  Double_t ptlower[11]={0.15,0.29,0.49,0.75,0.9,1.9,3.5,6.5, 9.,19.,27.};
  Double_t ptupper[11]={0.16,0.31,0.51,0.85,1.1,2.1,4.5,7.5,11.,21.,33.};

  for(Int_t i=0; i<11; i++) {
    if(pt>ptlower[i] && pt<ptupper[i]) {
      fCountsPerPtBin[i]++;
      return kTRUE;
    }
  }
  return kFALSE;
  //return kTRUE;
}
//---------------------------------------------------------------------------
Int_t AliAnalysisTaskITSTrackingCheck::MakeITSflag(AliESDtrack *track) const {
  //
  // ITSflag takes the value 0 if the track has no cluster assigned in the SPDs, 
  // 1 (2) if one cluster is assigned in SPD1(2), 3 if two clusters are present. 
  // Then the same adding 10,20 or 30 for SDD and 100,200 or 300 for SSD
  //
  Int_t iITSflag=0;
  if(track->HasPointOnITSLayer(0))             iITSflag+=1;
  if(track->HasPointOnITSLayer(fNITSLayers-5)) iITSflag+=2;
  if(track->HasPointOnITSLayer(fNITSLayers-4)) iITSflag+=10;
  if(track->HasPointOnITSLayer(fNITSLayers-3)) iITSflag+=20;
  if(track->HasPointOnITSLayer(fNITSLayers-2)) iITSflag+=100;
  if(track->HasPointOnITSLayer(fNITSLayers-1)) iITSflag+=200;
  
  // number of associated ITS clusters
  iITSflag += 1000*track->GetNcls(0);
  
  // number of associated TPC clusters
  iITSflag += 100000*track->GetNcls(1);
  
  // if MC info and is available
  // write the number of ITS clusters produced by this track
  Int_t nITSclsMC=0;
  if(fReadMC && fReadRPLabels) {
    nITSclsMC = NumberOfITSClustersMC(TMath::Abs(track->GetLabel()));
    if(nITSclsMC>=0) iITSflag += 10000*nITSclsMC;    
    // flag fake tracks
    if(track->GetLabel()<0) iITSflag *= -1;
  }

  return iITSflag;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskITSTrackingCheck::IsSelectedCentrality() const
{
  //
  // check if events is in the required multiplicity range
  //

  const AliMultiplicity *alimult = fESD->GetMultiplicity();
  //  Int_t ntrklets=1;
  Int_t nclsSPDouter=0;
  if(alimult) {
    //    ntrklets = alimult->GetNumberOfTracklets();
    nclsSPDouter = alimult->GetNumberOfITSClusters(1);
  }

  if(nclsSPDouter<fMinMult || nclsSPDouter>fMaxMult) return kFALSE;


  return kTRUE;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskITSTrackingCheck::ConditionSPD(AliESDtrack *track) const
{
  //
  // check if track has points in all inner barrel layers
  //


  if(NPointsInnerBarrel(track)==(fNITSLayers-4)) return kTRUE;
  return kFALSE;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskITSTrackingCheck::ConditionSPDone(AliESDtrack *track) const
{
  //
  // check if track has points in at least one inner barrel layer
  //


  if(NPointsInnerBarrel(track)>=1) return kTRUE;
  return kFALSE;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskITSTrackingCheck::ConditionSPDtwo(AliESDtrack *track) const
{
  //
  // check if track has points in at least two inner barrel layers
  //


  if(NPointsInnerBarrel(track)>=2) return kTRUE;
  return kFALSE;
}
//---------------------------------------------------------------------------
Int_t AliAnalysisTaskITSTrackingCheck::NPointsInnerBarrel(AliESDtrack *track) const
{
  //
  // N of points in inner barrel layers
  //
  Int_t n=0;

  for(Int_t i=0; i<(fNITSLayers-4); i++) {
    if(track->HasPointOnITSLayer(i)) n++;
  }

  return n;
}
