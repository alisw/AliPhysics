#include <chrono>
#include <ctime>
#include <ratio>
#include <regex>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliNeutralTrackParam.h"
#include "AliESDtrackCuts.h"
#include "AliVertexerTracks.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "DCAFitterN.h"
#include "AliAnalysisTaskHFFindJets.h"
#include "AliLog.h"

#include <TH1F.h>
#include <TSystem.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TObjString.h>


/**************************************************************************
 * Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//
//*************************************************************************

ClassImp(AliAnalysisTaskHFFindJets)
//______________________________________________________________________________
AliAnalysisTaskHFFindJets::AliAnalysisTaskHFFindJets():AliAnalysisTaskSE("HFFindJets"),
	fOutput{nullptr},
	fHistNEvents{nullptr},
	fHistTrackStatus{nullptr},
	fHistPtAllTracks{nullptr},
	fHistPtSelTracks{nullptr},
	fHistTglAllTracks{nullptr},
	fHistTglSelTracks{nullptr},
	fHistEtaAllTracks{nullptr},
	fHistEtaSelTracks2prong{nullptr},
	fHistEtaSelTracks3prong{nullptr},
	fHistImpParAllTracks{nullptr},
	fHistImpParSelTracks2prong{nullptr},
	fHistImpParSelTracks3prong{nullptr},
	fHistITSmapAllTracks{nullptr},
	fHistITSmapSelTracks{nullptr},
	fHistPrimVertX{nullptr},
	fHistPrimVertY{nullptr},
	fHistPrimVertZ{nullptr},
	fHist2ProngVertX{nullptr},
	fHist2ProngVertY{nullptr},
	fHist2ProngVertZ{nullptr},
	fHist3ProngVertX{nullptr},
	fHist3ProngVertY{nullptr},
	fHist3ProngVertZ{nullptr},
	fHistDist12LcpKpi{nullptr},
	fHistInvMassD0{nullptr},
	fHistPtD0{nullptr},
	fHistYPtD0{nullptr},
	fHistPtD0Dau0{nullptr},
	fHistPtD0Dau1{nullptr},
	fHistImpParD0Dau0{nullptr},
	fHistImpParD0Dau1{nullptr},
	fHistd0Timesd0{nullptr},
	fHistCosPointD0{nullptr},
	fHistDecLenD0{nullptr},
	fHistDecLenXYD0{nullptr},
	fHistImpParErrD0Dau{nullptr},
	fHistDecLenErrD0{nullptr},
	fHistDecLenXYErrD0{nullptr},
	fHistCovMatPrimVXX2Prong{nullptr},
	fHistCovMatSecVXX2Prong{nullptr},
	fHistD0SignalVertX{nullptr},
	fHistD0SignalVertY{nullptr},
	fHistD0SignalVertZ{nullptr},
	fHistInvMassD0Signal{nullptr},
	fHistInvMassD0Refl{nullptr},
	fHistInvMassJpsi{nullptr},
	fHistPtJpsi{nullptr},
	fHistPtJpsiDau0{nullptr},
	fHistPtJpsiDau1{nullptr},
	fHistImpParJpsiDau0{nullptr},
	fHistImpParJpsiDau1{nullptr},
	fHistd0Timesd0Jpsi{nullptr},
	fHistCosPointJpsi{nullptr},
	fHistDecLenJpsi{nullptr},
	fHistDecLenXYJpsi{nullptr},
	fHistDecLenErrJpsi{nullptr},
	fHistDecLenXYErrJpsi{nullptr},
	fHistJpsiSignalVertX{nullptr},
	fHistJpsiSignalVertY{nullptr},
	fHistJpsiSignalVertZ{nullptr},
	fHistInvMassJpsiSignal{nullptr},
	fHistInvMassDplus{nullptr},
	fHistInvMassDplusSignal{nullptr},
	fHistPtDplus{nullptr},
	fHistYPtDplus{nullptr},
	fHistPtDplusDau0{nullptr},
	fHistPtDplusDau1{nullptr},
	fHistPtDplusDau2{nullptr},
	fHistImpParDplusDau0{nullptr},
	fHistImpParDplusDau1{nullptr},
	fHistImpParDplusDau2{nullptr},
	fHistDecLenDplus{nullptr},
	fHistDecLenXYDplus{nullptr},
	fHistNormDecLenXYDplus{nullptr},
	fHistImpParErrDplusDau{nullptr},
	fHistDecLenErrDplus{nullptr},
	fHistDecLenXYErrDplus{nullptr},
	fHistCosPointDplus{nullptr},
	fHistCosPointXYDplus{nullptr},
	fHistImpParXYDplus{nullptr},
	fHistNormIPDplus{nullptr},
	fHistoSumSqImpParDplusDau{nullptr},
	fHistCovMatPrimVXX3Prong{nullptr},
	fHistCovMatSecVXX3Prong{nullptr},
	fHistInvMassDs{nullptr},
	fHistInvMassDsSignal{nullptr},
	fHistInvMassDsRefl{nullptr},
	fHistPtDs{nullptr},
	fHistYPtDs{nullptr},
	fHistDecLenDs{nullptr},
	fHistCosPointDs{nullptr},
	fHistInvMassLc{nullptr},
	fHistPtLc{nullptr},
	fHistYPtLc{nullptr},
	fHistPtLcDau0{nullptr},
	fHistPtLcDau1{nullptr},
	fHistPtLcDau2{nullptr},
	fHistDecLenLc{nullptr},
	fHistCosPointLc{nullptr},
	fHistInvMassK0s{nullptr},
	fHistInvMassLcK0sp{nullptr},
	fHistPtLcK0sp{nullptr},
	fHistCPUTimeTrackVsNTracks{nullptr},
	fHistCPUTimeCandVsNTracks{nullptr},
	fHistWallTimeTrackVsNTracks{nullptr},
	fHistWallTimeCandVsNTracks{nullptr},
	fReadMC(kFALSE),
	fUsePhysSel(kTRUE),
	fTriggerMask(AliVEvent::kAny),
	fSelectOnCentrality(kFALSE),
	fMinCentrality(-1.),
	fMaxCentrality(110.),
	fCentrEstimator("V0M"),
	fCutOnSPDVsTrackVtx(kFALSE),
	fMaxZVert(999.),
	fDo3Prong(kFALSE),
	fMaxDecVertRadius2(8),
	fMassDzero(0.),
	fMassJpsi(0.),
	fMassDplus(0.),
	fMassDs(0.),
	fMassLambdaC(0.),
	fMassK0s(0.),
	fSecVertexerAlgo(0),
	fVertexerTracks{nullptr},
	fO2Vertexer2Prong{},
	fO2Vertexer3Prong{},
	fVertexerPropagateToPCA(true),
	fVertexerMaxR(200.),
	fVertexerMaxDZIni(4.),
	fVertexerMinParamChange(1.e-3),
	fVertexerMinRelChi2Change(0.9),
	fVertexerUseAbsDCA(true),
	fTrackCuts2pr{nullptr},
	fTrackCuts3pr{nullptr},
	fTrackCutsBach{nullptr},
	fTrackCutsV0Dau{nullptr},
	fMaxTracksToProcess(9999999),
	fNPtBinsSingleTrack(6),
	fNPtBinsDzero(25),
	fPtWithoutVtxToll(0.1),
	fMinPtDzero(0.),
	fMaxPtDzero(9999.),
	fMinPtDplus(0.),
	fMaxPtDplus(9999.),
	fMinPtJpsi(0.),
	fMaxPtJpsi(9999.),
	fCandidateCutLevel(1),
	fSelectD0(1),
	fSelectD0bar(1),
	fMinPt3Prong(0.),
	fMaxRapidityCand(-999.),
	fNPtBinsDplus(12),
	fNPtBinsJpsi(9),
	fNPtBinsLc(10),
	fSelectDplus(1),
	fSelectJpsi(1),
	fSelectLcpKpi(1),
	fFindVertexForCascades(kFALSE),
	fMinPtV0(0.),
	fMinCosPointV0(0.),
	fCutOnK0sMass(0.1),
	fEnableCPUTimeCheck(kFALSE),
	fCountTimeInMilliseconds(kFALSE),
	hpt_nocuts{nullptr},
	htgl_nocuts{nullptr},
	hpt_cuts{nullptr},
	hdcatoprimxy_cuts{nullptr},
	htgl_cuts{nullptr},
	hvx{nullptr},
	hvy{nullptr},
	hvz{nullptr},
	hvx3{nullptr},
	hvy3{nullptr},
	hvz3{nullptr},
	hitsmap{nullptr},
	hvertexx{nullptr},
	hvertexy{nullptr},
	hvertexz{nullptr},
	hdecayxyz{nullptr},
	hdecayxy{nullptr},
	hmass0{nullptr},
	hmassP{nullptr},
	hptD0{nullptr},
	hptprong0{nullptr},
	hptprong1{nullptr},
	hd0{nullptr},
	hd0d0{nullptr},
	hImpParErr{nullptr},
	hDecLenErr{nullptr},
	hDecLenXYErr{nullptr},
	hCovPVXX{nullptr},
	hCovSVXX{nullptr},
	hjetpt{nullptr},
	hjetE{nullptr},
	hjetpx{nullptr},
	hjetpy{nullptr},
	hjetpz{nullptr}, 
	hjetphi{nullptr},
	hjetrap{nullptr},
	hjetconstituents{nullptr},
	hjetzg{nullptr},
	hjetrg{nullptr},
	hjetnsd{nullptr}

{
  //
  InitDefault();

  for(Int_t i=0; i<5; i++){
    fHistPtGenPrompt[i]=0x0;
    fHistPtGenFeeddw[i]=0x0;
    fHistPtGenLimAccPrompt[i]=0x0;
    fHistPtGenLimAccFeeddw[i]=0x0;
    fHistPtRecoGenPtPrompt[i]=0x0;
    fHistPtRecoGenPtFeeddw[i]=0x0;
    fHistPtRecoPrompt[i]=0x0;
    fHistPtRecoFeeddw[i]=0x0;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskHFFindJets::~AliAnalysisTaskHFFindJets(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;

  if(fOutput && !fOutput->IsOwner()){
	delete fOutput;
	delete fHistNEvents;
	delete hpt_nocuts;
	delete htgl_nocuts;
	delete hpt_cuts;
	delete hdcatoprimxy_cuts;
	delete htgl_cuts;
	delete hvx;
	delete hvy;
	delete hvz;
	delete hvx3;
	delete hvy3;
	delete hvz3;
	delete hitsmap;
	delete hvertexx;
	delete hvertexy;
	delete hvertexz;
	delete hdecayxyz;
	delete hdecayxy;
	delete hmass0;
	delete hmassP;
	delete hptD0;
	delete hptprong0;
	delete hptprong1;
	delete hd0;
	delete hd0d0;
	delete hImpParErr;
	delete hDecLenErr;
	delete hDecLenXYErr;
	delete hCovPVXX;
	delete hCovSVXX;
	delete hjetpt;
	delete hjetE;
	delete hjetpx;
	delete hjetpy;
	delete hjetpz; 
	delete hjetphi;
	delete hjetrap;
	delete hjetconstituents;
	delete hjetzg;
	delete hjetrg;
	delete hjetnsd;
  }
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::InitDefault() {
	/*kbitDplus = 0;
	kbitDs = 1;
	kbitLc = 2;
   
	fMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
	fMassDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
	fMassDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
	fMassLambdaC = TDatabasePDG::Instance()->GetParticle(4122)->Mass();*/
	
    //m     dca   cost* ptk  ptpi  d0k            d0pi         d0d0     cosp cosxy normdxy
	/*Double_t defaultCuts[npTBins][nCutVars] = {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   // pt<0.5
                                     {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   // 0.5<pt<1
                                     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  // 1<pt<1.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  // 1.5<pt<2 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  // 2<pt<2.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  // 2.5<pt<3 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  // 3<pt<3.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  // 3.5<pt<4
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 4<pt<4.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 4.5<pt<5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 5<pt<5.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 5.5<pt<6 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 6<pt<6.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 6.5<pt<7 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   // 7<pt<7.5
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   // 7.5<pt<8 
                                     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   // 8<pt<9 
                                     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   // 9<pt<10 
                                     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   // 10<pt<12 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0.},   // 12<pt<16 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 16<pt<20 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 20<pt<24 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 24<pt<36 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 36<pt<50 
                                     {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0.}}; // pt>50  ;
   for(Int_t ib=0; ib<npTBins; ib++){
		for(Int_t jc=0; jc<nCutVars; jc++){
			fCuts[ib][jc]=defaultCuts[ib][jc];
		}
	}  */
	
	
	/// initialization with default values

  fMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  fMassJpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  fMassDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  fMassDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  fMassLambdaC = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  fMassK0s = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  fTrackCuts2pr = new AliESDtrackCuts("AliESDtrackCuts", "default");
  fTrackCuts2pr->SetPtRange(0., 1.e10);
  // fTrackCuts->SetEtaRange(-0.8, +0.8);
  fTrackCuts2pr->SetMinNClustersTPC(50);
  fTrackCuts2pr->SetRequireITSRefit(kTRUE);
  fTrackCuts2pr->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
  // fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // fTrackCuts->SetMaxDCAToVertexZ(3.2);
  // fTrackCuts->SetMaxDCAToVertexXY(2.4);
  // fTrackCuts->SetDCAToVertex2D(kTRUE);

  // no possibility to apply pT-dependent DCA cut in AliESDtrackCuts --> done by hand
  fTrackCuts2pr->SetMaxDCAToVertexXY(1000.);
  fTrackCuts2pr->SetMinDCAToVertexXY(0.);

  fTrackCuts3pr = new AliESDtrackCuts("AliESDtrackCuts", "default3p");
  fTrackCuts3pr->SetPtRange(0., 1.e10);
  // fTrackCuts->SetEtaRange(-0.8, +0.8);
  fTrackCuts3pr->SetMinNClustersTPC(50);
  fTrackCuts3pr->SetRequireITSRefit(kTRUE);
  fTrackCuts3pr->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);

  // no possibility to apply pT-dependent DCA cut in AliESDtrackCuts --> done by hand
  fTrackCuts3pr->SetMaxDCAToVertexXY(1000.);
  fTrackCuts3pr->SetMinDCAToVertexXY(0.);

  fTrackCutsBach = new AliESDtrackCuts("AliESDtrackCuts", "defaultbach");
  fTrackCutsBach->SetPtRange(0.3, 1.e10);
  fTrackCutsBach->SetEtaRange(-0.8, +0.8);
  fTrackCutsBach->SetMinDCAToVertexXYPtDep("0.*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  fTrackCutsBach->SetMaxDCAToVertexXY(1.);
  
  fTrackCutsV0Dau = new AliESDtrackCuts("AliESDtrackCuts", "defaultV0dau");
  fTrackCutsV0Dau->SetPtRange(0.05, 1.e10);
  fTrackCutsV0Dau->SetEtaRange(-1.1, +1.1);
  fTrackCutsV0Dau->SetMinNCrossedRowsTPC(50);
  fTrackCutsV0Dau->SetRequireTPCRefit(kTRUE);
  
  fNPtBinsSingleTrack = 6;
  Double_t defaultPtBinsSingleTrack[7] = {0., 0.5, 1., 1.5, 2., 3., 1000.};
  for (Int_t ib = 0; ib < fNPtBinsSingleTrack + 1; ib++)
      fPtBinLimsSingleTrack[ib] = defaultPtBinsSingleTrack[ib];

  Double_t defaultSingleTrackCuts[6][kNCutVarsSingleTrack] =
    {{0.0025, 1000.},  /* 0   < pt < 0.5    */
     {0.0025, 1000.},  /* 0.5 < pt < 1      */
     {0.0025, 1000.},  /* 1   < pt < 1.5    */
     {0.0025, 1000.},  /* 1.5 < pt < 2      */
     {0.0000, 1000.},  /* 2   < pt < 3      */
     {0.0000, 1000.}}; /* 3   < pt < 1000   */

  for (Int_t ib = 0; ib < fNPtBinsSingleTrack; ib++) {
    for (Int_t jc = 0; jc < kNCutVarsSingleTrack; jc++) {
      fSingleTrackCuts2Prong[ib][jc] = defaultSingleTrackCuts[ib][jc];
      fSingleTrackCuts3Prong[ib][jc] = defaultSingleTrackCuts[ib][jc];
    }
  }

  Double_t defaultPtBins2Prongs[3] = {0., 5., 1000.};
  Double_t defaultPtBins3Prongs[3] = {0., 5., 1000.};
  for (Int_t ib = 0; ib < 2 + 1; ib++) {
    fPtBinLimsDzeroSkims[ib] = defaultPtBins2Prongs[ib];
    fPtBinLimsJpsiSkims[ib] = defaultPtBins2Prongs[ib];
    fPtBinLimsDplusSkims[ib] = defaultPtBins3Prongs[ib];
    fPtBinLimsDsSkims[ib] = defaultPtBins3Prongs[ib];
    fPtBinLimsLcSkims[ib] = defaultPtBins3Prongs[ib];
    fPtBinLimsXicSkims[ib] = defaultPtBins3Prongs[ib];
  }

  Double_t defaultDzeroSkimCuts[2][kNCutVars2Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultJpsiSkimCuts[2][kNCutVars2Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultDplusSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultDsSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultLcSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultXicSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  for(Int_t ib=0; ib<fNPtBins; ib++){
    for(Int_t jc=0; jc<kNCutVars2Prong; jc++){
      fDzeroSkimCuts[ib][jc]=defaultDzeroSkimCuts[ib][jc];
      fJpsiSkimCuts[ib][jc]=defaultJpsiSkimCuts[ib][jc];
      fDplusSkimCuts[ib][jc]=defaultDplusSkimCuts[ib][jc];
      fDsSkimCuts[ib][jc]=defaultDsSkimCuts[ib][jc];
      fLcSkimCuts[ib][jc]=defaultLcSkimCuts[ib][jc];
      fXicSkimCuts[ib][jc]=defaultXicSkimCuts[ib][jc];
    }
  }

  fNPtBinsDzero=25;
  Double_t defaultPtBins[26] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 16.0, 20.0, 24.0, 36.0, 50.0, 100.0};
  for(Int_t ib=0; ib<fNPtBinsDzero+1; ib++) fPtBinLimsDzero[ib]=defaultPtBins[ib];

  Double_t defaultD0Cuts[25][kNCutVarsDzero] =
    {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   /* pt<0.5*/
     {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   /* 0.5<pt<1*/
     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  /* 1<pt<1.5 */
     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  /* 1.5<pt<2 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  /* 2<pt<2.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  /* 2.5<pt<3 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  /* 3<pt<3.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  /* 3.5<pt<4 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 4<pt<4.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 4.5<pt<5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 5<pt<5.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 5.5<pt<6 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 6<pt<6.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 6.5<pt<7 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   /* 7<pt<7.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   /* 7.5<pt<8 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 8<pt<9 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 9<pt<10 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 10<pt<12 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0.},   /* 12<pt<16 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 16<pt<20 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 20<pt<24 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 24<pt<36 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 36<pt<50 */
     {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0.}}; /* pt>50 */
  for(Int_t ib=0; ib<fNPtBinsDzero; ib++){
   for(Int_t jc=0; jc<kNCutVarsDzero; jc++){
     fDzeroCuts[ib][jc]=defaultD0Cuts[ib][jc];
   }
  }

  fNPtBinsDplus = 12;
  Double_t defaultPtBinsDplus[13] = {1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36.};
  for (Int_t ib = 0; ib < fNPtBinsDplus + 1; ib++)
    fPtBinLimsDplus[ib] = defaultPtBinsDplus[ib];

  Double_t defaultDplusCuts[12][kNCutVarsDplus] =
    {{0.2, 0.3, 0.3, 0.07, 6., 0.96, 0.985, 2.5},  /* 1<pt<2   */
     {0.2, 0.3, 0.3, 0.07, 5., 0.96, 0.985, 2.5},  /* 2<pt<3   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.980, 2.5},  /* 3<pt<4   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 4<pt<5   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 5<pt<6   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 6<pt<7   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 7<pt<8   */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 8<pt<10  */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 10<pt<12 */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 12<pt<16 */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 16<pt<24 */
     {0.2, 0.3, 0.3, 0.20, 5., 0.94, 0.000, 2.5}}; /* 24<pt<36 */
  for(Int_t ib=0; ib<fNPtBinsDplus; ib++){
   for(Int_t jc=0; jc<kNCutVarsDplus; jc++){
     fDplusCuts[ib][jc]=defaultDplusCuts[ib][jc];
   }
  }

  fNPtBinsLc = 10;
  Double_t defaultPtBinsLc[11] = {0., 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 24.0, 36.0};
  for (Int_t ib = 0; ib < fNPtBinsLc + 1; ib++)
    fPtBinLimsLc[ib] = defaultPtBinsLc[ib];

  Double_t defaultLcCuts[10][kNCutVarsLc] =
    {{0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* pt<1*/
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 1<pt<2*/
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 2<pt<3 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 3<pt<4 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 4<pt<5 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 5<pt<6 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 6<pt<8 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 8<pt<12 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 12<pt<24 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.}}; /* 24<pt<36 */
  for (Int_t ib = 0; ib < fNPtBinsLc; ib++) {
    for (Int_t jc = 0; jc < kNCutVarsLc; jc++) {
      fLcCuts[ib][jc] = defaultLcCuts[ib][jc];
    }
  }

  fNPtBinsJpsi = 9;
  Double_t defaultPtBinsJpsi[10] = {0, 0.5, 1., 2., 3., 4., 5., 7., 10., 15.};
  for (Int_t ib = 0; ib < fNPtBinsJpsi + 1; ib++)
      fPtBinLimsJpsi[ib] = defaultPtBinsJpsi[ib];

  Double_t defaultJpsiCuts[9][kNCutVarsJpsi] =
    {{0.5, 0.2, 0.4, 1},  /* pt<0.5    */
    {0.5, 0.2, 0.4,  1},  /* 0.5<pt<1   */
    {0.5, 0.2, 0.4,  1},  /* 1<pt<2   */
    {0.5, 0.2, 0.4,  1},  /* 2<pt<3   */
    {0.5, 0.2, 0.4,  1},  /* 3<pt<4   */
    {0.5, 0.2, 0.4,  1},  /* 4<pt<5   */
    {0.5, 0.2, 0.4,  1},  /* 5<pt<7   */
    {0.5, 0.2, 0.4,  1},  /* 7<pt<10  */
    {0.5, 0.2, 0.4,  1}}; /* 10<pt<15 */

  for (Int_t ib = 0; ib < fNPtBinsJpsi; ib++) {
    for (Int_t jc = 0; jc < kNCutVarsJpsi; jc++) {
      fJpsiCuts[ib][jc] = defaultJpsiCuts[ib][jc];
    }
  }                              
                                                        	
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::ReadJson() {
	std::string value = GetJsonString("dpl-config_std.json", "aod-file");
	printf("%s\n", value);
	int i3p = GetJsonInteger("dpl-config_std.json", "do3prong");
	printf("%d\n", i3p);
	int tin = GetJsonInteger("dpl-config_std.json", "triggerindex");
	printf("%d\n", tin);
	float minpt = GetJsonFloat("dpl-config_std.json", "ptmintrack");
	printf("%f\n", minpt);
	float dcatoprimxymin = GetJsonFloat("dpl-config_std.json", "dcatoprimxymin");
	printf("%f\n", dcatoprimxymin);
	bool doit = GetJsonBool("dpl-config_std.json", "b_propdca");
	printf("%d\n", doit);
	float chi2 = GetJsonFloat("dpl-config_std.json", "d_minrelchi2change");
	printf("%f\n", chi2);
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::UserCreateOutputObjects() {  
	// create output histos
	fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos1");
    
    hpt_nocuts = new TH1F("hpt_nocuts", " ; pt tracks (#GeV) ; Entries", 100, 0, 10.);
    fOutput->Add(hpt_nocuts);
	htgl_nocuts = new TH1F("htgl_nocuts", "tgl tracks (#GeV)", 100, -5., 5.);
	fOutput->Add(htgl_nocuts);
	hpt_cuts = new TH1F("hpt_cuts", " ; pt tracks (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hpt_cuts);
	hdcatoprimxy_cuts = new TH1F("hdcatoprimxy_cuts", "dca xy to prim. vtx (cm)", 100, -1.0, 1.0);
	fOutput->Add(hdcatoprimxy_cuts);
	htgl_cuts = new TH1F("htgl_cuts", "tgl tracks (#GeV)", 100, -5., 5.);
	fOutput->Add(htgl_cuts);
	hvx = new TH1F("hvx", " Secondary vertex ; X vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvx);
	hvy = new TH1F("hvy", " Secondary vertex ; Y vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvy);
	hvz = new TH1F("hvz", " Secondary vertex ; Z vertex (cm) ; Entries", 1000, -20.0, 20.0);
	fOutput->Add(hvz);
	hvx3 = new TH1F("hvx3", " Secondary vertex 3prong ; X vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvx3);
	hvy3 = new TH1F("hvy3", " Secondary vertex 3prong ; Y vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvy3);
	hvz3 = new TH1F("hvz3", " Secondary vertex 3prong ; Z vertex (cm) ; Entries", 1000, -20.0, 20.0);
	fOutput->Add(hvz3);
	hitsmap = new TH1F("hitsmap", "hitsmap_cuts", 64, -0.5, 63.5);
	fOutput->Add(hitsmap);

	hvertexx = new TH1F("hvertexx", " Primary vertex ; X vertex (cm) ; Entries", 100, -0.5, 0.5);
	fOutput->Add(hvertexx);
	hvertexy = new TH1F("hvertexy", " Primary vertex ; Y vertex (cm) ; Entries", 100, -0.5, 0.5);
	fOutput->Add(hvertexy);
	hvertexz = new TH1F("hvertexz", " Primary vertex ; Z vertex (cm) ; Entries", 100, -20.0, 20.0);
	fOutput->Add(hvertexz);

	hdecayxyz = new TH1F("hdecayxyz", "hdecayxyz", 200, 0., 2.0);
	fOutput->Add(hdecayxyz);
	hdecayxy = new TH1F("hdecayxy", "hdecayxy", 200, 0., 2.0);
	fOutput->Add(hdecayxy);
	hmass0 = new TH1F("hmass0", "; Inv Mass (GeV/c^{2})", 500, 0, 5.0);
	fOutput->Add(hmass0);
	hmassP = new TH1F("hmassP", "; Inv Mass (GeV/c^{2})", 500, 1.6, 2.1);
	fOutput->Add(hmassP);
	hptD0 = new TH1F("hptD0", " ; pt D0 (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hptD0);
	hptprong0 = new TH1F("hptprong0", " ; pt prong0 (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hptprong0);
	hptprong1 = new TH1F("hptprong1", " ; pt prong1 (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hptprong1);
	hd0 = new TH1F("hd0", "dca xy to prim. vertex (cm)", 100, -1.0, 1.0);
	fOutput->Add(hd0);
	hd0d0 = new TH1F("hd0d0", "product of dca xy to prim. vertex (cm^{2})", 500, -1.0, 1.0);
	fOutput->Add(hd0d0);
	hImpParErr = new TH1F("hImpParErr", "impact parameter error", 100, -1.0, 1.0);
	fOutput->Add(hImpParErr);
	hDecLenErr = new TH1F("hDecLenErr", "decay length error", 100, 0., 1.0);
	fOutput->Add(hDecLenErr);
	hDecLenXYErr = new TH1F("hDecLenXYErr", "decay length XY error", 100, 0., 1.0);
	fOutput->Add(hDecLenXYErr);
	hCovPVXX = new TH1F("hCovPVXX", "XX element of PV cov. matrix", 100, 0., 1.0e-4);
	fOutput->Add(hCovPVXX);
	hCovSVXX = new TH1F("hCovSVXX", "XX element of SV cov. matrix", 100, 0., 0.2);
	fOutput->Add(hCovSVXX);

	hjetpt = new TH1F("hjetpt", " ; pt jet (#GeV) ; Entries", 100, 0., 100.);
	fOutput->Add(hjetpt);
	hjetE = new TH1F("hjetE", " ; jet E ; Entries", 100, 0., 100.);
	fOutput->Add(hjetE);
	hjetpx = new TH1F("hjetpx", " ; px jet (#GeV) ; Entries", 100, 0., 100.);
	fOutput->Add(hjetpx);
	hjetpy = new TH1F("hjetpy", " ; py jet (#GeV) ; Entries", 100, 0., 100.);
	fOutput->Add(hjetpy);
	hjetpz = new TH1F("hjetpz", " ; pz jet (#GeV) ; Entries", 100, 0., 100.);
	fOutput->Add(hjetpz);
	hjetphi = new TH1F("hjetphi", " ; jet phi ; Entries", 100, 0., 100.);
	fOutput->Add(hjetphi);
	hjetrap = new TH1F("hjetrap", " ; jet rapidity ; Entries", 100, 0., 100.);
	fOutput->Add(hjetrap);
	hjetconstituents = new TH1F("hjetconstituents", " ; jet constituents ; Entries", 100, 0., 100.);
	fOutput->Add(hjetconstituents);
	hjetzg = new TH1F("hjetzg", " ; jet zg ; Entries", 100, 0., 100.);
	fOutput->Add(hjetzg);
	hjetrg = new TH1F("hjetrrg", " ; jet rg ; Entries", 100, 0., 100.);
	fOutput->Add(hjetrg);
	hjetnsd = new TH1F("hjetnsd", " ; jet nsd ; Entries", 100, 0., 100.);
	fOutput->Add(hjetnsd);

	PostData(1,fOutput);
}

//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHFFindJets::SingleTrkCutsSimple(AliESDtrack* trk, Int_t minclutpc, int ptmintrack, double dcatoprimxymin, AliESDVertex* fV1, Double_t fBzkG) {
	Int_t status = trk->GetStatus();
	bool sel_track = status & AliESDtrack::kITSrefit && (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1));
	sel_track = sel_track && trk->GetNcls(1) >= minclutpc;
	sel_track = sel_track && trk->Pt() > ptmintrack;
	AliExternalTrackParam* track = (AliExternalTrackParam*)trk;
	double b[2];
	double bCov[3];
	track->PropagateToDCA(fV1, fBzkG, 100., b, bCov);
	sel_track = sel_track && abs(b[0]) > dcatoprimxymin;
	return sel_track;
}

//_______________________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::TwoProngSelectionCuts(AliAODRecoDecayHF2Prong* cand, Double_t candpTMin, Double_t candpTMax) {
	bool isD0 = true;
	bool isD0bar = true;
	Double_t candpT = cand->Pt();
	if (candpT < candpTMin || candpT >= candpTMax) return 0;
	Int_t pTBin = GetPtBin(candpT, fPtBinLimsDzeroSkims, kMaxNPtBins2ProngsSkims);
	if (pTBin==-1) return 0;
	if (cand->Prodd0d0() > fDzeroCuts[pTBin][7]) return 0;
	if (cand->CosPointingAngle() < fDzeroCuts[pTBin][8]) return 0;
	if (cand->CosPointingAngleXY() < fDzeroCuts[pTBin][9]) return 0;
	if (cand->NormalizedDecayLengthXY() < fDzeroCuts[pTBin][10]) return 0;
	Double_t decayLengthCut = TMath::Min((cand->P() * 0.0066) + 0.01, 0.06);
	if (TMath::Abs(cand->Normalizedd0Prong(0)) < 0.5 || TMath::Abs(cand->Normalizedd0Prong(1)) < 0.5) return 0;
	if (cand->DecayLength() * cand->DecayLength() < decayLengthCut * decayLengthCut) return 0;
	// if (cand->NormalizedDecayLength() * cand->NormalizedDecayLength() < 1.0) return 0;
	if (TMath::Abs(cand->InvMassD0()-fMassDzero) > fDzeroCuts[pTBin][0] ) isD0=false;
	if (TMath::Abs(cand->InvMassD0bar()-fMassDzero) > fDzeroCuts[pTBin][0] ) isD0bar=false;
	if (!isD0 && !isD0bar) return 0;

	if (cand->Pt2Prong(0) < fDzeroCuts[pTBin][4]*fDzeroCuts[pTBin][4] || cand->Pt2Prong(1) < fDzeroCuts[pTBin][3]*fDzeroCuts[pTBin][3] ) isD0=false;
	if (cand->Pt2Prong(0) < fDzeroCuts[pTBin][3]*fDzeroCuts[pTBin][3] || cand->Pt2Prong(1) < fDzeroCuts[pTBin][4]*fDzeroCuts[pTBin][4] ) isD0bar=false;
	if (!isD0 && !isD0bar) return 0;

	if (TMath::Abs(cand->Getd0Prong(0)) > fDzeroCuts[pTBin][6] || TMath::Abs(cand->Getd0Prong(1)) > fDzeroCuts[pTBin][5] ) isD0=false;
	if (TMath::Abs(cand->Getd0Prong(0)) > fDzeroCuts[pTBin][5] || TMath::Abs(cand->Getd0Prong(1)) > fDzeroCuts[pTBin][6] ) isD0bar=false;
	if (!isD0 && !isD0bar) return 0;

	Double_t cosThetaStarD0,cosThetaStarD0bar;
	cand->CosThetaStarD0(cosThetaStarD0,cosThetaStarD0bar);
	if (TMath::Abs(cosThetaStarD0) > fDzeroCuts[pTBin][2] ) isD0=false;
	if (TMath::Abs(cosThetaStarD0bar) > fDzeroCuts[pTBin][2] ) isD0bar=false;
	if (!isD0 && !isD0bar) return 0;

	Int_t returnValue=0;
	if(isD0) returnValue+=1;
	if(isD0bar) returnValue+=2;
	return returnValue;
}

//_______________________________________________________________________________________
/*AliESDVertex* AliAnalysisTaskHFFindJets::ReconstructSecondaryVertex(AliVertexerTracks* vt, TObjArray* trkArray, AliESDVertex* primvtx, double rmax) {
	vt->SetVtxStart(primvtx);

	AliESDVertex* trkv = (AliESDVertex*)vt->VertexForSelectedESDTracks(trkArray);
	if (trkv->GetNContributors() != trkArray->GetEntriesFast())
		return 0x0;
	Double_t vertRadius2 = trkv->GetX() * trkv->GetX() + trkv->GetY() * trkv->GetY();
	if(vertRadius2>rmax*rmax) return 0x0;
	return trkv;
}*/

//_______________________________________________________________________________________
AliESDVertex* AliAnalysisTaskHFFindJets::ReconstructSecondaryVertex(TObjArray* trkArray, AliESDVertex* primvtx)
{

  AliESDVertex* trkv =0x0;
  // printf("------\n");
  // for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
  //   AliExternalTrackParam *tp=(AliExternalTrackParam *)trkArray->At(jt);
  //   printf("%f ",tp->Pt());
  // }
  // printf("\n");
  if(fSecVertexerAlgo==0){
    fVertexerTracks->SetVtxStart(primvtx);
    trkv = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(trkArray);
    if (trkv->GetNContributors() != trkArray->GetEntriesFast()) return 0x0;
  }else if(fSecVertexerAlgo==1){
    o2::track::TrackParCov* o2Track[3] = {nullptr};
    for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
      o2Track[jt] = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)trkArray->At(jt));
    }
    Int_t nVert=0;
    Double_t vertCoord[3];
    Double_t vertCov[6];
    Double_t vertChi2=-999.;
    if(trkArray->GetEntriesFast()==2){
      nVert=fO2Vertexer2Prong.process(*o2Track[0], *o2Track[1]);
      if(nVert){
        fO2Vertexer2Prong.propagateTracksToVertex();
        auto vertPos = fO2Vertexer2Prong.getPCACandidate();
        auto vertCMat = fO2Vertexer2Prong.calcPCACovMatrix().Array();
        for(Int_t ic=0; ic<3; ic++) vertCoord[ic]=vertPos[ic];
        for(Int_t ic=0; ic<6; ic++) vertCov[ic]=vertCMat[ic];
        vertChi2 = fO2Vertexer2Prong.getChi2AtPCACandidate();
      }
    }else if(trkArray->GetEntriesFast()==3){
      nVert=fO2Vertexer3Prong.process(*o2Track[0], *o2Track[1], *o2Track[2]);
      if(nVert){
        fO2Vertexer3Prong.propagateTracksToVertex();
        auto vertPos = fO2Vertexer3Prong.getPCACandidate();
        auto vertCMat = fO2Vertexer3Prong.calcPCACovMatrix().Array();
        for(Int_t ic=0; ic<3; ic++) vertCoord[ic]=vertPos[ic];
        for(Int_t ic=0; ic<6; ic++) vertCov[ic]=vertCMat[ic];
        vertChi2 = fO2Vertexer3Prong.getChi2AtPCACandidate();
      }
    }
    if(nVert) trkv = new AliESDVertex(vertCoord,vertCov,vertChi2,trkArray->GetEntriesFast());
    //   else printf("nVert = %d\n",nVert);
  }else if(fSecVertexerAlgo==2){
    KFParticle  dMesonVert;
    double posmom[6],cov[21];
    for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
      AliExternalTrackParam* trpar=(AliExternalTrackParam*)trkArray->At(jt);
      trpar->GetXYZ(posmom);
      trpar->GetPxPyPz(posmom+3);
      trpar->GetCovarianceXYZPxPyPz(cov);
      KFParticle trKFpar;
      trKFpar.Create(posmom,cov,trpar->Charge(),0.13957); // mass of the pion for the time being
      dMesonVert.AddDaughter(trKFpar);
    }
    Double_t vertChi2 = dMesonVert.GetChi2() / dMesonVert.GetNDF();
    Double_t vertCoord[3]={dMesonVert.X(),dMesonVert.Y(),dMesonVert.Z()};
    Double_t vertCov[6];
    vertCov[0]=dMesonVert.GetCovariance(0,0);
    vertCov[1]=dMesonVert.GetCovariance(0,1);
    vertCov[2]=dMesonVert.GetCovariance(1,1);
    vertCov[3]=dMesonVert.GetCovariance(0,2);
    vertCov[4]=dMesonVert.GetCovariance(1,2);
    vertCov[5]=dMesonVert.GetCovariance(2,2);
    trkv = new AliESDVertex(vertCoord,vertCov,vertChi2,trkArray->GetEntriesFast());
  }
  if(!trkv) return 0x0;
  Double_t vertRadius2 = trkv->GetX() * trkv->GetX() + trkv->GetY() * trkv->GetY();
  if(vertRadius2>fMaxDecVertRadius2) return 0x0;
  //  trkv->Print("all");
  //  printf("=============\n");
  return trkv;
}

//_______________________________________________________________________________________
AliAODVertex* AliAnalysisTaskHFFindJets::ConvertToAODVertex(AliESDVertex* trkv)
{
  Double_t pos[3], cov[6], chi2perNDF;
  trkv->GetXYZ(pos);       // position
  trkv->GetCovMatrix(cov); // covariance matrix
  chi2perNDF = trkv->GetChi2toNDF();
  //  double dispersion = trkv->GetDispersion();
  AliAODVertex* vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, trkv->GetNContributors());
  vertexAOD->SetNContributors(trkv->GetNContributors());
  return vertexAOD;
}

//_______________________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3)
{

  Int_t retval = (1 << kbitDplus) + (1 << kbitDs) + (1 << kbitLc);
  Double_t momentum[3];
  Double_t px[3], py[3], pz[3];
  for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
    AliESDtrack* track = (AliESDtrack*)trkArray->UncheckedAt(iTrack);
    track->GetPxPyPz(momentum);
    px[iTrack] = momentum[0];
    py[iTrack] = momentum[1];
    pz[iTrack] = momentum[2];
  }
  UInt_t pdg3[3];
  Int_t nprongs = 3;
  rd4massCalc3->SetPxPyPzProngs(nprongs, px, py, pz);
  Double_t ptCand = rd4massCalc3->Pt() + fPtWithoutVtxToll;
  Int_t iPtBinDplus = GetPtBin(ptCand, fPtBinLimsDplusSkims, kMaxNPtBins3ProngsSkims);
  if(iPtBinDplus < 0) {
    retval &= ~(1 << kbitDplus);
  }
  Int_t iPtBinDs = GetPtBin(ptCand, fPtBinLimsDsSkims, kMaxNPtBins3ProngsSkims);
  if(iPtBinDs < 0) {
    retval &= ~(1 << kbitDs);
  }
  Int_t iPtBinLc = GetPtBin(ptCand, fPtBinLimsLcSkims, kMaxNPtBins3ProngsSkims);
  if(iPtBinLc < 0) {
    retval &= ~(1 << kbitLc);
  }
  Double_t minv2;
  Double_t lolim=fDplusSkimCuts[iPtBinDplus][0];
  Double_t hilim=fDplusSkimCuts[iPtBinDplus][1];
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if ((retval & (1 << kbitDplus)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
    retval &= ~(1 << kbitDplus);
  lolim=fDsSkimCuts[iPtBinDs][0];
  hilim=fDsSkimCuts[iPtBinDs][1];
  Bool_t isSelDs[2] = {true, true}; 
  for (Int_t ih = 0; ih < 2; ih++) {
    Int_t k = ih * 2;
    pdg3[k] = 321;
    pdg3[1] = 321;
    pdg3[2 - k] = 211;
    minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
    if ((retval & (1 << kbitDs)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
      isSelDs[ih] = false;
  }
  if(!isSelDs[0] && !isSelDs[1])
   retval &= ~(1 << kbitDs);
  lolim=fLcSkimCuts[iPtBinLc][0];
  hilim=fLcSkimCuts[iPtBinLc][1];
  Bool_t isSelLc[2] = {true, true}; 
  pdg3[0] = 2212;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if ((retval & (1 << kbitLc)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
    isSelLc[0] = false;
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 2212;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if ((retval & (1 << kbitLc)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
    isSelLc[1] = false;
  if(!isSelLc[0] && !isSelLc[1])
   retval &= ~(1 << kbitLc);

  return retval;
}

//______________________________________________________________________________
AliAODRecoDecayHF2Prong* AliAnalysisTaskHFFindJets::Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)twoTrackArray->UncheckedAt(1);

  Double_t px[2], py[2], pz[2], d0[2], d0err[2];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, bzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, bzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];

  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, bzkG, xdummy, ydummy);

  AliAODRecoDecayHF2Prong* the2Prong = new AliAODRecoDecayHF2Prong(0x0, px, py, pz, d0, d0err, dcap1n1);
  AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  the2Prong->SetOwnSecondaryVtx(ownsecv);
  return the2Prong;
}

//______________________________________________________________________________
AliAODRecoDecayHF3Prong* AliAnalysisTaskHFFindJets::Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)threeTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)threeTrackArray->UncheckedAt(1);
  AliESDtrack* track_2 = (AliESDtrack*)threeTrackArray->UncheckedAt(2);

  Double_t px[3], py[3], pz[3], d0[3], d0err[3];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, bzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, bzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];
  GetTrackMomentumAtSecVert(track_2, secVert, momentum, bzkG);
  px[2] = momentum[0];
  py[2] = momentum[1];
  pz[2] = momentum[2];
  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);
  track_2->GetImpactParameters(d0z0f, covd0z0f);
  d0[2] = d0z0f[0];
  d0err[2] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, bzkG, xdummy, ydummy);
  float dcap2n1 = track_2->GetDCA(track_1, bzkG, xdummy, ydummy);
  float dcap1p2 = track_0->GetDCA(track_2, bzkG, xdummy, ydummy);
  Double_t dca[3] = {dcap1n1, dcap2n1, dcap1p2};
  Double_t dispersion = 0;
  Double_t dist12 = 0.;
  Double_t dist23 = 0.;
  Short_t charge = (Short_t)(track_0->Charge() + track_1->Charge() + track_2->Charge());

  // construct the candidate passing a NULL pointer for the secondary vertex to avoid creation of TRef
  AliAODRecoDecayHF3Prong* the3Prong = new AliAODRecoDecayHF3Prong(0x0, px, py, pz, d0, d0err, dca, dispersion, dist12, dist23, charge);
  // add a pointer to the secondary vertex via SetOwnSecondaryVtx (no TRef created)
  AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  the3Prong->SetOwnSecondaryVtx(ownsecv);
  return the3Prong;
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::InitFromJson(TString filename){
  /// read configuration from json file
  if (filename != "" && gSystem->Exec(Form("ls %s > /dev/null", filename.Data())) == 0) {
    printf("------Read configuration from JSON file------\n");

    std::string triggerMaskFromJSON = GetJsonString(filename.Data(), "triggerClassName");

    if (triggerMask.find(triggerMaskFromJSON) != triggerMask.end()) {
      fUsePhysSel = kTRUE;
      fTriggerMask = triggerMask[triggerMaskFromJSON];
    }

    Double_t ptmintrack2 = GetJsonFloat(filename.Data(), "pTMinTrack2Prong");
    printf("Min pt track (2 prong)= %f\n", ptmintrack2);
    if(ptmintrack2>0) fTrackCuts2pr->SetPtRange(ptmintrack2, 1.e10);
    Double_t ptmintrack3 = GetJsonFloat(filename.Data(), "pTMinTrack3Prong");
    printf("Min pt track (3 prong)= %f\n", ptmintrack3);
    if(ptmintrack3>0) fTrackCuts3pr->SetPtRange(ptmintrack3, 1.e10);
    Int_t do3Prongs = GetJsonInteger(filename.Data(), "do3prong");
    printf("do3prong     = %d\n", do3Prongs);
    if(do3Prongs>0) fDo3Prong=kTRUE;
    Int_t selectD0 = GetJsonInteger(filename.Data(), "d_selectionFlagD0");
    printf("d_selectionFlagD0 = %d\n",selectD0);
    if(selectD0>=0) fSelectD0=selectD0;
    Int_t selectD0bar = GetJsonInteger(filename.Data(), "d_selectionFlagD0bar");
    printf("d_selectionFlagD0bar = %d\n",selectD0bar);
    if(selectD0>=0) fSelectD0bar=selectD0bar;
    Int_t selectDplus = GetJsonInteger(filename.Data(), "d_selectionFlagDPlus");
    printf("d_selectionFlagDplus = %d\n",selectDplus);
    if(selectDplus>=0) fSelectDplus=selectDplus;
    Int_t selectJpsi = GetJsonInteger(filename.Data(), "d_selectionFlagJpsi");
    printf("d_selectionFlagJpsi = %d\n",selectJpsi);
    if(selectJpsi>=0) fSelectJpsi=selectJpsi;
    Int_t selectLcpKpi = GetJsonInteger(filename.Data(), "d_selectionFlagLc");
    printf("d_selectionFlagLc = %d\n", selectLcpKpi);
    if (selectLcpKpi >= 0) fSelectLcpKpi=selectLcpKpi;
    Int_t minncluTPC = GetJsonInteger(filename.Data(), "tpcNClsFound");
    if(minncluTPC>0) printf("minncluTPC   = %d\n", minncluTPC);
    fTrackCuts2pr->SetMinNClustersTPC(minncluTPC);
    fTrackCuts3pr->SetMinNClustersTPC(minncluTPC);

    int nptbinlimsSingleTrack = 0;
    float* ptbinsSingleTrack = GetJsonArray(filename.Data(),"pTBinsTrack",nptbinlimsSingleTrack);
    int npt2Prong = 0, npt3Prong = 0, nc2Prong = 0, nc3Prong = 0;
    float** cutsSingleTrack2Prong = GetJsonMatrix(filename.Data(),"cutsTrack2Prong",npt2Prong,nc2Prong);
    float** cutsSingleTrack3Prong = GetJsonMatrix(filename.Data(),"cutsTrack3Prong",npt3Prong,nc3Prong);
    if((nptbinlimsSingleTrack-1 != npt3Prong) || (nptbinlimsSingleTrack-1 != npt2Prong))
      AliFatal("Number of pT bins in JSON for single track cuts of 2-prong and 3-prongs not consistent, please check it");

    for (Int_t ib = 0; ib < nptbinlimsSingleTrack; ib++) {
      fPtBinLimsSingleTrack[ib] = ptbinsSingleTrack[ib];
    }
    for (Int_t ib = 0; ib < nptbinlimsSingleTrack-1; ib++) {
      for (Int_t jc = 0; jc < nc2Prong; jc++) {
        fSingleTrackCuts2Prong[ib][jc] = cutsSingleTrack2Prong[ib][jc];
        AliWarning(Form("2prong %d, %d, %f", ib, jc, cutsSingleTrack2Prong[ib][jc]));
      }
      for (Int_t jc = 0; jc < nc3Prong; jc++) {
        fSingleTrackCuts3Prong[ib][jc] = cutsSingleTrack3Prong[ib][jc];
        AliWarning(Form("3prong %d %d, %d, %f", nc2Prong, ib, jc, cutsSingleTrack2Prong[ib][jc]));
      }
    }
    Double_t etamax2 = GetJsonFloat(filename.Data(), "etaMax2Prong");
    printf("Max eta  (2 prong) = %f\n", etamax2);
    if(etamax2>0) fTrackCuts2pr->SetEtaRange(-etamax2, +etamax2);
    Double_t etamax3 = GetJsonFloat(filename.Data(), "etaMax3Prong");
    printf("Max eta  (3 prong) = %f\n", etamax3);
    if(etamax3>0) fTrackCuts3pr->SetEtaRange(-etamax3, +etamax3);

    // vertexer parameters
    printf("--- DCAFitterN parameters ---\n");
    Int_t b_propdca = GetJsonBool(filename.Data(),"propToDCA");
    if(b_propdca==1){
      fVertexerPropagateToPCA = true;
      printf("propdca = %d\n",fVertexerPropagateToPCA);
    }else if(b_propdca==0){
      fVertexerPropagateToPCA = false;
      printf("propdca = %d\n",fVertexerPropagateToPCA);
    }
    Double_t d_maxr = GetJsonFloat(filename.Data(), "maxRad");
    if(d_maxr>0){
      fMaxDecVertRadius2=d_maxr*d_maxr;
      fVertexerMaxR=d_maxr;
      printf("maxr = %f\n",fVertexerMaxR);
    }
    Double_t d_maxdzini = GetJsonFloat(filename.Data(), "maxDZIni");
    if(d_maxdzini>0){
      fVertexerMaxDZIni=d_maxdzini;
      printf("maxdzini = %f\n",fVertexerMaxDZIni);
    }
    Double_t d_minparamchange = GetJsonFloat(filename.Data(), "minParamChange");
    if(d_minparamchange>0){
      fVertexerMinParamChange=d_minparamchange;
      printf("minparamchange = %f\n",fVertexerMinParamChange);
    }
    Double_t d_minrelchi2change = GetJsonFloat(filename.Data(), "minRelChi2Change");
    if(d_minrelchi2change){
      fVertexerMinRelChi2Change=d_minrelchi2change;
      printf("minrelchi2change = %f\n",fVertexerMinRelChi2Change);
    }
    printf("----------------\n");
    Double_t ptMinCand = GetJsonFloat(filename.Data(), "d_pTCandMin");
    printf("Min pt Dzero cand = %f\n", ptMinCand);
    if(ptMinCand>=0.) fMinPtDzero=ptMinCand;
    Double_t ptMaxCand = GetJsonFloat(filename.Data(), "d_pTCandMax");
    printf("Max pt Dzero cand = %f\n", ptMaxCand);
    if(ptMaxCand>=0. && ptMaxCand>=fMinPtDzero) fMaxPtDzero=ptMaxCand;
    Double_t ptMinCand3 = GetJsonFloat(filename.Data(), "cut3ProngPtCandMin");
    printf("Min pt 3-prong cand = %f\n", ptMinCand3);
    if( ptMinCand3>=0.) fMinPt3Prong=ptMinCand3;

    // Selections used in the skimming
    printf("------- CANDIDATE SELECTIONS FOR SKIMMING -------\n");

    Double_t ptTol = GetJsonFloat(filename.Data(), "pTTolerance");
    if(ptTol > 0)
      fPtWithoutVtxToll = ptTol;

    int nptbinlimsDzeroSkims = 0, ncDzeroSkims = 0, nptDzeroSkims = 0;
    float* ptbinlimsDzeroSkims = GetJsonArray(filename.Data(),"pTBinsD0ToPiK",nptbinlimsDzeroSkims);
    float** cutsDzeroSkims = GetJsonMatrix(filename.Data(),"cutsD0ToPiK",nptDzeroSkims,ncDzeroSkims);
    if(nptbinlimsDzeroSkims-1 != nptDzeroSkims)
      AliFatal("Number of pT bins in JSON for Dzero at skims level not consistent, please check it");

    int nptbinlimsJpsiSkims = 0, ncJpsiSkims = 0, nptJpsiSkims = 0;
    float* ptbinlimsJpsiSkims = GetJsonArray(filename.Data(),"pTBinsJpsiToEE",nptbinlimsJpsiSkims);
    float** cutsJpsiSkims = GetJsonMatrix(filename.Data(),"cutsJpsiToEE",nptJpsiSkims,ncJpsiSkims);
    if(nptbinlimsJpsiSkims-1 != nptJpsiSkims)
      AliFatal("Number of pT bins in JSON for J/psi at skims level not consistent, please check it");

    int nptbinlimsDplusSkims = 0, ncDplusSkims = 0, nptDplusSkims = 0;
    float* ptbinlimsDplusSkims = GetJsonArray(filename.Data(),"pTBinsDPlusToPiKPi",nptbinlimsDplusSkims);
    float** cutsDplusSkims = GetJsonMatrix(filename.Data(),"cutsDPlusToPiKPi",nptDplusSkims,ncDplusSkims);
    if(nptbinlimsDplusSkims-1 != nptDplusSkims)
      AliFatal("Number of pT bins in JSON for Dplus at skims level not consistent, please check it");

    int nptbinlimsDsSkims = 0, ncDsSkims = 0, nptDsSkims = 0;
    float* ptbinlimsDsSkims = GetJsonArray(filename.Data(),"pTBinsDsToPiKK",nptbinlimsDsSkims);
    float** cutsDsSkims = GetJsonMatrix(filename.Data(),"cutsDsToPiKK",nptDsSkims,ncDsSkims);
    if(nptbinlimsDsSkims-1 != nptDsSkims)
      AliFatal("Number of pT bins in JSON for Ds at skims level not consistent, please check it");

    int nptbinlimsLcSkims = 0, ncLcSkims = 0, nptLcSkims = 0;
    float* ptbinlimsLcSkims = GetJsonArray(filename.Data(),"pTBinsLcToPKPi",nptbinlimsLcSkims);
    float** cutsLcSkims = GetJsonMatrix(filename.Data(),"cutsLcToPKPi",nptLcSkims,ncLcSkims);
    if(nptbinlimsLcSkims-1 != nptLcSkims)
      AliFatal("Number of pT bins in JSON for Lc at skims level not consistent, please check it");

    int nptbinlimsXicSkims = 0, ncXicSkims = 0, nptXicSkims = 0;
    float* ptbinlimsXicSkims = GetJsonArray(filename.Data(),"pTBinsXicToPKPi",nptbinlimsXicSkims);
    float** cutsXicSkims = GetJsonMatrix(filename.Data(),"cutsXicToPKPi",nptXicSkims,ncXicSkims);
    if(nptbinlimsXicSkims-1 != nptXicSkims)
      AliFatal("Number of pT bins in JSON for Xic at skims level not consistent, please check it");

    for (Int_t ib = 0; ib < nptbinlimsDzeroSkims; ib++) {
      fPtBinLimsDzeroSkims[ib] = ptbinlimsDzeroSkims[ib];
    }
    for (Int_t ib = 0; ib < nptDzeroSkims; ib++) {
      for (Int_t jc = 0; jc < ncDzeroSkims; jc++) {
        fDzeroSkimCuts[ib][jc] = cutsDzeroSkims[ib][jc];
      }
      AliWarning(Form("Dzero cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; d0xd0  < %f\n", fPtBinLimsDzeroSkims[ib], fPtBinLimsDzeroSkims[ib+1], fDzeroSkimCuts[ib][0], fDzeroSkimCuts[ib][1], fDzeroSkimCuts[ib][2], fDzeroSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsJpsiSkims; ib++) {
      fPtBinLimsJpsiSkims[ib] = ptbinlimsJpsiSkims[ib];
    }
    for (Int_t ib = 0; ib < nptJpsiSkims; ib++) {
      for (Int_t jc = 0; jc < ncJpsiSkims; jc++) {
        fJpsiSkimCuts[ib][jc] = cutsJpsiSkims[ib][jc];
      }
      AliWarning(Form("J/psi cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; d0xd0  < %f\n", fPtBinLimsJpsiSkims[ib], fPtBinLimsJpsiSkims[ib+1], fJpsiSkimCuts[ib][0], fJpsiSkimCuts[ib][1], fJpsiSkimCuts[ib][2], fJpsiSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsDplusSkims; ib++) {
      fPtBinLimsDplusSkims[ib] = ptbinlimsDplusSkims[ib];
    }
    for (Int_t ib = 0; ib < nptDplusSkims; ib++) {
      for (Int_t jc = 0; jc < ncDplusSkims; jc++) {
        fDplusSkimCuts[ib][jc] = cutsDplusSkims[ib][jc];
      }
      AliWarning(Form("Dplus cuts: %f < pt < %f  ; %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsDplusSkims[ib], fPtBinLimsDplusSkims[ib+1], fDplusSkimCuts[ib][0], fDplusSkimCuts[ib][1], fDplusSkimCuts[ib][2], fDplusSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsDsSkims; ib++) {
      fPtBinLimsDsSkims[ib] = ptbinlimsDsSkims[ib];
    }
    for (Int_t ib = 0; ib < nptDsSkims; ib++) {
      for (Int_t jc = 0; jc < ncDsSkims; jc++) {
        fDsSkimCuts[ib][jc] = cutsDsSkims[ib][jc];
      }
      AliWarning(Form("Ds cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsDsSkims[ib], fPtBinLimsDsSkims[ib+1], fDsSkimCuts[ib][0], fDsSkimCuts[ib][1], fDsSkimCuts[ib][2], fDsSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsLcSkims; ib++) {
      fPtBinLimsLcSkims[ib] = ptbinlimsLcSkims[ib];
    }
    for (Int_t ib = 0; ib < nptLcSkims; ib++) {
      for (Int_t jc = 0; jc < ncLcSkims; jc++) {
        fLcSkimCuts[ib][jc] = cutsLcSkims[ib][jc];
      }
      AliWarning(Form("Lc cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsLcSkims[ib], fPtBinLimsLcSkims[ib+1], fLcSkimCuts[ib][0], fLcSkimCuts[ib][1], fLcSkimCuts[ib][2], fLcSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsXicSkims; ib++) {
      fPtBinLimsXicSkims[ib] = ptbinlimsXicSkims[ib];
    }
    for (Int_t ib = 0; ib < nptXicSkims; ib++) {
      for (Int_t jc = 0; jc < ncXicSkims; jc++) {
        fXicSkimCuts[ib][jc] = cutsXicSkims[ib][jc];
      }
      AliWarning(Form("Xic cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsXicSkims[ib], fPtBinLimsXicSkims[ib+1], fXicSkimCuts[ib][0], fXicSkimCuts[ib][1], fXicSkimCuts[ib][2], fXicSkimCuts[ib][3]));
    }

    
    Double_t cutcpaV0 = GetJsonFloat(filename.Data(), "cosPAV0");
    if(cutcpaV0>0){
      fMinCosPointV0=cutcpaV0;
      printf("cosPAV0 cut = %f\n",fMinCosPointV0);
    }

    Double_t cutinvmV0 = GetJsonFloat(filename.Data(), "cutInvMassV0");
    if(cutinvmV0>0){
      fCutOnK0sMass=cutinvmV0;
      printf("invmass V0 cut = %f\n",fCutOnK0sMass);
    }
    
    
    printf("---------------------------------------------\n");

    printf("---- TEST READOUT OF ARRAYS ----\n");
    int nptbinlims;
    float* ptbins = GetJsonArray(filename.Data(),"pTBins",nptbinlims);
    printf("%d ptbins. Limits: ",nptbinlims-1);
    for(int j=0; j<nptbinlims; j++) printf("%.1f ",ptbins[j]);
    printf("\n");
    int npt,nc;
    float** cuts = GetJsonMatrix(filename.Data(),"D0_to_pi_K_cuts",npt,nc);
    printf("D0 2D cuts: %d pt bins, %d cut variables:\n",npt,nc);
    for(int j=0; j<npt; j++){
      for(int k=0; k<nc; k++){
        printf("%.2f ",cuts[j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }else{
    AliError(Form("Json configuration file %s not found\n",filename.Data()));
  }
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::MakeJetFinding(AliESDEvent *esd) {
	Double_t d03[3] = {0., 0., 0.};
	AliAODRecoDecay* rd4massCalc3 = new AliAODRecoDecay(0x0, 3, 1, d03);
	Double_t covMatrix[6];
  
	TString trClass = esd->GetFiredTriggerClasses(); 
	if (triggerstring != "" && !trClass.Contains(triggerstring)) return;
	printf("      Fired Trigger Classes %s\n", trClass.Data());

	Int_t maxTracksToProcess = 9999999; /// temporary to limit the time duration of tests
	Int_t totTracks = TMath::Min(maxTracksToProcess, esd->GetNumberOfTracks());
	AliESDVertex* primvtx = (AliESDVertex*)esd->GetPrimaryVertex();
	TString title = primvtx->GetTitle();
	if (primvtx->IsFromVertexer3D() || primvtx->IsFromVertexerZ()) return;
	if (primvtx->GetNContributors() < 2) return;
	hvertexx->Fill(primvtx->GetX());
	hvertexy->Fill(primvtx->GetY());
	hvertexz->Fill(primvtx->GetZ());
	AliAODVertex *vertexAODp = ConvertToAODVertex(primvtx);
	if (triggerstring != "" && !trClass.Contains(triggerstring)) return;
	Double_t fBzkG = (Double_t)esd->GetMagneticField();

	// Apply single track cuts and flag them
	UChar_t* status = new UChar_t[totTracks];
	for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
		status[iTrack] = 0;
		AliESDtrack* track = esd->GetTrack(iTrack);
		hpt_nocuts->Fill(track->Pt());
		htgl_nocuts->Fill(track->GetTgl());
		if (SingleTrkCutsSimple(track, minncluTPC, ptmintrack, dcatoprimxymin, primvtx, fBzkG))
			status[iTrack] = 1; //FIXME
	}

	TObjArray* twoTrackArray = new TObjArray(2);
	TObjArray* threeTrackArray = new TObjArray(3);
	TObjArray *twoTrackArrayCasc = new TObjArray(2);//
	AliESDVertex* primVtxTrk = (AliESDVertex*)esd->GetPrimaryVertex();//

	AliVertexerTracks* vt = new AliVertexerTracks(fBzkG);

	Double_t mom0[3], mom1[3], mom2[3];
	for (Int_t iPosTrack_0 = 0; iPosTrack_0 < totTracks; iPosTrack_0++) {
		AliESDtrack* track_p0 = esd->GetTrack(iPosTrack_0);
		track_p0->GetPxPyPz(mom0);
		if (status[iPosTrack_0] == 0) continue;
		AliExternalTrackParam* trackext = (AliExternalTrackParam*)track_p0;
		double b[2];
		double bCov[3];
		trackext->PropagateToDCA(primvtx, fBzkG, 100., b, bCov);
		hpt_cuts->Fill(track_p0->Pt());
		hdcatoprimxy_cuts->Fill(b[0]);
		htgl_cuts->Fill(track_p0->GetTgl());
		hitsmap->Fill(track_p0->GetITSClusterMap());
		if (track_p0->Charge() < 0) continue;

		for (Int_t iNegTrack_0 = 0; iNegTrack_0 < totTracks; iNegTrack_0++) {
		AliESDtrack* track_n0 = esd->GetTrack(iNegTrack_0);
		track_n0->GetPxPyPz(mom1);
		if (track_n0->Charge() > 0) continue;
		if (status[iNegTrack_0] == 0) continue;

		twoTrackArray->AddAt(track_p0, 0);
		twoTrackArray->AddAt(track_n0, 1);
		AliESDVertex* trkv = ReconstructSecondaryVertex(twoTrackArrayCasc, primVtxTrk);
		if (trkv == 0x0) {
			twoTrackArray->Clear();
			continue;
		}

		hvx->Fill(trkv->GetX());
		hvy->Fill(trkv->GetY());
		hvz->Fill(trkv->GetZ());
		double deltax = trkv->GetX() - primvtx->GetX();
		double deltay = trkv->GetY() - primvtx->GetY();
		double deltaz = trkv->GetZ() - primvtx->GetZ();
		double decaylength = TMath::Sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
		double decaylengthxy = TMath::Sqrt(deltax * deltax + deltay * deltay);

		AliAODVertex* vertexAOD = ConvertToAODVertex(trkv);
		delete trkv;
		AliAODRecoDecayHF2Prong* the2Prong = Make2Prong(twoTrackArray, vertexAOD, fBzkG);
		the2Prong->SetOwnPrimaryVtx(vertexAODp);

		Int_t twoProngSelection = 3;
		if (selectD0 + selectD0bar > 0) twoProngSelection = TwoProngSelectionCuts(the2Prong, candpTMin, candpTMax);
		Double_t m0 = the2Prong->InvMassD0();
		Double_t m0b = the2Prong->InvMassD0bar();
		if (twoProngSelection > 0) {
		if (selectD0 == 0 || twoProngSelection == 1 || twoProngSelection == 3) hmass0->Fill(m0);
        if (selectD0bar == 0 || twoProngSelection == 2 || twoProngSelection == 3) hmass0->Fill(m0b);
        hdecayxyz->Fill(decaylength);
        hdecayxy->Fill(decaylengthxy);
        hptD0->Fill(the2Prong->Pt());
        hptprong0->Fill(the2Prong->PtProng(0));
        hptprong1->Fill(the2Prong->PtProng(1));
        hd0->Fill(the2Prong->Getd0Prong(0));
        hd0->Fill(the2Prong->Getd0Prong(1));
        hd0d0->Fill(the2Prong->Prodd0d0());
        hImpParErr->Fill(the2Prong->Getd0errProng(0));
        hImpParErr->Fill(the2Prong->Getd0errProng(1));
        hDecLenErr->Fill(the2Prong->DecayLengthError());
        hDecLenXYErr->Fill(the2Prong->DecayLengthXYError());
        the2Prong->GetPrimaryVtx()->GetCovMatrix(covMatrix);
        hCovPVXX->Fill(covMatrix[0]);
        the2Prong->GetSecondaryVtx()->GetCovMatrix(covMatrix);
        hCovSVXX->Fill(covMatrix[0]);
		
		AliFJWrapper *fFastJetWrapper;
		fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
		fFastJetWrapper->Clear();
		fFastJetWrapper->SetR(0.4); 
		fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
		fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::E_scheme);
		fFastJetWrapper->SetStrategy(fastjet::Strategy::Best);
		fFastJetWrapper->SetGhostArea(0.005); 
		fFastJetWrapper->SetAreaType(fastjet::AreaType::passive_area);

		bool isHFJet=false;
		fFastJetWrapper->Clear();
		for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
			AliESDtrack* track = esd->GetTrack(iTrack);
			if (track->Pt() >= 0.15 && TMath::Abs(track->Eta()) < 0.9){ 
				if (iTrack==iNegTrack_0 || iTrack==iPosTrack_0) continue;
				fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack+2);
			}
		}
		fFastJetWrapper->AddInputVector(the2Prong->Px(), the2Prong->Py(), the2Prong->Pz(), the2Prong->ED0(),1);

		fFastJetWrapper->Run();
		std::vector<fastjet::PseudoJet> jets = fFastJetWrapper->GetInclusiveJets();
		for (Int_t ijet=0; ijet<jets.size(); ijet++) {
			isHFJet=false;
			fastjet::PseudoJet jet = jets[ijet];
			if (jet.pt() < 0.15 || jet.perp() >= 1000.0 || TMath::Abs(jet.eta()) >= 0.5) continue;
			std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));
			for (Int_t iconstituent=0; iconstituent<constituents.size(); iconstituent++) {
				if (constituents[iconstituent].user_index()==1) isHFJet=true;
				break;
			}
			if(!isHFJet) continue;
			hjetpt->Fill(jet.pt());
			hjetE->Fill(jet.E());
			hjetpx->Fill(jet.px());
			hjetpy->Fill(jet.py());
			hjetpz->Fill(jet.pz());
			hjetphi->Fill(jet.phi());
			hjetrap->Fill(jet.rap());
			hjetconstituents->Fill(constituents.size());
		
			fastjet::JetDefinition subJetDef(fastjet::JetAlgorithm::cambridge_algorithm , 0.4*2.5,fastjet::RecombinationScheme::E_scheme, fastjet::Best);
			try{
				fastjet::ClusterSequence reclusterSeq(constituents, subJetDef);
				std::vector<fastjet::PseudoJet> reclusteredJet =  reclusterSeq.inclusive_jets(0.0);
				reclusteredJet = sorted_by_pt(reclusteredJet);
         
				fastjet::PseudoJet daughterSubJet = reclusteredJet[0];
				fastjet::PseudoJet parentSubJet1; 
				fastjet::PseudoJet parentSubJet2;

				Float_t zg=-1.0,rg=-1.0;
				Int_t nsd=0;
				bool softDropped=false;
				bool isHFSubJet=false;
				std::vector<fastjet::PseudoJet> constituentsSubJet;
				while(daughterSubJet.has_parents(parentSubJet1,parentSubJet2)){
					isHFSubJet=false;
					constituentsSubJet=parentSubJet1.constituents();
					for (Int_t iconstituent=0; iconstituent<constituentsSubJet.size(); iconstituent++){
						if (constituentsSubJet[iconstituent].user_index()==1) isHFSubJet=true;
					}
					if (!isHFSubJet) std::swap(parentSubJet1,parentSubJet2);
					zg=parentSubJet2.perp()/(parentSubJet1.perp()+parentSubJet2.perp());
					rg=parentSubJet1.delta_R(parentSubJet2);
	
					if (zg >= 0.1*TMath::Power(rg/0.4,0.0)){
						if(!softDropped){
							hjetzg->Fill(zg);
							hjetrg->Fill(rg);
							softDropped=true;
						}
						nsd++;
					}
					daughterSubJet=parentSubJet1;
				}
				hjetnsd->Fill(nsd);
			}catch (fastjet::Error) {}
		}
		delete fFastJetWrapper;
		}
		delete the2Prong;
		delete vertexAOD;
   
		if (do3Prongs) {
			for (Int_t iPosTrack_1 = iPosTrack_0 + 1; iPosTrack_1 < totTracks; iPosTrack_1++) {
				AliESDtrack* track_p1 = esd->GetTrack(iPosTrack_1);
				if (!track_p1) continue;
				if (track_p1->Charge() < 0) continue;
				track_p1->GetPxPyPz(mom2);
				if (status[iPosTrack_1] == 0) continue;
				// order tracks according to charge: +-+
				threeTrackArray->AddAt(track_p0, 0);
				threeTrackArray->AddAt(track_n0, 1);
				threeTrackArray->AddAt(track_p1, 2);
				Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
				if (massSel == 0) {
					threeTrackArray->Clear();
					continue;
				}
				//AliESDVertex* trkv3 = ReconstructSecondaryVertex(vt, threeTrackArray, primvtx, d_maxr);
				AliESDVertex* trkv3 = ReconstructSecondaryVertex(twoTrackArrayCasc, primVtxTrk);
				if (trkv3 == 0x0) {
					threeTrackArray->Clear();
					continue;
				}
				AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
				AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, fBzkG);
				//  the3Prong->SetOwnPrimaryVtx(vertexAODp);
				if (massSel & (1 << kbitDplus)) {
					Double_t mp = the3Prong->InvMassDplus();
					hmassP->Fill(mp);
					hvx3->Fill(trkv3->GetX());
					hvy3->Fill(trkv3->GetY());
					hvz3->Fill(trkv3->GetZ());
				}
				delete trkv3;
				delete the3Prong;
				delete vertexAOD3;
				threeTrackArray->Clear();
			}
			
			for (Int_t iNegTrack_1 = iNegTrack_0 + 1; iNegTrack_1 < totTracks; iNegTrack_1++) {
				AliESDtrack* track_n1 = esd->GetTrack(iNegTrack_1);
				if (!track_n1) continue;
				if (track_n1->Charge() > 0) continue;
				track_n1->GetPxPyPz(mom2);
				if (status[iNegTrack_1] == 0) continue;
				// order tracks according to charge: -+-
				threeTrackArray->AddAt(track_n0, 0);
				threeTrackArray->AddAt(track_p0, 1);
				threeTrackArray->AddAt(track_n1, 2);
				Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
				if (massSel == 0) {
					threeTrackArray->Clear();
					continue;
				}
				AliESDVertex* trkv3 = ReconstructSecondaryVertex(twoTrackArrayCasc, primVtxTrk);
				if (trkv3 == 0x0) {
					threeTrackArray->Clear();
					continue;
				}
				AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
				AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, fBzkG);
				//  the3Prong->SetOwnPrimaryVtx(vertexAODp);
				if (massSel & (1 << kbitDplus)) {
					Double_t mp = the3Prong->InvMassDplus();
					hmassP->Fill(mp);
				}
				delete trkv3;
				delete the3Prong;
				delete vertexAOD3;
				threeTrackArray->Clear();
			}
		}
        twoTrackArray->Clear();
		}
	}
    delete[] status;
    delete vt;
    delete twoTrackArray;
    delete threeTrackArray;
    delete rd4massCalc3; 
}

//_______________________________________________________________________________________
/*Int_t AliAnalysisTaskHFFindJets::GetpTBin(Double_t candpT){
	Double_t pTBins[npTBins + 1] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 16.0, 20.0, 24.0, 36.0, 50.0, 100.0};
	for (Int_t i = 0; i < npTBins; i++) {
		if (candpT >= pTBins[i] && candpT < pTBins[i + 1]) {
			return i;
		}
	}
	return -1;
}*/


//______________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::GetPtBin(Double_t ptCand, Double_t* ptBinLims, Double_t nPtBins)
{
  for (Int_t i = 0; i < nPtBins; i++) {
    if (ptCand>=ptBinLims[i] && ptCand<ptBinLims[i+1]){
      return i;
    }
  }
  return -1;
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::UserExec(Option_t *) {
	
	AliESDEvent *esd = (AliESDEvent*) (InputEvent()); //looping over events
	if(!esd) {
		printf("AliAnalysisTaskHFFindJets:::UserExec(): bad ESD\n");
		return;
	}
	
	Double_t centr=-1;
	AliMultSelection* mulSel=0x0;
	if(fSelectOnCentrality){
		mulSel = (AliMultSelection*)esd->FindListObject("MultSelection");
		if(mulSel){
		centr=mulSel->GetMultiplicityPercentile(fCentrEstimator.Data());
		}
	}
  
	if(fUsePhysSel){
		Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
		if(!isPhysSel) return;
	}
 
	if(doJetFinding){
		MakeJetFinding(esd);
	}
  
	PostData(1,fOutput);
}

//______________________________________________________________________________
std::string AliAnalysisTaskHFFindJets::GetJsonString(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  char* value=0x0;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      value=strtok(line, ":");
      value=strtok(NULL, ":");
      break;
    }
  }
  fclose(fj);
  TString sValue = value;
  sValue.ReplaceAll("\"", "");
  sValue.ReplaceAll("\n", "");
  sValue.ReplaceAll("\t", "");
  sValue.ReplaceAll(" ", "");
  return std::string(sValue.Data());
}
//______________________________________________________________________________
int AliAnalysisTaskHFFindJets::GetJsonBool(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  int value=-1;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      if(temp.Contains("true")) value=1;
      if(temp.Contains("false")) value=0;
      break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
int AliAnalysisTaskHFFindJets::GetJsonInteger(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  int value=-999;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      value=temp.Atoi();
      break;
    }
  }
  fclose(fj);
  return value;
}
//______________________________________________________________________________
float AliAnalysisTaskHFFindJets::GetJsonFloat(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float value=-999.;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      value=temp.Atof();
      break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
float* AliAnalysisTaskHFFindJets::GetJsonArray(const char* jsonFileName, const char* key, int &size){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float* arrVals=0x0;
  while(!feof(fj)){
    fgets(line,500,fj);
    if(strstr(line,key)){
      TString full="";
      while(!feof(fj)){
        fgets(line,500,fj);
        int len = strlen(line);
        if(line[len-1]=='\n') line[len-1]=0;
        full.Append(line);
        if(strstr(line,"}")) break;
      }
      full.ReplaceAll("\"values\":","");
      full.ReplaceAll(" ","");
      full.ReplaceAll("},","");
      full.ReplaceAll("}","");
      TObjArray* arrStr=full.Tokenize(",");
      size=arrStr->GetEntriesFast();
      arrVals=new float[size];
      for(int j=0; j<size; j++){
        TObjString* sss=(TObjString*)arrStr->At(j);
        TString strval=sss->GetString();
        strval.ReplaceAll("[","");
        strval.ReplaceAll("]","");
        strval.ReplaceAll("\"","");
        arrVals[j]=strval.Atof();
      }
      arrStr->Delete();
      delete arrStr;
    }
  }
  return arrVals;
}

//______________________________________________________________________________
float** AliAnalysisTaskHFFindJets::GetJsonMatrix(const char* jsonFileName, const char* key, int &size1, int &size2){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float** arrVals=0x0;
  while(!feof(fj)){
    fgets(line,500,fj);
    if(strstr(line,key)){
      TString full="";
      while(!feof(fj)){
        fgets(line,500,fj);
        int len = strlen(line);
        if(line[len-1]=='\n') line[len-1]=0;
        full.Append(line);
        if(strstr(line,"}")) break;
      }
      full.ReplaceAll("\"values\":","");
      full.ReplaceAll(" ","");
      full.ReplaceAll("},","");
      full.ReplaceAll("}","");
      TObjArray* rowArrStr=full.Tokenize("]");
      size1=rowArrStr->GetEntriesFast();
      arrVals=new float*[size1];
      for(int j=0; j<size1; j++){
        TObjString* rowStr=(TObjString*)rowArrStr->At(j);
        TString rowStrVal=rowStr->GetString();
        TObjArray* arrStr=rowStrVal.Tokenize(",");
        size2=arrStr->GetEntriesFast();
        arrVals[j]=new float[size2];
        for(int k=0; k<size2; k++){
          TObjString* sss=(TObjString*)arrStr->At(k);
          TString strval=sss->GetString();
          strval.ReplaceAll(",[","");
          strval.ReplaceAll("[","");
          strval.ReplaceAll("]","");
          strval.ReplaceAll("\"","");
          arrVals[j][k]=strval.Atof();
        }
      }
    }
  }
  return arrVals;
}

//start
Int_t AliAnalysisTaskHFFindJets::MatchToMC(AliAODRecoDecay* rd, Int_t pdgabs, AliMCEvent* mcEvent,
                                                 Int_t ndgCk, const TObjArray *trkArray, const Int_t *pdgDg) const {return 0;}
                                                
Int_t AliAnalysisTaskHFFindJets::LcSelectionCuts(
    AliAODRecoDecayHF3Prong *cand) {return 0;}
    
Int_t AliAnalysisTaskHFFindJets::JpsiSelectionCuts(AliAODRecoDecayHF2Prong* cand,AliESDtrack* trk_p,AliESDtrack* trk_n,AliESDVertex* primvtx,float bzkG)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::DplusSelectionCuts(AliAODRecoDecayHF3Prong* cand, Double_t bzkG)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::DzeroSelectionCuts(AliAODRecoDecayHF2Prong* cand)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::LcSkimCuts(AliAODRecoDecayHF3Prong* cand)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::DsSkimCuts(AliAODRecoDecayHF3Prong* cand)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::DplusSkimCuts(AliAODRecoDecayHF3Prong* cand)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::JpsiSkimCuts(AliAODRecoDecayHF2Prong* cand)
{return 0;}

Int_t AliAnalysisTaskHFFindJets::DzeroSkimCuts(AliAODRecoDecayHF2Prong* cand)
{return 0;}

Bool_t AliAnalysisTaskHFFindJets::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{return 0;}

AliAODRecoCascadeHF* AliAnalysisTaskHFFindJets::MakeCascade(TObjArray *twoTrackArray, AliAODVertex* secVert, Double_t bzkG){return 0;}


Int_t AliAnalysisTaskHFFindJets::SelectInvMassAndPt2prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc2)
{return 0;}

Bool_t AliAnalysisTaskHFFindJets::SelectV0(AliESDv0 *v0, AliESDVertex* primvtx)
{return 0;}

//Int_t AliAnalysisTaskHFSimpleVertices::SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG, Double_t d0track[2]){}

Bool_t AliAnalysisTaskHFFindJets::GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG)
{
  /// fast calculation (no covariance matrix treatment) of track momentum at secondary vertex

  Double_t alpha = tr->GetAlpha();
  Double_t sn = TMath::Sin(alpha), cs = TMath::Cos(alpha);
  Double_t x = tr->GetX(), y = tr->GetParameter()[0], snp = tr->GetParameter()[2];
  Double_t xv = secVert->GetX() * cs + secVert->GetY() * sn;
  Double_t yv = -secVert->GetX() * sn + secVert->GetY() * cs;
  x -= xv;
  y -= yv;
  Double_t crv = tr->GetC(bzkG);
  if (TMath::Abs(bzkG) < 0.000001)
    crv = 0.;
  double csp = TMath::Sqrt((1. - snp) * (1. + snp));

  Double_t tgfv = -(crv * x - snp) / (crv * y + csp);
  cs = 1. / TMath::Sqrt(1 + tgfv * tgfv);
  sn = cs < 1. ? tgfv * cs : 0.;

  x = xv * cs + yv * sn;
  Double_t alpNew = alpha + TMath::ASin(sn);
  Double_t ca = TMath::Cos(alpNew - alpha), sa = TMath::Sin(alpNew - alpha);
  Double_t p2 = tr->GetSnp();
  Double_t xNew = tr->GetX() * ca + tr->GetY() * sa;
  Double_t p2New = p2 * ca - TMath::Sqrt((1. - p2) * (1. + p2)) * sa;
  momentum[0] = tr->GetSigned1Pt();
  momentum[1] = p2New * (x - xNew) * tr->GetC(bzkG);
  momentum[2] = tr->GetTgl();
  Bool_t retCode = tr->Local2GlobalMomentum(momentum, alpNew);
  return retCode;
}

void AliAnalysisTaskHFFindJets::ProcessTriplet(TObjArray* threeTrackArray, AliAODRecoDecay* rd4massCalc3, AliESDVertex* primVtxTrk, AliAODVertex *vertexAODp, float bzkG, double dist12, AliMCEvent* mcEvent){
}

Int_t AliAnalysisTaskHFFindJets::GetPtBinSingleTrack(Double_t ptTrack)
{ return 0;}







//end

//______________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG, Double_t d0track[2])
{
  Double_t covd0track[3];
  if (!trk->PropagateToDCA(primVert, bzkG, kVeryBig, d0track, covd0track)) return kFALSE;
  trk->RelateToVertex(primVert, bzkG, kVeryBig);
  Int_t retCode=0;

  // test first impact parameter
  Int_t iPtTrack = GetPtBinSingleTrack(trk->Pt());
  if(TMath::Abs(d0track[0]) >= fSingleTrackCuts2Prong[iPtTrack][0] && TMath::Abs(d0track[0]) <= fSingleTrackCuts2Prong[iPtTrack][1])
    retCode+=1;
  if(TMath::Abs(d0track[0]) >= fSingleTrackCuts3Prong[iPtTrack][0] && TMath::Abs(d0track[0]) <= fSingleTrackCuts3Prong[iPtTrack][1])
    retCode+=2;

  // test other cuts if impact-parameter tested
  Int_t retCodeCurrent=retCode;
  if((retCodeCurrent == 1 || retCodeCurrent == 3) && !fTrackCuts2pr->AcceptTrack(trk)) retCode-=1;
  if(retCodeCurrent >= 2 && !fTrackCuts3pr->AcceptTrack(trk)) retCode-=2;

  //test cuts for bachelor track
  if(fTrackCutsBach->AcceptTrack(trk)) retCode+=4;
  
  return retCode;
}

//______________________________________________________________________________
/*AliESDVertex* AliAnalysisTaskHFFindJets::ReconstructSecondaryVertex(TObjArray* trkArray, AliESDVertex* primvtx)
{

  AliESDVertex* trkv =0x0;
  // printf("------\n");
  // for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
  //   AliExternalTrackParam *tp=(AliExternalTrackParam *)trkArray->At(jt);
  //   printf("%f ",tp->Pt());
  // }
  // printf("\n");
  if(fSecVertexerAlgo==0){
    fVertexerTracks->SetVtxStart(primvtx);
    trkv = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(trkArray);
    if (trkv->GetNContributors() != trkArray->GetEntriesFast()) return 0x0;
  }else if(fSecVertexerAlgo==1){
    o2::track::TrackParCov* o2Track[3] = {nullptr};
    for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
      o2Track[jt] = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)trkArray->At(jt));
    }
    Int_t nVert=0;
    Double_t vertCoord[3];
    Double_t vertCov[6];
    Double_t vertChi2=-999.;
    if(trkArray->GetEntriesFast()==2){
      nVert=fO2Vertexer2Prong.process(*o2Track[0], *o2Track[1]);
      if(nVert){
        fO2Vertexer2Prong.propagateTracksToVertex();
        auto vertPos = fO2Vertexer2Prong.getPCACandidate();
        auto vertCMat = fO2Vertexer2Prong.calcPCACovMatrix().Array();
        for(Int_t ic=0; ic<3; ic++) vertCoord[ic]=vertPos[ic];
        for(Int_t ic=0; ic<6; ic++) vertCov[ic]=vertCMat[ic];
        vertChi2 = fO2Vertexer2Prong.getChi2AtPCACandidate();
      }
    }else if(trkArray->GetEntriesFast()==3){
      nVert=fO2Vertexer3Prong.process(*o2Track[0], *o2Track[1], *o2Track[2]);
      if(nVert){
        fO2Vertexer3Prong.propagateTracksToVertex();
        auto vertPos = fO2Vertexer3Prong.getPCACandidate();
        auto vertCMat = fO2Vertexer3Prong.calcPCACovMatrix().Array();
        for(Int_t ic=0; ic<3; ic++) vertCoord[ic]=vertPos[ic];
        for(Int_t ic=0; ic<6; ic++) vertCov[ic]=vertCMat[ic];
        vertChi2 = fO2Vertexer3Prong.getChi2AtPCACandidate();
      }
    }
    if(nVert) trkv = new AliESDVertex(vertCoord,vertCov,vertChi2,trkArray->GetEntriesFast());
    //   else printf("nVert = %d\n",nVert);
  }else if(fSecVertexerAlgo==2){
    KFParticle  dMesonVert;
    double posmom[6],cov[21];
    for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
      AliExternalTrackParam* trpar=(AliExternalTrackParam*)trkArray->At(jt);
      trpar->GetXYZ(posmom);
      trpar->GetPxPyPz(posmom+3);
      trpar->GetCovarianceXYZPxPyPz(cov);
      KFParticle trKFpar;
      trKFpar.Create(posmom,cov,trpar->Charge(),0.13957); // mass of the pion for the time being
      dMesonVert.AddDaughter(trKFpar);
    }
    Double_t vertChi2 = dMesonVert.GetChi2() / dMesonVert.GetNDF();
    Double_t vertCoord[3]={dMesonVert.X(),dMesonVert.Y(),dMesonVert.Z()};
    Double_t vertCov[6];
    vertCov[0]=dMesonVert.GetCovariance(0,0);
    vertCov[1]=dMesonVert.GetCovariance(0,1);
    vertCov[2]=dMesonVert.GetCovariance(1,1);
    vertCov[3]=dMesonVert.GetCovariance(0,2);
    vertCov[4]=dMesonVert.GetCovariance(1,2);
    vertCov[5]=dMesonVert.GetCovariance(2,2);
    trkv = new AliESDVertex(vertCoord,vertCov,vertChi2,trkArray->GetEntriesFast());
  }
  if(!trkv) return 0x0;
  Double_t vertRadius2 = trkv->GetX() * trkv->GetX() + trkv->GetY() * trkv->GetY();
  if(vertRadius2>fMaxDecVertRadius2) return 0x0;
  //  trkv->Print("all");
  //  printf("=============\n");
  return trkv;
}*/

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::Terminate(Option_t */*option*/) 
{
	// Terminate analysis
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutput) {
		printf("ERROR: fOutput not available\n");
		return;
	}
	//fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
	//printf("AliAnalysisTaskHFSimpleVertices::Terminate --- Number of events: read = %.0f\n",fHistNEvents->GetBinContent(1));
	return;
}


