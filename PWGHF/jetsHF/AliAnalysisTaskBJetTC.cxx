#include "TRandom.h"
#include "TList.h"
#include "AliParticleContainer.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "TLorentzVector.h"
#include <TFile.h>
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include <THnSparse.h>
#include <TH3D.h>

#include "AliConvEventCuts.h"
#include "AliV0ReaderV1.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"

#include "AliAODMCHeader.h"
#include "AliJetContainer.h"
#include "AliPicoTrack.h"
#include "TMath.h"
#include "AliAnalysisTaskBJetTC.h"
#include "AliExternalTrackParam.h"
#include "AliVertexerTracks.h"
#include "AliHFJetsTagging.h"
#include "AliPIDResponse.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliRDHFJetsCutsVertex.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliGenPythiaEventHeader.h"
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;
using std::vector;
ClassImp(AliAnalysisTaskBJetTC)

const Double_t AliAnalysisTaskBJetTC::fgkMassPion   = 0.13957;
const Double_t AliAnalysisTaskBJetTC::fgkMassKshort = 0.497614;
const Double_t AliAnalysisTaskBJetTC::fgkMassProton = 0.938272;
const Double_t AliAnalysisTaskBJetTC::fgkMassLambda = 1.11568;

// ######################################################################################## CONSTRUCTORS
AliAnalysisTaskBJetTC::AliAnalysisTaskBJetTC(): AliAnalysisTaskEmcalJet("AliAnalysisTaskBJetTC", kTRUE),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),   fReaderGammas(NULL),
fHFJetUtils(0x0),
fRespoPID(0x0),
fPtHardThreshold(0.0),
fPythiaEventWeight(1.0),
fDoImprovedDCACut(kTRUE),
fVertexConstraint(kFALSE),
fThresholdIP(0.008),
fDoDeltaPtWithSignal(kFALSE),
fDiamond(0x0),
fVertexer(0x0),
fDoJetProbabilityAnalysis(kFALSE),
fDoCharmFractions(kFALSE),
fUsePartonDef(kTRUE),
fUseIPs(kFALSE),
fDoJetMass(kFALSE),
fDoSVEnergyFraction(kFALSE),
fDoPtRelAnalysis(0),
fDoSelectionPtRel(0),
//Bjet Cuts
fTCMinTrackPt(0.5),
fTCMinClusTPC(80),
fTCMinHitsITS(2),
fTCMaxChi2pNDF(5.),
fTCMaxIPxy(1.),
fTCMaxIPz(5.),
fTCMaxDecayLength(5),
fTCMaxDCATrackJet(0.07),
fMaxFactorPtHardJet(10.0),
fhistInclusiveJetCuts(0x0),
fhistbJetCuts(0x0),
fhistcJetCuts(0x0),
fhistlfJetCuts(0x0),
//____
fh1dEventRejectionRDHFCuts(0x0),
fh1dVertexZ(0x0),
fh1dVertexZAccepted(0x0),
fh1dVertexR(0x0),
fh1dVertexRAccepted(0x0),
fh1dTracksAccepeted(0x0),
fh1dTracksImpParXY(0x0),
fh1dTracksImpParXYZ(0x0),
fh1dTracksImpParXYSignificance(0x0),
fh1dTracksImpParXYZSignificance(0x0),
fh1dTracksImpParXYTruth(0x0),
fh1dTracksImpParXYZTruth(0x0),
fh1dTracksImpParXYResidualTruth(0x0),
fh1dTracksImpParXYZResidualTruth(0x0),
fh2dVertexChi2NDFNESDTracks(0x0),
fh1dJetGenPt(0x0),
fh1dJetGenPtUnidentified(0x0),
fh1dJetGenPtudsg(0x0),
fh1dJetGenPtc(0x0),
fh1dJetGenPtb(0x0),
fh1dJetRecPt(0x0),
fh1dPhotonPt(0x0),
fh2dKshortMassVsPt(0x0),
fh2dLamdaMassVsPt(0x0),
fh2dAnLamdaMassVsPt(0x0),
fh2dKshortMassVsPtReal(0x0),
fh2dLamdaMassVsPtReal(0x0),
fh2dAnLamdaMassVsPtReal(0x0),
fh2dKshortRecPtVsGenPt(0x0),
fh2dLamdaRecPtVsGenPt(0x0),
fh2dAnLamdaRecPtVsGenPt(0x0),
fh1dKshortPtMC(0x0),
fh1dLamdaPtMC(0x0),
fh1dAnLamdaPtMC(0x0),
fh2dKshortPtVsJetPtMC(0x0),
fh2dLamdaPtVsJetPtMC(0x0),
fh2dAnLamdaPtVsJetPtMC(0x0),
fh1dJetRecPtAccepted(0x0),
fhnV0InJetK0s(0x0),
fhnV0InJetLambda(0x0),
fhnV0InJetALambda(0x0),
fh1dJetRecPtAcceptedunCorr(0x0),
f2histRhoVsDeltaPt(0x0),
f2histRhoVsDeltaPtFirst(0x0),
f2histRhoVsDeltaPtSecond(0x0),
f2histRhoVsDeltaPtThird(0x0),
f2histRhoVsDeltaPtWithSignal(0x0),
f2histRhoVsDeltaPtWithSignalFirst(0x0),
f2histRhoVsDeltaPtWithSignalSecond(0x0),
f2histRhoVsDeltaPtWithSignalThird(0x0),
fRandom(new TRandom3(0)),
fh1dJetRecEtaPhiAccepted(0x0),
fh1dJetRecPtUnidentified(0x0),
fh1dJetRecPtudsg(0x0),
fh1dJetRecPtc(0x0),
fh1dJetRecPtb(0x0),
fh1dJetRecPtUnidentifiedAccepted(0x0),
fh1dJetRecPtudsgAccepted(0x0),
fh1dJetRecPtcAccepted(0x0),
fh1dJetRecPtbAccepted(0x0),
fDoTaggedDRM(kFALSE),
fh2dJetGenPtVsJetRecPt(0x0),
fh2dJetGenPtVsJetRecPtFirst(0x0),
fh2dJetGenPtVsJetRecPtSecond(0x0),
fh2dJetGenPtVsJetRecPtThird(0x0),
fh2dJetGenPtVsJetRecPtb(0x0),
fh2dJetGenPtVsJetRecPtc(0x0),
fh2dJetGenPtVsJetRecPtudsg(0x0),
//PtRel
fhistPtRelEvents(0x0),
fCaloClusters(0x0),
fhistPtRelVsJetPt(0x0),
fhistLepIPVsJetPt(0x0),
fhistPtRelVsJetPtUnidentified(0x0),
fhistPtRelVsJetPtudsg(0x0),
fhistPtRelVsJetPtc(0x0),
fhistPtRelVsJetPtb(0x0),
fhistLepIPVsJetPtUnidentified(0x0),
fhistLepIPVsJetPtudsg(0x0),
fhistLepIPVsJetPtc(0x0),
fhistLepIPVsJetPtb(0x0),
fHistMcEopEle(0x0),
fHistMcEopHad(0x0),
fTPCnsigMcEle(0x0),
fTPCnsigMcHad(0x0),
fhistPtRelVsJetPtTaggedFirst(0x0),
fhistLepIPVsJetPtTaggedFirst(0x0),
fhistPtRelVsJetPtTaggedUnidentifiedFirst(0x0),
fhistLepIPVsJetPtTaggedUnidentifiedFirst(0x0),
fhistPtRelVsJetPtTaggedudsgFirst(0x0),
fhistLepIPVsJetPtTaggedudsgFirst(0x0),
fhistPtRelVsJetPtTaggedcFirst(0x0),
fhistLepIPVsJetPtTaggedcFirst(0x0),
fhistPtRelVsJetPtTaggedbFirst(0x0),
fhistLepIPVsJetPtTaggedbFirst(0x0),
fhistPtRelVsJetPtTaggedSecond(0x0),
fhistLepIPVsJetPtTaggedSecond(0x0),
fhistPtRelVsJetPtTaggedUnidentifiedSecond(0x0),
fhistLepIPVsJetPtTaggedUnidentifiedSecond(0x0),
fhistPtRelVsJetPtTaggedudsgSecond(0x0),
fhistLepIPVsJetPtTaggedudsgSecond(0x0),
fhistPtRelVsJetPtTaggedcSecond(0x0),
fhistLepIPVsJetPtTaggedcSecond(0x0),
fhistPtRelVsJetPtTaggedbSecond(0x0),
fhistLepIPVsJetPtTaggedbSecond(0x0),
fhistPtRelVsJetPtTaggedThird(0x0),
fhistLepIPVsJetPtTaggedThird(0x0),
fhistPtRelVsJetPtTaggedUnidentifiedThird(0x0),
fhistLepIPVsJetPtTaggedUnidentifiedThird(0x0),
fhistPtRelVsJetPtTaggedudsgThird(0x0),
fhistLepIPVsJetPtTaggedudsgThird(0x0),
fhistPtRelVsJetPtTaggedcThird(0x0),
fhistLepIPVsJetPtTaggedcThird(0x0),
fhistPtRelVsJetPtTaggedbThird(0x0),
fhistLepIPVsJetPtTaggedbThird(0x0),
//__
fh2dJetSignedImpParXY(0x0),
fh2dJetSignedImpParXYUnidentified(0x0),
fh2dJetSignedImpParXYudsg(0x0),
fh2dJetSignedImpParXYb(0x0),
fh2dJetSignedImpParXYc(0x0),
fh2dJetSignedImpParXYSignificance(0x0),
fh2dJetSignedImpParXYSignificanceUnidentified(0x0),
fh2dJetSignedImpParXYSignificanceudsg(0x0),
fh2dJetSignedImpParXYSignificanceb(0x0),
fh2dJetSignedImpParXYSignificancec(0x0),
fh2dJetSignedImpParXYZ(0x0),
fh2dJetSignedImpParXYZUnidentified(0x0),
fh2dJetSignedImpParXYZudsg(0x0),
fh2dJetSignedImpParXYZb(0x0),
fh2dJetSignedImpParXYZc(0x0),
fh2dJetSignedImpParXYZSignificance(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentified(0x0),
fh2dJetSignedImpParXYZSignificanceudsg(0x0),
fh2dJetSignedImpParXYZSignificanceb(0x0),
fh2dJetSignedImpParXYZSignificancec(0x0),
fh2dJetSignedImpParXYFirst(0x0),
fh2dJetSignedImpParXYUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYudsgFirst(0x0),
fh2dJetSignedImpParXYbFirst(0x0),
fh2dJetSignedImpParXYcFirst(0x0),
fh2dJetSignedImpParXYSignificanceFirst(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYSignificanceudsgFirst(0x0),
fh2dJetSignedImpParXYSignificancebFirst(0x0),
fh2dJetSignedImpParXYSignificancecFirst(0x0),
fh2dJetSignedImpParXYZFirst(0x0),
fh2dJetSignedImpParXYZUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYZudsgFirst(0x0),
fh2dJetSignedImpParXYZbFirst(0x0),
fh2dJetSignedImpParXYZcFirst(0x0),
fh2dJetSignedImpParXYZSignificanceFirst(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYZSignificanceudsgFirst(0x0),
fh2dJetSignedImpParXYZSignificancebFirst(0x0),
fh2dJetSignedImpParXYZSignificancecFirst(0x0),
fh2dJetSignedImpParXYSecond(0x0),
fh2dJetSignedImpParXYUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYudsgSecond(0x0),
fh2dJetSignedImpParXYbSecond(0x0),
fh2dJetSignedImpParXYcSecond(0x0),
fh2dJetSignedImpParXYSignificanceSecond(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYSignificanceudsgSecond(0x0),
fh2dJetSignedImpParXYSignificancebSecond(0x0),
fh2dJetSignedImpParXYSignificancecSecond(0x0),
fh2dJetSignedImpParXYZSecond(0x0),
fh2dJetSignedImpParXYZUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYZudsgSecond(0x0),
fh2dJetSignedImpParXYZbSecond(0x0),
fh2dJetSignedImpParXYZcSecond(0x0),
fh2dJetSignedImpParXYZSignificanceSecond(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYZSignificanceudsgSecond(0x0),
fh2dJetSignedImpParXYZSignificancebSecond(0x0),
fh2dJetSignedImpParXYZSignificancecSecond(0x0),
fh2dJetSignedImpParXYThird(0x0),
fh2dJetSignedImpParXYUnidentifiedThird(0x0),
fh2dJetSignedImpParXYudsgThird(0x0),
fh2dJetSignedImpParXYbThird(0x0),
fh2dJetSignedImpParXYcThird(0x0),
fh2dJetSignedImpParXYSignificanceThird(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedThird(0x0),
fh2dJetSignedImpParXYSignificanceudsgThird(0x0),
fh2dJetSignedImpParXYSignificancebThird(0x0),
fh2dJetSignedImpParXYSignificancecThird(0x0),
fh2dJetSignedImpParXYZThird(0x0),
fh2dJetSignedImpParXYZUnidentifiedThird(0x0),
fh2dJetSignedImpParXYZudsgThird(0x0),
fh2dJetSignedImpParXYZbThird(0x0),
fh2dJetSignedImpParXYZcThird(0x0),
fh2dJetSignedImpParXYZSignificanceThird(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedThird(0x0),
fh2dJetSignedImpParXYZSignificanceudsgThird(0x0),
fh2dJetSignedImpParXYZSignificancebThird(0x0),
fh2dJetSignedImpParXYZSignificancecThird(0x0),
//forth
fDoForthIP(kFALSE),
fh2dJetSignedImpParXYForth(0x0),
fh2dJetSignedImpParXYudsgForth(0x0),
fh2dJetSignedImpParXYbForth(0x0),
fh2dJetSignedImpParXYcForth(0x0),
fh2dJetSignedImpParXYSignificanceForth(0x0),
fh2dJetSignedImpParXYSignificanceudsgForth(0x0),
fh2dJetSignedImpParXYSignificancebForth(0x0),
fh2dJetSignedImpParXYSignificancecForth(0x0),
//Jet Probabilty
fh2dJetSignedImpParXY_Class1(0x0),
fh2dJetSignedImpParXYSignificance_Class1(0x0),
fh2dJetSignedImpParXYSignificanceb_Class1(0x0),
fh2dJetSignedImpParXYSignificancec_Class1(0x0),
fh2dJetSignedImpParXYSignificancelf_Class1(0x0),
fh2dJetSignedImpParXYZ_Class1(0x0),
fh2dJetSignedImpParXYZSignificance_Class1(0x0),
fh2dJetSignedImpParXY_Class2(0x0),
fh2dJetSignedImpParXYSignificance_Class2(0x0),
fh2dJetSignedImpParXYSignificanceb_Class2(0x0),
fh2dJetSignedImpParXYSignificancec_Class2(0x0),
fh2dJetSignedImpParXYSignificancelf_Class2(0x0),
fh2dJetSignedImpParXYZ_Class2(0x0),
fh2dJetSignedImpParXYZSignificance_Class2(0x0),
fh2dJetSignedImpParXY_Class3(0x0),
fh2dJetSignedImpParXYSignificance_Class3(0x0),
fh2dJetSignedImpParXYSignificanceb_Class3(0x0),
fh2dJetSignedImpParXYSignificancec_Class3(0x0),
fh2dJetSignedImpParXYSignificancelf_Class3(0x0),
fh2dJetSignedImpParXYZ_Class3(0x0),
fh2dJetSignedImpParXYZSignificance_Class3(0x0),
fh2dJetSignedImpParXY_Class4(0x0),
fh2dJetSignedImpParXYSignificance_Class4(0x0),
fh2dJetSignedImpParXYSignificanceb_Class4(0x0),
fh2dJetSignedImpParXYSignificancec_Class4(0x0),
fh2dJetSignedImpParXYSignificancelf_Class4(0x0),
fh2dJetSignedImpParXYZ_Class4(0x0),
fh2dJetSignedImpParXYZSignificance_Class4(0x0),
//Jet Mass
fhistJetMass(0x0),
fhistJetMass_Unidentified(0x0),
fhistJetMass_udsg(0x0),
fhistJetMass_c(0x0),
fhistJetMass_b(0x0),
fhistJetMassFirst(0x0),
fhistJetMassSecond(0x0),
fhistJetMassThird(0x0),
fhistJetMass_UnidentifiedFirst(0x0),
fhistJetMass_udsgFirst(0x0),
fhistJetMass_cFirst(0x0),
fhistJetMass_bFirst(0x0),
fhistJetMass_UnidentifiedSecond(0x0),
fhistJetMass_udsgSecond(0x0),
fhistJetMass_cSecond(0x0),
fhistJetMass_bSecond(0x0),
fhistJetMass_UnidentifiedThird(0x0),
fhistJetMass_udsgThird(0x0),
fhistJetMass_cThird(0x0),
fhistJetMass_bThird(0x0),
//Secondary vertex energy fraction
fhistSVEnergyFraction(0x0),
fhistSVEnergyFraction_Unidentified(0x0),
fhistSVEnergyFraction_udsg(0x0),
fhistSVEnergyFraction_c(0x0),
fhistSVEnergyFraction_b(0x0),
fhistSVEnergyFractionFirst(0x0),
fhistSVEnergyFraction_UnidentifiedFirst(0x0),
fhistSVEnergyFraction_udsgFirst(0x0),
fhistSVEnergyFraction_cFirst(0x0),
fhistSVEnergyFraction_bFirst(0x0),
fhistSVEnergyFractionSecond(0x0),
fhistSVEnergyFraction_UnidentifiedSecond(0x0),
fhistSVEnergyFraction_udsgSecond(0x0),
fhistSVEnergyFraction_cSecond(0x0),
fhistSVEnergyFraction_bSecond(0x0),
fhistSVEnergyFractionThird(0x0),
fhistSVEnergyFraction_UnidentifiedThird(0x0),
fhistSVEnergyFraction_udsgThird(0x0),
fhistSVEnergyFraction_cThird(0x0),
fhistSVEnergyFraction_bThird(0x0),
fhistSVnProngs(0x0),
fhistSVnProngs_Unidentified(0x0),
fhistSVnProngs_udsg(0x0),
fhistSVnProngs_c(0x0),
fhistSVnProngs_b(0x0),
//Jet Probability
fhistJetProbability(0x0),
fhistJetProbability_Unidentified(0x0),
fhistJetProbability_udsg(0x0),
fhistJetProbability_c(0x0),
fhistJetProbability_b(0x0),
fhistJetProbabilityLog(0x0),
fhistJetProbability_UnidentifiedLog(0x0),
fhistJetProbability_udsgLog(0x0),
fhistJetProbability_cLog(0x0),
fhistJetProbability_cLog_D0(0x0),
fhistJetProbability_cLog_Dp(0x0),
fhistJetProbability_cLog_Ds(0x0),
fhistJetProbability_cLog_Lc(0x0),
fhistJetProbability_bLog(0x0),
fhistJetProbabilityLogFirst(0x0),
fhistJetProbability_UnidentifiedLogFirst(0x0),
fhistJetProbability_udsgLogFirst(0x0),
fhistJetProbability_cLogFirst(0x0),
fhistJetProbability_cLogFirst_D0(0x0),
fhistJetProbability_cLogFirst_Dp(0x0),
fhistJetProbability_cLogFirst_Ds(0x0),
fhistJetProbability_cLogFirst_Lc(0x0),
fhistJetProbability_bLogFirst(0x0),
fhistJetProbabilityLogSecond(0x0),
fhistJetProbability_UnidentifiedLogSecond(0x0),
fhistJetProbability_udsgLogSecond(0x0),
fhistJetProbability_cLogSecond(0x0),
fhistJetProbability_cLogSecond_D0(0x0),
fhistJetProbability_cLogSecond_Dp(0x0),
fhistJetProbability_cLogSecond_Ds(0x0),
fhistJetProbability_cLogSecond_Lc(0x0),
fhistJetProbability_bLogSecond(0x0),
fhistJetProbabilityLogThird(0x0),
fhistJetProbability_UnidentifiedLogThird(0x0),
fhistJetProbability_udsgLogThird(0x0),
fhistJetProbability_cLogThird(0x0),
fhistJetProbability_cLogThird_D0(0x0),
fhistJetProbability_cLogThird_Dp(0x0),
fhistJetProbability_cLogThird_Ds(0x0),
fhistJetProbability_cLogThird_Lc(0x0),
fhistJetProbability_bLogThird(0x0),
fhistJetProbabilityLogSVHE(0x0),
fhistJetProbability_UnidentifiedLogSVHE(0x0),
fhistJetProbability_udsgLogSVHE(0x0),
fhistJetProbability_cLogSVHE(0x0),
fhistJetProbability_bLogSVHE(0x0),
fhistJetProbabilityLogSVHP(0x0),
fhistJetProbability_UnidentifiedLogSVHP(0x0),
fhistJetProbability_udsgLogSVHP(0x0),
fhistJetProbability_cLogSVHP(0x0),
fhistJetProbability_bLogSVHP(0x0),
fMinTrackProb(0.0),
//__________
// V0Reconstruction
fh1V0CounterCentK0s(0x0),
fh1V0CounterCentLambda(0x0),
fh1V0CounterCentALambda(0x0),
fbTPCRefit(0),
fbRejectKinks(0),
fbFindableClusters(0),
fdCutNCrossedRowsTPCMin(-1),
fdCutCrossedRowsOverFindMin(-1),
fdCutCrossedRowsOverFindMax(-1),
fdCutPtDaughterMin(-1),
fdCutDCAToPrimVtxMin(-1),
fdCutDCADaughtersMax(-1),
fdCutEtaDaughterMax(-1),
fdCutNSigmadEdxMax(-1),
fdPtProtonPIDMax(-1),
fbOnFly(0),
fdCutCPAKMin(-1),
fdCutCPALMin(-1),
fdCutRadiusDecayMin(-1),
fdCutRadiusDecayMax(-1),
fdCutEtaV0Max(-1),
fdCutRapV0Max(-1),
fdCutNTauKMax(-1),
fdCutNTauLMax(-1),
fbCutArmPod(0),
fbCutCross(0),
fApplyV0Rec(kFALSE),
fApplyV0RejectionAll(kFALSE),
//______
fMCArray(0x0),
fUseCorrPt(kTRUE),
fUsePicoTracks(kTRUE),
fEnableV0GammaRejection(0),
fV0CandidateArray(0x0),
fUtils(new AliAnalysisUtils()),
fJetContainerMC(0x0),
fJetContainerData(0x0),
fAODIn(0x0),
fPrimaryVertex(0x0),
//SV Analysis
fDoSVAnalysis(kFALSE),
fDoTrackCountingAnalysis(kTRUE),
fVtxTagger3Prong(0x0),
fVtxTagger2Prong(0x0),
fjetCuts3Prong(0x0),
fjetCuts2Prong(0x0),
fTrackArray(0x0),
fEsdTrackCuts(0x0),
fInvariantMass(0.),
fJetMass(0.),
fDispersion(0.),
fDecayLength(0.),
fLxySign(0.),
fJetPt(0.),
fJetFlavor(0),
fValJetProb(-1),
fLogJetProb(-1.),
fCalcDCATruth(kFALSE),
fDecayVertex(0x0),
fHistSV2Prong(0x0),
fHistSV2ProngUnidentified(0x0),
fHistSV2Prongb(0x0),
fHistSV2Prongc(0x0),
fHistSV2Pronglf(0x0),
fHistDispersion2Prong(0x0),
fHistDispersion2ProngUnidentified(0x0),
fHistDispersion2Prongb(0x0),
fHistDispersion2Prongc(0x0),
fHistDispersion2Pronglf(0x0),
fHistSV3Prong(0x0),
fHistSV3ProngUnidentified(0x0),
fHistSV3Prongb(0x0),
fHistSV3Prongc(0x0),
fHistSV3Pronglf(0x0),
fHistDispersion3Prong(0x0),
fHistDispersion3ProngUnidentified(0x0),
fHistDispersion3Prongb(0x0),
fHistDispersion3Prongc(0x0),
fHistDispersion3Pronglf(0x0)
{
	SetMakeGeneralHistograms(kTRUE);

	for(int i=0; i<7; i++){
		fResolutionFunction[i]=0x0;
		fResolutionFunctionb[i]=0x0;
		fResolutionFunctionc[i]=0x0;
		fResolutionFunctionlf[i]=0x0;
	}
}
// ######################################################################################## CONSTRUCTORS
AliAnalysisTaskBJetTC::AliAnalysisTaskBJetTC(const char *name): AliAnalysisTaskEmcalJet(name, kTRUE),
  		fV0Reader(NULL),
  		fV0ReaderName("V0ReaderV1"),  fReaderGammas(NULL),
		fHFJetUtils(0x0),
		fRespoPID(0x0),
		fPtHardThreshold(0.0),
		fPythiaEventWeight(1.0),
		fDoImprovedDCACut(kTRUE),
		fVertexConstraint(kFALSE),
		fThresholdIP(0.008),
		fDoDeltaPtWithSignal(kFALSE),
		fDiamond(0x0),
		fVertexer(0x0),
		//Bjet Cuts
		fTCMinTrackPt(0.5),
		fTCMinClusTPC(80),
		fTCMinHitsITS(2),
		fTCMaxChi2pNDF(5.),
		fTCMaxIPxy(1.),
		fTCMaxIPz(5.),
		fTCMaxDecayLength(5),
		fTCMaxDCATrackJet(0.07),
		fMaxFactorPtHardJet(10.0),
		fhistInclusiveJetCuts(0x0),
		fhistbJetCuts(0x0),
		fhistcJetCuts(0x0),
		fhistlfJetCuts(0x0),
		//____
		fh1dPhotonPt(0x0),
		fh2dKshortMassVsPt(0x0),
		fh2dLamdaMassVsPt(0x0),
		fh2dAnLamdaMassVsPt(0x0),
		fh2dKshortMassVsPtReal(0x0),
		fh2dLamdaMassVsPtReal(0x0),
		fh2dAnLamdaMassVsPtReal(0x0),
		fh2dKshortRecPtVsGenPt(0x0),
		fh2dLamdaRecPtVsGenPt(0x0),
		fh2dAnLamdaRecPtVsGenPt(0x0),
		fh1dKshortPtMC(0x0),
		fh1dLamdaPtMC(0x0),
		fh1dAnLamdaPtMC(0x0),
		fh2dKshortPtVsJetPtMC(0x0),
		fh2dLamdaPtVsJetPtMC(0x0),
		fh2dAnLamdaPtVsJetPtMC(0x0),
		fhnV0InJetK0s(0x0),
		fhnV0InJetLambda(0x0),
		fhnV0InJetALambda(0x0),
		fDoJetProbabilityAnalysis(kFALSE),
		fDoCharmFractions(kFALSE),
		fUsePartonDef(kTRUE),
		fUseIPs(kFALSE),
		fDoJetMass(kFALSE),
		fDoSVEnergyFraction(kFALSE),
		fDoPtRelAnalysis(0),
		fDoSelectionPtRel(0),
		fh1dEventRejectionRDHFCuts(0x0),
		fh1dVertexZ(0x0),
		fh1dVertexZAccepted(0x0),
		fh1dVertexR(0x0),
		fh1dVertexRAccepted(0x0),
		fh1dTracksAccepeted(0x0),
		fh1dTracksImpParXY(0x0),
		fh1dTracksImpParXYZ(0x0),
		fh1dTracksImpParXYSignificance(0x0),
		fh1dTracksImpParXYZSignificance(0x0),
		fh1dTracksImpParXYTruth(0x0),
		fh1dTracksImpParXYZTruth(0x0),
		fh1dTracksImpParXYResidualTruth(0x0),
		fh1dTracksImpParXYZResidualTruth(0x0),
		fh2dVertexChi2NDFNESDTracks(0x0),
		fh1dJetGenPt(0x0),
		fh1dJetGenPtUnidentified(0x0),
		fh1dJetGenPtudsg(0x0),
		fh1dJetGenPtc(0x0),
		fh1dJetGenPtb(0x0),
		fh1dJetRecPt(0x0),
		fh1dJetRecPtAccepted(0x0),
		fh1dJetRecPtAcceptedunCorr(0x0),
 		fRandom(new TRandom3(0)),
		f2histRhoVsDeltaPt(0x0),
		f2histRhoVsDeltaPtFirst(0x0),
		f2histRhoVsDeltaPtSecond(0x0),
		f2histRhoVsDeltaPtThird(0x0),
		f2histRhoVsDeltaPtWithSignal(0x0),
		f2histRhoVsDeltaPtWithSignalFirst(0x0),
		f2histRhoVsDeltaPtWithSignalSecond(0x0),
		f2histRhoVsDeltaPtWithSignalThird(0x0),
		fh1dJetRecEtaPhiAccepted(0x0),
		fh1dJetRecPtUnidentified(0x0),
		fh1dJetRecPtudsg(0x0),
		fh1dJetRecPtc(0x0),
		fh1dJetRecPtb(0x0),
		fh1dJetRecPtUnidentifiedAccepted(0x0),
		fh1dJetRecPtudsgAccepted(0x0),
		fh1dJetRecPtcAccepted(0x0),
		fh1dJetRecPtbAccepted(0x0),
		fDoTaggedDRM(kFALSE),
		fh2dJetGenPtVsJetRecPt(0x0),
		fh2dJetGenPtVsJetRecPtFirst(0x0),
		fh2dJetGenPtVsJetRecPtSecond(0x0),
		fh2dJetGenPtVsJetRecPtThird(0x0),
		fh2dJetGenPtVsJetRecPtb(0x0),
		fh2dJetGenPtVsJetRecPtc(0x0),
		fh2dJetGenPtVsJetRecPtudsg(0x0),
		//PtRel
		fhistPtRelEvents(0x0),
		fCaloClusters(0x0),
		fhistPtRelVsJetPt(0x0),
		fhistLepIPVsJetPt(0x0),
		fhistPtRelVsJetPtUnidentified(0x0),
		fhistPtRelVsJetPtudsg(0x0),
		fhistPtRelVsJetPtc(0x0),
		fhistPtRelVsJetPtb(0x0),
		fhistLepIPVsJetPtUnidentified(0x0),
		fhistLepIPVsJetPtudsg(0x0),
		fhistLepIPVsJetPtc(0x0),
		fhistLepIPVsJetPtb(0x0),
		fHistMcEopEle(0x0),
		fHistMcEopHad(0x0),
		fTPCnsigMcEle(0x0),
		fTPCnsigMcHad(0x0),
		fhistPtRelVsJetPtTaggedFirst(0x0),
		fhistLepIPVsJetPtTaggedFirst(0x0),
		fhistPtRelVsJetPtTaggedUnidentifiedFirst(0x0),
		fhistLepIPVsJetPtTaggedUnidentifiedFirst(0x0),
		fhistPtRelVsJetPtTaggedudsgFirst(0x0),
		fhistLepIPVsJetPtTaggedudsgFirst(0x0),
		fhistPtRelVsJetPtTaggedcFirst(0x0),
		fhistLepIPVsJetPtTaggedcFirst(0x0),
		fhistPtRelVsJetPtTaggedbFirst(0x0),
		fhistLepIPVsJetPtTaggedbFirst(0x0),
		fhistPtRelVsJetPtTaggedSecond(0x0),
		fhistLepIPVsJetPtTaggedSecond(0x0),
		fhistPtRelVsJetPtTaggedUnidentifiedSecond(0x0),
		fhistLepIPVsJetPtTaggedUnidentifiedSecond(0x0),
		fhistPtRelVsJetPtTaggedudsgSecond(0x0),
		fhistLepIPVsJetPtTaggedudsgSecond(0x0),
		fhistPtRelVsJetPtTaggedcSecond(0x0),
		fhistLepIPVsJetPtTaggedcSecond(0x0),
		fhistPtRelVsJetPtTaggedbSecond(0x0),
		fhistLepIPVsJetPtTaggedbSecond(0x0),
		fhistPtRelVsJetPtTaggedThird(0x0),
		fhistLepIPVsJetPtTaggedThird(0x0),
		fhistPtRelVsJetPtTaggedUnidentifiedThird(0x0),
		fhistLepIPVsJetPtTaggedUnidentifiedThird(0x0),
		fhistPtRelVsJetPtTaggedudsgThird(0x0),
		fhistLepIPVsJetPtTaggedudsgThird(0x0),
		fhistPtRelVsJetPtTaggedcThird(0x0),
		fhistLepIPVsJetPtTaggedcThird(0x0),
		fhistPtRelVsJetPtTaggedbThird(0x0),
		fhistLepIPVsJetPtTaggedbThird(0x0),
		//___
		fh2dJetSignedImpParXY(0x0),
		fh2dJetSignedImpParXYUnidentified(0x0),
		fh2dJetSignedImpParXYudsg(0x0),
		fh2dJetSignedImpParXYb(0x0),
		fh2dJetSignedImpParXYc(0x0),
		fh2dJetSignedImpParXYSignificance(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentified(0x0),
		fh2dJetSignedImpParXYSignificanceudsg(0x0),
		fh2dJetSignedImpParXYSignificanceb(0x0),
		fh2dJetSignedImpParXYSignificancec(0x0),
		fh2dJetSignedImpParXYZ(0x0),
		fh2dJetSignedImpParXYZUnidentified(0x0),
		fh2dJetSignedImpParXYZudsg(0x0),
		fh2dJetSignedImpParXYZb(0x0),
		fh2dJetSignedImpParXYZc(0x0),
		fh2dJetSignedImpParXYZSignificance(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentified(0x0),
		fh2dJetSignedImpParXYZSignificanceudsg(0x0),
		fh2dJetSignedImpParXYZSignificanceb(0x0),
		fh2dJetSignedImpParXYZSignificancec(0x0),
		fh2dJetSignedImpParXYFirst(0x0),
		fh2dJetSignedImpParXYUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYudsgFirst(0x0),
		fh2dJetSignedImpParXYbFirst(0x0),
		fh2dJetSignedImpParXYcFirst(0x0),
		fh2dJetSignedImpParXYSignificanceFirst(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYSignificanceudsgFirst(0x0),
		fh2dJetSignedImpParXYSignificancebFirst(0x0),
		fh2dJetSignedImpParXYSignificancecFirst(0x0),
		fh2dJetSignedImpParXYZFirst(0x0),
		fh2dJetSignedImpParXYZUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYZudsgFirst(0x0),
		fh2dJetSignedImpParXYZbFirst(0x0),
		fh2dJetSignedImpParXYZcFirst(0x0),
		fh2dJetSignedImpParXYZSignificanceFirst(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgFirst(0x0),
		fh2dJetSignedImpParXYZSignificancebFirst(0x0),
		fh2dJetSignedImpParXYZSignificancecFirst(0x0),
		fh2dJetSignedImpParXYSecond(0x0),
		fh2dJetSignedImpParXYUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYudsgSecond(0x0),
		fh2dJetSignedImpParXYbSecond(0x0),
		fh2dJetSignedImpParXYcSecond(0x0),
		fh2dJetSignedImpParXYSignificanceSecond(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYSignificanceudsgSecond(0x0),
		fh2dJetSignedImpParXYSignificancebSecond(0x0),
		fh2dJetSignedImpParXYSignificancecSecond(0x0),
		fh2dJetSignedImpParXYZSecond(0x0),
		fh2dJetSignedImpParXYZUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYZudsgSecond(0x0),
		fh2dJetSignedImpParXYZbSecond(0x0),
		fh2dJetSignedImpParXYZcSecond(0x0),
		fh2dJetSignedImpParXYZSignificanceSecond(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgSecond(0x0),
		fh2dJetSignedImpParXYZSignificancebSecond(0x0),
		fh2dJetSignedImpParXYZSignificancecSecond(0x0),
		fh2dJetSignedImpParXYThird(0x0),
		fh2dJetSignedImpParXYUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYudsgThird(0x0),
		fh2dJetSignedImpParXYbThird(0x0),
		fh2dJetSignedImpParXYcThird(0x0),
		fh2dJetSignedImpParXYSignificanceThird(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYSignificanceudsgThird(0x0),
		fh2dJetSignedImpParXYSignificancebThird(0x0),
		fh2dJetSignedImpParXYSignificancecThird(0x0),
		fh2dJetSignedImpParXYZThird(0x0),
		fh2dJetSignedImpParXYZUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYZudsgThird(0x0),
		fh2dJetSignedImpParXYZbThird(0x0),
		fh2dJetSignedImpParXYZcThird(0x0),
		fh2dJetSignedImpParXYZSignificanceThird(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgThird(0x0),
		fh2dJetSignedImpParXYZSignificancebThird(0x0),
		fh2dJetSignedImpParXYZSignificancecThird(0x0),
		//forth
		fDoForthIP(kFALSE),
		fh2dJetSignedImpParXYForth(0x0),
		fh2dJetSignedImpParXYudsgForth(0x0),
		fh2dJetSignedImpParXYbForth(0x0),
		fh2dJetSignedImpParXYcForth(0x0),
		fh2dJetSignedImpParXYSignificanceForth(0x0),
		fh2dJetSignedImpParXYSignificanceudsgForth(0x0),
		fh2dJetSignedImpParXYSignificancebForth(0x0),
		fh2dJetSignedImpParXYSignificancecForth(0x0),
		//Jet Mass
		fhistJetMass(0x0),
		fhistJetMass_Unidentified(0x0),
		fhistJetMass_udsg(0x0),
		fhistJetMass_c(0x0),
		fhistJetMass_b(0x0),
		fhistJetMassFirst(0x0),
		fhistJetMassSecond(0x0),
		fhistJetMassThird(0x0),
		fhistJetMass_UnidentifiedFirst(0x0),
		fhistJetMass_udsgFirst(0x0),
		fhistJetMass_cFirst(0x0),
		fhistJetMass_bFirst(0x0),
		fhistJetMass_UnidentifiedSecond(0x0),
		fhistJetMass_udsgSecond(0x0),
		fhistJetMass_cSecond(0x0),
		fhistJetMass_bSecond(0x0),
		fhistJetMass_UnidentifiedThird(0x0),
		fhistJetMass_udsgThird(0x0),
		fhistJetMass_cThird(0x0),
		fhistJetMass_bThird(0x0),
		//Secondary vertex energy fraction
		fhistSVEnergyFraction(0x0),
		fhistSVEnergyFraction_Unidentified(0x0),
		fhistSVEnergyFraction_udsg(0x0),
		fhistSVEnergyFraction_c(0x0),
		fhistSVEnergyFraction_b(0x0),
		fhistSVEnergyFractionFirst(0x0),
		fhistSVEnergyFraction_UnidentifiedFirst(0x0),
		fhistSVEnergyFraction_udsgFirst(0x0),
		fhistSVEnergyFraction_cFirst(0x0),
		fhistSVEnergyFraction_bFirst(0x0),
		fhistSVEnergyFractionSecond(0x0),
		fhistSVEnergyFraction_UnidentifiedSecond(0x0),
		fhistSVEnergyFraction_udsgSecond(0x0),
		fhistSVEnergyFraction_cSecond(0x0),
		fhistSVEnergyFraction_bSecond(0x0),
		fhistSVEnergyFractionThird(0x0),
		fhistSVEnergyFraction_UnidentifiedThird(0x0),
		fhistSVEnergyFraction_udsgThird(0x0),
		fhistSVEnergyFraction_cThird(0x0),
		fhistSVEnergyFraction_bThird(0x0),
		fhistSVnProngs(0x0),
		fhistSVnProngs_Unidentified(0x0),
		fhistSVnProngs_udsg(0x0),
		fhistSVnProngs_c(0x0),
		fhistSVnProngs_b(0x0),
		//Jet Probability
		fh2dJetSignedImpParXY_Class1(0x0),
		fh2dJetSignedImpParXYSignificance_Class1(0x0),
		fh2dJetSignedImpParXYSignificanceb_Class1(0x0),
		fh2dJetSignedImpParXYSignificancec_Class1(0x0),
		fh2dJetSignedImpParXYSignificancelf_Class1(0x0),
		fh2dJetSignedImpParXYZ_Class1(0x0),
		fh2dJetSignedImpParXYZSignificance_Class1(0x0),
		fh2dJetSignedImpParXY_Class2(0x0),
		fh2dJetSignedImpParXYSignificance_Class2(0x0),
		fh2dJetSignedImpParXYSignificanceb_Class2(0x0),
		fh2dJetSignedImpParXYSignificancec_Class2(0x0),
		fh2dJetSignedImpParXYSignificancelf_Class2(0x0),
		fh2dJetSignedImpParXYZ_Class2(0x0),
		fh2dJetSignedImpParXYZSignificance_Class2(0x0),
		fh2dJetSignedImpParXY_Class3(0x0),
		fh2dJetSignedImpParXYSignificance_Class3(0x0),
		fh2dJetSignedImpParXYSignificanceb_Class3(0x0),
		fh2dJetSignedImpParXYSignificancec_Class3(0x0),
		fh2dJetSignedImpParXYSignificancelf_Class3(0x0),
		fh2dJetSignedImpParXYZ_Class3(0x0),
		fh2dJetSignedImpParXYZSignificance_Class3(0x0),
		fh2dJetSignedImpParXY_Class4(0x0),
		fh2dJetSignedImpParXYSignificance_Class4(0x0),
		fh2dJetSignedImpParXYSignificanceb_Class4(0x0),
		fh2dJetSignedImpParXYSignificancec_Class4(0x0),
		fh2dJetSignedImpParXYSignificancelf_Class4(0x0),
		fh2dJetSignedImpParXYZ_Class4(0x0),
		fh2dJetSignedImpParXYZSignificance_Class4(0x0),
		fhistJetProbability(0x0),
		fhistJetProbability_Unidentified(0x0),
		fhistJetProbability_udsg(0x0),
		fhistJetProbability_c(0x0),
		fhistJetProbability_b(0x0),
		fhistJetProbabilityLog(0x0),
		fhistJetProbability_UnidentifiedLog(0x0),
		fhistJetProbability_udsgLog(0x0),
		fhistJetProbability_cLog(0x0),
		fhistJetProbability_cLog_D0(0x0),
		fhistJetProbability_cLog_Dp(0x0),
		fhistJetProbability_cLog_Ds(0x0),
		fhistJetProbability_cLog_Lc(0x0),
		fhistJetProbability_bLog(0x0),
		fhistJetProbabilityLogFirst(0x0),
		fhistJetProbability_UnidentifiedLogFirst(0x0),
		fhistJetProbability_udsgLogFirst(0x0),
		fhistJetProbability_cLogFirst(0x0),
		fhistJetProbability_cLogFirst_D0(0x0),
		fhistJetProbability_cLogFirst_Dp(0x0),
		fhistJetProbability_cLogFirst_Ds(0x0),
		fhistJetProbability_cLogFirst_Lc(0x0),
		fhistJetProbability_bLogFirst(0x0),
		fhistJetProbabilityLogSecond(0x0),
		fhistJetProbability_UnidentifiedLogSecond(0x0),
		fhistJetProbability_udsgLogSecond(0x0),
		fhistJetProbability_cLogSecond(0x0),
		fhistJetProbability_cLogSecond_D0(0x0),
		fhistJetProbability_cLogSecond_Dp(0x0),
		fhistJetProbability_cLogSecond_Ds(0x0),
		fhistJetProbability_cLogSecond_Lc(0x0),
		fhistJetProbability_bLogSecond(0x0),
		fhistJetProbabilityLogThird(0x0),
		fhistJetProbability_UnidentifiedLogThird(0x0),
		fhistJetProbability_udsgLogThird(0x0),
		fhistJetProbability_cLogThird(0x0),
		fhistJetProbability_cLogThird_D0(0x0),
		fhistJetProbability_cLogThird_Dp(0x0),
		fhistJetProbability_cLogThird_Ds(0x0),
		fhistJetProbability_cLogThird_Lc(0x0),
		fhistJetProbability_bLogThird(0x0),
		fhistJetProbabilityLogSVHE(0x0),
		fhistJetProbability_UnidentifiedLogSVHE(0x0),
		fhistJetProbability_udsgLogSVHE(0x0),
		fhistJetProbability_cLogSVHE(0x0),
		fhistJetProbability_bLogSVHE(0x0),
		fhistJetProbabilityLogSVHP(0x0),
		fhistJetProbability_UnidentifiedLogSVHP(0x0),
		fhistJetProbability_udsgLogSVHP(0x0),
		fhistJetProbability_cLogSVHP(0x0),
		fhistJetProbability_bLogSVHP(0x0),
		fMinTrackProb(0.0),
		//__________V0 Reconstruction
		fh1V0CounterCentK0s(0x0),
		fh1V0CounterCentLambda(0x0),
		fh1V0CounterCentALambda(0x0),
		fbTPCRefit(0),
		fbRejectKinks(0),
		fbFindableClusters(0),
		fdCutNCrossedRowsTPCMin(-1),
		fdCutCrossedRowsOverFindMin(-1),
		fdCutCrossedRowsOverFindMax(-1),
		fdCutPtDaughterMin(-1),
		fdCutDCAToPrimVtxMin(-1),
		fdCutDCADaughtersMax(-1),
		fdCutEtaDaughterMax(-1),
		fdCutNSigmadEdxMax(-1),
		fdPtProtonPIDMax(-1),
		fbOnFly(0),
		fdCutCPAKMin(-1),
		fdCutCPALMin(-1),
		fdCutRadiusDecayMin(-1),
		fdCutRadiusDecayMax(-1),
		fdCutEtaV0Max(-1),
		fdCutRapV0Max(-1),
		fdCutNTauKMax(-1),
		fdCutNTauLMax(-1),
		fbCutArmPod(0),
		fbCutCross(0),
		fApplyV0Rec(kFALSE),
		fApplyV0RejectionAll(kFALSE),
		//_________
		fMCArray(0x0),
		fUseCorrPt(kTRUE),
		fUsePicoTracks(kTRUE),
		fEnableV0GammaRejection(0),
		fV0CandidateArray(0x0),
		fUtils(new AliAnalysisUtils()),
		fJetContainerMC(0x0),
		fJetContainerData(0x0),
		fAODIn(0x0),
		fPrimaryVertex(0x0),
		//SV Analysis
		fDoSVAnalysis(kFALSE),
		fDoTrackCountingAnalysis(kTRUE),
		fVtxTagger3Prong(0x0),
		fVtxTagger2Prong(0x0),
		fjetCuts3Prong(0x0),
		fjetCuts2Prong(0x0),
		fTrackArray(0x0),
		fEsdTrackCuts(0x0),
		fInvariantMass(0.),
		fJetMass(0.),
		fDispersion(0.),
		fDecayLength(0.),
		fLxySign(0.),
		fJetPt(0.),
		fJetFlavor(0),
		fValJetProb(-1),
		fLogJetProb(-1.),
		fCalcDCATruth(kFALSE),
		fDecayVertex(0x0),
		fHistSV2Prong(0x0),
		fHistSV2ProngUnidentified(0x0),
		fHistSV2Prongb(0x0),
		fHistSV2Prongc(0x0),
		fHistSV2Pronglf(0x0),
		fHistDispersion2Prong(0x0),
		fHistDispersion2ProngUnidentified(0x0),
		fHistDispersion2Prongb(0x0),
		fHistDispersion2Prongc(0x0),
		fHistDispersion2Pronglf(0x0),
		fHistSV3Prong(0x0),
		fHistSV3ProngUnidentified(0x0),
		fHistSV3Prongb(0x0),
		fHistSV3Prongc(0x0),
		fHistSV3Pronglf(0x0),
		fHistDispersion3Prong(0x0),
		fHistDispersion3ProngUnidentified(0x0),
		fHistDispersion3Prongb(0x0),
		fHistDispersion3Prongc(0x0),
		fHistDispersion3Pronglf(0x0)
{
	SetMakeGeneralHistograms(kTRUE);

	for(int i=0; i<7; i++){
		fResolutionFunction[i]=0x0;
		fResolutionFunctionb[i]=0x0;
		fResolutionFunctionc[i]=0x0;
		fResolutionFunctionlf[i]=0x0;
	}
}
//#######################################
AliAnalysisTaskBJetTC::~AliAnalysisTaskBJetTC()
{
	//Destructor
	delete fCaloClusters;
	delete fReaderGammas;
	delete fV0Reader;
	delete fV0CandidateArray;
	delete fMCArray;
	delete fOutput;
	delete fJetContainerMC;
	delete fJetContainerData;
	delete fAODIn;
	delete fPrimaryVertex;
	delete fDecayVertex;
	delete fVtxTagger3Prong;
	delete fVtxTagger2Prong;
	delete fjetCuts3Prong;
	delete fjetCuts2Prong;
  	delete fTrackArray;
	delete fEsdTrackCuts;
	delete fHFJetUtils;
	delete fRespoPID;
	delete fUtils;
	delete fRandom;
	delete fVertexer;
	delete fDiamond;
}
// #################################################################################
Bool_t AliAnalysisTaskBJetTC::Notify()
{

  if(!fIsPythia) return kTRUE;

  AliAnalysisTaskEmcal::UserNotify();
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  /*TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;
  Int_t   pthbin   = 0;

  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    
    if(PythiaInfoFromFile(curfile->GetName(),xsection,ftrials,pthbin) ){
	fPythiaEventWeight = xsection/ftrials;
	cout<<"This is the weighting : "<<fPythiaEventWeight<<endl;
    }
  }*/	

  fPythiaEventWeight=1;
  return kTRUE;
}
// ########################################################################################  Main Loop
Bool_t AliAnalysisTaskBJetTC::Run()
{

	fMCArray = NULL;


	if(fIsPythia){
  		fJetContainerMC = static_cast<AliJetContainer*>(fJetCollArray.At(1));
		fMCArray= dynamic_cast<TClonesArray*>(fAODIn->FindListObject(AliAODMCParticle::StdBranchName()));
	}


	fVertexer = new AliVertexerTracks(fAODIn->GetMagneticField());
	fVertexer->SetITSMode();
	fVertexer->SetMinClusters(3);
	fVertexer->SetConstraintOn();

	if(fVertexConstraint) {
		Float_t diamondcovxy[3];
		fAODIn->GetDiamondCovXY(diamondcovxy);
		Double_t pos[3]={fAODIn->GetDiamondX(),fAODIn->GetDiamondY(),0.};
		Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
		fDiamond = new AliESDVertex(pos,cov,1.,1);
		fVertexer->SetVtxStart(fDiamond);
	}

  	Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex

	if(fApplyV0RejectionAll){

		AliAODMCHeader* headerMC = 0; // MC header

		headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
		if(!headerMC)
		{
	     	  AliError("No MC header found!");
	      	  return kFALSE;
	    	}
	       // get position of the MC primary vertex
	       dPrimVtxMCX = headerMC->GetVtxX();
	       dPrimVtxMCY = headerMC->GetVtxY();
	       dPrimVtxMCZ = headerMC->GetVtxZ();

  		AliAODMCParticle *pAOD = 0;
		AliEmcalJet * jetMC  = 0x0;
		double fJetPt=0;

	  	for (Int_t i=0; i<fMCArray->GetEntriesFast(); i++) {

	    		pAOD = dynamic_cast<AliAODMCParticle*>(fMCArray->At(i));
		        if (!pAOD) continue;
	
	    		/*Bool_t bPri = kFALSE;
	    		if (pAOD) bPri = pAOD->IsPrimary();


	    		Bool_t bPhy = kFALSE;
	    		if (pAOD) bPhy =   pAOD->IsPhysicalPrimary();

	    		if ((!bPri) && (!bPhy)) { pAOD=0; continue; }*/

			// Get the distance between the production point of the MC V0 particle and the primary vertex
		        Double_t dx = dPrimVtxMCX - pAOD->Xv();
		        Double_t dy = dPrimVtxMCY - pAOD->Yv();
		        Double_t dz = dPrimVtxMCZ - pAOD->Zv();
		        Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
		        Bool_t bV0MCIsPrimaryDist = (dDistPrimary < 0.01); // Is close enough to be considered primary-like?

		        // Select only primary-like MC V0 particles
		        if(!bV0MCIsPrimaryDist)
			  continue;

	    		Int_t id = 0;
	    		if (pAOD) id = pAOD->GetPdgCode();

	    		Bool_t bV0 = ((id==3122) || (id==-3122) || (id==310));
	    		if (!bV0) { pAOD=0; continue; }

	    		if (TMath::Abs(pAOD->Eta()) > 0.8) { pAOD=0; continue; }
		
			if(pAOD->Pt() < 0.15) continue;

			 if(id==310) {fh1dKshortPtMC->Fill(pAOD->Pt(), fPythiaEventWeight);}
			 if(id==3122) { fh1dLamdaPtMC->Fill(pAOD->Pt(), fPythiaEventWeight);}
			 if(id==-3122) { fh1dAnLamdaPtMC->Fill(pAOD->Pt(), fPythiaEventWeight);}
			 
			  fJetContainerMC->ResetCurrentID();

			  while ((jetMC = fJetContainerMC->GetNextAcceptJet()))
			  {

				fJetPt= jetMC->Pt();

				if(!(fJetContainerMC->GetRhoParameter() == 0x0))
				{
					fJetPt = fJetPt - fJetContainerMC->GetRhoVal() * jetMC->Area();
				}

				//if(!(fJetCutsHF->IsJetSelected(jetMC))) continue;

				if(fJetPt < 5.) continue;

				if(IsParticleInCone(pAOD, jetMC, 0.4)) // If good jet in event, find out whether V0 is in that jet
				{
				    if(id==310) {fh2dKshortPtVsJetPtMC->Fill(pAOD->Pt(), fJetPt, fPythiaEventWeight);}
				    if(id==3122) { fh2dLamdaPtVsJetPtMC->Fill(pAOD->Pt(), fJetPt, fPythiaEventWeight);}
				    if(id==-3122) { fh2dAnLamdaPtVsJetPtMC->Fill(pAOD->Pt(), fJetPt, fPythiaEventWeight);}
				    break;
				}
				

			 }


	    	}
		jetMC=NULL;
		delete jetMC;
		pAOD=NULL;
		delete pAOD;
		headerMC=NULL;
		delete headerMC;
	}

	if(fEnableV0GammaRejection){

	  	fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

		AliAODConversionPhoton* PhotonCandidate = 0x0;
		// Loop over Photon Candidates allocated by ReaderV1
		for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){

		    PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
		    if(!PhotonCandidate) continue;

		    TVector3 vV0;
	    	    if (PhotonCandidate) vV0.SetXYZ(PhotonCandidate->GetPx(), PhotonCandidate->GetPy(), PhotonCandidate->GetPz());

		    if(IsV0InJet(vV0,5.)) {fh1dPhotonPt->Fill(PhotonCandidate->GetPhotonPt(), fPythiaEventWeight);}

		}
		if(fIsPythia && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
		    fV0Reader->RelabelAODs(kTRUE);
		}
		PhotonCandidate=NULL;
		delete PhotonCandidate;
	}

	//Main loop over AOD tracks with filterbit 4  || or ESD tracks with filter
	int nTracksInEvent = 0 ;

	nTracksInEvent  = fAODIn ->GetNumberOfTracks() ;

	AliAODTrack * trackAOD =  NULL;

	if(fCalcDCATruth){
		for(int itrack= 0; itrack<nTracksInEvent;++itrack)
		{
		
			trackAOD = (AliAODTrack*)fAODIn->GetTrack(itrack);
			if (!trackAOD) continue;

			fh1dTracksAccepeted->SetBinContent(1,fh1dTracksAccepeted->GetBinContent(1)+1);

			if(!IsTrackAccepted(trackAOD)) {
				fh1dTracksAccepeted->SetBinContent(3,fh1dTracksAccepeted->GetBinContent(3)+1);
				continue;
			}

			fh1dTracksAccepeted->SetBinContent(2,fh1dTracksAccepeted->GetBinContent(2)+1);
			//Calculate impact parameters and fill histograms
			double dca[2] = {-99999,-99999};
			double cov[3] = {-99999,-99999,-99999};

			if (!CalculateTrackImpactParameter(trackAOD,dca,cov)) continue;

			fh1dTracksImpParXY->Fill(GetValImpactParameter(kXY,dca,cov),fPythiaEventWeight);
			fh1dTracksImpParXYZ->Fill(GetValImpactParameter(kXYZ,dca,cov),fPythiaEventWeight);
			fh1dTracksImpParXYSignificance->Fill(GetValImpactParameter(kXYSig,dca,cov),fPythiaEventWeight);
			fh1dTracksImpParXYZSignificance->Fill(GetValImpactParameter(kXYZSig,dca,cov),fPythiaEventWeight);

			if(fIsPythia){
			
				double dcaMC[2] = {-99999,-99999};
				double covMC[3] = {-99999,-99999,-99999};

				if(!CalculateTrackImpactParameterTruth(trackAOD,dcaMC,covMC)) continue;

				fh1dTracksImpParXYTruth->Fill(GetValImpactParameter(kXY,dcaMC,covMC),fPythiaEventWeight);
				fh1dTracksImpParXYZTruth->Fill(GetValImpactParameter(kXYZ,dcaMC,covMC),fPythiaEventWeight);
				// Fill residual plots
				double residualxy = TMath::Abs(GetValImpactParameter(kXY,dca,cov)) - TMath::Abs(GetValImpactParameter(kXY,dcaMC,covMC));
				residualxy /= TMath::Sqrt(cov[0]);
				fh1dTracksImpParXYResidualTruth->Fill(residualxy,fPythiaEventWeight);
				double residualxyz = TMath::Abs(GetValImpactParameter(kXYZ,dca,cov)) - TMath::Abs(GetValImpactParameter(kXYZ,dcaMC,covMC));
				residualxyz /= 	GetValImpactParameter(kXYZSigmaOnly,dca,cov);
				fh1dTracksImpParXYZResidualTruth->Fill(residualxyz,fPythiaEventWeight);
			}
		}
		trackAOD = NULL;
		delete trackAOD;
	}

	// Main part jet analysis
	//preparation
	fJetContainerData = static_cast<AliJetContainer*>(fJetCollArray.At(0));

	TString JetContName = fJetContainerData->GetName();
	if(JetContName.Contains("PicoTracks")) fUsePicoTracks = kTRUE;
	else fUsePicoTracks = kFALSE;

	Double_t randomConePt = GetDeltaPtRandomCone();
	Double_t randomConePtWithSignal = 0.0;
	f2histRhoVsDeltaPt->Fill(randomConePt, fJetContainerData->GetRhoVal(), fPythiaEventWeight);

	if(fDoDeltaPtWithSignal){
		randomConePtWithSignal = GetDeltaPtRandomConeWithSignal();
		f2histRhoVsDeltaPtWithSignal->Fill(randomConePtWithSignal, fJetContainerData->GetRhoVal(), fPythiaEventWeight);
	}

	if(fApplyV0Rec) SelectV0CandidateVIT();


	if(fIsPythia)
	{
		// SetContainer
		AliEmcalJet * jetgen  = 0x0;
		AliAODMCParticle* partonAOD = NULL;

		if(!MatchJetsGeometricDefault()) cout << "Error running jet matching!" << endl;
		fJetContainerMC->ResetCurrentID();

		// Fill gen. level jet histograms
		while ((jetgen = fJetContainerMC->GetNextAcceptJet()))
		{
			if (!jetgen) continue;
			Int_t MCJetflavour =0;
			Int_t partonpdg=0;
			
			if(fUsePartonDef){
			  partonAOD = fHFJetUtils->IsMCJetParton(fMCArray, jetgen, 0.4);
			  if(!(partonAOD)) MCJetflavour =0;
			  else
			  {
				partonpdg = abs(partonAOD->PdgCode());
				if(partonpdg==1||partonpdg==2||partonpdg==3||partonpdg==21 )MCJetflavour=1;
				else if(partonpdg==4)MCJetflavour=2;
				else if(partonpdg==5)MCJetflavour=3;
			  }
			}else {
			  partonAOD = fHFJetUtils->IsMCJetMeson(fMCArray, jetgen, 0.4);
			  if(!(partonAOD)) MCJetflavour =0;
			  else
			  {
				partonpdg = abs(partonAOD->PdgCode());
				if(fHFJetUtils->IsBMeson(partonpdg)) MCJetflavour=3;
				else if(fHFJetUtils->IsDMeson(partonpdg)) MCJetflavour=2;
				else MCJetflavour=1;
			  }
			}

			double genpt = jetgen->Pt();
			if(!(fJetContainerMC->GetRhoParameter() == 0x0)){
				 genpt = genpt - fJetContainerMC->GetRhoVal() * jetgen->Area();
			}

			fh1dJetGenPt->Fill(genpt,fPythiaEventWeight);
			if(MCJetflavour ==0)
				fh1dJetGenPtUnidentified->Fill(genpt,fPythiaEventWeight);
			else if(MCJetflavour ==1)
				fh1dJetGenPtudsg->Fill(genpt,fPythiaEventWeight);
			else if(MCJetflavour ==2)
				fh1dJetGenPtc->Fill(genpt,fPythiaEventWeight);
			else if(MCJetflavour ==3)
				fh1dJetGenPtb->Fill(genpt,fPythiaEventWeight);
		}
		jetgen = 0x0;
		delete jetgen;

		partonAOD=NULL;
		delete partonAOD;
	}

	// loop rec level jets
	AliEmcalJet * jetrec  = 0x0;
	AliEmcalJet * jetmatched  = 0x0;
	fJetContainerData->ResetCurrentID();
	fJetPt=0;
	double jetptmc=0;

	//########################## Electron Enriched Sample
	Bool_t PtRelSample = kTRUE; 

	if(fDoPtRelAnalysis) fCaloClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));

	if(fDoPtRelAnalysis && fResolutionFunction[0] && fDoSelectionPtRel){

		Bool_t ElecJet(0), TagJet(0);
		double IPxy[2] = {-99999,-99999};
		double Cov[3] = {-99999,-99999,-99999};
		double IPsign=0;
		double DCA_Track_Jet(0), bDecayLength(0);
		AliEmcalJet* selJet = 0x0;

		PtRelSample = kFALSE;

		AliAODTrack* trackAOD = 0x0;

		while ((selJet = fJetContainerData->GetNextAcceptJet()))
		{
			//if(!(fJetCutsHF->IsJetSelected(selJet))) continue;

			Int_t ntracks = (Int_t)selJet->GetNumberOfTracks();

			for(Int_t itrack = 0; itrack < ntracks; ++itrack)
			{
				if(fUsePicoTracks) trackAOD = (AliAODTrack*)((AliPicoTrack*)selJet->Track(itrack))->GetTrack();
				else trackAOD = (AliAODTrack*)((fJetContainerData->GetParticleContainer())->GetParticle(selJet->TrackAt(itrack)));
				
				if(!trackAOD) 	continue;

				if(IsElectronHF(trackAOD)) {
					ElecJet=kTRUE;
					break;
				}
			}

              		Double_t value = CalculateJetProb(selJet, 0);
			if(value==0) value = 1e-5;
			Double_t LogValue = -1*TMath::Log(value);

			if(LogValue>=5) TagJet=kTRUE;

			PtRelSample = (TagJet && ElecJet);
			if(PtRelSample) break;
		}
		if(PtRelSample) fhistPtRelEvents->Fill(0.5);

		selJet = 0x0;
		trackAOD = NULL;
		delete selJet; delete trackAOD;

		fJetContainerData->ResetCurrentID();	
	}
	//######################

	AliAODMCParticle* partonAOD = NULL;

	std::vector<double> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;

	Bool_t TaggedFirst(0), TaggedSecond(0), TaggedThird(0);

	while ((jetrec = fJetContainerData->GetNextJet()))
	{
		//	Printf("%s:%i",__FUNCTION__,__LINE__);
		fJetPt= jetrec->Pt();
		//	Printf("%s:%i",__FUNCTION__,__LINE__);

		if(!(fJetContainerData->GetRhoParameter() == 0x0) && fUseCorrPt)
		{
			fJetPt = fJetPt - fJetContainerData->GetRhoVal() * jetrec->Area();
		}

		// make inclusive signed imp. parameter constituent histograms
		Int_t ntracks = (Int_t)jetrec->GetNumberOfTracks();

		double dca[2] = {-99999,-99999};
		double cov[3] = {-99999,-99999,-99999};
		double sign=0;
		Double_t ElePtRel[10]={0};
		Double_t EleIP[10]={0};
		Double_t EleID[10]={0};
		Int_t ElecNum=0;

		//Printf("%s:%i",__FUNCTION__,__LINE__);

		fJetFlavor =0;
		Int_t partonpdg=0;
		if(fIsPythia){
			jetmatched = 0x0;
			jetmatched =jetrec->MatchedJet();

			if(jetmatched){
			  if(fUsePartonDef){
				partonAOD = fHFJetUtils->IsMCJetParton(fMCArray, jetmatched, 0.4);

				if((!partonAOD)) fJetFlavor =0;
				else{
					partonpdg = abs(partonAOD->PdgCode());
					if(partonpdg==1||partonpdg==2||partonpdg==3||partonpdg==21 )fJetFlavor=1;
					else if(partonpdg==4)fJetFlavor=2;
					else if(partonpdg==5)fJetFlavor=3;
				}
			  }else {
				partonAOD = fHFJetUtils->IsMCJetMeson(fMCArray, jetmatched, 0.4);
				if(!(partonAOD)) fJetFlavor =0;
				else
				{
					partonpdg = abs(partonAOD->PdgCode());
					if(fHFJetUtils->IsBMeson(partonpdg)) fJetFlavor=3;
					else if(fHFJetUtils->IsDMeson(partonpdg)) fJetFlavor=2;
					else fJetFlavor=1;
				}
			}
			}
		}
		//	Printf("%s:%i",__FUNCTION__,__LINE__);

			fh1dJetRecPt->Fill(jetrec->Pt(),fPythiaEventWeight);
			if(fIsPythia){			
				  if(fJetFlavor==0) fh1dJetRecPtUnidentified->Fill(fJetPt,fPythiaEventWeight);
				  else if(fJetFlavor==1)fh1dJetRecPtudsg->Fill(fJetPt,fPythiaEventWeight);
				  else if(fJetFlavor==2)fh1dJetRecPtc->Fill(fJetPt,fPythiaEventWeight);
				  else if(fJetFlavor==3)fh1dJetRecPtb->Fill(fJetPt,fPythiaEventWeight);
				
			}
			
			UInt_t rejectionReason = 0;
			if(!(fJetContainerData->AcceptJet(jetrec, rejectionReason))) continue;

			fJetMass = jetrec->M();

			//if(fJetPt < 5.0) continue;

			fh1dJetRecEtaPhiAccepted->Fill(jetrec->Eta(),jetrec->Phi(),fPythiaEventWeight);
			fh1dJetRecPtAcceptedunCorr->Fill(jetrec->Pt(),fPythiaEventWeight);
			fh1dJetRecPtAccepted->Fill(fJetPt,fPythiaEventWeight);
			if(fIsPythia){
				if (jetrec->MatchedJet()) {
				  double genpt = jetrec->MatchedJet()->Pt();
				  if(!(fJetContainerMC->GetRhoParameter() == 0x0)){
					genpt = genpt - fJetContainerMC->GetRhoVal() * jetrec->MatchedJet()->Area();
				  }
				  fh2dJetGenPtVsJetRecPt->Fill(fJetPt,genpt,fPythiaEventWeight);

				  if(fJetFlavor==0){ fh1dJetRecPtUnidentifiedAccepted->Fill(fJetPt,fPythiaEventWeight);}
				  else if(fJetFlavor==1){fh1dJetRecPtudsgAccepted->Fill(fJetPt,fPythiaEventWeight); fh2dJetGenPtVsJetRecPtudsg->Fill(fJetPt,genpt,fPythiaEventWeight);}
				  else if(fJetFlavor==2){fh1dJetRecPtcAccepted->Fill(fJetPt,fPythiaEventWeight); fh2dJetGenPtVsJetRecPtc->Fill(fJetPt,genpt,fPythiaEventWeight);}
				  else if(fJetFlavor==3){fh1dJetRecPtbAccepted->Fill(fJetPt,fPythiaEventWeight); fh2dJetGenPtVsJetRecPtb->Fill(fJetPt,genpt,fPythiaEventWeight);}
				}
			}
			

			Double_t Mass = jetrec->M();

			if(fDoJetMass){

			      	fhistJetMass->Fill(fJetPt,Mass ,fPythiaEventWeight);
			    		
			  	if(fIsPythia){
			      	    switch(fJetFlavor)
					{
					case 0:
						fhistJetMass_Unidentified->Fill(fJetPt,jetrec->M(),fPythiaEventWeight);
						break;
					case 1:
						fhistJetMass_udsg->Fill(fJetPt,Mass,fPythiaEventWeight);
						break;
					case 2:
						fhistJetMass_c->Fill(fJetPt,Mass,fPythiaEventWeight);
						break;
					case 3:
						fhistJetMass_b->Fill(fJetPt,Mass,fPythiaEventWeight);
						break;
					default:
						break;
					}
			    	}
			}


			fValJetProb = -1.;
			fLogJetProb = -1.;


			if(fDoJetProbabilityAnalysis && fResolutionFunction[0]){

              			fValJetProb = CalculateJetProb(jetrec, fJetFlavor);
				if(fValJetProb>0){
					fLogJetProb = -1*TMath::Log(fValJetProb);

		      			fhistJetProbability->Fill(fJetPt,fValJetProb,fPythiaEventWeight);
		      			fhistJetProbabilityLog->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
		    			
		  			if(fIsPythia){
		      			    switch(fJetFlavor)
		        			{
		        			case 0:
		          				fhistJetProbability_Unidentified->Fill(fJetPt,fValJetProb,fPythiaEventWeight);
		          				fhistJetProbability_UnidentifiedLog->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
		          				break;
		        			case 1:
		          				fhistJetProbability_udsg->Fill(fJetPt,fValJetProb,fPythiaEventWeight);
		          				fhistJetProbability_udsgLog->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
		          				break;
		        			case 2:
		          				fhistJetProbability_c->Fill(fJetPt,fValJetProb,fPythiaEventWeight);
		          				fhistJetProbability_cLog->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
							if(fDoCharmFractions){
							  switch(TMath::Abs(partonpdg)){
							    case 421:
							      fhistJetProbability_cLog_D0->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
							      break;
							    case 411:
							      fhistJetProbability_cLog_Dp->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
							      break;
							    case 431:
							      fhistJetProbability_cLog_Ds->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
							      break;
							    case 4122:
							      fhistJetProbability_cLog_Lc->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
							      break;
							    default:
							      break;
							  }
							}
		          				break;
		        			case 3:
				  			fhistJetProbability_b->Fill(fJetPt,fValJetProb,fPythiaEventWeight);
				  			fhistJetProbability_bLog->Fill(fJetPt,fLogJetProb,fPythiaEventWeight);
				  			break;
		        			default:
		          				break;
		        			}
		    			}
				}
        		}


			Double_t SVEnergy=0.;
			Double_t TracksEnergy=0.;
			Double_t EnergyFraction=-1;

			if(fDoSVEnergyFraction){

			  Double_t vtxPos[3]   = {fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ()};
			  Double_t covMatrix[6] = {0};
			  fPrimaryVertex->GetCovarianceMatrix(covMatrix);
			  AliESDVertex* esdVtx = new AliESDVertex(vtxPos, covMatrix, fPrimaryVertex->GetChi2(), fPrimaryVertex->GetNContributors());

			  // 3 Prong Vertex
			  TClonesArray* secVertexArrProng = 0;
			  vector<pair <Double_t, Int_t>> arrDispersionProng;
			  arrDispersionProng.reserve(5);

			  secVertexArrProng = new TClonesArray("AliAODVertex");
			  Int_t nDauRejCountProng = 0;
			  Int_t nVtx6Prong = FindVertices6Prong(jetrec,
						                 static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
						                 fAODIn,
						                 esdVtx,
						                 fAODIn->GetMagneticField(),
						                 secVertexArrProng,
						                 nDauRejCountProng);

			  if(nVtx6Prong > 0)
			  {

				  Int_t MaxSVindex=-1;
				  Double_t MaxSVLxy=0.;
				  Double_t MaxSVLxyS=0.;

				  for(Int_t iv=0; iv<secVertexArrProng->GetEntriesFast(); iv++)
				  {
				    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(iv));

				    // Calculate vtx distance
				    Double_t effX = secVtx->GetX() - esdVtx->GetX();
				    Double_t effY = secVtx->GetY() - esdVtx->GetY();

				    // signed length
				    Double_t decLenXY  = fPrimaryVertex->DistanceXYToVertex(secVtx);
				    Double_t jetP[3]; jetrec->PxPyPz(jetP);
				    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
				    if (signLxy < 0.) decLenXY *= -1.;

				    Double_t errdecLenXY = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

				    if(decLenXY > MaxSVLxy){
					MaxSVLxy = decLenXY;
					MaxSVindex = iv;
					MaxSVLxyS = decLenXY / errdecLenXY;
				     }
   
				  }

				  if(MaxSVindex>=0){
				    	AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(MaxSVindex));

					for(int iDaugh=0; iDaugh<secVtx->GetNDaughters(); iDaugh++){
						AliAODTrack* Daught = (AliAODTrack*)secVtx->GetDaughter(iDaugh);
						SVEnergy+= Daught->E();
					}

					fhistSVnProngs->Fill(fJetPt, 6);

		  			if(fIsPythia){
		      			    switch(fJetFlavor)
		        			{
		        			case 0:
		          				fhistSVnProngs_Unidentified->Fill(fJetPt,6);
		          				break;
		        			case 1:
		          				fhistSVnProngs_udsg->Fill(fJetPt,6);
		          				break;
		        			case 2:
		          				fhistSVnProngs_c->Fill(fJetPt,6);
		          				break;
		        			case 3:
				  			fhistSVnProngs_b->Fill(fJetPt,6);
				  			break;
		        			default:
		          				break;
		        			}
		    			}
				  }

				    secVertexArrProng->Clear();
			  	    delete secVertexArrProng;
			  }else{
				    secVertexArrProng->Clear();
				    delete secVertexArrProng;

				  secVertexArrProng = new TClonesArray("AliAODVertex");
				  Int_t nDauRejCountProng = 0;
				  Int_t nVtx5Prong = FindVertices5Prong(jetrec,
								         static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
								         fAODIn,
								         esdVtx,
								         fAODIn->GetMagneticField(),
								         secVertexArrProng,
								         nDauRejCountProng);

				  if(nVtx5Prong > 0)
				  {

					  Int_t MaxSVindex=-1;
					  Double_t MaxSVLxy=0.;
					  Double_t MaxSVLxyS=0.;

					  for(Int_t iv=0; iv<secVertexArrProng->GetEntriesFast(); iv++)
					  {
					    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(iv));

					    // Calculate vtx distance
					    Double_t effX = secVtx->GetX() - esdVtx->GetX();
					    Double_t effY = secVtx->GetY() - esdVtx->GetY();

					    // signed length
					    Double_t decLenXY  = fPrimaryVertex->DistanceXYToVertex(secVtx);
					    Double_t jetP[3]; jetrec->PxPyPz(jetP);
					    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
					    if (signLxy < 0.) decLenXY *= -1.;

					    Double_t errdecLenXY = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

					    if(decLenXY > MaxSVLxy){
						MaxSVLxy = decLenXY;
						MaxSVindex = iv;
						MaxSVLxyS = decLenXY / errdecLenXY;
					     }
	   
					  }

					  if(MaxSVindex>=0){
					    	AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(MaxSVindex));

						for(int iDaugh=0; iDaugh<secVtx->GetNDaughters(); iDaugh++){
							AliAODTrack* Daught = (AliAODTrack*)secVtx->GetDaughter(iDaugh);
							SVEnergy+= Daught->E();
						}

						fhistSVnProngs->Fill(fJetPt, 5);

			  			if(fIsPythia){
			      			    switch(fJetFlavor)
							{
							case 0:
				  				fhistSVnProngs_Unidentified->Fill(fJetPt,5);
				  				break;
							case 1:
				  				fhistSVnProngs_udsg->Fill(fJetPt,5);
				  				break;
							case 2:
				  				fhistSVnProngs_c->Fill(fJetPt,5);
				  				break;
							case 3:
					  			fhistSVnProngs_b->Fill(fJetPt,5);
					  			break;
							default:
				  				break;
							}
			    			}
					  }

					    secVertexArrProng->Clear();
				  	    delete secVertexArrProng;
				  }else{
					    secVertexArrProng->Clear();
					    delete secVertexArrProng;

					    secVertexArrProng = new TClonesArray("AliAODVertex");
					    Int_t nDauRejCountProng = 0;
					    Int_t nVtx4Prong = FindVertices4Prong(jetrec,
											 static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
											 fAODIn,
											 esdVtx,
											 fAODIn->GetMagneticField(),
											 secVertexArrProng,
											 nDauRejCountProng);

					  if(nVtx4Prong > 0)
					  {

						  Int_t MaxSVindex=-1;
						  Double_t MaxSVLxy=0.;
						  Double_t MaxSVLxyS=0.;

						  for(Int_t iv=0; iv<secVertexArrProng->GetEntriesFast(); iv++)
						  {
						    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(iv));

						    // Calculate vtx distance
						    Double_t effX = secVtx->GetX() - esdVtx->GetX();
						    Double_t effY = secVtx->GetY() - esdVtx->GetY();

						    // signed length
						    Double_t decLenXY  = fPrimaryVertex->DistanceXYToVertex(secVtx);
						    Double_t jetP[3]; jetrec->PxPyPz(jetP);
						    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
						    if (signLxy < 0.) decLenXY *= -1.;

						    Double_t errdecLenXY = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

						    if(decLenXY > MaxSVLxy){
							MaxSVLxy = decLenXY;
							MaxSVindex = iv;
							MaxSVLxyS = decLenXY / errdecLenXY;
						     }
		   
						  }

						  if(MaxSVindex>=0){
						    	AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(MaxSVindex));

							for(int iDaugh=0; iDaugh<secVtx->GetNDaughters(); iDaugh++){
								AliAODTrack* Daught = (AliAODTrack*)secVtx->GetDaughter(iDaugh);
								SVEnergy+= Daught->E();
							}
							fhistSVnProngs->Fill(fJetPt, 4);

				  			if(fIsPythia){
				      			    switch(fJetFlavor)
								{
								case 0:
					  				fhistSVnProngs_Unidentified->Fill(fJetPt,4);
					  				break;
								case 1:
					  				fhistSVnProngs_udsg->Fill(fJetPt,4);
					  				break;
								case 2:
					  				fhistSVnProngs_c->Fill(fJetPt,4);
					  				break;
								case 3:
						  			fhistSVnProngs_b->Fill(fJetPt,4);
						  			break;
								default:
					  				break;
								}
				    			}
						  }

						    secVertexArrProng->Clear();
					  	    delete secVertexArrProng;
					  }else{
						    secVertexArrProng->Clear();
					  	    delete secVertexArrProng;

						    secVertexArrProng = new TClonesArray("AliAODVertex");
						    Int_t nDauRejCountProng = 0;

						    Int_t nVtx3Prong = fVtxTagger3Prong->FindVertices(jetrec,
											 static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
											 fAODIn,
											 esdVtx,
											 fAODIn->GetMagneticField(),
											 secVertexArrProng,
											 0,
											 arrDispersionProng,
											 nDauRejCountProng);

						  if(nVtx3Prong > 0)
						  {

							  Int_t MaxSVindex=-1;
							  Double_t MaxSVLxy=0.;
							  Double_t MaxSVLxyS=0.;

							  for(Int_t iv=0; iv<secVertexArrProng->GetEntriesFast(); iv++)
							  {
							    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(iv));

							    // Calculate vtx distance
							    Double_t effX = secVtx->GetX() - esdVtx->GetX();
							    Double_t effY = secVtx->GetY() - esdVtx->GetY();

							    // signed length
							    Double_t decLenXY  = fPrimaryVertex->DistanceXYToVertex(secVtx);
							    Double_t jetP[3]; jetrec->PxPyPz(jetP);
							    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
							    if (signLxy < 0.) decLenXY *= -1.;

							    Double_t errdecLenXY = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

							    if(decLenXY > MaxSVLxy){
								MaxSVLxy = decLenXY;
								MaxSVindex = iv;
								MaxSVLxyS = decLenXY / errdecLenXY;
							     }
			   
							  }

							  if(MaxSVindex>=0){
							    	AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(MaxSVindex));

								for(int iDaugh=0; iDaugh<secVtx->GetNDaughters(); iDaugh++){
									AliAODTrack* Daught = (AliAODTrack*)secVtx->GetDaughter(iDaugh);
									SVEnergy+= Daught->E();
								}

								fhistSVnProngs->Fill(fJetPt, 3);

					  			if(fIsPythia){
					      			    switch(fJetFlavor)
									{
									case 0:
						  				fhistSVnProngs_Unidentified->Fill(fJetPt,3);
						  				break;
									case 1:
						  				fhistSVnProngs_udsg->Fill(fJetPt,3);
						  				break;
									case 2:
						  				fhistSVnProngs_c->Fill(fJetPt,3);
						  				break;
									case 3:
							  			fhistSVnProngs_b->Fill(fJetPt,3);
							  			break;
									default:
						  				break;
									}
					    			}
							  }

							    secVertexArrProng->Clear();
						  	    delete secVertexArrProng;
						}else{
							    secVertexArrProng->Clear();
						  	    delete secVertexArrProng;

							    secVertexArrProng = new TClonesArray("AliAODVertex");
							    Int_t nDauRejCountProng = 0;

							    Int_t nVtx2Prong = fVtxTagger2Prong->FindVertices(jetrec,
												 static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
												 fAODIn,
												 esdVtx,
												 fAODIn->GetMagneticField(),
												 secVertexArrProng,
												 0,
												 arrDispersionProng,
												 nDauRejCountProng);

							  if(nVtx2Prong > 0)
							  {

								  Int_t MaxSVindex=-1;
								  Double_t MaxSVLxy=0.;
								  Double_t MaxSVLxyS=0.;

								  for(Int_t iv=0; iv<secVertexArrProng->GetEntriesFast(); iv++)
								  {
								    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(iv));

								    // Calculate vtx distance
								    Double_t effX = secVtx->GetX() - esdVtx->GetX();
								    Double_t effY = secVtx->GetY() - esdVtx->GetY();

								    // signed length
								    Double_t decLenXY  = fPrimaryVertex->DistanceXYToVertex(secVtx);
								    Double_t jetP[3]; jetrec->PxPyPz(jetP);
								    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
								    if (signLxy < 0.) decLenXY *= -1.;

								    Double_t errdecLenXY = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

								    if(decLenXY > MaxSVLxy){
									MaxSVLxy = decLenXY;
									MaxSVindex = iv;
									MaxSVLxyS = decLenXY / errdecLenXY;
								     }
				   
								  }

								  if(MaxSVindex>=0){
								    	AliAODVertex* secVtx = (AliAODVertex*)(secVertexArrProng->UncheckedAt(MaxSVindex));

									for(int iDaugh=0; iDaugh<secVtx->GetNDaughters(); iDaugh++){
										AliAODTrack* Daught = (AliAODTrack*)secVtx->GetDaughter(iDaugh);
										SVEnergy+= Daught->E();
									}

									fhistSVnProngs->Fill(fJetPt, 2);

						  			if(fIsPythia){
						      			    switch(fJetFlavor)
										{
										case 0:
							  				fhistSVnProngs_Unidentified->Fill(fJetPt,2);
							  				break;
										case 1:
							  				fhistSVnProngs_udsg->Fill(fJetPt,2);
							  				break;
										case 2:
							  				fhistSVnProngs_c->Fill(fJetPt,2);
							  				break;
										case 3:
								  			fhistSVnProngs_b->Fill(fJetPt,2);
								  			break;
										default:
							  				break;
										}
						    			}

							 	  }		

							    secVertexArrProng->Clear();
						  	    delete secVertexArrProng;
							}else{
							    secVertexArrProng->Clear();
						  	    delete secVertexArrProng;
							} // 2 prong
						} // 3 prong
					  }// 4 prong
				  }// 5 prong
			     }// 6 prong

			}//SV fE

			for(Int_t itrack = 0; itrack < ntracks; ++itrack)
			{
				double dcatrackjet =999;
				double lineardecaylenth = 999;

				if(fUsePicoTracks) trackAOD = (AliAODTrack*)((AliPicoTrack*)jetrec->Track(itrack))->GetTrack();
				else trackAOD = (AliAODTrack*)((fJetContainerData->GetParticleContainer())->GetParticle(jetrec->TrackAt(itrack)));
				
				if(!trackAOD) 	continue;

				TracksEnergy+=trackAOD->E();

				if (fDoJetProbabilityAnalysis && !fResolutionFunction[0]) 
					FillResolutionFunctionHists(trackAOD,jetrec,fJetFlavor);

				if(fDoPtRelAnalysis && PtRelSample){
					Double_t PtRel=0;
					double LepIP[2] = {-99999,-99999};
					double COV[3] = {-99999,-99999,-99999};
					Bool_t hasIP=0;
				
					if(IsElectronHF(trackAOD)){
					   PtRel = GetPtRel(trackAOD, jetrec, kFALSE);
					   if(PtRel>2) PtRel=1.99;
					   //if(CalculateTrackImpactParameter(trackAOD,LepIP,COV)) hasIP=kTRUE;
					   fhistPtRelVsJetPt->Fill(fJetPt, PtRel, fPythiaEventWeight);
					   if(hasIP) fhistLepIPVsJetPt->Fill(fJetPt, GetValImpactParameter(kXY,LepIP,COV), fPythiaEventWeight);
					   EleID[ElecNum]=itrack;
					   ElePtRel[ElecNum]=PtRel;
					   if(hasIP) EleIP[ElecNum]=GetValImpactParameter(kXY,LepIP,COV);
					   else EleIP[ElecNum]=-999;
					   if(fIsPythia){
						if(fJetFlavor==0){fhistPtRelVsJetPtUnidentified->Fill(fJetPt, PtRel,fPythiaEventWeight);}
					  	else if(fJetFlavor==1){fhistPtRelVsJetPtudsg->Fill(fJetPt, PtRel,fPythiaEventWeight);}
					  	else if(fJetFlavor==2){fhistPtRelVsJetPtc->Fill(fJetPt, PtRel,fPythiaEventWeight);}
					  	else if(fJetFlavor==3){fhistPtRelVsJetPtb->Fill(fJetPt, PtRel,fPythiaEventWeight);}
						if(hasIP){	     
						   if(fJetFlavor==0){fhistLepIPVsJetPtUnidentified->Fill(fJetPt, GetValImpactParameter(kXY,LepIP,COV),fPythiaEventWeight);}
					  	   else if(fJetFlavor==1){fhistLepIPVsJetPtudsg->Fill(fJetPt, GetValImpactParameter(kXY,LepIP,COV),fPythiaEventWeight);}
					  	   else if(fJetFlavor==2){fhistLepIPVsJetPtc->Fill(fJetPt, GetValImpactParameter(kXY,LepIP,COV),fPythiaEventWeight);}
					  	   else if(fJetFlavor==3){fhistLepIPVsJetPtb->Fill(fJetPt, GetValImpactParameter(kXY,LepIP,COV),fPythiaEventWeight);}
						}
					   }
					   ElecNum++;
					}
				}

				//if(!fDoTrackCountingAnalysis) continue;

				if(fIsPythia && fApplyV0RejectionAll){
					AliAODMCParticle* pMC = GetMCTrack(trackAOD);
					Int_t PDGcode = pMC->GetPdgCode();
					if( PDGcode==  211|| PDGcode==  -211|| PDGcode==  2212|| PDGcode==  -2212){
					Int_t pMotherLabel = 	pMC->GetMother();
					if(pMotherLabel>2){
					  AliAODMCParticle* mcpartMother = dynamic_cast<AliAODMCParticle*>(fMCArray->At(pMotherLabel));
					  if(mcpartMother){

						// Get the distance between the production point of the MC V0 particle and the primary vertex
						Double_t dx = dPrimVtxMCX - mcpartMother->Xv();
						Double_t dy = dPrimVtxMCY - mcpartMother->Yv();
						Double_t dz = dPrimVtxMCZ - mcpartMother->Zv();
						Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
						Bool_t bV0MCIsPrimaryDist = (dDistPrimary < 0.01); // Is close enough to be considered primary-like?

						// Select only primary-like MC V0 particles
						if(bV0MCIsPrimaryDist){

						    Int_t maPdgcode = mcpartMother->GetPdgCode();

						    if(TMath::Abs(maPdgcode) == 310){
						    	continue;
						    }
						    if(TMath::Abs(maPdgcode) == 3122){
							continue;
						    }
						}
					  }
					}
					}
				}

				if(fEnableV0GammaRejection) {
					if(IsV0PhotonFromBeamPipeDaughter(trackAOD))
						continue;
				}


				if(fApplyV0Rec){
					if(IsV0Daughter(trackAOD))
						continue;
				}

				if (!IsTrackAcceptedBJetCuts(trackAOD, fJetFlavor)) continue;
          			
				Double_t d[2]={0};
				Double_t covM[3]={0};

				if(!CalculateJetSignedTrackImpactParameter(trackAOD,jetrec,dca,cov,sign,dcatrackjet,lineardecaylenth)) continue;

				//Select only tracks with dca rphi <1 cm and dca z < 5 cm
				// linear decay length < 5 cm
				// dca track to jet < 0.07 cm

				if(abs(dca[0])>fTCMaxIPxy) continue; 	FillCandidateJet(10, fJetFlavor);
				if(abs(dca[1])>fTCMaxIPz) continue;	FillCandidateJet(11, fJetFlavor);
				if(lineardecaylenth > fTCMaxDecayLength) continue;	FillCandidateJet(12, fJetFlavor);
				if(dcatrackjet > fTCMaxDCATrackJet) continue;	FillCandidateJet(13, fJetFlavor);

				double cursImParXY =TMath::Abs(GetValImpactParameter(kXY,dca,cov))*sign;
				double cursImParXYZ =TMath::Abs(GetValImpactParameter(kXYZ,dca,cov))*sign;
				double cursImParXYSig =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
				double cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;

				if (cursImParXY > fThresholdIP) FillCandidateJet(14, fJetFlavor);

				fh2dJetSignedImpParXY->Fill(fJetPt,cursImParXY,fPythiaEventWeight);
				fh2dJetSignedImpParXYZ->Fill(fJetPt,cursImParXYZ,fPythiaEventWeight);
				fh2dJetSignedImpParXYSignificance->Fill(fJetPt,cursImParXYSig,fPythiaEventWeight);
				fh2dJetSignedImpParXYZSignificance->Fill(fJetPt,cursImParXYZSig,fPythiaEventWeight);

				if(!fDoTrackCountingAnalysis) continue;

				if(fIsPythia)
				{
					if(fJetFlavor ==0){
						fh2dJetSignedImpParXYUnidentified->Fill(fJetPt,cursImParXY,fPythiaEventWeight);
						fh2dJetSignedImpParXYZUnidentified->Fill(fJetPt,cursImParXYZ,fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceUnidentified->Fill(fJetPt,cursImParXYSig,fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceUnidentified->Fill(fJetPt,cursImParXYZSig,fPythiaEventWeight);
					}
					else if(fJetFlavor ==1){
						fh2dJetSignedImpParXYudsg->Fill(fJetPt,cursImParXY,fPythiaEventWeight);
						fh2dJetSignedImpParXYZudsg->Fill(fJetPt,cursImParXYZ,fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceudsg->Fill(fJetPt,cursImParXYSig,fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceudsg->Fill(fJetPt,cursImParXYZSig,fPythiaEventWeight);
					}
					else if(fJetFlavor ==2){
						fh2dJetSignedImpParXYc->Fill(fJetPt,cursImParXY,fPythiaEventWeight);
						fh2dJetSignedImpParXYZc->Fill(fJetPt,cursImParXYZ,fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancec->Fill(fJetPt,cursImParXYSig,fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancec->Fill(fJetPt,cursImParXYZSig,fPythiaEventWeight);
					}
					else if(fJetFlavor ==3){
						fh2dJetSignedImpParXYb->Fill(fJetPt,cursImParXY,fPythiaEventWeight);
						fh2dJetSignedImpParXYZb->Fill(fJetPt,cursImParXYZ,fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceb->Fill(fJetPt,cursImParXYSig,fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceb->Fill(fJetPt,cursImParXYZSig,fPythiaEventWeight);
					}
				}

					sImpParXY.push_back(cursImParXY);
					sImpParXYZ.push_back(cursImParXYZ);
					sImpParXYSig.push_back(cursImParXYSig);
					sImpParXYZSig.push_back(cursImParXYZSig);
			}// end of track loop

		EnergyFraction=SVEnergy/TracksEnergy;

		if(fDoSVEnergyFraction && EnergyFraction>0){

		      	fhistSVEnergyFraction->Fill(fJetPt,EnergyFraction ,fPythiaEventWeight);
		    		
		  	if(fIsPythia){
		      	    switch(fJetFlavor)
				{
				case 0:
					fhistSVEnergyFraction_Unidentified->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
					break;
				case 1:
					fhistSVEnergyFraction_udsg->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
					break;
				case 2:
					fhistSVEnergyFraction_c->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
					break;
				case 3:
					fhistSVEnergyFraction_b->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
					break;
				default:
					break;
				}
		    	}
		}

		if(fDoTrackCountingAnalysis){

			std::sort(sImpParXY.begin(),sImpParXY.end(), std::greater<double>());
			std::sort(sImpParXYZ.begin(),sImpParXYZ.end(), std::greater<double>());
			std::sort(sImpParXYSig.begin(),sImpParXYSig.end(), std::greater<double>());
			std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(), std::greater<double>());	

			std::vector<double> DefaultDiscriminator;

			if(fUseIPs)
				DefaultDiscriminator = sImpParXYSig;
			else
				DefaultDiscriminator = sImpParXY;

			//Ordered n=1,2,3 sip
			if (DefaultDiscriminator.size()>0){
				fh2dJetSignedImpParXYFirst->Fill(fJetPt,sImpParXY.at(0),fPythiaEventWeight);
				fh2dJetSignedImpParXYZFirst->Fill(fJetPt,sImpParXYZ.at(0),fPythiaEventWeight);
				fh2dJetSignedImpParXYSignificanceFirst->Fill(fJetPt,sImpParXYSig.at(0),fPythiaEventWeight);
				fh2dJetSignedImpParXYZSignificanceFirst->Fill(fJetPt,sImpParXYZSig.at(0),fPythiaEventWeight);

				if(DefaultDiscriminator.at(0) > fThresholdIP)
					TaggedFirst = kTRUE;

				if(fDoPtRelAnalysis && PtRelSample && DefaultDiscriminator.at(0) > fThresholdIP){
					for(int e=0; e<ntracks; e++){
					   if(ElePtRel[e]==0) break;
					   fhistPtRelVsJetPtTaggedFirst->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
					   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedFirst->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistPtRelVsJetPtTaggedUnidentifiedFirst->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedUnidentifiedFirst->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistPtRelVsJetPtTaggedudsgFirst->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedudsgFirst->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistPtRelVsJetPtTaggedcFirst->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedcFirst->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==3){
							   fhistPtRelVsJetPtTaggedbFirst->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedbFirst->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
					   }
					}
				}

				if(fDoJetProbabilityAnalysis && fResolutionFunction[0] && fValJetProb >= 0 && DefaultDiscriminator.at(0) >= fThresholdIP){
					   fhistJetProbabilityLogFirst->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistJetProbability_UnidentifiedLogFirst->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistJetProbability_udsgLogFirst->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistJetProbability_cLogFirst->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
							   if(fDoCharmFractions){
							     switch(TMath::Abs(partonpdg)){
							       case 421:
							         fhistJetProbability_cLogFirst_D0->Fill(fJetPt,fLogJetProb);
							         break;
							       case 411:
							         fhistJetProbability_cLogFirst_Dp->Fill(fJetPt,fLogJetProb);
							         break;
							       case 431:
							         fhistJetProbability_cLogFirst_Ds->Fill(fJetPt,fLogJetProb);
							         break;
							       case 4122:
							         fhistJetProbability_cLogFirst_Lc->Fill(fJetPt,fLogJetProb);
							         break;
							       default:
							         break;
							     }
							   }
						}
						else if(fJetFlavor ==3){
							   fhistJetProbability_bLogFirst->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
					   }
				}

				if(fDoJetMass && DefaultDiscriminator.at(0) >= fThresholdIP){
					   fhistJetMassFirst->Fill(fJetPt,Mass ,fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistJetMass_UnidentifiedFirst->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistJetMass_udsgFirst->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistJetMass_cFirst->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==3){
							   fhistJetMass_bFirst->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
					   }
				}

				if(fDoSVEnergyFraction && DefaultDiscriminator.at(0) >= fThresholdIP && EnergyFraction>0){

				      	fhistSVEnergyFractionFirst->Fill(fJetPt,EnergyFraction ,fPythiaEventWeight);
				    		
				  	if(fIsPythia){
				      	    switch(fJetFlavor)
						{
						case 0:
							fhistSVEnergyFraction_UnidentifiedFirst->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 1:
							fhistSVEnergyFraction_udsgFirst->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 2:
							fhistSVEnergyFraction_cFirst->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 3:
							fhistSVEnergyFraction_bFirst->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						default:
							break;
						}
				    	}
				}



				if(fIsPythia){
					
					if(fJetFlavor ==0){
						fh2dJetSignedImpParXYUnidentifiedFirst->Fill(fJetPt,sImpParXY.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZUnidentifiedFirst->Fill(fJetPt,sImpParXYZ.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceUnidentifiedFirst->Fill(fJetPt,sImpParXYSig.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst->Fill(fJetPt,sImpParXYZSig.at(0),fPythiaEventWeight);
					}
					else if(fJetFlavor ==1){
						fh2dJetSignedImpParXYudsgFirst->Fill(fJetPt,sImpParXY.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZudsgFirst->Fill(fJetPt,sImpParXYZ.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceudsgFirst->Fill(fJetPt,sImpParXYSig.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceudsgFirst->Fill(fJetPt,sImpParXYZSig.at(0),fPythiaEventWeight);						
					}
					else if(fJetFlavor ==2){
						fh2dJetSignedImpParXYcFirst->Fill(fJetPt,sImpParXY.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZcFirst->Fill(fJetPt,sImpParXYZ.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancecFirst->Fill(fJetPt,sImpParXYSig.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancecFirst->Fill(fJetPt,sImpParXYZSig.at(0),fPythiaEventWeight);						
					}
					else if(fJetFlavor ==3){
						fh2dJetSignedImpParXYbFirst->Fill(fJetPt,sImpParXY.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZbFirst->Fill(fJetPt,sImpParXYZ.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancebFirst->Fill(fJetPt,sImpParXYSig.at(0),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancebFirst->Fill(fJetPt,sImpParXYZSig.at(0),fPythiaEventWeight);
					}
				}
			}//N=1

			//Second largest
			if (DefaultDiscriminator.size()>1)
			{

				fh2dJetSignedImpParXYSecond->Fill(fJetPt,sImpParXY.at(1),fPythiaEventWeight);
				fh2dJetSignedImpParXYZSecond->Fill(fJetPt,sImpParXYZ.at(1),fPythiaEventWeight);
				fh2dJetSignedImpParXYSignificanceSecond->Fill(fJetPt,sImpParXYSig.at(1),fPythiaEventWeight);
				fh2dJetSignedImpParXYZSignificanceSecond->Fill(fJetPt,sImpParXYZSig.at(1),fPythiaEventWeight);

				if(DefaultDiscriminator.at(1) > fThresholdIP)
					TaggedSecond = kTRUE;

				if(fDoPtRelAnalysis && PtRelSample && DefaultDiscriminator.at(1) >= fThresholdIP){
					for(int e=0; e<ntracks; e++){
					   if(ElePtRel[e]==0) break;
					   fhistPtRelVsJetPtTaggedSecond->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
					   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedSecond->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistPtRelVsJetPtTaggedUnidentifiedSecond->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedUnidentifiedSecond->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistPtRelVsJetPtTaggedudsgSecond->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedudsgSecond->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistPtRelVsJetPtTaggedcSecond->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedcSecond->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==3){
							   fhistPtRelVsJetPtTaggedbSecond->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedbFirst->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
					   }
					}
				}

				if(fDoJetProbabilityAnalysis && fResolutionFunction[0] && fValJetProb >= 0 && DefaultDiscriminator.at(1) >= fThresholdIP){
					   fhistJetProbabilityLogSecond->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistJetProbability_UnidentifiedLogSecond->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistJetProbability_udsgLogSecond->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistJetProbability_cLogSecond->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
							   if(fDoCharmFractions){
							     switch(TMath::Abs(partonpdg)){
							       case 421:
							         fhistJetProbability_cLogSecond_D0->Fill(fJetPt,fLogJetProb);
							         break;
							       case 411:
							         fhistJetProbability_cLogSecond_Dp->Fill(fJetPt,fLogJetProb);
							         break;
							       case 431:
							         fhistJetProbability_cLogSecond_Ds->Fill(fJetPt,fLogJetProb);
							         break;
							       case 4122:
							         fhistJetProbability_cLogSecond_Lc->Fill(fJetPt,fLogJetProb);
							         break;
							       default:
							         break;
							     }
							   }
						}
						else if(fJetFlavor ==3){
							   fhistJetProbability_bLogSecond->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
					   }
				}


				if(fDoJetMass && DefaultDiscriminator.at(1) >= fThresholdIP){
					   fhistJetMassSecond->Fill(fJetPt,Mass ,fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistJetMass_UnidentifiedSecond->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistJetMass_udsgSecond->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistJetMass_cSecond->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==3){
							   fhistJetMass_bSecond->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
					   }
				}

				if(fDoSVEnergyFraction && DefaultDiscriminator.at(1) >= fThresholdIP && EnergyFraction>0){

				      	fhistSVEnergyFractionSecond->Fill(fJetPt,EnergyFraction ,fPythiaEventWeight);
				    		
				  	if(fIsPythia){
				      	    switch(fJetFlavor)
						{
						case 0:
							fhistSVEnergyFraction_UnidentifiedSecond->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 1:
							fhistSVEnergyFraction_udsgSecond->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 2:
							fhistSVEnergyFraction_cSecond->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 3:
							fhistSVEnergyFraction_bSecond->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						default:
							break;
						}
				    	}
				}

				if(fIsPythia){
					
					if(fJetFlavor ==0){
						fh2dJetSignedImpParXYUnidentifiedSecond->Fill(fJetPt,sImpParXY.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZUnidentifiedSecond->Fill(fJetPt,sImpParXYZ.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceUnidentifiedSecond->Fill(fJetPt,sImpParXYSig.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond->Fill(fJetPt,sImpParXYZSig.at(1),fPythiaEventWeight);
					}
					else if(fJetFlavor ==1){
						fh2dJetSignedImpParXYudsgSecond->Fill(fJetPt,sImpParXY.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZudsgSecond->Fill(fJetPt,sImpParXYZ.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceudsgSecond->Fill(fJetPt,sImpParXYSig.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceudsgSecond->Fill(fJetPt,sImpParXYZSig.at(1),fPythiaEventWeight);
					}
					else if(fJetFlavor ==2){
						fh2dJetSignedImpParXYcSecond->Fill(fJetPt,sImpParXY.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZcSecond->Fill(fJetPt,sImpParXYZ.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancecSecond->Fill(fJetPt,sImpParXYSig.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancecSecond->Fill(fJetPt,sImpParXYZSig.at(1),fPythiaEventWeight);
					}
					else if(fJetFlavor ==3)
					{
						fh2dJetSignedImpParXYbSecond->Fill(fJetPt,sImpParXY.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZbSecond->Fill(fJetPt,sImpParXYZ.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancebSecond->Fill(fJetPt,sImpParXYSig.at(1),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancebSecond->Fill(fJetPt,sImpParXYZSig.at(1),fPythiaEventWeight);
					}
				}
			}//N=2
			//Third largest

			if (DefaultDiscriminator.size()>2)
			{
				fh2dJetSignedImpParXYThird->Fill(fJetPt,sImpParXY.at(2),fPythiaEventWeight);
				fh2dJetSignedImpParXYZThird->Fill(fJetPt,sImpParXYZ.at(2),fPythiaEventWeight);
				fh2dJetSignedImpParXYSignificanceThird->Fill(fJetPt,sImpParXYSig.at(2),fPythiaEventWeight);
				fh2dJetSignedImpParXYZSignificanceThird->Fill(fJetPt,sImpParXYZSig.at(2),fPythiaEventWeight);

				if(DefaultDiscriminator.at(2) > fThresholdIP)
					TaggedThird = kTRUE;

				if(fDoPtRelAnalysis && PtRelSample && DefaultDiscriminator.at(2) >= fThresholdIP){
					for(int e=0; e<ntracks; e++){
					   if(ElePtRel[e]==0) break;
					   fhistPtRelVsJetPtTaggedThird->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
					   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedThird->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistPtRelVsJetPtTaggedUnidentifiedThird->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedUnidentifiedThird->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistPtRelVsJetPtTaggedudsgThird->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedudsgThird->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistPtRelVsJetPtTaggedcThird->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedcThird->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
						else if(fJetFlavor ==3){
							   fhistPtRelVsJetPtTaggedbThird->Fill(fJetPt, ElePtRel[e],fPythiaEventWeight);
							   if(EleIP[e]!=-999) fhistLepIPVsJetPtTaggedbThird->Fill(fJetPt, EleIP[e],fPythiaEventWeight);
						}
					   }
					}
				}

				if(fDoJetProbabilityAnalysis && fResolutionFunction[0] && fValJetProb >= 0 && DefaultDiscriminator.at(2) >= fThresholdIP){
					   fhistJetProbabilityLogThird->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistJetProbability_UnidentifiedLogThird->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistJetProbability_udsgLogThird->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistJetProbability_cLogThird->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
							   if(fDoCharmFractions){
							     switch(TMath::Abs(partonpdg)){
							       case 421:
							         fhistJetProbability_cLogThird_D0->Fill(fJetPt,fLogJetProb);
							         break;
							       case 411:
							         fhistJetProbability_cLogThird_Dp->Fill(fJetPt,fLogJetProb);
							         break;
							       case 431:
							         fhistJetProbability_cLogThird_Ds->Fill(fJetPt,fLogJetProb);
							         break;
							       case 4122:
							         fhistJetProbability_cLogThird_Lc->Fill(fJetPt,fLogJetProb);
							         break;
							       default:
							         break;
							     }
							   }
						}
						else if(fJetFlavor ==3){
							   fhistJetProbability_bLogThird->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
						}
					   }
				}

				if(fDoJetMass && DefaultDiscriminator.at(2) >= fThresholdIP){
					   fhistJetMassThird->Fill(fJetPt,Mass ,fPythiaEventWeight);
					   if(fIsPythia){
						if(fJetFlavor ==0){
							   fhistJetMass_UnidentifiedThird->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==1){
							   fhistJetMass_udsgThird->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==2){
							   fhistJetMass_cThird->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
						else if(fJetFlavor ==3){
							   fhistJetMass_bThird->Fill(fJetPt,Mass ,fPythiaEventWeight);
						}
					   }
				}

				if(fDoSVEnergyFraction && DefaultDiscriminator.at(2) >= fThresholdIP && EnergyFraction>0){

				      	fhistSVEnergyFractionThird->Fill(fJetPt,EnergyFraction ,fPythiaEventWeight);
				    		
				  	if(fIsPythia){
				      	    switch(fJetFlavor)
						{
						case 0:
							fhistSVEnergyFraction_UnidentifiedThird->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 1:
							fhistSVEnergyFraction_udsgThird->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 2:
							fhistSVEnergyFraction_cThird->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						case 3:
							fhistSVEnergyFraction_bThird->Fill(fJetPt,EnergyFraction,fPythiaEventWeight);
							break;
						default:
							break;
						}
				    	}
				}

				if(fIsPythia){

					if(fJetFlavor ==0){
						fh2dJetSignedImpParXYUnidentifiedThird->Fill(fJetPt,sImpParXY.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZUnidentifiedThird->Fill(fJetPt,sImpParXYZ.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceUnidentifiedThird->Fill(fJetPt,sImpParXYSig.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedThird->Fill(fJetPt,sImpParXYZSig.at(2),fPythiaEventWeight);
					}
					else if(fJetFlavor ==1){
						fh2dJetSignedImpParXYudsgThird->Fill(fJetPt,sImpParXY.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZudsgThird->Fill(fJetPt,sImpParXYZ.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceudsgThird->Fill(fJetPt,sImpParXYSig.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificanceudsgThird->Fill(fJetPt,sImpParXYZSig.at(2),fPythiaEventWeight);
					}
					else if(fJetFlavor ==2){
						fh2dJetSignedImpParXYcThird->Fill(fJetPt,sImpParXY.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZcThird->Fill(fJetPt,sImpParXYZ.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancecThird->Fill(fJetPt,sImpParXYSig.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancecThird->Fill(fJetPt,sImpParXYZSig.at(2),fPythiaEventWeight);
					}
					else if(fJetFlavor ==3){
						fh2dJetSignedImpParXYbThird->Fill(fJetPt,sImpParXY.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZbThird->Fill(fJetPt,sImpParXYZ.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancebThird->Fill(fJetPt,sImpParXYSig.at(2),fPythiaEventWeight);
						fh2dJetSignedImpParXYZSignificancebThird->Fill(fJetPt,sImpParXYZSig.at(2),fPythiaEventWeight);
					}
				}
			}//N=3


			//Forth largest

			if (DefaultDiscriminator.size()>3 && fDoForthIP)
			{
				fh2dJetSignedImpParXYForth->Fill(fJetPt,sImpParXY.at(3),fPythiaEventWeight);
				fh2dJetSignedImpParXYSignificanceForth->Fill(fJetPt,sImpParXYSig.at(3),fPythiaEventWeight);

				if(fIsPythia){

					if(fJetFlavor ==1){
						fh2dJetSignedImpParXYudsgForth->Fill(fJetPt,sImpParXY.at(3),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificanceudsgForth->Fill(fJetPt,sImpParXYSig.at(3),fPythiaEventWeight);
					}
					else if(fJetFlavor ==2){
						fh2dJetSignedImpParXYcForth->Fill(fJetPt,sImpParXY.at(3),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancecForth->Fill(fJetPt,sImpParXYSig.at(3),fPythiaEventWeight);
					}
					else if(fJetFlavor ==3){
						fh2dJetSignedImpParXYbForth->Fill(fJetPt,sImpParXY.at(3),fPythiaEventWeight);
						fh2dJetSignedImpParXYSignificancebForth->Fill(fJetPt,sImpParXYSig.at(3),fPythiaEventWeight);
					}
				}
			}//N=4

			if(fIsPythia && fDoTaggedDRM){

				if (jetrec->MatchedJet()) {
					  double genpt = jetrec->MatchedJet()->Pt();
					  if(!(fJetContainerMC->GetRhoParameter() == 0x0)){
						genpt = genpt - fJetContainerMC->GetRhoVal() * jetrec->MatchedJet()->Area();
					  }

					if (DefaultDiscriminator.size()>0){
					  if(DefaultDiscriminator.at(0) >= fThresholdIP)  fh2dJetGenPtVsJetRecPtFirst ->Fill(fJetPt,genpt,fPythiaEventWeight);
					}

					if (DefaultDiscriminator.size()>1){
					  if(DefaultDiscriminator.at(1) >= fThresholdIP)  fh2dJetGenPtVsJetRecPtSecond->Fill(fJetPt,genpt,fPythiaEventWeight);
					}

					if (DefaultDiscriminator.size()>2){
					  if(DefaultDiscriminator.at(2) >= fThresholdIP)  fh2dJetGenPtVsJetRecPtThird ->Fill(fJetPt,genpt,fPythiaEventWeight);
					}
				}
			}

			sImpParXY.clear();
			sImpParXYZ.clear();
			sImpParXYSig.clear();
			sImpParXYZSig.clear();
			DefaultDiscriminator.clear();
		}//end track counting

		// Secondary Vertex Tagger, for testing the DATA driven approach
		if(fDoSVAnalysis){

			  Double_t vtxPos[3]   = {fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ()};
			  Double_t covMatrix[6] = {0};
			  fPrimaryVertex->GetCovarianceMatrix(covMatrix);
			  AliESDVertex* esdVtx = new AliESDVertex(vtxPos, covMatrix, fPrimaryVertex->GetChi2(), fPrimaryVertex->GetNContributors());

			  // 3 Prong Vertex
			  TClonesArray* secVertexArr3Prong = 0;
			  vector<pair <Double_t, Int_t>> arrDispersion3Prong;
			  arrDispersion3Prong.reserve(5);

			  secVertexArr3Prong = new TClonesArray("AliAODVertex");
			  Int_t nDauRejCount3Prong = 0;
			  Int_t nVtx3Prong = fVtxTagger3Prong->FindVertices(jetrec,
						                 static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
						                 fAODIn,
						                 esdVtx,
						                 fAODIn->GetMagneticField(),
						                 secVertexArr3Prong,
						                 0,
						                 arrDispersion3Prong,
						                 nDauRejCount3Prong);

			  if(nVtx3Prong > 0)
			  {
				Double_t *decLenXY        = new Double_t[nVtx3Prong];
    				Double_t *errdecLenXY     = new Double_t[nVtx3Prong];
			    	Double_t *sigdecLenXY     = new Double_t[nVtx3Prong];
			    	Double_t *invMasses       = new Double_t[nVtx3Prong];
			    	Double_t *sigmavertex     = new Double_t[nVtx3Prong];

    				Int_t    *idxLxy = new Int_t[nVtx3Prong];

				  for(Int_t iv=0; iv<secVertexArr3Prong->GetEntriesFast(); iv++)
				  {
				    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArr3Prong->UncheckedAt(iv));

				    // Calculate vtx distance
				    Double_t effX = secVtx->GetX() - esdVtx->GetX();
				    Double_t effY = secVtx->GetY() - esdVtx->GetY();
				    //Double_t effZ = secVtx->GetZ() - esdVtx->GetZ();

				    // ##### Vertex properties
				    // vertex dispersion
				    sigmavertex[iv] = arrDispersion3Prong[iv].first;

				    // invariant mass
				    invMasses[iv] = fVtxTagger3Prong->GetVertexInvariantMass(secVtx);

				    // signed length
				    decLenXY[iv]  = fPrimaryVertex->DistanceXYToVertex(secVtx);
				    Double_t jetP[3]; jetrec->PxPyPz(jetP);
				    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
				    if (signLxy < 0.) decLenXY[iv] *= -1.;

				    errdecLenXY[iv] = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

				    sigdecLenXY[iv] = decLenXY[iv]/errdecLenXY[iv];

				    //cout << Form("(%3.3f, %3.3f, %3.3f), r=%3.3f, chi2=%3.3f, dispersion=%3.3f", effX, effY, effZ, vtxDistance, TMath::Abs(vtx->GetChi2perNDF()), dispersion) << endl;
				  }

				    TMath::Sort(nVtx3Prong, decLenXY, idxLxy);

				      	fInvariantMass  = invMasses[idxLxy[0]];
				      	fDispersion  = sigmavertex[idxLxy[0]];
				      	fDecayLength  = decLenXY[idxLxy[0]];
				      	fLxySign = sigdecLenXY[idxLxy[0]];

					fHistDispersion3Prong->Fill(fJetPt,fDispersion);
					if(fIsPythia){
						if(fJetFlavor==0) fHistDispersion3ProngUnidentified->Fill(fJetPt,fDispersion);
						if(fJetFlavor==1) fHistDispersion3Pronglf->Fill(fJetPt,fDispersion);
						if(fJetFlavor==2) fHistDispersion3Prongc->Fill(fJetPt,fDispersion);
						if(fJetFlavor==3) fHistDispersion3Prongb->Fill(fJetPt,fDispersion);
					}
					if(fDispersion < 0.04){
	  					fHistSV3Prong->Fill(fJetPt,fLxySign,fInvariantMass);
						if(fIsPythia){
							if(fJetFlavor==0) fHistSV3ProngUnidentified->Fill(fJetPt,fLxySign,fInvariantMass);
							if(fJetFlavor==1) fHistSV3Pronglf->Fill(fJetPt,fLxySign,fInvariantMass);
							if(fJetFlavor==2) fHistSV3Prongc->Fill(fJetPt,fLxySign,fInvariantMass);
							if(fJetFlavor==3) fHistSV3Prongb->Fill(fJetPt,fLxySign,fInvariantMass);
						}
						if(fDoJetProbabilityAnalysis && fResolutionFunction[0] && fValJetProb >= 0 && fLxySign >= 5.){
							   fhistJetProbabilityLogSVHP->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
							   if(fIsPythia){
								if(fJetFlavor ==0){
									   fhistJetProbability_UnidentifiedLogSVHP->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
								else if(fJetFlavor ==1){
									   fhistJetProbability_udsgLogSVHP->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
								else if(fJetFlavor ==2){
									   fhistJetProbability_cLogSVHP->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
								else if(fJetFlavor ==3){
									   fhistJetProbability_bLogSVHP->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
							   }
						}
					}

				    secVertexArr3Prong->Clear();
			  	    delete secVertexArr3Prong;
				    delete[] idxLxy;
				    delete[] decLenXY;
				    delete[] invMasses;
				    delete[] errdecLenXY;
				    delete[] sigdecLenXY;
				    delete[] sigmavertex;
			  }else{
			    secVertexArr3Prong->Clear();
			    delete secVertexArr3Prong;
			  }

			// 2 Prong Vertex
			  TClonesArray* secVertexArr2Prong = 0;
			  vector<pair <Double_t, Int_t>> arrDispersion2Prong;
			  arrDispersion2Prong.reserve(5);

			  secVertexArr2Prong = new TClonesArray("AliAODVertex");
			  Int_t nDauRejCount2Prong = 0;
			  Int_t nVtx2Prong = fVtxTagger2Prong->FindVertices(jetrec,
						                 static_cast<AliParticleContainer*>(fParticleCollArray.At(0))->GetArray(),
						                 fAODIn,
						                 esdVtx,
						                 fAODIn->GetMagneticField(),
						                 secVertexArr2Prong,
						                 0,
						                 arrDispersion2Prong,
						                 nDauRejCount2Prong);


			  if(nVtx2Prong > 0)
			  {
				Double_t *decLenXY        = new Double_t[nVtx2Prong];
    				Double_t *errdecLenXY     = new Double_t[nVtx2Prong];
			    	Double_t *sigdecLenXY     = new Double_t[nVtx2Prong];
			    	Double_t *invMasses       = new Double_t[nVtx2Prong];
			    	Double_t *sigmavertex     = new Double_t[nVtx2Prong];

    				Int_t    *idxLxy = new Int_t[nVtx2Prong];

				  for(Int_t iv=0; iv<secVertexArr2Prong->GetEntriesFast(); iv++)
				  {
				    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArr2Prong->UncheckedAt(iv));

				    // Calculate vtx distance
				    Double_t effX = secVtx->GetX() - esdVtx->GetX();
				    Double_t effY = secVtx->GetY() - esdVtx->GetY();
				    //Double_t effZ = secVtx->GetZ() - esdVtx->GetZ();

				    // ##### Vertex properties
				    // vertex dispersion
				    sigmavertex[iv] = arrDispersion2Prong[iv].first;

				    // invariant mass
				    invMasses[iv] = fVtxTagger2Prong->GetVertexInvariantMass(secVtx);

				    // signed length
				    decLenXY[iv]  = fPrimaryVertex->DistanceXYToVertex(secVtx);
				    Double_t jetP[3]; jetrec->PxPyPz(jetP);
				    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
				    if (signLxy < 0.) decLenXY[iv] *= -1.;

				    errdecLenXY[iv] = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

				    sigdecLenXY[iv] = decLenXY[iv]/errdecLenXY[iv];

				    //cout << Form("(%3.3f, %3.3f, %3.3f), r=%3.3f, chi2=%3.3f, dispersion=%3.3f", effX, effY, effZ, vtxDistance, TMath::Abs(vtx->GetChi2perNDF()), dispersion) << endl;
				  }

				    TMath::Sort(nVtx2Prong, decLenXY, idxLxy);

					fInvariantMass  = invMasses[idxLxy[0]];
				      	fDispersion  = sigmavertex[idxLxy[0]];
				      	fDecayLength  = decLenXY[idxLxy[0]];
				      	fLxySign = sigdecLenXY[idxLxy[0]];

					fHistDispersion2Prong->Fill(fJetPt,fDispersion);
					if(fIsPythia){
						if(fJetFlavor==0) fHistDispersion2ProngUnidentified->Fill(fJetPt,fDispersion);
						if(fJetFlavor==1) fHistDispersion2Pronglf->Fill(fJetPt,fDispersion);
						if(fJetFlavor==2) fHistDispersion2Prongc->Fill(fJetPt,fDispersion);
						if(fJetFlavor==3) fHistDispersion2Prongb->Fill(fJetPt,fDispersion);
					}
					if(fDispersion < 0.04){
	  					fHistSV2Prong->Fill(fJetPt,fLxySign,fInvariantMass);

						if(fIsPythia){
							if(fJetFlavor==0) fHistSV2ProngUnidentified->Fill(fJetPt,fLxySign,fInvariantMass);
							if(fJetFlavor==1) fHistSV2Pronglf->Fill(fJetPt,fLxySign,fInvariantMass);
							if(fJetFlavor==2) fHistSV2Prongc->Fill(fJetPt,fLxySign,fInvariantMass);
							if(fJetFlavor==3) fHistSV2Prongb->Fill(fJetPt,fLxySign,fInvariantMass);
						}
						if(fDoJetProbabilityAnalysis && fResolutionFunction[0] && fValJetProb >= 0 && fLxySign >= 5.){
							   fhistJetProbabilityLogSVHE->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
							   if(fIsPythia){
								if(fJetFlavor ==0){
									   fhistJetProbability_UnidentifiedLogSVHE->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
								else if(fJetFlavor ==1){
									   fhistJetProbability_udsgLogSVHE->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
								else if(fJetFlavor ==2){
									   fhistJetProbability_cLogSVHE->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
								else if(fJetFlavor ==3){
									   fhistJetProbability_bLogSVHE->Fill(fJetPt, fLogJetProb,fPythiaEventWeight);
								}
							   }
						}
					}

			  	    secVertexArr2Prong->Clear();
			  	    delete secVertexArr2Prong;
				    delete[] idxLxy;
				    delete[] decLenXY;
				    delete[] invMasses;
				    delete[] errdecLenXY;
				    delete[] sigdecLenXY;
				    delete[] sigmavertex;
			  }else{
			    secVertexArr2Prong->Clear();
			    delete secVertexArr2Prong;
			  }

			  delete esdVtx;

		}//end SV

		fJetPt = 0.;
		fJetMass=0.;
		trackAOD = 0x0;
		delete trackAOD;

	}//End jet loop
	jetrec = 0x0; jetmatched=0x0;
	delete jetmatched; delete jetrec;

	if(TaggedFirst){
		f2histRhoVsDeltaPtFirst->Fill( randomConePt , fJetContainerData->GetRhoVal(), fPythiaEventWeight);
		if(fDoDeltaPtWithSignal) f2histRhoVsDeltaPtWithSignalFirst->Fill( randomConePtWithSignal , fJetContainerData->GetRhoVal(), fPythiaEventWeight);
	}
	if(TaggedSecond){
		f2histRhoVsDeltaPtSecond->Fill( randomConePt , fJetContainerData->GetRhoVal(), fPythiaEventWeight);
		if(fDoDeltaPtWithSignal) f2histRhoVsDeltaPtWithSignalSecond->Fill( randomConePtWithSignal , fJetContainerData->GetRhoVal(), fPythiaEventWeight);
	}
	if(TaggedThird){
		f2histRhoVsDeltaPtThird->Fill( randomConePt , fJetContainerData->GetRhoVal(), fPythiaEventWeight);
		if(fDoDeltaPtWithSignal) f2histRhoVsDeltaPtWithSignalThird->Fill( randomConePtWithSignal , fJetContainerData->GetRhoVal(), fPythiaEventWeight);
	}


	if(fEnableV0GammaRejection){
		if( fIsPythia > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
		    fV0Reader->RelabelAODs(kFALSE);
		}
	}

	return kTRUE;
}
// ######################################################################################## JEt Probability Function
Double_t AliAnalysisTaskBJetTC::CalculateTrackProb(Double_t significance, Int_t trclass, Int_t jetFlavor)
{
  Double_t trackprob = 0;
  //switch resolution function based on track pt;
  if(TMath::Abs(significance) >100) significance =99.9; //Limit to function definition range
  if(fIsPythia){
 	 switch(jetFlavor)
	 {
		case 0:
  		      trackprob = fResolutionFunctionlf[trclass]->Integral(-100,-TMath::Abs(significance))/fResolutionFunctionlf[trclass]->Integral(-100,0);
		      break;
		case 1:
  		      trackprob = fResolutionFunctionlf[trclass]->Integral(-100,-TMath::Abs(significance))/fResolutionFunctionlf[trclass]->Integral(-100,0);
		      break;
		case 2:
  		      trackprob = fResolutionFunctionc[trclass]->Integral(-100,-TMath::Abs(significance))/fResolutionFunctionc[trclass]->Integral(-100,0);
			break;
		case 3:
  		      trackprob = fResolutionFunctionb[trclass]->Integral(-100,-TMath::Abs(significance))/fResolutionFunctionb[trclass]->Integral(-100,0);
			break;
		default:
			break;
	  }
  }else
  	trackprob = fResolutionFunction[trclass]->Integral(-100,-TMath::Abs(significance))/fResolutionFunction[trclass]->Integral(-100,0);

  if(fMinTrackProb)  trackprob=TMath::Max(trackprob,fMinTrackProb);
  return trackprob;
}
// ######################################################################################## Jet Probability Function
Double_t AliAnalysisTaskBJetTC::CalculateJetProb(AliEmcalJet *jet, Int_t jetFlavor)
{
  if(!jet) return -1;
  Double_t JetProb = -1;
  //Loop over all tracks calculate P(s) for all accepted later add looser cuts also
  Int_t ntracks = (Int_t)jet->GetNumberOfTracks();
  Double_t TrackProb = 1;
  Double_t curps=-1;
  Int_t trackcounter =0;
  for(Int_t itrack = 0; itrack < ntracks; ++itrack)
    {
      AliAODTrack* trackV = 0x0;
      if(fUsePicoTracks) trackV = (AliAODTrack*)((AliPicoTrack*)jet->Track(itrack))->GetTrack();
      else trackV = (AliAODTrack*)((fJetContainerData->GetParticleContainer())->GetParticle(jet->TrackAt(itrack)));

      if(!trackV) continue;

      //class selection
      Int_t QualityClass=0;
      double dca[2] = {0};
      double cov[3] = {0};
      double sign = 0;

      if(!IsTrackAcceptedQuality(trackV, jet, QualityClass, dca, cov, sign)) continue;
      if(sign<0) continue;//only take positive IP tracks
      curps =CalculateTrackProb(TMath::Abs(GetValImpactParameter(kXYSig,dca,cov)), QualityClass, jetFlavor);
      TrackProb*=TMath::Abs(curps);
      trackcounter++;
    }

  if(trackcounter<2) return -1;

  Double_t sumPS =0;
  bool chan=false;
  for(Int_t j=0;j<trackcounter;++j){
      double val = TMath::Power(-1 * TMath::Log(TrackProb),j)/TMath::Factorial(j);
      sumPS += val;
      chan=true;
    }
  if(!chan)return -1;
  JetProb =sumPS *TrackProb;
  return JetProb;
}
//###############################################################################################################

void AliAnalysisTaskBJetTC::FillResolutionFunctionHists(AliAODTrack * track,AliEmcalJet * jet, Int_t jetFlavor)
{
 
   Int_t QualityClass=0;
   double imp[2] = {0};
   double cov[3] = {0};
   double sign = 0;
 
   if(!IsTrackAcceptedQuality(track, jet, QualityClass, imp, cov, sign)) return;

  Double_t weight=1;//Weight of the tracks DATA/MC

	switch(QualityClass){
	case 0:	//Chi2perNDF()>2
          	fh2dJetSignedImpParXY_Class1->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
          	fh2dJetSignedImpParXYSignificance_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
          	fh2dJetSignedImpParXYZ_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
          	fh2dJetSignedImpParXYZSignificance_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class1->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		       }
		}
		break;
	case 1: //Pt<2 nITShits=2
		fh2dJetSignedImpParXY_Class2->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYSignificance_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZ_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZSignificance_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		        	}
		}
		break;
	case 3: //Pt<2 nITShits=3
		fh2dJetSignedImpParXY_Class3->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYSignificance_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZ_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZSignificance_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		        	}
		}
		break;
	case 5: //Pt<2 nITShits=4
		fh2dJetSignedImpParXY_Class4->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYSignificance_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZ_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZSignificance_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		        	}
		}
		break;
	case 2: //Pt>2 nITShits=2
		fh2dJetSignedImpParXY_Class2->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYSignificance_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZ_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZSignificance_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class2->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		        	}
		}
		break;
	case 4: //Pt>2 nITShits=3
		fh2dJetSignedImpParXY_Class3->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYSignificance_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZ_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZSignificance_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class3->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		        	}
		}
		break;
	case 6: //Pt>2 nITShits=4
		fh2dJetSignedImpParXY_Class4->Fill(track->Pt() , TMath::Abs(GetValImpactParameter(kXY, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYSignificance_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZ_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZ, imp,cov))*sign,weight);
		fh2dJetSignedImpParXYZSignificance_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYZSig, imp,cov))*sign,weight);
		if(fIsPythia){
		    	switch(jetFlavor)
			{
		        	case 0:
		  			fh2dJetSignedImpParXYSignificancelf_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 1:
		  			fh2dJetSignedImpParXYSignificancelf_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 2:
		  			fh2dJetSignedImpParXYSignificancec_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
		          		break;
		        	case 3:
		  			fh2dJetSignedImpParXYSignificanceb_Class4->Fill(track->Pt(),TMath::Abs(GetValImpactParameter(kXYSig, imp,cov))*sign,weight);
					break;
		        	default:
		          		break;
		        	}
		}
		break;

	default:
		break;
	}

  return;
}

// ######################################################################################## Event Selection
Bool_t AliAnalysisTaskBJetTC::IsEventSelected()	{

	fAODIn = NULL;

	Int_t WhyRejected =0;
	ULong_t RejectionBits=0;

	fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());

	fPrimaryVertex = dynamic_cast<AliAODVertex*>(fAODIn->GetPrimaryVertex());

		if(fAODIn && fPrimaryVertex && fPrimaryVertex->GetNContributors()>0){
			fh1dVertexZ->Fill(fPrimaryVertex->GetZ());
			double vtxx = fPrimaryVertex->GetX();
			double vtxy = fPrimaryVertex->GetY();
			fh1dVertexR->Fill(vtxx,vtxy);
		}else return kFALSE;

		if(!(IsSelected(WhyRejected,RejectionBits)))
		{
			fh1dEventRejectionRDHFCuts->SetBinContent(2,fh1dEventRejectionRDHFCuts->GetBinContent(2)+1);
			if(RejectionBits&(1<<kPhysicsSelection))fh1dEventRejectionRDHFCuts->SetBinContent(3,fh1dEventRejectionRDHFCuts->GetBinContent(3)+1);
			else if(RejectionBits&(1<<kNoVertex))fh1dEventRejectionRDHFCuts->SetBinContent(6,fh1dEventRejectionRDHFCuts->GetBinContent(6)+1);
			else if(RejectionBits&(1<<kNoVertexTracks))fh1dEventRejectionRDHFCuts->SetBinContent(11,fh1dEventRejectionRDHFCuts->GetBinContent(11)+1);
			else if(RejectionBits&(1<<kNoContributors))fh1dEventRejectionRDHFCuts->SetBinContent(12,fh1dEventRejectionRDHFCuts->GetBinContent(12)+1);
			else if(RejectionBits&(1<<kTooFewVtxContrib))fh1dEventRejectionRDHFCuts->SetBinContent(7,fh1dEventRejectionRDHFCuts->GetBinContent(7)+1);
			else if(RejectionBits&(1<<kVertexZContrib))fh1dEventRejectionRDHFCuts->SetBinContent(15,fh1dEventRejectionRDHFCuts->GetBinContent(15)+1);
			else if(RejectionBits&(1<<kVertexZResolution))fh1dEventRejectionRDHFCuts->SetBinContent(14,fh1dEventRejectionRDHFCuts->GetBinContent(14)+1);
			else if(RejectionBits&(1<<kDeltaVertexZ))fh1dEventRejectionRDHFCuts->SetBinContent(13,fh1dEventRejectionRDHFCuts->GetBinContent(13)+1);
			else if(RejectionBits&(1<<kZVtxOutFid))fh1dEventRejectionRDHFCuts->SetBinContent(5,fh1dEventRejectionRDHFCuts->GetBinContent(5)+1);
			else if(RejectionBits&(1<<kOutsideCentrality))fh1dEventRejectionRDHFCuts->SetBinContent(4,fh1dEventRejectionRDHFCuts->GetBinContent(4)+1);
			else if(RejectionBits&(1<<kNotSelTrigger))fh1dEventRejectionRDHFCuts->SetBinContent(8,fh1dEventRejectionRDHFCuts->GetBinContent(8)+1);
			else if(RejectionBits&(1<<kSPDClusterCut))fh1dEventRejectionRDHFCuts->SetBinContent(9,fh1dEventRejectionRDHFCuts->GetBinContent(9)+1);
			else if(RejectionBits&(1<<kMVPileup))fh1dEventRejectionRDHFCuts->SetBinContent(10,fh1dEventRejectionRDHFCuts->GetBinContent(10)+1);
			else if(RejectionBits&(1<<kSelPtHardBin))fh1dEventRejectionRDHFCuts->SetBinContent(16,fh1dEventRejectionRDHFCuts->GetBinContent(16)+1);

			return kFALSE;
		}else {
			fh1dEventRejectionRDHFCuts->SetBinContent(1,fh1dEventRejectionRDHFCuts->GetBinContent(1)+1);
			fh1dVertexZAccepted->Fill(fPrimaryVertex->GetZ());
			double vtxx = fPrimaryVertex->GetX();
			double vtxy = fPrimaryVertex->GetY();
			fh1dVertexRAccepted->Fill(vtxx,vtxy);
			fh2dVertexChi2NDFNESDTracks->Fill(fPrimaryVertex->GetChi2perNDF(),fAODIn->GetNumberOfESDTracks());
			return kTRUE;
		}
	return kFALSE;
}
//############################################################################################################### Event Selection
Bool_t AliAnalysisTaskBJetTC::IsSelected(Int_t &WhyRejected,ULong_t &RejectionBits){
	WhyRejected =0;
	Int_t fMinVtxContr=1;
	Bool_t accept=kTRUE;
	Int_t fMinContrPileup = 5;
	Float_t fMinDzPileup = 0.8;
	Double_t fMaxVtxZ = 10;
	RejectionBits=0;
	//Physics Selection Cut


	    Bool_t isSelected = kTRUE;/*(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);*/
	    if(!isSelected) {
	      if(accept) WhyRejected=7;
	      RejectionBits+=1<<kPhysicsSelection;
	      accept=kFALSE;
	    }

	    // vertex requirements
	    if(!fPrimaryVertex || fPrimaryVertex->GetNContributors() ==0){
	      accept=kFALSE;
	      if(!fPrimaryVertex)
	    	   RejectionBits+=1<<kNoVertex;
	      else {
	    	    RejectionBits+=1<<kNoContributors;
	      }
	    }else{
	    	const AliVVertex* trkVtx = dynamic_cast<const AliVVertex*>(fAODIn->GetPrimaryVertex());
	    	const AliVVertex* spdVtx = dynamic_cast<const AliVVertex*>(fAODIn->GetPrimaryVertexSPD());
	    	TString vtxTtl = trkVtx->GetTitle();
	    	if(!vtxTtl.Contains("VertexerTracks"))
	    	{
	    		  accept=kFALSE;
	    		  RejectionBits+=1<<kNoVertexTracks;

	    	}
		if(vtxTtl.Contains("WithConstraint")) fVertexConstraint = kTRUE;
	    	if(trkVtx->GetNContributors()<1)	{
	    	    		 accept=kFALSE;
	    	    		 RejectionBits+=1<<kTooFewVtxContrib;
	    	    	}
	    	Double_t cov[6] = { 0 };
	    	spdVtx->GetCovarianceMatrix(cov);
	    	Double_t zRes = TMath::Sqrt(cov[5]);
	    	if(spdVtx->IsFromVertexerZ() && (zRes > 0.25)) {
	    		  accept=kFALSE;
	    		  RejectionBits+=1<<kVertexZResolution;
	    	}

	    	if((TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ()) > 0.5)) {
	  		  accept=kFALSE;
	  		  RejectionBits+=1<<kDeltaVertexZ;
	    	}
	    	if(spdVtx->GetNContributors() <1) {
	    		 accept=kFALSE;
	    		RejectionBits+=1<<kVertexZContrib;
	    	}

	    	 if(TMath::Abs(trkVtx->GetZ())>=fMaxVtxZ) {
	    		        RejectionBits+=1<<kZVtxOutFid;
	    		        if(accept) WhyRejected=6;
	    		        accept=kFALSE;
	    	 }
	    }

	        Int_t cutc=(Int_t)fMinContrPileup;
	        Double_t cutz=(Double_t)fMinDzPileup;
	        /*if(fAODIn->IsPileupFromSPD(5, 0.8, 3.0, 2.0, 5.0)) {
	          if(accept) WhyRejected=1;
	          RejectionBits+=1<<kPileupSPD;
	          accept=kFALSE;
	        }*/
	//Special out-of bunch pileup cuts
	    	// SPD Cluster vs Tracklet plot to estimate pileup effect
	    	/*Int_t nClustersLayer0 = fAODIn->GetNumberOfITSClusters(0);
	    	Int_t nClustersLayer1 = fAODIn->GetNumberOfITSClusters(1);
	    	Int_t nTracklets = fAODIn->GetMultiplicity()->GetNumberOfTracklets();
			if(nClustersLayer0 + nClustersLayer1 > 65 + 4 * nTracklets){
				 accept=kFALSE;
				RejectionBits+=1<<kSPDClusterCut;
			}*/

    		Double_t  minContributors=5, minChi2=5., minWeiZDiff=15, checkPlpFromDifferentBC=kFALSE; 
		fUtils->SetMinPlpContribMV(minContributors);
   		fUtils->SetMaxPlpChi2MV(minChi2);
   		fUtils->SetMinWDistMV(minWeiZDiff);
   		fUtils->SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC);

		if(fUtils->IsPileUpMV(fAODIn)){
			accept=kFALSE;
			RejectionBits+=1<<kMVPileup;
		}

	if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin)  {
    		RejectionBits+=1<<kSelPtHardBin;
    		accept=kFALSE;
  	}
	if (fPtHardThreshold > 0){
		if(fPtHard > fPtHardThreshold){
			RejectionBits+=1<<kSelPtHardBin;
	    		accept=kFALSE;
		}
	}

	if(fIsPythia){
		if(!CheckMCOutliers()) {
			 	RejectionBits+=1<<kSelPtHardBin;
				accept=kFALSE;
		}
		if (fPtHardAndTrackPtFactor > 0.) {
			AliParticleContainer* mcpartcont = dynamic_cast<AliParticleContainer*>(fParticleCollArray.At(0));
			if ((Bool_t)mcpartcont) {
			      for (auto mctrack : mcpartcont->all()) {// Not cuts applied ; use accept for cuts
				Float_t trackpt = mctrack->Pt();
				if (trackpt > (fPtHardAndTrackPtFactor * fPtHard) ) {
				  AliInfo(Form("Reject : track %2.2f, factor %2.2f, ptHard %f", trackpt, fPtHardAndTrackPtFactor, fPtHard));
				  RejectionBits+=1<<kSelPtHardBin;
				  accept=kFALSE;
				}
			      }
			}
			mcpartcont=NULL; delete mcpartcont;
		}
	}	


return accept;
}
// ######################################################################################## Init histograms
void AliAnalysisTaskBJetTC::UserCreateOutputObjects(){

	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TString BjetCuts[15] = {
    "all"/*0*/,
    "FilterBit 4"/*1*/,
    "p_{T} cut"/*2*/,
    "#eta cut"/*3*/,
    "TPC refit"/*4*/,
    "ITS refit"/*5*/,
    "SPD hits"/*6*/,
    "TPC clusters"/*7*/,
    "ITS clusters"/*8*/,
    "Chi2/NDF"/*9*/,
    "IP_{xy}"/*10*/,
    "IP_{z}"/*11*/,
    "DecayLength"/*12*/,
    "DCA Track-Jet"/*13*/,
    "+IP_{xy}"/*14*/
  };

  if(fDoSVAnalysis || fDoSVEnergyFraction){
	  fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
	  fEsdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	  fEsdTrackCuts->SetMinNClustersTPC(80);
	  fEsdTrackCuts->SetMaxChi2PerClusterTPC(4);
	  fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
	  fEsdTrackCuts->SetRequireITSRefit(kTRUE);
	  fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	  fEsdTrackCuts->SetMinDCAToVertexXY(0.008);
	  fEsdTrackCuts->SetEtaRange(-0.9, 0.9);
	  fEsdTrackCuts->SetPtRange(1., 1.e10);


	  fjetCuts3Prong = new AliRDHFJetsCutsVertex("jetCuts3Prong");
	  fjetCuts3Prong->AddTrackCuts(fEsdTrackCuts);
	  fjetCuts3Prong->SetNprongs(3);
	  fjetCuts3Prong->SetMinPtHardestTrack(1.);

	  fjetCuts2Prong = new AliRDHFJetsCutsVertex("jetCuts2Prong");
	  fjetCuts2Prong->AddTrackCuts(fEsdTrackCuts);
	  fjetCuts2Prong->SetNprongs(2);
	  fjetCuts2Prong->SetMinPtHardestTrack(1.);

	  fVtxTagger3Prong = new AliHFJetsTaggingVertex();
	  fVtxTagger3Prong->SetCuts(fjetCuts3Prong);

	  fVtxTagger2Prong = new AliHFJetsTaggingVertex();
	  fVtxTagger2Prong->SetCuts(fjetCuts2Prong);

  	  fTrackArray = new TObjArray();
  }

  if(fDoImprovedDCACut){
	fDecayVertex = new AliAnalysisTaskWeakDecayVertexer();
	fDecayVertex->SetDoImprovedDCAV0DauPropagation(kTRUE);
  }

  if(fEnableV0GammaRejection){

	  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

	  if(fV0Reader)
	    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
	      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
		fOutput->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

	  if(fV0Reader)
	    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
	      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
		fOutput->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	  
	  if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
	    if (fV0Reader->GetV0FindingEfficiencyHistograms())
	      fOutput->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

	  if(fV0Reader && fV0Reader->GetProduceImpactParamHistograms())fOutput->Add(fV0Reader->GetImpactParamHistograms());
  }

  if(fApplyV0Rec){

	  Bool_t fbIsPbPb = kFALSE;

	  printf("-------------------------------------------------------\n");
	  if(fbTPCRefit) printf("TPC refit for daughter tracks\n");
	  if(fbRejectKinks) printf("reject kink-like production vertices of daughter tracks\n");
	  if(fbFindableClusters) printf("require positive number of findable clusters\n");
	  if(fdCutNCrossedRowsTPCMin > 0.) printf("min number of crossed TPC rows: %g\n", fdCutNCrossedRowsTPCMin);
	  if(fdCutCrossedRowsOverFindMin > 0.) printf("min ratio crossed rows / findable clusters: %g\n", fdCutCrossedRowsOverFindMin);
	  if(fdCutCrossedRowsOverFindMax > 0.) printf("max ratio crossed rows / findable clusters: %g\n", fdCutCrossedRowsOverFindMax);
	  if(fdCutPtDaughterMin > 0.) printf("min pt of daughter tracks [GeV/c]: %g\n", fdCutPtDaughterMin);
	  if(fdCutDCAToPrimVtxMin > 0.) printf("min DCA of daughters to the prim vtx [cm]: %g\n", fdCutDCAToPrimVtxMin);
	  if(fdCutDCADaughtersMax > 0.) printf("max DCA between daughters [sigma of TPC tracking]: %g\n", fdCutDCADaughtersMax);
	  if(fdCutEtaDaughterMax > 0.) printf("max |eta| of daughter tracks: %g\n", fdCutEtaDaughterMax);
	  if(fdCutNSigmadEdxMax > 0. && (!fbIsPbPb || (fbIsPbPb && fdPtProtonPIDMax > 0.))) printf("max |Delta(dE/dx)| in the TPC [sigma dE/dx]: %g\n", fdCutNSigmadEdxMax);
	  if(fdCutNSigmadEdxMax > 0. && fbIsPbPb && fdPtProtonPIDMax > 0.) printf("max pt of proton for applying PID cut [GeV/c]: %g\n", fdPtProtonPIDMax);
	  printf("V0 reconstruction method: %s\n", fbOnFly ? "on-the-fly" : "offline");
	  if(fdCutCPAKMin > 0.) printf("min CPA, K0S: %g\n", fdCutCPAKMin);
	  if(fdCutCPALMin > 0.) printf("min CPA, (A)Lambda: %g\n", fdCutCPALMin);
	  if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.) printf("R of the decay vertex [cm]: %g-%g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
	  if(fdCutEtaV0Max > 0.) printf("max |eta| of V0: %g\n", fdCutEtaV0Max);
	  if(fdCutRapV0Max > 0.) printf("max |y| of V0: %g\n", fdCutRapV0Max);
	  if(fdCutNTauKMax > 0.) printf("max proper lifetime, K0S [tau]: %g\n", fdCutNTauKMax);
	  if(fdCutNTauLMax > 0.) printf("max proper lifetime, (A)Lambda [tau]: %g\n", fdCutNTauLMax);
	  if(fbCutArmPod) printf("Armenteros-Podolanski cut for K0S\n");
	  if(fbCutCross) printf("cross-contamination cut\n");
	  printf("-------------------------------------------------------\n");

	  if(fV0CandidateArray == NULL){
      		fV0CandidateArray = new TClonesArray("AliAODv0",100);
    	  }
          fV0CandidateArray->Delete();//Reset the TClonesArray
  }

	const Int_t nBins2dSignificance =400;
	const Int_t nBins3dSignificance =250;
	const Int_t nBins2d=500;
	const Int_t nBins3d =250;

	if(fIsPythia) fHFJetUtils = new AliHFJetsTagging("fHFJetUtils");

	if (!fOutput) fOutput = new AliEmcalList();
	fOutput->SetOwner(kTRUE);

	if(fApplyV0Rec || fEnableV0GammaRejection || fDoPtRelAnalysis){
		AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
		AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
		fRespoPID = inputHandler->GetPIDResponse();
	}

	
	  // binning in jets
	  const Int_t iNDimInJC = 4;
	  Int_t binsKInJC[iNDimInJC] = {200, 200, 200, 200};
	  Double_t xminKInJC[iNDimInJC] = {0.35, 0., -1., 0.};
	  Double_t xmaxKInJC[iNDimInJC] = {0.65, 50., 1., 200.};
	  Int_t binsLInJC[iNDimInJC] = {200, 200, 200, 200};
	  Double_t xminLInJC[iNDimInJC] = {1.05, 0., -1., 0.};
	  Double_t xmaxLInJC[iNDimInJC] = {1.25, 50., 1., 200.};


	//Event selection histograms
	fh1dEventRejectionRDHFCuts = new TH1D("fh1dEventRejectionRDHFCuts;reason;count","Rejection reasons",16,0,16);
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(1,"Accepted");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(2,"Rejected");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(3,"DueToPhysicsSelection");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(4,"DueCentralitySelection");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(5,"ZVertexOutsideFiducialRegion");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(6,"IsEventRejectedDueToNotRecoVertex");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(7,"IsEventRejectedDueToVertexContributors");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(8,"DueToTrigger");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(9,"DueToSPDTrackletClusterCut");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(10,"DueToMVPileup");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(11,"NoVertexTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(12,"NoContributorsVertexTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(13,"DeltaVertexZSPDTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(14,"ZVertexResolution");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(15,"VertexZContributors");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(16,"SelPtHardBin");

	fhistInclusiveJetCuts = new TH1D("fhistInclusiveJetCuts", "Number of Inclusive jets after cuts", 15, 0, 15);
	if(fIsPythia){
		fhistbJetCuts = new TH1D("fhistbJetCuts", "Number of b jets after cuts", 15, 0, 15);
		fhistcJetCuts = new TH1D("fhistcJetCuts", "Number of c jets after cuts", 15, 0, 15);
		fhistlfJetCuts = new TH1D("fhistlfJetCuts", "Number of lf jets after cuts", 15, 0, 15);
	}

	for(Int_t j = 0; j < 15; j++){
      		fhistInclusiveJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
		if(fIsPythia){
	      		fhistbJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
	      		fhistcJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
	      		fhistlfJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
		}
	}	


	if(fApplyV0Rec){

		  // labels for stages of V0 selection
		  TString categV0[fgkiNCategV0] = {
		    "all"/*0*/,
		    "mass range"/*1*/,
		    "rec. method"/*2*/,
		    "tracks TPC"/*3*/,
		    "track pt"/*4*/,
		    "DCA prim v"/*5*/,
		    "DCA daughters"/*6*/,
		    "CPA"/*7*/,
		    "volume"/*8*/,
		    "track #it{#eta}"/*9*/,
		    "V0 #it{y} & #it{#eta}"/*10*/,
		    "lifetime"/*11*/,
		    "PID"/*12*/,
		    "Arm.-Pod."/*13*/,
		    "cross-cont."/*14*/,
		    "inclusive"/*15*/,
		    "in jet event"/*16*/,
		    "in jet"/*17*/
		  };

		fh2dKshortMassVsPt = new TH2D("fh2dKshortMassVsPt","KShort Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,0.35, 0.65);

		fh2dLamdaMassVsPt = new TH2D("fh2dLamdaMassVsPt","Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);

		fh2dAnLamdaMassVsPt = new TH2D("fh2dAnLamdaMassVsPt","Anti Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);

		if(fIsPythia){
			fh2dKshortMassVsPtReal = new TH2D("fh2dKshortMassVsPtReal","Real KShort Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,0.35, 0.65);

			fh2dLamdaMassVsPtReal = new TH2D("fh2dLamdaMassVsPtReal","Real Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);

			fh2dAnLamdaMassVsPtReal = new TH2D("fh2dAnLamdaMassVsPtReal","Real Anti Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);

			fh2dKshortRecPtVsGenPt = new TH2D("fh2dKshortRecPtVsGenPt","Real KShort Rec Pt Vs Gen Pt;p_{T,Rec} (GeV/c);p_{T,Gen} (GeV/c)",200,0,50,200,0,50);

			fh2dLamdaRecPtVsGenPt = new TH2D("fh2dLamdaRecPtVsGenPt","Real Lamda Rec Pt Vs Gen Pt;p_{T,Rec} (GeV/c);p_{T,Gen} (GeV/c)",200,0,50,200,0,50);

			fh2dAnLamdaRecPtVsGenPt = new TH2D("fh2dAnLamdaRecPtVsGenPt","Real Anti Lamda Rec Pt Vs Gen Pt;p_{T,Rec} (GeV/c);p_{T,Gen} (GeV/c)",200,0,50,200,0,50);
		}

		fhnV0InJetK0s = new THnSparseD("fhnV0InJetK0s", "K0s: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);

		fhnV0InJetLambda = new THnSparseD("fhnV0InJetLambda", "Lambda: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);

		fhnV0InJetALambda = new THnSparseD("fhnV0InJetALambda", "ALambda: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);

		fh1V0CounterCentK0s = new TH1D(Form("fh1V0CounterCentK0s_%d", 0), Form("Number of K0s candidates after cuts, cent %s;cut;counts","0-100%"), fgkiNCategV0, 0, fgkiNCategV0);
	    	fh1V0CounterCentLambda = new TH1D(Form("fh1V0CounterCentLambda_%d", 0), Form("Number of Lambda candidates after cuts, cent %s;cut;counts", "0-100%"), fgkiNCategV0, 0, fgkiNCategV0);
	    	fh1V0CounterCentALambda = new TH1D(Form("fh1V0CounterCentALambda_%d", 0), Form("Number of ALambda candidates after cuts, cent %s;cut;counts", "0-100%"), fgkiNCategV0, 0, fgkiNCategV0);

		fhnV0InJetK0s = new THnSparseD("fhnV0InJetK0s", "K0s: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);

		fhnV0InJetLambda = new THnSparseD("fhnV0InJetLambda", "Lambda: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);

		fhnV0InJetALambda = new THnSparseD("fhnV0InJetALambda", "ALambda: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);

		for(Int_t j = 0; j < fgkiNCategV0; j++)
	    	{
	      		fh1V0CounterCentK0s->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
	      		fh1V0CounterCentLambda->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
	      		fh1V0CounterCentALambda->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
	    	}
	}

	if(fEnableV0GammaRejection) fh1dPhotonPt = new TH1D("fh1dPhotonPt","Photon Mass Vs Pt;p_{T,photon} (GeV/c)",200,0,20);

	if(fApplyV0RejectionAll){
		fh1dKshortPtMC = new TH1D("fh1dKshortPtMC","KShort Pt MC;p_{T} (GeV/c)",200,0,50);

		fh1dLamdaPtMC = new TH1D("fh1dLamdaPtMC","Lamda Pt MC;p_{T} (GeV/c)",200,0,50);

		fh1dAnLamdaPtMC = new TH1D("fh1dAnLamdaPtMC","Anti Lamda Pt MC;p_{T} (GeV/c)",200,0,50);

		fh2dKshortPtVsJetPtMC = new TH2D("fh2dKshortPtVsJetPtMC","KShort Pt Vs Jet Pt MC;p_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/#it{c})",200,0,50,200,0, 200);

		fh2dLamdaPtVsJetPtMC = new TH2D("fh2dLamdaPtVsJetPtMC","Lamda Pt Vs Jet Pt MC;p_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/#it{c})",200,0,50,200,0,200);

		fh2dAnLamdaPtVsJetPtMC = new TH2D("fh2dAnLamdaPtVsJetPtMC","Anti Lamda Pt Vs Jet Pt MC;p_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/#it{c})",200,0,50,200,0,200);
	}


	//Vertex Z before and after
	fh1dVertexZ = new TH1D("fh1dVertexZ","Vertex Z before Event selection;primary vertex z (cm);count",500,-30,30);
	fh1dVertexZAccepted = new TH1D("fh1dVertexZAccepted","Vertex Z after Event selection;primary vertex z (cm);count",500,-30,30);
	//Vertex R before and after
	fh1dVertexR = new TH2D("fh1dVertexR","Vertex R before Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);
	fh1dVertexRAccepted = new TH2D("fh1dVertexRAccepted","Vertex R after Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);

	// Vertex Chi2/NDF vs TPC track multiplicity(ESD tracks)
	fh2dVertexChi2NDFNESDTracks = new TH2D("fh1dVertexChi2NDFNESDTracks","Vertex Chi2/NDF vs # tracks ESD;vertex #chi^{2}/NDF;# tracks esd",200,0,10,500,0,500);
	// AOD tracks accepted
	fh1dTracksAccepeted = new TH1D("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);
	fh1dTracksAccepeted->GetXaxis()->SetBinLabel(1,"total");
	fh1dTracksAccepeted->GetXaxis()->SetBinLabel(2,"accepted");
	fh1dTracksAccepeted->GetXaxis()->SetBinLabel(3,"rejected");
	// Tracks impact parameter histograms
	fh1dTracksImpParXY = new TH1D("fh1dTracksImpParXY","2D Impact Paramter ;impact parameter xy (cm);a.u.",nBins2d,-1.,1.);
	fh1dTracksImpParXYZ = new TH1D("fh1dTracksImpParXYZ","3d imp. parameter ;impact parameter 3d (cm);a.u.",nBins3d,0,1.);
	fh1dTracksImpParXYSignificance = new TH1D("fh1dTracksImpParXYSignificance","2D Impact Paramter ;impact parameter xy significance;a.u.",nBins2dSignificance,-100,100.);
	fh1dTracksImpParXYZSignificance = new TH1D("fh1dTracksImpParXYZSignificance","3d imp. parameter ;impact parameter 3d significance;a.u.",nBins3dSignificance/2,0.,100.);

	fh1dJetRecEtaPhiAccepted = new TH2D("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",200,-1.0,1.0,200,0.,TMath::TwoPi());

	//Jet Mass
	if (fDoJetMass){
		fhistJetMass = new TH2D("fhistJetMass","fJetMass;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
		if(fDoTrackCountingAnalysis){
			fhistJetMassFirst = new TH2D("fhistJetMassFirst","fJetMass N=1;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
			fhistJetMassSecond = new TH2D("fhistJetMassSecond","fJetMass N=2;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
			fhistJetMassThird = new TH2D("fhistJetMassThird","fJetMass N=3;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
		}

		if(fIsPythia){
			fhistJetMass_Unidentified = new TH2D("fhistJetMass_Unidentified","fhistJetMass_Unidentified;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
			fhistJetMass_udsg = new TH2D("fhistJetMass_udsg","fhistJetMass_udsg;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
			fhistJetMass_c = new TH2D("fhistJetMass_c","fhistJetMass_c;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
			fhistJetMass_b = new TH2D("fhistJetMass_b","fhistJetMass_b;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);

			if(fDoTrackCountingAnalysis){

				fhistJetMass_UnidentifiedFirst = new TH2D("fhistJetMass_UnidentifiedFirst","fhistJetMass_UnidentifiedFirst;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_udsgFirst = new TH2D("fhistJetMass_udsgFirst","fhistJetMass_udsgFirst;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_cFirst = new TH2D("fhistJetMass_cFirst","fhistJetMass_cFirst;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_bFirst = new TH2D("fhistJetMass_bFirst","fhistJetMass_bFirst;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);

				fhistJetMass_UnidentifiedSecond = new TH2D("fhistJetMass_UnidentifiedSecond","fhistJetMass_UnidentifiedSecond;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_udsgSecond = new TH2D("fhistJetMass_udsgSecond","fhistJetMass_udsgSecond;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_cSecond = new TH2D("fhistJetMass_cSecond","fhistJetMass_cSecond;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_bSecond = new TH2D("fhistJetMass_bSecond","fhistJetMass_bSecond;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);

				fhistJetMass_UnidentifiedThird = new TH2D("fhistJetMass_UnidentifiedThird","fhistJetMass_UnidentifiedThird;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_udsgThird = new TH2D("fhistJetMass_udsgThird","fhistJetMass_udsgThird;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_cThird = new TH2D("fhistJetMass_cThird","fhistJetMass_cThird;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);
				fhistJetMass_bThird = new TH2D("fhistJetMass_bThird","fhistJetMass_bThird;#it{p}_{T,jet} (GeV/#it{c});Mass (GeV/c^{2})",250,0,250,1000,0,25);

			}
		}

	}

	//Secondary vertex energy fraction
	if(fDoSVEnergyFraction){
		fhistSVEnergyFraction = new TH2D("fhistSVEnergyFraction","SV Energy Fraction;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);

		fhistSVnProngs = new TH2D("fhistSVnProngs", "number of tracks in the SV;#it{p}_{T,jet} (GeV/#it{c});nTracks SV", 250,0,250,5,2,7);

		if(fDoTrackCountingAnalysis){
			fhistSVEnergyFractionFirst = new TH2D("fhistSVEnergyFractionFirst","SV Energy Fraction N=1;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
			fhistSVEnergyFractionSecond = new TH2D("fhistSVEnergyFractionSecond","SV Energy Fraction N=2;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
			fhistSVEnergyFractionThird = new TH2D("fhistSVEnergyFractionThird","SV Energy Fraction N=3;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
		}

		if(fIsPythia){
			fhistSVEnergyFraction_Unidentified = new TH2D("fhistSVEnergyFraction_Unidentified","SV Energy Fraction Undef;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
			fhistSVEnergyFraction_udsg = new TH2D("fhistSVEnergyFraction_udsg","SV Energy Fraction lf-jet;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
			fhistSVEnergyFraction_c    = new TH2D("fhistSVEnergyFraction_c","SV Energy Fraction c-jet;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
			fhistSVEnergyFraction_b   = new TH2D("fhistSVEnergyFraction_b","SV Energy Fraction b-jet;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);

			fhistSVnProngs_Unidentified = new TH2D("fhistSVnProngs_Unidentified", "number of tracks in the SV Unidentified;#it{p}_{T,jet} (GeV/#it{c});nTracks SV", 250,0,250,5,2,7);
			fhistSVnProngs_udsg = new TH2D("fhistSVnProngs_udsg", "number of tracks in the SV lf-jet;#it{p}_{T,jet} (GeV/#it{c});nTracks SV", 250,0,250,5,2,7);

			fhistSVnProngs_b = new TH2D("fhistSVnProngs_b", "number of tracks in the SV b-jet;#it{p}_{T,jet} (GeV/#it{c});nTracks SV", 250,0,250,5,2,7);

			fhistSVnProngs_c = new TH2D("fhistSVnProngs_c", "number of tracks in the SV c-jet;#it{p}_{T,jet} (GeV/#it{c});nTracks SV", 250,0,250,5,2,7);


			if(fDoTrackCountingAnalysis){

				fhistSVEnergyFraction_UnidentifiedFirst = new TH2D("fhistSVEnergyFraction_UnidentifiedFirst","SV Energy Fraction Undef N=1;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_udsgFirst = new TH2D("fhistSVEnergyFraction_udsgFirst","SV Energy Fraction lf-jet N=1;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_cFirst    = new TH2D("fhistSVEnergyFraction_cFirst","SV Energy Fraction c-jet N=1;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_bFirst   = new TH2D("fhistSVEnergyFraction_bFirst","SV Energy Fraction b-jet N=1;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);

				fhistSVEnergyFraction_UnidentifiedSecond = new TH2D("fhistSVEnergyFraction_UnidentifiedSecond","SV Energy Fraction Undef N=2;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_udsgSecond = new TH2D("fhistSVEnergyFraction_udsgSecond","SV Energy Fraction lf-jet N=2;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_cSecond    = new TH2D("fhistSVEnergyFraction_cSecond","SV Energy Fraction c-jet N=2;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_bSecond   = new TH2D("fhistSVEnergyFraction_bSecond","SV Energy Fraction b-jet N=2;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);

				fhistSVEnergyFraction_UnidentifiedThird = new TH2D("fhistSVEnergyFraction_UnidentifiedThird","SV Energy Fraction Undef N=3;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_udsgThird = new TH2D("fhistSVEnergyFraction_udsgThird","SV Energy Fraction lf-jet N=3;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_cThird    = new TH2D("fhistSVEnergyFraction_cThird","SV Energy Fraction c-jet N=3;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);
				fhistSVEnergyFraction_bThird   = new TH2D("fhistSVEnergyFraction_bThird","SV Energy Fraction b-jet N=3;#it{p}_{T,jet} (GeV/#it{c}); f_{E}",250,0,250,500,0,1);

			}

		}

	}

	if(fDoJetProbabilityAnalysis){
		if(!fResolutionFunction[0]){
			fh2dJetSignedImpParXY_Class1 = new TH2D("fh2dJetSignedImpParXY_Class1","Tracks with chi2/NDF>2 IP_{xy};#it{p}_{T} (GeV/#it{c}); IP_{xy} (cm)",200,0,100,2000,-1,1);
			fh2dJetSignedImpParXYSignificance_Class1 = new TH2D("fh2dJetSignedImpParXYSignificance_Class1","Tracks with chi2/NDF>2 sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
			fh2dJetSignedImpParXYZ_Class1 = new TH2D("fh2dJetSignedImpParXYZ_Class1","Tracks with chi2/NDF>2 IP_{xyz};#it{p}_{T} (GeV/#it{c}); IP3D (cm)",200,0,100,2000,-2,2);
			fh2dJetSignedImpParXYZSignificance_Class1 = new TH2D("fh2dJetSignedImpParXYZSignificance_Class1","Tracks with chi2/NDF>2 sIP_{xyz};#it{p}_{T} (GeV/#it{c}); sIP3D",200,0,100,2000,-100,100);

			fh2dJetSignedImpParXY_Class2 = new TH2D("fh2dJetSignedImpParXY_Class2","Tracks with chi2/NDF<2 and 2 ITS hits IP_{xy};#it{p}_{T} (GeV/#it{c}); IP_{xy} (cm)",200,0,100,2000,-1,1);
			fh2dJetSignedImpParXYSignificance_Class2 = new TH2D("fh2dJetSignedImpParXYSignificance_Class2","Tracks with chi2/NDF<2 and 2 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
			fh2dJetSignedImpParXYZ_Class2 = new TH2D("fh2dJetSignedImpParXYZ_Class2","Tracks with chi2/NDF<2 and 2 ITS hits IP_{xyz};#it{p}_{T} (GeV/#it{c}); IP3D (cm)",200,0,100,2000,-2,2);
			fh2dJetSignedImpParXYZSignificance_Class2 = new TH2D("fh2dJetSignedImpParXYZSignificance_Class2","Tracks with chi2/NDF<2 and 2 ITS hits sIP_{xyz};#it{p}_{T} (GeV/#it{c}); sIP3D",200,0,100,2000,-100,100);

			fh2dJetSignedImpParXY_Class3 = new TH2D("fh2dJetSignedImpParXY_Class3","Tracks with chi2/NDF<2 and 3 ITS hits IP_{xy};#it{p}_{T} (GeV/#it{c}); IP_{xy} (cm)",200,0,100,2000,-1,1);
			fh2dJetSignedImpParXYSignificance_Class3 = new TH2D("fh2dJetSignedImpParXYSignificance_Class3","Tracks with chi2/NDF<2 and 3 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
			fh2dJetSignedImpParXYZ_Class3 = new TH2D("fh2dJetSignedImpParXYZ_Class3","Tracks with chi2/NDF<2 and 3 ITS hits IP_{xyz};#it{p}_{T} (GeV/#it{c}); IP3D (cm)",200,0,100,2000,-2,2);
			fh2dJetSignedImpParXYZSignificance_Class3 = new TH2D("fh2dJetSignedImpParXYZSignificance_Class3","Tracks with chi2/NDF<2 and 3 ITS hits sIP_{xyz};#it{p}_{T} (GeV/#it{c}); sIP3D",200,0,100,2000,-100,100);

			fh2dJetSignedImpParXY_Class4 = new TH2D("fh2dJetSignedImpParXY_Class4","Tracks with chi2/NDF<2 and 4 ITS hits IP_{xy};#it{p}_{T} (GeV/#it{c}); IP_{xy} (cm)",200,0,100,2000,-1,1);
			fh2dJetSignedImpParXYSignificance_Class4 = new TH2D("fh2dJetSignedImpParXYSignificance_Class4","Tracks with chi2/NDF<2 and 4 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
			fh2dJetSignedImpParXYZ_Class4 = new TH2D("fh2dJetSignedImpParXYZ_Class4","Tracks with chi2/NDF<2 and 4 ITS hits IP_{xyz};#it{p}_{T} (GeV/#it{c}); IP3D (cm)",200,0,100,2000,-2,2);
			fh2dJetSignedImpParXYZSignificance_Class4 = new TH2D("fh2dJetSignedImpParXYZSignificance_Class4","Tracks with chi2/NDF<2.5 and 4 ITS hits sIP_{xyz};#it{p}_{T} (GeV/#it{c}); sIP3D",200,0,100,2000,-100,100);

			if(fIsPythia){

				fh2dJetSignedImpParXYSignificanceb_Class1 = new TH2D("fh2dJetSignedImpParXYSignificanceb_Class1","Tracks with chi2/NDF>2 sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancec_Class1 = new TH2D("fh2dJetSignedImpParXYSignificancec_Class1","Tracks with chi2/NDF>2 sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancelf_Class1 = new TH2D("fh2dJetSignedImpParXYSignificancelf_Class1","Tracks with chi2/NDF>2 sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);

				fh2dJetSignedImpParXYSignificanceb_Class2 = new TH2D("fh2dJetSignedImpParXYSignificanceb_Class2","Tracks with chi2/NDF<2 and 2 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancec_Class2 = new TH2D("fh2dJetSignedImpParXYSignificancec_Class2","Tracks with chi2/NDF<2 and 2 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancelf_Class2 = new TH2D("fh2dJetSignedImpParXYSignificancelf_Class2","Tracks with chi2/NDF<2 and 2 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);

				fh2dJetSignedImpParXYSignificanceb_Class3 = new TH2D("fh2dJetSignedImpParXYSignificanceb_Class3","Tracks with chi2/NDF<2 and 3 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancec_Class3 = new TH2D("fh2dJetSignedImpParXYSignificancec_Class3","Tracks with chi2/NDF<2 and 3 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancelf_Class3 = new TH2D("fh2dJetSignedImpParXYSignificancelf_Class3","Tracks with chi2/NDF<2 and 3 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);

				fh2dJetSignedImpParXYSignificanceb_Class4 = new TH2D("fh2dJetSignedImpParXYSignificanceb_Class4","Tracks with chi2/NDF<2 and 4 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancec_Class4 = new TH2D("fh2dJetSignedImpParXYSignificancec_Class4","Tracks with chi2/NDF<2 and 4 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);
				fh2dJetSignedImpParXYSignificancelf_Class4 = new TH2D("fh2dJetSignedImpParXYSignificancelf_Class4","Tracks with chi2/NDF<2 and 4 ITS hits sIP_{xy};#it{p}_{T} (GeV/#it{c}); sIP_{xy}",200,0,100,2000,-100,100);

			}

		}else{	
			fhistJetProbability = new TH2D("fhistJetProbability","JetProbability;#it{p}_{T,jet} (GeV/#it{c});JP",250,0,250,1000,0,1);
			fhistJetProbabilityLog = new TH2D("fhistJetProbabilityLog","JetProbability Logarithmic;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

			if(fDoTrackCountingAnalysis){
				fhistJetProbabilityLogFirst = new TH2D("fhistJetProbabilityLogFirst","JetProbability Logarithmic;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				fhistJetProbabilityLogSecond = new TH2D("fhistJetProbabilityLogSecond","JetProbability Logarithmic;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				fhistJetProbabilityLogThird = new TH2D("fhistJetProbabilityLogThird","JetProbability Logarithmic;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
			}
			if(fDoSVAnalysis){
				fhistJetProbabilityLogSVHE = new TH2D("fhistJetProbabilityLogSVHE","JetProbability Logarithmic;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				fhistJetProbabilityLogSVHP = new TH2D("fhistJetProbabilityLogSVHP","JetProbability Logarithmic;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
			}

			if(fIsPythia){
				fhistJetProbability_Unidentified = new TH2D("fhistJetProbability_Unidentified","JetProbability_Unidentified;#it{p}_{T,jet} (GeV/#it{c});JP",250,0,250,1000,0,1);
				fhistJetProbability_udsg = new TH2D("fhistJetProbability_udsg","JetProbability_udsg;#it{p}_{T,jet} (GeV/#it{c});JP",250,0,250,1000,0,1);
				fhistJetProbability_c = new TH2D("fhistJetProbability_c","JetProbability_c;#it{p}_{T,jet} (GeV/#it{c});JP",250,0,250,1000,0,1);
				fhistJetProbability_b = new TH2D("fhistJetProbability_b","JetProbability_b;#it{p}_{T,jet} (GeV/#it{c});JP",250,0,250,1000,0,1);

				fhistJetProbability_UnidentifiedLog = new TH2D("fhistJetProbability_UnidentifiedLog","JetProbability_Unidentified;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				fhistJetProbability_udsgLog = new TH2D("fhistJetProbability_udsgLog","JetProbability_udsg;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				fhistJetProbability_cLog = new TH2D("fhistJetProbability_cLog","JetProbability_c;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				fhistJetProbability_bLog = new TH2D("fhistJetProbability_bLog","JetProbability_b;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

				if(fDoTrackCountingAnalysis){
					fhistJetProbability_UnidentifiedLogFirst = new TH2D("fhistJetProbability_UnidentifiedLogFirst","JetProbability_Unidentified N=1 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_udsgLogFirst = new TH2D("fhistJetProbability_udsgLogFirst","JetProbability_udsg N=1 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_cLogFirst = new TH2D("fhistJetProbability_cLogFirst","JetProbability_c N=1 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_bLogFirst = new TH2D("fhistJetProbability_bLogFirst","JetProbability_b N=1 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					fhistJetProbability_UnidentifiedLogSecond = new TH2D("fhistJetProbability_UnidentifiedLogSecond","JetProbability_Unidentified N=2 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_udsgLogSecond = new TH2D("fhistJetProbability_udsgLogSecond","JetProbability_udsg N=2 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_cLogSecond = new TH2D("fhistJetProbability_cLogSecond","JetProbability_c N=2 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_bLogSecond = new TH2D("fhistJetProbability_bLogSecond","JetProbability_b N=2 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					fhistJetProbability_UnidentifiedLogThird = new TH2D("fhistJetProbability_UnidentifiedLogThird","JetProbability_Unidentified N=3 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_udsgLogThird = new TH2D("fhistJetProbability_udsgLogThird","JetProbability_udsg N=3 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_cLogThird = new TH2D("fhistJetProbability_cLogThird","JetProbability_c N=3 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_bLogThird = new TH2D("fhistJetProbability_bLogThird","JetProbability_b N=3 Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					if(fDoCharmFractions){

					  fhistJetProbability_cLog_D0 = new TH2D("fhistJetProbability_cLog_D0","JetProbability_c D0;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLog_Dp = new TH2D("fhistJetProbability_cLog_Dp","JetProbability_c Dp;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLog_Ds = new TH2D("fhistJetProbability_cLog_Ds","JetProbability_c Ds;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLog_Lc = new TH2D("fhistJetProbability_cLog_Lc","JetProbability_c Lc;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					  fhistJetProbability_cLogFirst_D0 = new TH2D("fhistJetProbability_cLogFirst_D0","JetProbability_c D0;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogFirst_Dp = new TH2D("fhistJetProbability_cLogFirst_Dp","JetProbability_c Dp;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogFirst_Ds = new TH2D("fhistJetProbability_cLogFirst_Ds","JetProbability_c Ds;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogFirst_Lc = new TH2D("fhistJetProbability_cLogFirst_Lc","JetProbability_c Lc;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					  fhistJetProbability_cLogSecond_D0 = new TH2D("fhistJetProbability_cLogSecond_D0","JetProbability_c D0;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogSecond_Dp = new TH2D("fhistJetProbability_cLogSecond_Dp","JetProbability_c Dp;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogSecond_Ds = new TH2D("fhistJetProbability_cLogSecond_Ds","JetProbability_c Ds;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogSecond_Lc = new TH2D("fhistJetProbability_cLogSecond_Lc","JetProbability_c Lc;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					  fhistJetProbability_cLogThird_D0 = new TH2D("fhistJetProbability_cLogThird_D0","JetProbability_c D0;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogThird_Dp = new TH2D("fhistJetProbability_cLogThird_Dp","JetProbability_c Dp;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogThird_Ds = new TH2D("fhistJetProbability_cLogThird_Ds","JetProbability_c Ds;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					  fhistJetProbability_cLogThird_Lc = new TH2D("fhistJetProbability_cLogThird_Lc","JetProbability_c Lc;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					}
				}
				if(fDoSVAnalysis){
					fhistJetProbability_UnidentifiedLogSVHE = new TH2D("fhistJetProbability_UnidentifiedLogSVHE","JetProbability_Unidentified SVHE Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_udsgLogSVHE = new TH2D("fhistJetProbability_udsgLogSVHE","JetProbability_udsg SVHE Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_cLogSVHE = new TH2D("fhistJetProbability_cLogSVHE","JetProbability_c SVHE Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_bLogSVHE = new TH2D("fhistJetProbability_bLogSVHE","JetProbability_b SVHE Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);

					fhistJetProbability_UnidentifiedLogSVHP = new TH2D("fhistJetProbability_UnidentifiedLogSVHP","JetProbability_Unidentified SVHP Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_udsgLogSVHP = new TH2D("fhistJetProbability_udsgLogSVHP","JetProbability_udsg SVHP Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_cLogSVHP = new TH2D("fhistJetProbability_cLogSVHP","JetProbability_c SVHP Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
					fhistJetProbability_bLogSVHP = new TH2D("fhistJetProbability_bLogSVHP","JetProbability_b SVHP Tagged;#it{p}_{T,jet} (GeV/#it{c});-ln(JP)",250,0,250,375,0,30);
				}
			}
		}
	}


	fh1dJetRecPtAcceptedunCorr = new TH1D("fh1dJetRecPtAcceptedunCorr","Rec Jet Pt uncorrected;#it{p}_{T,jet} (GeV/#it{c})",250 ,0, 250);


	f2histRhoVsDeltaPt = new TH2D("f2histRhoVsDeltaPt","Rho Vs Delta Pt;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
	f2histRhoVsDeltaPtFirst = new TH2D("f2histRhoVsDeltaPtFirst","Rho Vs Delta Pt N=1 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
	f2histRhoVsDeltaPtSecond = new TH2D("f2histRhoVsDeltaPtSecond","Rho Vs Delta Pt N=2 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
	f2histRhoVsDeltaPtThird = new TH2D("f2histRhoVsDeltaPtThird","Rho Vs Delta Pt N=3 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);

	if(fDoDeltaPtWithSignal){
		f2histRhoVsDeltaPtWithSignal = new TH2D("f2histRhoVsDeltaPtWithSignal","Rho Vs Delta Pt;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
		f2histRhoVsDeltaPtWithSignalFirst = new TH2D("f2histRhoVsDeltaPtWithSignalFirst","Rho Vs Delta Pt N=1 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
		f2histRhoVsDeltaPtWithSignalSecond = new TH2D("f2histRhoVsDeltaPtWithSignalSecond","Rho Vs Delta Pt N=2 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
		f2histRhoVsDeltaPtWithSignalThird = new TH2D("f2histRhoVsDeltaPtWithSignalThird","Rho Vs Delta Pt N=3 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)",170,-20,150,30,0,30);
	}

        


	if (fIsPythia){

		fh1dTracksImpParXYTruth = new TH1D("fh1dTracksImpParXYTruth","True: 2D Impact Paramter ;impact parameter xy (cm);a.u.",nBins2d,-1.,1.);
		fh1dTracksImpParXYZTruth = new TH1D("fh1dTracksImpParXYZTruth","True: 3d imp. parameter ;impact parameter 3d (cm);a.u.",nBins3d,0,1.);
		fh1dTracksImpParXYResidualTruth  = new TH1D ("fh1dTracksImpParXYResidualTruth","Residual 2D Impact Paramter; #frac{|DCA_{xy}| - |DCA^{Truth}_{xy}|}{#sigma_{xy}} (N#sigma);a.u.",1000,-5,5);
		fh1dTracksImpParXYZResidualTruth  = new TH1D ("fh1dTracksImpParXYZResidualTruth","Residual 3d imp. parameter; #frac{|DCA_{xyz}| - |DCA^{Truth}_{xyz}|}{#sigma_{xyz}} (N#sigma);a.u.",1000,-5,5);
		fh1dJetGenPt = new TH1D("fh1dJetGenPt","generator level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetGenPtUnidentified = new TH1D("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetGenPtudsg = new TH1D("fh1dJetGenPtudsg","generator level udsg jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetGenPtc = new TH1D("fh1dJetGenPtc","generator level c jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetGenPtb = new TH1D("fh1dJetGenPtb","generator level b jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);

		fh2dJetGenPtVsJetRecPt = new TH2D("fh2dJetGenPtVsJetRecPt","detector momentum response;rec pt;gen pt",500,0,250,500,0,250);

		fh2dJetGenPtVsJetRecPtb = new TH2D("fh2dJetGenPtVsJetRecPtb","detector momentum response;rec pt;gen pt",500,0,250,500,0,250);

		fh2dJetGenPtVsJetRecPtc = new TH2D("fh2dJetGenPtVsJetRecPtc","detector momentum response;rec pt;gen pt",500,0,250,500,0,250);

		fh2dJetGenPtVsJetRecPtudsg = new TH2D("fh2dJetGenPtVsJetRecPtudsg","detector momentum response;rec pt;gen pt",500,0,250,500,0,250);

		if(fDoTaggedDRM){

			fh2dJetGenPtVsJetRecPtFirst = new TH2D("fh2dJetGenPtVsJetRecPtFirst","detector momentum response N=1;rec pt;gen pt",500,0,250,500,0,250);
			fh2dJetGenPtVsJetRecPtSecond = new TH2D("fh2dJetGenPtVsJetRecPtSecond","detector momentum response N=2 ;rec pt;gen pt",500,0,250,500,0,250);
			fh2dJetGenPtVsJetRecPtThird = new TH2D("fh2dJetGenPtVsJetRecPtThird","detector momentum response N=3;rec pt;gen pt",500,0,250,500,0,250);

		}

		//Track Counting Analysis
		if(fDoTrackCountingAnalysis){

			fh2dJetSignedImpParXYUnidentified = new TH2D("fh2dJetSignedImpParXYUnidentified","fh2dJetSignedImpParXYUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZUnidentified = new TH2D("fh2dJetSignedImpParXYZUnidentified","fh2dJetSignedImpParXYZUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceUnidentified = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentified","fh2dJetSignedImpParXYSignificanceUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceUnidentified = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentified","fh2dJetSignedImpParXYZSignificanceUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYudsg = new TH2D("fh2dJetSignedImpParXYudsg","fh2dJetSignedImpParXYudsg;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZudsg  = new TH2D("fh2dJetSignedImpParXYZudsg","fh2dJetSignedImpParXYZudsg;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceudsg  = new TH2D("fh2dJetSignedImpParXYSignificanceudsg","fh2dJetSignedImpParXYSignificanceudsg;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceudsg  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsg","fh2dJetSignedImpParXYZSignificanceudsg;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYc= new TH2D("fh2dJetSignedImpParXYc","fh2dJetSignedImpParXYc;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZc  = new TH2D("fh2dJetSignedImpParXYZc","fh2dJetSignedImpParXYZc;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancec  = new TH2D("fh2dJetSignedImpParXYSignificancec","fh2dJetSignedImpParXYSignificancec;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancec  = new TH2D("fh2dJetSignedImpParXYZSignificancec","fh2dJetSignedImpParXYZSignificancec;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYb= new TH2D("fh2dJetSignedImpParXYb","fh2dJetSignedImpParXYb;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZb  = new TH2D("fh2dJetSignedImpParXYZb","fh2dJetSignedImpParXYZb;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceb  = new TH2D("fh2dJetSignedImpParXYSignificanceb","fh2dJetSignedImpParXYSignificanceb;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceb  = new TH2D("fh2dJetSignedImpParXYZSignificanceb","fh2dJetSignedImpParXYZSignificanceb;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
			//N=1
			fh2dJetSignedImpParXYUnidentifiedFirst= new TH2D("fh2dJetSignedImpParXYUnidentifiedFirst","fh2dJetSignedImpParXYUnidentifiedFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYZUnidentifiedFirst","fh2dJetSignedImpParXYZUnidentifiedFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedFirst","fh2dJetSignedImpParXYSignificanceUnidentifiedFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst","fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYudsgFirst = new TH2D("fh2dJetSignedImpParXYudsgFirst","fh2dJetSignedImpParXYudsgFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZudsgFirst  = new TH2D("fh2dJetSignedImpParXYZudsgFirst","fh2dJetSignedImpParXYZudsgFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceudsgFirst  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgFirst","fh2dJetSignedImpParXYSignificanceudsgFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceudsgFirst  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgFirst","fh2dJetSignedImpParXYZSignificanceudsgFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYcFirst= new TH2D("fh2dJetSignedImpParXYcFirst","fh2dJetSignedImpParXYcFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZcFirst  = new TH2D("fh2dJetSignedImpParXYZcFirst","fh2dJetSignedImpParXYZcFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancecFirst  = new TH2D("fh2dJetSignedImpParXYSignificancecFirst","fh2dJetSignedImpParXYSignificancecFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancecFirst  = new TH2D("fh2dJetSignedImpParXYZSignificancecFirst","fh2dJetSignedImpParXYZSignificancecFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYbFirst= new TH2D("fh2dJetSignedImpParXYbFirst","fh2dJetSignedImpParXYbFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZbFirst  = new TH2D("fh2dJetSignedImpParXYZbFirst","fh2dJetSignedImpParXYZbFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancebFirst  = new TH2D("fh2dJetSignedImpParXYSignificancebFirst","fh2dJetSignedImpParXYSignificancebFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancebFirst  = new TH2D("fh2dJetSignedImpParXYZSignificancebFirst","fh2dJetSignedImpParXYZSignificancebFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			//N=2
			fh2dJetSignedImpParXYUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYUnidentifiedSecond","fh2dJetSignedImpParXYUnidentifiedSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYZUnidentifiedSecond","fh2dJetSignedImpParXYZUnidentifiedSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedSecond","fh2dJetSignedImpParXYSignificanceUnidentifiedSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond","fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYudsgSecond = new TH2D("fh2dJetSignedImpParXYudsgSecond","fh2dJetSignedImpParXYudsgSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZudsgSecond  = new TH2D("fh2dJetSignedImpParXYZudsgSecond","fh2dJetSignedImpParXYZudsgSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceudsgSecond  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgSecond","fh2dJetSignedImpParXYSignificanceudsgSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceudsgSecond  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgSecond","fh2dJetSignedImpParXYZSignificanceudsgSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYcSecond= new TH2D("fh2dJetSignedImpParXYcSecond","fh2dJetSignedImpParXYcSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZcSecond  = new TH2D("fh2dJetSignedImpParXYZcSecond","fh2dJetSignedImpParXYZcSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancecSecond  = new TH2D("fh2dJetSignedImpParXYSignificancecSecond","fh2dJetSignedImpParXYSignificancecSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancecSecond  = new TH2D("fh2dJetSignedImpParXYZSignificancecSecond","fh2dJetSignedImpParXYZSignificancecSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYbSecond= new TH2D("fh2dJetSignedImpParXYbSecond","fh2dJetSignedImpParXYbSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZbSecond  = new TH2D("fh2dJetSignedImpParXYZbSecond","fh2dJetSignedImpParXYZbSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancebSecond  = new TH2D("fh2dJetSignedImpParXYSignificancebSecond","fh2dJetSignedImpParXYSignificancebSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancebSecond  = new TH2D("fh2dJetSignedImpParXYZSignificancebSecond","fh2dJetSignedImpParXYZSignificancebSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
			//N=3
			fh2dJetSignedImpParXYUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYUnidentifiedThird","fh2dJetSignedImpParXYUnidentifiedThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYZUnidentifiedThird","fh2dJetSignedImpParXYZUnidentifiedThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedThird","fh2dJetSignedImpParXYSignificanceUnidentifiedThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedThird","fh2dJetSignedImpParXYZSignificanceUnidentifiedThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYudsgThird = new TH2D("fh2dJetSignedImpParXYudsgThird","fh2dJetSignedImpParXYudsgThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZudsgThird  = new TH2D("fh2dJetSignedImpParXYZudsgThird","fh2dJetSignedImpParXYZudsgThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificanceudsgThird  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgThird","fh2dJetSignedImpParXYSignificanceudsgThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificanceudsgThird  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgThird","fh2dJetSignedImpParXYZSignificanceudsgThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYcThird= new TH2D("fh2dJetSignedImpParXYcThird","fh2dJetSignedImpParXYcThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZcThird  = new TH2D("fh2dJetSignedImpParXYZcThird","fh2dJetSignedImpParXYZcThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancecThird  = new TH2D("fh2dJetSignedImpParXYSignificancecThird","fh2dJetSignedImpParXYSignificancecThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancecThird  = new TH2D("fh2dJetSignedImpParXYZSignificancecThird","fh2dJetSignedImpParXYZSignificancecThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

			fh2dJetSignedImpParXYbThird= new TH2D("fh2dJetSignedImpParXYbThird","fh2dJetSignedImpParXYbThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fh2dJetSignedImpParXYZbThird  = new TH2D("fh2dJetSignedImpParXYZbThird","fh2dJetSignedImpParXYZbThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
			fh2dJetSignedImpParXYSignificancebThird  = new TH2D("fh2dJetSignedImpParXYSignificancebThird","fh2dJetSignedImpParXYSignificancebThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			fh2dJetSignedImpParXYZSignificancebThird  = new TH2D("fh2dJetSignedImpParXYZSignificancebThird","fh2dJetSignedImpParXYZSignificancebThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);


			//N=4 Optional
			if(fDoForthIP){		

				fh2dJetSignedImpParXYudsgForth = new TH2D("fh2dJetSignedImpParXYudsgForth","fh2dJetSignedImpParXYudsgForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
				fh2dJetSignedImpParXYSignificanceudsgForth  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgForth","fh2dJetSignedImpParXYSignificanceudsgForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			
				fh2dJetSignedImpParXYcForth= new TH2D("fh2dJetSignedImpParXYcForth","fh2dJetSignedImpParXYcForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
				fh2dJetSignedImpParXYSignificancecForth  = new TH2D("fh2dJetSignedImpParXYSignificancecForth","fh2dJetSignedImpParXYSignificancecForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
			
				fh2dJetSignedImpParXYbForth= new TH2D("fh2dJetSignedImpParXYbForth","fh2dJetSignedImpParXYbForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
				fh2dJetSignedImpParXYSignificancebForth  = new TH2D("fh2dJetSignedImpParXYSignificancebForth","fh2dJetSignedImpParXYSignificancebForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
					
			}
		}
		

	}
	// Jet histograms
	fh1dJetRecPt = new TH1D("fh1dJetRecPt","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
	fh1dJetRecPtAccepted = new TH1D("fh1dJetRecPtAccepted","accepted detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);

	if(fIsPythia){
		fh1dJetRecPtUnidentified = new TH1D("fh1dJetRecPtUnidentified","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtudsg = new TH1D("fh1dJetRecPtudsg","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtc = new TH1D("fh1dJetRecPtc","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtb = new TH1D("fh1dJetRecPtb","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtUnidentifiedAccepted = new TH1D("fh1dJetRecPtUnidentifiedAccepted","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtudsgAccepted = new TH1D("fh1dJetRecPtudsgAccepted","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtcAccepted= new TH1D("fh1dJetRecPtcAccepted","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
		fh1dJetRecPtbAccepted = new TH1D("fh1dJetRecPtbAccepted","detector level jets;#it{p}_{T} (GeV/#it{c}); count",500,0,250);
	}

	//PtRel
	if(fDoPtRelAnalysis){

		fhistPtRelEvents = new TH1D("fhistPtRelEvents","Number of PtRel Events", 1,0,1);
		fhistPtRelVsJetPt = new TH2D("fhistPtRelVsJetPt","Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
		fhistLepIPVsJetPt = new TH2D("fhistLepIPVsJetPt","Electron 2D IP Vs Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fHistMcEopEle = new TH2D("fHistMcEopEle","Real Electrons E/P Vs Track p_{T};p_{T,track} (GeV/c); E/P;a.u.",200,0.,50,200,0.4,2);
		fHistMcEopHad = new TH2D("fHistMcEopHad","Hadrons E/P Vs Track p_{T};p_{T,track} (GeV/c); E/P;a.u.",200,0.,50,200,0.4,2);
		fTPCnsigMcEle = new TH2D("fTPCnsigMcEle","Real Electrons NsigmaTPC Vs Track p_{T};p_{T,track} (GeV/c); N#sigma_{TPC};a.u.",200,0.,50,200,-10,6);
		fTPCnsigMcHad = new TH2D("fTPCnsigMcHad","Hadrons NsigmaTPC Vs Track p_{T};p_{T,track} (GeV/c); N#sigma_{TPC};a.u.",200,0.,50,200,-10,6);

		if(fDoTrackCountingAnalysis){
			fhistPtRelVsJetPtTaggedFirst = new TH2D("fhistPtRelVsJetPtTaggedFirst","N=1 Tagged Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			fhistLepIPVsJetPtTaggedFirst = new TH2D("fhistLepIPVsJetPtTaggedFirst","Electron 2D IP Vs N=1 Tagged Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fhistPtRelVsJetPtTaggedSecond = new TH2D("fhistPtRelVsJetPtTaggedSecond","N=2 Tagged Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			fhistLepIPVsJetPtTaggedSecond  = new TH2D("fhistLepIPVsJetPtTaggedSecond","Electron 2D IP Vs N=2 Tagged Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			fhistPtRelVsJetPtTaggedThird = new TH2D("fhistPtRelVsJetPtTaggedThird","N=3 Tagged Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			fhistLepIPVsJetPtTaggedThird  = new TH2D("fhistLepIPVsJetPtTaggedThird","Electron 2D IP Vs N=3 Tagged Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		}

		if(fIsPythia){
		   fhistPtRelVsJetPtUnidentified = new TH2D("fhistPtRelVsJetPtUnidentified","Unidentified Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
		   fhistPtRelVsJetPtudsg = new TH2D("fhistPtRelVsJetPtudsg","lf-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
		   fhistPtRelVsJetPtc = new TH2D("fhistPtRelVsJetPtc","c-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
		   fhistPtRelVsJetPtb = new TH2D("fhistPtRelVsJetPtb","b-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
		   fhistLepIPVsJetPtUnidentified = new TH2D("fhistLepIPVsJetPtUnidentified","Electron 2D IP Vs Unidentified Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		   fhistLepIPVsJetPtudsg = new TH2D("fhistLepIPVsJetPtudsg","Electron 2D IP Vs lf Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		   fhistLepIPVsJetPtc = new TH2D("fhistLepIPVsJetPtc","Electron 2D IP Vs c-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		   fhistLepIPVsJetPtb = new TH2D("fhistLepIPVsJetPtb","Electron 2D IP Vs b-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);

		   if(fDoTrackCountingAnalysis){
			   //N=1
			   fhistPtRelVsJetPtTaggedUnidentifiedFirst = new TH2D("fhistPtRelVsJetPtTaggedUnidentifiedFirst","Unidentified N=1 Tagged Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedUnidentifiedFirst = new TH2D("fhistLepIPVsJetPtTaggedUnidentifiedFirst","Electron 2D IP Vs Unidentified N=1 Tagged Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedudsgFirst = new TH2D("fhistPtRelVsJetPtTaggedudsgFirst","N=1 Tagged lf-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedudsgFirst = new TH2D("fhistLepIPVsJetPtTaggedudsgFirst","Electron 2D IP Vs N=1 Tagged lf Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedcFirst = new TH2D("fhistPtRelVsJetPtTaggedcFirst","N=1 Tagged c-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedcFirst = new TH2D("fhistLepIPVsJetPtTaggedcFirst","Electron 2D IP Vs N=1 Tagged c-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedbFirst = new TH2D("fhistPtRelVsJetPtTaggedbFirst","N=1 Tagged b-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedbFirst = new TH2D("fhistLepIPVsJetPtTaggedbFirst","Electron 2D IP Vs N=1 Tagged b-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   //N=2
			   fhistPtRelVsJetPtTaggedUnidentifiedSecond = new TH2D("fhistPtRelVsJetPtTaggedUnidentifiedSecond","Unidentified N=2 Tagged Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedUnidentifiedSecond = new TH2D("fhistLepIPVsJetPtTaggedUnidentifiedSecond","Electron 2D IP Vs Unidentified N=2 Tagged Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedudsgSecond = new TH2D("fhistPtRelVsJetPtTaggedudsgSecond","N=2 Tagged lf-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedudsgSecond = new TH2D("fhistLepIPVsJetPtTaggedudsgSecond","Electron 2D IP Vs N=2 Tagged lf Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedcSecond = new TH2D("fhistPtRelVsJetPtTaggedcSecond","N=2 Tagged c-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedcSecond = new TH2D("fhistLepIPVsJetPtTaggedcSecond","Electron 2D IP Vs N=2 Tagged c-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedbSecond = new TH2D("fhistPtRelVsJetPtTaggedbSecond","N=2 Tagged b-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedbSecond = new TH2D("fhistLepIPVsJetPtTaggedbSecond","Electron 2D IP Vs N=2 Tagged b-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   //N=3
			   fhistPtRelVsJetPtTaggedUnidentifiedThird = new TH2D("fhistPtRelVsJetPtTaggedUnidentifiedThird","Unidentified N=3 Tagged Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedUnidentifiedThird = new TH2D("fhistLepIPVsJetPtTaggedUnidentifiedThird","Electron 2D IP Vs Unidentified N=3 Tagged Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedudsgThird = new TH2D("fhistPtRelVsJetPtTaggedudsgThird","N=3 Tagged lf-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedudsgThird = new TH2D("fhistLepIPVsJetPtTaggedudsgThird","Electron 2D IP Vs N=3 Tagged lf Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedcThird = new TH2D("fhistPtRelVsJetPtTaggedcThird","N=3 Tagged c-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedcThird = new TH2D("fhistLepIPVsJetPtTaggedcThird","Electron 2D IP Vs N=3 Tagged c-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
			   fhistPtRelVsJetPtTaggedbThird = new TH2D("fhistPtRelVsJetPtTaggedbThird","N=3 Tagged b-Jet p_{T} Vs Electron p_{T}^{Rel};#it{p}_{T,jet} (GeV/#it{c}); p_{T}^{Rel} (GeV/c);a.u.",500,0.,250,100,0,2);
			   fhistLepIPVsJetPtTaggedbThird = new TH2D("fhistLepIPVsJetPtTaggedbThird","Electron 2D IP Vs N=3 Tagged b-Jet Pt;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		   }
		}
	}

	fh2dJetSignedImpParXY = new TH2D("fh2dJetSignedImpParXY","fh2dJetSignedImpParXY;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);
	fh2dJetSignedImpParXYZ = new TH2D("fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYZ;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
	fh2dJetSignedImpParXYSignificance = new TH2D("fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYSignificance;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
	fh2dJetSignedImpParXYZSignificance = new TH2D("fh2dJetSignedImpParXYZSignificance","fh2dJetSignedImpParXYZSignificance;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

	if(fDoTrackCountingAnalysis){
		//N=1
		fh2dJetSignedImpParXYFirst = new TH2D("fh2dJetSignedImpParXYFirst","fh2dJetSignedImpParXYFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);


		fh2dJetSignedImpParXYZFirst = new TH2D("fh2dJetSignedImpParXYZFirst","fh2dJetSignedImpParXYZFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceFirst = new TH2D("fh2dJetSignedImpParXYSignificanceFirst","fh2dJetSignedImpParXYSignificanceFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceFirst = new TH2D("fh2dJetSignedImpParXYZSignificanceFirst","fh2dJetSignedImpParXYZSignificanceFirst;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=2
		fh2dJetSignedImpParXYSecond = new TH2D("fh2dJetSignedImpParXYSecond","fh2dJetSignedImpParXYSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);


		fh2dJetSignedImpParXYZSecond = new TH2D("fh2dJetSignedImpParXYZSecond","fh2dJetSignedImpParXYZSecond;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceSecond = new TH2D("fh2dJetSignedImpParXYSignificanceSecond","fh2dJetSignedImpParXYSignificanceSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceSecond = new TH2D("fh2dJetSignedImpParXYZSignificanceSecond","fh2dJetSignedImpParXYZSignificanceThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=3

		fh2dJetSignedImpParXYThird = new TH2D("fh2dJetSignedImpParXYThird","fh2dJetSignedImpParXYThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);	

		fh2dJetSignedImpParXYZThird = new TH2D("fh2dJetSignedImpParXYZThird","fh2dJetSignedImpParXYZThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceThird = new TH2D("fh2dJetSignedImpParXYSignificanceThird","fh2dJetSignedImpParXYSignificanceThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceThird = new TH2D("fh2dJetSignedImpParXYZSignificanceThird","fh2dJetSignedImpParXYZSignificanceThird;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		if(fDoForthIP){
			fh2dJetSignedImpParXYForth = new TH2D("fh2dJetSignedImpParXYForth","fh2dJetSignedImpParXYForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.",500,0.,250,nBins2d,-1,1);	
			fh2dJetSignedImpParXYSignificanceForth = new TH2D("fh2dJetSignedImpParXYSignificanceForth","fh2dJetSignedImpParXYSignificanceForth;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		}
	

	}

	if(fDoSVAnalysis){
		fHistSV2Prong = new TH3D("fHistSV2Prong","Secondary vertex 2Prong;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);
		fHistSV3Prong = new TH3D("fHistSV3Prong","Secondary vertex 3Prong;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);

		fHistDispersion2Prong = new TH2D("fHistDispersion2Prong","Secondary vertex Dispersion 2Prong;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",80,0.,250., 100,0,0.5);
		fHistDispersion3Prong = new TH2D("fHistDispersion3Prong","Secondary vertex Dispersion 3Prong;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",80,0.,250., 100,0,0.5);

		if(fIsPythia){

			fHistSV2ProngUnidentified = new TH3D("fHistSV2ProngUnidentified","Secondary vertex 2Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 400,0,80, 100,0,10);
			fHistSV3ProngUnidentified = new TH3D("fHistSV3ProngUnidentified","Secondary vertex 3Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 400,0,80, 100,0,10);

			fHistSV2Prongb = new TH3D("fHistSV2Prongb","Secondary vertex 2Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);
			fHistSV3Prongb = new TH3D("fHistSV3Prongb","Secondary vertex 3Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);

			fHistSV2Prongc = new TH3D("fHistSV2Prongc","Secondary vertex 2Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);
			fHistSV3Prongc = new TH3D("fHistSV3Prongc","Secondary vertex 3Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);

			fHistSV2Pronglf = new TH3D("fHistSV2Pronglf","Secondary vertex 2Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);
			fHistSV3Pronglf = new TH3D("fHistSV3Pronglf","Secondary vertex 3Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)",500,0.,250., 160,0,80, 100,0,10);

			fHistDispersion2ProngUnidentified = new TH2D("fHistDispersion2ProngUnidentified","Secondary vertex Dispersion 2Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);
			fHistDispersion3ProngUnidentified = new TH2D("fHistDispersion3ProngUnidentified","Secondary vertex Dispersion 3Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);

			fHistDispersion2Prongb = new TH2D("fHistDispersion2Prongb","Secondary vertex Dispersion 2Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);
			fHistDispersion3Prongb = new TH2D("fHistDispersion3Prongb","Secondary vertex Dispersion 3Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);

			fHistDispersion2Prongc = new TH2D("fHistDispersion2Prongc","Secondary vertex Dispersion 2Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);
			fHistDispersion3Prongc = new TH2D("fHistDispersion3Prongc","Secondary vertex Dispersion 3Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);

			fHistDispersion2Pronglf = new TH2D("fHistDispersion2Pronglf","Secondary vertex Dispersion 2Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);
			fHistDispersion3Pronglf = new TH2D("fHistDispersion3Pronglf","Secondary vertex Dispersion 3Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion",500,0.,250., 100,0,0.5);
		}
	}


	//Add to output list
	fOutput->Add(fh1dEventRejectionRDHFCuts);
	fOutput->Add(fh1dVertexZ);
	fOutput->Add(fh1dVertexZAccepted);
	fOutput->Add(fh1dVertexR);
	fOutput->Add(fh1dVertexRAccepted);
	fOutput->Add(fh1dTracksAccepeted);
	fOutput->Add(fh1dTracksImpParXY);
	fOutput->Add(fh1dTracksImpParXYZ);
	fOutput->Add(fh1dTracksImpParXYSignificance);
	fOutput->Add(fh1dTracksImpParXYZSignificance);
	fOutput->Add(fh2dVertexChi2NDFNESDTracks);

	fOutput->Add(fhistInclusiveJetCuts);
	if(fIsPythia){
		fOutput->Add(fhistbJetCuts);
		fOutput->Add(fhistcJetCuts);
		fOutput->Add(fhistlfJetCuts);
	}

	if(fApplyV0Rec){
		fOutput->Add(fh2dKshortMassVsPt);
		fOutput->Add(fh2dLamdaMassVsPt);
		fOutput->Add(fh2dAnLamdaMassVsPt);

		if(fIsPythia){
			fOutput->Add(fh2dKshortMassVsPtReal);
			fOutput->Add(fh2dLamdaMassVsPtReal);
			fOutput->Add(fh2dAnLamdaMassVsPtReal);

			fOutput->Add(fh2dKshortRecPtVsGenPt);
			fOutput->Add(fh2dLamdaRecPtVsGenPt);
			fOutput->Add(fh2dAnLamdaRecPtVsGenPt);
		}

		fOutput->Add(fhnV0InJetK0s);
		fOutput->Add(fhnV0InJetLambda);
		fOutput->Add(fhnV0InJetALambda);

		fOutput->Add(fh1V0CounterCentK0s);
		fOutput->Add(fh1V0CounterCentLambda);
		fOutput->Add(fh1V0CounterCentALambda);
	}

	if(fApplyV0RejectionAll){
		fOutput->Add(fh1dKshortPtMC);
		fOutput->Add(fh1dLamdaPtMC);
		fOutput->Add(fh1dAnLamdaPtMC);

		fOutput->Add(fh2dKshortPtVsJetPtMC);
		fOutput->Add(fh2dLamdaPtVsJetPtMC);
		fOutput->Add(fh2dAnLamdaPtVsJetPtMC);
	}

	if(fEnableV0GammaRejection) fOutput->Add(fh1dPhotonPt);

	//PtRel
	if(fDoPtRelAnalysis){
		fOutput->Add(fhistPtRelEvents);
		fOutput->Add(fhistPtRelVsJetPt);
		fOutput->Add(fhistLepIPVsJetPt);
		if(fDoTrackCountingAnalysis){
			fOutput->Add(fhistPtRelVsJetPtTaggedFirst);
			fOutput->Add(fhistLepIPVsJetPtTaggedFirst);
			fOutput->Add(fhistPtRelVsJetPtTaggedSecond);
			fOutput->Add(fhistLepIPVsJetPtTaggedSecond);
			fOutput->Add(fhistPtRelVsJetPtTaggedThird);
			fOutput->Add(fhistLepIPVsJetPtTaggedThird);
		}
		if(fIsPythia){
			fOutput->Add(fhistPtRelVsJetPtUnidentified);
			fOutput->Add(fhistPtRelVsJetPtudsg);
			fOutput->Add(fhistPtRelVsJetPtc);
			fOutput->Add(fhistPtRelVsJetPtb);
			fOutput->Add(fhistLepIPVsJetPtUnidentified);
			fOutput->Add(fhistLepIPVsJetPtudsg);
			fOutput->Add(fhistLepIPVsJetPtc);
			fOutput->Add(fhistLepIPVsJetPtb);
			if(fDoTrackCountingAnalysis){
				fOutput->Add(fhistPtRelVsJetPtTaggedUnidentifiedFirst);
				fOutput->Add(fhistLepIPVsJetPtTaggedUnidentifiedFirst);
				fOutput->Add(fhistPtRelVsJetPtTaggedudsgFirst);
				fOutput->Add(fhistLepIPVsJetPtTaggedudsgFirst);
				fOutput->Add(fhistPtRelVsJetPtTaggedcFirst);
				fOutput->Add(fhistLepIPVsJetPtTaggedcFirst);
				fOutput->Add(fhistPtRelVsJetPtTaggedbFirst);
				fOutput->Add(fhistLepIPVsJetPtTaggedbFirst);
				fOutput->Add(fhistPtRelVsJetPtTaggedUnidentifiedSecond);
				fOutput->Add(fhistLepIPVsJetPtTaggedUnidentifiedSecond);
				fOutput->Add(fhistPtRelVsJetPtTaggedudsgSecond);
				fOutput->Add(fhistLepIPVsJetPtTaggedudsgSecond);
				fOutput->Add(fhistPtRelVsJetPtTaggedcSecond);
				fOutput->Add(fhistLepIPVsJetPtTaggedcSecond);
				fOutput->Add(fhistPtRelVsJetPtTaggedbSecond);
				fOutput->Add(fhistLepIPVsJetPtTaggedbSecond);
				fOutput->Add(fhistPtRelVsJetPtTaggedUnidentifiedThird);
				fOutput->Add(fhistLepIPVsJetPtTaggedUnidentifiedThird);
				fOutput->Add(fhistPtRelVsJetPtTaggedudsgThird);
				fOutput->Add(fhistLepIPVsJetPtTaggedudsgThird);
				fOutput->Add(fhistPtRelVsJetPtTaggedcThird);
				fOutput->Add(fhistLepIPVsJetPtTaggedcThird);
				fOutput->Add(fhistPtRelVsJetPtTaggedbThird);
				fOutput->Add(fhistLepIPVsJetPtTaggedbThird);
			}
		}
		fOutput->Add(fHistMcEopEle);
		fOutput->Add(fHistMcEopHad);
		fOutput->Add(fTPCnsigMcEle);
		fOutput->Add(fTPCnsigMcHad);
	}
	//____

	fOutput->Add(fh1dJetRecPtAcceptedunCorr);

	fOutput->Add(f2histRhoVsDeltaPt);
	fOutput->Add(f2histRhoVsDeltaPtFirst);
	fOutput->Add(f2histRhoVsDeltaPtSecond);
	fOutput->Add(f2histRhoVsDeltaPtThird);

	if(fDoDeltaPtWithSignal){
		fOutput->Add(f2histRhoVsDeltaPtWithSignal);
		fOutput->Add(f2histRhoVsDeltaPtWithSignalFirst);
		fOutput->Add(f2histRhoVsDeltaPtWithSignalSecond);
		fOutput->Add(f2histRhoVsDeltaPtWithSignalThird);
	}

	//Jet Mass
	if(fDoJetMass){
		fOutput->Add(fhistJetMass);
		if(fDoTrackCountingAnalysis){
			fOutput->Add(fhistJetMassFirst);
			fOutput->Add(fhistJetMassSecond);
			fOutput->Add(fhistJetMassThird);
		}
		if(fIsPythia){
			fOutput->Add(fhistJetMass_Unidentified);
			fOutput->Add(fhistJetMass_udsg);
			fOutput->Add(fhistJetMass_c);
			fOutput->Add(fhistJetMass_b);

			if(fDoTrackCountingAnalysis){
				fOutput->Add(fhistJetMass_UnidentifiedFirst);
				fOutput->Add(fhistJetMass_udsgFirst);
				fOutput->Add(fhistJetMass_cFirst);
				fOutput->Add(fhistJetMass_bFirst);
				fOutput->Add(fhistJetMass_UnidentifiedSecond);
				fOutput->Add(fhistJetMass_udsgSecond);
				fOutput->Add(fhistJetMass_cSecond);
				fOutput->Add(fhistJetMass_bSecond);
				fOutput->Add(fhistJetMass_UnidentifiedThird);
				fOutput->Add(fhistJetMass_udsgThird);
				fOutput->Add(fhistJetMass_cThird);
				fOutput->Add(fhistJetMass_bThird);
			}
		}
	}
	
	//Energy Fraction carried by the SV
	if(fDoSVEnergyFraction){
		fOutput->Add(fhistSVEnergyFraction);

		fOutput->Add(fhistSVnProngs);

		if(fDoTrackCountingAnalysis){
			fOutput->Add(fhistSVEnergyFractionFirst);
			fOutput->Add(fhistSVEnergyFractionSecond);
			fOutput->Add(fhistSVEnergyFractionThird);
		}
		if(fIsPythia){
			fOutput->Add(fhistSVEnergyFraction_Unidentified);
			fOutput->Add(fhistSVEnergyFraction_udsg);
			fOutput->Add(fhistSVEnergyFraction_c);
			fOutput->Add(fhistSVEnergyFraction_b);

			fOutput->Add(fhistSVnProngs_Unidentified);
			fOutput->Add(fhistSVnProngs_udsg);
			fOutput->Add(fhistSVnProngs_c);
			fOutput->Add(fhistSVnProngs_b);

			if(fDoTrackCountingAnalysis){
				fOutput->Add(fhistSVEnergyFraction_UnidentifiedFirst);
				fOutput->Add(fhistSVEnergyFraction_udsgFirst);
				fOutput->Add(fhistSVEnergyFraction_cFirst);
				fOutput->Add(fhistSVEnergyFraction_bFirst);
				fOutput->Add(fhistSVEnergyFraction_UnidentifiedSecond);
				fOutput->Add(fhistSVEnergyFraction_udsgSecond);
				fOutput->Add(fhistSVEnergyFraction_cSecond);
				fOutput->Add(fhistSVEnergyFraction_bSecond);
				fOutput->Add(fhistSVEnergyFraction_UnidentifiedThird);
				fOutput->Add(fhistSVEnergyFraction_udsgThird);
				fOutput->Add(fhistSVEnergyFraction_cThird);
				fOutput->Add(fhistSVEnergyFraction_bThird);
			}
		}
	}


	//JetProbability
	if(fDoJetProbabilityAnalysis){
		if(!fResolutionFunction[0]){
			fOutput->Add(fh2dJetSignedImpParXY_Class1);
			fOutput->Add(fh2dJetSignedImpParXYSignificance_Class1);
			fOutput->Add(fh2dJetSignedImpParXYZ_Class1);
			fOutput->Add(fh2dJetSignedImpParXYZSignificance_Class1);
			fOutput->Add(fh2dJetSignedImpParXY_Class2);
			fOutput->Add(fh2dJetSignedImpParXYSignificance_Class2);
			fOutput->Add(fh2dJetSignedImpParXYZ_Class2);
			fOutput->Add(fh2dJetSignedImpParXYZSignificance_Class2);
			fOutput->Add(fh2dJetSignedImpParXY_Class3);
			fOutput->Add(fh2dJetSignedImpParXYSignificance_Class3);
			fOutput->Add(fh2dJetSignedImpParXYZ_Class3);
			fOutput->Add(fh2dJetSignedImpParXYZSignificance_Class3);
			fOutput->Add(fh2dJetSignedImpParXY_Class4);
			fOutput->Add(fh2dJetSignedImpParXYSignificance_Class4);
			fOutput->Add(fh2dJetSignedImpParXYZ_Class4);
			fOutput->Add(fh2dJetSignedImpParXYZSignificance_Class4);
			if(fIsPythia){
				fOutput->Add(fh2dJetSignedImpParXYSignificanceb_Class1);
				fOutput->Add(fh2dJetSignedImpParXYSignificancec_Class1);
				fOutput->Add(fh2dJetSignedImpParXYSignificancelf_Class1);

				fOutput->Add(fh2dJetSignedImpParXYSignificanceb_Class2);
				fOutput->Add(fh2dJetSignedImpParXYSignificancec_Class2);
				fOutput->Add(fh2dJetSignedImpParXYSignificancelf_Class2);

				fOutput->Add(fh2dJetSignedImpParXYSignificanceb_Class3);
				fOutput->Add(fh2dJetSignedImpParXYSignificancec_Class3);
				fOutput->Add(fh2dJetSignedImpParXYSignificancelf_Class3);

				fOutput->Add(fh2dJetSignedImpParXYSignificanceb_Class4);
				fOutput->Add(fh2dJetSignedImpParXYSignificancec_Class4);
				fOutput->Add(fh2dJetSignedImpParXYSignificancelf_Class4);
			}

		}else{
			fOutput->Add(fhistJetProbability);
			fOutput->Add(fhistJetProbabilityLog);

			if(fDoTrackCountingAnalysis){
				fOutput->Add(fhistJetProbabilityLogFirst);
				fOutput->Add(fhistJetProbabilityLogSecond);
				fOutput->Add(fhistJetProbabilityLogThird);
			}
			if(fDoSVAnalysis){
				fOutput->Add(fhistJetProbabilityLogSVHE);
				fOutput->Add(fhistJetProbabilityLogSVHP);
			}

			if(fIsPythia){
				fOutput->Add(fhistJetProbability_Unidentified);
				fOutput->Add(fhistJetProbability_udsg);
				fOutput->Add(fhistJetProbability_c);
				fOutput->Add(fhistJetProbability_b);

				fOutput->Add(fhistJetProbability_UnidentifiedLog);
				fOutput->Add(fhistJetProbability_udsgLog);
				fOutput->Add(fhistJetProbability_cLog);
				if(fDoCharmFractions){
				  fOutput->Add(fhistJetProbability_cLog_D0);
				  fOutput->Add(fhistJetProbability_cLog_Dp);
				  fOutput->Add(fhistJetProbability_cLog_Ds);
				  fOutput->Add(fhistJetProbability_cLog_Lc);
				}
				fOutput->Add(fhistJetProbability_bLog);

				if(fDoTrackCountingAnalysis){
					fOutput->Add(fhistJetProbability_UnidentifiedLogFirst);
					fOutput->Add(fhistJetProbability_udsgLogFirst);
					fOutput->Add(fhistJetProbability_cLogFirst);
					if(fDoCharmFractions){
					  fOutput->Add(fhistJetProbability_cLogFirst_D0);
					  fOutput->Add(fhistJetProbability_cLogFirst_Dp);
					  fOutput->Add(fhistJetProbability_cLogFirst_Ds);
					  fOutput->Add(fhistJetProbability_cLogFirst_Lc);
					}
					fOutput->Add(fhistJetProbability_bLogFirst);

					fOutput->Add(fhistJetProbability_UnidentifiedLogSecond);
					fOutput->Add(fhistJetProbability_udsgLogSecond);
					fOutput->Add(fhistJetProbability_cLogSecond);
					if(fDoCharmFractions){
					  fOutput->Add(fhistJetProbability_cLogSecond_D0);
					  fOutput->Add(fhistJetProbability_cLogSecond_Dp);
					  fOutput->Add(fhistJetProbability_cLogSecond_Ds);
					  fOutput->Add(fhistJetProbability_cLogSecond_Lc);
					}
					fOutput->Add(fhistJetProbability_bLogSecond);

					fOutput->Add(fhistJetProbability_UnidentifiedLogThird);
					fOutput->Add(fhistJetProbability_udsgLogThird);
					fOutput->Add(fhistJetProbability_cLogThird);
					if(fDoCharmFractions){
					  fOutput->Add(fhistJetProbability_cLogThird_D0);
					  fOutput->Add(fhistJetProbability_cLogThird_Dp);
					  fOutput->Add(fhistJetProbability_cLogThird_Ds);
					  fOutput->Add(fhistJetProbability_cLogThird_Lc);
					}
					fOutput->Add(fhistJetProbability_bLogThird);
				}
				if(fDoSVAnalysis){
					fOutput->Add(fhistJetProbability_UnidentifiedLogSVHE);
					fOutput->Add(fhistJetProbability_udsgLogSVHE);
					fOutput->Add(fhistJetProbability_cLogSVHE);
					fOutput->Add(fhistJetProbability_bLogSVHE);

					fOutput->Add(fhistJetProbability_UnidentifiedLogSVHP);
					fOutput->Add(fhistJetProbability_udsgLogSVHP);
					fOutput->Add(fhistJetProbability_cLogSVHP);
					fOutput->Add(fhistJetProbability_bLogSVHP);
				}
			}
		}
	}
	//_______

	if(fIsPythia){

		fOutput->Add(fh1dTracksImpParXYTruth);
		fOutput->Add(fh1dTracksImpParXYZTruth);
		fOutput->Add(fh1dTracksImpParXYResidualTruth);
		fOutput->Add(fh1dTracksImpParXYZResidualTruth);
		fOutput->Add(fh1dJetGenPt);
		fOutput->Add(fh1dJetGenPtUnidentified);
		fOutput->Add(fh1dJetGenPtudsg);
		fOutput->Add(fh1dJetGenPtc);
		fOutput->Add(fh1dJetGenPtb);
		fOutput->Add(fh2dJetGenPtVsJetRecPt);
		fOutput->Add(fh2dJetGenPtVsJetRecPtb);
		fOutput->Add(fh2dJetGenPtVsJetRecPtc);
		fOutput->Add(fh2dJetGenPtVsJetRecPtudsg);
		if(fDoTaggedDRM){

			fOutput->Add(fh2dJetGenPtVsJetRecPtFirst);
			fOutput->Add(fh2dJetGenPtVsJetRecPtSecond);
			fOutput->Add(fh2dJetGenPtVsJetRecPtThird);

		}

		if(fDoTrackCountingAnalysis){
			fOutput->Add(fh2dJetSignedImpParXYUnidentified);
			fOutput->Add(fh2dJetSignedImpParXYZUnidentified);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentified);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentified);
			fOutput->Add(fh2dJetSignedImpParXYudsg);
			fOutput->Add(fh2dJetSignedImpParXYZudsg);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceudsg);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsg);
			fOutput->Add(fh2dJetSignedImpParXYc);
			fOutput->Add(fh2dJetSignedImpParXYZc);
			fOutput->Add(fh2dJetSignedImpParXYSignificancec);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancec);
			fOutput->Add(fh2dJetSignedImpParXYb);
			fOutput->Add(fh2dJetSignedImpParXYZb);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceb);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceb);

			fOutput->Add(fh2dJetSignedImpParXYUnidentifiedFirst);
			fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedFirst);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedFirst);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst);
			fOutput->Add(fh2dJetSignedImpParXYudsgFirst);
			fOutput->Add(fh2dJetSignedImpParXYZudsgFirst);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgFirst);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgFirst);
			fOutput->Add(fh2dJetSignedImpParXYcFirst);
			fOutput->Add(fh2dJetSignedImpParXYZcFirst);
			fOutput->Add(fh2dJetSignedImpParXYSignificancecFirst);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancecFirst);
			fOutput->Add(fh2dJetSignedImpParXYbFirst);
			fOutput->Add(fh2dJetSignedImpParXYZbFirst);
			fOutput->Add(fh2dJetSignedImpParXYSignificancebFirst);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancebFirst);

			fOutput->Add(fh2dJetSignedImpParXYUnidentifiedSecond);
			fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedSecond);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedSecond);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond);
			fOutput->Add(fh2dJetSignedImpParXYudsgSecond);
			fOutput->Add(fh2dJetSignedImpParXYZudsgSecond);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgSecond);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgSecond);
			fOutput->Add(fh2dJetSignedImpParXYcSecond);
			fOutput->Add(fh2dJetSignedImpParXYZcSecond);
			fOutput->Add(fh2dJetSignedImpParXYSignificancecSecond);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancecSecond);
			fOutput->Add(fh2dJetSignedImpParXYbSecond);
			fOutput->Add(fh2dJetSignedImpParXYZbSecond);
			fOutput->Add(fh2dJetSignedImpParXYSignificancebSecond);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancebSecond);

			fOutput->Add(fh2dJetSignedImpParXYUnidentifiedThird);
			fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedThird);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedThird);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedThird);
			fOutput->Add(fh2dJetSignedImpParXYudsgThird);
			fOutput->Add(fh2dJetSignedImpParXYZudsgThird);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgThird);
			fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgThird);
			fOutput->Add(fh2dJetSignedImpParXYcThird);
			fOutput->Add(fh2dJetSignedImpParXYZcThird);
			fOutput->Add(fh2dJetSignedImpParXYSignificancecThird);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancecThird);
			fOutput->Add(fh2dJetSignedImpParXYbThird);
			fOutput->Add(fh2dJetSignedImpParXYZbThird);
			fOutput->Add(fh2dJetSignedImpParXYSignificancebThird);
			fOutput->Add(fh2dJetSignedImpParXYZSignificancebThird);

			if(fDoForthIP){
				fOutput->Add(fh2dJetSignedImpParXYudsgForth);
				fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgForth);
				fOutput->Add(fh2dJetSignedImpParXYcForth);
				fOutput->Add(fh2dJetSignedImpParXYSignificancecForth);
				fOutput->Add(fh2dJetSignedImpParXYbForth);
				fOutput->Add(fh2dJetSignedImpParXYSignificancebForth);
			}
		}


	}
	fOutput->Add(fh1dJetRecPt);
	fOutput->Add(fh1dJetRecPtAccepted);
	fOutput->Add(fh1dJetRecEtaPhiAccepted);

	if(fIsPythia){

		fOutput->Add(fh1dJetRecPtUnidentified);
		fOutput->Add(fh1dJetRecPtudsg);
		fOutput->Add(fh1dJetRecPtc);
		fOutput->Add(fh1dJetRecPtb);
		fOutput->Add(fh1dJetRecPtUnidentifiedAccepted);
		fOutput->Add(fh1dJetRecPtudsgAccepted);
		fOutput->Add(fh1dJetRecPtcAccepted);
		fOutput->Add(fh1dJetRecPtbAccepted);
	}

	fOutput->Add(fh2dJetSignedImpParXY);
	fOutput->Add(fh2dJetSignedImpParXYZ);
	fOutput->Add(fh2dJetSignedImpParXYSignificance);
	fOutput->Add(fh2dJetSignedImpParXYZSignificance);

	if(fDoTrackCountingAnalysis){
		fOutput->Add(fh2dJetSignedImpParXYFirst);
		fOutput->Add(fh2dJetSignedImpParXYZFirst);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceFirst);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceFirst);
		fOutput->Add(fh2dJetSignedImpParXYSecond);
		fOutput->Add(fh2dJetSignedImpParXYZSecond);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceSecond);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceSecond);
		fOutput->Add(fh2dJetSignedImpParXYThird);
		fOutput->Add(fh2dJetSignedImpParXYZThird);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceThird);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceThird);
		if(fDoForthIP){
			fOutput->Add(fh2dJetSignedImpParXYForth);
			fOutput->Add(fh2dJetSignedImpParXYSignificanceForth);
		}
	}

	if(fDoSVAnalysis){
		fOutput->Add(fHistDispersion2Prong);
		fOutput->Add(fHistSV2Prong);

		if(fIsPythia){
			fOutput->Add(fHistDispersion2ProngUnidentified);
			fOutput->Add(fHistDispersion2Prongb);
			fOutput->Add(fHistDispersion2Prongc);
			fOutput->Add(fHistDispersion2Pronglf);

			fOutput->Add(fHistSV2ProngUnidentified);
			fOutput->Add(fHistSV2Prongb);
			fOutput->Add(fHistSV2Prongc);
			fOutput->Add(fHistSV2Pronglf);
		}


		fOutput->Add(fHistDispersion3Prong);
		fOutput->Add(fHistSV3Prong);

		if(fIsPythia){
			fOutput->Add(fHistDispersion3ProngUnidentified);
			fOutput->Add(fHistDispersion3Prongb);
			fOutput->Add(fHistDispersion3Prongc);
			fOutput->Add(fHistDispersion3Pronglf);

			fOutput->Add(fHistSV3ProngUnidentified);
			fOutput->Add(fHistSV3Prongb);
			fOutput->Add(fHistSV3Prongc);
			fOutput->Add(fHistSV3Pronglf);
		}

	}

	TIter next(fOutput);
	while (TObject *obj = next.Next()){
		if(obj->IsA() == TH1D::Class() || obj->IsA() == TH2D::Class()){
			((TH1*)obj)->Sumw2();
		}
	}
	PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}
// ######################################################################################## Calculate impact parameters
Bool_t AliAnalysisTaskBJetTC::CalculateTrackImpactParameter(AliAODTrack * track,double *impar, double * cov)
{
	AliAODVertex *vtxAODNew=0x0;
	AliESDVertex *vtxESDNew =0x0;
	Bool_t recalculate = kFALSE;
	if( fPrimaryVertex->GetNContributors() < 30){
		recalculate=kTRUE;
		Int_t skipped[1] = {-1};
		Int_t id = (Int_t)track->GetID();
		if(id<0) return kFALSE;
		skipped[0] = id;
		fVertexer->SetSkipTracks(1,skipped);
		vtxESDNew = fVertexer->FindPrimaryVertex(fAODIn);
		if(!vtxESDNew) return kFALSE;
		if(vtxESDNew->GetNContributors()<=0) {
			delete vtxESDNew; vtxESDNew=NULL;
			return kFALSE;
		}
		// convert to AliAODVertex
		Double_t pos[3],cova[6],chi2perNDF;
		vtxESDNew->GetXYZ(pos); // position
		vtxESDNew->GetCovMatrix(cova); //covariance matrix
		chi2perNDF = vtxESDNew->GetChi2toNDF();
		delete vtxESDNew; vtxESDNew=NULL;
		vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);
	}
	// Calculate Impact Parameters
	AliExternalTrackParam etp; etp.CopyFromVTrack(track);
	if(etp.PropagateToDCA(vtxAODNew,fAODIn->GetMagneticField(),3.,impar,cov))
	{
		if(recalculate)
			delete vtxAODNew;
		return kTRUE;
	}
	else{
		if(recalculate)
			delete vtxAODNew;
		return kFALSE;

	}
}
// ######################################################################################## Calculate impact parameters based on MC event vertex and MC particle information (no special mass treatment)
Bool_t AliAnalysisTaskBJetTC::CalculateTrackImpactParameterTruth(AliAODTrack * track,double *impar, double * cov)
{
	AliAODMCParticle *pMC = 0x0;
	AliAODMCHeader* mcheader = dynamic_cast<AliAODMCHeader*>(fAODIn->FindListObject(AliAODMCHeader::StdBranchName()));
	if (!mcheader) return kFALSE;

	if (!fMCArray) return kFALSE;
	if(track->GetLabel()>-1)
		pMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(track->GetLabel()));
	if (!pMC) return kFALSE;

	Double_t pos[3]={0,0,0};
	mcheader->GetVertex(pos);

	Double_t cova[6]={0,0,0,0,0,0};
	Double_t chi2perNDF =0;
	AliAODVertex *vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);
	double xpart[3] = {0,0,0};
	pMC->XvYvZv(xpart);
	double ppart[3] = {0,0,0};
	pMC->PxPyPz(ppart);
	double cv[21]={0} ;
	AliExternalTrackParam trackparam(xpart,ppart,cv,(TMath::Sign((Short_t)1,(Short_t)pMC->Charge())));
	if(trackparam.PropagateToDCA(vtxAODNew,fAODIn->GetMagneticField(),3.,impar,cov))
	{
		pMC=NULL; delete pMC;
		delete vtxAODNew;
		return kTRUE;
	}
	else{
		pMC=NULL; delete pMC;
		delete vtxAODNew;
		return kFALSE;

	}

	return kFALSE;
}

// ######################################################################################## Calculate signed  impact parameters

Bool_t AliAnalysisTaskBJetTC::CalculateJetSignedTrackImpactParameter(AliAODTrack * track,AliEmcalJet * jet ,double *impar, double * cov, double &sign, double &dcajetrack, double &lineardecaylength){

	Int_t skipped[1] = {-1};
	Int_t id = (Int_t)track->GetID();
	if(id<0) return kFALSE;
	skipped[0] = id;
	fVertexer->SetSkipTracks(1,skipped);
	AliESDVertex *vtxESDNew = fVertexer->FindPrimaryVertex(fAODIn);
	if(!vtxESDNew) return kFALSE;
	if(vtxESDNew->GetNContributors()<=0) {
		vtxESDNew=NULL; delete vtxESDNew; 
		return kFALSE;
	}
	// convert to AliAODVertex
	Double_t pos[3],cova[6],chi2perNDF;
	vtxESDNew->GetXYZ(pos); // position
	vtxESDNew->GetCovMatrix(cova); //covariance matrix
	chi2perNDF = vtxESDNew->GetChi2toNDF();

	// Calculate Impact Parameters
	AliExternalTrackParam etp; etp.CopyFromVTrack(track);
	if(etp.PropagateToDCA(vtxESDNew,fAODIn->GetMagneticField(),3.,impar,cov))
	{
		//Calculate Sign
		Double_t posdcatrack[3]= {0.,0.,0.};
		etp.GetXYZ(posdcatrack);
		Double_t ipvector3[3] = { posdcatrack[0] - pos[0], posdcatrack[1] - pos[1], posdcatrack[2] - pos[2] };
		sign =TMath::Sign(1.,ipvector3[0]*jet->Px() +ipvector3[1]*jet->Py()+ipvector3[2]*jet->Pz() );

		// Calculate decay legnth and track jet DCA against new vertex
		Double_t bpos[3] = { 0,0,0 };
		vtxESDNew->GetXYZ(bpos);
		Double_t bpxpypz[3] = { jet->Px(), jet->Py(), jet->Pz() };
		Double_t bcv[21] = { 0 };
		AliExternalTrackParam bjetparam(bpos, bpxpypz, bcv, (Short_t)0);
		Double_t xa = 0., xb = 0.;
		if(!fDoImprovedDCACut){
			bjetparam.GetDCA(&etp, fAODIn->GetMagneticField(), xa, xb);
		}else{
			fDecayVertex->GetDCAV0Dau(&bjetparam, &etp, xa, xb, fAODIn->GetMagneticField() );
		}
		Double_t xyz[3] = { 0., 0., 0. };
		Double_t xyzb[3] = { 0., 0., 0. };
		bjetparam.GetXYZAt(xa, fAODIn->GetMagneticField(), xyz);
		etp.GetXYZAt(xb, fAODIn->GetMagneticField(), xyzb);
		double  bdecaylength =
				TMath::Sqrt(
						(bpos[0] - xyzb[0]) * (bpos[0] - xyzb[0]) +
						(bpos[1] - xyzb[1]) * (bpos[1] - xyzb[1]) +
						(bpos[2] - xyzb[2]) * (bpos[2] - xyzb[2]));
		dcajetrack =
				TMath::Sqrt(
						(xyzb[0] - xyz[0]) * (xyzb[0] - xyz[0]) +
						(xyzb[1] - xyz[1]) * (xyzb[1] - xyz[1]) +
						(xyzb[2] - xyz[2]) * (xyzb[2] - xyz[2]));
		if(bdecaylength>0) lineardecaylength=bdecaylength;
		delete vtxESDNew;
		return kTRUE;
	}
	else{
		delete vtxESDNew;
		return kFALSE;

	}
}
// ######################################################################################## Post-process ImpPar
Double_t AliAnalysisTaskBJetTC::GetValImpactParameter(TTypeImpPar type,double *impar, double * cov)
{
	double result =-999;
	double dFdx = 0;
	double dFdy = 0;

	switch(type){
	case kXY:
		result = impar[0];
		break;
	case kXYSig:
		result = impar[0]/TMath::Sqrt(cov[0]);
		break;
	case kXYZ:
		result = TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
		break;
	case kXYZSig:
		result =  TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
		dFdx = 2*impar[0]/result;
		dFdy = 2*impar[1]/result;
		result /=TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
		break;
	case kXYZSigmaOnly:
		dFdx = 2*impar[0]/result;
		dFdy = 2*impar[1]/result;
		result =TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
		break;

	default:
		break;
	}
	return result;
}
// ########################################################################################Track Selection
Bool_t AliAnalysisTaskBJetTC::IsTrackAccepted(AliAODTrack* track){
        if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4)))return kFALSE;
	if(track->Pt() < fTCMinTrackPt)return kFALSE; //0.5
	if(TMath::Abs( track->Eta() ) >0.9) return kFALSE;
	ULong_t status = track->GetStatus();
	if(!(status & AliAODTrack::kTPCrefit))return kFALSE;
	if(!(status & AliAODTrack::kITSrefit))return kFALSE;
	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) return kFALSE;
	
	if(track->GetNcls(1)<fTCMinClusTPC) return kFALSE;//80
	Int_t SPDSSDHits = track->HasPointOnITSLayer(0) + track->HasPointOnITSLayer(1) + track->HasPointOnITSLayer(4) + track->HasPointOnITSLayer(5);
	if(SPDSSDHits<fTCMinHitsITS) return kFALSE;//2 in case of FAST and CENT_woSDD and 3 if there was an SDD
	if(track->Chi2perNDF()>=fTCMaxChi2pNDF) return kFALSE;//5
	return kTRUE;
}
// ########################################################################################Track Selection
Bool_t AliAnalysisTaskBJetTC::IsTrackAcceptedBJetCuts(AliAODTrack* track, Int_t jetFlavour){

	Int_t iCutIndex = 0; // indicator of current selection step

	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

        if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4)))return kFALSE;
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	if(track->Pt() < fTCMinTrackPt)return kFALSE; //0.5
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	if(TMath::Abs( track->Eta() ) >0.9) return kFALSE;
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	ULong_t status = track->GetStatus();
	if(!(status & AliAODTrack::kTPCrefit))return kFALSE;
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	if(!(status & AliAODTrack::kITSrefit))return kFALSE;
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) return kFALSE;
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;
	
	if(track->GetNcls(1)<fTCMinClusTPC) return kFALSE;//80
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	Int_t SPDSSDHits = track->HasPointOnITSLayer(0) + track->HasPointOnITSLayer(1) + track->HasPointOnITSLayer(4) + track->HasPointOnITSLayer(5);
	if(SPDSSDHits<fTCMinHitsITS) return kFALSE;//2
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	if(track->Chi2perNDF()>=fTCMaxChi2pNDF) return kFALSE;//5
	FillCandidateJet(iCutIndex, jetFlavour);
	iCutIndex++;

	return kTRUE;
}
//____________________________________________________
void AliAnalysisTaskBJetTC::FillCandidateJet(Int_t CutIndex, Int_t JetFlavor){

	fhistInclusiveJetCuts->Fill(CutIndex);
	if(fIsPythia){
		if(JetFlavor==3)fhistbJetCuts->Fill(CutIndex);
		if(JetFlavor==2)fhistcJetCuts->Fill(CutIndex);
		if(JetFlavor==1)fhistlfJetCuts->Fill(CutIndex);
	}

}
// ########################################################################################Track Selection
Bool_t AliAnalysisTaskBJetTC::IsTrackAcceptedFidutial(AliAODTrack* track){

	if(track->Pt() < 0.15)return kFALSE;
	if(TMath::Abs( track->Eta() ) >0.9) return kFALSE;

	return kTRUE;
}
// ########################################################################################Track Selection
Bool_t AliAnalysisTaskBJetTC::IsTrackAcceptedQuality(AliAODTrack* track ,AliEmcalJet* Jet, Int_t &QualityClass, double *imp, double * cov, double &sign){

  if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4))) return kFALSE;
  if(track->Pt() < fTCMinTrackPt)return kFALSE;//0.5
  if(TMath::Abs( track->Eta() ) >0.9) return kFALSE;
  ULong_t status = track->GetStatus();
  if(!(status & AliAODTrack::kTPCrefit))return kFALSE;
  if(!(status & AliAODTrack::kITSrefit))return kFALSE;
  if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) return kFALSE;
  if(track->GetNcls(1)<fTCMinClusTPC) return kFALSE;//80
  Int_t SPDSSDHits = track->HasPointOnITSLayer(0) + track->HasPointOnITSLayer(1) + track->HasPointOnITSLayer(4) + track->HasPointOnITSLayer(5);
  if(SPDSSDHits<fTCMinHitsITS) return kFALSE;//2
  //if(track->GetNcls(0)<fTCMinHitsITS) return kFALSE;//2
  if(track->Chi2perNDF()>=fTCMaxChi2pNDF) return kFALSE;//5

  double dcaTrackJet =0,lindeclen =0 ;

  if(!CalculateJetSignedTrackImpactParameter(track,Jet,imp,cov,sign,dcaTrackJet,lindeclen)) return kFALSE;

  if(abs(imp[0])>fTCMaxIPxy) return kFALSE;//1cm
  if(abs(imp[1])>fTCMaxIPz) return kFALSE;//5cm
  if(lindeclen > fTCMaxDecayLength) return kFALSE;//5cm
  if (dcaTrackJet > fTCMaxDCATrackJet) return kFALSE;//0.07cm

  if(track->Chi2perNDF()>=2) QualityClass=0;

  else if(track->Chi2perNDF()<2){
	if(track->Pt() < 2.){
		if(SPDSSDHits==2) QualityClass=1;
		else if(SPDSSDHits==3) QualityClass=3;
		else if(SPDSSDHits==4) QualityClass=5;
	}else if(track->Pt() >= 2.){
		if(SPDSSDHits==2) QualityClass=2;
		else if(SPDSSDHits==3) QualityClass=4;
		else if(SPDSSDHits==4) QualityClass=6;
	}
   }
  
  return kTRUE;
}
// ######################################################################################## Jet matching 1/4
Bool_t AliAnalysisTaskBJetTC::MatchJetsGeometricDefault()
{
	double matchingpar1 =0.25;
	double matchingpar2 =0.25;
	if (!fJetContainerData || !fJetContainerData->GetArray() || !fJetContainerMC || !fJetContainerMC->GetArray()) return kFALSE;
	DoJetLoop();
	AliEmcalJet* jet1 = 0;
	AliEmcalJet *jet2 = 0;

	fJetContainerData->ResetCurrentID();
	while ((jet1 = fJetContainerData->GetNextJet())) {
		jet2 = jet1->ClosestJet();
		if (!jet2) continue;
		if (jet2->ClosestJet() != jet1) continue;
		if (jet1->ClosestJetDistance() > matchingpar1 || jet2->ClosestJetDistance() > matchingpar2) continue;
		// Matched jet found
		jet1->SetMatchedToClosest(1);
		jet2->SetMatchedToClosest(1);
	}
	jet1=NULL; delete jet1;
	jet2=NULL; delete jet2;
	return kTRUE;
}
// ######################################################################################## Jet matching 2/4
void AliAnalysisTaskBJetTC::DoJetLoop()
{
	// Do the jet loop.
	double minjetpt =1.;

	AliEmcalJet* jet1 = 0;
	AliEmcalJet* jet2 = 0;
	fJetContainerMC->ResetCurrentID();
	while ((jet2 = fJetContainerMC->GetNextJet())) jet2->ResetMatching();
	fJetContainerData->ResetCurrentID();
	while ((jet1 = fJetContainerData->GetNextJet())) {
		jet1->ResetMatching();
		if (jet1->MCPt() < minjetpt) continue;
		fJetContainerMC->ResetCurrentID();
		while ((jet2 = fJetContainerMC->GetNextJet())) {
			SetMatchingLevel(jet1, jet2, 1);
		} // jet2 loop
	} // jet1 loop

	jet1=NULL; jet2=NULL;
	delete jet1; delete jet2;
}
// ######################################################################################## Jet matching 3/4
void AliAnalysisTaskBJetTC::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, int matching)
{
	Double_t d1 = -1;
	Double_t d2 = -1;

	switch (matching) {
	case 1:
		GetGeometricalMatchingLevel(jet1,jet2,d1);
		d2 = d1;
		break;
	default:
		break;
	}
	if (d1 >= 0) {

		if (d1 < jet1->ClosestJetDistance()) {
			jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
			jet1->SetClosestJet(jet2, d1);
		}
		else if (d1 < jet1->SecondClosestJetDistance()) {
			jet1->SetSecondClosestJet(jet2, d1);
		}
	}
	if (d2 >= 0) {
		if (d2 < jet2->ClosestJetDistance()) {
			jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
			jet2->SetClosestJet(jet1, d2);
		}
		else if (d2 < jet2->SecondClosestJetDistance()) {
			jet2->SetSecondClosestJet(jet1, d2);
		}
	}
}

// ######################################################################################## Jet matching 4/4
void AliAnalysisTaskBJetTC::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
	Double_t deta = jet2->Eta() - jet1->Eta();
	Double_t dphi = jet2->Phi() - jet1->Phi();
	dphi = TVector2::Phi_mpi_pi(dphi);
	d = TMath::Sqrt(deta * deta + dphi * dphi);
}

AliAODMCParticle* AliAnalysisTaskBJetTC::GetMCTrack( const AliAODTrack* _track)
{
	//
	// return MC track
	//
	if(!fIsPythia) return NULL;
	if(!fMCArray) { AliError("No fMCArray"); return NULL;}
	Int_t nStack = fMCArray->GetEntriesFast();
	Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
	if(label > nStack) return NULL;
	AliAODMCParticle *mctrack =  dynamic_cast<AliAODMCParticle *>(fMCArray->At(label));
	return mctrack;
}

//
Bool_t AliAnalysisTaskBJetTC::IsV0PhotonFromBeamPipeDaughter(const AliAODTrack* track)
{
	if(!track)return kFALSE;
	int posid = -1;
	int negid = -1;
	int trackid = -1;

	fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

	AliAODConversionPhoton* PhotonCandidate = 0;
	// Loop over Photon Candidates allocated by ReaderV1
	for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){

		    PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
		    if(!PhotonCandidate) continue;

		    posid = PhotonCandidate->GetLabel1();
		    negid = PhotonCandidate->GetLabel2();
		    trackid = track->GetID();
		    
		    if(posid == trackid || negid == trackid) { return kTRUE;}

	}


	if(fIsPythia > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){

		    fV0Reader->RelabelAODs(kTRUE);
	}


	return kFALSE;
}
//////////////////////////////////////////////////
Bool_t AliAnalysisTaskBJetTC::IsV0Daughter(const AliAODTrack* track)
{
	if(!track)return kFALSE;
	AliAODv0* v0aod = 0x0;
	int posid = -1;
	int negid = -1;
	int trackid = -1;


	if(!fV0CandidateArray) {cout<<"No V0 Candidates \n"; return kFALSE;}
	
	for(int i = 0; i < fV0CandidateArray->GetEntriesFast(); ++i) {

		v0aod = dynamic_cast<AliAODv0*>(fV0CandidateArray->At(i));
		if(!v0aod) continue;

		posid = v0aod->GetPosID();
		negid = v0aod->GetNegID();
		trackid = track->GetID();

		if(posid == trackid || negid == trackid) {return kTRUE; }
	}



	return kFALSE;
}
//=============================================================================
Bool_t AliAnalysisTaskBJetTC::SelectV0CandidateVIT()
{
  AliAODv0* v0 = 0; // pointer to V0 candidates
  Double_t dMassV0K0s = 0; // invariant mass of the K0s candidate
  Double_t dMassV0Lambda = 0; // invariant mass of the Lambda candidate
  Double_t dMassV0ALambda = 0; // invariant mass of the Lambda candidate
  Int_t iNV0CandTot = 0; // counter of all V0 candidates at the beginning
  Int_t iNV0CandK0s = 0; // counter of K0s candidates at the end
  Int_t iNV0CandLambda = 0; // counter of Lambda candidates at the end
  Int_t iNV0CandALambda = 0; // counter of Lambda candidates at the end

  fV0CandidateArray->Delete();//Reset the TClonesArray



// axis: K0S invariant mass
const Double_t fgkdMassK0sMin = 0.35; // [GeV/c^2]
const Double_t fgkdMassK0sMax = 0.65; // [GeV/c^2]

// axis: Lambda invariant mass
const Double_t fgkdMassLambdaMin = 1.05; // [GeV/c^2]
const Double_t fgkdMassLambdaMax = 1.25; // [GeV/c^2]

  Bool_t bPrintCuts = 0; // print out which cuts are applied

  // Other cuts
  Double_t dNSigmaMassMax = 3.; // [sigma m] max difference between candidate mass and real particle mass (used only for mass peak method of signal extraction)
  Double_t dDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)

  // Mean lifetime
  Double_t dCTauK0s = 2.6844; // [cm] c*tau of K0S
  Double_t dCTauLambda = 7.89; // [cm] c*tau of Lambda

  // particle masses from PDG
  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

  // PDG codes of used particles
  Int_t iPdgCodePion = 211;
  Int_t iPdgCodeProton = 2212;
  Int_t iPdgCodeK0s = 310;
  Int_t iPdgCodeLambda = 3122;

  // Loading primary vertex info
  Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
  fPrimaryVertex->GetXYZ(dPrimVtxPos);

  Bool_t fbIsPbPb = kFALSE;

  Int_t iCentIndex=0;

  Int_t iNV0s = fAODIn->GetNumberOfV0s(); // get the number of V0 candidates
  if(!iNV0s)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No V0s found in event");
  }

  for(Int_t iV0 = 0; iV0 < iNV0s; iV0++)
  {
    v0 = fAODIn->GetV0(iV0); // get next candidate from the list in AOD
    if(!v0)
      continue;

    iNV0CandTot++;

    // Initialization of status indicators
    Bool_t bIsCandidateK0s = kTRUE; // candidate for K0s
    Bool_t bIsCandidateLambda = kTRUE; // candidate for Lambda
    Bool_t bIsCandidateALambda = kTRUE; // candidate for anti-Lambda
    Bool_t bIsInPeakK0s = kFALSE; // candidate within the K0s mass peak
    Bool_t bIsInPeakLambda = kFALSE; // candidate within the Lambda mass peak
    Bool_t bIsInPeakALambda = kFALSE; // candidate within the anti-Lambda mass peak
    Bool_t bIsInConeJet = kFALSE; // candidate within the jet cones
    Bool_t bIsInConePerp = kFALSE; // candidate within a perpendicular cone
    Bool_t bIsInConeRnd = kFALSE; // candidate within the random cone
    Bool_t bIsInConeMed = kFALSE; // candidate within the median-cluster cone
    Bool_t bIsOutsideCones = kFALSE; // candidate outside excluded cones

    // Invariant mass calculation
    dMassV0K0s = v0->MassK0Short();
    dMassV0Lambda = v0->MassLambda();
    dMassV0ALambda = v0->MassAntiLambda();

    Int_t iCutIndex = 0; // indicator of current selection step
    // 0
    // All V0 candidates
    FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    Double_t dPtV0 = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0

    // Sigma of the mass peak window
    Double_t dMassPeakWindowK0s = dNSigmaMassMax * ( 0.0044 + 0.0004 * (dPtV0 - 1.) );
    Double_t dMassPeakWindowLambda = dNSigmaMassMax * ( 0.0023 + 0.00034 * (dPtV0 - 1.) );
    if(!fbIsPbPb) // p-p
    {
      dMassPeakWindowK0s = 0.010; // LF p-p
      dMassPeakWindowLambda = 0.005; // LF p-p
    }

    // Invariant mass peak selection
    if(TMath::Abs(dMassV0K0s - dMassPDGK0s) < dMassPeakWindowK0s)
      bIsInPeakK0s = kTRUE;
    if(TMath::Abs(dMassV0Lambda - dMassPDGLambda) < dMassPeakWindowLambda)
      bIsInPeakLambda = kTRUE;
    if(TMath::Abs(dMassV0ALambda - dMassPDGLambda) < dMassPeakWindowLambda)
      bIsInPeakALambda = kTRUE;


    // Skip candidates outside the histogram range
    if((dMassV0K0s < fgkdMassK0sMin) || (dMassV0K0s >= fgkdMassK0sMax))
      bIsCandidateK0s = kFALSE;
    if((dMassV0Lambda < fgkdMassLambdaMin) || (dMassV0Lambda >= fgkdMassLambdaMax))
      bIsCandidateLambda = kFALSE;
    if((dMassV0ALambda < fgkdMassLambdaMin) || (dMassV0ALambda >= fgkdMassLambdaMax))
      bIsCandidateALambda = kFALSE;
    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

    // Retrieving all relevant properties of the V0 candidate
    Bool_t bOnFlyStatus = v0->GetOnFlyStatus(); // online (on fly) reconstructed vs offline reconstructed
    const AliAODTrack* trackPos = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
    const AliAODTrack* trackNeg = (AliAODTrack*)v0->GetDaughter(1); // negative daughter track
    Double_t dPtDaughterPos = trackPos->Pt(); // transverse momentum of a daughter track calculated as if primary, != v0->PtProng(0)
    Double_t dPtDaughterNeg = trackNeg->Pt(); // != v0->PtProng(1)
    Double_t dNRowsPos = trackPos->GetTPCClusterInfo(2, 1); // crossed TPC pad rows of a daughter track
    Double_t dNRowsNeg = trackNeg->GetTPCClusterInfo(2, 1);
    Double_t dFindablePos = Double_t(trackPos->GetTPCNclsF()); // Findable clusters
    Double_t dFindableNeg = Double_t(trackNeg->GetTPCNclsF());
    Double_t dDCAToPrimVtxPos = TMath::Abs(v0->DcaPosToPrimVertex()); // dca of a daughter to the primary vertex
    Double_t dDCAToPrimVtxNeg = TMath::Abs(v0->DcaNegToPrimVertex());
    Double_t dDCADaughters = v0->DcaV0Daughters(); // dca between daughters
    Double_t dCPA = v0->CosPointingAngle(fPrimaryVertex); // cosine of the pointing angle
    Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
//      Double_t dSecVtxPos[3] = {v0->DecayVertexV0X(),v0->DecayVertexV0Y(),v0->DecayVertexV0Z()}; // V0 vertex position
    v0->GetSecondaryVtx(dSecVtxPos);
    Double_t dRadiusDecay = TMath::Sqrt(dSecVtxPos[0] * dSecVtxPos[0] + dSecVtxPos[1] * dSecVtxPos[1]); // distance of the V0 vertex from the z-axis
    Double_t dEtaDaughterPos = trackPos->Eta(); // pseudorapidity of a daughter track calculated as if primary, != v0->EtaProng(0)
    Double_t dEtaDaughterNeg = trackNeg->Eta(); // != v0->EtaProng(1);
    Double_t dRapK0s = v0->RapK0Short(); // rapidity calculated for K0s assumption
    Double_t dRapLambda = v0->RapLambda(); // rapidity calculated for Lambda assumption
    Double_t dEtaV0 = v0->Eta(); // V0 pseudorapidity
    Double_t dPhiV0 = v0->Phi(); // V0 azimuth
    Double_t dDecayPath[3];
    for(Int_t iPos = 0; iPos < 3; iPos++)
      dDecayPath[iPos] = dSecVtxPos[iPos] - dPrimVtxPos[iPos]; // vector of the V0 path
    Double_t dDecLen = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1] + dDecayPath[2] * dDecayPath[2]); // path length L
    Double_t dDecLen2D = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1]); // transverse path length R
    Double_t dLOverP = dDecLen / v0->P(); // L/p
    Double_t dROverPt = dDecLen2D / dPtV0; // R/pT
    Double_t dMLOverPK0s = dMassPDGK0s * dLOverP; // m*L/p = c*(proper lifetime)
//      Double_t dMLOverPLambda = dMassPDGLambda*dLOverP; // m*L/p
    Double_t dMROverPtK0s = dMassPDGK0s * dROverPt; // m*R/pT
    Double_t dMROverPtLambda = dMassPDGLambda * dROverPt; // m*R/pT
    Double_t dNSigmaPosPion   = (fRespoPID ? TMath::Abs(fRespoPID->NumberOfSigmasTPC(trackPos, AliPID::kPion)) : 0.); // difference between measured and expected signal of the dE/dx in the TPC
    Double_t dNSigmaPosProton = (fRespoPID ? TMath::Abs(fRespoPID->NumberOfSigmasTPC(trackPos, AliPID::kProton)) : 0.);
    Double_t dNSigmaNegPion   = (fRespoPID ? TMath::Abs(fRespoPID->NumberOfSigmasTPC(trackNeg, AliPID::kPion)) : 0.);
    Double_t dNSigmaNegProton = (fRespoPID ? TMath::Abs(fRespoPID->NumberOfSigmasTPC(trackNeg, AliPID::kProton)) : 0.);
    Double_t dAlpha = v0->AlphaV0(); // Armenteros-Podolanski alpha
    Double_t dPtArm = v0->PtArmV0(); // Armenteros-Podolanski pT
    AliAODVertex* prodVtxDaughterPos = (AliAODVertex*)(trackPos->GetProdVertex()); // production vertex of the positive daughter track
    Char_t cTypeVtxProdPos = prodVtxDaughterPos->GetType(); // type of the production vertex
    AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*)(trackNeg->GetProdVertex()); // production vertex of the negative daughter track
    Char_t cTypeVtxProdNeg = prodVtxDaughterNeg->GetType(); // type of the production vertex

    //===== Start of reconstruction cutting =====

    // 1
    // All V0 candidates
    FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    // Start of global cuts
    // 2
    // Reconstruction method
    if(bPrintCuts) printf("Rec: Applying cut: Reconstruction method: %s\n", (fbOnFly ? "on-the-fly" : "offline"));
    if(bOnFlyStatus != fbOnFly)
      continue;
    FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    // 3
    // Tracks TPC OK
    if(bPrintCuts) printf("Rec: Applying cut: Correct charge of daughters\n");
    if(!trackNeg || !trackPos)
      continue;
    if(trackNeg->Charge() == trackPos->Charge()) // daughters have different charge?
      continue;
    if(trackNeg->Charge() != -1) // daughters have expected charge?
      continue;
    if(trackPos->Charge() != 1) // daughters have expected charge?
      continue;

    if(fbTPCRefit)
    {
      if(bPrintCuts) printf("Rec: Applying cut: TPC refit\n");
      if(!trackNeg->IsOn(AliAODTrack::kTPCrefit)) // TPC refit is ON?
        continue;
      if(!trackPos->IsOn(AliAODTrack::kTPCrefit))
        continue;
    }

    if(fbRejectKinks)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Type of production vertex of daughter: No kinks\n");
      if(cTypeVtxProdNeg == AliAODVertex::kKink) // kink daughter rejection
        continue;
      if(cTypeVtxProdPos == AliAODVertex::kKink)
        continue;
    }

    if(fbFindableClusters)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Positive number of findable clusters\n");
      if(dFindableNeg <= 0.)
        continue;
      if(dFindablePos <= 0.)
        continue;
    }

    if(fdCutNCrossedRowsTPCMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Number of TPC rows >= %g\n", fdCutNCrossedRowsTPCMin);
      if(dNRowsNeg < fdCutNCrossedRowsTPCMin) // Crossed TPC padrows
        continue;
      if(dNRowsPos < fdCutNCrossedRowsTPCMin)
        continue;
    }

    if(fdCutCrossedRowsOverFindMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: rows/findable >= %g\n", fdCutCrossedRowsOverFindMin);
      if(dNRowsNeg / dFindableNeg < fdCutCrossedRowsOverFindMin)
        continue;
      if(dNRowsPos / dFindablePos < fdCutCrossedRowsOverFindMin)
        continue;
    }

    if(fdCutCrossedRowsOverFindMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: rows/findable <= %g\n", fdCutCrossedRowsOverFindMax);
      if(dNRowsNeg / dFindableNeg > fdCutCrossedRowsOverFindMax)
        continue;
      if(dNRowsPos / dFindablePos > fdCutCrossedRowsOverFindMax)
        continue;
    }

        FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    // 4
    // Daughters: transverse momentum cut
    if(fdCutPtDaughterMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter pt >= %g\n", fdCutPtDaughterMin);
      if((dPtDaughterNeg < fdCutPtDaughterMin) || (dPtDaughterPos < fdCutPtDaughterMin))
        continue;
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;

    // 5
    // Daughters: Impact parameter of daughters to prim vtx
    if(fdCutDCAToPrimVtxMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter DCA to prim vtx >= %g\n", fdCutDCAToPrimVtxMin);
      if((dDCAToPrimVtxNeg < fdCutDCAToPrimVtxMin) || (dDCAToPrimVtxPos < fdCutDCAToPrimVtxMin))
        continue;
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;

    // 6
    // Daughters: DCA
    if(fdCutDCADaughtersMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: DCA between daughters <= %g\n", fdCutDCADaughtersMax);
      if(dDCADaughters > fdCutDCADaughtersMax)
        continue;
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;

    // 7
    // V0: Cosine of the pointing angle
    if(fdCutCPAKMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: CPA >= %g (K)\n", fdCutCPAKMin);
      if(dCPA < fdCutCPAKMin)
        bIsCandidateK0s = kFALSE;
    }
    if(fdCutCPALMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: CPA >= %g (L, AL)\n", fdCutCPALMin);
      if(dCPA < fdCutCPALMin)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutCPAKMin > 0. || fdCutCPALMin > 0.)
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    // 8
    // V0: Fiducial volume
    if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Decay radius >= %g, <= %g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
      if((dRadiusDecay < fdCutRadiusDecayMin) || (dRadiusDecay > fdCutRadiusDecayMax))
        continue;
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;

    // 9
    // Daughters: pseudorapidity cut
    if(fdCutEtaDaughterMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter |eta| < %g\n", fdCutEtaDaughterMax);
      if((TMath::Abs(dEtaDaughterNeg) > fdCutEtaDaughterMax) || (TMath::Abs(dEtaDaughterPos) > fdCutEtaDaughterMax))
        continue;
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;
    // End of global cuts

    // Start of particle-dependent cuts
    // 10
    // V0: pseudorapidity cut & rapidity cut
    if(fdCutEtaV0Max > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |eta| < %g\n", fdCutEtaV0Max);
      if(TMath::Abs(dEtaV0) > fdCutEtaV0Max)
      {
        bIsCandidateK0s = kFALSE;
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutRapV0Max > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |y| < %g\n", fdCutRapV0Max);
      if(TMath::Abs(dRapK0s) > fdCutRapV0Max)
        bIsCandidateK0s = kFALSE;
      if(TMath::Abs(dRapLambda) > fdCutRapV0Max)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutEtaV0Max > 0. || fdCutRapV0Max > 0.)
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    // 11
    // Lifetime cut
    if(fdCutNTauKMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime < %g (K)\n", fdCutNTauKMax);
      if(dMROverPtK0s > fdCutNTauKMax * dCTauK0s)
        bIsCandidateK0s = kFALSE;
    }
    if(fdCutNTauLMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime < %g (L, AL)\n", fdCutNTauLMax);
      if(dMROverPtLambda > fdCutNTauLMax * dCTauLambda)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutNTauKMax > 0. || fdCutNTauLMax > 0.)
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    iCutIndex++;

    // 12
    // Daughter PID
    if(fdCutNSigmadEdxMax > 0.)
    {
      if(fdPtProtonPIDMax > 0.) // Pb-Pb
      {
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (proton below %g GeV/c) < %g\n", fdPtProtonPIDMax, fdCutNSigmadEdxMax);
        if((dPtDaughterPos < fdPtProtonPIDMax) && (dNSigmaPosProton > fdCutNSigmadEdxMax)) // p+
          bIsCandidateLambda = kFALSE;
        if((dPtDaughterNeg < fdPtProtonPIDMax) && (dNSigmaNegProton > fdCutNSigmadEdxMax)) // p-
          bIsCandidateALambda = kFALSE;
      }
      else // p-p
      {
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (both daughters): < %g\n", fdCutNSigmadEdxMax);
        if(dNSigmaPosPion > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // pi+, pi-
          bIsCandidateK0s = kFALSE;
        if(dNSigmaPosProton > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // p+, pi-
          bIsCandidateLambda = kFALSE;
        if(dNSigmaNegProton > fdCutNSigmadEdxMax || dNSigmaPosPion > fdCutNSigmadEdxMax) // p-, pi+
          bIsCandidateALambda = kFALSE;
      }
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;


    // 13
    // Armenteros-Podolanski cut
    if(fbCutArmPod)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Armenteros-Podolanski (K0S) pT > %g * |alpha|\n", 0.2);
      if(dPtArm < TMath::Abs(0.2 * dAlpha))
        bIsCandidateK0s = kFALSE;
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;

    // 14
    // Cross-contamination
    /*
    if(bIsInPeakK0s)
    {
      if(bIsCandidateLambda) // Lambda candidates in K0s peak, excluded from Lambda candidates by CC cut
        fh2CCLambda->Fill(dMassV0Lambda, dPtV0);
    }
    if(bIsInPeakLambda)
    {
      if(bIsCandidateK0s) // K0s candidates in Lambda peak, excluded from K0s candidates by CC cut
        fh2CCK0s->Fill(dMassV0K0s, dPtV0);
    }
    */
    if(fbCutCross)
    {
      if(bIsInPeakK0s)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
      if(bIsInPeakLambda)
      {
        bIsCandidateK0s = kFALSE;
      }
      if(bIsInPeakALambda)
      {
        bIsCandidateK0s = kFALSE;
      }
          FillCandidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
    }
    iCutIndex++;
    // End of particle-dependent cuts

  if(bIsCandidateK0s) {fh2dKshortMassVsPt->Fill(dPtV0, dMassV0K0s,fPythiaEventWeight); }
  if(bIsCandidateLambda) {fh2dLamdaMassVsPt->Fill(dPtV0, dMassV0Lambda,fPythiaEventWeight); }
  if(bIsCandidateALambda) {fh2dAnLamdaMassVsPt->Fill(dPtV0, dMassV0ALambda,fPythiaEventWeight); }

  AliEmcalJet * jetrec  = 0x0;
  double fJetPt=0;
 
  fJetContainerData->ResetCurrentID();

  while ((jetrec = fJetContainerData->GetNextAcceptJet()))
  {

	fJetPt= jetrec->Pt();

	if(!(fJetContainerData->GetRhoParameter() == 0x0))
	{
		fJetPt = fJetPt - fJetContainerData->GetRhoVal() * jetrec->Area();
	}

	//if(!(fJetCutsHF->IsJetSelected(jetrec))) continue;

	if(fJetPt < 5.) continue;

	// make inclusive signed imp. parameter constituent histograms

	if(bIsCandidateK0s && IsParticleInCone(v0, jetrec, 0.4)) {

        	Double_t valueKInJC[4] = {dMassV0K0s, dPtV0, dEtaV0, fJetPt};
        	fhnV0InJetK0s->Fill(valueKInJC, fPythiaEventWeight);

	}
	if(bIsCandidateLambda && IsParticleInCone(v0, jetrec, 0.4)) { 

        	Double_t valueLInJC[4] = {dMassV0Lambda, dPtV0, dEtaV0, fJetPt};
        	fhnV0InJetLambda->Fill(valueLInJC, fPythiaEventWeight);

	}
	if(bIsCandidateALambda && IsParticleInCone(v0, jetrec, 0.4)) {

		Double_t valueLInJC[4] = {dMassV0ALambda, dPtV0, dEtaV0, fJetPt};
		fhnV0InJetALambda->Fill(valueLInJC, fPythiaEventWeight);

	}

   }
  fJetContainerData->ResetCurrentID();
  jetrec=NULL; delete jetrec;

  if(bIsCandidateK0s || bIsCandidateLambda || bIsCandidateALambda)
	 new((*fV0CandidateArray)[fV0CandidateArray->GetEntriesFast()]) AliAODv0(*v0);


  if(fIsPythia){

     AliAODMCHeader* headerMC = 0; // MC header
     Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex

     headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
     if(!headerMC)
     {
	   AliError("No MC header found!");
	   return kFALSE;
     }
    // get position of the MC primary vertex
    dPrimVtxMCX = headerMC->GetVtxX();
    dPrimVtxMCY = headerMC->GetVtxY();
    dPrimVtxMCZ = headerMC->GetVtxZ();

    Int_t iNTracksMC = fMCArray->GetEntriesFast();
      // Associate selected candidates only
//          if ( !(bIsCandidateK0s && bIsInPeakK0s) && !(bIsCandidateLambda && bIsInPeakLambda) ) // signal candidates
      if(!(bIsCandidateK0s) && !(bIsCandidateLambda)  && !(bIsCandidateALambda)) // chosen candidates with any mass
        continue;

      // Get MC labels of reconstructed daughter tracks
      Int_t iLabelPos = TMath::Abs(trackPos->GetLabel());
      Int_t iLabelNeg = TMath::Abs(trackNeg->GetLabel());

      // Make sure MC daughters are in the array range
      if((iLabelNeg < 0) || (iLabelNeg >= iNTracksMC) || (iLabelPos < 0) || (iLabelPos >= iNTracksMC))
        continue;

      // Get MC particles corresponding to reconstructed daughter tracks
      AliAODMCParticle* particleMCDaughterNeg = (AliAODMCParticle*)GetMCTrack(trackNeg);
      AliAODMCParticle* particleMCDaughterPos = (AliAODMCParticle*)GetMCTrack(trackPos);
      if(!particleMCDaughterNeg || !particleMCDaughterPos)
        continue;

      // Make sure MC daughter particles are not physical primary
      //if((particleMCDaughterNeg->IsPhysicalPrimary()) || (particleMCDaughterPos->IsPhysicalPrimary()))
      //  continue;

      // Get identities of MC daughter particles
      Int_t iPdgCodeDaughterPos = particleMCDaughterPos->GetPdgCode();
      Int_t iPdgCodeDaughterNeg = particleMCDaughterNeg->GetPdgCode();

      // Get index of the mother particle for each MC daughter particle
      Int_t iIndexMotherPos = particleMCDaughterPos->GetMother();
      Int_t iIndexMotherNeg = particleMCDaughterNeg->GetMother();

      if((iIndexMotherNeg < 0) || (iIndexMotherNeg >= iNTracksMC) || (iIndexMotherPos < 0) || (iIndexMotherPos >= iNTracksMC))
        continue;

      // Check whether MC daughter particles have the same mother
      if(iIndexMotherNeg != iIndexMotherPos)
        continue;

      // Get the MC mother particle of both MC daughter particles
      AliAODMCParticle* particleMCMother = (AliAODMCParticle*)fMCArray->At(iIndexMotherPos);
      if(!particleMCMother)
        continue;

      // Get identity of the MC mother particle
      Int_t iPdgCodeMother = particleMCMother->GetPdgCode();

      // Skip not interesting particles
      if((iPdgCodeMother != iPdgCodeK0s) && (TMath::Abs(iPdgCodeMother) != iPdgCodeLambda))
        continue;

      // Check identity of the MC mother particle and the decay channel
      // Is MC mother particle K0S?
      Bool_t bV0MCIsK0s = ((iPdgCodeMother == iPdgCodeK0s) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodePion));
      // Is MC mother particle Lambda?
      Bool_t bV0MCIsLambda = ((iPdgCodeMother == +iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
      // Is MC mother particle anti-Lambda?
      Bool_t bV0MCIsALambda = ((iPdgCodeMother == -iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));

      Double_t dPtV0Gen = particleMCMother->Pt();
      Double_t dRapV0Gen = particleMCMother->Y();
      Double_t dEtaV0Gen = particleMCMother->Eta();
//      Double_t dPhiV0Gen = particleMCMother->Phi();

      // V0 pseudorapidity cut applied on generated particles
      if(fdCutEtaV0Max > 0.)
      {
        if(bPrintCuts) printf("Rec->Gen: Applying cut: V0 |eta|: < %g\n", fdCutEtaV0Max);
        if((TMath::Abs(dEtaV0Gen) > fdCutEtaV0Max))
          continue;
      }
      // V0 rapidity cut applied on generated particles
      if(fdCutRapV0Max > 0.)
      {
        if(bPrintCuts) printf("Rec->Gen: Applying cut: V0 |y|: < %g\n", fdCutRapV0Max);
        if((TMath::Abs(dRapV0Gen) > fdCutRapV0Max))
          continue;
      }

      // Select only particles from a specific generator
      //if(!IsFromGoodGenerator(iIndexMotherPos))
      //  continue;

      // Is MC mother particle physical primary? Attention!! Definition of IsPhysicalPrimary may change!!

      // Get the distance between production point of the MC mother particle and the primary vertex
      Double_t dx = dPrimVtxMCX - particleMCMother->Xv();
      Double_t dy = dPrimVtxMCY - particleMCMother->Yv();
      Double_t dz = dPrimVtxMCZ - particleMCMother->Zv();
      Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      Bool_t bV0MCIsPrimaryDist = (dDistPrimary < dDistPrimaryMax); // Is close enough to be considered primary-like?

      // K0s
//          if (bIsCandidateK0s && bIsInPeakK0s) // selected candidates in peak
      if(bIsCandidateK0s) // selected candidates with any mass
      {
        if(bV0MCIsK0s && bV0MCIsPrimaryDist) // well reconstructed candidates
        {
          fh2dKshortMassVsPtReal->Fill(dPtV0, dMassV0K0s, fPythiaEventWeight);
	  fh2dKshortRecPtVsGenPt->Fill(dPtV0, dPtV0Gen, fPythiaEventWeight);
        }
      }

      // Lambda
//          if (bIsCandidateLambda && bIsInPeakLambda) // selected candidates in peak
      if(bIsCandidateLambda) // selected candidates with any mass
      {
        if(bV0MCIsLambda && bV0MCIsPrimaryDist) // well reconstructed candidates
        {
          fh2dLamdaMassVsPtReal->Fill(dPtV0, dMassV0Lambda, fPythiaEventWeight);
          fh2dLamdaRecPtVsGenPt->Fill(dPtV0, dPtV0Gen, fPythiaEventWeight);
        }
      }

      // anti-Lambda
//          if (bIsCandidateALambda && bIsInPeakALambda) // selected candidates in peak
      if(bIsCandidateALambda) // selected candidates with any mass
      {
        if(bV0MCIsALambda && bV0MCIsPrimaryDist) // well reconstructed candidates
        {
          fh2dAnLamdaMassVsPtReal->Fill(dPtV0, dMassV0ALambda, fPythiaEventWeight);
          fh2dAnLamdaRecPtVsGenPt->Fill(dPtV0, dPtV0Gen, fPythiaEventWeight);
        }
      }

    }
	

 }
  
  return kTRUE;
}
//==========================================
void AliAnalysisTaskBJetTC::FillCandidates(Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut/*cut index*/)
{
  if(isK)
  {
    fh1V0CounterCentK0s->Fill(iCut);
  }
  if(isL)
  {
    fh1V0CounterCentLambda->Fill(iCut);
  }
  if(isAL)
  {
    fh1V0CounterCentALambda->Fill(iCut);
  }
}
//=====================================================================
Bool_t AliAnalysisTaskBJetTC::IsV0InJet(TVector3 vV0, Double_t dJetPtMin)
{  
  if (!fJetContainerData) return kFALSE;

  TVector3 vJet;
  Double_t dJetRadius = fJetContainerData->GetJetRadius();
  fJetContainerData->ResetCurrentID();
  AliEmcalJet *pJet = fJetContainerData->GetNextAcceptJet(); while (pJet) {
    Double_t dPt = fJetContainerData->GetJetPtCorr(fJetContainerData->GetCurrentID());
    if (dPt<dJetPtMin) { pJet = fJetContainerData->GetNextAcceptJet(); continue; }

    vJet.SetPtEtaPhi(dPt, pJet->Eta(), pJet->Phi());
    if (vJet.DeltaR(vV0)<dJetRadius) return kTRUE;
    pJet = fJetContainerData->GetNextAcceptJet();
  }
  pJet=NULL;
  delete pJet;

  return kFALSE;
}
//=====================================================================
Bool_t AliAnalysisTaskBJetTC::IsParticleInCone(const AliVParticle* part, const AliEmcalJet* jet, Double_t dRMax) const
{
// decides whether a particle is inside a jet cone
  if(!part || !jet)
    return kFALSE;

  TVector3 vecMom2(jet->Px(), jet->Py(), jet->Pz());
  TVector3 vecMom1(part->Px(), part->Py(), part->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  if(dR <= dRMax) // momentum vectors of part1 and part2 are closer than dRMax
    return kTRUE;
  return kFALSE;
}

//=====================================================================
Bool_t AliAnalysisTaskBJetTC::IsElectronHF( AliAODTrack* track){

        Double_t pid_ele = 0.0;

        // get track information
        Double_t pt = track->Pt(); 
        Double_t eta = track->Eta(); 
        Double_t d0z0[2]={-999,-999}, cov[3];

          if(!track->PropagateToDCA(fPrimaryVertex, fAODIn->GetMagneticField(), 20., d0z0, cov)) return kFALSE;

	//Don't forget to reject Kink daughters :(

	if(pt<0.5) return kFALSE;
        if(fabs(eta)>0.9) return kFALSE;
        if(fabs(d0z0[0])>2.4) return kFALSE;
        if(fabs(d0z0[1])>3.2) return kFALSE;
        if(track->GetTPCNcls() < 80)  return kFALSE;
        if(track->GetITSNcls() < 2) return kFALSE;   // AOD track level

        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return kFALSE;    // kAny
        if((!(track->GetStatus()&AliESDtrack::kITSrefit)|| (!(track->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;

	AliAODVertex* ProdVertex = (AliAODVertex*)(track->GetProdVertex()); // production vertex of the track
        Char_t cTypeVtxProd = ProdVertex->GetType(); // type of the production vertex

 	if(cTypeVtxProd == AliAODVertex::kKink) return kFALSE; // kink daughter rejection 

	// Get TPC nSigma
        Double_t fTPCnSigma=-999;
        fTPCnSigma = fRespoPID->NumberOfSigmasTPC(track, AliPID::kElectron);

	// Get TOF nSigma
        Double_t fTOFnSigma=-999;
        fTOFnSigma = fRespoPID->NumberOfSigmasTOF(track, AliPID::kElectron);

	if(pt>0.5 && pt<2.5){
	     if(TMath::Abs(fTOFnSigma) > 3)
		  return kFALSE;
	}

	if(fIsPythia){
		AliAODMCParticle* ParticleMC = GetMCTrack(track);
		if(abs(ParticleMC->GetPdgCode())==11) pid_ele = 1.0;
	}

	if(pid_ele==1.0)
            fTPCnsigMcEle->Fill(track->Pt(),fTPCnSigma);
        else
            fTPCnsigMcHad->Fill(track->Pt(),fTPCnSigma);
	

        if(fTPCnSigma<-0.5 || fTPCnSigma>3.)return kFALSE;  //++++++++

	Bool_t EmcalAccepted = kFALSE;

	///////////////////////////
        //Track matching to EMCAL//
        ///////////////////////////

	if(pt>6 ){

            if(!track->IsEMCAL()) return kFALSE;
            Int_t EMCalIndex = -1;
            EMCalIndex = track->GetEMCALcluster();
            if(EMCalIndex < 0) return kFALSE;

            AliVCluster *clustMatch=0x0;
            //clustMatch = (AliVCluster*)fAODIn->GetCaloCluster(EMCalIndex);
            clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters->At(EMCalIndex));

            if(clustMatch && clustMatch->IsEMCAL())
            {
                Double_t fPhiDiff = -999, fEtaDiff = -999;
                GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);

                if(TMath::Abs(fPhiDiff) > 0.05 || TMath::Abs(fEtaDiff)> 0.05) return kFALSE;
            
            
                /////////////////////////////////////////////
                //Properties of tracks matched to the EMCAL//
                /////////////////////////////////////////////           
                Double_t clustMatchE = clustMatch->E();

                //EMCAL EID info
                Double_t eop = -1.0;
                Double_t m02 = -99999,m20 = -99999;
                if(track->P()>0)eop = clustMatchE/track->P();
                m02 =clustMatch->GetM02();
                m20 =clustMatch->GetM20();
            
                if(pid_ele==1.0)
                    fHistMcEopEle->Fill(track->Pt(),eop);
                else
                    fHistMcEopHad->Fill(track->Pt(),eop);
            
                if(!(eop>0.8 && eop<1.2) ) return kFALSE;

            }
	EmcalAccepted = kTRUE;

	}

	if(!EmcalAccepted && pt>6.0) return kFALSE;

	return kTRUE;

}
//++++++++++++++++++++++++++++++++++
Double_t AliAnalysisTaskBJetTC::GetPtRel(AliAODTrack* lep, AliEmcalJet* jet, Bool_t addLepToJet){

    TVector3 LepJetAxis;
    LepJetAxis.SetXYZ((addLepToJet) ? lep->Px()+jet->Px() : jet->Px(), (addLepToJet) ? lep->Py()+jet->Py() : jet->Py(), (addLepToJet) ? lep->Pz()+jet->Pz() : jet->Pz());

    TVector3 Lepton;
    Lepton.SetXYZ(lep->Px(), lep->Py(), lep->Pz());
 
    Double_t PtRel=0;
    PtRel = Lepton.Perp(LepJetAxis);

    return PtRel;

}
//________________________________________________________________________
void AliAnalysisTaskBJetTC::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
    // Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface
    
    phidiff = 999;
    etadiff = 999;
    
    if (!t||!v) return;
    
    Double_t veta = t->GetTrackEtaOnEMCal();
    Double_t vphi = t->GetTrackPhiOnEMCal();
    
    Float_t pos[3] = {0};
    v->GetPosition(pos);
    TVector3 cpos(pos);
    Double_t ceta     = cpos.Eta();
    Double_t cphi     = cpos.Phi();
    etadiff=veta-ceta;
    phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
////////////////////////////////////////////////////////////////////////////////
Double_t AliAnalysisTaskBJetTC::GetDeltaPtRandomCone()
{

	Double_t deltaPt = -1000.;
	AliParticleContainer* partcont = 0x0;
	partcont = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
	Double_t jetradius = fJetContainerData->GetJetRadius();
	Double_t minEta = -0.5;
	Double_t maxEta = 0.5;
	Double_t tmpRandConeEta = -999;
	Double_t tmpRandConePhi = -999;
	Double_t tmpConePt = -1.;

	AliEmcalJet* LeadingJet = NULL;
	AliEmcalJet* SubLeadingJet = NULL;
	Double_t LJeta = 999;
	Double_t LJphi = 999;

	Double_t SLJeta = 999;
	Double_t SLJphi = 999;

	if (fJetContainerData){

		Float_t maxJetPts[] = { 0,  0};

		fJetContainerData->ResetCurrentID();

  		AliEmcalJet * jet  = 0x0;

		while ((jet = fJetContainerData->GetNextAcceptJet())){

			if (!jet)  continue;

			if (jet->Pt() > maxJetPts[0]) {
				maxJetPts[1] = maxJetPts[0];
				SubLeadingJet = LeadingJet;
				maxJetPts[0] = jet->Pt();
				LeadingJet = jet;
			} else if (jet->Pt() > maxJetPts[1]) {
				maxJetPts[1] = jet->Pt();
				SubLeadingJet = jet;
			}
		}

		jet=NULL; delete jet;

		if(LeadingJet){
			LJeta = LeadingJet->Eta();
			LJphi = LeadingJet->Phi();
		}
		if(SubLeadingJet){
			SLJeta = SubLeadingJet->Eta();
			SLJphi = SubLeadingJet->Phi();
		}
	}

	LeadingJet=NULL; SubLeadingJet=NULL;
	delete LeadingJet; delete SubLeadingJet;

  	Double_t dLJ = 0;
  	Double_t dSLJ = 0;
  	Int_t repeats = 0;

	do {
	    tmpRandConeEta = minEta + fRandom->Rndm() * (maxEta - minEta);
	    tmpRandConePhi = fRandom->Rndm() * TMath::TwoPi();
	    dLJ = TMath::Sqrt((LJeta - tmpRandConeEta) * (LJeta - tmpRandConeEta) + (LJphi - tmpRandConePhi) * (LJphi - tmpRandConePhi));
	    dSLJ = TMath::Sqrt((SLJeta - tmpRandConeEta) * (SLJeta - tmpRandConeEta) + (SLJphi - tmpRandConePhi) * (SLJphi - tmpRandConePhi));
	    repeats++;

	  } while (dLJ < 0.45 || dSLJ < 0.45);

	AliVTrack* tmpTrack = 0x0;
	AliAODTrack* trackAOD = 0x0;

	for(Int_t i = 0; i < partcont->GetNAcceptedParticles(); i++) {

		if(!partcont->GetParticle(i)) continue;
		tmpTrack = static_cast<AliVTrack*>(partcont->GetParticle(i));
		trackAOD = (AliAODTrack*)partcont->GetParticle(i);
		if(!((trackAOD)->TestFilterBit(1 << 4)) && !((trackAOD)->TestFilterBit(1 << 9)) ) continue;

		if(fabs(tmpTrack->Eta()) > 0.9) continue;

		if(tmpTrack->Pt() < 0.15) continue;

		if(sqrt((tmpTrack->Eta() - tmpRandConeEta) * (tmpTrack->Eta() - tmpRandConeEta) +
				TVector2::Phi_mpi_pi((tmpTrack->Phi() - tmpRandConePhi)) *
				TVector2::Phi_mpi_pi((tmpTrack->Phi() - tmpRandConePhi))) < jetradius) {
			tmpConePt += tmpTrack->Pt();
		}
	}
	tmpTrack=NULL; trackAOD=NULL;
	delete tmpTrack; delete trackAOD;

	partcont=NULL;
	delete partcont;

	if(tmpConePt > 0) {
		deltaPt = tmpConePt - jetradius * jetradius * TMath::Pi() * fJetContainerData->GetRhoVal();
		return deltaPt;
	}
	return deltaPt;
}
//_________________________________________________________________________
Double_t AliAnalysisTaskBJetTC::GetDeltaPtRandomConeWithSignal()
{

	Double_t deltaPt = -1000.;
	AliParticleContainer* partcont = 0x0;
	partcont = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
	Double_t jetradius = fJetContainerData->GetJetRadius();
	Double_t minEta = -0.5;
	Double_t maxEta = 0.5;
	Double_t tmpRandConeEta = minEta + fRandom->Rndm() * (maxEta - minEta);
	Double_t tmpRandConePhi = fRandom->Rndm() * TMath::TwoPi();
	Double_t tmpConePt = -1.;

	for(Int_t i = 0; i < partcont->GetNAcceptedParticles(); i++) {

		if(!partcont->GetParticle(i)) continue;
		AliVTrack* tmpTrack = static_cast<AliVTrack*>(partcont->GetParticle(i));
		AliAODTrack* trackAOD = (AliAODTrack*)partcont->GetParticle(i);
		if(!((trackAOD)->TestFilterBit(1 << 4)) && !((trackAOD)->TestFilterBit(1 << 9)) ) continue;

		if(fabs(tmpTrack->Eta()) > 0.9) continue;

		if(tmpTrack->Pt() < 0.15) continue;

		if(sqrt((tmpTrack->Eta() - tmpRandConeEta) * (tmpTrack->Eta() - tmpRandConeEta) +
				TVector2::Phi_mpi_pi((tmpTrack->Phi() - tmpRandConePhi)) *
				TVector2::Phi_mpi_pi((tmpTrack->Phi() - tmpRandConePhi))) < jetradius) {
			tmpConePt += tmpTrack->Pt();
		}
	}

	if(tmpConePt > 0) {
		deltaPt = tmpConePt - jetradius * jetradius * TMath::Pi() * fJetContainerData->GetRhoVal();
		return deltaPt;
	}
	return deltaPt;
}
//_____________________________________________________________________________________
Int_t AliAnalysisTaskBJetTC::FindVertices6Prong(const AliEmcalJet* jet,
                                           TClonesArray*      fTrackArrayIn,
                                           AliAODEvent*       aodEvent,
                                           AliESDVertex*      primaryESDVertex,
                                           Double_t           magZkG,
                                           TClonesArray*      arrayVtxHF,
                                           Int_t&             nDauRejCount)
{

  Int_t nSecndVxtHF = 0;

  arrayVtxHF->Clear();

  Double_t vtxRes = 0.;

  Int_t nTrksInJet = jet->GetNumberOfTracks();
  AliDebugF(6, MSGINFO("nTrksInJet = %d \n"), nTrksInJet);
  if (nTrksInJet < 6) {
    AliDebug(2, MSGWARNING("Cannot find vertices w/ only one track"));
    return -3;
  }

  //make array of ESD tracks, then needed for fTrackArray
  vctr_pair_int_esdTrk vecESDTrks;
  vecESDTrks.reserve(nTrksInJet);
 
  for (Int_t j = 0; j < nTrksInJet; ++j) {
    AliAODTrack* jTrk   = ((AliAODTrack*)jet->TrackAt(j, fTrackArrayIn));
    if (!jTrk) {
      AliWarningF(MSGWARNING("Track in Jet with index %d/%d not found. Total number of AODtracks %d"),
                  j, nTrksInJet, fTrackArrayIn->GetEntries());
      continue;
    }
    //utilize dynamic cast and then check pointer
    Int_t jTrkID = jTrk->GetID();
    if (jTrkID < 0) {
      AliDebugF(6, MSGINFO("Track with index < 0 %d"), jTrkID);
      continue;
    }

    if (!fjetCuts3Prong->IsDaughterSelected(jTrk, primaryESDVertex, fEsdTrackCuts)){
      nDauRejCount++;
      continue;
    }
    
    AliESDtrack* tmpESDtrk = new AliESDtrack(jTrk);
    vecESDTrks.push_back(make_pair(j, tmpESDtrk));
  }

  Int_t nGoodTrks = (Int_t)vecESDTrks.size();
  if (nGoodTrks < 6) {
    AliDebugF(6, MSGDEBUG("Number of good tracks = %d"), nGoodTrks);
    return -4;
  }

  Int_t up = nGoodTrks - 5;
  Int_t nVtxContributorsBelongToV0 = 0;
  for (Int_t it1 = 0; it1 < up; ++it1) {

    Int_t        jTrkID_1 = (vecESDTrks.at(it1)).first;
    AliESDtrack* esdTrk_1 = (vecESDTrks.at(it1)).second;

    fTrackArray->Clear();
    fTrackArray->AddAt(esdTrk_1, 0);

    for (Int_t it2 = it1 + 1; it2 < up + 1; ++it2) {

      Int_t        jTrkID_2 = (vecESDTrks.at(it2)).first;
      AliESDtrack* esdTrk_2 = (vecESDTrks.at(it2)).second;

      fTrackArray->AddAt(esdTrk_2, 1);


        for (Int_t it3 = it2 + 1; it3 < up + 2; ++it3) {

          Int_t        jTrkID_3 = (vecESDTrks.at(it3)).first;
          AliESDtrack* esdTrk_3 = (vecESDTrks.at(it3)).second;

          fTrackArray->AddAt(esdTrk_3, 2);

	  for(Int_t it4= it3 + 1; it4 < up+3; ++it4){

		  Int_t        jTrkID_4 = (vecESDTrks.at(it4)).first;
		  AliESDtrack* esdTrk_4 = (vecESDTrks.at(it4)).second;

		  fTrackArray->AddAt(esdTrk_4, 3);

		  for(Int_t it5= it4 + 1; it5 < nGoodTrks; ++it5){

			  Int_t        jTrkID_5 = (vecESDTrks.at(it5)).first;
			  AliESDtrack* esdTrk_5 = (vecESDTrks.at(it5)).second;

			  fTrackArray->AddAt(esdTrk_5, 4);

			  for(Int_t it6= it5 + 1; it6 < nGoodTrks; ++it6){

				  Int_t        jTrkID_6 = (vecESDTrks.at(it6)).first;
				  AliESDtrack* esdTrk_6 = (vecESDTrks.at(it6)).second;

				  fTrackArray->AddAt(esdTrk_6, 5);

				  AliAODVertex* secAODVertex = fVtxTagger3Prong->ReconstructSecondaryVertex(fTrackArray, primaryESDVertex, magZkG, vtxRes); 

				  if (secAODVertex) {
				    AliAODTrack* aodTrk_1 = (AliAODTrack*)jet->TrackAt(jTrkID_1, fTrackArrayIn);
				    AliAODTrack* aodTrk_2 = (AliAODTrack*)jet->TrackAt(jTrkID_2, fTrackArrayIn);
				    AliAODTrack* aodTrk_3 = (AliAODTrack*)jet->TrackAt(jTrkID_3, fTrackArrayIn);
				    AliAODTrack* aodTrk_4 = (AliAODTrack*)jet->TrackAt(jTrkID_4, fTrackArrayIn);
				    AliAODTrack* aodTrk_5 = (AliAODTrack*)jet->TrackAt(jTrkID_5, fTrackArrayIn);
				    AliAODTrack* aodTrk_6 = (AliAODTrack*)jet->TrackAt(jTrkID_6, fTrackArrayIn);

				    secAODVertex->AddDaughter(aodTrk_1);
				    secAODVertex->AddDaughter(aodTrk_2);
				    secAODVertex->AddDaughter(aodTrk_3);
				    secAODVertex->AddDaughter(aodTrk_4);
				    secAODVertex->AddDaughter(aodTrk_5);
				    secAODVertex->AddDaughter(aodTrk_6);

				    if (!fjetCuts3Prong->IsVertexSelected(secAODVertex, aodEvent, magZkG, vtxRes))
				      continue;

				    new ((* arrayVtxHF)[nSecndVxtHF]) AliAODVertex(* secAODVertex);
				    nSecndVxtHF++;
				  } // end if (vert)
			}// end for it6
		} //end for it5
	   } //end for it4
        } // end for it3
    } // end for it2
  } // end for it1

    //cout<<"This is a 6prong SV\n";

  fTrackArray->Clear();

  for (vctr_pair_int_esdTrk::iterator it = vecESDTrks.begin(); it != vecESDTrks.end(); ++it) {
      AliESDtrack* lESDtrk = (* it).second;
      if (lESDtrk)
        delete lESDtrk;
  }

  return nSecndVxtHF;
}
//_____________________________________________________________________________________
Int_t AliAnalysisTaskBJetTC::FindVertices5Prong(const AliEmcalJet* jet,
                                           TClonesArray*      fTrackArrayIn,
                                           AliAODEvent*       aodEvent,
                                           AliESDVertex*      primaryESDVertex,
                                           Double_t           magZkG,
                                           TClonesArray*      arrayVtxHF,
                                           Int_t&             nDauRejCount)
{

  Int_t nSecndVxtHF = 0;

  arrayVtxHF->Clear();

  Double_t vtxRes = 0.;

  Int_t nTrksInJet = jet->GetNumberOfTracks();
  AliDebugF(6, MSGINFO("nTrksInJet = %d \n"), nTrksInJet);
  if (nTrksInJet < 5) {
    AliDebug(2, MSGWARNING("Cannot find vertices w/ only one track"));
    return -3;
  }

  //make array of ESD tracks, then needed for fTrackArray
  vctr_pair_int_esdTrk vecESDTrks;
  vecESDTrks.reserve(nTrksInJet);
 
  for (Int_t j = 0; j < nTrksInJet; ++j) {
    AliAODTrack* jTrk   = ((AliAODTrack*)jet->TrackAt(j, fTrackArrayIn));
    if (!jTrk) {
      AliWarningF(MSGWARNING("Track in Jet with index %d/%d not found. Total number of AODtracks %d"),
                  j, nTrksInJet, fTrackArrayIn->GetEntries());
      continue;
    }
    //utilize dynamic cast and then check pointer
    Int_t jTrkID = jTrk->GetID();
    if (jTrkID < 0) {
      AliDebugF(6, MSGINFO("Track with index < 0 %d"), jTrkID);
      continue;
    }

    if (!fjetCuts3Prong->IsDaughterSelected(jTrk, primaryESDVertex, fEsdTrackCuts)){
      nDauRejCount++;
      continue;
    }
    
    AliESDtrack* tmpESDtrk = new AliESDtrack(jTrk);
    vecESDTrks.push_back(make_pair(j, tmpESDtrk));
  }

  Int_t nGoodTrks = (Int_t)vecESDTrks.size();
  if (nGoodTrks < 5) {
    AliDebugF(6, MSGDEBUG("Number of good tracks = %d"), nGoodTrks);
    return -4;
  }

  Int_t up = nGoodTrks - 4;
  Int_t nVtxContributorsBelongToV0 = 0;
  for (Int_t it1 = 0; it1 < up; ++it1) {

    Int_t        jTrkID_1 = (vecESDTrks.at(it1)).first;
    AliESDtrack* esdTrk_1 = (vecESDTrks.at(it1)).second;

    fTrackArray->Clear();
    fTrackArray->AddAt(esdTrk_1, 0);

    for (Int_t it2 = it1 + 1; it2 < up + 1; ++it2) {

      Int_t        jTrkID_2 = (vecESDTrks.at(it2)).first;
      AliESDtrack* esdTrk_2 = (vecESDTrks.at(it2)).second;

      fTrackArray->AddAt(esdTrk_2, 1);


        for (Int_t it3 = it2 + 1; it3 < up + 2; ++it3) {

          Int_t        jTrkID_3 = (vecESDTrks.at(it3)).first;
          AliESDtrack* esdTrk_3 = (vecESDTrks.at(it3)).second;

          fTrackArray->AddAt(esdTrk_3, 2);

	  for(Int_t it4= it3 + 1; it4 < up+3; ++it4){

		  Int_t        jTrkID_4 = (vecESDTrks.at(it4)).first;
		  AliESDtrack* esdTrk_4 = (vecESDTrks.at(it4)).second;

		  fTrackArray->AddAt(esdTrk_4, 3);

		  for(Int_t it5= it4 + 1; it5 < nGoodTrks; ++it5){

			  Int_t        jTrkID_5 = (vecESDTrks.at(it5)).first;
			  AliESDtrack* esdTrk_5 = (vecESDTrks.at(it5)).second;

			  fTrackArray->AddAt(esdTrk_5, 4);

			  AliAODVertex* secAODVertex = fVtxTagger3Prong->ReconstructSecondaryVertex(fTrackArray, primaryESDVertex, magZkG, vtxRes); 

			  if (secAODVertex) {
			    AliAODTrack* aodTrk_1 = (AliAODTrack*)jet->TrackAt(jTrkID_1, fTrackArrayIn);
			    AliAODTrack* aodTrk_2 = (AliAODTrack*)jet->TrackAt(jTrkID_2, fTrackArrayIn);
			    AliAODTrack* aodTrk_3 = (AliAODTrack*)jet->TrackAt(jTrkID_3, fTrackArrayIn);
			    AliAODTrack* aodTrk_4 = (AliAODTrack*)jet->TrackAt(jTrkID_4, fTrackArrayIn);
			    AliAODTrack* aodTrk_5 = (AliAODTrack*)jet->TrackAt(jTrkID_5, fTrackArrayIn);

			    secAODVertex->AddDaughter(aodTrk_1);
			    secAODVertex->AddDaughter(aodTrk_2);
			    secAODVertex->AddDaughter(aodTrk_3);
			    secAODVertex->AddDaughter(aodTrk_4);
			    secAODVertex->AddDaughter(aodTrk_5);

			    if (!fjetCuts3Prong->IsVertexSelected(secAODVertex, aodEvent, magZkG, vtxRes))
			      continue;

			    new ((* arrayVtxHF)[nSecndVxtHF]) AliAODVertex(* secAODVertex);
			    nSecndVxtHF++;
			  } // end if (vert)
		} //end for it5
	   } //end for it4
        } // end for it3
    } // end for it2
  } // end for it1

    //cout<<"This is a 5prong SV\n";

  fTrackArray->Clear();

  for (vctr_pair_int_esdTrk::iterator it = vecESDTrks.begin(); it != vecESDTrks.end(); ++it) {
      AliESDtrack* lESDtrk = (* it).second;
      if (lESDtrk)
        delete lESDtrk;
  }

  return nSecndVxtHF;
}
//_____________________________________________________________________________________
Int_t AliAnalysisTaskBJetTC::FindVertices4Prong(const AliEmcalJet* jet,
                                           TClonesArray*      fTrackArrayIn,
                                           AliAODEvent*       aodEvent,
                                           AliESDVertex*      primaryESDVertex,
                                           Double_t           magZkG,
                                           TClonesArray*      arrayVtxHF,
                                           Int_t&             nDauRejCount)
{

  Int_t nSecndVxtHF = 0;

  arrayVtxHF->Clear();

  Double_t vtxRes = 0.;

  Int_t nTrksInJet = jet->GetNumberOfTracks();
  AliDebugF(6, MSGINFO("nTrksInJet = %d \n"), nTrksInJet);
  if (nTrksInJet < 4) {
    AliDebug(2, MSGWARNING("Cannot find vertices w/ only one track"));
    return -3;
  }

  //make array of ESD tracks, then needed for fTrackArray

  vctr_pair_int_esdTrk vecESDTrks;
  vecESDTrks.reserve(nTrksInJet);
 
  for (Int_t j = 0; j < nTrksInJet; ++j) {
    AliAODTrack* jTrk   = ((AliAODTrack*)jet->TrackAt(j, fTrackArrayIn));
    if (!jTrk) {
      AliWarningF(MSGWARNING("Track in Jet with index %d/%d not found. Total number of AODtracks %d"),
                  j, nTrksInJet, fTrackArrayIn->GetEntries());
      continue;
    }
    //utilize dynamic cast and then check pointer
    Int_t jTrkID = jTrk->GetID();
    if (jTrkID < 0) {
      AliDebugF(6, MSGINFO("Track with index < 0 %d"), jTrkID);
      continue;
    }

    if (!fjetCuts3Prong->IsDaughterSelected(jTrk, primaryESDVertex, fEsdTrackCuts)){
      nDauRejCount++;
      continue;
    }
    
    AliESDtrack* tmpESDtrk = new AliESDtrack(jTrk);
    vecESDTrks.push_back(make_pair(j, tmpESDtrk));
  }

  Int_t nGoodTrks = (Int_t)vecESDTrks.size();
  if (nGoodTrks < 4) {
    AliDebugF(6, MSGDEBUG("Number of good tracks = %d"), nGoodTrks);
    return -4;
  }

  Int_t up = nGoodTrks - 3;
  Int_t nVtxContributorsBelongToV0 = 0;
  for (Int_t it1 = 0; it1 < up; ++it1) {

    Int_t        jTrkID_1 = (vecESDTrks.at(it1)).first;
    AliESDtrack* esdTrk_1 = (vecESDTrks.at(it1)).second;

    fTrackArray->Clear();
    fTrackArray->AddAt(esdTrk_1, 0);

    for (Int_t it2 = it1 + 1; it2 < up + 1; ++it2) {

      Int_t        jTrkID_2 = (vecESDTrks.at(it2)).first;
      AliESDtrack* esdTrk_2 = (vecESDTrks.at(it2)).second;

      fTrackArray->AddAt(esdTrk_2, 1);


        for (Int_t it3 = it2 + 1; it3 < up + 2; ++it3) {

          Int_t        jTrkID_3 = (vecESDTrks.at(it3)).first;
          AliESDtrack* esdTrk_3 = (vecESDTrks.at(it3)).second;

          fTrackArray->AddAt(esdTrk_3, 2);

	  for(Int_t it4= it3 + 1; it4 < up+3; ++it4){

		  Int_t        jTrkID_4 = (vecESDTrks.at(it4)).first;
		  AliESDtrack* esdTrk_4 = (vecESDTrks.at(it4)).second;

		  fTrackArray->AddAt(esdTrk_4, 3);

		  AliAODVertex* secAODVertex = fVtxTagger3Prong->ReconstructSecondaryVertex(fTrackArray, primaryESDVertex, magZkG, vtxRes); 

		  if (secAODVertex) {
		    AliAODTrack* aodTrk_1 = (AliAODTrack*)jet->TrackAt(jTrkID_1, fTrackArrayIn);
		    AliAODTrack* aodTrk_2 = (AliAODTrack*)jet->TrackAt(jTrkID_2, fTrackArrayIn);
		    AliAODTrack* aodTrk_3 = (AliAODTrack*)jet->TrackAt(jTrkID_3, fTrackArrayIn);
		    AliAODTrack* aodTrk_4 = (AliAODTrack*)jet->TrackAt(jTrkID_4, fTrackArrayIn);

		    secAODVertex->AddDaughter(aodTrk_1);
		    secAODVertex->AddDaughter(aodTrk_2);
		    secAODVertex->AddDaughter(aodTrk_3);
		    secAODVertex->AddDaughter(aodTrk_4);

		    if (!fjetCuts3Prong->IsVertexSelected(secAODVertex, aodEvent, magZkG, vtxRes))
		      continue;

		    new ((* arrayVtxHF)[nSecndVxtHF]) AliAODVertex(* secAODVertex);
		    nSecndVxtHF++;
		  } // end if (vert)
	   } //end for it4
        } // end for it3
    } // end for it2
  } // end for it1

    //cout<<"This is a 4prong SV\n";

  fTrackArray->Clear();

  for (vctr_pair_int_esdTrk::iterator it = vecESDTrks.begin(); it != vecESDTrks.end(); ++it) {
      AliESDtrack* lESDtrk = (* it).second;
      if (lESDtrk)
        delete lESDtrk;
  }

  return nSecndVxtHF;
}
//_________________________________________________________________________
void AliAnalysisTaskBJetTC::Terminate(Option_t *)
{  
  
}
