/*********************************************************************
 *                                                                   *
 * ConfigFemtoAnalysis.C - configuration macro for the femtoscopic   *
 * analysis, to be run in the analysis train.                        *
 * Assumed input data: large (>10M) sample of MC pp events           *
 * Inluded analysis:                                                 *
 *    - positive pion HBT, 3 kt bins, 1D+3D functions                *
 *    - negative pion HBT, 3 kt bins, 1D+3D functions                *
 *    - positive kaon HBT, 1 kt bin,  1D+3D functions                *
 *                                                                   *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)                         *
 *                                                                   *
 *********************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include "AliFemtoModelCorrFctn3DSpherical.h"
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  
  AliFemtoEventReaderESDChain* Reader=new AliFemtoEventReaderESDChain();
  Reader->SetConstrained(true);
  Reader->SetUseTPCOnly(false);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  int runPositivePions = 1;
  int runNegativePions = 1;
  int runPositiveKaons = 1;
  int runNegativeKaons = 1;
  int runPositiveNegativeKaons = 1;

  // Z vertex mixing settings
  int nZVertexBins = 4;
  Double_t minZVertex = -10.0;
  Double_t maxZVertex = 10.0;

  if (runPositivePions) {
    // *** Begin pion-pion (positive) analysis ***
    AliFemtoVertexMultAnalysis *anpip = new AliFemtoVertexMultAnalysis(nZVertexBins, minZVertex, maxZVertex, 1, 2, 200000);
    anpip->SetNumEventsToMix(10);
    anpip->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpip = new AliFemtoBasicEventCut();
    mecpip->SetEventMult(2,100000);
    mecpip->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpip = new AliFemtoESDTrackCut();
    dtcpip->SetPidProbPion(0.2,1.001);
    dtcpip->SetPidProbMuon(0.0,1.0);
    dtcpip->SetPidProbKaon(0.0,1.0);
    dtcpip->SetPidProbProton(0.0,1.0);
    dtcpip->SetCharge(1.0);
    dtcpip->SetPt(0.05,1.0);
    dtcpip->SetRapidity(-0.9,0.9);
    dtcpip->SetMass(PionMass);
    // Track quality cuts
    dtcpip->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //  dtcpip->SetStatus(AliESDtrack::kTPCrefit);
    dtcpip->SetminTPCncls(50);
    dtcpip->SetRemoveKinks(kTRUE);
    dtcpip->SetLabel(kFALSE);
    dtcpip->SetMaxITSChiNdof(2.5);
    dtcpip->SetMaxTPCChiNdof(3.0);
    dtcpip->SetMaxImpactXY(3.0);
    dtcpip->SetMaxImpactZ(3.0);

    // Track monitors
    AliFemtoCutMonitorParticleYPt *cutPassYPtpip = new AliFemtoCutMonitorParticleYPt("cutPasspip", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpip = new AliFemtoCutMonitorParticleYPt("cutFailpip", 0.13957);
    dtcpip->AddCutMonitor(cutPassYPtpip, cutFailYPtpip);

    AliFemtoCutMonitorEventMult *cutPassEvMpip = new AliFemtoCutMonitorEventMult("cutPasspip");
    AliFemtoCutMonitorEventMult *cutFailEvMpip = new AliFemtoCutMonitorEventMult("cutFailpip");
    mecpip->AddCutMonitor(cutPassEvMpip, cutFailEvMpip);

    AliFemtoCutMonitorEventVertex *cutPassEvVpip = new AliFemtoCutMonitorEventVertex("cutPasspip");
    AliFemtoCutMonitorEventVertex *cutFailEvVpip = new AliFemtoCutMonitorEventVertex("cutFailpip");
    mecpip->AddCutMonitor(cutPassEvVpip, cutFailEvVpip);

    // Pair cut
    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpip = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpip->SetShareQualityMax(0.0);
    sqpcpip->SetShareFractionMax(0.02);
    sqpcpip->SetRemoveSameLabel(kFALSE);
    sqpcpip->SetTPCEntranceSepMinimum(2.0);

    anpip->SetEventCut(mecpip);
    anpip->SetFirstParticleCut(dtcpip);
    anpip->SetSecondParticleCut(dtcpip);
    anpip->SetPairCut(sqpcpip);

    // Two-track quality monitoring
    AliFemtoShareQualityCorrFctn *csqqinvpip= new AliFemtoShareQualityCorrFctn("sqqinvcfpip",40,0.0,0.4);
    AliFemtoChi2CorrFctn *cchiqinvpip= new AliFemtoChi2CorrFctn("chicfpip",40,0.0,0.4);
    AliFemtoCorrFctnTPCNcls *cqtpcnclspip = new AliFemtoCorrFctnTPCNcls("cqtpcnclspip",40,0.0,0.4);

    // Intrdouce kT binning
    AliFemtoKTPairCut *ktpairkT1pip = new AliFemtoKTPairCut(0.05,0.27);
    AliFemtoKTPairCut *ktpairkT2pip = new AliFemtoKTPairCut(0.27,0.37);
    AliFemtoKTPairCut *ktpairkT3pip = new AliFemtoKTPairCut(0.37,1.0);

    // Purely experimental correlation function
    AliFemtoCorrFctnDirectYlm *cylmkT1pip = new AliFemtoCorrFctnDirectYlm("cylmkT1pip",3,80,0.0,0.8,1);
    cylmkT1pip->SetPairSelectionCut(ktpairkT1pip);
    anpip->AddCorrFctn(cylmkT1pip);
    
    AliFemtoCorrFctnDirectYlm *cylmkT2pip = new AliFemtoCorrFctnDirectYlm("cylmkT2pip",3,80,0.0,0.8,1);
    cylmkT2pip->SetPairSelectionCut(ktpairkT2pip);
    anpip->AddCorrFctn(cylmkT2pip);
    
    AliFemtoCorrFctnDirectYlm *cylmkT3pip = new AliFemtoCorrFctnDirectYlm("cylmkT3pip",3,80,0.0,0.8,1);
    cylmkT3pip->SetPairSelectionCut(ktpairkT3pip);
    anpip->AddCorrFctn(cylmkT3pip);

    AliFemtoQinvCorrFctn *cqinvkt1pip = new AliFemtoQinvCorrFctn("qinvcfkt1pip", 100,0.0,1.0);
    cqinvkt1pip->SetPairSelectionCut(ktpairkT1pip);
    anpip->AddCorrFctn(cqinvkt1pip);

    AliFemtoQinvCorrFctn *cqinvkt2pip = new AliFemtoQinvCorrFctn("qinvcfkt2pip", 100,0.0,1.0);
    cqinvkt2pip->SetPairSelectionCut(ktpairkT2pip);
    anpip->AddCorrFctn(cqinvkt2pip);

    AliFemtoQinvCorrFctn *cqinvkt3pip = new AliFemtoQinvCorrFctn("qinvcfkt3pip", 100,0.0,1.0);
    cqinvkt3pip->SetPairSelectionCut(ktpairkT3pip);
    anpip->AddCorrFctn(cqinvkt3pip);

    // Add correlation functions to the analysis 
    anpip->AddCorrFctn(csqqinvpip);
    anpip->AddCorrFctn(cchiqinvpip);
    anpip->AddCorrFctn(cqtpcnclspip);

    Manager->AddAnalysis(anpip);	

    // *** End pion-pion (positive) analysis
  }

  if (runNegativePions) {
    // *** Begin pion-pion (negative) analysis ***
    AliFemtoVertexMultAnalysis *anpim = new AliFemtoVertexMultAnalysis(nZVertexBins, minZVertex, maxZVertex, 1, 2, 200000);
    anpim->SetNumEventsToMix(10);
    anpim->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpim = new AliFemtoBasicEventCut();
    mecpim->SetEventMult(2,100000);
    mecpim->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpim = new AliFemtoESDTrackCut();
    dtcpim->SetPidProbPion(0.2,1.001);
    dtcpim->SetPidProbMuon(0.0,1.0);
    dtcpim->SetPidProbKaon(0.0,1.0);
    dtcpim->SetPidProbProton(0.0,1.0);
    dtcpim->SetCharge(1.0);
    dtcpim->SetPt(0.05,1.0);
    dtcpim->SetRapidity(-0.9,0.9);
    dtcpim->SetMass(PionMass);
    // Track quality cuts
    dtcpim->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //  dtcpim->SetStatus(AliESDtrack::kTPCrefit);
    dtcpim->SetminTPCncls(50);
    dtcpim->SetRemoveKinks(kTRUE);
    dtcpim->SetLabel(kFALSE);
    dtcpim->SetMaxITSChiNdof(2.5);
    dtcpim->SetMaxTPCChiNdof(3.0);
    dtcpim->SetMaxImpactXY(3.0);
    dtcpim->SetMaxImpactZ(3.0);

    // Track monitors
    AliFemtoCutMonitorParticleYPt *cutPassYPtpim = new AliFemtoCutMonitorParticleYPt("cutPasspim", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpim = new AliFemtoCutMonitorParticleYPt("cutFailpim", 0.13957);
    dtcpim->AddCutMonitor(cutPassYPtpim, cutFailYPtpim);

    AliFemtoCutMonitorEventMult *cutPassEvMpim = new AliFemtoCutMonitorEventMult("cutPasspim");
    AliFemtoCutMonitorEventMult *cutFailEvMpim = new AliFemtoCutMonitorEventMult("cutFailpim");
    mecpim->AddCutMonitor(cutPassEvMpim, cutFailEvMpim);

    AliFemtoCutMonitorEventVertex *cutPassEvVpim = new AliFemtoCutMonitorEventVertex("cutPasspim");
    AliFemtoCutMonitorEventVertex *cutFailEvVpim = new AliFemtoCutMonitorEventVertex("cutFailpim");
    mecpim->AddCutMonitor(cutPassEvVpim, cutFailEvVpim);

    // Pair cut
    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpim = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpim->SetShareQualityMax(0.0);
    sqpcpim->SetShareFractionMax(0.02);
    sqpcpim->SetRemoveSameLabel(kFALSE);
    sqpcpim->SetTPCEntranceSepMinimum(2.0);

    anpim->SetEventCut(mecpim);
    anpim->SetFirstParticleCut(dtcpim);
    anpim->SetSecondParticleCut(dtcpim);
    anpim->SetPairCut(sqpcpim);

    // Two-track quality monitoring
    AliFemtoShareQualityCorrFctn *csqqinvpim= new AliFemtoShareQualityCorrFctn("sqqinvcfpim",40,0.0,0.4);
    AliFemtoChi2CorrFctn *cchiqinvpim= new AliFemtoChi2CorrFctn("chicfpim",40,0.0,0.4);
    AliFemtoCorrFctnTPCNcls *cqtpcnclspim = new AliFemtoCorrFctnTPCNcls("cqtpcnclspim",40,0.0,0.4);

    // Intrdouce kT binning
    AliFemtoKTPairCut *ktpairkT1pim = new AliFemtoKTPairCut(0.05,0.27);
    AliFemtoKTPairCut *ktpairkT2pim = new AliFemtoKTPairCut(0.27,0.37);
    AliFemtoKTPairCut *ktpairkT3pim = new AliFemtoKTPairCut(0.37,1.0);

    // Purely experimental correlation function
    AliFemtoCorrFctnDirectYlm *cylmkT1pim = new AliFemtoCorrFctnDirectYlm("cylmkT1pim",3,80,0.0,0.8,1);
    cylmkT1pim->SetPairSelectionCut(ktpairkT1pim);
    anpim->AddCorrFctn(cylmkT1pim);
    
    AliFemtoCorrFctnDirectYlm *cylmkT2pim = new AliFemtoCorrFctnDirectYlm("cylmkT2pim",3,80,0.0,0.8,1);
    cylmkT2pim->SetPairSelectionCut(ktpairkT2pim);
    anpim->AddCorrFctn(cylmkT2pim);
    
    AliFemtoCorrFctnDirectYlm *cylmkT3pim = new AliFemtoCorrFctnDirectYlm("cylmkT3pim",3,80,0.0,0.8,1);
    cylmkT3pim->SetPairSelectionCut(ktpairkT3pim);
    anpim->AddCorrFctn(cylmkT3pim);

    AliFemtoQinvCorrFctn *cqinvkt1pim = new AliFemtoQinvCorrFctn("qinvcfkt1pim", 100,0.0,1.0);
    cqinvkt1pim->SetPairSelectionCut(ktpairkT1pim);
    anpim->AddCorrFctn(cqinvkt1pim);

    AliFemtoQinvCorrFctn *cqinvkt2pim = new AliFemtoQinvCorrFctn("qinvcfkt2pim", 100,0.0,1.0);
    cqinvkt2pim->SetPairSelectionCut(ktpairkT2pim);
    anpim->AddCorrFctn(cqinvkt2pim);

    AliFemtoQinvCorrFctn *cqinvkt3pim = new AliFemtoQinvCorrFctn("qinvcfkt3pim", 100,0.0,1.0);
    cqinvkt3pim->SetPairSelectionCut(ktpairkT3pim);
    anpim->AddCorrFctn(cqinvkt3pim);

    // Add correlation functions to the analysis 
    anpim->AddCorrFctn(csqqinvpim);
    anpim->AddCorrFctn(cchiqinvpim);
    anpim->AddCorrFctn(cqtpcnclspim);

    Manager->AddAnalysis(anpim);	

    // *** End pion-pion (negative) analysis
  }

  if (runPositiveKaons) {
    // *** Begin Kaon-Kaon (positive) analysis
    AliFemtoVertexMultAnalysis *ankp = new AliFemtoVertexMultAnalysis(nZVertexBins, minZVertex, maxZVertex, 1, 2, 20000);
    ankp->SetNumEventsToMix(5);
    ankp->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* meckp = new AliFemtoBasicEventCut();
    meckp->SetEventMult(1,100000);
    meckp->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtckp = new AliFemtoESDTrackCut();
    dtckp->SetPidProbKaon(0.7,1.001);
    dtckp->SetPidProbMuon(0.0,0.5);
    dtckp->SetPidProbPion(0.0,0.5);
    dtckp->SetPidProbProton(0.0,0.5);
    dtckp->SetCharge(1.0);
    dtckp->SetMostProbableKaon();
    dtckp->SetMomRangeTOFpidIs(0.6,10000.);
    dtckp->SetPt(0.15,2.0);
    dtckp->SetMass(KaonMass);
    // Track quality cuts
    dtckp->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    dtckp->SetminTPCncls(50);
    dtckp->SetRemoveKinks(kTRUE);
    dtckp->SetLabel(kFALSE);
    dtckp->SetMaxITSChiNdof(2.5);
    dtckp->SetMaxTPCChiNdof(3.0);
    dtckp->SetMaxImpactXY(3.0);
    dtckp->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtkp = new AliFemtoCutMonitorParticleYPt("cutPasskp", 0.493677);
    AliFemtoCutMonitorParticleYPt *cutFailYPtkp = new AliFemtoCutMonitorParticleYPt("cutFailkp", 0.493677);
    dtckp->AddCutMonitor(cutPassYPtkp, cutFailYPtkp);

    AliFemtoCutMonitorParticlePtPDG *cutPassPidkp = new AliFemtoCutMonitorParticlePtPDG("cutPasskp", 0.493677);
    AliFemtoCutMonitorParticlePtPDG *cutFailPidkp = new AliFemtoCutMonitorParticlePtPDG("cutFailkp", 0.493677);
    dtckp->AddCutMonitor(cutPassPidkp, cutFailPidkp);

    AliFemtoCutMonitorParticleMomRes *cutPassMRkp = new AliFemtoCutMonitorParticleMomRes("cutPasskp");
    AliFemtoCutMonitorParticleMomRes *cutFailMRkp = new AliFemtoCutMonitorParticleMomRes("cutFailkp");
    dtckp->AddCutMonitor(cutPassMRkp, cutFailMRkp);

    AliFemtoCutMonitorParticleVertPos *cutPassVPkp = new AliFemtoCutMonitorParticleVertPos("cutPasskp");
    AliFemtoCutMonitorParticleVertPos *cutFailVPkp = new AliFemtoCutMonitorParticleVertPos("cutFailkp");
    dtckp->AddCutMonitor(cutPassVPkp, cutFailVPkp);

    AliFemtoCutMonitorEventMult *cutPassEvMkp = new AliFemtoCutMonitorEventMult("cutPasskp");
    AliFemtoCutMonitorEventMult *cutFailEvMkp = new AliFemtoCutMonitorEventMult("cutFailkp");
    meckp->AddCutMonitor(cutPassEvMkp, cutFailEvMkp);

    AliFemtoCutMonitorEventVertex *cutPassEvVkp = new AliFemtoCutMonitorEventVertex("cutPasskp");
    AliFemtoCutMonitorEventVertex *cutFailEvVkp = new AliFemtoCutMonitorEventVertex("cutFailkp");
    meckp->AddCutMonitor(cutPassEvVkp, cutFailEvVkp);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpckp = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpckp->SetShareQualityMax(0.0);
    sqpckp->SetShareFractionMax(0.02);
    sqpckp->SetRemoveSameLabel(kFALSE);
    sqpckp->SetTPCEntranceSepMinimum(3.0);

    ankp->SetEventCut(meckp);
    ankp->SetFirstParticleCut(dtckp);
    ankp->SetSecondParticleCut(dtckp);
    ankp->SetPairCut(sqpckp);

    AliFemtoQinvCorrFctn *cqinvkp= new AliFemtoQinvCorrFctn("qinvcf",100,0.0,1.0);

    //###
    ankp->AddCorrFctn(cqinvkp);

    Manager->AddAnalysis(ankp);	  

    // *** End Kaon-Kaon (positive) analysis
  }

  if (runNegativeKaons) {
    // *** Begin Kaon-Kaon (negative) analysis
    AliFemtoVertexMultAnalysis *ankm = new AliFemtoVertexMultAnalysis(nZVertexBins, minZVertex, maxZVertex, 1, 2, 20000);
    ankm->SetNumEventsToMix(5);
    ankm->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* meckm = new AliFemtoBasicEventCut();
    meckm->SetEventMult(1,100000);
    meckm->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtckm = new AliFemtoESDTrackCut();
    dtckm->SetPidProbKaon(0.7,1.001);
    dtckm->SetPidProbMuon(0.0,0.5);
    dtckm->SetPidProbPion(0.0,0.5);
    dtckm->SetPidProbProton(0.0,0.5);
    dtckp->SetMomRangeTOFpidIs(0.6,10000.);
    dtckm->SetCharge(-1.0);
    dtckm->SetMostProbableKaon();
    dtckm->SetPt(0.15,2.0);
    dtckm->SetMass(KaonMass);
    // Track quality cuts
    dtckm->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    dtckm->SetminTPCncls(50);
    dtckm->SetRemoveKinks(kTRUE);
    dtckm->SetLabel(kFALSE);
    dtckm->SetMaxITSChiNdof(2.5);
    dtckm->SetMaxTPCChiNdof(3.0);
    dtckm->SetMaxImpactXY(3.0);
    dtckm->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtkm = new AliFemtoCutMonitorParticleYPt("cutPasskm", 0.493677);
    AliFemtoCutMonitorParticleYPt *cutFailYPtkm = new AliFemtoCutMonitorParticleYPt("cutFailkm", 0.493677);
    dtckm->AddCutMonitor(cutPassYPtkm, cutFailYPtkm);

    AliFemtoCutMonitorParticlePtPDG *cutPassPidkm = new AliFemtoCutMonitorParticlePtPDG("cutPasskm", 0.493677);
    AliFemtoCutMonitorParticlePtPDG *cutFailPidkm = new AliFemtoCutMonitorParticlePtPDG("cutFailkm", 0.493677);
    dtckm->AddCutMonitor(cutPassPidkm, cutFailPidkm);

    AliFemtoCutMonitorParticleMomRes *cutPassMRkm = new AliFemtoCutMonitorParticleMomRes("cutPasskm");
    AliFemtoCutMonitorParticleMomRes *cutFailMRkm = new AliFemtoCutMonitorParticleMomRes("cutFailkm");
    dtckm->AddCutMonitor(cutPassMRkm, cutFailMRkm);

    AliFemtoCutMonitorParticleVertPos *cutPassVPkm = new AliFemtoCutMonitorParticleVertPos("cutPasskm");
    AliFemtoCutMonitorParticleVertPos *cutFailVPkm = new AliFemtoCutMonitorParticleVertPos("cutFailkm");
    dtckm->AddCutMonitor(cutPassVPkm, cutFailVPkm);

    AliFemtoCutMonitorEventMult *cutPassEvMkm = new AliFemtoCutMonitorEventMult("cutPasskm");
    AliFemtoCutMonitorEventMult *cutFailEvMkm = new AliFemtoCutMonitorEventMult("cutFailkm");
    meckm->AddCutMonitor(cutPassEvMkm, cutFailEvMkm);

    AliFemtoCutMonitorEventVertex *cutPassEvVkm = new AliFemtoCutMonitorEventVertex("cutPasskm");
    AliFemtoCutMonitorEventVertex *cutFailEvVkm = new AliFemtoCutMonitorEventVertex("cutFailkm");
    meckm->AddCutMonitor(cutPassEvVkm, cutFailEvVkm);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpckm = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpckm->SetShareQualityMax(0.0);
    sqpckm->SetShareFractionMax(0.02);
    sqpckm->SetRemoveSameLabel(kFALSE);
    sqpckm->SetTPCEntranceSepMinimum(3.0);

    ankm->SetEventCut(meckm);
    ankm->SetFirstParticleCut(dtckm);
    ankm->SetSecondParticleCut(dtckm);
    ankm->SetPairCut(sqpckm);

    AliFemtoQinvCorrFctn *cqinvkm= new AliFemtoQinvCorrFctn("qinvcf",100,0.0,1.0);


    //###
    ankm->AddCorrFctn(cqinvkm);
 
    Manager->AddAnalysis(ankm);	  

    // *** End Kaon-Kaon (negative) analysis
  }

  if (runPositiveNegativeKaons) {
    // *** Begin Kaon+Kaon- analysis
    AliFemtoVertexMultAnalysis *ankpkm = new AliFemtoVertexMultAnalysis(nZVertexBins, minZVertex, maxZVertex, 1, 2, 20000);
    ankpkm->SetNumEventsToMix(5);
    ankpkm->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* meckpkm = new AliFemtoBasicEventCut();
    meckpkm->SetEventMult(1,100000);
    meckpkm->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtckp = new AliFemtoESDTrackCut();
    dtckp->SetPidProbKaon(0.7,1.001);
    dtckp->SetPidProbMuon(0.0,0.5);
    dtckp->SetPidProbPion(0.0,0.5);
    dtckp->SetPidProbProton(0.0,0.5);
    dtckp->SetMomRangeTOFpidIs(0.6,10000.);
    dtckp->SetCharge(1.0);
    dtckp->SetMostProbableKaon();
    dtckp->SetPt(0.15,2.0);
    dtckp->SetMass(KaonMass);
    // Track quality cuts
    dtckp->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    dtckp->SetminTPCncls(50);
    dtckp->SetRemoveKinks(kTRUE);
    dtckp->SetLabel(kFALSE);
    dtckp->SetMaxITSChiNdof(2.5);
    dtckp->SetMaxTPCChiNdof(3.0);
    dtckp->SetMaxImpactXY(3.0);
    dtckp->SetMaxImpactZ(3.0);

    AliFemtoESDTrackCut* dtckm = new AliFemtoESDTrackCut();
    dtckm->SetPidProbKaon(0.7,1.001);
    dtckm->SetPidProbMuon(0.0,0.5);
    dtckm->SetPidProbPion(0.0,0.5);
    dtckm->SetPidProbProton(0.0,0.5);
    dtckm->SetMomRangeTOFpidIs(0.6,10000.);
    dtckm->SetCharge(-1.0);
    dtckm->SetMostProbableKaon();
    dtckm->SetPt(0.15,2.0);
    dtckm->SetMass(KaonMass);
    // Track quality cuts
    dtckm->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    dtckm->SetminTPCncls(50);
    dtckm->SetRemoveKinks(kTRUE);
    dtckm->SetLabel(kFALSE);
    dtckm->SetMaxITSChiNdof(2.5);
    dtckm->SetMaxTPCChiNdof(3.0);
    dtckm->SetMaxImpactXY(3.0);
    dtckm->SetMaxImpactZ(3.0);


    AliFemtoCutMonitorParticleYPt *cutPassYPtkp = new AliFemtoCutMonitorParticleYPt("cutPasskp", 0.493677);
    AliFemtoCutMonitorParticleYPt *cutFailYPtkp = new AliFemtoCutMonitorParticleYPt("cutFailkp", 0.493677);
    dtckp->AddCutMonitor(cutPassYPtkp, cutFailYPtkp);

    AliFemtoCutMonitorParticleYPt *cutPassYPtkm = new AliFemtoCutMonitorParticleYPt("cutPasskm", 0.493677);
    AliFemtoCutMonitorParticleYPt *cutFailYPtkm = new AliFemtoCutMonitorParticleYPt("cutFailkm", 0.493677);
    dtckm->AddCutMonitor(cutPassYPtkm, cutFailYPtkm);

    AliFemtoCutMonitorParticlePtPDG *cutPassPidkp = new AliFemtoCutMonitorParticlePtPDG("cutPasskp", 0.493677);
    AliFemtoCutMonitorParticlePtPDG *cutFailPidkp = new AliFemtoCutMonitorParticlePtPDG("cutFailkp", 0.493677);
    dtckp->AddCutMonitor(cutPassPidkp, cutFailPidkp);

    AliFemtoCutMonitorParticlePtPDG *cutPassPidkm = new AliFemtoCutMonitorParticlePtPDG("cutPasskm", 0.493677);
    AliFemtoCutMonitorParticlePtPDG *cutFailPidkm = new AliFemtoCutMonitorParticlePtPDG("cutFailkm", 0.493677);
    dtckm->AddCutMonitor(cutPassPidkm, cutFailPidkm);
  
    AliFemtoCutMonitorParticleMomRes *cutPassMRkp = new AliFemtoCutMonitorParticleMomRes("cutPasskp");
    AliFemtoCutMonitorParticleMomRes *cutFailMRkp = new AliFemtoCutMonitorParticleMomRes("cutFailkp");
    dtckp->AddCutMonitor(cutPassMRkp, cutFailMRkp);

    AliFemtoCutMonitorParticleMomRes *cutPassMRkm = new AliFemtoCutMonitorParticleMomRes("cutPasskm");
    AliFemtoCutMonitorParticleMomRes *cutFailMRkm = new AliFemtoCutMonitorParticleMomRes("cutFailkm");
    dtckm->AddCutMonitor(cutPassMRkm, cutFailMRkm);

    AliFemtoCutMonitorParticleVertPos *cutPassVPkp = new AliFemtoCutMonitorParticleVertPos("cutPasskp");
    AliFemtoCutMonitorParticleVertPos *cutFailVPkp = new AliFemtoCutMonitorParticleVertPos("cutFailkp");
    dtckp->AddCutMonitor(cutPassVPkp, cutFailVPkp);

    AliFemtoCutMonitorParticleVertPos *cutPassVPkm = new AliFemtoCutMonitorParticleVertPos("cutPasskm");
    AliFemtoCutMonitorParticleVertPos *cutFailVPkm = new AliFemtoCutMonitorParticleVertPos("cutFailkm");
    dtckm->AddCutMonitor(cutPassVPkm, cutFailVPkm);

    AliFemtoCutMonitorEventMult *cutPassEvMkpkm = new AliFemtoCutMonitorEventMult("cutPasskpkm");
    AliFemtoCutMonitorEventMult *cutFailEvMkpkm = new AliFemtoCutMonitorEventMult("cutFailkpkm");
    meckpkm->AddCutMonitor(cutPassEvMkpkm, cutFailEvMkpkm);
 
    AliFemtoCutMonitorEventVertex *cutPassEvVkpkm = new AliFemtoCutMonitorEventVertex("cutPasskpkm");
    AliFemtoCutMonitorEventVertex *cutFailEvVkpkm = new AliFemtoCutMonitorEventVertex("cutFailkpkm");
    meckpkm->AddCutMonitor(cutPassEvVkpkm, cutFailEvVkpkm);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpckpkm = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpckpkm->SetShareQualityMax(0.0);
    sqpckpkm->SetShareFractionMax(0.02);
    sqpckpkm->SetRemoveSameLabel(kFALSE);
    sqpckpkm->SetTPCEntranceSepMinimum(3.0);

    ankpkm->SetEventCut(meckpkm);
    ankpkm->SetFirstParticleCut(dtckp);
    ankpkm->SetSecondParticleCut(dtckm);
    ankpkm->SetPairCut(sqpckpkm);

    AliFemtoQinvCorrFctn *cqinvkpkm= new AliFemtoQinvCorrFctn("qinvcf",100,0.0,1.0);

    //###
    ankpkm->AddCorrFctn(cqinvkpkm);
 
    Manager->AddAnalysis(ankpkm);	  

    // *** End Kaon+Kaon-  analysis
  }

  return Manager;
}                         
                      
