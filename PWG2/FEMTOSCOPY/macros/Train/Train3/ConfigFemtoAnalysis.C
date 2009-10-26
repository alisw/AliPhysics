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
  
  AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  Reader->SetConstrained(true);
  Reader->SetUseTPCOnly(false);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  int runPositivePions = 1;
  int runNegativePions = 1;
  int runPositiveKaons = 1;
  int runNegativeKaons = 1;
  int runPositiveNegativeKaons = 1;

  if (runPositivePions) {
    // *** Begin pion-pion (positive) analysis ***
    AliFemtoVertexMultAnalysis *anpip = new AliFemtoVertexMultAnalysis(3, -15.6, 15.6, 1, 2, 200000);
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
    dtcpip->SetPt(0.15,0.5);
    dtcpip->SetRapidity(-0.8,0.8);
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
    AliFemtoKTPairCut *ktpairkT1pip = new AliFemtoKTPairCut(0.1,0.27);
    AliFemtoKTPairCut *ktpairkT2pip = new AliFemtoKTPairCut(0.27,0.37);
    AliFemtoKTPairCut *ktpairkT3pip = new AliFemtoKTPairCut(0.37,0.52);

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

    // Setting up the model calculation
    // First create the freeze-out generator
    AliFemtoModelGausLCMSFreezeOutGenerator *tFreezepip = new AliFemtoModelGausLCMSFreezeOutGenerator();
    tFreezepip->SetSizeOut(1.8*TMath::Sqrt(2.0));                                                
    tFreezepip->SetSizeSide(1.3*TMath::Sqrt(2.0));                                               
    tFreezepip->SetSizeLong(1.6*TMath::Sqrt(2.0));                                                

    // And the weight generator                                                                    
    AliFemtoModelWeightGeneratorBasic *tWeightpip = new AliFemtoModelWeightGeneratorBasic();
    tWeightpip->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionPlus());

    // Create a manager that will connect it                                                                                                               
    AliFemtoModelManager *tModelManagerpip = new AliFemtoModelManager();
    tModelManagerpip->AcceptFreezeOutGenerator(tFreezepip);
    tModelManagerpip->AcceptWeightGenerator(tWeightpip);
    tModelManagerpip->CreateCopyHiddenInfo(kFALSE);

    // Model correlation functions
    AliFemtoModelCorrFctn *c1dpipip;
    AliFemtoModelBPLCMSCorrFctn *c3dsmallkt1pip;
    AliFemtoModelBPLCMSCorrFctn *c3dsmallkt2pip;
    AliFemtoModelBPLCMSCorrFctn *c3dsmallkt3pip;

    c1dpipip = new AliFemtoModelCorrFctn("c1dpipip",100,0.0,1.0);
    c1dpipip->ConnectToManager(tModelManagerpip);

    c3dsmallkt1pip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt1pip",40, -0.4, 0.4);
    c3dsmallkt1pip->SetSpecificPairCut(ktpairkT1pip);
    c3dsmallkt1pip->ConnectToManager(tModelManagerpip);

    c3dsmallkt2pip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt2pip",40, -0.4, 0.4);
    c3dsmallkt2pip->SetSpecificPairCut(ktpairkT2pip);
    c3dsmallkt2pip->ConnectToManager(tModelManagerpip);

    c3dsmallkt3pip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt3pip",40, -0.4, 0.4);
    c3dsmallkt3pip->SetSpecificPairCut(ktpairkT3pip);
    c3dsmallkt3pip->ConnectToManager(tModelManagerpip);

    AliFemtoModelCorrFctnDirectYlm *cmylmkt1pip = new AliFemtoModelCorrFctnDirectYlm("mcylmkt1pip",3,80,0.0,0.8,1);
    cmylmkt1pip->SetPairSelectionCut(ktpairkT1pip);
    cmylmkt1pip->ConnectToManager(tModelManagerpip);

    AliFemtoModelCorrFctnDirectYlm *cmylmkt2pip = new AliFemtoModelCorrFctnDirectYlm("mcylmkt2pip",3,80,0.0,0.8,1);
    cmylmkt2pip->SetPairSelectionCut(ktpairkT2pip);
    cmylmkt2pip->ConnectToManager(tModelManagerpip);

    AliFemtoModelCorrFctnDirectYlm *cmylmkt3pip = new AliFemtoModelCorrFctnDirectYlm("mcylmkt3pip",3,80,0.0,0.8,1);
    cmylmkt3pip->SetPairSelectionCut(ktpairkT3pip);
    cmylmkt3pip->ConnectToManager(tModelManagerpip);

    // Add correlation functions to the analysis 
    anpip->AddCorrFctn(csqqinvpip);
    anpip->AddCorrFctn(cchiqinvpip);
    anpip->AddCorrFctn(cqtpcnclspip);
    anpip->AddCorrFctn(c3dsmallkt1pip);
    anpip->AddCorrFctn(c3dsmallkt2pip);
    anpip->AddCorrFctn(c3dsmallkt3pip);
    anpip->AddCorrFctn(c1dpipip);
    anpip->AddCorrFctn(cmylmkt1pip);
    anpip->AddCorrFctn(cmylmkt2pip);
    anpip->AddCorrFctn(cmylmkt3pip);

    Manager->AddAnalysis(anpip);	

    // *** End pion-pion (positive) analysis
  }

  if (runNegativePions) {
    // *** Begin pion-pion (negative) analysis ***
    AliFemtoVertexMultAnalysis *anpim = new AliFemtoVertexMultAnalysis(3, -15.6, 15.6, 1, 2, 200000);
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
    dtcpim->SetPt(0.15,0.5);
    dtcpim->SetRapidity(-0.8,0.8);
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
    AliFemtoKTPairCut *ktpairkT1pim = new AliFemtoKTPairCut(0.1,0.27);
    AliFemtoKTPairCut *ktpairkT2pim = new AliFemtoKTPairCut(0.27,0.37);
    AliFemtoKTPairCut *ktpairkT3pim = new AliFemtoKTPairCut(0.37,0.52);

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

    // Setting up the model calculation
    // First create the freeze-out generator
    AliFemtoModelGausLCMSFreezeOutGenerator *tFreezepim = new AliFemtoModelGausLCMSFreezeOutGenerator();
    tFreezepim->SetSizeOut(1.8*TMath::Sqrt(2.0));                                                
    tFreezepim->SetSizeSide(1.3*TMath::Sqrt(2.0));                                               
    tFreezepim->SetSizeLong(1.6*TMath::Sqrt(2.0));                                                

    // And the weight generator                                                                    
    AliFemtoModelWeightGeneratorBasic *tWeightpim = new AliFemtoModelWeightGeneratorBasic();
    tWeightpim->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionPlus());

    // Create a manager that will connect it                                                                                                               
    AliFemtoModelManager *tModelManagerpim = new AliFemtoModelManager();
    tModelManagerpim->AcceptFreezeOutGenerator(tFreezepim);
    tModelManagerpim->AcceptWeightGenerator(tWeightpim);
    tModelManagerpim->CreateCopyHiddenInfo(kFALSE);

    // Model correlation functions
    AliFemtoModelCorrFctn *c1dpipim;
    AliFemtoModelBPLCMSCorrFctn *c3dsmallkt1pim;
    AliFemtoModelBPLCMSCorrFctn *c3dsmallkt2pim;
    AliFemtoModelBPLCMSCorrFctn *c3dsmallkt3pim;

    c1dpipim = new AliFemtoModelCorrFctn("c1dpipim",100,0.0,1.0);
    c1dpipim->ConnectToManager(tModelManagerpim);

    c3dsmallkt1pim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt1pim",40, -0.4, 0.4);
    c3dsmallkt1pim->SetSpecificPairCut(ktpairkT1pim);
    c3dsmallkt1pim->ConnectToManager(tModelManagerpim);

    c3dsmallkt2pim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt2pim",40, -0.4, 0.4);
    c3dsmallkt2pim->SetSpecificPairCut(ktpairkT2pim);
    c3dsmallkt2pim->ConnectToManager(tModelManagerpim);

    c3dsmallkt3pim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt3pim",40, -0.4, 0.4);
    c3dsmallkt3pim->SetSpecificPairCut(ktpairkT3pim);
    c3dsmallkt3pim->ConnectToManager(tModelManagerpim);

    AliFemtoModelCorrFctnDirectYlm *cmylmkt1pim = new AliFemtoModelCorrFctnDirectYlm("mcylmkt1pim",3,80,0.0,0.8,1);
    cmylmkt1pim->SetPairSelectionCut(ktpairkT1pim);
    cmylmkt1pim->ConnectToManager(tModelManagerpim);

    AliFemtoModelCorrFctnDirectYlm *cmylmkt2pim = new AliFemtoModelCorrFctnDirectYlm("mcylmkt2pim",3,80,0.0,0.8,1);
    cmylmkt2pim->SetPairSelectionCut(ktpairkT2pim);
    cmylmkt2pim->ConnectToManager(tModelManagerpim);

    AliFemtoModelCorrFctnDirectYlm *cmylmkt3pim = new AliFemtoModelCorrFctnDirectYlm("mcylmkt3pim",3,80,0.0,0.8,1);
    cmylmkt3pim->SetPairSelectionCut(ktpairkT3pim);
    cmylmkt3pim->ConnectToManager(tModelManagerpim);

    // Add correlation functions to the analysis 
    anpim->AddCorrFctn(csqqinvpim);
    anpim->AddCorrFctn(cchiqinvpim);
    anpim->AddCorrFctn(cqtpcnclspim);
    anpim->AddCorrFctn(c3dsmallkt1pim);
    anpim->AddCorrFctn(c3dsmallkt2pim);
    anpim->AddCorrFctn(c3dsmallkt3pim);
    anpim->AddCorrFctn(c1dpipim);
    anpim->AddCorrFctn(cmylmkt1pim);
    anpim->AddCorrFctn(cmylmkt2pim);
    anpim->AddCorrFctn(cmylmkt3pim);

    Manager->AddAnalysis(anpim);	

    // *** End pion-pion (negative) analysis
  }

  if (runPositiveKaons) {
    // *** Begin Kaon-Kaon (positive) analysis
    AliFemtoVertexMultAnalysis *ankp = new AliFemtoVertexMultAnalysis(18, -15.6, 15.6, 1, 2, 20000);
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

    AliFemtoModelBPLCMSCorrFctn *c3dsmallkp;
    AliFemtoModelCorrFctn *c1dpikp;

    // Setting up the model calculation
    // First create the freeze-out generator
    AliFemtoModelGausRinvFreezeOutGenerator *tFreezekp = new AliFemtoModelGausRinvFreezeOutGenerator();
    tFreezekp->SetSizeInv(1.8*TMath::Sqrt(2.0));
    tFreezekp->SetSelectPrimaryFromHidden(false);

    // And the weight generator                                                                    
    //   AliFemtoModelWeightGeneratorBasic *tWeightkp = new AliFemtoModelWeightGeneratorBasic();
    //   tWeightkp->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
    AliFemtoModelWeightGeneratorLednicky *tWeightkp = new AliFemtoModelWeightGeneratorLednicky();
    tWeightkp->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
    tWeightkp->SetCoulOn();
    tWeightkp->SetQuantumOn();
    tWeightkp->SetStrongOff();
    tWeightkp->Set3BodyOff();

    // Create a manager that will connect it  
    AliFemtoModelManager *tModelManagerkp = new AliFemtoModelManager();
    tModelManagerkp->AcceptFreezeOutGenerator(tFreezekp);
    tModelManagerkp->AcceptWeightGenerator(tWeightkp);
    tModelManagerkp->CreateCopyHiddenInfo(kFALSE);

    c3dsmallkp = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkp",30, 0.0, 0.6);
    c3dsmallkp->ConnectToManager(tModelManagerkp);

    c1dpikp = new AliFemtoModelCorrFctn("c1dpikp",100,0.0,1.0);
    c1dpikp->ConnectToManager(tModelManagerkp);

    //###
    ankp->AddCorrFctn(cqinvkp);
 
    ankp->AddCorrFctn(c3dsmallkp);
    ankp->AddCorrFctn(c1dpikp);

    Manager->AddAnalysis(ankp);	  

    // *** End Kaon-Kaon (positive) analysis
  }

  if (runNegativeKaons) {
    // *** Begin Kaon-Kaon (negative) analysis
    AliFemtoVertexMultAnalysis *ankm = new AliFemtoVertexMultAnalysis(18, -15.6, 15.6, 1, 2, 20000);
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

    AliFemtoModelBPLCMSCorrFctn *c3dsmallkm;
    AliFemtoModelCorrFctn *c1dpikm;

    // Setting up the model calculation
    // First create the freeze-out generator
    AliFemtoModelGausRinvFreezeOutGenerator *tFreezekm = new AliFemtoModelGausRinvFreezeOutGenerator();
    tFreezekm->SetSizeInv(1.8*TMath::Sqrt(2.0));
    tFreezekm->SetSelectPrimaryFromHidden(false);

    // And the weight generator                                                                    
    //   AliFemtoModelWeightGeneratorBasic *tWeightkp = new AliFemtoModelWeightGeneratorBasic();
    //   tWeightkp->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
    AliFemtoModelWeightGeneratorLednicky *tWeightkm = new AliFemtoModelWeightGeneratorLednicky();
    tWeightkm->SetPairType(AliFemtoModelWeightGenerator::KaonMinusKaonMinus());
    tWeightkm->SetCoulOn();
    tWeightkm->SetQuantumOn();
    tWeightkm->SetStrongOff();
    tWeightkm->Set3BodyOff();

    // Create a manager that will connect it  
    AliFemtoModelManager *tModelManagerkm = new AliFemtoModelManager();
    tModelManagerkm->AcceptFreezeOutGenerator(tFreezekp);
    tModelManagerkm->AcceptWeightGenerator(tWeightkp);
    tModelManagerkm->CreateCopyHiddenInfo(kFALSE);

    c3dsmallkm = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkm",30, 0.0, 0.6);
    c3dsmallkm->ConnectToManager(tModelManagerkm);

    c1dpikm = new AliFemtoModelCorrFctn("c1dpikm",100,0.0,1.0);
    c1dpikm->ConnectToManager(tModelManagerkm);

    //###
    ankm->AddCorrFctn(cqinvkm);
 
    ankm->AddCorrFctn(c3dsmallkm);
    ankm->AddCorrFctn(c1dpikm);

    Manager->AddAnalysis(ankm);	  

    // *** End Kaon-Kaon (negative) analysis
  }

  if (runPositiveNegativeKaons) {
    // *** Begin Kaon+Kaon- analysis
    AliFemtoVertexMultAnalysis *ankpkm = new AliFemtoVertexMultAnalysis(18, -15.6, 15.6, 1, 2, 20000);
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

    AliFemtoModelBPLCMSCorrFctn *c3dsmallkpkm;
    AliFemtoModelCorrFctn *c1dpikpkm;

    // Setting up the model calculation
    // First create the freeze-out generator
    AliFemtoModelGausRinvFreezeOutGenerator *tFreezekpkm = new AliFemtoModelGausRinvFreezeOutGenerator();
    tFreezekpkm->SetSizeInv(1.8*TMath::Sqrt(2.0));
    tFreezekpkm->SetSelectPrimaryFromHidden(false);

    // And the weight generator                                                                    
    AliFemtoModelWeightGeneratorLednicky *tWeightkpkm = new AliFemtoModelWeightGeneratorLednicky();
    tWeightkpkm->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
    tWeightkpkm->SetCoulOn();
    tWeightkpkm->SetQuantumOn();
    tWeightkpkm->SetStrongOff();
    tWeightkpkm->Set3BodyOff();

    // Create a manager that will connect it  
    AliFemtoModelManager *tModelManagerkpkm = new AliFemtoModelManager();
    tModelManagerkpkm->AcceptFreezeOutGenerator(tFreezekpkm);
    tModelManagerkpkm->AcceptWeightGenerator(tWeightkpkm);
    tModelManagerkpkm->CreateCopyHiddenInfo(kFALSE);

    c3dsmallkpkm = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkpkm",30, 0.0, 0.6);
    c3dsmallkpkm->ConnectToManager(tModelManagerkpkm);

    c1dpikpkm = new AliFemtoModelCorrFctn("c1dpikpkm",100,0.0,1.0);
    c1dpikpkm->ConnectToManager(tModelManagerkpkm);

    //###
    ankpkm->AddCorrFctn(cqinvkpkm);
 
    ankpkm->AddCorrFctn(c3dsmallkpkm);
    ankpkm->AddCorrFctn(c1dpikpkm);

    Manager->AddAnalysis(ankpkm);	  

    // *** End Kaon+Kaon-  analysis
  }
  //   if (runPositiveKaons) {
  //     // *** Begin Kaon-Kaon (positive) analysis
  //     AliFemtoVertexMultAnalysis *ankp = new AliFemtoVertexMultAnalysis(10, -15.6, 15.6, 1, 2, 20000);
  //     ankp->SetNumEventsToMix(10);
  //     ankp->SetMinSizePartCollection(2);

  //     AliFemtoBasicEventCut* meckp = new AliFemtoBasicEventCut();
  //     meckp->SetEventMult(1,100000);
  //     meckp->SetVertZPos(-1000,1000);
	
  //     AliFemtoESDTrackCut* dtckp = new AliFemtoESDTrackCut();
  //     dtckp->SetPidProbKaon(0.7,1.001);
  //     dtckp->SetPidProbMuon(0.0,0.5);
  //     dtckp->SetPidProbPion(0.0,0.5);
  //     dtckp->SetPidProbProton(0.0,0.5);
  //     dtckp->SetCharge(1.0);
  //     dtckp->SetMostProbableKaon();
  //     dtckp->SetPt(0.15,2.0);
  //     dtckp->SetMass(KaonMass);
  //     dtckp->SetRapidity(-0.8,0.8);
  //     // Track quality cuts
  //     dtckp->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  //     //dtckp->SetStatus(AliESDtrack::kTPCrefit);
  //     //  dtckp->SetminTPCclsF(95);
  //     dtckp->SetminTPCncls(50);
  //     dtckp->SetRemoveKinks(kTRUE);
  //     dtckp->SetLabel(kFALSE);
  //     dtckp->SetMaxITSChiNdof(2.5);
  //     dtckp->SetMaxTPCChiNdof(3.0);
  //     //  dtckp->SetMaxSigmaToVertex(3.0);
  //     dtckp->SetMaxImpactXY(3.0);
  //     dtckp->SetMaxImpactZ(3.0);

  //     AliFemtoCutMonitorParticleYPt *cutPassYPtkp = new AliFemtoCutMonitorParticleYPt("cutPasskp", 0.493677);
  //     AliFemtoCutMonitorParticleYPt *cutFailYPtkp = new AliFemtoCutMonitorParticleYPt("cutFailkp", 0.493677);
  //     dtckp->AddCutMonitor(cutPassYPtkp, cutFailYPtkp);

  //     AliFemtoCutMonitorParticlePtPDG *cutPassPidkp = new AliFemtoCutMonitorParticlePtPDG("cutPasskp", 0.493677);
  //     AliFemtoCutMonitorParticlePtPDG *cutFailPidkp = new AliFemtoCutMonitorParticlePtPDG("cutFailkp", 0.493677);
  //     dtckp->AddCutMonitor(cutPassPidkp, cutFailPidkp);

  //     AliFemtoCutMonitorParticleMomRes *cutPassMRkp = new AliFemtoCutMonitorParticleMomRes("cutPasskp");
  //     AliFemtoCutMonitorParticleMomRes *cutFailMRkp = new AliFemtoCutMonitorParticleMomRes("cutFailkp");
  //     dtckp->AddCutMonitor(cutPassMRkp, cutFailMRkp);

  //     AliFemtoCutMonitorParticleVertPos *cutPassVPkp = new AliFemtoCutMonitorParticleVertPos("cutPasskp");
  //     AliFemtoCutMonitorParticleVertPos *cutFailVPkp = new AliFemtoCutMonitorParticleVertPos("cutFailkp");
  //     dtckp->AddCutMonitor(cutPassVPkp, cutFailVPkp);

  //     AliFemtoCutMonitorEventMult *cutPassEvMkp = new AliFemtoCutMonitorEventMult("cutPasskp");
  //     AliFemtoCutMonitorEventMult *cutFailEvMkp = new AliFemtoCutMonitorEventMult("cutFailkp");
  //     meckp->AddCutMonitor(cutPassEvMkp, cutFailEvMkp);

  //     AliFemtoCutMonitorEventVertex *cutPassEvVkp = new AliFemtoCutMonitorEventVertex("cutPasskp");
  //     AliFemtoCutMonitorEventVertex *cutFailEvVkp = new AliFemtoCutMonitorEventVertex("cutFailkp");
  //     meckp->AddCutMonitor(cutPassEvVkp, cutFailEvVkp);

  //     AliFemtoShareQualityTPCEntranceSepPairCut *sqpckp = new AliFemtoShareQualityTPCEntranceSepPairCut();
  //     sqpckp->SetShareQualityMax(0.0);
  //     sqpckp->SetShareFractionMax(0.02);
  //     sqpckp->SetRemoveSameLabel(kFALSE);
  //     sqpckp->SetTPCEntranceSepMinimum(3.0);

  //     ankp->SetEventCut(meckp);
  //     ankp->SetFirstParticleCut(dtckp);
  //     ankp->SetSecondParticleCut(dtckp);
  //     ankp->SetPairCut(sqpckp);

  //     AliFemtoQinvCorrFctn *cqinvkp= new AliFemtoQinvCorrFctn("qinvcf",100,0.0,1.0);

  //     AliFemtoModelBPLCMSCorrFctn *c3dsmallkp;
  //     AliFemtoModelCorrFctn *c1dpikp;

  //     // Setting up the model calculation
  //     // First create the freeze-out generator
  //     AliFemtoModelGausRinvFreezeOutGenerator *tFreezekp = new AliFemtoModelGausRinvFreezeOutGenerator();
  //     tFreezekp->SetSizeInv(1.8*TMath::Sqrt(2.0));
  //     tFreezekp->SetSelectPrimaryFromHidden(false);

  //     // And the weight generator                                                                    
  //     //   AliFemtoModelWeightGeneratorBasic *tWeightkp = new AliFemtoModelWeightGeneratorBasic();
  //     //   tWeightkp->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
  //     AliFemtoModelWeightGeneratorLednicky *tWeightkp = new AliFemtoModelWeightGeneratorLednicky();
  //     tWeightkp->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
  //     tWeightkp->SetCoulOn();
  //     tWeightkp->SetQuantumOn();
  //     tWeightkp->SetStrongOff();
  //     tWeightkp->Set3BodyOff();

  //     // Create a manager that will connect it  
  //     AliFemtoModelManager *tModelManagerkp = new AliFemtoModelManager();
  //     tModelManagerkp->AcceptFreezeOutGenerator(tFreezekp);
  //     tModelManagerkp->AcceptWeightGenerator(tWeightkp);
  //     tModelManagerkp->CreateCopyHiddenInfo(kFALSE);

  //     c3dsmallkp = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkp",30, 0.0, 0.6);
  //     c3dsmallkp->ConnectToManager(tModelManagerkp);

  //     c1dpikp = new AliFemtoModelCorrFctn("c1dpikp",100,0.0,1.0);
  //     c1dpikp->ConnectToManager(tModelManagerkp);

  //     //###
  //     ankp->AddCorrFctn(cqinvkp);
 
  //     ankp->AddCorrFctn(c3dsmallkp);
  //     ankp->AddCorrFctn(c1dpikp);

  //     Manager->AddAnalysis(ankp);	  

  //     // *** End Kaon-Kaon (positive) analysis
  //   }

  //   if (runNegativeKaons) {
  //     // *** Begin Kaon-Kaon (negative) analysis
  //     AliFemtoVertexMultAnalysis *ankm = new AliFemtoVertexMultAnalysis(10, -15.6, 15.6, 1, 2, 20000);
  //     ankm->SetNumEventsToMix(10);
  //     ankm->SetMinSizePartCollection(2);

  //     AliFemtoBasicEventCut* meckm = new AliFemtoBasicEventCut();
  //     meckm->SetEventMult(1,100000);
  //     meckm->SetVertZPos(-1000,1000);
	
  //     AliFemtoESDTrackCut* dtckm = new AliFemtoESDTrackCut();
  //     dtckm->SetPidProbKaon(0.7,1.001);
  //     dtckm->SetPidProbMuon(0.0,0.5);
  //     dtckm->SetPidProbPion(0.0,0.5);
  //     dtckm->SetPidProbProton(0.0,0.5);
  //     dtckm->SetCharge(1.0);
  //     dtckm->SetMostProbableKaon();
  //     dtckm->SetPt(0.15,2.0);
  //     dtckm->SetMass(KaonMass);
  //     dtckm->SetRapidity(-0.8,0.8);
  //     // Track quality cuts
  //     dtckm->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  //     //dtckm->SetStatus(AliESDtrack::kTPCrefit);
  //     //  dtckm->SetminTPCclsF(95);
  //     dtckm->SetminTPCncls(50);
  //     dtckm->SetRemoveKinks(kTRUE);
  //     dtckm->SetLabel(kFALSE);
  //     dtckm->SetMaxITSChiNdof(2.5);
  //     dtckm->SetMaxTPCChiNdof(3.0);
  //     //  dtckm->SetMaxSigmaToVertex(3.0);
  //     dtckm->SetMaxImpactXY(3.0);
  //     dtckm->SetMaxImpactZ(3.0);

  //     AliFemtoCutMonitorParticleYPt *cutPassYPtkm = new AliFemtoCutMonitorParticleYPt("cutPasskm", 0.493677);
  //     AliFemtoCutMonitorParticleYPt *cutFailYPtkm = new AliFemtoCutMonitorParticleYPt("cutFailkm", 0.493677);
  //     dtckm->AddCutMonitor(cutPassYPtkm, cutFailYPtkm);

  //     AliFemtoCutMonitorParticlePtPDG *cutPassPidkm = new AliFemtoCutMonitorParticlePtPDG("cutPasskm", 0.493677);
  //     AliFemtoCutMonitorParticlePtPDG *cutFailPidkm = new AliFemtoCutMonitorParticlePtPDG("cutFailkm", 0.493677);
  //     dtckm->AddCutMonitor(cutPassPidkm, cutFailPidkm);

  //     AliFemtoCutMonitorParticleMomRes *cutPassMRkm = new AliFemtoCutMonitorParticleMomRes("cutPasskm");
  //     AliFemtoCutMonitorParticleMomRes *cutFailMRkm = new AliFemtoCutMonitorParticleMomRes("cutFailkm");
  //     dtckm->AddCutMonitor(cutPassMRkm, cutFailMRkm);

  //     AliFemtoCutMonitorParticleVertPos *cutPassVPkm = new AliFemtoCutMonitorParticleVertPos("cutPasskm");
  //     AliFemtoCutMonitorParticleVertPos *cutFailVPkm = new AliFemtoCutMonitorParticleVertPos("cutFailkm");
  //     dtckm->AddCutMonitor(cutPassVPkm, cutFailVPkm);

  //     AliFemtoCutMonitorEventMult *cutPassEvMkm = new AliFemtoCutMonitorEventMult("cutPasskm");
  //     AliFemtoCutMonitorEventMult *cutFailEvMkm = new AliFemtoCutMonitorEventMult("cutFailkm");
  //     meckm->AddCutMonitor(cutPassEvMkm, cutFailEvMkm);

  //     AliFemtoCutMonitorEventVertex *cutPassEvVkm = new AliFemtoCutMonitorEventVertex("cutPasskm");
  //     AliFemtoCutMonitorEventVertex *cutFailEvVkm = new AliFemtoCutMonitorEventVertex("cutFailkm");
  //     meckm->AddCutMonitor(cutPassEvVkm, cutFailEvVkm);

  //     AliFemtoShareQualityTPCEntranceSepPairCut *sqpckm = new AliFemtoShareQualityTPCEntranceSepPairCut();
  //     sqpckm->SetShareQualityMax(0.0);
  //     sqpckm->SetShareFractionMax(0.02);
  //     sqpckm->SetRemoveSameLabel(kFALSE);
  //     sqpckm->SetTPCEntranceSepMinimum(3.0);

  //     ankm->SetEventCut(meckm);
  //     ankm->SetFirstParticleCut(dtckm);
  //     ankm->SetSecondParticleCut(dtckm);
  //     ankm->SetPairCut(sqpckm);

  //     AliFemtoQinvCorrFctn *cqinvkm= new AliFemtoQinvCorrFctn("qinvcf",100,0.0,1.0);

  //     AliFemtoModelBPLCMSCorrFctn *c3dsmallkm;
  //     AliFemtoModelCorrFctn *c1dpikm;

  //     // Setting up the model calculation
  //     // First create the freeze-out generator
  //     AliFemtoModelGausRinvFreezeOutGenerator *tFreezekm = new AliFemtoModelGausRinvFreezeOutGenerator();
  //     tFreezekm->SetSizeInv(1.8*TMath::Sqrt(2.0));
  //     tFreezekm->SetSelectPrimaryFromHidden(false);

  //     // And the weight generator                                                                    
  //     //   AliFemtoModelWeightGeneratorBasic *tWeightkm = new AliFemtoModelWeightGeneratorBasic();
  //     //   tWeightkm->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
  //     AliFemtoModelWeightGeneratorLednicky *tWeightkm = new AliFemtoModelWeightGeneratorLednicky();
  //     tWeightkm->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
  //     tWeightkm->SetCoulOn();
  //     tWeightkm->SetQuantumOn();
  //     tWeightkm->SetStrongOff();
  //     tWeightkm->Set3BodyOff();

  //     // Create a manager that will connect it  
  //     AliFemtoModelManager *tModelManagerkm = new AliFemtoModelManager();
  //     tModelManagerkm->AcceptFreezeOutGenerator(tFreezekm);
  //     tModelManagerkm->AcceptWeightGenerator(tWeightkm);
  //     tModelManagerkm->CreateCopyHiddenInfo(kFALSE);

  //     c3dsmallkm = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkm",30, 0.0, 0.6);
  //     c3dsmallkm->ConnectToManager(tModelManagerkm);

  //     c1dpikm = new AliFemtoModelCorrFctn("c1dpikm",100,0.0,1.0);
  //     c1dpikm->ConnectToManager(tModelManagerkm);

  //     //###
  //     ankm->AddCorrFctn(cqinvkm);
 
  //     ankm->AddCorrFctn(c3dsmallkm);
  //     ankm->AddCorrFctn(c1dpikm);

  //     Manager->AddAnalysis(ankm);	  

  //     // *** End Kaon-Kaon (positive) analysis
  //   }

  return Manager;
}                         
                      
