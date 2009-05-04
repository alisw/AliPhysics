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

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  
  AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  Reader->SetConstrained(true);
  Reader->SetUseTPCOnly(false);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);
	
  // *** Begin pion-pion (positive) analysis ***

  AliFemtoVertexMultAnalysis *anpip = new AliFemtoVertexMultAnalysis(3, -15.6, 15.6, 5, 2, 200);
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

  AliFemtoCutMonitorParticleYPt *cutPassYPtpip = new AliFemtoCutMonitorParticleYPt("cutPasspip", 0.13957);
  AliFemtoCutMonitorParticleYPt *cutFailYPtpip = new AliFemtoCutMonitorParticleYPt("cutFailpip", 0.13957);
  dtcpip->AddCutMonitor(cutPassYPtpip, cutFailYPtpip);

  AliFemtoCutMonitorEventMult *cutPassEvMpip = new AliFemtoCutMonitorEventMult("cutPasspip");
  AliFemtoCutMonitorEventMult *cutFailEvMpip = new AliFemtoCutMonitorEventMult("cutFailpip");
  mecpip->AddCutMonitor(cutPassEvMpip, cutFailEvMpip);

  AliFemtoCutMonitorEventVertex *cutPassEvVpip = new AliFemtoCutMonitorEventVertex("cutPasspip");
  AliFemtoCutMonitorEventVertex *cutFailEvVpip = new AliFemtoCutMonitorEventVertex("cutFailpip");
  mecpip->AddCutMonitor(cutPassEvVpip, cutFailEvVpip);

  AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpip = new AliFemtoShareQualityTPCEntranceSepPairCut();
  sqpcpip->SetShareQualityMax(0.0);
  sqpcpip->SetShareFractionMax(0.02);
  sqpcpip->SetRemoveSameLabel(kFALSE);
  sqpcpip->SetTPCEntranceSepMinimum(2.0);

  anpip->SetEventCut(mecpip);
  anpip->SetFirstParticleCut(dtcpip);
  anpip->SetSecondParticleCut(dtcpip);
  anpip->SetPairCut(sqpcpip);

  AliFemtoCorrFctn3DSpherical *cqsphpip = new AliFemtoCorrFctn3DSpherical("cqsphpip",60,0.0,0.3, 12, 12);

  AliFemtoShareQualityCorrFctn *csqqinvpip= new AliFemtoShareQualityCorrFctn("sqqinvcfpip",40,0.0,0.4);

  // Intrdouce kT binning
  AliFemtoKTPairCut *ktpairkT1pip = new AliFemtoKTPairCut(0.1,0.27);
  AliFemtoKTPairCut *ktpairkT2pip = new AliFemtoKTPairCut(0.27,0.37);
  AliFemtoKTPairCut *ktpairkT3pip = new AliFemtoKTPairCut(0.37,0.52);

  AliFemtoCorrFctnDirectYlm *cylmkT1pip = new AliFemtoCorrFctnDirectYlm("cylmkT1pip",3,60,0.0,0.3);
  cylmkT1pip->SetPairSelectionCut(ktpairkT1pip);
  anpip->AddCorrFctn(cylmkT1pip);
    
  AliFemtoCorrFctnDirectYlm *cylmkT2pip = new AliFemtoCorrFctnDirectYlm("cylmkT2pip",3,60,0.0,0.3);
  cylmkT2pip->SetPairSelectionCut(ktpairkT2pip);
  anpip->AddCorrFctn(cylmkT2pip);
    
  AliFemtoCorrFctnDirectYlm *cylmkT3pip = new AliFemtoCorrFctnDirectYlm("cylmkT3pip",3,60,0.0,0.3);
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

  AliFemtoChi2CorrFctn *cchiqinvpip= new AliFemtoChi2CorrFctn("chicfpip",40,0.0,0.4);
  AliFemtoCorrFctnTPCNcls *cqtpcnclspip = new AliFemtoCorrFctnTPCNcls("cqtpcnclspip",40,0.0,0.4);

  AliFemtoModelBPLCMSCorrFctn *c3dsmallkt1pip;
  AliFemtoModelBPLCMSCorrFctn *c3dsmallkt2pip;
  AliFemtoModelBPLCMSCorrFctn *c3dsmallkt3pip;
  AliFemtoModelCorrFctn *c1dpipip;
  AliFemtoModelCorrFctn3DSpherical *c3dmsphpip;

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

  c3dsmallkt1pip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt2pip",30, 0.0, 0.3);
  c3dsmallkt1pip->SetSpecificPairCut(ktpairkT2pip);
  c3dsmallkt1pip->ConnectToManager(tModelManagerpip);

  c3dsmallkt2pip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt1pip",30, 0.0, 0.3);
  c3dsmallkt2pip->SetSpecificPairCut(ktpairkT1pip);
  c3dsmallkt2pip->ConnectToManager(tModelManagerpip);

  c3dsmallkt3pip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt3pip",30, 0.0, 0.3);
  c3dsmallkt3pip->SetSpecificPairCut(ktpairkT3pip);
  c3dsmallkt3pip->ConnectToManager(tModelManagerpip);

  c1dpipip = new AliFemtoModelCorrFctn("c1dpipip",100,0.0,1.0);
  c1dpipip->ConnectToManager(tModelManagerpip);

  c3dmsphpip = new AliFemtoModelCorrFctn3DSpherical("c3dmsphpip",60, 0.0, 0.3, 12,12);
  c3dmsphpip->ConnectToManager(tModelManagerpip);

  AliFemtoModelCorrFctnDirectYlm *cmylmkt1pip = new AliFemtoModelCorrFctnDirectYlm("mcylmkt1pip",3,60,0.0,0.3);
  cmylmkt1pip->SetPairSelectionCut(ktpairkT1pip);
  cmylmkt1pip->ConnectToManager(tModelManagerpip);
  anpip->AddCorrFctn(cmylmkt1pip);

  AliFemtoModelCorrFctnDirectYlm *cmylmkt2pip = new AliFemtoModelCorrFctnDirectYlm("mcylmkt2pip",3,60,0.0,0.3);
  cmylmkt2pip->SetPairSelectionCut(ktpairkT2pip);
  cmylmkt2pip->ConnectToManager(tModelManagerpip);
  anpip->AddCorrFctn(cmylmkt2pip);

  AliFemtoModelCorrFctnDirectYlm *cmylmkt3pip = new AliFemtoModelCorrFctnDirectYlm("mcylmkt3pip",3,60,0.0,0.3);
  cmylmkt3pip->SetPairSelectionCut(ktpairkT3pip);
  cmylmkt3pip->ConnectToManager(tModelManagerpip);
  anpip->AddCorrFctn(cmylmkt3pip);

  //###
  anpip->AddCorrFctn(cqsphpip);
 
  anpip->AddCorrFctn(csqqinvpip);
  anpip->AddCorrFctn(cchiqinvpip);
  anpip->AddCorrFctn(cqtpcnclspip);
  anpip->AddCorrFctn(c3dsmallkt1pip);
  anpip->AddCorrFctn(c3dsmallkt2pip);
  anpip->AddCorrFctn(c3dsmallkt3pip);
  anpip->AddCorrFctn(c1dpipip);
  anpip->AddCorrFctn(c3dmsphpip);

  Manager->AddAnalysis(anpip);	

  // *** End pion-pion (positive) analysis

  // *** Begin pion-pion (negative) analysis ***

  AliFemtoVertexMultAnalysis *anpim = new AliFemtoVertexMultAnalysis(3, -15.6, 15.6, 5, 2, 200);
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

  AliFemtoCutMonitorParticleYPt *cutPassYPtpim = new AliFemtoCutMonitorParticleYPt("cutPasspim", 0.13957);
  AliFemtoCutMonitorParticleYPt *cutFailYPtpim = new AliFemtoCutMonitorParticleYPt("cutFailpim", 0.13957);
  dtcpim->AddCutMonitor(cutPassYPtpim, cutFailYPtpim);

  AliFemtoCutMonitorEventMult *cutPassEvMpim = new AliFemtoCutMonitorEventMult("cutPasspim");
  AliFemtoCutMonitorEventMult *cutFailEvMpim = new AliFemtoCutMonitorEventMult("cutFailpim");
  mecpim->AddCutMonitor(cutPassEvMpim, cutFailEvMpim);

  AliFemtoCutMonitorEventVertex *cutPassEvVpim = new AliFemtoCutMonitorEventVertex("cutPasspim");
  AliFemtoCutMonitorEventVertex *cutFailEvVpim = new AliFemtoCutMonitorEventVertex("cutFailpim");
  mecpim->AddCutMonitor(cutPassEvVpim, cutFailEvVpim);

  AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpim = new AliFemtoShareQualityTPCEntranceSepPairCut();
  sqpcpim->SetShareQualityMax(0.0);
  sqpcpim->SetShareFractionMax(0.02);
  sqpcpim->SetRemoveSameLabel(kFALSE);
  sqpcpim->SetTPCEntranceSepMinimum(2.0);

  anpim->SetEventCut(mecpim);
  anpim->SetFirstParticleCut(dtcpim);
  anpim->SetSecondParticleCut(dtcpim);
  anpim->SetPairCut(sqpcpim);

  AliFemtoCorrFctn3DSpherical *cqsphpim = new AliFemtoCorrFctn3DSpherical("cqsphpim",60,0.0,0.3, 12, 12);

  AliFemtoShareQualityCorrFctn *csqqinvpim= new AliFemtoShareQualityCorrFctn("sqqinvcfpim",40,0.0,0.4);

  // Intrdouce kT binning
  AliFemtoKTPairCut *ktpairkT1pim = new AliFemtoKTPairCut(0.1,0.27);
  AliFemtoKTPairCut *ktpairkT2pim = new AliFemtoKTPairCut(0.27,0.37);
  AliFemtoKTPairCut *ktpairkT3pim = new AliFemtoKTPairCut(0.37,0.52);

  AliFemtoCorrFctnDirectYlm *cylmkT1pim = new AliFemtoCorrFctnDirectYlm("cylmkT1pim",3,60,0.0,0.3);
  cylmkT1pim->SetPairSelectionCut(ktpairkT1pim);
  anpim->AddCorrFctn(cylmkT1pim);
    
  AliFemtoCorrFctnDirectYlm *cylmkT2pim = new AliFemtoCorrFctnDirectYlm("cylmkT2pim",3,60,0.0,0.3);
  cylmkT2pim->SetPairSelectionCut(ktpairkT2pim);
  anpim->AddCorrFctn(cylmkT2pim);
    
  AliFemtoCorrFctnDirectYlm *cylmkT3pim = new AliFemtoCorrFctnDirectYlm("cylmkT3pim",3,60,0.0,0.3);
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

  AliFemtoChi2CorrFctn *cchiqinvpim= new AliFemtoChi2CorrFctn("chicfpim",40,0.0,0.4);
  AliFemtoCorrFctnTPCNcls *cqtpcnclspim = new AliFemtoCorrFctnTPCNcls("cqtpcnclspim",40,0.0,0.4);

  AliFemtoModelBPLCMSCorrFctn *c3dsmallkt1pim;
  AliFemtoModelBPLCMSCorrFctn *c3dsmallkt2pim;
  AliFemtoModelBPLCMSCorrFctn *c3dsmallkt3pim;
  AliFemtoModelCorrFctn *c1dpimip;
  AliFemtoModelCorrFctn3DSpherical *c3dmsphpim;

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

  c3dsmallkt1pim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt2pim",30, 0.0, 0.3);
  c3dsmallkt1pim->SetSpecificPairCut(ktpairkT2pim);
  c3dsmallkt1pim->ConnectToManager(tModelManagerpim);

  c3dsmallkt2pim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt1pim",30, 0.0, 0.3);
  c3dsmallkt2pim->SetSpecificPairCut(ktpairkT1pim);
  c3dsmallkt2pim->ConnectToManager(tModelManagerpim);

  c3dsmallkt3pim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallkt3pim",30, 0.0, 0.3);
  c3dsmallkt3pim->SetSpecificPairCut(ktpairkT3pim);
  c3dsmallkt3pim->ConnectToManager(tModelManagerpim);

  c1dpimip = new AliFemtoModelCorrFctn("c1dpimip",100,0.0,1.0);
  c1dpimip->ConnectToManager(tModelManagerpim);

  c3dmsphpim = new AliFemtoModelCorrFctn3DSpherical("c3dmsphpim",60, 0.0, 0.3, 12,12);
  c3dmsphpim->ConnectToManager(tModelManagerpim);

  AliFemtoModelCorrFctnDirectYlm *cmylmkt1pim = new AliFemtoModelCorrFctnDirectYlm("mcylmkt1pim",3,60,0.0,0.3);
  cmylmkt1pim->SetPairSelectionCut(ktpairkT1pim);
  cmylmkt1pim->ConnectToManager(tModelManagerpim);
  anpim->AddCorrFctn(cmylmkt1pim);

  AliFemtoModelCorrFctnDirectYlm *cmylmkt2pim = new AliFemtoModelCorrFctnDirectYlm("mcylmkt2pim",3,60,0.0,0.3);
  cmylmkt2pim->SetPairSelectionCut(ktpairkT2pim);
  cmylmkt2pim->ConnectToManager(tModelManagerpim);
  anpim->AddCorrFctn(cmylmkt2pim);

  AliFemtoModelCorrFctnDirectYlm *cmylmkt3pim = new AliFemtoModelCorrFctnDirectYlm("mcylmkt3pim",3,60,0.0,0.3);
  cmylmkt3pim->SetPairSelectionCut(ktpairkT3pim);
  cmylmkt3pim->ConnectToManager(tModelManagerpim);
  anpim->AddCorrFctn(cmylmkt3pim);

  //###
  anpim->AddCorrFctn(cqsphpim);
 
  anpim->AddCorrFctn(csqqinvpim);
  anpim->AddCorrFctn(cchiqinvpim);
  anpim->AddCorrFctn(cqtpcnclspim);
  anpim->AddCorrFctn(c3dsmallkt1pim);
  anpim->AddCorrFctn(c3dsmallkt2pim);
  anpim->AddCorrFctn(c3dsmallkt3pim);
  anpim->AddCorrFctn(c1dpimip);
  anpim->AddCorrFctn(c3dmsphpim);

  Manager->AddAnalysis(anpim);	

  // *** End pion-pion (negative) analysis

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
  dtckp->SetPt(0.15,2.0);
  dtckp->SetMass(KaonMass);
  // Track quality cuts
  dtckp->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  //dtckp->SetStatus(AliESDtrack::kTPCrefit);
  //  dtckp->SetminTPCclsF(95);
  dtckp->SetminTPCncls(50);
  dtckp->SetRemoveKinks(kTRUE);
  dtckp->SetLabel(kFALSE);
  dtckp->SetMaxITSChiNdof(2.5);
  dtckp->SetMaxTPCChiNdof(3.0);
  //  dtckp->SetMaxSigmaToVertex(3.0);
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
  AliFemtoModelWeightGeneratorBasic *tWeightkp = new AliFemtoModelWeightGeneratorBasic();
  tWeightkp->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());

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

  return Manager;
}                         
                      
