#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
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

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  
//   AliFemtoEventReaderESDChain* Reader=new AliFemtoEventReaderESDChain();
//   Reader->SetConstrained(true);
//   Reader->SetReadTPCInner(true);
  AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  Reader->SetConstrained(true);
  //  Reader->SetReadTPCInner(true);
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
  //  dtcpip->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  dtcpip->SetStatus(AliESDtrack::kTPCrefit);
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

//   AliFemtoCutMonitorParticleMomRes *cutPassMR = new AliFemtoCutMonitorParticleMomRes("cutPass");
//   AliFemtoCutMonitorParticleMomRes *cutFailMR = new AliFemtoCutMonitorParticleMomRes("cutFail");
//   dtcpip->AddCutMonitor(cutPassMR, cutFailMR);

//   AliFemtoCutMonitorParticleVertPos *cutPassVP = new AliFemtoCutMonitorParticleVertPos("VPcutPass");
//   AliFemtoCutMonitorParticleVertPos *cutFailVP = new AliFemtoCutMonitorParticleVertPos("VPcutFail");
//   dtcpip->AddCutMonitor(cutPassVP, cutFailVP);

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

  AliFemtoQinvCorrFctn *cqinvpip= new AliFemtoQinvCorrFctn("qinvcfpip", 100,0.0,1.0);
  AliFemtoCorrFctn3DSpherical *cqsphpip = new AliFemtoCorrFctn3DSpherical("cqsphpip",60,0.0,0.3, 12, 12);

  AliFemtoShareQualityCorrFctn *csqqinvpip= new AliFemtoShareQualityCorrFctn("sqqinvcfpip",40,0.0,0.4);
  AliFemtoCorrFctnDirectYlm *cylmpip = new AliFemtoCorrFctnDirectYlm("cylmpip",2,60,0.0,0.3);
    
  AliFemtoChi2CorrFctn *cchiqinvpip= new AliFemtoChi2CorrFctn("chicfpip",40,0.0,0.4);
  AliFemtoCorrFctnTPCNcls *cqtpcnclspip = new AliFemtoCorrFctnTPCNcls("cqtpcnclspip",40,0.0,0.4);

  AliFemtoModelBPLCMSCorrFctn *c3dsmallpip;
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

  c3dsmallpip = new AliFemtoModelBPLCMSCorrFctn("c3dsmallpip",30, 0.0, 0.3);
  c3dsmallpip->ConnectToManager(tModelManagerpip);

  c1dpipip = new AliFemtoModelCorrFctn("c1dpipip",100,0.0,1.0);
  c1dpipip->ConnectToManager(tModelManagerpip);

  c3dmsphpip = new AliFemtoModelCorrFctn3DSpherical("c3dmsphpip",60, 0.0, 0.3, 12,12);
  c3dmsphpip->ConnectToManager(tModelManagerpip);

  AliFemtoModelCorrFctnDirectYlm *cmylmpip = new AliFemtoModelCorrFctnDirectYlm("mcylmpip",2,60,0.0,0.3);
  cmylmpip->ConnectToManager(tModelManagerpip);

  //###
  anpip->AddCorrFctn(cqinvpip);
  anpip->AddCorrFctn(cqsphpip);
  anpip->AddCorrFctn(cylmpip);
 
  anpip->AddCorrFctn(csqqinvpip);
  anpip->AddCorrFctn(cchiqinvpip);
  anpip->AddCorrFctn(cqtpcnclspip);
  anpip->AddCorrFctn(c3dsmallpip);
  anpip->AddCorrFctn(c1dpipip);
  anpip->AddCorrFctn(c3dmsphpip);
  anpip->AddCorrFctn(cmylmpip);

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
  dtcpim->SetCharge(-1.0);
  dtcpim->SetPt(0.15,0.5);
  dtcpim->SetMass(PionMass);
  // Track quality cuts
  //  dtcpim->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  dtcpim->SetStatus(AliESDtrack::kTPCrefit);
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

//   AliFemtoCutMonitorParticleMomRes *cutPassMR = new AliFemtoCutMonitorParticleMomRes("cutPass");
//   AliFemtoCutMonitorParticleMomRes *cutFailMR = new AliFemtoCutMonitorParticleMomRes("cutFail");
//   dtcpim->AddCutMonitor(cutPassMR, cutFailMR);

//   AliFemtoCutMonitorParticleVertPos *cutPassVP = new AliFemtoCutMonitorParticleVertPos("VPcutPass");
//   AliFemtoCutMonitorParticleVertPos *cutFailVP = new AliFemtoCutMonitorParticleVertPos("VPcutFail");
//   dtcpim->AddCutMonitor(cutPassVP, cutFailVP);

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

  AliFemtoQinvCorrFctn *cqinvpim= new AliFemtoQinvCorrFctn("qinvcfpim", 100,0.0,1.0);
  AliFemtoCorrFctn3DSpherical *cqsphpim = new AliFemtoCorrFctn3DSpherical("cqsphpim",60,0.0,0.3, 12, 12);

  AliFemtoShareQualityCorrFctn *csqqinvpim= new AliFemtoShareQualityCorrFctn("sqqinvcfpim",40,0.0,0.4);
  AliFemtoCorrFctnDirectYlm *cylmpim = new AliFemtoCorrFctnDirectYlm("cylmpim",2,60,0.0,0.3);
    
  AliFemtoChi2CorrFctn *cchiqinvpim= new AliFemtoChi2CorrFctn("chicfpim",40,0.0,0.4);
  AliFemtoCorrFctnTPCNcls *cqtpcnclspim = new AliFemtoCorrFctnTPCNcls("cqtpcnclspim",40,0.0,0.4);

  AliFemtoModelBPLCMSCorrFctn *c3dsmallpim;
  AliFemtoModelCorrFctn *c1dpipim;
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

  c3dsmallpim = new AliFemtoModelBPLCMSCorrFctn("c3dsmallpim",30, 0.0, 0.3);
  c3dsmallpim->ConnectToManager(tModelManagerpim);

  c1dpipim = new AliFemtoModelCorrFctn("c1dpipim",100,0.0,1.0);
  c1dpipim->ConnectToManager(tModelManagerpim);

  c3dmsphpim = new AliFemtoModelCorrFctn3DSpherical("c3dmsphpim",60, 0.0, 0.3, 12,12);
  c3dmsphpim->ConnectToManager(tModelManagerpim);

  AliFemtoModelCorrFctnDirectYlm *cmylmpim = new AliFemtoModelCorrFctnDirectYlm("mcylmpim",2,60,0.0,0.3);
  cmylmpim->ConnectToManager(tModelManagerpim);

  //###
  anpim->AddCorrFctn(cqinvpim);
  anpim->AddCorrFctn(cqsphpim);
  anpim->AddCorrFctn(cylmpim);
 
  anpim->AddCorrFctn(csqqinvpim);
  anpim->AddCorrFctn(cchiqinvpim);
  anpim->AddCorrFctn(cqtpcnclspim);
  anpim->AddCorrFctn(c3dsmallpim);
  anpim->AddCorrFctn(c1dpipim);
  anpim->AddCorrFctn(c3dmsphpim);
  anpim->AddCorrFctn(cmylmpim);

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
  //  dtckp->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  dtckp->SetStatus(AliESDtrack::kTPCrefit);
  //  dtckp->SetminTPCclsF(95);
  dtckp->SetminTPCncls(50);
  dtckp->SetRemoveKinks(kTRUE);
  dtckp->SetLabel(kFALSE);
  dtckp->SetMaxITSChiNdof(2.5);
  dtckp->SetMaxTPCChiNdof(3.0);
  dtckp->SetMaxSigmaToVertex(3.0);

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
                      
