/*********************************************************************
 *                                                                   *
 * ConfigFemtoAnalysis.C - configuration macro for the femtoscopic   *
 * analysis, meant as a QA process for two-particle effects          *
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
#include "AliFemtoCorrFctnDEtaDPhi.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
	
  // Switches for QA analyses
  int runPositiveTPCQA = 1;
  int runNegativeTPCQA = 1;

  int runktdep = 1;
  int run3d = 0;
  
  int runtype = 0; // Types 0 - global, 1 - ITS only
  int isrealdata = 0;

  AliFemtoEventReaderESDChain* Reader=new AliFemtoEventReaderESDChain();
  Reader->SetConstrained(true);
  Reader->SetUseTPCOnly(false);
  if (runtype == 0)
    Reader->SetReadTrackType(AliFemtoEventReaderESDChain::kGlobal);
  else if (runtype == 1)
    Reader->SetReadTrackType(AliFemtoEventReaderESDChain::kITSOnly);
  if (isrealdata)
    Reader->SetUsePhysicsSelection(kTRUE);
  else
    Reader->SetUsePhysicsSelection(kFALSE);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin pion-pion (positive) analysis ***
  if (runPositiveTPCQA) {
    AliFemtoVertexMultAnalysis *anpiptpc = new AliFemtoVertexMultAnalysis(6, -10.0, 10.0, 1, 0, 10000);
    anpiptpc->SetNumEventsToMix(10);
    anpiptpc->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpiptpc = new AliFemtoBasicEventCut();
    mecpiptpc->SetEventMult(0,100000);
    mecpiptpc->SetVertZPos(-1000,1000);

    AliFemtoCutMonitorEventMult *cutPassEvMpiptpc = new AliFemtoCutMonitorEventMult("cutPasspiptpc");
    AliFemtoCutMonitorEventMult *cutFailEvMpiptpc = new AliFemtoCutMonitorEventMult("cutFailpiptpc");
    mecpiptpc->AddCutMonitor(cutPassEvMpiptpc, cutFailEvMpiptpc);

    AliFemtoCutMonitorEventVertex *cutPassEvVpiptpc = new AliFemtoCutMonitorEventVertex("cutPasspiptpc");
    AliFemtoCutMonitorEventVertex *cutFailEvVpiptpc = new AliFemtoCutMonitorEventVertex("cutFailpiptpc");
    mecpiptpc->AddCutMonitor(cutPassEvVpiptpc, cutFailEvVpiptpc);
	
    AliFemtoESDTrackCut* dtcpiptpc = new AliFemtoESDTrackCut();
    //     dtcpiptpc->SetPidProbPion(0.0,1.001);
    //     dtcpiptpc->SetPidProbMuon(0.0,1.0);
    //     dtcpiptpc->SetPidProbKaon(0.0,1.0);
    //     dtcpiptpc->SetPidProbProton(0.0,1.0);
    dtcpiptpc->SetCharge(1.0);
    dtcpiptpc->SetPt(0.05,10.0);
    dtcpiptpc->SetMass(PionMass);
    // Track quality cuts
    if (runtype == 0) {
      dtcpiptpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
      //    dtcpiptpc->SetStatus(AliESDtrack::kTPCrefit);
      //    dtcpiptpc->SetStatus(AliESDtrack::kITSrefit);
      dtcpiptpc->SetminTPCncls(50);
      dtcpiptpc->SetRemoveKinks(kTRUE);
      dtcpiptpc->SetLabel(kFALSE);
      //    dtcpiptpc->SetMaxITSChiNdof(6.0);
      dtcpiptpc->SetMaxTPCChiNdof(6.0);
      dtcpiptpc->SetMaxImpactXY(0.2);
      dtcpiptpc->SetMaxImpactZ(0.25);
      dtcpiptpc->SetMaxSigmaToVertex(6.0);
    }
    else if (runtype == 1) {
      //      dtcpiptpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
      //    dtcpiptpc->SetStatus(AliESDtrack::kTPCrefit);
      dtcpiptpc->SetStatus(AliESDtrack::kITSrefit);
      //      dtcpiptpc->SetminTPCncls(50);
      dtcpiptpc->SetRemoveKinks(kTRUE);
      dtcpiptpc->SetLabel(kFALSE);
      //    dtcpiptpc->SetMaxITSChiNdof(6.0);
      //      dtcpiptpc->SetMaxTPCChiNdof(6.0);
      dtcpiptpc->SetMaxImpactXY(0.2);
      dtcpiptpc->SetMaxImpactZ(0.25);
      //      dtcpiptpc->SetMaxSigmaToVertex(6.0);
    }

    AliFemtoCutMonitorParticleYPt *cutPassYPtpiptpc = new AliFemtoCutMonitorParticleYPt("cutPasspiptpc", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpiptpc = new AliFemtoCutMonitorParticleYPt("cutFailpiptpc", 0.13957);
    dtcpiptpc->AddCutMonitor(cutPassYPtpiptpc, cutFailYPtpiptpc);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpiptpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
    if (runtype == 0) {
      sqpcpiptpc->SetShareQualityMax(0.0);
      sqpcpiptpc->SetShareFractionMax(0.05);
      sqpcpiptpc->SetRemoveSameLabel(kFALSE);
      sqpcpiptpc->SetTPCEntranceSepMinimum(0.0);
    }
    else if (runtype == 1) {
      sqpcpiptpc->SetShareQualityMax(1.0);
      sqpcpiptpc->SetShareFractionMax(1.05);
      sqpcpiptpc->SetRemoveSameLabel(kFALSE);
      sqpcpiptpc->SetTPCEntranceSepMinimum(0.0);
    }

    anpiptpc->SetEventCut(mecpiptpc);
    anpiptpc->SetFirstParticleCut(dtcpiptpc);
    anpiptpc->SetSecondParticleCut(dtcpiptpc);
    anpiptpc->SetPairCut(sqpcpiptpc);
    
    //     AliFemtoShareQualityCorrFctn *csqqinvpiptpc= new AliFemtoShareQualityCorrFctn("sqqinvcfpiptpc",40,0.0,0.4);
    //     anpiptpc->AddCorrFctn(csqqinvpiptpc);

    AliFemtoCorrFctnDirectYlm *cylmpiptpc = new AliFemtoCorrFctnDirectYlm("cylmpiptpc",3,100,0.0,1.5,1);
    anpiptpc->AddCorrFctn(cylmpiptpc);
    
    AliFemtoQinvCorrFctn *cqinvpiptpc = new AliFemtoQinvCorrFctn("qinvcfpiptpc", 100,0.0,1.5);
    anpiptpc->AddCorrFctn(cqinvpiptpc);

    AliFemtoChi2CorrFctn *cchiqinvpiptpc= new AliFemtoChi2CorrFctn("chicfpiptpc",40,0.0,0.4);
    anpiptpc->AddCorrFctn(cchiqinvpiptpc);

    if (run3d) {
      AliFemtoBPLCMS3DCorrFctn *cq3dallpiptpc = new AliFemtoBPLCMS3DCorrFctn("cq3dallpiptpc",100,-1.5,1.5);
      anpiptpc->AddCorrFctn(cq3dallpiptpc);
    }

    if (runktdep) {
      // Intrdouce kT binning
      AliFemtoKTPairCut *ktpairkT1piptpc = new AliFemtoKTPairCut(0.05,0.2);
      AliFemtoKTPairCut *ktpairkT2piptpc = new AliFemtoKTPairCut(0.3,0.4);
      AliFemtoKTPairCut *ktpairkT3piptpc = new AliFemtoKTPairCut(0.4,0.55);
      AliFemtoKTPairCut *ktpairkT4piptpc = new AliFemtoKTPairCut(0.55,2.0);
      
      AliFemtoCorrFctnDirectYlm *cylmkT1piptpc = new AliFemtoCorrFctnDirectYlm("cylmkT1piptpc",3,100,0.0,1.5,1);
      cylmkT1piptpc->SetPairSelectionCut(ktpairkT1piptpc);
      anpiptpc->AddCorrFctn(cylmkT1piptpc);
      
      AliFemtoQinvCorrFctn *cqinvkT1piptpc = new AliFemtoQinvCorrFctn("qinvcfkT1piptpc", 100,0.0,1.5);
      cqinvkT1piptpc->SetPairSelectionCut(ktpairkT1piptpc);
      anpiptpc->AddCorrFctn(cqinvkT1piptpc);
      
      AliFemtoCorrFctnDirectYlm *cylmkT2piptpc = new AliFemtoCorrFctnDirectYlm("cylmkT2piptpc",3,100,0.0,1.5,1);
      cylmkT2piptpc->SetPairSelectionCut(ktpairkT2piptpc);
      anpiptpc->AddCorrFctn(cylmkT2piptpc);
      
      AliFemtoQinvCorrFctn *cqinvkT2piptpc = new AliFemtoQinvCorrFctn("qinvcfkT2piptpc", 100,0.0,1.5);
      cqinvkT2piptpc->SetPairSelectionCut(ktpairkT2piptpc);
      anpiptpc->AddCorrFctn(cqinvkT2piptpc);
      
      AliFemtoCorrFctnDirectYlm *cylmkT3piptpc = new AliFemtoCorrFctnDirectYlm("cylmkT3piptpc",3,100,0.0,1.5,1);
      cylmkT3piptpc->SetPairSelectionCut(ktpairkT3piptpc);
      anpiptpc->AddCorrFctn(cylmkT3piptpc);
      
      AliFemtoQinvCorrFctn *cqinvkT3piptpc = new AliFemtoQinvCorrFctn("qinvcfkT3piptpc", 100,0.0,1.5);
      cqinvkT3piptpc->SetPairSelectionCut(ktpairkT3piptpc);
      anpiptpc->AddCorrFctn(cqinvkT3piptpc);
      
      AliFemtoCorrFctnDirectYlm *cylmkT4piptpc = new AliFemtoCorrFctnDirectYlm("cylmkT4piptpc",3,100,0.0,1.5,1);
      cylmkT4piptpc->SetPairSelectionCut(ktpairkT4piptpc);
      anpiptpc->AddCorrFctn(cylmkT4piptpc);
      
      AliFemtoQinvCorrFctn *cqinvkT4piptpc = new AliFemtoQinvCorrFctn("qinvcfkT4piptpc", 100,0.0,1.5);
      cqinvkT4piptpc->SetPairSelectionCut(ktpairkT4piptpc);
      anpiptpc->AddCorrFctn(cqinvkT4piptpc);
    }
    //     AliFemtoCorrFctnTPCNcls *cqtpcnclspiptpc = new AliFemtoCorrFctnTPCNcls("cqtpcnclspiptpc",40,0.0,0.4);
    //     anpiptpc->AddCorrFctn(cqtpcnclspiptpc);

    //     AliFemtoCorrFctnDEtaDPhi *cdetadphipiptpc = new AliFemtoCorrFctnDEtaDPhi("cdetadphipiptpc", 18, 20);
    //     anpiptpc->AddCorrFctn(cdetadphipiptpc);

    Manager->AddAnalysis(anpiptpc);	
  }
  // *** End pion-pion (positive) analysis

  // *** Begin pion-pion (negative) analysis ***
  if (runNegativeTPCQA) {
    AliFemtoVertexMultAnalysis *anpimtpc = new AliFemtoVertexMultAnalysis(6, -10.0, 10.0, 1, 0, 10000);
    anpimtpc->SetNumEventsToMix(10);
    anpimtpc->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpimtpc = new AliFemtoBasicEventCut();
    mecpimtpc->SetEventMult(0,100000);
    mecpimtpc->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpimtpc = new AliFemtoESDTrackCut();
    //     dtcpimtpc->SetPidProbPion(0.0,1.001);
    //     dtcpimtpc->SetPidProbMuon(0.0,1.0);
    //     dtcpimtpc->SetPidProbKaon(0.0,1.0);
    //     dtcpimtpc->SetPidProbProton(0.0,1.0);
    dtcpimtpc->SetCharge(-1.0);
    dtcpimtpc->SetPt(0.05,10.0);
    dtcpimtpc->SetMass(PionMass);
    // Track quality cuts
    if (runtype == 0) {
      dtcpimtpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
      //    dtcpimtpc->SetStatus(AliESDtrack::kTPCrefit);
      //dtcpimtpc->SetStatus(AliESDtrack::kITSrefit);
      dtcpimtpc->SetminTPCncls(50);
      dtcpimtpc->SetRemoveKinks(kTRUE);
      dtcpimtpc->SetLabel(kFALSE);
      //     dtcpimtpc->SetMaxITSChiNdof(6.0);
      dtcpimtpc->SetMaxTPCChiNdof(6.0);
      dtcpimtpc->SetMaxImpactXY(0.2);
      dtcpimtpc->SetMaxImpactZ(0.25);
      dtcpimtpc->SetMaxSigmaToVertex(6.0);
    }
    else if (runtype == 1) {
      //      dtcpimtpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
      //    dtcpimtpc->SetStatus(AliESDtrack::kTPCrefit);
      dtcpimtpc->SetStatus(AliESDtrack::kITSrefit);
      //    dtcpimtpc->SetminTPCncls(50);
      dtcpimtpc->SetRemoveKinks(kTRUE);
      dtcpimtpc->SetLabel(kFALSE);
      //     dtcpimtpc->SetMaxITSChiNdof(6.0);
      //     dtcpimtpc->SetMaxTPCChiNdof(6.0);
      dtcpimtpc->SetMaxImpactXY(0.2);
      dtcpimtpc->SetMaxImpactZ(0.25);
      //      dtcpimtpc->SetMaxSigmaToVertex(6.0);
    }
    
    AliFemtoCutMonitorParticleYPt *cutPassYPtpimtpc = new AliFemtoCutMonitorParticleYPt("cutPasspimtpc", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpimtpc = new AliFemtoCutMonitorParticleYPt("cutFailpimtpc", 0.13957);
    dtcpimtpc->AddCutMonitor(cutPassYPtpimtpc, cutFailYPtpimtpc);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpimtpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
    if (runtype == 0) {
      sqpcpimtpc->SetShareQualityMax(0.0);
      sqpcpimtpc->SetShareFractionMax(0.05);
      sqpcpimtpc->SetRemoveSameLabel(kFALSE);
      sqpcpimtpc->SetTPCEntranceSepMinimum(0.0);
    }
    else if (runtype == 1) {
      sqpcpimtpc->SetShareQualityMax(1.0);
      sqpcpimtpc->SetShareFractionMax(1.05);
      sqpcpimtpc->SetRemoveSameLabel(kFALSE);
      sqpcpimtpc->SetTPCEntranceSepMinimum(0.0);
    }
 
    anpimtpc->SetEventCut(mecpimtpc);
    anpimtpc->SetFirstParticleCut(dtcpimtpc);
    anpimtpc->SetSecondParticleCut(dtcpimtpc);
    anpimtpc->SetPairCut(sqpcpimtpc);
    
    //     AliFemtoShareQualityCorrFctn *csqqinvpimtpc= new AliFemtoShareQualityCorrFctn("sqqinvcfpimtpc",40,0.0,0.4);
    //     anpimtpc->AddCorrFctn(csqqinvpimtpc);

    AliFemtoCorrFctnDirectYlm *cylmpimtpc = new AliFemtoCorrFctnDirectYlm("cylmpimtpc",3,100,0.0,1.5,1);
    anpimtpc->AddCorrFctn(cylmpimtpc);
    
    AliFemtoQinvCorrFctn *cqinvpimtpc = new AliFemtoQinvCorrFctn("qinvcfpimtpc", 100,0.0,1.5);
    anpimtpc->AddCorrFctn(cqinvpimtpc);

    AliFemtoChi2CorrFctn *cchiqinvpimtpc= new AliFemtoChi2CorrFctn("chicfpimtpc",40,0.0,0.4);
    anpimtpc->AddCorrFctn(cchiqinvpimtpc);

    if (run3d) {
      AliFemtoBPLCMS3DCorrFctn *cq3dallpimtpc = new AliFemtoBPLCMS3DCorrFctn("cq3dallpimtpc",100,-1.5,1.5);
      anpimtpc->AddCorrFctn(cq3dallpimtpc);
    }

    // Intrdouce kT binning
    if (runktdep) {
      // Intrdouce kT binning
      AliFemtoKTPairCut *ktpairkT1pimtpc = new AliFemtoKTPairCut(0.05,0.3);
      AliFemtoKTPairCut *ktpairkT2pimtpc = new AliFemtoKTPairCut(0.3,0.4);
      AliFemtoKTPairCut *ktpairkT3pimtpc = new AliFemtoKTPairCut(0.4,0.55);
      AliFemtoKTPairCut *ktpairkT4pimtpc = new AliFemtoKTPairCut(0.55,2.0);
      
      AliFemtoCorrFctnDirectYlm *cylmkT1pimtpc = new AliFemtoCorrFctnDirectYlm("cylmkT1pimtpc",3,100,0.0,1.5,1);
      cylmkT1pimtpc->SetPairSelectionCut(ktpairkT1pimtpc);
      anpimtpc->AddCorrFctn(cylmkT1pimtpc);
      
      AliFemtoQinvCorrFctn *cqinvkT1pimtpc = new AliFemtoQinvCorrFctn("qinvcfkT1pimtpc", 100,0.0,1.5);
      cqinvkT1pimtpc->SetPairSelectionCut(ktpairkT1pimtpc);
      anpimtpc->AddCorrFctn(cqinvkT1pimtpc);
      
      AliFemtoCorrFctnDirectYlm *cylmkT2pimtpc = new AliFemtoCorrFctnDirectYlm("cylmkT2pimtpc",3,100,0.0,1.5,1);
      cylmkT2pimtpc->SetPairSelectionCut(ktpairkT2pimtpc);
      anpimtpc->AddCorrFctn(cylmkT2pimtpc);
      
      AliFemtoQinvCorrFctn *cqinvkT2pimtpc = new AliFemtoQinvCorrFctn("qinvcfkT2pimtpc", 100,0.0,1.5);
      cqinvkT2pimtpc->SetPairSelectionCut(ktpairkT2pimtpc);
      anpimtpc->AddCorrFctn(cqinvkT2pimtpc);
      
      AliFemtoCorrFctnDirectYlm *cylmkT3pimtpc = new AliFemtoCorrFctnDirectYlm("cylmkT3pimtpc",3,100,0.0,1.5,1);
      cylmkT3pimtpc->SetPairSelectionCut(ktpairkT3pimtpc);
      anpimtpc->AddCorrFctn(cylmkT3pimtpc);
      
      AliFemtoQinvCorrFctn *cqinvkT3pimtpc = new AliFemtoQinvCorrFctn("qinvcfkT3pimtpc", 100,0.0,1.5);
      cqinvkT3pimtpc->SetPairSelectionCut(ktpairkT3pimtpc);
      anpimtpc->AddCorrFctn(cqinvkT3pimtpc);
      
      AliFemtoCorrFctnDirectYlm *cylmkT4pimtpc = new AliFemtoCorrFctnDirectYlm("cylmkT4pimtpc",3,100,0.0,1.5,1);
      cylmkT4pimtpc->SetPairSelectionCut(ktpairkT4pimtpc);
      anpimtpc->AddCorrFctn(cylmkT4pimtpc);
      
      AliFemtoQinvCorrFctn *cqinvkT4pimtpc = new AliFemtoQinvCorrFctn("qinvcfkT4pimtpc", 100,0.0,1.5);
      cqinvkT4pimtpc->SetPairSelectionCut(ktpairkT4pimtpc);
      anpimtpc->AddCorrFctn(cqinvkT4pimtpc);
    }
    //     AliFemtoCorrFctnTPCNcls *cqtpcnclspimtpc = new AliFemtoCorrFctnTPCNcls("cqtpcnclspimtpc",40,0.0,0.4);
    //     anpimtpc->AddCorrFctn(cqtpcnclspimtpc);

    //     AliFemtoCorrFctnDEtaDPhi *cdetadphipimtpc = new AliFemtoCorrFctnDEtaDPhi("cdetadphipimtpc", 18, 20);
    //     anpimtpc->AddCorrFctn(cdetadphipimtpc);

    Manager->AddAnalysis(anpimtpc);	
  }
  // *** End pion-pion (negative) analysis

  return Manager;
}                         
                      
