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
#include "AliFemtoModelCorrFctn3DLCMSSym.h"
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
AliFemtoManager* ConfigFemtoAnalysis(int doESD=1, int iCollision=1, int iCent=0) {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  
//   AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
//   Reader->SetConstrained(true);
//   Reader->SetUseTPCOnly(false);

  int multmin=0, multmax=0;
  double zmin, zmax;

  int ppranges[9] = { 1, 12, 17, 23, 30, 37, 45, 58, 150 };
  int aaranges[11] = { 0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900 };

  if (iCollision == 0) { // pp setup
    if ((iCent >=0) && (iCent<8)) {
      multmin = ppranges[iCent];
      multmax = ppranges[iCent+1];
    }
    zmin = -12.0;
    zmax = 12.0;
  }
  if (iCollision == 1) { // PbPb setup
    if ((iCent >=0) && (iCent<10)) {
      multmin = aaranges[iCent];
      multmax = aaranges[iCent+1];
    }
    zmin = -8.0;
    zmax = 8.0;
  }

  AliFemtoEventReaderESDChain* Reader=new AliFemtoEventReaderESDChain();
  Reader->SetConstrained(true);
  Reader->SetUseTPCOnly(false);

  if (iCollision == 0) 
    Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);
  if (iCollision == 1) 
    Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kV0Centrality);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);
	
  // Switches for QA analyses
  int runPiPlusStandard = 1;
  int runPiMinusStandard = 1;
  int runPositiveQA = 1;
  int runNegativeQA = 1;
  int runPositiveTPCQA = 1;
  int runNegativeTPCQA = 0;

  // *** First QA task - standard HBT analysis with all the cut on ***

  // *** Begin pion-pion (positive) analysis ***
  if (runPiPlusStandard) {
    AliFemtoVertexMultAnalysis *anpipstd = new AliFemtoVertexMultAnalysis(1, zmin, zmax, 5, multmin, multmax);
    anpipstd->SetNumEventsToMix(iCollision?3:10);
    anpipstd->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpipstd = new AliFemtoBasicEventCut();
    mecpipstd->SetEventMult(0,100000);
    mecpipstd->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpipstd = new AliFemtoESDTrackCut();
    dtcpipstd->SetPidProbPion(0.2,1.001);
    dtcpipstd->SetPidProbMuon(0.0,1.0);
    dtcpipstd->SetPidProbKaon(0.0,1.0);
    dtcpipstd->SetPidProbProton(0.0,1.0);
    dtcpipstd->SetCharge(1.0);
    dtcpipstd->SetPt(0.15,0.5);
    dtcpipstd->SetMass(PionMass);
    // Track quality cuts
    dtcpipstd->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //  dtcpipstd->SetStatus(AliESDtrack::kTPCrefit);
    dtcpipstd->SetminTPCncls(50);
    dtcpipstd->SetRemoveKinks(kTRUE);
    dtcpipstd->SetLabel(kFALSE);
    dtcpipstd->SetMaxITSChiNdof(2.5);
    dtcpipstd->SetMaxTPCChiNdof(3.0);
    dtcpipstd->SetMaxImpactXY(3.0);
    dtcpipstd->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtpipstd = new AliFemtoCutMonitorParticleYPt("cutPasspipstd", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpipstd = new AliFemtoCutMonitorParticleYPt("cutFailpipstd", 0.13957);
    dtcpipstd->AddCutMonitor(cutPassYPtpipstd, cutFailYPtpipstd);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpipstd = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpipstd->SetShareQualityMax(0.0);
    sqpcpipstd->SetShareFractionMax(0.02);
    sqpcpipstd->SetRemoveSameLabel(kFALSE);
    sqpcpipstd->SetTPCEntranceSepMinimum(2.0);

    anpipstd->SetEventCut(mecpipstd);
    anpipstd->SetFirstParticleCut(dtcpipstd);
    anpipstd->SetSecondParticleCut(dtcpipstd);
    anpipstd->SetPairCut(sqpcpipstd);
    
    AliFemtoShareQualityCorrFctn *csqqinvpipstd= new AliFemtoShareQualityCorrFctn("sqqinvcfpipstd",40,0.0,0.4);
    anpipstd->AddCorrFctn(csqqinvpipstd);

    AliFemtoCorrFctnDirectYlm *cylmpipstd = new AliFemtoCorrFctnDirectYlm("cylmpipstd",3,60,0.0,0.3,1);
    anpipstd->AddCorrFctn(cylmpipstd);
    
    AliFemtoQinvCorrFctn *cqinvpipstd = new AliFemtoQinvCorrFctn("qinvcfpipstd", 100,0.0,1.0);
    anpipstd->AddCorrFctn(cqinvpipstd);

    AliFemtoChi2CorrFctn *cchiqinvpipstd= new AliFemtoChi2CorrFctn("chicfpipstd",40,0.0,0.4);
    anpipstd->AddCorrFctn(cchiqinvpipstd);

    AliFemtoCorrFctnTPCNcls *cqtpcnclspipstd = new AliFemtoCorrFctnTPCNcls("cqtpcnclspipstd",40,0.0,0.4);
    anpipstd->AddCorrFctn(cqtpcnclspipstd);

    AliFemtoCorrFctn3DLCMSSym *cq3dlcmspipstd = new AliFemtoCorrFctn3DLCMSSym("cq3dlcmspipstd",60,0.3);
    anpipstd->AddCorrFctn(cq3dlcmspipstd);

    Manager->AddAnalysis(anpipstd);	
  }
  // *** End pion-pion (positive) analysis

  // *** Begin pion-pion (negative) analysis ***
  if (runPiMinusStandard) {
    AliFemtoVertexMultAnalysis *anpimstd = new AliFemtoVertexMultAnalysis(1, zmin, zmax, 5, multmin, multmax);
    anpimstd->SetNumEventsToMix(iCollision?3:10);
    anpimstd->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpimstd = new AliFemtoBasicEventCut();
    mecpimstd->SetEventMult(0,100000);
    mecpimstd->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpimstd = new AliFemtoESDTrackCut();
    dtcpimstd->SetPidProbPion(0.2,1.001);
    dtcpimstd->SetPidProbMuon(0.0,1.0);
    dtcpimstd->SetPidProbKaon(0.0,1.0);
    dtcpimstd->SetPidProbProton(0.0,1.0);
    dtcpimstd->SetCharge(-1.0);
    dtcpimstd->SetPt(0.15,0.5);
    dtcpimstd->SetMass(PionMass);
    // Track quality cuts
    dtcpimstd->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //  dtcpimstd->SetStatus(AliESDtrack::kTPCrefit);
    dtcpimstd->SetminTPCncls(50);
    dtcpimstd->SetRemoveKinks(kTRUE);
    dtcpimstd->SetLabel(kFALSE);
    dtcpimstd->SetMaxITSChiNdof(2.5);
    dtcpimstd->SetMaxTPCChiNdof(3.0);
    dtcpimstd->SetMaxImpactXY(3.0);
    dtcpimstd->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtpimstd = new AliFemtoCutMonitorParticleYPt("cutPasspimstd", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpimstd = new AliFemtoCutMonitorParticleYPt("cutFailpimstd", 0.13957);
    dtcpimstd->AddCutMonitor(cutPassYPtpimstd, cutFailYPtpimstd);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpimstd = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpimstd->SetShareQualityMax(0.0);
    sqpcpimstd->SetShareFractionMax(0.02);
    sqpcpimstd->SetRemoveSameLabel(kFALSE);
    sqpcpimstd->SetTPCEntranceSepMinimum(2.0);

    anpimstd->SetEventCut(mecpimstd);
    anpimstd->SetFirstParticleCut(dtcpimstd);
    anpimstd->SetSecondParticleCut(dtcpimstd);
    anpimstd->SetPairCut(sqpcpimstd);
    
    AliFemtoShareQualityCorrFctn *csqqinvpimstd= new AliFemtoShareQualityCorrFctn("sqqinvcfpimstd",40,0.0,0.4);
    anpimstd->AddCorrFctn(csqqinvpimstd);

    AliFemtoCorrFctnDirectYlm *cylmpimstd = new AliFemtoCorrFctnDirectYlm("cylmpimstd",3,60,0.0,0.3,1);
    anpimstd->AddCorrFctn(cylmpimstd);
    
    AliFemtoQinvCorrFctn *cqinvpimstd = new AliFemtoQinvCorrFctn("qinvcfpimstd", 100,0.0,1.0);
    anpimstd->AddCorrFctn(cqinvpimstd);

    AliFemtoChi2CorrFctn *cchiqinvpimstd= new AliFemtoChi2CorrFctn("chicfpimstd",40,0.0,0.4);
    anpimstd->AddCorrFctn(cchiqinvpimstd);

    AliFemtoCorrFctnTPCNcls *cqtpcnclspimstd = new AliFemtoCorrFctnTPCNcls("cqtpcnclspimstd",40,0.0,0.4);
    anpimstd->AddCorrFctn(cqtpcnclspimstd);

    AliFemtoCorrFctn3DLCMSSym *cq3dlcmspimstd = new AliFemtoCorrFctn3DLCMSSym("cq3dlcmspimstd",60,0.3);
    anpimstd->AddCorrFctn(cq3dlcmspimstd);

    Manager->AddAnalysis(anpimstd);	
  }
  // *** End pion-pion (negative) analysis

  // *** Second QA task - HBT analysis with all pair cuts off ***
  // *** Begin pion-pion (positive) analysis ***
  if (runPositiveQA) {
    AliFemtoVertexMultAnalysis *anpipnct = new AliFemtoVertexMultAnalysis(1, zmin, zmax, 5, multmin, multmax);
    anpipnct->SetNumEventsToMix(iCollision?3:10);
    anpipnct->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpipnct = new AliFemtoBasicEventCut();
    mecpipnct->SetEventMult(0,100000);
    mecpipnct->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpipnct = new AliFemtoESDTrackCut();
    dtcpipnct->SetPidProbPion(0.0,1.001);
    dtcpipnct->SetPidProbMuon(0.0,1.0);
    dtcpipnct->SetPidProbKaon(0.0,1.0);
    dtcpipnct->SetPidProbProton(0.0,1.0);
    dtcpipnct->SetCharge(1.0);
    dtcpipnct->SetPt(0.1,1.0);
    dtcpipnct->SetMass(PionMass);
    // Track quality cuts
    dtcpipnct->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //  dtcpipnct->SetStatus(AliESDtrack::kTPCrefit);
    dtcpipnct->SetminTPCncls(50);
    dtcpipnct->SetRemoveKinks(kTRUE);
    dtcpipnct->SetLabel(kFALSE);
    dtcpipnct->SetMaxITSChiNdof(6.0);
    dtcpipnct->SetMaxTPCChiNdof(6.0);
    dtcpipnct->SetMaxImpactXY(3.0);
    dtcpipnct->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtpipnct = new AliFemtoCutMonitorParticleYPt("cutPasspipnct", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpipnct = new AliFemtoCutMonitorParticleYPt("cutFailpipnct", 0.13957);
    dtcpipnct->AddCutMonitor(cutPassYPtpipnct, cutFailYPtpipnct);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpipnct = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpipnct->SetShareQualityMax(1.0);
    sqpcpipnct->SetShareFractionMax(1.0);
    sqpcpipnct->SetRemoveSameLabel(kFALSE);
    sqpcpipnct->SetTPCEntranceSepMinimum(0.0);

    anpipnct->SetEventCut(mecpipnct);
    anpipnct->SetFirstParticleCut(dtcpipnct);
    anpipnct->SetSecondParticleCut(dtcpipnct);
    anpipnct->SetPairCut(sqpcpipnct);
    
    AliFemtoShareQualityCorrFctn *csqqinvpipnct= new AliFemtoShareQualityCorrFctn("sqqinvcfpipnct",40,0.0,0.4);
    anpipnct->AddCorrFctn(csqqinvpipnct);

    AliFemtoCorrFctnDirectYlm *cylmpipnct = new AliFemtoCorrFctnDirectYlm("cylmpipnct",3,60,0.0,0.3,1);
    anpipnct->AddCorrFctn(cylmpipnct);
    
    AliFemtoQinvCorrFctn *cqinvpipnct = new AliFemtoQinvCorrFctn("qinvcfpipnct", 100,0.0,1.0);
    anpipnct->AddCorrFctn(cqinvpipnct);

    AliFemtoChi2CorrFctn *cchiqinvpipnct= new AliFemtoChi2CorrFctn("chicfpipnct",40,0.0,0.4);
    anpipnct->AddCorrFctn(cchiqinvpipnct);

    AliFemtoCorrFctnTPCNcls *cqtpcnclspipnct = new AliFemtoCorrFctnTPCNcls("cqtpcnclspipnct",40,0.0,0.4);
    anpipnct->AddCorrFctn(cqtpcnclspipnct);

    Manager->AddAnalysis(anpipnct);	
  }
  // *** End pion-pion (positive) analysis

  // *** Begin pion-pion (negative) analysis ***
  if (runNegativeQA) {
    AliFemtoVertexMultAnalysis *anpimnct = new AliFemtoVertexMultAnalysis(1, zmin, zmax, 5, multmin, multmax);
    anpimnct->SetNumEventsToMix(iCollision?3:10);
    anpimnct->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpimnct = new AliFemtoBasicEventCut();
    mecpimnct->SetEventMult(0,100000);
    mecpimnct->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpimnct = new AliFemtoESDTrackCut();
    dtcpimnct->SetPidProbPion(0.0,1.001);
    dtcpimnct->SetPidProbMuon(0.0,1.0);
    dtcpimnct->SetPidProbKaon(0.0,1.0);
    dtcpimnct->SetPidProbProton(0.0,1.0);
    dtcpimnct->SetCharge(-1.0);
    dtcpimnct->SetPt(0.1,1.0);
    dtcpimnct->SetMass(PionMass);
    // Track quality cuts
    dtcpimnct->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //  dtcpimnct->SetStatus(AliESDtrack::kTPCrefit);
    dtcpimnct->SetminTPCncls(50);
    dtcpimnct->SetRemoveKinks(kTRUE);
    dtcpimnct->SetLabel(kFALSE);
    dtcpimnct->SetMaxITSChiNdof(6.0);
    dtcpimnct->SetMaxTPCChiNdof(6.0);
    dtcpimnct->SetMaxImpactXY(3.0);
    dtcpimnct->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtpimnct = new AliFemtoCutMonitorParticleYPt("cutPasspimnct", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpimnct = new AliFemtoCutMonitorParticleYPt("cutFailpimnct", 0.13957);
    dtcpimnct->AddCutMonitor(cutPassYPtpimnct, cutFailYPtpimnct);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpimnct = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpimnct->SetShareQualityMax(1.0);
    sqpcpimnct->SetShareFractionMax(1.0);
    sqpcpimnct->SetRemoveSameLabel(kFALSE);
    sqpcpimnct->SetTPCEntranceSepMinimum(0.0);

    anpimnct->SetEventCut(mecpimnct);
    anpimnct->SetFirstParticleCut(dtcpimnct);
    anpimnct->SetSecondParticleCut(dtcpimnct);
    anpimnct->SetPairCut(sqpcpimnct);
    
    AliFemtoShareQualityCorrFctn *csqqinvpimnct= new AliFemtoShareQualityCorrFctn("sqqinvcfpimnct",40,0.0,0.4);
    anpimnct->AddCorrFctn(csqqinvpimnct);

    AliFemtoCorrFctnDirectYlm *cylmpimnct = new AliFemtoCorrFctnDirectYlm("cylmpimnct",3,60,0.0,0.3,1);
    anpimnct->AddCorrFctn(cylmpimnct);
    
    AliFemtoQinvCorrFctn *cqinvpimnct = new AliFemtoQinvCorrFctn("qinvcfpimnct", 100,0.0,1.0);
    anpimnct->AddCorrFctn(cqinvpimnct);

    AliFemtoChi2CorrFctn *cchiqinvpimnct= new AliFemtoChi2CorrFctn("chicfpimnct",40,0.0,0.4);
    anpimnct->AddCorrFctn(cchiqinvpimnct);

    AliFemtoCorrFctnTPCNcls *cqtpcnclspimnct = new AliFemtoCorrFctnTPCNcls("cqtpcnclspimnct",40,0.0,0.4);
    anpimnct->AddCorrFctn(cqtpcnclspimnct);

    Manager->AddAnalysis(anpimnct);	
  }
  // *** End pion-pion (negative) analysis

  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin pion-pion (positive) analysis ***
  if (runPositiveTPCQA) {
    AliFemtoVertexMultAnalysis *anpiptpc = new AliFemtoVertexMultAnalysis(1, zmin, zmax, 5, multmin, multmax);
    anpiptpc->SetNumEventsToMix(iCollision?3:10);
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
    dtcpiptpc->SetCharge(0.0);
    dtcpiptpc->SetPt(0.01,1000.0);
    dtcpiptpc->SetMass(PionMass);
    // Track quality cuts
    //dtcpiptpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    dtcpiptpc->SetStatus(AliESDtrack::kTPCrefit);
    dtcpiptpc->SetminTPCncls(50);
    dtcpiptpc->SetRemoveKinks(kTRUE);
    dtcpiptpc->SetLabel(kFALSE);
//     dtcpiptpc->SetMaxITSChiNdof(6.0);
    dtcpiptpc->SetMaxTPCChiNdof(4.0);
//     dtcpiptpc->SetMaxImpactXY(3.0);
//     dtcpiptpc->SetMaxImpactZ(3.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtpiptpc = new AliFemtoCutMonitorParticleYPt("cutPasspiptpc", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpiptpc = new AliFemtoCutMonitorParticleYPt("cutFailpiptpc", 0.13957);
    dtcpiptpc->AddCutMonitor(cutPassYPtpiptpc, cutFailYPtpiptpc);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpiptpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpiptpc->SetShareQualityMax(1.0);
    sqpcpiptpc->SetShareFractionMax(1.0);
    sqpcpiptpc->SetRemoveSameLabel(kFALSE);
    sqpcpiptpc->SetTPCEntranceSepMinimum(0.0);

    anpiptpc->SetEventCut(mecpiptpc);
    anpiptpc->SetFirstParticleCut(dtcpiptpc);
    anpiptpc->SetSecondParticleCut(dtcpiptpc);
    anpiptpc->SetPairCut(sqpcpiptpc);
    
    AliFemtoShareQualityCorrFctn *csqqinvpiptpc= new AliFemtoShareQualityCorrFctn("sqqinvcfpiptpc",40,0.0,0.4);
    anpiptpc->AddCorrFctn(csqqinvpiptpc);

//     AliFemtoCorrFctnDirectYlm *cylmpiptpc = new AliFemtoCorrFctnDirectYlm("cylmpiptpc",3,60,0.0,0.3,1);
//     anpiptpc->AddCorrFctn(cylmpiptpc);
    
    AliFemtoQinvCorrFctn *cqinvpiptpc = new AliFemtoQinvCorrFctn("qinvcfpiptpc", 100,0.0,1.0);
    anpiptpc->AddCorrFctn(cqinvpiptpc);

    AliFemtoChi2CorrFctn *cchiqinvpiptpc= new AliFemtoChi2CorrFctn("chicfpiptpc",40,0.0,0.4);
    anpiptpc->AddCorrFctn(cchiqinvpiptpc);

    AliFemtoCorrFctnTPCNcls *cqtpcnclspiptpc = new AliFemtoCorrFctnTPCNcls("cqtpcnclspiptpc",40,0.0,0.4);
    anpiptpc->AddCorrFctn(cqtpcnclspiptpc);

    Manager->AddAnalysis(anpiptpc);	
  }
  // *** End pion-pion (positive) analysis

  // *** Begin pion-pion (negative) analysis ***
  if (runNegativeTPCQA) {
    AliFemtoVertexMultAnalysis *anpimtpc = new AliFemtoVertexMultAnalysis(1, zmin, zmax, 5, multmin, multmax);
    anpimtpc->SetNumEventsToMix(iCollision?3:10);
    anpimtpc->SetMinSizePartCollection(2);

    AliFemtoBasicEventCut* mecpimtpc = new AliFemtoBasicEventCut();
    mecpimtpc->SetEventMult(0,100000);
    mecpimtpc->SetVertZPos(-1000,1000);
	
    AliFemtoESDTrackCut* dtcpimtpc = new AliFemtoESDTrackCut();
//     dtcpimtpc->SetPidProbPion(0.0,1.001);
//     dtcpimtpc->SetPidProbMuon(0.0,1.0);
//     dtcpimtpc->SetPidProbKaon(0.0,1.0);
//     dtcpimtpc->SetPidProbProton(0.0,1.0);
    dtcpimtpc->SetCharge(0.0);
    dtcpimtpc->SetPt(0.001,100.0);
    dtcpimtpc->SetMass(PionMass);
    // Track quality cuts
    //dtcpimtpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
    //    dtcpimtpc->SetStatus(AliESDtrack::kTPCrefit);
//     dtcpimtpc->SetminTPCncls(50);
    dtcpimtpc->SetRemoveKinks(kTRUE);
    dtcpimtpc->SetLabel(kFALSE);
//     dtcpimtpc->SetMaxITSChiNdof(6.0);
//     dtcpimtpc->SetMaxTPCChiNdof(6.0);
    dtcpimtpc->SetMaxImpactXY(1.0);
    dtcpimtpc->SetMaxImpactZ(1.0);

    AliFemtoCutMonitorParticleYPt *cutPassYPtpimtpc = new AliFemtoCutMonitorParticleYPt("cutPasspimtpc", 0.13957);
    AliFemtoCutMonitorParticleYPt *cutFailYPtpimtpc = new AliFemtoCutMonitorParticleYPt("cutFailpimtpc", 0.13957);
    dtcpimtpc->AddCutMonitor(cutPassYPtpimtpc, cutFailYPtpimtpc);

    AliFemtoShareQualityTPCEntranceSepPairCut *sqpcpimtpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
    sqpcpimtpc->SetShareQualityMax(1.0);
    sqpcpimtpc->SetShareFractionMax(1.0);
    sqpcpimtpc->SetRemoveSameLabel(kFALSE);
    sqpcpimtpc->SetTPCEntranceSepMinimum(0.0);

    anpimtpc->SetEventCut(mecpimtpc);
    anpimtpc->SetFirstParticleCut(dtcpimtpc);
    anpimtpc->SetSecondParticleCut(dtcpimtpc);
    anpimtpc->SetPairCut(sqpcpimtpc);
    
    AliFemtoShareQualityCorrFctn *csqqinvpimtpc= new AliFemtoShareQualityCorrFctn("sqqinvcfpimtpc",40,0.0,0.4);
    anpimtpc->AddCorrFctn(csqqinvpimtpc);

    AliFemtoCorrFctnDirectYlm *cylmpimtpc = new AliFemtoCorrFctnDirectYlm("cylmpimtpc",3,60,0.0,0.3,1);
    anpimtpc->AddCorrFctn(cylmpimtpc);
    
    AliFemtoQinvCorrFctn *cqinvpimtpc = new AliFemtoQinvCorrFctn("qinvcfpimtpc", 100,0.0,1.0);
    anpimtpc->AddCorrFctn(cqinvpimtpc);

    AliFemtoChi2CorrFctn *cchiqinvpimtpc= new AliFemtoChi2CorrFctn("chicfpimtpc",40,0.0,0.4);
    anpimtpc->AddCorrFctn(cchiqinvpimtpc);

    AliFemtoCorrFctnTPCNcls *cqtpcnclspimtpc = new AliFemtoCorrFctnTPCNcls("cqtpcnclspimtpc",40,0.0,0.4);
    anpimtpc->AddCorrFctn(cqtpcnclspimtpc);

    Manager->AddAnalysis(anpimtpc);	
  }
  // *** End pion-pion (negative) analysis



  return Manager;
}                         
                      
