
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
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
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
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoCutMonitorCollections.h"
#endif
 
//_  6 fm  pure kaons _(0.2-0.4)______________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  double ProtonMass = 0.938272013;

  //multiplicity bins
  int runmults[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int multbins[11] = {0.01, 200000, 400, 600, 900, 950, 500, 600, 700, 800, 900};

  int runch[3] = {1, 1, 0};
  const char *chrgs[3] = { "kk", "akak", "kak" };

  int runktdep = 1;
 // double ktrng[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
  double ktrng[4] = {0.2, 0.4, 0.6, 1.0};

  int run3d = 0; // Do 3D cartesian analysis?
  int runshlcms = 1;

  int nbinssh = 100;
  double shqmax = 0.5;


  AliFemtoEventReaderESDChainKine *Reader = new AliFemtoEventReaderESDChainKine();
  Reader->SetKaonAnalysis(kTRUE);
  //Reader->SetUseMultiplicity(AliFemtoEventReaderESDChainKine::kCentrality);
  //Reader->SetReadTrackType(AliFemtoEventReaderESDChainKine::kGlobal);

  // AliFemtoEventReaderKinematicsChain *Reader = new AliFemtoEventReaderKinematicsChain();
  // Reader->SetUseMultiplicity(AliFemtoEventReaderKinematicsChain::kGlobalCount);
  // // // Reader->SetReadTrackType(AliFemtoEventReaderESDChainKine::kGlobal);

  AliFemtoModelGausRinvFreezeOutGenerator *tFreeze = new AliFemtoModelGausRinvFreezeOutGenerator();
  tFreeze->SetSizeInv(5.0);

 // AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
  // tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
   
    AliFemtoModelWeightGeneratorBasic *tWeight = new AliFemtoModelWeightGeneratorBasic();
        tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
   
 //////// tWeight->SetPairType(5);
  //tWeight->SetPairType(1);

  AliFemtoModelManager *tModelManager = new AliFemtoModelManager();
  tModelManager->AcceptFreezeOutGenerator(tFreeze);
  tModelManager->AcceptWeightGenerator(tWeight);
  tModelManager->CreateCopyHiddenInfo(kFALSE);
  //tModelManager->CreateCopyHiddenInfo(kTRUE);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[30];
  AliFemtoBasicEventCut         *mecetaphitpc[30];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[30];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[30];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[30];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[30];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[30];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[30];
  AliFemtoESDTrackCut           *dtc1etaphitpc[30];
  AliFemtoESDTrackCut           *dtc2etaphitpc[30];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[30];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[30];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[30];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[30];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[30];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[30];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[30];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[30];
  //  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[30];
  AliFemtoShareQualityPairCut      *sqpcetaphitpc[30];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[30];
  AliFemtoModelCorrFctnDirectYlm     *cylmetaphitpc[30];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[30];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[30];
  AliFemtoKTPairCut             *ktpcuts[30*7];
  AliFemtoModelCorrFctnDirectYlm     *cylmkttpc[30*7];
  //AliFemtoModelQinvCorrFctn          *cqinvkttpc[30*7];
  AliFemtoModelCorrFctnSource          *cqinvkttpc[30*7];
//  AliFemtoModelCorrFctn          *cqinvkttpc[30*7];
  //AliFemtoModelCorrFctn3DLCMSSym     *cq3dlcmskttpc[30*7];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[30];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[30*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[30];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[30*10];


  // *** Begin pion-pion analysis ***
  int aniter = 0;

  for (int imult=0; imult<10; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<3; ichg++) {
        if (runch[ichg]) {
          aniter = ichg*10+imult;

          anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -8.0, 8.0, 4, multbins[imult], multbins[imult+1]);
          anetaphitpc[aniter]->SetNumEventsToMix(5);
          anetaphitpc[aniter]->SetMinSizePartCollection(1);

          mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
          mecetaphitpc[aniter]->SetEventMult(0,200000);
          mecetaphitpc[aniter]->SetVertZPos(-8,8);


          dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
          dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();


    
          if (ichg == 0) {
            dtc1etaphitpc[aniter]->SetCharge(1.0);
            dtc1etaphitpc[aniter]->SetPt(0.2,1.5);
          }
          else if (ichg == 1) {
            dtc1etaphitpc[aniter]->SetCharge(-1.0);
            dtc1etaphitpc[aniter]->SetPt(0.2,1.5);
          }
          else if (ichg == 2) {
            dtc1etaphitpc[aniter]->SetCharge(-1.0);
            dtc2etaphitpc[aniter]->SetCharge(1.0);
            dtc1etaphitpc[aniter]->SetPt(0.2,1.5);
            dtc2etaphitpc[aniter]->SetPt(0.2,1.5);
          }

          dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
          dtc1etaphitpc[aniter]->SetMass(KaonMass);
//test          dtc1etaphitpc[aniter]->SetMostProbableKaon();
//test          dtc1etaphitpc[aniter]->SetNsigma(3.0);
          //dtc1etaphitpc[aniter]->SetNsigma(2.0);
//test          dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
          //dtc1etaphitpc[aniter]->SetNsigmaTPConly(kTRUE);

          if (ichg == 2) {
            dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
            dtc2etaphitpc[aniter]->SetMass(KaonMass);
//test            dtc2etaphitpc[aniter]->SetMostProbableKaon();
//test            dtc2etaphitpc[aniter]->SetNsigma(3.0);
            //dtc2etaphitpc[aniter]->SetNsigma(2.0);
//test            dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
            //dtc2etaphitpc[aniter]->SetNsigmaTPConly(kTRUE);

          }


           dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
           dtc1etaphitpc[aniter]->SetminTPCncls(80);
           dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
           dtc1etaphitpc[aniter]->SetLabel(kFALSE);
           dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
           dtc1etaphitpc[aniter]->SetMaxImpactXY(2.4); // 2.4 0.1
           dtc1etaphitpc[aniter]->SetMaxImpactZ(3.2); // 2.0 0.1

           AliFemtoCutMonitorParticleYPt *cutPassYPtkp = new AliFemtoCutMonitorParticleYPt("cutPasskp", 0.493677);
           AliFemtoCutMonitorParticleYPt *cutFailYPtkp = new AliFemtoCutMonitorParticleYPt("cutFailkp", 0.493677);
           dtc1etaphitpc[aniter]->AddCutMonitor(cutPassYPtkp, cutFailYPtkp);

           AliFemtoCutMonitorParticlePtPDG *cutPassPidkp = new AliFemtoCutMonitorParticlePtPDG("cutPasskp", 0.493677);
           AliFemtoCutMonitorParticlePtPDG *cutFailPidkp = new AliFemtoCutMonitorParticlePtPDG("cutFailkp", 0.493677);
           dtc1etaphitpc[aniter]->AddCutMonitor(cutPassPidkp, cutFailPidkp);

           AliFemtoCutMonitorParticleMomRes *cutPassMRkp = new AliFemtoCutMonitorParticleMomRes("cutPasskp");
           AliFemtoCutMonitorParticleMomRes *cutFailMRkp = new AliFemtoCutMonitorParticleMomRes("cutFailkp");
           dtc1etaphitpc[aniter]->AddCutMonitor(cutPassMRkp, cutFailMRkp);


          sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
          //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
          sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
          sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
          sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);


          anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
          if (ichg == 2) {
            anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
            anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
          }
          else {
            anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
            anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
          }

          anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);


          //Qinv (without kT bins)
          cqinvkttpc[aniter] = new AliFemtoModelCorrFctnSource(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
          //cqinvkttpc[aniter] = new AliFemtoModelCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
        //  ktpcuts[0] = new AliFemtoKTPairCut(0.2,0.4);
         // cqinvkttpc[aniter]->SetPairSelectionCut(ktpcuts[0]);
          cqinvkttpc[aniter]->ConnectToManager(tModelManager);
          anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);


          if (runktdep) {
            int ktm;
            for (int ikt=0; ikt<3; ikt++) {
              ktm = aniter*3 + ikt;
              ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);

              cqinvkttpc[ktm] = new AliFemtoModelCorrFctnSource(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
              cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              cqinvkttpc[ktm]->ConnectToManager(tModelManager);
              anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

              //if (run3d) {
              //		cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,(imult>3)?((imult>6)?((imult>7)?0.6:0.4):0.25):0.15);
              //cq3dlcmskttpc[ktm] = new AliFemtoModelCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,0.5);
              //cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              //anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
              //}
            }
          }

          //cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
          //anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

          Manager->AddAnalysis(anetaphitpc[aniter]);
        }
      }
    }
  }
  // *** End pion-pion analysis

  return Manager;
}
