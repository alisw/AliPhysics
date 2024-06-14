
/*********************************************************************
 *                                                                   *
 * ConfigFemtoAnalysis.C - configuration macro for the femtoscopic   *
 * analysis, meant as a QA process for two-particle effects          *
 *                                                                   *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)                         *
 *********************************************************************
 *    K+K- in PbPb@2.76TeV                                           *
 * Update: Konstantin.Mikhaylov@cern.ch                              *
 *         KpmHBT16  0010  03-Mar-2016                               *
 *         TOF PID from 0.45 GeV/c                                   *
 *         SetNsigmaTPC400_450(1.0) and other set to standard        *
 *         no cut on mom and no on gamma                             *
 *         no cut on phi* eta (set to 0)                             *
 *         10 MeV/q-bin in CF                                        *
 *********************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoManager.h"
#include "AliFemtoCorrFctnDEtaDPhi.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoKpm45TrackCut.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoPairCutRadialDistanceKK.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#endif

//________________________________________________________________________
AliFemtoManager *ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  const int cMu = 8;
  const int cKt = 1; // 4

  //-------Single track cuts------------------------------------------------->
  double DCAxy = 0.3; // cm // our standard is 0.20 cm; super narrow was 0.015cm
  double DCAz = 0.3;  // cm // our standard is 0.15 cm;
  
  double PhiStarDifferenceMinimum = 0.0; //[radian]
  double EtaDifferenceMinimum = 0.0;     //[radian]
  
  int runmults[cMu] = {1, 1, 1, 1, 1, 1, 1, 1};
  int multbins[cMu + 1] = {0, 50, 100, 200, 300, 400, 500, 700, 900};

  int runch[2] = {1, 1}; // K+-
  const char *chrgs[2] = {"Kp", "Km"};
  double ktrng[cKt + 1] = {0.2, 1.2};

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(128);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);

  Reader->SetReadMC(kTRUE);
  Reader->SetDCAglobalTrack(kTRUE);
  Reader->SetKaonAnalysis(kTRUE);

  AliFemtoManager *Manager = new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis *anetaphitpc[20];
  AliFemtoBasicEventCut *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoKpm45TrackCut *dtc1etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePtPDG *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];

  AliFemtoPairCutRadialDistanceKK *sqpcetaphitpc[20]; // Dhevan's dphi* cut
  AliFemtoCorrFctnDEtaDPhi *cdedpetaphi[20];
  AliFemtoKTPairCut *ktpcuts[20];

  int aniter = 0;

  bool verbose = false;

  for (int imult = 0; imult < cMu /*4*/; imult++) {
    if (runmults[imult]) {
      for (int ichg = 0; ichg < 2 /*K+-*/; ichg++) { // one loop
        if (runch[ichg]) {
          aniter = ichg * cMu + imult;

          anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(4, -8.0, 8.0, 5, multbins[imult], multbins[imult + 1]);
          anetaphitpc[aniter]->SetNumEventsToMix(20);
          anetaphitpc[aniter]->SetMinSizePartCollection(1);
          anetaphitpc[aniter]->SetVerboseMode(verbose);

          mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
          mecetaphitpc[aniter]->SetEventMult(0, 100000);
          mecetaphitpc[aniter]->SetVertZPos(-8.0, 8.0);

          cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);

          cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);

          dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();

          if (ichg == 0)
            dtc1etaphitpc[aniter]->SetCharge(1.0);
          else if (ichg == 1)
            dtc1etaphitpc[aniter]->SetCharge(-1.0);

          dtc1etaphitpc[aniter]->SetPt(0.14, 1.5);
          dtc1etaphitpc[aniter]->SetEta(-0.8, 0.8);
          dtc1etaphitpc[aniter]->SetMass(KaonMass);

          dtc1etaphitpc[aniter]->SetMostProbableKaon();

          dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
          dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
          dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);

          dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
          dtc1etaphitpc[aniter]->SetminTPCncls(80);
          dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
          dtc1etaphitpc[aniter]->SetLabel(kFALSE);
          dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
          dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
          dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);

          cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
          cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
          dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

          cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
          cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
          dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);

          sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKK(); // Dhevan's dphi* cut

          sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
          sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
          sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
          sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);

          anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
          anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]); // K+
          anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]); // K-
          anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);

          ktpcuts[aniter] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt + 1]);

          cq3dlcmskttpc[aniter]->SetPairSelectionCut(ktpcuts[aniter]);
          cq3dlcmskttpc[aniter]->SetKaonPDG(kTRUE);

          cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult), 1.2);
          anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

          Manager->AddAnalysis(anetaphitpc[aniter]);
        }
      }
    }
  }

  return Manager;
}
