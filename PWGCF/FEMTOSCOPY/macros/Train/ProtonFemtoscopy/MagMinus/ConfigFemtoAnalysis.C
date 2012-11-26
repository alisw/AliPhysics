/*********************************************************************
 *                                                                   *
 * Configfemtoanalysis.C - configuration macro for the femtoscopic   *
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
#include "AliFemtoCorrFctnNonIdDR.h"
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
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

        double PionMass = 0.13956995;
        double KaonMass = 0.493677;
        double ProtonMass = 0.938272013;

        int runmults[10] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
        int multbins[11] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};

        int runch[3] = {1, 1, 1};
        const char *chrgs[3] = { "PP", "APAP", "PAP" };

        int runktdep = 1;
        double ktrng[5] = {0.1, 0.3, 0.6, 1.0, 10.0};

        int numOfMultBins = 10;
        int numOfChTypes = 3;
        int numOfkTbins = 4;

        int runqinv = 1;
        int runshlcms = 0; // 0:PRF(PAP), 1:LCMS(PP,APAP)

        int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
        int isrealdata = 1;

        //  int gammacut = 1;

        double shqmax = 2.0;
        int nbinssh = 100;

        AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
        Reader->SetFilterBit(7);
        Reader->SetCentralityPreSelection(0.001, 510);

        AliFemtoManager* Manager = new AliFemtoManager();
        Manager->SetEventReader(Reader);

        AliFemtoVertexMultAnalysis    *anetaphitpc[10*3];
        AliFemtoBasicEventCut         *mecetaphitpc[10*3];

        AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[10*3];
        AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[10*3];

        AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[10*3];
        AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[10*3];

        AliFemtoESDTrackCut           *dtc1etaphitpc[10*3];
        AliFemtoESDTrackCut           *dtc2etaphitpc[10*3];

        AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[10*3];
        AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[10*3];

        AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[10*3];
        AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[10*3];

        AliFemtoPairCutRadialDistance      *sqpcetaphitpc[10*3];

        AliFemtoKTPairCut             *ktpcuts[100*3*4];
        AliFemtoCorrFctnDirectYlm     *cylmkttpc[100*3];

        AliFemtoQinvCorrFctn          *cqinvkttpc[100*3*4];
        AliFemtoQinvCorrFctn          *cqinvtpc[100*3];

        AliFemtoCorrFctnNonIdDR       *ckstartpc[100*3];
        AliFemtoCorrFctnNonIdDR       *ckstarkttpc[100*3*4];

        // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
        // *** Begin pion-pion (positive) analysis ***
        int aniter = 0;

        for (int imult = 0; imult < numOfMultBins; imult++) {
                if (runmults[imult]) {

                        for (int ichg = 0; ichg < numOfChTypes; ichg++) {
                                if (runch[ichg]) {

                                        aniter = ichg * numOfMultBins + imult;

                                        anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 4, multbins[imult], multbins[imult+1]);
                                        anetaphitpc[aniter]->SetNumEventsToMix(10);
                                        anetaphitpc[aniter]->SetMinSizePartCollection(1);

                                        mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
                                        mecetaphitpc[aniter]->SetEventMult(0.001,100000);
                                        mecetaphitpc[aniter]->SetVertZPos(-8,8);

                                        if (isrealdata)
                                                mecetaphitpc[aniter]->SetAcceptOnlyPhysics(kTRUE);

                                        cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                                        cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                                        mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);

                                        // cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                                        // cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                                        // mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);

                                        dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
                                        dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();

                                        if (ichg == 0) { // pp
                                                dtc1etaphitpc[aniter]->SetCharge(1.0);
                                                dtc1etaphitpc[aniter]->SetPt(0.5,5.0);
                                                dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
                                                dtc1etaphitpc[aniter]->SetMass(ProtonMass);
                                                dtc1etaphitpc[aniter]->SetMostProbableProton();

                                        }

                                        else if (ichg == 1) { // apap
                                                dtc1etaphitpc[aniter]->SetCharge(-1.0);
                                                dtc2etaphitpc[aniter]->SetPt(0.2,5.0);
                                                dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
                                                dtc1etaphitpc[aniter]->SetMass(ProtonMass);
                                                dtc1etaphitpc[aniter]->SetMostProbableProton();


                                        }
                                        else if (ichg == 2) { // pap
                                                dtc1etaphitpc[aniter]->SetCharge(-1.0);
                                                dtc1etaphitpc[aniter]->SetPt(0.2,5.0);
                                                dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
                                                dtc1etaphitpc[aniter]->SetMass(ProtonMass);
                                                dtc1etaphitpc[aniter]->SetMostProbableProton();

                                                dtc2etaphitpc[aniter]->SetCharge(1.0);
                                                dtc2etaphitpc[aniter]->SetPt(0.5,5.0);
                                                dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
                                                dtc2etaphitpc[aniter]->SetMass(ProtonMass);
                                                dtc2etaphitpc[aniter]->SetMostProbableProton();
                                        }

                                        // Track quality cuts

                                        if (runtype == 0) {
                                                dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                                //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
                                                //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
                                                dtc1etaphitpc[aniter]->SetminTPCncls(80);
                                                dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                                dtc1etaphitpc[aniter]->SetLabel(kFALSE);
                                                //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                                dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                                dtc1etaphitpc[aniter]->SetMaxImpactXY(0.2);
                                                //            dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
                                                dtc1etaphitpc[aniter]->SetMaxImpactZ(0.15);
                                                //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                        }
                                        else if (runtype == 1) {
                                                //      dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                                //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
                                                //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit|AliESDtrack::kITSpureSA);
                                                //      dtc1etaphitpc[aniter]->SetminTPCncls(70);
                                                dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
                                                dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                                dtc1etaphitpc[aniter]->SetLabel(kFALSE);
                                                //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                                //      dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(6.0);
                                                dtc1etaphitpc[aniter]->SetMaxImpactXY(0.2);
                                                dtc1etaphitpc[aniter]->SetMaxImpactZ(0.25);
                                                //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                        }

                                        else if (runtype == 2) {

                                                dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
                                                dtc1etaphitpc[aniter]->SetminTPCncls(70);
                                                dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                                dtc1etaphitpc[aniter]->SetLabel(kFALSE);
                                                dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                                // dtc1etaphitpc[aniter]->SetMaxImpactXY(0.1); // 2.4 0.1
                                                dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
                                                dtc1etaphitpc[aniter]->SetMaxImpactZ(2.0); // 2.0 0.1

                                                if (ichg == 2) {

                                                        dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
                                                        dtc2etaphitpc[aniter]->SetminTPCncls(70);
                                                        dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                                        dtc2etaphitpc[aniter]->SetLabel(kFALSE);
                                                        dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                                        // dtc2etaphitpc[aniter]->SetMaxImpactXY(0.1); // 2.4 0.1
                                                        dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
                                                        dtc2etaphitpc[aniter]->SetMaxImpactZ(2.0); // 2.0 0.1
                                                }

                                        }

                                        cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
                                        cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),2);
                                        dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);



                                        cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),ProtonMass);
                                        cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),ProtonMass);
                                        dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);



                                        // sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
                                        sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();

                                        if (runtype == 0) {
                                                sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
                                                sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
                                                sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
                                                // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
                                                // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
                                                //	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
                                                //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
                                                //	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
                                        }
                                        else if (runtype == 1) {
                                                sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
                                                sqpcetaphitpc[aniter]->SetShareFractionMax(1.05);
                                                sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
                                                // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
                                                // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
                                                //	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
                                                //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
                                                //	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
                                        }
                                        else if (runtype == 2) {

                                                sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
                                                sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
                                                sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);

                                                //	    if (gammacut == 0) {
                                                //sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
                                                //sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
                                                //}
                                                //else if (gammacut == 1) {
                                                //sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
                                                //sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
                                                //}

                                                sqpcetaphitpc[aniter]->SetMagneticFieldSign(-1); // field1 -1, field3 +1
                                                sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.012); // 0.012 - pions, 0.017 - kaons, 0.018
                                                sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.015); // 0.017 - pions, 0.015 - kaons

                                        }

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


                                        if (ichg == 2) {
                                                ckstartpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("ckstar%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
                                                anetaphitpc[aniter]->AddCorrFctn(ckstartpc[aniter]);
                                        }
                                        else {
                                                cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
                                                anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
                                        }

                                        cylmkttpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),2,nbinssh, 0.0,shqmax,runshlcms);
                                        anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[aniter]);


                                        if (runktdep) {
                                                int ktm;
                                                for (int ikt=0; ikt<numOfkTbins; ikt++) {

                                                        ktm = aniter * numOfkTbins + ikt;
                                                        ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);

                                                        cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),2,nbinssh,0.0,shqmax,0);
                                                        cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                                                        anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);

                                                        if (ichg == 2) {
                                                                ckstarkttpc[ktm] = new AliFemtoCorrFctnNonIdDR(Form("ckstar%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
                                                                ckstarkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                                                                anetaphitpc[aniter]->AddCorrFctn(ckstarkttpc[ktm]);


                                                        }
                                                        else {
                                                                cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
                                                                cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                                                                anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
                                                        }

                                                }
                                        }


                                        Manager->AddAnalysis(anetaphitpc[aniter]);
                                }
                        }
                }
        }
        // *** End  analysis

        return Manager;
}
