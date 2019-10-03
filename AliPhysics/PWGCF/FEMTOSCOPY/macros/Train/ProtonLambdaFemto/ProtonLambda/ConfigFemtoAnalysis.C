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
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoCutMonitorCollections.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoPairOriginMonitor.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false,
                                     Bool_t fAvgSep=kTRUE, Bool_t fCosPAngle=kTRUE,
                                     Bool_t fDCAdaughters=kTRUE, Bool_t fDCAXYprotons=kTRUE,
                                     Bool_t fDCAZprotons=kTRUE, Bool_t fMaxDCAdaughters=kTRUE,
                                     Bool_t fMaxDCAV0=kTRUE, Bool_t fDecayLength=kTRUE,
                                     Bool_t fEtaDaughters=kTRUE, Bool_t fEtaProtons=kTRUE,
                                     Bool_t fEtaV0=kTRUE, Bool_t fpTdaughters=kTRUE,
                                     Bool_t fpTprotons=kTRUE, Bool_t fpTV0=kTRUE,
                                     Bool_t fInvariantMassK0s=kFALSE, Bool_t fSharedDaughter=kTRUE) {
    
    double PionMass = 0.13956995;
    double KaonMass = 0.493677;
    double ProtonMass = 0.938272013;
    double LambdaMass = 1.115683;
    
    double psi = TMath::Pi()/2.;
    double psid = TMath::Pi()/6.;
    
    int runepvzero[7] = {0, 0, 0, 0, 0, 0, 1};
    double epvzerobins[7] = {-psi, -psi+psid, -psi+2*psid, -psi+3*psid, -psi+4*psid, -psi+5*psid, -psi+6*psid};
    
    // Switches for QA analyses
    int runmults[10] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
    int multbins[11] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};
    int runch[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    const char *chrgs[10] = { "V0LL", "V0ALAL", "V0LAL", "V0PL", "V0APL", "V0PAL", "V0APAL","PP","PAP","APAP" };
    
    double ktrng[3] = {0.01, 0.7, 100};
    
    int numOfMultBins = 10;
    int numOfChTypes = 10;
    int numOfkTBins = 2;
    int numOfEPvzero = 7;
    
    bool performSharedDaughterCut = fSharedDaughter;
    bool enablePairMonitors = false;
    
    double shqmax = 1.0;
    int nbinssh = 100;
    
    double minK0Cut = 0.48;
    double maxK0Cut = 0.515;
    
    AliFemtoEventReaderAODMultSelection* Reader = new AliFemtoEventReaderAODMultSelection();
    Reader->SetFilterBit(7);
    Reader->SetReadV0(1); //Read V0
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    Reader->SetEPVZERO(kTRUE);
    Reader->SetCentralityFlattening(kTRUE);
    if(mcAnalysis) Reader->SetReadMC(kTRUE);

    AliFemtoManager* Manager=new AliFemtoManager();
    Manager->SetEventReader(Reader);
    
    AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
    AliFemtoModelWeightGeneratorBasic *tWeight = new AliFemtoModelWeightGeneratorBasic();
//    modelMgr->AcceptWeightGenerator(tWeight);
    
    AliFemtoVertexMultAnalysis    *anetaphitpc[320];
    AliFemtoBasicEventCut         *mecetaphitpc[320];
    AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[320];
    AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[320];
    AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[320];
    AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[320];
    AliFemtoCutMonitorCollections   *cutPassColletaphitpc[320];
    AliFemtoCutMonitorCollections   *cutFailColletaphitpc[320];
    AliFemtoV0TrackCut           *dtc1etaphitpc[320];
    AliFemtoV0TrackCut           *dtc2etaphitpc[320];
    AliFemtoESDTrackCut           *dtc3etaphitpc[320];
    AliFemtoESDTrackCut           *dtc4etaphitpc[320];
    AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[320];
    AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[320];
    AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[320];
    AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[320];
    AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[320];
    AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[320];
    AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[320];
    AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[320];
    AliFemtoCutMonitorV0          *cutPass1V0[320];
    AliFemtoCutMonitorV0          *cutFail1V0[320];
    AliFemtoCutMonitorV0          *cutPass2V0[320];
    AliFemtoCutMonitorV0          *cutFail2V0[320];
    AliFemtoV0PairCut             *sqp1cetaphitpc[320];
    AliFemtoV0TrackPairCut        *sqp2cetaphitpc[320];
    AliFemtoV0TrackPairCut        *sqp3cetaphitpc[320];
    AliFemtoV0TrackPairCut        *sqp4cetaphitpc[320];
    AliFemtoPairCutRadialDistance *sqp5cetaphitpc[320];
    AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[320];
    AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[320*6];
    AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[320];
    AliFemtoKTPairCut             *ktpcuts[320*6];
    AliFemtoCorrFctnDirectYlm     *cylmkttpc[320*6];
    AliFemtoQinvCorrFctn          *cqinvkttpc[320*6];
    AliFemtoAvgSepCorrFctn        *avgsepcorr[320*6];
    AliFemtoCorrFctnNonIdDR       *cnonidtpc[320*6];
    AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[320];
    AliFemtoShareQualityCorrFctn  *cqinvsqtpc[320];
    AliFemtoTPCInnerCorrFctn      *cqinvtitpc[320];
    AliFemtoQinvCorrFctn          *cqinvtpc[100*3];
    
    AliFemtoModelCorrFctn         *cQinvModel[320];
    
    // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
    // *** Begin pion-pion (positive) analysis ***
    int aniter = 0;
    
    for (int imult=0; imult<numOfMultBins; imult++) {
        if (runmults[imult]) {
            for (int ichg=0; ichg<numOfChTypes; ichg++) {
                if (runch[ichg]) {
                    for (int iepvzero = 0; iepvzero < numOfEPvzero; iepvzero++) {
                        if (runepvzero[iepvzero]) {
                            aniter = imult * numOfChTypes + ichg * numOfEPvzero + iepvzero;
                            
                            anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 4, multbins[imult], multbins[imult+1]);
                            anetaphitpc[aniter]->SetNumEventsToMix(10);
                            anetaphitpc[aniter]->SetMinSizePartCollection(1);
                            anetaphitpc[aniter]->SetVerboseMode(kFALSE);
                            
                            mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
                            mecetaphitpc[aniter]->SetEventMult(0,100000);
                            mecetaphitpc[aniter]->SetVertZPos(-8,8);
                            
                            if (iepvzero == (numOfEPvzero-1))
                                mecetaphitpc[aniter]->SetEPVZERO(epvzerobins[0],epvzerobins[numOfEPvzero-1]);
                            else
                                mecetaphitpc[aniter]->SetEPVZERO(epvzerobins[iepvzero],epvzerobins[iepvzero+1]);
                            
                            // cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%iPsi%d", chrgs[ichg], imult, iepvzero));
                            // cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%iPsi%d", chrgs[ichg], imult, iepvzero));
                            // mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
                            
                            // cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%iPsi%d", chrgs[ichg], imult, iepvzero));
                            // cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%iPsi%d", chrgs[ichg], imult, iepvzero));
                            // mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
                            
                            //V0 first particle cut -> Lambda ichg 0, 2, 6
                            dtc1etaphitpc[aniter] = new AliFemtoV0TrackCut();
                            dtc1etaphitpc[aniter]->SetMass(LambdaMass);
                            if(fEtaV0)
                                dtc1etaphitpc[aniter]->SetEta(0.8);
                            if(fpTV0)
                                dtc1etaphitpc[aniter]->SetPt(0.5, 5.0);
                            if(fEtaDaughters)
                                dtc1etaphitpc[aniter]->SetEtaDaughters(0.8);
                            if(fpTdaughters) {
                                dtc1etaphitpc[aniter]->SetPtPosDaughter(0.5, 4.0);
                                dtc1etaphitpc[aniter]->SetPtNegDaughter(0.16, 4.0);
                            }
                            dtc1etaphitpc[aniter]->SetTPCnclsDaughters(80);
                            dtc1etaphitpc[aniter]->SetNdofDaughters(4.0);
                            dtc1etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit);
                            dtc1etaphitpc[aniter]->SetOnFlyStatus(kFALSE);
                            dtc1etaphitpc[aniter]->SetParticleType(0);
                            if(fDCAdaughters)
                                dtc1etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.1, 0.3);
                            if(fMaxDCAdaughters)
                                dtc1etaphitpc[aniter]->SetMaxDcaV0Daughters(0.4);
                            if(fMaxDCAV0)
                                dtc1etaphitpc[aniter]->SetMaxDcaV0(0.5);
                            if(fDecayLength)
                                dtc1etaphitpc[aniter]->SetMaxV0DecayLength(60.0);
                            if(fCosPAngle)
                                dtc1etaphitpc[aniter]->SetMaxCosPointingAngle(0.9993);
                            if(fInvariantMassK0s)
                                dtc1etaphitpc[aniter]->SetInvariantMassRejectK0s(minK0Cut, maxK0Cut);
                            dtc1etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038, LambdaMass+0.0043);
                            
                            //V0 second particle cut -> AntiLambda ichg 1, 3, 4, 5
                            dtc2etaphitpc[aniter] = new AliFemtoV0TrackCut();
                            dtc2etaphitpc[aniter]->SetMass(LambdaMass);
                            if(fEtaV0)
                                dtc2etaphitpc[aniter]->SetEta(0.8);
                            if(fpTV0)
                                dtc2etaphitpc[aniter]->SetPt(0.5, 5.0);
                            if(fEtaDaughters)
                                dtc2etaphitpc[aniter]->SetEtaDaughters(0.8);
                            if(fpTdaughters) {
                                dtc2etaphitpc[aniter]->SetPtPosDaughter(0.16, 4.0);
                                dtc2etaphitpc[aniter]->SetPtNegDaughter(0.3, 4.0);
                            }
                            dtc2etaphitpc[aniter]->SetTPCnclsDaughters(80);
                            dtc2etaphitpc[aniter]->SetNdofDaughters(4.0);
                            dtc2etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit);
                            dtc2etaphitpc[aniter]->SetOnFlyStatus(kFALSE);
                            dtc2etaphitpc[aniter]->SetParticleType(1);
                            if(fMaxDCAdaughters)
                                dtc2etaphitpc[aniter]->SetMaxDcaV0Daughters(0.4);
                            if(fMaxDCAV0)
                                dtc2etaphitpc[aniter]->SetMaxDcaV0(0.5);
                            if(fDCAdaughters)
                                dtc2etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.3, 0.1);
                            if(fCosPAngle)
                                dtc2etaphitpc[aniter]->SetMaxCosPointingAngle(0.9993);
                            if(fDecayLength)
                                dtc2etaphitpc[aniter]->SetMaxV0DecayLength(60.0);
                            if(fInvariantMassK0s)
                                dtc2etaphitpc[aniter]->SetInvariantMassRejectK0s(minK0Cut, maxK0Cut);
                            dtc2etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0036, LambdaMass+0.0041);
                            
                            //ESD first particle cut -> Proton 3, 5; AntiProton 4, 6, 7, 8
                            dtc3etaphitpc[aniter] = new AliFemtoESDTrackCut();
                            dtc3etaphitpc[aniter]->SetMostProbableProton();
                            dtc3etaphitpc[aniter]->SetMass(ProtonMass);
                            dtc3etaphitpc[aniter]->SetCharge(1.0);
                            if(fpTprotons)
                                dtc3etaphitpc[aniter]->SetPt(0.7, 4.0);
                            if(fEtaProtons)
                                dtc3etaphitpc[aniter]->SetEta(-0.8, 0.8);
                            
                            // Track quality cuts
                            dtc3etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
                            dtc3etaphitpc[aniter]->SetminTPCncls(80);
                            dtc3etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                            dtc3etaphitpc[aniter]->SetLabel(kFALSE);
                            dtc3etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                            if(fDCAXYprotons)
                                dtc3etaphitpc[aniter]->SetMaxImpactXY(2.8);
                            if(fDCAZprotons)
                                dtc3etaphitpc[aniter]->SetMaxImpactZ(3.2);
                            dtc3etaphitpc[aniter]->SetNsigma(3.0);
                            dtc3etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
                            
                            //ESD first particle cut -> Proton 3, 5; AntiProton 4, 6, 8, 9
                            dtc4etaphitpc[aniter] = new AliFemtoESDTrackCut();
                            dtc4etaphitpc[aniter]->SetMostProbableProton();
                            dtc4etaphitpc[aniter]->SetMass(ProtonMass);
                            dtc4etaphitpc[aniter]->SetCharge(-1.0);
                            if(fpTprotons)
                                dtc4etaphitpc[aniter]->SetPt(0.7, 5.0);
                            if(fEtaProtons)
                                dtc4etaphitpc[aniter]->SetEta(-0.8, 0.8);
                            
                            // Track quality cuts
                            dtc4etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
                            dtc4etaphitpc[aniter]->SetminTPCncls(80);
                            dtc4etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                            dtc4etaphitpc[aniter]->SetLabel(kFALSE);
                            dtc4etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                            if(fDCAXYprotons)
                                dtc4etaphitpc[aniter]->SetMaxImpactXY(2.8);
                            if(fDCAZprotons)
                                dtc4etaphitpc[aniter]->SetMaxImpactZ(3.2);
                            dtc4etaphitpc[aniter]->SetNsigma(3.0);
                            dtc4etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
                            
                            // // //V0 monitor
                            cutPass1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass1%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero));
                            cutFail1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail1%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero));
                            dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1V0[aniter], cutFail1V0[aniter]);
                            
                            cutPass2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass2%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero));
                            cutFail2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail2%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero));
                            dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2V0[aniter], cutFail2V0[aniter]);
                            
                            // cutPass1V0Origin[aniter] = new AliFemtoCutMonitorV0OriginDependent(Form("cutPass1%stpcM%i", chrgs[ichg], imult));
                            // cutFail1V0Origin[aniter] = new AliFemtoCutMonitorV0OriginDependent(Form("cutFail1%stpcM%i", chrgs[ichg], imult));
                            // dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1V0Origin[aniter], cutFail1V0Origin[aniter]);
                            
                            // cutPass2V0Origin[aniter] = new AliFemtoCutMonitorV0OriginDependent(Form("cutPass2%stpcM%i", chrgs[ichg], imult));
                            // cutFail2V0Origin[aniter] = new AliFemtoCutMonitorV0OriginDependent(Form("cutFail2%stpcM%i", chrgs[ichg], imult));
                            // dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2V0Origin[aniter], cutFail2V0Origin[aniter]);
                            
                            // cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("trackCutPass1%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero),2);
                            // cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("trackCutFail1%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero),2);
                            // dtc3etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
                            
                            // cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("trackCutPass2%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero),2);
                            // cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("trackCutFail2%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero),2);
                            // // //if(ichg == 4 || ichg == 5 || ichg == 6 || ichg == 7 || ichg == 8 || ichg == 9)
                            // dtc4etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
                            
                            
                            sqp1cetaphitpc[aniter] = new AliFemtoV0PairCut();
                            sqp1cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                            sqp1cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                            sqp1cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                            if(fAvgSep) {
                                sqp1cetaphitpc[aniter]->SetMinAvgSeparation(0, 5 + 0.5*fAvgSep); //proton-pion+
                                sqp1cetaphitpc[aniter]->SetMinAvgSeparation(1, 5 + 0.5*fAvgSep); //proton-antiproton
                                sqp1cetaphitpc[aniter]->SetMinAvgSeparation(2, 0); //pion- - pion+
                                sqp1cetaphitpc[aniter]->SetMinAvgSeparation(3, 5 + 0.5*fAvgSep); //antiproton - pion-
                            }
                            sqp2cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-proton
                            sqp2cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
                            sqp2cetaphitpc[aniter]->SetShareFractionMax(0.05);
                            sqp2cetaphitpc[aniter]->SetTPCOnly(kTRUE);
                            sqp2cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                            sqp2cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                            sqp2cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                            sqp2cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton
                            if(fAvgSep) {
                                sqp2cetaphitpc[aniter]->SetMinAvgSeparation(0, 5); //0 - track-pos, 1 - track-neg
                                sqp2cetaphitpc[aniter]->SetMinAvgSeparation(1, 5);
                            }
                            sqp3cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //antilambda-antiproton
                            sqp3cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
                            sqp3cetaphitpc[aniter]->SetShareFractionMax(0.05);
                            sqp3cetaphitpc[aniter]->SetTPCOnly(kTRUE);
                            sqp3cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                            sqp3cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                            sqp3cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                            sqp3cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
                            if(fAvgSep) {
                                sqp3cetaphitpc[aniter]->SetMinAvgSeparation(0, 5); //0 - track-pos, 1 - track-neg
                                sqp3cetaphitpc[aniter]->SetMinAvgSeparation(1, 5);
                            }
                            
                            sqp4cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-antiproton, antilambda-proton
                            sqp4cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
                            sqp4cetaphitpc[aniter]->SetShareFractionMax(0.05);
                            sqp4cetaphitpc[aniter]->SetTPCOnly(kTRUE);
                            sqp4cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                            sqp4cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                            sqp4cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
                            //SetMinAvgSeparation w if'ach ponizej
                            
                            sqp5cetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
                            sqp5cetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.012);
                            sqp5cetaphitpc[aniter]->SetEtaDifferenceMinimum(0.017);
                            sqp5cetaphitpc[aniter]->SetShareQualityMax(1.0);
                            sqp5cetaphitpc[aniter]->SetShareFractionMax(0.05);
                            sqp5cetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
                            sqp5cetaphitpc[aniter]->SetMaxEEMinv(0.002);
                            sqp5cetaphitpc[aniter]->SetMaxThetaDiff(0.008);
                            sqp5cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
                            sqp5cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
                            sqp5cetaphitpc[aniter]->SetAvgsepMinimum(5.0);
                            
                            avgsepcorr[aniter] = new AliFemtoAvgSepCorrFctn(Form("Avgsep%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero),5000,0,500);
                            
                            // pairOriginPass[aniter] = new AliFemtoPairOriginMonitor(Form("Pass%stpcM%i", chrgs[ichg], imult));
                            // pairOriginFail[aniter] = new AliFemtoPairOriginMonitor(Form("Fail%stpcM%i", chrgs[ichg], imult));
                            
                            // sqp1cetaphitpc[aniter]->AddCutMonitor(pairOriginPass[aniter], pairOriginFail[aniter]); //lambda systems
                            // sqp2cetaphitpc[aniter]->AddCutMonitor(pairOriginPass[aniter], pairOriginFail[aniter]); //lambda-proton systems
                            // sqp3cetaphitpc[aniter]->AddCutMonitor(pairOriginPass[aniter], pairOriginFail[aniter]); //(anti-)lambda-(anti-)proton mixed systems
                            // sqp4cetaphitpc[aniter]->AddCutMonitor(pairOriginPass[aniter], pairOriginFail[aniter]); //antilambda-antiprotons
                            
                            anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                            anetaphitpc[aniter]->SetEnablePairMonitors(enablePairMonitors);
                            
                            if(ichg == 0) //V0LL
                            {
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
                            }
                            else if(ichg == 1) //V0ALAL
                            {
                                anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
                            }
                            else if(ichg == 2) //VOLAL
                            {
                                anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
                            }
                            else if(ichg == 3) //V0PL
                            {
                                anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp2cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
                            }
                            else if(ichg == 4) //V0APL
                            {
                                if(fAvgSep) {
                                    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(0, 5);  // antiproton - proton
                                    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(1, 8); // antiproton - pion-
                                }
                                anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
                            }
                            else if(ichg == 5) //V0PAL
                            {
                                if(fAvgSep) {
                                    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(0, 8); // proton - pion+
                                    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(1, 5); // proton - antiproton-
                                }
                                anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
                            }
                            else if(ichg == 6) //V0APAL
                            {
                                anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp3cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
                            }
                            else if(ichg ==7) //PP
                            {
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp5cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
                            }
                            else if(ichg ==8) //PAP
                            {
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp5cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
                            }
                            else if(ichg ==9) //APAP
                            {
                                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                                anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
                                anetaphitpc[aniter]->SetPairCut(sqp5cetaphitpc[aniter]);
                                avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
                            }
                            
                            if(ichg==3 || ichg==4 || ichg==5 || ichg==6 || ichg==8) { //PL, APL, PAL, APAL, PAP
                                cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero), nbinssh, 0.0,shqmax);
                                anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
                            }
                            else {
                                cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%iPsi%i", chrgs[ichg], imult, iepvzero),nbinssh,0.0,shqmax);
                                anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
                            }
                            
                            anetaphitpc[aniter]->AddCorrFctn(avgsepcorr[aniter]);
                            
                            if(mcAnalysis)
                            {
                                cQinvModel[aniter] = new AliFemtoModelCorrFctn(Form("cQinv_Model_%s_M%i", chrgs[ichg],imult), 400, 0, 2);
                                cQinvModel[aniter]->ConnectToManager(modelMgr);
                                anetaphitpc[aniter]->AddCorrFctn(cQinvModel[aniter]);
                            }
                            Manager->AddAnalysis(anetaphitpc[aniter]);
                        }
                    }
                }
            }
        }
    }
    // *** End of analysis
    return Manager;
}
