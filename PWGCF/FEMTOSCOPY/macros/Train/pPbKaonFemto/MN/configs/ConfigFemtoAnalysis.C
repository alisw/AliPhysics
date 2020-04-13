
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
//#include "AliFemtoKKTrackCut.h"
#include "AliFemtoKpm45TrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
//#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistanceKK.h"
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
#include "AliFemtoCutMonitorParticlePtPDG.h"
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
#include "AliFemtoKTPairCut.h"
#include "AliFemtoCutMonitorCollections.h"
#include "AliFemtoCorrFctnNonIdDR.h"



#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {
    
    double PionMass = 0.13956995;
    double KaonMass = 0.493677;
    
    //multiplicity bins
    int runmults[3] = {1, 1, 1};
    int multbins[4] = {0.01, 200, 400, 900};
    
    double PhiStarDifferenceMinimum=0.;//0.02; //[radian]
    double EtaDifferenceMinimum=0.;//0.02; //[radian]
    
    int runch[2] = {1, 1}; // Why?
    const char *chrgs[2] = { "Kp", "Km"};
    
    int runktdep = 1;
    double ktrng[3] = {0.2, 0.5, 1.0};
    
    int run3d = 0; // Do 3D cartesian analysis?
    //int runshlcms = 1;
    int runshlcms = 0;
    double shqmax;
    //int nbinssh = 200;
    int nbinssh = 100;
    
    //if (runshlcms) shqmax = 2.0;
    if (runshlcms) shqmax = 0.25;
    else shqmax = 2.0;
    
    AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
    Reader->SetUseMultiplicity(AliFemtoEventReaderAODChain::kCentrality);
    Reader->SetFilterBit(7);
    Reader->SetCentralityPreSelection(0.0001, 900);
    Reader->SetReadMC(kTRUE);
    Reader->SetKaonAnalysis(kTRUE);
    
    Reader->SetDCAglobalTrack(kTRUE);
    //Reader->SetpA2013(kTRUE);
    
    AliFemtoModelGausLCMSFreezeOutGenerator *tFreeze = new AliFemtoModelGausLCMSFreezeOutGenerator();
    tFreeze->SetSizeOut(1.7*TMath::Sqrt(1.5));
    tFreeze->SetSizeSide(1.7*TMath::Sqrt(1.5));
    tFreeze->SetSizeLong(1.7*TMath::Sqrt(1.5));
    
    AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
    tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
    tWeight->SetCoulOff();
    //tWeight->SetCoulOn();
    tWeight->SetQuantumOff();
    tWeight->SetStrongOff();
    //tWeight->SetStrongOn();
    tWeight->Set3BodyOff();
    
    //tWeight->SetKpKmModelType(14,1);//(Martin-f0, Achasov2-a0,0->phi off;1->phi on)
    //No phi:
    tWeight->SetKpKmModelType(14,0);//(Martin-f0, Achasov2-a0,0->phi off;1->phi on)
    
    AliFemtoModelManager *tModelManager = new AliFemtoModelManager();
    tModelManager->AcceptFreezeOutGenerator(tFreeze);
    tModelManager->AcceptWeightGenerator(tWeight);
    tModelManager->CreateCopyHiddenInfo(kTRUE);
    
    AliFemtoManager* Manager=new AliFemtoManager();
    Manager->SetEventReader(Reader);
    
    AliFemtoVertexMultAnalysis    *anetaphitpc[20];
    AliFemtoBasicEventCut         *mecetaphitpc[20];
    AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
    AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
    AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
    AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
    AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
    AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
    AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
    AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
    //AliFemtoKKTrackCut           *dtc1etaphitpc[20];
    //AliFemtoKKTrackCut           *dtc2etaphitpc[20];
    //AliFemtoESDTrackCut           *dtc1etaphitpc[20];
    //AliFemtoESDTrackCut           *dtc2etaphitpc[20];
    AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
    AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
    AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
    AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
    AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
    AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
    AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
    AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
    //AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
    AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[20];
    //AliFemtoShareQualityPairCut      *sqpcetaphitpc[20];
    //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];
    AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
    AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
    AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
    AliFemtoCorrFctnGammaMonitor  *cgamma[20*10];
    //AliFemtoKTPairCut             *ktpcuts[20*7];
    //AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*7];
    //AliFemtoQinvCorrFctn          *cqinvkttpc[20*7];
    //AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*7];
    AliFemtoKTPairCut             *ktpcuts[20*8];
    AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
    AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
    AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
    AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
    AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
    AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
    AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
    AliFemtoModelCorrFctn         *cqinvkttpcmodel[20*8];
    AliFemtoCorrFctnNonIdDR       *cfdourat[20*10];


    
    // *** Begin pion-pion analysis ***
    int aniter = 0;
    int ichg=0;
    for (int imult=0; imult<3; imult++) {
        if (runmults[imult]) {
            if (runch[ichg]) {
                aniter = ichg*3+imult;
                
                anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
                anetaphitpc[aniter]->SetNumEventsToMix(30);
                anetaphitpc[aniter]->SetMinSizePartCollection(1);
                anetaphitpc[aniter]->SetVerboseMode(kFALSE); //why?
                
                mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
                mecetaphitpc[aniter]->SetEventMult(0,10000);
                mecetaphitpc[aniter]->SetVertZPos(-10,10);
                
                cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
                
                cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
                
                cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
                cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
                mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
                
                
                //-----------------------1 particle-------------------------------------------<
                dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();
                dtc1etaphitpc[aniter]->SetCharge(1.0);
                dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
                dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
                //PID method
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
                //------------------- November 2013 ----------------------------------->
                
                dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                
                dtc1etaphitpc[aniter]->SetLabel(kFALSE);
                //----------------------2particle----------------------< KR
                // dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();
                dtc2etaphitpc[aniter]=new AliFemtoKpm45TrackCut();
                dtc2etaphitpc[aniter]->SetCharge(-1.0);
                dtc2etaphitpc[aniter]->SetPt(0.14,1.5);
                dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
                //PID method
                dtc2etaphitpc[aniter]->SetMass(KaonMass);
                dtc2etaphitpc[aniter]->SetMostProbableKaon();
                //dtc2etaphitpc[aniter]->SetPIDMethod(AliFemtoESDTrackCut::kContour);
                //------------------- November 2013 -----------------------------------<
                // new cuts to remove electron (do not take into analysis if 400<p<500)
                dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
                dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
                dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(2.0);
                dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
                dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
                // new cuts are stronger, better separation of pion in TOF
                // when momentum is greater then 800 MeV/c
                dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
                dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
                dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
                //------------------- November 2013 ----------------------------------->
                //Track quality cuts
                
                dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                
                dtc2etaphitpc[aniter]->SetLabel(kFALSE);
                
                
                cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
                cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
                dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
                
                cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg+1], imult), 0.493677); //ichg+1 --> ichg=1 --> charge = -1
                cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg+1], imult), 0.493677);
                
                dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);
                /*****************************************************/
                
                cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
                cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
                dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
                
                cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[ichg+1], imult),1);
                cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[ichg+1], imult),1);
                dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
                
                //sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
                
                
                sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKK();
                sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
                sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
                sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);
                
                
                
                /*****************************************************/
                anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
                anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]); //druga czastka
                anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
                /*****************************************************/
                //Qinv (without kT bins)
                cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
                anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);
                
                //3D cartesian (without kT bins)
                if(run3d){
                    cq3dlcmskttpc[aniter] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%i", chrgs[ichg], imult),100,0.5);
                    anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[aniter]);
                }
                
                /*****************************************************/
                if (runktdep) {
                    int ktm;
                    for (int ikt=0; ikt<2; ikt++) {
                        ktm = aniter*2 + ikt;
                        ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
                        /*****************************************************/
                        cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
                        cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                        anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
                        
                        cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
                        cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                        anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);
                        
                        // model--------------
                        
                        cqinvkttpcmodel[ktm] = new AliFemtoModelCorrFctn(Form("cqinvModelKK%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);
                        
                        cqinvkttpcmodel[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
                        //cqinvkttpcmodel[ktm]->SetKaonPDG(kTRUE);//Special MC analysis for K selected by PDG code -->
                        cqinvkttpcmodel[ktm]->SetKaonPDG(kFALSE);//w/o special MC analysis
                        
                        cqinvkttpcmodel[ktm]->ConnectToManager(tModelManager);
                        anetaphitpc[aniter]->AddCorrFctn(cqinvkttpcmodel[ktm]);// add CF histos
                        
                        
                        
                        
                        cgamma[aniter] = new AliFemtoCorrFctnGammaMonitor(Form("cgammaM%ikT%i", imult, ikt),200,200);
                        anetaphitpc[aniter]->AddCorrFctn(cgamma[aniter]);
                        
                        
                        /*****************************************************/
                        //              cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
                        //              cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                        //              cqinvinnertpc[ktm]->SetRadius(1.6);
                        //              anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);
                        //              cgamma[aniter] = new AliFemtoCorrFctnGammaMonitor(Form("cgammaM%ikT%i", imult, ikt),200,200);
                        //              anetaphitpc[aniter]->AddCorrFctn(cgamma[aniter]);
                        
//                        cfdourat[ktm] = new AliFemtoCorrFctnNonIdDR(Form("cfKstr%stpcM%ikT%i", chrgs[ichg], imult, ikt), nbinssh,0.0,shqmax);//AliFemtoCorrFctnNonIdDR(char* title, const int& nbins, const float& QinvLo, const float& QinvHi)
//                        cfdourat[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
//                        anetaphitpc[aniter]->AddCorrFctn(cfdourat[ktm]);// add CF histos
//
//                        if (run3d) {
//                            cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,0.5);
//                            cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
//                            anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
//                        }
                    }
                }
                
                // cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
                // anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
                
                Manager->AddAnalysis(anetaphitpc[aniter]);
            }
        }
    }
    // *** End pion-pion analysis
    
    return Manager;
}

