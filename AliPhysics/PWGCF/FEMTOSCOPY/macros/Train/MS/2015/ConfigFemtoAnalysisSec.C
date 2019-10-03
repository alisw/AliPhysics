

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
#include "AliFemtoPairCutRadialDistanceLM.h"
#include "AliFemtoPairCutRadialDistanceKK.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
//#include "AliFemtoAnalysisAzimuthalPbPb3rdOrder.h"
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
AliFemtoManager* ConfigFemtoAnalysis(int CentL=10, int CentH=20, int kTRange=4) {

    double PionMass = 0.13956995;
    double KaonMass = 0.493677;




    if (CentL==0 && CentH==5) {
        int runmults[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    }
    if (CentL==5 && CentH==10) {
        int runmults[10] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    }
    if (CentL==10 && CentH==20) {
        int runmults[10] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
    }
    if (CentL==20 && CentH==30) {
        int runmults[10] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
    }
    if (CentL==30 && CentH==40) {
        int runmults[10] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
    }
    if (CentL==40 && CentH==50 ) {
        int runmults[10] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    }



    int runktdeprange[4]={1,1,1,1};


    if(kTRange==1){
        int runktdeprange[4]={1,0,0,0};
    }
    if(kTRange==2){
        int runktdeprange[4]={1,1,0,0};
    }
    if(kTRange==3){
        int runktdeprange[4]={1,1,1,0};
    }
    if(kTRange==4){
        int runktdeprange[4]={1,1,1,1};
    }
    char *nametest="sysstudy";
    double PionMass = 0.13956995;
    double KaonMass = 0.493677;

    int multbins[11] = {0, 50, 100, 200, 300, 400, 500, 600, 800, 900, 990};

    int runch[2] = {1, 1};
    const char *chrgs[2] = { "pip", "pim" };

    //int rundetadphistar[4]={1,1,1,1};
    //const char *detadphistar[4]={"nocuts","original","test3","test4"};
    int rundetadphistar[5]={1,0,0,0,0};
    //const char *detadphistar[5]={"nocuts","test1","test2","test3","test4"};
    const char *detadphistar[5]={"test1a","test2a","test3a","nocuts","b"};


    int runktdep = 1;
    double ktrng[5] = {0.2, 0.3, 0.4, 0.5, 0.7};
    //int phirange[7] = {-15, 15, 45, 75, 105, 135, 165};
    int phirange[10] = {-15, 5, 25,45,65,85,105,125,145,165};
    int phirange2[11] = {-9, 9, 27, 45, 63, 81, 99, 117, 136, 153, 171};// Moe


    int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
    int isrealdata = 0;

    AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODMultSelection();
    Reader->SetFilterBit(7);
    Reader->SetEPVZERO(kFALSE); //keep this line it fills the whole histogram of event plane
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    Reader->SetReadV0(0);
    Reader->SetCentralityFlattening(kFALSE);
    AliFemtoManager* Manager=new AliFemtoManager();
    Manager->SetEventReader(Reader);

    AliFemtoAnalysisAzimuthalPbPb2Order		*ana[2][8][2];				//[charge][mult][cuts]
    AliFemtoBasicEventCut				*mecetaphitpc[2][8][2];	//[charge][mult][cuts]
    AliFemtoCutMonitorEventMult			*cutPassEvMetaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoCutMonitorEventMult			*cutFailEvMetaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoCutMonitorEventVertex		*cutPassEvVetaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoCutMonitorEventVertex		*cutFailEvVetaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoESDTrackCut					*dtc1etaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoESDTrackCut					*dtc2etaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[2][8][2];//[charge][mult][cuts]
    AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpcRD[2][8][2];//[charge][mult][cuts]


    //AliFemtoKTPairCut             *ktpaircut[2][10][5][4][9];		//[charge][mult][cuts][ktrange][phibins]
    //AliFemtoBPLCMS3DCorrFctn     *cq3dlcmskttpc[2][10][5][4][9];	//[charge][mult][cuts][ktrange][phibins]
    AliFemtoKTPairCut             *ktpaircut[2][8][2][4][10];		//[charge][mult][cuts][ktrange][phibins]
    AliFemtoBPLCMS3DCorrFctn     *cq3dlcmskttpc[2][8][2][4][10];	//[charge][mult][cuts][ktrange][phibins]

    for (int icuts=0; icuts<2; icuts++) {
        if (rundetadphistar[icuts]) {
            for (int imult=0; imult<8; imult++) {
                if (runmults[imult]) {
                    for (int ichg=0; ichg<2; ichg++) {
                        if (runch[ichg]) {

                            ana[ichg][imult][icuts] = new AliFemtoAnalysisAzimuthalPbPb2Order(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1],10);
                            //ana[ichg][imult][icuts] = new AliFemtoAnalysisAzimuthalPbPb2Order(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1],6);

                            ana[ichg][imult][icuts]->SetNumEventsToMix(3);
                            ana[ichg][imult][icuts]->SetMinSizePartCollection(4);
                            ana[ichg][imult][icuts]->SetEPhistname(Form("hist%i%i_%s",ichg,imult,detadphistar[icuts]));

                            mecetaphitpc[ichg][imult][icuts]= new AliFemtoBasicEventCut();
                            //mecetaphitpc[ichg][imult][icuts]->SetEventMult(10,100000);  //was 0.001
                            mecetaphitpc[ichg][imult][icuts]->SetEventMult(10,100000);  //was 0.001

                            mecetaphitpc[ichg][imult][icuts]->SetVertZPos(-8,8);

                            cutPassEvMetaphitpc[ichg][imult][icuts] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i_%s", chrgs[ichg], imult,detadphistar[icuts]));
                            cutFailEvMetaphitpc[ichg][imult][icuts] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i_%s", chrgs[ichg], imult,detadphistar[icuts]));
                            mecetaphitpc[ichg][imult][icuts]->AddCutMonitor(cutPassEvMetaphitpc[ichg][imult][icuts], cutFailEvMetaphitpc[ichg][imult][icuts]);

                            cutPassEvVetaphitpc[ichg][imult][icuts] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i_%s", chrgs[ichg], imult,detadphistar[icuts]));
                            cutFailEvVetaphitpc[ichg][imult][icuts] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i_%s", chrgs[ichg], imult,detadphistar[icuts]));
                            mecetaphitpc[ichg][imult][icuts]->AddCutMonitor(cutPassEvVetaphitpc[ichg][imult][icuts], cutFailEvVetaphitpc[ichg][imult][icuts]);

                            dtc1etaphitpc[ichg][imult][icuts] = new AliFemtoESDTrackCut();

                            if (ichg == 0)
                                dtc1etaphitpc[ichg][imult][icuts]->SetCharge(1.0);
                            else if (ichg == 1)
                                dtc1etaphitpc[ichg][imult][icuts]->SetCharge(-1.0);

                            dtc1etaphitpc[ichg][imult][icuts]->SetPt(0.15,1.5);
                            dtc1etaphitpc[ichg][imult][icuts]->SetEta(-1,1);
                            dtc1etaphitpc[ichg][imult][icuts]->SetRapidity(-0.8,0.8);
                            dtc1etaphitpc[ichg][imult][icuts]->SetMass(PionMass);

                            dtc1etaphitpc[ichg][imult][icuts]->SetMostProbablePion();

                            dtc1etaphitpc[ichg][imult][icuts]->SetStatus(AliESDtrack::kTPCin);
                            dtc1etaphitpc[ichg][imult][icuts]->SetminTPCncls(80);//was 80
                            dtc1etaphitpc[ichg][imult][icuts]->SetRemoveKinks(kTRUE);
                            dtc1etaphitpc[ichg][imult][icuts]->SetLabel(kFALSE);
                            dtc1etaphitpc[ichg][imult][icuts]->SetMaxTPCChiNdof(4.0);
                            dtc1etaphitpc[ichg][imult][icuts]->SetMaxImpactXY(2.4);//was 0.20
                            dtc1etaphitpc[ichg][imult][icuts]->SetMaxImpactZ(3.0);//was 0.15

                            sqpcetaphitpc[ichg][imult][icuts] = new AliFemtoPairCutRadialDistanceKK();
                            sqpcetaphitpcRD[ichg][imult][icuts] = new AliFemtoPairCutRadialDistanceKK();

                            sqpcetaphitpc[ichg][imult][icuts]->SetShareQualityMax(1.0);
                            sqpcetaphitpc[ichg][imult][icuts]->SetShareFractionMax(1.0); //was 1
                            sqpcetaphitpc[ichg][imult][icuts]->SetRemoveSameLabel(kFALSE);
                            sqpcetaphitpc[ichg][imult][icuts]->SetMinimumRadius(1.6);


                            sqpcetaphitpcRD[ichg][imult][icuts]->SetShareQualityMax(1.0);
                            sqpcetaphitpcRD[ichg][imult][icuts]->SetShareFractionMax(0.05);
                            sqpcetaphitpcRD[ichg][imult][icuts]->SetRemoveSameLabel(kFALSE);
                            sqpcetaphitpcRD[ichg][imult][icuts]->SetMinimumRadius(1.6);
                            if (icuts==0) // ( test1)
                            {
                                /*sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0); //0.012,0.018,0.024,0.048
                                 sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0);     //0.017,0.0225,0.034,0.068
                                 sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //was 0
                                 sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);//was 0
                                 */
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.017); //0.012,0.018,0.024,0.048
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.015);     //0.017,0.0225,0.034,0.068
                                sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //was 0
                                sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);//was 0


                            }
                            if (icuts==1) // (test2)
                            {
                                /*sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.015); //0.012,0.018,0.024,0.048
                                 sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.017);     //0.017,0.0225,0.034,0.068
                                 sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //0.012,0.018,0.024,0.048
                                 sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);     //0.017,0.0225,0.034,0.068
                                 */
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.04); //0.012,0.018,0.024,0.048
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.02);     //0.017,0.0225,0.034,0.068
                                sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //was 0
                                sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);//was 0

                            }
                            if (icuts==2) // (test3)
                            {
                                /*sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.02); //0.012,0.018,0.024,0.048
                                 sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.04);     //0.017,0.0225,0.034,0.068
                                 sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //0.012,0.018,0.024,0.048
                                 sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);     //0.017,0.0225,0.034,0.068
                                 */
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.06); //0.012,0.018,0.024,0.048
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.02);     //0.017,0.0225,0.034,0.068
                                sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //was 0
                                sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);//was 0

                            }
                            if (icuts==3) //test4 deta=0.01, dphistar=0.04
                            {
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //0.012,0.018,0.024,0.048
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);     //0.017,0.0225,0.034,0.068
                                sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //0.012,0.018,0.024,0.048
                                sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);     //0.017,0.0225,0.034,0.068

                            }
                            if (icuts==4) //test4 deta=0.01, dphistar=0.04
                            {
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.045); //0.012,0.018,0.024,0.048
                                sqpcetaphitpcRD[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.02);     //0.017,0.0225,0.034,0.068
                                sqpcetaphitpc[ichg][imult][icuts]->SetPhiStarDifferenceMinimum(0.0); //0.012,0.018,0.024,0.048
                                sqpcetaphitpc[ichg][imult][icuts]->SetEtaDifferenceMinimum(0.0);     //0.017,0.0225,0.034,0.068

                            }

                            ana[ichg][imult][icuts]->SetEventCut(mecetaphitpc[ichg][imult][icuts]);
                            ana[ichg][imult][icuts]->SetFirstParticleCut(dtc1etaphitpc[ichg][imult][icuts]);
                            ana[ichg][imult][icuts]->SetSecondParticleCut(dtc1etaphitpc[ichg][imult][icuts]);
                            ana[ichg][imult][icuts]->SetPairCut(sqpcetaphitpc[ichg][imult][icuts]);
                            ana[ichg][imult][icuts]->SetPairCutRD(sqpcetaphitpcRD[ichg][imult][icuts]);


                            for (int ikt=0; ikt<4; ikt++){
                                if (runktdeprange[ikt]) {
                                    for (int iphi=0; iphi<10; iphi++){


                                        ktpaircut[ichg][imult][icuts][ikt][iphi] = new AliFemtoKTPairCut(ktrng[ikt],ktrng[ikt+1]);
                                        ktpaircut[ichg][imult][icuts][ikt][iphi]->SetPhiRange(phirange2[iphi],phirange2[iphi+1]);

                                        cq3dlcmskttpc[ichg][imult][icuts][ikt][iphi] = new AliFemtoBPLCMS3DCorrFctn(Form("%s_cq3d%imult%ikT%iRP%i_%s", nametest,ichg, imult, ikt, iphi,detadphistar[icuts]),96,-0.12,0.12);//was 30 bins
                                        cq3dlcmskttpc[ichg][imult][icuts][ikt][iphi]->SetPairSelectionCut(ktpaircut[ichg][imult][icuts][ikt][iphi]);
                                        ana[ichg][imult][icuts]->AddCorrFctn(cq3dlcmskttpc[ichg][imult][icuts][ikt][iphi]);



                                    } //end of phi
                                }//end of kt case
                            }//end of kt
                            Manager->AddAnalysis(ana[ichg][imult][icuts]);
                        }//end runcharge
                    }//end charge
                }//end runmult
            }  //end mult
        } //end rundeltaetadeltaphistar
    } //end deltaetap
    return Manager;
}
