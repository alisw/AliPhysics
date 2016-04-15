

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
AliFemtoManager* ConfigFemtoAnalysis(int CentL=20, int CentH=30, int kTRange=4) {
    
    double PionMass = 0.13956995;
    double KaonMass = 0.493677;
    
    
    
    int runmults[6] = {1, 0, 0, 0, 0, 0};

    if (CentL==0 && CentH==5) {
        int runmults[6] = {1, 0, 0, 0, 0, 0};
    }
    if (CentL==5 && CentH==10) {
        int runmults[6] = {0, 1, 0, 0, 0, 0};
    }
    if (CentL==10 && CentH==20) {
        int runmults[6] = {0, 0, 1, 0, 0, 0};
    }
    if (CentL==20 && CentH==30) {
        int runmults[6] = {0, 0, 0, 1, 0, 0};
    }
    if (CentL==30 && CentH==40) {
        int runmults[6] = {0, 0, 0, 0, 1, 0};
    }
    if (CentL==40 && CentH==50 ) {
        int runmults[6] = {0, 0, 0, 0, 0, 1};
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
    
    int multbins[7] = {0, 50, 100, 200, 300, 400, 500};
    
    int runch[2] = {1, 1};
    const char *chrgs[2] = { "pip", "pim" };


    int runktdep = 1;
    double ktrng[5] = {0.2, 0.3, 0.4, 0.5, 0.7};
    //int phirange[7] = {-15, 15, 45, 75, 105, 135, 165};
    int phirange2[11] = {-9, 9, 27, 45, 63, 81, 99, 117, 136, 153, 171};// Moe
    
    int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
    int isrealdata = 0;
    
    AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
    Reader->SetFilterBit(7);
    Reader->SetEPVZERO(kFALSE); //keep this line it fills the whole histogram of event plane
    //Reader->SetCentralityPreSelection(0.001, 110);//was 0.001
    Reader->SetCentralityPreSelection(0.001, 510);//was 0.001
    
    
    AliFemtoManager* Manager=new AliFemtoManager();
    Manager->SetEventReader(Reader);
    
    /*AliFemtoAnalysisAzimuthalPbPb		*ana[2][10];				//[charge][mult]
     AliFemtoBasicEventCut				*mecetaphitpc[2][10];	//[charge][mult]
     AliFemtoCutMonitorEventMult			*cutPassEvMetaphitpc[2][10];//[charge][mult]
     AliFemtoCutMonitorEventMult			*cutFailEvMetaphitpc[2][10];//[charge][mult]
     AliFemtoCutMonitorEventVertex		*cutPassEvVetaphitpc[2][10];//[charge][mult]
     AliFemtoCutMonitorEventVertex		*cutFailEvVetaphitpc[2][10];//[charge][mult]
     AliFemtoESDTrackCut					*dtc1etaphitpc[2][10];//[charge][mult]
     AliFemtoESDTrackCut					*dtc2etaphitpc[2][10];//[charge][mult]
     AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpc[2][10];//[charge][mult]
     AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpcRD[2][10];//[charge][mult]
     AliFemtoKTPairCut             *ktpaircut[2][10][4][9];		//[charge][mult][ktrange][phibins]
     AliFemtoBPLCMS3DCorrFctn     *cq3dlcmskttpc[2][10][4][9];	//[charge][mult][ktrange][phibins]
     */
    
    //including delta eta delta phi star cuts
    AliFemtoAnalysisAzimuthalPbPb2Order		*ana[2][6];				//[charge][mult][cuts]
    AliFemtoBasicEventCut				*mecetaphitpc[2][6];	//[charge][mult][cuts]
    AliFemtoESDTrackCut					*dtc1etaphitpc[2][6];//[charge][mult][cuts]
    AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[2][6];//[charge][mult][cuts]
    AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpcRD[2][6];//[charge][mult][cuts]
    AliFemtoKTPairCut             *ktpaircut[2][6][4][10];		//[charge][mult][cuts][ktrange][phibins]
    AliFemtoBPLCMS3DCorrFctn     *cq3dlcmskttpc[2][6][4][10];	//[charge][mult][cuts][ktrange][phibins]
    
    for (int imult=0; imult<6; imult++) {
        if (runmults[imult]) {
            for (int ichg=0; ichg<2; ichg++) {
                if (runch[ichg]) {
                    
                    ana[ichg][imult] = new AliFemtoAnalysisAzimuthalPbPb2Order(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1],10);
                    
                    ana[ichg][imult]->SetNumEventsToMix(3);
                    ana[ichg][imult]->SetMinSizePartCollection(4);
                    ana[ichg][imult]->SetEPhistname(Form("hist%i%i",ichg,imult));
                    
                    mecetaphitpc[ichg][imult]= new AliFemtoBasicEventCut();
                    //mecetaphitpc[ichg][imult][icuts]->SetEventMult(10,100000);  //was 0.001
                    mecetaphitpc[ichg][imult]->SetEventMult(0.01,100000);  //was 0.001
                    mecetaphitpc[ichg][imult]->SetVertZPos(-8,8);
                    
                    dtc1etaphitpc[ichg][imult]= new AliFemtoESDTrackCut();
                    
                    if (ichg == 0)
                        dtc1etaphitpc[ichg][imult]->SetCharge(1.0);
                    else if (ichg == 1)
                        dtc1etaphitpc[ichg][imult]->SetCharge(-1.0);
                    
                    dtc1etaphitpc[ichg][imult]->SetPt(0.15,1.5);
                    dtc1etaphitpc[ichg][imult]->SetEta(-1,1);
                    dtc1etaphitpc[ichg][imult]->SetRapidity(-0.8,0.8);
                    dtc1etaphitpc[ichg][imult]->SetMass(PionMass);
                    
                    dtc1etaphitpc[ichg][imult]->SetMostProbablePion();
                    
                    dtc1etaphitpc[ichg][imult]->SetStatus(AliESDtrack::kTPCin);
                    dtc1etaphitpc[ichg][imult]->SetminTPCncls(80);//was 80
                    dtc1etaphitpc[ichg][imult]->SetRemoveKinks(kTRUE);
                    dtc1etaphitpc[ichg][imult]->SetLabel(kFALSE);
                    dtc1etaphitpc[ichg][imult]->SetMaxTPCChiNdof(4.0);
                    dtc1etaphitpc[ichg][imult]->SetMaxImpactXY(2.4);//was 0.20
                    dtc1etaphitpc[ichg][imult]->SetMaxImpactZ(3.0);//was 0.15
                    
                    sqpcetaphitpc[ichg][imult]= new AliFemtoPairCutRadialDistanceKK();
                    sqpcetaphitpcRD[ichg][imult] = new AliFemtoPairCutRadialDistanceKK();
                    
                    sqpcetaphitpc[ichg][imult]->SetShareQualityMax(1.0);
                    sqpcetaphitpc[ichg][imult]->SetShareFractionMax(1.0); //was 1
                    sqpcetaphitpc[ichg][imult]->SetRemoveSameLabel(kFALSE);
                    sqpcetaphitpc[ichg][imult]->SetMinimumRadius(1.6);
                    
                    
                    sqpcetaphitpcRD[ichg][imult]->SetShareQualityMax(1.0);
                    sqpcetaphitpcRD[ichg][imult]->SetShareFractionMax(0.05);
                    sqpcetaphitpcRD[ichg][imult]->SetRemoveSameLabel(kFALSE);
                    sqpcetaphitpcRD[ichg][imult]->SetMinimumRadius(1.6);
                    
                    sqpcetaphitpcRD[ichg][imult]->SetPhiStarDifferenceMinimum(0.017); //0.012,0.018,0.024,0.048
                    sqpcetaphitpcRD[ichg][imult]->SetEtaDifferenceMinimum(0.015);     //0.017,0.0225,0.034,0.068
                    sqpcetaphitpc[ichg][imult]->SetPhiStarDifferenceMinimum(0.0); //was 0
                    sqpcetaphitpc[ichg][imult]->SetEtaDifferenceMinimum(0.0);//was 0
                    
                    
                    
                    
                    ana[ichg][imult]->SetEventCut(mecetaphitpc[ichg][imult]);
                    ana[ichg][imult]->SetFirstParticleCut(dtc1etaphitpc[ichg][imult]);
                    ana[ichg][imult]->SetSecondParticleCut(dtc1etaphitpc[ichg][imult]);
                    ana[ichg][imult]->SetPairCut(sqpcetaphitpc[ichg][imult]);
                    ana[ichg][imult]->SetPairCutRD(sqpcetaphitpcRD[ichg][imult]);
                    
                    
                    for (int ikt=0; ikt<4; ikt++){
                        if (runktdeprange[ikt]) {
                            for (int iphi=0; iphi<10; iphi++){
                                
                                
                                ktpaircut[ichg][imult][ikt][iphi] = new AliFemtoKTPairCut(ktrng[ikt],ktrng[ikt+1]);
                                ktpaircut[ichg][imult][ikt][iphi]->SetPhiRange(phirange2[iphi],phirange2[iphi+1]);
                                
                                cq3dlcmskttpc[ichg][imult][ikt][iphi] = new AliFemtoBPLCMS3DCorrFctn(Form("%s_cq3d%imult%ikT%iRP%i", nametest,ichg, imult, ikt, iphi),56,-0.14,0.14);//was 30 bins
                                cq3dlcmskttpc[ichg][imult][ikt][iphi]->SetPairSelectionCut(ktpaircut[ichg][imult][ikt][iphi]);
                                ana[ichg][imult]->AddCorrFctn(cq3dlcmskttpc[ichg][imult][ikt][iphi]);
                                
                                
                                
                            } //end of phi
                        }//end of kt case
                    }//end of kt
                    Manager->AddAnalysis(ana[ichg][imult]);
                }//end runcharge
            }//end charge
        }//end runmult
    }  //end mult
    return Manager;
}

