/*enum EtaWindowSets
{
 en_eta_standart,
 en_eta_forward_keep_small_and_distant,
};*/

//#include "AliLRCProcess.h"
//#include "AliEventShape.h"
//#include "etaWindowsConfigEnumerator.h"

enum etaWindowsConfigurations
{
    en_etaWinConfig_standart = 0,
    en_etaWinConfig_empty,
    en_etaWinConfig_standart_extended,
    en_etaWinConfig_fixedFwd,
    en_etaWinConfig_fixedBkwd,
    en_etaWinConfig_phiSectors,
    //en_etaWinConfig_phiSectorsWithFixedFwdWindow,
    //en_etaWinConfig_phiSectorsWithFixedBkwdWindow,
    en_etaWinConfig_phiSectorsWithFixedBkwdWindow_eta0_4,
    en_etaWinConfig_twistedPhiWindows,
    en_etaWinConfig_longRangeEtaWindows,
    en_etaWinConfig_ridge_eta0_2_step_0_2_phi_8

};



AliESDtrackCuts* createAliLRCcuts(TString mode)//char* mode)
{

    if(mode=="Global2_TPC_clusters_20")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(20);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }

    if(mode=="Global2_TPC_clusters_40")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(40);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }


    if(mode=="Global2_TPC_clusters_100")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(100);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }
    /*
    if(mode=="Global2_TPC_clusters_110")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(110);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit();
        cuts->SetRequireITSRefit();
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }*/


    if(mode=="Global2_TPC_clusters_120")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(120);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }

    if(mode=="Global2_DCA_0.1")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.1);
        cuts->SetMaxDCAToVertexZ(0.1);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }

    /*
    if(mode=="Global2_DCA_0.3")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.3);
        cuts->SetMaxDCAToVertexZ(0.3);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }*/
    
    if(mode=="Global2_DCA_2.0")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(2.);
        cuts->SetMaxDCAToVertexZ(2.);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }
    
    if(mode=="Global2_min_ITS_clusters_0")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(0);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }

    if(mode=="Global2_min_ITS_clusters_4")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(4);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }

    if(mode=="Global2_no_refits")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kFALSE);
        cuts->SetRequireITSRefit(kFALSE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }


    if(mode=="Global2")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        //cuts->SetMinNCrossedRowsTPC(80);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }
    
    if(mode=="Global2_noPtCut")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(80);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.01,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }
    /*
    if(mode=="Global2_old_AI")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(70);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit();
        cuts->SetRequireITSRefit();
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(0.5);
        cuts->SetMaxDCAToVertexZ(0.5);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.8,1.8);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }*/
    
    
    if(mode=="Global2_softCuts")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(40);
        cuts->SetMinNClustersITS(0);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kFALSE);
        cuts->SetRequireITSRefit(kFALSE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(3.);
        cuts->SetMaxDCAToVertexZ(3.);
        cuts->SetPtRange(0.1,100);
        cuts->SetEtaRange(-2.,2.);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }
    
    
    if(mode=="StandardITSTPCTrackCuts2010")
    {
        AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);


        return cuts;
    }
    
    if(mode=="StandardITSTPCTrackCuts2010no")
    {
        AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);


        return cuts;
    }
    
    if(mode=="StandardITSTPCTrackCuts2010trueAndCrossRows")
    {
        AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, 1);


        return cuts;
    }
    
    if(mode=="StandardTPCOnlyTrackCuts")
    {
        AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();


        return cuts;
    }

    if(mode=="StandardITSSATrackCuts2010")
    {
        AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardITSSATrackCuts2010();


        return cuts;
    }




    /*if(mode=="Global2_eta_in_1")
    {
        AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
        cuts->SetMinNClustersTPC(70);
        cuts->SetMinNClustersITS(2);
        cuts->SetMaxChi2PerClusterTPC(4.0);
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetAcceptKinkDaughters(kFALSE);
        cuts->SetMaxDCAToVertexXY(2.);
        cuts->SetMaxDCAToVertexZ(2.);
        cuts->SetPtRange(0.15,100);
        cuts->SetEtaRange(-1.,1.);
        cuts->SaveHistograms("trackCuts");

        return cuts;
    }*/
    
}

//===========================================================================
AliAnalysisTaskLRC* createLRCtaskSkeleton(char* name="Task_LRC", Bool_t RunKine=kFALSE)
{
    AliAnalysisTaskLRC *taskLRC = new AliAnalysisTaskLRC(name,RunKine);
    taskLRC->SetMinPtLimit(0.2);
    taskLRC->SetMaxPtLimit(3.5);
    taskLRC->SetCheckForVtxPosition(kTRUE);
    taskLRC->SetVtxDiamond(0.4,0.4,5.0);
    return taskLRC;
}
//===========================================================================


void tuneEtaPhiWindows( AliAnalysisTaskLRC* taskLRC
                        , int nEtaWindows
                        , int nPhiWindows
                        , double etaWinWidth
                        , double etaWindowStep
                        , int ridgeStudy = 0
                        , double phiWindowWidthByHand = -1
        )
{
    //create processors, tune Eta
    for(int i = 0; i < nPhiWindows; i++)
    {
        //double etaWinWidth = 0.2;
        //double etaWindowStep = 0.2;//1;
        for ( int etaId = 0; etaId < nEtaWindows; etaId++ )
        {
            if ( ridgeStudy == 0 )
            {
                double winEtaBegin  = etaWindowStep * (etaId)-etaWinWidth/2;
                double winEtaEnd    = etaWindowStep * (etaId)-etaWinWidth/2 + etaWinWidth;
                taskLRC->AddLRCProcess( new AliLRCProcess(-winEtaEnd,-winEtaBegin,winEtaBegin,winEtaEnd));

//                const int includeSameEtaGapWindows = 1; //simply add same eta window processors
//                if ( includeSameEtaGapWindows )
//                {

//                    if ( etaId == 0 ) //put same eta gap windows
//                        taskLRC->AddLRCProcess( new AliLRCProcess(-etaWinWidth/2,etaWinWidth/2,-etaWinWidth/2,etaWinWidth/2));
//                    else // continue with symmetrical
//                    {
//                        double winEtaBegin  = etaWindowStep * (etaId - 1);
//                        double winEtaEnd    = etaWindowStep * (etaId - 1) + etaWinWidth;
//                        taskLRC->AddLRCProcess( new AliLRCProcess(-winEtaEnd,-winEtaBegin,winEtaBegin,winEtaEnd));
//                    }
//                }
//                else //usual simmetrical study without overlaping
//                {
//                    double winEtaBegin  = etaWindowStep * etaId;
//                    double winEtaEnd    = etaWindowStep * etaId + etaWinWidth;
//                    taskLRC->AddLRCProcess( new AliLRCProcess(-winEtaEnd,-winEtaBegin,winEtaBegin,winEtaEnd));
//                }
            }
            else if ( ridgeStudy == 1 )
            {
                double winEtaBegin0  = etaWindowStep * ( nEtaWindows/2 - 1 );
                double winEtaEnd0    = winEtaBegin0 + etaWinWidth;

                double winEtaBegin  = -winEtaBegin0 + etaWindowStep * etaId;
                double winEtaEnd    = -winEtaEnd0 + etaWindowStep * etaId;
                //double fwdFixedBegin  = etaWindowStep * ( nEtaWindows - 1 );
                //double fwdFixedEnd    = etaWindowStep * ( nEtaWindows - 1 ) + etaWinWidth;

                taskLRC->AddLRCProcess( new AliLRCProcess(winEtaEnd,winEtaBegin,winEtaBegin0,winEtaEnd0));
            }
            else if ( ridgeStudy == 51 ) //expanding windows study - win pair position more presice, close eta
            {
                //cout << "study win pair position" << endl;
                double shiftEta = etaId * 0.2;
                taskLRC->AddLRCProcess( new AliLRCProcess( -0.8 + shiftEta,-0.6+ shiftEta, -0.6 + shiftEta, -0.4 + shiftEta ));
            }
            else if ( ridgeStudy == 52 ) //expanding windows study - win pair position more presice, wide gap
            {
                //cout << "study win pair position" << endl;
                double shiftEta = etaId * 0.2;
                taskLRC->AddLRCProcess( new AliLRCProcess( -0.8 + shiftEta,-0.6+ shiftEta, -0.2 + shiftEta, -0.0 + shiftEta ));
            }
        }
    }
    //tune Phi
    //!!!taskLRC->SetNumberOfPhiSectors( nPhiWindows );
    double phiStep = 2 * TMath::Pi() / nPhiWindows;
    for ( Int_t sectorId = 0; sectorId < nPhiWindows; sectorId++ )
    {
        for ( Int_t i = nEtaWindows * sectorId; i < nEtaWindows * ( sectorId + 1 ); i++ )
        {
            double lFwdWinWidth = phiStep;
            double lBkwPhi_1 = phiStep * sectorId;
            double lBkwPhi_2 = phiStep * ( sectorId + 1 );

            if ( phiWindowWidthByHand > 0 ) // width by hand!
            {
                lFwdWinWidth = phiWindowWidthByHand;
                lBkwPhi_2 = phiStep * sectorId + phiWindowWidthByHand;
            }
            AliLRCBase* lrcBaseTmp = (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)));
            //            Double_t a;
            //            Double_t b;
            //            Double_t c;
            //            Double_t d;
            //            lrcBaseTmp->GetPhiWindows(a,b,c,d);
            lrcBaseTmp->SetForwardWindowPhi( 0, lFwdWinWidth );
            lrcBaseTmp->SetBackwardWindowPhi( lBkwPhi_1, lBkwPhi_2 );
        }
    }


//    taskLRC->SetNumberOfPhiSectors( nPhiWindows );
//    double phiStep = 2 * TMath::Pi() / nPhiWindows;
//    for ( Int_t sectorId = 0; sectorId < nPhiWindows; sectorId++ )
//    {
//        for ( Int_t i = nEtaWindows * sectorId; i < nEtaWindows * ( sectorId + 1 ); i++ )
//        {
//            double lBkwPhi_1 = phiStep * sectorId;
//            double lBkwPhi_2 = phiStep * ( sectorId + 1 );

//            (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetForwardWindowPhi( 0, phiStep );
//            (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetBackwardWindowPhi( lBkwPhi_1, lBkwPhi_2 );
//        }
//    }
    //end of eta-phi windows settings
}



//===========================================================================
void addAliLRCProcessors(AliAnalysisTaskLRC* taskLRC
                         , Int_t windowsConfigurationSetId = en_etaWinConfig_standart
        , Int_t nPhiSectors = 1
        , Double_t gapAsPartOfPi = 0. // for twisted sectors
        )
{

    int nLRCprocessors = 0;

    if( windowsConfigurationSetId == en_etaWinConfig_standart/*0*/ ) //standart eta window set
    {
        nLRCprocessors = 11;
        for(int i = 0; i < nPhiSectors; i++)
        {
            //FB
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,0.0,0.2));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,0.0,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,0.0,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,0.0,0.8));
            //0.2 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,0.2,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,0.4,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,0.6,0.8));
            //0.4 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,0.2,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,0.4,0.8));
            //0.6 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,0.2,0.8));

            //FULL
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8,-0.8,0.8));
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_empty ) // add nothing
    {

    }
    else if( windowsConfigurationSetId == en_etaWinConfig_standart_extended/*0*/ ) //standart eta window set
    {
        nLRCprocessors = 8;//14;
        for(int i = 0; i < nPhiSectors; i++)
        {
            //FB
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,0.0,0.2));
//            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,0.0,0.4));
//            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,0.0,0.6));
//            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,0.0,0.8));
            //0.2 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.3,-0.1,0.1,0.3));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,0.2,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.5,-0.3,0.3,0.5));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,0.4,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.7,-0.5,0.5,0.7));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,0.6,0.8));
            //0.4 gap
//            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,0.2,0.6));
//            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,0.4,0.8));
            //0.6 gap
//            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,0.2,0.8));

            //FULL
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8,-0.8,0.8));
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_fixedFwd/*1*/ ) //fixing fwd win (to test independance from bkwd width)
    {
        nLRCprocessors = 13;
        for(int i = 0; i < nPhiSectors; i++)
        {
            cout << "fixed fwd win" << endl;
            double windowMin = -0.8;
            double windowMax = -0.6;
            //FB
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.0,0.2));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.0,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.0,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.0,0.8));
            //0.2 gap
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.1,0.3));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.2,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.3,0.5));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.4,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.5,0.7));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.6,0.8));
            //0.4 gap
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.2,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.4,0.8));
            //0.6 gap
            taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,0.2,0.8));

            //FULL
            //taskLRC->AddLRCProcess(new AliLRCProcess( windowMin, windowMax,-0.8,0.8));
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_fixedBkwd/*2*/ )  //fixing BKW win, and for phi-voroching with fixed window
    {
        nLRCprocessors = 13;
        for(int i = 0; i < nPhiSectors; i++)
        {
            cout << "fixed fwd win" << endl;
            double windowMin = 0.6;
            double windowMax = 0.8;
            //FB
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0, windowMin, windowMax));
            //0.2 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.3,-0.1, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.5,-0.3, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.7,-0.5, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6, windowMin, windowMax));
            //0.4 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2, windowMin, windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4, windowMin, windowMax));
            //0.6 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2, windowMin, windowMax));

            //FULL
            //taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8, windowMin, windowMax));
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_phiSectors/*3*/ )  //phi-sectors configs //TMP?..
    {
        cout << "fixed FWD win" << endl;
        double windowMin = -0.8;
        double windowMax = -0.6;
        nLRCprocessors = 11;//8;
        for(int i = 0; i < nPhiSectors; i++)
        {
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,0.6,0.8));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,0.4,0.6));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,0.2,0.4));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,0.0,0.2));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,-0.2,-0.0));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,-0.4,-0.2));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,-0.6,-0.4));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(windowMin,windowMax,-0.8,-0.6));
            
            //taskLRC->AddLRCProcess( new AliLRCProcess(-0.8,-0.4,0.4,0.8) );
            //taskLRC->AddLRCProcess( new AliLRCProcess(-0.4, 0.0,0.0,0.4) );
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,0.0,0.2));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,0.0,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,0.0,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,0.0,0.8));
            //0.2 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,0.2,0.4));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,0.4,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,0.6,0.8));
            //0.4 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,0.2,0.6));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,0.4,0.8));
            //0.6 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,0.2,0.8));
            //FULL
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8,-0.8,0.8));
        }
    }
    /*    else if( windowsConfigurationSetId == en_etaWinConfig_phiSectorsWithFixedBkwdWindow ) //4
    {
        nLRCprocessors = 10;
        for(int i = 0; i < nPhiSectors; i++)
        {
            cout << "fixed fwd win" << endl;
            double windowMin = 0.6;
            double windowMax = 0.8;
            //FB
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,windowMin,windowMax));
            //0.2 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,windowMin,windowMax));
            //0.4 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,windowMin,windowMax));
            //0.6 gap
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,windowMin,windowMax));
        }
    }*/
    else if( windowsConfigurationSetId == en_etaWinConfig_phiSectorsWithFixedBkwdWindow_eta0_4/*5*/ )  //fixing BKW win, and for phi-rotating with fixed window (0.4 phi-width)
    {
        cout << "fixed bkwd win" << endl;
        double windowMin = 0.4;
        double windowMax = 0.8;
        nLRCprocessors = 4;
        for(int i = 0; i < nPhiSectors; i++)
        {
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(0.0,0.4,windowMin,windowMax));
            taskLRC->AddLRCProcess(new AliLRCProcess(0.4,0.8,windowMin,windowMax));
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_twistedPhiWindows/*6*/ )  // 3.09.12 - twisted windows in phi
    {
        cout << "studying TWISTED phi windows..." << endl;
        nLRCprocessors = 18;//11;
        const int nPhiTwistedSectors = 3; //const number of twisted windows
        for(int i = 0; i < nPhiTwistedSectors; i++)
        {
            double etaWinWidth = 0.5; //this tuned for true MC, when we have wide eta
            double windowStep = 0.25;
            for ( int etaId = 0; etaId < nLRCprocessors; etaId++ )
            {
                double winEtaBegin  = windowStep * etaId; //etaWinWidth * etaId;
                double winEtaEnd    = windowStep * etaId + etaWinWidth; //etaWinWidth * ( etaId + windowOverlap);
                taskLRC->AddLRCProcess(new AliLRCProcess(-winEtaEnd,-winEtaBegin,winEtaBegin,winEtaEnd));
            }
            //            //FB
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,0.0,0.2));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,0.0,0.4));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,0.0,0.6));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,0.0,0.8));
            //            //0.2 gap
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,0.2,0.4));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,0.4,0.6));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,0.6,0.8));
            //            //0.4 gap
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,0.2,0.6));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,0.4,0.8));
            //            //0.6 gap
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,0.2,0.8));

            //            //FULL
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8,-0.8,0.8));
        }
        //tune phi for processors:
        taskLRC->SetNumberOfPhiSectors(6); // make rotation for 6 times!

        //double phiStep = 2 * TMath::Pi() / nPhiSectors;
        double gapInPhiBetweenWins = gapAsPartOfPi;
        for ( Int_t sectorId = 0; sectorId < nPhiTwistedSectors; sectorId++ )
        {
            for ( Int_t i = nLRCprocessors * sectorId; i < nLRCprocessors * ( sectorId + 1 ); i++ )
            {
                double lBkwPhi_1 = 0;   //phiStep * sectorId;
                double lBkwPhi_2 = 0;   //phiStep * ( sectorId + 1 );
                if ( sectorId == 0 ) // || phi-windows
                {
                    lBkwPhi_1 = 0 + gapInPhiBetweenWins;
                    lBkwPhi_2 = 2./3. * TMath::Pi() - gapInPhiBetweenWins;
                }
                else if ( sectorId == 1 ) // transverse phi-windows
                {
                    lBkwPhi_1 = 2./3. * TMath::Pi() + gapInPhiBetweenWins;
                    lBkwPhi_2 = TMath::Pi() - gapInPhiBetweenWins;
                }
                else if ( sectorId == 2 ) // opposite phi-windows
                {
                    lBkwPhi_1 = TMath::Pi() + gapInPhiBetweenWins;
                    lBkwPhi_2 = 5./3. * TMath::Pi() - gapInPhiBetweenWins;
                }

                (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetForwardWindowPhi( 0 + gapInPhiBetweenWins, 2./3. * TMath::Pi() - gapInPhiBetweenWins );
                (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetBackwardWindowPhi( lBkwPhi_1, lBkwPhi_2 );
                if ( sectorId == 1 ) // transverse phi-windows
                    (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetDoubleSidedBackwardWindowPhi( true );
            }
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_longRangeEtaWindows ) //for MC truth study
    {
        nLRCprocessors = 8;//7;//14;
        for(int i = 0; i < nPhiSectors; i++)
        {
            double etaWinWidth = 0.2;
            double windowStep = 0.1;//1;
            for ( int etaId = 0; etaId < nLRCprocessors; etaId++ )
            {
                double winEtaBegin  = windowStep * etaId;
                double winEtaEnd    = windowStep * etaId + etaWinWidth;
                //                double winEtaBegin = etaWinWidth * etaId;
                //                double winEtaEnd = etaWinWidth * ( etaId + 1);
                taskLRC->AddLRCProcess(new AliLRCProcess(-winEtaEnd,-winEtaBegin,winEtaBegin,winEtaEnd));
            }
        }
    }
    else if( windowsConfigurationSetId == en_etaWinConfig_ridge_eta0_2_step_0_2_phi_8 )
    {
        cout << "fixed bkwd win" << endl;
        double windowMin = 0.4;
        double windowMax = 0.8;
        double etaWinWidth = windowMax - windowMin;
        double windowStep = 0.4;//1;
        nLRCprocessors = 4;
        for(int i = 0; i < nPhiSectors; i++)
        {
            for ( int etaId = nLRCprocessors-1; etaId >=0 ; etaId-- )
            {
                double winEtaBegin  = windowStep * etaId;
                double winEtaEnd    = windowStep * etaId + etaWinWidth;
                //                double winEtaBegin = etaWinWidth * etaId;
                //                double winEtaEnd = etaWinWidth * ( etaId + 1);
                taskLRC->AddLRCProcess(new AliLRCProcess(windowMax-winEtaEnd,windowMax-winEtaBegin
                                                         ,windowMin,windowMax));
            }
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,windowMin,windowMax));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,windowMin,windowMax));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(0.0,0.4,windowMin,windowMax));
            //            taskLRC->AddLRCProcess(new AliLRCProcess(0.4,0.8,windowMin,windowMax));
        }
    }
    else if( windowsConfigurationSetId == 100 )  //event shape
    {
        cout << "event shape" << endl;
        //double windowMin = 0.4;
        //double windowMax = 0.8;
        nLRCprocessors = 1;
        //for(int i = 0; i < nPhiSectors; i++)
        //{
        taskLRC->AddLRCProcess(new AliEventShape(-0.8,0.8,-0.8,0.8));
        //taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,windowMin,windowMax));
        //taskLRC->AddLRCProcess(new AliLRCProcess(0.0,0.4,windowMin,windowMax));
        //taskLRC->AddLRCProcess(new AliLRCProcess(0.4,0.8,windowMin,windowMax));
        //}
    }

    if ( nPhiSectors > 1 ) //tune phi sectors
    {
        taskLRC->SetNumberOfPhiSectors(nPhiSectors);
        //int nLRCprocessors = 8;
        double phiStep = 2 * TMath::Pi() / nPhiSectors;
        for ( Int_t sectorId = 0; sectorId < nPhiSectors; sectorId++ )
        {
            for ( Int_t i = nLRCprocessors * sectorId; i < nLRCprocessors * ( sectorId + 1 ); i++ )
            {
                double lBkwPhi_1 = phiStep * sectorId;
                double lBkwPhi_2 = phiStep * ( sectorId + 1 );

                (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetForwardWindowPhi( 0, phiStep );
                (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetBackwardWindowPhi( lBkwPhi_1, lBkwPhi_2 );
            }
        }

    } // endif AddPhiWindows IA variant


    
}

void setHistPtRange(AliAnalysisTaskLRC* taskLRC,  Double_t LoPt,Double_t HiPt, Double_t PtStep = 0.005, Int_t xAxisRebinFactor = 1 )
{
    Int_t nPtBins = (Int_t) ((HiPt - LoPt) / PtStep ); // divide by "resolution" of Pt histogram to obtain number of bins
    for(Int_t i=0; i < taskLRC->GetListOfProcessors()->GetEntries(); i++)
    {
        (dynamic_cast<AliLRCBase*> (taskLRC->Proc(i)))->SetHistPtRange( LoPt, HiPt, nPtBins );
        (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetHistPtRangeForwardWindowRebinFactor( xAxisRebinFactor );
    }
}

void setHistMultRange(AliAnalysisTaskLRC* taskLRC,  Int_t whichWindow, Int_t LoMult,Int_t HiMult
                      ,Bool_t isStandartWindowSet = kTRUE, Int_t MultBins=0) // Sets range for Nch histos axis
{
    //whichWindow - 0=for both, 1=for fwd, 2=bwd
    if (isStandartWindowSet ) //usual case
        for(Int_t i=0; i < taskLRC->GetListOfProcessors()->GetEntries(); i++)
        {
            (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetHistMultRange( whichWindow, LoMult, HiMult );
        }
    else  //try to tune bins according to eta-width
    {
        double eta[4];
        AliLRCProcess* process;
        for(Int_t i=0; i < taskLRC->GetListOfProcessors()->GetEntries(); i++)
        {
            process = (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)));
            process->GetETAWindows(eta[0],eta[1],eta[2],eta[3]);
            //look at window width and correct NN mult bins!
            //standart is 70-270 for dEtaF+dEtaB == 0.4
            double dEtaF = fabs( eta[1] - eta[0] );
            double dEtaB = fabs( eta[3] - eta[2] );
            double multCoeff = ( dEtaF + dEtaB ) / 0.4;
            process->SetHistMultRange( whichWindow, LoMult*multCoeff, HiMult*multCoeff  );
            //Int_t lowM = LoMult;
            //Int_t hiM = HiMult;
            /*Int_t shift = ( HiMult  ) / 2.;

        if ( i == 1 || i == 7 || i == 8 )
            (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetHistMultRange( whichWindow, 2*shift - 1.5*shift, 2*shift + 1.5*shift );
        else if ( i == 2|| i == 9 )
            (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetHistMultRange( whichWindow, 3*shift - 1.*shift, 3*shift + 2.*shift );
        else if ( i == 3 )
            (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetHistMultRange( whichWindow, 4*shift - 0.5*shift, 4*shift + 2.5*shift );
        else 		//if ( i == 0 || i == 4 || i == 5 || i == 6 )
            (dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetHistMultRange( whichWindow, 0 , 2*shift );//LoMult, HiMult );
                */
        }
    }

}





//===========================================================================
void configureLRCtaskOutput(AliAnalysisTaskLRC* taskLRC/*,TString OutputRootFolder=":PWG2LRC"*/
                            , TString strPrefixName= "", TString strRunMode = "default"/*by IA*/)
{
    if(!taskLRC)
        return;
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        Error("configureLRCtaskOutput", "No analysis manager to connect to.");
        return ;
    }

    if (!mgr->GetInputEventHandler())
    {
        Error("AddTaskLRC", "This task requires an input event handler");
        return ;
    }

    //TString type = mgr->GetInputEventHandler()->GetDataType();

    TString outputFileName= mgr->GetCommonFileName();
    if( outputFileName == "" )
        outputFileName = "LRC." + strRunMode + ".root";
    TString strOutputRootFolder = ":PWGCFLRC_" + strPrefixName;
    outputFileName += strOutputRootFolder ;

    TString listOutName;
    listOutName = taskLRC->GetName();

    AliAnalysisDataContainer *cout_LRC = mgr->CreateContainer(listOutName, TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    mgr->ConnectInput(taskLRC, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskLRC, 1, cout_LRC);
//    if ( 1 )//taskLRC->GetIncludeEventTreeInOutput() )
//    {
//        AliAnalysisDataContainer *outLRCtree = mgr->CreateContainer("eventTreeContainer", TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
//        mgr->ConnectOutput(taskLRC, 2, outLRCtree);
//    }

    cout << "N of LRC Processors ="<< taskLRC->GetListOfProcessors()->GetEntries() <<"\n";
    for(Int_t i=0; i < taskLRC->GetListOfProcessors()->GetEntries(); i++)
    {
        mgr->ConnectOutput(taskLRC,taskLRC->Proc(i)->GetOutputSlotNumber(),mgr->CreateContainer(strPrefixName+"_"+((taskLRC->Proc(i)->GetShortDef()+"_")+=i),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName));
    }

}


