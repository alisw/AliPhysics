#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>


#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMath.h"

#include "TFile.h"
#include "TString.h"

#include "SupplementaryClasses.h"

#include <iostream>
#include <fstream>

#include <algorithm>


using namespace std;

const int REDUCTION_FACTOR = 50;//1;//50;//10; //to reduce memory consumption

BootstrapHistos::BootstrapHistos():
    hist_avF(0)
  , hist_avB(0)
  , hist_avFB(0)
  , hist_avF2(0)
  , hist_DFB(0)
  , hist_DF2(0)
{}

void BootstrapHistos::InitHistos( const int cW, const char *strType, float _cMin, float _cMax, int _etaW, int _phiW, int _ptW, double histRangeFB )
{
    //    cout << "InitHistos for BS - " << strType << endl;
    initWithNameTitle( &hist_avF  , strType, "avF"    , cW, _cMin, _cMax, _etaW, _phiW, _ptW, histRangeFB );
    initWithNameTitle( &hist_avB  , strType, "avB"    , cW, _cMin, _cMax, _etaW, _phiW, _ptW, histRangeFB );
    initWithNameTitle( &hist_avFB , strType, "avFB"   , cW, _cMin, _cMax, _etaW, _phiW, _ptW, histRangeFB*histRangeFB );
    initWithNameTitle( &hist_avF2 , strType, "avF2"   , cW, _cMin, _cMax, _etaW, _phiW, _ptW, histRangeFB*histRangeFB );
    initWithNameTitle( &hist_DFB  , strType, "DFB"    , cW, _cMin, _cMax, _etaW, _phiW, _ptW, histRangeFB/10 ); //*10 );
    initWithNameTitle( &hist_DF2  , strType, "DF2"    , cW, _cMin, _cMax, _etaW, _phiW, _ptW, histRangeFB/10 ); //*10 );
    initWithNameTitle( &hist_bCorr, strType, "bCorr"  , cW, _cMin, _cMax, _etaW, _phiW, _ptW, 4 );
    initWithNameTitle( &hist_C2   , strType, "C2"     , cW, _cMin, _cMax, _etaW, _phiW, _ptW, 4 );

    sum_bCorr = 0;
    sum_bCorr2 = 0;
    n_entries = 0;

    av_bCorr_from_outside = 0;
    sum_delta_bCorr2 = 0;
}

void BootstrapHistos::initWithNameTitle( TH1D** hist, const char *strType, const char *strName, const int cW, float _cMin, float _cMax, int _etaW, int _phiW, int _ptW, double histRange )
{
    TString strBSname  = Form("hist1D_%s_%s_cW%d_c%.1f-%.1f_etaW_%d_phiW_%d_ptW_%d", strType, strName, cW, _cMin, _cMax, _etaW, _phiW, _ptW);
    TString strBStitle = Form("hist1D_%s_%s_cW%d_c%.1f-%.1f_etaW_%d_phiW_%d_ptW_%d;bCorr;entries", strType, strName, cW, _cMin, _cMax, _etaW, _phiW, _ptW);
    *hist = new TH1D( strBSname, strBStitle, 8000/REDUCTION_FACTOR, -histRange, histRange );
}

void BootstrapHistos::FillHistos( const CorrCoeffInfo &corrInfo )
{
    //    cout << "FillHistos for BS" << endl;
    hist_avF  ->Fill( corrInfo.avF   );
    hist_avB  ->Fill( corrInfo.avB   );
    hist_avFB ->Fill( corrInfo.avFB  );
    hist_avF2 ->Fill( corrInfo.avF2  );
    hist_DFB  ->Fill( corrInfo.DFB   );
    hist_DF2  ->Fill( corrInfo.DF2   );
    hist_bCorr->Fill( corrInfo.bCorr );
    hist_C2   ->Fill( corrInfo.C2    );

    sum_bCorr += corrInfo.bCorr;
    sum_bCorr2 += corrInfo.bCorr * corrInfo.bCorr;
    n_entries++;

    double tmpDelta = (corrInfo.bCorr - av_bCorr_from_outside);
    //cout << tmpDelta << endl;
    sum_delta_bCorr2 += tmpDelta*tmpDelta;

}

void BootstrapHistos::WriteHistos()
{
    TString strDirName = Form( "dir_%s" , hist_bCorr->GetName() );
    gFile->mkdir( strDirName.Data() );
    gFile->cd( strDirName.Data() );

    hist_bCorr->Write();
    hist_C2   ->Write();

    hist_DFB -> Write();
    hist_DF2 -> Write();
    hist_avFB-> Write();
    hist_avF2-> Write();
    hist_avF -> Write();
    hist_avB -> Write();
    gFile->cd();
}

void BootstrapHistos::ResetHistos()
{
    hist_bCorr->Reset();
    hist_C2   ->Reset();

    hist_DFB -> Reset();
    hist_DF2 -> Reset();
    hist_avFB-> Reset();
    hist_avF2-> Reset();
    hist_avF -> Reset();
    hist_avB -> Reset();

    sum_bCorr = 0;
    sum_bCorr2 = 0;
    n_entries = 0;

    av_bCorr_from_outside = 0;
    sum_delta_bCorr2 = 0;
}

// ##### WinPair
WinPair::WinPair() :
    cBinMin(0)
  , cBinMax(100)
  , etaW(-1)
  , phiW(-1)
  , ptW(-1)

  , initHistos(true)

  , NN_nF(0)
  , NN_nB(0)
  , NN_nF_nB(0)
  , NN_nF2(0)
  , NN_nB2(0)
  , NN_Nevents(0)

  , PtPt_PtF(0)
  , PtPt_PtB(0)
  , PtPt_PtF_PtB(0)
  , PtPt_PtF2(0)
  , PtPt_PtB2(0)
  , PtPt_Nevents(0)

  , PtN_nF(0)
  , PtN_PtB(0)
  , PtN_nF_PtB(0)
  , PtN_nF2(0)
  , PtN_PtB2(0)
  , PtN_Nevents(0)

  , sumPtN_nF(0)
  , sumPtN_sumPtB(0)
  , sumPtN_nF_sumPtB(0)
  , sumPtN_nF2(0)
  , sumPtN_sumPtB2(0)
  , sumPtN_Nevents(0)


  , nRuns(-1)

  //  , NN_bCorr(-1000)
  //  , NN_C2(-1000)
  //  , PtPt_bCorr(-1000)
  //  , PtPt_C2(-1000)
  //  , PtN_bCorr(-1000)
  //  , PtN_C2(-1000)

  , hist2D_NN                (0)
  , hist2D_PtPt              (0)
  , hist2D_PtN               (0)
  , hist2D_sumPtN            (0)
  , hist1D_EstimatorEntries  (0)
  , hist1D_multDistrF        (0)
  , hist1D_multDistrB        (0)
  , hist1D_QA_PtF            (0)
  , hist1D_QA_PtB            (0)

  , doBootstrap(false)
  , SAMPLING_nF(0)
  , SAMPLING_nB(0)
  , SAMPLING_PtF(0)
  , SAMPLING_PtB(0)

  , hist1D_multDistr_RunByRun_F(0x0)
  , hist1D_multDistr_RunByRun_B(0x0)
  , hist1D_avPtDistr_RunByRun_F(0x0)
  , hist1D_avPtDistr_RunByRun_B(0x0)

{}
void WinPair::init(int cW, float _cMin, float _cMax, int _etaW, int _phiW, int _ptW, bool _initHistos, int MAX_N_EVENTS_FOR_BOOTSTRAP )
{
    cBinMin = _cMin;
    cBinMax = _cMax;
    etaW = _etaW;
    phiW = _phiW;
    ptW = _ptW;

    initHistos = _initHistos;

    if ( !initHistos )
        return;

    TString strWinDescr = Form("c%.1f-%.1f_etaW_%d_phiW_%d", cBinMin, cBinMax, etaW, phiW);
    // !!! modify if >1 pt bins!
    if ( _ptW > 1 )
        strWinDescr = Form("%s_ptW_%d", strWinDescr.Data(), ptW);

    TString strNN = Form("hist2D_NN_%s", strWinDescr.Data());
    hist2D_NN = new TH2D( strNN, strNN, 800/REDUCTION_FACTOR, -0.5, 799.5, 800/REDUCTION_FACTOR, -0.5, 799.5 );

    TString strPtPt = Form("hist2D_PtPt_%s", strWinDescr.Data());
//  !!! TMP COMMENTED:  hist2D_PtPt = new TH2D( strPtPt, strPtPt, 400/REDUCTION_FACTOR, 0, 2, 400/REDUCTION_FACTOR, 0, 2);
//    hist2D_PtPt = new TH2D( strPtPt, strPtPt, 60, 0.8, 3.8, 60, 0.8, 3.8);
//    hist2D_PtPt = new TH2D( strPtPt, strPtPt, 80, 0.8, 2.8, 80, 0.8, 2.8);
//    hist2D_PtPt = new TH2D( strPtPt, strPtPt, 160, 0.8, 4.8, 160, 0.8, 4.8);
    hist2D_PtPt = new TH2D( strPtPt, strPtPt, 180, 0.2, 2.0, 180, 0.2, 2.0);

    TString strPtN = Form("hist2D_PtN_%s", strWinDescr.Data());
    hist2D_PtN = new TH2D( strPtN, strPtN, 1000/REDUCTION_FACTOR, -0.5, 999.5, 400/REDUCTION_FACTOR, 0, 2);

    TString strSumPtN = Form("hist2D_sumPtN_%s", strWinDescr.Data());
    hist2D_sumPtN = new TH2D( strSumPtN, strSumPtN, 1000/REDUCTION_FACTOR, -0.5, 999.5, 400/REDUCTION_FACTOR, 0, 2);

    //QA histos:
    TString strEstPerc = Form("hist1D_EstimatorEntries_%s;percentile;entries", strWinDescr.Data());
    hist1D_EstimatorEntries = new TH1D( strEstPerc, strEstPerc, 20001/REDUCTION_FACTOR, -0.5, 400.5); //1000.5);

    TString strMultDistrF_name  = Form("hist1D_multDistrF_%s", strWinDescr.Data());
    TString strMultDistrF_title = Form("hist1D_multDistrF_%s;n tracks;n events", strWinDescr.Data());
    hist1D_multDistrF = new TH1D( strMultDistrF_name, strMultDistrF_title, 1001, -0.5, 1000.5);

    TString strMultDistrB_name  = Form("hist1D_multDistrB_%s", strWinDescr.Data());
    TString strMultDistrB_title = Form("hist1D_multDistrB_%s;n tracks;n events", strWinDescr.Data());
    hist1D_multDistrB = new TH1D( strMultDistrB_name, strMultDistrB_title, 1001, -0.5, 1000.5);

    TString strPtF_name = Form("hist1D_PtF_%s", strWinDescr.Data());
    TString strPtF_title = Form("hist1D_PtF_%s;#LTp_{T}#GT Forward;n events", strWinDescr.Data());
    hist1D_QA_PtF = new TH1D( strPtF_name, strPtF_title, 3100/*2002*/, -1.1, 4.1 );

    TString strPtB_name = Form("hist1D_PtB_%s", strWinDescr.Data());
    TString strPtB_title = Form("hist1D_PtB_%s;#LTp_{T}#GT Backward;n events", strWinDescr.Data());
    hist1D_QA_PtB = new TH1D( strPtB_name, strPtB_title, 3100/*2002*/, -1.1, 4.1 );

    //moved here to avoid seg. fault! (25.06.2016)
    histos_BS_NN.InitHistos( cW, "BS_NN", cBinMin, cBinMax, etaW, phiW, ptW, 1000 );
    histos_BS_PtPt.InitHistos( cW, "BS_PtPt", cBinMin, cBinMax, etaW, phiW, ptW );

    if ( MAX_N_EVENTS_FOR_BOOTSTRAP > 0 )
    {
        doBootstrap = true;

        SAMPLING_nF = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        SAMPLING_nB = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        SAMPLING_PtF = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];
        SAMPLING_PtB = new double[MAX_N_EVENTS_FOR_BOOTSTRAP];

        //        histos_BS_NN.InitHistos( cW, "BS_NN", strWinDescr.Data(), 1000 );
        //        histos_BS_PtPt.InitHistos( cW, "BS_PtPt", strWinDescr.Data() );

        //        histos_SS_NN.InitHistos( "SS_NN", strWinDescr.Data(), 1000 );
        //        histos_SS_PtPt.InitHistos( "SS_PtPt", strWinDescr.Data() );

        //        TString strBS_NN_name  = Form("hist1D_bCorr_BS_NN_%s", strWinDescr.Data());
        //        TString strBS_NN_title = Form("hist1D_bCorr_BS_NN_%s;bCorr;entries", strWinDescr.Data());
        //        hist1D_bCorr_BS_NN = new TH1D( strBS_NN_name, strBS_NN_title, 2000, -1, 1 );

        //        TString strBS_PtPt_name  = Form("hist1D_bCorr_BS_PtPt_%s", strWinDescr.Data());
        //        TString strBS_PtPt_title = Form("hist1D_bCorr_BS_PtPt_%s;bCorr;entries", strWinDescr.Data());
        //        hist1D_bCorr_BS_PtPt = new TH1D( strBS_PtPt_name, strBS_PtPt_title, 2000, -1, 1 );

    }

}

void WinPair::initRunByRunHistos( int _nRuns, const int *runListNumbers )
{
    nRuns = _nRuns;
    hist1D_multDistr_RunByRun_F = new TH1D*[nRuns];
    hist1D_multDistr_RunByRun_B = new TH1D*[nRuns];
    hist1D_avPtDistr_RunByRun_F = new TH1D*[nRuns];
    hist1D_avPtDistr_RunByRun_B = new TH1D*[nRuns];

    TString str_hist_name;
    for ( int r = 0; r < nRuns; r++ )
    {
        int runNumber = runListNumbers[r];
        TString str_hist_title  = Form("%d", runNumber);

        TString strWinDescr = Form("run_%d_c%.1f-%.1f_etaW_%d_phiW_%d", runNumber, cBinMin, cBinMax, etaW, phiW);
        // !!! modify if >1 pt bins!
        if ( ptW > 1 )
            strWinDescr = Form("%s_ptW_%d", strWinDescr.Data(), ptW);

        //mult
        str_hist_name  = Form("hist1D_multDistr_winF_%s", strWinDescr.Data());
        hist1D_multDistr_RunByRun_F[r] = new TH1D( str_hist_name, ";n tracks;n events", 1000, -0.5, 1000-0.5 );

        str_hist_name  = Form("hist1D_multDistr_winB_%s", strWinDescr.Data());
        hist1D_multDistr_RunByRun_B[r] = new TH1D( str_hist_name, ";n tracks;n events", 1000, -0.5, 1000-0.5 );

        hist1D_multDistr_RunByRun_F[r]->SetTitle( str_hist_title );
        hist1D_multDistr_RunByRun_B[r]->SetTitle( str_hist_title );

        //av pT
        str_hist_name  = Form("hist1D_avPtDistr_winF_%s", strWinDescr.Data());
        hist1D_avPtDistr_RunByRun_F[r] = new TH1D( str_hist_name, ";#LTp_{T}#GT;n events", 1000, 0, 2 );

        str_hist_name  = Form("hist1D_avPtDistr_winB_%s", strWinDescr.Data());
        hist1D_avPtDistr_RunByRun_B[r] = new TH1D( str_hist_name, ";#LTp_{T}#GT;n events", 1000, 0, 2 );

        hist1D_avPtDistr_RunByRun_F[r]->SetTitle( str_hist_title );
        hist1D_avPtDistr_RunByRun_B[r]->SetTitle( str_hist_title );
    }

}


//void WinPair::fill(Float_t cPerc, UShort_t nF, UShort_t nB, Float_t PtF, Float_t PtB, int treeId )
//void WinPair::fill(Float_t cPerc, Float_t nF, Float_t nB, Float_t PtF, Float_t PtB, int treeId )
void WinPair::fill(Double_t cPerc, Double_t nF, Double_t nB, Double_t PtF, Double_t PtB, int treeId )
{
    //check if we in centrality bin
    if ( cPerc < cBinMin || cPerc > cBinMax )
        return;

    if (initHistos)
    {
        //QA fills:
        hist1D_EstimatorEntries->Fill(cPerc);
        hist1D_multDistrF->Fill(nF);
        hist1D_multDistrB->Fill(nB);

        hist1D_QA_PtF->Fill( PtF );
        hist1D_QA_PtB->Fill( PtB );

        //QA run-by-run
        if ( hist1D_multDistr_RunByRun_F != 0x0 )
        {
            //F
            hist1D_multDistr_RunByRun_F[treeId]->Fill( nF);
            if ( nF > 0 )
                hist1D_avPtDistr_RunByRun_F[treeId]->Fill( PtF );

            //B
            hist1D_multDistr_RunByRun_B[treeId]->Fill( nB);
            if ( nB > 0 )
                hist1D_avPtDistr_RunByRun_B[treeId]->Fill( PtB );
        }
    }


    //NN
    if (initHistos)
        hist2D_NN->Fill( nF, nB );

    NN_nF       += nF;
    NN_nB       += nB;
    NN_nF_nB    += nF*nB;
    NN_nF2      += nF*nF;
    NN_nB2      += nB*nB;

    if ( doBootstrap )
    {
        SAMPLING_nF[NN_Nevents] = nF;
        SAMPLING_nB[NN_Nevents] = nB;
    }

    NN_Nevents++;


    //PtPt
    if ( nF > 0 && nB > 0 )
    {
        if (initHistos)
            hist2D_PtPt->Fill( PtF, PtB );

        PtPt_PtF += PtF;
        PtPt_PtB += PtB;
        PtPt_PtF_PtB += PtF*PtB;
        PtPt_PtF2 += PtF*PtF;
        PtPt_PtB2 += PtB*PtB;

        if ( doBootstrap )
        {
            SAMPLING_PtF[PtPt_Nevents] = PtF;
            SAMPLING_PtB[PtPt_Nevents] = PtB;
        }

        PtPt_Nevents++;
    }

    //PtN
    if ( nB > 0 )
    {
        if (initHistos)
            hist2D_PtN->Fill( nF, PtB );

        PtN_nF += nF;
        PtN_PtB += PtB;
        PtN_nF_PtB += nF*PtB;
        PtN_nF2 += nF*nF;
        PtN_PtB2 += PtB*PtB;
        PtN_Nevents++;
    }

    // sumPtN
    if ( nB > 0 )
    {
        if (initHistos)
            hist2D_sumPtN->Fill( nF, PtB*nB );

        sumPtN_nF += nF;
        sumPtN_sumPtB += PtB*nB;
        sumPtN_nF_sumPtB += nF*PtB*nB;
        sumPtN_nF2 += nF*nF;
        sumPtN_PtB2 += PtB*nB*PtB*nB;
        sumPtN_Nevents++;
    }

}

void WinPair::calcCorrCoeffs()
{
    corrInfo_NN    = _calc( NN_nF, NN_nB, NN_nF_nB, NN_nF2, NN_nB2, NN_Nevents, 0 );
    corrInfo_PtPt  = _calc( PtPt_PtF, PtPt_PtB, PtPt_PtF_PtB, PtPt_PtF2, PtPt_PtB2, PtPt_Nevents, 1 );
    corrInfo_PtN   = _calc( PtN_nF, PtN_PtB, PtN_nF_PtB, PtN_nF2, PtN_PtB2, PtN_Nevents, 2, true );

    corrInfo_sumPtN   = _calc( sumPtN_nF, sumPtN_sumPtB, sumPtN_nF_sumPtB, sumPtN_nF2, sumPtN_PtB2, sumPtN_Nevents, 3, true );
}

void WinPair::writeHistos()
{
    if ( !initHistos )
    {
        cout << ">>> !!! histos are not initialized!" << endl;
        return;
    }

    hist2D_NN->Write();
    hist2D_PtPt->Write();
    hist2D_PtN->Write();
    hist2D_sumPtN->Write();

    hist1D_multDistrF->Write();
    hist1D_multDistrB->Write();

    hist1D_QA_PtF->Write();
    hist1D_QA_PtB->Write();

    if ( nRuns > 0 ) //we have run-by-run histos
    {
        for ( int treeId = 0; treeId < nRuns; treeId++ )
        {
            hist1D_multDistr_RunByRun_F[treeId]->Write();
            hist1D_multDistr_RunByRun_B[treeId]->Write();
            hist1D_avPtDistr_RunByRun_F[treeId]->Write();
            hist1D_avPtDistr_RunByRun_B[treeId]->Write();
        }
    }
    //                    wins[cW][cBin][etaW][phiW].hist2D_NN->ProfileX()->Write();
    //                    wins[cW][cBin][etaW][phiW].hist2D_PtPt->ProfileX()->Write();
    //                    wins[cW][cBin][etaW][phiW].hist2D_PtN->ProfileX()->Write();

}

CorrCoeffInfo WinPair::_calc(const double &F, const double &B, const double &FB
                             , const double &F2, const double &B2, const int &nEvents
                             , const int corrType , const bool ifRel )
{
    // the code is independent of corrType,
    // until strongly intensive part!!

    CorrCoeffInfo corrInfo;

    if ( nEvents <= 0 )
        return corrInfo;
    double meanF     =  F   / nEvents;     //printf( "<F> = %f\n", meanF );
    double meanB     =  B   / nEvents;     //printf( "<B> = %f\n", meanB );
    double meanFB    =  FB  / nEvents;    //printf( "<FB> = %f\n", meanFB );
    double meanF2    =  F2  / nEvents;    //printf( "<FB> = %f\n", meanFB );
    double meanB2    =  B2  / nEvents;    //printf( "<FB> = %f\n", meanFB );

    //        printf( "nEvents = %d, ", nEvents );
    //        printf( "meanF = %f, ", meanF  );
    //        printf( "meanB = %f, ", meanB  );
    //        printf( "meanFB = %f, ",  meanFB );
    //        printf( "meanF2 = %f\n, ",  meanF2 );

    double numerator = meanFB - meanF * meanB;

    double denominator_bCorr = meanF2 - meanF*meanF;
//        double denominator_bCorr = meanB2 - meanB*meanB;
//        double denominator_bCorr = sqrt( (meanF2 - meanF*meanF)*(meanB2 - meanB*meanB) );

    double denominator_C2 = meanF * meanB;
    //    if (type == 0) // bCorr
    //        denominator = meanF2 - meanF*meanF;
    //    else if (type == 1) // C2
    //        denominator = meanF * meanB;

    //bcorr
    double bcorr = -1000;
    if ( denominator_bCorr != 0 )
        //    if ( fabs(denominator_bCorr) > 0.0001 )
    {
        bcorr = numerator / denominator_bCorr;
        if (ifRel)
            bcorr *= meanF/meanB;
    }
    //    else
    if ( fabs(denominator_bCorr) < 0.00001 )
        if(1)cout << "!!!!! WARNING denominator_bCorr=" << denominator_bCorr
                  << ", meanF2=" << meanF2 << ", meanF*meanF=" << meanF*meanF
                  << ", numerator = " << numerator
                  << ", meanFB=" << meanFB << ", meanF*meanB=" << meanF*meanB << endl;

    if ( bcorr < -10 )
    {
        cout << "!!!!! AHTUNG! bcorr<10! bcorr=" << bcorr;// << endl;
        cout << "  denominator_bCorr=" << denominator_bCorr
                          << ", meanF2=" << meanF2 << ", meanF*meanF=" << meanF*meanF
                          << ", numerator = " << numerator
                          << ", meanFB=" << meanFB << ", meanF*meanB=" << meanF*meanB << endl;
    }



    //C2
    double C2 = -1000;
    if ( denominator_C2 !=0 )
        C2 = numerator / denominator_C2;


    //strongly intensive:
    double omegaF = meanF2-meanF*meanF;
    double omegaB = meanB2-meanB*meanB;

    if ( corrType == 3 ) //sumPtN
    {
        //A=backward win (B), B=forward win (F)
        double delta = meanF*omegaB - meanB*omegaF;
        double sigma = meanF*omegaB + meanB*omegaF
                - 2 * ( meanFB - meanF*meanB );
        delta /= 1.;
        sigma /= 1.;
        //need to normalized them properly!
        // -> to be done in future!

    }



    //fill corrInfo data:
    corrInfo.avF = meanF;
    corrInfo.avB = meanB;
    corrInfo.avFB = meanFB;
    corrInfo.avF2 = meanF2;

    corrInfo.DFB = numerator;
    corrInfo.DF2 = denominator_bCorr;

    corrInfo.bCorr = bcorr;
    corrInfo.C2 = C2;

    return corrInfo;
}



void WinPair::performSubsampling(int corrType, int subsamplingType, int n_samples )
{
    // subsamplingType=0 - bootstrap
    // subsamplingType=1 - simple subsampling

    Long64_t _nDataEvents = 0;

    //corr type for subsampling:
    if ( corrType == 0 )
        _nDataEvents = NN_Nevents;
    else if ( corrType == 1 )
        _nDataEvents = PtPt_Nevents;
    else
        cout << ">>> WARNING: up to now - bootstrapping for NN and PtPt only!" << endl;
    //        else if ( corrType == 2 )
    //            _nDataEvents = PtN_Nevents;

    //    cout << ">>> PtPt_Nevents=" << PtPt_Nevents << endl;

    //for simple subsampling
    int subset_size = _nDataEvents / n_samples;
    int subset_begin    = 0;
    int subset_end      = 0;

    int _nEvents = (subsamplingType == 0 ? _nDataEvents : subset_size);
    //    cout << ">>> _nEvents = " << _nEvents << ", subset_size*n_samples = " << subset_size*n_samples << endl;

    //    int n_samples_fixed /*= n_samples */ = _nDataEvents / subset_size;
    int n_samples_fixed = (subsamplingType == 0 ? n_samples : _nDataEvents / subset_size );

    //    cout << ">>> n_samples_fixed = " << n_samples_fixed << ", _nDataEvents / subset_size = " << _nDataEvents / subset_size << endl;


    //shuffling array indeces:
    Long64_t *indexArr = 0x0;
    if ( subsamplingType == 1 ) //use shuffling for subsampling!
    {
        indexArr = new Long64_t[_nDataEvents];
        for ( int i = 0; i < _nDataEvents; i++ )
            indexArr[i] = i;
        random_shuffle(&indexArr[0], &indexArr[_nDataEvents]); // end must be n+1, not n!!!
    }


    for ( int t = 0; t < n_samples_fixed; t++ )
    {
        double BS_F = 0;
        double BS_B = 0;
        double BS_FB = 0;
        double BS_F2 = 0;
        double BS_B2 = 0;

        double *F;
        double *B;

        //                cout << ">>> start subsampling event loop..." << endl;
        if ( subsamplingType == 1 )
        {
            subset_begin    = t*subset_size;
            subset_end      = (t+1)*subset_size;
        }

        //                cout << ">>> _nEvents=" << _nEvents << endl;
        for ( int ev = 0; ev < _nEvents; ev++ )
        {
            //generate id for event:
            int id = -1;
            if ( subsamplingType == 0 ) // BS
                id = TMath::Nint( gRandom->Uniform( -0.5, _nDataEvents-0.5) );
            else if ( subsamplingType == 1 ) //simple subsampling
                //                id = TMath::Nint( gRandom->Uniform( subset_begin-0.5, subset_end-0.5 ) ); // ????????
                id = subset_begin+ev;

            int shuffledId = ( subsamplingType == 1 ? indexArr[id] : id ); //use shuffling for subsampling but not for bootstrap!
            if ( corrType == 0 )
            {
                F = &SAMPLING_nF[shuffledId];
                B = &SAMPLING_nB[shuffledId];
            }
            else if ( corrType == 1 )
            {
                F = &SAMPLING_PtF[shuffledId];
                B = &SAMPLING_PtB[shuffledId];
                //                        cout << ">>> *F = " << *F << ", *B = " << *B << ", ";
            }

            BS_F += *F;
            BS_B += *B;
            BS_FB += (*F)*(*B);
            BS_F2 += (*F)*(*F);
            BS_B2 += (*B)*(*B);

            //            cout << ">>> BS_F = " << BS_F << ", BS_B = " << BS_B
            //                << ">>> BS_FB = " << BS_FB << ", BS_F2 = " << BS_F2
            //                << ", _nDataEvents = " << _nDataEvents << endl;
            //            int a;
            //            cin >> a;
        }
        //        cout << ">>> _nDataEvents = " << _nDataEvents << endl;

        CorrCoeffInfo BS_corrInfo = _calc( BS_F, BS_B, BS_FB, BS_F2, BS_B2, _nEvents, corrType );
        //        double BS_bCorr = BS_corrInfo.bCorr;
        //                cout << ">>> BS_bCorr = " << BS_bCorr << endl;
        //double BS_C2    = _calc( 1, BS_F, BS_B, BS_FB, BS_F2, _nDataEvents );


        if ( corrType == 0 )
            //hist1D_bCorr_BS_NN->Fill( BS_bCorr );
            histos_BS_NN.FillHistos( BS_corrInfo );
        else if ( corrType == 1 )
            //            hist1D_bCorr_BS_PtPt->Fill( BS_bCorr );
            histos_BS_PtPt.FillHistos( BS_corrInfo );

    }

    if ( indexArr )
        delete [] indexArr;
    //    cout << ">>> END BS " << endl;
}

// ##### Event Mixing
CorrCoeffInfo WinPair::performEventMixing(int corrType) //, int subsamplingType, int n_samples )
{
    Long64_t _nDataEvents = 0;

    //corr type for subsampling:
    if ( corrType == 0 )
        _nDataEvents = NN_Nevents;
    else if ( corrType == 1 )
        _nDataEvents = PtPt_Nevents;
    else
        cout << ">>> WARNING: up to now - bootstrapping for NN and PtPt only!" << endl;
    //        else if ( corrType == 2 )
    //            _nDataEvents = PtN_Nevents;

    //    cout << ">>> PtPt_Nevents=" << PtPt_Nevents << endl;

    double MIX_F = 0;
    double MIX_B = 0;
    double MIX_FB = 0;
    double MIX_F2 = 0;
    double MIX_B2 = 0;

    double *F;
    double *B;


    // FLAG EVENT MIXING TYPE:
    bool flagEvMixType = 0; // 0 - prev. event, 1 - shuffled indeces

    //shuffling array indeces - for event mixing with random event
    Long64_t *indexArr = 0x0;
    if ( flagEvMixType == 1 ) //use shuffling for mixing with random event!
    {
        indexArr = new Long64_t[_nDataEvents];
        for ( int i = 0; i < _nDataEvents; i++ )
            indexArr[i] = i;
        random_shuffle(&indexArr[0], &indexArr[_nDataEvents]); // end must be n+1, not n!!!
    }


    for ( int ev = 0; ev < _nDataEvents; ev++ )
    {
        int ev1 = ev;
        int ev2 = ( flagEvMixType == 0 ? ev-1 : indexArr[ev] );
        while ( ev2 == ev1 )
        {
            cout << " ! warning: same event chosen, apply randomization..." << endl;
            ev2 = TMath::Nint( gRandom->Uniform( -0.5, _nDataEvents-0.5) );
        }

        //fix for the first event in case of ev2=ev-1: take 0 and 1 events for mixing
        if ( ev2 < 0 )
        {
            ev1 = 0;
            ev2 = 1;
        }

        if ( corrType == 0 )
        {
            F = &SAMPLING_nF[ev1];
            B = &SAMPLING_nB[ev2];
        }
        else if ( corrType == 1 )
        {
            F = &SAMPLING_PtF[ev1];
            B = &SAMPLING_PtB[ev2];
            //                        cout << ">>> *F = " << *F << ", *B = " << *B << ", ";
        }

        MIX_F += *F;
        MIX_B += *B;
        MIX_FB += (*F)*(*B);
        MIX_F2 += (*F)*(*F);
        MIX_B2 += (*B)*(*B);
    }
    CorrCoeffInfo MIX_corrInfo = _calc( MIX_F, MIX_B, MIX_FB, MIX_F2, MIX_B2, _nDataEvents, corrType );

    return MIX_corrInfo;
//    if ( corrType == 0 )
//        //hist1D_bCorr_BS_NN->Fill( BS_bCorr );
//        histos_BS_NN.FillHistos( BS_corrInfo );
//    else if ( corrType == 1 )
//        //            hist1D_bCorr_BS_PtPt->Fill( BS_bCorr );
//        histos_BS_PtPt.FillHistos( BS_corrInfo );

}



// ##### CentralityOccupancy

CentralityOccupancy::CentralityOccupancy() :
    cBinMin(0)
  , cBinMax(100)
  , nEventsV0M(0)
  , nEventsZDCZEM(0)
  , nEventsV0M_and_ZDCZEM(0)
{}
void CentralityOccupancy::fill( Double_t cPercV0M, Double_t cPercZDCvsZNA )
{
    //check if we in centrality bin
    bool isV0M = false;
    if ( cPercV0M > cBinMin && cPercV0M < cBinMax )
    {
        nEventsV0M++;
        isV0M = true;
    }

    bool isZDCvsZNA = false;
    if ( cPercZDCvsZNA > cBinMin && cPercZDCvsZNA < cBinMax )
    {
        nEventsZDCZEM++;
        isZDCvsZNA = true;
    }

    if ( isV0M && isZDCvsZNA )
        nEventsV0M_and_ZDCZEM++;
}




void GraphsCorrInfo::WriteGraphs()
{
    TString strDirName = Form( "dir_%s" , gr_bCorr->GetName() );
    gFile->mkdir( strDirName.Data() );
    gFile->cd( strDirName.Data() );

    gr_bCorr->Write();
    gr_C2   ->Write();

    gr_DFB -> Write();
    gr_DF2 -> Write();
    gr_avFB-> Write();
    gr_avF2-> Write();
    gr_avF -> Write();
    gr_avB -> Write();

    gFile->cd();

}
