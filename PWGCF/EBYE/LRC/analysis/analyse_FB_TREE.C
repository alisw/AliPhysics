#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>

#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include "TObject.h"
#include "TGraphErrors.h"
#include "TH1.h"
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
#include "TString.h"

//#include "TChain.h"
//#include "TFileCollection.h"

#include "SupplementaryClasses.h"


// !!! UNCOMMENT IF NEED UNFOLD->
//#include "/opt/RooUnfold/RooUnfold/src/RooUnfold.h"
// !!! UNCOMMENT IF NEED UNFOLD->
//#include "/opt/RooUnfold/RooUnfold/src/RooUnfoldResponse.h"

#include <iostream>
#include <fstream>
//#include <string>

using namespace std;

#include "utils.C"


// ########## event selection parameters
const float kVertexZcut = 8;
//    if ( isAnalysing502 )
//        kVertexZcut = 7;
const bool USE_VERTEXZ_CUT = 1;
const int flag_V0M_ZDC_CL1 = 0; // V0M, ZDCvsZEM, CL1, |V0M - CL1|<7.5

bool checkEventSelection( float &cEstimator
                          , const Bool_t evSel_isPileupSPD
                          , const Bool_t evSel_isPileupFromMV
                          , const Bool_t evSel_vtxSPD_isFromVertexer3D
                          , const Bool_t evSel_vtxSPD_isFromVertexerZ
                          , const Bool_t evSel_zRes_above_fMaxResol
                          , const Bool_t evSel_diff_vZ_SPD_tracks_above_05
                          , const Bool_t evSel_comb_cut_on_abs_and_nsigma_dist
                          , const Bool_t evSel_too_many_TPC_clusters

                          , const UShort_t evSel_multTrk_bit32
                          , const UShort_t evSel_multTrk_bit32_TOF

                          , const Float_t &brZEMvsZDC
                          , const Float_t &brV0M
                          , const Float_t &brCL1
                          , const Float_t &br_vertexZ
                          , int &nRejectedByPileup
                          )
{
    cEstimator = -1;
    
    // check isPileupSPD
    if(0)if ( evSel_isPileupSPD == 1 )
    {
        nRejectedByPileup++;
        return false;
    }

    // !!! try other checks!
    if(1)if ( evSel_zRes_above_fMaxResol || evSel_diff_vZ_SPD_tracks_above_05 || evSel_comb_cut_on_abs_and_nsigma_dist || evSel_too_many_TPC_clusters )
        return false;

    // !!! try cut by hand using TOF timing!!!
    if(1)if ( evSel_multTrk_bit32_TOF < 200./1400*(evSel_multTrk_bit32-100)
              || evSel_multTrk_bit32_TOF < 900./2000*(evSel_multTrk_bit32-400))
        return false;

    
    // vertex Z cut!!!
    if ( USE_VERTEXZ_CUT && fabs(br_vertexZ) > kVertexZcut
         //                  || fabs(br_vertexZ) < 5 // TEST INFLUENCE FROM TPC MEMBRANE!
//                  || br_vertexZ > 0 //> -1 //reject vertZ<0 ! -> keep vertZ>0
         )
        return false;
    
    //    cout << ">>>    vertexZ bool decision =" << (USE_VERTEXZ_CUT && fabs(br_vertexZ) > kVertexZcut
    //            && fabs(br_vertexZ) < 7.5) << endl;
    
    if ( flag_V0M_ZDC_CL1 == 0 )
    {
        if ( brV0M > 90 ) //V0M cut
            return false;
        cEstimator = brV0M;
    }
    else if ( flag_V0M_ZDC_CL1 == 1 )
    {
        if ( brZEMvsZDC > 50 ) //ZDCvsZEM cut
            return false;
        cEstimator = brZEMvsZDC;
    }
    else if ( flag_V0M_ZDC_CL1 == 2 )
    {
        if ( brCL1 > 90 ) //CL1 cut
            return false;
        cEstimator = brCL1;
    }
    else if ( flag_V0M_ZDC_CL1 == 3 ) //spec check!
    {
        if ( brV0M > 90 ) //V0M cut
            return false;
        if ( brCL1 > 90 ) //CL1 cut
            return false;
        if ( fabs( brV0M - brCL1 ) > 7.5 ) //|V0M - CL1|<7.5 cut from A.Dobrin
            return false;

        cEstimator = brV0M;
    }

    return true;
}



TTree *getTreeFromFile( TFile *myFile, int listIdByHand = -1 )
{
    int listId;
    TList *listKeys = myFile->GetListOfKeys();
    
    if ( listIdByHand >= 0 ) //set listId by hand!!!
    {
        listId = listIdByHand;
    }
    else //automatic search for LRC task list
    {
        listId = 0;//1;
        // !!! listId can be different for 2.76 and 5.02!
        //    if ( isAnalysing502 )
        //        listId = 1;
        
        TString listName = Form("%s",listKeys->At(listId)->GetName() );
        while ( 1 && !listName.Contains("PWGCFLRC") )
        {
            listId++;
            listName = Form("%s",listKeys->At(listId)->GetName() );
        }
        //return;
        //SET listId BY HAND:
        //    listId = 0;
        
        //return;
    }
    cout << "going into list: " << listKeys->At(listId)->GetName() << endl;
    myFile->cd( listKeys->At(listId)->GetName() );
    
    TList *listKeys2 = gDirectory->GetListOfKeys();
    TList *myTask = (TList*) gDirectory->Get( listKeys2->At(0)->GetName() );
    TTree *t1 = (TTree*)  myTask->FindObject( "treeFB" );
    
    return t1;
}


void get_run_list( int &nFiles, TString *&runListFullPath, int *&runListNumbers, const char *dirBase = "", const char *taskdirname=""
        , const char *begins="AnalysisResults", const char *ext=".root" )
{
    //    runList = new TString[200];
    int counter = 0;
    
    TString dirname = Form( "%s/%s", dirBase, taskdirname);
    
    TSystemDirectory dir( dirname.Data(), dirname.Data() );
    TList *files = dir.GetListOfFiles();
    
    dir.ls();
    
    cout << "Runlist: "  << endl;
    if (files)
    {
        TSystemFile *systFile;
        TString fname;
        TIter next(files);
        while ((systFile=(TSystemFile*)next()))
        {
            fname = systFile->GetName();
            if (!systFile->IsDirectory() && fname.BeginsWith(begins) && fname.EndsWith(ext))
            {
                //strip run number:
                TString runName = fname;
                runName.Remove(0,16);
                runName.Remove( runName.Length()-5, 5);
                runListNumbers[nFiles+counter] = runName.Atoi();
                
                //remember full path to file!!!
                runListFullPath[nFiles+counter] = Form( "%s/%s", dirname.Data(), fname.Data() ); //runName.Atoi();
                
                cout << runListFullPath[nFiles+counter] << endl;
                counter++;
            }
        }
    }
    cout << endl;
    nFiles += counter;
}




void analyse_FB_TREE( const char *dirBase
                      , int fileIdByHand = -1
        , int nEventsByHand = -1
        , int LIST_ID_BY_HAND = -1
        , int WHICH_BRANCH_BY_HAND = 0
        , int ptW = 0
        )
// !!! UNCOMMENT IF NEED UNFOLD
//        , RooUnfoldResponse *response_PtPt = 0x0 ) // TString inputFileName = "MergedOutput.root")
{
//    WHICH_BRANCH_BY_HAND = 0;//0;//2;
    cout << "starting analyse_FB_TREE..." << endl;
    
    //int LIST_ID_BY_HAND = -1;
//    LIST_ID_BY_HAND = 0;
    
    // !!! important for branch V0M below!
    bool isAnalysing502 = 1;
    if ( isAnalysing502 )
        LIST_ID_BY_HAND = 1; // by hand - in order to omit MultSelection list (id=0)!


    //    bool isAMPT_kine = 0;
//        LIST_ID_BY_HAND = 2;//4;
    
    
    //    const int nFiles = 24+21; //30; //24+22; //3;//10;
    //    TFile *myFiles[nFiles];
    
    // FULL STATISTICS
    // LHC10h
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all/MergedOutput.root" );
    //    myFiles[1] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFminus_all/MergedOutput.root" );
    
    //HIJING LHC11a10a_bis RECO level
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_28_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco_try5/MergedOutput.root" );
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_03_06_PbPb_MCAOD_LHC11a10a_bis_AOD162_TEST_Efficiency_kine_vs_reco/MergedOutput.root" );
    
    //AMPT AMPT_LHC13f3c_StrMelt_ON_rescatON_8_eta_wins KINE level
    //    isAMPT_kine = 1;
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_29_PbPb_AMPT_LHC13f3c_StrMelt_ON_rescatON_8_eta_wins/merged_results_with_scaled_factor.root" );
    
    
    
    
    //    const char *dirBase = "/Users/macbook/alice/aliceAnalysis/results";
    
    //LHC10h plus, minus:
    //    const char *taskdirname = "task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all";
    //    const char *taskdirname2 = "task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFminus_all";
    //HIJING:
    //const char *taskdirname = "task_2016_03_09_PbPb_MCAOD_LHC11a10a_bis_AOD162_Efficiency_kine_vs_reco_fixedCentrBins_try2";
    
    //LHC10h systematics study (March 2016)
    //        const char *taskdirname = "task_2016_03_21_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_SYSTEMATICS_filterBits_TPCclust";
    
    //LHC11a10a_bis_AOD162_kine_PLUS_reco_in_one_tree - TRY TO FILL RESPONSE MATRIX_FOR UNFOLDING:
    //    const char *taskdirname = "task_2016_04_28_PbPb_MCAOD_LHC11a10a_bis_AOD162_kine_PLUS_reco_in_one_tree";
    //LHC11a10a_bis_AOD162_kine_PLUS_reco_in_one_tree - TRY TO FILL RESPONSE MATRIX_FOR UNFOLDING:
    //    const char *taskdirname = "task_2016_04_29_PbPb_MCAOD_LHC12a11aANDb_AMPT_AOD157_kine_PLUS_reco_in_one_tree";
    
    //LHC12a1 AMPT c30_40 bit128 AOD081 kine_PLUS_reco - for kine/reco/primary/secondary studies
    //    const char *taskdirname = "task_2016_05_04_PbPb_AMPT_LHC12a1_c30_40_bit128_AOD081_kine_PLUS_reco_in_tree_SPARSE_try3";
    
    //LHC12a11abc... bit768_AOD157_WEAK_SEC
    //    const char *taskdirname = "task_2016_05_09_PbPb_AMPT_LHC12a11a_c0_5_bit768_AOD157_WEAK_SEC_in_tree";
    //    const char *taskdirname = "task_2016_05_09_PbPb_AMPT_LHC12a11d_c20_30_bit768_AOD157_WEAK_SEC_in_tree";
    //    const char *taskdirname = "task_2016_05_09_PbPb_AMPT_LHC12a11a_c0_5_bit768_AOD157_WEAK_SEC_in_tree";
    //    const char *taskdirname2 = "task_2016_05_09_PbPb_AMPT_LHC12a11b_c5_10_bit768_AOD157_WEAK_SEC_in_tree";
    
    //LHC12a11f for IneffByHand (1D hist) bit768_AOD157_WEAK_SEC
    //    const char *taskdirname = "task_2016_05_21_PbPb_AMPT_LHC12a11f_c40_50_bit768_AOD157_WEAK_SEC_in_tree_INEFF_1D_HIST";
    //    const char *taskdirname = "task_2016_05_22_PbPb_AMPT_LHC12a11f_c40_50_bit768_AOD157_WEAK_SEC_in_tree_INEFF_1D_HIST_FROM_HIJING";
    //    const char *taskdirname = "task_2016_05_24_PbPb_AMPT_LHC12a11f_c40_50_bit128_AOD157_WEAK_SEC_in_tree_INEFF_3D_HIJING_bit768";
    //    const char *taskdirname = "task_2016_05_25_PbPb_AMPT_LHC12a11f_c40_50_bit128_AOD157_WEAK_SEC_in_tree_INEFF_3D_HIJING_bit128";
    //    const char *taskdirname = "task_2016_05_29_PbPb_AMPT_LHC12a11f_c40_50_bit768_AOD157_WEAK_SEC_in_tree_INEFF_3D_HIJING_bit768";
    
    //LHC10h MAIN DATA June 2016, bit 768, also bit 128
//    const char *taskdirname = "task_2016_06_14_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_BIT768_FB_TREE_FINAL";
    //        const char *taskdirname = "task_2016_06_14_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_BIT128_FB_TREE_FINAL";
    
    //LHC10h MAIN DATA June 2016, bit 768, PHI 8 WINS
    //    const char *taskdirname = "task_2016_06_16_PbPb_Data_LHC10h_AOD160_eW_04_08_phi8_pt02_20_BIT768_FB_TREE_PHI_WINS";
    
    // LHC10h DATA June 2016, bit 768, PHI 32 WINS WITH PT BINNING, MOST CENTRAL (2% by V0M)
    //    OLD: not full stat! //const char *taskdirname = "task_2016_06_26_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt5_BIT768_PHI_PT_MANY_WINS";
//    const char *taskdirname = "task_2016_07_02_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt6_BIT768_PHI_PT_MANY_WINS";
    //    const char *taskdirname = "task_2016_07_05_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt6_BIT128_PHI_PT_MANY_WINS";
    // 36 PHI BINS 29.07.2016
    //    const char *taskdirname = "task_2016_07_29_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi36_pt4_BIT768_cPos0_5_20_40_cW2";
    // 20 PHI BINS 31.07.2016 (QA check)
    //    const char *taskdirname = "task_2016_07_30_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi20_pt1_08_10_BIT768_c0_2";
    //64 PHI BINS 31.07.2016
    //    const char *taskdirname = "task_2016_07_31_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi64_pt2_BIT768_c0_2";
    
    //23 PHI BINS 01.08.2016
    //    const char *taskdirname = "task_2016_08_01_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi23_pt2_BIT768_c0_2";
    
    //2.08.2016: DATA PbPb 32 PHI BINS, one pT bin 0.8-10.0 GeV/c, CENTRALITY 0-15% (V0M)
    //    const char *taskdirname = "task_2016_08_01_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt1_08_100_BIT768_c0_15";
    
    //2.08.2016: DATA pPb 32 PHI BINS, pT bins 0.8-10.0, 0.2-4.0, 0.2-0.8 GeV/c
    //    const char *taskdirname = "task_2016_08_02_FOR_LXPLUS_pPb_Data_LHC13c_AOD154_eW2_phi16_pt3_BIT768_V0M_c0_5_try2";
    //    const char *taskdirname = "task_2016_08_02_FOR_LXPLUS_pPb_Data_LHC13c_AOD154_eW1_phi29_pt2_BIT768_CL1_c0_6";
    
    //8.08.2016: DATA PbPb 32 PHI BINS, eW2_phi32_pt4_HIGHER_LOW_BOUNDARY_BIT768_c0_2% (V0M)
    //    const char *taskdirname = "task_2016_08_05_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt4_HIGHER_LOW_BOUNDARY_BIT768_c0_2";
    
    //9.08.2016: DATA PbPb 8 PHI BINS, c0_2% (V0M)
    //    const char *taskdirname = "task_2016_08_08_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW2_phi8_pt5_BIT768_c0_2";
    
    //3.08.2016: EXTRA eW (-0.8, 0)...(0, 0.8) AND CHARGE COMBINATIONS!
    //charges for phi bins 32
    //    const char *taskdirname = "task_2016_08_03_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW3_phi32_pt2_BIT768_c0_2_EXTRA_eW_AND_CHARGES/MFplus";
    //charges for phi bins 36
    //    const char *taskdirname = "task_2016_08_04_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW1_phi36_pt2_BIT768_c0_2_CHARGES/MFminus";
    
    //5.08: try phi in AMPT reco LHC12a11a (32 phi bins!!!)
    //    const char *taskdirname = "task_2016_08_04_FOR_LXPLUS_PbPb_AMPT_LHC12a11a_RECO_AOD157_bit768_c0_5";
    
    //3.08.2016: 36 phi bins, SHIFTED_PHI_BY_10_DEG
    //    const char *taskdirname = "task_2016_08_03_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW3_phi36_pt1_BIT768_c0_2_SHIFTED_PHI_BY_10_DEG";
    
    // !!! 1.08.2016: MF+/- test for PHI 32 BINS!
    //    const char *taskdirname = "task_2016_07_02_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt6_BIT768_PHI_PT_MANY_WINS/MFplus";
    //    const char *taskdirname = "task_2016_07_02_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt6_BIT768_PHI_PT_MANY_WINS/MFminus";
    // !!! 1.08.2016: analyse 9 largest runs
    //    const char *taskdirname = "task_2016_07_02_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt6_BIT768_PHI_PT_MANY_WINS/tmp_9_largest_runs";
    
    //20.08.2016: 32 phi bins, ETA WINS
    //    const char *taskdirname = "task_2016_08_18_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_eW4_phi32_pt1_BIT768_c0_2_MORE_ETA_WINS";
    
    
    //AMPT_LHC13f3c_ON_ON_eW2_phi32_pt5_BIT768_PHI_PT_MANY_WINS
    //    const char *taskdirname = "task_2016_06_27_PbPb_AMPT_LHC13f3c_ON_ON_eW2_phi32_pt5_BIT768_PHI_PT_MANY_WINS_TRY2";
    
    
    //SYSTEMATIC CHECKS: PT RES, DCA(19.06.2016):
    //    const char *taskdirname = "task_2016_06_17_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_BIT768_SYST_PT_RES_DCA";
    //    const char *taskdirname = "task_2016_06_20_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_BIT768_SYST_DCA_MORE";
    
    //HIJING FOR CORRECTIONS: task_2016_05_01_PbPb_MCAOD_LHC11a10a_bis_bit768_AOD162_kine_PLUS_reco_in_one_tree
    //    const char *taskdirname = "task_2016_05_01_PbPb_MCAOD_LHC11a10a_bis_bit768_AOD162_kine_PLUS_reco_in_one_tree";
//    const char *taskdirname = "task_2016_06_01_PbPb_HIJING_LHC11a10a_bis_bit768_AOD162_WEAK_SEC_in_tree_FINAL"; //RESULTS_SEPARATE_FILES";
    //    const char *taskdirname = "task_2016_06_01_PbPb_HIJING_LHC11a10a_bis_bit128_AOD162_WEAK_SEC_in_tree_FINAL"; //RESULTS_SEPARATE_FILES";
    
    // 08/06 AMPT LHC13f3c: task KINE + task with INEFF BY HAND
    // TMP: const char *taskdirname = "task_2016_06_01_PbPb_AMPT_LHC13f3c_StrMelt_ON_rescatON_eW8_phi1_pt02_20_FB_TREE_KINE_INEFF";
    //    const char *taskdirname = "task_2016_06_12_PbPb_AMPT_LHC13f3c_ON_ON_eW8_phi1_pt02_20_KINE_INEFF_cBinsEmul";
    //        const char *taskdirname = "task_2016_06_12_PbPb_AMPT_LHC13f3c_ON_ON_eW8_phi1_pt02_20_KINE_INEFF_cBinsEmul_NEW";
    
    
    //LHC10h plus, minus INEFF STUDY:
    //    const char *taskdirname  = "task_2016_03_10_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND_HIST3D/MFplus";
    //    const char *taskdirname2 = "task_2016_03_10_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND_HIST3D/MFminus";
    
    //LHC NEW HIJING RECO KINE:
    //    const char *taskdirname  = "task_2016_03_16_PbPb_MCAOD_LHC11a10a_bis_AOD162_Efficiency_kine_vs_reco_KINE_AND_RECO";
    
    
    
    
    
    // ######## SEPTEMBER 2016 ANALYSIS:
    
    // 7.09.2016: PbPb PHI_WINS 32 - HIJING 0-10% (V0M) - LHC11a10a_bis_bit768_AOD162_eW2_phi32_pt2
    //    const char *taskdirname = "task_2016_09_06_FOR_LXPLUS_PbPb_HIJING_LHC11a10a_bis_bit768_AOD162_eW2_phi32_pt2_BIT768_c0_10";
    //    const char *taskdirname = "task_2016_09_07_FOR_LXPLUS_PbPb_HIJING_LHC11a10a_bis_AOD162_eW2_phi32_pt2_c0_10_MC_KINE";
    //    const char *taskdirname = "task_2016_09_28_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_SYST_DCA_TPCclust";
    
    // HIJING DCA, TPC clust
    //    const char *taskdirname = "task_2016_09_28_FOR_LXPLUS_PbPb_HIJING_LHC11a10a_bis_bit768_AOD162_8eW_phi1_pt02_20_SYST_DCA_TPCclust";
    //    const char *taskdirname = "task_2016_10_12_FOR_LXPLUS_PbPb_HIJING_LHC11a10a_bis_bit768_AOD162_8eW_phi1_pt02_20_SYST_DCA_TPCclust";
    
    
    // 12.09.2016: PbPb PHI_WINS 32 - RUN2 Data High Intensity Runs
    //    const char *taskdirname = "task_2016_09_12_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_eW2_phi32_pt2_BIT32_c0_10_try2";
    //    const char *taskdirname = "task_2016_09_12_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_eW2_phi32_pt2_BIT128_c0_2_try2";
    //    const char *taskdirname = "task_2016_09_24_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_eW2_phi32_pt1_BIT128_c0_2_MANY_RUNS";
    
    // 13.09.2016: PbPb FULL PHI - RUN2 Data High Intensity Runs
    //    const char *taskdirname = "task_2016_09_13_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_8eW_phi1_pt02_20_BIT32_ETA_WINS";
    //    const char *taskdirname = "task_2016_09_14_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_8eW_phi1_pt02_08_and_08_10_BIT32";
    //    const char *taskdirname = "task_2016_09_21_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_08_and_08_10_BIT768";
    //    const char *taskdirname = "task_2016_09_24_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_8eW_phi1_pt02_20_BIT32_ETA_WINS_MANY_RUNS";
    //    const char *taskdirname = "task_2016_09_30_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_8eW_phi1_pt02_20_BIT32_ETA_WINS_NewAliPhys_v";
    //    const char *taskdirname = "task_2016_09_29_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_8eW_pt4_NewPtWins_BIT128_ETA_WINS_MANY_RUNS";
    
    // 27.09.2016: EPOS phi32 all centralities (studied 0-5% by multTPC)
    //    const char *taskdirname = "task_2016_09_25_FOR_LXPLUS_PbPb_EPOS_LHC16d2_502TeV_run245064_LHC15o_kine_1Mln_PHI_WINS32_allEvents";
    //    const char *taskdirname = "task_2016_09_25_FOR_LXPLUS_PbPb_EPOS_LHC16d2_502TeV_run245064_LHC15o_recoESD_1Mln_PHI_WINS32_allEvents";
    
    // 28.09.2016: EPOS FULL PHI all centralities
    //    const char *taskdirname = "task_2016_09_26_FOR_LXPLUS_PbPb_EPOS_LHC16d2_502TeV_run245064_LHC15o_kine_1Mln_ETA_WINS";
    //    const char *taskdirname = "task_2016_09_26_FOR_LXPLUS_PbPb_EPOS_LHC16d2_502TeV_run245064_LHC15o_recoESD_1Mln_ETA_WINS";
    
    
    // ######## NOVEMBER 2016 ANALYSIS:
//    const char *taskdirname = "task_2016_10_29_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_1eW_phi1_pt22bins_FOR_NN_CORRECTION";
//    const char *taskdirname = "task_2016_10_31_FOR_LXPLUS_PbPb_HIJING_LHC11a10a_bis_bit768_AOD162_1eW_phi1_pt9bins_FOR_NN_CORRECTION";
//    const char *taskdirname = "task_2016_11_02_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_1eW_phi1_pt9bins_FOR_NN_CORRECTION";
    // bit128 and 32!!! (to cross-check the results):
//    const char *taskdirname = "task_2017_01_19_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_1eW_phi1_pt9bins_FOR_NN_CORRECTION_bit128";
//    const char *taskdirname = "task_2017_01_19_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_1eW_phi1_pt9bins_FOR_NN_CORRECTION_bit32";

    // ######## NOVEMBER 2016 ANALYSIS WITH PID:
    // BEFORE ELLIPTIC sigma cut, fPIDnSigmaCut=2, WITH lMaxProb < 0.8,  ### filterBit768!  // WAS DONE WITH lMostProbablePIDdirty FOR FB_TREE!!! (~MISTAKE)
//    const char *taskdirname = "task_2016_11_08_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_1eW_phi1_pt4bins_PID_FIRST_TRY";
//    const char *taskdirname = "task_2016_11_08_FOR_LXPLUS_PbPb_Data_LHC10h_AOD160_1eW_02_08_phi1_pt2bins_PID_SECOND_TRY";

    // ######## DECEMBER 2016 ANALYSIS:
//    const char *taskdirname = "task_2016_12_01_PbPb_Data_LHC10h_AOD160_eW08_02_phi1_pt2_BIT32_CHECK_SUGGESTED_BY_MW";
//    const char *taskdirname = "task_2016_12_02_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_BIT32_CHECK_SUGGESTED_BY_MW";
//    const char *taskdirname = "task_2016_12_03_PbPb_HIJING_LHC11a10a_bis_bit768_AOD162_8eW_phi1_pt02_20_BIT32";
//5.02 TeV:
//    const char *taskdirname = "task_2016_12_04_HIJING_Run2_LHC16g1_AOD_8eW_phi1_pt4_BIT32";
//        const char *taskdirname = "task_2016_12_04_Data_LHC15o_HighInt_AOD_8eW_phi1_pt4_BIT32";
//    const char *taskdirname = "task_2016_12_06_Data_LHC15o_HighInt_AOD_1eW_phi1_pt4_BIT32";
//    const char *taskdirname = "task_2017_01_20_FOR_LXPLUS_PbPb_Data_LHC15o_HighInt_AOD_1eW_phi1_BIT768";
//    const char *taskdirname = "task_2017_01_22_PbPb_Data_LHC15o_LowInt_AOD_1eW_phi1_BIT32";
//    const char *taskdirname = "task_2017_01_24_Data_LHC15o_HighInt_AOD_1eW_phi1_pt1_BIT32_PHI_GAP_BY_HAND_17_25_try2";
//    const char *taskdirname = "task_2017_01_22_HIJING_LHC16g1_AOD_1eW_phi1_pt1_BIT32_PHI_GAP_BY_HAND_17_25";
//    const char *taskdirname = "task_2017_01_25_Data_LHC15o_HighInt_AOD_1eW_phi1_pt1_BIT768_DCA_TPC_cuts";
//    const char *taskdirname = "task_2017_01_26_HIJING_LHC16g1_AOD_1eW_phi1_pt1_BIT768_DCA_TPC_cuts";

    // MC with injected flow? (Run-1)
//    const char *taskdirname = "task_2017_01_27_TRY_MC_WITH_FLOW_LHC14b8a_AOD_1eW_phi1_pt1_BIT768";

    // March 1: test more evSel criteria (pileup, etc.)
//    const char *taskdirname = "task_2017_03_01_PbPb_Data_LHC15o_HighInt_AOD_1eW_phi1_pt1_BIT768_NewEvSelCuts_try2";
//    const char *taskdirname = "task_2017_03_01_PbPb_Data_LHC15o_HighInt_AOD_1eW_phi1_pt2_BIT768_NewEvSelCuts";
//    const char *taskdirname = "task_2017_03_25_PbPb_Data_LHC15o_HighInt_AOD_1eW_phi1_pt1_BIT768_PhysSelWithPileUpCuts";
//    const char *taskdirname = "task_2017_03_29_PbPb_Data_LHC15o_HighInt_AOD_eW2_phi32_pt08_50_BIT768_c0_2_NewEvSelCuts";
//    const char *taskdirname = "task_2017_04_01_PbPb_Data_LHC15o_HighInt_AOD_eW1_phi32_pt10_15_20_50_BIT768_c0_2_NewEvSelCuts";
    const char *taskdirname = "task_2017_04_02_PbPb_Data_LHC15o_HighInt_AOD_eW1_phi32_pt02_05_17_BIT768_c0_2_NewEvSelCuts";

    int nFiles = 0; // = 24+21; //30; //24+22; //3;//10;
    TString *runFullPathList = new TString[400];
    int *runListNumbers = new int[400];
    
    get_run_list( nFiles, runFullPathList, runListNumbers, dirBase, taskdirname, "AnalysisResults" );
    //    get_run_list( nFiles, runFullPathList, runListNumbers, dirBase, taskdirname2, "AnalysisResults" );
    
    if (0) //combined AMPT to increase statistics
    {
        const char *taskdirname3 = "task_2016_05_09_PbPb_AMPT_LHC12a11c_c10_20_bit768_AOD157_WEAK_SEC_in_tree";
        const char *taskdirname4 = "task_2016_05_09_PbPb_AMPT_LHC12a11d_c20_30_bit768_AOD157_WEAK_SEC_in_tree";
        const char *taskdirname5 = "task_2016_05_09_PbPb_AMPT_LHC12a11e_c30_40_bit768_AOD157_WEAK_SEC_in_tree";
        const char *taskdirname6 = "task_2016_05_09_PbPb_AMPT_LHC12a11f_c40_50_bit768_AOD157_WEAK_SEC_in_tree";
        get_run_list( nFiles, runFullPathList, runListNumbers, dirBase, taskdirname3, "AnalysisResults" );
        get_run_list( nFiles, runFullPathList, runListNumbers, dirBase, taskdirname4, "AnalysisResults" );
        get_run_list( nFiles, runFullPathList, runListNumbers, dirBase, taskdirname5, "AnalysisResults" );
        get_run_list( nFiles, runFullPathList, runListNumbers, dirBase, taskdirname6, "AnalysisResults" );
    }
    
    cout << ">>> nFiles = " << nFiles << endl;
    TFile *myFiles[400];
    
//        fileIdByHand = 0; // !!!!
    //    return;
    
    //    runListNumbers[0] = 138442;
    //    runFullPathList[0] = Form("/Users/macbook/alice/aliceAnalysis/results/task_2016_07_02_PbPb_Data_LHC10h_AOD160_eW2_phi32_pt6_BIT768_PHI_PT_MANY_WINS/AnalysisResults_000138442.root");
    
    if ( fileIdByHand >= 0 )
    {
        myFiles[0] = new TFile( runFullPathList[fileIdByHand] );
        runListNumbers[0] = runListNumbers[fileIdByHand];
        nFiles = 1;
    }
    else //usual case:
    {
        for (Int_t i=0; i<nFiles; i++)
            myFiles[i] = new TFile( runFullPathList[i] );
    }
    
    
    
    //        myFiles[i] = new TFile( Form("/Users/macbook/alice/aliceAnalysis/results/task_2016_03_10_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND_HIST3D/MFplus/AnalysisResults.000%d.root", runFullPathList[i] ) );
    //    myFiles[i] = new TFile( Form("/Users/macbook/alice/aliceAnalysis/results/task_2016_03_10_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND_HIST3D/MFminus/AnalysisResults_000%d.root", runFullPathList[i] ) );
    
    //    myFiles[0] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFplus_all/MergedOutput.root" );
    //    myFiles[1] = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_26_PbPb_Data_LHC10h_AOD160_8eW_phi1_pt02_20_FB_TREE_MFminus_all/MergedOutput.root" );
    
    
    
    // INEFF BY HAND
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_25_PbPb_Data_LHC10h_AOD160_3etaWins_phi1_pt02_20_FB_TREE_INEFF_BY_HAND/MergedOutput.root" );
    
    
    
    //LHC10h
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_2_76TeV_LHC10h_AOD0160_3etaWins_phi1_pt02_20_FB_TREE_blocks123/block1/AnalysisResults.139465.root" );
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_09_PbPb_Data_2_76TeV_LHC10h_AOD0160_3etaWins_phi1_pt02_20_FB_TREE_blocks123/MergedOutput.root" );
    
    //LHC11h
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoMinus/MergedOutput.root" );
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoPlus/MergedOutput.root" );
    
    //LHC11h Femto Plus Minus merged
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_23_PbPb_Data_2_76TeV_LHC11h_AOD145_3etaWins_phi1_pt02_20_FB_TREE_FemtoPlus/MergedResults/MergedOutput_Femto_Plus_Minus.root" );
    
    //LHC15o MF ++ --
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_22_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldMM/MergedOutput.root" );
    //        TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_21_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldPP/MergedOutput.root" );
    
    //LHC15o MF merged
    //    TFile *myFile = new TFile( "/Users/macbook/alice/aliceAnalysis/results/task_2016_02_22_PbPb_Data_5_02TeV_LHC15o_AOD_3etaWins_phi1_pt02_20_FB_TREE_fieldMM/MergedResults/MergedOutput_PP_and_MM.root" );
    
    
    
    //TMP!!! 13.06.2016:
    int nEvAMPTbyHand[10] =
    {
        //        230000,
        //        230000,
        //        460000,
        //        460000,
        //        460000,
        //        460000,
        //        460000,
        //        460000,
        //        460000,
        //        460000,
        395000,
        395000,
        790000,
        790000,
        790000,
        790000,
        790000,
        790000,
        790000,
        790000,
    };
    
    const int nTrees = nFiles;
    TTree *trees[nTrees];
    int nEvInTrees[nTrees];
    
    //get trees from files
    int nEvents = 0;
    for ( int i = 0; i < nTrees; i++ )
    {
        if ( !myFiles[i] )
        {
            cout << "No input file! i=" << i << endl;
            return;
        }
        myFiles[i]->ls();
        
        //        if (!isAMPT_kine)
        trees[i] = (TTree*) getTreeFromFile( myFiles[i], LIST_ID_BY_HAND );
        //        else
        //            trees[i] = (TTree*) myFiles[i]->Get( "t1" );
        
        
        nEvInTrees[i] = trees[i]->GetEntries();
        
        if (0)
            nEvInTrees[i] = nEvAMPTbyHand[i];
        cout << "nEvents in tree " << i << ": " << nEvInTrees[i] << endl;
        
        nEvents += nEvInTrees[i];
    }
    
    if ( nEventsByHand >= 0 )
        nEvents = nEventsByHand;
    
    //        nEvents = 1000000;
    //    nEvents = trees[0]->GetEntries() + trees[1]->GetEntries();
    //    nEvents = trees[0]->GetEntries();
    //    cout <<"nEvents = " << nEvents << endl;
    
    //    return;
    
    //if take 1 run: apply cut on number of events!
    if(0)if ( fileIdByHand >=0 )
    {
        if ( nEvents < 1e5 || nEvents > 6e5)
        {
            myFiles[0]->Close();
            return;
        }
    }
    
    // ########## Branches for event-info
    Bool_t fEvSel_isPileupSPD;
    Bool_t fEvSel_isPileupFromMV ;
    Bool_t fEvSel_vtxSPD_isFromVertexer3D         ;
    Bool_t fEvSel_vtxSPD_isFromVertexerZ          ;
    Bool_t fEvSel_zRes_above_fMaxResol            ;
    Bool_t fEvSel_diff_vZ_SPD_tracks_above_05     ;
    Bool_t fEvSel_comb_cut_on_abs_and_nsigma_dist ;
    Bool_t fEvSel_too_many_TPC_clusters           ;

    UShort_t fEvSel_multTrk_bit32           ;
    UShort_t fEvSel_multTrk_bit32_TOF           ;

    Float_t fDiff_vZ_SPD_tracks;
    Float_t fDiff_vZ_TPC_tracks;

    Float_t brZEMvsZDC = 0;
    Float_t brV0M = 0;
    Float_t brCL1 = 0;
    Float_t br_vertexZ = 0;
    
    
    //    cout << "TEST1" << endl;
    
    // ########## Number of eta-phi wins
    //    const int nEtaBr = 3;
    //        const int nEtaBr = 1;
    const int nEtaBr = 1;//8;//1; //8;//1;//8; // ANCHOR ANALYSIS
    //    const int nEtaBr = 1; // FOR EXAMPLE, FOR PHI WINS
    const int nPhiWins = 32;//1;//32;//8;//32;//8;//32;//29;//32;//20;//23;//32;//64;//32; //20;//1;//32;//8;
    const int nEtaIdSPEC = 0;//1;//1;//1; //to spec set exact match between br array and branch name in tree
    bool FLAG_PERCOLATING_WINS = false;//true;
    bool FLAG_ONE_WIDE_ETA_WIN = false;
    
//    const int nPtBr = 22;//1;
//    const int nPtBr = 9;
    const int nPtBr = 1;

    BranchFB br[nEtaBr][nPhiWins][nPtBr];
    
    // for unfolding:
    bool isUseSpecRecoWithKineForResponse = 0;
    BranchFB br_spec_reco[nEtaBr][nPhiWins][nPtBr];
    
    
    //    const double nbins = 30;//25;//24;
    //    const double nmin = -0.5;
    //    const double nmax = 300-0.5;
    
    //    const double nbins = 25;//24;
    //    const double nmax = 100;
    //    TH2D histFB_NN_truth("histFB_NN_truth","histFB_NN_truth", nbins, -0.5, nmax-0.5, nbins, -0.5, nmax-0.5 );
    //    TH2D histFB_NN_rec("histFB_NN_rec","histFB_NN_rec", nbins, -0.5, nmax-0.5, nbins, -0.5, nmax-0.5 );
    //    RooUnfoldResponse response_NN;
    //    response_NN.Setup( &histFB_NN_rec, &histFB_NN_truth );
    //    response_NN.SetName( "unfoldResponse_NN" );
    
    //    const double nbins_ptpt = 30;//24;
    //    const double min_ptpt = 0.5;
    //    const double max_ptpt = 0.65;
    
    const double nbins_ptpt = 32;
    const double min_ptpt = 0.4;
    const double max_ptpt = 0.8;
    
    
    // for UNFOLDING:
    //    tmp: // TH2D *histFB_PtPt_truth = new TH2D( "histFB_PtPt_truth","histFB_PtPt_truth", 40, 0.50, 0.70,     40, 0.50, 0.70 );
    //    tmp: // TH2D *histFB_PtPt_rec = new TH2D( "histFB_PtPt_rec","histFB_PtPt_rec", 30, 0.50, 0.65,   30, 0.50, 0.65 );
    // 2D-unfolding to recover FB correlations!
    //    TH2D *histFB_PtPt_truth = new TH2D( "histFB_PtPt_truth","histFB_PtPt_truth", nbins_ptpt, min_ptpt, max_ptpt, nbins_ptpt, min_ptpt, max_ptpt );
    //    TH2D *histFB_PtPt_rec = new TH2D( "histFB_PtPt_rec","histFB_PtPt_rec", nbins_ptpt, min_ptpt, max_ptpt, nbins_ptpt, min_ptpt, max_ptpt );
    // 1D-unfolding to recover Nch_truth!
    //    TH1D *histFB_PtPt_truth = new TH1D( "histFB_PtPt_truth","histFB_PtPt_truth", nbins_ptpt, min_ptpt, max_ptpt );
    //    TH1D *histFB_PtPt_rec   = new TH1D( "histFB_PtPt_rec","histFB_PtPt_rec", nbins_ptpt, min_ptpt, max_ptpt );
    //    TH1D *histFB_PtPt_truth = new TH1D( "histFB_NN_truth",  "histFB_NN_truth"   , 1001, -0.5, 1000.5 );
    //    TH1D *histFB_PtPt_rec   = new TH1D( "histFB_NN_rec",    "histFB_NN_rec"     , 1001, -0.5, 1000.5 );
    TH1D *histFB_PtPt_truth = new TH1D( "histFB_NN_truth",  "histFB_NN_truth"   , 11, -0.5, 10.5 ); //tmp to reduce memory usage
    TH1D *histFB_PtPt_rec   = new TH1D( "histFB_NN_rec",    "histFB_NN_rec"     , 11, -0.5, 10.5 ); //tmp to reduce memory usage

    //    RooUnfoldResponse response_PtPt;
    if (isUseSpecRecoWithKineForResponse)
    {
        // !!! UNCOMMENT IF NEED UNFOLD->   response_PtPt->Setup( histFB_PtPt_rec, histFB_PtPt_truth );
        // !!! UNCOMMENT IF NEED UNFOLD->   response_PtPt->SetName( "unfoldResponse_PtPt" );

        // unfolding in centarlity bins:
        // !!! UNCOMMENT IF NEED UNFOLD
        /*
        for (Int_t cBin = 0; cBin < 8; cBin++ ) // for each cBin!
        {
            response_PtPt[cBin].Setup( histFB_PtPt_rec, histFB_PtPt_truth );
            response_PtPt[cBin].SetName( Form( "unfoldResponse_PtPt_cBin%d", cBin ) );
        }
        */
    }
    
    
    //set branch addresses:
    for ( int iTree = 0; iTree < nTrees; iTree++ )
    {
        TTree *tr = trees[iTree];
        
        // begin of EvSel cuts
        tr->SetBranchAddress( "isPileupSPD",                     &fEvSel_isPileupSPD                    );
        tr->SetBranchAddress( "isPileupFromMV",                  &fEvSel_isPileupFromMV                 );

        tr->SetBranchAddress( "vtxSPD_isFromVertexer3D",         &fEvSel_vtxSPD_isFromVertexer3D         );
        tr->SetBranchAddress( "vtxSPD_isFromVertexerZ",          &fEvSel_vtxSPD_isFromVertexerZ          );
        tr->SetBranchAddress( "zRes_above_fMaxResol",            &fEvSel_zRes_above_fMaxResol            );
        tr->SetBranchAddress( "diff_vZ_SPD_tracks_above_05",     &fEvSel_diff_vZ_SPD_tracks_above_05     );
        tr->SetBranchAddress( "comb_cut_on_abs_and_nsigma_dist", &fEvSel_comb_cut_on_abs_and_nsigma_dist );
        tr->SetBranchAddress( "too_many_TPC_clusters",           &fEvSel_too_many_TPC_clusters           );

        tr->SetBranchAddress( "multTrk_bit32",      &fEvSel_multTrk_bit32 );
        tr->SetBranchAddress( "multTrk_bit32_TOF",  &fEvSel_multTrk_bit32_TOF );

        tr->SetBranchAddress( "diff_vZ_SPD_tracks", &fDiff_vZ_SPD_tracks );
        tr->SetBranchAddress( "diff_vZ_TPC_tracks", &fDiff_vZ_TPC_tracks );
        // end of EvSel cuts

        tr->SetBranchAddress( "vertexZ", &br_vertexZ );
        tr->SetBranchAddress( "centr_ZEMvsZDC", &brZEMvsZDC );
        tr->SetBranchAddress( "centr_CL1", &brCL1 );
        
        if ( !isAnalysing502 )
            tr->SetBranchAddress( "centr_V0M", &brV0M );
        else
            tr->SetBranchAddress( "centrV0M_NEW_MULT_SEL", &brV0M );
        
        // binding of wins branches:
        for ( int ptBr = 0; ptBr < nPtBr; ptBr++ )
        {
            ptW = 2;//1;//ptBr;//0; // !!!!!!!! be careful
//            ptW = ptBr; // !!!!!!!! be careful
            int ptW_F = ptW;
            int ptW_B = ptW;

            for ( int etaBr = 0; etaBr < nEtaBr; etaBr++ )
                for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                {
                    TString brNamePostfixF = Form("eta_%d_phi%d_pt_%d"
                                                  , etaBr+nEtaIdSPEC, phiW, ptW_F );
                    TString brNamePostfixB = Form("eta_%d_phi%d_pt_%d"
                                                  , etaBr+nEtaIdSPEC, phiW, ptW_B );

                    if( WHICH_BRANCH_BY_HAND == 0 )
                    {
                        tr->SetBranchAddress( Form("nF_%s", brNamePostfixF.Data() ),  &br[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_%s", brNamePostfixB.Data() ),  &br[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_%s", brNamePostfixF.Data() ), &br[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_%s", brNamePostfixB.Data() ), &br[etaBr][phiW][ptBr].PtB );
                    }
                    else if( WHICH_BRANCH_BY_HAND == 1 ) // TMP!!!! for manual operation!!!
                    {
                        tr->SetBranchAddress( Form("nF_spec_reco_%s", brNamePostfixF.Data() ),  &br[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_spec_reco_%s", brNamePostfixB.Data() ),  &br[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_spec_reco_%s", brNamePostfixF.Data() ), &br[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_spec_reco_%s", brNamePostfixB.Data() ), &br[etaBr][phiW][ptBr].PtB );
                    }
                    else if( WHICH_BRANCH_BY_HAND == 2 ) // TMP!!!! for manual operation!!!
                    {
                        tr->SetBranchAddress( Form("nF_physPrimary_%s", brNamePostfixF.Data() ),  &br[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_physPrimary_%s", brNamePostfixB.Data() ),  &br[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_physPrimary_%s", brNamePostfixF.Data() ), &br[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_physPrimary_%s", brNamePostfixB.Data() ), &br[etaBr][phiW][ptBr].PtB );
                    }
                    else if( WHICH_BRANCH_BY_HAND == 3 ) // TMP!!!! for manual operation!!!
                    {
                        tr->SetBranchAddress( Form("nF_nonPhysPrimary_%s", brNamePostfixF.Data() ),  &br[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_nonPhysPrimary_%s", brNamePostfixB.Data() ),  &br[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_nonPhysPrimary_%s", brNamePostfixF.Data() ), &br[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_nonPhysPrimary_%s", brNamePostfixB.Data() ), &br[etaBr][phiW][ptBr].PtB );
                    }
                    else if( WHICH_BRANCH_BY_HAND == 4 ) // TMP!!!! for manual operation!!!
                    {
                        tr->SetBranchAddress( Form("nF_fromWeakDecays_%s", brNamePostfixF.Data() ),  &br[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_fromWeakDecays_%s", brNamePostfixB.Data() ),  &br[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_fromWeakDecays_%s", brNamePostfixF.Data() ), &br[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_fromWeakDecays_%s", brNamePostfixB.Data() ), &br[etaBr][phiW][ptBr].PtB );
                    }
                    else if( WHICH_BRANCH_BY_HAND == 5 ) // TMP!!!! for manual operation!!!
                    {
                        tr->SetBranchAddress( Form("nF_fromMaterial_%s", brNamePostfixF.Data() ),  &br[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_fromMaterial_%s", brNamePostfixB.Data() ),  &br[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_fromMaterial_%s", brNamePostfixF.Data() ), &br[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_fromMaterial_%s", brNamePostfixB.Data() ), &br[etaBr][phiW][ptBr].PtB );
                    }



                    if ( isUseSpecRecoWithKineForResponse )
                    {
                        tr->SetBranchAddress( Form("nF_spec_reco_%s", brNamePostfixF.Data() ),  &br_spec_reco[etaBr][phiW][ptBr].nF );
                        tr->SetBranchAddress( Form("nB_spec_reco_%s", brNamePostfixB.Data() ),  &br_spec_reco[etaBr][phiW][ptBr].nB );
                        tr->SetBranchAddress( Form("PtF_spec_reco_%s", brNamePostfixF.Data() ), &br_spec_reco[etaBr][phiW][ptBr].PtF );
                        tr->SetBranchAddress( Form("PtB_spec_reco_%s", brNamePostfixB.Data() ), &br_spec_reco[etaBr][phiW][ptBr].PtB );
                        //                    tr->SetBranchAddress( Form("nF_physPrimary_%s", brNamePostfixF.Data() ),  &br_spec_reco[etaBr][phiW].nF );
                        //                    tr->SetBranchAddress( Form("nB_physPrimary_%s", brNamePostfixB.Data() ),  &br_spec_reco[etaBr][phiW].nB );
                        //                    tr->SetBranchAddress( Form("PtF_physPrimary_%s", brNamePostfixF.Data() ), &br_spec_reco[etaBr][phiW].PtF );
                        //                    tr->SetBranchAddress( Form("PtB_physPrimary_%s", brNamePostfixB.Data() ), &br_spec_reco[etaBr][phiW].PtB );
                    }
                }
        } // end of pt branches
    }
    
    // ##### QA pre-loop (for mult binning, etc.)
    
    TH1D *hist1D_QA_percentilesEstimator = new TH1D( "hist1D_QA_percentilesEstimator", "hist1D_QA_percentilesEstimator;percentile;entries", 302, -1, 301);
    TH1D *hist1D_QA_multALL = new TH1D( "hist1D_QA_multALL", "hist1D_QA_multALL;mult;entries", 3001, -0.5, 3000.5);
    
    TH2D *hist2D_ESTIMATOR_VS_multTPC = new TH2D( "hist2D_ESTIMATOR_VS_multTPC", "hist2D_ESTIMATOR_VS_multTPC;estimator;mult in TPC", 4080, -2, 100, 301, -0.5, 3000.5);
    
    TH2D *hist2D_QA_multBit32_vs_multWithTOFtiming = new TH2D( "hist2D_QA_multBit32_vs_multWithTOFtiming", "Pileup from out-of-bunch using TOF;mult (bit32);mult (bit32+TOF)",250,0,5000,250,0,2500);

    //QA mult in eta win IN EACH TREE (=run-by-run):
    TH1D *hist1D_multInWin[nTrees];
    TH1D *hist1D_vertexZ[nTrees];
    
    //23.03.2016: new more useful histos: mult distr in each win run-by-run:
    //    TH1D *hist1D_multDistrInWin[nTrees][nEtaWins];
    //    TH1D *hist1D_avPtDistrInWin[nTrees][nEtaWins];
    for ( int i = 0; i < nTrees; i++ )
    {
        TString strMultDistr_name = Form("hist1D_multDistr_run_%d", runListNumbers[i] );//, etaW );
        //            cout << "strMultDistr_name=" << strMultDistr_name << endl;
        hist1D_multInWin[i] = new TH1D( strMultDistr_name, ";etaWin;n tracks", 2*nEtaBr, -0.5, 2*nEtaBr-0.5 );
        
        TString strVertexZ_name = Form("hist1D_vertZdistr_run_%d", runListNumbers[i] );//, etaW );
        //            cout << "strVertexZ_name=" << strVertexZ_name << endl;
        hist1D_vertexZ[i] = new TH1D( strVertexZ_name, ";vertex Z, cm;n events", 300, 15, 15 );
        
        //23.03.2016: new more useful histos: mult distr in each win run-by-run:
        //        for ( int etaW = 0; etaW < nEtaBr; etaW++ )
        //        {
        //            strMultDistr_name = Form("hist1D_multDistr_run_%d_inWin%d", runListNumbers[i], etaW );
        //            hist1D_multDistrInWin[i][etaW] = new TH1D( strMultDistr_name, ";n tracks;n events", 400, -0.5, 400-0.5 );
        
        //            TString strAvPtDistr_name = Form("hist1D_avPtDistr_run_%d_inWin%d", runListNumbers[i], etaW );
        //            hist1D_avPtDistrInWin[i][etaW] = new TH1D( strAvPtDistr_name, ";#LTp_{T}#GT;n events", 400, 0, 2 );
        //        }
        
    }
    
    
    // ##### pre-loop over events for mult bins
    int nAccepted_PRE_LOOP = 0;
    int nRejectedByPileup_PRE_LOOP = 0;
    int treeId = 0;
    int sumEvPrevTrees = 0;//nEvInTrees[0];
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "pre-loop: getting " << (int)i << endl;
        
        
        if ( i >= sumEvPrevTrees + nEvInTrees[treeId] ) // it's time to go to next tree...
        {
            sumEvPrevTrees += nEvInTrees[treeId]; // sum of events in prev trees
            treeId++;
            cout << "check treeId = " << treeId << " brV0M = " << brV0M << endl;
        }
        trees[treeId]->GetEntry( i - sumEvPrevTrees );
        
        // ### event selection
        float cEstimator;
        bool isEventSelected = checkEventSelection(
                    cEstimator

                    , fEvSel_isPileupSPD
                    , fEvSel_isPileupFromMV
                    , fEvSel_vtxSPD_isFromVertexer3D
                    , fEvSel_vtxSPD_isFromVertexerZ
                    , fEvSel_zRes_above_fMaxResol
                    , fEvSel_diff_vZ_SPD_tracks_above_05
                    , fEvSel_comb_cut_on_abs_and_nsigma_dist
                    , fEvSel_too_many_TPC_clusters

                    , fEvSel_multTrk_bit32
                    , fEvSel_multTrk_bit32_TOF


                    , brZEMvsZDC
                    , brV0M
                    , brCL1
                    , br_vertexZ
                    , nRejectedByPileup_PRE_LOOP
                    );
        
        if ( !isEventSelected )
            continue;
        
        //                cout << ">>> br_vertexZ=" << br_vertexZ << endl;
        
        hist1D_QA_percentilesEstimator->Fill(cEstimator);
        
        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int ptBr = 0; ptBr < nPtBr; ptBr++ )
            for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            {
                if ( FLAG_ONE_WIDE_ETA_WIN )
                    multTPC += br[0][phiW][ptBr].nF; // + br[etaW][phiW].nB;
                else if ( !FLAG_PERCOLATING_WINS ) // usual case!
                    for ( int etaW = 0; etaW < nEtaBr; etaW++ )
                        multTPC += br[etaW][phiW][ptBr].nF + br[etaW][phiW][ptBr].nB;
                else
                    multTPC += (br[0][phiW][ptBr].nF + br[2][phiW][ptBr].nF) + (br[0][phiW][ptBr].nB + br[2][phiW][ptBr].nB);

                //fill QA mult in eta wins run-by-run
                for ( int etaW = 0; etaW < nEtaBr; etaW++ )
                {
                    //B
                    double binContent = hist1D_multInWin[treeId]->GetBinContent(etaW+1);
                    hist1D_multInWin[treeId]->SetBinContent(etaW+1, binContent + br[etaW][phiW][ptBr].nB);
                    //new:
                    //                hist1D_multDistrInWin[treeId][etaW]->Fill( br[etaW][phiW][ptBr].nB);
                    //                if ( br[etaW][phiW][ptBr].nB > 0 )
                    //                    hist1D_avPtDistrInWin[treeId][etaW]->Fill( br[etaW][phiW][ptBr].PtB );// / br[etaW][phiW][ptBr].nB);

                    //F
                    int etaWinMod = 2*nEtaBr-1-etaW;
                    binContent = hist1D_multInWin[treeId]->GetBinContent(etaWinMod+1);
                    hist1D_multInWin[treeId]->SetBinContent(etaWinMod+1, binContent + br[etaW][phiW][ptBr].nF);
                    //new:
                    //                hist1D_multDistrInWin[treeId][etaWinMod]->Fill( br[etaW][phiW].nF);
                    //                if ( br[etaW][phiW][ptBr].nF > 0 )
                    //                    hist1D_avPtDistrInWin[treeId][etaWinMod]->Fill( br[etaW][phiW][ptBr].PtF );// / br[etaW][phiW][ptBr].nF);

                }
            }
        
        // cut on zero particles in TPC!
        if ( multTPC < 1 )
            continue;
        
        hist1D_vertexZ[treeId]->Fill( br_vertexZ );
        
        //        cout << "multTPC=" << multTPC << endl;
        hist1D_QA_multALL->Fill(multTPC);
        
        hist2D_ESTIMATOR_VS_multTPC->Fill( cEstimator, multTPC );
        hist2D_QA_multBit32_vs_multWithTOFtiming->Fill( fEvSel_multTrk_bit32, fEvSel_multTrk_bit32_TOF );

        nAccepted_PRE_LOOP++;
    } // end of pre-loop
    
    
    cout << ">>> nAccepted_PRE_LOOP = " << nAccepted_PRE_LOOP << endl;
    cout << ">>> nRejectedByPileup_PRE_LOOP = " << nRejectedByPileup_PRE_LOOP << endl;
    
    
    //write run-by-run QA mult histos to file
    TFile *file_QA_mult_in_eta_wins_RunByRun = new TFile( "histos_QA_mult_in_eta_wins_run-by-run.root", "RECREATE" );
    for ( int i = 0; i < nTrees; i++ )
    {
        //        for ( int etaW = 0; etaW < 2*nEtaBr; etaW++ )
        hist1D_multInWin[i]->Write();
        hist1D_vertexZ[i]->Write();
        //        for ( int etaW = 0; etaW < nEtaBr; etaW++ )
        //        {
        //            hist1D_multDistrInWin[i][etaW]->Write();
        //            hist1D_avPtDistrInWin[i][etaW]->Write();
        //        }
    }
    file_QA_mult_in_eta_wins_RunByRun->Close();
    
    
    // ########## QA PLOTTING:
    
    TCanvas *canv_QA_multBit32_vs_multWithTOFtiming = new TCanvas("canv_QA_multBit32_vs_multWithTOFtiming","canv_QA_multBit32_vs_multWithTOFtiming",2,2,600,600 );
    tuneCanvas(canv_QA_multBit32_vs_multWithTOFtiming);
    canv_QA_multBit32_vs_multWithTOFtiming->SetLogz();
    tuneHist2D(hist2D_QA_multBit32_vs_multWithTOFtiming);
    hist2D_QA_multBit32_vs_multWithTOFtiming->DrawCopy("colz");


    //percentiles QA hist
    TCanvas *canv_estimatorPercentiles_QA_all = new TCanvas("canv_estimatorPercentiles_QA_all","canv_estimatorPercentiles_QA_all",0,0,700,600 );
    tuneCanvas(canv_estimatorPercentiles_QA_all);
    hist1D_QA_percentilesEstimator->DrawCopy();
    
    TCanvas *canv_hist1D_QA_multALL = new TCanvas("canv_hist1D_QA_multALL","canv_hist1D_QA_multALL",50,50,700,600 );
    tuneCanvas(canv_hist1D_QA_multALL);
    canv_hist1D_QA_multALL->SetLogy();
    hist1D_QA_multALL->DrawCopy();
    
    TCanvas *canv_ESTIMATOR_VS_multTPC = new TCanvas("canv_ESTIMATOR_VS_multTPC","canv_ESTIMATOR_VS_multTPC",10,10,800,800 );
    tuneCanvas(canv_ESTIMATOR_VS_multTPC);
    canv_ESTIMATOR_VS_multTPC->SetLogz();
    tuneHist2D(hist2D_ESTIMATOR_VS_multTPC);
    hist2D_ESTIMATOR_VS_multTPC->DrawCopy("colz");
    


    canv_ESTIMATOR_VS_multTPC->cd();

    // ##### FUNCTIONAL CUT FOR OUTLIERS:
    TF1 *fBorderToCutOutliersLower = 0x0;
    TF1 *fBorderToCutOutliersUpper = 0x0;
    
    // ### NEW BOUNDERS: Try mean+/-nsigma in slices:
    TGraphErrors *grMultBoundaryLower = new TGraphErrors;
    TGraphErrors *grMultBoundaryUpper = new TGraphErrors;
    if(1)
    {
        const int stepPercOnHist = 39;
        for ( Int_t perc = 0; perc < hist2D_ESTIMATOR_VS_multTPC->GetNbinsX(); perc+=1+stepPercOnHist )
        {
            TH1D* histSlice= hist2D_ESTIMATOR_VS_multTPC->ProjectionY( "tmpProj", perc+1,perc+1+stepPercOnHist);
            double sliceMean = histSlice->GetMean();
            double sliceSigma = histSlice->GetRMS();
            cout << "sliceMean=" << sliceMean << ", sliceSigma=" << sliceSigma << ", nEntries = " << histSlice->GetEntries() << endl;
            if ( sliceMean > 10 && hist2D_ESTIMATOR_VS_multTPC->GetEntries() > 10 )
            {
                double percBinCenter = hist2D_ESTIMATOR_VS_multTPC->GetXaxis()->GetBinCenter( perc+(1+stepPercOnHist)/2 );
                if(1) // "ANCHOR" boundaries
                {
                    grMultBoundaryLower->SetPoint( grMultBoundaryLower->GetN(), percBinCenter, sliceMean-4*sliceSigma-10 );
                    grMultBoundaryUpper->SetPoint( grMultBoundaryUpper->GetN(), percBinCenter, sliceMean+5*sliceSigma+40 );
                }
                else if(0)
                {
                    grMultBoundaryLower->SetPoint( grMultBoundaryLower->GetN(), percBinCenter, sliceMean-3*sliceSigma-10 );
                    grMultBoundaryUpper->SetPoint( grMultBoundaryUpper->GetN(), percBinCenter, sliceMean+3*sliceSigma+40 );
                }
                else if(0)
                {
                    grMultBoundaryLower->SetPoint( grMultBoundaryLower->GetN(), percBinCenter, sliceMean-2*sliceSigma-10 );
                    grMultBoundaryUpper->SetPoint( grMultBoundaryUpper->GetN(), percBinCenter, sliceMean+2*sliceSigma+40 );
                }
                else if(0)
                {
                    grMultBoundaryLower->SetPoint( grMultBoundaryLower->GetN(), percBinCenter, sliceMean-1*sliceSigma-10 );
                    grMultBoundaryUpper->SetPoint( grMultBoundaryUpper->GetN(), percBinCenter, sliceMean+1*sliceSigma+10 );
                }
                else if(0)
                {
                    grMultBoundaryLower->SetPoint( grMultBoundaryLower->GetN(), percBinCenter, sliceMean-7*sliceSigma-10 );
                    grMultBoundaryUpper->SetPoint( grMultBoundaryUpper->GetN(), percBinCenter, sliceMean+7*sliceSigma+10 );
                }
                else if(0) // BY HAND
                {
                    grMultBoundaryLower->SetPoint( grMultBoundaryLower->GetN(), percBinCenter, 1 );//100 );
                    grMultBoundaryUpper->SetPoint( grMultBoundaryUpper->GetN(), percBinCenter, 10000 );
                }

            }
            
        }
        
        // lower boundary
        grMultBoundaryLower->SetLineColor(kRed);
        grMultBoundaryLower->DrawClone("L");
        fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","[0]+[1]*exp([2]*x)",0,90);
        fBorderToCutOutliersLower->SetParameters( -100, 1500, -0.04 );
        fBorderToCutOutliersLower->SetLineColor( kMagenta);
        //        //TMP!!!! :
        //        fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","10",0,90);
        grMultBoundaryLower->Fit( fBorderToCutOutliersLower, "QN" );
        fBorderToCutOutliersLower->DrawClone( "same" );
        
        // upper boundary
        grMultBoundaryUpper->SetLineColor(kRed);
        grMultBoundaryUpper->DrawClone("L");
        if(1) // to fit real data!
        {
            fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","[0]+[1]*exp([2]*x)",0,90);
            fBorderToCutOutliersUpper->SetParameters( -100, 1500, -0.04 );
            fBorderToCutOutliersUpper->SetLineColor( kMagenta);
            grMultBoundaryUpper->Fit( fBorderToCutOutliersUpper, "QN" );
        }
        else
            fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","2950",0,90);
        
        
        fBorderToCutOutliersUpper->DrawClone( "same" );
        
        //return;
    }
    
    TProfile *prof_ESTIMATOR_VS_multTPC = hist2D_ESTIMATOR_VS_multTPC->ProfileX();
    prof_ESTIMATOR_VS_multTPC->SetLineColor(kBlue+1);
    //    prof_ESTIMATOR_VS_multTPC->DrawCopy("same");
    
    
    if (0) // OLD implementation of boundries BY HAND
    {
        if ( isAnalysing502 == 0 ) //=2.76 default analysis
        {
            fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","-120+1800*exp(-0.042*x)",0,90);
            fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","50+2450*exp(-0.038*x)",0,90);
            
            //tight cuts:
            //        fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-120+2000*exp(-0.042*x)",0,90);
            //        fBorderToCutOutliersUpper = new TF1("fBorderToCutOutliersUpper","-100+2400*exp(-0.034*x)",0,90);
            
            
            //for INEFF study:
            if ( 0 && LIST_ID_BY_HAND > 0 ) // RECREATE LOWER BOUND: to account for ineff loses in mult!
            {
                fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-120+[0]*exp(-0.042*x)",0,90);
                fBorderToCutOutliersLower->SetParameter(0, 1800-170*LIST_ID_BY_HAND );
            }
            //for TPC-only tracks
            if ( 0 && LIST_ID_BY_HAND > 0 ) // RECREATE UPPER BOUND
            {
                fBorderToCutOutliersUpper = new TF1("fBorderToCutOutliersUpper","0+2870*exp(-0.038*x)",0,90);
            }
            //for 272 bit
            if ( 0 && LIST_ID_BY_HAND > 0 ) // RECREATE UPPER BOUND
            {
                fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-100+1550*exp(-0.042*x)",0,90);
                fBorderToCutOutliersUpper = new TF1("fBorderToCutOutliersUpper","0+2500*exp(-0.038*x)",0,90);
            }
            //for TPC clusters 90
            if ( 0 && LIST_ID_BY_HAND > 0 ) // RECREATE UPPER BOUND
            {
                fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-120+1750*exp(-0.042*x)",0,90);
            }
            
            // FOR HIJING:
            if (0)
            {
                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","-120+1600*exp(-0.042*x)",0,90);
                fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","50+3700*exp(-0.03*x)",0,90);
            }
            
            // FOR AMPT spec c30_40:
            if (0)
            {
                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","-500+1000*exp(-0.042*x)",0,90);
                fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","50+3700*exp(-0.03*x)",0,90);
            }
            // JUST ALL EVENTS:
            if (1)
            {
                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","10",-10,100);
                fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","2980",-10,100);
            }
            // FOR AMPT: LHC12a11a-i (AOD157)
            if (0)
            {
                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","10",0,90);
                //            fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","1000",0,90);
                fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","10000",0,90);
                
                //            if ( WHICH_BRANCH_BY_HAND == 3 )
                //                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","8"/*"30"*/,0,90);
                //            if ( WHICH_BRANCH_BY_HAND == 4 ) //weak decays
                //                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","5"/*"15"*/,0,90);
                //            if ( WHICH_BRANCH_BY_HAND == 5 ) //material
                //                fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","5"/*"10"*/,0,90);
                
                if (0)
                {
                    if ( WHICH_BRANCH_BY_HAND == 3 )
                        fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","40",0,90);
                    if ( WHICH_BRANCH_BY_HAND == 4 ) //weak decays
                        fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","15",0,90);
                    if ( WHICH_BRANCH_BY_HAND == 5 ) //material
                        fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","15",0,90);
                }
                
            }
            
        }
        else
        {
            //        fBorderToCutOutliersLower = new TF1("fBorderToCutOutliersLower","-150+2050*exp(-0.042*x)",0,90);
            fBorderToCutOutliersLower = new TF1( "fBorderToCutOutliersLower","-100+2150*exp(-0.042*x)",0,90);
            fBorderToCutOutliersUpper = new TF1( "fBorderToCutOutliersUpper","120+3500*exp(-0.042*x)",0,90);
        }
    } // end of // OLD implementation of boundries BY HAND
    
    //    if (!isAMPT_kine)
    //    {
    fBorderToCutOutliersLower->SetLineColor(kRed+1);
    fBorderToCutOutliersLower->DrawCopy("same");
    
    fBorderToCutOutliersUpper->SetLineColor(kRed+2);
    fBorderToCutOutliersUpper->DrawCopy("same");
    //    }
    
    
    
    //        return;
    
    // ########## Centrality bins:
    double cLowerEdge = 0;

    //    const int nCW = 2; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 9, 18 }; //n centrality bins
    
    //    const int nCW = 3; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5, 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 9, 17, 35 }; //n centrality bins
    
    
    //    const int nCW = 2; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins
    
    //    const int nCW = 3; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 1.0 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5, 1.0 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35, 90 }; //n centrality bins
    
    //    const int nCW = 4; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 1.0, 0.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5, 1.0, 1.0 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35, 90, 90 }; //n centrality bins
    
    //    const int nCW = 5; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2.5, 1.0, 0.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5, 2.5,  2.5, 1.0, 1.0 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 17, 35, 36,  90, 90 }; //n centrality bins
    
    //USED IN FINAL ANALYSIS (width=10%):
    //    const int nCW = 1; //nCentrWidths
    //    const double cWidths[nCW] = { 10 }; //width of the centrality bins
    //    const double cStep[nCW] = { 5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 15 }; //n centrality bins
    
    
    //    const int nCW = 1; //nCentrWidths // AMPT c30_40
    //    const double cWidths[nCW] = { 60 }; //100 }; //width of the centrality bins
    //    const double cStep[nCW] = { 60 };   //{ 100 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 1 }; //n centrality bins
    
    
    //    const int nCW = 1; //nCentrWidths // AMPT c30_40
    //    const double cWidths[nCW] = { 60 }; //100 }; //width of the centrality bins
    //    const double cStep[nCW] = { 60 };   //{ 100 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 1 }; //n centrality bins
    
    // FOR AMPT: SPEC CENTRALITY TUNING 12.06.2016
    //        const int nCW = 1; //nCentrWidths // AMPT c30_40
    //        const double cWidths[nCW] = { 200 }; //100 }; //width of the centrality bins
    //        const double cStep[nCW] = { 200 };   //{ 100 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 1 }; //n centrality bins
    //        cLowerEdge = -10;
    
    
    //USED IN FINAL ANALYSIS (width=20%):
    //        const int nCW = 1; //nCentrWidths
    //        const double cWidths[nCW] = { 20 }; //width of the centrality bins
    //        const double cStep[nCW] = { 20 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 4 }; //10 }; //n centrality bins
    
    
    //USED IN FINAL ANALYSIS (width=10%):
//        const int nCW = 1; //nCentrWidths
//        const double cWidths[nCW] = { 10 }; //width of the centrality bins
//        const double cStep[nCW] = { 10 }; //centrality bins step
//        const int nCentrBins[nCW] = { 8 }; //10 }; //n centrality bins

    //        const int nCW = 20; //nCentrWidths
    //        const double cWidths[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }; //width of the centrality bins
    //        const double cStep[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }; //centrality bins step
    //        const int nCentrBins[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }; //n centrality bins
    
    
    
    //USED IN FINAL ANALYSIS (width=5%):
//        const int nCW = 1; //nCentrWidths
//        const double cWidths[nCW] = { 5 }; //width of the centrality bins
//        const double cStep[nCW] = { 5 }; //centrality bins step
//        const int nCentrBins[nCW] = { 16 }; //n centrality bins

    //USED IN FINAL ANALYSIS (width=2.5%):
    //    const int nCW = 1; //nCentrWidths
    //    const double cWidths[nCW] = { 2.5 }; //width of the centrality bins
    //    const double cStep[nCW] = { 2.5 }; //centrality bins step
    //    const int nCentrBins[nCW] = { 34 };//18 };//34 }; //n centrality bins
    
    //USED IN FINAL ANALYSIS (width=2%):
//    const int nCW = 1; //nCentrWidths
//    const double cWidths[nCW] = { 2 }; //width of the centrality bins
//    const double cStep[nCW] = { 2 }; //centrality bins step
//    const int nCentrBins[nCW] = { 40 };//18 };//34 }; //n centrality bins

    //USED IN FINAL ANALYSIS (width=1%):
//        const int nCW = 1; //nCentrWidths
//        const double cWidths[nCW] = { 1 }; //2 }; //width of the centrality bins
//        const double cStep[nCW] = { 1 }; //centrality bins step
//        const int nCentrBins[nCW] = { 80 };//18 };//34 }; //n centrality bins
    
    //USED IN FINAL ANALYSIS (width=0.5%):
//    const int nCW = 1; //nCentrWidths
//    const double cWidths[nCW] = { 0.5  };  //width of the centrality bins
//    const double cStep[nCW] = { 0.5  }; //centrality bins step
//    const int nCentrBins[nCW] = { 40 };//160 }; //42  };//18 };//34 }; //n centrality bins
//    cLowerEdge = 20;

//        const int nCW = 1; //nCentrWidths
//        const double cWidths[nCW] = { 0.5  };  //width of the centrality bins
//        const double cStep[nCW] = { 0.5  }; //centrality bins step
//        const int nCentrBins[nCW] = { 160 }; //n centrality bins


    //        const int nCW = 2; //nCentrWidths
    //        const double cWidths[nCW] = { 10, 5.001 }; //width of the centrality bins
    //        const double cStep[nCW] = { 5, 2.5 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 17, 35 }; //n centrality bins
    
    
    //USED IN FINAL ANALYSIS: many centralities in one run
    //        const int nCW = 4; //nCentrWidths
    //        const double cWidths[nCW] = { 10, 5, 2, 1 };  //width of the centrality bins
    //        const double cStep[nCW] = { 10, 5, 2, 1 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 8, 17, 42, 83 };//18 };//34 }; //n centrality bins
    
    //    // 9.06.2016: FINAL BINNING:
//    const int nCW = 3; //nCentrWidths
//    const double cWidths[nCW] = { 10, 5, 2  };  //width of the centrality bins
//    const double cStep[nCW] = { 10, 5, 2  }; //centrality bins step
//    const int nCentrBins[nCW] = { 8, 16, 40 }; //42  };//18 };//34 }; //n centrality bins
//    cLowerEdge = 0.0;

    const int nCW = 1; //nCentrWidths
    const double cWidths[nCW] = { 2  };  //width of the centrality bins
    const double cStep[nCW] = { 2  }; //centrality bins step
    const int nCentrBins[nCW] = { 1 }; //42  };//18 };//34 }; //n centrality bins
//    cLowerEdge = 7;


    // Sept. 2016 Anchor analysis
//            const int nCW = 5; //nCentrWidths
//            const double cWidths[nCW] = { 10, 5, 2, 1  };  //width of the centrality bins
//            const double cStep[nCW] = { 10, 5, 2, 1  }; //centrality bins step
//            const int nCentrBins[nCW] = { 8, 16, 40, 80 }; //42  };//18 };//34 }; //n centrality bins
//            const double cWidths[nCW] = { 10, 5, 2, 1, 0.5  };  //width of the centrality bins
//            const double cStep[nCW] = { 10, 5, 2, 1, 0.5  }; //centrality bins step
//            const int nCentrBins[nCW] = { 8, 16, 40, 80, 160 }; //42  };//18 };//34 }; //n centrality bins
//            const double cWidths[nCW] = { 10, 5, 2, 1, 0.5  };  //width of the centrality bins
//            const double cStep[nCW] = { 10, 5, 2, 1, 0.5  }; //centrality bins step
//            const int nCentrBins[nCW] = { 4, 9, 22, 44, 88 }; //42  };//18 };//34 }; //n centrality bins



//    const int nCW = 3; //STUDY MOST CENTRAL:
//    const double cWidths[nCW] = { 1, 0.5, 0.25  };  //width of the centrality bins
//    const double cStep[nCW] = { 1, 0.5, 0.25  }; //centrality bins step
//    const int nCentrBins[nCW] = { 5, 10, 20 }; //42  };//18 };//34 }; //n centrality bins
//    cLowerEdge = 0.;

    //            const int nCW = 5; //nCentrWidths
    //            const double cWidths[nCW] = { 10, 5, 2, 1, 0.5  };  //width of the centrality bins
    //            const double cStep[nCW] = { 10, 5, 2, 1, 0.5  }; //centrality bins step
    //            const int nCentrBins[nCW] = { 8, 16, 40, 80, 160 }; //42  };//18 };//34 }; //n centrality bins

    
    //        const int nCW = 2; //nCentrWidths
    //        const double cWidths[nCW] = { 10, 5 };  //width of the centrality bins
    //        const double cStep[nCW] = { 10, 5  }; //centrality bins step
    //        const int nCentrBins[nCW] = { 8, 16 };  //n centrality bins
    
    // 13.06.2016: BINNING FOR AMPT MERGED BY HAND:
    //        const int nCW = 2; //nCentrWidths
    //        const double cWidths[nCW] = { 10, 5 };  //width of the centrality bins
    //        const double cStep[nCW] = { 10, 5 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 10, 20 };  //n centrality bins
    
    // 28.06.2016: most centr 0-2% with 32 phi-bins:
    //    cLowerEdge = 1; // !!!!!!!!!!
    
    //    const int nCW = 1;//1; //nCentrWidths
    //    const double cWidths[] = { 5 };  //width of the centrality bins
    //    const double cStep[] = { 5 }; //centrality bins step
    //    const int nCentrBins[] = { 1 };  //n centrality bins
    //    cLowerEdge = 0;
    //        const int nCW = 1;//1; //nCentrWidths
    //        const double cWidths[] = { 100 };  //width of the centrality bins
    //        const double cStep[] = { 100 }; //centrality bins step
    //        const int nCentrBins[] = { 1 };  //n centrality bins
    //        cLowerEdge = 0;
    
    
    //USED IN EPOS PHI32 ANALYSIS (width=5%) => 100% divided by 20!
    //        const int nCW = 1; //nCentrWidths
    //        const double cWidths[nCW] = { 5 }; //width of the centrality bins
    //        const double cStep[nCW] = { 5 }; //centrality bins step
    //        const int nCentrBins[nCW] = { 20 }; //n centrality bins
    
    //USED IN EPOS FULL PHI ANALYSIS => divide 100% range into bins! (by MultTPC)
    //    const int nCW = 4; //nCentrWidths
    //    const double cWidths[nCW] = { 10, 5, 2, 1  };  //width of the centrality bins
    //    const double cStep[nCW] = { 10, 5, 2, 1  }; //centrality bins step
    //    const int nCentrBins[nCW] = { 10, 20, 50, 100 }; //42  };//18 };//34 }; //n centrality bins
    
    //TRY NN CORRECTION - USE FEW CENTRALITY WINDOWS (width=5%):
//    const int nCW = 1; //nCentrWidths
//    const double cWidths[nCW] = { 5 }; //width of the centrality bins
//    const double cStep[nCW] = { 10 }; //centrality bins step
//    const int nCentrBins[nCW] = { 4 }; //n centrality bins

//        const int nCW = 1; //try HIJING ZDCvsZEM 10-45%
//        const double cWidths[nCW] = { 5 }; //width of the centrality bins
//        const double cStep[nCW] = { 5 }; //centrality bins step
//        const int nCentrBins[nCW] = { 7 }; //n centrality bins
//        cLowerEdge = 10;

    
    // Split QA mult hist into quantiles: FOR MULT CLASSES
    double **multBounds = new double*[nCW];
    double **multBinCenters = new double*[nCW]; //by mean of 1D-histograms in bins!
    
    double **boundsMin = new double*[nCW];
    double **boundsMax = new double*[nCW];
    
    
    TCanvas *canvColoredMultClasses[nCW];
    for ( int cW = 0; cW < nCW; cW++ )
    {
        int nCBinsForQuant = nCentrBins[cW];//+1;
        cout << "###### Quantiles for estimator: n bins = " << nCBinsForQuant << endl;
        
        multBounds[cW] = new double[nCBinsForQuant]; // array to contain the quantiles
        getQuantiles(hist1D_QA_multALL, nCBinsForQuant, multBounds[cW]);
        
        // prepare canvas with colored bins
        TString canvName = Form("canv_mult_withClasses_byMultTPC_%d_bins_cW%d", nCBinsForQuant, cW );
        canvColoredMultClasses[cW] = new TCanvas( canvName, canvName,300,100,800,600 );
        tuneCanvas(canvColoredMultClasses[cW]);
        
        // draw mult classes, get mult bin centers:
        multBinCenters[cW] = new double[nCBinsForQuant]; // array to contain means of 1D-histograms in bins
        drawCanvasWithClasses(hist1D_QA_multALL//, //Form("byMultTPC_%d_bins_cW%d", nCBinsForQuant, cW )
                              , nCBinsForQuant, multBounds[cW], multBinCenters[cW]);//, canvColoredClasses[cW] );
        
        //        cout << "test canvas name: " << canvColoredMultClasses[cW]->GetName() << endl;
        
        int nCBins = nCentrBins[cW];
        boundsMin[cW] = new double[nCBins];
        boundsMax[cW] = new double[nCBins];
        // for "overlapping" multBins:
        //        rearrangeBoundaries(nCBins, multBounds[cW], boundsMin[cW], boundsMax[cW] );
        
        for ( int bin = 0; bin < nCBins; bin++ )
        {
            boundsMin[cW][bin] = ( bin > 0 ? multBounds[cW][bin-1] : 0);
            boundsMax[cW][bin] = multBounds[cW][bin];
            cout << "bounds for bin=" << bin << ": " << boundsMin[cW][bin] << " " << boundsMax[cW][bin]
                    << ", multBinCenter=" << multBinCenters[cW][bin] << endl;
        }
        
        
        // for "overlapping" multBins:
        //        for ( int bin = 0; bin < nCBins; bin++ )
        //            multBinCenters[cW][bin] = (multBinCenters[cW][bin]+multBinCenters[cW][bin+1])/2;
        
    }
    
    //    int nCentrBinsMult = 10;
    //    cout << "nCentrBins=" << nCentrBinsMult << endl;
    //    double *estBounds = new double[nCentrBinsMult]; // array to contain the quantiles
    //    getQuantiles(hist1D_QA_multALL, nCentrBinsMult, estBounds);
    //    drawCanvasWithClasses( hist1D_QA_multALL, "byMultTPC", nCentrBinsMult, estBounds );
    
    
    
    //    return;
    
    

    
    
    // ########## Select and initiate windows:
    //    const int nEtaWins = 3;
    //    const int howMany = 1; //how many windows to merge from branches
    //    const int nEtaWins = 2;
    //    const int howMany = 4; //how many windows to merge from branches
    //    const int nEtaWins = 8;//1;
    
    // #### 7 eta-wins of 0.2 size with 0.1 step (full phi):
    //    const int nEtaWins = 7;//1;
    //    const int howMany = 2;//4;//2;//1;//4; //how many windows to merge from branches
    //    const int etaStep = 1;  // step for eta win (in terms of number of branches)

    // #### 13 eta-wins of 0.2 size with 0.1 step (full phi) - FOR TRANSLATIONAL INVARIANCE
    //    const int nEtaWins = 13;//7;//1;
    //    const int howMany = 1;//2;//4;//2;//1;//4; //how many windows to merge from branches
    //    const int etaStep = 1;  // step for eta win (in terms of number of branches)

    // #### 1 eta-win of 0.4 size (full phi):
    //        const int nEtaWins = 1;//2;
    //        const int howMany = 4;//1; //how many windows to merge from branches
    //        const int etaStep = 1;  // step for eta win (in terms of number of branches)
    
    // #### 1 eta-win of 0.4 size (full phi) - FOR NN CORRECTIONS ANALYSIS
    //    const int nEtaWins = 1;//2;
    //    const int howMany = 1;//1; //how many windows to merge from branches
    //    const int etaStep = 1;  // step for eta win (in terms of number of branches)

    // #### PHI 32 bins study: 2 eta-wins of 0.4 and 1.6 size (which of the two - is set up by nEtaIdSPEC)
    const int nEtaWins = 1;//7;//1;//2;
    const int howMany = 1;//4;//2;//1;//4;//1; //how many windows to merge from branches
    const int etaStep = 1;  // step for eta win (in terms of number of branches)

    //return;
    
    const int maxNCentrBins = 1;//40;//160;//80;//100;//1;//100;//50; //TMath::MaxElement(nCW, &nCentrBins);
    //    WinPair wins[nCW][maxNCentrBins][nEtaWins][nPhiWins];
    WinPair wins[nCW][maxNCentrBins][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr]; // try SEPARATE WINS FOR PHI ROTATIONS (29.06.2016)
    WinPair winsTmpForCorrNN[nCW][maxNCentrBins][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr]; // try SEPARATE WINS FOR PHI ROTATIONS (29.06.2016)
    CentralityOccupancy cOccupancy[nCW][maxNCentrBins];
    
    bool useSameWinVarForCorrectionNN = false;//true;//false; // for correction NN!

    bool useCentrPercOrMult = 0;//1; // 0 - perc bins, 1-mult bins
    bool doBootstrap = 0;//1;
    bool doBS_in_only_1_cBin = 0;//1;
    bool flag_fill_run_by_run_histos_in_wins = 0; //from 23.03.2016
    bool doEventMixing = 0; // put =1 ONLY WHEN doBootstrap==1 (because we keep events)

    bool doInitHistos = false;
    
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {
            
            float cBinMin = 0;
            float cBinMax = 0;
            
            if ( useCentrPercOrMult==0 ) //bins according to centrality percentiles
            {
                cBinMin = cLowerEdge + cStep[cW] * cBin;
                cBinMax = cLowerEdge + cWidths[cW] + cStep[cW] * cBin;
            }
            else //bins according to mult bins
            {
                //                cBinMin = ( cBin==0 ? 0 : multBounds[cW][cBin-1] );
                //                cBinMax = multBounds[cW][cBin];
                cBinMin = boundsMin[cW][cBin];
                cBinMax = boundsMax[cW][cBin];
            }
            
            cOccupancy[cW][cBin].cBinMin = cBinMin;
            cOccupancy[cW][cBin].cBinMax = cBinMax;
            

            for ( int ptBr = 0; ptBr < nPtBr*nPtBr; ptBr++ )
            {
                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                    for ( int phiW = 0; phiW < nPhiWins*nPhiWins; phiW++ )
                    {
                        int _ptW = nPtBr>1 ? ptBr : -1;
                        // TMP for corrections for NN!
                        winsTmpForCorrNN[cW][cBin][etaW][phiW][ptBr].init(cW, cBinMin, cBinMax, etaW, phiW, _ptW, doInitHistos );

                        if ( !doBootstrap || ( doBootstrap && doBS_in_only_1_cBin && cBin != 0 ) ) // CHECK!!!
                            wins[cW][cBin][etaW][phiW][ptBr].init(cW, cBinMin, cBinMax, etaW, phiW, _ptW, doInitHistos );
                        else
                        {
                            int nEventsForSubsampling = /*nPhiWins* */nAccepted_PRE_LOOP;
                            wins[cW][cBin][etaW][phiW][ptBr].init(cW, cBinMin, cBinMax, etaW, phiW, _ptW , doInitHistos, nEventsForSubsampling );
                        }

                        if( flag_fill_run_by_run_histos_in_wins )
                            wins[cW][cBin][etaW][phiW][ptBr].initRunByRunHistos( nTrees, runListNumbers );
                    }
            }
        }
    
    
    
    
    
    //    return;
    
    // 12.09.2016:  prepare array of e-by-e graphs
    const int nGraphsEbyE = 40;
    int counterGraphsEbyE = 0;
    //    TGraphErrors *grEbyE_nB_inPhiBins[nGraphsEbyE];
    //    TGraphErrors *grEbyE_avPtB_inPhiBins[nGraphsEbyE];
    TH1D *grEbyE_nB_inPhiBins[nGraphsEbyE];
    TH1D *grEbyE_avPtB_inPhiBins[nGraphsEbyE];
    for (int i=0; i < nGraphsEbyE; i++)
    {
        //        grEbyE_nB_inPhiBins[i] = new TGraphErrors;
        //        grEbyE_avPtB_inPhiBins[i] = new TGraphErrors;
        grEbyE_nB_inPhiBins[i] = new TH1D( Form( "grEbyE_nB_ev%d", i ), ";#varphi bin;value", nPhiWins, -0.5, nPhiWins-0.5 );
        grEbyE_avPtB_inPhiBins[i] = new TH1D( Form( "grEbyE_avPtB_ev%d", i ), ";#varphi bin;value", nPhiWins, -0.5, nPhiWins-0.5 );
        
    }
    
    
    // ##### main loop over events
    int nAccepted = 0;
    int nRejectedByPileup = 0;
    treeId = 0;
    sumEvPrevTrees = 0; //nEvInTrees[0];
    for (Long64_t i=0; i < nEvents; i++)
    {
        if ( i % 100000 == 0 )
            cout << "getting " << (int)i << endl;
        //                cout <<"getting " << (int)i << "\r"; cout.flush();
        
        if ( i >= sumEvPrevTrees + nEvInTrees[treeId] ) // it's time to go to next tree...
        {
            sumEvPrevTrees += nEvInTrees[treeId]; // sum of events in prev trees
            treeId++;
            cout << "check treeId = " << treeId << " brV0M = " << brV0M << endl;
        }
        trees[treeId]->GetEntry( i - sumEvPrevTrees );
        
        //        cout << "check i = " << i - (treeId>0 ? nEvInTrees[treeId-1] : 0) << endl;
        
        
        //        if ( i >= nEvInTrees[treeId] )
        //        {
        //            treeId++;
        ////            if ( i < 5 )
        //                cout << "check treeId = " << treeId << " brV0M = " << brV0M << endl;
        //        }
        //        trees[treeId]->GetEntry( i - (treeId>0 ? nEvInTrees[treeId-1] : 0) );
        
        //        cout << "check i = " << i - (treeId>0 ? nEvInTrees[treeId-1] : 0) << endl;
        
        //        if ( i < nEvInTrees[0] )
        //        {
        //            trees[0]->GetEntry( i );
        //            if ( i < 5 )
        //                cout << "check i = " << i << " brV0M = " << brV0M << endl;
        //        }
        //        else
        //        {
        //            trees[1]->GetEntry( i-nEvInTrees[0] );
        //            if ( i < nEvInTrees[0]+5 )
        //                cout << "check i = " << i << " brV0M = " << brV0M << endl;
        //        }
        
        
        // ### event selection
        float cEstimator;
        bool isEventSelected = checkEventSelection(
                    cEstimator

                    , fEvSel_isPileupSPD
                    , fEvSel_isPileupFromMV
                    , fEvSel_vtxSPD_isFromVertexer3D
                    , fEvSel_vtxSPD_isFromVertexerZ
                    , fEvSel_zRes_above_fMaxResol
                    , fEvSel_diff_vZ_SPD_tracks_above_05
                    , fEvSel_comb_cut_on_abs_and_nsigma_dist
                    , fEvSel_too_many_TPC_clusters

                    , fEvSel_multTrk_bit32
                    , fEvSel_multTrk_bit32_TOF

                    , brZEMvsZDC
                    , brV0M
                    , brCL1
                    , br_vertexZ
                    , nRejectedByPileup
                    );
        
        if ( !isEventSelected )
        {
            //            cout << "TEST: cEstimator=" << cEstimator << endl;
            continue;
        }
        
        //calc mult in whole TPC using wins:
        int multTPC = 0;
        for ( int ptBr = 0; ptBr < nPtBr; ptBr++ )
            for ( int phiW = 0; phiW < nPhiWins; phiW++ )
            {
                if ( FLAG_ONE_WIDE_ETA_WIN )
                    multTPC += br[0][phiW][ptBr].nF; // + br[etaW][phiW][ptBr].nB;
                else if ( !FLAG_PERCOLATING_WINS ) // usual case!
                    for ( int etaW = 0; etaW < nEtaBr; etaW++ )
                        multTPC += br[etaW][phiW][ptBr].nF + br[etaW][phiW][ptBr].nB;
                else
                    multTPC += (br[0][phiW][ptBr].nF + br[2][phiW][ptBr].nF) + (br[0][phiW][ptBr].nB + br[2][phiW][ptBr].nB);
            }
        
        // cut on zero particles in TPC!
        if ( multTPC < 1 )
            continue;
        
        
        // !!!! test cut by line on Perc_vs_mult plot:
        //        if ( !isAMPT_kine )
        if ( multTPC < fBorderToCutOutliersLower->Eval(cEstimator)
             || multTPC > fBorderToCutOutliersUpper->Eval(cEstimator) )
            continue;
        
        //BY HAND!!! FOR EPOS MOST CENTRAL STUDY!
        //        if ( multTPC < 231 )//300 ) // for SEP ETA
        //        if ( multTPC < 487 ) //640 )  // for FULL ETA
        //            continue;
        
        
        //assign "centrality" for this event: either cPerc or multTPC
        float centrValue = ( useCentrPercOrMult==0 ? cEstimator : multTPC );

        // TMP!!! 18.11.2016: check low-centr V0M cut to figure out jumping point at c=0:
//        if ( cEstimator < 0.1 )
//            continue;
        
        //prepare vars
        UShort_t _nF, _nB;
        Float_t _nF_PtF, _nB_PtB, _PtF, _PtB;
        
        UShort_t _nF1, _nB1;
        Float_t _nF_PtF1, _nB_PtB1, _PtF1, _PtB1;
        
        //
        UShort_t _nFj; // spec for correction of NN!
        Float_t _nFj_PtFj, _PtFj; // spec for correction of NN! // ... but also keep PtFj for completeness

        //fill occupancy info
        for ( int cW = 0; cW < nCW; cW++ )
            for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
                cOccupancy[cW][cBin].fill(brV0M, brZEMvsZDC);
        
        //fill wins
        for ( int ptFwd = 0; ptFwd < nPtBr; ptFwd++ ) // ptFwd branch
            for ( int ptBwd = 0; ptBwd < nPtBr; ptBwd++ ) // ptBwd branch
            {
                for ( int phiR = 0; phiR < nPhiWins; phiR++ ) // "rotations" of FB phi-win pair
                {
                    for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                        //            for ( int phiW = 0; phiW < nPhiWins/2; phiW++ ) // 13.08.2016: study with mirrored backward windows
                    {
                        int phiF = phiR;
                        int phiB;

                        int phiFadd, phiBmir;

                        if(1) // BASIC STUDY
                        {
                            phiB = phiR+phiW;
                            if ( phiB >= nPhiWins )
                                phiB -= nPhiWins;
                        }
                        // 13.08.2016: study with mirrored backward windows
                        else
                        {
                            // F: add neighbour window to forward
                            phiFadd = phiF+1;
                            if ( phiFadd >= nPhiWins )
                                phiFadd -= nPhiWins;

                            // B
                            phiB = phiR+1+phiW;
                            if ( phiB >= nPhiWins )
                                phiB -= nPhiWins;

                            // Bmir
                            phiBmir = phiR-phiW;
                            if ( phiBmir < 0 )
                                phiBmir += nPhiWins;

                            // print (QA):
                            if ( nAccepted == 0 )
                                cout << "phiWins: F, Fadd, B, Bmir = " << phiF << ", " << phiFadd << ", " << phiB << ", " << phiBmir << endl;
                        }

                        //                    break;    // DON'T DO SAME FB PAIR!!! (even if F and B are swapped)

                        // TMP??? 31.10.2016: remember array with data in each eta win (n and pt)
                        double arr_n_inEta[2*nEtaBr];
                        double arr_Pt_inEta[2*nEtaBr];
                        if(0)for ( int specEtaId = 0; specEtaId < nEtaBr; specEtaId++ )
                        {
                            arr_n_inEta[specEtaId] = br[specEtaId][phiF][0].nF;
                            arr_n_inEta[2*nEtaBr-1-specEtaId] = br[specEtaId][phiB][0].nB;

                            arr_Pt_inEta[specEtaId] = br[specEtaId][phiF][0].nF * br[specEtaId][phiF][0].PtF;
                            arr_Pt_inEta[2*nEtaBr-1-specEtaId] = br[specEtaId][phiF][0].nB * br[specEtaId][phiF][0].PtB;
                        }

                        // QA printing:
                        if(0 && centrValue<5) for ( int specEtaId = 0; specEtaId < nEtaBr*2; specEtaId++ )
                        {
                            double nInBin = arr_n_inEta[specEtaId];
                            double ptSum = arr_Pt_inEta[specEtaId];
                            cout << "specEtaId = " << specEtaId
                                 << ", n in eta = " << nInBin
                                 << ", avPt in eta = " << (nInBin>0 ? ptSum/nInBin : 0) << endl;
                        }



                        for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                        {
                            // 27.02.16: fill wins "by hand" using info from branches:
                            _nF = 0;
                            _nB = 0;
                            _nF_PtF = 0;
                            _nB_PtB = 0;
                            _PtF = -1; // important to put -1!
                            _PtB = -1; // important to put -1!

                            if(useSameWinVarForCorrectionNN) // TMP!!! for correction NN!
                            {
                                _nFj = 0;
                                _nFj_PtFj = 0;
                                _PtFj = -1; // important to put -1!
                            }

                            //for spec_reco
                            if ( isUseSpecRecoWithKineForResponse )
                            {
                                _nF1 = 0;
                                _nB1 = 0;
                                _nF_PtF1 = 0;
                                _nB_PtB1 = 0;
                                _PtF1 = -1; // important to put -1!
                                _PtB1 = -1; // important to put -1!
                            }

                            //                    for ( int _eW = howMany*etaW; _eW < howMany*(etaW+1); _eW++ )
                            for ( int _eW = etaStep*etaW; _eW < etaStep*etaW+howMany; _eW++ )
                            {
                                // cout << "etaW=" << etaW << ", _eW=" << _eW << endl;
                                //                        int a;
                                //                        cin >> a;

                                //BASIC STUDY:
                                if(1)
                                {
                                    _nF  += br[_eW][phiF][ptFwd].nF;
                                    _nB  += br[_eW][phiB][ptBwd].nB;
                                    _nF_PtF += br[_eW][phiF][ptFwd].nF * br[_eW][phiF][ptFwd].PtF;
                                    _nB_PtB += br[_eW][phiB][ptBwd].nB * br[_eW][phiB][ptBwd].PtB;

                                    // TMP!!! for correction NN!
                                    if(useSameWinVarForCorrectionNN)
                                    {
                                        _nFj  += br[_eW][phiB][ptBwd].nF; // nF, but for ptBwd!!! (also phiB?)
                                        _nFj_PtFj += br[_eW][phiB][ptBwd].nF * br[_eW][phiB][ptBwd].PtF;
                                    }
                                }
                                // 13.08.2016: study with mirrored backward windows
                                else if (0)
                                {
                                    _nF  += br[_eW][phiF][ptFwd].nF + br[_eW][phiFadd][ptFwd].nF;
                                    _nB  += br[_eW][phiB][ptBwd].nB + br[_eW][phiBmir][ptBwd].nB;
                                    _nF_PtF += br[_eW][phiF][ptFwd].nF * br[_eW][phiF][ptFwd].PtF + br[_eW][phiFadd][ptFwd].nF * br[_eW][phiFadd][ptFwd].PtF;
                                    _nB_PtB += br[_eW][phiB][ptBwd].nB * br[_eW][phiB][ptBwd].PtB + br[_eW][phiBmir][ptBwd].nB * br[_eW][phiBmir][ptBwd].PtB;
                                }
                                // 31.10.2016: shift wins with gap=0 (check "translational invariance")
                                else
                                {
                                    // if howMany=1 and step=1:
                                    _nF  += arr_n_inEta[_eW] + arr_n_inEta[_eW+1];
                                    _nB  += arr_n_inEta[_eW+2] + arr_n_inEta[_eW+3];
                                    _nF_PtF += arr_Pt_inEta[_eW] + arr_Pt_inEta[_eW+1];
                                    _nB_PtB += arr_Pt_inEta[_eW+2] + arr_Pt_inEta[_eW+3];
                                }

//                                cout << "centrValue=" << centrValue << "  ptFwd=" << ptFwd << ", ptB=" << ptBwd << "    " << _nF << " " << _nB << endl;

                                // ONLY F or B! (used for study in SAME eta wins like (-0.8, -0.4) or (0.4, 0.8) )
                                //                        _nF  += br[_eW][phiF].nB;
                                //                        _nB  += br[_eW][phiB].nB;
                                //                        _nF_PtF += br[_eW][phiF].nB * br[_eW][phiF].PtB;
                                //                        _nB_PtB += br[_eW][phiB].nB * br[_eW][phiB].PtB;

                                // isUseSpecRecoWithKineForResponse
                                if ( isUseSpecRecoWithKineForResponse )
                                {
                                    _nF1  += br_spec_reco[_eW][phiF][ptFwd].nF;
                                    _nB1  += br_spec_reco[_eW][phiB][ptBwd].nB;
                                    _nF_PtF1 += br_spec_reco[_eW][phiF][ptFwd].nF * br_spec_reco[_eW][phiF][ptFwd].PtF;
                                    _nB_PtB1 += br_spec_reco[_eW][phiB][ptBwd].nB * br_spec_reco[_eW][phiB][ptBwd].PtB;
                                }
                            } // end of eW
                            //cout << "_nF=" << _nF << ", _nB=" << _nB << endl;

                            if ( _nF > 0 )
                                _PtF = _nF_PtF / _nF;
                            if ( _nB > 0 )
                                _PtB = _nB_PtB / _nB;

                            if(useSameWinVarForCorrectionNN) // TMP!!! for correction NN!
                            {
                                if ( _nFj > 0 )
                                    _PtFj = _nFj_PtFj / _nFj;
                            }


                            //cout << "_PtF=" << _PtF << ", _PtB=" << _PtB << endl;

                            // isUseSpecRecoWithKineForResponse
                            if ( isUseSpecRecoWithKineForResponse )
                            {
                                if ( _nF1 > 0 )
                                    _PtF1 = _nF_PtF1 / _nF1;
                                if ( _nB1 > 0 )
                                    _PtB1 = _nB_PtB1 / _nB1;
                            }

                            // loop over centrality bins in all binnings
                            for ( int cW = 0; cW < nCW; cW++ )
                            {
                                //cout << "   cW =" << cW << endl;
                                for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
                                {
                                    //for 0.2-0.8:
                                    // if ( _nF >= 38 && _nF < 52
                                    //      && _nB >= 38 && _nB < 52 )
                                    //for 0.8-10.0:
                                    // if ( _nF >= 15 && _nF < 26
                                    //      && _nB >= 15 && _nB < 26 )
                                    //                            if ( _nF >= 21
                                    //                                 && _nB >= 21 )
                                    //                                if ( _nF >= 2
                                    //                                     && _nB >= 2 )
                                    wins[cW][cBin][etaW][nPhiWins*phiR+phiW][ptFwd*nPtBr+ptBwd].fill( centrValue, _nF, _nB, _PtF, _PtB, treeId );
                                    if(useSameWinVarForCorrectionNN)
                                        winsTmpForCorrNN[cW][cBin][etaW][nPhiWins*phiR+phiW][ptFwd*nPtBr+ptBwd].fill( centrValue, _nF, _nFj, _PtF, _PtFj, treeId );


                                    // 12.09.2016: graphs e-by-e with avPt:
                                    if ( phiR == 0 )
                                    {
                                        if ( counterGraphsEbyE < nGraphsEbyE )
                                        {
                                            //                                    grEbyE_nB_inPhiBins[counterGraphsEbyE]->SetPoint( phiW, TMath::TwoPi()/nPhiWins*phiW, _PtB );
                                            //                                    grEbyE_avPtB_inPhiBins[counterGraphsEbyE]->SetPoint( phiW, TMath::TwoPi()/nPhiWins*phiW, _PtB );
                                            //                                    cout << "phiW = " << phiW << ", _nB=" << _nB << ", _PtB=" << _PtB << endl;

                                            grEbyE_nB_inPhiBins[counterGraphsEbyE]->SetBinContent( phiW+1, _nB );
                                            grEbyE_avPtB_inPhiBins[counterGraphsEbyE]->SetBinContent( phiW+1, _PtB );
                                        }
                                    }

                                    //fill response matrix
                                    if ( isUseSpecRecoWithKineForResponse )
                                        //&& (cBin == 0 && centrValue>40 && centrValue<60) ) //tmp!!! centr bin by hand!
                                    {
                                        if ( etaW==0 && centrValue > cWidths[0]*cBin && centrValue < cWidths[0]*(cBin+1) ) // FOR CENTR BINS 10% WIDTH!
                                        {
                                            if ( _nF > 0 && _nB > 0 )
                                            {
                                                // histFB_PtPt_truth->Fill( _PtF, _PtB );
                                                histFB_PtPt_truth->Fill( _nF );
                                                if ( _nF1 > 0 && _nB1 > 0 )
                                                {
                                                    // histFB_PtPt_rec->Fill( _PtF1, _PtB1 );
                                                    histFB_PtPt_rec->Fill( _nF1 );
                                                    // !!! UNCOMMENT IF NEED UNFOLD-> response_PtPt->Fill( _PtF1, _PtB1, _PtF, _PtB );
                                                    // !!! UNCOMMENT IF NEED UNFOLD 1D ->
                                                    //response_PtPt[cBin].Fill( _nF1, _nF );
                                                }
                                                // !!! UNCOMMENT IF NEED UNFOLD->
                                                //else
                                                // !!! UNCOMMENT IF NEED UNFOLD-> response_PtPt->Miss( _PtF, _PtB );
                                                // !!! UNCOMMENT IF NEED UNFOLD 1D ->
                                                //response_PtPt[cBin].Miss( _nF );
                                            }
                                        }
                                    }
                                } // end of loop over centrality bins in one binning
                            } // end of loop over centrality binnings
                        } // end of eta wins loop
                    } // end of phi wins loop
                } // end of "rotations" of FB phi-win pairs
            } // end of pt B branch
        
        
        // 12.09.2016 - increment counter for graphs e-by-e with avPt:
        counterGraphsEbyE++;
        
        
        //cout << " event " << i << " is analyzed." << endl;
        
        
        //... IA: po-normalnomu bylo tak:
        //                for ( int etaW = 0; etaW < nEtaBr; etaW++ )
        //                    for ( int phiW = 0; phiW < nPhiWins; phiW++ )
        //                        wins[cW][cBin][etaW][phiW].fill( centrValue, br[etaW][phiW].nF, br[etaW][phiW].nB, br[etaW][phiW].PtF, br[etaW][phiW].PtB );
        nAccepted++;
    } // end of events
    cout << "nAccepted = " << nAccepted << endl;
    cout << "nAccepted/nAll = " << (float)nAccepted/nEvents << endl;
    
    cout << ">>> main loop: nRejectedByPileup = " << nRejectedByPileup << endl;
    
    cout << "histFB_PtPt_rec mean F = " << histFB_PtPt_rec->GetMean(1) << endl;
    cout << "histFB_PtPt_rec mean B = " << histFB_PtPt_rec->GetMean(2) << endl;
    cout << "histFB_PtPt_truth mean F = " << histFB_PtPt_truth->GetMean(1) << endl;
    cout << "histFB_PtPt_truth mean B = " << histFB_PtPt_truth->GetMean(2) << endl;
    
    // ########## PREPARE OUTPUT ROOT-FILE:
    TString strOutFile;
    //    TFile *fileOutput = new TFile( "output_histos_graphs.root", "RECREATE" );
    //    strOutFile = Form( "output_histos_graphs_ineff%d.root", LIST_ID_BY_HAND );
    //        strOutFile = Form( "output_histos_graphs_ineff%d_whichBranch%d.root", LIST_ID_BY_HAND, WHICH_BRANCH_BY_HAND );
    //    strOutFile = Form( "output_histos_graphs_nEvents_%d.root", nEventsByHand );
    //    strOutFile = Form( "output_histos_graphs_run_%d.root", runListNumbers[0] );
    
    //    strOutFile = Form( "output_fileId_%d_ineff_%d.root", fileIdByHand, LIST_ID_BY_HAND );
    //    strOutFile = Form( "output_ptW_%d.root", ptW );
    strOutFile = Form( "output_listId_%d.root", LIST_ID_BY_HAND );
    
    
    if (0)
    {
        TString tmpNameAMPT[6] = {
            "gen.root"
            , "allRec.root"
            , "recPrim.root"
            , "nonPrim.root"
            , "secWeak.root"
            , "secMaterial.root"
        };
        //        strOutFile = tmpNameAMPT[ WHICH_BRANCH_BY_HAND ];
        strOutFile = Form( "ineffList%d_%s", LIST_ID_BY_HAND, tmpNameAMPT[ WHICH_BRANCH_BY_HAND ].Data() );
    }
    
    TFile *fileOutput = new TFile( strOutFile, "RECREATE" );
    
    
    if ( isUseSpecRecoWithKineForResponse )
    {
        histFB_PtPt_truth->Write();
        histFB_PtPt_rec->Write();
        //        response_PtPt->Write();
        // !!! UNCOMMENT IF NEED UNFOLD->
        //for (Int_t cBin = 0; cBin < 8; cBin++ ) // for each cBin!
        //    response_PtPt[cBin].Write();
    }
    
    
    
    
    
    // ########## SAVE HISTOS TO ROOT-FILE:
    if(doInitHistos)for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = /*nCentrBins[cW]-3*/ 0; cBin < nCentrBins[cW]; cBin++ )
            for ( int ptBr = 0; ptBr < nPtBr/* *nPtBr*/; ptBr++ ) // ?? should be nPtBr*nPtBr
                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                    // !!! for ( int phiW = 0; phiW < nPhiWins*nPhiWins; phiW++ )
                    for ( int phiW = 0; phiW < nPhiWins; phiW++ )
                    {
                        if(0)
                            wins[cW][cBin][etaW][phiW][ptBr].writeHistos();
                        // for QA: write only mult distr:
                        if(1)
                        {
                            wins[cW][cBin][etaW][phiW][ptBr].hist1D_multDistrF->Write();
                            wins[cW][cBin][etaW][phiW][ptBr].hist1D_multDistrB->Write();

                            wins[cW][cBin][etaW][phiW][ptBr].hist1D_QA_PtF->Write();
                            wins[cW][cBin][etaW][phiW][ptBr].hist1D_QA_PtB->Write();
                        }
                    }
    
    hist1D_QA_percentilesEstimator->Write();
    hist1D_QA_multALL->Write();
    hist2D_ESTIMATOR_VS_multTPC->Write();
    hist2D_ESTIMATOR_VS_multTPC->ProfileX()->Write();
    
    //    fileOutput->WriteObject(canv_ESTIMATOR_VS_multTPC);
    canv_ESTIMATOR_VS_multTPC->Write();
    canv_hist1D_QA_multALL->Write();
    hist2D_QA_multBit32_vs_multWithTOFtiming->Write();
    

    for ( int cW = 0; cW < nCW; cW++ )
    {
        //        cout << "test: writing canvas to file, cW=" << cW << endl;
        //        cout << "test canvas name: " << canvColoredMultClasses[cW]->GetName() << endl;
        canvColoredMultClasses[cW]->Write();
    }
    

    // ########## MAIN PLOTTING FOR CORRS:
    GraphsCorrInfo grNN[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    GraphsCorrInfo grPtPt[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    GraphsCorrInfo grPtN[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    
    GraphsCorrInfo grNN_fromMultF[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    GraphsCorrInfo grPtPt_fromMultF[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    GraphsCorrInfo grPtN_fromMultF[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    
    // tmp for correction of NN!
    GraphsCorrInfo grNNTmpForCorrNN[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];

    TGraphErrors *grFractEstByV0M[nCW];
    TGraphErrors *grFractEstByZDC[nCW];
    
    for ( int cW = 0; cW < nCW; cW++ )
    {
        grFractEstByV0M[cW] = new TGraphErrors;
        grFractEstByZDC[cW] = new TGraphErrors;
        for ( int ptW = 0; ptW < nPtBr*nPtBr; ptW++ )
            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
            {
                for ( int phiW = 0; phiW < nPhiWins*nPhiWins; phiW++ )
                {
                    grNN[cW][etaW][phiW][ptW].SetNames( "NN", cW, etaW, phiW, ptW );
                    grPtPt[cW][etaW][phiW][ptW].SetNames( "PtPt", cW, etaW, phiW, ptW );
                    grPtN[cW][etaW][phiW][ptW].SetNames( "PtN", cW, etaW, phiW, ptW );

                    grNN_fromMultF[cW][etaW][phiW][ptW].SetNames( "NN_fromMultF", cW, etaW, phiW, ptW );
                    grPtPt_fromMultF[cW][etaW][phiW][ptW].SetNames( "PtPt_fromMultF", cW, etaW, phiW, ptW );
                    grPtN_fromMultF[cW][etaW][phiW][ptW].SetNames( "PtN_fromMultF", cW, etaW, phiW, ptW );

                    if(useSameWinVarForCorrectionNN)
                        grNNTmpForCorrNN[cW][etaW][phiW][ptW].SetNames( "NN_SPEC_j", cW, etaW, phiW, ptW );
                }
            }
    }
    

    //calc (1) - occupancies in centr bins, (2) - corr coeffs
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
        {
            if ( cOccupancy[cW][cBin].nEventsV0M > 0 )
            {
                CentralityOccupancy *c = &cOccupancy[cW][cBin];
                float centr = c->cBinMin + (c->cBinMax - c->cBinMin)/2;
                float cRatio = 0;
                if (c->nEventsV0M>0)
                    cRatio = (float)c->nEventsV0M_and_ZDCZEM / c->nEventsV0M;
                grFractEstByV0M[cW]->SetPoint(grFractEstByV0M[cW]->GetN(), centr, cRatio);
                if (c->nEventsZDCZEM>0)
                    cRatio = (float)c->nEventsV0M_and_ZDCZEM / c->nEventsZDCZEM;
                grFractEstByZDC[cW]->SetPoint(grFractEstByZDC[cW]->GetN(), centr, cRatio);
            }
            for ( int ptW = 0; ptW < nPtBr*nPtBr; ptW++ )
                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                    for ( int phiW = 0; phiW < nPhiWins*nPhiWins; phiW++ )
                    {
                        WinPair *w = &wins[cW][cBin][etaW][phiW][ptW];
                        float centr = -1;
                        if ( useCentrPercOrMult == 0 )
                            centr = w->cBinMin + (w->cBinMax - w->cBinMin)/2;
                        else
                            centr = multBinCenters[cW][cBin];

                        double multF = -1;
                        if (doInitHistos) // to avoid seg fault if do not init histos
                            multF = w->hist1D_multDistrF->GetMean();

                        w->calcCorrCoeffs();
                        if(0)cout << "cMin=" << w->cBinMin << ", cMax=" << w->cBinMax << ", etaW=" << etaW
                                  << ", corrInfo_NN.bCorr= " << w->corrInfo_NN.bCorr
                                  << ", corrInfo_PtPt.bCorr= " << w->corrInfo_PtPt.bCorr
                                  << endl;


                        //fill graphs
                        grNN[cW][etaW][phiW][ptW].SetPoints( centr, &w->corrInfo_NN );
                        grPtPt[cW][etaW][phiW][ptW].SetPoints( centr, &w->corrInfo_PtPt );
                        grPtN[cW][etaW][phiW][ptW].SetPoints( centr, &w->corrInfo_PtN );

                        grNN_fromMultF[cW][etaW][phiW][ptW].SetPoints( multF, &w->corrInfo_NN );
                        grPtPt_fromMultF[cW][etaW][phiW][ptW].SetPoints( multF, &w->corrInfo_PtPt );
                        grPtN_fromMultF[cW][etaW][phiW][ptW].SetPoints( multF, &w->corrInfo_PtN );

                        if(useSameWinVarForCorrectionNN)
                        {
                            WinPair *wj = &winsTmpForCorrNN[cW][cBin][etaW][phiW][ptW];
                            wj->calcCorrCoeffs();
                            grNNTmpForCorrNN[cW][etaW][phiW][ptW].SetPoints( centr, &wj->corrInfo_NN );
                        }

                    }
        }



    // ########## BOOTSTRAPING
    GraphsCorrInfo gr_BS_NN[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    GraphsCorrInfo gr_BS_PtPt[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    
    GraphsCorrInfo gr_BS_NN_fromMultF[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    GraphsCorrInfo gr_BS_PtPt_fromMultF[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    
    // ########## Event Mixing
    TGraphErrors *gr_bCorr_MIX_NN[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    TGraphErrors *gr_bCorr_MIX_PtPt[nCW][nEtaWins][nPhiWins*nPhiWins][nPtBr*nPtBr];
    
    if ( doBootstrap )
    {
        const int kSubsamplingType = 1;//0; // 0 - bootstrap, 1 - simple subsampling
        cout << "Start subsampling...";
        if ( kSubsamplingType == 0 )
            cout << " (bootstrap)" << endl;
        else
            cout << " (simple subsampling)" << endl;
        
        TCanvas *canv_bootStrapPhiWins = new TCanvas("canv_bootStrapPhiWins","canv_bootStrapPhiWins",350,150,900,700 );
        
        for ( int cW = 0; cW < nCW; cW++ )
        {
            for ( int ptW = 0; ptW < nPtBr*nPtBr; ptW++ )
            {
                if(1)cout << "  entering ptW " << ptW << "..." << endl;

                for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                {
                    if(1)cout << "  entering etaW " << etaW << "..." << endl;

                    for ( int phiW = 0; phiW < nPhiWins*nPhiWins; phiW++ )
                    {
                        if(1)cout << "  entering phiW " << phiW << "..." << endl;
                        //                int etaW = 0; // TMP, FIXED FOR TESTS!

                        gr_BS_NN[cW][etaW][phiW][ptW].SetNames( "BS_NN", cW, etaW, phiW, ptW );
                        gr_BS_PtPt[cW][etaW][phiW][ptW].SetNames( "BS_PtPt", cW, etaW, phiW, ptW );

                        gr_BS_NN_fromMultF[cW][etaW][phiW][ptW].SetNames( "BS_NN_fromMultF", cW, etaW, phiW, ptW );
                        gr_BS_PtPt_fromMultF[cW][etaW][phiW][ptW].SetNames( "BS_PtPt_fromMultF", cW, etaW, phiW, ptW );

                        //for event mixing:
                        gr_bCorr_MIX_NN[cW][etaW][phiW][ptW] = new TGraphErrors;
                        gr_bCorr_MIX_PtPt[cW][etaW][phiW][ptW] = new TGraphErrors;

                        TString strWinDescr = Form( "gr_evMix_NN_bcorr_cW%d_etaW%d_phiW%d_ptW%d", cW, etaW, phiW, ptW );
                        gr_bCorr_MIX_NN[cW][etaW][phiW][ptW]->SetName( strWinDescr );
                        strWinDescr = Form( "gr_evMix_PtPt_bcorr_cW%d_etaW%d_phiW%d_ptW%d", cW, etaW, phiW, ptW );
                        gr_bCorr_MIX_PtPt[cW][etaW][phiW][ptW]->SetName( strWinDescr );

                        for ( int cBin = /*nCentrBins[cW]-3*/ 0; cBin < nCentrBins[cW]; cBin++ )
                        {
                            WinPair *w = &wins[cW][cBin][etaW][phiW][ptW];

                            float centr = -1;
                            if ( useCentrPercOrMult == 0 )
                                centr = w->cBinMin + (w->cBinMax - w->cBinMin)/2;
                            else
                                centr = multBinCenters[cW][cBin];

                            double multF = w->hist1D_multDistrF->GetMean();

                            int nSubsamples = 400;//2;//10;//400;// !!!! to speed-up for phi wins //1000;
                            //                        nSubsamples = 10*(cW+1);
                            //                        nSubsamples = 200+50*(cW+1);
                            double modFactorForSubsamplingError = 1.;
                            if ( kSubsamplingType == 1 )
                            {
                                //                            nSubsamples = 10*(cW+1);
                                //                            nSubsamples = 200+50*(cW+1);
                                modFactorForSubsamplingError = 1./sqrt(nSubsamples);
                            }

                            // DO BOOTSTRAP NN
                            if (1)
                            {
                                if(1)cout << "  NN: processing cBin " << cBin << "..." << endl;
                                if ( !doBS_in_only_1_cBin || ( doBS_in_only_1_cBin && cBin == 0 ) )
                                    w->performSubsampling( 0, kSubsamplingType, nSubsamples );

                                //from centr
                                gr_BS_NN[cW][etaW][phiW][ptW].SetPointsErrorsFromBootstrapHistos( centr, w->histos_BS_NN
                                                                                                  , doBS_in_only_1_cBin, &w->corrInfo_NN, modFactorForSubsamplingError );

                                //from multF
                                gr_BS_NN_fromMultF[cW][etaW][phiW][ptW].SetPointsErrorsFromBootstrapHistos( multF, w->histos_BS_NN
                                                                                                            , doBS_in_only_1_cBin, &w->corrInfo_NN, modFactorForSubsamplingError );

                                w->histos_BS_NN.WriteHistos();

                                //do mixing
                                if ( doEventMixing )
                                {
                                    if(1)cout << "  NN: do event mixing for cBin " << cBin << "..." << endl;
                                    CorrCoeffInfo corrInfoEvMix = w->performEventMixing( 0 );
                                    int nPoints = gr_bCorr_MIX_NN[cW][etaW][phiW][ptW]->GetN();
                                    gr_bCorr_MIX_NN[cW][etaW][phiW][ptW]->SetPoint(nPoints, centr, corrInfoEvMix.bCorr );
                                }
                            }

                            // DO BOOTSTRAP PtPt
                            if (1)
                            {
                                if(1)cout << "  PtPt: processing cBin " << cBin << "..." << endl;
                                if ( !doBS_in_only_1_cBin || ( doBS_in_only_1_cBin && cBin == 0 ) )
                                    w->performSubsampling( 1, kSubsamplingType, nSubsamples );

                                //from centr
                                gr_BS_PtPt[cW][etaW][phiW][ptW].SetPointsErrorsFromBootstrapHistos( centr, w->histos_BS_PtPt, doBS_in_only_1_cBin, &w->corrInfo_PtPt, modFactorForSubsamplingError );

                                //from multF
                                gr_BS_PtPt_fromMultF[cW][etaW][phiW][ptW].SetPointsErrorsFromBootstrapHistos( multF, w->histos_BS_PtPt, doBS_in_only_1_cBin, &w->corrInfo_PtPt, modFactorForSubsamplingError );

                                //QA drawing for BS PtPt:
                                w->histos_BS_PtPt.hist_bCorr->SetLineColor(kOrange - 5 + cW);

                                if (cBin==0)
                                    w->histos_BS_PtPt.hist_bCorr->DrawCopy();
                                else
                                    w->histos_BS_PtPt.hist_bCorr->DrawCopy("same");

                                // !!! TMP!!! cout bootstrap error:
                                double absError = w->histos_BS_PtPt.hist_bCorr->GetRMS();
                                if ( kSubsamplingType == 1 )
                                    absError *= modFactorForSubsamplingError;

                                if(0)cout << " TMP!!! >>> subsampling value, abs. and stat. error (relative) = "
                                          << w->histos_BS_PtPt.hist_bCorr->GetMean() << " "
                                          << absError << " "
                                          << absError/w->histos_BS_PtPt.hist_bCorr->GetMean() << endl;

                                if(0)cout << nSubsamples << " "
                                          << absError << endl;

                                //                            w->histos_BS_PtPt.hist_bCorr->Write();
                                w->histos_BS_PtPt.WriteHistos();

                                //do mixing
                                if ( doEventMixing )
                                {
                                    if(1)cout << "  PtPt: do event mixing for cBin " << cBin << "..." << endl;
                                    CorrCoeffInfo corrInfoEvMix = w->performEventMixing( 1 );
                                    int nPoints = gr_bCorr_MIX_PtPt[cW][etaW][phiW][ptW]->GetN();
                                    gr_bCorr_MIX_PtPt[cW][etaW][phiW][ptW]->SetPoint(nPoints, centr, corrInfoEvMix.bCorr );
                                }

                            }
                        } // end of cBin
                        // #### write BS graphs
                        //BS NN graph
                        gr_BS_NN[cW][etaW][phiW][ptW].gr_bCorr->SetLineColor(kMagenta);
                        gr_BS_NN[cW][etaW][phiW][ptW].gr_bCorr->SetMarkerColor(kMagenta);
                        gr_BS_NN[cW][etaW][phiW][ptW].gr_bCorr->SetMarkerStyle(24);
                        gr_BS_NN[cW][etaW][phiW][ptW].WriteGraphs();

                        gr_BS_NN_fromMultF[cW][etaW][phiW][ptW].WriteGraphs();

                        //BS PtPt graph
                        gr_BS_PtPt[cW][etaW][phiW][ptW].gr_bCorr->SetLineColor(kMagenta);
                        gr_BS_PtPt[cW][etaW][phiW][ptW].gr_bCorr->SetMarkerColor(kMagenta);
                        gr_BS_PtPt[cW][etaW][phiW][ptW].gr_bCorr->SetMarkerStyle(24);
                        gr_BS_PtPt[cW][etaW][phiW][ptW].WriteGraphs();

                        gr_BS_PtPt_fromMultF[cW][etaW][phiW][ptW].WriteGraphs();

                        //write event mixing graphs:
                        if ( doEventMixing && doBootstrap )
                        {
                            gr_bCorr_MIX_NN[cW][etaW][phiW][ptW]->Write();
                            gr_bCorr_MIX_PtPt[cW][etaW][phiW][ptW]->Write();
                        }

                    } // end of phiW
                } // end of etaW
            } // end of ptW
        } // end of cW
        
        //NN
        TCanvas *canv_GrCoeff_BS_NN = new TCanvas("canv_GrCoeff_BS_NN","canv_GrCoeff_BS_NN",80,150,900,700 );
        tuneCanvas(canv_GrCoeff_BS_NN);
        tuneGraphAxisLabels(gr_BS_NN[0][0][0][0].gr_bCorr);
        gr_BS_NN[0][0][0][0].gr_bCorr->DrawClone("APL");
        
        
        //PtPt
        TCanvas *canv_GrCoeff_BS_PtPt = new TCanvas("canv_GrCoeff_BS_PtPt","canv_GrCoeff_BS_PtPt",80,150,900,700 );
        tuneCanvas(canv_GrCoeff_BS_PtPt);
        //        grC2->Draw("APL");
        
        //        grFromFit2D->SetLineColor(kRed);
        //        grFromFit2D->DrawClone("PL");
        
        tuneGraphAxisLabels(gr_BS_PtPt[0][0][0][0].gr_bCorr);
        gr_BS_PtPt[0][0][0][0].gr_bCorr->DrawClone("APL");
    } // end of bootstrap
    
    
    // ########### draw graphs ###########
    int colors[] = { kBlack, kBlack, kBlue, kRed, kMagenta };
    int markers[] = { 20, 24, 5, 2, 21 };
    
    const int ptId = 0;
    const int etaId = 0;
    const int phiId = 0;
    

    // ########## NN
    TCanvas *canv_grNN = new TCanvas("canv_grNN","canv_grNN",20,50,700,600 );
    tuneCanvas(canv_grNN);
    grNN[0][etaId][phiId][ptId].gr_bCorr->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels( grNN[0][etaId][phiId][ptId].gr_bCorr );
    
    
    //centr
    for ( int cW = 0; cW < nCW; cW++ )
    {
        int marker = (cW < 5 ? markers[cW] : 20 );
        int color = (cW < 5 ? colors[cW] : kGray+1 );
        drawGraph(grNN[cW][etaId][phiId][ptId].gr_bCorr, marker, color, cW == 0 ? "AP" : "P");
    }
    
    if ( doBootstrap )
    {
        gr_BS_NN[0][0][0][0].gr_bCorr->Draw("P");
    }

    
    grNN[0][etaId][phiId][ptId].gr_bCorr->SetMinimum( 0 );
    
    TLegend *leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grNN[cW][etaId][phiId][ptId].gr_bCorr, Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();
    
    TLatex *tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);
    
    tex = new TLatex(0.4,0.89, "NN");
    drawTex(tex, 0.045);
    
    
    TString strPostfix;
    
    if ( flag_V0M_ZDC_CL1 == 0 )
        strPostfix = Form("V0M.eps");
    else if ( flag_V0M_ZDC_CL1 == 1 )
        strPostfix = Form("ZDCZEM.eps");
    else if ( flag_V0M_ZDC_CL1 == 2 )
        strPostfix = Form("CL1.eps");
    
    canv_grNN->SaveAs( Form("output/NN_%s", strPostfix.Data() ) );
    
    // ########## PtPt
    TCanvas *canv_grPtPt = new TCanvas("canv_grPtPt","canv_grPtPt",250,50,700,600 );
    tuneCanvas(canv_grPtPt);
    grPtPt[0][etaId][phiId][ptId].gr_bCorr->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels( grPtPt[0][etaId][phiId][ptId].gr_bCorr );
    
    
    //centr
    for ( int cW = 0; cW < nCW; cW++ )
    {
        int marker = (cW < 5 ? markers[cW] : 20 );
        int color = (cW < 5 ? colors[cW] : kGray+1 );
        drawGraph(grPtPt[cW][etaId][phiId][ptId].gr_bCorr, marker, color, cW == 0 ? "AP" : "P");
    }
    
    if ( doBootstrap )
    {
        gr_BS_PtPt[0][0][0][0].gr_bCorr->Draw("P");
    }
    
    grPtPt[0][etaId][phiId][ptId].gr_bCorr->SetMinimum( 0 );
    
    leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grPtPt[cW][etaId][phiId][ptId].gr_bCorr, Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();
    
    tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);
    
    tex = new TLatex(0.4,0.89, "PtPt");
    drawTex(tex, 0.045);
    
    
    canv_grNN->SaveAs( Form("output/PtPt_%s", strPostfix.Data() ) );
    
    
    
    // ########## PtN
    TCanvas *canv_grPtN = new TCanvas("canv_grPtN","canv_grPtN",450,50,700,600 );
    tuneCanvas(canv_grPtN);
    grPtN[0][etaId][phiId][ptId].gr_bCorr->SetTitle( ";centrality percentile;b_{corr}" ); //C_{2}
    tuneGraphAxisLabels( grPtN[0][etaId][phiId][ptId].gr_bCorr );
    
    
    //centr 10
    for ( int cW = 0; cW < nCW; cW++ )
    {
        int marker = (cW < 5 ? markers[cW] : 20 );
        int color = (cW < 5 ? colors[cW] : kGray+1 );
        drawGraph(grPtN[cW][etaId][phiId][ptId].gr_bCorr, marker, color, cW == 0 ? "AP" : "P");
    }
    
    grPtN[0][etaId][phiId][ptId].gr_bCorr->SetMinimum( 0 );
    leg = new TLegend(0.65,0.65,0.999,0.95);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grPtN[cW][etaId][phiId][ptId].gr_bCorr, Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();
    
    tex = new TLatex(0.7,0.4, "#eta_{gap}=0.8, #delta#eta=0.4");
    drawTex(tex, 0.045);
    
    tex = new TLatex(0.4,0.89, "PtN");
    drawTex(tex, 0.045);
    
    
    canv_grPtN->SaveAs( Form("output/PtN_%s", strPostfix.Data() ) );
    
    
    
    TCanvas *canv_grPtN_2D = new TCanvas("canv_grPtN_2D","canv_grPtN_2D",450,50,700,600 );
    tuneCanvas(canv_grPtN_2D);
    
    
    //    wins[0][0][0][0].hist2D_PtN->DrawCopy();
    if (doInitHistos)
        wins[0][0][0][0][0].hist2D_PtN->ProfileX()->DrawCopy();

    // ########## CENTR ESTIMATOR EVENT RATIO:
    TCanvas *canv_grCentrRatio = new TCanvas("canv_grCentrRatio","canv_grCentrRatio",450,250,700,600 );
    tuneCanvas(canv_grCentrRatio);
    grFractEstByV0M[0]->SetTitle(";centrality percentile;ratio");
    tuneGraphAxisLabels(grFractEstByV0M[0]);
    
    //centr
    for ( int cW = 0; cW < nCW; cW++ )
    {
        int marker = (cW < 5 ? markers[cW] : 20 );
        int color = (cW < 5 ? colors[cW] : kGray+1 );
        drawGraph(grFractEstByV0M[cW], marker, color, cW == 0 ? "AP" : "P");
    }
    
    leg = new TLegend(0.65,0.65,0.999,0.95,"ratio #frac{V0M-and-ZEMZDC}{V0M}");
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    
    for ( int cW = 0; cW < nCW; cW++ )
        leg->AddEntry(grFractEstByV0M[cW], Form("class width %.1f", cWidths[cW]), "p");
    leg->Draw();
    //    leg->SetHeader("ratio #frac{V0M-and-ZEMZDC}/{V0M}");
    
    canv_grCentrRatio->SaveAs("output/ratio_V0M-and-ZEMZDC_by_V0M.eps");
    
    
    if(doInitHistos)
    {
        // CENTR ESTIMATOR PERCENTILES QA:
        TCanvas *canv_estimatorPercentiles_QA = new TCanvas("canv_estimatorPercentiles_QA","canv_estimatorPercentiles_QA",50,350,700,600 );
        for ( int cBin = 0; cBin < nCentrBins[0]; cBin++ )
        {
            TH1D *h = wins[0][cBin][0][0][0].hist1D_EstimatorEntries;
            h->SetLineColor(kOrange-9+cBin);
            if ( cBin == 0 )
                h->DrawCopy();
            else
                h->DrawCopy("same");
        }

        // MULT F IN CENTR CLASSES QA:
        TCanvas *canv_mult_F_in_centr_QA = new TCanvas("canv_mult_F_in_centr_QA","canv_mult_F_in_centr_QA",50,400,700,600 );
        gPad->SetLogy();
        for ( int cBin = 0; cBin < nCentrBins[0]; cBin++ )
        {
            TH1D *h = wins[0][cBin][0][0][0].hist1D_multDistrF;
            h->SetLineColor(kOrange-9+cBin);

            if ( cBin == 0 )
            {
                h->SetLineColor(kRed);
                h->GetYaxis()->SetRangeUser(1,100000);
            }

            if ( cBin == 0 )
                h->DrawCopy();
            else
                h->DrawCopy("same");
        }


        // Check entries in centrality bins:
        for ( int cW = 0; cW < nCW; cW++ )
        {
            cout << " ###### cW = " << cW << endl;
            for ( int cBin = 0; cBin < nCentrBins[cW]; cBin++ )
            {
                TH1D *h = wins[cW][cBin][0][0][0].hist1D_multDistrF;
                cout << "entries in cBin " << cBin << " = " << h->GetEntries() << endl;
            }
        }
    }
    
    
    
    // ###### write graphs to output file
    for ( int cW = 0; cW < nCW; cW++ )
        for ( int ptW = 0; ptW < nPtBr*nPtBr; ptW++ )
            for ( int etaW = 0; etaW < nEtaWins; etaW++ )
                for ( int phiW = 0; phiW < nPhiWins*nPhiWins; phiW++ )
                {
                    grNN[cW][etaW][phiW][ptW].WriteGraphs();
                    grPtPt[cW][etaW][phiW][ptW].WriteGraphs();
                    grPtN[cW][etaW][phiW][ptW].WriteGraphs();

                    grNN_fromMultF[cW][etaW][phiW][ptW].WriteGraphs();
                    grPtPt_fromMultF[cW][etaW][phiW][ptW].WriteGraphs();
                    grPtN_fromMultF[cW][etaW][phiW][ptW].WriteGraphs();

                    if(useSameWinVarForCorrectionNN)
                        grNNTmpForCorrNN[cW][etaW][phiW][ptW].WriteGraphs();
                }
    
    canv_grCentrRatio->Write();
    
    fileOutput->Close();
    
    
    // WRITE E-BY-E GRAPHS TO FILE
    if(0)
    {
        TFile *fileGrEbyE_avPt = new TFile( "out_grEbyE_avPt.root", "RECREATE" );
        for (int i=0; i < nGraphsEbyE; i++)
        {
            //        grEbyE_avPtB_inPhiBins[i]->SetName( Form( "grEbyE_avPtB_ev%d", i ) );
            grEbyE_nB_inPhiBins[i]->Write();
            grEbyE_avPtB_inPhiBins[i]->Write();
        }
        
        fileGrEbyE_avPt->Close();
    }
    
    //    for (Int_t i=0; i<nFiles; i++)
    //        myFiles[i]->Close();
    
    //            gROOT->ProcessLine( ".q");
    
    return;
}



