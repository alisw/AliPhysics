/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskCorPIDTOFdeut
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
//#include <iostream>
//#include <fstream>
#include <cmath>
//#include <iomanip>
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TMath.h"
//#include "TGraphErrors.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCorPIDTOFdeut.h"
#include "AliPIDResponse.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
//#include "AliEmcalTrackSelection.h"
//#include "AliEmcalTrackSelectionAOD.h"

//#include "TFile.h"

//#include "AliMultSelection.h"


// particle identifications in ALICE
//
//    enum EParticleType
//
//    kElectron =  0;
//    kMuon     =  1;
//    kPion     =  2;
//    kKaon     =  3;
//    kProton   =  4;
//    kDeuteron =  5;
//    kTriton   =  6;
//    kHe3      =  7;
//    kAlpha    =  8;
//    kPhoton   =  9;
//    kPi0      = 10;
//    kNeutron  = 11;
//    kKaon0    = 12;
//    kEleCon   = 13;
//    kUnknown  = 14;




//class AliAnalysisTaskCorPIDTOFdeut;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

//using namespace BSchaefer_devel;

//ofstream fout;
//Int_t NTracks            = 0;
//Int_t NGoodTracks        = 0;
//Int_t NAssociatedTracks  = 0;
//Int_t DeuteronCandidates = 0;


// SetCutGeoNcrNcl  cut out 3 cm on the sides of the TPC



//Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]
//TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);

ClassImp(AliAnalysisTaskCorPIDTOFdeut) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFdeut::AliAnalysisTaskCorPIDTOFdeut() : AliAnalysisTaskSE(), 
fAOD(0), fOutputList(0), fPIDResponse(0),

    fHistPt(0),                 //  1
    cent_ntracks(0),            //  2
    
    m2_pos(0),                  //  3
    m2_neg(0),                  //  4
    beta_p_pos(0),              //  5
    beta_p_neg(0),              //  6
    deltat_p_pos(0),            //  7
    deltat_p_neg(0),            //  8
    m2_pos_cut(0),              //  9
    m2_neg_cut(0),              // 10
    beta_p_pos_cut(0),          // 11
    beta_p_neg_cut(0),          // 12
    deltat_p_pos_cut(0),        // 13
    deltat_p_neg_cut(0),        // 14
    m2_pos_cut_T(0),            // 15
    m2_neg_cut_T(0),            // 16
    deut_per_event(0),          // 17

    dphi_pt_deut_C(0),          // 18
    dphi_kt_deut_C(0),          // 19
    dphi_et_deut_C(0),          // 20
    dphi_en_deut_C(0),          // 21

    deut_pt_count_C(0),         // 22
    deut_kt_count_C(0),         // 23
    deut_et_count_C(0),         // 24
    deut_en_count_C(0),         // 25

    dphi_pt_deut_P(0),          // 26
    dphi_kt_deut_P(0),          // 27
    dphi_et_deut_P(0),          // 28
    dphi_en_deut_P(0),          // 29

    deut_pt_count_P(0),         // 30
    deut_kt_count_P(0),         // 31
    deut_et_count_P(0),         // 32
    deut_en_count_P(0),         // 33
    
    centrality_dist(0),         // 34
    multiplicity_dist(0),       // 35
    centrality_raw(0),          // 36
    multiplicity_raw(0),        // 37

    deut_dphi_deta_p1_0510(0),  // 38
    deut_dphi_deta_p1_1020(0),  // 39
    deut_dphi_deta_p1_2030(0),  // 40
    deut_dphi_deta_p1_3040(0),  // 41
    deut_dphi_deta_p1_4050(0),  // 42

    deut_dphi_deta_p2_0510(0),  // 43
    deut_dphi_deta_p2_1020(0),  // 44
    deut_dphi_deta_p2_2030(0),  // 45
    deut_dphi_deta_p2_3040(0),  // 46
    deut_dphi_deta_p2_4050(0),  // 47
    
    deut_dphi_deta_p3_0510(0),  // 48
    deut_dphi_deta_p3_1020(0),  // 49
    deut_dphi_deta_p3_2030(0),  // 50
    deut_dphi_deta_p3_3040(0),  // 51
    deut_dphi_deta_p3_4050(0),  // 52

    deut_cor_radius(0)          // 53
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFdeut::AliAnalysisTaskCorPIDTOFdeut(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0),

    fHistPt(0),                 //  1
    cent_ntracks(0),            //  2
    
    m2_pos(0),                  //  3
    m2_neg(0),                  //  4
    beta_p_pos(0),              //  5
    beta_p_neg(0),              //  6
    deltat_p_pos(0),            //  7
    deltat_p_neg(0),            //  8
    m2_pos_cut(0),              //  9
    m2_neg_cut(0),              // 10
    beta_p_pos_cut(0),          // 11
    beta_p_neg_cut(0),          // 12
    deltat_p_pos_cut(0),        // 13
    deltat_p_neg_cut(0),        // 14
    m2_pos_cut_T(0),            // 15
    m2_neg_cut_T(0),            // 16
    deut_per_event(0),          // 17

    dphi_pt_deut_C(0),          // 18
    dphi_kt_deut_C(0),          // 19
    dphi_et_deut_C(0),          // 20
    dphi_en_deut_C(0),          // 21

    deut_pt_count_C(0),         // 22
    deut_kt_count_C(0),         // 23
    deut_et_count_C(0),         // 24
    deut_en_count_C(0),         // 25

    dphi_pt_deut_P(0),          // 26
    dphi_kt_deut_P(0),          // 27
    dphi_et_deut_P(0),          // 28
    dphi_en_deut_P(0),          // 29

    deut_pt_count_P(0),         // 30
    deut_kt_count_P(0),         // 31
    deut_et_count_P(0),         // 32
    deut_en_count_P(0),         // 33									       

    centrality_dist(0),         // 34
    multiplicity_dist(0),       // 35
    centrality_raw(0),          // 36
    multiplicity_raw(0),        // 37

    deut_dphi_deta_p1_0510(0),  // 38
    deut_dphi_deta_p1_1020(0),  // 39
    deut_dphi_deta_p1_2030(0),  // 40
    deut_dphi_deta_p1_3040(0),  // 41
    deut_dphi_deta_p1_4050(0),  // 42

    deut_dphi_deta_p2_0510(0),  // 43
    deut_dphi_deta_p2_1020(0),  // 44
    deut_dphi_deta_p2_2030(0),  // 45
    deut_dphi_deta_p2_3040(0),  // 46
    deut_dphi_deta_p2_4050(0),  // 47
    
    deut_dphi_deta_p3_0510(0),  // 48
    deut_dphi_deta_p3_1020(0),  // 49
    deut_dphi_deta_p3_2030(0),  // 50
    deut_dphi_deta_p3_3040(0),  // 51
    deut_dphi_deta_p3_4050(0),  // 52

    deut_cor_radius(0)          // 53
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFdeut::~AliAnalysisTaskCorPIDTOFdeut()
{
    // destructor
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdeut::UserCreateOutputObjects()
{

// old //   3.43396 0.0133493 0.0593701 -0.125677 0.0646736 0.148637 3.59796 -0.00619808 -0.107274 0.034074 0.0359503 0.0131823 
//3.09901 0.0446792 0.469684 0.0973401 0.0424328 -0.0685078 2.91151 0.0730297 0.664744 -0.115955 0.0709238 0.15557 


    deut_curves[0][0][0] = 3.09901;     // pos deut mean curve
    deut_curves[0][0][1] = 0.0446792;
    deut_curves[0][0][2] = 0.469684;

    deut_curves[0][1][0] = 0.0973401;  // pos deut sigma curve
    deut_curves[0][1][1] = 0.0424328;
    deut_curves[0][1][2] = -0.0685078;

    deut_curves[1][0][0] = 2.91151;     // neg deut mean curve
    deut_curves[1][0][1] = 0.0730297;
    deut_curves[1][0][2] = -0.664744;

    deut_curves[1][1][0] = -0.115955;  // neg deut sigma curve
    deut_curves[1][1][1] = 0.0709238;
    deut_curves[1][1][2] = 0.15557;

    
//    pi = TMath::Pi();
//    fout.open("output.txt");


//  Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]

//  TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);
//    fit_deut_curve->SetParNames("a", "b x", "c/#sqrt(x)");
	
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    fHistPt                = new TH1F("fHistPt",                "Pt()",                 1000,     0,     10);                              //  1
    cent_ntracks           = new TH2F("cent_ntracks",           "cent_ntracks",          100,     0,    100,     100,       0,     800);   //  2
    
    m2_pos                 = new TH2F("m2_pos",                 "m2_pos",                320,   0.0,    8.0,    2400,    -1.0,     7.0);   //  3
    m2_neg                 = new TH2F("m2_neg",                 "m2_neg",                320,   0.0,    8.0,    2400,    -1.0,     7.0);   //  4
    beta_p_pos             = new TH2F("beta_p_pos",             "beta_p_pos",            780,   0.0,    8.0,    3000,     0.1,     1.1);   //  5
    beta_p_neg             = new TH2F("beta_p_neg",             "beta_p_neg",            780,   0.0,    8.0,    3000,     0.1,     1.1);   //  6
    deltat_p_pos           = new TH2F("deltat_p_pos",           "deltat_p_pos",        27000,   0.3,    3.0,    1000,    -1.0,     9.0);   //  7
    deltat_p_neg           = new TH2F("deltat_p_neg",           "deltat_p_neg",        27000,   0.3,    3.0,    1000,    -1.0,     9.0);   //  8

    m2_pos_cut             = new TH2F("m2_pos_cut",             "m2_pos_cut",            320,   0.0,    8.0,    2400,    -1.0,     7.0);   //  9
    m2_neg_cut             = new TH2F("m2_neg_cut",             "m2_neg_cut",            320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 10
    beta_p_pos_cut         = new TH2F("beta_p_pos_cut",         "beta_p_pos_cut",        780,   0.0,    8.0,    3000,     0.1,     1.1);   // 11
    beta_p_neg_cut         = new TH2F("beta_p_neg_cut",         "beta_p_neg_cut",        780,   0.0,    8.0,    3000,     0.1,     1.1);   // 12
    deltat_p_pos_cut       = new TH2F("deltat_p_pos_cut",       "deltat_p_pos_cut",    27000,   0.3,    3.0,    1000,    -1.0,     9.0);   // 13
    deltat_p_neg_cut       = new TH2F("deltat_p_neg_cut",       "deltat_p_neg_cut",    27000,   0.3,    3.0,    1000,    -1.0,     9.0);   // 14
    m2_pos_cut_T           = new TH2F("m2_pos_cut_T",           "m2_pos_cut_T",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 15
    m2_neg_cut_T           = new TH2F("m2_neg_cut_T",           "m2_neg_cut_T",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 16
    deut_per_event         = new TH1I("deut_per_event",         "deut_per_event",          5,     0,     5);                               // 17
    
    dphi_pt_deut_C         = new TH2F("dphi_pt_deut_C",         "dphi_pt_deut_C",        160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 18 // dphi, pt
    dphi_kt_deut_C         = new TH2F("dphi_kt_deut_C",         "dphi_kt_deut_C",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 19 // dphi, kt
    dphi_et_deut_C         = new TH2F("dphi_et_deut_C",         "dphi_et_deut_C",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 20 // dphi, et
    dphi_en_deut_C         = new TH2F("dphi_en_deut_C",         "dphi_en_deut_C",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 21 // dphi, en

    deut_pt_count_C        = new TH1F("deut_pt_count_C",        "deut_pt_count_C",       160,  1.00,   5.00);                              // 22
    deut_kt_count_C        = new TH1F("deut_kt_count_C",        "deut_kt_count_C",       239,  1.00,   5.78);                              // 23
    deut_et_count_C        = new TH1F("deut_et_count_C",        "deut_et_count_C",       239,  1.00,   5.78);                              // 24
    deut_en_count_C        = new TH1F("deut_en_count_C",        "deut_en_count_C",       239,  1.00,   5.78);                              // 25

    dphi_pt_deut_P         = new TH2F("dphi_pt_deut_P",         "dphi_pt_deut_P",        160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 26 // dphi, pt
    dphi_kt_deut_P         = new TH2F("dphi_kt_deut_P",         "dphi_kt_deut_P",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 27 // dphi, kt
    dphi_et_deut_P         = new TH2F("dphi_et_deut_P",         "dphi_et_deut_P",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 28 // dphi, et
    dphi_en_deut_P         = new TH2F("dphi_en_deut_P",         "dphi_en_deut_P",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 29 // dphi, en

    deut_pt_count_P        = new TH1F("deut_pt_count_P",        "deut_pt_count_P",       160,  1.00,   5.00);                              // 30
    deut_kt_count_P        = new TH1F("deut_kt_count_P",        "deut_kt_count_P",       239,  1.00,   5.78);                              // 31
    deut_et_count_P        = new TH1F("deut_et_count_P",        "deut_et_count_P",       239,  1.00,   5.78);                              // 32
    deut_en_count_P        = new TH1F("deut_en_count_P",        "deut_en_count_P",       239,  1.00,   5.78);                              // 33
    
    centrality_dist        = new TH1F("centrality_dist",        "centrality_dist",       100,     0,    100);                              // 34
    multiplicity_dist      = new TH1F("multiplicity_dist",      "multiplicity_dist",     100,     0,    800);                              // 35
    centrality_raw         = new TH1F("centrality_raw",         "centrality_raw",        100,     0,    100);                              // 36
    multiplicity_raw       = new TH1F("multiplicity_raw",       "multiplicity_raw",      100,     0,    800);                              // 37
    
    deut_dphi_deta_p1_0510 = new TH2F("deut_dphi_deta_p1_0510", "deut_dphi_deta_p1_0510", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 38  // dphi    ,    deta
    deut_dphi_deta_p1_1020 = new TH2F("deut_dphi_deta_p1_1020", "deut_dphi_deta_p1_1020", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 39
    deut_dphi_deta_p1_2030 = new TH2F("deut_dphi_deta_p1_2030", "deut_dphi_deta_p1_2030", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 40
    deut_dphi_deta_p1_3040 = new TH2F("deut_dphi_deta_p1_3040", "deut_dphi_deta_p1_3040", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 41
    deut_dphi_deta_p1_4050 = new TH2F("deut_dphi_deta_p1_4050", "deut_dphi_deta_p1_4050", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 42

    deut_dphi_deta_p2_0510 = new TH2F("deut_dphi_deta_p2_0510", "deut_dphi_deta_p2_0510", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 43// dphi    ,    deta
    deut_dphi_deta_p2_1020 = new TH2F("deut_dphi_deta_p2_1020", "deut_dphi_deta_p2_1020", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 44
    deut_dphi_deta_p2_2030 = new TH2F("deut_dphi_deta_p2_2030", "deut_dphi_deta_p2_2030", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 45
    deut_dphi_deta_p2_3040 = new TH2F("deut_dphi_deta_p2_3040", "deut_dphi_deta_p2_3040", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 46
    deut_dphi_deta_p2_4050 = new TH2F("deut_dphi_deta_p2_4050", "deut_dphi_deta_p2_4050", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 47
        
    deut_dphi_deta_p3_0510 = new TH2F("deut_dphi_deta_p3_0510", "deut_dphi_deta_p3_0510", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 48// dphi    ,    deta
    deut_dphi_deta_p3_1020 = new TH2F("deut_dphi_deta_p3_1020", "deut_dphi_deta_p3_1020", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 49
    deut_dphi_deta_p3_2030 = new TH2F("deut_dphi_deta_p3_2030", "deut_dphi_deta_p3_2030", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 50
    deut_dphi_deta_p3_3040 = new TH2F("deut_dphi_deta_p3_3040", "deut_dphi_deta_p3_3040", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 51
    deut_dphi_deta_p3_4050 = new TH2F("deut_dphi_deta_p3_4050", "deut_dphi_deta_p3_4050", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 52

    deut_cor_radius        = new TH2F("deut_cor_radius",        "deut_cor_radius",       160,  1.00,   5.00,     325,   -3.53,    3.53);   // 53

	
    fOutputList->Add(fHistPt);                //  1                                                             // objects added to output file
    fOutputList->Add(cent_ntracks);           //  2
    
    fOutputList->Add(m2_pos);                 //  3
    fOutputList->Add(m2_neg);                 //  4 
    fOutputList->Add(beta_p_pos);             //  5
    fOutputList->Add(beta_p_neg);             //  6
    fOutputList->Add(deltat_p_pos);           //  7
    fOutputList->Add(deltat_p_neg);           //  8

    fOutputList->Add(m2_pos_cut);             //  9
    fOutputList->Add(m2_neg_cut);             // 10
    fOutputList->Add(beta_p_pos_cut);         // 11
    fOutputList->Add(beta_p_neg_cut);         // 12
    fOutputList->Add(deltat_p_pos_cut);       // 13
    fOutputList->Add(deltat_p_neg_cut);       // 14
    fOutputList->Add(m2_pos_cut_T);           // 15
    fOutputList->Add(m2_neg_cut_T);           // 16
    fOutputList->Add(deut_per_event);         // 17
    
    fOutputList->Add(dphi_pt_deut_C);         // 18
    fOutputList->Add(dphi_kt_deut_C);         // 19
    fOutputList->Add(dphi_et_deut_C);         // 20
    fOutputList->Add(dphi_en_deut_C);         // 21
    
    fOutputList->Add(deut_pt_count_C);        // 22
    fOutputList->Add(deut_kt_count_C);        // 23
    fOutputList->Add(deut_et_count_C);        // 24
    fOutputList->Add(deut_en_count_C);        // 25
    
    fOutputList->Add(dphi_pt_deut_P);         // 26
    fOutputList->Add(dphi_kt_deut_P);         // 27
    fOutputList->Add(dphi_et_deut_P);         // 28
    fOutputList->Add(dphi_en_deut_P);         // 29
    
    fOutputList->Add(deut_pt_count_P);        // 30
    fOutputList->Add(deut_kt_count_P);        // 31
    fOutputList->Add(deut_et_count_P);        // 32
    fOutputList->Add(deut_en_count_P);        // 33
    
    fOutputList->Add(centrality_dist);        // 34
    fOutputList->Add(multiplicity_dist);      // 35
    fOutputList->Add(centrality_raw);         // 36
    fOutputList->Add(multiplicity_raw);       // 37
    
    fOutputList->Add(deut_dphi_deta_p1_0510); // 38
    fOutputList->Add(deut_dphi_deta_p1_1020); // 39
    fOutputList->Add(deut_dphi_deta_p1_2030); // 40
    fOutputList->Add(deut_dphi_deta_p1_3040); // 41
    fOutputList->Add(deut_dphi_deta_p1_4050); // 42
    
    fOutputList->Add(deut_dphi_deta_p2_0510); // 43
    fOutputList->Add(deut_dphi_deta_p2_1020); // 44
    fOutputList->Add(deut_dphi_deta_p2_2030); // 45
    fOutputList->Add(deut_dphi_deta_p2_3040); // 46
    fOutputList->Add(deut_dphi_deta_p2_4050); // 47

    fOutputList->Add(deut_dphi_deta_p3_0510); // 48
    fOutputList->Add(deut_dphi_deta_p3_1020); // 49
    fOutputList->Add(deut_dphi_deta_p3_2030); // 50
    fOutputList->Add(deut_dphi_deta_p3_3040); // 51
    fOutputList->Add(deut_dphi_deta_p3_4050); // 52

    fOutputList->Add(deut_cor_radius);        // 53
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdeut::UserExec(Option_t *)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    Int_t iTracks(fAOD->GetNumberOfTracks());


    Float_t centrality = 1.0;
//    centrality = fAOD->GetCentrality()->GetCentralityPercentile("VOM");
//  cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("VOM"), iTracks);

        
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    centrality = fAOD->GetCentrality()->GetCentralityPercentile("V0M");

    centrality_raw->Fill(centrality);
    multiplicity_raw->Fill(iTracks);

    
    int NDeut              = 0;
    int NDeut_pos          = 0;
    int NDeut_neg          = 0;
    int deut_flag_T        = 0;
    int deut_track_num_T   = -9999;;
    int associated_tracks  = 0;



	
    struct track_node        // linked list, for efficient, dynamically allocated memory use
    {
	int track_num;
	track_node *next;
    };
    track_node *root;        // don't change this one, or you will lose the list in memory
    track_node *conductor;   // points to each node as it traverses the list
    root = new track_node;
    root->next = NULL;
    conductor = root;



    // loop over all these tracks
    for(Int_t i(0); i < iTracks; i++)
    {

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                                                       // if we failed, skip this track
	if(!(track->IsHybridGlobalConstrainedGlobal())){   continue;   }

//	if(!track->IsPrimaryCandidate())   continue;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	AliExternalTrackParam trkparam;
	trkparam.CopyFromVTrack(track);
	
	track_node *new_track = new track_node;
	new_track->track_num = i;
	new_track->next = NULL;
	conductor->next = new_track;
	conductor = conductor->next;

	
	short isGlobal    = 1;
	short isModGlobal = 0;  // is it suitable for deuteron candidate
	
	if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal))){   isGlobal = 0;   }
	
	Double_t nsigmaTPC = 999.0;
	Double_t nsigmaTOF = 999.0;

	if(isGlobal == 1)  
	{
	    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	    Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);     /* && trk->IsOn(AliESDtrack::kTPCpid)*/;
	    Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);

	    if(tpcIsOk  && tofIsOk){   isModGlobal = 1;   }
	    else                   {   isModGlobal = 0;   }
		
	}

	    
	Double_t mom            = track->P();
	float    pt             = track->Pt();
	Short_t  charge         = track->Charge();

	Double_t sigma_min      = -999.0;
	int      id             = 14;
	Double_t deut_mean      = 0.0;
	Double_t deut_sigma     = 0.0;
	Double_t cut_width      = 1.0;

	Double_t phi = track->Phi();
	Double_t eta = track->Eta();


	if(fabs(eta > 0.8))   continue;

	   
	if(phi <  -TMath::PiOver2())    phi = phi + TMath::TwoPi();			if(phi <  -TMath::PiOver2())    phi = phi + TMath::TwoPi();
	if(phi > 3*TMath::PiOver2())    phi = phi - TMath::TwoPi();			if(phi > 3*TMath::PiOver2())    phi = phi - TMath::TwoPi();


	fHistPt->Fill(pt);


	// markc
	if(isModGlobal == 1)
	{
	    Double_t m2tof = get_mass_squared(track);
	    Float_t deltat = tof_minus_tpion(track);

//    kElectron =  0;
//    kMuon     =  1;
//    kPion     =  2;
//    kKaon     =  3;
//    kProton   =  4;
//    kDeuteron =  5;
//    kTriton   =  6;
//    kHe3      =  7;
//    kAlpha    =  8;
//    kPhoton   =  9;
//    kPi0      = 10;
//    kNeutron  = 11;
//    kKaon0    = 12;
//    kEleCon   = 13;
//    kUnknown  = 14;  

	    
      ////  markA

	    Float_t beta_kaon = 0.0;
	    Float_t beta      = 0.0;

	    beta      = Beta(track);
	    beta_kaon = mom / sqrt(0.2437169 + pow(mom,2));  // mass-squared of kaon is 0.2437169
	    
	    if(charge > 0)
	    {
		m2_pos->Fill(        pt, m2tof);
		beta_p_pos->Fill(   mom, Beta(track));
		deltat_p_pos->Fill(  pt, deltat);

//		if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340  &&  isGlobal == 1)
		if(                                                                      isGlobal == 1)
		{
//		    if(    (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
//		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
//		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
//		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
//			   )
//		    {
			AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
			if(nsigmaTPC <= 2.0)
			{
			    m2_pos_cut->Fill(        pt, m2tof);
			    beta_p_pos_cut->Fill(   mom, Beta(track));
			    deltat_p_pos_cut->Fill(  pt, deltat);


			    if(pt >= 1.4  &&  pt < 4.4)
			    {
				for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
				deut_mean = fit_deut_curve->Eval(mom);
				for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
				deut_sigma = fit_deut_curve->Eval(mom);

				if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
				{
				    m2_pos_cut_T->Fill(mom,m2tof);
				    deut_track_num_T = i;
				    NDeut++;
				    NDeut_pos++;
				    deut_flag_T = 1;
				}		    
			    }
			}
		 // }
		}
	    }
	    
	    else if(charge < 0)
	    {
		m2_neg->Fill(        pt, m2tof);
		beta_p_neg->Fill(   mom, Beta(track));
		deltat_p_neg->Fill(  pt, deltat);
		
//		if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340  &&  isGlobal == 1)
		if(                                                                      isGlobal == 1)
		{
//		    if(    (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
//		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
//		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
//		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
//			   )
//		    {
			AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
			if(nsigmaTPC <= 2.0)
			{	
			    m2_neg_cut->Fill(        pt, m2tof);
			    beta_p_neg_cut->Fill(   mom, Beta(track));
			    deltat_p_neg_cut->Fill(  pt, deltat);

			    if(pt >= 1.4  &&  pt < 4.4)
			    {
				for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
				deut_mean = fit_deut_curve->Eval(mom);
				for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
				deut_sigma = fit_deut_curve->Eval(mom);
				
				if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
				{
				    m2_neg_cut_T->Fill(mom,m2tof);
				    deut_track_num_T = i;
				    NDeut++;
				    NDeut_neg++;
				    deut_flag_T = 1;
				}
			    }  
			}
//		    }
		}
	    }
	}    //   track has both TPC and TOF PID, tracks are o.k.
    }        //   end of track loop


    deut_per_event->Fill(NDeut);

    if(deut_flag_T > 0)
    {
//	DeuteronCandidates++;
	int j = deut_track_num_T;

	centrality_dist->Fill(centrality);
	multiplicity_dist->Fill(iTracks);
	
	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t deut_phi = track->Phi();
	    Double_t deut_mom = track->P();
	    Short_t  charge   = track->Charge();
	    Double_t deut_eta = track->Eta();

	    Double_t deut_pt  = track->Pt();
	    Double_t deut_kt  = 0.00;
	    Double_t deut_et  = 0.00;
	    Double_t deut_en  = 0.00;
	    
	    deut_kt = sqrt(pow(deut_pt, 2) + pow(1.87561,2)) - 1.87561;
	    deut_et = sqrt(pow(deut_pt, 2) + pow(1.87561,2));
	    deut_en = sqrt(pow(deut_mom,2) + pow(1.87561,2));

	    
	    int asso = 0;
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;
		    
		    if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
		    {
			AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
			if(!track){}                                                               // if we failed, skip this track
			else
			{
			    Double_t mom            = track->P();
			    Double_t phi            = track->Phi();
			    Double_t pt             = track->Pt();
			    Double_t eta            = track->Eta();
			    
			    if(pt > 2.0  &&  pt < 5.0)
			    {
				Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();

				if(centrality <= 20)
				{
				    dphi_pt_deut_C->Fill(deut_pt, dphi);
				    dphi_kt_deut_C->Fill(deut_kt, dphi);
				    dphi_et_deut_C->Fill(deut_et, dphi);
				    dphi_en_deut_C->Fill(deut_en, dphi);
				}
				else if(centrality > 60)
				{
				    dphi_pt_deut_P->Fill(deut_pt, dphi);
				    dphi_kt_deut_P->Fill(deut_kt, dphi);
				    dphi_et_deut_P->Fill(deut_et, dphi);
				    dphi_en_deut_P->Fill(deut_en, dphi);
				}

				dphi = deut_phi - phi;
				if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
				if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
				if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
				if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
				deut_cor_radius->Fill(deut_pt, TMath::Sign(1.0,dphi)*sqrt(pow(dphi,2) + pow(deut_eta - eta,2)));			    
				
				asso++;
			    }

			    if(deut_pt >= 1.4  &&  deut_pt < 2.4)
			    {
				Double_t eta  = track->Eta();	 Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p1_0510->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1_1020->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p1_2030->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 3.0  &&  mom < 4.0){   deut_dphi_deta_p1_3040->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 4.0  &&  mom < 5.0){   deut_dphi_deta_p1_4050->Fill(dphi, deut_eta - eta);    }
			    }
			    else if(deut_pt >= 2.4  &&  deut_pt < 3.4)
			    {
				Double_t eta  = track->Eta();	 Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p2_0510->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p2_1020->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2_2030->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 3.0  &&  mom < 4.0){   deut_dphi_deta_p2_3040->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 4.0  &&  mom < 5.0){   deut_dphi_deta_p2_4050->Fill(dphi, deut_eta - eta);    }
			    }
			    else if(deut_pt >= 3.4  &&  deut_pt < 4.4)
			    {
				Double_t eta  = track->Eta();	 Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p3_0510->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p3_1020->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p3_2030->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 3.0  &&  mom < 4.0){   deut_dphi_deta_p3_3040->Fill(dphi, deut_eta - eta);    }
				else if(mom >= 4.0  &&  mom < 5.0){   deut_dphi_deta_p3_4050->Fill(dphi, deut_eta - eta);    }
			    }
			}
		    }
		}
	    
		int k = conductor->track_num;
		conductor = conductor->next;
		if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t mom            = track->P();
			Double_t phi            = track->Phi();
			Double_t eta            = track->Eta();
			Double_t pt             = track->Pt();
			
			if(pt > 2.0  &&  pt < 5.0)
			{
			    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();

			    if(centrality <= 20)
			    {
				dphi_pt_deut_C->Fill(deut_pt, dphi);
				dphi_kt_deut_C->Fill(deut_kt, dphi);
				dphi_et_deut_C->Fill(deut_et, dphi);
				dphi_en_deut_C->Fill(deut_en, dphi);
			    }
			    else if(centrality > 60)
			    {
				dphi_pt_deut_P->Fill(deut_pt, dphi);
				dphi_kt_deut_P->Fill(deut_kt, dphi);
				dphi_et_deut_P->Fill(deut_et, dphi);
				dphi_en_deut_P->Fill(deut_en, dphi);
			    }
			    
			    dphi = deut_phi - phi;
			    if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
			    if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
			    if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
			    if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
			    deut_cor_radius->Fill(deut_pt, TMath::Sign(1.0,dphi)*sqrt(pow(dphi,2) + pow(deut_eta - eta,2)));
			    asso++;			    
			}

			if(deut_pt >= 1.4  &&  deut_pt < 2.4)
			{
			    Double_t eta  = track->Eta();    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p1_0510->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1_1020->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p1_2030->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 3.0  &&  mom < 4.0){   deut_dphi_deta_p1_3040->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 4.0  &&  mom < 5.0){   deut_dphi_deta_p1_4050->Fill(dphi, deut_eta - eta);    }
			}
			else if(deut_pt >= 2.4  &&  deut_pt < 3.4)
			{
			    Double_t eta  = track->Eta();    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p2_0510->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p2_1020->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2_2030->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 3.0  &&  mom < 4.0){   deut_dphi_deta_p2_3040->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 4.0  &&  mom < 5.0){   deut_dphi_deta_p2_4050->Fill(dphi, deut_eta - eta);    }
			}
			else if(deut_pt >= 3.4  &&  deut_pt < 4.4)
			{
			    Double_t eta  = track->Eta();    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();			    
			    if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p3_0510->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p3_1020->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p3_2030->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 3.0  &&  mom < 4.0){   deut_dphi_deta_p3_3040->Fill(dphi, deut_eta - eta);    }
			    else if(mom >= 4.0  &&  mom < 5.0){   deut_dphi_deta_p3_4050->Fill(dphi, deut_eta - eta);    }
			}			
		    }
		}
	    }

	    if(asso > 0)
	    {
		if(centrality <= 20.0)
		{
		    deut_pt_count_C->Fill(deut_pt);
		    deut_kt_count_C->Fill(deut_kt);
		    deut_et_count_C->Fill(deut_et);
		    deut_en_count_C->Fill(deut_en);
		}
		else if(centrality > 60)
		{
		    deut_pt_count_P->Fill(deut_pt);
		    deut_kt_count_P->Fill(deut_kt);
		    deut_et_count_P->Fill(deut_et);
		    deut_en_count_P->Fill(deut_en);
		}
	    }
	}
    }  // end of deuteron correlation PRIMARY LOOP

                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file

}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdeut::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFdeut::Beta(AliAODTrack *track)
{
    Double_t startTime     = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());         //in ps
    Double_t stoptime      = track->GetTOFsignal();
    Double_t c             = TMath::C()*1.E-9;                                                              // m/ns
    Double_t length        = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
    stoptime -= startTime;      
    Double_t scaleStopTime = stoptime*1E-3;          
    scaleStopTime          = scaleStopTime*c;
    return length/scaleStopTime;
}



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFdeut::tof_minus_tpion(AliAODTrack *track)
{
    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
    Double_t stop_time      = track->GetTOFsignal();
    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
    Double_t time_of_flight = stop_time - start_time;
    time_of_flight          = time_of_flight * 0.001;                                                       // convert from ps to ns
    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
    length = length * 100;                                                                                  // convert to cm
    Double_t pion_time      = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion);        // // take a close look at this function -BRENNAN SCHAEFER-
    pion_time               = pion_time * 0.001;
    return time_of_flight - pion_time;
}



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFdeut::get_mass_squared(AliAODTrack *track)
{
    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
    Double_t stop_time      = track->GetTOFsignal();
    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
    Double_t tof            = stop_time - start_time;
    tof                     = tof * 0.001;                                                                  // convert from ps to ns
    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
    length                  = length * 100;                                                                 // convert to cm
    Double_t mom            = track->P();
    Double_t m2      = 0.0;

    Double_t c2 = 29.9792458;
    
    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);
    return m2;

}
