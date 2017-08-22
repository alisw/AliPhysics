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

/* AliAnaysisTaskCorPIDTOFhadr
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
#include "AliAnalysisTaskCorPIDTOFhadr.h"
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
//    kHadreron =  5;
//    kTriton   =  6;
//    kHe3      =  7;
//    kAlpha    =  8;
//    kPhoton   =  9;
//    kPi0      = 10;
//    kNeutron  = 11;
//    kKaon0    = 12;
//    kEleCon   = 13;
//    kUnknown  = 14;




//class AliAnalysisTaskCorPIDTOFhadr;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

//using namespace BSchaefer_devel;

//ofstream fout;
//Int_t NTracks            = 0;
//Int_t NGoodTracks        = 0;
//Int_t NAssociatedTracks  = 0;
//Int_t HadreronCandidates = 0;


// SetCutGeoNcrNcl  cut out 3 cm on the sides of the TPC



//Double_t hadr_curves[2][2][3];  // [charge][mean,sigma][par]
//TF1 *fit_hadr_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);

ClassImp(AliAnalysisTaskCorPIDTOFhadr) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFhadr::AliAnalysisTaskCorPIDTOFhadr() : AliAnalysisTaskSE(), 
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
    hadr_per_event(0),          // 17

    dphi_pt_hadr_C(0),          // 18
    dphi_kt_hadr_C(0),          // 19
    dphi_et_hadr_C(0),          // 20
    dphi_en_hadr_C(0),          // 21

    hadr_pt_count_C(0),         // 22
    hadr_kt_count_C(0),         // 23
    hadr_et_count_C(0),         // 24
    hadr_en_count_C(0),         // 25

    dphi_pt_hadr_P(0),          // 26
    dphi_kt_hadr_P(0),          // 27
    dphi_et_hadr_P(0),          // 28
    dphi_en_hadr_P(0),          // 29

    hadr_pt_count_P(0),         // 30
    hadr_kt_count_P(0),         // 31
    hadr_et_count_P(0),         // 32
    hadr_en_count_P(0),         // 33
    
    centrality_dist(0),         // 34
    multiplicity_dist(0),       // 35
    centrality_raw(0),          // 36
    multiplicity_raw(0),        // 37

    hadr_dphi_deta_p1_0510(0),  // 38
    hadr_dphi_deta_p1_1020(0),  // 39
    hadr_dphi_deta_p1_2030(0),  // 40
    hadr_dphi_deta_p1_3040(0),  // 41
    hadr_dphi_deta_p1_4050(0),  // 42

    hadr_dphi_deta_p2_0510(0),  // 43
    hadr_dphi_deta_p2_1020(0),  // 44
    hadr_dphi_deta_p2_2030(0),  // 45
    hadr_dphi_deta_p2_3040(0),  // 46
    hadr_dphi_deta_p2_4050(0),  // 47
    
    hadr_dphi_deta_p3_0510(0),  // 48
    hadr_dphi_deta_p3_1020(0),  // 49
    hadr_dphi_deta_p3_2030(0),  // 50
    hadr_dphi_deta_p3_3040(0),  // 51
    hadr_dphi_deta_p3_4050(0),  // 52

    hadr_cor_radius(0)          // 53
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFhadr::AliAnalysisTaskCorPIDTOFhadr(const char* name) : AliAnalysisTaskSE(name),
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
    hadr_per_event(0),          // 17

    dphi_pt_hadr_C(0),          // 18
    dphi_kt_hadr_C(0),          // 19
    dphi_et_hadr_C(0),          // 20
    dphi_en_hadr_C(0),          // 21

    hadr_pt_count_C(0),         // 22
    hadr_kt_count_C(0),         // 23
    hadr_et_count_C(0),         // 24
    hadr_en_count_C(0),         // 25

    dphi_pt_hadr_P(0),          // 26
    dphi_kt_hadr_P(0),          // 27
    dphi_et_hadr_P(0),          // 28
    dphi_en_hadr_P(0),          // 29

    hadr_pt_count_P(0),         // 30
    hadr_kt_count_P(0),         // 31
    hadr_et_count_P(0),         // 32
    hadr_en_count_P(0),         // 33									       

    centrality_dist(0),         // 34
    multiplicity_dist(0),       // 35
    centrality_raw(0),          // 36
    multiplicity_raw(0),        // 37

    hadr_dphi_deta_p1_0510(0),  // 38
    hadr_dphi_deta_p1_1020(0),  // 39
    hadr_dphi_deta_p1_2030(0),  // 40
    hadr_dphi_deta_p1_3040(0),  // 41
    hadr_dphi_deta_p1_4050(0),  // 42

    hadr_dphi_deta_p2_0510(0),  // 43
    hadr_dphi_deta_p2_1020(0),  // 44
    hadr_dphi_deta_p2_2030(0),  // 45
    hadr_dphi_deta_p2_3040(0),  // 46
    hadr_dphi_deta_p2_4050(0),  // 47
    
    hadr_dphi_deta_p3_0510(0),  // 48
    hadr_dphi_deta_p3_1020(0),  // 49
    hadr_dphi_deta_p3_2030(0),  // 50
    hadr_dphi_deta_p3_3040(0),  // 51
    hadr_dphi_deta_p3_4050(0),  // 52

    hadr_cor_radius(0)          // 53
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFhadr::~AliAnalysisTaskCorPIDTOFhadr()
{
    // destructor
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFhadr::UserCreateOutputObjects()
{

// old //   3.43396 0.0133493 0.0593701 -0.125677 0.0646736 0.148637 3.59796 -0.00619808 -0.107274 0.034074 0.0359503 0.0131823 
//3.09901 0.0446792 0.469684 0.0973401 0.0424328 -0.0685078 2.91151 0.0730297 0.664744 -0.115955 0.0709238 0.15557 


    
//    pi = TMath::Pi();
//    fout.open("output.txt");


//  Double_t hadr_curves[2][2][3];  // [charge][mean,sigma][par]

//  TF1 *fit_hadr_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);
//    fit_hadr_curve->SetParNames("a", "b x", "c/#sqrt(x)");
	
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
    hadr_per_event         = new TH1I("hadr_per_event",         "hadr_per_event",          5,     0,     5);                               // 17
    
    dphi_pt_hadr_C         = new TH2F("dphi_pt_hadr_C",         "dphi_pt_hadr_C",        160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 18 // dphi, pt
    dphi_kt_hadr_C         = new TH2F("dphi_kt_hadr_C",         "dphi_kt_hadr_C",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 19 // dphi, kt
    dphi_et_hadr_C         = new TH2F("dphi_et_hadr_C",         "dphi_et_hadr_C",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 20 // dphi, et
    dphi_en_hadr_C         = new TH2F("dphi_en_hadr_C",         "dphi_en_hadr_C",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 21 // dphi, en

    hadr_pt_count_C        = new TH1F("hadr_pt_count_C",        "hadr_pt_count_C",       160,  1.00,   5.00);                              // 22
    hadr_kt_count_C        = new TH1F("hadr_kt_count_C",        "hadr_kt_count_C",       239,  1.00,   5.78);                              // 23
    hadr_et_count_C        = new TH1F("hadr_et_count_C",        "hadr_et_count_C",       239,  1.00,   5.78);                              // 24
    hadr_en_count_C        = new TH1F("hadr_en_count_C",        "hadr_en_count_C",       239,  1.00,   5.78);                              // 25

    dphi_pt_hadr_P         = new TH2F("dphi_pt_hadr_P",         "dphi_pt_hadr_P",        160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 26 // dphi, pt
    dphi_kt_hadr_P         = new TH2F("dphi_kt_hadr_P",         "dphi_kt_hadr_P",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 27 // dphi, kt
    dphi_et_hadr_P         = new TH2F("dphi_et_hadr_P",         "dphi_et_hadr_P",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 28 // dphi, et
    dphi_en_hadr_P         = new TH2F("dphi_en_hadr_P",         "dphi_en_hadr_P",        239,  1.00,   5.78,     300, -1.6708,  4.8124);   // 29 // dphi, en

    hadr_pt_count_P        = new TH1F("hadr_pt_count_P",        "hadr_pt_count_P",       160,  1.00,   5.00);                              // 30
    hadr_kt_count_P        = new TH1F("hadr_kt_count_P",        "hadr_kt_count_P",       239,  1.00,   5.78);                              // 31
    hadr_et_count_P        = new TH1F("hadr_et_count_P",        "hadr_et_count_P",       239,  1.00,   5.78);                              // 32
    hadr_en_count_P        = new TH1F("hadr_en_count_P",        "hadr_en_count_P",       239,  1.00,   5.78);                              // 33
    
    centrality_dist        = new TH1F("centrality_dist",        "centrality_dist",       100,     0,    100);                              // 34
    multiplicity_dist      = new TH1F("multiplicity_dist",      "multiplicity_dist",     100,     0,    800);                              // 35
    centrality_raw         = new TH1F("centrality_raw",         "centrality_raw",        100,     0,    100);                              // 36
    multiplicity_raw       = new TH1F("multiplicity_raw",       "multiplicity_raw",      100,     0,    800);                              // 37
    
    hadr_dphi_deta_p1_0510 = new TH2F("hadr_dphi_deta_p1_0510", "hadr_dphi_deta_p1_0510", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 38  // dphi    ,    deta
    hadr_dphi_deta_p1_1020 = new TH2F("hadr_dphi_deta_p1_1020", "hadr_dphi_deta_p1_1020", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 39
    hadr_dphi_deta_p1_2030 = new TH2F("hadr_dphi_deta_p1_2030", "hadr_dphi_deta_p1_2030", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 40
    hadr_dphi_deta_p1_3040 = new TH2F("hadr_dphi_deta_p1_3040", "hadr_dphi_deta_p1_3040", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 41
    hadr_dphi_deta_p1_4050 = new TH2F("hadr_dphi_deta_p1_4050", "hadr_dphi_deta_p1_4050", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 42

    hadr_dphi_deta_p2_0510 = new TH2F("hadr_dphi_deta_p2_0510", "hadr_dphi_deta_p2_0510", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 43// dphi    ,    deta
    hadr_dphi_deta_p2_1020 = new TH2F("hadr_dphi_deta_p2_1020", "hadr_dphi_deta_p2_1020", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 44
    hadr_dphi_deta_p2_2030 = new TH2F("hadr_dphi_deta_p2_2030", "hadr_dphi_deta_p2_2030", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 45
    hadr_dphi_deta_p2_3040 = new TH2F("hadr_dphi_deta_p2_3040", "hadr_dphi_deta_p2_3040", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 46
    hadr_dphi_deta_p2_4050 = new TH2F("hadr_dphi_deta_p2_4050", "hadr_dphi_deta_p2_4050", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 47
        
    hadr_dphi_deta_p3_0510 = new TH2F("hadr_dphi_deta_p3_0510", "hadr_dphi_deta_p3_0510", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 48// dphi    ,    deta
    hadr_dphi_deta_p3_1020 = new TH2F("hadr_dphi_deta_p3_1020", "hadr_dphi_deta_p3_1020", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 49
    hadr_dphi_deta_p3_2030 = new TH2F("hadr_dphi_deta_p3_2030", "hadr_dphi_deta_p3_2030", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 50
    hadr_dphi_deta_p3_3040 = new TH2F("hadr_dphi_deta_p3_3040", "hadr_dphi_deta_p3_3040", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 51
    hadr_dphi_deta_p3_4050 = new TH2F("hadr_dphi_deta_p3_4050", "hadr_dphi_deta_p3_4050", 60, -1.6708,  4.8124,   46,   -2.30,    2.30);   // 52

    hadr_cor_radius        = new TH2F("hadr_cor_radius",        "hadr_cor_radius",       160,  1.00,   5.00,     325,   -3.53,    3.53);   // 53

	
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
    fOutputList->Add(hadr_per_event);         // 17
    
    fOutputList->Add(dphi_pt_hadr_C);         // 18
    fOutputList->Add(dphi_kt_hadr_C);         // 19
    fOutputList->Add(dphi_et_hadr_C);         // 20
    fOutputList->Add(dphi_en_hadr_C);         // 21
    
    fOutputList->Add(hadr_pt_count_C);        // 22
    fOutputList->Add(hadr_kt_count_C);        // 23
    fOutputList->Add(hadr_et_count_C);        // 24
    fOutputList->Add(hadr_en_count_C);        // 25
    
    fOutputList->Add(dphi_pt_hadr_P);         // 26
    fOutputList->Add(dphi_kt_hadr_P);         // 27
    fOutputList->Add(dphi_et_hadr_P);         // 28
    fOutputList->Add(dphi_en_hadr_P);         // 29
    
    fOutputList->Add(hadr_pt_count_P);        // 30
    fOutputList->Add(hadr_kt_count_P);        // 31
    fOutputList->Add(hadr_et_count_P);        // 32
    fOutputList->Add(hadr_en_count_P);        // 33
    
    fOutputList->Add(centrality_dist);        // 34
    fOutputList->Add(multiplicity_dist);      // 35
    fOutputList->Add(centrality_raw);         // 36
    fOutputList->Add(multiplicity_raw);       // 37
    
    fOutputList->Add(hadr_dphi_deta_p1_0510); // 38
    fOutputList->Add(hadr_dphi_deta_p1_1020); // 39
    fOutputList->Add(hadr_dphi_deta_p1_2030); // 40
    fOutputList->Add(hadr_dphi_deta_p1_3040); // 41
    fOutputList->Add(hadr_dphi_deta_p1_4050); // 42
    
    fOutputList->Add(hadr_dphi_deta_p2_0510); // 43
    fOutputList->Add(hadr_dphi_deta_p2_1020); // 44
    fOutputList->Add(hadr_dphi_deta_p2_2030); // 45
    fOutputList->Add(hadr_dphi_deta_p2_3040); // 46
    fOutputList->Add(hadr_dphi_deta_p2_4050); // 47

    fOutputList->Add(hadr_dphi_deta_p3_0510); // 48
    fOutputList->Add(hadr_dphi_deta_p3_1020); // 49
    fOutputList->Add(hadr_dphi_deta_p3_2030); // 50
    fOutputList->Add(hadr_dphi_deta_p3_3040); // 51
    fOutputList->Add(hadr_dphi_deta_p3_4050); // 52

    fOutputList->Add(hadr_cor_radius);        // 53
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFhadr::UserExec(Option_t *)
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

    
    int NHadr              = 0;
    int NHadr_pos          = 0;
    int NHadr_neg          = 0;
    int hadr_flag_T        = 0;
    int hadr_track_num_T   = -9999;;
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
	short isModGlobal = 0;  // is it suitable for hadreron candidate
	
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
	Double_t hadr_mean      = 0.0;
	Double_t hadr_sigma     = 0.0;
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
//    kHadreron =  5;
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
				m2_pos_cut_T->Fill(mom,m2tof);
				hadr_track_num_T = i;
				NHadr++;
				NHadr_pos++;
				hadr_flag_T = 1;
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
				m2_neg_cut_T->Fill(mom,m2tof);
				hadr_track_num_T = i;
				NHadr++;
				NHadr_neg++;
				hadr_flag_T = 1;
			    }  
			}
//		    }
		}
	    }
	}    //   track has both TPC and TOF PID, tracks are o.k.
    }        //   end of track loop


    hadr_per_event->Fill(NHadr);

    if(hadr_flag_T > 0)
    {
//	HadreronCandidates++;
	int j = hadr_track_num_T;

	centrality_dist->Fill(centrality);
	multiplicity_dist->Fill(iTracks);
	
	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t hadr_phi = track->Phi();
	    Double_t hadr_mom = track->P();
	    Short_t  charge   = track->Charge();
	    Double_t hadr_eta = track->Eta();

	    Double_t hadr_pt  = track->Pt();
	    Double_t hadr_kt  = 0.00;
	    Double_t hadr_et  = 0.00;
	    Double_t hadr_en  = 0.00;
	    
	    hadr_kt = sqrt(pow(hadr_pt, 2) + pow(1.87561,2)) - 1.87561;
	    hadr_et = sqrt(pow(hadr_pt, 2) + pow(1.87561,2));
	    hadr_en = sqrt(pow(hadr_mom,2) + pow(1.87561,2));

	    
	    int asso = 0;
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;
		    
		    if(k != j)  // if the k-th track is not the j-th track which is the hadreron of interest
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
				Double_t dphi = hadr_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();

				if(centrality <= 20)
				{
				    dphi_pt_hadr_C->Fill(hadr_pt, dphi);	    dphi_kt_hadr_C->Fill(hadr_kt, dphi);
				    dphi_et_hadr_C->Fill(hadr_et, dphi);	    dphi_en_hadr_C->Fill(hadr_en, dphi);
				}
				else if(centrality > 60)
				{
				    dphi_pt_hadr_P->Fill(hadr_pt, dphi);	    dphi_kt_hadr_P->Fill(hadr_kt, dphi);
				    dphi_et_hadr_P->Fill(hadr_et, dphi);	    dphi_en_hadr_P->Fill(hadr_en, dphi);
				}

				dphi = hadr_phi - phi;
				if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
				if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
				if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
				if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
				hadr_cor_radius->Fill(hadr_pt, TMath::Sign(1.0,dphi)*sqrt(pow(dphi,2) + pow(hadr_eta - eta,2)));			    
				
				asso++;
			    }

			    if(hadr_pt >= 1.4  &&  hadr_pt < 2.4)
			    {
				Double_t eta  = track->Eta();	 Double_t dphi = hadr_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(     mom >= 0.5  &&  mom < 1.0){   hadr_dphi_deta_p1_0510->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 1.0  &&  mom < 2.0){   hadr_dphi_deta_p1_1020->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 2.0  &&  mom < 3.0){   hadr_dphi_deta_p1_2030->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 3.0  &&  mom < 4.0){   hadr_dphi_deta_p1_3040->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 4.0  &&  mom < 5.0){   hadr_dphi_deta_p1_4050->Fill(dphi, hadr_eta - eta);    }
			    }
			    else if(hadr_pt >= 2.4  &&  hadr_pt < 3.4)
			    {
				Double_t eta  = track->Eta();	 Double_t dphi = hadr_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(     mom >= 0.5  &&  mom < 1.0){   hadr_dphi_deta_p2_0510->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 1.0  &&  mom < 2.0){   hadr_dphi_deta_p2_1020->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 2.0  &&  mom < 3.0){   hadr_dphi_deta_p2_2030->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 3.0  &&  mom < 4.0){   hadr_dphi_deta_p2_3040->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 4.0  &&  mom < 5.0){   hadr_dphi_deta_p2_4050->Fill(dphi, hadr_eta - eta);    }
			    }
			    else if(hadr_pt >= 3.4  &&  hadr_pt < 4.4)
			    {
				Double_t eta  = track->Eta();	 Double_t dphi = hadr_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(     mom >= 0.5  &&  mom < 1.0){   hadr_dphi_deta_p3_0510->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 1.0  &&  mom < 2.0){   hadr_dphi_deta_p3_1020->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 2.0  &&  mom < 3.0){   hadr_dphi_deta_p3_2030->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 3.0  &&  mom < 4.0){   hadr_dphi_deta_p3_3040->Fill(dphi, hadr_eta - eta);    }
				else if(mom >= 4.0  &&  mom < 5.0){   hadr_dphi_deta_p3_4050->Fill(dphi, hadr_eta - eta);    }
			    }
			}
		    }
		}
	    
		int k = conductor->track_num;
		conductor = conductor->next;
		if(k != j)  // if the k-th track is not the j-th track which is the hadreron of interest
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
			    Double_t dphi = hadr_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();

			    if(centrality <= 20)
			    {
				dphi_pt_hadr_C->Fill(hadr_pt, dphi);		dphi_kt_hadr_C->Fill(hadr_kt, dphi);
				dphi_et_hadr_C->Fill(hadr_et, dphi);		dphi_en_hadr_C->Fill(hadr_en, dphi);
			    }
			    else if(centrality > 60)
			    {
				dphi_pt_hadr_P->Fill(hadr_pt, dphi);		dphi_kt_hadr_P->Fill(hadr_kt, dphi);
				dphi_et_hadr_P->Fill(hadr_et, dphi);		dphi_en_hadr_P->Fill(hadr_en, dphi);
			    }
			    
			    dphi = hadr_phi - phi;
			    if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
			    if(dphi < -TMath::Pi())    dphi = dphi + TMath::TwoPi();
			    if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
			    if(dphi >  TMath::Pi())    dphi = dphi - TMath::TwoPi();
			    hadr_cor_radius->Fill(hadr_pt, TMath::Sign(1.0,dphi)*sqrt(pow(dphi,2) + pow(hadr_eta - eta,2)));
			    asso++;			    
			}

			if(hadr_pt >= 1.4  &&  hadr_pt < 2.4)
			{
			    Double_t eta  = track->Eta();    Double_t dphi = hadr_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(     mom >= 0.5  &&  mom < 1.0){   hadr_dphi_deta_p1_0510->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 1.0  &&  mom < 2.0){   hadr_dphi_deta_p1_1020->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 2.0  &&  mom < 3.0){   hadr_dphi_deta_p1_2030->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 3.0  &&  mom < 4.0){   hadr_dphi_deta_p1_3040->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 4.0  &&  mom < 5.0){   hadr_dphi_deta_p1_4050->Fill(dphi, hadr_eta - eta);    }
			}
			else if(hadr_pt >= 2.4  &&  hadr_pt < 3.4)
			{
			    Double_t eta  = track->Eta();    Double_t dphi = hadr_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(     mom >= 0.5  &&  mom < 1.0){   hadr_dphi_deta_p2_0510->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 1.0  &&  mom < 2.0){   hadr_dphi_deta_p2_1020->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 2.0  &&  mom < 3.0){   hadr_dphi_deta_p2_2030->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 3.0  &&  mom < 4.0){   hadr_dphi_deta_p2_3040->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 4.0  &&  mom < 5.0){   hadr_dphi_deta_p2_4050->Fill(dphi, hadr_eta - eta);    }
			}
			else if(hadr_pt >= 3.4  &&  hadr_pt < 4.4)
			{
			    Double_t eta  = track->Eta();    Double_t dphi = hadr_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();			    
			    if(     mom >= 0.5  &&  mom < 1.0){   hadr_dphi_deta_p3_0510->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 1.0  &&  mom < 2.0){   hadr_dphi_deta_p3_1020->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 2.0  &&  mom < 3.0){   hadr_dphi_deta_p3_2030->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 3.0  &&  mom < 4.0){   hadr_dphi_deta_p3_3040->Fill(dphi, hadr_eta - eta);    }
			    else if(mom >= 4.0  &&  mom < 5.0){   hadr_dphi_deta_p3_4050->Fill(dphi, hadr_eta - eta);    }
			}			
		    }
		}
	    }

	    if(asso > 0)
	    {
		if(centrality <= 20.0)
		{
		    hadr_pt_count_C->Fill(hadr_pt);	    hadr_kt_count_C->Fill(hadr_kt);
		    hadr_et_count_C->Fill(hadr_et);	    hadr_en_count_C->Fill(hadr_en);
		}
		else if(centrality > 60)
		{
		    hadr_pt_count_P->Fill(hadr_pt);	    hadr_kt_count_P->Fill(hadr_kt);
		    hadr_et_count_P->Fill(hadr_et);	    hadr_en_count_P->Fill(hadr_en);
		}
	    }
	}
    }  // end of hadreron correlation PRIMARY LOOP

                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file

}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFhadr::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFhadr::Beta(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFhadr::tof_minus_tpion(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFhadr::get_mass_squared(AliAODTrack *track)
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
