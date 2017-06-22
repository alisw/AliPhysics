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

/* AliAnaysisTaskCorPIDTOFQA
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




//class AliAnalysisTaskCorPIDTOFhadr;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

//using namespace BSchaefer_devel;

//ofstream fout;
//Int_t NTracks            = 0;
//Int_t NGoodTracks        = 0;
//Int_t NAssociatedTracks  = 0;
//Int_t DeuteronCandidates = 0;


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
    dedx_p_pos(0),              //  9
    dedx_p_neg(0),              // 10
    dedx_deltat_pos(0),         // 11
    dedx_deltat_neg(0),         // 12
    dedx_p_deltat_pos(0),       // 13
    dedx_p_deltat_neg(0),       // 14

    m2_pos_cut(0),              // 15
    m2_neg_cut(0),              // 16
    beta_p_pos_cut(0),          // 17
    beta_p_neg_cut(0),          // 18
    deltat_p_pos_cut(0),        // 19
    deltat_p_neg_cut(0),        // 20

    dedx_p_pos_cut(0),          // 21
    dedx_p_neg_cut(0),          // 22

    dedx_deltat_pos_cut(0),     // 23
    dedx_deltat_neg_cut(0),     // 24
    dedx_p_deltat_pos_cut(0),   // 25
    dedx_p_deltat_neg_cut(0),   // 26




    hadr_dphi_T(0),             // 27
    hadr_dphi_pos_T(0),         // 28
    hadr_dphi_neg_T(0),         // 29
//    hadr_dphi_A(0),             // 30
//    hadr_dphi_pos_A(0),         // 31
//    hadr_dphi_neg_A(0),         // 32
//    hadr_dphi_B(0),             // 33
//    hadr_dphi_pos_B(0),         // 34
//    hadr_dphi_neg_B(0),         // 35

    hadr_per_event(0),          // 36
    hadr_per_event_pos(0),      // 37
    hadr_per_event_neg(0),      // 38
    
    m2_pos_cut_T(0),            // 39
//    m2_pos_cut_A(0),            // 40
//    m2_pos_cut_B(0),            // 41
    m2_neg_cut_T(0),            // 42
//    m2_neg_cut_A(0),            // 43
//    m2_neg_cut_B(0),            // 44

    dphi_et_hadr_T(0),         // 45
//    dphi_et_hadr_A(0),         // 46
//    dphi_et_hadr_B(0),         // 47

    deltat_channel(0),           // 48

    hadr_et_count(0),          // 49
    hadr_e_count(0),           // 50
    dphi_e_hadr_T(0),          // 51

    hadr_pt_count(0),           // 52
    asso_phi_hist(0),           // 53
    hadr_phi_hist(0)            // 54  
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
    dedx_p_pos(0),              //  9
    dedx_p_neg(0),              // 10
    dedx_deltat_pos(0),         // 11
    dedx_deltat_neg(0),         // 12
    dedx_p_deltat_pos(0),       // 13
    dedx_p_deltat_neg(0),       // 14

    m2_pos_cut(0),              // 15
    m2_neg_cut(0),              // 16
    beta_p_pos_cut(0),          // 17
    beta_p_neg_cut(0),          // 18
    deltat_p_pos_cut(0),        // 19
    deltat_p_neg_cut(0),        // 20

    dedx_p_pos_cut(0),          // 21
    dedx_p_neg_cut(0),          // 22

    dedx_deltat_pos_cut(0),     // 23
    dedx_deltat_neg_cut(0),     // 24
    dedx_p_deltat_pos_cut(0),   // 25
    dedx_p_deltat_neg_cut(0),   // 26


    hadr_dphi_T(0),             // 27
    hadr_dphi_pos_T(0),         // 28
    hadr_dphi_neg_T(0),         // 29
//    hadr_dphi_A(0),             // 30
//    hadr_dphi_pos_A(0),         // 31
//    hadr_dphi_neg_A(0),         // 32
//    hadr_dphi_B(0),             // 33
//    hadr_dphi_pos_B(0),         // 34
//    hadr_dphi_neg_B(0),         // 35

    hadr_per_event(0),          // 36
    hadr_per_event_pos(0),      // 37
    hadr_per_event_neg(0),      // 38
    
    m2_pos_cut_T(0),            // 39
//    m2_pos_cut_A(0),            // 40
//    m2_pos_cut_B(0),            // 41
    m2_neg_cut_T(0),            // 42
//    m2_neg_cut_A(0),            // 43
//    m2_neg_cut_B(0),            // 44

    dphi_et_hadr_T(0),         // 45
//    dphi_et_hadr_A(0),         // 46
//    dphi_et_hadr_B(0),         // 47

    deltat_channel(0),           // 48

    hadr_et_count(0),          // 49
    hadr_e_count(0),           // 50
    dphi_e_hadr_T(0),          // 51

    hadr_pt_count(0),           // 52
    asso_phi_hist(0),           // 53
    hadr_phi_hist(0)            // 54  
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

//    pi = TMath::Pi();
//    fout.open("output.txt");


//  Double_t prot_curves[2][2][3];  // [charge][mean,sigma][par]

//  TF1 *fit_prot_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);
//    fit_prot_curve->SetParNames("a", "b x", "c/#sqrt(x)");
	
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    fHistPt               = new TH1F("fHistPt",               "Pt()",                 1000,     0,     10);                              //  1
    cent_ntracks          = new TH2F("cent_ntracks",          "cent_ntracks",          100,     0,    100,     100,       0,     800);   //  2
    
    m2_pos                = new TH2F("m2_pos",                "m2_pos",                320,   0.0,    8.0,    2400,    -1.0,     7.0);   //  3
    m2_neg                = new TH2F("m2_neg",                "m2_neg",                320,   0.0,    8.0,    2400,    -1.0,     7.0);   //  4
    beta_p_pos            = new TH2F("beta_p_pos",            "beta_p_pos",            780,   0.0,    8.0,    3000,     0.1,     1.1);   //  5
    beta_p_neg            = new TH2F("beta_p_neg",            "beta_p_neg",            780,   0.0,    8.0,    3000,     0.1,     1.1);   //  6
    deltat_p_pos          = new TH2F("deltat_p_pos",          "deltat_p_pos",        27000,   0.3,    3.0,    1000,    -1.0,     9.0);   //  7
    deltat_p_neg          = new TH2F("deltat_p_neg",          "deltat_p_neg",        27000,   0.3,    3.0,    1000,    -1.0,     9.0);   //  8
    dedx_p_pos            = new TH2F("dedx_p_pos",            "dedx_p_pos",            780,   0.2,    8.0,     250,     0.0,   500.0);   //  9
    dedx_p_neg            = new TH2F("dedx_p_neg",            "dedx_p_neg",            780,   0.2,    8.0,     250,     0.0,   500.0);   // 10
    dedx_deltat_pos       = new TH2F("dedx_deltat_pos",       "dedx_deltat_pos",      2100,  -1.0,   20.0,     250,     0.0,   500.0);   // 11
    dedx_deltat_neg       = new TH2F("dedx_deltat_neg",       "dedx_deltat_neg",      2100,  -1.0,   20.0,     250,     0.0,   500.0);   // 12
    dedx_p_deltat_pos     = new TH3F("dedx_p_deltat_pos",     "dedx_p_deltat_pos",     250, 0.0, 500.0, 156, 0.2, 8.0, 420, -1.0, 20.0); // 13 //// (dedx, mom, deltat)
    dedx_p_deltat_neg     = new TH3F("dedx_p_deltat_neg",     "dedx_p_deltat_neg",     250, 0.0, 500.0, 156, 0.2, 8.0, 420, -1.0, 20.0); // 14 //// (dedx, mom, deltat)
    m2_pos_cut            = new TH2F("m2_pos_cut",            "m2_pos_cut",            320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 15
    m2_neg_cut            = new TH2F("m2_neg_cut",            "m2_neg_cut",            320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 16
    beta_p_pos_cut        = new TH2F("beta_p_pos_cut",        "beta_p_pos_cut",        780,   0.0,    8.0,    3000,     0.1,     1.1);   // 17
    beta_p_neg_cut        = new TH2F("beta_p_neg_cut",        "beta_p_neg_cut",        780,   0.0,    8.0,    3000,     0.1,     1.1);   // 18
    deltat_p_pos_cut      = new TH2F("deltat_p_pos_cut",      "deltat_p_pos_cut",    27000,   0.3,    3.0,    1000,    -1.0,     9.0);   // 19
    deltat_p_neg_cut      = new TH2F("deltat_p_neg_cut",      "deltat_p_neg_cut",    27000,   0.3,    3.0,    1000,    -1.0,     9.0);   // 20
    dedx_p_pos_cut        = new TH2F("dedx_p_pos_cut",        "dedx_p_pos_cut",        780,   0.2,    8.0,     250,     0.0,   500.0);   // 21
    dedx_p_neg_cut        = new TH2F("dedx_p_neg_cut",        "dedx_p_neg_cut",        780,   0.2,    8.0,     250,     0.0,   500.0);   // 22
    dedx_deltat_pos_cut   = new TH2F("dedx_deltat_pos_cut",   "dedx_deltat_pos_cut",  2100,  -1.0,   20.0,     250,     0.0,   500.0);   // 23
    dedx_deltat_neg_cut   = new TH2F("dedx_deltat_neg_cut",   "dedx_deltat_neg_cut",  2100,  -1.0,   20.0,     250,     0.0,   500.0);   // 24
    dedx_p_deltat_pos_cut = new TH3F("dedx_p_deltat_pos_cut", "dedx_p_deltat_pos_cut", 250, 0.0, 500.0, 156, 0.2, 8.0, 420, -1.0, 20.0); // 25 //// (dedx, mom, deltat)
    dedx_p_deltat_neg_cut = new TH3F("dedx_p_deltat_neg_cut", "dedx_p_deltat_neg_cut", 250, 0.0, 500.0, 156, 0.2, 8.0, 420, -1.0, 20.0); // 26 //// (dedx, mom, deltat)

    hadr_dphi_T           = new TH2F("hadr_dphi_T",           "hadr_dphi_T",           160,   1.0,    5.0,     300, -1.6708,  4.8124);   // 27
    hadr_dphi_pos_T       = new TH2F("hadr_dphi_pos_T",       "hadr_dphi_pos_T",       160,   1.0,    5.0,     300, -1.6708,  4.8124);   // 28
    hadr_dphi_neg_T       = new TH2F("hadr_dphi_neg_T",       "hadr_dphi_neg_T",       160,   1.0,    5.0,     300, -1.6708,  4.8124);   // 29    
//    hadr_dphi_A           = new TH2F("hadr_dphi_A",           "hadr_dphi_A",            14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 30
//    hadr_dphi_pos_A       = new TH2F("hadr_dphi_pos_A",       "hadr_dphi_pos_A",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 31
//    hadr_dphi_neg_A       = new TH2F("hadr_dphi_neg_A",       "hadr_dphi_neg_A",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 32    
//    hadr_dphi_B           = new TH2F("hadr_dphi_B",           "hadr_dphi_B",            14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 33
//    hadr_dphi_pos_B       = new TH2F("hadr_dphi_pos_B",       "hadr_dphi_pos_B",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 34
//    hadr_dphi_neg_B       = new TH2F("hadr_dphi_neg_B",       "hadr_dphi_neg_B",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 35
    hadr_per_event        = new TH1I("hadr_per_event",        "hadr_per_event",         20,     0,    20);                               // 36
    hadr_per_event_pos    = new TH1I("hadr_per_event_pos",    "hadr_per_event_pos",     20,     0,    20);                               // 37
    hadr_per_event_neg    = new TH1I("hadr_per_event_neg",    "hadr_per_event_neg",     20,     0,    20);                               // 38
    m2_pos_cut_T          = new TH2F("m2_pos_cut_T",          "m2_pos_cut_T",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 39
//    m2_pos_cut_A          = new TH2F("m2_pos_cut_A",          "m2_pos_cut_A",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 40
//    m2_pos_cut_B          = new TH2F("m2_pos_cut_B",          "m2_pos_cut_B",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 41
    m2_neg_cut_T          = new TH2F("m2_neg_cut_T",          "m2_neg_cut_T",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 42
//    m2_neg_cut_A          = new TH2F("m2_neg_cut_A",          "m2_neg_cut_A",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 43
//    m2_neg_cut_B          = new TH2F("m2_neg_cut_B",          "m2_neg_cut_B",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 44
//    dphi_et_hadr_T       = new TH2F("dphi_et_hadr_T",       "dphi_et_hadr_T",        46, 0.445,  2.745,      60, -1.6708,  4.8124);   // 45 // dphi, ket
    dphi_et_hadr_T        = new TH2F("dphi_et_hadr_T",        "dphi_et_hadr_T",        189,  1.00,   4.78,      60, -1.6708,  4.8124);   // 45 // dphi, et
//    dphi_et_hadr_A       = new TH2F("dphi_et_hadr_A",       "dphi_et_hadr_A",        46, 0.445,  2.745,      60, -1.6708,  4.8124);   // 46 // dphi, ket
//    dphi_et_hadr_B       = new TH2F("dphi_et_hadr_B",       "dphi_et_hadr_B",        46, 0.445,  2.745,      60, -1.6708,  4.8124);   // 47 // dphi, ket

    deltat_channel        = new TProfile("deltat_channel",    "deltat_channel",       4800,     0,   4800,    -0.6,     0.6);            // 48 // tof channel and deltat




    hadr_et_count         = new TH1F("hadr_et_count",         "hadr_et_count",         189,  1.00,   4.78);                              // 49
    hadr_e_count          = new TH1F("hadr_e_count",          "hadr_e_count",          189,  1.00,   4.78);                              // 50

    dphi_e_hadr_T         = new TH2F("dphi_e_hadr_T",         "dphi_e_hadr_T",         189,  1.00,   4.78,      60, -1.6708,  4.8124);   // 51 // dphi, e
    hadr_pt_count         = new TH1F("hadr_pt_count",         "hadr_pt_count",         160,   1.0,    5.0);                              // 52
  
    asso_phi_hist         = new TH1F("asso_phi_hist",         "asso_phi_hist",         480, -0.1, 6.383185);                             // 53
    hadr_phi_hist         = new TH1F("hadr_phi_hist",         "hadr_phi_hist",         480, -0.1, 6.383185);                             // 54
    
    
//    path_length             = new TH1F("path_length",               "path_length",             200, 382.0,  398.0);
//    ttof                    = new TH1F("ttof",                      "time of flight",          500,  12.0,   20.0);
//    alpha                   = new TH1F("alpha",                     "alpha",                  1000,     0,    4.0);
//    theta                   = new TH1F("theta",                     "theta",                   100,  -7.0,    7.0);
//    dca                     = new TH2F("dca",                       "dca",                     100,  -5.0,    5.0,     100,   -5.0,    5.0);
//    generic                 = new TH1F("generic",                   "generic",                 4, -1,3);
    //markc


    
    
    fOutputList->Add(fHistPt);                //  1                                                             // objects added to output file
    fOutputList->Add(cent_ntracks);           //  2
    
    fOutputList->Add(m2_pos);                 //  3
    fOutputList->Add(m2_neg);                 //  4 
    fOutputList->Add(beta_p_pos);             //  5
    fOutputList->Add(beta_p_neg);             //  6
    fOutputList->Add(deltat_p_pos);           //  7
    fOutputList->Add(deltat_p_neg);           //  8
    fOutputList->Add(dedx_p_pos);             //  9
    fOutputList->Add(dedx_p_neg);             // 10
    fOutputList->Add(dedx_deltat_pos);        // 11
    fOutputList->Add(dedx_deltat_neg);        // 12
    fOutputList->Add(dedx_p_deltat_pos);      // 13
    fOutputList->Add(dedx_p_deltat_neg);      // 14
    
    fOutputList->Add(m2_pos_cut);             // 15
    fOutputList->Add(m2_neg_cut);             // 16
    fOutputList->Add(beta_p_pos_cut);         // 17
    fOutputList->Add(beta_p_neg_cut);         // 18
    fOutputList->Add(deltat_p_pos_cut);       // 19
    fOutputList->Add(deltat_p_neg_cut);       // 20
    fOutputList->Add(dedx_p_pos_cut);         // 21
    fOutputList->Add(dedx_p_neg_cut);         // 22
    fOutputList->Add(dedx_deltat_pos_cut);    // 23   
    fOutputList->Add(dedx_deltat_neg_cut);    // 24
    fOutputList->Add(dedx_p_deltat_pos_cut);  // 25  
    fOutputList->Add(dedx_p_deltat_neg_cut);  // 26
    fOutputList->Add(hadr_dphi_T);            // 27
    fOutputList->Add(hadr_dphi_pos_T);        // 28
    fOutputList->Add(hadr_dphi_neg_T);        // 29
//    fOutputList->Add(hadr_dphi_A);            // 30
//    fOutputList->Add(hadr_dphi_pos_A);        // 31
//    fOutputList->Add(hadr_dphi_neg_A);        // 32
//    fOutputList->Add(hadr_dphi_B);            // 33
//    fOutputList->Add(hadr_dphi_pos_B);        // 34
//    fOutputList->Add(hadr_dphi_neg_B);        // 35
    fOutputList->Add(hadr_per_event);         // 36
    fOutputList->Add(hadr_per_event_pos);     // 37
    fOutputList->Add(hadr_per_event_neg);     // 38
    fOutputList->Add(m2_pos_cut_T);           // 39
//    fOutputList->Add(m2_pos_cut_A);           // 40
//    fOutputList->Add(m2_pos_cut_B);	      // 41
    fOutputList->Add(m2_neg_cut_T);           // 42
//    fOutputList->Add(m2_neg_cut_A);           // 43
//    fOutputList->Add(m2_neg_cut_B);           // 44
    fOutputList->Add(dphi_et_hadr_T);        // 45
//    fOutputList->Add(dphi_et_hadr_A);        // 46
//    fOutputList->Add(dphi_et_hadr_B);        // 47
    
    fOutputList->Add(deltat_channel);         // 48


    fOutputList->Add(hadr_et_count);          // 49
    fOutputList->Add(hadr_e_count);           // 50
    fOutputList->Add(dphi_e_hadr_T);          // 51
    fOutputList->Add(hadr_pt_count);          // 52
    fOutputList->Add(asso_phi_hist);          // 53
    fOutputList->Add(hadr_phi_hist);          // 54
    
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


    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////



    
    int NHadr                = 0;
    int NHadr_pos            = 0;
    int NHadr_neg            = 0;
    int hadr_flag_T          = 0;
    int hadr_track_num_T     = -9999;;
    int associated_tracks    = 0;



	
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
//	NTracks++;


        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                                                       // if we failed, skip this track
	if(!(track->IsHybridGlobalConstrainedGlobal())){   continue;   }

//	if(!track->IsPrimaryCandidate())   continue;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//	AliExternalTrackParam trkparam;
//	trkparam.CopyFromVTrack(track);
	
//	if(trkparam.GetX() > 3.) continue; // only valid for propagation inside the beam pipe
		
	// impact parameters
//	Double_t b[2];
//	Double_t bCov[3];
//	if(!trkparam.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 10000., b, bCov))   continue;
	
//	Double_t dcaxy = b[0];
//	Double_t dcaz  = b[1];

//	if(fabs(dcaxy) > 3.0)   continue;
//	if(fabs(dcaz)  > 2.0)   continue;
//	dca->Fill(dcaxy, dcaz);


//	UChar_t map = track->GetITSClusterMap();
//	Int_t nITSpoints = 0;
//	for(Int_t j=2; j<6; j++)
//	{
//		if(map&(1<<j)) ++nITSpoints;
//	}
	

//	if(nITSpoints < 2)      continue;
	
//	Short_t nTPCclusters     = track->GetTPCsignalN();
//	if(nTPCclusters<70)     continue;                          // No of TPC Clusters



	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

//	if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { continue; }

//	if(!fTrackSelection->IsTrackAccepted(track)){   continue;   }  // from markus' code
//	if(!(track->TestFilterBit(AliAODTrack::klsHybridTPCCG)) ) { continue; }
//	DCA_XY->Fill(track->XAtDCA(), track->YAtDCA());

//	NGoodTracks++;
//	associated_tracks++;
	track_node *new_track = new track_node;
	new_track->track_num = i;
	new_track->next = NULL;
	conductor->next = new_track;
	conductor = conductor->next;

	
	short isGlobal    = 1;
	short isModGlobal = 0;  // is it suitable for proteron candidate
	
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



//	if(tpcIsOk  &&  tofIsOk)
//	if(tofIsOk)
//	{

	    
	Double_t mom            = track->P();
//	Double_t energy         = track->E();
//	float    mass           = track->M();
	float    pt             = track->Pt();
//	Double_t tpc_mom        = track->GetTPCmomentum();
	Short_t  charge         = track->Charge();
//	Double_t eta            = track->Eta();

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
//	if(pt >= 0.3  &&  pt < 5)   alpha->Fill(1/pt);

//	if(pt >= 1.0  &&  pt < 1.05)   cout<<fPIDResponse->GetTOFResponse().GetTOFchannel(track)<<endl;

	// markc
	if(isModGlobal == 1)
	{
	    Double_t m2tof = get_mass_squared(track);
	    Float_t dedx   = track->GetTPCsignal();
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
		dedx_p_pos->Fill(   mom, dedx);
		if(mom >= 2.0  &&  mom < 3.0)    dedx_deltat_pos->Fill(deltat, dedx);
		dedx_p_deltat_pos->Fill(dedx, mom, deltat);
//		dedx_p_m2_pos->Fill(dedx, mom, m2tof);
		deltat_p_pos->Fill(  pt, deltat);

		
		if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340  &&  isGlobal == 1)
		{
		    if(    (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
			   )
		    {
			AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
			if(nsigmaTPC <= 2.0)
			{
			    m2_pos_cut->Fill(        pt, m2tof);
			    beta_p_pos_cut->Fill(   mom, Beta(track));
			    deltat_p_pos_cut->Fill(  pt, deltat);
			    dedx_p_pos_cut->Fill(   mom, dedx);
			    if(mom >= 2.0  &&  mom < 3.0)    dedx_deltat_pos_cut->Fill(deltat, dedx);
			    dedx_p_deltat_pos_cut->Fill(dedx, mom, deltat);
//			    dedx_p_m2_pos_cut->Fill(dedx, mom, m2tof);
			}
		    }

		    Int_t channel           = 0;
		    channel                 = fPIDResponse->GetTOFResponse().GetTOFchannel(track);
//		    Double_t xy_length      = 0.00;
//		    xy_length               = length * TMath::Sin(track_theta);

		    if(channel > -1  &&  channel < 4801)
		    {
			if(pt >= 0.995  &&  pt < 1.005)
			{
			    if(deltat >= -0.6  &&  deltat < 0.6)
			    {
				deltat_channel->Fill(channel, deltat);
			    }
			}
		    }

		    
		    if(mom >= 1.0  &&  mom < 4.0)
		    {
			m2_pos_cut_T->Fill(mom,m2tof);
			hadr_track_num_T = i;
			NHadr++;
			NHadr_pos++;
			hadr_flag_T = 1;
		    }
		}
	    }
	    
	    else if(charge < 0)
	    {
		m2_neg->Fill(      pt,m2tof);
		beta_p_neg->Fill(   mom, Beta(track));
		dedx_p_neg->Fill(   mom, dedx);
		if(mom >= 2.0  &&  mom < 3.0)    dedx_deltat_neg->Fill(deltat, dedx);
		dedx_p_deltat_neg->Fill(dedx, mom, deltat);
		deltat_p_neg->Fill(  pt, deltat);

		
		if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340  &&  isGlobal == 1)
		{
		    if(    (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
			   )
		    {
			AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
			if(nsigmaTPC <= 2.0)
			{	
			    m2_neg_cut->Fill(      pt, m2tof);
			    beta_p_neg_cut->Fill(   mom, Beta(track));
			    deltat_p_neg_cut->Fill(  pt, deltat);
			    dedx_p_neg_cut->Fill(   mom, dedx);
			    if(mom >= 2.0  &&  mom < 3.0)    dedx_deltat_neg_cut->Fill(deltat, dedx);
			    dedx_p_deltat_neg_cut->Fill(dedx, mom, deltat);
			}
		    }

		    if(mom >= 1.0  &&  mom < 4.0)
		    {
			m2_neg_cut_T->Fill(mom,m2tof);
			hadr_track_num_T = i;
		     // hadr_track_num_neg_T = i;
			NHadr++;
			NHadr_neg++;
			hadr_flag_T = 1;
		    }
		}
	    }
	}    //   track has both TPC and TOF PID, tracks are o.k.
    }        //   end of track loop


    hadr_per_event->Fill(NHadr);
    hadr_per_event_pos->Fill(NHadr_pos);
    hadr_per_event_neg->Fill(NHadr_neg);

    if(hadr_flag_T > 0)
    {
//	HadreronCandidates++;
	int j = hadr_track_num_T;

	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t hadr_phi = track->Phi();
	    Double_t hadr_pt  = track->Pt();
	    Double_t hadr_mom = track->P();
	    Short_t  charge   = track->Charge();
	    Double_t hadr_eta = track->Eta();

	    Double_t hadr_et = 0.0;
	    Double_t hadr_e   = 0.0;
	    
	    hadr_et = sqrt(pow(hadr_pt, 2) + pow(0.2902614,2));// - 0.2902614;
	    hadr_e  = sqrt(pow(hadr_mom,2) + pow(0.2902614,2));// - 0.2902614;
	    

	    hadr_phi_hist->Fill(hadr_phi);

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
//		    NAssociatedTracks++;
			AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
			if(!track){}                                                               // if we failed, skip this track
			else
			{
			    Double_t mom            = track->P();
		 	    Double_t phi            = track->Phi();
			    

			    if(mom > 2.0  &&  mom < 5.0)
			    {
				Double_t dphi = hadr_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
				hadr_dphi_T->Fill(hadr_mom, dphi);
				
				if(charge>0) {   hadr_dphi_pos_T->Fill(hadr_mom, dphi);   }
				if(charge<0) {   hadr_dphi_neg_T->Fill(hadr_mom, dphi);   }
				dphi_et_hadr_T->Fill(hadr_et, dphi);
				dphi_e_hadr_T->Fill(hadr_e, dphi);
				asso_phi_hist->Fill(phi);
				asso++;
			    }
			}
		    }
		}
	    
		int k = conductor->track_num;
		conductor = conductor->next;
//	    NAssociatedTracks++;
		if(k != j)  // if the k-th track is not the j-th track which is the proteron of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t mom            = track->P();
			Double_t phi            = track->Phi();
			
			if(mom > 2.0  &&  mom < 5.0)
			{
			    Double_t dphi = hadr_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
			    hadr_dphi_T->Fill(hadr_mom, dphi);
			    if(charge>0) {   hadr_dphi_pos_T->Fill(hadr_mom, dphi);   }
			    if(charge<0) {   hadr_dphi_neg_T->Fill(hadr_mom, dphi);   }
			    dphi_et_hadr_T->Fill(hadr_et, dphi);
			    dphi_e_hadr_T->Fill(hadr_e, dphi);
			    asso_phi_hist->Fill(phi);
			    asso++;
			}
		    }
		}
	    }

	    if(asso > 0)
	    {
		hadr_et_count->Fill(hadr_et);
		hadr_e_count->Fill(hadr_e);
		hadr_pt_count->Fill(hadr_pt);
	    }
	}
    }  // end of proteron correlation PRIMARY LOOP


     
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
//  Double_t start_time     = fPIDResponse->GetTOFResponse().GetTimeZero() * 1e-12;
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
//  Double_t start_time     = fPIDResponse->GetTOFResponse().GetTimeZero() * 1e-12;                         // lower resolution than previous line
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
