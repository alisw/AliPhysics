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
#include "AliAnalysisTaskCorPIDTOFQA.h"
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




//class AliAnalysisTaskCorPIDTOFQA;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

//using namespace BSchaefer_devel;

//ofstream fout;
Int_t NEvents            = 0;
//Int_t NTracks            = 0;
//Int_t NGoodTracks        = 0;
//Int_t NAssociatedTracks  = 0;
//Int_t DeuteronCandidates = 0;


// SetCutGeoNcrNcl  cut out 3 cm on the sides of the TPC



//Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]
//TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);

ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
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




    deut_dphi_T(0),             // 27
    deut_dphi_pos_T(0),         // 28
    deut_dphi_neg_T(0),         // 29
    deut_dphi_A(0),             // 30
    deut_dphi_pos_A(0),         // 31
    deut_dphi_neg_A(0),         // 32
    deut_dphi_B(0),             // 33
    deut_dphi_pos_B(0),         // 34
    deut_dphi_neg_B(0),         // 35

    deut_per_event(0),          // 36
    deut_per_event_pos(0),      // 37
    deut_per_event_neg(0),      // 38
    
    m2_pos_cut_T(0),            // 39
    m2_pos_cut_A(0),            // 40
    m2_pos_cut_B(0),            // 41
    m2_neg_cut_T(0),            // 42
    m2_neg_cut_A(0),            // 43
    m2_neg_cut_B(0),            // 44

    dphi_ket_deut_T(0),         // 45
    dphi_ket_deut_A(0),         // 46
    dphi_ket_deut_B(0),         // 47

    deltat_channel(0)           // 48
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
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


    deut_dphi_T(0),             // 27
    deut_dphi_pos_T(0),         // 28
    deut_dphi_neg_T(0),         // 29
    deut_dphi_A(0),             // 30
    deut_dphi_pos_A(0),         // 31
    deut_dphi_neg_A(0),         // 32
    deut_dphi_B(0),             // 33
    deut_dphi_pos_B(0),         // 34
    deut_dphi_neg_B(0),         // 35

    deut_per_event(0),          // 36
    deut_per_event_pos(0),      // 37
    deut_per_event_neg(0),      // 38
    
    m2_pos_cut_T(0),            // 39
    m2_pos_cut_A(0),            // 40
    m2_pos_cut_B(0),            // 41
    m2_neg_cut_T(0),            // 42
    m2_neg_cut_A(0),            // 43
    m2_neg_cut_B(0),            // 44

    dphi_ket_deut_T(0),         // 45
    dphi_ket_deut_A(0),         // 46
    dphi_ket_deut_B(0),         // 47

    deltat_channel(0)           // 48
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::~AliAnalysisTaskCorPIDTOFQA()
{
    // destructor
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserCreateOutputObjects()
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

    deut_dphi_T           = new TH2F("deut_dphi_T",           "deut_dphi_T",            14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 27
    deut_dphi_pos_T       = new TH2F("deut_dphi_pos_T",       "deut_dphi_pos_T",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 28
    deut_dphi_neg_T       = new TH2F("deut_dphi_neg_T",       "deut_dphi_neg_T",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 29    
    deut_dphi_A           = new TH2F("deut_dphi_A",           "deut_dphi_A",            14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 30
    deut_dphi_pos_A       = new TH2F("deut_dphi_pos_A",       "deut_dphi_pos_A",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 31
    deut_dphi_neg_A       = new TH2F("deut_dphi_neg_A",       "deut_dphi_neg_A",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 32    
    deut_dphi_B           = new TH2F("deut_dphi_B",           "deut_dphi_B",            14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 33
    deut_dphi_pos_B       = new TH2F("deut_dphi_pos_B",       "deut_dphi_pos_B",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 34
    deut_dphi_neg_B       = new TH2F("deut_dphi_neg_B",       "deut_dphi_neg_B",        14,   1.0,    4.5,     300, -1.6708,  4.8124);   // 35
    deut_per_event        = new TH1I("deut_per_event",        "deut_per_event",          5,     0,     5);                               // 36
    deut_per_event_pos    = new TH1I("deut_per_event_pos",    "deut_per_event_pos",      5,     0,     5);                               // 37
    deut_per_event_neg    = new TH1I("deut_per_event_neg",    "deut_per_event_neg",      5,     0,     5);                               // 38
    m2_pos_cut_T          = new TH2F("m2_pos_cut_T",          "m2_pos_cut_T",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 39
    m2_pos_cut_A          = new TH2F("m2_pos_cut_A",          "m2_pos_cut_A",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 40
    m2_pos_cut_B          = new TH2F("m2_pos_cut_B",          "m2_pos_cut_B",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 41
    m2_neg_cut_T          = new TH2F("m2_neg_cut_T",          "m2_neg_cut_T",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 42
    m2_neg_cut_A          = new TH2F("m2_neg_cut_A",          "m2_neg_cut_A",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 43
    m2_neg_cut_B          = new TH2F("m2_neg_cut_B",          "m2_neg_cut_B",          320,   0.0,    8.0,    2400,    -1.0,     7.0);   // 44
    dphi_ket_deut_T       = new TH2F("dphi_ket_deut_T",       "dphi_ket_deut_T",        46, 0.445,  2.745,      60, -1.6708,  4.8124);   // 45 // dphi, ket
    dphi_ket_deut_A       = new TH2F("dphi_ket_deut_A",       "dphi_ket_deut_A",        46, 0.445,  2.745,      60, -1.6708,  4.8124);   // 46 // dphi, ket
    dphi_ket_deut_B       = new TH2F("dphi_ket_deut_B",       "dphi_ket_deut_B",        46, 0.445,  2.745,      60, -1.6708,  4.8124);   // 47 // dphi, ket

    deltat_channel        = new TProfile("deltat_channel",    "deltat_channel",       4800,     0,   4800,    -0.6,     0.6);            // 48 // tof channel and deltat
    
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
    fOutputList->Add(deut_dphi_T);            // 27
    fOutputList->Add(deut_dphi_pos_T);        // 28
    fOutputList->Add(deut_dphi_neg_T);        // 29
    fOutputList->Add(deut_dphi_A);            // 30
    fOutputList->Add(deut_dphi_pos_A);        // 31
    fOutputList->Add(deut_dphi_neg_A);        // 32
    fOutputList->Add(deut_dphi_B);            // 33
    fOutputList->Add(deut_dphi_pos_B);        // 34
    fOutputList->Add(deut_dphi_neg_B);        // 35
    fOutputList->Add(deut_per_event);         // 36
    fOutputList->Add(deut_per_event_pos);     // 37
    fOutputList->Add(deut_per_event_neg);     // 38
    fOutputList->Add(m2_pos_cut_T);           // 39
    fOutputList->Add(m2_pos_cut_A);           // 40
    fOutputList->Add(m2_pos_cut_B);	      // 41
    fOutputList->Add(m2_neg_cut_T);           // 42
    fOutputList->Add(m2_neg_cut_A);           // 43
    fOutputList->Add(m2_neg_cut_B);           // 44
    fOutputList->Add(dphi_ket_deut_T);        // 45
    fOutputList->Add(dphi_ket_deut_A);        // 46
    fOutputList->Add(dphi_ket_deut_B);        // 47
    
    fOutputList->Add(deltat_channel);         // 48
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserExec(Option_t *)
{

//    	Double_t pi = 3.1415926535897932384626434;
    NEvents++;

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    Int_t iTracks(fAOD->GetNumberOfTracks());


    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////



    
    int NDeut           = 0;
    int NDeut_pos       = 0;
    int NDeut_neg       = 0;
    int deut_flag_T          = 0;
    int deut_flag_A          = 0;
    int deut_flag_B          = 0;
    int deut_track_num_T     = -9999;;
    int deut_track_num_A     = -9999;;
    int deut_track_num_B     = -9999;;
    int associated_tracks   = 0;



	
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

	if(!track->IsPrimaryCandidate())   continue;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	AliExternalTrackParam trkparam;
	trkparam.CopyFromVTrack(track);
	
	if(trkparam.GetX() > 3.) continue; // only valid for propagation inside the beam pipe
		
	// impact parameters
	Double_t b[2];
	Double_t bCov[3];
	if(!trkparam.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 10000., b, bCov))   continue;
	
	Double_t dcaxy = b[0];
	Double_t dcaz  = b[1];

	if(fabs(dcaxy) > 3.0)   continue;
	if(fabs(dcaz)  > 2.0)   continue;
//	dca->Fill(dcaxy, dcaz);


	UChar_t map = track->GetITSClusterMap();
	Int_t nITSpoints = 0;
	for(Int_t j=2; j<6; j++)
	{
		if(map&(1<<j)) ++nITSpoints;
	}
	

	if(nITSpoints < 2)      continue;
	
	Short_t nTPCclusters     = track->GetTPCsignalN();
	if(nTPCclusters<70)     continue;                          // No of TPC Clusters



	
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
	Double_t deut_mean      = 0.0;
	Double_t deut_sigma     = 0.0;
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

		    
		    if(mom >= 1.5  &&  mom < 4.4)
		    {
			for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
			deut_mean = fit_deut_curve->Eval(mom);
			for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
			deut_sigma = fit_deut_curve->Eval(mom);

			if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
			{
			    m2_pos_cut_T->Fill(mom,m2tof);
			    deut_track_num_T = i;
			    // deut_track_num_pos_T = i;
			    NDeut++;
			    NDeut_pos++;
			    deut_flag_T = 1;
			}
			else if(m2tof < deut_mean -1 + cut_width * deut_sigma  &&   m2tof > deut_mean -1 - cut_width * deut_sigma)
			{
			    m2_pos_cut_A->Fill(mom,m2tof);
			    deut_flag_A = 1;
			    deut_track_num_A = i;
			}
			else if(m2tof < deut_mean +1 + cut_width * deut_sigma  &&   m2tof > deut_mean +1 - cut_width * deut_sigma)
			{
			    m2_pos_cut_B->Fill(mom,m2tof);
			    deut_flag_B = 1;
			    deut_track_num_B = i;
			}		    
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
//		dedx_p_m2_neg->Fill(dedx, mom, m2tof);
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

		    if(mom >= 1.5  &&  mom < 4.4)
		    {
			for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
			deut_mean = fit_deut_curve->Eval(mom);
			for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
			deut_sigma = fit_deut_curve->Eval(mom);

			if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
			{
			    m2_neg_cut_T->Fill(mom,m2tof);
			    deut_track_num_T = i;
			    // deut_track_num_neg_T = i;
			    NDeut++;
			    NDeut_neg++;
			    deut_flag_T = 1;
			}
			else if(m2tof < deut_mean -1 + cut_width * deut_sigma  &&   m2tof > deut_mean -1 - cut_width * deut_sigma)
			{
			    m2_neg_cut_A->Fill(mom,m2tof);
			    deut_flag_A = 1;
			    deut_track_num_A = i;
			}
			else if(m2tof < deut_mean +1 + cut_width * deut_sigma  &&   m2tof > deut_mean +1 - cut_width * deut_sigma)
			{
			    m2_neg_cut_B->Fill(mom,m2tof);
			    deut_flag_B = 1;
			    deut_track_num_B = i;
			}		    
		    }
		    
		    
		}
	    }
	}    //   track has both TPC and TOF PID, tracks are o.k.
    }        //   end of track loop


    deut_per_event->Fill(NDeut);
    deut_per_event_pos->Fill(NDeut_pos);
    deut_per_event_neg->Fill(NDeut_neg);

    if(deut_flag_T > 0)
    {
//	DeuteronCandidates++;
	int j = deut_track_num_T;

	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t deut_phi = track->Phi();
	    Double_t deut_mom = track->P();
	    Short_t  charge   = track->Charge();
	    Double_t deut_eta = track->Eta();

	    Double_t deut_ket = 0.0;
	    
	    deut_ket = sqrt(pow(deut_mom,2) + pow(2.224,2)) - 2.224;
	    
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;
		    
		    if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
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
				Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
				deut_dphi_T->Fill(deut_mom, dphi);
				
				if(charge>0) {   deut_dphi_pos_T->Fill(deut_mom, dphi);   }
				if(charge<0) {   deut_dphi_neg_T->Fill(deut_mom, dphi);   }
				dphi_ket_deut_T->Fill(deut_ket, dphi);
				
			    }
			}
		    }
		}
	    
		int k = conductor->track_num;
		conductor = conductor->next;
//	    NAssociatedTracks++;
		if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t mom            = track->P();
			Double_t phi            = track->Phi();
			
			if(mom > 2.0  &&  mom < 5.0)
			{
			    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
			    deut_dphi_T->Fill(deut_mom, dphi);
			    if(charge>0) {   deut_dphi_pos_T->Fill(deut_mom, dphi);   }
			    if(charge<0) {   deut_dphi_neg_T->Fill(deut_mom, dphi);   }
			    dphi_ket_deut_T->Fill(deut_ket, dphi);
			}
//			if(deut_mom >= 3.25  &&  deut_mom < 4.25)
//			{
//			    Double_t eta            = track->Eta();			    
//			    Double_t dphi = deut_phi - phi;
//			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
//			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
//			    if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p0510->Fill(dphi, deut_eta - eta);    }
//			    else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1020->Fill(dphi, deut_eta - eta);    }
//			    else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2030->Fill(dphi, deut_eta - eta);    }
//			    if(     mom >= 1.0  &&  mom < 5.0){   deut_dphi_deta_p1050->Fill(dphi, deut_eta - eta);    }
//			}
		    }
		}
	    }
	}
    }  // end of deuteron correlation PRIMARY LOOP


    if(deut_flag_A > 0)
    {
//	DeuteronCandidates++;
	int j = deut_track_num_A;

	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t deut_phi = track->Phi();
	    Double_t deut_mom = track->P();
	    Short_t  charge   = track->Charge();
	    Double_t deut_eta = track->Eta();

	    Double_t deut_ket = 0.0;
	    
	    deut_ket = sqrt(pow(deut_mom,2) + pow(2.224,2)) - 2.224;
	    
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;
		    
		    if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
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
				Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
				deut_dphi_A->Fill(deut_mom, dphi);
				
				if(charge>0) {   deut_dphi_pos_A->Fill(deut_mom, dphi);   }
				if(charge<0) {   deut_dphi_neg_A->Fill(deut_mom, dphi);   }
				dphi_ket_deut_A->Fill(deut_ket, dphi);
				
			    }
//			    if(deut_mom >= 3.25  &&  deut_mom < 4.25)
//			    {
//				Double_t eta  = track->Eta();
//				Double_t dphi = deut_phi - phi;
//				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
//				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
//				if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p0510->Fill(dphi, deut_eta - eta);    }
//				else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1020->Fill(dphi, deut_eta - eta);    }
//				else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2030->Fill(dphi, deut_eta - eta);    }
//				if(     mom >= 1.0  &&  mom < 5.0){   deut_dphi_deta_p1050->Fill(dphi, deut_eta - eta);    }
//			    }
			}
		    }
		}
	    
		int k = conductor->track_num;
		conductor = conductor->next;
//	    NAssociatedTracks++;
		if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t mom            = track->P();
			Double_t phi            = track->Phi();
			
			if(mom > 2.0  &&  mom < 5.0)
			{
			    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
			    deut_dphi_A->Fill(deut_mom, dphi);
			    if(charge>0) {   deut_dphi_pos_A->Fill(deut_mom, dphi);   }
			    if(charge<0) {   deut_dphi_neg_A->Fill(deut_mom, dphi);   }
			    dphi_ket_deut_A->Fill(deut_ket, dphi);
			}
//			if(deut_mom >= 3.25  &&  deut_mom < 4.25)
//			{
//			    Double_t eta            = track->Eta();			    
//			    Double_t dphi = deut_phi - phi;
//			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
//			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
//			    if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p0510->Fill(dphi, deut_eta - eta);    }
//			    else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1020->Fill(dphi, deut_eta - eta);    }
//			    else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2030->Fill(dphi, deut_eta - eta);    }
//			    if(     mom >= 1.0  &&  mom < 5.0){   deut_dphi_deta_p1050->Fill(dphi, deut_eta - eta);    }
//			}
		    }
		}
	    }
	}
    }  // end of deuteron correlation A LOOP


    if(deut_flag_B > 0)
    {
//	DeuteronCandidates++;
	int j = deut_track_num_B;

	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t deut_phi = track->Phi();
	    Double_t deut_mom = track->P();
	    Short_t  charge   = track->Charge();
	    Double_t deut_eta = track->Eta();

	    Double_t deut_ket = 0.0;
	    
	    deut_ket = sqrt(pow(deut_mom,2) + pow(2.224,2)) - 2.224;
	    
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;
		    
		    if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
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
				Double_t dphi = deut_phi - phi;
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
				deut_dphi_B->Fill(deut_mom, dphi);
				
				if(charge>0) {   deut_dphi_pos_B->Fill(deut_mom, dphi);   }
				if(charge<0) {   deut_dphi_neg_B->Fill(deut_mom, dphi);   }
				dphi_ket_deut_B->Fill(deut_ket, dphi);
				
			    }
//			    if(deut_mom >= 3.25  &&  deut_mom < 4.25)
//			    {
//				Double_t eta  = track->Eta();
//				Double_t dphi = deut_phi - phi;
//				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
//				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
//				if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p0510->Fill(dphi, deut_eta - eta);    }
//				else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1020->Fill(dphi, deut_eta - eta);    }
//				else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2030->Fill(dphi, deut_eta - eta);    }
//				if(     mom >= 1.0  &&  mom < 5.0){   deut_dphi_deta_p1050->Fill(dphi, deut_eta - eta);    }
//			    }
			}
		    }
		}
	    
		int k = conductor->track_num;
		conductor = conductor->next;
//	    NAssociatedTracks++;
		if(k != j)  // if the k-th track is not the j-th track which is the deuteron of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t mom            = track->P();
			Double_t phi            = track->Phi();
			
			if(mom > 2.0  &&  mom < 5.0)
			{
			    Double_t dphi = deut_phi - phi;
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
			    deut_dphi_B->Fill(deut_mom, dphi);
			    if(charge>0) {   deut_dphi_pos_B->Fill(deut_mom, dphi);   }
			    if(charge<0) {   deut_dphi_neg_B->Fill(deut_mom, dphi);   }
			    dphi_ket_deut_B->Fill(deut_ket, dphi);
			}
//			if(deut_mom >= 3.25  &&  deut_mom < 4.25)
//			{
//			    Double_t eta            = track->Eta();			    
//			    Double_t dphi = deut_phi - phi;
//			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//			    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
//			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
//			    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
			    
//			    if(     mom >= 0.5  &&  mom < 1.0){   deut_dphi_deta_p0510->Fill(dphi, deut_eta - eta);    }
//			    else if(mom >= 1.0  &&  mom < 2.0){   deut_dphi_deta_p1020->Fill(dphi, deut_eta - eta);    }
//			    else if(mom >= 2.0  &&  mom < 3.0){   deut_dphi_deta_p2030->Fill(dphi, deut_eta - eta);    }
//			    if(     mom >= 1.0  &&  mom < 5.0){   deut_dphi_deta_p1050->Fill(dphi, deut_eta - eta);    }
//			}
		    }
		}
	    }
	}
    }  // end of deuteron correlation PRIMARY LOOP

     
                                                     // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file

//    NEvents++;

}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
    cout<<endl<<endl;
//    cout<<"total number deuteron candidates:  "<<DeuteronCandidates<<endl;
    cout<<"total number of events analyzed:   "<<NEvents<<endl;
//    cout<<"total number of tracks analyzed:   "<<NTracks<<endl;
//    cout<<"total number of good tracks:       "<<NGoodTracks<<endl;
//    cout<<"tracks from events with deut: "<<NAssociatedTracks<<endl;
    cout<<endl<<endl;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFQA::Beta(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFQA::tof_minus_tpion(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFQA::get_mass_squared(AliAODTrack *track)
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
