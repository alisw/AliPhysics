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

/* AliAnaysisTaskCorPIDTOFdiprot
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
//#include <iostream>
#include <fstream>
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
#include "AliAnalysisTaskCorPIDTOFdiprot.h"
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




using namespace std;            // std namespace: so you can do things like 'cout'
//using namespace BSchaefer_devel;

//ofstream file_output("output.txt");


ClassImp(AliAnalysisTaskCorPIDTOFdiprot) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFdiprot::AliAnalysisTaskCorPIDTOFdiprot() : AliAnalysisTaskSE(), 
fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),

    fHistPt(0),                    //  1
    cent_ntracks(0),               //  2
    
    m2_pt_pos(0),                  //  3
    m2_pt_neg(0),                  //  4
    beta_p_pos(0),                 //  5
    beta_p_neg(0),                 //  6

    deltat_pt_pos(0),              //  7
    deltat_pt_neg(0),              //  8

    m2_pt_pos_cut(0),              //  9
    m2_pt_neg_cut(0),              // 10
    beta_p_pos_cut(0),             // 11
    beta_p_neg_cut(0),             // 12
    deltat_pt_pos_cut(0),          // 13
    deltat_pt_neg_cut(0),          // 14

    prot_per_event(0),             // 15
//  prot_per_event_pos(0),         // 16
//  prot_per_event_neg(0),         // 17
    
    m2_pt_pos_cut_T(0),            // 18
    m2_pt_neg_cut_T(0),            // 19

    prot_phi_pt_pos(0),            // 20
    prot_phi_pt_neg(0),            // 21
    
    prot_p0_pt_pos_pos(0),         // 22
    prot_p0_pt_pos_neg(0),         // 23
    prot_p0_pt_neg_neg(0),         // 24
    
    di_prot_phi_pt(0),             // 25
    trig_05_phi_pt(0),             // 26
    trig_08_phi_pt(0),             // 27
    
    di_prot_dphi_p0_pos_pos_05(0), // 28
    di_prot_dphi_p0_pos_neg_05(0), // 29
    di_prot_dphi_p0_neg_neg_05(0), // 30

    di_prot_dphi_p0_pos_pos_08(0), // 31
    di_prot_dphi_p0_pos_neg_08(0), // 32
    di_prot_dphi_p0_neg_neg_08(0), // 33

    di_prot_dphi_r02_pos_pos_05(0),// 34
    di_prot_dphi_r02_pos_neg_05(0),// 35
    di_prot_dphi_r02_neg_neg_05(0),// 36

    di_prot_dphi_r02_pos_pos_08(0),// 37
    di_prot_dphi_r02_pos_neg_08(0),// 38
    di_prot_dphi_r02_neg_neg_08(0),// 39


    di_prot_dphi_r04_pos_pos_05(0),// 40
    di_prot_dphi_r04_pos_neg_05(0),// 41
    di_prot_dphi_r04_neg_neg_05(0),// 42

    di_prot_dphi_r04_pos_pos_08(0),// 43
    di_prot_dphi_r04_pos_neg_08(0),// 44
    di_prot_dphi_r04_neg_neg_08(0),// 45
    
    primary_vertex_z(0),           // 46
    primary_vertex_z_cut1(0),      // 47
    primary_vertex_z_cut2(0)       // 48

//    track_cor_radius_pt(0),        // 34
//    track_cor_radius_pt_cut(0)     // 35
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFdiprot::AliAnalysisTaskCorPIDTOFdiprot(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),
									   

    fHistPt(0),                    //  1
    cent_ntracks(0),               //  2
    
    m2_pt_pos(0),                  //  3
    m2_pt_neg(0),                  //  4
    beta_p_pos(0),                 //  5
    beta_p_neg(0),                 //  6

    deltat_pt_pos(0),              //  7
    deltat_pt_neg(0),              //  8

    m2_pt_pos_cut(0),              //  9
    m2_pt_neg_cut(0),              // 10
    beta_p_pos_cut(0),             // 11
    beta_p_neg_cut(0),             // 12
    deltat_pt_pos_cut(0),          // 13
    deltat_pt_neg_cut(0),          // 14

    prot_per_event(0),             // 15
//  prot_per_event_pos(0),         // 16
//  prot_per_event_neg(0),         // 17
    
    m2_pt_pos_cut_T(0),            // 18
    m2_pt_neg_cut_T(0),            // 19

    prot_phi_pt_pos(0),            // 20
    prot_phi_pt_neg(0),            // 21

    prot_p0_pt_pos_pos(0),         // 22
    prot_p0_pt_pos_neg(0),         // 23
    prot_p0_pt_neg_neg(0),         // 24
    
    di_prot_phi_pt(0),             // 25
    trig_05_phi_pt(0),             // 26
    trig_08_phi_pt(0),             // 27
    
    di_prot_dphi_p0_pos_pos_05(0), // 28
    di_prot_dphi_p0_pos_neg_05(0), // 29
    di_prot_dphi_p0_neg_neg_05(0), // 30

    di_prot_dphi_p0_pos_pos_08(0), // 31
    di_prot_dphi_p0_pos_neg_08(0), // 32
    di_prot_dphi_p0_neg_neg_08(0), // 33

    di_prot_dphi_r02_pos_pos_05(0),// 34
    di_prot_dphi_r02_pos_neg_05(0),// 35
    di_prot_dphi_r02_neg_neg_05(0),// 36

    di_prot_dphi_r02_pos_pos_08(0),// 37
    di_prot_dphi_r02_pos_neg_08(0),// 38
    di_prot_dphi_r02_neg_neg_08(0),// 39


    di_prot_dphi_r04_pos_pos_05(0),// 40
    di_prot_dphi_r04_pos_neg_05(0),// 41
    di_prot_dphi_r04_neg_neg_05(0),// 42

    di_prot_dphi_r04_pos_pos_08(0),// 43
    di_prot_dphi_r04_pos_neg_08(0),// 44
    di_prot_dphi_r04_neg_neg_08(0),// 45										   
//    track_cor_radius_pt(0),        // 34
//    track_cor_radius_pt_cut(0)     // 35

    primary_vertex_z(0),           // 46
    primary_vertex_z_cut1(0),      // 47
    primary_vertex_z_cut2(0)       // 48									       
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFdiprot::~AliAnalysisTaskCorPIDTOFdiprot()
{
    // destructor
    if(fAnalysisUtils) delete fAnalysisUtils;
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdiprot::UserCreateOutputObjects()
{
    prot_curves[0][0][0] = 0.763625;      // pos prot mean curve
    prot_curves[0][0][1] = 0.0194772;
    prot_curves[0][0][2] = 0.0903159;
    prot_curves[0][0][3] =-0.000109114;
    
    prot_curves[0][1][0] =-0.129833;      // pos prot sigma curve
    prot_curves[0][1][1] = 0.0625516;
    prot_curves[0][1][2] = 0.10124;
    prot_curves[0][1][3] = 0.000109926;

    prot_curves[1][0][0] = 0.778006;      // neg prot mean curve
    prot_curves[1][0][1] = 0.0183448;
    prot_curves[1][0][2] = 0.0754383;
    prot_curves[1][0][3] =-0.00010526;

    prot_curves[1][1][0] =-0.117123;      // neg prot sigma curve
    prot_curves[1][1][1] = 0.0597632;
    prot_curves[1][1][2] = 0.0921575;
    prot_curves[1][1][3] = 0.000118412;


    fAnalysisUtils = new AliAnalysisUtils;
    
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    Double_t pt_binning[2001];


    Float_t moving_marker = 0.10;
    for(int i=0; i<1202; i++)
    {
	pt_binning[i] = moving_marker;
	moving_marker = moving_marker + pt_binning[i] * 0.005;
    }
    
    
     
    fHistPt                    = new TH1F("fHistPt",                    "Pt()",                        800,       pt_binning);                              //  1
    cent_ntracks               = new TH2F("cent_ntracks",               "cent_ntracks",                100,        0,    100,     100,       0,     800);   //  2
    
    m2_pt_pos                  = new TH2F("m2_pt_pos",                  "m2_pt_pos",                   800,       pt_binning,    2400,    -1.0,     7.0);   //  3
    m2_pt_neg                  = new TH2F("m2_pt_neg",                  "m2_pt_neg",                   800,       pt_binning,    2400,    -1.0,     7.0);   //  4
    beta_p_pos                 = new TH2F("beta_p_pos",                 "beta_p_pos",                  800,       pt_binning,    3000,     0.1,     1.1);   //  5
    beta_p_neg                 = new TH2F("beta_p_neg",                 "beta_p_neg",                  800,       pt_binning,    3000,     0.1,     1.1);   //  6
    deltat_pt_pos              = new TH2F("deltat_pt_pos",              "deltat_pt_pos",               800,       pt_binning,   10000,    -1.0,    99.0);   //  7
    deltat_pt_neg              = new TH2F("deltat_pt_neg",              "deltat_pt_neg",               800,       pt_binning,   10000,    -1.0,    99.0);   //  8

    m2_pt_pos_cut              = new TH2F("m2_pt_pos_cut",              "m2_pt_pos_cut",               800,       pt_binning,    2400,    -1.0,     7.0);   //  9
    m2_pt_neg_cut              = new TH2F("m2_pt_neg_cut",              "m2_pt_neg_cut",               800,       pt_binning,    2400,    -1.0,     7.0);   // 10
    beta_p_pos_cut             = new TH2F("beta_p_pos_cut",             "beta_p_pos_cut",              800,       pt_binning,    3000,     0.1,     1.1);   // 11
    beta_p_neg_cut             = new TH2F("beta_p_neg_cut",             "beta_p_neg_cut",              800,       pt_binning,    3000,     0.1,     1.1);   // 12
    deltat_pt_pos_cut          = new TH2F("deltat_pt_pos_cut",          "deltat_pt_pos_cut",           800,       pt_binning,   10000,    -1.0,    99.0);   // 13
    deltat_pt_neg_cut          = new TH2F("deltat_pt_neg_cut",          "deltat_pt_neg_cut",           800,       pt_binning,   10000,    -1.0,    99.0);   // 14


    prot_per_event             = new TH1I("prot_per_event",             "prot_per_event",               12,        0,     12);                              // 15
//  prot_per_event_pos         = new TH1I("prot_per_event_pos",         "prot_per_event_pos",           12,        0,     12);                              // 16
//  prot_per_event_neg         = new TH1I("prot_per_event_neg",         "prot_per_event_neg",           12,        0,     12);                              // 17
    m2_pt_pos_cut_T            = new TH2F("m2_pt_pos_cut_T",            "m2_pt_pos_cut_T",             800,       pt_binning,    2400,    -1.0,     7.0);   // 18
    m2_pt_neg_cut_T            = new TH2F("m2_pt_neg_cut_T",            "m2_pt_neg_cut_T",             800,       pt_binning,    2400,    -1.0,     7.0);   // 19

    prot_phi_pt_pos            = new TH2F("prot_phi_pt_pos",            "prot_phi_pt_pos",             800,       pt_binning,     300, -1.6708,  4.8124);   // 20
    prot_phi_pt_neg            = new TH2F("prot_phi_pt_neg",            "prot_phi_pt_neg",             800,       pt_binning,     300, -1.6708,  4.8124);   // 21
    
    prot_p0_pt_pos_pos         = new TH2F("prot_p0_pt_pos_pos",         "prot_p0_pt_pos_pos",          800,       pt_binning,      76,    0.00,     3.8);   // 22
    prot_p0_pt_pos_neg         = new TH2F("prot_p0_pt_pos_neg",         "prot_p0_pt_pos_neg",          800,       pt_binning,      76,    0.00,     3.8);   // 23
    prot_p0_pt_neg_neg         = new TH2F("prot_p0_pt_neg_neg",         "prot_p0_pt_neg_neg",          800,       pt_binning,      76,    0.00,     3.8);   // 24

    di_prot_phi_pt             = new TH2F("di_prot_phi_pt",             "di_prot_phi_pt",              900,       pt_binning,     300, -1.6708,  4.8124);   // 25
    trig_05_phi_pt             = new TH2F("trig_05_phi_pt",             "trig_05_phi_pt",             1200,       pt_binning,     300, -1.6708,  4.8124);   // 26
    trig_08_phi_pt             = new TH2F("trig_08_phi_pt",             "trig_08_phi_pt",             1200,       pt_binning,     300, -1.6708,  4.8124);   // 27
    
    di_prot_dphi_p0_pos_pos_05 = new TH2F("di_prot_dphi_p0_pos_pos_05", "di_prot_dphi_p0_pos_pos_05", 1100,       pt_binning,     300, -1.6708,  4.8124);   // 28
    di_prot_dphi_p0_pos_neg_05 = new TH2F("di_prot_dphi_p0_pos_neg_05", "di_prot_dphi_p0_pos_neg_05", 1100,       pt_binning,     300, -1.6708,  4.8124);   // 29
    di_prot_dphi_p0_neg_neg_05 = new TH2F("di_prot_dphi_p0_neg_neg_05", "di_prot_dphi_p0_neg_neg_05", 1100,       pt_binning,     300, -1.6708,  4.8124);   // 30
    
    di_prot_dphi_p0_pos_pos_08 = new TH2F("di_prot_dphi_p0_pos_pos_08", "di_prot_dphi_p0_pos_pos_08", 1100,       pt_binning,     300, -1.6708,  4.8124);   // 31
    di_prot_dphi_p0_pos_neg_08 = new TH2F("di_prot_dphi_p0_pos_neg_08", "di_prot_dphi_p0_pos_neg_08", 1100,       pt_binning,     300, -1.6708,  4.8124);   // 32
    di_prot_dphi_p0_neg_neg_08 = new TH2F("di_prot_dphi_p0_neg_neg_08", "di_prot_dphi_p0_neg_neg_08", 1100,       pt_binning,     300, -1.6708,  4.8124);   // 33
    
    di_prot_dphi_r02_pos_pos_05 = new TH1F("di_prot_dphi_r02_pos_pos_05", "di_prot_dphi_r02_pos_pos_05",       300, -1.6708,  4.8124);   // 34
    di_prot_dphi_r02_pos_neg_05 = new TH1F("di_prot_dphi_r02_pos_neg_05", "di_prot_dphi_r02_pos_neg_05",       300, -1.6708,  4.8124);   // 35
    di_prot_dphi_r02_neg_neg_05 = new TH1F("di_prot_dphi_r02_neg_neg_05", "di_prot_dphi_r02_neg_neg_05",       300, -1.6708,  4.8124);   // 36
    
    di_prot_dphi_r02_pos_pos_08 = new TH1F("di_prot_dphi_r02_pos_pos_08", "di_prot_dphi_r02_pos_pos_08",       300, -1.6708,  4.8124);   // 37
    di_prot_dphi_r02_pos_neg_08 = new TH1F("di_prot_dphi_r02_pos_neg_08", "di_prot_dphi_r02_pos_neg_08",       300, -1.6708,  4.8124);   // 38
    di_prot_dphi_r02_neg_neg_08 = new TH1F("di_prot_dphi_r02_neg_neg_08", "di_prot_dphi_r02_neg_neg_08",       300, -1.6708,  4.8124);   // 39

    di_prot_dphi_r04_pos_pos_05 = new TH1F("di_prot_dphi_r04_pos_pos_05", "di_prot_dphi_r04_pos_pos_05",       300, -1.6708,  4.8124);   // 40
    di_prot_dphi_r04_pos_neg_05 = new TH1F("di_prot_dphi_r04_pos_neg_05", "di_prot_dphi_r04_pos_neg_05",       300, -1.6708,  4.8124);   // 41
    di_prot_dphi_r04_neg_neg_05 = new TH1F("di_prot_dphi_r04_neg_neg_05", "di_prot_dphi_r04_neg_neg_05",       300, -1.6708,  4.8124);   // 42
    
    di_prot_dphi_r04_pos_pos_08 = new TH1F("di_prot_dphi_r04_pos_pos_08", "di_prot_dphi_r04_pos_pos_08",       300, -1.6708,  4.8124);   // 43
    di_prot_dphi_r04_pos_neg_08 = new TH1F("di_prot_dphi_r04_pos_neg_08", "di_prot_dphi_r04_pos_neg_08",       300, -1.6708,  4.8124);   // 44
    di_prot_dphi_r04_neg_neg_08 = new TH1F("di_prot_dphi_r04_neg_neg_08", "di_prot_dphi_r04_neg_neg_08",       300, -1.6708,  4.8124);   // 45
   
    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            400,  -20.0,   20.0);                              // 46
    primary_vertex_z_cut1      = new TH1F("primary_vertex_z_cut1",      "primary_vertex_z_cut1",       400,  -20.0,   20.0);                              // 47
    primary_vertex_z_cut2      = new TH1F("primary_vertex_z_cut2",      "primary_vertex_z_cut2",       400,  -20.0,   20.0);                              // 48

//    track_cor_radius_pt        = new TH2F("track_cor_radius_pt",        "track_cor_radius_pt",         900,       pt_binning,     325,   -3.53,    3.53);   // 34
//    track_cor_radius_pt_cut    = new TH2F("track_cor_radius_pt_cut",    "track_cor_radius_pt_cut",     900,       pt_binning,     325,   -3.53,    3.53);   // 35
    
    // objects added to output file



    fOutputList->Add(fHistPt);                     //  1
    fOutputList->Add(cent_ntracks);                //  2
    
    fOutputList->Add(m2_pt_pos);                   //  3
    fOutputList->Add(m2_pt_neg);                   //  4 
    fOutputList->Add(beta_p_pos);                  //  5
    fOutputList->Add(beta_p_neg);                  //  6
    fOutputList->Add(deltat_pt_pos);               //  7
    fOutputList->Add(deltat_pt_neg);               //  8

    fOutputList->Add(m2_pt_pos_cut);               //  9
    fOutputList->Add(m2_pt_neg_cut);               // 10
    fOutputList->Add(beta_p_pos_cut);              // 11
    fOutputList->Add(beta_p_neg_cut);              // 12
    fOutputList->Add(deltat_pt_pos_cut);           // 13
    fOutputList->Add(deltat_pt_neg_cut);           // 14

    fOutputList->Add(prot_per_event);              // 15
//  fOutputList->Add(prot_per_event_pos);          // 16
//  fOutputList->Add(prot_per_event_neg);          // 17
    fOutputList->Add(m2_pt_pos_cut_T);             // 18
    fOutputList->Add(m2_pt_neg_cut_T);             // 19

    fOutputList->Add(prot_phi_pt_pos);             // 20
    fOutputList->Add(prot_phi_pt_neg);             // 21
    
    fOutputList->Add(prot_p0_pt_pos_pos);          // 22
    fOutputList->Add(prot_p0_pt_pos_neg);          // 23
    fOutputList->Add(prot_p0_pt_neg_neg);          // 24

    fOutputList->Add(di_prot_phi_pt);              // 25
    fOutputList->Add(trig_05_phi_pt);              // 26
    fOutputList->Add(trig_08_phi_pt);              // 27

    fOutputList->Add(di_prot_dphi_p0_pos_pos_05);  // 28
    fOutputList->Add(di_prot_dphi_p0_pos_neg_05);  // 29
    fOutputList->Add(di_prot_dphi_p0_neg_neg_05);  // 30
    
    fOutputList->Add(di_prot_dphi_p0_pos_pos_08);  // 31
    fOutputList->Add(di_prot_dphi_p0_pos_neg_08);  // 32
    fOutputList->Add(di_prot_dphi_p0_neg_neg_08);  // 33

    
    fOutputList->Add(di_prot_dphi_r02_pos_pos_05); // 34
    fOutputList->Add(di_prot_dphi_r02_pos_neg_05); // 35
    fOutputList->Add(di_prot_dphi_r02_neg_neg_05); // 36
    
    fOutputList->Add(di_prot_dphi_r02_pos_pos_08); // 37
    fOutputList->Add(di_prot_dphi_r02_pos_neg_08); // 38
    fOutputList->Add(di_prot_dphi_r02_neg_neg_08); // 39

    fOutputList->Add(di_prot_dphi_r04_pos_pos_05); // 40
    fOutputList->Add(di_prot_dphi_r04_pos_neg_05); // 41
    fOutputList->Add(di_prot_dphi_r04_neg_neg_05); // 42
    
    fOutputList->Add(di_prot_dphi_r04_pos_pos_08); // 43
    fOutputList->Add(di_prot_dphi_r04_pos_neg_08); // 44
    fOutputList->Add(di_prot_dphi_r04_neg_neg_08); // 45


    
    fOutputList->Add(primary_vertex_z);            // 46
    fOutputList->Add(primary_vertex_z_cut1);       // 47
    fOutputList->Add(primary_vertex_z_cut2);       // 48
    
//    fOutputList->Add(track_cor_radius_pt);         // 34
//    fOutputList->Add(track_cor_radius_pt_cut);     // 35
       
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdiprot::UserExec(Option_t *)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    Int_t iTracks(fAOD->GetNumberOfTracks());


    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0A"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    const AliAODVertex *primVertex = fAOD->GetPrimaryVertex();
    Double_t pv = primVertex->GetZ();                             primary_vertex_z->Fill(pv);
    
    if(!fAnalysisUtils->IsVertexSelected2013pA(fAOD)) return;     primary_vertex_z_cut1->Fill(pv);

    if(fAnalysisUtils->IsPileUpSPD(fAOD)) return;                 primary_vertex_z_cut2->Fill(pv);
    

    int prot_track_num[20];    
    int prot_count              = 0;
    
    int trig_05_track_num[20];
    int trig_05_track_count     = 0;

    int trig_08_track_num[20];
    int trig_08_track_count     = 0;

    Float_t pio2 = TMath::PiOver2();
    Float_t twopi  = TMath::TwoPi();
    // markA

    // loop over all these tracks
    //
    // pull out track numbers for high-pt triggers and also protons
    
    for(Int_t i(0); i < iTracks; i++)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track)                                                                      {    continue;    }
	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }
	Float_t eta = track->Eta();	if(TMath::Abs(eta) > 0.9)                       {    continue;    }
	if(!track->IsPrimaryCandidate())                                                {    continue;    }
	Double_t nsigmaTPC = 999.0;	Double_t nsigmaTOF = 999.0;
	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);     /* && trk->IsOn(AliESDtrack::kTPCpid)*/;


	Float_t pt            = track->Pt();
	fHistPt->Fill(pt);
	
	if(pt >= 5.0)	{    trig_05_track_num[trig_05_track_count] = i;	    trig_05_track_count++;	}
	if(pt >= 8.0)	{    trig_08_track_num[trig_08_track_count] = i;	    trig_08_track_count++;	}

	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
	if(!tpcIsOk  ||  !tofIsOk)                                                      {    continue;    }

	
	float   mom           = track->P();

	if(mom > 4.0)                                                                   {    continue;    }
	
	Short_t charge        = track->Charge();
	Float_t sigma_min     = -999.0;
	Float_t prot_mean     = 0.0;
	Float_t prot_sigma    = 0.0;

	Float_t phi           = track->Phi();
	   
	if(phi <  -pio2){    phi = phi + twopi; }	if(phi <  -pio2){    phi = phi + twopi;   }
	if(phi > 3*pio2){    phi = phi - twopi;	}       if(phi > 3*pio2){    phi = phi - twopi;   }





	Float_t m2tof  = get_mass_squared(track);
	Float_t deltat = tof_minus_tpion(track);
	Float_t dedx   = track->GetTPCsignal();
	
	Float_t beta      = 0.0;
	beta      = Beta(track);

	    
	if(charge > 0)
	{
	    m2_pt_pos      ->Fill(pt,  m2tof);
	    beta_p_pos     ->Fill(mom, Beta(track));
	    deltat_pt_pos  ->Fill(pt,  deltat);

	    if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340)
	    {
		if(        (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
		    )
		{
//		    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
//		    if(nsigmaTPC <= 2.0)
//		    {
			m2_pt_pos_cut    ->Fill(pt,  m2tof);
			beta_p_pos_cut   ->Fill(mom, Beta(track));
			deltat_pt_pos_cut->Fill(pt,  deltat);
			
			if(mom >= 0.5  &&  mom < 4.0)
			{
			    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
			    prot_mean = fit_prot_curve->Eval(mom);
			    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
			    prot_sigma = fit_prot_curve->Eval(mom);
			    
			    if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
			    {
				prot_track_num[prot_count] = i;
				prot_count++;
				m2_pt_pos_cut_T->Fill(pt,m2tof);

				Float_t prot_phi = track->Phi();
				if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }	if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }
				if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }	if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }
				prot_phi_pt_pos->Fill(pt, prot_phi);
			    }		    
			}			    
//		    }
		}
	    }
	}
	    
	else if(charge < 0)
	{
	    m2_pt_neg      ->Fill(pt,  m2tof);
	    beta_p_neg     ->Fill(mom, Beta(track));
	    deltat_pt_neg  ->Fill(pt,  deltat);
		
	    if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340)
	    {
		if(        (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
		    )
		{
//		    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
//		    if(nsigmaTPC <= 2.0)
//		    {	
			m2_pt_neg_cut    ->Fill(pt, m2tof);
			beta_p_neg_cut   ->Fill(mom, Beta(track));
			deltat_pt_neg_cut->Fill(pt, deltat);


			if(mom >= 0.5  &&  mom < 4.0)
			{
			    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[1][0][w]);   }
			    prot_mean = fit_prot_curve->Eval(mom);
			    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[1][1][w]);   }
			    prot_sigma = fit_prot_curve->Eval(mom);
			    
			    if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
			    {
				prot_track_num[prot_count] = i;
				prot_count++;
				m2_pt_neg_cut_T->Fill(pt,m2tof);
				
				Float_t prot_phi = track->Phi();
				if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }	if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }
				if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }	if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }
				prot_phi_pt_neg->Fill(pt, prot_phi);
			    }		    
			}
//		    }
		}
	    }
	}    //   end of neg charge if statement
    }        //   end of track loop


    prot_per_event->Fill(prot_count);
//    prot_per_event_pos->Fill(NProt_pos);
//    prot_per_event_neg->Fill(NProt_neg);


//    cout<<endl<<endl;
    // markA
    
    if(prot_count >= 2  &&  trig_05_track_count > 0)
    {
//	cout<<"trig count: "<<trig_05_track_count<<" associate proton count: "<<prot_count<<" ";
//	cout<<trig_05_track_count<<","<<prot_count<<" ";
	for(int i=0; i<trig_05_track_count; i++)  // trigger loop
	{


	    
	    int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));

	    Float_t phi_H       = trackH->Phi();
	    Float_t px_H        = trackH->Px();
	    Float_t py_H        = trackH->Py();
	    Float_t pz_H        = trackH->Pz();
	    Float_t pt_H        = trackH->Pt();
//	    Float_t eta_H       = trackH->Eta();
	    
	    if(phi_H <  -pio2){    phi_H = phi_H + twopi;   }	    if(phi_H <  -pio2){    phi_H = phi_H + twopi;   }
	    if(phi_H > 3*pio2){    phi_H = phi_H - twopi;   }	    if(phi_H > 3*pio2){    phi_H = phi_H - twopi;   }
	    trig_05_phi_pt->Fill(pt_H, phi_H);


	    for(int j=0; j<prot_count-1; j++)  // five protons means 4+3+2+1 di-protons
	    {
		int A               = prot_track_num[j];
//		cout<<"A"<<A<<" ";
		AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		
		Float_t px_A        = trackA->Px();
		Float_t py_A        = trackA->Py();
		Float_t pz_A        = trackA->Pz();
		Float_t phi_A       = trackA->Phi();
		Float_t eta_A       = trackA->Eta();
		Short_t charge_A    = trackA->Charge();

		for(int k=j+1; k<prot_count; k++)
		{
		    int B               = prot_track_num[k];
//		    cout<<"B"<<B<<" ";
		    AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));
		    
		    Float_t px_B        = trackB->Px();
		    Float_t py_B        = trackB->Py();
		    Float_t pz_B        = trackB->Pz();
		    Float_t phi_B       = trackB->Phi();
		    Float_t eta_B       = trackB->Eta();
		    Short_t charge_B    = trackB->Charge();
//		    cout<<i<<j<<k<<" ";

		    Float_t P0 = 0.00;	    P0 = pow(px_A - px_B,2)  +  pow(py_A - py_B,2)  +  pow(pz_A - pz_B,2);


		    Float_t Sx = px_A + px_B;
		    Float_t Sy = py_A + py_B;
		    Float_t PT = sqrt(pow(px_A + px_B,2) + pow(py_A + py_B,2));
		    
		    Float_t Sphi = atan2(Sy,Sx);
		    if(Sphi <  -pio2){    Sphi = Sphi + twopi;   }    if(Sphi <  -pio2){    Sphi = Sphi + twopi;   }
		    if(Sphi > 3*pio2){    Sphi = Sphi - twopi;   }    if(Sphi > 3*pio2){    Sphi = Sphi - twopi;   }
		    di_prot_phi_pt->Fill(PT, Sphi);

		    Float_t Sdphi = Sphi - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }


		    if(charge_A > 0)
		    {
			if     (charge_B > 0)
			{
			    di_prot_dphi_p0_pos_pos_05->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_pos_pos_05->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_pos_pos_05->Fill(Sdphi);	}
			}
			else if(charge_B < 0)
			{
			    di_prot_dphi_p0_pos_neg_05->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_pos_neg_05->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_pos_neg_05->Fill(Sdphi);	}			    
			}

		    }
		    else if(charge_A < 0)
		    {
			if     (charge_B > 0)
			{
			    di_prot_dphi_p0_pos_neg_05->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_pos_neg_05->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_pos_neg_05->Fill(Sdphi);	}			    
			}
			else if(charge_B < 0)
			{
			    di_prot_dphi_p0_neg_neg_05->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_neg_neg_05->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_neg_neg_05->Fill(Sdphi);	}  
			}
		    }

		    if(charge_A > 0  &&  i == 0)
		    {
			if     (charge_B > 0){    prot_p0_pt_pos_pos->Fill(PT, P0);    	}
			else if(charge_B < 0){    prot_p0_pt_pos_neg->Fill(PT, P0);    	}
		    }
		    else if(charge_A < 0  &&  i == 0)
		    {
			if     (charge_B > 0){    prot_p0_pt_pos_neg->Fill(PT, P0);     }
			else if(charge_B < 0){    prot_p0_pt_neg_neg->Fill(PT, P0);    	}
		    }
		    
		}	
	    }	    
	}
//	cout<<endl;
    }




    if(prot_count >= 2  &&  trig_08_track_count > 0)
    {
//	cout<<"trig count: "<<trig_08_track_count<<" associate proton count: "<<prot_count<<" ";
//	cout<<trig_08_track_count<<"88,"<<prot_count<<" ";
	for(int i=0; i<trig_08_track_count; i++)  // trigger loop
	{


//	    file_output<<trig_08_track_count<<endl;


	    
	    int H               = trig_08_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));

	    Float_t phi_H       = trackH->Phi();
	    Float_t px_H        = trackH->Px();
	    Float_t py_H        = trackH->Py();
	    Float_t pz_H        = trackH->Pz();
	    Float_t pt_H        = trackH->Pt();
	    
	    if(phi_H <  -pio2){    phi_H = phi_H + twopi;   }	    if(phi_H <  -pio2){    phi_H = phi_H + twopi;   }
	    if(phi_H > 3*pio2){    phi_H = phi_H - twopi;   }	    if(phi_H > 3*pio2){    phi_H = phi_H - twopi;   }
	    trig_08_phi_pt->Fill(pt_H, phi_H);


	    for(int j=0; j<prot_count-1; j++)  // five protons means 4+3+2+1 di-protons
	    {
		
		int A               = prot_track_num[j];
		AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		
		Float_t px_A        = trackA->Px();
		Float_t py_A        = trackA->Py();
		Float_t pz_A        = trackA->Pz();
		Float_t phi_A       = trackA->Phi();
		Float_t eta_A       = trackA->Eta();
		Short_t charge_A    = trackA->Charge();

		for(int k=j+1; k<prot_count; k++)
		{
		    int B               = prot_track_num[k];
		    AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));
		    
		    Float_t px_B        = trackB->Px();
		    Float_t py_B        = trackB->Py();
		    Float_t pz_B        = trackB->Pz();
		    Float_t phi_B       = trackB->Phi();
		    Float_t eta_B       = trackB->Eta();
		    Short_t charge_B    = trackB->Charge();
//		    cout<<i<<j<<k<<" ";

		    Float_t P0 = 0.00;	    P0 = pow(px_A - px_B,2)  +  pow(py_A - py_B,2)  +  pow(pz_A - pz_B,2);


		    Float_t Sx = px_A + px_B;
		    Float_t Sy = py_A + py_B;
		    Float_t PT = sqrt(pow(px_A + px_B,2) + pow(py_A + py_B,2));
		    
		    Float_t Sphi = atan2(Sy,Sx);
		    if(Sphi <  -pio2){    Sphi = Sphi + twopi;   }    if(Sphi <  -pio2){    Sphi = Sphi + twopi;   }
		    if(Sphi > 3*pio2){    Sphi = Sphi - twopi;   }    if(Sphi > 3*pio2){    Sphi = Sphi - twopi;   }


		    Float_t Sdphi = Sphi - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }


		    if(charge_A > 0)
		    {
			if     (charge_B > 0)
			{
			    di_prot_dphi_p0_pos_pos_08->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_pos_pos_08->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_pos_pos_08->Fill(Sdphi);	}
			}
			else if(charge_B < 0)
			{
			    di_prot_dphi_p0_pos_neg_08->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_pos_neg_08->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_pos_neg_08->Fill(Sdphi);	}
			}
		    }
		    else if(charge_A < 0)
		    {
			if     (charge_B > 0)
			{
			    di_prot_dphi_p0_pos_neg_08->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_pos_neg_08->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_pos_neg_08->Fill(Sdphi);	}
			}
			else if(charge_B < 0)
			{
			    di_prot_dphi_p0_neg_neg_08->Fill(P0, Sdphi);
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.2){    di_prot_dphi_r02_neg_neg_08->Fill(Sdphi);	}
			    if(sqrt(pow(eta_A - eta_B,2) + pow(phi_A - phi_B,2)) < 0.4){    di_prot_dphi_r04_neg_neg_08->Fill(Sdphi);	}
			}
		    }

		    
		}
	
	    }

	    
	}

//	cout<<endl;
    }
//    cout<<endl<<endl;





	
/*
	int A = prot_track_num_1;
	int B = prot_track_num_2;
	
	AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));

	if(!trackA  ||  !trackB  ||  !trackH) {}
	else
	{
	    Float_t prot_phi_A      = trackA->Phi();
	    Float_t prot_pt_A       = trackA->Pt();
	    Float_t prot_mom_A      = trackA->P();
	    Short_t prot_charge_A   = trackA->Charge();
	    Float_t prot_eta_A      = trackA->Eta();

	    Float_t prot_et_A       = 0.0;
	    Float_t prot_e_A        = 0.0;
	    Float_t prot_px_A       = trackA->Px();
	    Float_t prot_py_A       = trackA->Py();
	    Float_t prot_pz_A       = trackA->Pz();
	    
	    prot_et_A = sqrt(pow(prot_pt_A ,2) + pow(0.93827,2));
	    prot_e_A  = sqrt(pow(prot_mom_A,2) + pow(0.93827,2));
	    	    
	    if(prot_phi_A <  -pio2)    prot_phi_A = prot_phi_A + twopi;
	    if(prot_phi_A <  -pio2)    prot_phi_A = prot_phi_A + twopi;
	    if(prot_phi_A > 3*pio2)    prot_phi_A = prot_phi_A - twopi;
	    if(prot_phi_A > 3*pio2)    prot_phi_A = prot_phi_A - twopi;
	    prot_phi_pt_A->Fill(prot_phi_A);
	    
	    
	    Float_t prot_phi_B      = trackB->Phi();
	    Float_t prot_pt_B       = trackB->Pt();
	    Float_t prot_mom_B      = trackB->P();
	    Short_t  prot_charge_B   = trackB->Charge();
	    Float_t prot_eta_B      = trackB->Eta();

	    Float_t prot_et_B       = 0.0;
	    Float_t prot_e_B        = 0.0;
	    Float_t prot_px_B       = trackB->Px();
	    Float_t prot_py_B       = trackB->Py();
	    Float_t prot_pz_B       = trackB->Pz();
	    
	    prot_et_B = sqrt(pow(prot_pt_B ,2) + pow(0.93827,2));
	    prot_e_B  = sqrt(pow(prot_mom_B,2) + pow(0.93827,2));

	    if(prot_phi_B <  -pio2)    prot_phi_B = prot_phi_B + twopi;
	    if(prot_phi_B <  -pio2)    prot_phi_B = prot_phi_B + twopi;
	    if(prot_phi_B > 3*pio2)    prot_phi_B = prot_phi_B - twopi;
	    if(prot_phi_B > 3*pio2)    prot_phi_B = prot_phi_B - twopi;
	    prot_phi_pt_B->Fill(prot_phi_B);
	    	    


	    
	    Float_t P0 = 0.00;	    P0 = pow(prot_px_A - prot_px_B,2)  +  pow(prot_py_A - prot_py_B,2)  +  pow(prot_pz_A - prot_pz_B,2);

	    Float_t Q  = 0.00;	    Q  = sqrt(P0);
	    Float_t Sx = prot_px_A + prot_px_B;
	    Float_t Sy = prot_py_A + prot_py_B;

	    Float_t Sphi = atan2(Sy,Sx);
	    if(Sphi <  -pio2)    Sphi = Sphi + twopi;
	    if(Sphi <  -pio2)    Sphi = Sphi + twopi;
	    if(Sphi > 3*pio2)    Sphi = Sphi - twopi;
	    if(Sphi > 3*pio2)    Sphi = Sphi - twopi;
				
	    di_prot_phi_pt->Fill(Sphi);

	    Float_t Sdphi = Sphi - phi_H;
	    if(Sdphi <  -pio2)    Sdphi = Sdphi + twopi;
	    if(Sdphi <  -pio2)    Sdphi = Sdphi + twopi;
	    if(Sdphi > 3*pio2)    Sdphi = Sdphi - twopi;
	    if(Sdphi > 3*pio2)    Sdphi = Sdphi - twopi;



	    

	}

*/
    


//    file_output.close();
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file


}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdiprot::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFdiprot::Beta(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFdiprot::tof_minus_tpion(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFdiprot::get_mass_squared(AliAODTrack *track)
{
    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
    Double_t stop_time      = track->GetTOFsignal();
    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
    Double_t tof            = stop_time - start_time;
    tof                     = tof * 0.001;                                                                  // convert from ps to ns
    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
    length                  = length * 100;                                                                 // convert to cm
    Double_t mom            = track->P();
    Double_t m2             = 0.0;
    Double_t c2             = 29.9792458;
    
    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);
    return m2;

}
