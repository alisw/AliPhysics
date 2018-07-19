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
#include <bitset>
#include <iomanip>
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TMath.h"
//#include "TGraphErrors.h"
#include "TList.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCorPIDTOFdeut.h"
//#include "AliPIDResponse.h"
//#include "AliAnalysisUtils.h"

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

//  event trigger selection
//  if(datajets)
//  {
//      if(!(fInputHandler->IsEventSelected() & fTriggerSelectionBits)) return false;
//      if(fTriggerSelectionString.Length())
//      {
//	  if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerSelectionString)) return false;
//      }
//  }



using namespace std;            // std namespace: so you can do things like 'cout'
//using namespace BSchaefer_devel;

//ofstream file_output("output.txt");


ClassImp(AliAnalysisTaskCorPIDTOFdeut) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFdeut::AliAnalysisTaskCorPIDTOFdeut() : AliAnalysisTaskSE(), 
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

    m2_pt_pos_cut_T(0),            // 15
    m2_pt_neg_cut_T(0),            // 16
    m2_pt_pos_cut_G(0),            // 17
    m2_pt_neg_cut_G(0),            // 18
    m2_pt_pos_cut_A(0),            // 19
    m2_pt_neg_cut_A(0),            // 20
    m2_pt_pos_cut_B(0),            // 21
    m2_pt_neg_cut_B(0),            // 22
    
    deut_phi_pt_pos_T(0),          // 23
    deut_phi_pt_neg_T(0),          // 24
    deut_phi_pt_pos_A(0),          // 25
    deut_phi_pt_neg_A(0),          // 26
    deut_phi_pt_pos_B(0),          // 27
    deut_phi_pt_neg_B(0),          // 28

    deut_per_event(0),             // 29
    trig_03_per_event(0),          // 30
    trig_05_per_event(0),          // 31
    trig_08_per_event(0),          // 32
    
    trig_03_phi_pt_pos(0),         // 33
    trig_03_phi_pt_neg(0),         // 34    
    trig_05_phi_pt_pos(0),         // 35
    trig_05_phi_pt_neg(0),         // 36
    trig_08_phi_pt_pos(0),         // 37
    trig_08_phi_pt_neg(0),         // 38
    
    tof_phi_eta_pos(0),            // 39
    tof_phi_eta_neg(0),            // 40
    tof_phi_eta_pos_deut(0),       // 41
    tof_phi_eta_neg_deut(0),       // 42

    deut_pt_compare_pos(0),        // 43
    deut_pt_compare_neg(0),        // 44
    tpc_sector_fraction(0),        // 45
    primary_vertex_z(0),           // 46
    primary_vertex_z_cut(0),       // 47
    
    deut_dphi_pt_pos_pos_03_T(0),  // 48
    deut_dphi_pt_pos_neg_03_T(0),  // 49
    deut_dphi_pt_neg_neg_03_T(0),  // 50
    deut_dphi_pt_pos_pos_05_T(0),  // 51
    deut_dphi_pt_pos_neg_05_T(0),  // 52
    deut_dphi_pt_neg_neg_05_T(0),  // 53
    deut_dphi_pt_pos_pos_08_T(0),  // 54
    deut_dphi_pt_pos_neg_08_T(0),  // 55
    deut_dphi_pt_neg_neg_08_T(0),  // 56

    deut_dphi_pt_pos_pos_03_A(0),  // 57
    deut_dphi_pt_pos_neg_03_A(0),  // 58
    deut_dphi_pt_neg_neg_03_A(0),  // 59
    deut_dphi_pt_pos_pos_05_A(0),  // 60
    deut_dphi_pt_pos_neg_05_A(0),  // 61
    deut_dphi_pt_neg_neg_05_A(0),  // 62
    deut_dphi_pt_pos_pos_08_A(0),  // 63
    deut_dphi_pt_pos_neg_08_A(0),  // 64
    deut_dphi_pt_neg_neg_08_A(0),  // 65
 
    deut_dphi_pt_pos_pos_03_B(0),  // 66
    deut_dphi_pt_pos_neg_03_B(0),  // 67
    deut_dphi_pt_neg_neg_03_B(0),  // 68
    deut_dphi_pt_pos_pos_05_B(0),  // 69
    deut_dphi_pt_pos_neg_05_B(0),  // 70
    deut_dphi_pt_neg_neg_05_B(0),  // 71
    deut_dphi_pt_pos_pos_08_B(0),  // 72
    deut_dphi_pt_pos_neg_08_B(0),  // 73
    deut_dphi_pt_neg_neg_08_B(0),  // 74

    DCAxy_pos(0),                  // 75
    DCAxy_neg(0),                  // 76
    DCAz_pos(0),                   // 77
    DCAz_neg(0),                   // 78

    m2_pt_pos_fine(0),             // 79
    m2_pt_neg_fine(0),             // 80

    m2_pt_pos_cut_fine(0),         // 81
    m2_pt_neg_cut_fine(0),         // 82

    m2_pt_pos_cut_T_fine(0),       // 83
    m2_pt_neg_cut_T_fine(0)        // 84    
//    deut_eta_vs_phi_pos(0),        // 50
//    deut_eta_vs_phi_neg(0)         // 51

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFdeut::AliAnalysisTaskCorPIDTOFdeut(const char* name) : AliAnalysisTaskSE(name),
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

    m2_pt_pos_cut_T(0),            // 15
    m2_pt_neg_cut_T(0),            // 16
    m2_pt_pos_cut_G(0),            // 17
    m2_pt_neg_cut_G(0),            // 18
    m2_pt_pos_cut_A(0),            // 19
    m2_pt_neg_cut_A(0),            // 20
    m2_pt_pos_cut_B(0),            // 21
    m2_pt_neg_cut_B(0),            // 22
    
    deut_phi_pt_pos_T(0),          // 23
    deut_phi_pt_neg_T(0),          // 24
    deut_phi_pt_pos_A(0),          // 25
    deut_phi_pt_neg_A(0),          // 26
    deut_phi_pt_pos_B(0),          // 27
    deut_phi_pt_neg_B(0),          // 28

    deut_per_event(0),             // 29
    trig_03_per_event(0),          // 30
    trig_05_per_event(0),          // 31
    trig_08_per_event(0),          // 32
    
    trig_03_phi_pt_pos(0),         // 33
    trig_03_phi_pt_neg(0),         // 34    
    trig_05_phi_pt_pos(0),         // 35
    trig_05_phi_pt_neg(0),         // 36
    trig_08_phi_pt_pos(0),         // 37
    trig_08_phi_pt_neg(0),         // 38
    
    tof_phi_eta_pos(0),            // 39
    tof_phi_eta_neg(0),            // 40
    tof_phi_eta_pos_deut(0),       // 41
    tof_phi_eta_neg_deut(0),       // 42

    deut_pt_compare_pos(0),        // 43
    deut_pt_compare_neg(0),        // 44
    tpc_sector_fraction(0),        // 45
    primary_vertex_z(0),           // 46
    primary_vertex_z_cut(0),       // 47
    
    deut_dphi_pt_pos_pos_03_T(0),  // 48
    deut_dphi_pt_pos_neg_03_T(0),  // 49
    deut_dphi_pt_neg_neg_03_T(0),  // 50
    deut_dphi_pt_pos_pos_05_T(0),  // 51
    deut_dphi_pt_pos_neg_05_T(0),  // 52
    deut_dphi_pt_neg_neg_05_T(0),  // 53
    deut_dphi_pt_pos_pos_08_T(0),  // 54
    deut_dphi_pt_pos_neg_08_T(0),  // 55
    deut_dphi_pt_neg_neg_08_T(0),  // 56

    deut_dphi_pt_pos_pos_03_A(0),  // 57
    deut_dphi_pt_pos_neg_03_A(0),  // 58
    deut_dphi_pt_neg_neg_03_A(0),  // 59
    deut_dphi_pt_pos_pos_05_A(0),  // 60
    deut_dphi_pt_pos_neg_05_A(0),  // 61
    deut_dphi_pt_neg_neg_05_A(0),  // 62
    deut_dphi_pt_pos_pos_08_A(0),  // 63
    deut_dphi_pt_pos_neg_08_A(0),  // 64
    deut_dphi_pt_neg_neg_08_A(0),  // 65
 
    deut_dphi_pt_pos_pos_03_B(0),  // 66
    deut_dphi_pt_pos_neg_03_B(0),  // 67
    deut_dphi_pt_neg_neg_03_B(0),  // 68
    deut_dphi_pt_pos_pos_05_B(0),  // 69
    deut_dphi_pt_pos_neg_05_B(0),  // 70
    deut_dphi_pt_neg_neg_05_B(0),  // 71
    deut_dphi_pt_pos_pos_08_B(0),  // 72
    deut_dphi_pt_pos_neg_08_B(0),  // 73
    deut_dphi_pt_neg_neg_08_B(0),   // 74									   

    DCAxy_pos(0),                  // 75
    DCAxy_neg(0),                  // 76
    DCAz_pos(0),                   // 77
    DCAz_neg(0),                   // 78

    m2_pt_pos_fine(0),             // 79
    m2_pt_neg_fine(0),             // 80

    m2_pt_pos_cut_fine(0),         // 81
    m2_pt_neg_cut_fine(0),         // 82

    m2_pt_pos_cut_T_fine(0),       // 83
    m2_pt_neg_cut_T_fine(0)        // 84

//    deut_eta_vs_phi_pos(0),        // 50
//    deut_eta_vs_phi_neg(0)         // 51
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFdeut::~AliAnalysisTaskCorPIDTOFdeut()
{
    // destructor
    if(fAnalysisUtils) delete fAnalysisUtils;
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdeut::UserCreateOutputObjects()
{
    // 2.88465 0.0761582 0.709281 0.124386 0.017642 -0.0316078 2.65738 0.115151 0.918566 0.0986592 0.0187545 0.00346519    // pp 2016 untriggered

    
    
    deut_curves[0][0][0] = 2.88465;     // pos deut mean curve
    deut_curves[0][0][1] = 0.0761582;
    deut_curves[0][0][2] = 0.709281;

    deut_curves[0][1][0] = 0.124386;  // pos deut sigma curve
    deut_curves[0][1][1] = 0.017642;
    deut_curves[0][1][2] = -0.0316078;

    deut_curves[1][0][0] = 2.65738;     // neg deut mean curve
    deut_curves[1][0][1] = 0.115151;
    deut_curves[1][0][2] = 0.918566;

    deut_curves[1][1][0] = 0.0986592;  // neg deut sigma curve
    deut_curves[1][1][1] = 0.0187545;
    deut_curves[1][1][2] = 0.00346519;

    fAnalysisUtils = new AliAnalysisUtils;

    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

/*
    Double_t pt_binning[2001];
    Float_t moving_marker = 0.10;
    for(int i=0; i<1202; i++)
    {
	pt_binning[i] = moving_marker;
	moving_marker = moving_marker + pt_binning[i] * 0.005;
    }
*/  
    float lower = -pio2;
    float upper = 3.0*pio2;

    Double_t pt_binning[2001];
    Float_t moving_marker = 0.10;
    for(int i=0; i<1202; i++)
    {
	pt_binning[i] = moving_marker;
	moving_marker = moving_marker + pt_binning[i] * 0.005;
    }
    
    
    fHistPt                    = new TH1F("fHistPt",                    "Pt()",                         50,  0.0, 5.0);                              //  1
    cent_ntracks               = new TH2F("cent_ntracks",               "cent_ntracks",                100,    0, 100,     100,       0,     800);   //  2
    
    m2_pt_pos                  = new TH2F("m2_pt_pos",                  "m2_pt_pos",                    50,  0.0, 5.0,    2400,    -1.0,     7.0);   //  3
    m2_pt_neg                  = new TH2F("m2_pt_neg",                  "m2_pt_neg",                    50,  0.0, 5.0,    2400,    -1.0,     7.0);   //  4
    beta_p_pos                 = new TH2F("beta_p_pos",                 "beta_p_pos",                   50,  0.0, 5.0,    3000,     0.1,     1.1);   //  5
    beta_p_neg                 = new TH2F("beta_p_neg",                 "beta_p_neg",                   50,  0.0, 5.0,    3000,     0.1,     1.1);   //  6
    deltat_pt_pos              = new TH2F("deltat_pt_pos",              "deltat_pt_pos",                50,  0.0, 5.0,    2500,    -1.0,    24.0);   //  7
    deltat_pt_neg              = new TH2F("deltat_pt_neg",              "deltat_pt_neg",                50,  0.0, 5.0,    2500,    -1.0,    24.0);   //  8
    m2_pt_pos_cut              = new TH2F("m2_pt_pos_cut",              "m2_pt_pos_cut",                50,  0.0, 5.0,    2400,    -1.0,     7.0);   //  9
    m2_pt_neg_cut              = new TH2F("m2_pt_neg_cut",              "m2_pt_neg_cut",                50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 10
    beta_p_pos_cut             = new TH2F("beta_p_pos_cut",             "beta_p_pos_cut",               50,  0.0, 5.0,    3000,     0.1,     1.1);   // 11
    beta_p_neg_cut             = new TH2F("beta_p_neg_cut",             "beta_p_neg_cut",               50,  0.0, 5.0,    3000,     0.1,     1.1);   // 12
    deltat_pt_pos_cut          = new TH2F("deltat_pt_pos_cut",          "deltat_pt_pos_cut",            50,  0.0, 5.0,    2500,    -1.0,    24.0);   // 13
    deltat_pt_neg_cut          = new TH2F("deltat_pt_neg_cut",          "deltat_pt_neg_cut",            50,  0.0, 5.0,    2500,    -1.0,    24.0);   // 14

    m2_pt_pos_cut_T            = new TH2F("m2_pt_pos_cut_T",            "m2_pt_pos_cut_T",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 15
    m2_pt_neg_cut_T            = new TH2F("m2_pt_neg_cut_T",            "m2_pt_neg_cut_T",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 16
    m2_pt_pos_cut_G            = new TH2F("m2_pt_pos_cut_G",            "m2_pt_pos_cut_G",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 17
    m2_pt_neg_cut_G            = new TH2F("m2_pt_neg_cut_G",            "m2_pt_neg_cut_G",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 18
    m2_pt_pos_cut_A            = new TH2F("m2_pt_pos_cut_A",            "m2_pt_pos_cut_A",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 19
    m2_pt_neg_cut_A            = new TH2F("m2_pt_neg_cut_A",            "m2_pt_neg_cut_A",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 20
    m2_pt_pos_cut_B            = new TH2F("m2_pt_pos_cut_B",            "m2_pt_pos_cut_B",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 21
    m2_pt_neg_cut_B            = new TH2F("m2_pt_neg_cut_B",            "m2_pt_neg_cut_B",              50,  0.0, 5.0,    2400,    -1.0,     7.0);   // 22

    deut_phi_pt_pos_T          = new TH2F("deut_phi_pt_pos_T",          "deut_phi_pt_pos_T",            50,  0.0, 5.0,     288,   lower,   upper);   // 23
    deut_phi_pt_neg_T          = new TH2F("deut_phi_pt_neg_T",          "deut_phi_pt_neg_T",            50,  0.0, 5.0,     288,   lower,   upper);   // 24
    deut_phi_pt_pos_A          = new TH2F("deut_phi_pt_pos_A",          "deut_phi_pt_pos_A",            50,  0.0, 5.0,     288,   lower,   upper);   // 25
    deut_phi_pt_neg_A          = new TH2F("deut_phi_pt_neg_A",          "deut_phi_pt_neg_A",            50,  0.0, 5.0,     288,   lower,   upper);   // 26
    deut_phi_pt_pos_B          = new TH2F("deut_phi_pt_pos_B",          "deut_phi_pt_pos_B",            50,  0.0, 5.0,     288,   lower,   upper);   // 27
    deut_phi_pt_neg_B          = new TH2F("deut_phi_pt_neg_B",          "deut_phi_pt_neg_B",            50,  0.0, 5.0,     288,   lower,   upper);   // 28

    deut_per_event             = new TH1I("deut_per_event",             "deut_per_event",               12,    0,  12);                              // 29
    trig_03_per_event          = new TH1I("trig_03_per_event",          "trig_03_per_event",            20,    0,  20);                              // 30
    trig_05_per_event          = new TH1I("trig_05_per_event",          "trig_05_per_event",            20,    0,  20);                              // 31
    trig_08_per_event          = new TH1I("trig_08_per_event",          "trig_08_per_event",            20,    0,  20);                              // 32

    trig_03_phi_pt_pos         = new TH2F("trig_03_phi_pt_pos",         "trig_03_phi_pt_pos",          170,  3.0, 20.0,    288,   lower,   upper);   // 33
    trig_03_phi_pt_neg         = new TH2F("trig_03_phi_pt_neg",         "trig_03_phi_pt_neg",          170,  3.0, 20.0,    288,   lower,   upper);   // 34
    trig_05_phi_pt_pos         = new TH2F("trig_05_phi_pt_pos",         "trig_05_phi_pt_pos",          170,  3.0, 20.0,    288,   lower,   upper);   // 35
    trig_05_phi_pt_neg         = new TH2F("trig_05_phi_pt_neg",         "trig_05_phi_pt_neg",          170,  3.0, 20.0,    288,   lower,   upper);   // 36
    trig_08_phi_pt_pos         = new TH2F("trig_08_phi_pt_pos",         "trig_08_phi_pt_pos",          170,  3.0, 20.0,    288,   lower,   upper);   // 37
    trig_08_phi_pt_neg         = new TH2F("trig_08_phi_pt_neg",         "trig_08_phi_pt_neg",          170,  3.0, 20.0,    288,   lower,   upper);   // 38

    tof_phi_eta_pos            = new TH2F("tof_phi_eta_pos",            "tof_phi_eta_pos",             190,lower,upper,     60,   -0.92,    0.92);   // 39
    tof_phi_eta_neg            = new TH2F("tof_phi_eta_neg",            "tof_phi_eta_neg",             190,lower,upper,     60,   -0.92,    0.92);   // 40
    tof_phi_eta_pos_deut       = new TH2F("tof_phi_eta_pos_deut",       "tof_phi_eta_pos_deut",        190,lower,upper,     60,   -0.92,    0.92);   // 41
    tof_phi_eta_neg_deut       = new TH2F("tof_phi_eta_neg_deut",       "tof_phi_eta_neg_deut",        190,lower,upper,     60,   -0.92,    0.92);   // 42

    deut_pt_compare_pos        = new TH1F("deut_pt_compare_pos",        "deut_pt_compare_pos",         100, -0.2,  0.2);                             // 43
    deut_pt_compare_neg        = new TH1F("deut_pt_compare_neg",        "deut_pt_compare_neg",         100, -0.2,  0.2);                             // 44
    tpc_sector_fraction        = new TH1F("tpc_sector_fraction",        "tpc_sector_fraction",         100,  0.0,  1.0);                             // 45
    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            400,-20.0, 20.0);                             // 46
    primary_vertex_z_cut       = new TH1F("primary_vertex_z_cut",       "primary_vertex_z_cut",        400,-20.0, 20.0);                             // 47

    deut_dphi_pt_pos_pos_03_T  = new TH2F("deut_dphi_pt_pos_pos_03_T",  "deut_dphi_pt_pos_pos_03_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 48
    deut_dphi_pt_pos_neg_03_T  = new TH2F("deut_dphi_pt_pos_neg_03_T",  "deut_dphi_pt_pos_neg_03_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 49
    deut_dphi_pt_neg_neg_03_T  = new TH2F("deut_dphi_pt_neg_neg_03_T",  "deut_dphi_pt_neg_neg_03_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 50
    deut_dphi_pt_pos_pos_05_T  = new TH2F("deut_dphi_pt_pos_pos_05_T",  "deut_dphi_pt_pos_pos_05_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 51
    deut_dphi_pt_pos_neg_05_T  = new TH2F("deut_dphi_pt_pos_neg_05_T",  "deut_dphi_pt_pos_neg_05_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 52
    deut_dphi_pt_neg_neg_05_T  = new TH2F("deut_dphi_pt_neg_neg_05_T",  "deut_dphi_pt_neg_neg_05_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 53
    deut_dphi_pt_pos_pos_08_T  = new TH2F("deut_dphi_pt_pos_pos_08_T",  "deut_dphi_pt_pos_pos_08_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 54
    deut_dphi_pt_pos_neg_08_T  = new TH2F("deut_dphi_pt_pos_neg_08_T",  "deut_dphi_pt_pos_neg_08_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 55
    deut_dphi_pt_neg_neg_08_T  = new TH2F("deut_dphi_pt_neg_neg_08_T",  "deut_dphi_pt_neg_neg_08_T",    50,  0.0, 5.0,     288,   lower,   upper);   // 56
    
    deut_dphi_pt_pos_pos_03_A  = new TH2F("deut_dphi_pt_pos_pos_03_A",  "deut_dphi_pt_pos_pos_03_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 57
    deut_dphi_pt_pos_neg_03_A  = new TH2F("deut_dphi_pt_pos_neg_03_A",  "deut_dphi_pt_pos_neg_03_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 58
    deut_dphi_pt_neg_neg_03_A  = new TH2F("deut_dphi_pt_neg_neg_03_A",  "deut_dphi_pt_neg_neg_03_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 59
    deut_dphi_pt_pos_pos_05_A  = new TH2F("deut_dphi_pt_pos_pos_05_A",  "deut_dphi_pt_pos_pos_05_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 60
    deut_dphi_pt_pos_neg_05_A  = new TH2F("deut_dphi_pt_pos_neg_05_A",  "deut_dphi_pt_pos_neg_05_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 61
    deut_dphi_pt_neg_neg_05_A  = new TH2F("deut_dphi_pt_neg_neg_05_A",  "deut_dphi_pt_neg_neg_05_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 62
    deut_dphi_pt_pos_pos_08_A  = new TH2F("deut_dphi_pt_pos_pos_08_A",  "deut_dphi_pt_pos_pos_08_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 63
    deut_dphi_pt_pos_neg_08_A  = new TH2F("deut_dphi_pt_pos_neg_08_A",  "deut_dphi_pt_pos_neg_08_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 64
    deut_dphi_pt_neg_neg_08_A  = new TH2F("deut_dphi_pt_neg_neg_08_A",  "deut_dphi_pt_neg_neg_08_A",    50,  0.0, 5.0,     288,   lower,   upper);   // 65
    
    deut_dphi_pt_pos_pos_03_B  = new TH2F("deut_dphi_pt_pos_pos_03_B",  "deut_dphi_pt_pos_pos_03_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 66
    deut_dphi_pt_pos_neg_03_B  = new TH2F("deut_dphi_pt_pos_neg_03_B",  "deut_dphi_pt_pos_neg_03_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 67
    deut_dphi_pt_neg_neg_03_B  = new TH2F("deut_dphi_pt_neg_neg_03_B",  "deut_dphi_pt_neg_neg_03_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 68
    deut_dphi_pt_pos_pos_05_B  = new TH2F("deut_dphi_pt_pos_pos_05_B",  "deut_dphi_pt_pos_pos_05_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 69
    deut_dphi_pt_pos_neg_05_B  = new TH2F("deut_dphi_pt_pos_neg_05_B",  "deut_dphi_pt_pos_neg_05_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 70
    deut_dphi_pt_neg_neg_05_B  = new TH2F("deut_dphi_pt_neg_neg_05_B",  "deut_dphi_pt_neg_neg_05_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 71
    deut_dphi_pt_pos_pos_08_B  = new TH2F("deut_dphi_pt_pos_pos_08_B",  "deut_dphi_pt_pos_pos_08_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 72
    deut_dphi_pt_pos_neg_08_B  = new TH2F("deut_dphi_pt_pos_neg_08_B",  "deut_dphi_pt_pos_neg_08_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 73
    deut_dphi_pt_neg_neg_08_B  = new TH2F("deut_dphi_pt_neg_neg_08_B",  "deut_dphi_pt_neg_neg_08_B",    50,  0.0, 5.0,     288,   lower,   upper);   // 74


    DCAxy_pos                  = new TH1F("DCAxy_pos",                  "DCAxy_pos",                   100, -3.0, 3.0);                              // 75
    DCAxy_neg                  = new TH1F("DCAxy_neg",                  "DCAxy_neg",                   100, -3.0, 3.0);                              // 76
    DCAz_pos                   = new TH1F("DCAz_pos",                   "DCAz_pos",                    100, -3.0, 3.0);                              // 77
    DCAz_neg                   = new TH1F("DCAz_neg",                   "DCAz_neg",                    100, -3.0, 3.0);                              // 78

    m2_pt_pos_fine             = new TH2F("m2_pt_pos_fine",             "m2_pt_pos_fine",              800,pt_binning,    2400,    -1.0,     7.0);   // 79
    m2_pt_neg_fine             = new TH2F("m2_pt_neg_fine",             "m2_pt_neg_fine",              800,pt_binning,    2400,    -1.0,     7.0);   // 80
    m2_pt_pos_cut_fine         = new TH2F("m2_pt_pos_cut_fine",         "m2_pt_pos_cut_fine",          800,pt_binning,    2400,    -1.0,     7.0);   // 79
    m2_pt_neg_cut_fine         = new TH2F("m2_pt_neg_cut_fine",         "m2_pt_neg_cut_fine",          800,pt_binning,    2400,    -1.0,     7.0);   // 80
    m2_pt_pos_cut_T_fine       = new TH2F("m2_pt_pos_cut_T_fine",       "m2_pt_pos_cut_T_fine",        800,pt_binning,    2400,    -1.0,     7.0);   // 79
    m2_pt_neg_cut_T_fine       = new TH2F("m2_pt_neg_cut_T_fine",       "m2_pt_neg_cut_T_fine",        800,pt_binning,    2400,    -1.0,     7.0);   // 80    
//  deut_eta_vs_phi_pos        = new TH2F("deut_eta_vs_phi_pos",        "deut_eta_vs_phi_pos",         300,    lower,  upper,     150,   -0.82,    0.82);   // 50
//  deut_eta_vs_phi_neg        = new TH2F("deut_eta_vs_phi_neg",        "deut_eta_vs_phi_neg",         300,    lower,  upper,     150,   -0.82,    0.82);   // 50

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
    
    fOutputList->Add(m2_pt_pos_cut_T);             // 15
    fOutputList->Add(m2_pt_neg_cut_T);             // 16
    fOutputList->Add(m2_pt_pos_cut_G);             // 17
    fOutputList->Add(m2_pt_neg_cut_G);             // 18
    fOutputList->Add(m2_pt_pos_cut_A);             // 19
    fOutputList->Add(m2_pt_neg_cut_A);             // 20
    fOutputList->Add(m2_pt_pos_cut_B);             // 21
    fOutputList->Add(m2_pt_neg_cut_B);             // 22
    
    fOutputList->Add(deut_phi_pt_pos_T);           // 23
    fOutputList->Add(deut_phi_pt_neg_T);           // 24
    fOutputList->Add(deut_phi_pt_pos_A);           // 25
    fOutputList->Add(deut_phi_pt_neg_A);           // 26
    fOutputList->Add(deut_phi_pt_pos_B);           // 27
    fOutputList->Add(deut_phi_pt_neg_B);           // 28
    
    fOutputList->Add(deut_per_event);              // 29
    fOutputList->Add(trig_03_per_event);           // 30
    fOutputList->Add(trig_05_per_event);           // 31
    fOutputList->Add(trig_08_per_event);           // 32
    
    fOutputList->Add(trig_03_phi_pt_pos);          // 33
    fOutputList->Add(trig_03_phi_pt_neg);          // 34
    fOutputList->Add(trig_05_phi_pt_pos);          // 35
    fOutputList->Add(trig_05_phi_pt_neg);          // 36
    fOutputList->Add(trig_08_phi_pt_pos);          // 37
    fOutputList->Add(trig_08_phi_pt_neg);          // 38

    fOutputList->Add(tof_phi_eta_pos);             // 39
    fOutputList->Add(tof_phi_eta_neg);             // 40
    fOutputList->Add(tof_phi_eta_pos_deut);        // 41
    fOutputList->Add(tof_phi_eta_neg_deut);        // 42

    fOutputList->Add(deut_pt_compare_pos);         // 43
    fOutputList->Add(deut_pt_compare_neg);         // 44
    fOutputList->Add(tpc_sector_fraction);         // 45
    fOutputList->Add(primary_vertex_z);            // 46
    fOutputList->Add(primary_vertex_z_cut);        // 47

    fOutputList->Add(deut_dphi_pt_pos_pos_03_T);   // 48
    fOutputList->Add(deut_dphi_pt_pos_neg_03_T);   // 49
    fOutputList->Add(deut_dphi_pt_neg_neg_03_T);   // 50
    fOutputList->Add(deut_dphi_pt_pos_pos_05_T);   // 51
    fOutputList->Add(deut_dphi_pt_pos_neg_05_T);   // 52
    fOutputList->Add(deut_dphi_pt_neg_neg_05_T);   // 53
    fOutputList->Add(deut_dphi_pt_pos_pos_08_T);   // 54
    fOutputList->Add(deut_dphi_pt_pos_neg_08_T);   // 55
    fOutputList->Add(deut_dphi_pt_neg_neg_08_T);   // 56

    fOutputList->Add(deut_dphi_pt_pos_pos_03_A);   // 57
    fOutputList->Add(deut_dphi_pt_pos_neg_03_A);   // 58
    fOutputList->Add(deut_dphi_pt_neg_neg_03_A);   // 59
    fOutputList->Add(deut_dphi_pt_pos_pos_05_A);   // 60
    fOutputList->Add(deut_dphi_pt_pos_neg_05_A);   // 61
    fOutputList->Add(deut_dphi_pt_neg_neg_05_A);   // 62
    fOutputList->Add(deut_dphi_pt_pos_pos_08_A);   // 63
    fOutputList->Add(deut_dphi_pt_pos_neg_08_A);   // 64
    fOutputList->Add(deut_dphi_pt_neg_neg_08_A);   // 65

    fOutputList->Add(deut_dphi_pt_pos_pos_03_B);   // 66
    fOutputList->Add(deut_dphi_pt_pos_neg_03_B);   // 67
    fOutputList->Add(deut_dphi_pt_neg_neg_03_B);   // 68
    fOutputList->Add(deut_dphi_pt_pos_pos_05_B);   // 69
    fOutputList->Add(deut_dphi_pt_pos_neg_05_B);   // 70
    fOutputList->Add(deut_dphi_pt_neg_neg_05_B);   // 71
    fOutputList->Add(deut_dphi_pt_pos_pos_08_B);   // 72
    fOutputList->Add(deut_dphi_pt_pos_neg_08_B);   // 73
    fOutputList->Add(deut_dphi_pt_neg_neg_08_B);   // 74

    fOutputList->Add(DCAxy_pos);                   // 75
    fOutputList->Add(DCAxy_neg);                   // 76
    fOutputList->Add(DCAz_pos);                    // 77
    fOutputList->Add(DCAz_neg);                    // 78

    fOutputList->Add(m2_pt_pos_fine);              // 79
    fOutputList->Add(m2_pt_neg_fine);              // 80
    fOutputList->Add(m2_pt_pos_cut_fine);          // 81
    fOutputList->Add(m2_pt_neg_cut_fine);          // 82
    fOutputList->Add(m2_pt_pos_cut_T_fine);        // 83
    fOutputList->Add(m2_pt_neg_cut_T_fine);        // 84
    
//  fOutputList->Add(trigger_selection);           // 40

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

//    fPIDResponse->SetRecoPass(1);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFdeut::UserExec(Option_t *)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

//  int x=5;
//  std::cout<<std::bitset<32>(x)<<std::endl;
    
//
//    ULong64_t mask    = AliVEvent::kEMCEJE;
//    ULong64_t anded   = trigger & mask;
//
//    else                      cout<<"false  ";
//    cout<<endl;
//    cout<<"Trigger mask is: "<<trigger<<"   "<<mask<<"   "<<anded<<endl;

//    cout<<std::bitset<64>(trigger)<<"  &  "<<std::bitset<64>(mask)<<"  =  "<<std::bitset<64>(anded)<<endl;
    
//    task->SetTriggerMask(AliVEvent::kEMCEGA);
//    cout<<"Trigger is: "<<fAOD->GetFiredTriggerClasses()<<endl<<endl;
//    cout<<"Trigger is: "<<fAOD->GetEventType()<<endl<<endl;

//    if(TString(fAOD->GetFiredTriggerClasses()).Contains("CINT7"))//   return;
//    {
//	ULong64_t trigger = fAOD->GetTriggerMask();
//	ULong64_t mask    = AliVEvent::kMB;
//	ULong64_t anded   = trigger & mask;
//	cout<<std::bitset<64>(trigger)<<"  &  "<<std::bitset<64>(mask)<<"  =  "<<std::bitset<64>(anded)<<endl;

//	if(trigger & mask ? 1:0)           ;
//	else                         return;

//	cout<<std::bitset<64>(trigger)<<endl;
	    
//	ULong64_t iTrig = 1;
//	for(int i=0; i<64; i++)
//	{
//	    int i = iTrig;
//	    long j = iTrig;
//	    cout<<iTrig<<" ";
//	    if(i==15  ||  i==31){	                                                            }
//	    else	        {
//	    if(trigger & iTrig ? 1:0)   trigger_selection->Fill(i);//	    }
//	    ULong64_t anded   = trigger & iTrig;
//	    cout<<std::bitset<64>(trigger)<<" & "<<std::bitset<64>(iTrig)<<" = "<<std::bitset<64>(anded)<<endl;
//	    iTrig = iTrig * 2;
//	}
//	cout<<endl;
//	cout<<endl<<endl<<endl;
//    }
    
    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7)   ;
    else return;                                                                // count 7208

    const AliAODVertex *primVertex = fAOD->GetPrimaryVertex();
    Double_t pv = primVertex->GetZ();
    primary_vertex_z->Fill(pv);
    
    if(!fAnalysisUtils->IsVertexSelected2013pA(fAOD)) return;
    
//    if(TMath::Abs(pv) > 10.) return;  // handled in IsVertexSelected2013pA
    primary_vertex_z_cut->Fill(pv);                                             // count 6983

    if(fAnalysisUtils->IsPileUpSPD(fAOD)) return;                               // count  !-187
    
    Int_t iTracks(fAOD->GetNumberOfTracks());
    cent_ntracks->Fill(50.0, iTracks);


    
    
/*        
    if(!TString(fAOD->GetFiredTriggerClasses()).Contains("DG1")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("DJ1")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("EG1")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("EJ1")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("DG2")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("DJ2")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("EG2")  &&
       !TString(fAOD->GetFiredTriggerClasses()).Contains("EJ2")
      )                                                                         return;
*/

    
    
//    cout<<fAOD->GetFiredTriggerClasses()<<endl<<endl;



    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

//    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0A"), iTracks);

    ///////////////////////////////////////////////////////////////////////////////////////////////////



    int deut_track_num_T[100];    
    int deut_count_T            = 0;

    int deut_track_num_A[100];    
    int deut_count_A            = 0;
    
    int deut_track_num_B[100];    
    int deut_count_B            = 0;
    
    int trig_03_track_num[100];
    int trig_03_track_count     = 0;
    
    int trig_05_track_num[100];
    int trig_05_track_count     = 0;

    int trig_08_track_num[100];
    int trig_08_track_count     = 0;


    // loop over all these tracks
    //
    // pull out track numbers for high-pt triggers and also deutons

    Float_t max_03_pt = 0.0;    Int_t   max_03    = 0;
    Float_t max_05_pt = 0.0;    Int_t   max_05    = 0;
    Float_t max_08_pt = 0.0;    Int_t   max_08    = 0;

    for(Int_t i(0); i < iTracks; i++)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track)                                                                      {    continue;    }

	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }
//	if(!track->IsGlobalConstrained())                                               {    continue;    }
//	if(!track->IsTPCConstrained())                                                  {    continue;    }
//	if(!track->IsPrimaryCandidate())                                                {    continue;    }
	
	Float_t pt            = track->Pt();
//	if(pt < 0.2)                                                                    {    continue;    }
	if(pt < 0.95)                                                                   {    continue;    }  // makes for fast processings, low momentum play no part of analysis anyway

	Float_t dedx   = track->GetTPCsignal();
	if(dedx > 1000)                                                                 {    continue;    }

	Float_t eta = track->Eta();	if(TMath::Abs(eta) > 0.9)                       {    continue;    }




	
	
	Float_t phi           = track->Phi();
	if(phi <  -pio2){    phi = phi + twopi; }	if(phi <  -pio2){    phi = phi + twopi;   }
	if(phi > 3*pio2){    phi = phi - twopi;	}       if(phi > 3*pio2){    phi = phi - twopi;   }
	

	Float_t fraction = 0.0;
	fraction = track->GetTPCFoundFraction();
	tpc_sector_fraction->Fill(fraction);
	
	fHistPt->Fill(pt);

	// markA
	
	if(pt >= 3.0)	{   trig_03_track_num[trig_03_track_count] = i;	  if(pt>max_03_pt){ max_03_pt = pt; max_03 = trig_03_track_count;   }   trig_03_track_count++;   }
	if(pt >= 5.0)	{   trig_05_track_num[trig_05_track_count] = i;	  if(pt>max_05_pt){ max_05_pt = pt; max_05 = trig_05_track_count;   }   trig_05_track_count++;   }
	if(pt >= 8.0)	{   trig_08_track_num[trig_08_track_count] = i;	  if(pt>max_08_pt){ max_08_pt = pt; max_08 = trig_08_track_count;   }   trig_08_track_count++;   }

	Short_t charge        = track->Charge();
	if(pt >= 3.0)	{   if(charge > 0)   trig_03_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_03_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 5.0)	{   if(charge > 0)   trig_05_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_05_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 8.0)	{   if(charge > 0)   trig_08_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_08_phi_pt_neg->Fill(pt, phi);	}


//	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }
	

	
	Double_t nsigmaTPC = 999.0;	Double_t nsigmaTOF = 999.0;
	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);     /* && trk->IsOn(AliESDtrack::kTPCpid)*/;
	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
	if(!tpcIsOk)	                                                                {    continue;    }
	if(!tofIsOk)	                                                                {    continue;    }




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	// Kalweit and Sharma cuts
	AliExternalTrackParam trkparam;
	trkparam.CopyFromVTrack(track);
	if(trkparam.GetX() > 3.)                                                        {    continue;    }   // only valid for propagation inside the beam pipe

//	dca->Fill(dcaxy, dcaz);


	Short_t nTPCclusters     = track->GetTPCsignalN();
	if(nTPCclusters<70)                                                             {    continue;    }     // cut a

	
	Float_t fraction = 0.0;
	fraction = track->GetTPCFoundFraction();
	tpc_sector_fraction->Fill(fraction);
	if(fraction < 0.8)                                                              {    continue;    }     // cut b
//	if(this->GetNTPCXRowsOverFindable(trk)< 0.8 ) return kFALSE;

	
	fHistPt->Fill(pt);
	
	UChar_t map = track->GetITSClusterMap();
	Int_t nITSpoints = 0;
	for(Int_t j=2; j<6; j++)
	{
		if(map&(1<<j)) ++nITSpoints;
	}
	if(nITSpoints < 2)                                                              {    continue;    }     // 
*/



/*
//	if(fAOD->IsKinkDaughter(track))                                                 {    continue;    }     // cut d
	AliAODVertex* vtx = track->GetProdVertex();
	if(vtx == 0)                                                                    {    continue;    }     // cut d
	if(vtx->GetType() == AliAODVertex::kKink)                                       {    continue;    }     // cut d
	


	if(!track->IsOn(AliAODTrack::kTPCrefit))                                        {    continue;    }     // cut e
	if(!track->IsOn(AliAODTrack::kITSrefit))                                        {    continue;    }     // cut f


//	if(track->GetITSchi2PerCluster()>36)                                            {    continue;    }     // cut k
	
//	Double_t b[2];
//	Double_t bCov[3];
//	if(!trkparam.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 10000., b, bCov))    continue;
//	Double_t dcaxy = b[0];
//	Double_t dcaz  = b[1];
//	if(fabs(dcaxy) > 3.0)                                                           {    continue;    }
//	if(fabs(dcaz)  > 2.0)                                                           {    continue;    }     // cut h


	Double_t dca[2];          // 0 - dcar, 1 - dcaz -> output parameter
	Double_t cov[3];          // covariance matrix -> output parameter
	Double_t kmaxd = 10000.0;    // something big
	AliExternalTrackParam param;
	param.CopyFromVTrack(track);
	if(!param.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), kmaxd, dca, cov)) continue; // don’t accept tracks for which the DCA fit failed.
	if(fabs(dca[0]) > 3.0)                                                          {    continue;    }
	if(fabs(dca[1]) > 2.0)                                                          {    continue;    }	    
*/

//	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }

	Double_t dca[2];          // 0 - dcar, 1 - dcaz -> output parameter
	Double_t cov[3];          // covariance matrix -> output parameter
	Double_t kmaxd = 10000.0;    // something big
	AliExternalTrackParam param;
	param.CopyFromVTrack(track);
	if(!param.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), kmaxd, dca, cov)) continue; // don’t accept tracks for which the DCA fit failed.
//	if(fabs(dca[0]) > 3.0)                                                          {    continue;    }
//	if(fabs(dca[1]) > 2.0)                                                          {    continue;    }
	
	int is_global = 0;
	if(!track->IsGlobalConstrained())                                               {    is_global=1; }
	
	Float_t deltat        = tof_minus_tpion(track);
	float   mom           = track->P();

	Float_t deut_mean     = 0.0;
	Float_t deut_sigma    = 0.0;
	Float_t m2tof         = get_mass_squared(track);
	Float_t beta          = 0.0;
	beta                  = Beta(track);

	    
	if(charge > 0)
	{
	    m2_pt_pos      ->Fill(pt,  m2tof);
	    m2_pt_pos_fine ->Fill(pt,  m2tof);
	    beta_p_pos     ->Fill(mom, Beta(track));
	    deltat_pt_pos  ->Fill(pt,  deltat);
	    if(pt >= 1.0  &&  pt < 2.0)   tof_phi_eta_pos->Fill(phi, eta);

	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < 3.0)
	    {
		m2_pt_pos_cut     ->Fill(pt,  m2tof);
		m2_pt_pos_cut_fine->Fill(pt,  m2tof);
		beta_p_pos_cut    ->Fill(mom, Beta(track));
		deltat_pt_pos_cut ->Fill(pt,  deltat);
			
		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
		    deut_sigma = fit_deut_curve->Eval(pt);

		    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_T[deut_count_T] = i;
			deut_count_T++;
			m2_pt_pos_cut_T->Fill(pt,m2tof);
			m2_pt_pos_cut_T_fine->Fill(pt,m2tof);
			DCAxy_pos->Fill(dca[0]);
			DCAz_pos->Fill(dca[1]);
			if(pt >= 1.00  &&  pt < 2.00){   tof_phi_eta_pos_deut->Fill(phi, eta);                     }
			if(pt >= 1.00  &&  pt < 1.25){   deut_pt_compare_pos->Fill(pt - get_deut_tof_pt(track));   }
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_T->Fill(pt, deut_phi);
		    }

		    if(m2tof < deut_mean + 3.0 * deut_sigma  &&   m2tof > deut_mean - 3.0 * deut_sigma)
		    {
			if(is_global == 1)           {   m2_pt_pos_cut_G->Fill(pt,m2tof);                          }
		    }
		    
		    if(m2tof +1.2< deut_mean + cut_width * deut_sigma  &&   m2tof +1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_A[deut_count_A] = i;
			deut_count_A++;
			m2_pt_pos_cut_A->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_A->Fill(pt, deut_phi);
		    }
		    if(m2tof -1.2< deut_mean + cut_width * deut_sigma  &&   m2tof -1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_B[deut_count_A] = i;
			deut_count_B++;
			m2_pt_pos_cut_B->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_B->Fill(pt, deut_phi);
		    }		    
		}
	    }
	}    //   end of pos charge if statement
	    
	else if(charge < 0)
	{
	    m2_pt_neg      ->Fill(pt,  m2tof);
	    m2_pt_neg_fine ->Fill(pt,  m2tof);
	    beta_p_neg     ->Fill(mom, Beta(track));
	    deltat_pt_neg  ->Fill(pt,  deltat);
	    if(pt >= 1.0  &&  pt < 2.0)   tof_phi_eta_neg->Fill(phi, eta);

	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < 3.0)
	    {
		m2_pt_neg_cut     ->Fill(pt, m2tof);
		m2_pt_neg_cut_fine->Fill(pt, m2tof);
		beta_p_neg_cut    ->Fill(mom, Beta(track));
		deltat_pt_neg_cut ->Fill(pt, deltat);

		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][1][w]);   }
		    deut_sigma = fit_deut_curve->Eval(pt);

		    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_T[deut_count_T] = i;
			deut_count_T++;
			m2_pt_neg_cut_T->Fill(pt,m2tof);
			m2_pt_neg_cut_T_fine->Fill(pt,m2tof);
			DCAxy_neg->Fill(dca[0]);
			DCAz_neg->Fill(dca[1]);
			if(pt >= 1.00  &&  pt < 2.00){   tof_phi_eta_neg_deut->Fill(phi, eta);                     }
			if(pt >= 1.00  &&  pt < 1.25){   deut_pt_compare_neg->Fill(pt - get_deut_tof_pt(track));   }	
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_T->Fill(pt, deut_phi);
		    }

		    
		    if(m2tof < deut_mean + 3.0 * deut_sigma  &&   m2tof > deut_mean - 3.0 * deut_sigma)
		    {
			if(is_global == 1)           {   m2_pt_neg_cut_G->Fill(pt,m2tof);                          }
		    }
		    
		    if(m2tof +1.2< deut_mean + cut_width * deut_sigma  &&   m2tof +1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_A[deut_count_A] = i;
			deut_count_A++;
			m2_pt_neg_cut_A->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_A->Fill(pt, deut_phi);
		    }
		    if(m2tof -1.2< deut_mean + cut_width * deut_sigma  &&   m2tof -1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_B[deut_count_B] = i;
			deut_count_B++;
			m2_pt_neg_cut_B->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_B->Fill(pt, deut_phi);
		    }
		}
	    }
	}    //   end of neg charge if statement
    }        //   end of track loop


    deut_per_event->Fill(deut_count_T);
    if(trig_03_track_count > 0)    {	trig_03_per_event->Fill(trig_03_track_count);    }
    if(trig_05_track_count > 0)    {	trig_05_per_event->Fill(trig_05_track_count);    }
    if(trig_08_track_count > 0)    {	trig_08_per_event->Fill(trig_08_track_count);    }



    
    // markB
    int begin = 0;
    if(deut_count_T > 0  &&  trig_03_track_count > 0  &&  do_lead_only == 1)    {	begin = max_03;   trig_03_track_count = max_03 + 1;    }
    if(deut_count_T > 0  &&  trig_03_track_count > 0)
    {
	for(int i=begin; i<trig_03_track_count; i++)  // trigger loop
	{
	    int H               = trig_03_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    
	    for(int j=0; j<deut_count_T; j++)
	    {
		int A               = deut_track_num_T[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_03_T->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_03_T->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_03_T->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_03_T->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }

    

    begin = 0;
    if(deut_count_T > 0  &&  trig_05_track_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_track_count = max_05 + 1;    }
    if(deut_count_T > 0  &&  trig_05_track_count > 0)
    {
	for(int i=begin; i<trig_05_track_count; i++)  // trigger loop
	{
	    int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_T; j++)
	    {
		int A               = deut_track_num_T[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_T->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_T->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }


    begin = 0;
    if(deut_count_T > 0  &&  trig_08_track_count > 0  &&  do_lead_only == 1)    {   begin = max_08;   trig_08_track_count = max_08 + 1;    }
    if(deut_count_T > 0  &&  trig_08_track_count > 0)
    {
	for(int i=begin; i<trig_08_track_count; i++)  // trigger loop
	{
	    int H               = trig_08_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_T; j++)
	    {
		int A               = deut_track_num_T[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_08_T->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_08_T->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_08_T->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_08_T->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_A > 0  &&  trig_03_track_count > 0  &&  do_lead_only == 1)    {	begin = max_03;   trig_03_track_count = max_03 + 1;    }
    if(deut_count_A > 0  &&  trig_03_track_count > 0)
    {
	for(int i=begin; i<trig_03_track_count; i++)  // trigger loop
	{
	    int H               = trig_03_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_A; j++)
	    {
		int A               = deut_track_num_A[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_03_A->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_03_A->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_03_A->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_03_A->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }

    

    begin = 0;
    if(deut_count_A > 0  &&  trig_05_track_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_track_count = max_05 + 1;    }
    if(deut_count_A > 0  &&  trig_05_track_count > 0)
    {
	for(int i=begin; i<trig_05_track_count; i++)  // trigger loop
	{
	    int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_A; j++)
	    {
		int A               = deut_track_num_A[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_A->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_A->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_A->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_A->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }


    begin = 0;
    if(deut_count_A > 0  &&  trig_08_track_count > 0  &&  do_lead_only == 1)    {   begin = max_08;   trig_08_track_count = max_08 + 1;    }
    if(deut_count_A > 0  &&  trig_08_track_count > 0)
    {
	for(int i=begin; i<trig_08_track_count; i++)  // trigger loop
	{
	    int H               = trig_08_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_A; j++)
	    {
		int A               = deut_track_num_A[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_08_A->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_08_A->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_08_A->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_08_A->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }    

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_B > 0  &&  trig_03_track_count > 0  &&  do_lead_only == 1)    {	begin = max_03;   trig_03_track_count = max_03 + 1;    }
    if(deut_count_B > 0  &&  trig_03_track_count > 0)
    {
	for(int i=begin; i<trig_03_track_count; i++)  // trigger loop
	{
	    int H               = trig_03_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_B; j++)
	    {
		int A               = deut_track_num_B[j];
		if(A != H)
		{
		    if(AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A)))
//		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_03_B->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_03_B->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_03_B->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_03_B->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }

    

    begin = 0;
    if(deut_count_B > 0  &&  trig_05_track_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_track_count = max_05 + 1;    }
    if(deut_count_B > 0  &&  trig_05_track_count > 0)
    {
	for(int i=begin; i<trig_05_track_count; i++)  // trigger loop
	{
	    int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_B; j++)
	    {
		int A               = deut_track_num_B[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_B->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_B->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_B->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_B->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }


    begin = 0;
    if(deut_count_B > 0  &&  trig_08_track_count > 0  &&  do_lead_only == 1)    {   begin = max_08;   trig_08_track_count = max_08 + 1;    }
    if(deut_count_B > 0  &&  trig_08_track_count > 0)
    {
	for(int i=begin; i<trig_08_track_count; i++)  // trigger loop
	{
	    int H               = trig_08_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	    for(int j=0; j<deut_count_B; j++)
	    {
		int A               = deut_track_num_B[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    Float_t phi_A       = trackA->Phi();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_08_B->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_08_B->Fill(pt_A, Sdphi);	}    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_08_B->Fill(pt_A, Sdphi);	}
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_08_B->Fill(pt_A, Sdphi);	}    }
		    }
    }   }   }   }
    
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
    Double_t m2             = 0.0;
    Double_t c2             = 29.9792458;
    
    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);
    return m2;

}



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFdeut::get_deut_tof_pt(AliAODTrack *track)
{
    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
    Double_t stop_time      = track->GetTOFsignal();
    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
    Double_t tof            = stop_time - start_time;
    tof                     = tof * 0.001;                                                                  // convert from ps to ns
    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
    length                  = length * 100;                                                                 // convert to cm
//  Double_t mom            = track->P();
//  Double_t m2             = 0.0;
    Double_t c2             = 29.9792458;

    Double_t deut_mass      = 1.8756;                                                                       // GeV/(c^2)
    
//    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);

    Double_t mom            = 0.00;

    mom                     = deut_mass / sqrt(
	                                          (tof*tof * c2 * c2 / (length * length) - 1)
                                              );

    Double_t pt             = 0.00;
    pt                      = sqrt(pow(mom,2) - pow(track->Pz(),2));
    return pt;

}




//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFdeut::get_deut_tof_p(AliAODTrack *track)
{
    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
    Double_t stop_time      = track->GetTOFsignal();
    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
    Double_t tof            = stop_time - start_time;
    tof                     = tof * 0.001;                                                                  // convert from ps to ns
    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
    length                  = length * 100;                                                                 // convert to cm
//  Double_t mom            = track->P();
//  Double_t m2             = 0.0;
    Double_t c2             = 29.9792458;

    Double_t deut_mass      = 1.8756;                                                                       // GeV/(c^2)
    
//    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);

    Double_t mom            = 0.00;

    mom                     = deut_mass / sqrt(
	                                          (tof*tof * c2 * c2 / (length * length) - 1)
                                              );
    return mom;

}
