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
#include "AliAnalysisTaskCorPIDTOFQA.h"
#include "AliAODHandler.h"
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


ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),


    fHistPt(0),                    //  1
    cent_ntracks(0),               //  2
    
//  m2_pt_pos(0),                  //  3
//  m2_pt_neg(0),                  //  4
    m2_pt_pos_TPC(0),              //  5
    m2_pt_neg_TPC(0),              //  6
    m2_pt_pos_cut_T(0),            //  7
    m2_pt_neg_cut_T(0),            //  8
    m2_pt_pos_cut_G(0),            //  9
    m2_pt_neg_cut_G(0),            // 10
    m2_pt_pos_cut_A(0),            // 11
    m2_pt_neg_cut_A(0),            // 12
    m2_pt_pos_cut_B(0),            // 13
    m2_pt_neg_cut_B(0),            // 14

    m2_pt_pos_cut_T_prot(0),       //  7a
    m2_pt_neg_cut_T_prot(0),       //  8a
	
//  m2_pt_pos_cut_with_trig_04(0), // 15
//  m2_pt_neg_cut_with_trig_04(0), // 16
    m2_pt_pos_cut_with_trig_05(0), // 17
    m2_pt_neg_cut_with_trig_05(0), // 18
//  m2_pt_pos_cut_with_trig_06(0), // 19
//  m2_pt_neg_cut_with_trig_06(0), // 20

    deut_phi_pt_pos_T(0),          // 21
    deut_phi_pt_neg_T(0),          // 22

    prot_phi_pt_pos_T(0),          // 21a
    prot_phi_pt_neg_T(0),          // 22a
    
    deut_phi_pt_pos_A(0),          // 23
    deut_phi_pt_neg_A(0),          // 24
    deut_phi_pt_pos_B(0),          // 25
    deut_phi_pt_neg_B(0),          // 26

    deut_per_event(0),             // 27
    prot_per_event(0),             // 27a
//  trig_04_per_event(0),          // 28
    trig_05_per_event(0),          // 29
//  trig_06_per_event(0),          // 30
    
//  trig_04_phi_pt_pos(0),         // 31
//  trig_04_phi_pt_neg(0),         // 32    
    trig_05_phi_pt_pos(0),         // 33
    trig_05_phi_pt_neg(0),         // 34
//  trig_06_phi_pt_pos(0),         // 35
//  trig_06_phi_pt_neg(0),         // 36
    
    tof_phi_eta_pos(0),            // 37
    tof_phi_eta_neg(0),            // 38
    tof_phi_eta_pos_deut(0),       // 39
    tof_phi_eta_neg_deut(0),       // 40

    deut_pt_compare_pos(0),        // 41
    deut_pt_compare_neg(0),        // 42
    tpc_sector_fraction(0),        // 43
    primary_vertex_z(0),           // 44
    primary_vertex_z_cut1(0),      // 45
    primary_vertex_z_cut2(0),      // 46
	
//  deut_dphi_pt_pos_pos_04_T(0),  // 47
//  deut_dphi_pt_pos_neg_04_T(0),  // 48
//  deut_dphi_pt_neg_neg_04_T(0),  // 49
    deut_dphi_pt_pos_pos_05_T(0),  // 50
    deut_dphi_pt_pos_neg_05_T(0),  // 51
    deut_dphi_pt_neg_neg_05_T(0),  // 52

//  deut_dphi_pt_pos_pos_06_T(0),  // 53
//  deut_dphi_pt_pos_neg_06_T(0),  // 54
//  deut_dphi_pt_neg_neg_06_T(0),  // 55
    
    prot_dphi_pt_pos_pos_05_T(0),  // 50a
    prot_dphi_pt_pos_neg_05_T(0),  // 51a
    prot_dphi_pt_neg_neg_05_T(0),  // 52a
    
//  deut_dphi_pt_pos_pos_04_A(0),  // 56
//  deut_dphi_pt_pos_neg_04_A(0),  // 57
//  deut_dphi_pt_neg_neg_04_A(0),  // 58
    deut_dphi_pt_pos_pos_05_A(0),  // 59
    deut_dphi_pt_pos_neg_05_A(0),  // 60
    deut_dphi_pt_neg_neg_05_A(0),  // 61
//  deut_dphi_pt_pos_pos_06_A(0),  // 62
//  deut_dphi_pt_pos_neg_06_A(0),  // 63
//  deut_dphi_pt_neg_neg_06_A(0),  // 64
 
//  deut_dphi_pt_pos_pos_04_B(0),  // 65
//  deut_dphi_pt_pos_neg_04_B(0),  // 66
//  deut_dphi_pt_neg_neg_04_B(0),  // 67
    deut_dphi_pt_pos_pos_05_B(0),  // 68
    deut_dphi_pt_pos_neg_05_B(0),  // 69
    deut_dphi_pt_neg_neg_05_B(0),  // 70
//  deut_dphi_pt_pos_pos_06_B(0),  // 71
//  deut_dphi_pt_pos_neg_06_B(0),  // 72
//  deut_dphi_pt_neg_neg_06_B(0),  // 73

    DCAxy_pos(0),                  // 74
    DCAxy_neg(0),                  // 75
    DCAz_pos(0),                   // 76
    DCAz_neg(0),                   // 77
    
    m2_pt_pos_fine(0),             // 78
    m2_pt_neg_fine(0),             // 79
    m2_pt_pos_TPC_fine(0),         // 80
    m2_pt_neg_TPC_fine(0),         // 81 
    m2_pt_pos_cut_T_fine(0),       // 82
    m2_pt_neg_cut_T_fine(0),       // 83
    m2_pt_pos_cut_A_fine(0),       // 84
    m2_pt_neg_cut_A_fine(0),       // 85
    m2_pt_pos_cut_B_fine(0),       // 86
    m2_pt_neg_cut_B_fine(0),       // 87

    deut_m2_dphi_04_pt1(0),        // 88
    deut_m2_dphi_04_pt2(0),        // 89
    deut_m2_dphi_04_pt3(0),        // 90
    deut_m2_dphi_04_pt4(0),        // 91
    
    m2_pt_pos_TPC_prot_fine(0),    // 92
    m2_pt_neg_TPC_prot_fine(0),    // 93

    m2_pt_pos_cut_T_prot_fine(0),  // 94
    m2_pt_neg_cut_T_prot_fine(0)   // 95
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),

    fHistPt(0),                    //  1
    cent_ntracks(0),               //  2
    
//  m2_pt_pos(0),                  //  3
//  m2_pt_neg(0),                  //  4
    m2_pt_pos_TPC(0),              //  5
    m2_pt_neg_TPC(0),              //  6
    m2_pt_pos_cut_T(0),            //  7
    m2_pt_neg_cut_T(0),            //  8
    m2_pt_pos_cut_G(0),            //  9
    m2_pt_neg_cut_G(0),            // 10
    m2_pt_pos_cut_A(0),            // 11
    m2_pt_neg_cut_A(0),            // 12
    m2_pt_pos_cut_B(0),            // 13
    m2_pt_neg_cut_B(0),            // 14

    m2_pt_pos_cut_T_prot(0),       //  7a
    m2_pt_neg_cut_T_prot(0),       //  8a
	
//  m2_pt_pos_cut_with_trig_04(0), // 15
//  m2_pt_neg_cut_with_trig_04(0), // 16
    m2_pt_pos_cut_with_trig_05(0), // 17
    m2_pt_neg_cut_with_trig_05(0), // 18
//  m2_pt_pos_cut_with_trig_06(0), // 19
//  m2_pt_neg_cut_with_trig_06(0), // 20

    deut_phi_pt_pos_T(0),          // 21
    deut_phi_pt_neg_T(0),          // 22

    prot_phi_pt_pos_T(0),          // 21a
    prot_phi_pt_neg_T(0),          // 22a
    
    deut_phi_pt_pos_A(0),          // 23
    deut_phi_pt_neg_A(0),          // 24
    deut_phi_pt_pos_B(0),          // 25
    deut_phi_pt_neg_B(0),          // 26

    deut_per_event(0),             // 27
    prot_per_event(0),             // 27a
//  trig_04_per_event(0),          // 28
    trig_05_per_event(0),          // 29
//  trig_06_per_event(0),          // 30
    
//  trig_04_phi_pt_pos(0),         // 31
//  trig_04_phi_pt_neg(0),         // 32    
    trig_05_phi_pt_pos(0),         // 33
    trig_05_phi_pt_neg(0),         // 34
//  trig_06_phi_pt_pos(0),         // 35
//  trig_06_phi_pt_neg(0),         // 36
    
    tof_phi_eta_pos(0),            // 37
    tof_phi_eta_neg(0),            // 38
    tof_phi_eta_pos_deut(0),       // 39
    tof_phi_eta_neg_deut(0),       // 40

    deut_pt_compare_pos(0),        // 41
    deut_pt_compare_neg(0),        // 42
    tpc_sector_fraction(0),        // 43
    primary_vertex_z(0),           // 44
    primary_vertex_z_cut1(0),      // 45
    primary_vertex_z_cut2(0),      // 46
	
//  deut_dphi_pt_pos_pos_04_T(0),  // 47
//  deut_dphi_pt_pos_neg_04_T(0),  // 48
//  deut_dphi_pt_neg_neg_04_T(0),  // 49
    deut_dphi_pt_pos_pos_05_T(0),  // 50
    deut_dphi_pt_pos_neg_05_T(0),  // 51
    deut_dphi_pt_neg_neg_05_T(0),  // 52

//  deut_dphi_pt_pos_pos_06_T(0),  // 53
//  deut_dphi_pt_pos_neg_06_T(0),  // 54
//  deut_dphi_pt_neg_neg_06_T(0),  // 55
    
    prot_dphi_pt_pos_pos_05_T(0),  // 50a
    prot_dphi_pt_pos_neg_05_T(0),  // 51a
    prot_dphi_pt_neg_neg_05_T(0),  // 52a
    
//  deut_dphi_pt_pos_pos_04_A(0),  // 56
//  deut_dphi_pt_pos_neg_04_A(0),  // 57
//  deut_dphi_pt_neg_neg_04_A(0),  // 58
    deut_dphi_pt_pos_pos_05_A(0),  // 59
    deut_dphi_pt_pos_neg_05_A(0),  // 60
    deut_dphi_pt_neg_neg_05_A(0),  // 61
//  deut_dphi_pt_pos_pos_06_A(0),  // 62
//  deut_dphi_pt_pos_neg_06_A(0),  // 63
//  deut_dphi_pt_neg_neg_06_A(0),  // 64
 
//  deut_dphi_pt_pos_pos_04_B(0),  // 65
//  deut_dphi_pt_pos_neg_04_B(0),  // 66
//  deut_dphi_pt_neg_neg_04_B(0),  // 67
    deut_dphi_pt_pos_pos_05_B(0),  // 68
    deut_dphi_pt_pos_neg_05_B(0),  // 69
    deut_dphi_pt_neg_neg_05_B(0),  // 70
//  deut_dphi_pt_pos_pos_06_B(0),  // 71
//  deut_dphi_pt_pos_neg_06_B(0),  // 72
//  deut_dphi_pt_neg_neg_06_B(0),  // 73

    DCAxy_pos(0),                  // 74
    DCAxy_neg(0),                  // 75
    DCAz_pos(0),                   // 76
    DCAz_neg(0),                   // 77
    
    m2_pt_pos_fine(0),             // 78
    m2_pt_neg_fine(0),             // 79
    m2_pt_pos_TPC_fine(0),         // 80
    m2_pt_neg_TPC_fine(0),         // 81 
    m2_pt_pos_cut_T_fine(0),       // 82
    m2_pt_neg_cut_T_fine(0),       // 83
    m2_pt_pos_cut_A_fine(0),       // 84
    m2_pt_neg_cut_A_fine(0),       // 85
    m2_pt_pos_cut_B_fine(0),       // 86
    m2_pt_neg_cut_B_fine(0),       // 87

    deut_m2_dphi_04_pt1(0),        // 88
    deut_m2_dphi_04_pt2(0),        // 89
    deut_m2_dphi_04_pt3(0),        // 90
    deut_m2_dphi_04_pt4(0),        // 91
    
    m2_pt_pos_TPC_prot_fine(0),    // 92
    m2_pt_neg_TPC_prot_fine(0),    // 93

    m2_pt_pos_cut_T_prot_fine(0),  // 94
    m2_pt_neg_cut_T_prot_fine(0)   // 95
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

//    cout<<endl<<endl<<endl<<name<<endl<<endl<<endl;

//    int i = 0;
    run_mode = atoi(name);
    cout<<endl<<endl<<endl<<run_mode<<endl<<endl<<endl;

    if(run_mode == 1)   cut_width = 3;
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::~AliAnalysisTaskCorPIDTOFQA()
{
    // destructor
    if(fAnalysisUtils) delete fAnalysisUtils;
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserCreateOutputObjects()
{
    // 2.88465 0.0761582 0.709281 0.124386 0.017642 -0.0316078 2.65738 0.115151 0.918566 0.0986592 0.0187545 0.00346519    // pp 2016 untriggered

    //3.38661   0.0499084  0.200357   0.00351171 -0.0019421  0.101265
    //0.668755 -0.313122  -0.639249   0.0563519   0          0.359363
    //3.45988   0.0545625  0.0827224 -0.0253765   0.00486567 0.15201
    //0.636145 -0.288126  -0.632198   0.0523836   0          0.368737

    deut_curves[0][0][0] = 3.38661;     // pos deut mean curve
    deut_curves[0][0][1] = 0.0499084;
    deut_curves[0][0][2] = 0.200357;
    deut_curves[0][0][3] = 0.00351171;
    deut_curves[0][0][4] =-0.0019421;
    deut_curves[0][0][5] = 0.101265;    

    deut_curves[0][1][0] = 0.668755;    // pos deut sigma curve
    deut_curves[0][1][1] =-0.313122;
    deut_curves[0][1][2] =-0.639249;
    deut_curves[0][1][3] = 0.0563519;
    deut_curves[0][1][4] = 0.0;
    deut_curves[0][1][5] = 0.359363;

    deut_curves[1][0][0] = 3.45988;     // neg deut mean curve
    deut_curves[1][0][1] = 0.0545625;
    deut_curves[1][0][2] = 0.0827224;
    deut_curves[1][0][3] =-0.0253765;
    deut_curves[1][0][4] = 0.00486567;
    deut_curves[1][0][5] = 0.15201;

    deut_curves[1][1][0] = 0.636145;   // neg deut sigma curve
    deut_curves[1][1][1] =-0.288126;
    deut_curves[1][1][2] =-0.632198;
    deut_curves[1][1][3] = 0.0523836;
    deut_curves[1][1][4] = 0.0;
    deut_curves[1][1][5] = 0.368737;


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

    
    /*
    deut_curves[0][0][0] = 2.88465;     // pos deut mean curve
    deut_curves[0][0][1] = 0.0761582;
    deut_curves[0][0][2] = 0.709281;

    deut_curves[0][1][0] = 0.124386;    // pos deut sigma curve
    deut_curves[0][1][1] = 0.017642;
    deut_curves[0][1][2] = -0.0316078;

    deut_curves[1][0][0] = 2.65738;     // neg deut mean curve
    deut_curves[1][0][1] = 0.115151;
    deut_curves[1][0][2] = 0.918566;

    deut_curves[1][1][0] = 0.0986592;   // neg deut sigma curve
    deut_curves[1][1][1] = 0.0187545;
    deut_curves[1][1][2] = 0.00346519;
    */

    
//2.6836 0.09 0.92872 -0.445429 0.121074 0.452984 2.69467 0.09 0.92872 -0.121833 0.06763 0.179748  // p+Pb 2013 untriggered

/*
    deut_curves[0][0][0] =  2.6836;     // pos deut mean curve
    deut_curves[0][0][1] =  0.09;
    deut_curves[0][0][2] =  0.92872;

    deut_curves[0][1][0] = -0.445429;   // pos deut sigma curve
    deut_curves[0][1][1] =  0.121074;
    deut_curves[0][1][2] =  0.452984;

    deut_curves[1][0][0] =  2.69467;    // neg deut mean curve
    deut_curves[1][0][1] =  0.09;
    deut_curves[1][0][2] =  0.92872;

    deut_curves[1][1][0] = -0.121833;   // neg deut sigma curve
    deut_curves[1][1][1] =  0.06763;
    deut_curves[1][1][2] =  0.179748;
*/  
    
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


    Double_t coarse_binning[19]  = {0.8, 0.9,  1.0, 1.1,  1.2, 1.35,  1.5,  1.65,  1.8, 2.0, 2.2, 2.4, 2.7,  3.0, 3.3,  3.6, 4.0, 4.5,  5.0};
    Double_t coarse_binning2[19] = {0.4, 0.45, 0.5, 0.55, 0.6, 0.675, 0.75, 0.825, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.65, 1.8, 2.0, 2.25, 2.5};
    
    fHistPt                    = new TH1F("fHistPt",                    "Pt()",                         18, coarse_binning);                              //  1
    cent_ntracks               = new TH2F("cent_ntracks",               "cent_ntracks",                100,      0,    100,     100,       0,     800);   //  2
    
//  m2_pt_pos                  = new TH2F("m2_pt_pos",                  "m2_pt_pos",                    18, coarse_binning,    2400,    -1.0,     7.0);   //  3
//  m2_pt_neg                  = new TH2F("m2_pt_neg",                  "m2_pt_neg",                    18, coarse_binning,    2400,    -1.0,     7.0);   //  4
    m2_pt_pos_TPC              = new TH2F("m2_pt_pos_TPC",              "m2_pt_pos_TPC",                18, coarse_binning,    2400,    -1.0,     7.0);   //  5
    m2_pt_neg_TPC              = new TH2F("m2_pt_neg_TPC",              "m2_pt_neg_TPC",                18, coarse_binning,    2400,    -1.0,     7.0);   //  6
    m2_pt_pos_cut_T            = new TH2F("m2_pt_pos_cut_T",            "m2_pt_pos_cut_T",              18, coarse_binning,    2400,    -1.0,     7.0);   //  7
    m2_pt_neg_cut_T            = new TH2F("m2_pt_neg_cut_T",            "m2_pt_neg_cut_T",              18, coarse_binning,    2400,    -1.0,     7.0);   //  8
    m2_pt_pos_cut_G            = new TH2F("m2_pt_pos_cut_G",            "m2_pt_pos_cut_G",              18, coarse_binning,    2400,    -1.0,     7.0);   //  9
    m2_pt_neg_cut_G            = new TH2F("m2_pt_neg_cut_G",            "m2_pt_neg_cut_G",              18, coarse_binning,    2400,    -1.0,     7.0);   // 10
    m2_pt_pos_cut_A            = new TH2F("m2_pt_pos_cut_A",            "m2_pt_pos_cut_A",              18, coarse_binning,    2400,    -1.0,     7.0);   // 11
    m2_pt_neg_cut_A            = new TH2F("m2_pt_neg_cut_A",            "m2_pt_neg_cut_A",              18, coarse_binning,    2400,    -1.0,     7.0);   // 12
    m2_pt_pos_cut_B            = new TH2F("m2_pt_pos_cut_B",            "m2_pt_pos_cut_B",              18, coarse_binning,    2400,    -1.0,     7.0);   // 13
    m2_pt_neg_cut_B            = new TH2F("m2_pt_neg_cut_B",            "m2_pt_neg_cut_B",              18, coarse_binning,    2400,    -1.0,     7.0);   // 14

    m2_pt_pos_cut_T_prot       = new TH2F("m2_pt_pos_cut_T_prot",       "m2_pt_pos_cut_T_prot",         18,coarse_binning2,    2400,    -1.0,     7.0);   //  7a
    m2_pt_neg_cut_T_prot       = new TH2F("m2_pt_neg_cut_T_prot",       "m2_pt_neg_cut_T_prot",         18,coarse_binning2,    2400,    -1.0,     7.0);   //  8a
    
//  m2_pt_pos_cut_with_trig_04 = new TH2F("m2_pt_pos_cut_with_trig_04", "m2_pt_pos_cut_with_trig_04",   18, coarse_binning,    2400,    -1.0,     7.0);   // 15
//  m2_pt_neg_cut_with_trig_04 = new TH2F("m2_pt_neg_cut_with_trig_04", "m2_pt_neg_cut_with_trig_04",   18, coarse_binning,    2400,    -1.0,     7.0);   // 16
    m2_pt_pos_cut_with_trig_05 = new TH2F("m2_pt_pos_cut_with_trig_05", "m2_pt_pos_cut_with_trig_05",   18, coarse_binning,    2400,    -1.0,     7.0);   // 17
    m2_pt_neg_cut_with_trig_05 = new TH2F("m2_pt_neg_cut_with_trig_05", "m2_pt_neg_cut_with_trig_05",   18, coarse_binning,    2400,    -1.0,     7.0);   // 18
//  m2_pt_pos_cut_with_trig_06 = new TH2F("m2_pt_pos_cut_with_trig_06", "m2_pt_pos_cut_with_trig_06",   18, coarse_binning,    2400,    -1.0,     7.0);   // 19
//  m2_pt_neg_cut_with_trig_06 = new TH2F("m2_pt_neg_cut_with_trig_06", "m2_pt_neg_cut_with_trig_06",   18, coarse_binning,    2400,    -1.0,     7.0);   // 20
    
    deut_phi_pt_pos_T          = new TH2F("deut_phi_pt_pos_T",          "deut_phi_pt_pos_T",            18, coarse_binning,     288,   lower,   upper);   // 21
    deut_phi_pt_neg_T          = new TH2F("deut_phi_pt_neg_T",          "deut_phi_pt_neg_T",            18, coarse_binning,     288,   lower,   upper);   // 22
    deut_phi_pt_pos_A          = new TH2F("deut_phi_pt_pos_A",          "deut_phi_pt_pos_A",            18, coarse_binning,     288,   lower,   upper);   // 23
    deut_phi_pt_neg_A          = new TH2F("deut_phi_pt_neg_A",          "deut_phi_pt_neg_A",            18, coarse_binning,     288,   lower,   upper);   // 24
    deut_phi_pt_pos_B          = new TH2F("deut_phi_pt_pos_B",          "deut_phi_pt_pos_B",            18, coarse_binning,     288,   lower,   upper);   // 25
    deut_phi_pt_neg_B          = new TH2F("deut_phi_pt_neg_B",          "deut_phi_pt_neg_B",            18, coarse_binning,     288,   lower,   upper);   // 26
    prot_phi_pt_pos_T          = new TH2F("prot_phi_pt_pos_T",          "prot_phi_pt_pos_T",            18,coarse_binning2,     288,   lower,   upper);   // 21a
    prot_phi_pt_neg_T          = new TH2F("prot_phi_pt_neg_T",          "prot_phi_pt_neg_T",            18,coarse_binning2,     288,   lower,   upper);   // 22a
    
    deut_per_event             = new TH1I("deut_per_event",             "deut_per_event",                4,      0,      4);                              // 27
    prot_per_event             = new TH1I("prot_per_event",             "prot_per_event",               20,      0,     20);                              // 27
//  trig_04_per_event          = new TH1I("trig_04_per_event",          "trig_04_per_event",            10,      0,     10);                              // 28
    trig_05_per_event          = new TH1I("trig_05_per_event",          "trig_05_per_event",            10,      0,     10);                              // 29
//  trig_06_per_event          = new TH1I("trig_06_per_event",          "trig_06_per_event",            10,      0,     10);                              // 30

//  trig_04_phi_pt_pos         = new TH2F("trig_04_phi_pt_pos",         "trig_04_phi_pt_pos",          170,    3.0,   20.0,     288,   lower,   upper);   // 31
//  trig_04_phi_pt_neg         = new TH2F("trig_04_phi_pt_neg",         "trig_04_phi_pt_neg",          170,    3.0,   20.0,     288,   lower,   upper);   // 32
    trig_05_phi_pt_pos         = new TH2F("trig_05_phi_pt_pos",         "trig_05_phi_pt_pos",          170,    3.0,   20.0,     288,   lower,   upper);   // 33
    trig_05_phi_pt_neg         = new TH2F("trig_05_phi_pt_neg",         "trig_05_phi_pt_neg",          170,    3.0,   20.0,     288,   lower,   upper);   // 34
//  trig_06_phi_pt_pos         = new TH2F("trig_06_phi_pt_pos",         "trig_06_phi_pt_pos",          170,    3.0,   20.0,     288,   lower,   upper);   // 35
//  trig_06_phi_pt_neg         = new TH2F("trig_06_phi_pt_neg",         "trig_06_phi_pt_neg",          170,    3.0,   20.0,     288,   lower,   upper);   // 36

    tof_phi_eta_pos            = new TH2F("tof_phi_eta_pos",            "tof_phi_eta_pos",             190,  lower,  upper,      60,   -0.92,    0.92);   // 37
    tof_phi_eta_neg            = new TH2F("tof_phi_eta_neg",            "tof_phi_eta_neg",             190,  lower,  upper,      60,   -0.92,    0.92);   // 38
    tof_phi_eta_pos_deut       = new TH2F("tof_phi_eta_pos_deut",       "tof_phi_eta_pos_deut",        190,  lower,  upper,      60,   -0.92,    0.92);   // 39
    tof_phi_eta_neg_deut       = new TH2F("tof_phi_eta_neg_deut",       "tof_phi_eta_neg_deut",        190,  lower,  upper,      60,   -0.92,    0.92);   // 40

    deut_pt_compare_pos        = new TH1F("deut_pt_compare_pos",        "deut_pt_compare_pos",         100,   -0.2,    0.2);                              // 41
    deut_pt_compare_neg        = new TH1F("deut_pt_compare_neg",        "deut_pt_compare_neg",         100,   -0.2,    0.2);                              // 42
    tpc_sector_fraction        = new TH1F("tpc_sector_fraction",        "tpc_sector_fraction",         100,    0.0,    1.0);                              // 43
    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            400,  -20.0,   20.0);                              // 44
    primary_vertex_z_cut1      = new TH1F("primary_vertex_z_cut1",      "primary_vertex_z_cut1",       400,  -20.0,   20.0);                              // 45
    primary_vertex_z_cut2      = new TH1F("primary_vertex_z_cut2",      "primary_vertex_z_cut2",       400,  -20.0,   20.0);                              // 46

//  deut_dphi_pt_pos_pos_04_T  = new TH2F("deut_dphi_pt_pos_pos_04_T",  "deut_dphi_pt_pos_pos_04_T",    18, coarse_binning,     288,   lower,   upper);   // 47
//  deut_dphi_pt_pos_neg_04_T  = new TH2F("deut_dphi_pt_pos_neg_04_T",  "deut_dphi_pt_pos_neg_04_T",    18, coarse_binning,     288,   lower,   upper);   // 48
//  deut_dphi_pt_neg_neg_04_T  = new TH2F("deut_dphi_pt_neg_neg_04_T",  "deut_dphi_pt_neg_neg_04_T",    18, coarse_binning,     288,   lower,   upper);   // 49
    deut_dphi_pt_pos_pos_05_T  = new TH2F("deut_dphi_pt_pos_pos_05_T",  "deut_dphi_pt_pos_pos_05_T",    18, coarse_binning,     288,   lower,   upper);   // 50
    deut_dphi_pt_pos_neg_05_T  = new TH2F("deut_dphi_pt_pos_neg_05_T",  "deut_dphi_pt_pos_neg_05_T",    18, coarse_binning,     288,   lower,   upper);   // 51
    deut_dphi_pt_neg_neg_05_T  = new TH2F("deut_dphi_pt_neg_neg_05_T",  "deut_dphi_pt_neg_neg_05_T",    18, coarse_binning,     288,   lower,   upper);   // 52
//  deut_dphi_pt_pos_pos_06_T  = new TH2F("deut_dphi_pt_pos_pos_06_T",  "deut_dphi_pt_pos_pos_06_T",    18, coarse_binning,     288,   lower,   upper);   // 53
//  deut_dphi_pt_pos_neg_06_T  = new TH2F("deut_dphi_pt_pos_neg_06_T",  "deut_dphi_pt_pos_neg_06_T",    18, coarse_binning,     288,   lower,   upper);   // 54
//  deut_dphi_pt_neg_neg_06_T  = new TH2F("deut_dphi_pt_neg_neg_06_T",  "deut_dphi_pt_neg_neg_06_T",    18, coarse_binning,     288,   lower,   upper);   // 55
    prot_dphi_pt_pos_pos_05_T  = new TH2F("prot_dphi_pt_pos_pos_05_T",  "prot_dphi_pt_pos_pos_05_T",    18,coarse_binning2,     288,   lower,   upper);   // 50a
    prot_dphi_pt_pos_neg_05_T  = new TH2F("prot_dphi_pt_pos_neg_05_T",  "prot_dphi_pt_pos_neg_05_T",    18,coarse_binning2,     288,   lower,   upper);   // 51a
    prot_dphi_pt_neg_neg_05_T  = new TH2F("prot_dphi_pt_neg_neg_05_T",  "prot_dphi_pt_neg_neg_05_T",    18,coarse_binning2,     288,   lower,   upper);   // 52a
    
    
//  deut_dphi_pt_pos_pos_04_A  = new TH2F("deut_dphi_pt_pos_pos_04_A",  "deut_dphi_pt_pos_pos_04_A",    18, coarse_binning,     288,   lower,   upper);   // 56
//  deut_dphi_pt_pos_neg_04_A  = new TH2F("deut_dphi_pt_pos_neg_04_A",  "deut_dphi_pt_pos_neg_04_A",    18, coarse_binning,     288,   lower,   upper);   // 57
//  deut_dphi_pt_neg_neg_04_A  = new TH2F("deut_dphi_pt_neg_neg_04_A",  "deut_dphi_pt_neg_neg_04_A",    18, coarse_binning,     288,   lower,   upper);   // 58
    deut_dphi_pt_pos_pos_05_A  = new TH2F("deut_dphi_pt_pos_pos_05_A",  "deut_dphi_pt_pos_pos_05_A",    18, coarse_binning,     288,   lower,   upper);   // 59
    deut_dphi_pt_pos_neg_05_A  = new TH2F("deut_dphi_pt_pos_neg_05_A",  "deut_dphi_pt_pos_neg_05_A",    18, coarse_binning,     288,   lower,   upper);   // 60
    deut_dphi_pt_neg_neg_05_A  = new TH2F("deut_dphi_pt_neg_neg_05_A",  "deut_dphi_pt_neg_neg_05_A",    18, coarse_binning,     288,   lower,   upper);   // 61
//  deut_dphi_pt_pos_pos_06_A  = new TH2F("deut_dphi_pt_pos_pos_06_A",  "deut_dphi_pt_pos_pos_06_A",    18, coarse_binning,     288,   lower,   upper);   // 62
//  deut_dphi_pt_pos_neg_06_A  = new TH2F("deut_dphi_pt_pos_neg_06_A",  "deut_dphi_pt_pos_neg_06_A",    18, coarse_binning,     288,   lower,   upper);   // 63
//  deut_dphi_pt_neg_neg_06_A  = new TH2F("deut_dphi_pt_neg_neg_06_A",  "deut_dphi_pt_neg_neg_06_A",    18, coarse_binning,     288,   lower,   upper);   // 64
    
//  deut_dphi_pt_pos_pos_04_B  = new TH2F("deut_dphi_pt_pos_pos_04_B",  "deut_dphi_pt_pos_pos_04_B",    18, coarse_binning,     288,   lower,   upper);   // 65
//  deut_dphi_pt_pos_neg_04_B  = new TH2F("deut_dphi_pt_pos_neg_04_B",  "deut_dphi_pt_pos_neg_04_B",    18, coarse_binning,     288,   lower,   upper);   // 66
//  deut_dphi_pt_neg_neg_04_B  = new TH2F("deut_dphi_pt_neg_neg_04_B",  "deut_dphi_pt_neg_neg_04_B",    18, coarse_binning,     288,   lower,   upper);   // 67
    deut_dphi_pt_pos_pos_05_B  = new TH2F("deut_dphi_pt_pos_pos_05_B",  "deut_dphi_pt_pos_pos_05_B",    18, coarse_binning,     288,   lower,   upper);   // 68
    deut_dphi_pt_pos_neg_05_B  = new TH2F("deut_dphi_pt_pos_neg_05_B",  "deut_dphi_pt_pos_neg_05_B",    18, coarse_binning,     288,   lower,   upper);   // 69
    deut_dphi_pt_neg_neg_05_B  = new TH2F("deut_dphi_pt_neg_neg_05_B",  "deut_dphi_pt_neg_neg_05_B",    18, coarse_binning,     288,   lower,   upper);   // 70
//  deut_dphi_pt_pos_pos_06_B  = new TH2F("deut_dphi_pt_pos_pos_06_B",  "deut_dphi_pt_pos_pos_06_B",    18, coarse_binning,     288,   lower,   upper);   // 71
//  deut_dphi_pt_pos_neg_06_B  = new TH2F("deut_dphi_pt_pos_neg_06_B",  "deut_dphi_pt_pos_neg_06_B",    18, coarse_binning,     288,   lower,   upper);   // 72
//  deut_dphi_pt_neg_neg_06_B  = new TH2F("deut_dphi_pt_neg_neg_06_B",  "deut_dphi_pt_neg_neg_06_B",    18, coarse_binning,     288,   lower,   upper);   // 73

    DCAxy_pos                  = new TH1F("DCAxy_pos",                  "DCAxy_pos",                   100,   -1.2,    1.2);                              // 74
    DCAxy_neg                  = new TH1F("DCAxy_neg",                  "DCAxy_neg",                   100,   -1.2,    1.2);                              // 75
    DCAz_pos                   = new TH1F("DCAz_pos",                   "DCAz_pos",                    100,   -1.2,    1.2);                              // 76
    DCAz_neg                   = new TH1F("DCAz_neg",                   "DCAz_neg",                    100,   -1.2,    1.2);                              // 77

    Double_t fine_binning[2001];

    Float_t moving_marker = 0.10;
    for(int i=0; i<1602; i++)
    {
	fine_binning[i] = moving_marker;
	moving_marker = moving_marker + fine_binning[i] * 0.005;
    }
    // display plots

    m2_pt_pos_fine             = new TH2F("m2_pt_pos_fine",             "m2_pt_pos_fine",              800,   fine_binning,    2400,    -1.0,     7.0);   // 78
    m2_pt_neg_fine             = new TH2F("m2_pt_neg_fine",             "m2_pt_neg_fine",              800,   fine_binning,    2400,    -1.0,     7.0);   // 79
    m2_pt_pos_TPC_fine         = new TH2F("m2_pt_pos_TPC_fine",         "m2_pt_pos_TPC_fine",          800,   fine_binning,    2400,    -1.0,     7.0);   // 80
    m2_pt_neg_TPC_fine         = new TH2F("m2_pt_neg_TPC_fine",         "m2_pt_neg_TPC_fine",          800,   fine_binning,    2400,    -1.0,     7.0);   // 81
    m2_pt_pos_cut_T_fine       = new TH2F("m2_pt_pos_cut_T_fine",       "m2_pt_pos_cut_T_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 82
    m2_pt_neg_cut_T_fine       = new TH2F("m2_pt_neg_cut_T_fine",       "m2_pt_neg_cut_T_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 83
    m2_pt_pos_cut_A_fine       = new TH2F("m2_pt_pos_cut_A_fine",       "m2_pt_pos_cut_A_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 84
    m2_pt_neg_cut_A_fine       = new TH2F("m2_pt_neg_cut_A_fine",       "m2_pt_neg_cut_A_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 85
    m2_pt_pos_cut_B_fine       = new TH2F("m2_pt_pos_cut_B_fine",       "m2_pt_pos_cut_B_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 86
    m2_pt_neg_cut_B_fine       = new TH2F("m2_pt_neg_cut_B_fine",       "m2_pt_neg_cut_B_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 87


    deut_m2_dphi_04_pt1        = new TH2F("deut_m2_dphi_04_pt1",        "deut_m2_dphi_04_pt1",          75,    2.5,    5.0,      36,   lower,   upper);   // 88
    deut_m2_dphi_04_pt2        = new TH2F("deut_m2_dphi_04_pt2",        "deut_m2_dphi_04_pt2",          75,    2.5,    5.0,      36,   lower,   upper);   // 89
    deut_m2_dphi_04_pt3        = new TH2F("deut_m2_dphi_04_pt3",        "deut_m2_dphi_04_pt3",          75,    2.5,    5.0,      36,   lower,   upper);   // 90
    deut_m2_dphi_04_pt4        = new TH2F("deut_m2_dphi_04_pt4",        "deut_m2_dphi_04_pt4",          75,    2.5,    5.0,      36,   lower,   upper);   // 91
    
    m2_pt_pos_TPC_prot_fine    = new TH2F("m2_pt_pos_TPC_prot_fine",    "m2_pt_pos_TPC_prot_fine",     800,   fine_binning,    2400,    -1.0,     7.0);   // 92
    m2_pt_neg_TPC_prot_fine    = new TH2F("m2_pt_neg_TPC_prot_fine",    "m2_pt_neg_TPC_prot_fine",     800,   fine_binning,    2400,    -1.0,     7.0);   // 93

    m2_pt_pos_cut_T_prot_fine  = new TH2F("m2_pt_pos_cut_T_prot_fine",  "m2_pt_pos_cut_T_prot_fine",   800,   fine_binning,    2400,    -1.0,     7.0);   // 94
    m2_pt_neg_cut_T_prot_fine  = new TH2F("m2_pt_neg_cut_T_prot_fine",  "m2_pt_neg_cut_T_prot_fine",   800,   fine_binning,    2400,    -1.0,     7.0);   // 95


    
  // objects added to output file

    fOutputList->Add(fHistPt);                      //  1
    fOutputList->Add(cent_ntracks);                 //  2
    
//  fOutputList->Add(m2_pt_pos);                    //  3
//  fOutputList->Add(m2_pt_neg);                    //  4
    fOutputList->Add(m2_pt_pos_TPC);                //  5
    fOutputList->Add(m2_pt_neg_TPC);                //  6
    fOutputList->Add(m2_pt_pos_cut_T);              //  7
    fOutputList->Add(m2_pt_neg_cut_T);              //  8
    fOutputList->Add(m2_pt_pos_cut_G);              //  9
    fOutputList->Add(m2_pt_neg_cut_G);              // 10
    fOutputList->Add(m2_pt_pos_cut_A);              // 11
    fOutputList->Add(m2_pt_neg_cut_A);              // 12
    fOutputList->Add(m2_pt_pos_cut_B);              // 13
    fOutputList->Add(m2_pt_neg_cut_B);              // 14

    if(run_mode == 0)
    {
	fOutputList->Add(m2_pt_pos_cut_T_prot);  //  7a
	fOutputList->Add(m2_pt_neg_cut_T_prot);  //  8a
    }
    
//  fOutputList->Add(m2_pt_pos_cut_with_trig_04);   // 15
//  fOutputList->Add(m2_pt_neg_cut_with_trig_04);   // 16
    fOutputList->Add(m2_pt_pos_cut_with_trig_05);   // 17
    fOutputList->Add(m2_pt_neg_cut_with_trig_05);   // 18
//  fOutputList->Add(m2_pt_pos_cut_with_trig_06);   // 19
//  fOutputList->Add(m2_pt_neg_cut_with_trig_06);   // 20

    fOutputList->Add(deut_phi_pt_pos_T);            // 21
    fOutputList->Add(deut_phi_pt_neg_T);            // 22

    if(run_mode == 0)
    {
	fOutputList->Add(prot_phi_pt_pos_T);            // 21a
	fOutputList->Add(prot_phi_pt_neg_T);            // 22a
    }
    
    fOutputList->Add(deut_phi_pt_pos_A);            // 23
    fOutputList->Add(deut_phi_pt_neg_A);            // 24
    fOutputList->Add(deut_phi_pt_pos_B);            // 25
    fOutputList->Add(deut_phi_pt_neg_B);            // 26

    fOutputList->Add(deut_per_event);               // 27
    if(run_mode == 0)
    {
	fOutputList->Add(prot_per_event);               // 27a
    }
    
//  fOutputList->Add(trig_04_per_event);            // 28
    fOutputList->Add(trig_05_per_event);            // 29
//  fOutputList->Add(trig_06_per_event);            // 30

//  fOutputList->Add(trig_04_phi_pt_pos);           // 31
//  fOutputList->Add(trig_04_phi_pt_neg);           // 32    
    fOutputList->Add(trig_05_phi_pt_pos);           // 33
    fOutputList->Add(trig_05_phi_pt_neg);           // 34
//  fOutputList->Add(trig_06_phi_pt_pos);           // 35
//  fOutputList->Add(trig_06_phi_pt_neg);           // 36

    fOutputList->Add(tof_phi_eta_pos);              // 37
    fOutputList->Add(tof_phi_eta_neg);              // 38
    fOutputList->Add(tof_phi_eta_pos_deut);         // 39
    fOutputList->Add(tof_phi_eta_neg_deut);         // 40

    fOutputList->Add(deut_pt_compare_pos);          // 41
    fOutputList->Add(deut_pt_compare_neg);          // 42
    fOutputList->Add(tpc_sector_fraction);          // 43
    fOutputList->Add(primary_vertex_z);             // 44
    fOutputList->Add(primary_vertex_z_cut1);        // 45
    fOutputList->Add(primary_vertex_z_cut2);        // 46
	
//  fOutputList->Add(deut_dphi_pt_pos_pos_04_T);    // 47
//  fOutputList->Add(deut_dphi_pt_pos_neg_04_T);    // 48
//  fOutputList->Add(deut_dphi_pt_neg_neg_04_T);    // 49
    fOutputList->Add(deut_dphi_pt_pos_pos_05_T);    // 50
    fOutputList->Add(deut_dphi_pt_pos_neg_05_T);    // 51
    fOutputList->Add(deut_dphi_pt_neg_neg_05_T);    // 52
//  fOutputList->Add(deut_dphi_pt_pos_pos_06_T);    // 53
//  fOutputList->Add(deut_dphi_pt_pos_neg_06_T);    // 54
//  fOutputList->Add(deut_dphi_pt_neg_neg_06_T);    // 55

    if(run_mode == 0)
    {
	fOutputList->Add(prot_dphi_pt_pos_pos_05_T);    // 50
	fOutputList->Add(prot_dphi_pt_pos_neg_05_T);    // 51
	fOutputList->Add(prot_dphi_pt_neg_neg_05_T);    // 52
    }
//  fOutputList->Add(deut_dphi_pt_pos_pos_04_A);    // 56
//  fOutputList->Add(deut_dphi_pt_pos_neg_04_A);    // 57
//  fOutputList->Add(deut_dphi_pt_neg_neg_04_A);    // 58
    fOutputList->Add(deut_dphi_pt_pos_pos_05_A);    // 59
    fOutputList->Add(deut_dphi_pt_pos_neg_05_A);    // 60
    fOutputList->Add(deut_dphi_pt_neg_neg_05_A);    // 61
//  fOutputList->Add(deut_dphi_pt_pos_pos_06_A);    // 62
//  fOutputList->Add(deut_dphi_pt_pos_neg_06_A);    // 63
//  fOutputList->Add(deut_dphi_pt_neg_neg_06_A);    // 64
 
//  fOutputList->Add(deut_dphi_pt_pos_pos_04_B);    // 65
//  fOutputList->Add(deut_dphi_pt_pos_neg_04_B);    // 66
//  fOutputList->Add(deut_dphi_pt_neg_neg_04_B);    // 67
    fOutputList->Add(deut_dphi_pt_pos_pos_05_B);    // 68
    fOutputList->Add(deut_dphi_pt_pos_neg_05_B);    // 69
    fOutputList->Add(deut_dphi_pt_neg_neg_05_B);    // 70
//  fOutputList->Add(deut_dphi_pt_pos_pos_06_B);    // 71
//  fOutputList->Add(deut_dphi_pt_pos_neg_06_B);    // 72
//  fOutputList->Add(deut_dphi_pt_neg_neg_06_B);    // 73

    fOutputList->Add(DCAxy_pos);                    // 74
    fOutputList->Add(DCAxy_neg);                    // 75
    fOutputList->Add(DCAz_pos);                     // 76
    fOutputList->Add(DCAz_neg);                     // 77

    fOutputList->Add(m2_pt_pos_fine);               // 78
    fOutputList->Add(m2_pt_neg_fine);               // 79
    fOutputList->Add(m2_pt_pos_TPC_fine);           // 80
    fOutputList->Add(m2_pt_neg_TPC_fine);           // 81 
    fOutputList->Add(m2_pt_pos_cut_T_fine);         // 82
    fOutputList->Add(m2_pt_neg_cut_T_fine);         // 83
    fOutputList->Add(m2_pt_pos_cut_A_fine);         // 84
    fOutputList->Add(m2_pt_neg_cut_A_fine);         // 85
    fOutputList->Add(m2_pt_pos_cut_B_fine);         // 86
    fOutputList->Add(m2_pt_neg_cut_B_fine);         // 87

    fOutputList->Add(deut_m2_dphi_04_pt1);          // 88
    fOutputList->Add(deut_m2_dphi_04_pt2);          // 89
    fOutputList->Add(deut_m2_dphi_04_pt3);          // 90
    fOutputList->Add(deut_m2_dphi_04_pt4);          // 91

    if(run_mode == 0)
    {
    fOutputList->Add(m2_pt_pos_TPC_prot_fine);      // 92
    fOutputList->Add(m2_pt_neg_TPC_prot_fine);      // 93

    fOutputList->Add(m2_pt_pos_cut_T_prot_fine);    // 94
    fOutputList->Add(m2_pt_neg_cut_T_prot_fine);    // 95
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

//    fPIDResponse->SetRecoPass(1);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserExec(Option_t *)
{


//    cout<<endl<<endl<<"         trial         "<<endl<<endl;
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7)   ;
    else return;

    const AliAODVertex *primVertex = fAOD->GetPrimaryVertex();
    Double_t pv = primVertex->GetZ();                             primary_vertex_z->Fill(pv);
    
    if(!fAnalysisUtils->IsVertexSelected2013pA(fAOD)) return;     primary_vertex_z_cut1->Fill(pv);

    if(fAnalysisUtils->IsPileUpSPD(fAOD)) return;                 primary_vertex_z_cut2->Fill(pv);
    
    Int_t iTracks(fAOD->GetNumberOfTracks());
    cent_ntracks->Fill(50.0, iTracks);

    int wide_cut = 0;

    int deut_track_num_T[200];     int deut_count_T            = 0;
    int deut_track_num_A[200];     int deut_count_A            = 0;
    int deut_track_num_B[200];     int deut_count_B            = 0;
    int deut_track_num_W[200];     int deut_count_W            = 0;
    int prot_track_num_T[500];     int prot_count_T            = 0;
    
//  int trig_04_track_num[200];    int trig_04_track_count     = 0;
    int trig_05_track_num[200];    int trig_05_track_count     = 0;
//  int trig_06_track_num[200];    int trig_06_track_count     = 0;

    // loop over all these tracks
    //
    // pull out track numbers for high-pt triggers and also deutons

//  Float_t max_04_pt = 0.0;    Int_t   max_04    = 0;
    Float_t max_05_pt = 0.0;    Int_t   max_05    = 0;
//  Float_t max_06_pt = 0.0;    Int_t   max_06    = 0;

    for(Int_t i(0); i < iTracks; i++)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track)                                                                      {    continue;    }


	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }   // cut 1
//	if(!track->IsGlobalConstrained())                                               {    continue;    }
//	if(!track->IsTPCConstrained())                                                  {    continue;    }
//	if(!track->IsPrimaryCandidate())                                                {    continue;    }
	
	Float_t pt   = track->Pt();	 //     if(pt < 0.95)                           {    continue;    }   // cut 2   // makes for fast processings, low momentum play no part of analysis anyway
        Float_t dedx = track->GetTPCsignal();	if(dedx > 1000)                         {    continue;    }   //
	Float_t eta  = track->Eta();	        if(TMath::Abs(eta) > 0.9)               {    continue;    }


	int is_global = 0;
	if(!track->IsGlobalConstrained())                                               {    is_global=1; }

	if(run_mode == 6  && track->IsGlobalConstrained())                              {    continue;    }


	Float_t phi           = track->Phi();
	if(phi <  -pio2){    phi = phi + twopi; }	if(phi <  -pio2){    phi = phi + twopi;   }
	if(phi > 3*pio2){    phi = phi - twopi;	}       if(phi > 3*pio2){    phi = phi - twopi;   }


	Float_t fraction = 0.0;
	fraction = track->GetTPCFoundFraction();
	if(run_mode == 3  &&  fraction < 0.90)                                          {    continue;    }  /////////////// OPTIONAL FOR SYSTEMATIC STUDY /////////////////
	tpc_sector_fraction->Fill(fraction);
	
	fHistPt->Fill(pt);

	// markA
	
//	if(pt >= 4.0  &&  pt < 6.0)	{   trig_04_track_num[trig_04_track_count] = i;	  if(pt>max_04_pt){ max_04_pt = pt; max_04 = trig_04_track_count;   }   trig_04_track_count++;   }
	if(pt >= 5.0              )	{   trig_05_track_num[trig_05_track_count] = i;	  if(pt>max_05_pt){ max_05_pt = pt; max_05 = trig_05_track_count;   }   trig_05_track_count++;   }
//	if(pt >= 6.0              )	{   trig_06_track_num[trig_06_track_count] = i;	  if(pt>max_06_pt){ max_06_pt = pt; max_06 = trig_06_track_count;   }   trig_06_track_count++;   }

	Short_t charge        = track->Charge();
//	if(pt >= 4.0  &&  pt < 6.0)	{   if(charge > 0)   trig_04_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_04_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 5.0              )	{   if(charge > 0)   trig_05_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_05_phi_pt_neg->Fill(pt, phi);	}
//	if(pt >= 6.0              )	{   if(charge > 0)   trig_06_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_06_phi_pt_neg->Fill(pt, phi);	}


//	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }
	

	
	Double_t nsigmaTPC = 999.0;	Double_t nsigmaTOF = 999.0;
	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);     /* && trk->IsOn(AliESDtrack::kTPCpid)*/;
	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
	if(!tpcIsOk)	                                                                {    continue;    }
	if(!tofIsOk)	                                                                {    continue;    }



	Double_t dca[2];          // 0 - dcar, 1 - dcaz -> output parameter
	Double_t cov[3];          // covariance matrix -> output parameter
	Double_t kmaxd = 10000.0;    // something big
	AliExternalTrackParam param;
	param.CopyFromVTrack(track);
	if(!param.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), kmaxd, dca, cov)) continue; // don’t accept tracks for which the DCA fit failed.
/////////  SYSTEMATIC STUDIES  /////////
	if(fabs(dca[0]) > 0.5)                                                          {    continue;    }      // originally not used, but 0.5 expected   'r' component (confirmed)
	if(fabs(dca[1]) > 1.0)                                                          {    continue;    }      // originally not used, but 1.0 expected   'z' component (cofirmed)
////////////////////////////////////////
	if(run_mode == 4  &&  fabs(dca[0]) > 0.1)                                       {    continue;    }      // used for systematic uncertainty calculation
	if(run_mode == 4  &&  fabs(dca[1]) > 0.1)                                       {    continue;    }
	
//	Float_t deltat        = tof_minus_tpion(track);
	float   mom           = track->P();

	Float_t deut_mean     = 0.0;
	Float_t deut_sigma    = 0.0;
	Float_t prot_mean     = 0.0;
	Float_t prot_sigma    = 0.0;
	Float_t m2tof         = get_mass_squared(track);
//	Float_t beta          = 0.0;
//	beta                  = Beta(track);

	Float_t TPC_PID       = 3.0;

	if(run_mode == 2){   TPC_PID = 2.0;   }
	
	if(charge > 0)
	{
//	    m2_pt_pos     ->Fill(pt, m2tof);
	    m2_pt_pos_fine->Fill(pt, m2tof);
	    if(pt >= 1.0  &&  pt < 2.0)   tof_phi_eta_pos->Fill(phi, eta);

	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < TPC_PID)
	    {
		m2_pt_pos_TPC     ->Fill(pt,  m2tof);
		m2_pt_pos_TPC_fine->Fill(pt,  m2tof);
			
		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
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
		    if(m2tof >= 2.5  &&  m2tof < 5.0)
		    {
			deut_track_num_W[deut_count_W] = i;
			deut_count_W++;
		    }
		    
		    if(run_mode != 5  &&  m2tof +1.2< deut_mean + cut_width * deut_sigma  &&   m2tof +1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_A[deut_count_A] = i;
			deut_count_A++;
			m2_pt_pos_cut_A->Fill(pt,m2tof);
			m2_pt_pos_cut_A_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_A->Fill(pt, deut_phi);
		    }
		    if(run_mode != 5  &&  m2tof -1.2< deut_mean + cut_width * deut_sigma  &&   m2tof -1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_B[deut_count_B] = i;
			deut_count_B++;
			m2_pt_pos_cut_B->Fill(pt,m2tof);
			m2_pt_pos_cut_B_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_B->Fill(pt, deut_phi);
		    }
		    if(run_mode == 5  &&  m2tof >=2.5  &&  m2tof < 3.0)
		    {
			deut_track_num_A[deut_count_A] = i;
			deut_count_A++;
			m2_pt_pos_cut_A->Fill(pt,m2tof);
			m2_pt_pos_cut_A_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_A->Fill(pt, deut_phi);
		    }
		    if(run_mode == 5  &&  m2tof >=4.0  &&  m2tof < 5.0)
		    {
			deut_track_num_B[deut_count_B] = i;
			deut_count_B++;
			m2_pt_pos_cut_B->Fill(pt,m2tof);
			m2_pt_pos_cut_B_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_B->Fill(pt, deut_phi);
		    }

		    if(m2tof >= 1.7  &&  m2tof < 5.5)    {   wide_cut++;   }
		}
	    }

	    if(run_mode == 0)
	    {
	    Double_t nSigmaTPCProt = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)4);  // 4 = proton
	    if(TMath::Abs(nSigmaTPCProt) < TPC_PID)
	    {
		m2_pt_pos_TPC_prot_fine->Fill(pt,  m2tof);
			
		if(pt >= 0.5  &&  pt < 2.5)
		{
		    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
		    prot_mean = fit_prot_curve->Eval(pt);
		    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
		    prot_sigma = fit_prot_curve->Eval(pt);

		    if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
		    {
			prot_track_num_T[prot_count_T] = i;
			prot_count_T++;
			m2_pt_pos_cut_T_prot->Fill(pt,m2tof);
			m2_pt_pos_cut_T_prot_fine->Fill(pt,m2tof);

			Float_t prot_phi = track->Phi();
			if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }	if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }
			if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }	if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }
			prot_phi_pt_pos_T->Fill(pt, prot_phi);
		    }
		}
	    }
	    }
	    
	}    //   end of pos charge if statement
	    
	else if(charge < 0)
	{
//	    m2_pt_neg     ->Fill(pt, m2tof);
	    m2_pt_neg_fine->Fill(pt, m2tof);
	    if(pt >= 1.0  &&  pt < 2.0)   tof_phi_eta_neg->Fill(phi, eta);

	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < TPC_PID)
	    {
		m2_pt_neg_TPC     ->Fill(pt, m2tof);
		m2_pt_neg_TPC_fine->Fill(pt, m2tof);

		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][1][w]);   }
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

		    
		    if(m2tof < deut_mean + 3.0 * deut_sigma  &&   m2tof > deut_mean - 3.0 * deut_sigma)  // used for comparison allowing use of eff+acc cuts
		    {
			if(is_global == 1)           {   m2_pt_neg_cut_G->Fill(pt,m2tof);                          }
		    }
		    if(m2tof >= 2.5  &&  m2tof < 5.0)
		    {
			deut_track_num_W[deut_count_W] = i;
			deut_count_W++;
		    }

		    if(run_mode != 5  &&  m2tof +1.2< deut_mean + cut_width * deut_sigma  &&   m2tof +1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_A[deut_count_A] = i;
			deut_count_A++;
			m2_pt_neg_cut_A->Fill(pt,m2tof);
			m2_pt_neg_cut_A_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_A->Fill(pt, deut_phi);
		    }
		    if(run_mode != 5  &&  m2tof -1.2< deut_mean + cut_width * deut_sigma  &&   m2tof -1.2> deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_B[deut_count_B] = i;
			deut_count_B++;
			m2_pt_neg_cut_B->Fill(pt,m2tof);
			m2_pt_neg_cut_B_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_B->Fill(pt, deut_phi);
		    }
		    if(run_mode == 5  &&  m2tof >=2.5  &&  m2tof < 3.0)
		    {
			deut_track_num_A[deut_count_A] = i;
			deut_count_A++;
			m2_pt_neg_cut_A->Fill(pt,m2tof);
			m2_pt_neg_cut_A_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_A->Fill(pt, deut_phi);
		    }
		    if(run_mode == 5  &&  m2tof >=4.0  &&  m2tof < 5.0)
		    {
			deut_track_num_B[deut_count_B] = i;
			deut_count_B++;
			m2_pt_neg_cut_B->Fill(pt,m2tof);
			m2_pt_neg_cut_B_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_B->Fill(pt, deut_phi);
		    }

		    if(m2tof >= 1.7  &&  m2tof < 5.5)    {   wide_cut++;   }
		}
	    }

	    if(run_mode == 0)
	    {
	    Double_t nSigmaTPCProt = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)4);  // 4 = proton
	    if(TMath::Abs(nSigmaTPCProt) < TPC_PID)
	    {
		m2_pt_neg_TPC_prot_fine->Fill(pt,  m2tof);
			
		if(pt >= 0.4  &&  pt < 2.5)
		{
		    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[1][0][w]);   }
		    prot_mean = fit_prot_curve->Eval(pt);
		    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[1][1][w]);   }
		    prot_sigma = fit_prot_curve->Eval(pt);

		    if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
		    {
			prot_track_num_T[prot_count_T] = i;
			prot_count_T++;
			m2_pt_neg_cut_T_prot->Fill(pt,m2tof);
			m2_pt_neg_cut_T_prot_fine->Fill(pt,m2tof);

			Float_t prot_phi = track->Phi();
			if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }	if(prot_phi <  -pio2){    prot_phi = prot_phi + twopi;   }
			if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }	if(prot_phi > 3*pio2){    prot_phi = prot_phi - twopi;   }
			prot_phi_pt_neg_T->Fill(pt, prot_phi);
		    }
		}
	    }
	    }
	    
	}    //   end of neg charge if statement
    }        //   end of track loop


    deut_per_event->Fill(deut_count_T);
    prot_per_event->Fill(prot_count_T);
//  if(trig_04_track_count > 0)    {	trig_04_per_event->Fill(trig_04_track_count);    }
    if(trig_05_track_count > 0)    {	trig_05_per_event->Fill(trig_05_track_count);    }
//  if(trig_06_track_count > 0)    {	trig_06_per_event->Fill(trig_06_track_count);    }

    int begin = 0;
    /*
    if(deut_count_W > 0  &&  trig_04_track_count > 0)
    {
	// check if the >3.0 GeV trigger coincices with the deuteron
	int matcher         = 0;
	int unique_trigger  = 0;
	for(int i=begin; i<trig_04_track_count; i++)  // trigger loop
	{
	    matcher = 0;
	    int H = trig_04_track_num[i];
	    for(int j=0; j<deut_count_W; j++)
	    {
		int A = deut_track_num_W[j];
		if(A == H)  matcher = 1;
	    }
	    if(matcher == 0){   unique_trigger = 1;   }
	}

	if(unique_trigger == 1)
	{
	    for(int j=0; j<deut_count_W; j++)
	    {
		int A = deut_track_num_W[j];
		AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		if(trackA)
		{
		    Float_t m2tof       = get_mass_squared(trackA);		    
		    Short_t charge_A    = trackA->Charge();
		    Float_t pt_A        = trackA->Pt();
		    if     (charge_A > 0){   m2_pt_pos_cut_with_trig_04->Fill(pt_A, m2tof);	}
		    else if(charge_A < 0){   m2_pt_neg_cut_with_trig_04->Fill(pt_A, m2tof);	}
		}
	    }
	}
    }
    */
    if(deut_count_W > 0  &&  trig_05_track_count > 0)
    {
	for(int j=0; j<deut_count_W; j++)
	{
	    int A               = deut_track_num_W[j];
	    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	    if(trackA)
	    {
		Float_t m2tof       = get_mass_squared(trackA);		    
		Short_t charge_A    = trackA->Charge();
		Float_t pt_A        = trackA->Pt();
		if     (charge_A > 0){   m2_pt_pos_cut_with_trig_05->Fill(pt_A, m2tof);	}
		else if(charge_A < 0){   m2_pt_neg_cut_with_trig_05->Fill(pt_A, m2tof);	}
    }   }   }

    /*
    if(deut_count_W > 0  &&  trig_06_track_count > 0)
    {
	for(int j=0; j<deut_count_W; j++)
	{
	    int A               = deut_track_num_W[j];
	    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	    if(trackA)
	    {
		Float_t m2tof       = get_mass_squared(trackA);		    
		Short_t charge_A    = trackA->Charge();
		Float_t pt_A        = trackA->Pt();
		if     (charge_A > 0){   m2_pt_pos_cut_with_trig_06->Fill(pt_A, m2tof);	}
		else if(charge_A < 0){   m2_pt_neg_cut_with_trig_06->Fill(pt_A, m2tof);	}
    }   }   }    
    */

    // mark
////////////////////////////////////////////////////////////////////////////////////////////////////    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    if(deut_count_W > 0  &&  trig_04_track_count > 0)
    {
	for(int i=begin; i<trig_04_track_count; i++)  // trigger loop
	{
	    int H               = trig_04_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t phi_H       = trackH->Phi();	    
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    
	    for(int j=0; j<deut_count_W; j++)
	    {
		int A               = deut_track_num_W[j];
		if(A != H)
		{
		    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {
			Float_t pt_A        = trackA->Pt();
			Float_t phi_A       = trackA->Phi();
			if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
			if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
			Float_t Sdphi       =  phi_A - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }

			Float_t m2tof_A     = get_mass_squared(trackA);
			if     (pt_A >= 1.0  &&  pt_A < 1.5  &&  TMath::Abs(Sdphi) < 0.698){   deut_m2_dphi_04_pt1->Fill(m2tof_A, Sdphi);   }
			else if(pt_A >= 1.5  &&  pt_A < 2.0  &&  TMath::Abs(Sdphi) < 0.698){   deut_m2_dphi_04_pt2->Fill(m2tof_A, Sdphi);   }
			else if(pt_A >= 2.0  &&  pt_A < 3.0  &&  TMath::Abs(Sdphi) < 0.698){   deut_m2_dphi_04_pt3->Fill(m2tof_A, Sdphi);   }
			else if(pt_A >= 3.0  &&  pt_A < 4.0  &&  TMath::Abs(Sdphi) < 0.698){   deut_m2_dphi_04_pt4->Fill(m2tof_A, Sdphi);   }			    
    }   }   }   }   }
*/			  



    ////////////////////////////////////////////////////////////////////////////////////////////////////    
////////////////////////////////////////////////////////////////////////////////////////////////////

    /*
    begin = 0;
    if(deut_count_T > 0  &&  trig_04_track_count > 0  &&  do_lead_only == 1)    {	begin = max_04;   trig_04_track_count = max_04 + 1;    }
    if(deut_count_T > 0  &&  trig_04_track_count > 0)
    {
	for(int i=begin; i<trig_04_track_count; i++)  // trigger loop
	{
	    int H               = trig_04_track_num[i];
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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_04_T->Fill(pt_A, Sdphi);      }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_04_T->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_04_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_04_T->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }
    */
    

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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_T->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }

    if(run_mode == 0)
    {
    begin = 0;
    if(prot_count_T > 0  &&  trig_05_track_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_track_count = max_05 + 1;    }
    if(prot_count_T > 0  &&  trig_05_track_count > 0)
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

	    for(int j=0; j<prot_count_T; j++)
	    {
		int A               = prot_track_num_T[j];
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
		    if     (charge_H > 0){  if     (charge_A > 0){    prot_dphi_pt_pos_pos_05_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    prot_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    prot_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    prot_dphi_pt_neg_neg_05_T->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }
    }
/*
    begin = 0;
    if(deut_count_T > 0  &&  trig_06_track_count > 0  &&  do_lead_only == 1)    {   begin = max_06;   trig_06_track_count = max_06 + 1;    }
    if(deut_count_T > 0  &&  trig_06_track_count > 0)
    {
	for(int i=begin; i<trig_06_track_count; i++)  // trigger loop
	{
	    int H               = trig_06_track_num[i];
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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_06_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_06_T->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_06_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_06_T->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    begin = 0;
    if(deut_count_A > 0  &&  trig_04_track_count > 0  &&  do_lead_only == 1)    {	begin = max_04;   trig_04_track_count = max_04 + 1;    }
    if(deut_count_A > 0  &&  trig_04_track_count > 0)
    {
	for(int i=begin; i<trig_04_track_count; i++)  // trigger loop
	{
	    int H               = trig_04_track_num[i];
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
		    if(     charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_04_A->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_04_A->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_04_A->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_04_A->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }
*/
    

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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_A->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_A->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_A->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_A->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }

/*
    begin = 0;
    if(deut_count_A > 0  &&  trig_06_track_count > 0  &&  do_lead_only == 1)    {   begin = max_06;   trig_06_track_count = max_06 + 1;    }
    if(deut_count_A > 0  &&  trig_06_track_count > 0)
    {
	for(int i=begin; i<trig_06_track_count; i++)  // trigger loop
	{
	    int H               = trig_06_track_num[i];
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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_06_A->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_06_A->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_06_A->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_06_A->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }    
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    begin = 0;
    if(deut_count_B > 0  &&  trig_04_track_count > 0  &&  do_lead_only == 1)    {	begin = max_04;   trig_04_track_count = max_04 + 1;    }
    if(deut_count_B > 0  &&  trig_04_track_count > 0)
    {
	for(int i=begin; i<trig_04_track_count; i++)  // trigger loop
	{
	    int H               = trig_04_track_num[i];
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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_04_B->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_04_B->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_04_B->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_04_B->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }
*/
    

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
		    if(charge_H > 0)	 {  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_B->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_B->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_B->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_B->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }

/*
    begin = 0;
    if(deut_count_B > 0  &&  trig_06_track_count > 0  &&  do_lead_only == 1)    {   begin = max_06;   trig_06_track_count = max_06 + 1;    }
    if(deut_count_B > 0  &&  trig_06_track_count > 0)
    {
	for(int i=begin; i<trig_06_track_count; i++)  // trigger loop
	{
	    int H               = trig_06_track_num[i];
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
		    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_06_B->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_pos_neg_06_B->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_06_B->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    deut_dphi_pt_neg_neg_06_B->Fill(pt_A, Sdphi);	 }    }
    }   }   }   }   }
*/
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file







    
    AliAODHandler *oh = (AliAODHandler*)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if (oh)
      oh->SetFillAOD(kFALSE);

//  if ((deut_count>=1)  &&  (trig_04_track_count>0)  &&  oh)
//  {
    if ((wide_cut >= 1) && oh)
    {
	oh->SetFillAOD(kTRUE);
	AliAODEvent *eout = dynamic_cast<AliAODEvent*>(oh->GetAOD());
	AliAODEvent *evin = dynamic_cast<AliAODEvent*>(InputEvent());
	TTree *tout = oh->GetTree();
	if (tout)
	{
	    TList *lout = tout->GetUserInfo();
	    if (lout->FindObject("alirootVersion")==0)
	    {
		TList *lin = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetUserInfo();
		for (Int_t jj=0;jj<lin->GetEntries()-1;++jj)
		{ 
		    lout->Add(lin->At(jj)->Clone(lin->At(jj)->GetName()));
		}
	    }
	}
	
	if (1) {   AliAODHeader    *out =            (AliAODHeader*)eout->GetHeader();           AliAODHeader    *in = (AliAODHeader*)evin->GetHeader();  	    *out = *in;                  }
	if (1) {   AliTOFHeader    *out = const_cast<AliTOFHeader*>(eout->GetTOFHeader()); const AliTOFHeader    *in =                evin->GetTOFHeader();	    *out = *in;                  }
	if (1) {   AliAODVZERO     *out =                           eout->GetVZEROData();        AliAODVZERO     *in =                evin->GetVZEROData();	    *out = *in;                  }
	if (1) {   AliAODTZERO     *out =                           eout->GetTZEROData();        AliAODTZERO     *in =                evin->GetTZEROData(); 	    *out = *in;                  }
	if (1) {   TClonesArray    *out =                           eout->GetTracks();	         TClonesArray    *in =                evin->GetTracks();	new (out) TClonesArray(*in);     }
	if (1) {   TClonesArray    *out =                           eout->GetVertices();         TClonesArray    *in =                evin->GetVertices();      new (out) TClonesArray(*in);     }
	if (1) {   TClonesArray    *out =                           eout->GetCaloClusters();     TClonesArray    *in =                evin->GetCaloClusters();  new (out) TClonesArray(*in);     }
	if (1) {   AliAODCaloCells *out =                           eout->GetEMCALCells();       AliAODCaloCells *in =                evin->GetEMCALCells();    new (out) AliAODCaloCells(*in);  }
    }


    
}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
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
Double_t AliAnalysisTaskCorPIDTOFQA::get_deut_tof_pt(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFQA::get_deut_tof_p(AliAODTrack *track)
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
