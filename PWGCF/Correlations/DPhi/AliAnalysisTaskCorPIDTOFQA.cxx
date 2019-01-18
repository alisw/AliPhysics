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
#include "AliPIDResponse.h"
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

int      eta_limiting = 1;
Double_t cut_width    = 2.0;
short    do_lead_only = 0;     // 0 = all,  1 = only leading

ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),


    deut_eta_pos_pt1(0),
    deut_eta_neg_pt1(0),
    deut_eta_pos_pt2(0),
    deut_eta_neg_pt2(0),
    deut_eta_pos_pt3(0),
    deut_eta_neg_pt3(0),
    deut_eta_pos_pt4(0),
    deut_eta_neg_pt4(0),
    deut_eta_pos_pt5(0),
    deut_eta_neg_pt5(0),
    
    hadr_eta_pos(0),
    hadr_eta_neg(0),

    deut_hadr_dphi36_deta20_pt1(0),
    deut_hadr_dphi36_deta20_pt2(0),
    deut_hadr_dphi36_deta20_pt3(0),
    deut_hadr_dphi36_deta20_pt4(0),
    deut_hadr_dphi36_deta20_pt5(0),

    deut_hadr_dphi36_deta21_pt1(0),
    deut_hadr_dphi36_deta21_pt2(0),
    deut_hadr_dphi36_deta21_pt3(0),
    deut_hadr_dphi36_deta21_pt4(0),
    deut_hadr_dphi36_deta21_pt5(0),
    
    multiplicity_w_deut_1(0),          //  n1
    multiplicity_w_deut_2(0),          //  n2
    multiplicity_w_deut_3(0),          //  n3
    multiplicity_w_deut_4(0),          //  n4
    multiplicity_w_deut_5(0),          //  n5
    multiplicity_w_trig_and_deut_1(0), //  n6
    multiplicity_w_trig_and_deut_2(0), //  n7
    multiplicity_w_trig_and_deut_3(0), //  n8
    multiplicity_w_trig_and_deut_4(0), //  n9
    multiplicity_w_trig_and_deut_5(0), //  n10

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

    m2_pt_pos_cut_sa0(0),          // 11
    m2_pt_neg_cut_sa0(0),          // 12
    m2_pt_pos_cut_sa1(0),          // 13
    m2_pt_neg_cut_sa1(0),          // 14
    m2_pt_pos_cut_sa2(0),          // 15
    m2_pt_neg_cut_sa2(0),          // 16
    m2_pt_pos_cut_sa3(0),          // 17
    m2_pt_neg_cut_sa3(0),          // 18
    m2_pt_pos_cut_sa4(0),          // 19
    m2_pt_neg_cut_sa4(0),          // 20
    m2_pt_pos_cut_sa5(0),          // 21
    m2_pt_neg_cut_sa5(0),          // 22
    m2_pt_pos_cut_sa6(0),          // 19
    m2_pt_neg_cut_sa6(0),          // 20
    m2_pt_pos_cut_sa7(0),          // 21
    m2_pt_neg_cut_sa7(0),          // 22
    
    m2_pt_pos_cut_sb0(0),          // 23
    m2_pt_neg_cut_sb0(0),          // 24
    m2_pt_pos_cut_sb1(0),          // 25
    m2_pt_neg_cut_sb1(0),          // 26
    m2_pt_pos_cut_sb2(0),          // 27
    m2_pt_neg_cut_sb2(0),          // 28
    m2_pt_pos_cut_sb3(0),          // 29
    m2_pt_neg_cut_sb3(0),          // 30
    m2_pt_pos_cut_sb4(0),          // 31
    m2_pt_neg_cut_sb4(0),          // 32
    m2_pt_pos_cut_sb5(0),          // 33
    m2_pt_neg_cut_sb5(0),          // 34
    m2_pt_pos_cut_sb6(0),          // 19
    m2_pt_neg_cut_sb6(0),          // 20
    m2_pt_pos_cut_sb7(0),          // 21
    m2_pt_neg_cut_sb7(0),          // 22

    m2_pt_pos_cut_T_prot(0),       // 35
    m2_pt_neg_cut_T_prot(0),       // 36
	
    m2_pt_pos_cut_with_trig_05(0), // 37
    m2_pt_neg_cut_with_trig_05(0), // 38

    deut_phi_pt_pos_T(0),          // 39
    deut_phi_pt_neg_T(0),          // 40

    prot_phi_pt_pos_T(0),          // 41
    prot_phi_pt_neg_T(0),          // 42

    deut_phi_pt_pos_sa0(0),        // 43
    deut_phi_pt_neg_sa0(0),        // 44
    deut_phi_pt_pos_sa1(0),        // 45
    deut_phi_pt_neg_sa1(0),        // 46
    deut_phi_pt_pos_sa2(0),        // 47
    deut_phi_pt_neg_sa2(0),        // 48
    deut_phi_pt_pos_sa3(0),        // 49
    deut_phi_pt_neg_sa3(0),        // 50
    deut_phi_pt_pos_sa4(0),        // 51
    deut_phi_pt_neg_sa4(0),        // 52
    deut_phi_pt_pos_sa5(0),        // 53
    deut_phi_pt_neg_sa5(0),        // 54
    deut_phi_pt_pos_sa6(0),        // 51
    deut_phi_pt_neg_sa6(0),        // 52
    deut_phi_pt_pos_sa7(0),        // 53
    deut_phi_pt_neg_sa7(0),        // 54
    
    
    deut_phi_pt_pos_sb0(0),        // 55
    deut_phi_pt_neg_sb0(0),        // 56
    deut_phi_pt_pos_sb1(0),        // 57
    deut_phi_pt_neg_sb1(0),        // 58
    deut_phi_pt_pos_sb2(0),        // 59
    deut_phi_pt_neg_sb2(0),        // 60
    deut_phi_pt_pos_sb3(0),        // 61
    deut_phi_pt_neg_sb3(0),        // 62
    deut_phi_pt_pos_sb4(0),        // 63
    deut_phi_pt_neg_sb4(0),        // 64
    deut_phi_pt_pos_sb5(0),        // 65
    deut_phi_pt_neg_sb5(0),        // 66
    deut_phi_pt_pos_sb6(0),        // 63
    deut_phi_pt_neg_sb6(0),        // 64
    deut_phi_pt_pos_sb7(0),        // 65
    deut_phi_pt_neg_sb7(0),        // 66

    hadr_per_event(0),
    deut_per_event(0),             // 67
    prot_per_event(0),             // 68
    deut_per_event_with_trig(0),
    trig_05_per_event(0),          // 69
    trig_05_phi_pt_pos(0),         // 70
    trig_05_phi_pt_neg(0),         // 71
    
    tof_phi_eta_pos(0),            // 72
    tof_phi_eta_neg(0),            // 73
    tof_phi_eta_pos_deut(0),       // 74
    tof_phi_eta_neg_deut(0),       // 75

    deut_pt_compare_pos(0),        // 76
    deut_pt_compare_neg(0),        // 77
    tpc_sector_fraction(0),        // 78
    primary_vertex_z(0),           // 79
    primary_vertex_z_cut1(0),      // 80
    primary_vertex_z_cut2(0),      // 81
	
    deut_dphi_pt_pos_pos_05_T(0),  // 82
    deut_dphi_pt_pos_neg_05_T(0),  // 83
    deut_dphi_pt_neg_neg_05_T(0),  // 84
    
    prot_dphi_pt_pos_pos_05_T(0),  // 85
    prot_dphi_pt_pos_neg_05_T(0),  // 86
    prot_dphi_pt_neg_neg_05_T(0),  // 87

    deut_dphi_pt_pos_pos_05_sa0(0),// 88
    deut_dphi_pt_pos_neg_05_sa0(0),// 89
    deut_dphi_pt_neg_neg_05_sa0(0),// 90
    deut_dphi_pt_pos_pos_05_sa1(0),// 91
    deut_dphi_pt_pos_neg_05_sa1(0),// 92
    deut_dphi_pt_neg_neg_05_sa1(0),// 93
    deut_dphi_pt_pos_pos_05_sa2(0),// 94
    deut_dphi_pt_pos_neg_05_sa2(0),// 95
    deut_dphi_pt_neg_neg_05_sa2(0),// 96
    deut_dphi_pt_pos_pos_05_sa3(0),// 97
    deut_dphi_pt_pos_neg_05_sa3(0),// 98
    deut_dphi_pt_neg_neg_05_sa3(0),// 99
    deut_dphi_pt_pos_pos_05_sa4(0),//100
    deut_dphi_pt_pos_neg_05_sa4(0),//101
    deut_dphi_pt_neg_neg_05_sa4(0),//102
    deut_dphi_pt_pos_pos_05_sa5(0),//103
    deut_dphi_pt_pos_neg_05_sa5(0),//104
    deut_dphi_pt_neg_neg_05_sa5(0),//105
    deut_dphi_pt_pos_pos_05_sa6(0),//100
    deut_dphi_pt_pos_neg_05_sa6(0),//101
    deut_dphi_pt_neg_neg_05_sa6(0),//102
    deut_dphi_pt_pos_pos_05_sa7(0),//103
    deut_dphi_pt_pos_neg_05_sa7(0),//104
    deut_dphi_pt_neg_neg_05_sa7(0),//105    
    
    deut_dphi_pt_pos_pos_05_sb0(0),//106
    deut_dphi_pt_pos_neg_05_sb0(0),//107
    deut_dphi_pt_neg_neg_05_sb0(0),//108
    deut_dphi_pt_pos_pos_05_sb1(0),//109
    deut_dphi_pt_pos_neg_05_sb1(0),//110
    deut_dphi_pt_neg_neg_05_sb1(0),//111
    deut_dphi_pt_pos_pos_05_sb2(0),//112
    deut_dphi_pt_pos_neg_05_sb2(0),//113
    deut_dphi_pt_neg_neg_05_sb2(0),//114
    deut_dphi_pt_pos_pos_05_sb3(0),//115
    deut_dphi_pt_pos_neg_05_sb3(0),//116
    deut_dphi_pt_neg_neg_05_sb3(0),//117
    deut_dphi_pt_pos_pos_05_sb4(0),//118
    deut_dphi_pt_pos_neg_05_sb4(0),//119
    deut_dphi_pt_neg_neg_05_sb4(0),//120
    deut_dphi_pt_pos_pos_05_sb5(0),//121
    deut_dphi_pt_pos_neg_05_sb5(0),//122
    deut_dphi_pt_neg_neg_05_sb5(0),//123
    deut_dphi_pt_pos_pos_05_sb6(0),//118
    deut_dphi_pt_pos_neg_05_sb6(0),//119
    deut_dphi_pt_neg_neg_05_sb6(0),//120
    deut_dphi_pt_pos_pos_05_sb7(0),//121
    deut_dphi_pt_pos_neg_05_sb7(0),//122
    deut_dphi_pt_neg_neg_05_sb7(0),//123

    DCAxy_pos(0),                  //124
    DCAxy_neg(0),                  //125
    DCAz_pos(0),                   //126
    DCAz_neg(0),                   //127
    
    m2_pt_pos_fine(0),             //128
    m2_pt_neg_fine(0),             //129
    m2_pt_pos_TPC_fine(0),         //130
    m2_pt_neg_TPC_fine(0),         //131 
    m2_pt_pos_cut_T_fine(0),       //132
    m2_pt_neg_cut_T_fine(0),       //133

    m2_pt_pos_cut_sb0_fine(0),     //134
    m2_pt_neg_cut_sb0_fine(0),     //135
    m2_pt_pos_cut_sb1_fine(0),     //136
    m2_pt_neg_cut_sb1_fine(0),     //137
    m2_pt_pos_cut_sb2_fine(0),     //138
    m2_pt_neg_cut_sb2_fine(0),     //139
    m2_pt_pos_cut_sb3_fine(0),     //140
    m2_pt_neg_cut_sb3_fine(0),     //141
    m2_pt_pos_cut_sb4_fine(0),     //142
    m2_pt_neg_cut_sb4_fine(0),     //143
    m2_pt_pos_cut_sb5_fine(0),     //144
    m2_pt_neg_cut_sb5_fine(0),     //145
    m2_pt_pos_cut_sb6_fine(0),     //142
    m2_pt_neg_cut_sb6_fine(0),     //143
    m2_pt_pos_cut_sb7_fine(0),     //144
    m2_pt_neg_cut_sb7_fine(0),     //145
    
    m2_pt_pos_cut_sa0_fine(0),     //146
    m2_pt_neg_cut_sa0_fine(0),     //147
    m2_pt_pos_cut_sa1_fine(0),     //148
    m2_pt_neg_cut_sa1_fine(0),     //149
    m2_pt_pos_cut_sa2_fine(0),     //150
    m2_pt_neg_cut_sa2_fine(0),     //151
    m2_pt_pos_cut_sa3_fine(0),     //152
    m2_pt_neg_cut_sa3_fine(0),     //153
    m2_pt_pos_cut_sa4_fine(0),     //154
    m2_pt_neg_cut_sa4_fine(0),     //155
    m2_pt_pos_cut_sa5_fine(0),     //156
    m2_pt_neg_cut_sa5_fine(0),     //157
    m2_pt_pos_cut_sa6_fine(0),     //154
    m2_pt_neg_cut_sa6_fine(0),     //155
    m2_pt_pos_cut_sa7_fine(0),     //156
    m2_pt_neg_cut_sa7_fine(0),     //157

    m2_pt_pos_TPC_prot_fine(0),    //158
    m2_pt_neg_TPC_prot_fine(0),    //159

    m2_pt_pos_cut_T_prot_fine(0),  //160
    m2_pt_neg_cut_T_prot_fine(0),  //161

    deut_deut_ptA_ptB_same(0),
    deut_deut_ptA_ptB_diff(0),
    deut_deut_dpt_dcharge(0),
    deut_deut_dphi_deta_same(0),
    deut_deut_dphi_deta_diff(0),
    deut_deut_inv_mass(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),


    deut_eta_pos_pt1(0),
    deut_eta_neg_pt1(0),
    deut_eta_pos_pt2(0),
    deut_eta_neg_pt2(0),
    deut_eta_pos_pt3(0),
    deut_eta_neg_pt3(0),
    deut_eta_pos_pt4(0),
    deut_eta_neg_pt4(0),
    deut_eta_pos_pt5(0),
    deut_eta_neg_pt5(0),

    hadr_eta_pos(0),
    hadr_eta_neg(0),

    deut_hadr_dphi36_deta20_pt1(0),
    deut_hadr_dphi36_deta20_pt2(0),
    deut_hadr_dphi36_deta20_pt3(0),
    deut_hadr_dphi36_deta20_pt4(0),
    deut_hadr_dphi36_deta20_pt5(0),

    deut_hadr_dphi36_deta21_pt1(0),
    deut_hadr_dphi36_deta21_pt2(0),
    deut_hadr_dphi36_deta21_pt3(0),
    deut_hadr_dphi36_deta21_pt4(0),
    deut_hadr_dphi36_deta21_pt5(0),

    multiplicity_w_deut_1(0),          //  n1
    multiplicity_w_deut_2(0),          //  n2
    multiplicity_w_deut_3(0),          //  n3
    multiplicity_w_deut_4(0),          //  n4
    multiplicity_w_deut_5(0),          //  n5
    multiplicity_w_trig_and_deut_1(0), //  n6
    multiplicity_w_trig_and_deut_2(0), //  n7
    multiplicity_w_trig_and_deut_3(0), //  n8
    multiplicity_w_trig_and_deut_4(0), //  n9
    multiplicity_w_trig_and_deut_5(0), //  n10

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

    m2_pt_pos_cut_sa0(0),          // 11
    m2_pt_neg_cut_sa0(0),          // 12
    m2_pt_pos_cut_sa1(0),          // 13
    m2_pt_neg_cut_sa1(0),          // 14
    m2_pt_pos_cut_sa2(0),          // 15
    m2_pt_neg_cut_sa2(0),          // 16
    m2_pt_pos_cut_sa3(0),          // 17
    m2_pt_neg_cut_sa3(0),          // 18
    m2_pt_pos_cut_sa4(0),          // 19
    m2_pt_neg_cut_sa4(0),          // 20
    m2_pt_pos_cut_sa5(0),          // 21
    m2_pt_neg_cut_sa5(0),          // 22
    m2_pt_pos_cut_sa6(0),          // 19
    m2_pt_neg_cut_sa6(0),          // 20
    m2_pt_pos_cut_sa7(0),          // 21
    m2_pt_neg_cut_sa7(0),          // 22
									   
    m2_pt_pos_cut_sb0(0),          // 23
    m2_pt_neg_cut_sb0(0),          // 24
    m2_pt_pos_cut_sb1(0),          // 25
    m2_pt_neg_cut_sb1(0),          // 26
    m2_pt_pos_cut_sb2(0),          // 27
    m2_pt_neg_cut_sb2(0),          // 28
    m2_pt_pos_cut_sb3(0),          // 29
    m2_pt_neg_cut_sb3(0),          // 30
    m2_pt_pos_cut_sb4(0),          // 31
    m2_pt_neg_cut_sb4(0),          // 32
    m2_pt_pos_cut_sb5(0),          // 33
    m2_pt_neg_cut_sb5(0),          // 34
    m2_pt_pos_cut_sb6(0),          // 31
    m2_pt_neg_cut_sb6(0),          // 32
    m2_pt_pos_cut_sb7(0),          // 33
    m2_pt_neg_cut_sb7(0),          // 34
		
    m2_pt_pos_cut_T_prot(0),       // 35
    m2_pt_neg_cut_T_prot(0),       // 36
	
    m2_pt_pos_cut_with_trig_05(0), // 37
    m2_pt_neg_cut_with_trig_05(0), // 38

    deut_phi_pt_pos_T(0),          // 39
    deut_phi_pt_neg_T(0),          // 40

    prot_phi_pt_pos_T(0),          // 41
    prot_phi_pt_neg_T(0),          // 42
		
    deut_phi_pt_pos_sa0(0),        // 43
    deut_phi_pt_neg_sa0(0),        // 44
    deut_phi_pt_pos_sa1(0),        // 45
    deut_phi_pt_neg_sa1(0),        // 46
    deut_phi_pt_pos_sa2(0),        // 47
    deut_phi_pt_neg_sa2(0),        // 48
    deut_phi_pt_pos_sa3(0),        // 49
    deut_phi_pt_neg_sa3(0),        // 50
    deut_phi_pt_pos_sa4(0),        // 51
    deut_phi_pt_neg_sa4(0),        // 52
    deut_phi_pt_pos_sa5(0),        // 53
    deut_phi_pt_neg_sa5(0),        // 54
    deut_phi_pt_pos_sa6(0),        // 51
    deut_phi_pt_neg_sa6(0),        // 52
    deut_phi_pt_pos_sa7(0),        // 53
    deut_phi_pt_neg_sa7(0),        // 54
									   
    deut_phi_pt_pos_sb0(0),        // 55
    deut_phi_pt_neg_sb0(0),        // 56
    deut_phi_pt_pos_sb1(0),        // 57
    deut_phi_pt_neg_sb1(0),        // 58
    deut_phi_pt_pos_sb2(0),        // 59
    deut_phi_pt_neg_sb2(0),        // 60
    deut_phi_pt_pos_sb3(0),        // 61
    deut_phi_pt_neg_sb3(0),        // 62
    deut_phi_pt_pos_sb4(0),        // 63
    deut_phi_pt_neg_sb4(0),        // 64
    deut_phi_pt_pos_sb5(0),        // 65
    deut_phi_pt_neg_sb5(0),        // 66
    deut_phi_pt_pos_sb6(0),        // 51
    deut_phi_pt_neg_sb6(0),        // 52
    deut_phi_pt_pos_sb7(0),        // 53
    deut_phi_pt_neg_sb7(0),        // 54

    hadr_per_event(0),
    deut_per_event(0),             // 67
    prot_per_event(0),             // 68
    deut_per_event_with_trig(0),
    trig_05_per_event(0),          // 69
    trig_05_phi_pt_pos(0),         // 70
    trig_05_phi_pt_neg(0),         // 71
    
    tof_phi_eta_pos(0),            // 72
    tof_phi_eta_neg(0),            // 73
    tof_phi_eta_pos_deut(0),       // 74
    tof_phi_eta_neg_deut(0),       // 75

    deut_pt_compare_pos(0),        // 76
    deut_pt_compare_neg(0),        // 77
    tpc_sector_fraction(0),        // 78
    primary_vertex_z(0),           // 79
    primary_vertex_z_cut1(0),      // 80
    primary_vertex_z_cut2(0),      // 81
	
    deut_dphi_pt_pos_pos_05_T(0),  // 82
    deut_dphi_pt_pos_neg_05_T(0),  // 83
    deut_dphi_pt_neg_neg_05_T(0),  // 84
    
    prot_dphi_pt_pos_pos_05_T(0),  // 85
    prot_dphi_pt_pos_neg_05_T(0),  // 86
    prot_dphi_pt_neg_neg_05_T(0),  // 87

    deut_dphi_pt_pos_pos_05_sa0(0),// 88
    deut_dphi_pt_pos_neg_05_sa0(0),// 89
    deut_dphi_pt_neg_neg_05_sa0(0),// 90
    deut_dphi_pt_pos_pos_05_sa1(0),// 91
    deut_dphi_pt_pos_neg_05_sa1(0),// 92
    deut_dphi_pt_neg_neg_05_sa1(0),// 93
    deut_dphi_pt_pos_pos_05_sa2(0),// 94
    deut_dphi_pt_pos_neg_05_sa2(0),// 95
    deut_dphi_pt_neg_neg_05_sa2(0),// 96
    deut_dphi_pt_pos_pos_05_sa3(0),// 97
    deut_dphi_pt_pos_neg_05_sa3(0),// 98
    deut_dphi_pt_neg_neg_05_sa3(0),// 99
    deut_dphi_pt_pos_pos_05_sa4(0),//100
    deut_dphi_pt_pos_neg_05_sa4(0),//101
    deut_dphi_pt_neg_neg_05_sa4(0),//102
    deut_dphi_pt_pos_pos_05_sa5(0),//103
    deut_dphi_pt_pos_neg_05_sa5(0),//104
    deut_dphi_pt_neg_neg_05_sa5(0),//105
    deut_dphi_pt_pos_pos_05_sa6(0),//100
    deut_dphi_pt_pos_neg_05_sa6(0),//101
    deut_dphi_pt_neg_neg_05_sa6(0),//102
    deut_dphi_pt_pos_pos_05_sa7(0),//103
    deut_dphi_pt_pos_neg_05_sa7(0),//104
    deut_dphi_pt_neg_neg_05_sa7(0),//105
									   
    deut_dphi_pt_pos_pos_05_sb0(0),//106
    deut_dphi_pt_pos_neg_05_sb0(0),//107
    deut_dphi_pt_neg_neg_05_sb0(0),//108
    deut_dphi_pt_pos_pos_05_sb1(0),//109
    deut_dphi_pt_pos_neg_05_sb1(0),//110
    deut_dphi_pt_neg_neg_05_sb1(0),//111
    deut_dphi_pt_pos_pos_05_sb2(0),//112
    deut_dphi_pt_pos_neg_05_sb2(0),//113
    deut_dphi_pt_neg_neg_05_sb2(0),//114
    deut_dphi_pt_pos_pos_05_sb3(0),//115
    deut_dphi_pt_pos_neg_05_sb3(0),//116
    deut_dphi_pt_neg_neg_05_sb3(0),//117
    deut_dphi_pt_pos_pos_05_sb4(0),//118
    deut_dphi_pt_pos_neg_05_sb4(0),//119
    deut_dphi_pt_neg_neg_05_sb4(0),//120
    deut_dphi_pt_pos_pos_05_sb5(0),//121
    deut_dphi_pt_pos_neg_05_sb5(0),//122
    deut_dphi_pt_neg_neg_05_sb5(0),//123
    deut_dphi_pt_pos_pos_05_sb6(0),//118
    deut_dphi_pt_pos_neg_05_sb6(0),//119
    deut_dphi_pt_neg_neg_05_sb6(0),//120
    deut_dphi_pt_pos_pos_05_sb7(0),//121
    deut_dphi_pt_pos_neg_05_sb7(0),//122
    deut_dphi_pt_neg_neg_05_sb7(0),//123

    DCAxy_pos(0),                  //124
    DCAxy_neg(0),                  //125
    DCAz_pos(0),                   //126
    DCAz_neg(0),                   //127
    
    m2_pt_pos_fine(0),             //128
    m2_pt_neg_fine(0),             //129
    m2_pt_pos_TPC_fine(0),         //130
    m2_pt_neg_TPC_fine(0),         //131 
    m2_pt_pos_cut_T_fine(0),       //132
    m2_pt_neg_cut_T_fine(0),       //133

    m2_pt_pos_cut_sb0_fine(0),     //134
    m2_pt_neg_cut_sb0_fine(0),     //135
    m2_pt_pos_cut_sb1_fine(0),     //136
    m2_pt_neg_cut_sb1_fine(0),     //137
    m2_pt_pos_cut_sb2_fine(0),     //138
    m2_pt_neg_cut_sb2_fine(0),     //139
    m2_pt_pos_cut_sb3_fine(0),     //140
    m2_pt_neg_cut_sb3_fine(0),     //141
    m2_pt_pos_cut_sb4_fine(0),     //142
    m2_pt_neg_cut_sb4_fine(0),     //143
    m2_pt_pos_cut_sb5_fine(0),     //144
    m2_pt_neg_cut_sb5_fine(0),     //145
    m2_pt_pos_cut_sb6_fine(0),     //142
    m2_pt_neg_cut_sb6_fine(0),     //143
    m2_pt_pos_cut_sb7_fine(0),     //144
    m2_pt_neg_cut_sb7_fine(0),     //145
									   
    m2_pt_pos_cut_sa0_fine(0),     //146
    m2_pt_neg_cut_sa0_fine(0),     //147
    m2_pt_pos_cut_sa1_fine(0),     //148
    m2_pt_neg_cut_sa1_fine(0),     //149
    m2_pt_pos_cut_sa2_fine(0),     //150
    m2_pt_neg_cut_sa2_fine(0),     //151
    m2_pt_pos_cut_sa3_fine(0),     //152
    m2_pt_neg_cut_sa3_fine(0),     //153
    m2_pt_pos_cut_sa4_fine(0),     //154
    m2_pt_neg_cut_sa4_fine(0),     //155
    m2_pt_pos_cut_sa5_fine(0),     //156
    m2_pt_neg_cut_sa5_fine(0),     //157
    m2_pt_pos_cut_sa6_fine(0),     //154
    m2_pt_neg_cut_sa6_fine(0),     //155
    m2_pt_pos_cut_sa7_fine(0),     //156
    m2_pt_neg_cut_sa7_fine(0),     //157

    m2_pt_pos_TPC_prot_fine(0),    //158
    m2_pt_neg_TPC_prot_fine(0),    //159

    m2_pt_pos_cut_T_prot_fine(0),  //160
    m2_pt_neg_cut_T_prot_fine(0),   //161

    deut_deut_ptA_ptB_same(0),
    deut_deut_ptA_ptB_diff(0),
    deut_deut_dpt_dcharge(0),
    deut_deut_dphi_deta_same(0),
    deut_deut_dphi_deta_diff(0),
    deut_deut_inv_mass(0)									   
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

//    cout<<endl<<endl<<endl<<name<<endl<<endl<<endl;

//    int i = 0;
    run_mode = atoi(name);
    cout<<endl<<run_mode<<endl;

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



    deut_eta_pos_pt1               = new TH1F("deut_eta_pos_pt1",      "Deut eta, 1.00 < pt <1.35", 40, -1.0, 1.0); 
    deut_eta_neg_pt1               = new TH1F("deut_eta_neg_pt1", "Anti-Deut eta, 1.00 < pt <1.35", 40, -1.0, 1.0); 
    deut_eta_pos_pt2               = new TH1F("deut_eta_pos_pt2",      "Deut eta, 1.35 < pt <1.80", 40, -1.0, 1.0); 
    deut_eta_neg_pt2               = new TH1F("deut_eta_neg_pt2", "Anti-Deut eta, 1.35 < pt <1.80", 40, -1.0, 1.0);
    deut_eta_pos_pt3               = new TH1F("deut_eta_pos_pt3",      "Deut eta, 1.80 < pt <2.40", 40, -1.0, 1.0); 
    deut_eta_neg_pt3               = new TH1F("deut_eta_neg_pt3", "Anti-Deut eta, 1.80 < pt <2.40", 40, -1.0, 1.0);
    deut_eta_pos_pt4               = new TH1F("deut_eta_pos_pt4",      "Deut eta, 2.40 < pt <3.00", 40, -1.0, 1.0); 
    deut_eta_neg_pt4               = new TH1F("deut_eta_neg_pt4", "Anti-Deut eta, 2.40 < pt <3.00", 40, -1.0, 1.0);
    deut_eta_pos_pt5               = new TH1F("deut_eta_pos_pt5",      "Deut eta, 3.00 < pt <4.00", 40, -1.0, 1.0); 
    deut_eta_neg_pt5               = new TH1F("deut_eta_neg_pt5", "Anti-Deut eta, 3.00 < pt <4.00", 40, -1.0, 1.0);
    
    hadr_eta_pos               = new TH1F("hadr_eta_pos", "Trig Pos Hadr eta,  pt > 5.0", 40, -1.0, 1.0); 
    hadr_eta_neg               = new TH1F("hadr_eta_neg", "Trig Neg Hadr eta,  pt > 5.0", 40, -1.0, 1.0); 

    deut_hadr_dphi36_deta20_pt1    = new TH2F("deut_hadr_dphi36_deta20_pt1", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:1.0-1.35) deta-dphi", 36,   lower,   upper, 20, -1.8, 1.8);
    deut_hadr_dphi36_deta20_pt2    = new TH2F("deut_hadr_dphi36_deta20_pt2", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:1.35-1.8) deta-dphi", 36,   lower,   upper, 20, -1.8, 1.8);
    deut_hadr_dphi36_deta20_pt3    = new TH2F("deut_hadr_dphi36_deta20_pt3", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:1.8-2.4) deta-dphi",  36,   lower,   upper, 20, -1.8, 1.8);
    deut_hadr_dphi36_deta20_pt4    = new TH2F("deut_hadr_dphi36_deta20_pt4", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:2.4-3.0) deta-dphi",  36,   lower,   upper, 20, -1.8, 1.8);
    deut_hadr_dphi36_deta20_pt5    = new TH2F("deut_hadr_dphi36_deta20_pt5", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:3.0-4.0) deta-dphi",  36,   lower,   upper, 20, -1.8, 1.8);

    deut_hadr_dphi36_deta21_pt1    = new TH2F("deut_hadr_dphi36_deta21_pt1", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:1.0-1.35) deta-dphi", 36,   lower,   upper, 21, -1.8, 1.8);
    deut_hadr_dphi36_deta21_pt2    = new TH2F("deut_hadr_dphi36_deta21_pt2", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:1.35-1.8) deta-dphi", 36,   lower,   upper, 21, -1.8, 1.8);
    deut_hadr_dphi36_deta21_pt3    = new TH2F("deut_hadr_dphi36_deta21_pt3", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:1.8-2.4) deta-dphi",  36,   lower,   upper, 21, -1.8, 1.8);
    deut_hadr_dphi36_deta21_pt4    = new TH2F("deut_hadr_dphi36_deta21_pt4", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:2.4-3.0) deta-dphi",  36,   lower,   upper, 21, -1.8, 1.8);
    deut_hadr_dphi36_deta21_pt5    = new TH2F("deut_hadr_dphi36_deta21_pt5", "Trig Hadr (pt>5.0) to (anti-)Deut.(pt:3.0-4.0) deta-dphi",  36,   lower,   upper, 21, -1.8, 1.8);

    
    multiplicity_w_deut_1      = new TH1I("multiplicity_w_deut_1",      "multiplicity_w_deut_1",       150,     0,   150);                                //  n1
    multiplicity_w_deut_2      = new TH1I("multiplicity_w_deut_2",      "multiplicity_w_deut_2",       150,     0,   150);                                //  n2
    multiplicity_w_deut_3      = new TH1I("multiplicity_w_deut_3",      "multiplicity_w_deut_3",       150,     0,   150);                                //  n3
    multiplicity_w_deut_4      = new TH1I("multiplicity_w_deut_4",      "multiplicity_w_deut_4",       150,     0,   150);                                //  n4
    multiplicity_w_deut_5      = new TH1I("multiplicity_w_deut_5",      "multiplicity_w_deut_5",       150,     0,   150);                                //  n5

    multiplicity_w_trig_and_deut_1 = new TH1I("multiplicity_w_trig_and_deut_1", "multiplicity_w_trig_and_deut_1", 150,     0,   150);                     //  n6
    multiplicity_w_trig_and_deut_2 = new TH1I("multiplicity_w_trig_and_deut_2", "multiplicity_w_trig_and_deut_2", 150,     0,   150);                     //  n7
    multiplicity_w_trig_and_deut_3 = new TH1I("multiplicity_w_trig_and_deut_3", "multiplicity_w_trig_and_deut_3", 150,     0,   150);                     //  n8
    multiplicity_w_trig_and_deut_4 = new TH1I("multiplicity_w_trig_and_deut_4", "multiplicity_w_trig_and_deut_4", 150,     0,   150);                     //  n9
    multiplicity_w_trig_and_deut_5 = new TH1I("multiplicity_w_trig_and_deut_5", "multiplicity_w_trig_and_deut_5", 150,     0,   150);                     //  n10
    
    
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

    m2_pt_pos_cut_sa0          = new TH2F("m2_pt_pos_cut_sa0",          "m2_pt_pos_cut_sa0",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa0          = new TH2F("m2_pt_neg_cut_sa0",          "m2_pt_neg_cut_sa0",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa1          = new TH2F("m2_pt_pos_cut_sa1",          "m2_pt_pos_cut_sa1",            18, coarse_binning,    2400,    -1.0,     7.0);   // 15
    m2_pt_neg_cut_sa1          = new TH2F("m2_pt_neg_cut_sa1",          "m2_pt_neg_cut_sa1",            18, coarse_binning,    2400,    -1.0,     7.0);   // 16
    m2_pt_pos_cut_sa2          = new TH2F("m2_pt_pos_cut_sa2",          "m2_pt_pos_cut_sa2",            18, coarse_binning,    2400,    -1.0,     7.0);   // 17
    m2_pt_neg_cut_sa2          = new TH2F("m2_pt_neg_cut_sa2",          "m2_pt_neg_cut_sa2",            18, coarse_binning,    2400,    -1.0,     7.0);   // 18
    m2_pt_pos_cut_sa3          = new TH2F("m2_pt_pos_cut_sa3",          "m2_pt_pos_cut_sa3",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa3          = new TH2F("m2_pt_neg_cut_sa3",          "m2_pt_neg_cut_sa3",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa4          = new TH2F("m2_pt_pos_cut_sa4",          "m2_pt_pos_cut_sa4",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa4          = new TH2F("m2_pt_neg_cut_sa4",          "m2_pt_neg_cut_sa4",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa5          = new TH2F("m2_pt_pos_cut_sa5",          "m2_pt_pos_cut_sa5",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa5          = new TH2F("m2_pt_neg_cut_sa5",          "m2_pt_neg_cut_sa5",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa6          = new TH2F("m2_pt_pos_cut_sa6",          "m2_pt_pos_cut_sa6",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa6          = new TH2F("m2_pt_neg_cut_sa6",          "m2_pt_neg_cut_sa6",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa7          = new TH2F("m2_pt_pos_cut_sa7",          "m2_pt_pos_cut_sa7",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa7          = new TH2F("m2_pt_neg_cut_sa7",          "m2_pt_neg_cut_sa7",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    
    m2_pt_pos_cut_sb0          = new TH2F("m2_pt_pos_cut_sb0",          "m2_pt_pos_cut_sb0",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb0          = new TH2F("m2_pt_neg_cut_sb0",          "m2_pt_neg_cut_sb0",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb1          = new TH2F("m2_pt_pos_cut_sb1",          "m2_pt_pos_cut_sb1",            18, coarse_binning,    2400,    -1.0,     7.0);   // 15
    m2_pt_neg_cut_sb1          = new TH2F("m2_pt_neg_cut_sb1",          "m2_pt_neg_cut_sb1",            18, coarse_binning,    2400,    -1.0,     7.0);   // 16
    m2_pt_pos_cut_sb2          = new TH2F("m2_pt_pos_cut_sb2",          "m2_pt_pos_cut_sb2",            18, coarse_binning,    2400,    -1.0,     7.0);   // 17
    m2_pt_neg_cut_sb2          = new TH2F("m2_pt_neg_cut_sb2",          "m2_pt_neg_cut_sb2",            18, coarse_binning,    2400,    -1.0,     7.0);   // 18
    m2_pt_pos_cut_sb3          = new TH2F("m2_pt_pos_cut_sb3",          "m2_pt_pos_cut_sb3",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb3          = new TH2F("m2_pt_neg_cut_sb3",          "m2_pt_neg_cut_sb3",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb4          = new TH2F("m2_pt_pos_cut_sb4",          "m2_pt_pos_cut_sb4",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb4          = new TH2F("m2_pt_neg_cut_sb4",          "m2_pt_neg_cut_sb4",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb5          = new TH2F("m2_pt_pos_cut_sb5",          "m2_pt_pos_cut_sb5",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb5          = new TH2F("m2_pt_neg_cut_sb5",          "m2_pt_neg_cut_sb5",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb6          = new TH2F("m2_pt_pos_cut_sb6",          "m2_pt_pos_cut_sb6",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb6          = new TH2F("m2_pt_neg_cut_sb6",          "m2_pt_neg_cut_sb6",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb7          = new TH2F("m2_pt_pos_cut_sb7",          "m2_pt_pos_cut_sb7",            18, coarse_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb7          = new TH2F("m2_pt_neg_cut_sb7",          "m2_pt_neg_cut_sb7",            18, coarse_binning,    2400,    -1.0,     7.0);   //


    m2_pt_pos_cut_T_prot       = new TH2F("m2_pt_pos_cut_T_prot",       "m2_pt_pos_cut_T_prot",         18,coarse_binning2,    2400,    -1.0,     7.0);   // 19
    m2_pt_neg_cut_T_prot       = new TH2F("m2_pt_neg_cut_T_prot",       "m2_pt_neg_cut_T_prot",         18,coarse_binning2,    2400,    -1.0,     7.0);   // 20
    
    m2_pt_pos_cut_with_trig_05 = new TH2F("m2_pt_pos_cut_with_trig_05", "m2_pt_pos_cut_with_trig_05",   18, coarse_binning,    2400,    -1.0,     7.0);   // 21
    m2_pt_neg_cut_with_trig_05 = new TH2F("m2_pt_neg_cut_with_trig_05", "m2_pt_neg_cut_with_trig_05",   18, coarse_binning,    2400,    -1.0,     7.0);   // 22
    
    deut_phi_pt_pos_T          = new TH2F("deut_phi_pt_pos_T",          "deut_phi_pt_pos_T",            18, coarse_binning,     288,   lower,   upper);   // 23
    deut_phi_pt_neg_T          = new TH2F("deut_phi_pt_neg_T",          "deut_phi_pt_neg_T",            18, coarse_binning,     288,   lower,   upper);   // 24
    prot_phi_pt_pos_T          = new TH2F("prot_phi_pt_pos_T",          "prot_phi_pt_pos_T",            18,coarse_binning2,     288,   lower,   upper);   // 25
    prot_phi_pt_neg_T          = new TH2F("prot_phi_pt_neg_T",          "prot_phi_pt_neg_T",            18,coarse_binning2,     288,   lower,   upper);   // 26

    deut_phi_pt_pos_sa0        = new TH2F("deut_phi_pt_pos_sa0",        "deut_phi_pt_pos_sa0",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sa0        = new TH2F("deut_phi_pt_neg_sa0",        "deut_phi_pt_neg_sa0",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sa1        = new TH2F("deut_phi_pt_pos_sa1",        "deut_phi_pt_pos_sa1",          18, coarse_binning,     288,   lower,   upper);   // 31
    deut_phi_pt_neg_sa1        = new TH2F("deut_phi_pt_neg_sa1",        "deut_phi_pt_neg_sa1",          18, coarse_binning,     288,   lower,   upper);   // 32
    deut_phi_pt_pos_sa2        = new TH2F("deut_phi_pt_pos_sa2",        "deut_phi_pt_pos_sa2",          18, coarse_binning,     288,   lower,   upper);   // 33
    deut_phi_pt_neg_sa2        = new TH2F("deut_phi_pt_neg_sa2",        "deut_phi_pt_neg_sa2",          18, coarse_binning,     288,   lower,   upper);   // 34
    deut_phi_pt_pos_sa3        = new TH2F("deut_phi_pt_pos_sa3",        "deut_phi_pt_pos_sa3",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sa3        = new TH2F("deut_phi_pt_neg_sa3",        "deut_phi_pt_neg_sa3",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sa4        = new TH2F("deut_phi_pt_pos_sa4",        "deut_phi_pt_pos_sa4",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sa4        = new TH2F("deut_phi_pt_neg_sa4",        "deut_phi_pt_neg_sa4",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sa5        = new TH2F("deut_phi_pt_pos_sa5",        "deut_phi_pt_pos_sa5",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sa5        = new TH2F("deut_phi_pt_neg_sa5",        "deut_phi_pt_neg_sa5",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sa6        = new TH2F("deut_phi_pt_pos_sa6",        "deut_phi_pt_pos_sa6",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sa6        = new TH2F("deut_phi_pt_neg_sa6",        "deut_phi_pt_neg_sa6",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sa7        = new TH2F("deut_phi_pt_pos_sa7",        "deut_phi_pt_pos_sa7",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sa7        = new TH2F("deut_phi_pt_neg_sa7",        "deut_phi_pt_neg_sa7",          18, coarse_binning,     288,   lower,   upper);   //
    
    deut_phi_pt_pos_sb0        = new TH2F("deut_phi_pt_pos_sb0",        "deut_phi_pt_pos_sb0",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sb0        = new TH2F("deut_phi_pt_neg_sb0",        "deut_phi_pt_neg_sb0",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sb1        = new TH2F("deut_phi_pt_pos_sb1",        "deut_phi_pt_pos_sb1",          18, coarse_binning,     288,   lower,   upper);   // 31
    deut_phi_pt_neg_sb1        = new TH2F("deut_phi_pt_neg_sb1",        "deut_phi_pt_neg_sb1",          18, coarse_binning,     288,   lower,   upper);   // 32
    deut_phi_pt_pos_sb2        = new TH2F("deut_phi_pt_pos_sb2",        "deut_phi_pt_pos_sb2",          18, coarse_binning,     288,   lower,   upper);   // 33
    deut_phi_pt_neg_sb2        = new TH2F("deut_phi_pt_neg_sb2",        "deut_phi_pt_neg_sb2",          18, coarse_binning,     288,   lower,   upper);   // 34
    deut_phi_pt_pos_sb3        = new TH2F("deut_phi_pt_pos_sb3",        "deut_phi_pt_pos_sb3",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sb3        = new TH2F("deut_phi_pt_neg_sb3",        "deut_phi_pt_neg_sb3",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sb4        = new TH2F("deut_phi_pt_pos_sb4",        "deut_phi_pt_pos_sb4",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sb4        = new TH2F("deut_phi_pt_neg_sb4",        "deut_phi_pt_neg_sb4",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sb5        = new TH2F("deut_phi_pt_pos_sb5",        "deut_phi_pt_pos_sb5",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sb5        = new TH2F("deut_phi_pt_neg_sb5",        "deut_phi_pt_neg_sb5",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sb6        = new TH2F("deut_phi_pt_pos_sb6",        "deut_phi_pt_pos_sb6",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sb6        = new TH2F("deut_phi_pt_neg_sb6",        "deut_phi_pt_neg_sb6",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_pos_sb7        = new TH2F("deut_phi_pt_pos_sb7",        "deut_phi_pt_pos_sb7",          18, coarse_binning,     288,   lower,   upper);   //
    deut_phi_pt_neg_sb7        = new TH2F("deut_phi_pt_neg_sb7",        "deut_phi_pt_neg_sb7",          18, coarse_binning,     288,   lower,   upper);   //

    hadr_per_event             = new TH1I("hadr_per_event",             "hadr_per_event",              150,      0,    150);                              // 35
    deut_per_event             = new TH1I("deut_per_event",             "deut_per_event",                4,      0,      4);                              // 35
    prot_per_event             = new TH1I("prot_per_event",             "prot_per_event",               20,      0,     20);                              // 36
    deut_per_event_with_trig   = new TH1I("deut_per_event_with_trig",   "deut_per_event_with_trig",      4,      0,      4);
    trig_05_per_event          = new TH1I("trig_05_per_event",          "trig_05_per_event",            10,      0,     10);                              // 37
    trig_05_phi_pt_pos         = new TH2F("trig_05_phi_pt_pos",         "trig_05_phi_pt_pos",          170,    3.0,   20.0,     288,   lower,   upper);   // 38
    trig_05_phi_pt_neg         = new TH2F("trig_05_phi_pt_neg",         "trig_05_phi_pt_neg",          170,    3.0,   20.0,     288,   lower,   upper);   // 39

    tof_phi_eta_pos            = new TH2F("tof_phi_eta_pos",            "tof_phi_eta_pos",             190,  lower,  upper,      60,   -0.92,    0.92);   // 40
    tof_phi_eta_neg            = new TH2F("tof_phi_eta_neg",            "tof_phi_eta_neg",             190,  lower,  upper,      60,   -0.92,    0.92);   // 41
    tof_phi_eta_pos_deut       = new TH2F("tof_phi_eta_pos_deut",       "tof_phi_eta_pos_deut",        190,  lower,  upper,      60,   -0.92,    0.92);   // 42
    tof_phi_eta_neg_deut       = new TH2F("tof_phi_eta_neg_deut",       "tof_phi_eta_neg_deut",        190,  lower,  upper,      60,   -0.92,    0.92);   // 43

    deut_pt_compare_pos        = new TH1F("deut_pt_compare_pos",        "deut_pt_compare_pos",         100,   -0.2,    0.2);                              // 44
    deut_pt_compare_neg        = new TH1F("deut_pt_compare_neg",        "deut_pt_compare_neg",         100,   -0.2,    0.2);                              // 45
    tpc_sector_fraction        = new TH1F("tpc_sector_fraction",        "tpc_sector_fraction",         100,    0.0,    1.0);                              // 46
    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            400,  -20.0,   20.0);                              // 47
    primary_vertex_z_cut1      = new TH1F("primary_vertex_z_cut1",      "primary_vertex_z_cut1",       400,  -20.0,   20.0);                              // 48
    primary_vertex_z_cut2      = new TH1F("primary_vertex_z_cut2",      "primary_vertex_z_cut2",       400,  -20.0,   20.0);                              // 49

    deut_dphi_pt_pos_pos_05_T  = new TH2F("deut_dphi_pt_pos_pos_05_T",  "deut_dphi_pt_pos_pos_05_T",          18, coarse_binning,     288,   lower,   upper);   // 50
    deut_dphi_pt_pos_neg_05_T  = new TH2F("deut_dphi_pt_pos_neg_05_T",  "deut_dphi_pt_pos_neg_05_T",          18, coarse_binning,     288,   lower,   upper);   // 51
    deut_dphi_pt_neg_neg_05_T  = new TH2F("deut_dphi_pt_neg_neg_05_T",  "deut_dphi_pt_neg_neg_05_T",          18, coarse_binning,     288,   lower,   upper);   // 52

    prot_dphi_pt_pos_pos_05_T  = new TH2F("prot_dphi_pt_pos_pos_05_T",  "prot_dphi_pt_pos_pos_05_T",          18,coarse_binning2,     288,   lower,   upper);   // 53
    prot_dphi_pt_pos_neg_05_T  = new TH2F("prot_dphi_pt_pos_neg_05_T",  "prot_dphi_pt_pos_neg_05_T",          18,coarse_binning2,     288,   lower,   upper);   // 54
    prot_dphi_pt_neg_neg_05_T  = new TH2F("prot_dphi_pt_neg_neg_05_T",  "prot_dphi_pt_neg_neg_05_T",          18,coarse_binning2,     288,   lower,   upper);   // 55
    
    deut_dphi_pt_pos_pos_05_sa0  = new TH2F("deut_dphi_pt_pos_pos_05_sa0",  "deut_dphi_pt_pos_pos_05_sa0",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sa0  = new TH2F("deut_dphi_pt_pos_neg_05_sa0",  "deut_dphi_pt_pos_neg_05_sa0",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sa0  = new TH2F("deut_dphi_pt_neg_neg_05_sa0",  "deut_dphi_pt_neg_neg_05_sa0",    18, coarse_binning,     288,   lower,   upper);   //
    
    deut_dphi_pt_pos_pos_05_sa1  = new TH2F("deut_dphi_pt_pos_pos_05_sa1",  "deut_dphi_pt_pos_pos_05_sa1",    18, coarse_binning,     288,   lower,   upper);   // 62
    deut_dphi_pt_pos_neg_05_sa1  = new TH2F("deut_dphi_pt_pos_neg_05_sa1",  "deut_dphi_pt_pos_neg_05_sa1",    18, coarse_binning,     288,   lower,   upper);   // 63
    deut_dphi_pt_neg_neg_05_sa1  = new TH2F("deut_dphi_pt_neg_neg_05_sa1",  "deut_dphi_pt_neg_neg_05_sa1",    18, coarse_binning,     288,   lower,   upper);   // 64
    
    deut_dphi_pt_pos_pos_05_sa2  = new TH2F("deut_dphi_pt_pos_pos_05_sa2",  "deut_dphi_pt_pos_pos_05_sa2",    18, coarse_binning,     288,   lower,   upper);   // 65
    deut_dphi_pt_pos_neg_05_sa2  = new TH2F("deut_dphi_pt_pos_neg_05_sa2",  "deut_dphi_pt_pos_neg_05_sa2",    18, coarse_binning,     288,   lower,   upper);   // 66
    deut_dphi_pt_neg_neg_05_sa2  = new TH2F("deut_dphi_pt_neg_neg_05_sa2",  "deut_dphi_pt_neg_neg_05_sa2",    18, coarse_binning,     288,   lower,   upper);   // 67

    deut_dphi_pt_pos_pos_05_sa3  = new TH2F("deut_dphi_pt_pos_pos_05_sa3",  "deut_dphi_pt_pos_pos_05_sa3",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sa3  = new TH2F("deut_dphi_pt_pos_neg_05_sa3",  "deut_dphi_pt_pos_neg_05_sa3",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sa3  = new TH2F("deut_dphi_pt_neg_neg_05_sa3",  "deut_dphi_pt_neg_neg_05_sa3",    18, coarse_binning,     288,   lower,   upper);   //
    
    deut_dphi_pt_pos_pos_05_sa4  = new TH2F("deut_dphi_pt_pos_pos_05_sa4",  "deut_dphi_pt_pos_pos_05_sa4",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sa4  = new TH2F("deut_dphi_pt_pos_neg_05_sa4",  "deut_dphi_pt_pos_neg_05_sa4",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sa4  = new TH2F("deut_dphi_pt_neg_neg_05_sa4",  "deut_dphi_pt_neg_neg_05_sa4",    18, coarse_binning,     288,   lower,   upper);   //

    deut_dphi_pt_pos_pos_05_sa5  = new TH2F("deut_dphi_pt_pos_pos_05_sa5",  "deut_dphi_pt_pos_pos_05_sa5",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sa5  = new TH2F("deut_dphi_pt_pos_neg_05_sa5",  "deut_dphi_pt_pos_neg_05_sa5",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sa5  = new TH2F("deut_dphi_pt_neg_neg_05_sa5",  "deut_dphi_pt_neg_neg_05_sa5",    18, coarse_binning,     288,   lower,   upper);   //

    deut_dphi_pt_pos_pos_05_sa6  = new TH2F("deut_dphi_pt_pos_pos_05_sa6",  "deut_dphi_pt_pos_pos_05_sa6",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sa6  = new TH2F("deut_dphi_pt_pos_neg_05_sa6",  "deut_dphi_pt_pos_neg_05_sa6",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sa6  = new TH2F("deut_dphi_pt_neg_neg_05_sa6",  "deut_dphi_pt_neg_neg_05_sa6",    18, coarse_binning,     288,   lower,   upper);   //

    deut_dphi_pt_pos_pos_05_sa7  = new TH2F("deut_dphi_pt_pos_pos_05_sa7",  "deut_dphi_pt_pos_pos_05_sa7",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sa7  = new TH2F("deut_dphi_pt_pos_neg_05_sa7",  "deut_dphi_pt_pos_neg_05_sa7",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sa7  = new TH2F("deut_dphi_pt_neg_neg_05_sa7",  "deut_dphi_pt_neg_neg_05_sa7",    18, coarse_binning,     288,   lower,   upper);   //


    
    deut_dphi_pt_pos_pos_05_sb0  = new TH2F("deut_dphi_pt_pos_pos_05_sb0",  "deut_dphi_pt_pos_pos_05_sb0",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sb0  = new TH2F("deut_dphi_pt_pos_neg_05_sb0",  "deut_dphi_pt_pos_neg_05_sb0",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sb0  = new TH2F("deut_dphi_pt_neg_neg_05_sb0",  "deut_dphi_pt_neg_neg_05_sb0",    18, coarse_binning,     288,   lower,   upper);   //
    
    deut_dphi_pt_pos_pos_05_sb1  = new TH2F("deut_dphi_pt_pos_pos_05_sb1",  "deut_dphi_pt_pos_pos_05_sb1",    18, coarse_binning,     288,   lower,   upper);   // 62
    deut_dphi_pt_pos_neg_05_sb1  = new TH2F("deut_dphi_pt_pos_neg_05_sb1",  "deut_dphi_pt_pos_neg_05_sb1",    18, coarse_binning,     288,   lower,   upper);   // 63
    deut_dphi_pt_neg_neg_05_sb1  = new TH2F("deut_dphi_pt_neg_neg_05_sb1",  "deut_dphi_pt_neg_neg_05_sb1",    18, coarse_binning,     288,   lower,   upper);   // 64
    
    deut_dphi_pt_pos_pos_05_sb2  = new TH2F("deut_dphi_pt_pos_pos_05_sb2",  "deut_dphi_pt_pos_pos_05_sb2",    18, coarse_binning,     288,   lower,   upper);   // 65
    deut_dphi_pt_pos_neg_05_sb2  = new TH2F("deut_dphi_pt_pos_neg_05_sb2",  "deut_dphi_pt_pos_neg_05_sb2",    18, coarse_binning,     288,   lower,   upper);   // 66
    deut_dphi_pt_neg_neg_05_sb2  = new TH2F("deut_dphi_pt_neg_neg_05_sb2",  "deut_dphi_pt_neg_neg_05_sb2",    18, coarse_binning,     288,   lower,   upper);   // 67

    deut_dphi_pt_pos_pos_05_sb3  = new TH2F("deut_dphi_pt_pos_pos_05_sb3",  "deut_dphi_pt_pos_pos_05_sb3",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sb3  = new TH2F("deut_dphi_pt_pos_neg_05_sb3",  "deut_dphi_pt_pos_neg_05_sb3",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sb3  = new TH2F("deut_dphi_pt_neg_neg_05_sb3",  "deut_dphi_pt_neg_neg_05_sb3",    18, coarse_binning,     288,   lower,   upper);   //
    
    deut_dphi_pt_pos_pos_05_sb4  = new TH2F("deut_dphi_pt_pos_pos_05_sb4",  "deut_dphi_pt_pos_pos_05_sb4",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sb4  = new TH2F("deut_dphi_pt_pos_neg_05_sb4",  "deut_dphi_pt_pos_neg_05_sb4",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sb4  = new TH2F("deut_dphi_pt_neg_neg_05_sb4",  "deut_dphi_pt_neg_neg_05_sb4",    18, coarse_binning,     288,   lower,   upper);   //

    deut_dphi_pt_pos_pos_05_sb5  = new TH2F("deut_dphi_pt_pos_pos_05_sb5",  "deut_dphi_pt_pos_pos_05_sb5",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sb5  = new TH2F("deut_dphi_pt_pos_neg_05_sb5",  "deut_dphi_pt_pos_neg_05_sb5",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sb5  = new TH2F("deut_dphi_pt_neg_neg_05_sb5",  "deut_dphi_pt_neg_neg_05_sb5",    18, coarse_binning,     288,   lower,   upper);   //

    deut_dphi_pt_pos_pos_05_sb6  = new TH2F("deut_dphi_pt_pos_pos_05_sb6",  "deut_dphi_pt_pos_pos_05_sb6",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sb6  = new TH2F("deut_dphi_pt_pos_neg_05_sb6",  "deut_dphi_pt_pos_neg_05_sb6",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sb6  = new TH2F("deut_dphi_pt_neg_neg_05_sb6",  "deut_dphi_pt_neg_neg_05_sb6",    18, coarse_binning,     288,   lower,   upper);   //

    deut_dphi_pt_pos_pos_05_sb7  = new TH2F("deut_dphi_pt_pos_pos_05_sb7",  "deut_dphi_pt_pos_pos_05_sb7",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_pos_neg_05_sb7  = new TH2F("deut_dphi_pt_pos_neg_05_sb7",  "deut_dphi_pt_pos_neg_05_sb7",    18, coarse_binning,     288,   lower,   upper);   //
    deut_dphi_pt_neg_neg_05_sb7  = new TH2F("deut_dphi_pt_neg_neg_05_sb7",  "deut_dphi_pt_neg_neg_05_sb7",    18, coarse_binning,     288,   lower,   upper);   //
    
    
    DCAxy_pos                  = new TH1F("DCAxy_pos",                  "DCAxy_pos",                   100,   -1.2,    1.2);                              // 68
    DCAxy_neg                  = new TH1F("DCAxy_neg",                  "DCAxy_neg",                   100,   -1.2,    1.2);                              // 69
    DCAz_pos                   = new TH1F("DCAz_pos",                   "DCAz_pos",                    100,   -1.2,    1.2);                              // 70
    DCAz_neg                   = new TH1F("DCAz_neg",                   "DCAz_neg",                    100,   -1.2,    1.2);                              // 71

    Double_t fine_binning[2001];

    Float_t moving_marker = 0.10;
    for(int i=0; i<1602; i++)
    {
	fine_binning[i] = moving_marker;
	moving_marker = moving_marker + fine_binning[i] * 0.005;
    }
    // display plots

    m2_pt_pos_fine             = new TH2F("m2_pt_pos_fine",             "m2_pt_pos_fine",              800,   fine_binning,    2400,    -1.0,     7.0);   // 72
    m2_pt_neg_fine             = new TH2F("m2_pt_neg_fine",             "m2_pt_neg_fine",              800,   fine_binning,    2400,    -1.0,     7.0);   // 73
    m2_pt_pos_TPC_fine         = new TH2F("m2_pt_pos_TPC_fine",         "m2_pt_pos_TPC_fine",          800,   fine_binning,    2400,    -1.0,     7.0);   // 74
    m2_pt_neg_TPC_fine         = new TH2F("m2_pt_neg_TPC_fine",         "m2_pt_neg_TPC_fine",          800,   fine_binning,    2400,    -1.0,     7.0);   // 75
    m2_pt_pos_cut_T_fine       = new TH2F("m2_pt_pos_cut_T_fine",       "m2_pt_pos_cut_T_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 76
    m2_pt_neg_cut_T_fine       = new TH2F("m2_pt_neg_cut_T_fine",       "m2_pt_neg_cut_T_fine",        800,   fine_binning,    2400,    -1.0,     7.0);   // 77

    m2_pt_pos_cut_sa0_fine     = new TH2F("m2_pt_pos_cut_sa0_fine",     "m2_pt_pos_cut_sa0_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa0_fine     = new TH2F("m2_pt_neg_cut_sa0_fine",     "m2_pt_neg_cut_sa0_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //    
    m2_pt_pos_cut_sa1_fine     = new TH2F("m2_pt_pos_cut_sa1_fine",     "m2_pt_pos_cut_sa1_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 82
    m2_pt_neg_cut_sa1_fine     = new TH2F("m2_pt_neg_cut_sa1_fine",     "m2_pt_neg_cut_sa1_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 83
    m2_pt_pos_cut_sa2_fine     = new TH2F("m2_pt_pos_cut_sa2_fine",     "m2_pt_pos_cut_sa2_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 84
    m2_pt_neg_cut_sa2_fine     = new TH2F("m2_pt_neg_cut_sa2_fine",     "m2_pt_neg_cut_sa2_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 85
    m2_pt_pos_cut_sa3_fine     = new TH2F("m2_pt_pos_cut_sa3_fine",     "m2_pt_pos_cut_sa3_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa3_fine     = new TH2F("m2_pt_neg_cut_sa3_fine",     "m2_pt_neg_cut_sa3_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa4_fine     = new TH2F("m2_pt_pos_cut_sa4_fine",     "m2_pt_pos_cut_sa4_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa4_fine     = new TH2F("m2_pt_neg_cut_sa4_fine",     "m2_pt_neg_cut_sa4_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa5_fine     = new TH2F("m2_pt_pos_cut_sa5_fine",     "m2_pt_pos_cut_sa5_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa5_fine     = new TH2F("m2_pt_neg_cut_sa5_fine",     "m2_pt_neg_cut_sa5_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa6_fine     = new TH2F("m2_pt_pos_cut_sa6_fine",     "m2_pt_pos_cut_sa6_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa6_fine     = new TH2F("m2_pt_neg_cut_sa6_fine",     "m2_pt_neg_cut_sa6_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sa7_fine     = new TH2F("m2_pt_pos_cut_sa7_fine",     "m2_pt_pos_cut_sa7_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sa7_fine     = new TH2F("m2_pt_neg_cut_sa7_fine",     "m2_pt_neg_cut_sa7_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //

    
    m2_pt_pos_cut_sb0_fine     = new TH2F("m2_pt_pos_cut_sb0_fine",     "m2_pt_pos_cut_sb0_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb0_fine     = new TH2F("m2_pt_neg_cut_sb0_fine",     "m2_pt_neg_cut_sb0_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //    
    m2_pt_pos_cut_sb1_fine     = new TH2F("m2_pt_pos_cut_sb1_fine",     "m2_pt_pos_cut_sb1_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 82
    m2_pt_neg_cut_sb1_fine     = new TH2F("m2_pt_neg_cut_sb1_fine",     "m2_pt_neg_cut_sb1_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 83
    m2_pt_pos_cut_sb2_fine     = new TH2F("m2_pt_pos_cut_sb2_fine",     "m2_pt_pos_cut_sb2_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 84
    m2_pt_neg_cut_sb2_fine     = new TH2F("m2_pt_neg_cut_sb2_fine",     "m2_pt_neg_cut_sb2_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   // 85
    m2_pt_pos_cut_sb3_fine     = new TH2F("m2_pt_pos_cut_sb3_fine",     "m2_pt_pos_cut_sb3_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb3_fine     = new TH2F("m2_pt_neg_cut_sb3_fine",     "m2_pt_neg_cut_sb3_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb4_fine     = new TH2F("m2_pt_pos_cut_sb4_fine",     "m2_pt_pos_cut_sb4_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb4_fine     = new TH2F("m2_pt_neg_cut_sb4_fine",     "m2_pt_neg_cut_sb4_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb5_fine     = new TH2F("m2_pt_pos_cut_sb5_fine",     "m2_pt_pos_cut_sb5_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb5_fine     = new TH2F("m2_pt_neg_cut_sb5_fine",     "m2_pt_neg_cut_sb5_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb6_fine     = new TH2F("m2_pt_pos_cut_sb6_fine",     "m2_pt_pos_cut_sb6_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb6_fine     = new TH2F("m2_pt_neg_cut_sb6_fine",     "m2_pt_neg_cut_sb6_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_pos_cut_sb7_fine     = new TH2F("m2_pt_pos_cut_sb7_fine",     "m2_pt_pos_cut_sb7_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //
    m2_pt_neg_cut_sb7_fine     = new TH2F("m2_pt_neg_cut_sb7_fine",     "m2_pt_neg_cut_sb7_fine",      800,   fine_binning,    2400,    -1.0,     7.0);   //

    
    m2_pt_pos_TPC_prot_fine    = new TH2F("m2_pt_pos_TPC_prot_fine",    "m2_pt_pos_TPC_prot_fine",     800,   fine_binning,    2400,    -1.0,     7.0);   // 86
    m2_pt_neg_TPC_prot_fine    = new TH2F("m2_pt_neg_TPC_prot_fine",    "m2_pt_neg_TPC_prot_fine",     800,   fine_binning,    2400,    -1.0,     7.0);   // 87

    m2_pt_pos_cut_T_prot_fine  = new TH2F("m2_pt_pos_cut_T_prot_fine",  "m2_pt_pos_cut_T_prot_fine",   800,   fine_binning,    2400,    -1.0,     7.0);   // 88
    m2_pt_neg_cut_T_prot_fine  = new TH2F("m2_pt_neg_cut_T_prot_fine",  "m2_pt_neg_cut_T_prot_fine",   800,   fine_binning,    2400,    -1.0,     7.0);   // 89


    deut_deut_dphi_deta_same   = new TH2F("deut_deut_dphi_deta_same",   "deut_deut_dphi_deta_same",    6, lower, upper,     5, -1.8, 1.8);
    deut_deut_dphi_deta_diff   = new TH2F("deut_deut_dphi_deta_diff",   "deut_deut_dphi_deta_diff",    6, lower, upper,     5, -1.8, 1.8);
    deut_deut_dpt_dcharge      = new TH2F("deut_deut_dpt_dcharge",      "deut_deut_dpt_dcharge",      40, 0.0, 8.0,         3,  0.0, 3.0);
    deut_deut_ptA_ptB_same     = new TH2F("deut_deut_ptA_ptB_same",     "deut_deut_ptA_ptB_same",     22, 0.0, 4.4,        22,  0.0, 4.4);
    deut_deut_ptA_ptB_diff     = new TH2F("deut_deut_ptA_ptB_diff",     "deut_deut_ptA_ptB_diff",     22, 0.0, 4.4,        22,  0.0, 4.4);

    deut_deut_inv_mass         = new TH1F("deut_deut_inv_mass",         "deut_deut_inv_mass",        120, 0.0,12.0);
    
  // objects added to output file

    fOutputList->Add(deut_eta_pos_pt1);
    fOutputList->Add(deut_eta_neg_pt1);
    fOutputList->Add(deut_eta_pos_pt2);
    fOutputList->Add(deut_eta_neg_pt2);
    fOutputList->Add(deut_eta_pos_pt3);
    fOutputList->Add(deut_eta_neg_pt3);
    fOutputList->Add(deut_eta_pos_pt4);
    fOutputList->Add(deut_eta_neg_pt4);
    fOutputList->Add(deut_eta_pos_pt5);
    fOutputList->Add(deut_eta_neg_pt5);
    
    fOutputList->Add(hadr_eta_pos);
    fOutputList->Add(hadr_eta_neg);

    fOutputList->Add(deut_hadr_dphi36_deta20_pt1);
    fOutputList->Add(deut_hadr_dphi36_deta20_pt2);
    fOutputList->Add(deut_hadr_dphi36_deta20_pt3);
    fOutputList->Add(deut_hadr_dphi36_deta20_pt4);
    fOutputList->Add(deut_hadr_dphi36_deta20_pt5);

    fOutputList->Add(deut_hadr_dphi36_deta21_pt1);
    fOutputList->Add(deut_hadr_dphi36_deta21_pt2);
    fOutputList->Add(deut_hadr_dphi36_deta21_pt3);
    fOutputList->Add(deut_hadr_dphi36_deta21_pt4);
    fOutputList->Add(deut_hadr_dphi36_deta21_pt5);

    
    fOutputList->Add(multiplicity_w_deut_1);           //  n1
    fOutputList->Add(multiplicity_w_deut_2);           //  n2
    fOutputList->Add(multiplicity_w_deut_3);           //  n3
    fOutputList->Add(multiplicity_w_deut_4);           //  n4
    fOutputList->Add(multiplicity_w_deut_5);           //  n5
    fOutputList->Add(multiplicity_w_trig_and_deut_1);  //  n6
    fOutputList->Add(multiplicity_w_trig_and_deut_2);  //  n7
    fOutputList->Add(multiplicity_w_trig_and_deut_3);  //  n8
    fOutputList->Add(multiplicity_w_trig_and_deut_4);  //  n9
    fOutputList->Add(multiplicity_w_trig_and_deut_5);  //  n10
    
    fOutputList->Add(deut_deut_dphi_deta_same);
    fOutputList->Add(deut_deut_dphi_deta_diff);
    fOutputList->Add(deut_deut_dpt_dcharge);
    fOutputList->Add(deut_deut_ptA_ptB_same);
    fOutputList->Add(deut_deut_ptA_ptB_diff);
    fOutputList->Add(deut_deut_inv_mass);
    
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

    fOutputList->Add(m2_pt_pos_cut_sa0);            //
    fOutputList->Add(m2_pt_neg_cut_sa0);            //
    fOutputList->Add(m2_pt_pos_cut_sa1);            // 15
    fOutputList->Add(m2_pt_neg_cut_sa1);            // 16
    fOutputList->Add(m2_pt_pos_cut_sa2);            // 17
    fOutputList->Add(m2_pt_neg_cut_sa2);            // 18
    fOutputList->Add(m2_pt_pos_cut_sa3);            //
    fOutputList->Add(m2_pt_neg_cut_sa3);            //
    fOutputList->Add(m2_pt_pos_cut_sa4);            //
    fOutputList->Add(m2_pt_neg_cut_sa4);            //
    fOutputList->Add(m2_pt_pos_cut_sa5);            //
    fOutputList->Add(m2_pt_neg_cut_sa5);            //
    fOutputList->Add(m2_pt_pos_cut_sa6);            //
    fOutputList->Add(m2_pt_neg_cut_sa6);            //
    fOutputList->Add(m2_pt_pos_cut_sa7);            //
    fOutputList->Add(m2_pt_neg_cut_sa7);            //
    
    fOutputList->Add(m2_pt_pos_cut_sb0);            //
    fOutputList->Add(m2_pt_neg_cut_sb0);            //
    fOutputList->Add(m2_pt_pos_cut_sb1);            // 15
    fOutputList->Add(m2_pt_neg_cut_sb1);            // 16
    fOutputList->Add(m2_pt_pos_cut_sb2);            // 17
    fOutputList->Add(m2_pt_neg_cut_sb2);            // 18
    fOutputList->Add(m2_pt_pos_cut_sb3);            //
    fOutputList->Add(m2_pt_neg_cut_sb3);            //
    fOutputList->Add(m2_pt_pos_cut_sb4);            //
    fOutputList->Add(m2_pt_neg_cut_sb4);            //
    fOutputList->Add(m2_pt_pos_cut_sb5);            //
    fOutputList->Add(m2_pt_neg_cut_sb5);            //
    fOutputList->Add(m2_pt_pos_cut_sb6);            //
    fOutputList->Add(m2_pt_neg_cut_sb6);            //
    fOutputList->Add(m2_pt_pos_cut_sb7);            //
    fOutputList->Add(m2_pt_neg_cut_sb7);            //


    fOutputList->Add(m2_pt_pos_cut_T_prot);         // 19
    fOutputList->Add(m2_pt_neg_cut_T_prot);         // 20
    
    fOutputList->Add(m2_pt_pos_cut_with_trig_05);   // 21
    fOutputList->Add(m2_pt_neg_cut_with_trig_05);   // 22

    fOutputList->Add(deut_phi_pt_pos_T);            // 23
    fOutputList->Add(deut_phi_pt_neg_T);            // 24
    fOutputList->Add(prot_phi_pt_pos_T);            // 25
    fOutputList->Add(prot_phi_pt_neg_T);            // 26

    fOutputList->Add(hadr_per_event);               // 35
    fOutputList->Add(deut_per_event);               // 35
    fOutputList->Add(deut_per_event_with_trig);

    fOutputList->Add(deut_phi_pt_pos_sb0);          //
    fOutputList->Add(deut_phi_pt_neg_sb0);          //
    fOutputList->Add(deut_phi_pt_pos_sb1);          // 31
    fOutputList->Add(deut_phi_pt_neg_sb1);          // 32
    fOutputList->Add(deut_phi_pt_pos_sb2);          // 33
    fOutputList->Add(deut_phi_pt_neg_sb2);          // 34
    fOutputList->Add(deut_phi_pt_pos_sb3);          //
    fOutputList->Add(deut_phi_pt_neg_sb3);          //
    fOutputList->Add(deut_phi_pt_pos_sb4);          //
    fOutputList->Add(deut_phi_pt_neg_sb4);          //
    fOutputList->Add(deut_phi_pt_pos_sb5);          //
    fOutputList->Add(deut_phi_pt_neg_sb5);          //
    fOutputList->Add(deut_phi_pt_pos_sb6);          //
    fOutputList->Add(deut_phi_pt_neg_sb6);          //
    fOutputList->Add(deut_phi_pt_pos_sb7);          //
    fOutputList->Add(deut_phi_pt_neg_sb7);          //
    
    fOutputList->Add(deut_phi_pt_pos_sa0);          //
    fOutputList->Add(deut_phi_pt_neg_sa0);          //
    fOutputList->Add(deut_phi_pt_pos_sa1);          // 31
    fOutputList->Add(deut_phi_pt_neg_sa1);          // 32
    fOutputList->Add(deut_phi_pt_pos_sa2);          // 33
    fOutputList->Add(deut_phi_pt_neg_sa2);          // 34
    fOutputList->Add(deut_phi_pt_pos_sa3);          //
    fOutputList->Add(deut_phi_pt_neg_sa3);          //
    fOutputList->Add(deut_phi_pt_pos_sa4);          //
    fOutputList->Add(deut_phi_pt_neg_sa4);          //
    fOutputList->Add(deut_phi_pt_pos_sa5);          //
    fOutputList->Add(deut_phi_pt_neg_sa5);          //
    fOutputList->Add(deut_phi_pt_pos_sa6);          //
    fOutputList->Add(deut_phi_pt_neg_sa6);          //
    fOutputList->Add(deut_phi_pt_pos_sa7);          //
    fOutputList->Add(deut_phi_pt_neg_sa7);          //

    fOutputList->Add(prot_per_event);               // 36
    
    
    fOutputList->Add(trig_05_per_event);            // 37
   
    fOutputList->Add(trig_05_phi_pt_pos);           // 38
    fOutputList->Add(trig_05_phi_pt_neg);           // 39


    fOutputList->Add(tof_phi_eta_pos);              // 40
    fOutputList->Add(tof_phi_eta_neg);              // 41
    fOutputList->Add(tof_phi_eta_pos_deut);         // 42
    fOutputList->Add(tof_phi_eta_neg_deut);         // 43

    fOutputList->Add(deut_pt_compare_pos);          // 44
    fOutputList->Add(deut_pt_compare_neg);          // 45
    fOutputList->Add(tpc_sector_fraction);          // 46
    fOutputList->Add(primary_vertex_z);             // 47
    fOutputList->Add(primary_vertex_z_cut1);        // 48
    fOutputList->Add(primary_vertex_z_cut2);        // 49

    fOutputList->Add(deut_dphi_pt_pos_pos_05_T);    // 50
    fOutputList->Add(deut_dphi_pt_pos_neg_05_T);    // 51
    fOutputList->Add(deut_dphi_pt_neg_neg_05_T);    // 52

    fOutputList->Add(prot_dphi_pt_pos_pos_05_T);    // 53
    fOutputList->Add(prot_dphi_pt_pos_neg_05_T);    // 54
    fOutputList->Add(prot_dphi_pt_neg_neg_05_T);    // 55

    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa0);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa0);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa0);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa1);  // 62
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa1);  // 63
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa1);  // 64
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa2);  // 65
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa2);  // 66
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa2);  // 67
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa3);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa3);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa3);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa4);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa4);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa4);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa5);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa5);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa5);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa6);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa6);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa6);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sa7);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sa7);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sa7);  //
    
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb0);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb0);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb0);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb1);  // 62
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb1);  // 63
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb1);  // 64
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb2);  // 65
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb2);  // 66
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb2);  // 67
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb3);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb3);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb3);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb4);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb4);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb4);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb5);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb5);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb5);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb6);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb6);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb6);  //
    fOutputList->Add(deut_dphi_pt_pos_pos_05_sb7);  //
    fOutputList->Add(deut_dphi_pt_pos_neg_05_sb7);  //
    fOutputList->Add(deut_dphi_pt_neg_neg_05_sb7);  //    

    fOutputList->Add(DCAxy_pos);                    // 68
    fOutputList->Add(DCAxy_neg);                    // 69
    fOutputList->Add(DCAz_pos);                     // 70
    fOutputList->Add(DCAz_neg);                     // 71

    fOutputList->Add(m2_pt_pos_fine);               // 72
    fOutputList->Add(m2_pt_neg_fine);               // 73
    fOutputList->Add(m2_pt_pos_TPC_fine);           // 74
    fOutputList->Add(m2_pt_neg_TPC_fine);           // 75 
    fOutputList->Add(m2_pt_pos_cut_T_fine);         // 76
    fOutputList->Add(m2_pt_neg_cut_T_fine);         // 77

    fOutputList->Add(m2_pt_pos_cut_sa0_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sa0_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sa1_fine);       // 82
    fOutputList->Add(m2_pt_neg_cut_sa1_fine);       // 83
    fOutputList->Add(m2_pt_pos_cut_sa2_fine);       // 84
    fOutputList->Add(m2_pt_neg_cut_sa2_fine);       // 85
    fOutputList->Add(m2_pt_pos_cut_sa3_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sa3_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sa4_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sa4_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sa5_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sa5_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sa6_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sa6_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sa7_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sa7_fine);       //
    
    fOutputList->Add(m2_pt_pos_cut_sb0_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sb0_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sb1_fine);       // 82
    fOutputList->Add(m2_pt_neg_cut_sb1_fine);       // 83
    fOutputList->Add(m2_pt_pos_cut_sb2_fine);       // 84
    fOutputList->Add(m2_pt_neg_cut_sb2_fine);       // 85
    fOutputList->Add(m2_pt_pos_cut_sb3_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sb3_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sb4_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sb4_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sb5_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sb5_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sb6_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sb6_fine);       //
    fOutputList->Add(m2_pt_pos_cut_sb7_fine);       //
    fOutputList->Add(m2_pt_neg_cut_sb7_fine);       //

    fOutputList->Add(m2_pt_pos_TPC_prot_fine);      // 86
    fOutputList->Add(m2_pt_neg_TPC_prot_fine);      // 87

    fOutputList->Add(m2_pt_pos_cut_T_prot_fine);    // 88
    fOutputList->Add(m2_pt_neg_cut_T_prot_fine);    // 89



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

    int deut_count_T_1 = 0;
    int deut_count_T_2 = 0;
    int deut_count_T_3 = 0;
    int deut_count_T_4 = 0;
    int deut_count_T_5 = 0;

    int hybrid_count   = 0;
	
    int deut_track_num_T[200];     int deut_count_T            = 0;

    int deut_track_num_sa0[200];   int deut_count_sa0          = 0;
    int deut_track_num_sa1[200];   int deut_count_sa1          = 0;
    int deut_track_num_sa2[200];   int deut_count_sa2          = 0;
    int deut_track_num_sa3[200];   int deut_count_sa3          = 0;
    int deut_track_num_sa4[200];   int deut_count_sa4          = 0;
    int deut_track_num_sa5[200];   int deut_count_sa5          = 0;
    int deut_track_num_sa6[200];   int deut_count_sa6          = 0;
    int deut_track_num_sa7[200];   int deut_count_sa7          = 0;
    
    int deut_track_num_sb0[200];   int deut_count_sb0          = 0;
    int deut_track_num_sb1[200];   int deut_count_sb1          = 0;
    int deut_track_num_sb2[200];   int deut_count_sb2          = 0;
    int deut_track_num_sb3[200];   int deut_count_sb3          = 0;
    int deut_track_num_sb4[200];   int deut_count_sb4          = 0;
    int deut_track_num_sb5[200];   int deut_count_sb5          = 0;
    int deut_track_num_sb6[200];   int deut_count_sb6          = 0;
    int deut_track_num_sb7[200];   int deut_count_sb7          = 0;

    int deut_track_num_W[200];     int deut_count_W            = 0;
    int prot_track_num_T[500];     int prot_count_T            = 0;
    
    int trig_05_track_num[200];    int trig_05_count     = 0;

    // loop over all these tracks
    //
    // pull out track numbers for high-pt triggers and also deutons

    Float_t max_05_pt = 0.0;    Int_t   max_05    = 0;

    for(Int_t i(0); i < iTracks; i++)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track)                                                                      {    continue;    }

	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }   // cut 1
	hybrid_count++;
//	if(!track->IsGlobalConstrained())                                               {    continue;    }
//	if(!track->IsTPCConstrained())                                                  {    continue;    }
//	if(!track->IsPrimaryCandidate())                                                {    continue;    }
	Float_t pt   = track->Pt();	 //     if(pt < 0.95)                           {    continue;    }   // cut 2   // makes for fast processings, low momentum play no part of analysis anyway
        Float_t dedx = track->GetTPCsignal();	if(dedx > 1000)                         {    continue;    }   //
	Float_t eta  = track->Eta();	        if(TMath::Abs(eta) > 0.9)               {    continue;    }


	int is_global = 0;
	if(!track->IsGlobalConstrained())                                               {    is_global=1; }




	Float_t phi           = track->Phi();
	if(phi <  -pio2){    phi = phi + twopi; }	if(phi <  -pio2){    phi = phi + twopi;   }
	if(phi > 3*pio2){    phi = phi - twopi;	}       if(phi > 3*pio2){    phi = phi - twopi;   }


	Float_t fraction = 0.0;
	fraction = track->GetTPCFoundFraction();
	if(run_mode == 3  &&  fraction < 0.90)                                          {    continue;    }  /////////////// OPTIONAL FOR SYSTEMATIC STUDY /////////////////
	tpc_sector_fraction->Fill(fraction);
	
	fHistPt->Fill(pt);

	// markA
	
	if(pt >= 5.0              )	{   trig_05_track_num[trig_05_count] = i;	  if(pt>max_05_pt){ max_05_pt = pt; max_05 = trig_05_count;   }   trig_05_count++;   }

	Short_t charge        = track->Charge();
	if(pt >= 5.0              )	{   if(charge > 0)   trig_05_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_05_phi_pt_neg->Fill(pt, phi);	}


	
	if(run_mode == 6  && track->IsGlobalConstrained())                              {    continue;    }

	
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

	if(run_mode == 1){   cut_width = 3;   }
	if(run_mode == 2){   TPC_PID = 2.0;   }

	
	if(charge > 0)
	{
//	    m2_pt_pos     ->Fill(pt, m2tof);
	    m2_pt_pos_fine->Fill(pt, m2tof);
	    if(pt >= 1.0  &&  pt < 2.0)   tof_phi_eta_pos->Fill(phi, eta);

	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < TPC_PID)  // standard, good
//	    if(nSigmaTPCDeut > 0.0 - TPC_PID  &&  nSigmaTPCDeut < 0.0)
	    {
		m2_pt_pos_TPC     ->Fill(pt,  m2tof);
		m2_pt_pos_TPC_fine->Fill(pt,  m2tof);
			
		if(pt >= 1.0  &&  pt < 5.4)
		{
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
		    deut_sigma = fit_deut_curve->Eval(pt);

		    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_T[deut_count_T] = i;
			deut_count_T++;
			if     (pt >= 1.00  &&  pt < 1.35)   deut_count_T_1++;
			else if(pt >= 1.35  &&  pt < 1.80)   deut_count_T_2++;
			else if(pt >= 1.80  &&  pt < 2.40)   deut_count_T_3++;
			else if(pt >= 2.40  &&  pt < 3.00)   deut_count_T_4++;
			else if(pt >= 3.00  &&  pt < 4.00)   deut_count_T_5++;
			
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
			if(is_global == 1){   m2_pt_pos_cut_G->Fill(pt,m2tof);   }
		    }
		    if(m2tof >= 2.5  &&  m2tof < 5.0)
		    {
			deut_track_num_W[deut_count_W] = i;
			deut_count_W++;
		    }


		    if(m2tof < deut_mean - 2.0 * deut_sigma  &&   m2tof > deut_mean - 2.5 * deut_sigma)
		    {
			deut_track_num_sa0[deut_count_sa0] = i;			deut_count_sa0++;
			m2_pt_pos_cut_sa0->Fill(pt,m2tof);			m2_pt_pos_cut_sa0_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa0->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 2.5 * deut_sigma  &&   m2tof > deut_mean - 3.0 * deut_sigma)
		    {
			deut_track_num_sa1[deut_count_sa1] = i;			deut_count_sa1++;
			m2_pt_pos_cut_sa1->Fill(pt,m2tof);			m2_pt_pos_cut_sa1_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa1->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 3.0 * deut_sigma  &&   m2tof > deut_mean - 3.5 * deut_sigma)
		    {
			deut_track_num_sa2[deut_count_sa2] = i;			deut_count_sa2++;
			m2_pt_pos_cut_sa2->Fill(pt,m2tof);			m2_pt_pos_cut_sa2_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa2->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 3.5 * deut_sigma  &&   m2tof > deut_mean - 4.0 * deut_sigma)
		    {
			deut_track_num_sa3[deut_count_sa3] = i;			deut_count_sa3++;
			m2_pt_pos_cut_sa3->Fill(pt,m2tof);			m2_pt_pos_cut_sa3_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa3->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 4.0 * deut_sigma  &&   m2tof > deut_mean - 4.5 * deut_sigma)
		    {
			deut_track_num_sa4[deut_count_sa4] = i;			deut_count_sa4++;
			m2_pt_pos_cut_sa4->Fill(pt,m2tof);			m2_pt_pos_cut_sa4_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa4->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 4.5 * deut_sigma  &&   m2tof > deut_mean - 5.0 * deut_sigma)
		    {
			deut_track_num_sa5[deut_count_sa5] = i;			deut_count_sa5++;
			m2_pt_pos_cut_sa5->Fill(pt,m2tof);			m2_pt_pos_cut_sa5_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa5->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 5.0 * deut_sigma  &&   m2tof > deut_mean - 5.5 * deut_sigma)
		    {
			deut_track_num_sa6[deut_count_sa6] = i;			deut_count_sa6++;
			m2_pt_pos_cut_sa6->Fill(pt,m2tof);			m2_pt_pos_cut_sa6_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa6->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 5.5 * deut_sigma  &&   m2tof > deut_mean - 6.0 * deut_sigma)
		    {
			deut_track_num_sa7[deut_count_sa7] = i;			deut_count_sa7++;
			m2_pt_pos_cut_sa7->Fill(pt,m2tof);			m2_pt_pos_cut_sa7_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sa7->Fill(pt, deut_phi);
		    }


		    if(m2tof > deut_mean + 2.0 * deut_sigma  &&   m2tof < deut_mean + 2.5 * deut_sigma)
		    {
			deut_track_num_sb0[deut_count_sb0] = i;			deut_count_sb0++;
			m2_pt_pos_cut_sb0->Fill(pt,m2tof);			m2_pt_pos_cut_sb0_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb0->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 2.5 * deut_sigma  &&   m2tof < deut_mean + 3.0 * deut_sigma)
		    {
			deut_track_num_sb1[deut_count_sb1] = i;			deut_count_sb1++;
			m2_pt_pos_cut_sb1->Fill(pt,m2tof);			m2_pt_pos_cut_sb1_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb1->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 3.0 * deut_sigma  &&   m2tof < deut_mean + 3.5 * deut_sigma)
		    {
			deut_track_num_sb2[deut_count_sb2] = i;			deut_count_sb2++;
			m2_pt_pos_cut_sb2->Fill(pt,m2tof);			m2_pt_pos_cut_sb2_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb2->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 3.5 * deut_sigma  &&   m2tof < deut_mean + 4.0 * deut_sigma)
		    {
			deut_track_num_sb3[deut_count_sb3] = i;			deut_count_sb3++;
			m2_pt_pos_cut_sb3->Fill(pt,m2tof);			m2_pt_pos_cut_sb3_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb3->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 4.0 * deut_sigma  &&   m2tof < deut_mean + 4.5 * deut_sigma)
		    {
			deut_track_num_sb4[deut_count_sb4] = i;			deut_count_sb4++;
			m2_pt_pos_cut_sb4->Fill(pt,m2tof);			m2_pt_pos_cut_sb4_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb4->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 4.5 * deut_sigma  &&   m2tof < deut_mean + 5.0 * deut_sigma)
		    {
			deut_track_num_sb5[deut_count_sb5] = i;			deut_count_sb5++;
			m2_pt_pos_cut_sb5->Fill(pt,m2tof);			m2_pt_pos_cut_sb5_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb5->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 5.0 * deut_sigma  &&   m2tof < deut_mean + 5.5 * deut_sigma)
		    {
			deut_track_num_sb6[deut_count_sb6] = i;			deut_count_sb6++;
			m2_pt_pos_cut_sb6->Fill(pt,m2tof);			m2_pt_pos_cut_sb6_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb6->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 5.5 * deut_sigma  &&   m2tof < deut_mean + 6.0 * deut_sigma)
		    {
			deut_track_num_sb7[deut_count_sb7] = i;			deut_count_sb7++;
			m2_pt_pos_cut_sb7->Fill(pt,m2tof);			m2_pt_pos_cut_sb7_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_pos_sb7->Fill(pt, deut_phi);
		    }

		    
		    if(m2tof >= 1.7  &&  m2tof < 5.5)    {   wide_cut++;   }
		}
	    }

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
	}    //   end of pos charge if statement
	    
	else if(charge < 0)
	{
	    m2_pt_neg_fine->Fill(pt, m2tof);
	    if(pt >= 1.0  &&  pt < 2.0)   tof_phi_eta_neg->Fill(phi, eta);

	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < TPC_PID)  // standard, good
//	    if(nSigmaTPCDeut < TPC_PID  &&  nSigmaTPCDeut >= 0.0)
//	    if(nSigmaTPCDeut > 0.0 - TPC_PID  &&  nSigmaTPCDeut < 0.0)
	    {
		m2_pt_neg_TPC     ->Fill(pt, m2tof);
		m2_pt_neg_TPC_fine->Fill(pt, m2tof);

		if(pt >= 1.0  &&  pt < 5.4)
		{
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<6; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][1][w]);   }
		    deut_sigma = fit_deut_curve->Eval(pt);

		    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
		    {
			deut_track_num_T[deut_count_T] = i;
			deut_count_T++;
			if     (pt >= 1.00  &&  pt < 1.35)   deut_count_T_1++;
			else if(pt >= 1.35  &&  pt < 1.80)   deut_count_T_2++;
			else if(pt >= 1.80  &&  pt < 2.40)   deut_count_T_3++;
			else if(pt >= 2.40  &&  pt < 3.00)   deut_count_T_4++;
			else if(pt >= 3.00  &&  pt < 4.00)   deut_count_T_5++;
			
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
		    



		    if(m2tof < deut_mean - 2.0 * deut_sigma  &&   m2tof > deut_mean - 2.5 * deut_sigma)
		    {
			deut_track_num_sa0[deut_count_sa0] = i;			deut_count_sa0++;
			m2_pt_neg_cut_sa0->Fill(pt,m2tof);			m2_pt_neg_cut_sa0_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa0->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 2.5 * deut_sigma  &&   m2tof > deut_mean - 3.0 * deut_sigma)
		    {
			deut_track_num_sa1[deut_count_sa1] = i;			deut_count_sa1++;
			m2_pt_neg_cut_sa1->Fill(pt,m2tof);			m2_pt_neg_cut_sa1_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa1->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 3.0 * deut_sigma  &&   m2tof > deut_mean - 3.5 * deut_sigma)
		    {
			deut_track_num_sa2[deut_count_sa2] = i;			deut_count_sa2++;
			m2_pt_neg_cut_sa2->Fill(pt,m2tof);			m2_pt_neg_cut_sa2_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa2->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 3.5 * deut_sigma  &&   m2tof > deut_mean - 4.0 * deut_sigma)
		    {
			deut_track_num_sa3[deut_count_sa3] = i;			deut_count_sa3++;
			m2_pt_neg_cut_sa3->Fill(pt,m2tof);			m2_pt_neg_cut_sa3_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa3->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 4.0 * deut_sigma  &&   m2tof > deut_mean - 4.5 * deut_sigma)
		    {
			deut_track_num_sa4[deut_count_sa4] = i;			deut_count_sa4++;
			m2_pt_neg_cut_sa4->Fill(pt,m2tof);			m2_pt_neg_cut_sa4_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa4->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 4.5 * deut_sigma  &&   m2tof > deut_mean - 5.0 * deut_sigma)
		    {
			deut_track_num_sa5[deut_count_sa5] = i;			deut_count_sa5++;
			m2_pt_neg_cut_sa5->Fill(pt,m2tof);			m2_pt_neg_cut_sa5_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa5->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 5.0 * deut_sigma  &&   m2tof > deut_mean - 5.5 * deut_sigma)
		    {
			deut_track_num_sa6[deut_count_sa6] = i;			deut_count_sa6++;
			m2_pt_neg_cut_sa6->Fill(pt,m2tof);			m2_pt_neg_cut_sa6_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa6->Fill(pt, deut_phi);
		    }
		    else if(m2tof < deut_mean - 5.5 * deut_sigma  &&   m2tof > deut_mean - 6.0 * deut_sigma)
		    {
			deut_track_num_sa7[deut_count_sa7] = i;			deut_count_sa7++;
			m2_pt_neg_cut_sa7->Fill(pt,m2tof);			m2_pt_neg_cut_sa7_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sa7->Fill(pt, deut_phi);
		    }


		    if(m2tof > deut_mean + 2.0 * deut_sigma  &&   m2tof < deut_mean + 2.5 * deut_sigma)
		    {
			deut_track_num_sb0[deut_count_sb0] = i;			deut_count_sb0++;
			m2_pt_neg_cut_sb0->Fill(pt,m2tof);			m2_pt_neg_cut_sb0_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb0->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 2.5 * deut_sigma  &&   m2tof < deut_mean + 3.0 * deut_sigma)
		    {
			deut_track_num_sb1[deut_count_sb1] = i;			deut_count_sb1++;
			m2_pt_neg_cut_sb1->Fill(pt,m2tof);			m2_pt_neg_cut_sb1_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb1->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 3.0 * deut_sigma  &&   m2tof < deut_mean + 3.5 * deut_sigma)
		    {
			deut_track_num_sb2[deut_count_sb2] = i;			deut_count_sb2++;
			m2_pt_neg_cut_sb2->Fill(pt,m2tof);			m2_pt_neg_cut_sb2_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb2->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 3.5 * deut_sigma  &&   m2tof < deut_mean + 4.0 * deut_sigma)
		    {
			deut_track_num_sb3[deut_count_sb3] = i;			deut_count_sb3++;
			m2_pt_neg_cut_sb3->Fill(pt,m2tof);			m2_pt_neg_cut_sb3_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb3->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 4.0 * deut_sigma  &&   m2tof < deut_mean + 4.5 * deut_sigma)
		    {
			deut_track_num_sb4[deut_count_sb4] = i;			deut_count_sb4++;
			m2_pt_neg_cut_sb4->Fill(pt,m2tof);			m2_pt_neg_cut_sb4_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb4->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 4.5 * deut_sigma  &&   m2tof < deut_mean + 5.0 * deut_sigma)
		    {
			deut_track_num_sb5[deut_count_sb5] = i;			deut_count_sb5++;
			m2_pt_neg_cut_sb5->Fill(pt,m2tof);			m2_pt_neg_cut_sb5_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb5->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 5.0 * deut_sigma  &&   m2tof < deut_mean + 5.5 * deut_sigma)
		    {
			deut_track_num_sb6[deut_count_sb6] = i;			deut_count_sb6++;
			m2_pt_neg_cut_sb6->Fill(pt,m2tof);			m2_pt_neg_cut_sb6_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb6->Fill(pt, deut_phi);
		    }
		    else if(m2tof > deut_mean + 5.5 * deut_sigma  &&   m2tof < deut_mean + 6.0 * deut_sigma)
		    {
			deut_track_num_sb7[deut_count_sb7] = i;			deut_count_sb7++;
			m2_pt_neg_cut_sb7->Fill(pt,m2tof);			m2_pt_neg_cut_sb7_fine->Fill(pt,m2tof);
			Float_t deut_phi = track->Phi();
			if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }	if(deut_phi <  -pio2){    deut_phi = deut_phi + twopi;   }
			if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }	if(deut_phi > 3*pio2){    deut_phi = deut_phi - twopi;   }
			deut_phi_pt_neg_sb7->Fill(pt, deut_phi);
		    }


		    if(m2tof >= 1.7  &&  m2tof < 5.5)    {   wide_cut++;   }
		}
	    }


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
	}    //   end of neg charge if statement
    }        //   end of track loop

    
    if(deut_count_T_1 > 0)                       {   multiplicity_w_deut_1->Fill(hybrid_count);	           }
    if(deut_count_T_2 > 0)                       {   multiplicity_w_deut_2->Fill(hybrid_count);	           }
    if(deut_count_T_3 > 0)                       {   multiplicity_w_deut_3->Fill(hybrid_count);	           }
    if(deut_count_T_4 > 0)                       {   multiplicity_w_deut_4->Fill(hybrid_count);	           }
    if(deut_count_T_5 > 0)                       {   multiplicity_w_deut_5->Fill(hybrid_count);	           }
    
    if(trig_05_count > 0  &&  deut_count_T_1 > 0) {   multiplicity_w_trig_and_deut_1->Fill(hybrid_count);   }
    if(trig_05_count > 0  &&  deut_count_T_2 > 0) {   multiplicity_w_trig_and_deut_2->Fill(hybrid_count);   }
    if(trig_05_count > 0  &&  deut_count_T_3 > 0) {   multiplicity_w_trig_and_deut_3->Fill(hybrid_count);   }
    if(trig_05_count > 0  &&  deut_count_T_4 > 0) {   multiplicity_w_trig_and_deut_4->Fill(hybrid_count);   }
    if(trig_05_count > 0  &&  deut_count_T_5 > 0) {   multiplicity_w_trig_and_deut_5->Fill(hybrid_count);   }
	
    
    hadr_per_event->Fill(hybrid_count);
    deut_per_event->Fill(deut_count_T);
    prot_per_event->Fill(prot_count_T);

    if(trig_05_count > 0)    {       deut_per_event_with_trig->Fill(deut_count_T);   }
    trig_05_per_event->Fill(trig_05_count); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(deut_count_T == 2)
    {
	int A = deut_track_num_T[0];
	int B = deut_track_num_T[1];
	
	AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));

	if(trackA)
	{
	    if(trackB)
	    {
		Short_t charge_A    = trackA->Charge();
		Float_t pt_A        = trackA->Pt();
		Float_t phi_A       = trackA->Phi();
		Float_t eta_A       = trackA->Eta();
		
		Short_t charge_B    = trackB->Charge();
		Float_t pt_B        = trackB->Pt();
		Float_t phi_B       = trackB->Phi();
		Float_t eta_B       = trackB->Eta();

		Float_t Sdphi       =  phi_A - phi_B;
		if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }


		Float_t dcharge = charge_A - charge_B;		if(dcharge < 0.0)  dcharge = -1.0 * dcharge;

		if(dcharge < 1){   if(pt_A < 2.0  && pt_B < 2.0){   deut_deut_dphi_deta_same->Fill(Sdphi, eta_A - eta_B);   }		deut_deut_ptA_ptB_same->Fill(pt_A, pt_B);   }
		else           {   if(pt_A < 2.0  && pt_B < 2.0){   deut_deut_dphi_deta_diff->Fill(Sdphi, eta_A - eta_B);   }		deut_deut_ptA_ptB_diff->Fill(pt_A, pt_B);   }
		
		Float_t dpt     = pt_A - pt_B;		        if(dpt     < 0.0)  dpt     = -1.0 * dpt;
		deut_deut_dpt_dcharge->Fill(dpt, dcharge);

		Float_t mass = 1.8756;
		Float_t p_A  = trackA->P();		Float_t px_A = trackA->Px();		Float_t py_A = trackA->Py();		Float_t pz_A = trackA->Pz();
		Float_t p_B  = trackB->P();		Float_t px_B = trackB->Px();		Float_t py_B = trackB->Py();		Float_t pz_B = trackB->Pz();

		
		deut_deut_inv_mass->Fill(sqrt(2*(mass*mass + sqrt(mass*mass+p_A*p_A)*sqrt(mass*mass+p_B*p_B)-(px_A*px_B + py_A*py_B + pz_A*pz_B)) ));

	    }
	}
    }	

    int begin = 0;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_T > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_T > 0  &&  trig_05_count > 0)
    {
	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{
	    int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }

	         if(charge_H > 0)   hadr_eta_pos->Fill(eta_H);
	    else if(charge_H < 0)   hadr_eta_neg->Fill(eta_H);
	    
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
			Float_t eta_A       = trackA->Eta();
			if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
			if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
			Float_t Sdphi       =  phi_A - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			    if     (charge_H > 0){  if     (charge_A > 0){    deut_dphi_pt_pos_pos_05_T->Fill(pt_A, Sdphi);	 }
				else if(charge_A < 0){    deut_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }    }
			    else if(charge_H < 0){  if     (charge_A > 0){    deut_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }
				else if(charge_A < 0){    deut_dphi_pt_neg_neg_05_T->Fill(pt_A, Sdphi);	 }    }
			}

			if(pt_A > 1.00  &&  pt_A < 1.35)
			{
			    Float_t deta = 0.0;
			    deta = eta_H - eta_A;
			    deut_hadr_dphi36_deta20_pt1->Fill(Sdphi, deta);
			    deut_hadr_dphi36_deta21_pt1->Fill(Sdphi, deta);
			}

			else if(pt_A > 1.35  &&  pt_A < 1.80)
			{
			    Float_t deta = 0.0;
			    deta = eta_H - eta_A;
			    deut_hadr_dphi36_deta20_pt2->Fill(Sdphi, deta);
			    deut_hadr_dphi36_deta21_pt2->Fill(Sdphi, deta);
			}

			else if(pt_A > 1.80  &&  pt_A < 2.40)
			{
			    Float_t deta = 0.0;
			    deta = eta_H - eta_A;
			    deut_hadr_dphi36_deta20_pt3->Fill(Sdphi, deta);
			    deut_hadr_dphi36_deta21_pt3->Fill(Sdphi, deta);
			}

			else if(pt_A > 2.40  &&  pt_A < 3.00)
			{
			    Float_t deta = 0.0;
			    deta = eta_H - eta_A;
			    deut_hadr_dphi36_deta20_pt4->Fill(Sdphi, deta);
			    deut_hadr_dphi36_deta21_pt4->Fill(Sdphi, deta);
			}

			else if(pt_A > 3.00  &&  pt_A < 4.00)
			{
			    Float_t deta = 0.0;
			    deta = eta_H - eta_A;
			    deut_hadr_dphi36_deta20_pt5->Fill(Sdphi, deta);
			    deut_hadr_dphi36_deta21_pt5->Fill(Sdphi, deta);
			}
		    }
		}
	    }
	}
//	cout<<deut_count_T<<endl;
	for(int j=0; j<deut_count_T; j++)
	{
	    int A               = deut_track_num_T[j];
	    AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	    if(trackA)
	    {
		Short_t charge_A    = trackA->Charge();
		Float_t eta_A       = trackA->Eta();
		Float_t pt_A        = trackA->Pt();

		if(pt_A > 1.00  &&  pt_A < 1.35)
		{
		         if(charge_A > 0)   deut_eta_pos_pt1->Fill(eta_A);
		    else if(charge_A < 0)   deut_eta_neg_pt1->Fill(eta_A); 
		}
		else if(pt_A > 1.35  &&  pt_A < 1.80)
		{
		         if(charge_A > 0)   deut_eta_pos_pt2->Fill(eta_A);
		    else if(charge_A < 0)   deut_eta_neg_pt2->Fill(eta_A); 
		}
		else if(pt_A > 1.80  &&  pt_A < 2.40)
		{
		         if(charge_A > 0)   deut_eta_pos_pt3->Fill(eta_A);
		    else if(charge_A < 0)   deut_eta_neg_pt3->Fill(eta_A); 
		}
		else if(pt_A > 2.40  &&  pt_A < 3.00)
		{
		         if(charge_A > 0)   deut_eta_pos_pt4->Fill(eta_A);
		    else if(charge_A < 0)   deut_eta_neg_pt4->Fill(eta_A); 
		}
		else if(pt_A > 3.00  &&  pt_A < 4.00)
		{
		         if(charge_A > 0)   deut_eta_pos_pt5->Fill(eta_A);
		    else if(charge_A < 0)   deut_eta_neg_pt5->Fill(eta_A); 
		}
		
	    }
	}
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(prot_count_T > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(prot_count_T > 0  &&  trig_05_count > 0)
    {
	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{
	    int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
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
		    Float_t eta_A       = trackA->Eta();
		    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }    if(phi_A <  -pio2){    phi_A = phi_A + twopi;   }
		    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }    if(phi_A > 3*pio2){    phi_A = phi_A - twopi;   }
		    Float_t Sdphi       =  phi_A - phi_H;
		    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
		    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
		    if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
		    {
		    if     (charge_H > 0){  if     (charge_A > 0){    prot_dphi_pt_pos_pos_05_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    prot_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }    }
		    else if(charge_H < 0){  if     (charge_A > 0){    prot_dphi_pt_pos_neg_05_T->Fill(pt_A, Sdphi);	 }
			                    else if(charge_A < 0){    prot_dphi_pt_neg_neg_05_T->Fill(pt_A, Sdphi);	 }    }
		    }
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    begin = 0;
    if(deut_count_sa0 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa0 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa0; j++)
	    {	int A               = deut_track_num_sa0[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa0    = trackA->Charge();	    Float_t pt_sa0      = trackA->Pt();
			Float_t phi_sa0       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa0 <  -pio2){    phi_sa0 = phi_sa0 + twopi;   }    if(phi_sa0 <  -pio2){    phi_sa0 = phi_sa0 + twopi;   }
			if(phi_sa0 > 3*pio2){    phi_sa0 = phi_sa0 - twopi;   }    if(phi_sa0 > 3*pio2){    phi_sa0 = phi_sa0 - twopi;   }
			Float_t Sdphi       =  phi_sa0 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa0 > 0){    deut_dphi_pt_pos_pos_05_sa0->Fill(pt_sa0, Sdphi);	 }
			                        else if(charge_sa0 < 0){    deut_dphi_pt_pos_neg_05_sa0->Fill(pt_sa0, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa0 > 0){    deut_dphi_pt_pos_neg_05_sa0->Fill(pt_sa0, Sdphi);	 }
			                        else if(charge_sa0 < 0){    deut_dphi_pt_neg_neg_05_sa0->Fill(pt_sa0, Sdphi);	 }    }
			}
    }   }   }   }   }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa1 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa1 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa1; j++)
	    {	int A               = deut_track_num_sa1[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa1    = trackA->Charge();	    Float_t pt_sa1      = trackA->Pt();
			Float_t phi_sa1       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa1 <  -pio2){    phi_sa1 = phi_sa1 + twopi;   }    if(phi_sa1 <  -pio2){    phi_sa1 = phi_sa1 + twopi;   }
			if(phi_sa1 > 3*pio2){    phi_sa1 = phi_sa1 - twopi;   }    if(phi_sa1 > 3*pio2){    phi_sa1 = phi_sa1 - twopi;   }
			Float_t Sdphi       =  phi_sa1 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa1 > 0){    deut_dphi_pt_pos_pos_05_sa1->Fill(pt_sa1, Sdphi);	 }
			                        else if(charge_sa1 < 0){    deut_dphi_pt_pos_neg_05_sa1->Fill(pt_sa1, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa1 > 0){    deut_dphi_pt_pos_neg_05_sa1->Fill(pt_sa1, Sdphi);	 }
			                        else if(charge_sa1 < 0){    deut_dphi_pt_neg_neg_05_sa1->Fill(pt_sa1, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa2 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa2 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa2; j++)
	    {	int A               = deut_track_num_sa2[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa2    = trackA->Charge();	    Float_t pt_sa2      = trackA->Pt();
			Float_t phi_sa2       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa2 <  -pio2){    phi_sa2 = phi_sa2 + twopi;   }    if(phi_sa2 <  -pio2){    phi_sa2 = phi_sa2 + twopi;   }
			if(phi_sa2 > 3*pio2){    phi_sa2 = phi_sa2 - twopi;   }    if(phi_sa2 > 3*pio2){    phi_sa2 = phi_sa2 - twopi;   }
			Float_t Sdphi       =  phi_sa2 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa2 > 0){    deut_dphi_pt_pos_pos_05_sa2->Fill(pt_sa2, Sdphi);	 }
			                        else if(charge_sa2 < 0){    deut_dphi_pt_pos_neg_05_sa2->Fill(pt_sa2, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa2 > 0){    deut_dphi_pt_pos_neg_05_sa2->Fill(pt_sa2, Sdphi);	 }
			                        else if(charge_sa2 < 0){    deut_dphi_pt_neg_neg_05_sa2->Fill(pt_sa2, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa3 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa3 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa3; j++)
	    {	int A               = deut_track_num_sa3[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa3    = trackA->Charge();	    Float_t pt_sa3      = trackA->Pt();
			Float_t phi_sa3       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa3 <  -pio2){    phi_sa3 = phi_sa3 + twopi;   }    if(phi_sa3 <  -pio2){    phi_sa3 = phi_sa3 + twopi;   }
			if(phi_sa3 > 3*pio2){    phi_sa3 = phi_sa3 - twopi;   }    if(phi_sa3 > 3*pio2){    phi_sa3 = phi_sa3 - twopi;   }
			Float_t Sdphi       =  phi_sa3 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa3 > 0){    deut_dphi_pt_pos_pos_05_sa3->Fill(pt_sa3, Sdphi);	 }
			                        else if(charge_sa3 < 0){    deut_dphi_pt_pos_neg_05_sa3->Fill(pt_sa3, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa3 > 0){    deut_dphi_pt_pos_neg_05_sa3->Fill(pt_sa3, Sdphi);	 }
			                        else if(charge_sa3 < 0){    deut_dphi_pt_neg_neg_05_sa3->Fill(pt_sa3, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa4 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa4 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa4; j++)
	    {	int A               = deut_track_num_sa4[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa4    = trackA->Charge();	    Float_t pt_sa4      = trackA->Pt();
			Float_t phi_sa4       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa4 <  -pio2){    phi_sa4 = phi_sa4 + twopi;   }    if(phi_sa4 <  -pio2){    phi_sa4 = phi_sa4 + twopi;   }
			if(phi_sa4 > 3*pio2){    phi_sa4 = phi_sa4 - twopi;   }    if(phi_sa4 > 3*pio2){    phi_sa4 = phi_sa4 - twopi;   }
			Float_t Sdphi       =  phi_sa4 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa4 > 0){    deut_dphi_pt_pos_pos_05_sa4->Fill(pt_sa4, Sdphi);	 }
			                        else if(charge_sa4 < 0){    deut_dphi_pt_pos_neg_05_sa4->Fill(pt_sa4, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa4 > 0){    deut_dphi_pt_pos_neg_05_sa4->Fill(pt_sa4, Sdphi);	 }
			                        else if(charge_sa4 < 0){    deut_dphi_pt_neg_neg_05_sa4->Fill(pt_sa4, Sdphi);	 }    }
			}
    }   }   }   }   }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa5 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa5 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa5; j++)
	    {	int A               = deut_track_num_sa5[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa5    = trackA->Charge();	    Float_t pt_sa5      = trackA->Pt();
			Float_t phi_sa5       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa5 <  -pio2){    phi_sa5 = phi_sa5 + twopi;   }    if(phi_sa5 <  -pio2){    phi_sa5 = phi_sa5 + twopi;   }
			if(phi_sa5 > 3*pio2){    phi_sa5 = phi_sa5 - twopi;   }    if(phi_sa5 > 3*pio2){    phi_sa5 = phi_sa5 - twopi;   }
			Float_t Sdphi       =  phi_sa5 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa5 > 0){    deut_dphi_pt_pos_pos_05_sa5->Fill(pt_sa5, Sdphi);	 }
			                        else if(charge_sa5 < 0){    deut_dphi_pt_pos_neg_05_sa5->Fill(pt_sa5, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa5 > 0){    deut_dphi_pt_pos_neg_05_sa5->Fill(pt_sa5, Sdphi);	 }
			                        else if(charge_sa5 < 0){    deut_dphi_pt_neg_neg_05_sa5->Fill(pt_sa5, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa6 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa6 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa6; j++)
	    {	int A               = deut_track_num_sa6[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa6    = trackA->Charge();	    Float_t pt_sa6      = trackA->Pt();
			Float_t phi_sa6       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa6 <  -pio2){    phi_sa6 = phi_sa6 + twopi;   }    if(phi_sa6 <  -pio2){    phi_sa6 = phi_sa6 + twopi;   }
			if(phi_sa6 > 3*pio2){    phi_sa6 = phi_sa6 - twopi;   }    if(phi_sa6 > 3*pio2){    phi_sa6 = phi_sa6 - twopi;   }
			Float_t Sdphi       =  phi_sa6 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa6 > 0){    deut_dphi_pt_pos_pos_05_sa6->Fill(pt_sa6, Sdphi);	 }
			                        else if(charge_sa6 < 0){    deut_dphi_pt_pos_neg_05_sa6->Fill(pt_sa6, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa6 > 0){    deut_dphi_pt_pos_neg_05_sa6->Fill(pt_sa6, Sdphi);	 }
			                        else if(charge_sa6 < 0){    deut_dphi_pt_neg_neg_05_sa6->Fill(pt_sa6, Sdphi);	 }    }
			}
    }   }   }   }   }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    begin = 0;
    if(deut_count_sa7 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sa7 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sa7; j++)
	    {	int A               = deut_track_num_sa7[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sa7    = trackA->Charge();	    Float_t pt_sa7      = trackA->Pt();
			Float_t phi_sa7       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sa7 <  -pio2){    phi_sa7 = phi_sa7 + twopi;   }    if(phi_sa7 <  -pio2){    phi_sa7 = phi_sa7 + twopi;   }
			if(phi_sa7 > 3*pio2){    phi_sa7 = phi_sa7 - twopi;   }    if(phi_sa7 > 3*pio2){    phi_sa7 = phi_sa7 - twopi;   }
			Float_t Sdphi       =  phi_sa7 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sa7 > 0){    deut_dphi_pt_pos_pos_05_sa7->Fill(pt_sa7, Sdphi);	 }
			                        else if(charge_sa7 < 0){    deut_dphi_pt_pos_neg_05_sa7->Fill(pt_sa7, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sa7 > 0){    deut_dphi_pt_pos_neg_05_sa7->Fill(pt_sa7, Sdphi);	 }
			                        else if(charge_sa7 < 0){    deut_dphi_pt_neg_neg_05_sa7->Fill(pt_sa7, Sdphi);	 }    }
			}
    }   }   }   }   }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    begin = 0;
    if(deut_count_sb0 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb0 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb0; j++)
	    {	int A               = deut_track_num_sb0[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb0    = trackA->Charge();	    Float_t pt_sb0      = trackA->Pt();
			Float_t phi_sb0       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb0 <  -pio2){    phi_sb0 = phi_sb0 + twopi;   }    if(phi_sb0 <  -pio2){    phi_sb0 = phi_sb0 + twopi;   }
			if(phi_sb0 > 3*pio2){    phi_sb0 = phi_sb0 - twopi;   }    if(phi_sb0 > 3*pio2){    phi_sb0 = phi_sb0 - twopi;   }
			Float_t Sdphi       =  phi_sb0 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb0 > 0){    deut_dphi_pt_pos_pos_05_sb0->Fill(pt_sb0, Sdphi);	 }
			                        else if(charge_sb0 < 0){    deut_dphi_pt_pos_neg_05_sb0->Fill(pt_sb0, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb0 > 0){    deut_dphi_pt_pos_neg_05_sb0->Fill(pt_sb0, Sdphi);	 }
			                        else if(charge_sb0 < 0){    deut_dphi_pt_neg_neg_05_sb0->Fill(pt_sb0, Sdphi);	 }    }
			}
    }   }   }   }   }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    begin = 0;
    if(deut_count_sb1 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb1 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb1; j++)
	    {	int A               = deut_track_num_sb1[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb1    = trackA->Charge();	    Float_t pt_sb1      = trackA->Pt();
			Float_t phi_sb1       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb1 <  -pio2){    phi_sb1 = phi_sb1 + twopi;   }    if(phi_sb1 <  -pio2){    phi_sb1 = phi_sb1 + twopi;   }
			if(phi_sb1 > 3*pio2){    phi_sb1 = phi_sb1 - twopi;   }    if(phi_sb1 > 3*pio2){    phi_sb1 = phi_sb1 - twopi;   }
			Float_t Sdphi       =  phi_sb1 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb1 > 0){    deut_dphi_pt_pos_pos_05_sb1->Fill(pt_sb1, Sdphi);	 }
			                        else if(charge_sb1 < 0){    deut_dphi_pt_pos_neg_05_sb1->Fill(pt_sb1, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb1 > 0){    deut_dphi_pt_pos_neg_05_sb1->Fill(pt_sb1, Sdphi);	 }
			                        else if(charge_sb1 < 0){    deut_dphi_pt_neg_neg_05_sb1->Fill(pt_sb1, Sdphi);	 }    }
			}
    }   }   }   }   }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    begin = 0;
    if(deut_count_sb2 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb2 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb2; j++)
	    {	int A               = deut_track_num_sb2[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb2    = trackA->Charge();	    Float_t pt_sb2      = trackA->Pt();
			Float_t phi_sb2       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb2 <  -pio2){    phi_sb2 = phi_sb2 + twopi;   }    if(phi_sb2 <  -pio2){    phi_sb2 = phi_sb2 + twopi;   }
			if(phi_sb2 > 3*pio2){    phi_sb2 = phi_sb2 - twopi;   }    if(phi_sb2 > 3*pio2){    phi_sb2 = phi_sb2 - twopi;   }
			Float_t Sdphi       =  phi_sb2 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb2 > 0){    deut_dphi_pt_pos_pos_05_sb2->Fill(pt_sb2, Sdphi);	 }
			                        else if(charge_sb2 < 0){    deut_dphi_pt_pos_neg_05_sb2->Fill(pt_sb2, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb2 > 0){    deut_dphi_pt_pos_neg_05_sb2->Fill(pt_sb2, Sdphi);	 }
			                        else if(charge_sb2 < 0){    deut_dphi_pt_neg_neg_05_sb2->Fill(pt_sb2, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_sb3 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb3 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb3; j++)
	    {	int A               = deut_track_num_sb3[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb3    = trackA->Charge();	    Float_t pt_sb3      = trackA->Pt();
			Float_t phi_sb3       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb3 <  -pio2){    phi_sb3 = phi_sb3 + twopi;   }    if(phi_sb3 <  -pio2){    phi_sb3 = phi_sb3 + twopi;   }
			if(phi_sb3 > 3*pio2){    phi_sb3 = phi_sb3 - twopi;   }    if(phi_sb3 > 3*pio2){    phi_sb3 = phi_sb3 - twopi;   }
			Float_t Sdphi       =  phi_sb3 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb3 > 0){    deut_dphi_pt_pos_pos_05_sb3->Fill(pt_sb3, Sdphi);	 }
			                        else if(charge_sb3 < 0){    deut_dphi_pt_pos_neg_05_sb3->Fill(pt_sb3, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb3 > 0){    deut_dphi_pt_pos_neg_05_sb3->Fill(pt_sb3, Sdphi);	 }
			                        else if(charge_sb3 < 0){    deut_dphi_pt_neg_neg_05_sb3->Fill(pt_sb3, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_sb4 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb4 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb4; j++)
	    {	int A               = deut_track_num_sb4[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb4    = trackA->Charge();	    Float_t pt_sb4      = trackA->Pt();
			Float_t phi_sb4       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb4 <  -pio2){    phi_sb4 = phi_sb4 + twopi;   }    if(phi_sb4 <  -pio2){    phi_sb4 = phi_sb4 + twopi;   }
			if(phi_sb4 > 3*pio2){    phi_sb4 = phi_sb4 - twopi;   }    if(phi_sb4 > 3*pio2){    phi_sb4 = phi_sb4 - twopi;   }
			Float_t Sdphi       =  phi_sb4 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb4 > 0){    deut_dphi_pt_pos_pos_05_sb4->Fill(pt_sb4, Sdphi);	 }
			                        else if(charge_sb4 < 0){    deut_dphi_pt_pos_neg_05_sb4->Fill(pt_sb4, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb4 > 0){    deut_dphi_pt_pos_neg_05_sb4->Fill(pt_sb4, Sdphi);	 }
			                        else if(charge_sb4 < 0){    deut_dphi_pt_neg_neg_05_sb4->Fill(pt_sb4, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_sb5 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb5 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb5; j++)
	    {	int A               = deut_track_num_sb5[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb5    = trackA->Charge();	    Float_t pt_sb5      = trackA->Pt();
			Float_t phi_sb5       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb5 <  -pio2){    phi_sb5 = phi_sb5 + twopi;   }    if(phi_sb5 <  -pio2){    phi_sb5 = phi_sb5 + twopi;   }
			if(phi_sb5 > 3*pio2){    phi_sb5 = phi_sb5 - twopi;   }    if(phi_sb5 > 3*pio2){    phi_sb5 = phi_sb5 - twopi;   }
			Float_t Sdphi       =  phi_sb5 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb5 > 0){    deut_dphi_pt_pos_pos_05_sb5->Fill(pt_sb5, Sdphi);	 }
			                        else if(charge_sb5 < 0){    deut_dphi_pt_pos_neg_05_sb5->Fill(pt_sb5, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb5 > 0){    deut_dphi_pt_pos_neg_05_sb5->Fill(pt_sb5, Sdphi);	 }
			                        else if(charge_sb5 < 0){    deut_dphi_pt_neg_neg_05_sb5->Fill(pt_sb5, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_sb6 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb6 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb6; j++)
	    {	int A               = deut_track_num_sb6[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb6    = trackA->Charge();	    Float_t pt_sb6      = trackA->Pt();
			Float_t phi_sb6       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb6 <  -pio2){    phi_sb6 = phi_sb6 + twopi;   }    if(phi_sb6 <  -pio2){    phi_sb6 = phi_sb6 + twopi;   }
			if(phi_sb6 > 3*pio2){    phi_sb6 = phi_sb6 - twopi;   }    if(phi_sb6 > 3*pio2){    phi_sb6 = phi_sb6 - twopi;   }
			Float_t Sdphi       =  phi_sb6 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb6 > 0){    deut_dphi_pt_pos_pos_05_sb6->Fill(pt_sb6, Sdphi);	 }
			                        else if(charge_sb6 < 0){    deut_dphi_pt_pos_neg_05_sb6->Fill(pt_sb6, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb6 > 0){    deut_dphi_pt_pos_neg_05_sb6->Fill(pt_sb6, Sdphi);	 }
			                        else if(charge_sb6 < 0){    deut_dphi_pt_neg_neg_05_sb6->Fill(pt_sb6, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    begin = 0;
    if(deut_count_sb7 > 0  &&  trig_05_count > 0  &&  do_lead_only == 1)    {	begin = max_05;   trig_05_count = max_05 + 1;    }
    if(deut_count_sb7 > 0  &&  trig_05_count > 0)
    {	for(int i=begin; i<trig_05_count; i++)  // trigger loop
	{   int H               = trig_05_track_num[i];
	    AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	    Float_t pt_H        = trackH->Pt();
	    Short_t charge_H    = trackH->Charge();
	    Float_t phi_H       = trackH->Phi();
	    Float_t eta_H       = trackH->Eta();
	    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }    if(phi_H <  -pio2){  phi_H = phi_H + twopi; }
	    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }    if(phi_H > 3*pio2){  phi_H = phi_H - twopi; }
	    for(int j=0; j<deut_count_sb7; j++)
	    {	int A               = deut_track_num_sb7[j];
		if(A != H)
		{   AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
		    if(trackA)
		    {	Short_t charge_sb7    = trackA->Charge();	    Float_t pt_sb7      = trackA->Pt();
			Float_t phi_sb7       = trackA->Phi();		    Float_t eta_A       = trackA->Eta();
			if(phi_sb7 <  -pio2){    phi_sb7 = phi_sb7 + twopi;   }    if(phi_sb7 <  -pio2){    phi_sb7 = phi_sb7 + twopi;   }
			if(phi_sb7 > 3*pio2){    phi_sb7 = phi_sb7 - twopi;   }    if(phi_sb7 > 3*pio2){    phi_sb7 = phi_sb7 - twopi;   }
			Float_t Sdphi       =  phi_sb7 - phi_H;
			if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }
			if(fabs(eta_H - eta_A) < 1.4  ||  eta_limiting == 0)
			{
			if     (charge_H > 0){  if     (charge_sb7 > 0){    deut_dphi_pt_pos_pos_05_sb7->Fill(pt_sb7, Sdphi);	 }
			                        else if(charge_sb7 < 0){    deut_dphi_pt_pos_neg_05_sb7->Fill(pt_sb7, Sdphi);	 }    }
			else if(charge_H < 0){  if     (charge_sb7 > 0){    deut_dphi_pt_pos_neg_05_sb7->Fill(pt_sb7, Sdphi);	 }
			                        else if(charge_sb7 < 0){    deut_dphi_pt_neg_neg_05_sb7->Fill(pt_sb7, Sdphi);	 }    }
			}
    }   }   }   }   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file






    if(run_mode == 0)
    {
    
	AliAODHandler *oh = (AliAODHandler*)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
	if (oh)
	    oh->SetFillAOD(kFALSE);

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
	    if (1) {   AliAODTZERO     *out =                           eout->GetTZEROData();        AliAODTZERO     *in =                evin->GetTZEROData();     *out = *in;                  }
	    if (1) {   TClonesArray    *out =                           eout->GetTracks();	     TClonesArray    *in =                evin->GetTracks();	    new (out) TClonesArray(*in);     }
	    if (1) {   TClonesArray    *out =                           eout->GetVertices();         TClonesArray    *in =                evin->GetVertices();      new (out) TClonesArray(*in);     }
	    if (1) {   TClonesArray    *out =                           eout->GetCaloClusters();     TClonesArray    *in =                evin->GetCaloClusters();  new (out) TClonesArray(*in);     }
	    if (1) {   AliAODCaloCells *out =                           eout->GetEMCALCells();       AliAODCaloCells *in =                evin->GetEMCALCells();    new (out) AliAODCaloCells(*in);  }
	}
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
