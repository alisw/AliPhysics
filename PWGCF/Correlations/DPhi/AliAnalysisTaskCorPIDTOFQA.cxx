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
    
    high_per_event_05(0),       //  3a
    high_per_event_10(0),       //  3b
    high_per_event_pos(0),      //  4
    high_per_event_neg(0),      //  5


    dphi_pt_deut_05(0),         //  6
    dphi_et_deut_05(0),         //  7
    dphi_en_deut_05(0),         //  8
    dphi_kt_deut_05(0),         //  9
    
    dphi_pt_prot_05(0),         // 10
    dphi_et_prot_05(0),         // 11
    dphi_en_prot_05(0),         // 12
    dphi_kt_prot_05(0),         // 13

    dphi_pt_hadr_05(0),         // 14
    dphi_et_hadr_05(0),         // 15
    dphi_en_hadr_05(0),         // 16
    dphi_kt_hadr_05(0),         // 17
    
//  high_pt_count_05(0),        // 18
//  high_et_count_05(0),        // 19
//  high_en_count_05(0),        // 20
//  high_kt_count_05(0),        // 21
    
    deut_phi_hist_05(0),        // 22
    prot_phi_hist_05(0),        // 23
    hadr_phi_hist_05(0),        // 24
    high_phi_hist_05(0),        // 25

    dphi_pt_deut_10(0),         // 26
    dphi_et_deut_10(0),         // 27
    dphi_en_deut_10(0),         // 28
    dphi_kt_deut_10(0),         // 29
    
    dphi_pt_prot_10(0),         // 30
    dphi_et_prot_10(0),         // 31
    dphi_en_prot_10(0),         // 32
    dphi_kt_prot_10(0),         // 33

    dphi_pt_hadr_10(0),         // 34
    dphi_et_hadr_10(0),         // 35
    dphi_en_hadr_10(0),         // 36
    dphi_kt_hadr_10(0),         // 37
    
//  high_pt_count_10(0),        // 38
//  high_et_count_10(0),        // 39
//  high_en_count_10(0),        // 40
//  high_kt_count_10(0),        // 41
    
    deut_phi_hist_10(0),        // 42
    prot_phi_hist_10(0),        // 43
    hadr_phi_hist_10(0),        // 44
    high_phi_hist_10(0),        // 45    
    
    m2_cut(0)                   // 46

    
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0),
									   
    fHistPt(0),                 //  1
    cent_ntracks(0),            //  2
    
    high_per_event_05(0),       //  3a
    high_per_event_10(0),       //  3b									       
    high_per_event_pos(0),      //  4
    high_per_event_neg(0),      //  5

    dphi_pt_deut_05(0),         //  6
    dphi_et_deut_05(0),         //  7
    dphi_en_deut_05(0),         //  8
    dphi_kt_deut_05(0),         //  9
    
    dphi_pt_prot_05(0),         // 10
    dphi_et_prot_05(0),         // 11
    dphi_en_prot_05(0),         // 12
    dphi_kt_prot_05(0),         // 13

    dphi_pt_hadr_05(0),         // 14
    dphi_et_hadr_05(0),         // 15
    dphi_en_hadr_05(0),         // 16
    dphi_kt_hadr_05(0),         // 17
    
//  high_pt_count_05(0),        // 18
//  high_et_count_05(0),        // 19
//  high_en_count_05(0),        // 20
//  high_kt_count_05(0),        // 21
    
    deut_phi_hist_05(0),        // 22
    prot_phi_hist_05(0),        // 23
    hadr_phi_hist_05(0),        // 24
    high_phi_hist_05(0),        // 25

    dphi_pt_deut_10(0),         // 26
    dphi_et_deut_10(0),         // 27
    dphi_en_deut_10(0),         // 28
    dphi_kt_deut_10(0),         // 29
    
    dphi_pt_prot_10(0),         // 30
    dphi_et_prot_10(0),         // 31
    dphi_en_prot_10(0),         // 32
    dphi_kt_prot_10(0),         // 33

    dphi_pt_hadr_10(0),         // 34
    dphi_et_hadr_10(0),         // 35
    dphi_en_hadr_10(0),         // 36
    dphi_kt_hadr_10(0),         // 37
    
//  high_pt_count_10(0),        // 38
//  high_et_count_10(0),        // 39
//  high_en_count_10(0),        // 40
//  high_kt_count_10(0),        // 41
    
    deut_phi_hist_10(0),        // 42
    prot_phi_hist_10(0),        // 43
    hadr_phi_hist_10(0),        // 44
    high_phi_hist_10(0),        // 45    
    
    m2_cut(0)                   // 46

									       
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


    
    prot_curves[0][0][0] = 0.763625;     // pos prot mean curve
    prot_curves[0][0][1] = 0.0194772;
    prot_curves[0][0][2] = 0.0903159;
    prot_curves[0][0][3] =-0.000109114;
    
    prot_curves[0][1][0] =-0.129833;  // pos prot sigma curve
    prot_curves[0][1][1] = 0.0625516;
    prot_curves[0][1][2] = 0.10124;
    prot_curves[0][1][3] = 0.000109926;

    prot_curves[1][0][0] = 0.778006;     // neg prot mean curve
    prot_curves[1][0][1] = 0.0183448;
    prot_curves[1][0][2] = 0.0754383;
    prot_curves[1][0][3] =-0.00010526;

    prot_curves[1][1][0] =-0.117123;  // neg prot sigma curve
    prot_curves[1][1][1] = 0.0597632;
    prot_curves[1][1][2] = 0.0921575;
    prot_curves[1][1][3] = 0.000118412;

	
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    fHistPt               = new TH1F("fHistPt",               "Pt()",                 1000,     0,     10);                              //  1
    cent_ntracks          = new TH2F("cent_ntracks",          "cent_ntracks",          100,     0,    100,     100,       0,     800);   //  2


    
    high_per_event_05     = new TH1I("high_per_event_05",     "high_per_event_05",      50,     0,     50);                              //  3a
    high_per_event_10     = new TH1I("high_per_event_10",     "high_per_event_10",      50,     0,     50);                              //  3b
    high_per_event_pos    = new TH1I("high_per_event_pos",    "high_per_event_pos",     50,     0,     50);                              //  4
    high_per_event_neg    = new TH1I("high_per_event_neg",    "high_per_event_neg",     50,     0,     50);                              //  5



    dphi_pt_deut_05       = new TH2F("dphi_pt_deut_05",       "dphi_pt_deut_05",       160,  1.00,   5.00,     300, -1.6708,  4.8124);   //  6
    dphi_et_deut_05       = new TH2F("dphi_et_deut_05",       "dphi_et_deut_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   //  7
    dphi_en_deut_05       = new TH2F("dphi_en_deut_05",       "dphi_en_deut_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   //  8
    dphi_kt_deut_05       = new TH2F("dphi_kt_deut_05",       "dphi_kt_deut_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   //  9

    dphi_pt_prot_05       = new TH2F("dphi_pt_prot_05",       "dphi_pt_prot_05",       160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 10
    dphi_et_prot_05       = new TH2F("dphi_et_prot_05",       "dphi_et_prot_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 11
    dphi_en_prot_05       = new TH2F("dphi_en_prot_05",       "dphi_en_prot_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 12
    dphi_kt_prot_05       = new TH2F("dphi_kt_prot_05",       "dphi_kt_prot_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 13

    dphi_pt_hadr_05       = new TH2F("dphi_pt_hadr_05",       "dphi_pt_hadr_05",       160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 14
    dphi_et_hadr_05       = new TH2F("dphi_et_hadr_05",       "dphi_et_hadr_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 15
    dphi_en_hadr_05       = new TH2F("dphi_en_hadr_05",       "dphi_en_hadr_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 16
    dphi_kt_hadr_05       = new TH2F("dphi_kt_hadr_05",       "dphi_kt_hadr_05",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 17

//  high_pt_count_05      = new TH1F("high_pt_count_05",      "high_pt_count_05",      160,  1.00,   5.00);                              // 18
//  high_et_count_05      = new TH1F("high_et_count_05",      "high_et_count_05",      160,  1.00,   4.78);                              // 19
//  high_en_count_05      = new TH1F("high_en_count_05",      "high_en_count_05",      160,  1.00,   4.78);                              // 20
//  high_kt_count_05      = new TH1F("high_kt_count_05",      "high_kt_count_05",      160,  1.00,   4.78);                              // 21

    deut_phi_hist_05      = new TH1F("deut_phi_hist_05",      "deut_phi_hist_05",      480, -0.10, 6.383185);                            // 22
    prot_phi_hist_05      = new TH1F("prot_phi_hist_05",      "prot_phi_hist_05",      480, -0.10, 6.383185);                            // 23
    hadr_phi_hist_05      = new TH1F("hadr_phi_hist_05",      "hadr_phi_hist_05",      480, -0.10, 6.383185);                            // 24
    high_phi_hist_05      = new TH1F("high_phi_hist_05",      "high_phi_hist_05",      480, -0.10, 6.383185);                            // 25
    

    dphi_pt_deut_10       = new TH2F("dphi_pt_deut_10",       "dphi_pt_deut_10",       160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 26
    dphi_et_deut_10       = new TH2F("dphi_et_deut_10",       "dphi_et_deut_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 27
    dphi_en_deut_10       = new TH2F("dphi_en_deut_10",       "dphi_en_deut_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 28
    dphi_kt_deut_10       = new TH2F("dphi_kt_deut_10",       "dphi_kt_deut_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 29

    dphi_pt_prot_10       = new TH2F("dphi_pt_prot_10",       "dphi_pt_prot_10",       160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 30
    dphi_et_prot_10       = new TH2F("dphi_et_prot_10",       "dphi_et_prot_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 31
    dphi_en_prot_10       = new TH2F("dphi_en_prot_10",       "dphi_en_prot_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 32
    dphi_kt_prot_10       = new TH2F("dphi_kt_prot_10",       "dphi_kt_prot_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 33

    dphi_pt_hadr_10       = new TH2F("dphi_pt_hadr_10",       "dphi_pt_hadr_10",       160,  1.00,   5.00,     300, -1.6708,  4.8124);   // 34
    dphi_et_hadr_10       = new TH2F("dphi_et_hadr_10",       "dphi_et_hadr_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 35
    dphi_en_hadr_10       = new TH2F("dphi_en_hadr_10",       "dphi_en_hadr_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 36
    dphi_kt_hadr_10       = new TH2F("dphi_kt_hadr_10",       "dphi_kt_hadr_10",       189,  1.00,   4.78,     300, -1.6708,  4.8124);   // 37

//  high_pt_count_10      = new TH1F("high_pt_count_10",      "high_pt_count_10",      160,  1.00,   5.00);                              // 38
//  high_et_count_10      = new TH1F("high_et_count_10",      "high_et_count_10",      160,  1.00,   4.78);                              // 39
//  high_en_count_10      = new TH1F("high_en_count_10",      "high_en_count_10",      160,  1.00,   4.78);                              // 40
//  high_kt_count_10      = new TH1F("high_kt_count_10",      "high_kt_count_10",      160,  1.00,   4.78);                              // 41

    deut_phi_hist_10      = new TH1F("deut_phi_hist_10",      "deut_phi_hist_10",      480, -0.10, 6.383185);                            // 42
    prot_phi_hist_10      = new TH1F("prot_phi_hist_10",      "prot_phi_hist_10",      480, -0.10, 6.383185);                            // 43
    hadr_phi_hist_10      = new TH1F("hadr_phi_hist_10",      "hadr_phi_hist_10",      480, -0.10, 6.383185);                            // 44
    high_phi_hist_10      = new TH1F("high_phi_hist_10",      "high_phi_hist_10",      480, -0.10, 6.383185);                            // 45
    
    m2_cut                = new TH2F("m2_cut",                "m2_cut",                200,  1.00, 5.00,      400,    -1.00,    7.00);   // 46


    
    
    fOutputList->Add(fHistPt);                //  1                                                             // objects added to output file
    fOutputList->Add(cent_ntracks);           //  2
    

    fOutputList->Add(high_per_event_05);      //  3a
    fOutputList->Add(high_per_event_10);      //  3b
    fOutputList->Add(high_per_event_pos);     //  4
    fOutputList->Add(high_per_event_neg);     //  5

    fOutputList->Add(dphi_pt_deut_05);        //  6
    fOutputList->Add(dphi_et_deut_05);        //  7
    fOutputList->Add(dphi_en_deut_05);        //  8
    fOutputList->Add(dphi_kt_deut_05);        //  9
    
    fOutputList->Add(dphi_pt_prot_05);        // 10
    fOutputList->Add(dphi_et_prot_05);        // 11
    fOutputList->Add(dphi_en_prot_05);        // 12
    fOutputList->Add(dphi_kt_prot_05);        // 13

    fOutputList->Add(dphi_pt_hadr_05);        // 14
    fOutputList->Add(dphi_et_hadr_05);        // 15
    fOutputList->Add(dphi_en_hadr_05);        // 16
    fOutputList->Add(dphi_kt_hadr_05);        // 17
    
//  fOutputList->Add(high_pt_count_05);       // 18
//  fOutputList->Add(high_et_count_05);       // 19
//  fOutputList->Add(high_en_count_05);       // 20
//  fOutputList->Add(high_kt_count_05);       // 21
    
    fOutputList->Add(deut_phi_hist_05);       // 22
    fOutputList->Add(prot_phi_hist_05);       // 23
    fOutputList->Add(hadr_phi_hist_05);       // 24
    fOutputList->Add(high_phi_hist_05);       // 25

    fOutputList->Add(dphi_pt_deut_10);        // 26
    fOutputList->Add(dphi_et_deut_10);        // 27
    fOutputList->Add(dphi_en_deut_10);        // 28
    fOutputList->Add(dphi_kt_deut_10);        // 29
    
    fOutputList->Add(dphi_pt_prot_10);        // 30
    fOutputList->Add(dphi_et_prot_10);        // 31
    fOutputList->Add(dphi_en_prot_10);        // 32
    fOutputList->Add(dphi_kt_prot_10);        // 33

    fOutputList->Add(dphi_pt_hadr_10);        // 34
    fOutputList->Add(dphi_et_hadr_10);        // 35
    fOutputList->Add(dphi_en_hadr_10);        // 36
    fOutputList->Add(dphi_kt_hadr_10);        // 37
    
//  fOutputList->Add(high_pt_count_10);       // 38
//  fOutputList->Add(high_et_count_10);       // 39
//  fOutputList->Add(high_en_count_10);       // 40
//  fOutputList->Add(high_kt_count_10);       // 41
    
    fOutputList->Add(deut_phi_hist_10);       // 42
    fOutputList->Add(prot_phi_hist_10);       // 43
    fOutputList->Add(hadr_phi_hist_10);       // 44
    fOutputList->Add(high_phi_hist_10);       // 45

    fOutputList->Add(m2_cut);                 // 46

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

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    Int_t iTracks(fAOD->GetNumberOfTracks());
    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////



     
    int NHigh_05             = 0;
    int NHigh_pos            = 0;
    int NHigh_neg            = 0;
    int high_flag_05         = 0;

    int high_track_num_05[100];

    int NHigh_10             = 0;
    int high_flag_10         = 0;

    int high_track_num_10[100];
    

	
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	track_node *new_track = new track_node;
	new_track->track_num = i;
	new_track->next = NULL;
	conductor->next = new_track;
	conductor = conductor->next;

	Double_t eta = track->Eta();                     	if(fabs(eta > 0.8))   continue;
	
	float    pt             = track->Pt();
	Short_t  charge         = track->Charge();


	fHistPt->Fill(pt);
 

	    
	if(charge > 0)
	{
	    if(pt > 5.0)
	    {
		high_track_num_05[NHigh_05] = i;
		NHigh_05++;
		NHigh_pos++;
		high_flag_05 = 1;
	    }
	    if(pt > 10.0)
	    {
		high_track_num_10[NHigh_10] = i;
		NHigh_10++;
		high_flag_10 = 1;
	    }
	    
	}
	    
	else if(charge < 0)
	{
	    if(pt > 5.0)
	    {
		high_track_num_05[NHigh_05] = i;
		NHigh_05++;
		NHigh_neg++;
		high_flag_05 = 1;
	    }
	    if(pt > 10.0)
	    {
		high_track_num_10[NHigh_10] = i;
		NHigh_10++;
		high_flag_10 = 1;
	    }
	}
    }        //   end of track loop


    high_per_event_05->Fill(NHigh_05);
    high_per_event_10->Fill(NHigh_10);
    high_per_event_pos->Fill(NHigh_pos);
    high_per_event_neg->Fill(NHigh_neg);



    
    for(int ii=0; ii < high_flag_05; ii++)
    {
	int j = high_track_num_05[ii];

	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t high_phi = track->Phi();
	    
	    high_phi_hist_05->Fill(high_phi);
	    
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;

		    short match = 0;
		    for(int kk=0; kk < high_flag_05; kk++)
		    {
			if(k == high_track_num_05[kk])   match = 1;
		    }

		    
		    if(match == 0)  // if the k-th track is not the j-th track which is the trigger particle of interest
		    {
			AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
			if(!track){}                                                               // if we failed, skip this track
			else
			{
			    Double_t pt             = track->Pt();

			    
			    if(pt > 1.1  &&  pt < 4.4)
			    {
				short isGlobal    = 1;
				short isModGlobal = 0;  // is it suitable for deuteron candidate
	
				if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal))){   isGlobal = 0;   }
	
				Double_t nsigmaTPC = 999.0;
				Double_t nsigmaTOF = 999.0;

				if(isGlobal == 1)  
				{
				    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
				    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
				    Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);  
				    Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);

				    if(tpcIsOk  && tofIsOk){   isModGlobal = 1;   }
				    else                   {   isModGlobal = 0;   }
				}

				if(isModGlobal == 1)
				{
				    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
				    if(nsigmaTPC <= 2.0)
				    {					
					Double_t deut_mean      = 0.00;
					Double_t deut_sigma     = 0.00;
					Double_t cut_width      = 1.00;

					Double_t mom   = track->P();
					Double_t m2tof = 0.00;
					m2tof = get_mass_squared(track);
					
					for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
					deut_mean = fit_deut_curve->Eval(mom);
					for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
					deut_sigma = fit_deut_curve->Eval(mom);
					
					if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
					{
//					    Double_t mom            = track->P();
					    Double_t phi            = track->Phi();

					    Double_t et = 0.00;
					    Double_t en = 0.00;
					    Double_t kt = 0.00;
	    
					    et = sqrt(pow(pt ,2) + pow(1.87561,2));
					    en = sqrt(pow(mom,2) + pow(1.87561,2));
					    kt = sqrt(pow(pt ,2) + pow(1.87561,2)) - 1.87561;
				
					    Double_t dphi = 0.00;
					    dphi = high_phi - phi;
				
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
					    dphi_pt_deut_05->Fill(pt, dphi);
					    dphi_et_deut_05->Fill(et, dphi);
					    dphi_en_deut_05->Fill(en, dphi);
					    dphi_kt_deut_05->Fill(kt, dphi);
					    
					    deut_phi_hist_05->Fill(phi);

					    m2_cut->Fill(pt, m2tof);
					}
				    }
				    statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
				    if(nsigmaTPC <= 2.0)
				    {					
					Double_t prot_mean      = 0.00;
					Double_t prot_sigma     = 0.00;
					Double_t cut_width      = 1.00;
				    
					Double_t mom            = track->P();
					Double_t m2tof          = 0.00;
					
					m2tof = get_mass_squared(track);
					
					for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
					prot_mean = fit_prot_curve->Eval(mom);
					for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
					prot_sigma = fit_prot_curve->Eval(mom);
					
					if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
					{
					    Double_t phi            = track->Phi();
					    
					    Double_t et = 0.00;
					    Double_t en = 0.00;
					    Double_t kt = 0.00;
					    
					    et = sqrt(pow(pt ,2) + pow(0.93827,2));
					    en = sqrt(pow(mom,2) + pow(0.93827,2));
					    kt = sqrt(pow(pt ,2) + pow(0.93827,2)) - 0.93827;
					    
					    Double_t dphi = 0.00;
					    dphi = high_phi - phi;
					    
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    
					    dphi_pt_prot_05->Fill(pt, dphi);
					    dphi_et_prot_05->Fill(et, dphi);
					    dphi_en_prot_05->Fill(en, dphi);
					    dphi_kt_prot_05->Fill(kt, dphi);
					    
					    prot_phi_hist_05->Fill(phi);

					    m2_cut->Fill(pt, m2tof);
					}
				    }

				    Double_t phi = track->Phi();
				    Double_t mom = track->P();
				    Double_t et  = 0.00;
				    Double_t en  = 0.00;
				    Double_t kt  = 0.00;

				    et = sqrt(pow(pt ,2) + pow(0.29026,2));
				    en = sqrt(pow(mom,2) + pow(0.29026,2));
				    kt = sqrt(pow(pt ,2) + pow(0.29026,2)) - 0.29026;
					    
				    Double_t dphi = 0.00;
				    dphi = high_phi - phi;
					    
				    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    
				    dphi_pt_hadr_05->Fill(pt, dphi);
				    dphi_et_hadr_05->Fill(et, dphi);
				    dphi_en_hadr_05->Fill(en, dphi);
				    dphi_kt_hadr_05->Fill(kt, dphi);
				    
				    hadr_phi_hist_05->Fill(phi);
				}
			    }
			}
		    }
		}
	    
	    
		int k = conductor->track_num;
		conductor = conductor->next;
	
		short match = 0;
		for(int kk=0; kk < high_flag_05; kk++)
		{
		    if(k == high_track_num_05[kk])   match = 1;
		}

		    
		if(match == 0)  // if the k-th track is not the j-th track which is the trigger particle of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t pt             = track->Pt();

			    
			if(pt > 1.1  &&  pt < 4.4)
			{
			    short isGlobal    = 1;
			    short isModGlobal = 0;  // is it suitable for deuteron candidate
	
			    if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal))){   isGlobal = 0;   }
	
			    Double_t nsigmaTPC = 999.0;
			    Double_t nsigmaTOF = 999.0;

			    if(isGlobal == 1)  
			    {
				AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
				AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
				Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);  
				Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
				
				if(tpcIsOk  && tofIsOk){   isModGlobal = 1;   }
				else                   {   isModGlobal = 0;   }
			    }

			    if(isModGlobal == 1)
			    {
				AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
				if(nsigmaTPC <= 2.0)
				{					
				    Double_t deut_mean      = 0.0;
				    Double_t deut_sigma     = 0.0;
				    Double_t cut_width      = 1.0;
				    
				    Double_t mom   = track->P();
				    Double_t m2tof = 0.00;
				    m2tof = get_mass_squared(track);
					
				    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
				    deut_mean = fit_deut_curve->Eval(mom);
				    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
				    deut_sigma = fit_deut_curve->Eval(mom);
					
				    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
				    {
					Double_t phi            = track->Phi();

					Double_t et = 0.00;
					Double_t en = 0.00;
					Double_t kt = 0.00;
	    
					et = sqrt(pow(pt ,2) + pow(1.87561,2));
					en = sqrt(pow(mom,2) + pow(1.87561,2));
					kt = sqrt(pow(pt ,2) + pow(1.87561,2)) - 1.87561;
				
					Double_t dphi = 0.00;
					dphi = high_phi - phi;
				
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
					dphi_pt_deut_05->Fill(pt, dphi);
					dphi_et_deut_05->Fill(et, dphi);
					dphi_en_deut_05->Fill(en, dphi);
					dphi_kt_deut_05->Fill(kt, dphi);

					deut_phi_hist_05->Fill(phi);

					m2_cut->Fill(pt, m2tof);
				    }
				}

				statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
				if(nsigmaTPC <= 2.0)
				{					
				    Double_t prot_mean      = 0.0;
				    Double_t prot_sigma     = 0.0;
				    Double_t cut_width      = 1.0;
				    
				    Double_t mom   = track->P();
				    Double_t m2tof = 0.00;
				    m2tof = get_mass_squared(track);
					
				    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
				    prot_mean = fit_prot_curve->Eval(mom);
				    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
				    prot_sigma = fit_prot_curve->Eval(mom);
					
				    if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
				    {
					Double_t phi            = track->Phi();

					Double_t et = 0.00;
					Double_t en = 0.00;
					Double_t kt = 0.00;
	
					et = sqrt(pow(pt ,2) + pow(0.93827,2));
					en = sqrt(pow(mom,2) + pow(0.93827,2));
					kt = sqrt(pow(pt ,2) + pow(0.93827,2)) - 0.93827;
				
					Double_t dphi = 0.00;
					dphi = high_phi - phi;
				
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
					dphi_pt_prot_05->Fill(pt, dphi);
					dphi_et_prot_05->Fill(et, dphi);
					dphi_en_prot_05->Fill(en, dphi);
					dphi_kt_prot_05->Fill(kt, dphi);

					prot_phi_hist_05->Fill(phi);

					m2_cut->Fill(pt, m2tof);
				    }
				}

				
				Double_t phi = track->Phi();
				Double_t mom = track->P();
				Double_t et = 0.00;
				Double_t en = 0.00;
				Double_t kt = 0.00;

				et = sqrt(pow(pt ,2) + pow(0.29026,2));
				en = sqrt(pow(mom,2) + pow(0.29026,2));
				kt = sqrt(pow(pt ,2) + pow(0.29026,2)) - 0.29026;
					    
				Double_t dphi = 0.00;
				dphi = high_phi - phi;
					    
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    
				dphi_pt_hadr_05->Fill(pt, dphi);
				dphi_et_hadr_05->Fill(et, dphi);
				dphi_en_hadr_05->Fill(en, dphi);
				dphi_kt_hadr_05->Fill(kt, dphi);
				    
				hadr_phi_hist_05->Fill(phi);
			    }
			}
		    }
		}
	    }
	}
    }  // end of high-pt correlation PRIMARY LOOP

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int ii=0; ii < high_flag_10; ii++)
    {
	int j = high_track_num_10[ii];

	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
	if(!track) {}
	else
	{
	    Double_t high_phi = track->Phi();
	    
	    high_phi_hist_10->Fill(high_phi);
	    
	    conductor = root->next;
	    if (conductor != 0)
	    {
		while(conductor->next != 0)
		{
		    int k = conductor->track_num;
		    conductor = conductor->next;

		    short match = 0;
		    for(int kk=0; kk < high_flag_10; kk++)
		    {
			if(k == high_track_num_10[kk])   match = 1;
		    }

		    
		    if(match == 0)  // if the k-th track is not the j-th track which is the trigger particle of interest
		    {
			AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
			if(!track){}                                                               // if we failed, skip this track
			else
			{
			    Double_t pt             = track->Pt();

			    
			    if(pt > 1.1  &&  pt < 4.4)
			    {
				short isGlobal    = 1;
				short isModGlobal = 0;  // is it suitable for deuteron candidate
	
				if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal))){   isGlobal = 0;   }
	
				Double_t nsigmaTPC = 999.0;
				Double_t nsigmaTOF = 999.0;

				if(isGlobal == 1)  
				{
				    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
				    AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
				    Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);  
				    Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);

				    if(tpcIsOk  && tofIsOk){   isModGlobal = 1;   }
				    else                   {   isModGlobal = 0;   }
				}

				if(isModGlobal == 1)
				{
				    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
				    if(nsigmaTPC <= 2.0)
				    {					
					Double_t deut_mean      = 0.00;
					Double_t deut_sigma     = 0.00;
					Double_t cut_width      = 1.00;

					Double_t mom   = track->P();
					Double_t m2tof = 0.00;
					m2tof = get_mass_squared(track);
					
					for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
					deut_mean = fit_deut_curve->Eval(mom);
					for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
					deut_sigma = fit_deut_curve->Eval(mom);
					
					if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
					{
//					    Double_t mom            = track->P();
					    Double_t phi            = track->Phi();

					    Double_t et = 0.00;
					    Double_t en = 0.00;
					    Double_t kt = 0.00;
	    
					    et = sqrt(pow(pt ,2) + pow(1.87561,2));
					    en = sqrt(pow(mom,2) + pow(1.87561,2));
					    kt = sqrt(pow(pt ,2) + pow(1.87561,2)) - 1.87561;
				
					    Double_t dphi = 0.00;
					    dphi = high_phi - phi;
				
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
					    dphi_pt_deut_10->Fill(pt, dphi);
					    dphi_et_deut_10->Fill(et, dphi);
					    dphi_en_deut_10->Fill(en, dphi);
					    dphi_kt_deut_10->Fill(kt, dphi);
					    
					    deut_phi_hist_10->Fill(phi);
					}
				    }
				    statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
				    if(nsigmaTPC <= 2.0)
				    {					
					Double_t prot_mean      = 0.00;
					Double_t prot_sigma     = 0.00;
					Double_t cut_width      = 1.00;
				    
					Double_t mom            = track->P();
					Double_t m2tof          = 0.00;
					
					m2tof = get_mass_squared(track);
					
					for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
					prot_mean = fit_prot_curve->Eval(mom);
					for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
					prot_sigma = fit_prot_curve->Eval(mom);
					
					if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
					{
					    Double_t phi            = track->Phi();
					    
					    Double_t et = 0.00;
					    Double_t en = 0.00;
					    Double_t kt = 0.00;
					    
					    et = sqrt(pow(pt ,2) + pow(0.93827,2));
					    en = sqrt(pow(mom,2) + pow(0.93827,2));
					    kt = sqrt(pow(pt ,2) + pow(0.93827,2)) - 0.93827;
					    
					    Double_t dphi = 0.00;
					    dphi = high_phi - phi;
					    
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    
					    dphi_pt_prot_10->Fill(pt, dphi);
					    dphi_et_prot_10->Fill(et, dphi);
					    dphi_en_prot_10->Fill(en, dphi);
					    dphi_kt_prot_10->Fill(kt, dphi);
					    
					    prot_phi_hist_10->Fill(phi);
					}
				    }

				    Double_t phi            = track->Phi();
				    Double_t mom = track->P();
				    Double_t et = 0.00;
				    Double_t en = 0.00;
				    Double_t kt = 0.00;

				    et = sqrt(pow(pt ,2) + pow(0.29026,2));
				    en = sqrt(pow(mom,2) + pow(0.29026,2));
				    kt = sqrt(pow(pt ,2) + pow(0.29026,2)) - 0.29026;
					    
				    Double_t dphi = 0.00;
				    dphi = high_phi - phi;
					    
				    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				    if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				    if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    
				    dphi_pt_hadr_10->Fill(pt, dphi);
				    dphi_et_hadr_10->Fill(et, dphi);
				    dphi_en_hadr_10->Fill(en, dphi);
				    dphi_kt_hadr_10->Fill(kt, dphi);
				    
				    hadr_phi_hist_10->Fill(phi);
				}
			    }
			}
		    }
		}
	    
	    
		int k = conductor->track_num;
		conductor = conductor->next;
	
		short match = 0;
		for(int kk=0; kk < high_flag_10; kk++)
		{
		    if(k == high_track_num_10[kk])   match = 1;
		}

		    
		if(match == 0)  // if the k-th track is not the j-th track which is the trigger particle of interest
		{
		    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(k));         // get a track (type AliAODTrack) from the event
		    if(!track){}                                                               // if we failed, skip this track
		    else
		    {
			Double_t pt             = track->Pt();

			    
			if(pt > 1.1  &&  pt < 4.4)
			{
			    short isGlobal    = 1;
			    short isModGlobal = 0;  // is it suitable for deuteron candidate
	
			    if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal))){   isGlobal = 0;   }
	
			    Double_t nsigmaTPC = 999.0;
			    Double_t nsigmaTOF = 999.0;

			    if(isGlobal == 1)  
			    {
				AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
				AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
				Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);  
				Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
				
				if(tpcIsOk  && tofIsOk){   isModGlobal = 1;   }
				else                   {   isModGlobal = 0;   }
			    }

			    if(isModGlobal == 1)
			    {
				AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 5, nsigmaTPC);
				if(nsigmaTPC <= 2.0)
				{					
				    Double_t deut_mean      = 0.0;
				    Double_t deut_sigma     = 0.0;
				    Double_t cut_width      = 1.0;
				    
				    Double_t mom   = track->P();
				    Double_t m2tof = 0.00;
				    m2tof = get_mass_squared(track);
					
				    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
				    deut_mean = fit_deut_curve->Eval(mom);
				    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
				    deut_sigma = fit_deut_curve->Eval(mom);
					
				    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
				    {
					Double_t phi            = track->Phi();

					Double_t et = 0.00;
					Double_t en = 0.00;
					Double_t kt = 0.00;
	    
					et = sqrt(pow(pt ,2) + pow(1.87561,2));
					en = sqrt(pow(mom,2) + pow(1.87561,2));
					kt = sqrt(pow(pt ,2) + pow(1.87561,2)) - 1.87561;
				
					Double_t dphi = 0.00;
					dphi = high_phi - phi;
				
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
					dphi_pt_deut_10->Fill(pt, dphi);
					dphi_et_deut_10->Fill(et, dphi);
					dphi_en_deut_10->Fill(en, dphi);
					dphi_kt_deut_10->Fill(kt, dphi);

					deut_phi_hist_10->Fill(phi);
				    }
				}

				statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
				if(nsigmaTPC <= 2.0)
				{					
				    Double_t prot_mean      = 0.0;
				    Double_t prot_sigma     = 0.0;
				    Double_t cut_width      = 1.0;
				    
				    Double_t mom   = track->P();
				    Double_t m2tof = 0.00;
				    m2tof = get_mass_squared(track);
					
				    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
				    prot_mean = fit_prot_curve->Eval(mom);
				    for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
				    prot_sigma = fit_prot_curve->Eval(mom);
					
				    if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
				    {
					Double_t phi            = track->Phi();

					Double_t et = 0.00;
					Double_t en = 0.00;
					Double_t kt = 0.00;
	
					et = sqrt(pow(pt ,2) + pow(0.93827,2));
					en = sqrt(pow(mom,2) + pow(0.93827,2));
					kt = sqrt(pow(pt ,2) + pow(0.93827,2)) - 0.93827;
				
					Double_t dphi = 0.00;
					dphi = high_phi - phi;
				
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				
					dphi_pt_prot_10->Fill(pt, dphi);
					dphi_et_prot_10->Fill(et, dphi);
					dphi_en_prot_10->Fill(en, dphi);
					dphi_kt_prot_10->Fill(kt, dphi);

					prot_phi_hist_10->Fill(phi);
				    }
				}

				
				Double_t phi = track->Phi();
				Double_t mom = track->P();
				Double_t et = 0.00;
				Double_t en = 0.00;
				Double_t kt = 0.00;

				et = sqrt(pow(pt ,2) + pow(0.29026,2));
				en = sqrt(pow(mom,2) + pow(0.29026,2));
				kt = sqrt(pow(pt ,2) + pow(0.29026,2)) - 0.29026;
					    
				Double_t dphi = 0.00;
				dphi = high_phi - phi;
					    
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi <  -TMath::PiOver2())    dphi = dphi + TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
				if(dphi > 3*TMath::PiOver2())    dphi = dphi - TMath::TwoPi();
					    
				dphi_pt_hadr_10->Fill(pt, dphi);
				dphi_et_hadr_10->Fill(et, dphi);
				dphi_en_hadr_10->Fill(en, dphi);
				dphi_kt_hadr_10->Fill(kt, dphi);
				    
				hadr_phi_hist_10->Fill(phi);
			    }
			}
		    }
		}
	    }
	}
    }
    
                                                     // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file

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
    Double_t m2      = 0.0;

    Double_t c2 = 29.9792458;
    
    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);
    return m2;

}
