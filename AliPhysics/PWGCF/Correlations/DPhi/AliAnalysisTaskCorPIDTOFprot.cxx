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
#include "AliAnalysisTaskCorPIDTOFprot.h"
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


ClassImp(AliAnalysisTaskCorPIDTOFprot) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFprot::AliAnalysisTaskCorPIDTOFprot() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),

    fHistPt(0),                    //  1
    cent_ntracks(0),               //  2
    
    m2_pt_pos_cut_T_prot(0),       // 35
    m2_pt_neg_cut_T_prot(0),       // 36

    prot_phi_pt_pos_T(0),          // 41
    prot_phi_pt_neg_T(0),          // 42

    prot_per_event(0),             // 68
 
    tof_phi_eta_pos(0),            // 72
    tof_phi_eta_neg(0),            // 73

    tpc_sector_fraction(0),        // 78
    primary_vertex_z(0),           // 79
    primary_vertex_z_cut1(0),      // 80
    primary_vertex_z_cut2(0),      // 81
	
    m2_pt_pos_fine(0),             //128
    m2_pt_neg_fine(0),             //129

    m2_pt_pos_TPC_prot_fine(0),    //158
    m2_pt_neg_TPC_prot_fine(0),    //159

    m2_pt_pos_cut_T_prot_fine(0),  //160
    m2_pt_neg_cut_T_prot_fine(0),  //161

    prot_prot_0509_same(0),
    prot_prot_0509_diff(0),
    prot_prot_1018_same(0),
    prot_prot_1018_diff(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFprot::AliAnalysisTaskCorPIDTOFprot(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),

    fHistPt(0),                    //  1
    cent_ntracks(0),               //  2
    
    m2_pt_pos_cut_T_prot(0),       // 35
    m2_pt_neg_cut_T_prot(0),       // 36

    prot_phi_pt_pos_T(0),          // 41
    prot_phi_pt_neg_T(0),          // 42

    prot_per_event(0),             // 68
 
    tof_phi_eta_pos(0),            // 72
    tof_phi_eta_neg(0),            // 73

    tpc_sector_fraction(0),        // 78
    primary_vertex_z(0),           // 79
    primary_vertex_z_cut1(0),      // 80
    primary_vertex_z_cut2(0),      // 81
	
    m2_pt_pos_fine(0),             //128
    m2_pt_neg_fine(0),             //129

    m2_pt_pos_TPC_prot_fine(0),    //158
    m2_pt_neg_TPC_prot_fine(0),    //159

    m2_pt_pos_cut_T_prot_fine(0),  //160
    m2_pt_neg_cut_T_prot_fine(0),  //161

    prot_prot_0509_same(0),
    prot_prot_0509_diff(0),
    prot_prot_1018_same(0),
    prot_prot_1018_diff(0)									   
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
AliAnalysisTaskCorPIDTOFprot::~AliAnalysisTaskCorPIDTOFprot()
{
    // destructor
    if(fAnalysisUtils) delete fAnalysisUtils;
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFprot::UserCreateOutputObjects()
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

    fHistPt                    = new TH1F("fHistPt",                    "Pt()",                         18, coarse_binning);                              //  1
    cent_ntracks               = new TH2F("cent_ntracks",               "cent_ntracks",                100,      0,    100,     100,       0,     800);   //  2

    m2_pt_pos_cut_T_prot       = new TH2F("m2_pt_pos_cut_T_prot",       "m2_pt_pos_cut_T_prot",         18,coarse_binning2,    2400,    -1.0,     7.0);   // 19
    m2_pt_neg_cut_T_prot       = new TH2F("m2_pt_neg_cut_T_prot",       "m2_pt_neg_cut_T_prot",         18,coarse_binning2,    2400,    -1.0,     7.0);   // 20
    
    prot_phi_pt_pos_T          = new TH2F("prot_phi_pt_pos_T",          "prot_phi_pt_pos_T",            18,coarse_binning2,     288,   lower,   upper);   // 25
    prot_phi_pt_neg_T          = new TH2F("prot_phi_pt_neg_T",          "prot_phi_pt_neg_T",            18,coarse_binning2,     288,   lower,   upper);   // 26

    prot_per_event             = new TH1I("prot_per_event",             "prot_per_event",               20,      0,     20);                              // 36

    tof_phi_eta_pos            = new TH2F("tof_phi_eta_pos",            "tof_phi_eta_pos",             190,  lower,  upper,      60,   -0.92,    0.92);   // 40
    tof_phi_eta_neg            = new TH2F("tof_phi_eta_neg",            "tof_phi_eta_neg",             190,  lower,  upper,      60,   -0.92,    0.92);   // 41

    tpc_sector_fraction        = new TH1F("tpc_sector_fraction",        "tpc_sector_fraction",         100,    0.0,    1.0);                              // 46
    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            400,  -20.0,   20.0);                              // 47
    primary_vertex_z_cut1      = new TH1F("primary_vertex_z_cut1",      "primary_vertex_z_cut1",       400,  -20.0,   20.0);                              // 48
    primary_vertex_z_cut2      = new TH1F("primary_vertex_z_cut2",      "primary_vertex_z_cut2",       400,  -20.0,   20.0);                              // 49

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

    m2_pt_pos_TPC_prot_fine    = new TH2F("m2_pt_pos_TPC_prot_fine",    "m2_pt_pos_TPC_prot_fine",     800,   fine_binning,    2400,    -1.0,     7.0);   // 86
    m2_pt_neg_TPC_prot_fine    = new TH2F("m2_pt_neg_TPC_prot_fine",    "m2_pt_neg_TPC_prot_fine",     800,   fine_binning,    2400,    -1.0,     7.0);   // 87

    m2_pt_pos_cut_T_prot_fine  = new TH2F("m2_pt_pos_cut_T_prot_fine",  "m2_pt_pos_cut_T_prot_fine",   800,   fine_binning,    2400,    -1.0,     7.0);   // 88
    m2_pt_neg_cut_T_prot_fine  = new TH2F("m2_pt_neg_cut_T_prot_fine",  "m2_pt_neg_cut_T_prot_fine",   800,   fine_binning,    2400,    -1.0,     7.0);   // 89

    prot_prot_0509_same        = new TH1F("prot_prot_0509_same",        "prot_prot_0509_same",     288,   lower,   upper);   //
    prot_prot_0509_diff        = new TH1F("prot_prot_0509_diff",        "prot_prot_0509_diff",     288,   lower,   upper);   //
    prot_prot_1018_same        = new TH1F("prot_prot_1018_same",        "prot_prot_1018_same",     288,   lower,   upper);   //
    prot_prot_1018_diff        = new TH1F("prot_prot_1018_diff",        "prot_prot_1018_diff",     288,   lower,   upper);   //

    
  // objects added to output file

    fOutputList->Add(prot_prot_0509_same);
    fOutputList->Add(prot_prot_0509_diff);
    fOutputList->Add(prot_prot_1018_same);
    fOutputList->Add(prot_prot_1018_diff);
    
    fOutputList->Add(fHistPt);                      //  1
    fOutputList->Add(cent_ntracks);                 //  2

    fOutputList->Add(m2_pt_pos_cut_T_prot);         // 19
    fOutputList->Add(m2_pt_neg_cut_T_prot);         // 20
    
    fOutputList->Add(prot_phi_pt_pos_T);            // 25
    fOutputList->Add(prot_phi_pt_neg_T);            // 26
    
    fOutputList->Add(prot_per_event);               // 36
    
    fOutputList->Add(tof_phi_eta_pos);              // 40
    fOutputList->Add(tof_phi_eta_neg);              // 41

    fOutputList->Add(tpc_sector_fraction);          // 46
    fOutputList->Add(primary_vertex_z);             // 47
    fOutputList->Add(primary_vertex_z_cut1);        // 48
    fOutputList->Add(primary_vertex_z_cut2);        // 49

    fOutputList->Add(m2_pt_pos_fine);               // 72
    fOutputList->Add(m2_pt_neg_fine);               // 73

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
void AliAnalysisTaskCorPIDTOFprot::UserExec(Option_t *)
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

    int hybrid_count   = 0;
	
    int prot_track_num_T[500];     int prot_count_T            = 0;
    
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
	
	Short_t charge        = track->Charge();
	
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
/*
	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < TPC_PID)
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

		    if(m2tof >= 2.5  &&  m2tof < 5.0)
		    {
			deut_track_num_W[deut_count_W] = i;
			deut_count_W++;
		    }

		    if(m2tof >= 1.7  &&  m2tof < 5.5)    {   wide_cut++;   }
		}
	    }
*/
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
/*
	    Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	    if(TMath::Abs(nSigmaTPCDeut) < TPC_PID)
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

		    if(m2tof >= 2.5  &&  m2tof < 5.0)
		    {
			deut_track_num_W[deut_count_W] = i;
			deut_count_W++;
		    }

		    if(m2tof >= 1.7  &&  m2tof < 5.5)    {   wide_cut++;   }
		}
	    }
*/

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


    prot_per_event->Fill(prot_count_T);

//    if(trig_05_count > 0)    {       deut_per_event_with_trig->Fill(deut_count_T);   }
//    trig_05_per_event->Fill(trig_05_count); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    if(prot_count_T == 2)

    // markB
    for(int i=0; i<prot_count_T-1; i++)
    {
	int A = prot_track_num_T[i];
	AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	if(trackA)
	{
	    Short_t charge_A    = trackA->Charge();
	    Float_t pt_A        = trackA->Pt();
	    Float_t phi_A       = trackA->Phi();
//	    Float_t eta_A       = trackA->Eta();

	    if(pt_A > 1.5  &&  pt_A < 2.0)
	    {
		for(int j=i+1; j<prot_count_T; j++)
		{
		    int B = prot_track_num_T[j];
	
		    AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));

		    if(trackB)
		    {
			Short_t charge_B    = trackB->Charge();
			Float_t pt_B        = trackB->Pt();
			Float_t phi_B       = trackB->Phi();
//			Float_t eta_B       = trackB->Eta();

			if(pt_B > 0.5  &&  pt_B < 0.9)
			{
			    Float_t Sdphi       =  phi_A - phi_B;
			    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }

			    Float_t dcharge = charge_A - charge_B;		if(dcharge < 0.0)  dcharge = -1.0 * dcharge;
			    
			    if(dcharge < 1){	prot_prot_0509_same->Fill(Sdphi);   }
			    else           {	prot_prot_0509_diff->Fill(Sdphi);   }
			}
			else if(pt_B > 1.0  &&  pt_B < 1.8)
			{
			    Float_t Sdphi       =  phi_A - phi_B;
			    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }

			    Float_t dcharge = charge_A - charge_B;		if(dcharge < 0.0)  dcharge = -1.0 * dcharge;
			    
			    if(dcharge < 1){	prot_prot_1018_same->Fill(Sdphi);   }
			    else           {	prot_prot_1018_diff->Fill(Sdphi);   }
			}
		    }
		}
	    }
	    else if(pt_A > 0.5  &&  pt_A < 0.9)
	    {
		for(int j=i+1; j<prot_count_T; j++)
		{
		    int B = prot_track_num_T[j];
	
		    AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));

		    if(trackB)
		    {
			Short_t charge_B    = trackB->Charge();
			Float_t pt_B        = trackB->Pt();
			Float_t phi_B       = trackB->Phi();
//			Float_t eta_B       = trackB->Eta();

			if(pt_B > 1.5  &&  pt_B < 2.0)
			{
			    Float_t Sdphi       =  phi_A - phi_B;
			    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }

			    Float_t dcharge = charge_A - charge_B;		if(dcharge < 0.0)  dcharge = -1.0 * dcharge;
			    
			    if(dcharge < 1){	prot_prot_0509_same->Fill(Sdphi);   }
			    else           {	prot_prot_0509_diff->Fill(Sdphi);   }
			}
		    }
		}
	    }
	    else if(pt_A > 1.0  &&  pt_A < 1.8)
	    {
		for(int j=i+1; j<prot_count_T; j++)
		{
		    int B = prot_track_num_T[j];
	
		    AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));

		    if(trackB)
		    {
			Short_t charge_B    = trackB->Charge();
			Float_t pt_B        = trackB->Pt();
			Float_t phi_B       = trackB->Phi();
//			Float_t eta_B       = trackB->Eta();

			if(pt_B > 1.5  &&  pt_B < 2.0)
			{
			    Float_t Sdphi       =  phi_A - phi_B;
			    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }    if(Sdphi <  -pio2){    Sdphi = Sdphi + twopi;   }
			    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }    if(Sdphi > 3*pio2){    Sdphi = Sdphi - twopi;   }

			    Float_t dcharge = charge_A - charge_B;		if(dcharge < 0.0)  dcharge = -1.0 * dcharge;
			    
			    if(dcharge < 1){	prot_prot_1018_same->Fill(Sdphi);   }
			    else           {	prot_prot_1018_diff->Fill(Sdphi);   }
			}
		    }
		}
	    }
	}
    }	


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file



}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFprot::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________



//_____________________________________________________________________________
Double_t AliAnalysisTaskCorPIDTOFprot::Beta(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFprot::tof_minus_tpion(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFprot::get_mass_squared(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFprot::get_deut_tof_pt(AliAODTrack *track)
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
Double_t AliAnalysisTaskCorPIDTOFprot::get_deut_tof_p(AliAODTrack *track)
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
