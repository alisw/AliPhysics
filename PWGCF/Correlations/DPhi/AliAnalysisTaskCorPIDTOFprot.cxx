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

/* AliAnaysisTaskCorPIDTOFprot
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
#include "AliAnalysisTaskCorPIDTOFprot.h"
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
// using namespace BSchaefer_devel;




ClassImp(AliAnalysisTaskCorPIDTOFprot) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFprot::AliAnalysisTaskCorPIDTOFprot() : AliAnalysisTaskSE(), 
fAOD(0), fOutputList(0), fPIDResponse(0),

    fHistPt(0),                   //  1
    cent_ntracks(0),              //  2
    
    m2_pos(0),                    //  3
    m2_neg(0),                    //  4
    beta_p_pos(0),                //  5
    beta_p_neg(0),                //  6

//  deltat_p_pos(0),              //  7
//  deltat_p_neg(0),              //  8

    m2_pos_cut(0),                //  9
    m2_neg_cut(0),                // 10
    beta_p_pos_cut(0),            // 11
    beta_p_neg_cut(0),            // 12
//  deltat_p_pos_cut(0),          // 13
//  deltat_p_neg_cut(0),          // 14

    prot_per_event(0),            // 15
    prot_per_event_pos(0),        // 16
    prot_per_event_neg(0),        // 17
    
    m2_pos_cut_T(0),              // 18
    m2_neg_cut_T(0),              // 19

    prot_phi_hist(0),             // 20

    pt_mom(0),                    // 21
    proton_pt_q_pos_pos(0),       // 22
    proton_pt_q_pos_neg(0),       // 23
    proton_pt_q_neg_neg(0),       // 24
    di_proton_phi(0),             // 25
    hi_05_phi(0),                 // 26
    hi_08_phi(0),                 // 27
    di_prot_q_dphi_pos_pos_05(0), // 28
    di_prot_q_dphi_pos_neg_05(0), // 29
    di_prot_q_dphi_neg_neg_05(0), // 30

    di_prot_q_dphi_pos_pos_08(0), // 31
    di_prot_q_dphi_pos_neg_08(0), // 32
    di_prot_q_dphi_neg_neg_08(0)  // 33

//    track_cor_radius(0),          // 34
//    track_cor_radius_cut(0)       // 35
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFprot::AliAnalysisTaskCorPIDTOFprot(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0),
									   

    fHistPt(0),                   //  1
    cent_ntracks(0),              //  2
    
    m2_pos(0),                    //  3
    m2_neg(0),                    //  4
    beta_p_pos(0),                //  5
    beta_p_neg(0),                //  6

//  deltat_p_pos(0),              //  7
//  deltat_p_neg(0),              //  8

    m2_pos_cut(0),                //  9
    m2_neg_cut(0),                // 10
    beta_p_pos_cut(0),            // 11
    beta_p_neg_cut(0),            // 12
//  deltat_p_pos_cut(0),          // 13
//  deltat_p_neg_cut(0),          // 14

    prot_per_event(0),            // 15
    prot_per_event_pos(0),        // 16
    prot_per_event_neg(0),        // 17
    
    m2_pos_cut_T(0),              // 18
    m2_neg_cut_T(0),              // 19

    prot_phi_hist(0),             // 20

    pt_mom(0),                    // 21
    proton_pt_q_pos_pos(0),       // 22
    proton_pt_q_pos_neg(0),       // 23
    proton_pt_q_neg_neg(0),       // 24
    di_proton_phi(0),             // 25
    hi_05_phi(0),                 // 26
    hi_08_phi(0),                 // 27
    di_prot_q_dphi_pos_pos_05(0), // 28
    di_prot_q_dphi_pos_neg_05(0), // 29
    di_prot_q_dphi_neg_neg_05(0), // 30

    di_prot_q_dphi_pos_pos_08(0), // 31
    di_prot_q_dphi_pos_neg_08(0), // 32
    di_prot_q_dphi_neg_neg_08(0)  // 33
									       
//    track_cor_radius(0),          // 34
//    track_cor_radius_cut(0)       // 35
									       
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFprot::~AliAnalysisTaskCorPIDTOFprot()
{
    // destructor
    if(fOutputList)
    {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFprot::UserCreateOutputObjects()
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



    
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    fHistPt               = new TH1F("fHistPt",               "Pt()",                 1000,         0,      10);                              //  1
    cent_ntracks          = new TH2F("cent_ntracks",          "cent_ntracks",          100,         0,     100,     100,       0,     800);   //  2
    
    m2_pos                = new TH2F("m2_pos",                "m2_pos",                320,       0.0,     8.0,    2400,    -1.0,     7.0);   //  3
    m2_neg                = new TH2F("m2_neg",                "m2_neg",                320,       0.0,     8.0,    2400,    -1.0,     7.0);   //  4
    beta_p_pos            = new TH2F("beta_p_pos",            "beta_p_pos",            780,       0.0,     8.0,    3000,     0.1,     1.1);   //  5
    beta_p_neg            = new TH2F("beta_p_neg",            "beta_p_neg",            780,       0.0,     8.0,    3000,     0.1,     1.1);   //  6
//  deltat_p_pos          = new TH2F("deltat_p_pos",          "deltat_p_pos",        27000,       0.3,     3.0,    1000,    -1.0,     9.0);   //  7
//  deltat_p_neg          = new TH2F("deltat_p_neg",          "deltat_p_neg",        27000,       0.3,     3.0,    1000,    -1.0,     9.0);   //  8

    m2_pos_cut            = new TH2F("m2_pos_cut",            "m2_pos_cut",            320,       0.0,     8.0,    2400,    -1.0,     7.0);   //  9
    m2_neg_cut            = new TH2F("m2_neg_cut",            "m2_neg_cut",            320,       0.0,     8.0,    2400,    -1.0,     7.0);   // 10
    beta_p_pos_cut        = new TH2F("beta_p_pos_cut",        "beta_p_pos_cut",        780,       0.0,     8.0,    3000,     0.1,     1.1);   // 11
    beta_p_neg_cut        = new TH2F("beta_p_neg_cut",        "beta_p_neg_cut",        780,       0.0,     8.0,    3000,     0.1,     1.1);   // 12
//  deltat_p_pos_cut      = new TH2F("deltat_p_pos_cut",      "deltat_p_pos_cut",    27000,       0.3,     3.0,    1000,    -1.0,     9.0);   // 13
//  deltat_p_neg_cut      = new TH2F("deltat_p_neg_cut",      "deltat_p_neg_cut",    27000,       0.3,     3.0,    1000,    -1.0,     9.0);   // 14


    prot_per_event        = new TH1I("prot_per_event",        "prot_per_event",          8,         0,      8);                               // 15
    prot_per_event_pos    = new TH1I("prot_per_event_pos",    "prot_per_event_pos",      8,         0,      8);                               // 16
    prot_per_event_neg    = new TH1I("prot_per_event_neg",    "prot_per_event_neg",      8,         0,      8);                               // 17
    m2_pos_cut_T          = new TH2F("m2_pos_cut_T",          "m2_pos_cut_T",          320,       0.0,     8.0,    2400,    -1.0,     7.0);   // 18
    m2_neg_cut_T          = new TH2F("m2_neg_cut_T",          "m2_neg_cut_T",          320,       0.0,     8.0,    2400,    -1.0,     7.0);   // 19

    prot_phi_hist         = new TH1F("prot_phi_hist",         "prot_phi_hist",         300,   -1.6708,  4.8124);                              // 20

    pt_mom                = new TH2F("pt_mom",                "pt_mom",                 92,       0.4,    5.00,    152,  0.40,   8.0);        // 21
    proton_pt_q_pos_pos   = new TH2F("proton_pt_q_pos_pos",   "proton_pt_q_pos_pos",    37,       0.4,    4.10,   76,  0.00,  3.8);           // 22
    proton_pt_q_pos_neg   = new TH2F("proton_pt_q_pos_neg",   "proton_pt_q_pos_neg",    37,       0.4,    4.10,   76,  0.00,  3.8);           // 23
    proton_pt_q_neg_neg   = new TH2F("proton_pt_q_neg_neg",   "proton_pt_q_neg_neg",    37,       0.4,    4.10,   76,  0.00,  3.8);           // 24

    di_proton_phi         = new TH1F("di_proton_phi",         "di_proton_phi",         300,   -1.6708,  4.8124);                              // 25
    hi_05_phi             = new TH1F("hi_05_phi",             "hi_05_phi",             300,   -1.6708,  4.8124);                              // 26
    hi_08_phi             = new TH1F("hi_08_phi",             "hi_08_phi",             300,   -1.6708,  4.8124);                              // 27
    
    di_prot_q_dphi_pos_pos_05 = new TH2F("di_prot_q_dphi_pos_pos_05", "di_prot_q_dphi_pos_pos_05", 76, 0.00, 3.8, 300, -1.6708, 4.8124);      // 28
    di_prot_q_dphi_pos_neg_05 = new TH2F("di_prot_q_dphi_pos_neg_05", "di_prot_q_dphi_pos_neg_05", 76, 0.00, 3.8, 300, -1.6708, 4.8124);      // 29
    di_prot_q_dphi_neg_neg_05 = new TH2F("di_prot_q_dphi_neg_neg_05", "di_prot_q_dphi_neg_neg_05", 76, 0.00, 3.8, 300, -1.6708, 4.8124);      // 30
    
    di_prot_q_dphi_pos_pos_08 = new TH2F("di_prot_q_dphi_pos_pos_08", "di_prot_q_dphi_pos_pos_08", 76, 0.00, 3.8, 300, -1.6708, 4.8124);      // 31
    di_prot_q_dphi_pos_neg_08 = new TH2F("di_prot_q_dphi_pos_neg_08", "di_prot_q_dphi_pos_neg_08", 76, 0.00, 3.8, 300, -1.6708, 4.8124);      // 32
    di_prot_q_dphi_neg_neg_08 = new TH2F("di_prot_q_dphi_neg_neg_08", "di_prot_q_dphi_neg_neg_08", 76, 0.00, 3.8, 300, -1.6708, 4.8124);      // 33
    
//    track_cor_radius      = new TH2F("track_cor_radius",        "track_cor_radius",    160,  1.00,   5.00,     325,   -3.53,    3.53);        // 34
//    track_cor_radius_cut  = new TH2F("track_cor_radius_cut",    "track_cor_radius_cut",160,  1.00,   5.00,     325,   -3.53,    3.53);        // 35
    
    // objects added to output file
    
    fOutputList->Add(fHistPt);                    //  1
    fOutputList->Add(cent_ntracks);               //  2
    
    fOutputList->Add(m2_pos);                     //  3
    fOutputList->Add(m2_neg);                     //  4 
    fOutputList->Add(beta_p_pos);                 //  5
    fOutputList->Add(beta_p_neg);                 //  6
//  fOutputList->Add(deltat_p_pos);               //  7
//  fOutputList->Add(deltat_p_neg);               //  8

    fOutputList->Add(m2_pos_cut);                 //  9
    fOutputList->Add(m2_neg_cut);                 // 10
    fOutputList->Add(beta_p_pos_cut);             // 11
    fOutputList->Add(beta_p_neg_cut);             // 12
//  fOutputList->Add(deltat_p_pos_cut);           // 13
//  fOutputList->Add(deltat_p_neg_cut);           // 14

    fOutputList->Add(prot_per_event);             // 15
    fOutputList->Add(prot_per_event_pos);         // 16
    fOutputList->Add(prot_per_event_neg);         // 17
    fOutputList->Add(m2_pos_cut_T);               // 18
    fOutputList->Add(m2_neg_cut_T);               // 19

    fOutputList->Add(prot_phi_hist);              // 20

    fOutputList->Add(pt_mom);                     // 21
    fOutputList->Add(proton_pt_q_pos_pos);        // 22
    fOutputList->Add(proton_pt_q_pos_neg);        // 23
    fOutputList->Add(proton_pt_q_neg_neg);        // 24

    fOutputList->Add(di_proton_phi);              // 25
    fOutputList->Add(hi_05_phi);                  // 26
    fOutputList->Add(hi_08_phi);                  // 27

    fOutputList->Add(di_prot_q_dphi_pos_pos_05);  // 28
    fOutputList->Add(di_prot_q_dphi_pos_neg_05);  // 29
    fOutputList->Add(di_prot_q_dphi_neg_neg_05);  // 30
    
    fOutputList->Add(di_prot_q_dphi_pos_pos_08);  // 31
    fOutputList->Add(di_prot_q_dphi_pos_neg_08);  // 32
    fOutputList->Add(di_prot_q_dphi_neg_neg_08);  // 33

//    fOutputList->Add(track_cor_radius);           // 34
//    fOutputList->Add(track_cor_radius_cut);       // 35
       
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFprot::UserExec(Option_t *)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    Int_t iTracks(fAOD->GetNumberOfTracks());


    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0A"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////



    
    int NProt                = 0;
    int NProt_pos            = 0;
    int NProt_neg            = 0;

    int prot_track_num_1    = -9999;
    int prot_track_num_2    = -9999;
    int associated_tracks    = 0;

    
    int high_05_track_num    = -9999;
    int high_05_track_count  = 0;
    Float_t high_05_track_pt = 0.0;

    int high_08_track_num    = -9999;
    int high_08_track_count  = 0;
    Float_t high_08_track_pt = 0.0;

    
    
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

    // markA

    // loop over all these tracks
    for(Int_t i(0); i < iTracks; i++)
    {

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                                                       // if we failed, skip this track
	if(!(track->IsHybridGlobalConstrainedGlobal())){   continue;   }

	if(!track->IsPrimaryCandidate())   continue;

	Double_t pt            = track->Pt();
	Short_t  yes_high_pt   = 0;
	
	if(pt >= 5.0)
	{
	    if(high_05_track_pt < pt)
	    {
		high_05_track_num = i;
		high_05_track_pt  = pt;	
	    }
	    high_05_track_count++;
	    yes_high_pt = 1;
	}

	if(pt >= 8.0)
	{
	    if(high_08_track_pt < pt)
	    {
		high_08_track_num = i;
		high_08_track_pt  = pt;	
	    }
	    high_08_track_count++;
	}

	if(yes_high_pt == 1) continue;
	
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

	    

//	Double_t energy         = track->E();
//	float    mass           = track->M();
	float    mom            = track->P();
//	Double_t tpc_mom        = track->GetTPCmomentum();
	Short_t  charge         = track->Charge();
//	Double_t eta            = track->Eta();

	Double_t sigma_min      = -999.0;
	int      id             = 14;
	Double_t prot_mean      = 0.0;
	Double_t prot_sigma     = 0.0;


	Double_t phi = track->Phi();
	Double_t eta = track->Eta();


	if(fabs(eta > 0.8))   continue;

	   
	if(phi <  -TMath::PiOver2())    phi = phi + TMath::TwoPi();			if(phi <  -TMath::PiOver2())    phi = phi + TMath::TwoPi();
	if(phi > 3*TMath::PiOver2())    phi = phi - TMath::TwoPi();			if(phi > 3*TMath::PiOver2())    phi = phi - TMath::TwoPi();


	fHistPt->Fill(pt);

	if(isModGlobal == 1)
	{
	    Double_t m2tof = get_mass_squared(track);
	    Float_t dedx   = track->GetTPCsignal();
	    Float_t deltat = tof_minus_tpion(track);

	    pt_mom->Fill(pt, mom);

	    Float_t beta_kaon = 0.0;
	    Float_t beta      = 0.0;

	    beta      = Beta(track);
	    beta_kaon = mom / sqrt(0.2437169 + pow(mom,2));  // mass-squared of kaon is 0.2437169
	    
	    if(charge > 0)
	    {
		m2_pos->Fill(        pt, m2tof);
		beta_p_pos->Fill(   mom, Beta(track));
//		deltat_p_pos->Fill(  pt, deltat);

		
		if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340  &&  isGlobal == 1)
		{
		    if(    (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
			   )
		    {
			AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
			if(nsigmaTPC <= 2.0)
			{
			    m2_pos_cut->Fill(        pt, m2tof);
			    beta_p_pos_cut->Fill(   mom, Beta(track));
//			    deltat_p_pos_cut->Fill(  pt, deltat);

			    if(mom >= 0.5  &&  mom < 4.0)
			    {
				for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][0][w]);   }
				prot_mean = fit_prot_curve->Eval(mom);
				for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[0][1][w]);   }
				prot_sigma = fit_prot_curve->Eval(mom);

				if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
				{
				    if(NProt == 0)
				    {
					prot_track_num_1 = i;
				    }
				    else if(NProt == 1)
				    {
					prot_track_num_2 = i;
				    }
				    m2_pos_cut_T->Fill(mom,m2tof);
				    NProt++;
				    NProt_pos++;
				}		    
			    }			    
			}
		    }
		}
	    }
	    
	    else if(charge < 0)
	    {
		m2_neg->Fill(      pt,m2tof);
		beta_p_neg->Fill(   mom, Beta(track));
//		deltat_p_neg->Fill(  pt, deltat);
		
		if(dedx > 7.91143*deltat+28.8714  &&  deltat > 0.07216*dedx-5.11340  &&  isGlobal == 1)
		{
		    if(    (1.0 <= deltat  &&  deltat < 6.0  &&  dedx <   9.6774*deltat+46.7742)
		       ||  (0.5 <= deltat  &&  deltat < 1.0  &&  dedx <  56.4516)
		       ||  (0.5 >  deltat  &&                    dedx <  56.4516)
		       ||  (6.0 <= deltat  &&                    dedx <  12.9032*deltat+27.4193)
			   )
		    {
			AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 4, nsigmaTPC);
			if(nsigmaTPC <= 2.0)
			{	
			    m2_neg_cut->Fill(      pt, m2tof);
			    beta_p_neg_cut->Fill(   mom, Beta(track));
//			    deltat_p_neg_cut->Fill(  pt, deltat);


			    if(mom >= 0.5  &&  mom < 4.0)
			    {
				for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[1][0][w]);   }
				prot_mean = fit_prot_curve->Eval(mom);
				for(int w=0; w<4; w++){   fit_prot_curve->SetParameter(w, prot_curves[1][1][w]);   }
				prot_sigma = fit_prot_curve->Eval(mom);
				
				if(m2tof < prot_mean + cut_width * prot_sigma  &&   m2tof > prot_mean - cut_width * prot_sigma)
				{
				    if(NProt == 0)
				    {
					prot_track_num_1 = i;
				    }
				    else if(NProt == 1)
				    {
					prot_track_num_2 = i;
				    }
				    m2_neg_cut_T->Fill(mom,m2tof);
				    NProt++;
				    NProt_neg++;
				}		    
			    }
			}
		    }
		}
	    }
	}    //   track has both TPC and TOF PID, tracks are o.k.
    }        //   end of track loop


    prot_per_event->Fill(NProt);
    prot_per_event_pos->Fill(NProt_pos);
    prot_per_event_neg->Fill(NProt_neg);


    // markA
    if(NProt == 2  &&  high_05_track_count > 0)
    {
//	ProteronCandidates++;
	int A = prot_track_num_1;
	int B = prot_track_num_2;
	int H = high_05_track_num;
	
	AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));
	AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	if(!trackA  ||  !trackB  ||  !trackH) {}
	else
	{
	    Double_t prot_phi_A      = trackA->Phi();
	    Double_t prot_pt_A       = trackA->Pt();
	    Double_t prot_mom_A      = trackA->P();
	    Short_t  prot_charge_A   = trackA->Charge();
	    Double_t prot_eta_A      = trackA->Eta();

	    Double_t prot_et_A       = 0.0;
	    Double_t prot_e_A        = 0.0;
	    Double_t prot_px_A       = trackA->Px();
	    Double_t prot_py_A       = trackA->Py();
	    Double_t prot_pz_A       = trackA->Pz();
	    
	    prot_et_A = sqrt(pow(prot_pt_A ,2) + pow(0.93827,2));
	    prot_e_A  = sqrt(pow(prot_mom_A,2) + pow(0.93827,2));

	    prot_phi_hist->Fill(prot_phi_A);
	    	    
	    
	    Double_t prot_phi_B      = trackB->Phi();
	    Double_t prot_pt_B       = trackB->Pt();
	    Double_t prot_mom_B      = trackB->P();
	    Short_t  prot_charge_B   = trackB->Charge();
	    Double_t prot_eta_B      = trackB->Eta();

	    Double_t prot_et_B       = 0.0;
	    Double_t prot_e_B        = 0.0;
	    Double_t prot_px_B       = trackB->Px();
	    Double_t prot_py_B       = trackB->Py();
	    Double_t prot_pz_B       = trackB->Pz();
	    
	    prot_et_B = sqrt(pow(prot_pt_B ,2) + pow(0.93827,2));
	    prot_e_B  = sqrt(pow(prot_mom_B,2) + pow(0.93827,2));

	    
	    Double_t phi_H      = trackH->Phi();
	    Double_t pt_H       = trackH->Pt();
	    Double_t mom_H      = trackH->P();
	    Short_t  charge_H   = trackH->Charge();
	    Double_t eta_H      = trackH->Eta();

	    Double_t et_H       = 0.0;
	    Double_t e_H        = 0.0;
	    Double_t px_H       = trackH->Px();
	    Double_t py_H       = trackH->Py();
	    Double_t pz_H       = trackH->Pz();
	    
	    et_H = sqrt(pow(pt_H ,2) + pow(0.93827,2));
	    e_H  = sqrt(pow(mom_H,2) + pow(0.93827,2));

	    
	    Double_t Q2 = 0.00;	    Q2 = pow(prot_px_A - prot_px_B,2)  +  pow(prot_py_A - prot_py_B,2)  +  pow(prot_pz_A - prot_pz_B,2);

	    Double_t Q  = 0.00;	    Q  = sqrt(Q2);
	    Double_t Sx = prot_px_A + prot_px_B;
	    Double_t Sy = prot_py_A + prot_py_B;

	    Double_t Sphi = atan2(Sy,Sx);
	    if(Sphi <  -TMath::PiOver2())    Sphi = Sphi + TMath::TwoPi();
	    if(Sphi <  -TMath::PiOver2())    Sphi = Sphi + TMath::TwoPi();
	    if(Sphi > 3*TMath::PiOver2())    Sphi = Sphi - TMath::TwoPi();
	    if(Sphi > 3*TMath::PiOver2())    Sphi = Sphi - TMath::TwoPi();
				
	    di_proton_phi->Fill(Sphi);

	    Double_t Sdphi = Sphi - phi_H;
	    if(Sdphi <  -TMath::PiOver2())    Sdphi = Sdphi + TMath::TwoPi();
	    if(Sdphi <  -TMath::PiOver2())    Sdphi = Sdphi + TMath::TwoPi();
	    if(Sdphi > 3*TMath::PiOver2())    Sdphi = Sdphi - TMath::TwoPi();
	    if(Sdphi > 3*TMath::PiOver2())    Sdphi = Sdphi - TMath::TwoPi();

	    if(phi_H <  -TMath::PiOver2())    phi_H = phi_H + TMath::TwoPi();
	    if(phi_H <  -TMath::PiOver2())    phi_H = phi_H + TMath::TwoPi();
	    if(phi_H > 3*TMath::PiOver2())    phi_H = phi_H - TMath::TwoPi();
	    if(phi_H > 3*TMath::PiOver2())    phi_H = phi_H - TMath::TwoPi();
	    hi_05_phi->Fill(phi_H);

	    
	    if(prot_charge_A > 0)
	    {
		if(prot_charge_B > 0)
		{
		    proton_pt_q_pos_pos->Fill(prot_pt_A, Q);
		    di_prot_q_dphi_pos_pos_05->Fill(Q, Sdphi);
		}
		else if(prot_charge_B < 0)
		{
		    proton_pt_q_pos_neg->Fill(prot_pt_A, Q);
		    di_prot_q_dphi_pos_neg_05->Fill(Q, Sdphi);
		}
	    }
	    else if(prot_charge_A < 0)
	    {
		if(prot_charge_B > 0)
		{
		    proton_pt_q_pos_neg->Fill(prot_pt_A, Q);
		    di_prot_q_dphi_pos_neg_05->Fill(Q, Sdphi);
		}
		else if(prot_charge_B < 0)
		{
		    proton_pt_q_neg_neg->Fill(prot_pt_A, Q);
		    di_prot_q_dphi_neg_neg_05->Fill(Q, Sdphi);
		}
	    }
	}
    }

    if(NProt == 2  &&  high_08_track_count > 0)
    {
//	ProteronCandidates++;
	int A = prot_track_num_1;
	int B = prot_track_num_2;
	int H = high_08_track_num;
	
	AliAODTrack* trackA = static_cast<AliAODTrack*>(fAOD->GetTrack(A));
	AliAODTrack* trackB = static_cast<AliAODTrack*>(fAOD->GetTrack(B));
	AliAODTrack* trackH = static_cast<AliAODTrack*>(fAOD->GetTrack(H));
	if(!trackA  ||  !trackB  ||  !trackH) {}
	else
	{
	    Double_t prot_phi_A      = trackA->Phi();
	    Double_t prot_pt_A       = trackA->Pt();
	    Double_t prot_mom_A      = trackA->P();
	    Short_t  prot_charge_A   = trackA->Charge();
	    Double_t prot_eta_A      = trackA->Eta();

	    Double_t prot_et_A       = 0.0;
	    Double_t prot_e_A        = 0.0;
	    Double_t prot_px_A       = trackA->Px();
	    Double_t prot_py_A       = trackA->Py();
	    Double_t prot_pz_A       = trackA->Pz();
	    
	    prot_et_A = sqrt(pow(prot_pt_A ,2) + pow(0.93827,2));
	    prot_e_A  = sqrt(pow(prot_mom_A,2) + pow(0.93827,2));

//	    prot_phi_hist->Fill(prot_phi_A);
	    	    
	    
	    Double_t prot_phi_B      = trackB->Phi();
	    Double_t prot_pt_B       = trackB->Pt();
	    Double_t prot_mom_B      = trackB->P();
	    Short_t  prot_charge_B   = trackB->Charge();
	    Double_t prot_eta_B      = trackB->Eta();

	    Double_t prot_et_B       = 0.0;
	    Double_t prot_e_B        = 0.0;
	    Double_t prot_px_B       = trackB->Px();
	    Double_t prot_py_B       = trackB->Py();
	    Double_t prot_pz_B       = trackB->Pz();
	    
	    prot_et_B = sqrt(pow(prot_pt_B ,2) + pow(0.93827,2));
	    prot_e_B  = sqrt(pow(prot_mom_B,2) + pow(0.93827,2));

	    
	    Double_t phi_H      = trackH->Phi();
	    Double_t pt_H       = trackH->Pt();
	    Double_t mom_H      = trackH->P();
	    Short_t  charge_H   = trackH->Charge();
	    Double_t eta_H      = trackH->Eta();

	    Double_t et_H       = 0.0;
	    Double_t e_H        = 0.0;
	    Double_t px_H       = trackH->Px();
	    Double_t py_H       = trackH->Py();
	    Double_t pz_H       = trackH->Pz();
	    
	    et_H = sqrt(pow(pt_H ,2) + pow(0.93827,2));
	    e_H  = sqrt(pow(mom_H,2) + pow(0.93827,2));

	    
	    Double_t Q2 = 0.00;	    Q2 = pow(prot_px_A - prot_px_B,2)  +  pow(prot_py_A - prot_py_B,2)  +  pow(prot_pz_A - prot_pz_B,2);
	    Double_t Q  = 0.00;	    Q  = sqrt(Q2);
	    
	    Double_t Sx = prot_px_A + prot_px_B;
	    Double_t Sy = prot_py_A + prot_py_B;

	    Double_t Sphi = atan2(Sy,Sx);
	    if(Sphi <  -TMath::PiOver2())    Sphi = Sphi + TMath::TwoPi();
	    if(Sphi <  -TMath::PiOver2())    Sphi = Sphi + TMath::TwoPi();
	    if(Sphi > 3*TMath::PiOver2())    Sphi = Sphi - TMath::TwoPi();
	    if(Sphi > 3*TMath::PiOver2())    Sphi = Sphi - TMath::TwoPi();
				
//	    di_proton_phi->Fill(Sphi);   // already filled under 05


	    Double_t Sdphi = Sphi - phi_H;
	    if(Sdphi <  -TMath::PiOver2())    Sdphi = Sdphi + TMath::TwoPi();
	    if(Sdphi <  -TMath::PiOver2())    Sdphi = Sdphi + TMath::TwoPi();
	    if(Sdphi > 3*TMath::PiOver2())    Sdphi = Sdphi - TMath::TwoPi();
	    if(Sdphi > 3*TMath::PiOver2())    Sdphi = Sdphi - TMath::TwoPi();


	    if(phi_H <  -TMath::PiOver2())    phi_H = phi_H + TMath::TwoPi();
	    if(phi_H <  -TMath::PiOver2())    phi_H = phi_H + TMath::TwoPi();
	    if(phi_H > 3*TMath::PiOver2())    phi_H = phi_H - TMath::TwoPi();
	    if(phi_H > 3*TMath::PiOver2())    phi_H = phi_H - TMath::TwoPi();

	    hi_08_phi->Fill(phi_H);

	    
	    
	    if(prot_charge_A > 0)
	    {
		if(prot_charge_B > 0)
		{
//		    proton_pt_q_pos_pos->Fill(prot_pt_A, Q);    // already filled under 05
		    di_prot_q_dphi_pos_pos_08->Fill(Q, Sdphi);
		}
		else if(prot_charge_B < 0)
		{
//		    proton_pt_q_pos_neg->Fill(prot_pt_A, Q);    // already filled under 05
		    di_prot_q_dphi_pos_neg_08->Fill(Q, Sdphi);
		}
	    }
	    else if(prot_charge_A < 0)
	    {
		if(prot_charge_B > 0)
		{
//		    proton_pt_q_pos_neg->Fill(prot_pt_A, Q);
		    di_prot_q_dphi_pos_neg_08->Fill(Q, Sdphi);
		}
		else if(prot_charge_B < 0)
		{
//		    proton_pt_q_neg_neg->Fill(prot_pt_A, Q);
		    di_prot_q_dphi_neg_neg_08->Fill(Q, Sdphi);
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
