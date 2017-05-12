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
fAOD(0), fOutputList(0), fPIDResponse(0), fHistPt(0),  cent_ntracks(0),
    
    m_squared_pos(0),       beta_vs_mom_pos(0),        deltat_vs_mom_pos(0),
    m_squared_neg(0),       beta_vs_mom_neg(0),        deltat_vs_mom_neg(0),
    dedx_vs_mom_pos(0),     dedx_vs_deltat_pos(0),     dedx_mom_deltat_pos(0),      //     dedx_mom_m2_pos(0),
    dedx_vs_mom_neg(0),     dedx_vs_deltat_neg(0),     dedx_mom_deltat_neg(0),      //     dedx_mom_m2_neg(0),

    m_squared_pos_cut(0),   beta_vs_mom_pos_cut(0),    deltat_vs_mom_pos_cut(0),
    m_squared_neg_cut(0),   beta_vs_mom_neg_cut(0),    deltat_vs_mom_neg_cut(0),
    dedx_vs_mom_pos_cut(0), dedx_vs_deltat_pos_cut(0), dedx_mom_deltat_pos_cut(0),  //     dedx_mom_m2_pos_cut(0),
    dedx_vs_mom_neg_cut(0), dedx_vs_deltat_neg_cut(0), dedx_mom_deltat_neg_cut(0)   //     dedx_mom_m2_neg_cut(0),

//    path_length(0), ttof(0), theta(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0), fHistPt(0),  cent_ntracks(0),

    m_squared_pos(0),       beta_vs_mom_pos(0),        deltat_vs_mom_pos(0),
    m_squared_neg(0),       beta_vs_mom_neg(0),        deltat_vs_mom_neg(0),
    dedx_vs_mom_pos(0),     dedx_vs_deltat_pos(0),     dedx_mom_deltat_pos(0),      //     dedx_mom_m2_pos(0),
    dedx_vs_mom_neg(0),     dedx_vs_deltat_neg(0),     dedx_mom_deltat_neg(0),      //     dedx_mom_m2_neg(0),

    m_squared_pos_cut(0),   beta_vs_mom_pos_cut(0),    deltat_vs_mom_pos_cut(0),
    m_squared_neg_cut(0),   beta_vs_mom_neg_cut(0),    deltat_vs_mom_neg_cut(0),
    dedx_vs_mom_pos_cut(0), dedx_vs_deltat_pos_cut(0), dedx_mom_deltat_pos_cut(0),  //     dedx_mom_m2_pos_cut(0),
    dedx_vs_mom_neg_cut(0), dedx_vs_deltat_neg_cut(0), dedx_mom_deltat_neg_cut(0)   //     dedx_mom_m2_neg_cut(0),

//    path_length(0), ttof(0), theta(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::~AliAnalysisTaskCorPIDTOFQA()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //

//  fTrackSelection = new AliEmcalTrackSelectionAOD;
//  fTrackSelection->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks, "lhc13c");
//    pi = TMath::Pi();
//    fout.open("output.txt");

// 3.38835 0.0203511 0.101903 -0.0487972 0.0513143 0.0815525 3.69323 -0.0206797 -0.200752 0.00644255 0.0391803 0.0399812
    
//  Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]

//  TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);
    fit_deut_curve->SetParNames("a", "b x", "c/#sqrt(x)");
	
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    fHistPt                 = new TH1F("fHistPt",                   "Pt()",                  10000,     0,     10);
    cent_ntracks            = new TH2F("cent_vs_ntracks",           "cent_vs_ntracks",         100,     0,    100,     100,      0,    800);
    
    m_squared_pos           = new TH2F("m_squared_pos",             "m_squared_pos",           320,   0.0,    8.0,    2400,   -1.0,    7.0);
    m_squared_neg           = new TH2F("m_squared_neg",             "m_squared_neg",           320,   0.0,    8.0,    2400,   -1.0,    7.0);
   
    beta_vs_mom_pos         = new TH2F("beta_vs_mom_pos",           "beta_vs_mom_pos",         780,   0.0,    8.0,    3000,    0.1,    1.1);
    beta_vs_mom_neg         = new TH2F("beta_vs_mom_neg",           "beta_vs_mom_neg",         780,   0.0,    8.0,    3000,    0.1,    1.1);
    deltat_vs_mom_pos       = new TH2F("deltat_vs_mom_pos",         "deltat_vs_mom_pos",     17000,   0.3,    2.0,    1500,   -2.5,   12.5);
    deltat_vs_mom_neg       = new TH2F("deltat_vs_mom_neg",         "deltat_vs_mom_neg",     17000,   0.3,    2.0,    1500,   -2.5,   12.5);

    dedx_vs_mom_pos         = new TH2F("dedx_vs_mom_pos",           "dedx_vs_mom_pos",         780,   0.2,    8.0,     250,    0.0,  500.0);
    dedx_vs_mom_neg         = new TH2F("dedx_vs_mom_neg",           "dedx_vs_mom_neg",         780,   0.2,    8.0,     250,    0.0,  500.0);
    dedx_vs_deltat_pos      = new TH2F("dedx_vs_deltat_pos",        "dedx_vs_deltat_pos",     2100,  -1.0,   20.0,     250,    0.0,  500.0);
    dedx_vs_deltat_neg      = new TH2F("dedx_vs_deltat_neg",        "dedx_vs_deltat_neg",     2100,  -1.0,   20.0,     250,    0.0,  500.0);
    dedx_mom_deltat_pos     = new TH3F("dedx_mom_deltat_pos",       "dedx_mom_deltat_pos",     250,   0.0,  500.0,     156,    0.2,    8.0,      420, -1.0, 20.0);  //// (dedx, mom, deltat)
    dedx_mom_deltat_neg     = new TH3F("dedx_mom_deltat_neg",       "dedx_mom_deltat_neg",     250,   0.0,  500.0,     156,    0.2,    8.0,      420, -1.0, 20.0);  //// (dedx, mom, deltat)

//  dedx_mom_m2_pos         = new TH3F("dedx_mom_m2_pos",           "dedx_mom_m2_pos",         250,   0.0,  500.0,     156,    0.2,    8.0,      240, -1.0,  7.0);  //// (dedx, mom, m2)
//  dedx_mom_m2_neg         = new TH3F("dedx_mom_m2_neg",           "dedx_mom_m2_neg",         250,   0.0,  500.0,     156,    0.2,    8.0,      240, -1.0,  7.0);  //// (dedx, mom, m2)
    
    m_squared_pos_cut       = new TH2F("m_squared_pos_cut",         "m_squared_pos_cut",       320,   0.0,    8.0,    2400,   -1.0,    7.0);
    m_squared_neg_cut       = new TH2F("m_squared_neg_cut",         "m_squared_neg_cut",       320,   0.0,    8.0,    2400,   -1.0,    7.0);
    beta_vs_mom_pos_cut     = new TH2F("beta_vs_mom_pos_cut",       "beta_vs_mom_pos_cut",     780,   0.0,    8.0,    3000,    0.1,    1.1);
    beta_vs_mom_neg_cut     = new TH2F("beta_vs_mom_neg_cut",       "beta_vs_mom_neg_cut",     780,   0.0,    8.0,    3000,    0.1,    1.1);
    deltat_vs_mom_pos_cut   = new TH2F("deltat_vs_mom_pos_cut",     "deltat_vs_mom_pos_cut", 17000,   0.3,    2.0,    1500,   -2.5,   12.5);
    deltat_vs_mom_neg_cut   = new TH2F("deltat_vs_mom_neg_cut",     "deltat_vs_mom_neg_cut", 17000,   0.3,    2.0,    1500,   -2.5,   12.5);

    dedx_vs_mom_pos_cut     = new TH2F("dedx_vs_mom_pos_cut",       "dedx_vs_mom_pos_cut",     780,   0.2,    8.0,     250,    0.0,  500.0);
    dedx_vs_mom_neg_cut     = new TH2F("dedx_vs_mom_neg_cut",       "dedx_vs_mom_neg_cut",     780,   0.2,    8.0,     250,    0.0,  500.0);
    dedx_vs_deltat_pos_cut  = new TH2F("dedx_vs_deltat_pos_cut",    "dedx_vs_deltat_pos_cut", 2100,  -1.0,   20.0,     250,    0.0,  500.0);
    dedx_vs_deltat_neg_cut  = new TH2F("dedx_vs_deltat_neg_cut",    "dedx_vs_deltat_neg_cut", 2100,  -1.0,   20.0,     250,    0.0,  500.0);
    dedx_mom_deltat_pos_cut = new TH3F("dedx_mom_deltat_pos_cut",   "dedx_mom_deltat_pos_cut", 250,   0.0,  500.0,     156,    0.2,    8.0,      420, -1.0, 20.0);  //// (dedx, mom, deltat)
    dedx_mom_deltat_neg_cut = new TH3F("dedx_mom_deltat_neg_cut",   "dedx_mom_deltat_neg_cut", 250,   0.0,  500.0,     156,    0.2,    8.0,      420, -1.0, 20.0);  //// (dedx, mom, deltat)

//    dedx_mom_m2_pos_cut     = new TH3F("dedx_mom_m2_pos_cut",       "dedx_mom_m2_pos_cut",     250,   0.0,  500.0,     156,    0.2,    8.0,      240, -1.0,  7.0);  //// (dedx, mom, m2)
//    dedx_mom_m2_neg_cut     = new TH3F("dedx_mom_m2_neg_cut",       "dedx_mom_m2_neg_cut",     250,   0.0,  500.0,     156,    0.2,    8.0,      240, -1.0,  7.0);  //// (dedx, mom, m2)

//    path_length             = new TH1F("path_length",               "path_length",             200, 382.0,  398.0);
//    ttof                    = new TH1F("ttof",                      "time of flight",          500,  12.0,   20.0);
//    alpha                   = new TH1F("alpha",                     "alpha",                  1000,     0,    4.0);
//    theta                   = new TH1F("theta",                     "theta",                   100,  -7.0,    7.0);
//    dca                     = new TH2F("dca",                       "dca",                     100,  -5.0,    5.0,     100,   -5.0,    5.0);
//    generic                 = new TH1F("generic",                   "generic",                 4, -1,3);
    //markc


//    fOutputList->Add(generic);
//    fOutputList->Add(dca);
//    fOutputList->Add(alpha);
//    fOutputList->Add(theta);
//    fOutputList->Add(path_length);
//    fOutputList->Add(ttof);
    
    fOutputList->Add(fHistPt);                                                                      // objects added to output file
    fOutputList->Add(cent_ntracks);
    
    fOutputList->Add(m_squared_pos);
    fOutputList->Add(m_squared_pos_cut);    
    fOutputList->Add(m_squared_neg);
    fOutputList->Add(m_squared_neg_cut);    
    fOutputList->Add(beta_vs_mom_pos);
    fOutputList->Add(beta_vs_mom_pos_cut);    
    fOutputList->Add(beta_vs_mom_neg);
    fOutputList->Add(beta_vs_mom_neg_cut);    
    fOutputList->Add(deltat_vs_mom_pos);
    fOutputList->Add(deltat_vs_mom_pos_cut);    
    fOutputList->Add(deltat_vs_mom_neg);
    fOutputList->Add(deltat_vs_mom_neg_cut);

    
    
    fOutputList->Add(dedx_vs_mom_pos);
    fOutputList->Add(dedx_vs_mom_pos_cut);
    fOutputList->Add(dedx_vs_mom_neg);
    fOutputList->Add(dedx_vs_mom_neg_cut);    
    fOutputList->Add(dedx_vs_deltat_pos);
    fOutputList->Add(dedx_vs_deltat_pos_cut);    
    fOutputList->Add(dedx_vs_deltat_neg);
    fOutputList->Add(dedx_vs_deltat_neg_cut);    
    fOutputList->Add(dedx_mom_deltat_pos);
    fOutputList->Add(dedx_mom_deltat_pos_cut);    
    fOutputList->Add(dedx_mom_deltat_neg);
    fOutputList->Add(dedx_mom_deltat_neg_cut);

//    fOutputList->Add(dedx_mom_m2_pos);
//    fOutputList->Add(dedx_mom_m2_pos_cut);    
//    fOutputList->Add(dedx_mom_m2_neg);
//    fOutputList->Add(dedx_mom_m2_neg_cut);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserExec(Option_t *)
{

//    	Double_t pi = 3.1415926535897932384626434;
    NEvents++;
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram



    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event



//  fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kTOF0);
//  fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kT00);
//  fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kBest0);

    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////

    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////


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
//	if(!(track->IsHybridGlobalConstrainedGlobal())){   continue;   }

//	generic->Fill(track->IsPrimaryCandidate());
	
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



	//generic->Fill(nPointsITS);

	
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
	    Double_t m2tof          = get_mass_squared(track);

	    if(fabs(pt-1.0) < 0.01)
	    {
		Double_t track_theta    = 0.0;
		track_theta             = track->Theta();
//		theta->Fill(track_theta);
		Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
		Double_t length         = 0.0;
		length                  = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
		length                  = length * 100;                                                                 // convert mm to c

		Int_t channel           = 0;
		channel                 = fPIDResponse->GetTOFResponse().GetTOFchannel(track);
//		Double_t xy_length      = 0.00;
//		xy_length               = length * TMath::Sin(track_theta);

		if(channel >=4330)
		{
		    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
		    Double_t stop_time      = track->GetTOFsignal();
		    Double_t time_of_flight = stop_time - start_time;
		    time_of_flight          = time_of_flight * 0.001;                                                       // convert ps to ns

//		    path_length->Fill(length);
//		    ttof->Fill(time_of_flight);
		}
	    }
	    
//	    Float_t mom_tpc = track->GetTPCmomentum();
	    Float_t dedx   = track->GetTPCsignal();
	    Float_t deltat = tof_minus_tpion(track);

//	    int pid_cut = 0;
//	    Double_t exptimes[20];
//	    track->GetIntegratedTimes(exptimes);
//	    Double_t time = track->GetTOFsignal();  //
//	    fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,(AliPID::EParticleType) 2, nsigmaTOF);  // always 1 // 2 = pion, 4 = proton
//	    dedx_exp_deut->Fill(fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kProton));  //
//	    if(NEvents % 200 == 0)   cout<<fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kProton)<<" ";  // 10,000 to 50,000
//	    track->GetIntegratedTimes(exptimes, AliPID::kProton);
//	    track->GetIntegratedTimes(exptimes, AliPID::kDeuteron);
//	    if(NEvents % 200 == 0)   cout<<(exptimes[5]-exptimes[4])/1000<<" ";
//	    if(NEvents % 200 == 0)   cout<<(exptimes[5])/1000<<" ";  // doesn't work
//	    if(mom > 0.8  &&  mom < 1.0)
//	    {
//		if(NEvents % 200 == 0)   cout<<fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kDeuteron)*1E-3<<" ";
//		if(NEvents % 200 == 0)   cout<<fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kProton)*1E-3<<endl;
//	    }
//	    dedx_exp_deut->Fill(exptimes[4]/1000);
	    
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
		m_squared_pos->Fill(      pt, m2tof);
		beta_vs_mom_pos->Fill(   mom, Beta(track));
		dedx_vs_mom_pos->Fill(   mom, dedx);
		if(mom >= 2.0  &&  mom < 3.0)    dedx_vs_deltat_pos->Fill(deltat, dedx);
		dedx_mom_deltat_pos->Fill(dedx, mom, deltat);
//		dedx_mom_m2_pos->Fill(dedx, mom, m2tof);
		deltat_vs_mom_pos->Fill(  pt, deltat);

		
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
			    m_squared_pos_cut->Fill(      pt, m2tof);
			    beta_vs_mom_pos_cut->Fill(   mom, Beta(track));
			    deltat_vs_mom_pos_cut->Fill(  pt, deltat);
			    dedx_vs_mom_pos_cut->Fill(   mom, dedx);
			    if(mom >= 2.0  &&  mom < 3.0)    dedx_vs_deltat_pos_cut->Fill(deltat, dedx);
			    dedx_mom_deltat_pos_cut->Fill(dedx, mom, deltat);
//			    dedx_mom_m2_pos_cut->Fill(dedx, mom, m2tof);
			}
		    }
		}
		
		//Float_t beta_prot = 0.0;
		//Float_t beta_deut = 0.0;
		//beta_prot = mom / sqrt(0.88035 + pow(mom,2));  // mass-squared of proton is 0.88035
		//beta_deut = mom / sqrt(3.52332 + pow(mom,2));  // mass-squared of proton is 3.52332
		//if(fabs(beta_prot - beta) > fabs(beta_deut - beta))
	    }
	    
	    else if(charge < 0)
	    {
		m_squared_neg->Fill(      pt,m2tof);
		beta_vs_mom_neg->Fill(   mom, Beta(track));
		dedx_vs_mom_neg->Fill(   mom, dedx);
		if(mom >= 2.0  &&  mom < 3.0)    dedx_vs_deltat_neg->Fill(deltat, dedx);
		dedx_mom_deltat_neg->Fill(dedx, mom, deltat);
//		dedx_mom_m2_neg->Fill(dedx, mom, m2tof);
		deltat_vs_mom_neg->Fill(  pt, deltat);

		
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
			    m_squared_neg_cut->Fill(      pt, m2tof);
			    beta_vs_mom_neg_cut->Fill(   mom, Beta(track));
			    deltat_vs_mom_neg_cut->Fill(  pt, deltat);
			    dedx_vs_mom_neg_cut->Fill(   mom, dedx);
			    if(mom >= 2.0  &&  mom < 3.0)    dedx_vs_deltat_neg_cut->Fill(deltat, dedx);
			    dedx_mom_deltat_neg_cut->Fill(dedx, mom, deltat);
//			    dedx_mom_m2_neg_cut->Fill(dedx, mom, m2tof);
			}
		    }
		}

		//Float_t beta_prot = 0.0;
		//Float_t beta_deut = 0.0;
		//beta_prot = mom / sqrt(0.88035 + pow(mom,2));  // mass-squared of proton is 0.88035
		//beta_deut = mom / sqrt(3.52332 + pow(mom,2));  // mass-squared of proton is 3.52332
		//if(fabs(beta_prot - beta) > fabs(beta_deut - beta))
	    }
	}    //   track has both TPC and TOF PID, tracks are o.k.
    }        //   end of track loop



     
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
//    cout<<"tracks from events with deuterons: "<<NAssociatedTracks<<endl;
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
    Double_t m_squared      = 0.0;

    Double_t c2 = 29.9792458;
    
    m_squared = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);
    return m_squared;

}
