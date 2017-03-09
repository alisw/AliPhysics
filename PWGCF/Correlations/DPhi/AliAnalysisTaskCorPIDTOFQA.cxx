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





//Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]
//TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);

ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
fAOD(0), fOutputList(0), fPIDResponse(0),
    fHistPt(0),  cent_ntracks(0), m_squared_pos(0), m_squared_pos_cut(0), m_squared_neg(0), m_squared_neg_cut(0),
    plength_vs_mom_pos(0), plength_vs_mom_neg(0), ttof_vs_mom_pos(0), ttof_vs_mom_neg(0), beta_vs_mom_pos(0), beta_vs_mom_neg(0), deltat_vs_mom_pos(0), deltat_vs_mom_neg(0),
    plength_vs_mom_pos_cut(0), plength_vs_mom_neg_cut(0), ttof_vs_mom_pos_cut(0), ttof_vs_mom_neg_cut(0), beta_vs_mom_pos_cut(0), beta_vs_mom_neg_cut(0), deltat_vs_mom_pos_cut(0), deltat_vs_mom_neg_cut(0)
    ,h_dedx_exp_deut(0),h_dedx_exp_prot(0), h_dedx(0), h_TPC_resolution_deut(0), h_TPC_resolution_prot(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0),
    fHistPt(0),  cent_ntracks(0), m_squared_pos(0), m_squared_pos_cut(0), m_squared_neg(0), m_squared_neg_cut(0),
    plength_vs_mom_pos(0), plength_vs_mom_neg(0), ttof_vs_mom_pos(0), ttof_vs_mom_neg(0), beta_vs_mom_pos(0), beta_vs_mom_neg(0), deltat_vs_mom_pos(0), deltat_vs_mom_neg(0),
    plength_vs_mom_pos_cut(0), plength_vs_mom_neg_cut(0), ttof_vs_mom_pos_cut(0), ttof_vs_mom_neg_cut(0), beta_vs_mom_pos_cut(0), beta_vs_mom_neg_cut(0), deltat_vs_mom_pos_cut(0), deltat_vs_mom_neg_cut(0)
    ,h_dedx_exp_deut(0),h_dedx_exp_prot(0), h_dedx(0), h_TPC_resolution_deut(0), h_TPC_resolution_prot(0)
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


    fHistPt                 = new TH1F("fHistPt",                   "Pt()",                  100,     0,    10);
    cent_ntracks            = new TH2F("cent_vs_ntracks",           "cent_vs_ntracks",       100,     0,   100,     100,     0,   800);
    m_squared_pos           = new TH2F("m_squared_pos",             "m_squared_pos",         320,   0.0,   8.0,    2400,  -1.0,   7.0);
    m_squared_pos_cut       = new TH2F("m_squared_pos_cut",         "m_squared_pos_cut",     320,   0.0,   8.0,    2400,  -1.0,   7.0);
    m_squared_neg           = new TH2F("m_squared_neg",             "m_squared_neg",         320,   0.0,   8.0,    2400,  -1.0,   7.0);
    m_squared_neg_cut       = new TH2F("m_squared_neg_cut",         "m_squared_neg_cut",     320,   0.0,   8.0,    2400,  -1.0,   7.0);

    plength_vs_mom_pos      = new TH2F("plength_vs_mom_pos",        "plength_vs_mom_pos",    780,   0.2,   8.0,    4000, 350.0,   750);  // (GeV/c 10MeV/bin,   cm  1mm/bin)
    plength_vs_mom_neg      = new TH2F("plength_vs_mom_neg",        "plength_vs_mom_neg",    780,   0.2,   8.0,    4000, 350.0,   750);  // (GeV/c 10MeV/bin,   cm  1mm/bin)
    ttof_vs_mom_pos         = new TH2F("ttof_vs_mom_pos",           "ttof_vs_mom_pos",       780,   0.2,   8.0,    4000,  10.0,  50.0);  // (GeV/c 10MeV/bin,   ns 10ps/bin)
    ttof_vs_mom_neg         = new TH2F("ttof_vs_mom_neg",           "ttof_vs_mom_neg",       780,   0.2,   8.0,    4000,  10.0,  50.0);  // (GeV/c 10MeV/bin,   ns 10ps/bin)    
    beta_vs_mom_pos         = new TH2F("beta_vs_mom_pos",           "beta_vs_mom_pos",       780,   0.0,   8.0,    3000,   0.1,   1.1);
    beta_vs_mom_neg         = new TH2F("beta_vs_mom_neg",           "beta_vs_mom_neg",       780,   0.0,   8.0,    3000,   0.1,   1.1);
    deltat_vs_mom_pos       = new TH2F("deltat_vs_mom_pos",         "deltat_vs_mom_pos",     780,   0.2,   8.0,    2100,  -1.0,  20.0);
    deltat_vs_mom_neg       = new TH2F("deltat_vs_mom_neg",         "deltat_vs_mom_neg",     780,   0.2,   8.0,    2100,  -1.0,  20.0);

    plength_vs_mom_pos_cut  = new TH2F("plength_vs_mom_pos_cut",    "plength_vs_mom_pos_cut",780,   0.2,   8.0,    4000, 350.0,   750);  // (GeV/c 10MeV/bin,   cm  1mm/bin)
    plength_vs_mom_neg_cut  = new TH2F("plength_vs_mom_neg_cut",    "plength_vs_mom_neg_cut",780,   0.2,   8.0,    4000, 350.0,   750);  // (GeV/c 10MeV/bin,   cm  1mm/bin)
    ttof_vs_mom_pos_cut     = new TH2F("ttof_vs_mom_pos_cut",       "ttof_vs_mom_pos_cut",   780,   0.2,   8.0,    4000,  10.0,  50.0);  // (GeV/c 10MeV/bin,   ns 10ps/bin)
    ttof_vs_mom_neg_cut     = new TH2F("ttof_vs_mom_neg_cut",       "ttof_vs_mom_neg_cut",   780,   0.2,   8.0,    4000,  10.0,  50.0);  // (GeV/c 10MeV/bin,   ns 10ps/bin)    
    beta_vs_mom_pos_cut     = new TH2F("beta_vs_mom_pos_cut",       "beta_vs_mom_pos_cut",   780,   0.0,   8.0,    3000,   0.1,   1.1);
    beta_vs_mom_neg_cut     = new TH2F("beta_vs_mom_neg_cut",       "beta_vs_mom_neg_cut",   780,   0.0,   8.0,    3000,   0.1,   1.1);
    deltat_vs_mom_pos_cut   = new TH2F("deltat_vs_mom_pos_cut",     "deltat_vs_mom_pos_cut", 780,   0.2,   8.0,    2100,  -1.0,  20.0);
    deltat_vs_mom_neg_cut   = new TH2F("deltat_vs_mom_neg_cut",     "deltat_vs_mom_neg_cut", 780,   0.2,   8.0,    2100,  -1.0,  20.0);
    h_dedx_exp_deut          = new TH1F("dedx_exp_deut",            "dedx_exp_deut",        400,   0.0, 1200.0);
    h_dedx_exp_prot          = new TH1F("dedx_exp_prot",            "dedx_exp_prot",        400,   0.0, 800.0);
    h_dedx                   = new TH1F("dedx",                     "dedx",                 400,   0.0, 200.0);
    h_TPC_resolution_deut    = new TH1F("TPC_resolution_deut",      "TPC_resolution_deut",  400,   0.0,100.0);
    h_TPC_resolution_prot    = new TH1F("TPC_resolution_prot",      "TPC_resolution_prot",  400,   0.0, 40.0);
    
    //Mark B
    fOutputList->Add(fHistPt);                                                                      // objects added to output file
    fOutputList->Add(cent_ntracks);
    fOutputList->Add(m_squared_pos);
    fOutputList->Add(m_squared_pos_cut);
    fOutputList->Add(m_squared_neg);    
    fOutputList->Add(m_squared_neg_cut);

    
    fOutputList->Add(plength_vs_mom_pos);
    fOutputList->Add(plength_vs_mom_pos_cut);
    fOutputList->Add(plength_vs_mom_neg);
    fOutputList->Add(plength_vs_mom_neg_cut);
    fOutputList->Add(ttof_vs_mom_pos);
    fOutputList->Add(ttof_vs_mom_pos_cut);    
    fOutputList->Add(ttof_vs_mom_neg);
    fOutputList->Add(ttof_vs_mom_neg_cut);    
    fOutputList->Add(beta_vs_mom_pos);
    fOutputList->Add(beta_vs_mom_pos_cut);    
    fOutputList->Add(beta_vs_mom_neg);
    fOutputList->Add(beta_vs_mom_neg_cut);    
    fOutputList->Add(deltat_vs_mom_pos);
    fOutputList->Add(deltat_vs_mom_pos_cut);    
    fOutputList->Add(deltat_vs_mom_neg);
    fOutputList->Add(deltat_vs_mom_neg_cut);
    fOutputList->Add(h_dedx_exp_deut);
    fOutputList->Add(h_dedx_exp_prot);
    fOutputList->Add(h_dedx);
    fOutputList->Add(h_TPC_resolution_deut);
    fOutputList->Add(h_TPC_resolution_prot);
    
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
	if(!(track->IsHybridGlobalConstrainedGlobal())){   continue;   }


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
	
	if(phi <  -TMath::PiOver2())    phi = phi + TMath::TwoPi();			if(phi <  -TMath::PiOver2())    phi = phi + TMath::TwoPi();
	if(phi > 3*TMath::PiOver2())    phi = phi - TMath::TwoPi();			if(phi > 3*TMath::PiOver2())    phi = phi - TMath::TwoPi();

	fHistPt->Fill(track->Pt());
	
	if(isModGlobal == 1)
	{

	    Double_t m2tof          = get_mass_squared(track);
	    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
	    Double_t stop_time      = track->GetTOFsignal();
	    Double_t time_of_flight = stop_time - start_time;
	    time_of_flight          = time_of_flight * 0.001;                                                       // convert ps to ns
	    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
	    Double_t length         = 0.0;
	    length                  = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
	    length                  = length * 100;                                                                 // convert mm to cm

	    Float_t mom_tpc = track->GetTPCmomentum();

	    int pid_cut = 0;
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
	
	    for(int iSpecies = 2; iSpecies < 5; iSpecies++)
	    {
		fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,(AliPID::EParticleType) iSpecies, nsigmaTOF);
		if(nsigmaTOF < 5)
		{
		    Float_t dedx = track->GetTPCsignal()-2;
		    Float_t dedx_exp_prot = fPIDResponse->GetTPCResponse().GetExpectedSignal(mom_tpc,(AliPID::EParticleType) iSpecies);
		    Float_t resolution_TPC_prot = fPIDResponse->GetTPCResponse().GetExpectedSigma(mom_tpc,track->GetTPCsignalN(),AliPID::kProton);
		    if(TMath::Abs(dedx - dedx_exp_prot) < resolution_TPC_prot * 5.0)
		    {
			pid_cut = 1;
		    }
		}
	    }
	    

//	    if(fPIDResponse->GetTOFMismatchProbability(track))   pid_cut = 0;  // didn't work very well
//	    else                                                 pid_cut = 1;
	    
	    ///// allow deuterons /////

	    
	    if(pid_cut == 0)
	    {
		Float_t time_exp_deut = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kDeuteron)*1E-3;
		Float_t time_exp_prot = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kProton)*1E-3;
	    
		Float_t dedx = track->GetTPCsignal()-2;  // -2 factor is from copying from algorithm from Francesco Noferini
		Float_t time = track->GetTOFsignal();
		Float_t dedx_exp_prot      = fPIDResponse->GetTPCResponse().GetExpectedSignal(mom_tpc,AliPID::kProton);
		Float_t dedx_exp_deut      = fPIDResponse->GetTPCResponse().GetExpectedSignal(mom_tpc,AliPID::kDeuteron);
		Float_t resolutionTPC_prot = fPIDResponse->GetTPCResponse().GetExpectedSigma(mom_tpc,track->GetTPCsignalN(),AliPID::kProton);
//		Float_t dedxExp_deut       = fPIDResponse->GetTPCResponse().GetExpectedSignal(mom_tpc,AliPID::kDeuteron);
//		dedx_exp_deut->Fill(resolutionTPC_prot);
		Float_t resolutionTOF_deut = fPIDResponse->GetTOFResponse().GetExpectedSigma(mom,time,AliPID::kDeuteron);
		Float_t resolutionTPC_deut = fPIDResponse->GetTPCResponse().GetExpectedSigma(mom_tpc,track->GetTPCsignalN(),AliPID::kDeuteron);		
		


//		    cout<<resolutionTPC_deut<<" ";
//		    Float_t resolutionTOF_deut = fPIDResponse->GetTOFResponse().GetExpectedSigma(mom,time_of_flight,AliPID::kDeuteron);

//		    cout<<resolutionTOF_deut<<endl;
		

		
		time = time * 0.001;

//		dedx_exp_deut->Fill(time_exp_deut);
//		dedx_exp_deut->Fill(resolutionTOF_deut);
//		if(time > time_exp_prot)   pid_cut = 1;                                                                                  // worked poorly, admitted lots of spurious tracks
////////////	if(TMath::Abs(time - time_exp_deut) < resolutionTOF_deut*0.001*5)   pid_cut = 1;                                         // works well, cutting off sides of deuteron

		// mark A


//		if(TMath::Abs(time - time_exp_deut) < resolutionTOF_deut*0.001*5)
		if(time > time_exp_prot  &&  dedx > (dedx_exp_prot + resolutionTPC_prot*5))
		{
		    h_dedx_exp_deut->Fill(dedx_exp_deut);
		    h_dedx_exp_prot->Fill(dedx_exp_prot);
		    h_dedx->Fill(dedx);
		    h_TPC_resolution_deut->Fill(resolutionTPC_deut);
		    h_TPC_resolution_prot->Fill(resolutionTPC_prot);
		    pid_cut = 1;
		}
		
//		if(NEvents%500 == 0)
//		{
		    
//		}
		
//		if(dedx>(dedx_exp_prot+resolutionTPC_prot*5))   pid_cut = 1;  // works, but limited pt range
//		if(dedx>(dedx_exp_prot+resolutionTPC_prot*5)  &&  TMath::Abs(time_exp_deut - time) < TMath::Abs(time_exp_prot - time))   pid_cut = 1;
//		if(TMath::Abs(time_exp_deut - time/1000) < TMath::Abs(time_exp_prot - time/1000))   pid_cut = 1;
//		if(exptimes[4] < time  &&  dedx > (dedx_exp_prot + resolutionTPC_prot*2)) pid_cut = 1;  // worked, to low momentum deuterons below about 2.5 GeV/c
//		if(exptimes[4]<time  &&  dedx>(dedxExp_deut-3*resolutionTPC_deut)  &&  dedx<(dedxExp_deut+3*resolutionTPC_deut)  &&  dedx>(dedx_exp_prot+resolutionTPC_prot*2))    pid_cut = 1;
	    }
	    
	    
	    if(charge > 0)
	    {
		m_squared_pos->Fill(mom,m2tof);
		plength_vs_mom_pos->Fill(mom, length);
		ttof_vs_mom_pos->Fill(   mom, time_of_flight);
		beta_vs_mom_pos->Fill(   mom, Beta(track));
		deltat_vs_mom_pos->Fill( mom, tof_minus_tpion(track));

		if(pid_cut > 0)
		{
		    m_squared_pos_cut->Fill(mom,m2tof);
		    plength_vs_mom_pos_cut->Fill(mom, length);
		    ttof_vs_mom_pos_cut->Fill(   mom, time_of_flight);
		    beta_vs_mom_pos_cut->Fill(   mom, Beta(track));
		    deltat_vs_mom_pos_cut->Fill( mom, tof_minus_tpion(track));
		}
	    }
	    else if(charge < 0)
	    {
		m_squared_neg->Fill(mom,m2tof);
		plength_vs_mom_neg->Fill(mom, length);
		ttof_vs_mom_neg->Fill(   mom, time_of_flight);
		beta_vs_mom_neg->Fill(   mom, Beta(track));
		deltat_vs_mom_neg->Fill( mom, tof_minus_tpion(track));

		if(pid_cut > 0)
		{
		    m_squared_neg_cut->Fill(mom,m2tof);
		    plength_vs_mom_neg_cut->Fill(mom, length);
		    ttof_vs_mom_neg_cut->Fill(   mom, time_of_flight);
		    beta_vs_mom_neg_cut->Fill(   mom, Beta(track));
		    deltat_vs_mom_neg_cut->Fill( mom, tof_minus_tpion(track));
		}
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
    
    m_squared = mom*mom * (tof*tof * c2 * c2 / (length * length) - 1);
    return m_squared;

}
