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
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
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




class AliAnalysisTaskCorPIDTOFQA;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'


//ofstream fout;
Int_t NEvents = 0;

ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fHistPt(0), /* fEMCal_cluster(0), fEMCal_phi(0), fEMCal_eta(0), fEMCal_eta_phi(0), fEMCal_E(0), fEMCal_P(0), fEMCal_E_over_P(0), fEMCal_E_P(0), fEMCal_M(0), pion_pt(0), kaon_pt(0), prot_pt(0), kaon_pion(0), prot_pion(0),*/ cent_ntracks(0), t_minus_tpion(0), m_squared(0), species_num(0), /* tof_lit_eta(0), tof_big_eta(0), TPCrows(0), DCA_XY(0),*/ dtof_dEdx(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
								     fAOD(0), fOutputList(0), fPIDResponse(0), fHistPt(0), /* fEMCal_cluster(0), fEMCal_phi(0), fEMCal_eta(0), fEMCal_eta_phi(0), fEMCal_E(0), fEMCal_P(0), fEMCal_E_over_P(0), fEMCal_E_P(0), fEMCal_M(0), pion_pt(0), kaon_pt(0), prot_pt(0), kaon_pion(0), prot_pion(0), */ cent_ntracks(0), t_minus_tpion(0), m_squared(0), species_num(0), /*tof_lit_eta(0), tof_big_eta(0), TPCrows(0), DCA_XY(0),*/ dtof_dEdx(0)
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

//    fout.open("output.txt");

    
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    // example of a histogram

    fHistPt             = new TH1F("fHistPt",           "Pt()",              100,     0,    10);
    cent_ntracks        = new TH2F("cent_vs_ntracks",   "cent_vs_ntracks",   100,     0,   100,     100,     0,   800);
    t_minus_tpion       = new TH1F("t_minus_tpion",     "t_minus_tpion",     900,  -2.0,   7.0);
    m_squared           = new TH2F("m_squared",         "m_squared",         320,   0.0,   8.0,     350,  -1.0,   6.0);
    species_num         = new TH1F("species_num",       "species_num",        14,     0,    14);
    beta_vs_p           = new TH2F("beta_vs_p",         "beta_vs_p",         320,   0.0,   8.0,     300,   0.1,   1.1);
    sigma_vs_p          = new TH2F("sigma_vs_p",        "sigma_vs_p",        320,   0.0,   8.0,     300,  -4.0,   4.0);
    dtof_dEdx           = new TH2F("dtof_dEdx",         "dtof_dEdx",         600, -30.0,  30.0,     900,  -2.0,   7.0);
    


    fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!

    fOutputList->Add(cent_ntracks);                                                                 //// added by Brennan
    fOutputList->Add(t_minus_tpion);                                                                //// added by Brennan
    fOutputList->Add(m_squared);                                                                    //// added by Brennan
    fOutputList->Add(species_num);                                                                  //// added by Brennan
    fOutputList->Add(beta_vs_p);                                                                    //// added by Brennan
    fOutputList->Add(sigma_vs_p);                                                                   //// added by Brennan
    fOutputList->Add(dtof_dEdx);                                                                    //// added by Brennan

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    AliAnalysisManager *man            = AliAnalysisManager::GetAnalysisManager();                  //// added by Brennan
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());     //// added by Brennan
    fPIDResponse                       = inputHandler->GetPIDResponse();                            //// added by Brennan

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    fESDpid->SetTOFResponse(fESD, AliESDpid::kBest_T0);  // found with 'grep'    
//    pidres->SetTOFResponse(ev,AliPIDResponse::kTOF_T0);


    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}







//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::UserExec(Option_t *)
{

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


//    fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kFILL_T0);  // set event-time method


    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event


// these made no difference
//  fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kTOF_T0);
//  fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kT0_T0);
//  fPIDResponse->SetTOFResponse(fAOD,AliPIDResponse::kBest_T0);


    
    
    //////////////////////////////////////// MULTIPLICITY PART ////////////////////////////////////////
/*    
    AliMultSelection *MultSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));

//    AliMultSelection *MultSelection=0;
    MultSelection=(AliMultSelection*) fAOD->FindListObject("MultSelection");

//  if(!MultSelection) return kTRUE;
//  if(MultSelection->IsEventSelected()) return kTRUE;
//  return kFALSE;
    
//  AliMultSelectionTask * MultSelection = AddTaskMultSelection(kFALSE); // user mode:    
//  use the default calibration for runs which have not yet been calibrated
//  MultSelection->SetUseDefaultCalib(kTRUE); // data
    
    //output to file to check if calibration existed in this dataset
    if (MultSelection) {   fout<<"MULTSELECTION OBJECT WORKED"<<endl;   }
    else               {   fout<<"MULTSELECTION OBJECT didn't WORK"<<endl;   }
*/  
    cent_ntracks->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"), iTracks);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////




    
//    cout<<fAOD->GetCentrality()->GetCentralityPercentile("V0M")<<" "<<iTracks<<endl;
    
    for(Int_t i(0); i < iTracks; i++)
    {                 // loop ove rall these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                                                       // if we failed, skip this track
	if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { continue; }
//	DCA_XY->Fill(track->XAtDCA(), track->YAtDCA());


	Double_t mom     = track->P();
	Double_t energy  = track->E();
	float    mass    = track->M();
	float    pt      = track->Pt();
	Double_t tpc_mom = track->GetTPCmomentum();


	fHistPt->Fill(track->Pt());                                                // plot the pt value of the track in a histogram
	/////////////////////////////////////////////////

	
	Double_t nsigmaTPC = 999.0;
	Double_t nsigmaTOF = 999.0;

	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);/* && trk->IsOn(AliESDtrack::kTPCpid)*/;
	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);

	Bool_t fHasTPCPID=kFALSE;
	Bool_t fHasTOFPID=kFALSE;


	Double_t m2tof = get_mass_squared(track);

	if(tpcIsOk  &&  tofIsOk)
	{
	    beta_vs_p->Fill(mom, Beta(track));
	    Double_t sigma_min = -999.0;
	    int      id        = 14;
	    
	    for(int i=0; i<15; i++)
	    {
		AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) i, nsigmaTPC);
//		AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType) i, nsigmaTOF);
//		AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) i, nsigmaTOF);

//		nsigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType) i);
		if(fabs(nsigmaTPC) < fabs(sigma_min))
		{
		    sigma_min = nsigmaTPC;
		    id = i;
		}
	    }
	    species_num->Fill(id);
	    sigma_vs_p->Fill(mom, sigma_min);
	    m_squared->Fill(mom,m2tof);
	    
	    if(mom > 2.90  &&  mom < 3.00)
	    {
		Double_t pion_dEdx = fPIDResponse->GetTPCResponse().GetExpectedSignal(tpc_mom,AliPID::kPion);
		dtof_dEdx->Fill(track->GetTPCsignal() - pion_dEdx, tof_minus_tpion(track));
	    }


	    if(mom > 0.95  &&  mom < 1.05)
	    {
		t_minus_tpion->Fill(tof_minus_tpion(track));
		
	    }	    
	}	
	

	
	
//	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
//	if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
	
    }                                                   // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file

    NEvents++;

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
    Double_t stoptime      = track->GetTOFsignal();
    Double_t c             = TMath::C()*1.E-9;// m/ns
    Float_t startTime      = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());//in ps
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
//    Double_t start_time     = fPIDResponse->GetTOFResponse().GetTimeZero() * 1e-12;
    Double_t stop_time      = track->GetTOFsignal();
    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
    Double_t time_of_flight = stop_time - start_time;
    time_of_flight          = time_of_flight * 0.001;                                                       // convert from ps to ns
    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
    length = length * 100;                                                                                  // convert to cm
    Double_t pion_time      = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion);        // doesn't work, returns essentially 0
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
