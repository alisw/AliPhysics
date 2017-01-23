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
Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]
TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);

ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fHistPt(0),  cent_ntracks(0), m_squared_pos_raw(0), m_squared_pos_cut(0), m_squared_pos(0), m_squared_neg_raw(0), m_squared_neg_cut(0), m_squared_neg(0), dtof_dEdx(0), mapped_ttof(0), m_squared_pos_deut(0), m_squared_neg_deut(0), plength_vs_mom(0), ttof_vs_mom(0), deltat_vs_mom(0), beta_vs_mom(0)

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
									   fAOD(0), fOutputList(0), fPIDResponse(0), fHistPt(0), cent_ntracks(0), m_squared_pos_raw(0), m_squared_pos_cut(0), m_squared_pos(0), m_squared_neg_raw(0),m_squared_neg_cut(0), m_squared_neg(0), dtof_dEdx(0), mapped_ttof(0), m_squared_pos_deut(0), m_squared_neg_deut(0), plength_vs_mom(0), ttof_vs_mom(0), deltat_vs_mom(0), beta_vs_mom(0)
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

// 3.38835 0.0203511 0.101903 -0.0487972 0.0513143 0.0815525 3.69323 -0.0206797 -0.200752 0.00644255 0.0391803 0.0399812
    
//  Double_t deut_curves[2][2][3];  // [charge][mean,sigma][par]

    deut_curves[0][0][0] = 3.38835;     // pos deut mean curve
    deut_curves[0][0][1] = 0.020351;
    deut_curves[0][0][2] = 0.101903;

    deut_curves[0][1][0] = -0.0487972;  // pos deut sigma curve
    deut_curves[0][1][1] = 0.0513143;
    deut_curves[0][1][2] = 0.0815525;

    deut_curves[1][0][0] = 3.69323;     // neg deut mean curve
    deut_curves[1][0][1] = -0.0206797;
    deut_curves[1][0][2] = -0.200752;

    deut_curves[1][1][0] = 0.00644255;  // neg deut sigma curve
    deut_curves[1][1][1] = 0.0391803;
    deut_curves[1][1][2] = 0.0399812;
    
    
//  TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) ",       1.1, 4.4);
    fit_deut_curve->SetParNames("a", "b x", "c/#sqrt(x)");
	
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    fHistPt             = new TH1F("fHistPt",           "Pt()",              100,     0,    10);
    cent_ntracks        = new TH2F("cent_vs_ntracks",   "cent_vs_ntracks",   100,     0,   100,     100,     0,   800);
    m_squared_pos       = new TH2F("m_squared_pos",     "m_squared_pos",     320,   0.0,   8.0,    1500,  -1.0,   5.0);
    m_squared_pos_raw   = new TH2F("m_squared_pos_raw", "m_squared_pos_raw", 320,   0.0,   8.0,    1500,  -1.0,   5.0);
    m_squared_pos_cut   = new TH2F("m_squared_pos_cut", "m_squared_pos_cut", 320,   0.0,   8.0,    1500,  -1.0,   5.0);
    m_squared_neg       = new TH2F("m_squared_neg",     "m_squared_neg",     320,   0.0,   8.0,    1500,  -1.0,   5.0);
    m_squared_neg_raw   = new TH2F("m_squared_neg_raw", "m_squared_neg_raw", 320,   0.0,   8.0,    1500,  -1.0,   5.0);
    m_squared_neg_cut   = new TH2F("m_squared_neg_cut", "m_squared_neg_cut", 320,   0.0,   8.0,    1500,  -1.0,   5.0);
    beta_vs_mom         = new TH2F("beta_vs_mom",       "beta_vs_mom",       780,   0.0,   8.0,    1500,   0.1,   1.1);
    dtof_dEdx           = new TH2F("dtof_dEdx",         "dtof_dEdx",         600, -30.0,  30.0,     900,  -2.0,   7.0);
    mapped_ttof         = new TH1F("mapped_ttof",       "mapped_ttof",       500,  -0.0,   2.0);
    m_squared_pos_deut  = new TH2F("m_squared_pos_deut","m_squared_pos_deut",320,   0.0,   8.0,    1500,  -1.0,   5.0);
    m_squared_neg_deut  = new TH2F("m_squared_neg_deut","m_squared_neg_deut",320,   0.0,   8.0,    1500,  -1.0,   5.0);
    plength_vs_mom      = new TH2F("plength_vs_mom",    "plength_vs_mom",    780,   0.2,   8.0,    4000, 350.0,   750);  // (GeV/c 10MeV/bin,   cm  1mm/bin)
    ttof_vs_mom         = new TH2F("ttof_vs_mom",       "ttof_vs_mom",       780,   0.2,   8.0,    4000,  10.0,  50.0);  // (GeV/c 10MeV/bin,   ns 10ps/bin)
    deltat_vs_mom       = new TH2F("deltat_vs_mom",     "deltat_vs_mom",     780,   0.2,   8.0,    2100,  -1.0,  20.0);

    fOutputList->Add(fHistPt);                                                                      // objects added to output file
    fOutputList->Add(cent_ntracks);
    fOutputList->Add(m_squared_pos);
    fOutputList->Add(m_squared_pos_raw);
    fOutputList->Add(m_squared_pos_cut);
    fOutputList->Add(m_squared_neg);    
    fOutputList->Add(m_squared_neg_raw);
    fOutputList->Add(m_squared_neg_cut);

    fOutputList->Add(dtof_dEdx);        
    fOutputList->Add(mapped_ttof);      
    fOutputList->Add(m_squared_pos_deut);
    fOutputList->Add(m_squared_neg_deut);
    fOutputList->Add(ttof_vs_mom);       
    fOutputList->Add(beta_vs_mom);       
    fOutputList->Add(plength_vs_mom);    
    fOutputList->Add(deltat_vs_mom);     
    
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




    

    // loop over all these tracks
    for(Int_t i(0); i < iTracks; i++)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                                                       // if we failed, skip this track
	if(!(track->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { continue; }
//	DCA_XY->Fill(track->XAtDCA(), track->YAtDCA());


	Double_t mom     = track->P();
	Double_t energy  = track->E();
	float    mass    = track->M();
	float    pt      = track->Pt();
	Double_t tpc_mom = track->GetTPCmomentum();
	Short_t  charge  = track->Charge();
	Double_t eta     = track->Eta();

	fHistPt->Fill(track->Pt());                                                // plot the pt value of the track in a histogram
	/////////////////////////////////////////////////

	
	Double_t nsigmaTPC = 999.0;
	Double_t nsigmaTOF = 999.0;

	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);/* && trk->IsOn(AliESDtrack::kTPCpid)*/;
	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);


	Double_t m2tof = get_mass_squared(track);

	if(tpcIsOk  &&  tofIsOk)
	{
	    beta_vs_mom->Fill(mom, Beta(track));
	    Double_t sigma_min = -999.0;
	    int      id        = 14;
	    
	    Double_t deut_mean  = 0.0;
	    Double_t deut_sigma = 0.0;
	    Double_t cut_width  = 2.0;
	    
	    if(charge > 0)
	    {
		m_squared_pos_raw->Fill(mom,m2tof);

		if(mom >= 1.1  &&  mom < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(mom);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }
		    deut_sigma = fit_deut_curve->Eval(mom);

		    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
		    {   m_squared_pos_cut->Fill(mom,m2tof);   }
		}
	    }
	    else if(charge < 0)
	    {
		m_squared_neg_raw->Fill(mom,m2tof);
		if(mom >= 1.1  &&  mom < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][0][w]);   }
		    deut_mean = fit_deut_curve->Eval(mom);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][1][w]);   }
		    deut_sigma = fit_deut_curve->Eval(mom);


		    if(m2tof < deut_mean + cut_width * deut_sigma  &&   m2tof > deut_mean - cut_width * deut_sigma)
		    {   m_squared_neg_cut->Fill(mom,m2tof);   }
		}
	    }


	    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,(AliPID::EParticleType) 2, nsigmaTPC);
	    if(fabs(nsigmaTPC) > 2.0)  // is it a pion?
	    {
		AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,(AliPID::EParticleType) 3, nsigmaTPC);
		if(fabs(nsigmaTPC) > 2.0)  // is it a kaon?
		{
		    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,(AliPID::EParticleType) 4, nsigmaTPC);
		    if(fabs(nsigmaTPC) > 2.0)  // is it a proton?
		    {
			AliPIDResponse::EDetPidStatus statusTPC=fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,(AliPID::EParticleType)1,nsigmaTPC);
			if(fabs(nsigmaTPC) > 2.0)  // is it a muon?
			{
			    AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,(AliPID::EParticleType) 0, nsigmaTPC);
			    if(fabs(nsigmaTPC) > 2.0)  // is it an electron?
			    {
				if(charge > 0)
				{
				    m_squared_pos->Fill(mom,m2tof);
				}
				else if(charge < 0)
				{
				    m_squared_neg->Fill(mom,m2tof);
				}
			    }
			}
		    }
		}
	    }


	    AliPIDResponse::EDetPidStatus statusTPC2 = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,(AliPID::EParticleType) 5, nsigmaTPC);
	    if(fabs(nsigmaTPC) < 3.0)  // is it a deuteron?
	    {
		if(charge > 0)
		{
		    m_squared_pos_deut->Fill(mom,m2tof);
		}
		else if(charge < 0)
		{
		    m_squared_neg_deut->Fill(mom,m2tof);
		}
	    }
			


	    
	    if(mom > 2.90  &&  mom < 3.00)
	    {
		Double_t pion_dEdx = fPIDResponse->GetTPCResponse().GetExpectedSignal(tpc_mom,AliPID::kPion);
		dtof_dEdx->Fill(track->GetTPCsignal() - pion_dEdx, tof_minus_tpion(track));
	    }


	    Double_t start_time     = fPIDResponse->GetTOFResponse().GetStartTime(track->P());                      // in ps
	    Double_t stop_time      = track->GetTOFsignal();
	    Double_t time_of_flight = stop_time - start_time;
	    time_of_flight          = time_of_flight * 0.001;                                                       // convert ps to ns
	    Double_t c              = TMath::C()*1.E-9;                                                             // in m/ns
	    Double_t length         = 0.0;
	    length                  = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
	    length                  = length * 100;                                                                 // convert mm to cm
	    plength_vs_mom->Fill(mom, length);
	    ttof_vs_mom->Fill(mom, time_of_flight);
	    deltat_vs_mom->Fill(mom, tof_minus_tpion(track));

	    
	    if(mom > 0.95  &&  mom < 1.05)
	    {
		Double_t tmtp = 0.0;
		tmtp = tof_minus_tpion(track);

		if(tmtp > 0.8  &&  tmtp < 2.4)  // select only kaons
		{
                                                      // convert from ps to ns
//		    Double_t length         = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion)*1E-3*c; // in meters
//		    length = length * 100;                                                                                  // convert to cm

//		    Double_t mapped = 0.0;
//		    mapped = m2tof*length / sqrt(time_of_flight*time_of_flight * c*c - length*length);
		    c = 29.97925;
		    mapped_ttof->Fill(0.493677*length / sqrt(time_of_flight*time_of_flight * c*c - length*length));
		    
		}
	    }	    
	}	
	

    }                                                   // continue until all the tracks are processed
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
