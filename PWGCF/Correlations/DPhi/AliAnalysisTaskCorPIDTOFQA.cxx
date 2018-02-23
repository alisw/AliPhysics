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

#include <fstream>
#include <cmath>

#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TList.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
//#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCorPIDTOFQA.h"
#include "AliPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAnalysisUtils.h"
//#include "AliEmcalTrackSelection.h"
//#include "AliEmcalTrackSelectionAOD.h"
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
//using namespace BSchaefer_devel;


ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),

    htree(0),
    primary_vertex_z(0),          //  E  1 (event)
    primary_vertex_z_cut(0),      //  E  2
    multiplicity_hybrid(0),
    multiplicity_global(0),
    deut_per_event(0),            //  E  3
    fHistPt(0),                   //  T  4 (track)
    m2_pt_pos(0),                 //  T  5
    m2_pt_neg(0),                 //  T  6
    m2_pt_pos_cut(0),             //  T  7
    m2_pt_neg_cut(0),             //  T  8
    m2_pt_pos_fine(0),            //  T  9
    m2_pt_neg_fine(0),            //  T 10
    m2_pt_pos_cut_fine(0),        //  T 11
    m2_pt_neg_cut_fine(0),        //  T 12
    trig_03_phi_pt_pos(0),         // T 13
    trig_03_phi_pt_neg(0),         // T 14  
    trig_05_phi_pt_pos(0),         // T 15
    trig_05_phi_pt_neg(0),         // T 16
    trig_08_phi_pt_pos(0),         // T 17
    trig_08_phi_pt_neg(0)          // T 18
    
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),
					 
    htree(0),
    primary_vertex_z(0),          //  E  1 (event)
    primary_vertex_z_cut(0),      //  E  2
    multiplicity_hybrid(0),
    multiplicity_global(0),
    deut_per_event(0),            //  E  3
    fHistPt(0),                   //  T  4 (track)
    m2_pt_pos(0),                 //  T  5
    m2_pt_neg(0),                 //  T  6
    m2_pt_pos_cut(0),             //  T  7
    m2_pt_neg_cut(0),             //  T  8
    m2_pt_pos_fine(0),            //  T  9
    m2_pt_neg_fine(0),            //  T 10
    m2_pt_pos_cut_fine(0),        //  T 11
    m2_pt_neg_cut_fine(0),        //  T 12
    trig_03_phi_pt_pos(0),         // T 13
    trig_03_phi_pt_neg(0),         // T 14  
    trig_05_phi_pt_pos(0),         // T 15
    trig_05_phi_pt_neg(0),         // T 16
    trig_08_phi_pt_pos(0),         // T 17
    trig_08_phi_pt_neg(0)          // T 18									   
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
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



    
    fAnalysisUtils = new AliAnalysisUtils;
    
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

    Double_t pt_binning[2001];

    Float_t moving_marker = 0.10;
    for(int i=0; i<1602; i++)
    {
	pt_binning[i] = moving_marker;
	moving_marker = moving_marker + pt_binning[i] * 0.005;
    }

    float lower = -pio2;
    float upper = 3.0*pio2;

    htree                      = new TTree("htree", "tree containing mixed event pt, phi, eta");
    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            100,    -50.0,   50.0);                              //  1
    primary_vertex_z_cut       = new TH1F("primary_vertex_z_cut",       "primary_vertex_z_cut",        100,    -50.0,   50.0);                              //  2
    multiplicity_hybrid        = new TH1I("multiplicity_hybrid",        "multiplicity_hybrid",         150,      0,    150);
    multiplicity_global        = new TH1I("multiplicity_global",        "multiplicity_global",         150,      0,    150);
    deut_per_event             = new TH1I("deut_per_event",             "deut_per_event",               12,        0,     12);                              //  3
    fHistPt                    = new TH1F("fHistPt",                    "Pt()",                       1000,       0.0, 20.0);                               //  4
    m2_pt_pos                  = new TH2F("m2_pt_pos",                  "m2_pt_pos",                    50,         0.0, 5.0,    2400,    -1.0,     7.0);   //  5
    m2_pt_neg                  = new TH2F("m2_pt_neg",                  "m2_pt_neg",                    50,         0.0, 5.0,    2400,    -1.0,     7.0);   //  6
    m2_pt_pos_cut              = new TH2F("m2_pt_pos_cut",              "m2_pt_pos_cut",                50,         0.0, 5.0,    2400,    -1.0,     7.0);   //  7
    m2_pt_neg_cut              = new TH2F("m2_pt_neg_cut",              "m2_pt_neg_cut",                50,         0.0, 5.0,    2400,    -1.0,     7.0);   //  8
    m2_pt_pos_fine             = new TH2F("m2_pt_pos_fine",             "m2_pt_pos_fine",              800,       pt_binning,    2400,    -1.0,     7.0);   //  9
    m2_pt_neg_fine             = new TH2F("m2_pt_neg_fine",             "m2_pt_neg_fine",              800,       pt_binning,    2400,    -1.0,     7.0);   // 10
    m2_pt_pos_cut_fine         = new TH2F("m2_pt_pos_cut_fine",         "m2_pt_pos_cut_fine",          800,       pt_binning,    2400,    -1.0,     7.0);   // 11
    m2_pt_neg_cut_fine         = new TH2F("m2_pt_neg_cut_fine",         "m2_pt_neg_cut_fine",          800,       pt_binning,    2400,    -1.0,     7.0);   // 12
    trig_03_phi_pt_pos         = new TH2F("trig_03_phi_pt_pos",         "trig_03_phi_pt_pos",          170,  3.0, 20.0,    288,   lower,   upper);          // 13
    trig_03_phi_pt_neg         = new TH2F("trig_03_phi_pt_neg",         "trig_03_phi_pt_neg",          170,  3.0, 20.0,    288,   lower,   upper);          // 14
    trig_05_phi_pt_pos         = new TH2F("trig_05_phi_pt_pos",         "trig_05_phi_pt_pos",          170,  3.0, 20.0,    288,   lower,   upper);          // 15
    trig_05_phi_pt_neg         = new TH2F("trig_05_phi_pt_neg",         "trig_05_phi_pt_neg",          170,  3.0, 20.0,    288,   lower,   upper);          // 16
    trig_08_phi_pt_pos         = new TH2F("trig_08_phi_pt_pos",         "trig_08_phi_pt_pos",          170,  3.0, 20.0,    288,   lower,   upper);          // 17
    trig_08_phi_pt_neg         = new TH2F("trig_08_phi_pt_neg",         "trig_08_phi_pt_neg",          170,  3.0, 20.0,    288,   lower,   upper);          // 18

    htree->Branch("mult",              &d_mult,                          "mult/B");         // E-1 event info
    htree->Branch("zvert",             &d_zvert,                        "zvert/B");         // E-2
    htree->Branch("assoc_count",       &d_assoc_count,            "assoc_count/S");         // E-3
    htree->Branch("assoc_pt",          &d_assoc_pt,     "assoc_pt[assoc_count]/F");         // T-1 track info
    htree->Branch("assoc_phi",         &d_assoc_phi,   "assoc_phi[assoc_count]/F");         // T-2
    htree->Branch("assoc_eta",         &d_assoc_eta,   "assoc_eta[assoc_count]/F");         // T-3
    htree->Branch("trigg_count",       &d_trigg_count,            "trigg_count/S");         // E-4
    htree->Branch("trigg_pt",          &d_trigg_pt,     "trigg_pt[trigg_count]/F");         // T-4 track info
    htree->Branch("trigg_phi",         &d_trigg_phi,   "trigg_phi[trigg_count]/F");         // T-5
    htree->Branch("trigg_eta",         &d_trigg_eta,   "trigg_eta[trigg_count]/F");         // T-6

    fOutputList->Add(htree);
    fOutputList->Add(primary_vertex_z);            //  1
    fOutputList->Add(primary_vertex_z_cut);        //  2
    fOutputList->Add(multiplicity_hybrid);
    fOutputList->Add(multiplicity_global);
    fOutputList->Add(deut_per_event);              //  3
    fOutputList->Add(fHistPt);                     //  4
    fOutputList->Add(m2_pt_pos);                   //  5
    fOutputList->Add(m2_pt_neg);                   //  6
    fOutputList->Add(m2_pt_pos_cut);               //  7
    fOutputList->Add(m2_pt_neg_cut);               //  8
    fOutputList->Add(m2_pt_pos_fine);              //  9
    fOutputList->Add(m2_pt_neg_fine);              // 10
    fOutputList->Add(m2_pt_pos_cut_fine);          // 11
    fOutputList->Add(m2_pt_neg_cut_fine);          // 12
    fOutputList->Add(trig_03_phi_pt_pos);          // 13
    fOutputList->Add(trig_03_phi_pt_neg);          // 14
    fOutputList->Add(trig_05_phi_pt_pos);          // 15
    fOutputList->Add(trig_05_phi_pt_neg);          // 16
    fOutputList->Add(trig_08_phi_pt_pos);          // 17
    fOutputList->Add(trig_08_phi_pt_neg);          // 18

    
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

    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7)   ;
    else return;

    const AliAODVertex *primVertex = fAOD->GetPrimaryVertex();
    Double_t pv = primVertex->GetZ();
    
    primary_vertex_z->Fill(pv);
//  if(fUtils−>IsFirstEventInChunk(fAOD))            return;  // not needed for new reconstructions
    if(!fAnalysisUtils->IsVertexSelected2013pA(fAOD)) return;
//  fAnalysisUtils−>SetCutOnZVertexSPD(kFALSE);
//  fAnalysisUtils->SetCutOnZVertexSPD(0);
    if(fAnalysisUtils->IsPileUpSPD(fAOD)) return;

    d_zvert = 0;
    primary_vertex_z_cut->Fill(pv);
    if(pv < 0.0){   d_zvert = 1;   }
//    else        {   d_zvert = 0;   }
    Int_t iTracks(fAOD->GetNumberOfTracks());

    int deut_count       = 0;
    int hybrid_count     = 0;
    int global_count     = 0;

    int trig_track_count = 0;
    Float_t max_pt       = 0.0;
    Int_t   max          = 0;

    short deut_candidate_event   = 0;
    
//  int trig_03_track_count     = 0;

    d_trigg_count        = 0;
    d_assoc_count        = 0;
    // loop over all these tracks
    //
    // pull out track numbers for high-pt triggers and also deutons

   
    for(Int_t i(0); i < iTracks; i++)
    {

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track)                                                                      {    continue;    }

	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }  // cut 1
	hybrid_count++;
	Float_t pt     = track->Pt();   	   if(pt   < 0.95)                      {    continue;    }  // cut 2
	Float_t dedx   = track->GetTPCsignal();    if(dedx > 1000)                      {    continue;    }  // cut 3
	Float_t eta    = track->Eta();	if(TMath::Abs(eta) > 0.90)                      {    continue;    }  // cut 4

	int is_global = 0;
	if(!track->IsGlobalConstrained())                                               {    is_global=1; }

	if(is_global==1)   global_count++;

	Float_t phi           = track->Phi();
	if(phi <  -pio2){    phi = phi + twopi; }	if(phi <  -pio2){    phi = phi + twopi;   }
	if(phi > 3*pio2){    phi = phi - twopi;	}       if(phi > 3*pio2){    phi = phi - twopi;   }

	if(pt >= 3.0)
	{
//	    trig_track_num[trig_track_count] = i;
//	    skimmer_trigger[trig_track_count][0] = pt;
//	    skimmer_trigger[trig_track_count][1] = phi;
//	    skimmer_trigger[trig_track_count][2] = eta;
	    if(pt>max_pt){   max_pt = pt; max = trig_track_count;   }
	    trig_track_count++;

	    d_trigg_pt[d_trigg_count]  = pt;
	    d_trigg_phi[d_trigg_count] = phi;
	    d_trigg_eta[d_trigg_count] = eta;
	    d_trigg_count++;
	    
	}

	

	Short_t charge        = track->Charge();
	if(pt >= 3.0)	{   if(charge > 0)   trig_03_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_03_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 5.0)	{   if(charge > 0)   trig_05_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_05_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 8.0)	{   if(charge > 0)   trig_08_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_08_phi_pt_neg->Fill(pt, phi);	}

	Double_t nsigmaTPC = 999.0;	Double_t nsigmaTOF = 999.0;
	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk); 
	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
	if(!tpcIsOk)	                                                                {    continue;    }
	fHistPt->Fill(pt);
	if(!tofIsOk)	                                                                {    continue;    }

	Float_t m2tof       = get_mass_squared(track);
	Float_t deut_mean   = 0.0;
	Float_t deut_sigma  = 0.0;


	Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron

	if(TMath::Abs(nSigmaTPCDeut) < 3.0)
	{
	    if(charge > 0)
	    {
		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }	    deut_mean  = fit_deut_curve->Eval(pt);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }	    deut_sigma = fit_deut_curve->Eval(pt);
		    if(m2tof < deut_mean + cut_width * deut_sigma  &&  m2tof > deut_mean - cut_width * deut_sigma){	deut_candidate_event = 1;    }
		}
	    }
	    else if(charge < 0)
	    {
		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][0][w]);   }	    deut_mean  = fit_deut_curve->Eval(pt);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][1][w]);   }	    deut_sigma = fit_deut_curve->Eval(pt);
		    if(m2tof < deut_mean + cut_width * deut_sigma  &&  m2tof > deut_mean - cut_width * deut_sigma){	deut_candidate_event = 1;    }
		}
	    }
	}
	
	
	if(pt >= 1.0  &&  pt < 4.4  &&  deut_candidate_event == 0  &&  m2tof < 1.7)
	{
	    d_assoc_pt[d_assoc_count]  = pt;
	    d_assoc_phi[d_assoc_count] = phi;
	    d_assoc_eta[d_assoc_count] = eta;
	    d_assoc_count++;
	}
    }        //   end of track loop


    deut_per_event->Fill(deut_count);
    multiplicity_hybrid->Fill(hybrid_count);
    multiplicity_global->Fill(global_count);

         if(hybrid_count <  16){   d_mult = 0;   }
    else if(hybrid_count <  23){   d_mult = 1;   }
    else if(hybrid_count <  31){   d_mult = 2;   }
    else if(hybrid_count <  42){   d_mult = 3;   }
    else if(hybrid_count < 149){   d_mult = 4;   }

    if((d_assoc_count > 0  ||  d_trigg_count > 0)  &&  deut_candidate_event == 0)
    {
	htree->Fill();
    }

    
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
    Double_t m2             = 0.0;
    Double_t c2             = 29.9792458;
    
    m2 = pow(mom,2) * (tof*tof * c2 * c2 / (length * length) - 1);
    return m2;

}
