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
//#include <iomanip>
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
//using namespace BSchaefer_devel;

//ofstream file("output.txt");
//const int multiplicity_cut       = 0;
//const int z_vertex_cut           = 0;


ClassImp(AliAnalysisTaskCorPIDTOFQA) // classimp: necessary for root

AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),

    
    primary_vertex_z(0),          //  E  1 (event)
    primary_vertex_z_cut(0),      //  E  2
    deut_per_event(0),            //  E  3

    multiplicity_hybrid(0),
    multiplicity_global(0),
    
    fHistPt(0),                   //  T  4 (track)
    
    trig_03_phi_pt_pos(0),         // T 13
    trig_03_phi_pt_neg(0),         // T 14  
    trig_05_phi_pt_pos(0),         // T 15
    trig_05_phi_pt_neg(0),         // T 16
    trig_08_phi_pt_pos(0),         // T 17
    trig_08_phi_pt_neg(0)          // T 18

//    associates(0),
//    triggers(0)
    
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorPIDTOFQA::AliAnalysisTaskCorPIDTOFQA(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fPIDResponse(0), fAnalysisUtils(0),


    primary_vertex_z(0),          //  E  1 (event)
    primary_vertex_z_cut(0),      //  E  2
    deut_per_event(0),            //  E  3

    multiplicity_hybrid(0),
    multiplicity_global(0),

    fHistPt(0),                   //  T  4 (track)
    
    trig_03_phi_pt_pos(0),        // T 13
    trig_03_phi_pt_neg(0),        // T 14  
    trig_05_phi_pt_pos(0),        // T 15
    trig_05_phi_pt_neg(0),        // T 16
    trig_08_phi_pt_pos(0),        // T 17
    trig_08_phi_pt_neg(0)         // T 18

//    associates(0),
//    triggers(0)

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

    primary_vertex_z           = new TH1F("primary_vertex_z",           "primary_vertex_z",            100, -50.0,  50.0);                               //  1
    primary_vertex_z_cut       = new TH1F("primary_vertex_z_cut",       "primary_vertex_z_cut",        100, -50.0,  50.0);                               //  2
    deut_per_event             = new TH1I("deut_per_event",             "deut_per_event",               12,     0,    12);                               //  3
    fHistPt                    = new TH1F("fHistPt",                    "Pt()",                       1000,   0.0,  20.0);                               //  4

    multiplicity_hybrid        = new TH1I("multiplicity_hybrid",        "multiplicity_hybrid",         150,     0,   150);
    multiplicity_global        = new TH1I("multiplicity_global",        "multiplicity_global",         150,     0,   150);
    
    trig_03_phi_pt_pos         = new TH2F("trig_03_phi_pt_pos",         "trig_03_phi_pt_pos",          170,   3.0,  20.0,   288, lower, upper);          // 13
    trig_03_phi_pt_neg         = new TH2F("trig_03_phi_pt_neg",         "trig_03_phi_pt_neg",          170,   3.0,  20.0,   288, lower, upper);          // 14
    trig_05_phi_pt_pos         = new TH2F("trig_05_phi_pt_pos",         "trig_05_phi_pt_pos",          170,   3.0,  20.0,   288, lower, upper);          // 15
    trig_05_phi_pt_neg         = new TH2F("trig_05_phi_pt_neg",         "trig_05_phi_pt_neg",          170,   3.0,  20.0,   288, lower, upper);          // 16
    trig_08_phi_pt_pos         = new TH2F("trig_08_phi_pt_pos",         "trig_08_phi_pt_pos",          170,   3.0,  20.0,   288, lower, upper);          // 17
    trig_08_phi_pt_neg         = new TH2F("trig_08_phi_pt_neg",         "trig_08_phi_pt_neg",          170,   3.0,  20.0,   288, lower, upper);          // 18

    associates                 = new std::ofstream("associates.txt");
    triggers                   = new std::ofstream("triggers.txt");



    
    char trigs[3]  = {'3', '5', '8'};
    char pts[4]    = {'1', '2', '3', '4'};
    char zverts[2] = {'+', '-'};
    char mults[5]  = {'1', '2', '3', '4', '5'};

    fOutputList->Add(primary_vertex_z);            //  1
    fOutputList->Add(primary_vertex_z_cut);        //  2
    fOutputList->Add(multiplicity_hybrid);
    fOutputList->Add(multiplicity_global);
    fOutputList->Add(deut_per_event);              //  3
    fOutputList->Add(fHistPt);                     //  4

    fOutputList->Add(trig_03_phi_pt_pos);          // 13
    fOutputList->Add(trig_03_phi_pt_neg);          // 14
    fOutputList->Add(trig_05_phi_pt_pos);          // 15
    fOutputList->Add(trig_05_phi_pt_neg);          // 16
    fOutputList->Add(trig_08_phi_pt_pos);          // 17
    fOutputList->Add(trig_08_phi_pt_neg);          // 18



//    file.open("output.txt");
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


    int i_trig  = 0;
    int i_pt    = 0;
    int i_zvert = 0;
    int i_mult  = 0;
    int i_phys  = 0;

    
    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7);
    else return;

    const AliAODVertex *primVertex = fAOD->GetPrimaryVertex();
    Double_t pv = primVertex->GetZ();
    
    primary_vertex_z->Fill(pv);

    if(pv < 0) i_zvert = 1;
//    if(fUtils−>IsFirstEventInChunk(fAOD))            return;  // not needed for new reconstructions

    if(!fAnalysisUtils->IsVertexSelected2013pA(fAOD)) return;
//    fAnalysisUtils−>SetCutOnZVertexSPD(kFALSE);
//    fAnalysisUtils->SetCutOnZVertexSPD(0);
    if(fAnalysisUtils->IsPileUpSPD(fAOD)) return;

    primary_vertex_z_cut->Fill(pv);
    


    
    Int_t iTracks(fAOD->GetNumberOfTracks());

    int non_deut_count   = 0;
    int hybrid_count     = 0;
    int global_count     = 0;

    int trigger_count    = 0;

//    int trig_track_num[200];
    int trig_track_count = 0;

    Float_t max_pt       = 0.0;
    Int_t   max          = 0;
    
    // loop over all these tracks
    //
    // pull out track numbers for high-pt triggers and also deutons



    Float_t skimmer_associate[500][3];
    Float_t skimmer_trigger[200][3];

    
    short deut_candidate_event   = 0;
    
    for(Int_t i(0); i < iTracks; i++)
    {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track)                                                                      {    continue;    }

	if(!(track->IsHybridGlobalConstrainedGlobal()))                                 {    continue;    }  // cut 1
	hybrid_count++;
	Float_t pt     = track->Pt();	           if(pt   < 0.95)                      {    continue;    }  // cut 2
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
	    skimmer_trigger[trig_track_count][0] = pt;
	    skimmer_trigger[trig_track_count][1] = phi;
	    skimmer_trigger[trig_track_count][2] = eta;
	    if(pt>max_pt){ max_pt = pt; max = trig_track_count;   }
	    trig_track_count++;
	}


	
	Short_t charge        = track->Charge();
	if(pt >= 3.0)	{   if(charge > 0)   trig_03_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_03_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 5.0)	{   if(charge > 0)   trig_05_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_05_phi_pt_neg->Fill(pt, phi);	}
	if(pt >= 8.0)	{   if(charge > 0)   trig_08_phi_pt_pos->Fill(pt, phi);	    else if(charge < 0)   trig_08_phi_pt_neg->Fill(pt, phi);	}
	
	Double_t nsigmaTPC = 999.0;	Double_t nsigmaTOF = 999.0;
	AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
	AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
	Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);     /* && trk->IsOn(AliESDtrack::kTPCpid)*/;
	Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
	if(!tpcIsOk)	                                                                {    continue;    }
	fHistPt->Fill(pt);
	if(!tofIsOk)	                                                                {    continue;    }

	Float_t m2tof          = get_mass_squared(track);
	Float_t deut_mean      = 0.0;
	Float_t deut_sigma     = 0.0;

	
	Double_t nSigmaTPCDeut = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);  // 5 = deuteron
	if(TMath::Abs(nSigmaTPCDeut) < 3.0)
	{
	    if(charge > 0)
	    {
		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][0][w]);   }	    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[0][1][w]);   }	    deut_sigma = fit_deut_curve->Eval(pt);
		    if(m2tof < deut_mean + cut_width * deut_sigma  &&  m2tof > deut_mean - cut_width * deut_sigma){	deut_candidate_event = 1;    }
		}
	    }
	    else if(charge < 0)
	    {
		if(pt >= 1.0  &&  pt < 4.4)
		{
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][0][w]);   }	    deut_mean = fit_deut_curve->Eval(pt);
		    for(int w=0; w<3; w++){   fit_deut_curve->SetParameter(w, deut_curves[1][1][w]);   }	    deut_sigma = fit_deut_curve->Eval(pt);
		    if(m2tof < deut_mean + cut_width * deut_sigma  &&  m2tof > deut_mean - cut_width * deut_sigma){	deut_candidate_event = 1;    }
		}
	    }
	}

	
	if(pt >= 1.0  &&  pt < 4.4  &&  deut_candidate_event == 0  &&  m2tof < 1.7)
	{
	    skimmer_associate[non_deut_count][0] = pt;
	    skimmer_associate[non_deut_count][1] = phi;
	    skimmer_associate[non_deut_count][2] = eta;
	    non_deut_count++;
	}		    
    }        //   end of track loop


    deut_per_event->Fill(non_deut_count);

    multiplicity_hybrid->Fill(hybrid_count);
    multiplicity_global->Fill(global_count);

         if(hybrid_count <  16) i_mult = 0;
    else if(hybrid_count <  23) i_mult = 1;
    else if(hybrid_count <  31) i_mult = 2;
    else if(hybrid_count <  42) i_mult = 3;
    else if(hybrid_count < 149) i_mult = 4;

    if(non_deut_count > 0  &&  deut_candidate_event == 0)
    {
	*associates<<i_mult<<" "<<i_zvert<<" "<<non_deut_count<<endl;
	for(int i=0; i<non_deut_count; i++)
	{
	    *associates<<skimmer_associate[i][0]<<" "<<skimmer_associate[i][1]<<" "<<skimmer_associate[i][2]<<endl;
	}
    }


    if(trig_track_count > 0  &&  deut_candidate_event == 0)
    {
	*triggers<<i_mult<<" "<<i_zvert<<" "<<trig_track_count<<endl;
	for(int i=0; i<trig_track_count; i++)
	{
	    *triggers<<skimmer_trigger[i][0]<<" "<<skimmer_trigger[i][1]<<" "<<skimmer_trigger[i][2]<<endl;
	}
    }



    


	 
    if(do_aod_skim == 1)
    {
	AliAODHandler *oh = (AliAODHandler*)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
	if (oh)    oh->SetFillAOD(kFALSE);

//  if ((non_deut_count>=1)  &&  (trig_03_track_count>0)  &&  oh)
//  { // markA
//    if ((non_deut_count >= 1) && oh  &&  i_mult == multiplicity_cut)
	if ((non_deut_count >= 1) && oh)
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
	
	    if (1) {   AliAODHeader    *out = (AliAODHeader*)eout->GetHeader();                      AliAODHeader    *in = (AliAODHeader*)evin->GetHeader();  	        *out = *in;                  }
	    if (1) {   AliTOFHeader    *out = const_cast<AliTOFHeader*>(eout->GetTOFHeader()); const AliTOFHeader    *in =                evin->GetTOFHeader();	        *out = *in;                  }
	    if (1) {   AliAODVZERO     *out =                eout->GetVZEROData();	             AliAODVZERO     *in =                evin->GetVZEROData();	        *out = *in;                  }
	    if (1) {   AliAODTZERO     *out =                eout->GetTZEROData();	             AliAODTZERO     *in =                evin->GetTZEROData();         *out = *in;                  }
	    if (1) {   TClonesArray    *out =                eout->GetTracks();	                     TClonesArray    *in =                evin->GetTracks();	    new (out) TClonesArray(*in);     }
	    if (1) {   TClonesArray    *out =                eout->GetVertices();	             TClonesArray    *in =                evin->GetVertices();      new (out) TClonesArray(*in);     }
	    if (1) {   TClonesArray    *out =                eout->GetCaloClusters();	             TClonesArray    *in =                evin->GetCaloClusters();  new (out) TClonesArray(*in);     }
	    if (1) {   AliAODCaloCells *out =                eout->GetEMCALCells();                  AliAODCaloCells *in =                evin->GetEMCALCells();    new (out) AliAODCaloCells(*in);  }
	}
    }
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}




//_____________________________________________________________________________
void AliAnalysisTaskCorPIDTOFQA::Terminate(Option_t *)
{
    associates->close();
    triggers->close();
//    file.close();
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
