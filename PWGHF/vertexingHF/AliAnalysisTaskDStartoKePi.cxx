#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include "TVector3.h"

//------------------------
#include "TH1D.h"  //!
#include "TH2D.h"  //!
#include "TProfile.h" //!
#include "TProfile2D.h" //!
#include "TCanvas.h" //!
#include "TF1.h" //!
#include "TLorentzVector.h"
//------------------------

//------------------------
#include "AliHelix.h"
#include "TSystem.h"
//------------------------

#include "AliTracker.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliKalmanTrack.h"

#include "AliTRDpadPlane.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"

#include "AliPIDResponse.h"
#include "AliPID.h"

#include "AliESDtrackCuts.h"

#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliESDRun.h"

#include "AliMultSelection.h"

//------------------------
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
//------------------------

#include "AliAnalysisTaskDStartoKePi.h"

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;


static const Int_t    N_D0_cuts               = 20;
static const Int_t    N_pT_bins               = 10;
static const Int_t    N_centrality_bins       = 10;
static const Int_t    N_ME_z_vertex_positions = 5;

static Int_t flag_plot_event = 0;
static TString HistName;

static const Double_t mass_array[16]         = {0.0,0.0,0.00051099892,0.00051099892,0.0,0.105658369,0.105658369,0.1349766,0.13957018,0.13957018,0.497648,0.493677,0.493677,0.93956536,0.93827203,0.93827203}; // in GeV/c^2  {empty,gamma,e+,e-,neutrino,mu+,mu-,pi0,pi+,pi-,K0long,K+,K-,n,p,anti-p}
static const Double_t pT_ranges[N_pT_bins+1] = {0.0,0.5,1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.0,10.0};
static const Double_t z_vertex_ranges[N_ME_z_vertex_positions+1] = {-20.0,-8.0,-3.0,3.0,8.0,20.0};
static const Double_t mass_D0 = 1.86484;
static const Double_t c_cm_per_ps = 0.0299792; // speed of light in cm per pico second
static const Int_t    N_pid_cuts = 5;
static const Int_t    N_cuts = 4;
static const Int_t    ME_event_buffer_size    = 4;
static TFile* dfile;

ClassImp(AliDStarEvent)
    ClassImp(AliDStarTrack)


    ClassImp(AliAnalysisTaskDStartoKePi)

    //________________________________________________________________________
    AliAnalysisTaskDStartoKePi::AliAnalysisTaskDStartoKePi(const char *name)
    : AliAnalysisTaskSE(name),
    fESD(0x0),
    fAOD(0x0),
    fDigitsInputFileName("TRD.FltDigits.root"), fDigitsInputFile(0),
    fDigitsOutputFileName(""), fDigitsOutputFile(0),
    fDigMan(0),fGeo(0),
    AS_Event(0),AS_Track(0),
    Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
    h_delta_dca(0x0),h_invariant_mass_DStar(0x0),h_invariant_mass_DStar_nom(0x0),h_invariant_mass_D0(0x0),h_delta_invariant_mass(0x0),h_multiplicities(0x0),
    TP_pt_bins(0x0),TP_dca_vs_p(0x0),h2D_mass2_vs_p(0x0),TP_counts_pid(0x0),h_pT(0x0),TP_TPC_chi2_vs_pT(0x0),h_runnumber(0x0),
    fListOfHistos(0x0),fTree(0x0), fPIDResponse(0), EsdTrackCuts(0)
{
    // Constructor

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());

    // Output slot #0 id reserved by the base class for AOD
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());

}

//_______________________________________________________________________
TFile* AliAnalysisTaskDStartoKePi::OpenDigitsFile(TString inputfile,
				   TString digfile,
				   TString opt)
{
    // we should check if we are reading ESDs or AODs - for now, only
    // ESDs are supported

    cout << "" << endl;
    cout << "In OpenDigitsFile" << endl;
    //cout << "Digits file name: " << digfile.Data() << endl;

    if(digfile == "")
    {
	cout << "WARNING: No TRD digits file available" << endl;
	return NULL;
    }

    // TGrid::Connect("alien")
    // construct the name of the digits file from the input file
    inputfile.ReplaceAll("AliESDs.root", digfile);
    //TString inputfile_LF = "alien:///alice/data/2016/LHC16q/000265525/pass1_CENT_wSDD/16000265525037.6203/TRD.FltDigits.root";
    TString inputfile_LF = inputfile;

    // open the file
    AliInfo( "opening digits file " + inputfile_LF + " with option \"" + opt + "\"");

    cout << "inputfile: " << inputfile_LF.Data() << endl;
    //TFile* dfile = new TFile(inputfile_LF, opt);
    //if(dfile) delete dfile;
    dfile = TFile::Open(inputfile_LF);
    cout << "After TRD digits file" << endl;
    cout << "" << endl;

    if(!dfile)
    {
	AliWarning("digits file '" + inputfile + "' cannot be opened");
    }

    return dfile;
}


//_______________________________________________________________________
Bool_t AliAnalysisTaskDStartoKePi::UserNotify()
{
    cout << "" << endl;
    cout << "In UserNotify" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;


    delete fDigitsInputFile;
    delete fDigitsOutputFile;

    cout << "Digits file pointers deleted" << endl;

    cout << "All pointers deleted" << endl;

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
	(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if ( ! esdH ) return kFALSE;
    if ( ! esdH->GetTree() ) return kFALSE;
    if ( ! esdH->GetTree()->GetCurrentFile() ) return kFALSE;

    TString fname = esdH->GetTree()->GetCurrentFile()->GetName();
    TString Tree_name = esdH->GetTree()->GetName();
    cout << "fname: " << fname << ", tree name: " << Tree_name <<  endl;
    FileStat_t file_stat;
    Int_t PathInfo = gSystem->GetPathInfo(fname.Data(),file_stat);
    cout << "PathInfo: " << PathInfo << endl;
    TFile* file = TFile::Open(fname.Data());
    cout << "Zombie: " << file->IsZombie() << ", header size: " << file->Sizeof() << ", FileBytesRead: " << file->GetFileBytesRead() << endl;

    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!inputHandler)
    {
	printf("WARNING: Inputhandler not available \n");
    }
    else
    {
	printf("Inputhandler available \n");

	fPIDResponse = inputHandler->GetPIDResponse();

        cout << "Got PID response" << endl;
    }

    fEventNoInFile = -1;
    N_good_events  = 0;


    fDigitsInputFile  = OpenDigitsFile(fname,fDigitsInputFileName,""); // <-
    EsdTrackCuts = new AliESDtrackCuts();

    EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
    EsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.52);
    EsdTrackCuts->AliESDtrackCuts::SetMinNClustersTPC(70); // 70, Automatically requires TPC refitted tracks?
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXY(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetPtRange(0.15,200.0); // 0.15, 200.0
    EsdTrackCuts->AliESDtrackCuts::SetEtaRange(-0.85,0.85); // 0.85
    EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE); // new
    EsdTrackCuts->AliESDtrackCuts::SetMaxFractionSharedTPCClusters(0.2); // new

    // create the digits manager
    cout << "" << endl;
    cout << "________________________________________________________________________" << endl;
    cout << "Created AliTRDdigitsManager" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;
    cout << "" << endl;

    // create a TRD geometry, needed for matching digits to tracks
    fGeo = new AliTRDgeometry;
    if(!fGeo)
    {
	AliFatal("cannot create geometry ");
    }

    //if(fDigMan) delete fDigMan;
    fDigMan = new AliTRDdigitsManager;
    fDigMan->CreateArrays();

    for(Int_t i_det = 0; i_det < 5; i_det++)
    {
	Int_t N_columns   = fDigMan->GetDigits(i_det)->GetNcol();
	cout << "i_det: " << i_det << ", N_columns: " << N_columns << endl;
    }


    cout << "End of UserNotify" << endl;
    return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskDStartoKePi::UserCreateOutputObjects()
{
    cout << "" << endl;
    cout << "In UserCreateOutputObjects" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;



    //-------------------------------------------------
    vec_track_info.resize(N_ME_z_vertex_positions);
    vec_two_prong_info.resize(N_ME_z_vertex_positions);
    vec_buffer_counter.resize(N_ME_z_vertex_positions);
    for(Int_t i_z_vertex = 0; i_z_vertex < N_ME_z_vertex_positions; i_z_vertex++)
    {
	vec_track_info[i_z_vertex].resize(N_centrality_bins);
        vec_two_prong_info[i_z_vertex].resize(N_centrality_bins);
	vec_buffer_counter[i_z_vertex].resize(N_centrality_bins);
	for(Int_t i_mult = 0; i_mult < N_centrality_bins; i_mult++)
	{
	    vec_buffer_counter[i_z_vertex][i_mult] = 0;
	    vec_track_info[i_z_vertex][i_mult].vec_data.resize(ME_event_buffer_size);
	    vec_two_prong_info[i_z_vertex][i_mult].vec_data.resize(ME_event_buffer_size);
	    for(Int_t i_buffer = 0; i_buffer < ME_event_buffer_size; i_buffer++)
	    {
		vec_track_info[i_z_vertex][i_mult].vec_data[i_buffer].resize(8); // pid: e-, pi-, K-, p-, e+, pi+, K+, p+
		vec_two_prong_info[i_z_vertex][i_mult].vec_data[i_buffer].resize(12); // reconstruction channels
		for(Int_t i_pid = 0; i_pid < 8; i_pid++)
		{
		    vec_track_info[i_z_vertex][i_mult].vec_data[i_buffer][i_pid].resize(17); // nSigma_dEdx, nSigma_TOF_e, px, py, pz, ITS refit, dca_xy, track_number + 6 for dir and base vectors + path dca, dca_z, ITS_cls
		}
		for(Int_t i_rec = 0; i_rec < 12; i_rec++)
		{
		    vec_two_prong_info[i_z_vertex][i_mult].vec_data[i_buffer][i_rec].resize(N_cuts*30); // cuts*(px, py, pz, M, track_id_A, track_id_B)
		}
	    }
	}
    }
    //-------------------------------------------------


    OpenFile(1);
    cout << "File opened" << endl;

    fListOfHistos = new TList();
    fListOfHistos ->SetOwner();

    h_delta_dca = new TH1D("h_delta_dca","h_delta_dca",1000,-0.05,1.0);
    h_delta_dca ->GetXaxis()->SetTitle("#Deltadca");
    fListOfHistos->Add(h_delta_dca);


    h_invariant_mass_DStar.resize(12); // K- + pi+, K+ + pi-, K- + e+, K+ + e-, same for ME, same for LS
    for(Int_t i_channel = 0 ; i_channel < 12; i_channel++)
    {
	h_invariant_mass_DStar[i_channel].resize(N_centrality_bins);
	for(Int_t i_centrality = 0; i_centrality < N_centrality_bins; i_centrality++)
	{
	    h_invariant_mass_DStar[i_channel][i_centrality].resize(N_pT_bins);
	    for(Int_t i_pT = 0; i_pT < h_invariant_mass_DStar[i_channel][i_centrality].size(); i_pT++)
	    {
		h_invariant_mass_DStar[i_channel][i_centrality][i_pT].resize(N_cuts);
		for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
		{
		    HistName = "h_invariant_mass_DStar_";
		    HistName += i_channel;
		    HistName += "_centr_";
		    HistName += i_centrality;
		    HistName += "_pT_";
		    HistName += i_pT;
		    HistName += "_cut_";
		    HistName += i_cut;
		    h_invariant_mass_DStar[i_channel][i_centrality][i_pT][i_cut] = new TH1D(HistName.Data(),HistName.Data(),1400,1.950,3.35);
		    fListOfHistos->Add(h_invariant_mass_DStar[i_channel][i_centrality][i_pT][i_cut]);
		}
	    }
	}
    }


    h_invariant_mass_DStar_nom.resize(12); // K- + pi+, K+ + pi-, K- + e+, K+ + e-, same for ME, same for LS
    for(Int_t i_channel = 0 ; i_channel < 12; i_channel++)
    {
	h_invariant_mass_DStar_nom[i_channel].resize(N_centrality_bins);
	for(Int_t i_centrality = 0; i_centrality < N_centrality_bins; i_centrality++)
	{
	    h_invariant_mass_DStar_nom[i_channel][i_centrality].resize(N_pT_bins);
	    for(Int_t i_pT = 0; i_pT < h_invariant_mass_DStar_nom[i_channel][i_centrality].size(); i_pT++)
	    {
		h_invariant_mass_DStar_nom[i_channel][i_centrality][i_pT].resize(N_cuts);
		for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
		{
		    HistName = "h_invariant_mass_DStar_nom_";
		    HistName += i_channel;
		    HistName += "_centr_";
		    HistName += i_centrality;
		    HistName += "_pT_";
		    HistName += i_pT;
		    HistName += "_cut_";
		    HistName += i_cut;
		    h_invariant_mass_DStar_nom[i_channel][i_centrality][i_pT][i_cut] = new TH1D(HistName.Data(),HistName.Data(),1400,1.950,3.35);
		    fListOfHistos->Add(h_invariant_mass_DStar_nom[i_channel][i_centrality][i_pT][i_cut]);
		}
	    }
	}
    }

    h_invariant_mass_D0.resize(12); // K- + pi+, K+ + pi-, K- + e+, K+ + e-, same for ME, same for LS
    for(Int_t i_channel = 0 ; i_channel < 12; i_channel++)
    {
	h_invariant_mass_D0[i_channel].resize(N_centrality_bins);
	for(Int_t i_centrality = 0; i_centrality < N_centrality_bins; i_centrality++)
	{
	    h_invariant_mass_D0[i_channel][i_centrality].resize(N_pT_bins);
	    for(Int_t i_pT = 0; i_pT < h_invariant_mass_D0[i_channel][i_centrality].size(); i_pT++)
	    {
		h_invariant_mass_D0[i_channel][i_centrality][i_pT].resize(N_cuts);
		for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
		{
		    HistName = "h_invariant_mass_D0_";
		    HistName += i_channel;
		    HistName += "_centr_";
		    HistName += i_centrality;
		    HistName += "_pT_";
		    HistName += i_pT;
		    HistName += "_cut_";
		    HistName += i_cut;
		    h_invariant_mass_D0[i_channel][i_centrality][i_pT][i_cut] = new TH1D(HistName.Data(),HistName.Data(),2400,0.450,2.850);
		    fListOfHistos->Add(h_invariant_mass_D0[i_channel][i_centrality][i_pT][i_cut]);
		}
	    }
	}
    }

    h_delta_invariant_mass.resize(12); // K- + pi+, K+ + pi-, K- + e+, K+ + e-, same for ME, same for LS
    for(Int_t i_channel = 0 ; i_channel < 12; i_channel++)
    {
	h_delta_invariant_mass[i_channel].resize(N_centrality_bins);
	for(Int_t i_centrality = 0; i_centrality < N_centrality_bins; i_centrality++)
	{
	    h_delta_invariant_mass[i_channel][i_centrality].resize(N_pT_bins);
	    for(Int_t i_pT = 0; i_pT < h_delta_invariant_mass[i_channel][i_centrality].size(); i_pT++)
	    {
		h_delta_invariant_mass[i_channel][i_centrality][i_pT].resize(N_cuts);
		for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
		{
		    HistName = "h_delta_invariant_mass_";
		    HistName += i_channel;
		    HistName += "_centr_";
		    HistName += i_centrality;
		    HistName += "_pT_";
		    HistName += i_pT;
		    HistName += "_cut_";
		    HistName += i_cut;
		    h_delta_invariant_mass[i_channel][i_centrality][i_pT][i_cut] = new TH1D(HistName.Data(),HistName.Data(),6000,0.0,1.5);
		    fListOfHistos->Add(h_delta_invariant_mass[i_channel][i_centrality][i_pT][i_cut]);
		}
	    }
	}
    }

    h_multiplicities.resize(11);
    for(Int_t i_mult = 0; i_mult < 11; i_mult++)
    {
	HistName = "h_multiplicities_";
	HistName += i_mult;
	h_multiplicities[i_mult] = new TH1D(HistName.Data(),HistName.Data(),100,0,100);
	fListOfHistos->Add(h_multiplicities[i_mult]);
    }

    TP_pt_bins = new TProfile("TP_pt_bins","TP_pt_bins",100,0,100);
    fListOfHistos->Add(TP_pt_bins);

    TP_dca_vs_p = new TProfile("TP_dca_vs_p","TP_dca_vs_p",100,0.0,5.0);
    fListOfHistos->Add(TP_dca_vs_p);

    h2D_mass2_vs_p.resize(8); // e-, pi-, K-, p-, e+, pi+, K+, p+
    for(Int_t i_pid = 0; i_pid < 8; i_pid++)
    {
	h2D_mass2_vs_p[i_pid].resize(N_pid_cuts);
	for(Int_t i_cut = 0; i_cut < N_pid_cuts; i_cut++)
	{
	    HistName = "h2D_mass2_vs_p_pid_";
	    HistName += i_pid;
	    HistName += "_cut_";
            HistName += i_cut;
	    h2D_mass2_vs_p[i_pid][i_cut] = new TH2D(HistName.Data(),HistName.Data(),40,0,10,2000,-0.5,3.5);
            fListOfHistos->Add(h2D_mass2_vs_p[i_pid][i_cut]);
	}
    }

    TP_counts_pid.resize(N_pT_bins);
    for(Int_t i_pT = 0; i_pT < N_pT_bins; i_pT++)
    {
	HistName = "TP_counts_pid_pT";
        HistName += i_pT;
	TP_counts_pid[i_pT] = new TProfile(HistName.Data(),HistName.Data(),8,0,8);
        fListOfHistos->Add(TP_counts_pid[i_pT]);
    }

    h_pT.resize(2); // -,+
    for(Int_t i_charge = 0; i_charge < 2; i_charge++)
    {
        h_pT[i_charge].resize(4);
	for(Int_t i_pid = 0; i_pid < 4; i_pid++)
	{
	    HistName = "h_pT_charge";
	    HistName += i_charge;
	    HistName += "_pid_";
	    HistName += i_pid;
	    h_pT[i_charge][i_pid] = new TH1D(HistName.Data(),HistName.Data(),500,0,20);
	    fListOfHistos->Add(h_pT[i_charge][i_pid]);
	}
    }

    TP_TPC_chi2_vs_pT.resize(2); // -,+
    for(Int_t i_charge = 0; i_charge < 2; i_charge++)
    {
	HistName = "TP_TPC_chi2_vs_pT_charge";
	HistName += i_charge;
	TP_TPC_chi2_vs_pT[i_charge] = new TProfile(HistName.Data(),HistName.Data(),500,0,20);
        fListOfHistos->Add(TP_TPC_chi2_vs_pT[i_charge]);
    }

    h_runnumber = new TH1D("h_runnumber","h_runnumber",2000,264000,266000);
    fListOfHistos ->Add(h_runnumber);


    OpenFile(2);

    AS_Event       = new AliDStarEvent();
    AS_Track       = new AliDStarTrack();
    Tree_AS_Event  = NULL;
    Tree_AS_Event  = new TTree("Tree_AS_Event" , "AS_Events" );
    Tree_AS_Event  ->Branch("Tree_AS_Event_branch"  , "AS_Event", AS_Event );

    PostData(1,fListOfHistos);
    PostData(2,Tree_AS_Event);


    for(Int_t i_pT = 0; i_pT < N_pT_bins+1; i_pT++)
    {
	TP_pt_bins ->SetBinContent(i_pT+1,pT_ranges[i_pT]);
    }

    cout << "PostData called" << endl;

}



//________________________________________________________________________
Bool_t AliAnalysisTaskDStartoKePi::NextEvent(Bool_t preload)
{
    fEventNoInFile++;

    if (fEventNoInFile != 0  &&  fEventNoInFile % 5 == 0)
	cout << "." << flush;
    if (fEventNoInFile != 0  &&  fEventNoInFile % 50 == 0)
    {
	cout << "Event: " << fEventNoInFile << " ==> Processing data (DStar Analysis) " << endl;
    }

    return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskDStartoKePi::UserExec(Option_t *)
{
    //cout << "" << endl;
    //cout << "Analysis started" << endl;
    //cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;


    //-----------------------------------------------------------------
    // IMPORTANT: call NextEvent() for book-keeping
    NextEvent();
    //if(fEventNoInFile > 4) return;
    //-----------------------------------------------------------------


    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    // prepare event data structures
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());

    if(fIsAOD) // aod
    {
	fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
	if(!fAOD)
	{
	    printf("ERROR: fAOD not available\n");
	    return;
	}
    }
    else // esd
    {
	fESD = dynamic_cast<AliESDEvent*>(InputEvent());
	if(!fESD)
	{
	    printf("ERROR: fESD not available\n");
	    return;
	}
    }


    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man)
    {
	//cout << "Got AliAnalysisManager" << endl;
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
	if(inputHandler)
	{
	    //cout << "Got AliInputEventHandler" << endl;
	    fPIDResponse = inputHandler->GetPIDResponse();
	}
    }
    //cout << "cent: " << fPIDResponse->GetCurrentCentrality() << endl;


    Int_t    N_tracks        = -1;
    Int_t    N_TRD_tracks    = -1;
    Int_t    N_TRD_tracklets = -1;
    Int_t    RunNum          = -1;
    Float_t  magF            = -100.0;
    Double_t T0zVertex       = -100.0;
    Double_t MeanBeamIntAA   = -100.0;
    AliCentrality* Centrality;
    TVector3 PrimVertex;
    if(fIsAOD) // AOD
    {
	N_tracks         = fAOD ->GetNumberOfTracks();
	N_TRD_tracks     = fAOD ->GetNumberOfTrdTracks();
	RunNum           = fAOD ->GetRunNumber();
	AliAODVertex*  PrimVertex_AOD   = fAOD ->GetPrimaryVertex();
	magF             = fAOD ->GetMagneticField();
	Centrality       = fAOD ->GetCentrality();

	AS_Event ->setTriggerWord(fAOD->GetFiredTriggerClasses());
        PrimVertex.SetXYZ(PrimVertex_AOD->GetX(),PrimVertex_AOD->GetY(),PrimVertex_AOD->GetZ());
    }
    if(!fIsAOD) // ESD
    {
	N_tracks         = fESD ->GetNumberOfTracks();
	N_TRD_tracks     = fESD ->GetNumberOfTrdTracks();
	N_TRD_tracklets  = fESD ->GetNumberOfTrdTracklets();
	RunNum           = fESD ->GetRunNumber();
	const AliESDVertex*  PrimVertex_ESD   = fESD ->GetPrimaryVertex();
	magF             = fESD ->GetMagneticField();
	T0zVertex        = fESD ->GetT0zVertex();
	Centrality       = fESD ->GetCentrality();
	MeanBeamIntAA    = fESD ->GetESDRun()->GetMeanIntensity(0,0);

	AS_Event ->setTriggerWord(fESD->GetFiredTriggerClasses());
        PrimVertex.SetXYZ(PrimVertex_ESD->GetX(),PrimVertex_ESD->GetY(),PrimVertex_ESD->GetZ());
    }

    TV3_prim_vertex.SetXYZ(PrimVertex.X(),PrimVertex.Y(),PrimVertex.Z());

    //AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");


    // Fill event information
    AS_Event ->clearTrackList();
    AS_Event ->setx(PrimVertex.X());
    AS_Event ->sety(PrimVertex.Y());
    AS_Event ->setz(PrimVertex.Z());
    AS_Event ->setid(RunNum);
    AS_Event ->setN_tracks(N_tracks);
    AS_Event ->setN_TRD_tracklets(N_TRD_tracklets);
    AS_Event ->setBeamIntAA(MeanBeamIntAA);
    AS_Event ->setT0zVertex(T0zVertex);

    if(MultSelection)
    {
	// V0MEq, V0AEq, V0CEq, SPDTracklets

	AS_Event ->setcent_class_ZNA(MultSelection->GetMultiplicityPercentile("ZNA"));
	AS_Event ->setcent_class_ZNC(MultSelection->GetMultiplicityPercentile("ZNC"));
	AS_Event ->setcent_class_V0A(MultSelection->GetMultiplicityPercentile("V0A"));
	AS_Event ->setcent_class_V0C(MultSelection->GetMultiplicityPercentile("V0C"));
	AS_Event ->setcent_class_V0M(MultSelection->GetMultiplicityPercentile("V0M"));
	AS_Event ->setcent_class_CL0(MultSelection->GetMultiplicityPercentile("CL0"));
	AS_Event ->setcent_class_CL1(MultSelection->GetMultiplicityPercentile("CL1"));
	AS_Event ->setcent_class_SPD(MultSelection->GetMultiplicityPercentile("SPDTracklets"));
	AS_Event ->setcent_class_V0MEq(MultSelection->GetMultiplicityPercentile("V0MEq"));
	AS_Event ->setcent_class_V0AEq(MultSelection->GetMultiplicityPercentile("V0AEq"));
	AS_Event ->setcent_class_V0CEq(MultSelection->GetMultiplicityPercentile("V0CEq"));
    }


    //-----------------------------------------------------------------
    //cout << "" << endl;
    //cout << "" << endl;
    //cout << "----------------------------------------------------------------------------------------" << endl;
    //cout << "Event number: " << fEventNoInFile << ", event number with TRD digits: " << N_good_events << ", SPC centr.: " << MultSelection->GetMultiplicityPercentile("SPDTracklets") << endl;
    //-----------------------------------------------------------------


    //Int_t centrality_bin = (Int_t)(MultSelection->GetMultiplicityPercentile("SPDTracklets")/((Double_t)N_centrality_bins));
    Int_t centrality_bin = (Int_t)(MultSelection->GetMultiplicityPercentile("V0M")/((Double_t)N_centrality_bins));

    // Determine z vertex bin
    Int_t z_bin = -1;
    for(Int_t i_z_vertex = 0; i_z_vertex < N_ME_z_vertex_positions; i_z_vertex++)
    {
	if(PrimVertex.Z() >= z_vertex_ranges[i_z_vertex] && PrimVertex.Z() < z_vertex_ranges[i_z_vertex+1])
	{
	    z_bin = i_z_vertex;
	    break;
	}
    }
    if(z_bin >= 0 && centrality_bin >= 0)
    {
	h_multiplicities[0]  ->Fill(MultSelection->GetMultiplicityPercentile("ZNA"));
	h_multiplicities[1]  ->Fill(MultSelection->GetMultiplicityPercentile("ZNC"));
	h_multiplicities[2]  ->Fill(MultSelection->GetMultiplicityPercentile("V0A"));
	h_multiplicities[3]  ->Fill(MultSelection->GetMultiplicityPercentile("V0C"));
	h_multiplicities[4]  ->Fill(MultSelection->GetMultiplicityPercentile("V0M"));
	h_multiplicities[5]  ->Fill(MultSelection->GetMultiplicityPercentile("CL0"));
	h_multiplicities[6]  ->Fill(MultSelection->GetMultiplicityPercentile("CL1"));
	h_multiplicities[7]  ->Fill(MultSelection->GetMultiplicityPercentile("SPDTracklets"));
	h_multiplicities[8]  ->Fill(MultSelection->GetMultiplicityPercentile("V0MEq"));
	h_multiplicities[9]  ->Fill(MultSelection->GetMultiplicityPercentile("V0AEq"));
	h_multiplicities[10] ->Fill(MultSelection->GetMultiplicityPercentile("V0CEq"));



	//-----------------------------------------------------------------
	if(N_tracks == 0)
	{
	    // Skip empty event
	    return;
	}

        Int_t buffer_bin = vec_buffer_counter[z_bin][centrality_bin]%ME_event_buffer_size; // actual event
        vec_buffer_counter[z_bin][centrality_bin]++;
	Int_t flag_do_ME = 0;
        Int_t buffer_bin_ME = -1;
	if(vec_buffer_counter[z_bin][centrality_bin] >= ME_event_buffer_size)
	{
	    flag_do_ME = 1; // ME ring buffer for this z and mult index is full
            buffer_bin_ME = vec_buffer_counter[z_bin][centrality_bin]%ME_event_buffer_size; // reference event for mixing (note that the vec_buffer_counter did a ++ in comparison to buffer_bin);
	}

	//printf("There are %d tracks in this event\n", N_tracks);
	//printf("There are %d TRD tracks in this event\n", N_TRD_tracks);
	//printf("There are %d TRD tracklets in this event\n", N_TRD_tracklets);

	Int_t N_good_tracks = 0;

	//-----------------------------------------------------------------
	// Track loop
	//cout << "" << endl;
	//cout << "-----------------------------------------------------------------" << endl;
	//cout << "Start matching " << N_tracks << " TPC tracks with" << TV3_TRD_hits_middle.size() << " TRD pads" << endl;
	N_good_tracks = 0;
	Int_t N_matched_TRD_hits_total = 0;

	// Store track information
        // centrality, z-vertex position, buffer
	// 8 PIDs, e-, pi-, K-, p-, e+, pi+, K+, p+ -> basic PID cuts, 3 sigma on TPC dE/dx + 3 sigma on TOF (only if available)
	// nSigma_dEdx, nSigma_TOF_e, px, py, pz, ITS refit, dca, track_number


	Int_t N_pids[8]   = {0,0,0,0,0,0,0,0};
	Int_t N_charge[2] = {0,0};
	//cout << "" << endl;
	//printf("New event: %d \n",fEventNoInFile);
        h_runnumber ->Fill(RunNum);
	for(Int_t iTracks = 0; iTracks < N_tracks; iTracks++)
	{
	    //---------------------------------------------------------------
	    // Gather track information

	    // We always want the track
	    AliESDtrack* track_ESD;
            AliAODTrack* track_AOD;

	    Double_t TRD_signal;
	    Double_t Track_pT;
	    Double_t Track_p;
	    Double_t p_vec[3];
	    Int_t    charge;
	    Double_t Track_phi;
	    Double_t Track_theta;
	    Double_t Track_eta;
	    Double_t TPC_chi2;
	    Double_t TPC_signal;
	    Double_t TOF_signal;
	    Double_t Track_length;
	    UShort_t N_TPC_cls;
	    Double_t TPC_crossed_rows;
            ULong_t  status;

	    // e = 0, pion = 1, kaon = 2, proton = 3
	    Double_t Track_PID[8] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};

	    Float_t track_xy_impact,track_z_impact;

	    Int_t N_ITS_cls = 0;

	    if(fIsAOD) // AOD
	    {
                track_AOD = (AliAODTrack*)fAOD->GetTrack(iTracks);
		if(!track_AOD)
		{
		    printf("ERROR: Could not receive track %d\n", iTracks);
		    continue;
		}

		track_AOD->GetPxPyPz(p_vec);
		TRD_signal       = track_AOD ->GetTRDsignal(); // truncated mean signal?
		Track_pT         = track_AOD ->Pt();
		Track_p          = track_AOD ->P();
		charge           = track_AOD ->Charge();
		Track_phi        = track_AOD ->Phi();
		Track_theta      = track_AOD ->Theta();
		Track_eta        = track_AOD ->Eta();
		TPC_chi2         = track_AOD ->GetTPCchi2();
		TPC_signal       = track_AOD ->GetTPCsignal(); // dE/dx?
		TOF_signal       = track_AOD ->GetTOFsignal(); // time-of-flight?
		Track_length     = track_AOD ->GetIntegratedLength();
		N_TPC_cls        = track_AOD ->GetTPCNcls();
		TPC_crossed_rows = track_AOD ->GetTPCCrossedRows();
		status           = track_AOD ->GetStatus();

		track_AOD ->GetImpactParameters(track_xy_impact,track_z_impact);

		for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
		{
		    if(track_AOD ->HasPointOnITSLayer(i_ITS_layer))
		    {
			N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
		    }
		}

		// nSigma TPC
		Track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track_AOD,AliPID::kElectron);
		Track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track_AOD,AliPID::kPion);
		Track_PID[2] = fPIDResponse->NumberOfSigmasTPC(track_AOD,AliPID::kKaon);
		Track_PID[3] = fPIDResponse->NumberOfSigmasTPC(track_AOD,AliPID::kProton);

		// nSigma TOF, -999 in case there is no TOF hit
		Track_PID[4] = fPIDResponse->NumberOfSigmasTOF(track_AOD,AliPID::kElectron);
		Track_PID[5] = fPIDResponse->NumberOfSigmasTOF(track_AOD,AliPID::kPion);
		Track_PID[6] = fPIDResponse->NumberOfSigmasTOF(track_AOD,AliPID::kKaon);
		Track_PID[7] = fPIDResponse->NumberOfSigmasTOF(track_AOD,AliPID::kProton);
	    }
	    else // ESD
	    {
		track_ESD = (AliESDtrack*)fESD->GetTrack(iTracks);
		if(!track_ESD)
		{
		    printf("ERROR: Could not receive track %d\n", iTracks);
		    continue;
		}

		if(!EsdTrackCuts->AcceptTrack(track_ESD)) continue;
		track_ESD->GetPxPyPz(p_vec);

		TRD_signal       = track_ESD ->GetTRDsignal(); // truncated mean signal?
		Track_pT         = track_ESD ->Pt();
		Track_p          = track_ESD ->P();
		charge           = track_ESD ->Charge();
		Track_phi        = track_ESD ->Phi();
		Track_theta      = track_ESD ->Theta();
		Track_eta        = track_ESD ->Eta();
		TPC_chi2         = track_ESD ->GetTPCchi2();
		TPC_signal       = track_ESD ->GetTPCsignal(); // dE/dx?
		TOF_signal       = track_ESD ->GetTOFsignal(); // time-of-flight?
		Track_length     = track_ESD ->GetIntegratedLength();
		N_TPC_cls        = track_ESD ->GetTPCNcls();
		TPC_crossed_rows = track_ESD ->GetTPCCrossedRows();
		status           = track_ESD ->GetStatus();

		track_ESD ->GetImpactParameters(track_xy_impact,track_z_impact);

		for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
		{
		    if(track_ESD ->HasPointOnITSLayer(i_ITS_layer))
		    {
			N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
		    }
		}

		// nSigma TPC
		Track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track_ESD,AliPID::kElectron);
		Track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track_ESD,AliPID::kPion);
		Track_PID[2] = fPIDResponse->NumberOfSigmasTPC(track_ESD,AliPID::kKaon);
		Track_PID[3] = fPIDResponse->NumberOfSigmasTPC(track_ESD,AliPID::kProton);

		// nSigma TOF, -999 in case there is no TOF hit
		Track_PID[4] = fPIDResponse->NumberOfSigmasTOF(track_ESD,AliPID::kElectron);
		Track_PID[5] = fPIDResponse->NumberOfSigmasTOF(track_ESD,AliPID::kPion);
		Track_PID[6] = fPIDResponse->NumberOfSigmasTOF(track_ESD,AliPID::kKaon);
		Track_PID[7] = fPIDResponse->NumberOfSigmasTOF(track_ESD,AliPID::kProton);

		//Double_t x_param,p_param[5]; // p are the track paramters, 0 = Y, 1 = Z, 2 = Snp, 3 = Tgl (p[2]/pt), 4 = Signed1Pt
		//track_ESD->GetExternalParameters(x_param,p_param); // Those are the parameters used for vertexing
	    }



	    Double_t mass_sq = -999.0;
            if(TOF_signal > 0.0)
            {
                Double_t beta_sq  = TMath::Power((Track_length/TOF_signal)/c_cm_per_ps,2.0);
		if(beta_sq > 0.0)
		{
		    mass_sq = Track_p*Track_p*((1.0/beta_sq) - 1.0);
		}
            }

	    Int_t i_charge = 0;
	    if(charge > 0) i_charge = 1;
            N_charge[i_charge]++;
	    //cout << "charge: " << charge << ", i_charge: " << i_charge << endl;

	    Int_t ITS_refit = 0;
	    Int_t TPC_refit = 0;
	    Int_t track_status = 0;
	    if(((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit))
	    {
		ITS_refit = 1;
		track_status |= 1 << 0; // setting bit 0 to 1
	    }
	    if(((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit))
	    {
		TPC_refit = 1;
		track_status |= 1 << 1; // setting bit 1 to 1
	    }

	    //printf("track number: %d, TPC_refit: %d \n",iTracks,TPC_refit);


	    Double_t track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);

	    Double_t TRD_ADC_bin_width = 100.0;


	    //----------------------------------------------------------------------------
	    //Helix
	    AliAnalysisTaskDStartoKePi::FillHelix(track_ESD,magF);
	    //printf("pT: %f, curvature: %f \n",Track_pT,1.0/aliHelix.fHelix[4]);

	    Double_t path_min    = -9999.0;
	    Double_t dca_min_xy  = -9999.0;
	    Double_t dca_min_z   = -9999.0;
	    Double_t helix_points_dca_A[3] = {-9999.0,-9999.0,-9999.0};
	    Double_t helix_points_dca_B[3] = {-9999.0,-9999.0,-9999.0};
            Double_t dca_min = -100.0;
	    if(ITS_refit && track_total_impact < 0.5) // CHANGE
	    {
		// Find helix path closest to primary vertex
		Double_t pathA = 0.0;   // usually more outer
		FindDCAHelixPoint(TV3_prim_vertex,aliHelix,pathA-5.0,pathA+5.0,path_min,dca_min_xy,dca_min_z);
		aliHelix.Evaluate(path_min,helix_points_dca_A);
                aliHelix.Evaluate(path_min+0.01,helix_points_dca_B);
		//printf("dca: %f, dca_min: %f \n",track_total_impact,dca_min);
		Double_t diff_length = TMath::Sqrt(TMath::Power(helix_points_dca_A[0]-helix_points_dca_B[0],2) + TMath::Power(helix_points_dca_A[1]-helix_points_dca_B[1],2) + TMath::Power(helix_points_dca_A[2]-helix_points_dca_B[2],2));
		//printf("point A: {%f, %f, %f}, pointB: {%f, %f, %f}, |diff|: %f \n",helix_points_dca_A[0],helix_points_dca_A[1],helix_points_dca_A[2],helix_points_dca_B[0],helix_points_dca_B[1],helix_points_dca_B[2],diff_length);
	    }
	    dca_min = TMath::Sqrt(dca_min_xy*dca_min_xy + dca_min_z*dca_min_z);
	    //----------------------------------------------------------------------------



	    //-------------------
	    // The following is only to fill PID histograms for QA
	    const Double_t n_sigma_dEdx_cuts[N_pid_cuts] = {1000.0,3.5,3.0,2.5,2.0};
	    for(Int_t i_cut = 0; i_cut < N_pid_cuts; i_cut++)
	    {
		Int_t PID_array_dEdx[4] = {0,0,0,0}; // e, pi, K, p: multiple PIDs per track possible
		if(fabs(Track_PID[0]) < n_sigma_dEdx_cuts[i_cut] // Electron
		  )
		{
		    PID_array_dEdx[0] = 1;
		}
		if(fabs(Track_PID[1]) < n_sigma_dEdx_cuts[i_cut] // Pion
		  )
		{
		    PID_array_dEdx[1] = 1;
		}
		if(fabs(Track_PID[2]) < n_sigma_dEdx_cuts[i_cut] // Kaon
		  )
		{
		    PID_array_dEdx[2] = 1;
		}
		if(fabs(Track_PID[3]) < n_sigma_dEdx_cuts[i_cut] // Proton
		  )
		{
		    PID_array_dEdx[3] = 1;
		}
		for(Int_t i_PID = 0; i_PID < 4; i_PID++) // e, pi, K, p
		{
		    if(PID_array_dEdx[i_PID])
		    {
			Int_t PID_ME_D0 = i_PID + 4*i_charge;
			if(mass_sq > -500.0)
			{
			    //cout << "i_cut: " << i_cut << ", PID_ME_D0: " << PID_ME_D0 << ", p: " << Track_p << ", mass_sq: " << mass_sq << endl;
			    h2D_mass2_vs_p[PID_ME_D0][i_cut] ->Fill(Track_p,mass_sq); // e-, pi-, K-, p-, e+, pi+, K+, p+
			}
		    }
		}
	    }
	    //-------------------


	    //-------------------
	    // Do multiple PID for each track with wide PID cuts
            Int_t PID_array[4] = {0,0,0,0}; // e, pi, K, p: multiple PIDs per track possible
	    if(
	       (!ITS_refit && track_total_impact < 1.0) ||
	       (ITS_refit && track_total_impact < 0.1)
               //(!ITS_refit && dca_min < 1.0) ||
	       //(ITS_refit && dca_min < 0.1)
	      )
	    {
                //-------------------
		if(fabs(Track_PID[0]) < 2.5 // Electron
		   && (fabs(Track_PID[4]) < 2.5 || Track_PID[4] < -100.0) // Use TOF information only if available
		  )
		{
		    PID_array[0] = 1;
		}
		if(fabs(Track_PID[1]) < 2.5 // Pion
		   && (fabs(Track_PID[5]) < 2.5 || Track_PID[5] < -100.0) // Use TOF information only if available
		  )
		{
		    PID_array[1] = 1;
		}
		if(fabs(Track_PID[2]) < 2.5 // Kaon
		   && (fabs(Track_PID[6]) < 2.5 || Track_PID[6] < -100.0) // Use TOF information only if available
		  )
		{
		    PID_array[2] = 1;
		}
		if(fabs(Track_PID[3]) < 2.5 // Proton
		   && (fabs(Track_PID[7]) < 2.5 || Track_PID[7] < -100.0) // Use TOF information only if available
		  )
		{
		    PID_array[3] = 1;
		}
		//-------------------
	    }
	    //-------------------


	    for(Int_t i_pid = 0; i_pid < 4; i_pid++) // e-, pi-, K-, p-, e+, pi+, K+, p+
	    {
		if(!PID_array[i_pid]) continue;
		h_pT[i_charge][i_pid]   ->Fill(Track_pT);
	    }
	    TP_TPC_chi2_vs_pT[i_charge] ->Fill(Track_pT,TPC_chi2);


	    //-------------------
	    // Store the track information
	    for(Int_t i_pid = 0; i_pid < 4; i_pid++) // e-, pi-, K-, p-, e+, pi+, K+, p+
	    {
		if(!PID_array[i_pid]) continue;
                Int_t i_pid_use = i_pid;
		if(i_charge == 1) i_pid_use = i_pid + 4;

		// nSigma_dEdx, nSigma_TOF_e, px, py, pz, ITS refit, dca, track_number
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][0].push_back(Track_PID[i_pid_use%4]); // nSigma_dEdx
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][1].push_back(Track_PID[i_pid_use%4 + 4]); // nSigma_TOF
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][2].push_back(p_vec[0]); // px
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][3].push_back(p_vec[1]); // py
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][4].push_back(p_vec[2]); // pz
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][5].push_back(ITS_refit); // ITS refit
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][6].push_back(dca_min_xy); // dca to primary vertex in xy
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][7].push_back(iTracks); // track number
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][8].push_back(helix_points_dca_A[0]);
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][9].push_back(helix_points_dca_A[1]);
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][10].push_back(helix_points_dca_A[2]);
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][11].push_back(helix_points_dca_B[0] - helix_points_dca_A[0]);
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][12].push_back(helix_points_dca_B[1] - helix_points_dca_A[1]);
                vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][13].push_back(helix_points_dca_B[2] - helix_points_dca_A[2]);
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][14].push_back(path_min);
		vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][15].push_back(dca_min_z); // dca to primary vertex in z
                vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin][i_pid_use][16].push_back(N_ITS_cls);
	    }
	    //-------------------


	    N_good_tracks++;

	} // End of TPC track loop



	//-----------------------------------------------------------------
	// Do combinatorics with vertexing
	//-----------------------------------------------------------------



	//-----------------------------------------------------------------
        // Do combinatorics
        //      0    1   2   3   4   5    6   7
        // PID: e-, pi-, K-, p-, e+, pi+, K+, p+

        // reconstruction channels: K- + pi+, K+ + pi-, K- + e+, K+ + e-, same for ME, same for LS
        Int_t channel_array[2][12] = // Mother -> PID_A + PID_B
        {
            {2,6,2,6,2,6,2,6,2,6,2,6}, // PID_A
            {5,1,4,0,5,1,4,0,1,5,0,4}  // PID_B
	};
	Int_t pion_electron_channel[12]  = {0,0,1,1,0,0,1,1,0,0,1,1}; // 0 = pion, 1 = electron
	Int_t SE_ME_flag_array[12]       = {0,0,0,0,1,1,1,1,0,0,0,0}; // 0 = SE, 1 = ME
	Int_t Three_prong_channel_AB[12] = {0,1,2,3,0,1,2,3,0,1,2,3}; // Always take the like sign combinations from the two prong (first four channels)
	Int_t Three_prong_pid_C[12]      = {5,1,5,1,5,1,5,1,1,5,1,5}; // First 4 are SE, next 4 are ME and last 4 are LS, but they will be mixed with the first four of the two prong decay

	for(Int_t i_channel = 0; i_channel < 12; i_channel++) // K- + pi+, K+ + pi-, K- + e+, K+ + e-, same for ME, same for LS
	{
	    if(!flag_do_ME && SE_ME_flag_array[i_channel]) continue;

	    std::vector< std::vector< std::vector<Int_t> > > vec_track_id;
	    std::vector< std::vector<TLorentzVector> > vec_TLV = Rec_two_prong_primary_decay(vec_track_info,channel_array[0][i_channel],channel_array[1][i_channel],
									      z_bin,centrality_bin,buffer_bin,buffer_bin_ME,SE_ME_flag_array[i_channel],
									      N_cuts,vec_track_id);

	    //if(i_channel == 0 || i_channel == 1) printf("i_channel: %d, N(two prong): %d \n",i_channel,vec_TLV.size());

	    for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
	    {
		for(Int_t i_rec = 0; i_rec < vec_TLV[i_cut].size(); i_rec++)
		{
		    Double_t pT_rec         = vec_TLV[i_cut][i_rec].Pt();
		    Double_t invariant_mass = vec_TLV[i_cut][i_rec].M();

		    Int_t pT_bin = -1;
		    for(Int_t i_pT = 0; i_pT < N_pT_bins+1; i_pT++)
		    {
			if(pT_rec > pT_ranges[i_pT] && pT_rec <= pT_ranges[i_pT+1])
			{
			    pT_bin = i_pT;
			    break;
			}
		    }
		    if(pT_bin >= 0)
		    {
			h_invariant_mass_D0[i_channel][centrality_bin][pT_bin][i_cut] ->Fill(invariant_mass);
                        //if(i_channel == 0 && centrality_bin == 0 && pT_bin_D0 == 0) printf("mass_D0: %f \n",invariant_mass_D0);
		    }

		    //printf("M: %f, channel: %d \n",invariant_mass,pion_electron_channel[i_channel]);

		    if(!(((invariant_mass > 1.8 && invariant_mass < 1.9) && pion_electron_channel[i_channel] == 0 ) ||
			 ((invariant_mass > 0.5 && invariant_mass < 1.9) && pion_electron_channel[i_channel] == 1 ))) continue; // cut on D0 invariant mass

		    // px, py, pz, track_id_A, track_id_B
		    vec_two_prong_info[z_bin][centrality_bin].vec_data[buffer_bin][i_channel][0+30*i_cut].push_back(vec_TLV[i_cut][i_rec].Px()); // px
		    vec_two_prong_info[z_bin][centrality_bin].vec_data[buffer_bin][i_channel][1+30*i_cut].push_back(vec_TLV[i_cut][i_rec].Py()); // py
		    vec_two_prong_info[z_bin][centrality_bin].vec_data[buffer_bin][i_channel][2+30*i_cut].push_back(vec_TLV[i_cut][i_rec].Pz()); // pz
		    vec_two_prong_info[z_bin][centrality_bin].vec_data[buffer_bin][i_channel][3+30*i_cut].push_back(vec_TLV[i_cut][i_rec].M()); // M
		    for(Int_t i_info_track = 0; i_info_track < 26; i_info_track++)
		    {
			vec_two_prong_info[z_bin][centrality_bin].vec_data[buffer_bin][i_channel][4+i_info_track+30*i_cut].push_back(vec_track_id[i_cut][i_info_track][i_rec]); // track_id_A
		    }
		    //cout << "track_id_A: " << vec_track_id[0][i_rec] << ", track_id_B: " << vec_track_id[1][i_rec] << endl;

		}
	    }

	    std::vector< std::vector<TLorentzVector> > vec_TLV_AB;
	    std::vector< std::vector<TLorentzVector> > vec_TLV_ABC = Rec_three_prong_primary_decay(vec_track_info,vec_two_prong_info,Three_prong_channel_AB[i_channel],
												   N_cuts,Three_prong_pid_C[i_channel],z_bin,centrality_bin,buffer_bin,
												   buffer_bin_ME,SE_ME_flag_array[i_channel],vec_TLV_AB);
	    for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
	    {
		for(Int_t i_rec = 0; i_rec < vec_TLV_ABC[i_cut].size(); i_rec++)
		{
		    Double_t pT_rec_ABC         = vec_TLV_ABC[i_cut][i_rec].Pt();
		    Double_t invariant_mass_ABC = vec_TLV_ABC[i_cut][i_rec].M();
		    TLorentzVector TLV_ABC_nom  = vec_TLV_ABC[i_cut][i_rec];
		    TLorentzVector TLV_C        = vec_TLV_ABC[i_cut][i_rec] - vec_TLV_AB[i_cut][i_rec];
		    TLorentzVector TLV_AB_nom;
                    TLV_AB_nom.SetXYZM(vec_TLV_AB[i_cut][i_rec].Px(),vec_TLV_AB[i_cut][i_rec].Py(),
                                       vec_TLV_AB[i_cut][i_rec].Pz(),mass_D0);
		    TLV_ABC_nom = TLV_AB_nom + TLV_C;
		    Double_t invariant_mass_ABC_nom = TLV_ABC_nom.M();

		    Double_t pT_rec_AB          = vec_TLV_AB[i_cut][i_rec].Pt();
		    Double_t invariant_mass_AB  = vec_TLV_AB[i_cut][i_rec].M();

		    Double_t delta_invariant_mass = invariant_mass_ABC - invariant_mass_AB;

		    Int_t pT_bin = -1;
		    for(Int_t i_pT = 0; i_pT < N_pT_bins+1; i_pT++)
		    {
			if(pT_rec_ABC > pT_ranges[i_pT] && pT_rec_ABC <= pT_ranges[i_pT+1])
			{
			    pT_bin = i_pT;
			    break;
			}
		    }
		    if(pT_bin >= 0)
		    {
			h_invariant_mass_DStar_nom[i_channel][centrality_bin][pT_bin][i_cut] ->Fill(invariant_mass_ABC_nom);
			h_invariant_mass_DStar[i_channel][centrality_bin][pT_bin][i_cut]     ->Fill(invariant_mass_ABC);
			h_delta_invariant_mass[i_channel][centrality_bin][pT_bin][i_cut]     ->Fill(delta_invariant_mass);
		    }
		}
	    }
	} // end of i_channel


	if(flag_do_ME) // clear the used buffers
	{
	    for(Int_t i_pid = 0; i_pid < 8; i_pid++) // e-, pi-, K-, p-, e+, pi+, K+, p+
	    {
		for(Int_t i_info = 0; i_info < 17; i_info++)
		{
		    vec_track_info[z_bin][centrality_bin].vec_data[buffer_bin_ME][i_pid][i_info].clear();
		}
	    }
	    for(Int_t i_channel = 0; i_channel < 12; i_channel++)
	    {
		for(Int_t i_info = 0; i_info < (N_cuts*30); i_info++)
		{
		    vec_two_prong_info[z_bin][centrality_bin].vec_data[buffer_bin_ME][i_channel][i_info].clear();
		}
	    }
	}
	//-----------------------------------------------------------------

	Tree_AS_Event ->Fill();

	N_good_events++;
    }
}



//________________________________________________________________________
void AliAnalysisTaskDStartoKePi::Terminate(Option_t *)
{
    cout << "In terminate" << endl;
}


//________________________________________________________________________
void AliAnalysisTaskDStartoKePi::FillHelix(AliESDtrack* track_in, Double_t magF_in)
{
    //-------------------
    // Get helix
    // Track parametrization:
    // https://www.physi.uni-heidelberg.de/~sma/alice/LukasLayer_bachelor.pdf
    Double_t x_param,p_param[5]; // p are the track paramters, 0 = Y, 1 = Z, 2 = Snp, 3 = Tgl (p[2]/pt), 4 = Signed1Pt
    track_in->GetExternalParameters(x_param,p_param); // Those are the parameters used for vertexing
    Double_t fX     = track_in->GetX();

    //-------------------
    // Correct way of filling aliHelix from
    // http://personalpages.to.infn.it/~puccio/htmldoc/src/AliHelix.cxx.html#PalT1E
    // line 52

    // CONCLUSION: AliTracker::GetBz() is not identical to GetC(magF_in), magnetic field is slightly different but it doesn't matter...

    Double_t fHelix_alt[9];
    Double_t alpha_alt,x_alt,cs_alt,sn_alt;
    //track_in->GetExternalParameters(x_alt,fHelix_alt); // Those are the parameters used for vertexing
    x_alt = x_param;

    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
	fHelix_alt[i_param] = p_param[i_param];
    }


    alpha_alt=track_in->GetAlpha();
    //
    //circle parameters
    //PH Sometimes fP4 and fHelix[4] are very big and the calculation
    //PH of the Sqrt cannot be done. To be investigated...

    // kB2C=-0.299792458e-3; // from /AliRoot/STEER/STEERBase/AliVParticle.h
    //Double_t kB2C_test =-0.299792458e-3;
    //Double_t GetC(Double_t b) const
    //{return fP[4]*b*kB2C;}

    //Double_t par4test = fHelix_alt[4]*magF_in*kB2C_test;
    fHelix_alt[4] = track_in->GetC(magF_in);
    //fHelix_alt[4] = par4test; // take the one with the magnetic field directly from the ESD file

    cs_alt = TMath::Cos(alpha_alt);
    sn_alt = TMath::Sin(alpha_alt);

    Double_t xc_alt, yc_alt, rc_alt;
    rc_alt  =  1/fHelix_alt[4];
    xc_alt  =  x_alt-fHelix_alt[2]*rc_alt;
    Double_t dummy = 1-(x_alt-xc_alt)*(x_alt-xc_alt)*fHelix_alt[4]*fHelix_alt[4];
    yc_alt  =  fHelix_alt[0]+TMath::Sqrt(dummy)/fHelix_alt[4];

    fHelix_alt[6] = xc_alt*cs_alt - yc_alt*sn_alt;
    fHelix_alt[7] = xc_alt*sn_alt + yc_alt*cs_alt;
    fHelix_alt[8] =  TMath::Abs(rc_alt);
    //
    //
    fHelix_alt[5]=x_alt*cs_alt - fHelix_alt[0]*sn_alt;            // x0
    fHelix_alt[0]=x_alt*sn_alt + fHelix_alt[0]*cs_alt;            // y0
    fHelix_alt[2]=TMath::ATan2(-(fHelix_alt[5]-fHelix_alt[6]),fHelix_alt[0]-fHelix_alt[7]); // phi0
    if (fHelix_alt[4]>0) fHelix_alt[2]-=TMath::Pi();
    fHelix_alt[5]   = fHelix_alt[6];
    fHelix_alt[0]   = fHelix_alt[7];

    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
	aliHelix.fHelix[i_param] = fHelix_alt[i_param];
    }
    //-------------------


    //const AliExternalTrackParam* inner_param = track_in->GetInnerParam();
    const AliExternalTrackParam* inner_param = track_in->GetOuterParam();
    //Double_t helix_param[6],b_inner; // p are the track paramters, 0 = Y, 1 = Z, 2 = Snp, 3 = Tgl (p[2]/pt), 4 = Signed1Pt
    //inner_param ->GetHelixParameters(helix_param,b_inner);
    const Double_t* helix_param = inner_param ->GetParameter();

    //printf("p1: {%f,%f}, p2: {%f,%f}, p3: {%f,%f}, p4: {%f,%f}, p5: {%f,%f}, p6: {%f,%f} \n",helix_param[0],fHelix_alt[0],helix_param[1],fHelix_alt[1],helix_param[2],fHelix_alt[2],helix_param[3],fHelix_alt[3],helix_param[4],fHelix_alt[4],helix_param[5],fHelix_alt[5]);
    //printf("p1: {%f,%f}, p2: {%f,%f}, p3: {%f,%f}, p4: {%f,%f}, p5: {%f,%f}, p6: {%f,%f} \n",helix_param[0],p_param[0],helix_param[1],p_param[1],helix_param[2],p_param[2],helix_param[3],p_param[3],helix_param[4],p_param[4]);

}


//________________________________________________________________________
void AliAnalysisTaskDStartoKePi::FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Double_t path_initA, Double_t path_initB,
						  Double_t &pathA, Double_t &dca_xy, Double_t &dca_z)
{
    // V1.0
    Double_t pA[2] = {path_initA,path_initB}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Double_t distarray[2];
    TVector3 testA;
    TVector3 vec_dist[2];
    for(Int_t r = 0; r < 2; r++)
    {
	Double_t helix_point[3];
	helixA.Evaluate(pA[r],helix_point);
	testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[r]
	distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Double_t scale = 1.0;
    Double_t flip  = 1.0; // checks if the minimization direction changed
    Double_t scale_length = 60.0;
    while(fabs(scale_length) > 0.0001 && loopcounter < 100) // stops when the length is too small
    {
	//cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
	//    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
	//    << ", d[1] = " << distarray[1] << ", flip = " << flip
	//    << ", scale_length = " << scale_length << endl;
	if(distarray[0] > distarray[1])
	{
	    if(loopcounter != 0)
	    {
		if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
		else scale = 0.7; // go on in this direction but only by the way * 0.7
	    }
	    scale_length = (pA[1]-pA[0])*scale; // the next length interval
	    pA[0]     = pA[1] + scale_length; // the new path

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[0],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[0]
	    distarray[0] = (testA-space_vec).Mag(); // new dca
            vec_dist[0]  = testA-space_vec;
	    flip = 1.0;
	}
	else
	{
	    if(loopcounter != 0)
	    {
		if(flip == -1.0) scale = 0.4;
		else scale = 0.7;
	    }
	    scale_length = (pA[0]-pA[1])*scale;
	    pA[1]     = pA[0] + scale_length;

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[1],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[1]
	    distarray[1] = (testA-space_vec).Mag();
            vec_dist[1]  = testA-space_vec;
	    flip = -1.0;
	}
	loopcounter++;
    }

    if(loopcounter >= 100) cout << "WARNING: FindDCAHelixPoint exceeded maximum of 100 loops" << endl;

    if(distarray[0] < distarray[1])
    {
	pathA = pA[0];
	//dcaAB = distarray[0];
	dca_xy = TMath::Sqrt(vec_dist[0].X()*vec_dist[0].X() + vec_dist[0].Y()*vec_dist[0].Y());
	dca_z  = vec_dist[0].Z();
    }
    else
    {
	pathA = pA[1];
	//dcaAB = distarray[1];
        dca_xy = TMath::Sqrt(vec_dist[1].X()*vec_dist[1].X() + vec_dist[1].Y()*vec_dist[1].Y());
        dca_z = vec_dist[1].Z();
    }
}


Double_t AliAnalysisTaskDStartoKePi::calcDeterminant(TVector3& v1,TVector3& v2,TVector3& v3)
{
  // calculating the Determinant of a 3 x 3 Matrix 
  // with the column vectors [v1, v2, v3]
  // using the RULE of SARRUS
  //
  // | v1(0)   v2(0)   v3(0) |      | v1(0) v2(0) v3(0)|v1(0) v2(0)  .
  // |                       |      |  \\     \\     X   |  /     /    . 
  // |                       |      |   \\     \\   / \\  | /     /     . 
  // |                       |      |    \\     \\ /   \\ |/     /      . 
  // |                       |      |     \\     X     \\/     /       . 
  // |                       |      |      \\   / \\    /\\    /        .  
  // |                       |      |       \\ /   \\  / |\\  /         . 
  // | v1(1)   v2(1)   v3(1) |   =  | v1(1) v2(1) v3(1)|v1(1) v2(1)  .
  // |                       |      |       / \\    /\\  | /\\          . 
  // |                       |      |      /   \\  /  \\ |/  \\         . 
  // |                       |      |     /     \\/    \\/    \\        . 
  // |                       |      |    /      /\\    /\\     \\       . 
  // |                       |      |   /      /  \\  / |\\     \\      .  
  // |                       |      |  /      /    \\/  | \\     \\     . 
  // | v1(2)   v2(2)   v3(2) |      | v1(2) v2(2) v3(2)| v1(2) v2(2) .  
  //                                 /      /     /  \\     \\     \\   .
  //                                                                
  //                                -      -     -    +     +     +  .

  return ( v1(0) * v2(1) * v3(2) 
	   + v2(0) * v3(1) * v1(2) 
	   + v3(0) * v1(1) * v2(2) 
	   - v3(0) * v2(1) * v1(2) 
	   - v1(0) * v3(1) * v2(2) 
	   - v2(0) * v1(1) * v3(2)); 
}


TVector3 AliAnalysisTaskDStartoKePi::calculatePointOfClosestApproach(TVector3 &base1, TVector3 &dir1,
								    TVector3 &base2, TVector3 &dir2)
{
  //  calculating point of closest approach
  //        
  //        from the equations of the straight lines of g and h 
  //        g: x1 = base1 + l * dir1 
  //        h: x2 = base2 + m * dir2 
  //        
  //        you can construct the following planes:
  //        
  //        E1: e1 = base1  +  a * dir1  +  b * (dir1 x dir2)
  //        E2: e2 = base2  +  s * dir2  +  t * (dir1 x dir2)
  //        
  //        now the intersection point of E1 with g2 = {P1} 
  //        and the intersection point of E2 with g1 = {P2}
  //        
  //        form the base points of the perpendicular to both straight lines.
  //        
  //        The point of closest approach is the middle point between P1 and P2: 
  //        
  //        vertex = (p2 - p1)/2
  // 
  //        E1 ^ g2:
  //
  //           e1 = x2
  //    -->    base1  +  a * dir1  +  b * (dir1 x dir2) = base2 + m * dir2 
  //    -->    base1 - base2 = m * dir2  -  a * dir1  -  b * (dir1 x dir2)       
  //                                          (m)
  //    -->    [ dir2, -dir1, -(dir1 x dir2)] (a) = base1 - base2        
  //                                          (b)
  //           
  //           using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D12 = det [dir2, -dir1, -(dir1 x dir2)] 
  //               = det [dir2,  dir1,  (dir1 x dir2)]
  //           
  //           Dm  = det [base1 - base2, -dir1, -(dir1 x dir2)]
  //               = det [base1 - base2,  dir1,  (dir1 x dir2)]
  //  
  //            m  = Dm/D12
  //           
  //           P1: p1 = x2(m)
  //                  = base2 + Dm/D12 * dir2
  //
  //        E2 ^ g1:
  //
  //           e2 = x1
  //    -->    base2  +  s * dir2  +  t * (dir1 x dir2) = base1 + l * dir1 
  //    -->    base2 - base1 = l * dir1  -  s * dir2  -  t * (dir1 x dir2)       
  //                                          (l)
  //    -->    [ dir1, -dir2, -(dir1 x dir2)] (s) = base2 - base1        
  //                                          (t)
  //           
  //           again using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D21 =  det [dir1, -dir2, -(dir1 x dir2)] 
  //               =  det [dir1,  dir2,  (dir1 x dir2)]
  //               = -det [dir2,  dir1,  (dir1 x dir2)]
  //               = -D12
  //           
  //           Dl  =  det [base2 - base1, -dir2, -(dir1 x dir2)]
  //               =  det [base2 - base1,  dir1,  (dir1 x dir2)]
  //               = -det [base1 - base2,  dir1,  (dir1 x dir2)]
  //
  //            l  =   Dl/D21
  //               = - Dl/D12
  //           
  //           P2: p2 = x1(m)
  //                  = base1 - Dl/D12 * dir1
  //           
  //           
  //           vertex = p1 + 1/2 * (p2 - p1)
  //                  = 1/2 * (p2 + p1)
  //                  = 1/2 *( (base1 + base2) +  1/D12 * ( Dm * dir2 - Dl * dir1) )
  //                      

  TVector3 Cross = dir1.Cross(dir2); // Cross product: dir1 x dir2

  // straight lines are either skew or have a Cross point
	      
  TVector3 diff = base1;
  diff-=base2; // Difference of two base vectors base1 - base2
		
  Double_t D;
  D =  calcDeterminant(dir2, dir1 ,Cross);

  if (!(fabs(D) > 0.))
    {
      ::Warning(":calculatePointOfClosestApproach","Dirs and Cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");
      return TVector3(-20000.,-20000.,-20000.);
    }

  Double_t Dm =  calcDeterminant(diff , dir1, Cross);
  Double_t Dl = -calcDeterminant(diff , dir2, Cross);

  TVector3 vertex;
  TVector3 dm;
  TVector3 dl;

  dm = dir2;
  dm *= Dm;

  dl = dir1;
  dl *= Dl;

  vertex = dm - dl;

  vertex *= ((1.)/D);

  vertex+=base1;
  vertex+=base2;
  vertex*=0.5;

  return TVector3(vertex);
}



TVector3 AliAnalysisTaskDStartoKePi::calculateCrossPoint(TVector3 &base1, TVector3 &dir1,
							TVector3 &base2, TVector3 &dir2)
{ 
  // calculating Cross point 
  // taking all three equations into account solving the overdetermined set of lin. equations
  // of 
  // base1 + l * dir2 =  base1 + m * dir2 
  //
  // set of lin. equations:
  //  
  //   base1(0) + l * dir1(0) = base2(0) + m * dir2(0) 
  //   base1(1) + l * dir1(1) = base2(1) + m * dir2(1)
  //   base1(2) + l * dir1(2) = base2(2) + m * dir2(2) this line is ignored
  //
  //   written in matrix form
  //
  //        l
  //   M * |   | = base2 - base1
  //       \\ m /
  //
  //   M is a 3x2 matrix
  //     
  // to solve multiply the equation by the transposed Matrix of M from the left: M 
  //     
  //  T      /  l \\                                                               .
  // M * M * |    | = M  * (base2 - base1)
  //         \\ -m /
  // MIND THE '-' of m
  //
  //     / dir1(0) dir2(0) \\                                                      .
  //     |                 |    T   / dir1(0) dir1(1) dir1(2) \\                   .
  // M = | dir1(1) dir2(1) |,  M  = |                         |
  //     |                 |        \\ dir2(0) dir2(1) dir2(2) /                   .
  //     \\ dir1(2) dir2(2) /                                    
  //
  //  T      / (dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2))   (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))  \\ .
  // M * M = |                                                                                                                |
  //         \\ (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))   (dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2))  /                        
  //
  //  T       / d1d1 d1d2 \\                           .
  // M  * M = |           |
  //          \\ d1d2 d2d2 /
  //
  // diff = base2 - base1
  //
  //  T           /  (dir1(0)*diff(0) + dir1(1)*diff(1) + dir1(2)*diff(2)) \\         .
  // M  * diff =  |                                                        |
  //              \\  (dir2(0)*diff(0) + dir2(1)*diff(1) + dir2(2)*diff(2)) /
  //
  //  T           /  d1diff  \\                                          .
  // M  * diff =  |          |
  //              \\  d2diff  /
  // 
  // now the new Matrix set is to be solved by CRAMER'S Rule:
  // 
  // / d1d1 d1d2 \\   /  l \\   /  d1diff \\                   .
  // |           | * |    | = |          |
  // \\ d1d2 d2d2 /   \\ -m /   \\  d2diff /
  //
  //     | d1d1 d1d2 |
  // D = |           | = d1d1*d2d2 - d1d2*d1d2;
  //     | d1d2 d2d2 |
  // 
  //     | d1diff d1d2 |
  // Dl= |              | = d1diff*d2d2 - d1d2*d2diff;
  //     | d2diff d2d2 |              
  //
  // l = Dl/D = l_Cross
  // 
  // vertex = base1 + l_Cross * dir1
  //

  Double_t d1d1 = dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2);
  Double_t d2d2 = dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2);
  Double_t d1d2 = dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2);
  
  Double_t D = d1d1*d2d2 - (d1d2*d1d2);
  
  if (!(fabs(D) > 0.))
    {
      ::Warning("calculateCrossPoint","Error while calculating Cross point ... eqns are lin. dependent:returning default Vertex (-20000,-20000,-20000)");
      return TVector3(-20000.,-20000.,-20000.);
    }

  Double_t d1diff = dir1(0)*(base2(0)-base1(0))+dir1(1)*(base2(1)-base1(1))+dir1(2)*(base2(2)-base1(2));
  Double_t d2diff = dir2(0)*(base2(0)-base1(0))+dir2(1)*(base2(1)-base1(1))+dir2(2)*(base2(2)-base1(2));

  Double_t Dlambda = d1diff*d2d2-d1d2*d2diff;
  
  Double_t lambda = Dlambda/D;
  
  TVector3 vertex;
  vertex += dir1;
  vertex *= lambda;
  vertex += base1;

  //cout << "Cross point calculated" << endl;
  return TVector3(vertex);

 // return TVector3(-20000.,-20000.,-20000.);
}



TVector3 AliAnalysisTaskDStartoKePi::calcVertexAnalytical(TVector3 &base1, TVector3 &dir1,
							 TVector3 &base2, TVector3 &dir2)
{
  // Calculates the Vertex of two straight lines define by the vectors base and dir
  //
  //      g: x1 = base1 + l * dir1 
  //      h: x2 = base2 + m * dir2 , where l,m are real numbers 
  //
  // 1. are g and h
  //       parallel / identical, i.e. are dir1 and dir2 linear dependent?
  //       
  //                                        /-                               
  //                                        |
  //                                        |   = 0    linear dependent, no unique solution, returning dummy  
  //      => Cross product : dir1 x dir2 = -|  
  //                                        |  != 0    linear independent
  //                                        |
  //                                        \\-         
  //
  // 2. are g and h 
  //       skew or do they have a Crossing point, i.e are dir1, dir2 and (base1 - base2) linear dependent ?
  //
  //                                                    /-                               
  //                                                    |
  //                                                    |   = 0    linear dependent
  //                                                    |          g and h are intersecting
  //                                                    |          calculating vertex as point of intersection
  //                                                    |
  //    => determinant: det[ dir1, dir2, base1-base2]= -|
  //                                                    |  != 0    linear independent
  //                                                    |          g and h are skew
  //                                                    |          calulating vertex as point of closest approach
  //                                                    |
  //                                                    \\-         
  //  
  // 3.
  //    (a) calculating intersection point
  //    (b) calculating point of closest approach



  // 1. exists a unique solution ?

  if ((dir1.Cross(dir2)).Mag()> 0.) // dir1 and dir2 linear independent
    {
      // straight lines are either skew or have a Cross point

      TVector3 diff = base1;
      diff-=base2; // Difference of two base vectors base1 - base2
      
      // 2. skew or intersecting ?
	
      if (fabs(calcDeterminant(dir2, dir1 ,diff))>0.) 
	{
	  // 3. (b) skew 
	  return TVector3(calculatePointOfClosestApproach(base1, dir1, base2, dir2));
	}
      else
	{
	  // 3. (a) intersection 
	  return TVector3(calculateCrossPoint(base1 ,dir1, base2 ,dir2));
	}
    }
  else
    {
      // dir1 and dir2 linear dependent -> g1 and g2 identical or parallel
      return TVector3(-10000000.,-10000000.,-10000000.);
    }
  return TVector3(-10000000.,-10000000.,-10000000.);
}

Double_t AliAnalysisTaskDStartoKePi::calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
									 TVector3 &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.Mag()>0))
    {
      return -1000000.;
    }
  
  TVector3 diff = base-point;

  TVector3 cross = dir.Cross(diff);
  
  return cross.Mag()/dir.Mag();
}

Double_t AliAnalysisTaskDStartoKePi::calculateMinimumDistance(TVector3 &base1, TVector3 &dir1,
							  TVector3 &base2, TVector3 &dir2)
{
  // calculates the minimum distance of two tracks given as parametric straights x = base + n * dir

  TVector3 cross = dir1.Cross(dir2);

  TVector3 ab = base1 - base2;

  if ( !( fabs(cross.Mag())>0.)) // dir1 || dir2
    {
      return calculateMinimumDistanceStraightToPoint(base1, dir1, base2);
    }
 
  return fabs(ab.Dot(cross)/cross.Mag());
}

TVector3 AliAnalysisTaskDStartoKePi::calc_straight_to_point_dca_vector(TVector3 base1, TVector3 dir1,TVector3 point)
{
    TVector3 diff_vec = base1 - point;
    if(dir1.Mag() == 0.0)
    {
	cout << "WARNING: directrion vector in calc_straight_to_point_dca_vector is 0" << endl;
        return base1;
    }
    dir1 *= 1.0/dir1.Mag(); // normalize direction vector
    Double_t diff_proj_length = diff_vec.Dot(dir1);
    TVector3 proj_vec = dir1;
    proj_vec *= diff_proj_length;
    TVector3 dca_vec = diff_vec - proj_vec;
    return dca_vec;
}



//----------------------------------------------------------------------------------------
std::vector< std::vector<TLorentzVector> > AliAnalysisTaskDStartoKePi::Rec_two_prong_primary_decay(std::vector< std::vector<vectors_t> > vec_track_info_in,
								       Int_t Pid_A_in, Int_t Pid_B_in, Int_t z_bin, Int_t centrality_bin, Int_t buffer_bin, Int_t buffer_bin_ME, Int_t SE_ME_mode,
								       Int_t N_cuts, std::vector< std::vector< std::vector<Int_t> > > &vec_track_id)
{
    Double_t mass_array[8] = {0.00051099892,0.13957018,0.493677,0.93827203,0.00051099892,0.13957018,0.493677,0.93827203}; // pid:  e-, pi-, K-, p-, e+, pi+, K+, p+

    std::vector< std::vector<TLorentzVector> > vec_TLV_AB;
    vec_TLV_AB.resize(N_cuts);
    vec_track_id.resize(N_cuts);
    for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
    {
	vec_track_id[i_cut].resize(26); //
    }

    // pid:  e-, pi-, K-, p-, e+, pi+, K+, p+
    // SE_ME_mode: 0 = same event, 1 = mixed event

    for(Int_t i_ME_PID_flip = 0; i_ME_PID_flip < 2; i_ME_PID_flip++) // Switch PID_A and PID_B for mixed event mode
    {
	Int_t Pid_A = Pid_A_in;
	Int_t Pid_B = Pid_B_in;

	if(SE_ME_mode == 0 && i_ME_PID_flip == 1) continue; // not needed for same event mode
	if(i_ME_PID_flip == 1) // switch A and B for second event mixing loop
	{
	    Pid_A = Pid_B_in;
	    Pid_B = Pid_A_in;
	}
	Double_t mass_A = mass_array[Pid_A];
	Double_t mass_B = mass_array[Pid_B];

	//printf("Pid_A_in: %d, Pid_B_in: %d, Pid_A: %d, Pid_B: %d \n",Pid_A_in,Pid_B_in,Pid_A,Pid_B);

	for(Int_t i_buffer_event = 0; i_buffer_event < ME_event_buffer_size; i_buffer_event++) // loop for event mixing
	{
	    Int_t ref_event = buffer_bin;
	    if(SE_ME_mode == 1) ref_event = buffer_bin_ME;

	    if(SE_ME_mode == 0 && i_buffer_event != ref_event) continue; // use only the actual event for same event mode
	    if(SE_ME_mode == 1 && i_buffer_event == ref_event) continue; // don't use the same event for both particles in mixed event mode

	    Int_t n_tracks_A = vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][0].size();
	    Int_t n_tracks_B = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][0].size();

	    for(Int_t i_track_A = 0; i_track_A < n_tracks_A; i_track_A++)
	    {
		Int_t track_id_A = (Int_t)vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][7][i_track_A];
		TLorentzVector TLV_A;
                TLV_A.SetXYZM(vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][2][i_track_A],
                              vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][3][i_track_A],
                              vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][4][i_track_A],
                              mass_A);
		Double_t nSigma_dEdx_A = vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][0][i_track_A];
		Double_t nSigma_TOF_A  = vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][1][i_track_A];
		Int_t    ITS_refit_A   = (Int_t)vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][5][i_track_A];
		Double_t dca_A_xy      = vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][6][i_track_A];
		Double_t dca_A_z       = vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][15][i_track_A];
                Double_t ITS_cls_A     = vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][16][i_track_A];

		if(dca_A_xy > 1.0) continue; // Don't use global tracks too far away from primary vertex
		if(fabs(nSigma_dEdx_A) > 2.5) continue;
		if(nSigma_TOF_A > -100 && fabs(nSigma_TOF_A) > 2.5) continue;

		for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
		{
		    Int_t cut_A = i_cut;
		    Int_t cut_B = i_cut;

		    if(cut_A == 1)
		    {
			if(Pid_A == 0 || Pid_A == 4) // electrons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			}
			if(Pid_A == 1 || Pid_A == 5 || Pid_A == 2 || Pid_A == 6) // pions, kaons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() > 0.7 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID for high momenta
			}
			if(Pid_A == 3 || Pid_A == 7) // protons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() > 1.2 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID for high momenta
			}
		    }
		    if(cut_A == 2)
		    {
			if(Pid_A == 0 || Pid_A == 4) // electrons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() < 0.22 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID
			    if(TLV_A.P() > 0.45 && TLV_A.P() < 0.6 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID
			    if(TLV_A.P() > 0.9 && TLV_A.P() < 1.1 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID
			}
			//if(Pid_A == 1 || Pid_A == 5 || Pid_A == 2 || Pid_A == 6) // pions, kaons
			//{
			//    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			//    if(TLV_A.P() > 0.7 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID for high momenta
			//}
			if(Pid_A == 3 || Pid_A == 7) // protons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() > 1.2 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID for high momenta
			}
		    }
		    if(cut_A == 3) // This is for topological reconstruction
		    {
			if(!ITS_refit_A) continue;
			if(Pid_A == 0 || Pid_A == 4) // electrons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() < 0.22 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID
			    if(TLV_A.P() > 0.45 && TLV_A.P() < 0.6 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID
			    if(TLV_A.P() > 0.9 && TLV_A.P() < 1.1 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID
			}
			if(Pid_A == 1 || Pid_A == 5 || Pid_A == 2 || Pid_A == 6) // pions, kaons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() > 0.7 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID for high momenta
			}
			if(Pid_A == 3 || Pid_A == 7) // protons
			{
			    if(fabs(nSigma_dEdx_A) > 2.5) continue;
			    if(TLV_A.P() > 1.2 && fabs(nSigma_TOF_A) > 2.5) continue; // enforce TOF PID for high momenta
			}
		    }

		    for(Int_t i_track_B = 0; i_track_B < n_tracks_B; i_track_B++)
		    {
			Int_t track_id_B = (Int_t)vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][7][i_track_B];
			if(SE_ME_mode == 0 && track_id_A == track_id_B)  continue; // same event mode
			TLorentzVector TLV_B;
                        TLV_B.SetXYZM(vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][2][i_track_B],
                                      vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][3][i_track_B],
                                      vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][4][i_track_B],
                                      mass_B);
			Double_t nSigma_dEdx_B = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][0][i_track_B];
			Double_t nSigma_TOF_B  = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][1][i_track_B];
			Double_t ITS_refit_B   = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][5][i_track_B];
			Double_t dca_B_xy      = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][6][i_track_B];
			Double_t dca_B_z       = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][15][i_track_B];
                        Double_t ITS_cls_B     = vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][16][i_track_B];

			if(dca_B_xy > 1.0) continue;
			if(fabs(nSigma_dEdx_B) > 2.5) continue;
			if(nSigma_TOF_B > -100 && fabs(nSigma_TOF_B) > 2.5) continue;
			if(cut_B == 1)
			{
			    if(Pid_B == 0 || Pid_B == 4) // electrons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
			    }
			    if(Pid_B == 1 || Pid_B == 5 || Pid_B == 2 || Pid_B == 6) // pions, kaons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() > 0.7 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			    if(Pid_B == 3 || Pid_B == 7) // protons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() > 1.2 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			}
			if(cut_B == 2)
			{
			    if(Pid_B == 0 || Pid_B == 4) // electrons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() < 0.22 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID
				if(TLV_B.P() > 0.45 && TLV_B.P() < 0.6 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID
				if(TLV_B.P() > 0.9 && TLV_B.P() < 1.1 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID
			    }
			    if(Pid_B == 1 || Pid_B == 5 || Pid_B == 2 || Pid_B == 6) // pions, kaons
			    {
				//if(fabs(nSigma_dEdx_B) > 2.5) continue;
                                if((TLV_A.P() > 0.7 && fabs(nSigma_TOF_A) > 2.5) && (TLV_B.P() > 0.7 && fabs(nSigma_TOF_B) > 2.5)) continue; // one of the two needs to have a very good PID
				//if(TLV_B.P() > 0.7 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			    if(Pid_B == 3 || Pid_B == 7) // protons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() > 1.2 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			}
			if(cut_B == 3) // This is for topological reconstruction
			{
			    if(!ITS_refit_B) continue;
			    if(Pid_B == 0 || Pid_B == 4) // electrons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() < 0.22 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID
				if(TLV_B.P() > 0.45 && TLV_B.P() < 0.6 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID
				if(TLV_B.P() > 0.9 && TLV_B.P() < 1.1 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID
			    }
			    if(Pid_B == 1 || Pid_B == 5 || Pid_B == 2 || Pid_B == 6) // pions, kaons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() > 0.7 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			    if(Pid_B == 3 || Pid_B == 7) // protons
			    {
				if(fabs(nSigma_dEdx_B) > 2.5) continue;
				if(TLV_B.P() > 1.2 && fabs(nSigma_TOF_B) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			}

                        //----------------------------------------
			// Do topological reconstruction
			Double_t dist_to_prim_xy       = -1.0;
			Double_t dist_to_prim_z        = -1.0;
			Double_t dca_to_secondary_A_xy = -1.0;
			Double_t dca_to_secondary_A_z  = -1.0;
			Double_t dca_to_primary_AB_xy  = -1.0;
                        Double_t dca_to_primary_AB_z   = -1.0;

			if(cut_A == 3 && cut_B == 3)
			{
			    TVector3 base_A, base_B, dir_A, dir_B;
			    base_A.SetXYZ(vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][8][i_track_A],
					  vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][9][i_track_A],
					  vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][10][i_track_A]
					 );
			    dir_A.SetXYZ(vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][11][i_track_A],
					 vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][12][i_track_A],
					 vec_track_info_in[z_bin][centrality_bin].vec_data[ref_event][Pid_A][13][i_track_A]
					);
			    base_B.SetXYZ(vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][8][i_track_B],
					  vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][9][i_track_B],
					  vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][10][i_track_B]
					 );
			    dir_B.SetXYZ(vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][11][i_track_B],
					 vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][12][i_track_B],
					 vec_track_info_in[z_bin][centrality_bin].vec_data[i_buffer_event][Pid_B][13][i_track_B]
                                        );

			    if(dir_A.Mag() <= 0.0 || dir_B.Mag() <= 0.0) continue;
                            dir_A *= 1.0/dir_A.Mag();
                            dir_B *= 1.0/dir_B.Mag();


			    // Calculate decay vertex, linear approximation is good enough for a few mum decay length
			    Double_t dca_track_lin  = calculateMinimumDistance(base_A,dir_A,base_B,dir_B); // dca of the two tracks
			    if(dca_track_lin > 0.1) continue; // ITS z pad size is 400 mum

                            // Calculate secondary decay vertex
                            TVector3 vector_AB_lin  = calcVertexAnalytical(base_A,dir_A,base_B,dir_B); // vertex of the two tracks

                            // Calculate dca of particle A to secondary decay vertex
			    TVector3 dca_vec_to_secondary_A = calc_straight_to_point_dca_vector(base_A,dir_A,vector_AB_lin);
			    dca_to_secondary_A_xy  = TMath::Sqrt(dca_vec_to_secondary_A.X()*dca_vec_to_secondary_A.X() + dca_vec_to_secondary_A.Y()*dca_vec_to_secondary_A.Y());
			    dca_to_secondary_A_z   = dca_vec_to_secondary_A.Z();
			    if(dca_to_secondary_A_xy > 0.025) continue; // dca between track A and secondary vertex
			    //printf("dca_xy: %f, dca_z: %f \n",dca_to_secondary_A_xy,dca_to_secondary_A_z);

			    // Calculate dca of particle B to secondary decay vertex
			    TVector3 dca_vec_to_secondary_B = calc_straight_to_point_dca_vector(base_B,dir_B,vector_AB_lin);      \
                            Double_t dca_to_secondary_B_xy  = TMath::Sqrt(dca_vec_to_secondary_B.X()*dca_vec_to_secondary_B.X() + dca_vec_to_secondary_B.Y()*dca_vec_to_secondary_B.Y());
                            Double_t dca_to_secondary_B_z   = dca_vec_to_secondary_B.Z();
                            if(dca_to_secondary_B_xy > 0.025) continue; // dca between track B and secondary vertex

                            // Calculate decay distance
                            TVector3  vector_AB_to_prim = vector_AB_lin;
			    vector_AB_to_prim -= TV3_prim_vertex;
                            Double_t dist_to_prim    = vector_AB_to_prim.Mag();
			    dist_to_prim_xy = TMath::Sqrt(vector_AB_to_prim.X()*vector_AB_to_prim.X() + vector_AB_to_prim.Y()*vector_AB_to_prim.Y());
                            if(dist_to_prim_xy < 0.008) continue; // distance to primary decay vertex in xy
			    dist_to_prim_z  = fabs(vector_AB_to_prim.Z());
                            //if(dist_to_prim_z < 0.008) continue; // distance to primary decay vertex in z

                            // Set new Lorentz vectors
                            Double_t p_A = TLV_A.P();
                            TLV_A.SetXYZM(dir_A.X()*p_A,dir_A.Y()*p_A,dir_A.Z()*p_A,mass_A);
                            Double_t p_B = TLV_B.P();
                            TLV_B.SetXYZM(dir_B.X()*p_B,dir_B.Y()*p_B,dir_B.Z()*p_B,mass_B);
			    //printf("p_A: %f, p_B: %f, vec_A: {%f, %f, %f}, vec_B: {%f, %f, %f}, mass_A: %f, mass_B: %f \n",p_A,p_B,dir_A.X()*p_A,dir_A.Y()*p_A,dir_A.Z()*p_A,dir_B.X()*p_B,dir_B.Y()*p_B,dir_B.Z()*p_B,mass_A,mass_B);

                            // Calculate dca of mother particle to primary vertex
                            TLorentzVector TLV_AB_dir;
                            TLV_AB_dir = TLV_A + TLV_B;
                            TVector3 vec_dir_AB;
			    vec_dir_AB.SetXYZ(TLV_AB_dir.Px(),TLV_AB_dir.Py(),TLV_AB_dir.Pz());
			    //printf("vecAB: {%f, %f, %f} \n",vec_dir_AB.Px(),vec_dir_AB.Py(),vec_dir_AB.Pz());
			    TVector3 dca_vec_to_primary_AB = calc_straight_to_point_dca_vector(vector_AB_lin,vec_dir_AB,TV3_prim_vertex);
			    dca_to_primary_AB_xy  = TMath::Sqrt(dca_vec_to_primary_AB.X()*dca_vec_to_primary_AB.X() + dca_vec_to_primary_AB.Y()*dca_vec_to_primary_AB.Y());
			    dca_to_primary_AB_z   = dca_vec_to_primary_AB.Z();
			    if(dca_to_primary_AB_xy > 0.025) continue; // mother particle does not result of a primary production

                            // Check if emission direction is the same as the direction of primary to secondary vertex
                            TVector3 vec_dir_primary_to_secondary = vector_AB_lin;
                            vec_dir_primary_to_secondary -= TV3_prim_vertex;
                            Double_t dot_product = vec_dir_primary_to_secondary.Dot(vec_dir_AB);
                            if(dot_product < 0.0) continue; // direction is not the same
			}
                        //----------------------------------------

			TLorentzVector TLV_AB;
			TLV_AB = TLV_A + TLV_B;
			vec_TLV_AB[i_cut].push_back(TLV_AB);
                        vec_track_id[i_cut][0].push_back(track_id_A);
                        vec_track_id[i_cut][1].push_back(track_id_B);
                        vec_track_id[i_cut][2].push_back(nSigma_dEdx_A);
                        vec_track_id[i_cut][3].push_back(nSigma_TOF_A);
                        vec_track_id[i_cut][4].push_back(dca_A_xy);
                        vec_track_id[i_cut][5].push_back(dca_A_z);
			vec_track_id[i_cut][6].push_back(TLV_A.Px());
			vec_track_id[i_cut][7].push_back(TLV_A.Py());
			vec_track_id[i_cut][8].push_back(TLV_A.Pz());
			vec_track_id[i_cut][9].push_back(TLV_A.M());
			vec_track_id[i_cut][10].push_back(ITS_cls_A);
			vec_track_id[i_cut][11].push_back(nSigma_dEdx_B);
                        vec_track_id[i_cut][12].push_back(nSigma_TOF_B);
                        vec_track_id[i_cut][13].push_back(dca_B_xy);
                        vec_track_id[i_cut][14].push_back(dca_B_z);
			vec_track_id[i_cut][15].push_back(TLV_B.Px());
			vec_track_id[i_cut][16].push_back(TLV_B.Py());
			vec_track_id[i_cut][17].push_back(TLV_B.Pz());
			vec_track_id[i_cut][18].push_back(TLV_B.M());
			vec_track_id[i_cut][19].push_back(ITS_cls_B);
			vec_track_id[i_cut][20].push_back(dist_to_prim_xy);
			vec_track_id[i_cut][21].push_back(dist_to_prim_z);
                        vec_track_id[i_cut][22].push_back(dca_to_secondary_A_xy*2.0);
			vec_track_id[i_cut][23].push_back(dca_to_secondary_A_z*2.0);
			vec_track_id[i_cut][24].push_back(dca_to_primary_AB_xy);
			vec_track_id[i_cut][25].push_back(dca_to_primary_AB_z);

			//printf("mass_A: %f, mass_B: %f, invariant_mass: %f, Pid_A: %d, Pid_B: %d \n",mass_A, mass_B,TLV_AB.M(),Pid_A,Pid_B);
		    }
		}
	    }
	}
    }

    return vec_TLV_AB;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
std::vector< std::vector<TLorentzVector> > AliAnalysisTaskDStartoKePi::Rec_three_prong_primary_decay(std::vector< std::vector<vectors_t> > vec_track_info_in,
                                                                         std::vector< std::vector<vectors_t> > vec_two_prong_info_in,
									 Int_t i_channel,Int_t N_cuts, Int_t Pid_C_in, Int_t z_bin, Int_t centrality_bin, Int_t buffer_bin, Int_t buffer_bin_ME, Int_t SE_ME_mode,
									 std::vector< std::vector<TLorentzVector> > &vec_TLV_AB
									)
{
    Double_t mass_array[8] = {0.00051099892,0.13957018,0.493677,0.93827203,0.00051099892,0.13957018,0.493677,0.93827203}; // pid:  e-, pi-, K-, p-, e+, pi+, K+, p+

    std::vector< std::vector<TLorentzVector> > vec_TLV_ABC;
    vec_TLV_ABC.resize(N_cuts);
    vec_TLV_AB.resize(N_cuts);

    Int_t    Pid_C  = Pid_C_in;
    Double_t mass_C = mass_array[Pid_C];

    // pid:  e-, pi-, K-, p-, e+, pi+, K+, p+
    // SE_ME_mode: 0 = same event, 1 = mixed event

    for(Int_t i_ME_PID_flip = 0; i_ME_PID_flip < 2; i_ME_PID_flip++) // Switch PID_AB and PID_C for mixed event mode
    {
	if(SE_ME_mode == 0 && i_ME_PID_flip == 1) continue; // not needed for same event mode
	for(Int_t i_buffer_event = 0; i_buffer_event < ME_event_buffer_size; i_buffer_event++) // loop for event mixing
	{
	    Int_t ref_event = buffer_bin;
	    if(SE_ME_mode == 1) ref_event = buffer_bin_ME; // mixed event mode

	    if(SE_ME_mode == 0 && i_buffer_event != ref_event) continue; // use only the actual event for same event mode
	    if(SE_ME_mode == 1 && i_buffer_event == ref_event) continue; // don't use the same event for both particles in mixed event mode

	    Int_t event_AB = i_buffer_event;
	    Int_t event_C  = ref_event;
	    if(i_ME_PID_flip == 1)
	    {
		event_AB = ref_event;
		event_C  = i_buffer_event;
	    }

            Int_t n_tracks_C  = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][0].size();
	    for(Int_t i_cut = 0; i_cut < N_cuts; i_cut++)
	    {
		Int_t n_tracks_AB = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][0+30*i_cut].size();

		//if(i_channel <= 1) printf("i_buffer_event: %d, i_channel: %d, n_track_AB: %d, n_track_C: %d \n",i_buffer_event, i_channel, n_tracks_AB, n_tracks_C);
		//if(SE_ME_mode == 1) cout << "centrality_bin: " << centrality_bin << ", i_ME_PID_flip: " << i_ME_PID_flip << ", n_tracks_AB: " << n_tracks_AB << ", n_tracks_C: " << n_tracks_C << endl;

		for(Int_t i_track_AB = 0; i_track_AB < n_tracks_AB; i_track_AB++)
		{
		    TLorentzVector TLV_AB;
		    TLV_AB.SetXYZM(vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][0+30*i_cut][i_track_AB],
				   vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][1+30*i_cut][i_track_AB],
				   vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][2+30*i_cut][i_track_AB],
				   vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][3+30*i_cut][i_track_AB]);


		    Int_t track_id_A = (Int_t)vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+0+30*i_cut][i_track_AB];
		    Int_t track_id_B = (Int_t)vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+1+30*i_cut][i_track_AB];

		    for(Int_t i_track_C = 0; i_track_C < n_tracks_C; i_track_C++)
		    {
			Int_t track_id_C = (Int_t)vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][7][i_track_C];
			if(SE_ME_mode == 0 && track_id_A == track_id_C) continue; // Don't use the same track twice in same event mode
			if(SE_ME_mode == 0 && track_id_B == track_id_C) continue;
			TLorentzVector TLV_C;

			Double_t nSigma_dEdx_C = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][0][i_track_C];
			Double_t nSigma_TOF_C  = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][1][i_track_C];
			Double_t ITS_refit_C   = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][5][i_track_C];
			Double_t dca_C_xy      = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][6][i_track_C];
			Double_t dca_C_z       = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][15][i_track_C];
                        Double_t ITS_cls_C     = vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][16][i_track_C];


			if(dca_C_xy > 1.0 || (ITS_refit_C == 1 && dca_C_xy > 0.1)) continue;
			if(fabs(nSigma_dEdx_C) > 2.5) continue;
			if(nSigma_TOF_C > -100 && fabs(nSigma_TOF_C) > 2.5) continue;

			TLV_C.SetXYZM(vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][2][i_track_C],
				      vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][3][i_track_C],
				      vec_track_info_in[z_bin][centrality_bin].vec_data[event_C][Pid_C][4][i_track_C],
				      mass_C);

			Int_t cut_C = i_cut;

			if(cut_C == 1)
			{
			    if(Pid_C == 0 || Pid_C == 4) // electrons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
			    }
			    if(Pid_C == 1 || Pid_C == 5 || Pid_C == 2 || Pid_C == 6) // pions, kaons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() > 0.7 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			    if(Pid_C == 3 || Pid_C == 7) // protons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() > 1.2 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			}
			if(cut_C == 2)
			{
			    if(Pid_C == 0 || Pid_C == 4) // electrons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() < 0.22 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID
				if(TLV_C.P() > 0.45 && TLV_C.P() < 0.6 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID
				if(TLV_C.P() > 0.9 && TLV_C.P() < 1.1 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID
			    }
			    if(Pid_C == 1 || Pid_C == 5 || Pid_C == 2 || Pid_C == 6) // pions, kaons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() > 0.7 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			    if(Pid_C == 3 || Pid_C == 7) // protons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() > 1.2 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			}
			if(cut_C == 3) // This is for topological reconstruction
			{
			    if(!ITS_refit_C) continue;
			    if(Pid_C == 0 || Pid_C == 4) // electrons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() < 0.22 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID
				if(TLV_C.P() > 0.45 && TLV_C.P() < 0.6 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID
				if(TLV_C.P() > 0.9 && TLV_C.P() < 1.1 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID
			    }
			    if(Pid_C == 1 || Pid_C == 5 || Pid_C == 2 || Pid_C == 6) // pions, kaons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() > 0.7 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			    if(Pid_C == 3 || Pid_C == 7) // protons
			    {
				if(fabs(nSigma_dEdx_C) > 2.5) continue;
				if(TLV_C.P() > 1.2 && fabs(nSigma_TOF_C) > 2.5) continue; // enforce TOF PID for high momenta
			    }
			}

			TLorentzVector TLV_ABC;
			TLV_ABC = TLV_AB + TLV_C;
			vec_TLV_ABC[i_cut].push_back(TLV_ABC);
			vec_TLV_AB[i_cut].push_back(TLV_AB);

			if(cut_C == 3)
			{
			    TLorentzVector TLV_A, TLV_B;
			    TLV_A.SetXYZM(vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+6+30*i_cut][i_track_AB],
					  vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+7+30*i_cut][i_track_AB],
					  vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+8+30*i_cut][i_track_AB],
                                          vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+9+30*i_cut][i_track_AB]
					 );
			    TLV_B.SetXYZM(vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+15+30*i_cut][i_track_AB],
					  vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+16+30*i_cut][i_track_AB],
					  vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+17+30*i_cut][i_track_AB],
                                          vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+18+30*i_cut][i_track_AB]
					 );

			    Double_t dist_to_prim_xy = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+20+30*i_cut][i_track_AB];
                            Double_t dist_to_prim_z  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+21+30*i_cut][i_track_AB];

			    Double_t dca_AB_xy = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+22+30*i_cut][i_track_AB];
			    Double_t dca_AB_z  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+23+30*i_cut][i_track_AB];

			    Double_t dca_AB_to_prim_xy = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+24+30*i_cut][i_track_AB];
                            Double_t dca_AB_to_prim_z  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+25+30*i_cut][i_track_AB];

			    Double_t dca_A_xy = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+4+30*i_cut][i_track_AB];
			    Double_t dca_A_z  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+5+30*i_cut][i_track_AB];

			    Double_t dca_B_xy = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+13+30*i_cut][i_track_AB];
                            Double_t dca_B_z  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+14+30*i_cut][i_track_AB];

                            Double_t nSigma_dEdx_A = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+2+30*i_cut][i_track_AB];
			    Double_t nSigma_TOF_A  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+3+30*i_cut][i_track_AB];

			    Double_t nSigma_dEdx_B = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+11+30*i_cut][i_track_AB];
                            Double_t nSigma_TOF_B  = vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+12+30*i_cut][i_track_AB];

			    UShort_t ITS_cls_A = (UShort_t)vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+10+30*i_cut][i_track_AB];
                            UShort_t ITS_cls_B = (UShort_t)vec_two_prong_info_in[z_bin][centrality_bin].vec_data[event_AB][i_channel][4+19+30*i_cut][i_track_AB];

			    AS_Track  = AS_Event ->createTrack();
			    AS_Track  ->set_channel(i_channel);
			    AS_Track  ->set_TLV_part_A(TLV_A);
			    AS_Track  ->set_TLV_part_B(TLV_B);
			    AS_Track  ->set_TLV_part_C(TLV_C);
			    AS_Track  ->set_dist_to_prim(dist_to_prim_xy, dist_to_prim_z);
			    AS_Track  ->set_dca_AB(dca_AB_xy,dca_AB_z);
			    AS_Track  ->set_dca_AB_to_prim(dca_AB_to_prim_xy,dca_AB_to_prim_z);
			    AS_Track  ->set_dca_A(dca_A_xy,dca_A_z);
			    AS_Track  ->set_dca_B(dca_B_xy,dca_B_z);
			    AS_Track  ->set_dca_C(dca_C_xy,dca_C_z);
			    AS_Track  ->set_Pid_A(nSigma_dEdx_A,nSigma_TOF_A);
			    AS_Track  ->set_Pid_B(nSigma_dEdx_B,nSigma_TOF_B);
			    AS_Track  ->set_Pid_C(nSigma_dEdx_C,nSigma_TOF_C);
			    AS_Track  ->set_NITScls_A(ITS_cls_A);
			    AS_Track  ->set_NITScls_B(ITS_cls_B);
			    AS_Track  ->set_NITScls_C(ITS_cls_C);
			}

		    }
		}
	    }
	}
    }

    return vec_TLV_ABC;
}
//----------------------------------------------------------------------------------------


