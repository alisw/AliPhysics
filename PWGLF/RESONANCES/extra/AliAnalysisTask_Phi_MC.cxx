// !TODO LIST

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TTree.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODHeader.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisTask_Phi_MC.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVCluster.h"

class AliAnalysisTask_Phi_MC;

using namespace std;

ClassImp(AliAnalysisTask_Phi_MC)

AliAnalysisTask_Phi_MC::AliAnalysisTask_Phi_MC() : AliAnalysisTaskSE(), kComputeSpherocity(0), kSpherocityPTWeight(0), kComputeRT(0), fIs_p_p(0), fIs_p_Pb(0), fIs_Pb_Pb(0), fRunName(0), fQC_Event_Enum_FLL(0), fMCD(0), AODMCTrackArray(0), fMultSelection(0), fCurrent_SPH(0), fCurrent_V0M(0), fCurrent_TRK(0), fCurrent_RT(0), fCurrent_Run(0), fTrueEventMask(0), tParticle(0), fTNParticle(0), lQCParticle(0), lAnalysisOutputList(0)  {
            }

//_____________________________________________________________________________

            AliAnalysisTask_Phi_MC::AliAnalysisTask_Phi_MC(const char* name) : AliAnalysisTaskSE(name), kComputeSpherocity(0), kSpherocityPTWeight(0), kComputeRT(0), fIs_p_p(0), fIs_p_Pb(0), fIs_Pb_Pb(0), fRunName(0), fQC_Event_Enum_FLL(0), fMCD(0), AODMCTrackArray(0), fMultSelection(0), fCurrent_SPH(0), fCurrent_V0M(0), fCurrent_TRK(0), fCurrent_RT(0), fCurrent_Run(0), fTrueEventMask(0), tParticle(0), fTNParticle(0), lQCParticle(0), lAnalysisOutputList(0) {
                // Define Input
                DefineInput(0, TChain::Class());

                // Define Output
                DefineOutput(1, TList::Class());
                DefineOutput(2, TList::Class());
                DefineOutput(3, TTree::Class());
            }

//_____________________________________________________________________________

            AliAnalysisTask_Phi_MC::~AliAnalysisTask_Phi_MC()                 {
                // Deleting TLists
                if  ( lQCParticle )    delete lQCParticle;
                if  ( lAnalysisOutputList ) delete lAnalysisOutputList;
                // Deleting TTrees
                if  ( tParticle )      delete tParticle;
            }

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::UserCreateOutputObjects()                  {
    //
    //  --- Analysis Output List
    lAnalysisOutputList =   new TList();
    lAnalysisOutputList ->  SetOwner( kTRUE );
    PostData ( 1, lAnalysisOutputList );
    //
    //  --- Particle QC Histograms List
    lQCParticle    =   new TList();
    lQCParticle    ->  SetOwner( kTRUE );
    //
    fQC_Event_Enum_FLL              = new TH1D("fQC_Event_Enum_FLL",    "Event Selection",                                      29, -0.5, 28.5);
    fQC_Event_Enum_FLL              ->  GetXaxis()  ->  SetTitle("");
    fQC_Event_Enum_FLL              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    uSetEventCountLabels(fQC_Event_Enum_FLL);
    lQCParticle    ->  Add(fQC_Event_Enum_FLL);
    //
    fQC_Event_Enum_E05              = new TH1D("fQC_Event_Enum_E05",    "Events with tracks in #eta < 0.5",                     5000,   0., 5000.);
    fQC_Event_Enum_E05              ->  GetXaxis()  ->  SetTitle("N Tracks");
    fQC_Event_Enum_E05              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    lQCParticle    ->  Add(fQC_Event_Enum_E05);
    //
    fQC_Event_Enum_E08              = new TH1D("fQC_Event_Enum_E08",    "Events with tracks in #eta < 0.8",                     5000,   0., 5000.);
    fQC_Event_Enum_E08              ->  GetXaxis()  ->  SetTitle("N Tracks");
    fQC_Event_Enum_E08              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    lQCParticle    ->  Add(fQC_Event_Enum_E08);
    //
    fQC_Event_Enum_E10              = new TH1D("fQC_Event_Enum_E10",    "Events with tracks in #eta < 1.0",                     5000,   0., 5000.);
    fQC_Event_Enum_E10              ->  GetXaxis()  ->  SetTitle("N Tracks");
    fQC_Event_Enum_E10              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    lQCParticle    ->  Add(fQC_Event_Enum_E10);
    //
    fQC_Event_Enum_V0A              = new TH1D("fQC_Event_Enum_V0A",    "Events with tracks in V0A acceptance",                 5000,   0., 5000.);
    fQC_Event_Enum_V0A              ->  GetXaxis()  ->  SetTitle("N Tracks");
    fQC_Event_Enum_V0A              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    lQCParticle    ->  Add(fQC_Event_Enum_V0A);
    //
    fQC_Event_Enum_V0A_05           = new TH2D("fQC_Event_Enum_V0A_05", "Events with tracks in V0A acceptance vs Mult Eta05",   5000,   0., 5000., 5000,   0., 5000.);
    fQC_Event_Enum_V0A_05           ->  GetXaxis()  ->  SetTitle("N Tracks V0A");
    fQC_Event_Enum_V0A_05           ->  GetYaxis()  ->  SetTitle("N Tracks |#eta| < 0.5");
    lQCParticle    ->  Add(fQC_Event_Enum_V0A_05);
    //
    fQC_Event_Enum_V0A_08           = new TH2D("fQC_Event_Enum_V0A_08", "Events with tracks in V0A acceptance vs Mult Eta08",   5000,   0., 5000., 5000,   0., 5000.);
    fQC_Event_Enum_V0A_08           ->  GetXaxis()  ->  SetTitle("N Tracks V0A");
    fQC_Event_Enum_V0A_08           ->  GetYaxis()  ->  SetTitle("N Tracks |#eta| < 0.8");
    lQCParticle    ->  Add(fQC_Event_Enum_V0A_08);
    //
    fQC_Event_Enum_V0A_10           = new TH2D("fQC_Event_Enum_V0A_10", "Events with tracks in V0A acceptance vs Mult Eta10",   5000,   0., 5000., 5000,   0., 5000.);
    fQC_Event_Enum_V0A_10           ->  GetXaxis()  ->  SetTitle("N Tracks V0A");
    fQC_Event_Enum_V0A_10           ->  GetYaxis()  ->  SetTitle("N Tracks |#eta| < 1.0");
    lQCParticle    ->  Add(fQC_Event_Enum_V0A_10);
    //
    fQC_Event_Enum_V0M              = new TH1D("fQC_Event_Enum_V0M",    "Events with tracks in V0M acceptance",                 5000,   0., 5000.);
    fQC_Event_Enum_V0M              ->  GetXaxis()  ->  SetTitle("N Tracks");
    fQC_Event_Enum_V0M              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    lQCParticle    ->  Add(fQC_Event_Enum_V0M);
    //
    fQC_Event_Enum_V0M_05           = new TH2D("fQC_Event_Enum_V0M_05", "Events with tracks in V0M acceptance vs Mult Eta05",   5000,   0., 5000., 5000,   0., 5000.);
    fQC_Event_Enum_V0M_05           ->  GetXaxis()  ->  SetTitle("N Tracks V0A");
    fQC_Event_Enum_V0M_05           ->  GetYaxis()  ->  SetTitle("N Tracks |#eta| < 0.5");
    lQCParticle    ->  Add(fQC_Event_Enum_V0M_05);
    //
    fQC_Event_Enum_V0M_08           = new TH2D("fQC_Event_Enum_V0M_08", "Events with tracks in V0M acceptance vs Mult Eta08",   5000,   0., 5000., 5000,   0., 5000.);
    fQC_Event_Enum_V0M_08           ->  GetXaxis()  ->  SetTitle("N Tracks V0A");
    fQC_Event_Enum_V0M_08           ->  GetYaxis()  ->  SetTitle("N Tracks |#eta| < 0.8");
    lQCParticle    ->  Add(fQC_Event_Enum_V0M_08);
    //
    fQC_Event_Enum_V0M_10           = new TH2D("fQC_Event_Enum_V0M_10", "Events with tracks in V0M acceptance vs Mult Eta10",   5000,   0., 5000., 5000,   0., 5000.);
    fQC_Event_Enum_V0M_10           ->  GetXaxis()  ->  SetTitle("N Tracks V0A");
    fQC_Event_Enum_V0M_10           ->  GetYaxis()  ->  SetTitle("N Tracks |#eta| < 1.0");
    lQCParticle    ->  Add(fQC_Event_Enum_V0M_10);
    PostData ( 2, lQCParticle );
    //
    //  --- Particle Tree
    OpenFile(3);
    //  --- --- Allocating tree and setting branches
    tParticle      =   new TTree   ( Form( "PhiCandidate_%s", fRunName.Data() ), "Data Tree for Phi Candidates" );
    tParticle      ->  Branch      ( "nPart",       &fTNParticle,       "fTNParticle/I");
    tParticle      ->  Branch      ( "Eta_10",      &fCurrent_E10,      "fCurrent_E10/F" );
    tParticle      ->  Branch      ( "Eta_08",      &fCurrent_E08,      "fCurrent_E08/F" );
    tParticle      ->  Branch      ( "Eta_05",      &fCurrent_E05,      "fCurrent_E05/F" );
    tParticle      ->  Branch      ( "V0A",         &fCurrent_V0A,      "fCurrent_V0A/F" );
    tParticle      ->  Branch      ( "V0M",         &fCurrent_V0M,      "fCurrent_V0M/F" );
    tParticle      ->  Branch      ( "Px",          &fTPx,              "fTPx[fTNParticle]/F" );
    tParticle      ->  Branch      ( "Py",          &fTPy,              "fTPy[fTNParticle]/F" );
    tParticle      ->  Branch      ( "Pz",          &fTPz,              "fTPz[fTNParticle]/F" );
    PostData(3, tParticle);
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::UserExec( Option_t* )                      {
    //  --- Check the Event is available and within requirements
    if ( !fIsEventCandidate() )             return;
    uIsEventMultiplicityAvailable();
    //
    for ( Int_t iTrack = 0; iTrack < fMCD->GetNumberOfTracks(); iTrack++ )    {
        //
        //  --- Recover current particle
        TParticle* fCurrent_Particle = static_cast<TParticle*>( fMCD->Particle(iTrack) );
        //
        //  Check it is a requested particle
        if ( !fCurrent_Particle )                               continue;
        if (  fCurrent_Particle->GetPdgCode() != kParticlePDG ) continue;
        //
        //  --- Save Kinematics
        fTPx[fTNParticle]  =   fCurrent_Particle->Px();
        fTPy[fTNParticle]  =   fCurrent_Particle->Py();
        fTPz[fTNParticle]  =   fCurrent_Particle->Pz();
        fTNParticle++;
    }
    // Saving output
    fFillTrees();
    fPostData();
}

//_____________________________________________________________________________

bool        AliAnalysisTask_Phi_MC::fIsEventCandidate ()                       {
    //
    //  Counting the Triggered events
    //
    fQC_Event_Enum_FLL->Fill("ALL",1);
    //
    //>-> Starting Mandatory Checks
    //>->____________________________________________________________________//
    //
    fMCD = MCEvent();
    //
    // Check the event is there
    if ( !MCEvent() )  {
        uFillEventEnumerate("fAOD-fMCD");
        fPostData();
        return false;
    }
    //
    // Recover and Check the MC tracks
    if ( !((MCEvent())->Stack()) )  {
        uFillEventEnumerate("TrackArray");
        fPostData();
        return false;
    }
    //  Event Counting
    uFillEventEnumerate("Accepted");
    //
    //  --- Zero all counters
    fTNParticle         = 0;
    fCurrent_E05        = 0;
    fCurrent_E08        = 0;
    fCurrent_E10        = 0;
    fCurrent_V0A        = 0;
    fCurrent_V0M        = 0;
    //
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTask_Phi_MC::uIsEventMultiplicityAvailable ( )           {
    //
    //  --- Loop on MC particles
    for ( Int_t iTrack = 0; iTrack < fMCD->GetNumberOfTracks(); iTrack++ )    {
        //
        //  --- Recover current particle
        TParticle* fCurrent_Particle = static_cast<TParticle*>( fMCD->Particle(iTrack) );
        //
        //  Check it is a requested particle
        if ( !fCurrent_Particle )                           continue;
        if ( !fCurrent_Particle->GetPDG())                  continue;
        if (  fCurrent_Particle->GetPDG()->Charge() == 0 )  continue;
        if ( !fMCD->IsPhysicalPrimary( iTrack ) )           continue;
        if (  fCurrent_Particle->Pt() < 0.15 )              continue;
        //
        //  --- Eta Nch Estimators
        //  --- --- #eta 0.5
        if (  fabs( fCurrent_Particle->Eta() ) < 0.5 )  fCurrent_E05++;
        //  --- --- #eta 0.8
        if (  fabs( fCurrent_Particle->Eta() ) < 0.8 )  fCurrent_E08++;
        //  --- --- #eta 1.0
        if (  fabs( fCurrent_Particle->Eta() ) < 1.0 )  fCurrent_E10++;
        //
        //  --- V0 Estimators
        //  --- --- V0A & M
        if ( (fCurrent_Particle->Eta() < +5.1) && (fCurrent_Particle->Eta() > +2.8) )   { fCurrent_V0A++; fCurrent_V0M++; }
        //  --- --- V0M
        if ( (fCurrent_Particle->Eta() < -1.7) && (fCurrent_Particle->Eta() > -3.7) )   fCurrent_V0M++;
    }
    //
    fQC_Event_Enum_E05      ->  Fill( fCurrent_E05 );
    fQC_Event_Enum_E08      ->  Fill( fCurrent_E08 );
    fQC_Event_Enum_E10      ->  Fill( fCurrent_E10 );
    fQC_Event_Enum_V0A      ->  Fill( fCurrent_V0A );
    fQC_Event_Enum_V0M      ->  Fill( fCurrent_V0M );
    fQC_Event_Enum_V0A_05   ->  Fill( fCurrent_V0A, fCurrent_E05 );
    fQC_Event_Enum_V0M_05   ->  Fill( fCurrent_V0M, fCurrent_E05 );
    fQC_Event_Enum_V0A_08   ->  Fill( fCurrent_V0A, fCurrent_E08 );
    fQC_Event_Enum_V0M_08   ->  Fill( fCurrent_V0M, fCurrent_E08 );
    fQC_Event_Enum_V0A_10   ->  Fill( fCurrent_V0A, fCurrent_E10 );
    fQC_Event_Enum_V0M_10   ->  Fill( fCurrent_V0M, fCurrent_E10 );
    //
    return true;
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::Terminate( Option_t* )                     {
}

//-----------------------------------------------------------------------------
//      UTILITY FUNCTIONS
//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::fPostData( )                               {
    // Post-data for TLists
    PostData(1, lAnalysisOutputList );
    PostData(2, lQCParticle );
    if ( fTNParticle > 0 ) PostData(3, tParticle );
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::fFillTrees( )                              {
    if ( fTNParticle > 0 ) tParticle->Fill();
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::uSetEventCountLabels ( TH1D * fEvCount )   {
    fEvCount->GetXaxis()->SetBinLabel(1,"ALL");
    fEvCount->GetXaxis()->SetBinLabel(2,"fAOD-fMCD");
    fEvCount->GetXaxis()->SetBinLabel(3,"TrackArray");
    fEvCount->GetXaxis()->SetBinLabel(4,"NoTrigger");
    fEvCount->GetXaxis()->SetBinLabel(5,"PIDResponse");
    fEvCount->GetXaxis()->SetBinLabel(6,"IncompleteDAQ");
    fEvCount->GetXaxis()->SetBinLabel(7,"NoSPDVTX");
    fEvCount->GetXaxis()->SetBinLabel(8,"TRK-SPD Mismatch");
    fEvCount->GetXaxis()->SetBinLabel(9,"VTX<Cut");
    fEvCount->GetXaxis()->SetBinLabel(10,"Accepted");
    fEvCount->GetXaxis()->SetBinLabel(11,"HasMult");
    fEvCount->GetXaxis()->SetBinLabel(12,"Pile-Up");
    fEvCount->GetXaxis()->SetBinLabel(13,"Pile-Up in Mult");
    fEvCount->GetXaxis()->SetBinLabel(14,"NoPhiCand");
    fEvCount->GetXaxis()->SetBinLabel(15,"OFPhiCand");
    fEvCount->GetXaxis()->SetBinLabel(16,"NoPhTCand");
    fEvCount->GetXaxis()->SetBinLabel(17,"OFPhTCand");
    fEvCount->GetXaxis()->SetBinLabel(18,"NoKaonCand");
    fEvCount->GetXaxis()->SetBinLabel(19,"OFKaonCand");
    fEvCount->GetXaxis()->SetBinLabel(20,"NoKaoTCand");
    fEvCount->GetXaxis()->SetBinLabel(21,"OFKaoTCand");
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::uFillEventEnumerate ( Int_t iIndex )       {
    for ( Int_t iFill = 1; iFill <= iIndex; iFill++ )   {
        fQC_Event_Enum_FLL->Fill(iFill);
    }
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::uFillEventEnumerate ( TString iBinName )   {
    for ( Int_t iBin = 1; iBin <= fQC_Event_Enum_FLL->GetNbinsX(); iBin++ )    {
        if ( strcmp(iBinName.Data(),fQC_Event_Enum_FLL->GetXaxis()->GetBinLabel(iBin)) == 0 )  {
            uFillEventEnumerate(iBin-1);
        }
    }
}

//_____________________________________________________________________________
