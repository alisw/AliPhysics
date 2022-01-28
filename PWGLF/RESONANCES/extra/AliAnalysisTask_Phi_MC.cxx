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
#include "AliAODInputHandler.h"
#include "AliMCVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODHeader.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisTask_Phi_MC.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTask_Phi_MC;

using namespace std;

ClassImp(AliAnalysisTask_Phi_MC)

            AliAnalysisTask_Phi_MC::AliAnalysisTask_Phi_MC() : AliAnalysisTaskSE(), kComputeSpherocity(0), kSpherocityPTWeight(0), kComputeRT(0), fIs_p_p(0), fIs_p_Pb(0), fIs_Pb_Pb(0), fRunName(0), fQC_Event_Enum_FLL(0)  {
            }

//_____________________________________________________________________________

            AliAnalysisTask_Phi_MC::AliAnalysisTask_Phi_MC(const char* name) : AliAnalysisTaskSE(name), kComputeSpherocity(0), kSpherocityPTWeight(0), kComputeRT(0), fIs_p_p(0), fIs_p_Pb(0), fIs_Pb_Pb(0), fRunName(0), fQC_Event_Enum_FLL(0)  {
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
                if  ( lQCParticle_0333 )    delete lQCParticle_0333;
                if  ( lAnalysisOutputList ) delete lAnalysisOutputList;
                // Deleting TTrees
                if  ( tParticle_0333 )      delete tParticle_0333;
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
    lQCParticle_0333    =   new TList();
    lQCParticle_0333    ->  SetOwner( kTRUE );
    //
    fQC_Event_Enum_FLL              = new TH1D("fQC_Event_Enum_FLL",       "Event Selection",                                  29, -0.5, 28.5);
    fQC_Event_Enum_FLL              ->  GetXaxis()  ->  SetTitle("");
    fQC_Event_Enum_FLL              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    uSetEventCountLabels(fQC_Event_Enum_FLL);
    lQCParticle_0333    ->  Add(fQC_Event_Enum_FLL);
    //
    PostData ( 2, lQCParticle_0333 );
    //
    //  --- Particle Tree
    OpenFile(3);
    //  --- --- Allocating tree and setting branches
    tParticle_0333      =   new TTree   ( Form( "PhiCandidate_%s", fRunName.Data() ), "Data Tree for Phi Candidates" );
    tParticle_0333      ->  Branch      ( "n0333",      &fTNParticle_0333,      "fTNParticle_0333/I");
    tParticle_0333      ->  Branch      ( "Px",         &fTPx,                  "fTPx[n0333]/F" );
    tParticle_0333      ->  Branch      ( "Py",         &fTPy,                  "fTPy[n0333]/F" );
    tParticle_0333      ->  Branch      ( "Pz",         &fTPz,                  "fTPz[n0333]/F" );
    PostData(3, tParticle_0333);
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::UserExec( Option_t* )                      {
    //  --- Check the Event is available and within requirements
    if ( !fIsEventCandidate() )    return;
    //
    Int_t   fNTracks    =   AODMCTrackArray->GetEntriesFast();
    for ( Int_t iTrack = 0; iTrack < fNTracks; iTrack++ )    {
        //
        //  --- Recover current particle
        AliAODMCParticle* fCurrent_Particle = static_cast<AliAODMCParticle*>( AODMCTrackArray->At(iTrack) );
        //
        //  Check it is a requested particle
        if ( !fCurrent_Particle )                       continue;
        if (  fCurrent_Particle->GetPdgCode() == 333 )  continue;
        //
        //  --- Save Kinematics
        fTPx[fTNParticle_0333]  =   fCurrent_Particle->Px();
        fTPy[fTNParticle_0333]  =   fCurrent_Particle->Py();
        fTPz[fTNParticle_0333]  =   fCurrent_Particle->Pz();
        fTNParticle_0333++;
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
    // Recovering MC Event
    fMCD = dynamic_cast<AliMCEvent*>(MCEvent());
    //
    // Check the event is there
    if ( !fMCD )  {
        uFillEventEnumerate("fAOD-fMCD");
        fPostData();
        return false;
    }
    //
    // Recover and Check the MC tracks
    AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !AODMCTrackArray )  {
        uFillEventEnumerate("TrackArray");
        fPostData();
        return false;
    }
    //  Event Counting
    uFillEventEnumerate("Accepted");
    //
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTask_Phi_MC::uIsEventMultiplicityAvailable ()           {
    return true;
    /*
    // Recovering Multiplicity information
    AliMultSelection   *fMultSelectio2  = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if ( !fMultSelection )  fCurrent_V0M    =   102.;
    // if pPb should be V0A
                    fCurrent_V0M    =   fMultSelectio2->GetMultiplicityPercentile("V0M");
    if ( fIs_p_Pb ) fCurrent_V0M    =   fMultSelectio2->GetMultiplicityPercentile("V0A");
    fCurrent_TRK    =   (dynamic_cast<AliAODHeader*>(fAOD->GetHeader()))->GetRefMultiplicityComb08();
    //
    //  INELGT0 Check
    if ( fCurrent_TRK != -1 && fCurrent_TRK != -2 ) {
        for ( Int_t iTRK = 0; iTRK < fAOD->GetMultiplicity()->GetNumberOfTracklets(); iTRK++ )  {
            if ( TMath::Abs(fAOD->GetMultiplicity()->GetEta(iTRK)) < 1.0 )  {
                fSetEventMask(0);
                break;
            }
        }
    }
    //
    if ( !fCheckMask(0) )       fCurrent_V0M = 104;
    //
    if ( fCurrent_V0M == -200 ) fCurrent_V0M = 154;
    if ( fCurrent_V0M == -201 ) fCurrent_V0M = 156;
    if ( fCurrent_V0M == -202 ) fCurrent_V0M = 158;
    if ( fCurrent_V0M == -203 ) fCurrent_V0M = 160;
    if ( fCurrent_V0M == -204 ) fCurrent_V0M = 162;
    if ( ( fCurrent_V0M >= 0 ) && ( fCurrent_V0M <= 100 ) )   {
        if ( fAOD->IsPileupFromSPD() )              fCurrent_V0M += 300;
        else if ( fAOD->IsPileupFromSPDInMultBins() )    fCurrent_V0M += 400;
    }
    //
    // Fill the QC on Multiplicity
    fQC_Event_Enum_V0M->Fill(fCurrent_V0M);
    fQC_Event_Enum_TRK->Fill(fCurrent_TRK);
    fQC_Event_Enum_V0T->Fill(fCurrent_V0M,fCurrent_TRK);
    //
    return true;
     */
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
    PostData(2, lQCParticle_0333 );
    PostData(3, tParticle_0333 );
}

//_____________________________________________________________________________

void        AliAnalysisTask_Phi_MC::fFillTrees( )                              {
    if ( fTNParticle_0333 > 0 ) tParticle_0333->Fill();
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
