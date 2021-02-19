// !TODO LIST
// TODO: You're all set!

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TTree.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODHeader.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskPhiCount.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskPhiCount;
class AliPIDResponse;
 
using namespace std;

ClassImp(AliAnalysisTaskPhiCount)

            AliAnalysisTaskPhiCount::AliAnalysisTaskPhiCount() : AliAnalysisTaskSE(), kMCbool(0), kPhibool(0), kKaonbool(0), kFilterBit(5), fAOD(0), fESD(0), fMCD(0), AODMCTrackArray(0), fPIDResponse(0), fPrimaryVertex(0), fCurrent_Track(0), fCurrent_Track_Charge(0), fCurrent_Track_Momentum(0),fCurrent_Track_TransMom(0), fCurrent_Track_Eta(0), fCurrent_Track_Phi(0), fCurrent_Track_DCAXY(0), fCurrent_Track_DCAZ(0), fIsTPCAvailable(0), fIsTOFAvailable(0), fBetaFromTOFSignal(0), fTPCSignal(0), fQCOutputList(0), fQC_Event_Enumerate(0), fQC_Event_Vertex_Fll(0), fQC_Event_Vertex_Cut(0), fQC_Event_Enum_Mult(0), fQC_Tracks_P_Momentum(0), fQC_Tracks_P_TransMom(0), fQC_Tracks_P_Eta(0), fQC_Tracks_P_Phi(0), fQC_Tracks_M_Momentum(0), fQC_Tracks_M_TransMom(0), fQC_Tracks_M_Eta(0), fQC_Tracks_M_Phi(0), fQC_Tracks_DCAXY_P(0), fQC_Tracks_DCAZ_P(0), fQC_Tracks_DCAXY_PT(0), fQC_Tracks_DCAZ_PT(0), fQC_Kaons_P_Momentum(0), fQC_Kaons_P_TransMom(0), fQC_Kaons_P_Eta(0), fQC_Kaons_P_Phi(0),fQC_Kaons_M_Momentum(0), fQC_Kaons_M_TransMom(0), fQC_Kaons_M_Eta(0), fQC_Kaons_M_Phi(0), fQC_Kaons_DCAXY_P(0), fQC_Kaons_DCAZ_P(0), fQC_Kaons_DCAXY_PT(0), fQC_Kaons_DCAZ_PT(0), fQC_PID_SignalTPC_P(0), fQC_PID_SignalTOF_P(0), fQC_PID_SignalTPC_PT(0), fQC_PID_SignalTOF_PT(0), fQC_Kaons_SigmaTPC_P(0), fQC_Kaons_SigmaTOF_P(0), fQC_Kaons_SigmaTPC_PT(0), fQC_Kaons_SigmaTOF_PT(0), fQC_Kaons_SigmaTOF_TPC(0), fQC_Kaons_SignalTPC_P(0), fQC_Kaons_SignalTOF_P(0), fQC_Kaons_SignalTPC_PT(0), fQC_Kaons_SignalTOF_PT(0), fMultiplicity(0), fPhiCandidate(0), fnPhi(0), fKaonCandidate(0), fnKaon(0), fPhiEfficiency(0), fnPhiTru(0), fKaonEfficiency(0), fnKaonTru(0), fAnalysisOutputList(0)   {
    
}

//_____________________________________________________________________________

            AliAnalysisTaskPhiCount::AliAnalysisTaskPhiCount(const char* name) : AliAnalysisTaskSE(name), kMCbool(0), kPhibool(0), kKaonbool(0), kFilterBit(5), fAOD(0), fESD(0), fMCD(0), AODMCTrackArray(0), fPIDResponse(0), fPrimaryVertex(0), fCurrent_Track(0), fCurrent_Track_Charge(0), fCurrent_Track_Momentum(0),fCurrent_Track_TransMom(0), fCurrent_Track_Eta(0), fCurrent_Track_Phi(0), fCurrent_Track_DCAXY(0), fCurrent_Track_DCAZ(0), fIsTPCAvailable(0), fIsTOFAvailable(0), fBetaFromTOFSignal(0), fTPCSignal(0), fQCOutputList(0), fQC_Event_Enumerate(0), fQC_Event_Vertex_Fll(0), fQC_Event_Vertex_Cut(0), fQC_Event_Enum_Mult(0), fQC_Tracks_P_Momentum(0), fQC_Tracks_P_TransMom(0), fQC_Tracks_P_Eta(0), fQC_Tracks_P_Phi(0), fQC_Tracks_M_Momentum(0), fQC_Tracks_M_TransMom(0), fQC_Tracks_M_Eta(0), fQC_Tracks_M_Phi(0), fQC_Tracks_DCAXY_P(0), fQC_Tracks_DCAZ_P(0), fQC_Tracks_DCAXY_PT(0), fQC_Tracks_DCAZ_PT(0), fQC_Kaons_P_Momentum(0), fQC_Kaons_P_TransMom(0), fQC_Kaons_P_Eta(0), fQC_Kaons_P_Phi(0),fQC_Kaons_M_Momentum(0), fQC_Kaons_M_TransMom(0), fQC_Kaons_M_Eta(0), fQC_Kaons_M_Phi(0), fQC_Kaons_DCAXY_P(0), fQC_Kaons_DCAZ_P(0), fQC_Kaons_DCAXY_PT(0), fQC_Kaons_DCAZ_PT(0), fQC_PID_SignalTPC_P(0), fQC_PID_SignalTOF_P(0), fQC_PID_SignalTPC_PT(0), fQC_PID_SignalTOF_PT(0), fQC_Kaons_SigmaTPC_P(0), fQC_Kaons_SigmaTOF_P(0), fQC_Kaons_SigmaTPC_PT(0), fQC_Kaons_SigmaTOF_PT(0), fQC_Kaons_SigmaTOF_TPC(0), fQC_Kaons_SignalTPC_P(0), fQC_Kaons_SignalTOF_P(0), fQC_Kaons_SignalTPC_PT(0), fQC_Kaons_SignalTOF_PT(0), fMultiplicity(0), fPhiCandidate(0), fnPhi(0), fKaonCandidate(0), fnKaon(0), fPhiEfficiency(0), fnPhiTru(0), fKaonEfficiency(0), fnKaonTru(0), fAnalysisOutputList(0)    {
    // Define Input
    DefineInput(0, TChain::Class());
    
    // Define Output
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TTree::Class());
    DefineOutput(4, TTree::Class());
    DefineOutput(5, TTree::Class());
    DefineOutput(6, TTree::Class());
}

//_____________________________________________________________________________

            AliAnalysisTaskPhiCount::~AliAnalysisTaskPhiCount()                 {
    // Deleting TLists
    if( fAnalysisOutputList )
    {
        delete fAnalysisOutputList;
    }
    if( fQCOutputList )
    {
        delete fQCOutputList;
    }
    if( fKaonCandidate )
    {
        delete fKaonCandidate;
    }
    if( fPhiCandidate )
    {
        delete fPhiCandidate;
    }
    if( fKaonEfficiency )
    {
        delete fKaonEfficiency;
    }
    if( fPhiEfficiency )
    {
        delete fPhiEfficiency;
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::UserCreateOutputObjects()                  {
    // Various utility Histograms TList initialisation
    fAnalysisOutputList     = new TList();
    fAnalysisOutputList     ->SetOwner(kTRUE);
    PostData(1, fAnalysisOutputList);
    
    // QC utility Histograms TList initialisation
    fQCOutputList   = new TList();
    fQCOutputList   ->SetOwner(kTRUE);
    
    //>->   Event
    //
    // TODO: Acc Tracks in mult bins
    //
    fQC_Event_Enumerate     = new TH1D("fQC_Event_Enumerate",   "Event Selection",                  17, 0.5, 17.5);
    fQC_Event_Vertex_Fll    = new TH1F("fQC_Event_Vertex_Fll",  "Collision Vertex (FULL)",          300, -15, 15);
    fQC_Event_Vertex_Cut    = new TH1F("fQC_Event_Vertex_Cut",  "Collision Vertex (CUTS)",          300, -15, 15);
    fQC_Event_Enum_Mult     = new TH1F("fQC_Event_Enum_Mult",   "Collision Vertex (CUTS)",          202, -1., 201.);
    fQCOutputList->Add(fQC_Event_Enumerate);
    fQCOutputList->Add(fQC_Event_Vertex_Fll);
    fQCOutputList->Add(fQC_Event_Vertex_Cut);
    fQCOutputList->Add(fQC_Event_Enum_Mult);
    
    //>->   Tracks
    fQC_Tracks_P_Momentum   = new TH1F("fQC_Tracks_P_Momentum",   "Acc. Tracks Momentum",               100, 0., 10.);
    fQC_Tracks_P_TransMom   = new TH1F("fQC_Tracks_P_TransMom",   "Acc. Tracks Trasnverse Momentum",    100, 0., 10.);
    fQC_Tracks_P_Eta        = new TH1F("fQC_Tracks_P_Eta",        "Acc. Tracks Eta",                    100, -1., 1.);
    fQC_Tracks_P_Phi        = new TH1F("fQC_Tracks_P_Phi",        "Acc. Tracks Phi",                    370, -5., 365.);
    fQC_Tracks_M_Momentum   = new TH1F("fQC_Tracks_M_Momentum",   "Acc. Tracks Momentum",               100, 0., 10.);
    fQC_Tracks_M_TransMom   = new TH1F("fQC_Tracks_M_TransMom",   "Acc. Tracks Trasnverse Momentum",    100, 0., 10.);
    fQC_Tracks_M_Eta        = new TH1F("fQC_Tracks_M_Eta",        "Acc. Tracks Eta",                    100, -1., 1.);
    fQC_Tracks_M_Phi        = new TH1F("fQC_Tracks_M_Phi",        "Acc. Tracks Phi",                    370, -5., 365.);
    fQC_Tracks_DCAXY_P      = new TH2F("fQC_Tracks_DCAXY_P",    "Acc. Tracks DCAXY",                    100, 0., 10., 400, -2., 2.);
    fQC_Tracks_DCAZ_P       = new TH2F("fQC_Tracks_DCAZ_P",     "Acc. Tracks DCAZ",                     100, 0., 10., 400, -2., 2.);
    fQC_Tracks_DCAXY_PT     = new TH2F("fQC_Tracks_DCAXY_PT",   "Acc. Tracks DCAXY",                    100, 0., 10., 400, -2., 2.);
    fQC_Tracks_DCAZ_PT      = new TH2F("fQC_Tracks_DCAZ_PT",    "Acc. Tracks DCAZ",                     100, 0., 10., 400, -2., 2.);
    fQCOutputList->Add(fQC_Tracks_P_Momentum);
    fQCOutputList->Add(fQC_Tracks_P_TransMom);
    fQCOutputList->Add(fQC_Tracks_P_Eta);
    fQCOutputList->Add(fQC_Tracks_P_Phi);
    fQCOutputList->Add(fQC_Tracks_M_Momentum);
    fQCOutputList->Add(fQC_Tracks_M_TransMom);
    fQCOutputList->Add(fQC_Tracks_M_Eta);
    fQCOutputList->Add(fQC_Tracks_M_Phi);
    fQCOutputList->Add(fQC_Tracks_DCAXY_P);
    fQCOutputList->Add(fQC_Tracks_DCAZ_P);
    fQCOutputList->Add(fQC_Tracks_DCAXY_PT);
    fQCOutputList->Add(fQC_Tracks_DCAZ_PT);
    
    //>->->     Kaons
    fQC_Kaons_P_Momentum    = new TH1F("fQC_Kaons_P_Momentum",   "Acc. Kaons Momentum",                 100, 0., 10.);
    fQC_Kaons_P_TransMom    = new TH1F("fQC_Kaons_P_TransMom",   "Acc. Kaons Trasnverse Momentum",      100, 0., 10.);
    fQC_Kaons_P_Eta         = new TH1F("fQC_Kaons_P_Eta",        "Acc. Kaons Eta",                      100, -1., 1.);
    fQC_Kaons_P_Phi         = new TH1F("fQC_Kaons_P_Phi",        "Acc. Kaons Phi",                      370, -5., 365.);
    fQC_Kaons_M_Momentum    = new TH1F("fQC_Kaons_M_Momentum",   "Acc. Kaons Momentum",                 100, 0., 10.);
    fQC_Kaons_M_TransMom    = new TH1F("fQC_Kaons_M_TransMom",   "Acc. Kaons Trasnverse Momentum",      100, 0., 10.);
    fQC_Kaons_M_Eta         = new TH1F("fQC_Kaons_M_Eta",        "Acc. Kaons Eta",                      100, -1., 1.);
    fQC_Kaons_M_Phi         = new TH1F("fQC_Kaons_M_Phi",        "Acc. Kaons Phi",                      370, -5., 365.);
    fQC_Kaons_DCAXY_P       = new TH2F("fQC_Kaons_DCAXY_P",    "Acc. Kaons DCAXY",                      100, 0., 10., 400, -2., 2.);
    fQC_Kaons_DCAZ_P        = new TH2F("fQC_Kaons_DCAZ_P",     "Acc. Kaons DCAZ",                       100, 0., 10., 400, -2., 2.);
    fQC_Kaons_DCAXY_PT      = new TH2F("fQC_Kaons_DCAXY_PT",   "Acc. Kaons DCAXY",                      100, 0., 10., 400, -2., 2.);
    fQC_Kaons_DCAZ_PT       = new TH2F("fQC_Kaons_DCAZ_PT",    "Acc. Kaons DCAZ",                       100, 0., 10., 400, -2., 2.);
    fQCOutputList->Add(fQC_Kaons_P_Momentum);
    fQCOutputList->Add(fQC_Kaons_P_TransMom);
    fQCOutputList->Add(fQC_Kaons_P_Eta);
    fQCOutputList->Add(fQC_Kaons_P_Phi);
    fQCOutputList->Add(fQC_Kaons_M_Momentum);
    fQCOutputList->Add(fQC_Kaons_M_TransMom);
    fQCOutputList->Add(fQC_Kaons_M_Eta);
    fQCOutputList->Add(fQC_Kaons_M_Phi);
    fQCOutputList->Add(fQC_Kaons_DCAXY_P);
    fQCOutputList->Add(fQC_Kaons_DCAZ_P);
    fQCOutputList->Add(fQC_Kaons_DCAXY_PT);
    fQCOutputList->Add(fQC_Kaons_DCAZ_PT);
    
    //>->   PID
    fQC_PID_SignalTPC_P     = new TH2F("fQC_PID_SignalTPC_P",  "TPC Response",                          400, 0.15, 10., 100, 0., 500.);
    fQC_PID_SignalTOF_P     = new TH2F("fQC_PID_SignalTOF_P",  "TOF Response",                          400, 0.15, 10., 120, 0., 1.2);
    fQC_PID_SignalTPC_PT    = new TH2F("fQC_PID_SignalTPC_PT", "TPC Response",                          400, 0.15, 10., 100, 0., 500.);
    fQC_PID_SignalTOF_PT    = new TH2F("fQC_PID_SignalTOF_PT", "TOF Response",                          400, 0.15, 10., 120, 0., 1.2);
    fQCOutputList->Add(fQC_PID_SignalTPC_P);
    fQCOutputList->Add(fQC_PID_SignalTPC_PT);
    fQCOutputList->Add(fQC_PID_SignalTOF_P);
    fQCOutputList->Add(fQC_PID_SignalTOF_PT);
    
    //>->->     Kaons
    fQC_Kaons_SigmaTPC_P    = new TH2F("fQC_Kaons_SigmaTPC_P",  "#sigma_{TPC}(K) Response",             100, 0., 10., 100, -10, 10);
    fQC_Kaons_SigmaTOF_P    = new TH2F("fQC_Kaons_SigmaTOF_P",  "#sigma_{TOF}(K) Response",             100, 0., 10., 100, -10, 10);
    fQC_Kaons_SigmaTPC_PT   = new TH2F("fQC_Kaons_SigmaTPC_PT", "#sigma_{TPC}(K) Response",             100, 0., 10., 100, -10, 10);
    fQC_Kaons_SigmaTOF_PT   = new TH2F("fQC_Kaons_SigmaTOF_PT", "#sigma_{TOF}(K) Response",             100, 0., 10., 100, -10, 10);
    fQC_Kaons_SigmaTOF_TPC  = new TH2F("fQC_Kaons_SigmaTOF_TPC","#sigma_{TPC}(K) vs #sigma_{TOF}(K) Response",  100, -10., 10., 100, -10, 10);
    fQC_Kaons_SignalTPC_P   = new TH2F("fQC_Kaons_SignalTPC_P",  "TPC Response",                        400, 0.15, 10., 100, 0., 500.);
    fQC_Kaons_SignalTOF_P   = new TH2F("fQC_Kaons_SignalTOF_P",  "TOF Response",                        400, 0.15, 10., 120, 0., 1.2);
    fQC_Kaons_SignalTPC_PT  = new TH2F("fQC_Kaons_SignalTPC_PT", "TPC Response",                        400, 0.15, 10., 100, 0., 500.);
    fQC_Kaons_SignalTOF_PT  = new TH2F("fQC_Kaons_SignalTOF_PT", "TOF Response",                        400, 0.15, 10., 120, 0., 1.2);
    fQCOutputList->Add(fQC_PID_SignalTPC_P);
    fQCOutputList->Add(fQC_PID_SignalTPC_PT);
    fQCOutputList->Add(fQC_PID_SignalTOF_P);
    fQCOutputList->Add(fQC_PID_SignalTOF_PT);
    fQCOutputList->Add(fQC_Kaons_SigmaTPC_P);
    fQCOutputList->Add(fQC_Kaons_SigmaTOF_P);
    fQCOutputList->Add(fQC_Kaons_SigmaTPC_PT);
    fQCOutputList->Add(fQC_Kaons_SigmaTOF_PT);
    fQCOutputList->Add(fQC_Kaons_SigmaTOF_TPC);
    
    PostData(2, fQCOutputList);
    
    // Where i is the slot that writes to a tree
    OpenFile(3);
    // PhiCandidate Tree Set-Up
    fPhiCandidate = new TTree   ("PhiCandidate",    "Data Tree for Phi Candidates");
    fPhiCandidate->Branch       ("Multiplicity",    &fMultiplicity,     "fMultiplicity/b");
    fPhiCandidate->Branch       ("nPhi",            &fnPhi,             "fnPhi/b");
    fPhiCandidate->Branch       ("Px",              &fPhiPx,            "fPhiPx[fnPhi]/F");
    fPhiCandidate->Branch       ("Py",              &fPhiPy,            "fPhiPy[fnPhi]/F");
    fPhiCandidate->Branch       ("Pz",              &fPhiPz,            "fPhiPz[fnPhi]/F");
    fPhiCandidate->Branch       ("InvMass",         &fInvMass,          "fInvMass[fnPhi]/F");
    fPhiCandidate->Branch       ("iKaon",           &fiKaon,            "fiKaon[fnPhi]/b");
    fPhiCandidate->Branch       ("jKaon",           &fjKaon,            "fjKaon[fnPhi]/b");
    if ( kMCbool )  {
        fPhiCandidate->Branch   ("Nature",          &fNature,           "fNature[fnPhi]/b");
    }
    
    if ( kPhibool )                 PostData(3, fPhiCandidate);
    
    // KaonCandidate Tree Set-Up
    fKaonCandidate = new TTree ("KaonCandidate",    "Data Tree for Kaon Candidates");
    fKaonCandidate->Branch     ("Multiplicity",     &fMultiplicity,     "fMultiplicity/b");
    fKaonCandidate->Branch     ("nKaon",            &fnKaon,            "fnKaon/b");
    fKaonCandidate->Branch     ("Px",               &fKaonPx,           "fKaonPx[fnKaon]/F");
    fKaonCandidate->Branch     ("Py",               &fKaonPy,           "fKaonPy[fnKaon]/F");
    fKaonCandidate->Branch     ("Pz",               &fKaonPz,           "fKaonPz[fnKaon]/F");
    fKaonCandidate->Branch     ("Charge",           &fCharge,           "fCharge[fnKaon]/B");
    fKaonCandidate->Branch     ("TOFSigma",         &fTOFSigma,         "fTOFSigma[fnKaon]/B");
    fKaonCandidate->Branch     ("TPCSigma",         &fTPCSigma,         "fTPCSigma[fnKaon]/B");
    
    if ( kKaonbool )                PostData(4, fKaonCandidate);

    fPhiEfficiency = new TTree  ("PhiEfficiency",   "MC Tree for Phi Efficiency");
    fPhiEfficiency->Branch      ("Multiplicity",    &fMultiplicity,     "fMultiplicity/b");
    fPhiEfficiency->Branch      ("nPhi",            &fnPhiTru,          "fnPhiTru/b");
    fPhiEfficiency->Branch      ("Px",              &fPhiTruPx,         "fPhiTruPx[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Py",              &fPhiTruPy,         "fPhiTruPy[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Pz",              &fPhiTruPz,         "fPhiTruPz[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Selection",       &fSelection,        "fSelection[fnPhiTru]/b");
    
    if ( kPhibool   &&  kMCbool )   PostData(5, fPhiEfficiency);
    
    fKaonEfficiency = new TTree ("KaonEfficiency",  "MC Tree for Kaon Efficiency");
    fKaonEfficiency->Branch     ("Multiplicity",    &fMultiplicity,     "fMultiplicity/b");
    
    if ( kKaonbool  &&  kMCbool )   PostData(6, fKaonEfficiency);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::UserExec( Option_t* )                      {
    if ( !fIsEventCandidate() )    return;
    
    // Setting utility variables
    Int_t           nTrack(fAOD->GetNumberOfTracks());
    TLorentzVector  fKaon1, fKaon2, fPhi;
    
    // Setting zero all counters and global variables
    fnPhi           =   0;
    fnPhiTru        =   0;
    fnPhiRec        =   0;
    fnKaon          =   0;
    fnPhiTru        =   0;
    
    //fSetZero();

    // Looping over tracks
    for ( Int_t iTrack(0); iTrack < nTrack; iTrack++ )
    {
        // Recovering Track
        fCurrent_Track  =   static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        
        // Check the track has due requirements
        if ( !fIsTrackCandidate() ) continue;
        
        // Fill the QC for Event Tracks
        fQC_TRK             ();
        
        // Fill the QC for PID
        fQC_PID             ();
        
        // Check the PID is present and within requirements
        if ( !fIsKaonCandidate() ) continue;
        
        // Fill the QC for Kaons Event Tracks
        fQC_TRK_Kaons       ();
        
        // Fill the QC for Kaons PID
        fQC_PID_Sel_Kaons   ();
        
        // Filling the Kaon Tree
        fKaonPx[fnKaon] =   fCurrent_Track->Px();
        fKaonPy[fnKaon] =   fCurrent_Track->Py();
        fKaonPz[fnKaon] =   fCurrent_Track->Pz();
        fCharge[fnKaon] =   static_cast<Char_t>(fCurrent_Track->GetSign());
        fKaonLabels[fnKaon] =   fCurrent_Track->GetLabel();
        fnKaon++;
    }
    
    //Coupling Kaons
    for ( Int_t iKaon(0); iKaon < fnKaon; iKaon++)
    {
        // Storing first Kaon kinematics and sign
        fKaon1.SetXYZM(fKaonPx[iKaon],fKaonPy[iKaon],fKaonPz[iKaon],.493677);
        
        for ( Int_t jKaon(iKaon+1); jKaon < fnKaon; jKaon++)
        {
            if ( fCharge[iKaon] == fCharge[jKaon] ) continue;
            
            // Storing second Kaon Kinematics and combining the two
            fKaon2.SetXYZM(fKaonPx[jKaon],fKaonPy[jKaon],fKaonPz[jKaon],.493677);
            
            // Check the Phi is a good candidate
            fPhi = (fKaon1 + fKaon2);
            if ( !fIsPhiCandidate(fPhi) ) continue;
            
            fPhiPx[fnPhi]       =   fPhi.Px();
            fPhiPy[fnPhi]       =   fPhi.Py();
            fPhiPz[fnPhi]       =   fPhi.Pz();
            fInvMass[fnPhi]     =   (fPhi).Mag();
            fiKaon[fnPhi]       =   iKaon;
            fjKaon[fnPhi]       =   jKaon;
            fNature[fnPhi]      =   0;
            if ( kMCbool && fIsCandidateTruPhi(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(fKaonLabels[iKaon])),static_cast<AliAODMCParticle*>(AODMCTrackArray->At(fKaonLabels[jKaon]))) )
            {
                fNature[fnPhi]          =   1;
                fPhiRecParticles[fnPhiRec] =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(fKaonLabels[iKaon]))->GetMother()));
                fnPhiRec++;
            }
            fnPhi++;
        }
    }
    
    // Loop over all primary MC particle
    if ( kMCbool )
    {
        Int_t           nTrack(AODMCTrackArray->GetEntriesFast());
        for ( Int_t iTrack(0); iTrack < nTrack; iTrack++ )
        {
            AliAODMCParticle* fPhiTru = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(iTrack));
            
            if ( !fIsPhi( fPhiTru ) ) continue;
            
            // Kinematics
            fPhiTruPx[fnPhiTru]        =   fPhiTru->Px();
            fPhiTruPy[fnPhiTru]        =   fPhiTru->Py();
            fPhiTruPz[fnPhiTru]        =   fPhiTru->Pz();
            
            // Options
            fSelection[fnPhiTru]       =   0;
            if( fIsPhiGen(fPhiTru) )    fSelection[fnPhiTru] = 1;
            if( fIsPhiRec(fPhiTru) )    fSelection[fnPhiTru] = 2;
            
            fnPhiTru++;
        }
    }
    
    // Saving output
    fFillTrees();
    fPostData();
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::Terminate( Option_t* )                     {

}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fPostData( )                               {
    // Post-data for TLists
    PostData(1, fAnalysisOutputList);
    PostData(2, fQCOutputList);
    
    // Post-data for TTrees
    if ( kPhibool )                 PostData(3, fPhiCandidate);
    if ( kKaonbool )                PostData(4, fKaonCandidate);
    if ( kPhibool   &&  kMCbool )   PostData(5, fPhiEfficiency);
    if ( kKaonbool  &&  kMCbool )   PostData(6, fKaonEfficiency);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fFillTrees( )                              {
    if ( fnPhi      != 0 )  {
        fPhiCandidate   ->Fill();
    }
    if ( fnPhi      == 0 )  {
        fQC_Event_Enumerate->Fill(10);
    }
    if ( fnPhi      >= 153 ){
        fQC_Event_Enumerate->Fill(11);
    }
    if ( fnPhiTru   != 0 )  {
        fPhiEfficiency  ->Fill();
    }
    if ( fnPhiTru   == 0 )  {
        fQC_Event_Enumerate->Fill(12);
    }
    if ( fnPhiTru   >= 153 ){
        fQC_Event_Enumerate->Fill(13);
    }
    if ( fnKaon     != 0 )  {
        fKaonCandidate   ->Fill();
    }
    if ( fnKaon     == 0 )  {
        fQC_Event_Enumerate->Fill(14);
    }
    if ( fnKaon     >= 153 ){
        fQC_Event_Enumerate->Fill(15);
    }
    if ( fnKaonTru  != 0 )  {
        fKaonEfficiency   ->Fill();
    }
    if ( fnKaonTru  == 0 )  {
        fQC_Event_Enumerate->Fill(16);
    }
    if ( fnKaonTru  >= 153 ){
        fQC_Event_Enumerate->Fill(17);
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fSetZero( )                                {
    fnPhi = 0;
    fnPhiTru = 0;
    fnKaon = 0;
    fnKaonTru = 0;
}

//_____________________________________________________________________________

Double_t   *AliAnalysisTaskPhiCount::fGetDCA( )                                 {
  
    // Create result
    Double_t   *fDCACoordinates     =   new Double_t[2];
    fDCACoordinates[0]              =   -999.;
    fDCACoordinates[1]              =   -999.;
    Double_t   *fDCACovariantMx     =   new Double_t[3];
    
    // Create utility
    AliExternalTrackParam fExternalTrack;
    fExternalTrack.CopyFromVTrack(fCurrent_Track);
  
    // Original comment:
    // * Taken from https://github.com/alisw/AliPhysics/blob/master/PWGLF/RESONANCES/extra/AliAnalysisTaskLambdaStar.cxx
    // // Propagation through the beam pipe would need a correction
    // // for material, I guess.
    //
    if  ( fExternalTrack.GetX() > 3.0 ) {
        
        return nullptr;
    }
    
    Bool_t  fSuccess    =   fExternalTrack.PropagateToDCA(fPrimaryVertex,fAOD->GetMagneticField(),10.,fDCACoordinates,fDCACovariantMx);
    
    if  ( !fSuccess ) return nullptr;
    
    return fDCACoordinates;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsEventCandidate ()                       {
    // Recovering Event Data
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fMCD = dynamic_cast<AliMCEvent*>(MCEvent());
    
    // Check the event is there
    if ( !fAOD || (!fMCD && kMCbool) )
    {
        fFillEventEnumerate(1);
        fPostData();
        return false;
    }
    
    // Recover the MC tracks
    if ( kMCbool )      AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !AODMCTrackArray && kMCbool )
    {
        fFillEventEnumerate(2);
        fPostData();
        return false;
    }
    
    // Define and Fetch PID with Manager
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man)
    {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
    }
    
    // Check the PID is present
    if ( !fPIDResponse )
    {
        fFillEventEnumerate(3);
        fPostData();
        return false;
    }
    
    // Check the Primary Vertex is there and is correct
    if ( !fIsPrimaryVertexCandidate() )   return false;
    
    // Check the Multiplicity is there and is correct
    if ( !fIsEventMultiplOK() )   return false;
    
    // Check the event is Pile-up from the SPD
    if ( fAOD->IsPileupFromSPDInMultBins() )
    {
        fFillEventEnumerate(8);
        fPostData();
        return false;
    }
    
    fFillEventEnumerate(9);
    
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsEventMultiplOK ()                       {
    // Recovering Multiplicity information
    AliMultSelection   *fMultSelection  = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    
    if ( !fMultSelection )    {
        fFillEventEnumerate(7);
        fPostData();
        return false;
    }
        
    fMultiplicity   =   static_cast<UChar_t>(fMultSelection->GetMultiplicityPercentile("V0M"));
    
    // Fill the QC on Multiplicity
    fQC_Event_Enum_Mult->Fill(fMultiplicity);
    
    return true;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fFillEventEnumerate ( Int_t iIndex )       {
    for ( Int_t iFill = 0; iFill <= iIndex; iFill++ )
    {
        fQC_Event_Enumerate->Fill(iFill);
    }
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPrimaryVertexCandidate ( )              {
    // Recovering Primary Vertex from General methods and SPD
    auto    PrimaryVertexSPD    = fAOD->GetPrimaryVertexSPD();
    auto    PrimaryVertexTRK    = fAOD->GetPrimaryVertex();
            fPrimaryVertex      = PrimaryVertexTRK;
    
    // Requires the vertex is reconstructed by the SPD
    if ( !PrimaryVertexSPD  ||  PrimaryVertexSPD->GetNContributors() < 1 )
    {
        fFillEventEnumerate(4);
        fPostData();
        return false;
    }

    // In lack of the general Method reconstructed Vertex, take the SPD reconstruction
    if ( !PrimaryVertexTRK  ||  PrimaryVertexTRK->GetNContributors() < 1 )
    {
        fPrimaryVertex = PrimaryVertexSPD;
    }
    
    // If Both are present and reconstructed discard event if the two mismatch for more than 0.5cm
    else
    {
        auto VertexZSPD = PrimaryVertexSPD->GetZ();
        auto VertexZTRK = PrimaryVertexTRK->GetZ();
        if ( std::fabs(VertexZSPD-VertexZTRK) > 0.5 )
        {
            fFillEventEnumerate(5);
            fPostData();
            return false;
        }
    }
    
    // Fill the Vertex Z position histogram
    fQC_Event_Vertex_Fll->Fill(fPrimaryVertex->GetZ());
    
    if ( std::fabs(fPrimaryVertex->GetZ()) > 10. )
    {
        fFillEventEnumerate(6);
        fPostData();
        return false;
    }
    
    // Fill the Vertex Z position histogram
    fQC_Event_Vertex_Cut->Fill(fPrimaryVertex->GetZ());
    
    fPostData();
    return  true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsTrackCandidate ( )                      {
    // Check the track is there and has proper kinematics
    if ( !fCurrent_Track    || !fCurrent_Track->TestFilterBit(kFilterBit) )                              return false;
    fCurrent_Track_TransMom = fCurrent_Track->Pt();
    fCurrent_Track_Eta      = fCurrent_Track->Eta();
    fCurrent_Track_Charge   = fCurrent_Track->Charge();
    if (  std::fabs(fCurrent_Track_TransMom) < 0.15 ||  std::fabs(fCurrent_Track_Eta) > 0.80 || std::fabs(fCurrent_Track_Charge) != 1 )  return false;
    fCurrent_Track_Momentum = fCurrent_Track->P();
    fCurrent_Track_Phi      = fCurrent_Track->Phi()*360/(2*TMath::Pi());
    auto fDCACoordinates    = fGetDCA( );
    fCurrent_Track_DCAXY    = fDCACoordinates[0];
    fCurrent_Track_DCAZ     = fDCACoordinates[1];
    return  true;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_TRK( )                                 {
    if ( fCurrent_Track_Charge > 0. )    {
        fQC_Tracks_P_Momentum       ->Fill(fCurrent_Track_Momentum);
        fQC_Tracks_P_TransMom       ->Fill(fCurrent_Track_TransMom);
        fQC_Tracks_P_Eta            ->Fill(fCurrent_Track_Eta);
        fQC_Tracks_P_Phi            ->Fill(fCurrent_Track_Phi);
    }
    if ( fCurrent_Track_Charge < 0. )    {
        fQC_Tracks_M_Momentum       ->Fill(fCurrent_Track_Momentum);
        fQC_Tracks_M_TransMom       ->Fill(fCurrent_Track_TransMom);
        fQC_Tracks_M_Eta            ->Fill(fCurrent_Track_Eta);
        fQC_Tracks_M_Phi            ->Fill(fCurrent_Track_Phi);
    }
    fQC_Tracks_DCAXY_P      ->Fill(fCurrent_Track_Momentum,fCurrent_Track_DCAXY);
    fQC_Tracks_DCAZ_P       ->Fill(fCurrent_Track_Momentum,fCurrent_Track_DCAZ);
    fQC_Tracks_DCAXY_PT     ->Fill(fCurrent_Track_TransMom,fCurrent_Track_DCAXY);
    fQC_Tracks_DCAZ_PT      ->Fill(fCurrent_Track_TransMom,fCurrent_Track_DCAZ);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_TRK_Kaons( )                           {
    if ( fCurrent_Track->Charge() > 0. )    {
        fQC_Kaons_P_Momentum        ->Fill(fCurrent_Track_Momentum);
        fQC_Kaons_P_TransMom        ->Fill(fCurrent_Track_TransMom);
        fQC_Kaons_P_Eta             ->Fill(fCurrent_Track_Eta);
        fQC_Kaons_P_Phi             ->Fill(fCurrent_Track_Phi);
    }
    if ( fCurrent_Track->Charge() < 0. )    {
        fQC_Kaons_M_Momentum        ->Fill(fCurrent_Track_Momentum);
        fQC_Kaons_M_TransMom        ->Fill(fCurrent_Track_TransMom);
        fQC_Kaons_M_Eta             ->Fill(fCurrent_Track_Eta);
        fQC_Kaons_M_Phi             ->Fill(fCurrent_Track_Phi);
    }
    fQC_Kaons_DCAXY_P      ->Fill(fCurrent_Track_Momentum,fCurrent_Track_DCAXY);
    fQC_Kaons_DCAZ_P       ->Fill(fCurrent_Track_Momentum,fCurrent_Track_DCAZ);
    fQC_Kaons_DCAXY_PT     ->Fill(fCurrent_Track_TransMom,fCurrent_Track_DCAXY);
    fQC_Kaons_DCAZ_PT      ->Fill(fCurrent_Track_TransMom,fCurrent_Track_DCAZ);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_PID ( )                                {
    fIsTPCAvailable     = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, fCurrent_Track) == AliPIDResponse::kDetPidOk);
    fIsTOFAvailable     = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, fCurrent_Track) == AliPIDResponse::kDetPidOk);
    fBetaFromTOFSignal  = fTOFBeta();
    fTPCSignal          = fCurrent_Track->GetTPCsignal();
    
    // PID Signal histograms
    if ( fIsTPCAvailable )    {
        fQC_PID_SignalTPC_P    ->  Fill(fCurrent_Track_Momentum,  fTPCSignal);
        fQC_PID_SignalTPC_PT   ->  Fill(fCurrent_Track_TransMom,  fTPCSignal);
    }
    if ( fIsTOFAvailable )    {
        fQC_PID_SignalTOF_P    ->  Fill(fCurrent_Track_Momentum,  fBetaFromTOFSignal);
        fQC_PID_SignalTOF_PT   ->  Fill(fCurrent_Track_TransMom,  fBetaFromTOFSignal);
    }
    
    // Sigma PID QC for all particles
    fQC_PID_Kaons();
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_PID_Kaons( )                           {
    auto fSigma_TOF     = (fPIDResponse->NumberOfSigmasTOF(fCurrent_Track,AliPID::kKaon));
    auto fSigma_TPC     = (fPIDResponse->NumberOfSigmasTPC(fCurrent_Track,AliPID::kKaon));
    
    // PID Sigmas histograms
    if ( fIsTPCAvailable )    {
        fQC_Kaons_SigmaTPC_P    ->  Fill(fCurrent_Track_Momentum,   fSigma_TPC);
        fQC_Kaons_SigmaTPC_PT   ->  Fill(fCurrent_Track_TransMom,   fSigma_TPC);
    }
    if ( fIsTOFAvailable )    {
        fQC_Kaons_SigmaTOF_P    ->  Fill(fCurrent_Track_Momentum,   fSigma_TOF);
        fQC_Kaons_SigmaTOF_PT   ->  Fill(fCurrent_Track_TransMom,   fSigma_TOF);
    }
    if ( fIsTPCAvailable && fIsTOFAvailable )
    {
        fQC_Kaons_SigmaTOF_TPC  ->  Fill(fSigma_TOF,            fSigma_TPC);
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_PID_Sel_Kaons( )                       {
    // PID Signal histograms
    if ( fIsTPCAvailable )    {
        fQC_Kaons_SignalTPC_P    ->  Fill(fCurrent_Track_Momentum,  fTPCSignal);
        fQC_Kaons_SignalTPC_PT   ->  Fill(fCurrent_Track_TransMom,  fTPCSignal);
    }
    if ( fIsTOFAvailable )    {
        fQC_Kaons_SignalTOF_P    ->  Fill(fCurrent_Track_Momentum,  fBetaFromTOFSignal);
        fQC_Kaons_SignalTOF_PT   ->  Fill(fCurrent_Track_TransMom,  fBetaFromTOFSignal);
    }
}

//_____________________________________________________________________________

Double_t    AliAnalysisTaskPhiCount::fTOFBeta( )                                {
    if ( !fIsTOFAvailable ) return -1.;
    Float_t fTrackTime(fCurrent_Track->GetTOFsignal()*1e-12);
    Float_t fTrackLngt(fCurrent_Track->GetIntegratedLength()*1e-2);
    if ( fTrackTime <= 0. || fTrackLngt <= 0. ) return -1.;
    return  fTrackLngt/(fTrackTime*TMath::C());
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsKaonCandidate ( )                       {
    auto fSigma_TOF     = std::fabs(fPIDResponse->NumberOfSigmasTOF(fCurrent_Track,AliPID::kKaon));
    auto fSigma_TPC     = std::fabs(fPIDResponse->NumberOfSigmasTPC(fCurrent_Track,AliPID::kKaon));
    
    //  CUSTOM
    if ( !fIsTPCAvailable || ( fIsTOFAvailable && fSigma_TOF >= 3 ) )       return false;
    if (  fIsTOFAvailable && fSigma_TPC > 5. )                              return false;
    if ( !fIsTOFAvailable && fSigma_TPC > 3. )                              return false;
    fTOFSigma[fnKaon]   =   (Char_t)(fSigma_TOF*10.);
    fTPCSigma[fnKaon]   =   (Char_t)(fSigma_TPC*10.);
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPhiCandidate ( TLorentzVector fPhi )    {
    if ( fPhi.Mag() < 0.75 || fPhi.Mag() > 1.25 ) return false;
    return true;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsCandidateTruPhi ( AliAODMCParticle* piKaon, AliAODMCParticle* pjKaon )
{
    if ( !piKaon || !pjKaon ) return false;
    auto    fiKaonMother    =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(piKaon->GetMother()));
    auto    fjKaonMother    =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(pjKaon->GetMother()));
    if  ( fiKaonMother->GetPdgCode() != fjKaonMother->GetPdgCode() || fjKaonMother->GetPdgCode() != 333 ) return false;
    return fiKaonMother == fjKaonMother;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhiGen ( AliAODMCParticle* particle )
{
    if ( particle->GetNDaughters() != 2 ) return false;
    auto const Dau1 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterFirst()));
    auto const Dau2 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterLast()));
    
    if ( !Dau1 || !Dau2 ) return false;
    return  ( ( Dau1->GetPdgCode() == -Dau2->GetPdgCode() ) && ( abs(Dau1->GetPdgCode()) == 321 ) );
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhiRec ( AliAODMCParticle* particle )
{
    // looping over all recordable phis
    for ( int iPhi = 0; iPhi < fnPhiRec; iPhi++ )
    {
        if ( fPhiRecParticles[iPhi] == particle ) return true;
    }
    return false;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhi ( AliAODMCParticle* particle )
{
    if ( !particle ) return false;
    if ( particle->GetPdgCode() == 333 ) return true;
    else return false;
}

//_____________________________________________________________________________
