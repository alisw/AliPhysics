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
#include "AliAnalysisTaskPhiCount.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskPhiCount;
class AliPIDResponse;
 
using namespace std;

ClassImp(AliAnalysisTaskPhiCount)

AliAnalysisTaskPhiCount::AliAnalysisTaskPhiCount() : AliAnalysisTaskSE(),
fAOD(0), fMCD(0), AODMCTrackArray(0), fKaonCandidate(0), fPhiCandidate(0), fKaonEfficiency(0), fPhiEfficiency(0), fAnalysisOutputList(0), fQCOutputList(0), fHistVertex0(0), fHistVertex1(0), fHistTPCPID0(0), fHistTPCPID1(0), fHistTPCPID2(0), fHistTOFPID0(0), fHistTOFPID1(0), fHistTOFPID2(0), fHistTOFPID3(0), fHistTPCPID3(0), fPIDResponse(0)
{
    
}

//_____________________________________________________________________________

AliAnalysisTaskPhiCount::AliAnalysisTaskPhiCount(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fMCD(0), AODMCTrackArray(0), fKaonCandidate(0), fPhiCandidate(0), fKaonEfficiency(0), fPhiEfficiency(0), fAnalysisOutputList(0), fQCOutputList(0), fHistVertex0(0), fHistVertex1(0), fHistTPCPID0(0), fHistTPCPID1(0), fHistTPCPID2(0), fHistTOFPID0(0), fHistTOFPID1(0), fHistTOFPID2(0), fHistTOFPID3(0), fHistTPCPID3(0), fPIDResponse(0)
{
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

AliAnalysisTaskPhiCount::~AliAnalysisTaskPhiCount()
{
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

void    AliAnalysisTaskPhiCount::UserCreateOutputObjects()
{
    // Various utility Histograms TList initialisation
    fAnalysisOutputList     = new TList();
    fAnalysisOutputList     ->SetOwner(kTRUE);
    fHistVertex0    = new TH1F("fHistVertex0", "Collision Vertex (FULL)", 300, -15, 15);
    fHistVertex1    = new TH1F("fHistVertex1", "Collision Vertex (CUTS)", 300, -15, 15);
    fHistTPCPID0    = new TH2F("fHistTPCPID0", "TPC Response (ALL)"     , 50, 0, 4, 100, 0, 400);
    fHistTPCPID1    = new TH2F("fHistTPCPID1", "TPC Response (Sel1)"    , 50, 0, 4, 100, 0, 400);
    fHistTPCPID2    = new TH2F("fHistTPCPID2", "TPC Response (Sel2)"    , 50, 0, 4, 100, 0, 400);
    fHistTOFPID0    = new TH2F("fHistTOFPID0", "TOF Response (ALL)"     , 50, 0, 4, 120, 0, 1.2);
    fHistTOFPID1    = new TH2F("fHistTOFPID1", "TOF Response (Sel1)"    , 50, 0, 4, 120, 0, 1.2);
    fHistTOFPID2    = new TH2F("fHistTOFPID2", "TOF Response (Sel2)"    , 50, 0, 4, 120, 0, 1.2);
    fHistEvntEff    = new TH1F("fHistEvntEff", "fHistEvntEff"           , 5,   0.5, 5.5);
    fAnalysisOutputList->Add(fHistEvntEff);
    fAnalysisOutputList->Add(fHistVertex0);
    fAnalysisOutputList->Add(fHistTPCPID0);
    fAnalysisOutputList->Add(fHistTOFPID0);
    fAnalysisOutputList->Add(fHistVertex1);
    fAnalysisOutputList->Add(fHistTPCPID1);
    fAnalysisOutputList->Add(fHistTOFPID1);
    fAnalysisOutputList->Add(fHistTPCPID2);
    fAnalysisOutputList->Add(fHistTOFPID2);
    
    PostData(1, fAnalysisOutputList);
    
    // QC utility Histograms TList initialisation
    fQCOutputList     = new TList();
    fQCOutputList     ->SetOwner(kTRUE);
    fHistTPCPID3    = new TH2F("fHistTPCPID3", "TPC Response (Sel3)"    , 50, 0.15, 4.15, 100, -10, 10);
    fHistTOFPID3    = new TH2F("fHistTOFPID3", "TOF Response (Sel3)"    , 50, 0.15, 4.15, 100, -10, 10);
    fQCOutputList->Add(fHistTPCPID3);
    fQCOutputList->Add(fHistTOFPID3);
    
    PostData(2, fQCOutputList);
    
    // PhiCandidate Tree Set-Up
    fPhiCandidate = new TTree   ("PhiCandidate",    "Data Tree for Phi Candidates");
    fPhiCandidate->Branch       ("Multiplicity",   &fMultiplicity,     "fMultiplicity/F");
    fPhiCandidate->Branch       ("Multiplicit2",   &fMultiplicit2,     "fMultiplicit2/F");
    fPhiCandidate->Branch       ("Multiplicit3",   &fMultiplicit3,     "fMultiplicit3/F");
    fPhiCandidate->Branch       ("nPhi",            &fnPhi,             "fnPhi/b");
    fPhiCandidate->Branch       ("Px",              &fPhiPx,            "fPhiPx[fnPhi]/F");
    fPhiCandidate->Branch       ("Py",              &fPhiPy,            "fPhiPy[fnPhi]/F");
    fPhiCandidate->Branch       ("Pz",              &fPhiPz,            "fPhiPz[fnPhi]/F");
    fPhiCandidate->Branch       ("InvMass",         &fInvMass,          "fInvMass[fnPhi]/F");
    fPhiCandidate->Branch       ("iKaon",           &fiKaon,            "fiKaon[fnPhi]/b");
    fPhiCandidate->Branch       ("jKaon",           &fjKaon,            "fjKaon[fnPhi]/b");
    
    if ( kPhibool )                 PostData(3, fPhiCandidate);
    
    // KaonCandidate Tree Set-Up
    fKaonCandidate = new TTree ("KaonCandidate",    "Data Tree for Kaon Candidates");
    fKaonCandidate->Branch     ("Multiplicity",     &fMultiplicity,     "fMultiplicity/F");
    fKaonCandidate->Branch     ("nKaon",            &fnKaon,            "fnKaon/b");
    fKaonCandidate->Branch     ("Px",               &fKaonPx,           "fKaonPx[fnKaon]/F");
    fKaonCandidate->Branch     ("Py",               &fKaonPy,           "fKaonPy[fnKaon]/F");
    fKaonCandidate->Branch     ("Pz",               &fKaonPz,           "fKaonPz[fnKaon]/F");
    fKaonCandidate->Branch     ("Charge",           &fCharge,           "fCharge[fnKaon]/B");
    fKaonCandidate->Branch     ("TOFSigma",         &fTOFSigma,         "fTOFSigma[fnKaon]/B");
    fKaonCandidate->Branch     ("TPCSigma",         &fTPCSigma,         "fTPCSigma[fnKaon]/B");
    
    if ( kKaonbool )                PostData(4, fKaonCandidate);

    fPhiEfficiency = new TTree  ("PhiEfficiency",   "MC Tree for Phi Efficiency");
    fPhiEfficiency->Branch      ("nPhi",            &fnPhiTru,          "fnPhiTru/b");
    fPhiEfficiency->Branch      ("Px",              &fPhiTruPx,         "fPhiTruPx[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Py",              &fPhiTruPy,         "fPhiTruPy[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Pz",              &fPhiTruPz,         "fPhiTruPz[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Selection",       &fSelection,        "fSelection[fnPhiTru]/b");
    
    if ( kPhibool   &&  kMCbool )   PostData(5, fPhiEfficiency);
    
    fKaonEfficiency = new TTree ("KaonEfficiency",   "MC Tree for Kaon Efficiency");
    
    if ( kKaonbool  &&  kMCbool )   PostData(6, fKaonEfficiency);
}

//_____________________________________________________________________________

void    AliAnalysisTaskPhiCount::fSetZero()
{
    //Setting all counters and global variables to zero
    fMultiplicity   =   0;
    fnPhi           =   0;
    fnPhiTru        =   0;
    fnKaon          =   0;
    fnPhiTru        =   0;
}

//_____________________________________________________________________________

void    AliAnalysisTaskPhiCount::fPostData( Bool_t fEventEfficiency = false )
{
    // Setting postdata options
    
    // Post-data for TLists
    PostData(1, fAnalysisOutputList);
    PostData(2, fQCOutputList);
    
    if ( !fEventEfficiency )
    {
        // Filling data for TTrees
        fPhiCandidate   ->  Fill();
        fKaonCandidate  ->  Fill();
        fPhiEfficiency  ->  Fill();
        fKaonEfficiency ->  Fill();
    
        // Post-data for TTrees
        if ( kPhibool )                 PostData(3, fPhiCandidate);
        if ( kKaonbool )                PostData(4, fKaonCandidate);
        if ( kPhibool   &&  kMCbool )   PostData(5, fPhiEfficiency);
        if ( kKaonbool  &&  kMCbool )   PostData(6, fKaonEfficiency);
    }
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPrimaryVertexCandidate ( AliAODEvent* fCurrent_Event )
{
    // Recovering Primary Vertex from General methods and SPD
    auto    PrimaryVertexSPD    = fAOD->GetPrimaryVertexSPD();
    auto    PrimaryVertexTRK    = fAOD->GetPrimaryVertex();
            fPrimaryVertex      = PrimaryVertexTRK;
    
    // Requires the vertex is reconstructed by the SPD
    if ( !PrimaryVertexSPD  ||  PrimaryVertexSPD->GetNContributors() < 1 )
    {
        fFillVtxHist(1);
        fPostData(kTRUE);
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
            fFillVtxHist(2);
            fPostData(kTRUE);
            return false;
        }
    }
    
    // Fill the Vertex Z position histogram
    fHistVertex0->Fill(fPrimaryVertex->GetZ());
    
    if ( std::fabs(fPrimaryVertex->GetZ()) > 10. )
    {
        fFillVtxHist(3);
        fPostData(kTRUE);
        return false;
    }
    
    // Fill the Vertex Z position histogram
    fHistVertex1->Fill(fPrimaryVertex->GetZ());
    
    // Check the event is Pile-up from the SPD
    if ( fCurrent_Event->IsPileupFromSPD() )
    {
        fFillVtxHist(4);
        fPostData(kTRUE);
        return false;
    }
    
    fFillVtxHist(5);
    return  true;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsTrackCandidate ( AliAODTrack* track )
{
    // Check the track is there and has proper kinematics
    if ( !track                         || !track->TestFilterBit(5)         ) return false;
    if (  std::fabs(track->Pt()) < 0.15 ||  std::fabs(track->Eta()) > 0.80  ) return false;
    return  true;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsKaonCandidate ( AliAODTrack* track )
{
    // Check the PID is present and within desired parameters
    if ( !fPIDResponse ) return false;
    auto fbTPC       = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk);
    auto fbTOF       = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk);
    auto ffSigTOF    = std::fabs(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
    auto ffSigTPC    = std::fabs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
    
    fFillPIDHist(track,0);
    fFillPIDHist(track,3);
    
    //  CUSTOM
    if ( !fbTPC || (fbTOF && ffSigTOF > 3) )      return false;
    if ( track->Pt() >= 0.28 &&  fbTOF && ffSigTPC > 5. )   return false;
    if ( track->Pt() >= 0.28 && !fbTOF && ffSigTPC > 3. )   return false;
    fFillPIDHist(track,1);
    if ( track->Pt() <  0.28  && track->Pt() >=  0.24  && ffSigTPC > 6. )   return false;
    if ( track->Pt() <  0.24  && track->Pt() >=  0.16  && ffSigTPC > 7. )   return false;
    if ( track->Pt() <  0.16  && track->Pt() >=  0.00  && ffSigTPC > 7.5 )  return false;
    fFillPIDHist(track,2);
    return true;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fSetKaonPID ( AliAODTrack* track )
{
    // Check the PID is present
    if ( !fPIDResponse ) return false;
    auto fbTPC       = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk);
    auto fbTOF       = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk);
    auto ffSigTOF    = (fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
    auto ffSigTPC    = (fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
    
    if ( fabs(ffSigTOF) <= 10 && fbTOF )    fTOFSigma[fnKaon]= static_cast<Char_t>(ffSigTOF*10);
    else                                    fTOFSigma[fnKaon]= static_cast<Char_t>(-127);
    if ( fabs(ffSigTPC) <= 10 && fbTPC )    fTPCSigma[fnKaon]= static_cast<Char_t>(ffSigTPC*10);
    else                                    fTPCSigma[fnKaon]= static_cast<Char_t>(-127);
    
    return  fIsKaonCandidate(track);
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsKaonTruPhi ( AliAODMCParticle* piKaon, AliAODMCParticle* pjKaon )
{
    if ( !piKaon || !pjKaon ) return false;
    if ( piKaon->GetMother() == pjKaon->GetMother() && (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(pjKaon->GetMother()))->GetPdgCode() == 333 ) ) return true;
    else return false;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhiCandidate ( TLorentzVector fPhi )
{
    if ( fPhi.Mag() < 0.75 || fPhi.Mag() > 1.25 ) return false;
    return true;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhiGen ( AliAODMCParticle* particle )
{
    if ( particle->GetNDaughters() != 2 ) return false;
    auto const Dau1 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterFirst()));
    auto const Dau2 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterLast()));
    
    if ( !Dau1 || !Dau2 ) return false;
    return  (
                ( particle->GetNDaughters() == 2 ) &&
                ( Dau1->GetPdgCode() == -Dau2->GetPdgCode() ) &&
                ( abs(Dau1->GetPdgCode()) == 321 )
             );
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhiRec ( AliAODMCParticle* particle )
{
    // To be recrodable, it must come from a K^+k^- decay
    if ( fSelection[fnPhi] != 1 ) return false;
    
    // Generating Daughter 1 and 2 Particles instances
    auto const Dau1 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterFirst()));
    auto const Dau2 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterLast()));
    
    Bool_t  fbDau1  =   false;
    Bool_t  fbDau2  =   false;

    // looping over all kaons
    for ( int iKaon = 0; iKaon < fnKaon; iKaon++ )
    {
        // recovering kaon label
        if ( fKaonLabels[iKaon] == Dau1->GetLabel() ) fbDau1 = true;
        if ( fKaonLabels[iKaon] == Dau2->GetLabel() ) fbDau2 = true;
    }
    return fbDau1 && fbDau2;
}

//_____________________________________________________________________________

bool    AliAnalysisTaskPhiCount::fIsPhi ( AliAODMCParticle* particle )
{
    if ( !particle ) return false;
    if ( particle->GetPdgCode() == 333 ) return true;
    else return false;
}

//_____________________________________________________________________________

void    AliAnalysisTaskPhiCount::fFillPIDHist ( AliAODTrack * track , Int_t iIndex )
{
    if ( !track ) return;
    if ( (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk) )
    {
        auto ffSigTPC    = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
        if ( iIndex == 0 ) fHistTPCPID0->Fill(track->P(), track->GetTPCsignal());
        if ( iIndex == 1 ) fHistTPCPID1->Fill(track->P(), track->GetTPCsignal());
        if ( iIndex == 2 ) fHistTPCPID2->Fill(track->P(), track->GetTPCsignal());
        if ( iIndex == 3 ) fHistTPCPID3->Fill(track->Pt(), ffSigTPC);
    }
    if ( (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk) )
    {
        auto ffSigTOF    = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
        
        Float_t fTrackLength(track->GetIntegratedLength()*1e-2);   // Track Length in cm
        if ( fTrackLength < 0. ) return;
        
        Float_t fTrackTime(track->GetTOFsignal()*1e-12);           // Track Time in s
        if ( fTrackTime < 0. ) return;
        
        // Track Beta
        Float_t fTrackBeta = fTrackLength/(fTrackTime*TMath::C());
        
        if ( iIndex == 0 ) fHistTOFPID0->Fill(track->P(), fTrackBeta);
        if ( iIndex == 1 ) fHistTOFPID1->Fill(track->P(), fTrackBeta);
        if ( iIndex == 2 ) fHistTOFPID2->Fill(track->P(), fTrackBeta);
        if ( iIndex == 3 ) fHistTOFPID3->Fill(track->Pt(), ffSigTOF);
    }
    return;
}

//_____________________________________________________________________________

void    AliAnalysisTaskPhiCount::fFillVtxHist ( Int_t iIndex )
{
    for ( Int_t iFill = 1; iFill <= iIndex; iFill++ )
    {
        fHistEvntEff->Fill(iFill);
    }
}

//_____________________________________________________________________________

void    AliAnalysisTaskPhiCount::UserExec(Option_t *)
{
    // Recovering Event Data
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fMCD = dynamic_cast<AliMCEvent*>(MCEvent());
    
    // Check the event is there
    if ( !fAOD ) return;
    if ( !fMCD && kMCbool ) return;
    
    // Recover the MC tracks
    if ( kMCbool )      AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !AODMCTrackArray && kMCbool )     return;
    
    // Setting utility variables
    Int_t           nTrack(fAOD->GetNumberOfTracks());
    TLorentzVector  fKaon1, fKaon2, fPhi;
    
    // Check the event is there and has a primary vertex with due requirements
    if ( !fIsPrimaryVertexCandidate(dynamic_cast<AliAODEvent*>(InputEvent())) ) return;
     
    // Define and Fetch PID with Manager
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man)
    {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
    }
    
    // Setting zero all counters and global variables
    fMultiplicity   =   0;
    fnPhi           =   0;
    fnPhiTru        =   0;
    fnKaon          =   0;
    fnPhiTru        =   0;
    
    fMultiplicity = ((AliAODHeader*)fAOD->GetHeader()) -> GetRefMultiplicityComb08();
    
    fMultUtil = new AliPPVsMultUtils();
    fMultiplicit2 = fMultUtil->GetMultiplicityPercentile(fAOD,"V0M",kFALSE);
    fMultiplicit3 = fMultUtil->GetStandardReferenceMultiplicity(fAOD,kFALSE);
    
    // Looping over tracks
    for ( Int_t iTrack(0); iTrack < nTrack; iTrack++ )
    {
        // Recovering Track
        auto    fCurrent_Track  =   static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        
        // Check the track has due requirements
        if ( !fIsTrackCandidate(fCurrent_Track) ) continue;
        
        // Check the PID is present and within requirements
        if ( !fSetKaonPID(fCurrent_Track) ) continue;
        
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
            if( fIsPhiGen(fPhiTru) )    fSelection[fnPhiTru]++;
            if( fIsPhiRec(fPhiTru) )    fSelection[fnPhiTru]++;
            
            fnPhiTru++;
        }
    }
    
    // Saving output
    fPostData();
    
}

//_____________________________________________________________________________

void    AliAnalysisTaskPhiCount::Terminate(Option_t *)
{

}

//_____________________________________________________________________________
