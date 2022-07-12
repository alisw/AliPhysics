/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisMultPt:
// Description: Analysis task to get multiplicity
// and pT distributions
// Author: Negin Alizadehvandchali
// (negin.alizadehvandchali@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TH3D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TString.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisMultPt.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliAODMCHeader.h"
#include "TCanvas.h"


using namespace std;            // std namespace: so you can do things like 'cout'
ClassImp(AliAnalysisMultPt)     // classimp: necessary for root
//___________________________________________________________________________________________

AliAnalysisMultPt::AliAnalysisMultPt() : AliAnalysisTaskSE(),
fAOD(0), fOutputList(0), MultHist(0), MultPtHist(0), MultPtHistRec(0), MultPtHistGen(0), MultHistRatio(0), fIsMC(0), fPtmin(0.2), fPtmax(5.0), fEtaMin(-0.8), fEtaMax(0.8), fBit(96), fPVzMax(8.0), fPVzMin(0.0), fChi2DoF(3), fTPCNCrossedRows(70), fIsRunFBOnly(0), fTPCNcls(70), fIsPileUpCuts(0), fPileUpLevel(2), fGenName("Hijing"), fEventCuts(0)
{
    

}
//___________________________________________________________________________________________

AliAnalysisMultPt::AliAnalysisMultPt(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), MultHist(0), MultPtHist(0), MultPtHistRec(0), MultPtHistGen(0), MultHistRatio(0), fIsMC(0), fPtmin(0.2), fPtmax(5.0), fEtaMin(-0.8), fEtaMax(0.8), fBit(96), fPVzMax(8.0), fPVzMin(0.0), fChi2DoF(3), fTPCNCrossedRows(70), fIsRunFBOnly(0), fTPCNcls(70), fIsPileUpCuts(0), fPileUpLevel(2), fGenName("Hijing"), fEventCuts(0)
{
    // Default constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 writes into a TList container
    DefineOutput(1, TList::Class());
}
//___________________________________________________________________________________________

AliAnalysisMultPt::~AliAnalysisMultPt()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}
//___________________________________________________________________________________________

void AliAnalysisMultPt::UserCreateOutputObjects()
    {
    // Initialize output list of containers
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    //______________________________________ Raw Data:
    MultHist = new TH1D("MultHist", "Mult", 4000, 0, 4000);
    fOutputList->Add(MultHist);
    MultHist->SetTitle("Multiplicity Distribution - Data");
    MultHist->SetXTitle("N_{ch}");
    MultHist->SetYTitle("Number of Events");
    MultHist->SetMarkerSize(1.2);
              
    MultPtHist = new TH2D("MultPtHist", "MultpT", 4000, 0, 4000, 50, 0, 5);
    fOutputList->Add(MultPtHist);
    MultPtHist->SetTitle("pT vs Multiplicity - Data");
    MultPtHist->SetXTitle("N_{ch}");
    MultPtHist->SetYTitle("pT");
    MultPtHist->SetMarkerSize(1.2);
        //______________________________________ Reconstructed:
    MultPtHistRec = new TH2D("MultPtHistRec", "MultPt-Rec", 4000, 0, 4000, 50, 0, 5);
    fOutputList->Add(MultPtHistRec);
    MultPtHistRec->SetTitle("pT vs Reconstructed Multiplicity");
    MultPtHistRec->SetXTitle("Reconstructed N_{ch}");
    MultPtHistRec->SetYTitle("pT");
    MultPtHistRec->SetMarkerSize(1.2);
    //______________________________________ Generated:
    MultPtHistGen = new TH2D("MultPtHistGen", "MultPt-Gen", 4000, 0, 4000, 50, 0, 5);
    fOutputList->Add(MultPtHistGen);
    MultPtHistGen->SetTitle("pT vs Generated Multiplicity");
    MultPtHistGen->SetXTitle("Generated N_{ch}");
    MultPtHistGen->SetYTitle("pT");
    MultPtHistGen->SetMarkerSize(1.2);
    //______________________________________ Ratio of Gen mult over Rec mult:
    MultHistRatio = new TH2D("MultHistRatio", "Ratio",  4000, 0, 4000, 4000, 0, 4000);
    fOutputList->Add(MultHistRatio);
    MultHistRatio->SetTitle("Generated Multiplicity vs Reconstructed");
    MultHistRatio->SetXTitle("Reconstructed N_{ch}");
    MultHistRatio->SetYTitle("Generated N_{ch}");
    MultHistRatio->SetMarkerSize(1.2);

    // add the list to our output file
    PostData(1, fOutputList);

    }
//___________________________________________________________________________________Raw Data
void AliAnalysisMultPt::BuildData()
{
    // get an event from the analysis manager:
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    // check if there actually is an event:
    if(!fAOD) return;
    
    if(fIsPileUpCuts){
        fEventCuts.fUseITSTPCCluCorrelationCut = fPileUpLevel;
        if (!fEventCuts.AcceptEvent(fAOD)) return;
    }
    //making a cut in pvz -8 to 8cm
    const AliAODVertex* PrimaryVertex = fAOD->GetVertex(0);
    if(!PrimaryVertex) return;
    Float_t PVz = PrimaryVertex->GetZ();
    if(fabs(PVz)>fPVzMax) return;
    if(fabs(PVz)<fPVzMin) return;

    int Mult=0;
    // loop over all these tracks:
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;
        if(!track->TestFilterBit(fBit)) continue;
        if(track->Charge()==0)continue;//only get charged tracks
        if(track->Eta() > fEtaMax || track->Eta() < fEtaMin) continue;//eta cut
        if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
        if(!fIsRunFBOnly){
        if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls , crossed rows and chi2/dof
      }
        Mult++;
    }
    //Number of events vs multiplicity
    MultHist->Fill(Mult);
    for( Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
        // get a track (type AliAODTrack) from the event:
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        // if we failed, skip this track:
        if(!track) continue;
        if(!track->TestFilterBit(fBit)) continue;
        if(!fIsRunFBOnly){
            if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls , crossed rows and chi2/dof
        }
        if(track->Charge()==0) continue;//only get charged tracks
        if(track->Eta() > fEtaMax || track->Eta() < fEtaMin) continue;//eta cut
        if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
        MultPtHist->Fill(Mult, track->Pt());
        }
    return;
}
//___________________________________________________________________________________MONTE-CARLO:
void AliAnalysisMultPt::BuildMC()
{
    // get an event from the analysis manager:
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    // check if there actually is an event:
    if(!fAOD) return;
    
    if(fIsPileUpCuts){
        fEventCuts.fUseITSTPCCluCorrelationCut = fPileUpLevel;
        if (!fEventCuts.AcceptEvent(fAOD)) return;
    }
    //making a cut in pvz -8 to 8cm
    const AliAODVertex* PrimaryVertex = fAOD->GetVertex(0);
    if(!PrimaryVertex) return;
    Float_t PVz = PrimaryVertex->GetZ();
    if(fabs(PVz)>fPVzMax) return;
    if(fabs(PVz)<fPVzMin) return;
    
    TClonesArray *stack =0;
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if(!stack) return;
    
    //______________________________________________________RECONSTRUCTED part:
    int MultRec=0; //Reconstructed Multiplicity
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;
        if(!track->TestFilterBit(fBit)) continue;
        int label = TMath::Abs(track->GetLabel());
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
        //mcTrack is the reconstructed track
        if (!mcTrack) continue;
        if (!mcTrack->IsPhysicalPrimary()) continue;
        Float_t etarec  = mcTrack->Eta(); //reconstructed Eta
        Float_t pTrec   = mcTrack->Pt(); //reconstructed pT
        Float_t ncl = track->GetTPCNcls();
        Float_t crossedrows = track->GetTPCNCrossedRows();
        Float_t chi2 = track->Chi2perNDF();
        if(abs(mcTrack->Charge())<=1)continue;//only get charged tracks
        if(etarec > fEtaMax || etarec < fEtaMin) continue;//eta cut
        if(pTrec  < fPtmin || pTrec > fPtmax) continue; //pt cut
        if(!fIsRunFBOnly){
          if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls, crossed rows and chi2/dof
        }
        MultRec++;
    }
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;
        if(!track->TestFilterBit(fBit)) continue;
        int label = TMath::Abs(track->GetLabel());
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
        //mcTrack is the reconstructed track
        if (!mcTrack) continue;
        if (!mcTrack->IsPhysicalPrimary()) continue;
        Float_t etarec  = mcTrack->Eta(); //reconstructed Eta
        Float_t pTrec   = mcTrack->Pt(); //reconstructed pT
        Float_t ncl = track->GetTPCNcls();
        Float_t crossedrows = track->GetTPCNCrossedRows();
        Float_t chi2 = track->Chi2perNDF();
        if(abs(mcTrack->Charge())<=1)continue;//only get charged tracks
        if(etarec > fEtaMax || etarec < fEtaMin) continue;//eta cut
        if(pTrec  < fPtmin || pTrec > fPtmax) continue; //pt cut
        if(!fIsRunFBOnly){
          if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls, crossed rows and chi2/dof
        }
        MultPtHistRec->Fill(MultRec, pTrec);
    }
    //______________________________________________________GENERATED part:
    int nMCTracks;
    if (!stack) nMCTracks = 0;
    else nMCTracks = stack->GetEntries();
    
    if(fIsPileUpCuts){
        AliAODMCHeader *mcHeader = 0;
        mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
        if(!mcHeader) {
            printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
            return;
        }
        Bool_t isPileupInGeneratedEvent = kFALSE;
        isPileupInGeneratedEvent = AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader, fGenName);
        if(isPileupInGeneratedEvent) return;
    }
    
    int MultGen  = 0;  //Generated Multiplicity
    for ( Int_t i(0); i < nMCTracks; i++) {
        AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
        if(!p1) continue;
        if(abs(p1->Charge())<=1)continue;//only get charged tracks
        if(!p1->IsPhysicalPrimary()) continue;
        Float_t etaMCgen  = p1->Eta();
        Float_t pTMCgen = p1->Pt();
        if(etaMCgen > fEtaMax || etaMCgen<fEtaMin ) continue;
        if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
        MultGen++;
    }
    for ( Int_t i(0); i < nMCTracks; i++) {
        AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
        if(!p1) continue;
        if(abs(p1->Charge())<=1)continue;//only get charged tracks
        if(!p1->IsPhysicalPrimary()) continue;
        Float_t etaMCgen  = p1->Eta();
        Float_t pTMCgen = p1->Pt();
        if(etaMCgen > fEtaMax || etaMCgen<fEtaMin ) continue;
        if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
        MultPtHistGen->Fill(MultGen, pTMCgen);
    }
    //Generated Mult vs Reconstructed Mult
    MultHistRatio->Fill(MultRec, MultGen);
    return;
}
//___________________________________________________________________________________
void AliAnalysisMultPt::UserExec(Option_t *)
{
    if(fIsMC) BuildMC();
    if(!fIsMC) BuildData();
    PostData(1, fOutputList);
    cout<<"Flag is: "<<fIsMC<<endl;
}
//________________________________________________________________________
void AliAnalysisMultPt::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
    











    
    
  

