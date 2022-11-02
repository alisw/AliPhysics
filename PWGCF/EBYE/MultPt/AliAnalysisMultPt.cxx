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
#include <AliHeader.h>
#include "TCanvas.h"


using namespace std;            // std namespace: so you can do things like 'cout'
ClassImp(AliAnalysisMultPt)     // classimp: necessary for root
//___________________________________________________________________________________________
AliAnalysisMultPt::AliAnalysisMultPt() : AliAnalysisTaskSE(),
    fAOD(0), fMC(0), fOutputList(0), MultHist(0), MultPtHist(0), MultPtHistRec(0), MultPtHistGen(0), MultHistRatio(0), hV0MVsPtGen(0), hV0MVsPtRec(0), fIsMC(0), fPtmin(0.2), fPtmax(5.0), fEtaMin(-0.8), fEtaMax(0.8), fBit(96), fPVzMax(10.0), fPVzMin(0.0), fEventCuts(0)
{


}
//___________________________________________________________________________________________
AliAnalysisMultPt::AliAnalysisMultPt(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fMC(0), fOutputList(0), MultHist(0), MultPtHist(0), MultPtHistRec(0), MultPtHistGen(0), MultHistRatio(0), hV0MVsPtGen(0), hV0MVsPtRec(0), fIsMC(0), fPtmin(0.2), fPtmax(5.0), fEtaMin(-0.8), fEtaMax(0.8), fBit(96), fPVzMax(10.0), fPVzMin(0.0), fEventCuts(0)
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

    const Int_t nV0Mbins = 9;
    Double_t V0Mbins[nV0Mbins+1] = {0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

    const Int_t nPtbins = 53;
    Double_t Ptbins[nPtbins+1] = {
        0.0, 0.2, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
        0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
        1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
        4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 18.0, 20.0 };

    hV0MVsPtGen = new TH2D("hV0MVsPtGen", "; V0M Percentiles;#it{p}_{T} (GeV/#it{c})", nV0Mbins, V0Mbins, nPtbins, Ptbins);
    fOutputList->Add(hV0MVsPtGen);

    hV0MVsPtRec = new TH2D("hV0MVsPtRec", "; V0M Percentiles;#it{p}_{T} (GeV/#it{c})", nV0Mbins, V0Mbins, nPtbins, Ptbins);
    fOutputList->Add(hV0MVsPtRec);

    fEventCuts.AddQAplotsToList(fOutputList);
    PostData(1, fOutputList);

}
//___________________________________________________________________________________Raw Data
void AliAnalysisMultPt::BuildData()
{
    // get an event from the analysis manager:
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    // check if there actually is an event:
    if(!fAOD) return;

    if (!fEventCuts.AcceptEvent(fAOD)) return;
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
        if(TMath::Abs(track->Charge()) < 0.1)continue;//only get charged tracks
        if(track->Eta() > fEtaMax || track->Eta() < fEtaMin) continue;//eta cut
        if(track->Pt() < fPtmin || track->Pt() > fPtmax) continue; //pt cut
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
        if(TMath::Abs(track->Charge()) < 0.1) continue;//only get charged tracks
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
    if(!fAOD) { return; }

    fMC = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMC) {
        Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,
                __LINE__);
        this->Dump();
        return;
    }

    if (!fEventCuts.AcceptEvent(fAOD)) { return; }

    Float_t Centrality = fEventCuts.GetCentrality();

    TClonesArray *stack = nullptr;
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) { return; }

    //______________________________________________________RECONSTRUCTED part:
    int MultRec = 0; //Reconstructed Multiplicity
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if (!track) continue;
        if (!track->TestFilterBit(fBit)) continue;
        int label = TMath::Abs(track->GetLabel());
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
        //mcTrack is the reconstructed track
        if (!mcTrack) continue;
        if (!mcTrack->IsPhysicalPrimary()) continue;
        Float_t etarec  = mcTrack->Eta(); //reconstructed Eta
        Float_t pTrec   = mcTrack->Pt(); //reconstructed pT
        if(TMath::Abs(mcTrack->Charge()) < 0.1)continue;//only get charged tracks
        if(etarec > fEtaMax || etarec < fEtaMin) continue;//eta cut
        if(pTrec  < fPtmin || pTrec > fPtmax) continue; //pt cut
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
        if(TMath::Abs(mcTrack->Charge()) < 0.1)continue;//only get charged tracks
        if(etarec > fEtaMax || etarec < fEtaMin) continue;//eta cut
        if(pTrec  < fPtmin || pTrec > fPtmax) continue; //pt cut
        MultPtHistRec->Fill(MultRec, pTrec);
        hV0MVsPtRec->Fill(Centrality,pTrec);
    }
    //______________________________________________________GENERATED part:

    int nMCTracks;
    if (!stack) nMCTracks = 0;
    else nMCTracks = stack->GetEntries();

    AliAODMCHeader *mcHeader = 0;
    mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
        printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
        return;
    }

    TClonesArray* fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
    if (!fMCArray){
        Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }

    int MyMultGen  = 0;  //Generated Multiplicity

    for (Int_t i= 0; i < fMC->GetNumberOfTracks(); ++i) {

        AliAODMCParticle *particle = (AliAODMCParticle*)fMC->GetTrack(i);
        if (!particle) { continue; }

        if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, mcHeader, stack)) { continue; }

        Double_t vz = particle->Zv();
        if (TMath::Abs(vz) > 10.0) { continue; }

        if (!particle->IsPhysicalPrimary()) { continue; }

        Float_t etaMCgen  = particle->Eta();
        Float_t pTMCgen = particle->Pt();

        if (etaMCgen > fEtaMax || etaMCgen<fEtaMin ) { continue; }

        if (particle->Pt() < fPtmin || particle->Pt() > fPtmax) { continue; }

        if (TMath::Abs(particle->Charge()) < 0.1) { continue; }

        MyMultGen++;

    }

    for (Int_t i= 0; i < fMC->GetNumberOfTracks(); ++i) {

        AliAODMCParticle *particle = (AliAODMCParticle*)fMC->GetTrack(i);
        if (!particle) { continue; }

        if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, mcHeader, fMCArray)) { continue; }

        Double_t vz = particle->Zv();
        if (TMath::Abs(vz) > 10.0) { continue; }

        if (!particle->IsPhysicalPrimary()) { continue; }

        Float_t etaMCgen  = particle->Eta();
        Float_t pTMCgen = particle->Pt();

        if (etaMCgen > fEtaMax || etaMCgen < fEtaMin ) continue;
        
        if (particle->Pt() < fPtmin || particle->Pt() > fPtmax) continue;

        if (TMath::Abs(particle->Charge()) < 0.1)
            continue;

        MultPtHistGen->Fill(MyMultGen, pTMCgen);
        hV0MVsPtGen->Fill(Centrality,pTMCgen);

    }
    //Generated Mult vs Reconstructed Mult
    MultHistRatio->Fill(MultRec, MyMultGen);
    return;
}
//___________________________________________________________________________________
void AliAnalysisMultPt::UserExec(Option_t *)
{
    if (fIsMC) { BuildMC(); }
    if (!fIsMC) { BuildData(); }
    PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisMultPt::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
