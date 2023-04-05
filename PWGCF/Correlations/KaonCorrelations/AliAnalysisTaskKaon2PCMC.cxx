/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskKaon2PCMC:
// Description: Analysis task to calculate Two-Particle Angular Correlation Functions of Neutral and Charged Kaons
// Author: Anjaly Sasikumar Menon
// (anjaly.sasikumar.menon@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliMultiInputEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisAlien.h"
#include "AliEventCuts.h"
#include <THnSparse.h>
#include "AliAnalysisTaskKaon2PCMC.h"


class AliAnalysisTaskKaon2PCMC;    // This analysis class
using namespace std;            // std namespace: so you can do things like 'cout'
ClassImp(AliAnalysisTaskKaon2PCMC) // classimp: necessary for root
const Double_t Pi = TMath::Pi();

AliAnalysisTaskKaon2PCMC::AliAnalysisTaskKaon2PCMC() : AliAnalysisTaskSE(),  
//fAOD(0), 
fOutputList(0),
fMCK0(0),
fMCKpos(0),
fMCKneg(0),
fHistKpKnMC(0),
fHistK0KchMC(0),
fHistGenMultiplicity(0) 

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PCMC::AliAnalysisTaskKaon2PCMC(const char* name) : AliAnalysisTaskSE(name),
//fAOD(0), 
fOutputList(0), 
fMCK0(0),
fMCKpos(0),
fMCKneg(0),
fHistKpKnMC(0),
fHistK0KchMC(0),
fHistGenMultiplicity(0)
{
    // constructor
    DefineInput(0, TChain::Class());    
    DefineOutput(1, TList::Class());    
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PCMC::~AliAnalysisTaskKaon2PCMC()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PCMC::UserCreateOutputObjects()
{
    // create output objects
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file.

    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);       

    //+++++++++++++++++++++ MC ++++++++++++++++++++++++++

    Int_t bins[4] = {100,32,32,100};
    Double_t min[4] = {0.2,0,-0.8,0.4};
    Double_t max[4] = {1.0,2*Pi,0.8,0.6};
    
    fMCK0 = new THnSparseF("fMCK0","fMCK0",4,bins,min,max);
    fMCK0->GetAxis(0)->SetTitle("p_{T} of K^{0}_{S}");
    fMCK0->GetAxis(1)->SetTitle("#phi");
    fMCK0->GetAxis(2)->SetTitle("#eta");
    fMCK0->GetAxis(3)->SetTitle("mass");

    fMCKpos = new THnSparseF("fMCKpos","fMCKpos",4,bins,min,max);
    fMCKpos->GetAxis(0)->SetTitle("p_{T} of K^{+}");
    fMCKpos->GetAxis(1)->SetTitle("#phi of K^{+}");
    fMCKpos->GetAxis(2)->SetTitle("#eta of K^{+}");
    fMCKpos->GetAxis(3)->SetTitle("mass of K^{+}");

    fMCKneg = new THnSparseF("fMCKneg","fMCKneg",4,bins,min,max);
    fMCKneg->GetAxis(0)->SetTitle("p_{T} of K^{-}");
    fMCKneg->GetAxis(1)->SetTitle("#phi of K^{-}");
    fMCKneg->GetAxis(2)->SetTitle("#eta of K^{-}");
    fMCKneg->GetAxis(3)->SetTitle("mass of K^{-}");

    fHistKpKnMC = new TH2F("fHistKpKnMC","K^{+}-K^{-} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKpKnMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKpKnMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistKpKnMC->SetOption("colz");

    fHistK0KchMC = new TH2F("fHistK0KchMC","K^{+}-K^{-} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistK0KchMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistK0KchMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistK0KchMC->SetOption("colz");

    fHistGenMultiplicity = new TH1D ("fHistGenMultiplicity","fHistGenMultiplicity",500,0,500);

    fOutputList->Add(fMCK0);
    fOutputList->Add(fMCKpos);
    fOutputList->Add(fMCKneg);
    fOutputList->Add(fHistKpKnMC);
    fOutputList->Add(fHistK0KchMC);
    fOutputList->Add(fHistGenMultiplicity);

    //fEventCuts.AddQAplotsToList(fOutputList);

    PostData(1, fOutputList);
    
}


Bool_t AliAnalysisTaskKaon2PCMC::SelectK0TracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectK0 =  mcPartPdg==310&& (isPhysPrim); 
    if (SelectK0) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PCMC::SelectKPosTracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKpos =  mcPartPdg==321&& (isPhysPrim); 
    if (SelectKpos) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PCMC::SelectKNegTracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKneg =  mcPartPdg==-321&& (isPhysPrim); 
    if (SelectKneg) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PCMC::SelectKchTracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKch =  mcPartPdg==321&& mcPartPdg==-321&& (isPhysPrim); 
    if (SelectKch) return kFALSE;
    return kTRUE;
}

void AliAnalysisTaskKaon2PCMC::UserExec(Option_t *)
{

Int_t nAcceptedParticles =0;
AliMCParticle *mcTrack = 0x0;
// /*
// anjaly

// */
fmcEvent  = dynamic_cast<AliMCEvent*> (MCEvent());
if(!fmcEvent){
    Printf("No MC particle branch found");
    return;
    }

AliVVertex * mcVertex = (AliVVertex*)fmcEvent->GetPrimaryVertex();
fPV[2] = mcVertex->GetZ();
if (TMath::Abs(fPV[2])>=10) return;

Int_t nMCTracks = fmcEvent->GetNumberOfTracks(); // MC Truth, Total number of tracks per event
AliVTrack *genTrackMix = 0x0;

for (Int_t i = 0; i < nMCTracks; i++){
    AliMCParticle *mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!mcTrack) {
        Error("ReadEventAODMC", "Could not receive particle %d", i);
        continue;
    }

    Double_t trackPseudorap = mcTrack->Eta();
    if( mcTrack->IsPhysicalPrimary()&&mcTrack->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8))) {
                nAcceptedParticles += 1;
            }
    }

cout << "number of accepted particles from MC tracks is"<< nAcceptedParticles << endl;
fHistGenMultiplicity->Fill(nAcceptedParticles);

AliMCParticle *mcMotherParticle = 0x0;
AliMCParticle* daughter0 = 0x0;
AliMCParticle* daughter1 = 0x0;
Bool_t SelectK0;
Bool_t SelectKpos;
Bool_t SelectKneg;

for (Int_t i = 0; i < nMCTracks; i++){
	mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!mcTrack) continue;
        
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();

    SelectK0 =  mcPartPdg==310&& (isPhysPrim);  // (8010 total. When we include only primary, 7546)-310 for neutral
    SelectKpos = mcPartPdg==321&& (isPhysPrim); // 321 is the code for positive kaons
    SelectKneg = mcPartPdg==-321&& (isPhysPrim); // 321 is the code for positive kaons

    Double_t TrackPt = mcTrack->Pt();
    Double_t TrackPhi = mcTrack->Phi();
    Double_t TrackEta = mcTrack->Eta();
    Double_t TrackMass = mcTrack->M();
    
    Double_t KaonVariables[4]= {TrackPt, TrackPhi, TrackEta, TrackMass};
    if(SelectK0) fMCK0->Fill(KaonVariables);
    if(SelectKpos) fMCKpos->Fill(KaonVariables);
    if(SelectKneg) fMCKneg->Fill(KaonVariables);
   
	Bool_t TrIsPrim = mcTrack->IsPhysicalPrimary();
	Bool_t TrCharge = (mcTrack->Charge())!=0;
    Short_t cha;
    if (mcTrack->Charge()>0) cha=1.;
    else if (mcTrack->Charge()<0) cha= -1.;
    else cha =0;

}


for (Int_t i = 0; i < nMCTracks; i++){
    AliMCParticle *mcTrack1 = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!SelectKPosTracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    Double_t trackPseudorap = mcTrack1->Eta();
    if(! (mcTrack1->IsPhysicalPrimary()&&mcTrack1->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = i+1; j < nMCTracks; j++){
        AliMCParticle *mcTrack2 = (AliMCParticle*)fmcEvent->GetTrack(j);
        if (!SelectKNegTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t trackPseudorap = mcTrack2->Eta();
        if(! (mcTrack2->IsPhysicalPrimary()&&mcTrack2->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistKpKnMC);
    }

}


for (Int_t i = 0; i < nMCTracks; i++){
    AliMCParticle *mcTrack1 = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!SelectK0TracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    Double_t trackPseudorap = mcTrack1->Eta();
    if(! (mcTrack1->IsPhysicalPrimary()&&mcTrack1->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = i+1; j < nMCTracks; j++){
        AliMCParticle *mcTrack2 = (AliMCParticle*)fmcEvent->GetTrack(j);
        if (!SelectKchTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t trackPseudorap = mcTrack2->Eta();
        if(! (mcTrack2->IsPhysicalPrimary()&&mcTrack2->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistK0KchMC);
    }

}

//fHistNEvents->Fill(1);
//fHistMult->Fill(iTracks);
//fHistCent->Fill(CentV0M);                          
PostData(1, fOutputList);
}

//****************** filling functions ***********************************************

void AliAnalysisTaskKaon2PCMC::Fill2DHist(Double_t DPhi, Double_t DEta, TH3F* hist, Double_t fWeight=1){
    DPhi = fabs(DPhi);
    DEta = fabs(DEta);
    if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
    hist->Fill(DPhi,DEta,fWeight);
    hist->Fill(DPhi,-DEta,fWeight);
    if (DPhi < 0.5*Pi ) {
        hist->Fill(-DPhi, DEta,fWeight);
        hist->Fill(-DPhi,-DEta,fWeight); 
    }
    else {
        hist->Fill(2*Pi-(DPhi), DEta,fWeight);
        hist->Fill(2*Pi-(DPhi),-DEta,fWeight);
    }
}

void AliAnalysisTaskKaon2PCMC::Fill2DHistMCTruth(Double_t DPhi, Double_t DEta, TH2F* hist){
    DPhi = fabs(DPhi);
    DEta = fabs(DEta);
    if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
    hist->Fill(DPhi,DEta);
    hist->Fill(DPhi,-DEta);
    if (DPhi < 0.5*Pi ) {
        hist->Fill(-DPhi, DEta);
        hist->Fill(-DPhi,-DEta); 
    }
    else {
        hist->Fill(2*Pi-(DPhi), DEta);
        hist->Fill(2*Pi-(DPhi),-DEta);
    }
}






//_____________________________________________________________________________
void AliAnalysisTaskKaon2PCMC::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________


