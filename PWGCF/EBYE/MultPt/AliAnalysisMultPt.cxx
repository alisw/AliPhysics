//////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisMultPt:
// Description: Analysis task get multiplicity
// and pT distributions
// Author: Negin Alizadehvandchali
// (negin.alizadehvandchali@cern.ch)
/////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TString.h"
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TBits.h"
#include "TRefArray.h"
#include "TAxis.h"
#include "TRandom.h"
#include "TList.h"
#include "TComplex.h"
#include "TMath.h"
#include <cmath>
#include <cstdlib>
#include "AliAODTrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisMultPt.h"



using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisMultPt)     // classimp: necessary for root
//________________________________________________________________________

AliAnalysisMultPt::AliAnalysisMultPt() : AliAnalysisTaskSE(),
fAOD(0), fOutputList(0), fHistMultPt(0), fHistMult(0), fHistMultPtMC(0), fHistMultRatio(0), fIsMC(0), fPtmin(0.2), fPtmax(5), fEta(0.8), fBit(128)

{
    

}

//________________________________________________________________________
AliAnalysisMultPt::AliAnalysisMultPt(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0), fHistMultPt(0), fHistMult(0), fHistMultPtMC(0), fHistMultRatio(0), fIsMC(0), fPtmin(0.2), fPtmax(5), fEta(0.8), fBit(128)


{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisMultPt::~AliAnalysisMultPt()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}


//________________________________________________________________________

void AliAnalysisMultPt::SetMCRead(Bool_t flag) { fIsMC = flag; }

void AliAnalysisMultPt::UserCreateOutputObjects()
    {


    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
        
    fHistMult = new TH1D("fHistMult", "Mult", 4000, 0, 4000);
    fOutputList->Add(fHistMult);
    fHistMult->SetStats(0);
    fHistMult->SetTitle("Multiplicity");
    fHistMult->SetXTitle("Multiplicity");
    fHistMult->SetMarkerSize(1.2);
              
    fHistMultPt = new TH2D("fHistMultPt", "MultpT", 4000, 0, 4000, 50, 0, 5);
    fOutputList->Add(fHistMultPt);
    fHistMultPt->SetStats(0);
    fHistMultPt->SetTitle("pT vs Multiplicity");
    fHistMultPt->SetXTitle("Multiplicity");
    fHistMultPt->SetYTitle("pT");
    fHistMultPt->SetLineColor(2);
    fHistMultPt->SetMarkerSize(1.2);

    fHistMultPtMC = new TH2D("fHistMultPtMC", "MultPt-MC", 4000, 0, 4000, 50, 0, 5);
    fOutputList->Add(fHistMultPtMC);
    fHistMultPtMC->SetStats(0);
    fHistMultPtMC->SetTitle("pT vs Generated Multiplicity");
    fHistMultPtMC->SetXTitle("Generated N_{ch}");
    fHistMultPtMC->SetYTitle("pT");
    fHistMultPtMC->SetLineColor(2);
    fHistMultPtMC->SetMarkerSize(1.2);
          
          
    fHistMultRatio = new TH2D("fHistMultRatio", "Ratio",  4000, 0, 4000, 4000, 0, 4000);
    fOutputList->Add(fHistMultRatio);
    fHistMultRatio->SetStats(0);
    fHistMultRatio->SetTitle("Generated Multiplicity vs Reconstructed");
    fHistMultRatio->SetXTitle("Reconstructed N_{ch}");
    fHistMultRatio->SetYTitle("Generated N_{ch}");
    fHistMultRatio->SetLineColor(2);
    fHistMultRatio->SetMarkerSize(1.2);


    // add the list to our output file
    PostData(1, fOutputList);
    //PostData(1, mixer);
    }

//________________________________________________________________________

void AliAnalysisMultPt::UserExec(Option_t *)
{

////////////////////////////////////////////////////////////////////////////////Track-based part (Either data or reconstructed MC):
    // get an event from the analysis manager:
    fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
    // check if there actually is an event:
    if(!fAOD) return;
    // ::Fatal("AliAnalysisTaskMultPt::UserExec", "No AOD event found, check the event handler.");
     //For charged particle
    Double_t Mult=0.0;
    
    // loop over all these tracks:
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
         // get a track (type AliAODTrack) from the event:
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        // if we failed, skip this track
        if(!track) continue;
        Int_t charge = track->Charge();
        if(TMath :: Abs(charge)>0) {
            if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue;
            if(fabs(track->Eta()) > fEta) continue;   //eta cut
            if(track->TestFilterBit(fBit)){
                Mult++;
            }
        }
    }
    //Number of events vs multiplicity
    fHistMult->Fill(Mult);

    for( Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {
         // get a track (type AliAODTrack) from the event:
         AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
         // if we failed, skip this track:
         if(!track) continue;
         Int_t charge = track->Charge();
        
         //Count the charge particles:
         if(TMath :: Abs(charge)>0) {
             if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue;
             if(fabs(track->Eta()) > fEta) continue;   //eta cut
             if(track->TestFilterBit(fBit)){
                 fHistMultPt->Fill(Mult, track->Pt());
             }
         }
        
    }
/////////////////////////////////////////////////////////////////////////////////// Generated part:
    if(fIsMC){
        TClonesArray *stack = 0;
        TList *lst = fAOD->GetList();
        stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
     
        int nMCTracks;
        if (!stack) nMCTracks = 0;
        else nMCTracks = stack->GetEntries();
    
        Double_t mNMC  = 0.0;
        for ( Int_t i(0); i < nMCTracks; i++) {
            AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
            if (!p1) continue;
            if(!p1->IsPhysicalPrimary()) continue;
            Int_t charge = p1->Charge();
            if(TMath :: Abs(charge)>0) {
                if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
                if(fabs(p1->Eta()) > fEta ) continue;
                mNMC++;
            }
            
        }
        
        
        for ( Int_t i(0); i < nMCTracks; i++) {
            AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
            if (!p1) continue;
            if(!p1->IsPhysicalPrimary()) continue;
            Int_t charge = p1->Charge();
            if(TMath :: Abs(charge)>0) {
                if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
                if(fabs(p1->Eta()) > fEta ) continue;
                fHistMultPtMC->Fill(mNMC, p1->Pt());
            }
        }
        //Generated Mult vs Reconstructed Mult
        fHistMultRatio->Fill(Mult, mNMC);
    }
    
    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisMultPt::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
    











    
    
  
