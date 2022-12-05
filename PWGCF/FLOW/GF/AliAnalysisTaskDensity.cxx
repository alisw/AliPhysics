/* AliAnaysisTasDensityCorrelation
 *
 *  Task for transverse momentum correlation with sub-events
 * 
*/

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

// root basics
#include "TChain.h"
#include "TMath.h"

// ali Basic
#include "AliAnalysisTaskDensity.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"

// monte carlo events
#include "AliAODMCParticle.h"

// plotting 
#include "TH1D.h"


class AliAnalysisTaskDensity;

ClassImp(AliAnalysisTaskDensity); // classimp: necessary for root

AliAnalysisTaskDensity::AliAnalysisTaskDensity() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fPtSampleList(0),
    PtSubContainer(0),
    rndGenerator(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskDensity::AliAnalysisTaskDensity(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), 
    fPtSampleList(0),
    PtSubContainer(0),
    rndGenerator(0)
{
    // Set cuts for default configuration
    SetDefaultCut();

    // constructor
    DefineInput(0, TChain::Class());   
    DefineOutput(1, TList::Class());    
}
//_____________________________________________________________________________
AliAnalysisTaskDensity::~AliAnalysisTaskDensity()
{
    if(fPtSampleList) {
        delete fPtSampleList;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::UserCreateOutputObjects()
{
    // Random number generator for bs samples
    rndGenerator = new TRandom();

    Int_t nbins = 100;
    Double_t xbinmax = 100.;

    PtSubContainer = new AliPtSubEventContainer*[10];
    for(int i(0); i < 10; i++){
        PtSubContainer[i] = new AliPtSubEventContainer(Form("ptcontsample%i",i),"ptcont", 6);
        PtSubContainer[i]->SetNamedPostfix(Form("_sample%i", i+1));
        PtSubContainer[i]->Initialize(nbins, 0., xbinmax); 
        PtSubContainer[i]->InitializeTwoSub(nbins, 0., xbinmax); 
        PtSubContainer[i]->InitializeThreeSub(nbins, 0., xbinmax);
    }

    fPtSampleList = new TList();     
    fPtSampleList->SetOwner(kTRUE);     

    for(int i(0); i < 10; i++){
        fPtSampleList->Add( new TList() ); ((TList*)fPtSampleList->At(i))->SetName(Form("sample%i", i+1));
        ((TList*)fPtSampleList->At(i))->Add( PtSubContainer[i]->GetTwoSubAnalysisList() );
        ((TList*)fPtSampleList->At(i))->Add( PtSubContainer[i]->GetThreeSubAnalysisList() );
        ((TList*)fPtSampleList->At(i))->Add( PtSubContainer[i]->GetCorrList() );
        ((TList*)((TList*)fPtSampleList->At(i))->At(2))->SetName("FullTPCCorrelation");
    }

    PostData(1, fPtSampleList);
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::UserExec(Option_t *)
{
    // MC Event
    if(fMode == "HIJING"){
        fMCEvent = MCEvent();
        if(fMCEvent) ProcessMCParticles();
    }

    // get an AOD event from input file
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());   
    if(!fAOD) return;                                  

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    }

    if(!(fAOD->GetNumberOfPileupVerticesTracks()<fMaxPileup)) return;
    if(!(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ())<fLimitVertexZ)) return;

    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(!multSelection) return; 
    fCentrality = multSelection->GetMultiplicityPercentile(fMultSelMethod);


    ClearWPCounter();   // restart and initiate containr for <m> calcualtions

    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) 
    {                 
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(AcceptAODTrack(track)) continue;             // track cuts

        FillWPCounter(wp, 1., track->Pt());
        
        // Fill pt-correlation for sub-events  
        if( -0.8 >  track->Eta() > -0.4) FillWPCounter(wpSubEvent[0], 1., track->Pt());     // -0.8 < eta < -0.4
        if(std::abs(track->Eta()) < 0.2)  FillWPCounter(wpSubEvent[1], 1., track->Pt());    // -0.2 < eta < 0.2
        if( 0.4 <   track->Eta() < 0.8)  FillWPCounter(wpSubEvent[2], 1., track->Pt());     //  0.4 < eta < 0.8 
    }                             
    // Fill sub-events
    const int rn = rndGenerator->Integer(10);      

    if(wp[0][0] < fMPar) return;
    PtSubContainer[rn]->FillRecursive(wp,fCentrality);
    
    if(wpSubEvent[0][0][0] < fMPar || wpSubEvent[2][0][0] < fMPar) return;
    PtSubContainer[rn]->FillTwoSubAnalsysis(wpSubEvent[0], wpSubEvent[2], fCentrality);

    if(wpSubEvent[0][0][0] < fMPar || wpSubEvent[1][0][0] < fMPar || wpSubEvent[2][0][0] < fMPar) return;
    PtSubContainer[rn]->FillThreeSubAnalsysis(wpSubEvent[0], wpSubEvent[1], wpSubEvent[2], fCentrality);



    PostData(1, fPtSampleList);
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::Terminate(Option_t *)
{

}
//_____________________________________________________________________________


void AliAnalysisTaskDensity::ProcessMCParticles()
{
    ClearWPCounter();       // clear and initiate dimension of wp counters

    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL) return;

    Double_t centrality(0);
    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fInputEvent->FindListObject("MultSelection"));
    
    if(!multSelection){ AliFatal("MC centrality not found!"); return; }
    centrality = multSelection->GetMultiplicityPercentile(fMultSelMethod);
    
    AliAODMCParticle* iParticle; 

    // loop over all primary MC particle
    for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
        iParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
        if (!iParticle) continue;
        if (!iParticle->IsPhysicalPrimary()) continue;

        double pt = iParticle->Pt();
        double eta = iParticle->Eta();
    } 

}
//_____________________________________________________________________________

void AliAnalysisTaskDensity::SetDefaultCut(){
    SetCorrelationOrder(6);
    SetFilterBit(96);
    SetMaxPileup(15000);
    SetMode("physics");
    SetMultSelectionMethod("V0M");
    SetPrimaryVertexZ(10.);
    SetPtRange(0.2, 3.0);
    SetTPCMinCls(70);
}


//_____________________________________________________________________________

bool AliAnalysisTaskDensity::AcceptAODTrack(AliAODTrack* track){
    // detector cut
    if(!track || !track->TestFilterBit(fFilterBit)) return kFALSE;         
    if(track->GetTPCNcls()>fTPCMinCls) return kFALSE;

    // physics cuts
    if(track->Pt() < fPtLow) return kFALSE;
    if(track->Pt() > fPtHigh) return kFALSE;
    if(std::abs(track->Eta()) > 0.8 ) return kFALSE;

    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskDensity::ClearWPCounter()
{
    wp.clear();
    wp.resize(fMPar+1,vector<Double_t>(fMPar+1));

    wpSubEvent.resize(3);
    for(int i(0); i < 3; i++){
        wpSubEvent[i].clear();
        wpSubEvent[i].resize(fMPar+1,vector<Double_t>(fMPar+1));
    }
}
void AliAnalysisTaskDensity::FillWPCounter(vector<vector<Double_t>> &inarr, double w, double p)
{
    for(int i=0;i<=fMPar;++i){
        for(int j=0;j<=fMPar;++j){
            inarr[i][j] += std::pow(w,i)*std::pow(p,j);
        }
    }
    return;
}
