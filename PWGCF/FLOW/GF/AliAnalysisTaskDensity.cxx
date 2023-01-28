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
#include "AliEventCuts.h"
#include "AliAODVertex.h"


// monte carlo events
#include "AliAODMCParticle.h"

// plotting 
#include "TH1D.h"


class AliAnalysisTaskDensity;
ClassImp(AliAnalysisTaskDensity);

AliAnalysisTaskDensity::AliAnalysisTaskDensity() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fPtSampleList(0),
    fQAList(0),
    fQATrackList(0),
    PtSubContainer(0),
    rndGenerator(0),
    fhDCAzDistribution(0),
    fhDCAxyDistribution(0),
    fhEtaDistribution(0),
    fhPtDistribution(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskDensity::AliAnalysisTaskDensity(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), 
    fPtSampleList(0),
    fQAList(0),
    fQATrackList(0),
    PtSubContainer(0),
    rndGenerator(0),
    fhDCAzDistribution(0),
    fhDCAxyDistribution(0),
    fhEtaDistribution(0),
    fhPtDistribution(0)
{
    // Set cuts for default configuration
    SetDefaultSettings();

    // constructor
    DefineInput(0, TChain::Class());   
    DefineOutput(1, TList::Class());    
    DefineOutput(2, TList::Class());   
    DefineOutput(3, TList::Class()); 
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
    printf("Pt correlation objects created!\n");
    PostData(1, fPtSampleList);
    
    printf("Creating QA event objects\n");
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAList,kTRUE);
    PostData(2,fQAList);
    printf("QA AliEventCuts objects created!\n");

    printf("Creating QA track objects\n");
    fQATrackList = new TList();
    fQATrackList->SetOwner(kTRUE);
    fhDCAzDistribution = new TH1D("DCAz distribution", "DCAz distribution", 600, -3, 3);
    fhDCAxyDistribution = new TH1D("DCAxy distribution", "DCAxy distribution", 400, -2, 2);
    fhEtaDistribution = new TH1D("eta distribution", "eta distribution", 200, -1, 1);
    fhPtDistribution = new TH1D("pt distributoon", "pt distribution", 500, 0., 5.);

    fQATrackList->Add(fhDCAzDistribution);
    fQATrackList->Add(fhDCAxyDistribution);
    fQATrackList->Add(fhEtaDistribution);
    fQATrackList->Add(fhPtDistribution);
    PostData(3,fQATrackList);
    printf("QA Tracks objects created!\n");

}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::UserExec(Option_t *)
{
    // MC Event
    if(fMode == "HIJING"){
        fMCEvent = MCEvent();
        if(fMCEvent) ProcessMCParticles();
        return;
    }

    // get an AOD event from input file
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());   
    if(!fAOD) return;                                  

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    }

    if(!AcceptAODEvent(fAOD)) return;

    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(!multSelection) { AliFatal("AOD Centrality not found!"); return; }
    fCentrality = multSelection->GetMultiplicityPercentile(fMultSelMethod);

    ClearWPCounter();   // restart and initiate containr for <m> calcualtions
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) 
    {                 
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!AcceptAODTrack(track)) continue;             // track cuts
        
        FillWPCounter(wp, 1., track->Pt());
        
        // Fill pt-correlation for sub-events  
        if( -0.8 >  track->Eta() > -0.4) FillWPCounter(wpSubEvent[0], 1., track->Pt());     // -0.8 < eta < -0.4
        if(std::abs(track->Eta()) < 0.2)  FillWPCounter(wpSubEvent[1], 1., track->Pt());    // -0.2 < eta < 0.2
        if( 0.4 <   track->Eta() < 0.8)  FillWPCounter(wpSubEvent[2], 1., track->Pt());     //  0.4 < eta < 0.8 
    }                             
    ProcessEventCorrelation(rndGenerator->Integer(10));
    PostData(1, fPtSampleList);
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::Terminate(Option_t *)
{

}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::ProcessEventCorrelation(Int_t id){
    if(wp[0][0] < fMPar) return;
    PtSubContainer[id]->FillRecursive(wp,fCentrality);
    
    if(wpSubEvent[0][0][0] < fMPar || wpSubEvent[2][0][0] < fMPar) return;
    PtSubContainer[id]->FillTwoSubAnalsysis(wpSubEvent[0], wpSubEvent[2], fCentrality);

    if(wpSubEvent[0][0][0] < fMPar || wpSubEvent[1][0][0] < fMPar || wpSubEvent[2][0][0] < fMPar) return;
    PtSubContainer[id]->FillThreeSubAnalsysis(wpSubEvent[0], wpSubEvent[1], wpSubEvent[2], fCentrality);
}




void AliAnalysisTaskDensity::ProcessMCParticles()
{
    ClearWPCounter();       // clear and initiate dimension of wp counters

    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL) return;

    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fInputEvent->FindListObject("MultSelection"));    
    if(!multSelection){ AliFatal("MC centrality not found!"); return; }
    fCentrality = multSelection->GetMultiplicityPercentile(fMultSelMethod);
    
    AliAODMCParticle* MCTrack; 

    for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
        MCTrack = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
        if (!MCTrack) continue;
        if (!MCTrack->IsPhysicalPrimary()) continue;      // loop over all primary MC particle only

        if(!AcceptMCTrack(MCTrack)) continue;             // track cuts
        FillWPCounter(wp, 1., MCTrack->Pt());

        // Fill pt-correlation for sub-events  
        if( -0.8 >  MCTrack->Eta() > -0.4) FillWPCounter(wpSubEvent[0], 1., MCTrack->Pt());     // -0.8 < eta < -0.4
        if(std::abs(MCTrack->Eta()) < 0.2)  FillWPCounter(wpSubEvent[1], 1., MCTrack->Pt());    // -0.2 < eta < 0.2
        if( 0.4 <   MCTrack->Eta() < 0.8)  FillWPCounter(wpSubEvent[2], 1., MCTrack->Pt());     //  0.4 < eta < 0.8 
    } 
    ProcessEventCorrelation(rndGenerator->Integer(10));
    PostData(1, fPtSampleList);
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::SetDefaultSettings(){
    SetCorrelationOrder(6);
    SetFilterBit(96);
    SetMode("physics");
    SetMultSelectionMethod("V0M");
    // physics selection
    SetPtRange(0.2, 3.0);
    // +track selection cut
    SetTPCMinCls(70);
    // +event selection cut
    SetMaxPileup(15000);
    //SetMaxPrimaryVz(10.);
}
void AliAnalysisTaskDensity::SetTightTrackCuts(UInt_t filterbit, UShort_t ncls, Double_t dcaz){
    SetFilterBit(filterbit);
    SetTPCMinCls(ncls);
    SetMaxDCAz(dcaz);
}
void AliAnalysisTaskDensity::SetTightEventCuts(Double_t zvertex, Int_t pileup){
    SetMaxPrimaryVz(zvertex);
    SetMaxPileup(pileup);
}



//_____________________________________________________________________________
bool AliAnalysisTaskDensity::AcceptAODEvent(AliAODEvent* event){
    if(!fEventCuts.AcceptEvent(event)) return kFALSE;
    //if(fMaxPileup < event->GetNumberOfPileupVerticesTracks()) return kFALSE;
    //if(fMaxPrimaryVertexZ < TMath::Abs(event->GetPrimaryVertex()->GetZ())) return kFALSE;
    return kTRUE;
}

//_____________________________________________________________________________
bool AliAnalysisTaskDensity::AcceptAODTrack(AliAODTrack* track){
    // detector cut
    if(!track || !track->TestFilterBit(fFilterBit)) return kFALSE;         
    if(track->GetTPCNcls()<fTPCMinCls) return kFALSE;
    // physics cuts
    if(track->Pt() < fPtLow) return kFALSE;
    if(track->Pt() > fPtHigh) return kFALSE;
    if(std::abs(track->Eta()) > 0.8 ) return kFALSE;
    FillTrackSelection(track);
    return kTRUE;
}
bool AliAnalysisTaskDensity::AcceptMCTrack(AliAODMCParticle* track){
    // physics cuts
    if(track->Pt() < fPtLow) return kFALSE;
    if(track->Pt() > fPtHigh) return kFALSE;
    if(std::abs(track->Eta()) > 0.8 ) return kFALSE;
    return kTRUE;
}


void AliAnalysisTaskDensity::FillTrackSelection(AliAODTrack* track){
    fhEtaDistribution->Fill(track->Eta());
    fhPtDistribution->Fill(track->Pt());

    //Propegate track to DCA
    const AliAODVertex* vtxp = dynamic_cast<const AliAODVertex*>(fAOD->GetPrimaryVertex());
    Double_t ltrackXYZ[] = {0.,0.,0.};
    track->GetXYZ(ltrackXYZ);
    if(ltrackXYZ && vtxp) {
        ltrackXYZ[0] = ltrackXYZ[0]-vtxp->GetX();
        ltrackXYZ[1] = ltrackXYZ[1]-vtxp->GetY();
        ltrackXYZ[2] = ltrackXYZ[2]-vtxp->GetZ();
        fhDCAxyDistribution->Fill(ltrackXYZ[0]+ltrackXYZ[1]);
        fhDCAzDistribution->Fill(ltrackXYZ[2]);
    }
    PostData(3,fQATrackList);
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
