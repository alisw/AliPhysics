/* AliAnaysisTasDensityCorrelation
 *
 *  Task for transverse momentum correlation with sub-events
 * 
*/
#include <cmath>
#include <vector>

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
#include "TH2D.h"

class AliAnalysisTaskDensity;
ClassImp(AliAnalysisTaskDensity);

AliAnalysisTaskDensity::AliAnalysisTaskDensity() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fPtSampleList(0),
    fInputListEfficiency(0),
    fQAEventList(0),
    fQATrackList(0),
    PtSubContainer(0),
    rndGenerator(0),
    fhDCAzDistribution(0),
    fhDCAxyDistribution(0),
    fhEtaDistribution(0),
    fhPtDistribution(0),
    fhCentSelected(0),
    fhPrimaryVzSelected(0),
    fhCrossedRowsTPC(0),
    fhChiPerTPCCls(0),
    fhChiPerITSCls(0),
    fhNofPileupSelected(0),
    fhSubEtaDistribution(0),
    fhMultSelected(0),
    fUseEffeciency(0),
    fSystFlag(0),
    fGFWSelection(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskDensity::AliAnalysisTaskDensity(const char* name, Bool_t bUseEff) : AliAnalysisTaskSE(name),
    fAOD(0), 
    fPtSampleList(0),
    fInputListEfficiency(0),
    fQAEventList(0),
    fQATrackList(0),
    PtSubContainer(0),
    rndGenerator(0),
    fhDCAzDistribution(0),
    fhDCAxyDistribution(0),
    fhEtaDistribution(0),
    fhPtDistribution(0),
    fhCentSelected(0),
    fhPrimaryVzSelected(0),
    fhCrossedRowsTPC(0),
    fhChiPerTPCCls(0),
    fhChiPerITSCls(0),
    fhNofPileupSelected(0),
    fhSubEtaDistribution(0),
    fhMultSelected(0),
    fUseEffeciency(bUseEff),
    fSystFlag(0),
    fGFWSelection(0)
{
    // Set cuts for default configuration
    SetDefaultSettings();
    // constructor
    DefineInput(0, TChain::Class());                // event chain
    if(bUseEff){DefineInput(1, TList::Class());};   // effeciency list    
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
    printf("\n\nCreating user output object!\n");
    if(!fGFWSelection) SetSystFlag(0);
    fGFWSelection->PrintSetup();
    fSystFlag = fGFWSelection->GetSystFlagIndex();
    if(fGFWSelection->GetSystFlagIndex() == 20) SetMultSelectionMethod("CL0");
    else if(fGFWSelection->GetSystFlagIndex() == 21) SetMultSelectionMethod("CL1");
    if(!fDCAxyFunctionalForm.IsNull()) { fGFWSelection->SetPtDepDCAXY(fDCAxyFunctionalForm); }
    
    // random number generator for bs samples
    rndGenerator = new TRandom();

    Int_t nbins = 100;
    Double_t xbinmax = 100.;
    PtSubContainer = new AliPtSubEventContainer*[10];
    for(int i(0); i < 10; i++){
        PtSubContainer[i] = new AliPtSubEventContainer(Form("ptcontsample%i",i),"ptcont", fMPar);
        PtSubContainer[i]->SetNamedPostfix(Form("_sample%i", i+1));
        PtSubContainer[i]->Initialize(nbins, 0., xbinmax); 
        PtSubContainer[i]->InitializeTwoSub(nbins, 0., xbinmax); 
        PtSubContainer[i]->InitializeThreeSub(nbins, 0., xbinmax);
    }
    fPtSampleList = new TList(); fPtSampleList->SetOwner(kTRUE);     
    for(int i(0); i < 10; i++){
        fPtSampleList->Add( new TList() ); ((TList*)fPtSampleList->At(i))->SetName(Form("sample%i", i+1));
        ((TList*)fPtSampleList->At(i))->Add( PtSubContainer[i]->GetTwoSubAnalysisList() );
        ((TList*)fPtSampleList->At(i))->Add( PtSubContainer[i]->GetThreeSubAnalysisList() );
        ((TList*)fPtSampleList->At(i))->Add( PtSubContainer[i]->GetCorrList() );
        ((TList*)((TList*)fPtSampleList->At(i))->At(2))->SetName("FullTPCCorrelation");
    }

    if(fUseEffeciency) {
        printf("Getting effeciency input list!\n");
        fInputListEfficiency = (TList*)GetInputData(1);
        if(!fInputListEfficiency) {AliFatal("Efficiency input list not loaded"); return;}
        if(fPtLow < 0.2 || 3.0 < fPtHigh) AliWarning("Efficiency loading -- pt cut can be out of range!");
        if(fSystFlag==0){
            fhEffeciencyPt = (TH1D*)fInputListEfficiency->FindObject("EffRescaled_Cent0");
        }
        else{
            fhEffeciencyPt = (TH1D*)fInputListEfficiency->FindObject(Form("EffRescaled_Cent0_SystFlag%i_", fSystFlag));
        }
        if(!fhEffeciencyPt){AliError(Form("Efficiency for pt not loaded with SystFlag: %i", fSystFlag)); return; }
    }
    printf("Pt correlation objects created!\n");
    PostData(1, fPtSampleList);
    
    printf("Creating QA event objects\n");
    fQAEventList = new TList(); fQAEventList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAEventList,kTRUE);
    fhNofPileupSelected = new TH1D("nofpileup",                 "nofpileup",                200, 0, 20000);
    fhPrimaryVzSelected = new TH1D("vtxz distribution",         "vtxz distribution",        300, -15., 15);
    fhCentSelected      = new TH2D("centrality distribution",   "centrality distribution",  200, 0., 100., 3., 0, 3.);
    fhMultSelected  = new TH2D("multiplicity distribution",     "multiplicity distribution",300, 0., 12000., 3., 0, 3.);
    const char* ntriggerbins[3] = {"Minbias", "Central", "Semi-Central"};

    for (int i(1);i<=3;i++) { fhCentSelected->GetYaxis()->SetBinLabel(i,ntriggerbins[i-1]); }
    for (int i(1);i<=3;i++) { fhMultSelected->GetYaxis()->SetBinLabel(i,ntriggerbins[i-1]); }

    fQAEventList->Add(fhNofPileupSelected);
    fQAEventList->Add(fhPrimaryVzSelected);
    fQAEventList->Add(fhCentSelected);
    fQAEventList->Add(fhMultSelected);
    PostData(2,fQAEventList);
    printf("QA Event objects created!\n");

    printf("Creating QA track objects\n");
    fQATrackList = new TList(); fQATrackList->SetOwner(kTRUE);
    fhCrossedRowsTPC    = new TH1D("CrossedTPC",         "crossed rows in TPC", 200, 0, 200);
    fhChiPerTPCCls      = new TH1D("ChiPerTPC",          "chi2 / tpc clusters", 200, 0, 5);
    fhChiPerITSCls      = new TH1D("ChiPerITS",          "chi2 / its clusters", 500, 0, 50);
    fhDCAzDistribution  = new TH1D("DCAz distribution",  "DCAz distribution",   600, -3, 3);
    fhDCAxyDistribution = new TH1D("DCAxy distribution", "DCAxy distribution",  400, -2, 2);
    fhEtaDistribution   = new TH1D("eta distribution",   "eta distribution",    200, -1, 1);
    fhPtDistribution    = new TH2D("pt distributoon",    "pt distribution",     600, 0., 10., 3., 0, 3.);
    fhSubEtaDistribution = new TH2D("sub eta distributoon",    "sub eta distributoon", 200, -1, 1, 3., 0, 3.);

    const char* nptbins[3] = {"Full", "Trig", "TrigWeighted"};
    const char* netabins[3] = {"Full", "Two-Sub", "Three-Sub"};
    for (int i(1);i<=3;i++) { fhPtDistribution->GetYaxis()->SetBinLabel(i,nptbins[i-1]); }
    for (int i(1);i<=3;i++) { fhSubEtaDistribution->GetYaxis()->SetBinLabel(i,netabins[i-1]); }
    fQATrackList->Add(fhCrossedRowsTPC);
    fQATrackList->Add(fhChiPerTPCCls);
    fQATrackList->Add(fhChiPerITSCls);
    fQATrackList->Add(fhDCAzDistribution);
    fQATrackList->Add(fhDCAxyDistribution);
    fQATrackList->Add(fhEtaDistribution);
    fQATrackList->Add(fhPtDistribution);
    fQATrackList->Add(fhSubEtaDistribution);
    PostData(3,fQATrackList);
    printf("QA Tracks objects created!\n\n");
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::UserExec(Option_t *)
{
    ClearWPCounter();   // restart and initiate containr for <m> calcualtions
    if(fMode == "HIJING"){
        fMCEvent = MCEvent();
        if(fMCEvent) ProcessMCParticles();
        return;
    }

    // get an AOD event from input file
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());   
    if(!fAOD) return;                                  

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man){ AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());}
    
    AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(!multSelection) { AliFatal("AOD Centrality not found!"); return; }
    fCentrality = multSelection->GetMultiplicityPercentile(fMultSelMethod);
    
    if(!CheckTrigger(fCentrality)) return;
    if(!AcceptAODEvent(fAOD)) return;
    Int_t NchSelected=0;
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) 
    {                 
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!AcceptAODTrack(track)) continue;
        ProcessTrack( WeightPt(track->Pt()), track->Pt(), track->Eta());
        NchSelected++;
    }                      
    FillEventQA(fCentrality, NchSelected);       
    ProcessEventCorrelation(rndGenerator->Integer(10));
    PostData(1, fPtSampleList);
}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::Terminate(Option_t *)
{

}
//_____________________________________________________________________________
void AliAnalysisTaskDensity::ProcessTrack(Double_t lweight, Double_t lpt, Double_t leta)
{
    FillWPCounter(wp, lweight, lpt);                                                            // regular [-0.8 < eta < 0.8]
    if( (-fAbsEta<leta) && (leta<-fEtaGap))  FillWPCounter(wpTwoSubEvent[0], lweight, lpt);     // two sub [-0.8 < eta < -0.4]
    if( (fEtaGap<leta)  && (leta<fAbsEta))   FillWPCounter(wpTwoSubEvent[1], lweight, lpt);     // two sub [0.4 < eta < 0.8 ]
    if( -fAbsEta<leta && leta<-0.4) FillWPCounter(wpThreeSubEvent[0], lweight, lpt);                        // three sub [-0.8 < eta < -0.4]
    if(std::abs(leta) < 0.2)        FillWPCounter(wpThreeSubEvent[1], lweight, lpt);                        // three sub [-0.2 < eta < 0.2]
    if( 0.4<leta && leta<fAbsEta)   FillWPCounter(wpThreeSubEvent[2], lweight, lpt);                        // three sub [0.4 < eta < 0.8 ]
    // Fill sub-event
    if( std::abs(leta)<fAbsEta) fhSubEtaDistribution->Fill(leta, "Full", 1);
    if( (-fAbsEta<leta) && (leta<-fEtaGap) || (fEtaGap<leta)  && (leta<fAbsEta)) fhSubEtaDistribution->Fill(leta, "Two-Sub", 1);
    if( (-fAbsEta<leta) && (leta<-0.4) || (std::abs(leta)<0.2) || (0.4<leta) && (leta<fAbsEta)) fhSubEtaDistribution->Fill(leta, "Three-Sub", 1);

}
void AliAnalysisTaskDensity::ProcessEventCorrelation(Int_t id)
{
    if(wp[0][0] < fMPar) return;
    PtSubContainer[id]->FillRecursive(wp,fCentrality);
    if(wpTwoSubEvent[0][0][0] < fMPar || wpTwoSubEvent[1][0][0] < fMPar) return;
    PtSubContainer[id]->FillTwoSubAnalsysis(wpTwoSubEvent[0], wpTwoSubEvent[1], fCentrality);
    if(wpThreeSubEvent[0][0][0] < fMPar || wpThreeSubEvent[1][0][0] < fMPar || wpThreeSubEvent[2][0][0] < fMPar) return;
    PtSubContainer[id]->FillThreeSubAnalsysis(wpThreeSubEvent[0], wpThreeSubEvent[1], wpThreeSubEvent[2], fCentrality);
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskDensity::ProcessMCParticles()
{
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
        if (!AcceptMCTrack(MCTrack)) continue;            // track cuts
        ProcessTrack( WeightPt(MCTrack->Pt()), MCTrack->Pt(), MCTrack->Eta());
    } 
    ProcessEventCorrelation(rndGenerator->Integer(10));
    PostData(1, fPtSampleList);
}
// --------------------------------------------------------------------------------------------
void AliAnalysisTaskDensity::SetDefaultSettings(){
    SetCorrelationOrder(8);
    SetEtaGap(0.1);         // gap to zero -> gap 0.1 -> |deta|=0.2
    SetMode("physics");
    SetMultSelectionMethod("V0M");
    // physics selection
    SetPtRange(0.2, 3.0);
    SetAbsEta(0.8);
    // detector selection
    SetMaxPileup(15000);    //  
}

//_____________________________________________________________________________
bool AliAnalysisTaskDensity::AcceptAODEvent(AliAODEvent* event){
    if(!fEventCuts.AcceptEvent(event)) return kFALSE;
    if(!fGFWSelection->AcceptVertex(event)) return kFALSE;
    // all triggers passed -> write to hist
    fhNofPileupSelected->Fill( event->GetNumberOfPileupVerticesTracks() );
    fhPrimaryVzSelected->Fill( event->GetPrimaryVertex()->GetZ() );
    PostData(2,fQAEventList);
    return kTRUE;
}

bool AliAnalysisTaskDensity::CheckTrigger(Double_t lCent){
    UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(!fSelMask) { return kFALSE; }; 
    if(fSelMask&AliVEvent::kINT7) { fhCentSelected->Fill(lCent, "Minbias", 1); return kTRUE; }; 
    if((fSelMask&AliVEvent::kCentral) && lCent>10) {return kFALSE; }; 
    if((fSelMask&AliVEvent::kSemiCentral) && (lCent<30 || lCent>50)) {return kFALSE; };
    return kTRUE;
}


void AliAnalysisTaskDensity::FillEventQA(Double_t lCent, Int_t lNch)
{
    UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(fSelMask&AliVEvent::kINT7){
        fhCentSelected->Fill(lCent, "Minbias", 1);
        fhMultSelected->Fill(lNch,  "Minbias", 1);
    }
    if(fSelMask&AliVEvent::kCentral){
        fhCentSelected->Fill(lCent, "Central", 1);
        fhMultSelected->Fill(lNch,  "Central", 1);
    }
    if(fSelMask&AliVEvent::kSemiCentral){
        fhCentSelected->Fill(lCent, "Semi-Central", 1);
        fhMultSelected->Fill(lNch,  "Semi-Central", 1);
    }
}


//_____________________________________________________________________________
bool AliAnalysisTaskDensity::AcceptAODTrack(AliAODTrack* track)
{
    if(!track || !track->TestFilterBit( fGFWSelection->fFilterBit )) return kFALSE;         
    fhPtDistribution->Fill(track->Pt(), "Full", 1);
    if(track->Pt() < fPtLow) return kFALSE;
    if(track->Pt() > fPtHigh) return kFALSE;
    if(std::abs(track->Eta()) > 0.8 ) return kFALSE;

    const AliAODVertex* vtxp = dynamic_cast<const AliAODVertex*>(fAOD->GetPrimaryVertex());
    Double_t ltrackXYZ[] = {0.,0.,0.};
    track->GetXYZ(ltrackXYZ);
    if(vtxp) {
        ltrackXYZ[0] = ltrackXYZ[0]-vtxp->GetX();
        ltrackXYZ[1] = ltrackXYZ[1]-vtxp->GetY();
        ltrackXYZ[2] = ltrackXYZ[2]-vtxp->GetZ();
    }
    // check if track passes trigger
    if(!fGFWSelection->AcceptTrack(track,fSystFlag==1?0:ltrackXYZ,0,kFALSE)) return kFALSE;
    // physics selection
    fhEtaDistribution->Fill(track->Eta());
    fhPtDistribution->Fill(track->Pt(), "Trig", 1);
    fhPtDistribution->Fill(track->Pt(), "TrigWeighted", WeightPt(track->Pt()));
    // detector selection
    fhDCAxyDistribution->Fill(ltrackXYZ[0]+ltrackXYZ[1]);
    fhDCAzDistribution->Fill(ltrackXYZ[2]);
    fhChiPerTPCCls->Fill(track->GetTPCchi2perCluster());
    fhChiPerITSCls->Fill(track->GetITSchi2()/(double)track->GetITSNcls());
    fhCrossedRowsTPC->Fill(track->GetTPCNclsF());

    PostData(3,fQATrackList);
    return kTRUE;
}
bool AliAnalysisTaskDensity::AcceptMCTrack(AliAODMCParticle* track)
{
    // physics cuts
    if(track->Pt() < fPtLow) return kFALSE;
    if(track->Pt() > fPtHigh) return kFALSE;
    if(std::abs(track->Eta()) > fAbsEta ) return kFALSE;
    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskDensity::ClearWPCounter()
{
    wp.clear(); wp.resize(fMPar+1,std::vector<Double_t>(fMPar+1));
    wpTwoSubEvent.resize(2);
    for(int i(0); i < 2; i++){
        wpTwoSubEvent[i].clear();
        wpTwoSubEvent[i].resize(fMPar+1,std::vector<Double_t>(fMPar+1));
    }
    wpThreeSubEvent.resize(3);
    for(int i(0); i < 3; i++){
        wpThreeSubEvent[i].clear();
        wpThreeSubEvent[i].resize(fMPar+1,std::vector<Double_t>(fMPar+1));
    }
}
void AliAnalysisTaskDensity::FillWPCounter(std::vector<std::vector<Double_t>> &inarr, Double_t w, Double_t p)
{
    for(int i=0;i<=fMPar;++i){
        for(int j=0;j<=fMPar;++j){
            inarr[i][j] += std::pow(w,i)*std::pow(p,j);
        }
    }
    return;
}
Double_t AliAnalysisTaskDensity::WeightPt(Double_t lpt){
    if(!fUseEffeciency) return 1.;
    if(fhEffeciencyPt->GetBinContent(fhEffeciencyPt->FindBin(lpt))==0) return 1.;   // do not apply weight if outside scope
    return 1./(fhEffeciencyPt->GetBinContent( fhEffeciencyPt->FindBin(lpt) ));
}