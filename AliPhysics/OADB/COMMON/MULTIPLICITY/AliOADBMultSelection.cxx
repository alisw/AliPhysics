#include "AliOADBMultSelection.h"
#include "TH1F.h"
#include "TList.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"
#include "TFolder.h"
#include "TObjString.h"
#include "TBrowser.h"
#include <TMap.h>
#include <TROOT.h>

ClassImp(AliOADBMultSelection);

//________________________________________________________________
//Constructors/Destructor
AliOADBMultSelection::AliOADBMultSelection() :
TNamed("multSel",""), fCalibList(0), fEventCuts(0), fSelection(0), fMap(0)
{
    // constructor
    // fCalibList = new TList();
    // fCalibList -> SetOwner (kTRUE) ;
}
//________________________________________________________________
AliOADBMultSelection::AliOADBMultSelection(const AliOADBMultSelection& o)
: TNamed(o),
fCalibList(0),
fEventCuts(0),
fSelection(0),
fMap(0)
{
    fCalibList = new TList();
    fCalibList->SetOwner (kTRUE);
    TIter next(o.fCalibList);
    TObject* obj = 0;
    while ((obj = next())) fCalibList->Add(obj->Clone());
    fSelection = new AliMultSelection(*o.fSelection);
    fEventCuts = new AliMultSelectionCuts(*o.fEventCuts);
}
//________________________________________________________________
AliOADBMultSelection::AliOADBMultSelection(const char * name, const char * title) :
TNamed(name, title), fCalibList(0), fEventCuts(0), fSelection(0), fMap(0)
{
    // constructor
    fCalibList = new TList();
    fCalibList -> SetOwner (kTRUE) ;
}
//________________________________________________________________
AliOADBMultSelection& AliOADBMultSelection::operator=(const AliOADBMultSelection& o)
{
    if (&o == this) return *this;
    SetName(o.GetName());
    SetTitle(o.GetTitle());
    if (fCalibList) {
        fCalibList->Delete();
        fCalibList = 0;
    }
    if (fMap) {
        delete fMap;
        fMap = 0;
    }
    fCalibList = new TList();
    fCalibList->SetOwner (kTRUE);
    TIter next(o.fCalibList);
    TObject* obj = 0;
    while ((obj = next())) fCalibList->Add(obj->Clone());
    SetMultSelection(new AliMultSelection(*o.fSelection));
    SetEventCuts(new AliMultSelectionCuts(*o.fEventCuts));
    return *this;
}
//________________________________________________________________
AliOADBMultSelection::~AliOADBMultSelection(){
    // Destructor
    if(fEventCuts)     delete fEventCuts;
    if(fSelection)     delete fSelection;
    
    //if( fCalibList) {
    //    fCalibList -> Delete();
    //}
}

//________________________________________________________________
void AliOADBMultSelection::Dissociate()
{
    TIter next(fCalibList);
    TH1*  hist = 0;
    while ((hist = static_cast<TH1*>(next()))) hist->SetDirectory(0);
}
//________________________________________________________________
void AliOADBMultSelection::AddCalibHisto (TH1F * var) {
    if (!var) return;
    if (!fCalibList) {
        //create if needed...
        fCalibList = new TList;
        fCalibList->SetOwner();
    }
    fCalibList->Add(var);
}
//________________________________________________________________
void AliOADBMultSelection::SetEventCuts(AliMultSelectionCuts* c)
{
    if (fEventCuts) delete fEventCuts;
    fEventCuts = c;
}
//________________________________________________________________
void AliOADBMultSelection::SetMultSelection(AliMultSelection* s)
{
    if (fSelection) delete fSelection;
    fSelection = s;
}
//________________________________________________________________
TH1F* AliOADBMultSelection::GetCalibHisto(Long_t iEst) const
{
    if (!fCalibList) return 0;
    if (iEst < 0 || iEst >= fCalibList->GetEntries()) return 0;
    return (TH1F*) fCalibList->At(iEst);
}
//________________________________________________________________
TH1F* AliOADBMultSelection::GetCalibHisto(const TString& lCalibHistoName) const
{
    if (!fCalibList) return 0;
    return ((TH1F*)fCalibList->FindObject(lCalibHistoName));
}
//________________________________________________________________
void AliOADBMultSelection::Print(Option_t* option) const
{
    Printf("%s: %s/%s", ClassName(), GetName(), GetTitle());
    gROOT->IncreaseDirLevel();
    if (fSelection) {
        gROOT->IndentLevel();
        fSelection->Print(option);
    }
    if (fEventCuts) {
        gROOT->IndentLevel();
        fEventCuts->Print(option);
    }
    if (fMap) {
        gROOT->IndentLevel();
        Printf("Mapping from estimator to histogram");
        gROOT->IncreaseDirLevel();
        TIter  next(fMap);
        TObject* key = 0;
        while ((key = static_cast<TPair*>(next()))) {
            TObject* value = fMap->GetValue(key);
            TString  k(key->GetName());
            TString  h(value ? value->GetName() : "?");
            gROOT->IndentLevel();
            Printf(" Estimator %s -> %s", k.Data(), h.Data());
        }
        gROOT->DecreaseDirLevel();
    }
    else if (fCalibList) {
        gROOT->IndentLevel();
        Printf("Calibration histograms");
        gROOT->IncreaseDirLevel();
        TIter next(fCalibList);
        TObject* o = 0;
        while ((o = next())) {
            gROOT->IndentLevel();
            o->Print(option);
        }
        gROOT->DecreaseDirLevel();  
    }
    
}
//________________________________________________________________
void AliOADBMultSelection::Browse(TBrowser *b)
{
    // Browse this object.
    // If b=0, there is no Browser: call TObject::Browse(0) instead.
    //         This means TObject::Inspect() will be invoked indirectly
    // FIXME: this should be implemented for all histograms
    // FIXME: who deletes the folders? Make sure there are no memory leaks here!
    if (b) {
        TFolder  * histoFolder = new TFolder ("histograms", "Calibration Histograms");
        TFolder  * cutsFolder  = new TFolder ("cuts"      , "Event Selection Cuts");
        TFolder  * estFolder   = new TFolder ("estimators", "Estimators");
        histoFolder->SetOwner();
        cutsFolder->SetOwner();
        estFolder->SetOwner();
        for(Int_t iEst=0; iEst<fCalibList->GetEntries(); iEst++) histoFolder->Add(GetCalibHisto(iEst));
        
        //Info for Event Selection: Printouts 
        cutsFolder->Add(new TObjString(Form("Vertex Z cut: %f", GetEventCuts()->GetVzCut())));
        cutsFolder->Add(new TObjString(Form("Physics Selection: %i", GetEventCuts()->GetTriggerCut())));
        cutsFolder->Add(new TObjString(Form("INEL > 0: %i", GetEventCuts()->GetINELgtZEROCut())));
        cutsFolder->Add(new TObjString(Form("Tracklets Vs Clusters: %i", GetEventCuts()->GetTrackletsVsClustersCut())));
        cutsFolder->Add(new TObjString(Form("Reject Pileup SPD (mult bins): %i", GetEventCuts()->GetRejectPileupInMultBinsCut())));
        cutsFolder->Add(new TObjString(Form("SPD and Track Vertex Consistency: %i", GetEventCuts()->GetVertexConsistencyCut())));
        cutsFolder->Add(new TObjString(Form("NContribs to PV > 0: %i", GetEventCuts()->GetNonZeroNContribs())));
        
        const long lNEstimators = fSelection->GetNEstimators();
        TFolder *lEst[lNEstimators];
        for(Int_t iEst=0; iEst<lNEstimators; iEst++){
            lEst[iEst] = new TFolder (Form("%s_Info",fSelection->GetEstimator(iEst)->GetName()), "Info");
            lEst[iEst]->SetOwner();
            lEst[iEst]->Add(new TObjString(Form("Estimator Name: %s", fSelection->GetEstimator(iEst)->GetName() ) ) );
            TString lThisDef = fSelection->GetEstimator(iEst)->GetDefinition();
            lEst[iEst]->Add(new TObjString(Form("Estimator Definition: %s", lThisDef.Data() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Is an Integer: %i", fSelection->GetEstimator(iEst)->IsInteger() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Current Value: %f", fSelection->GetEstimator(iEst)->GetValue() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Current Percentile: %f", fSelection->GetEstimator(iEst)->GetPercentile() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Mean Value: %f", fSelection->GetEstimator(iEst)->GetMean() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Is Anchored: %i", fSelection->GetEstimator(iEst)->GetUseAnchor() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Anchor Point (raw): %f", fSelection->GetEstimator(iEst)->GetAnchorPoint() ) ) );
            lEst[iEst]->Add(new TObjString(Form("Anchor Percentile: %f", fSelection->GetEstimator(iEst)->GetAnchorPercentile() ) ) );
            estFolder->Add(lEst[iEst]);
        }
        b->Add(histoFolder);
        b->Add(cutsFolder);
        b->Add(estFolder);
        
    }
    else
        TObject::Browse(b);
}
//________________________________________________________________
TH1F* AliOADBMultSelection::FindHisto(AliMultEstimator* e)
{
    if (!fMap) return 0;
    TPair* ret = static_cast<TPair*>(fMap->FindObject(e));
    if (!ret) return 0;
    return static_cast<TH1F*>(ret->Value());
}
//________________________________________________________________
void AliOADBMultSelection::Setup()
{
    if (fMap) {
        delete fMap;
        fMap = 0;
    }
    AliMultSelection* sel = GetMultSelection();
    if (!sel) return;
    
    fMap = new TMap;
    fMap->SetOwner(false);
    
    for(Long_t iEst=0; iEst<sel->GetNEstimators(); iEst++) {
        AliMultEstimator* e = sel->GetEstimator(iEst);
        if (!e) continue;
        
        TString name(Form("hCalib_%s", e->GetName()));
        TH1F*   h = GetCalibHisto(name);
        if (!h) continue;
        
        fMap->Add(e, h);
    }
}


