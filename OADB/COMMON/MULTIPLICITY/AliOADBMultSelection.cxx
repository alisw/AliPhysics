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


ClassImp(AliOADBMultSelection);
AliOADBMultSelection::AliOADBMultSelection() :
TNamed("multSel",""), fEventCuts(0), fSelection(0), fCalibList(0)
{
    // constructor
    fCalibList = new TList();
    fCalibList -> SetOwner (kTRUE) ;
}

AliOADBMultSelection::AliOADBMultSelection(const char * name, const char * title) :
TNamed(name, title), fEventCuts(0), fSelection(0), fCalibList(0)
{
    // constructor
    fCalibList = new TList();
    fCalibList -> SetOwner (kTRUE) ;
}

AliOADBMultSelection::~AliOADBMultSelection(){
    // Destructor
    if(fEventCuts)     delete fEventCuts;
    if(fSelection)     delete fSelection;
    
    if( fCalibList) {
        fCalibList -> Delete();
    }
}

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
        cutsFolder->Add(new TObjString(Form("Vz: %f", GetEventCuts()->GetVzCut())));
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
            estFolder->Add(lEst[iEst]);
        }
        b->Add(histoFolder);
        b->Add(cutsFolder);
        b->Add(estFolder);
        
    }
    else
        TObject::Browse(b);
}
