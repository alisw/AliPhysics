#ifndef ALIOADBMULTSELECTION_H
#define ALIOADBMULTSELECTION_H

#include <iostream>
#include "TNamed.h"
#include "TH1F.h"
#include "TList.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"
#include "TFolder.h"
#include "TFolder.h"
#include "TObjString.h"
#include "TBrowser.h"

using namespace std;

class TH1F;
class TList; 
class AliMultSelectionCuts;
class AliMultSelection;

class AliOADBMultSelection : public TNamed {
    
public:
    AliOADBMultSelection();
    AliOADBMultSelection(const char * name, const char * title = "OADB for multiplicity and centrality selection");
    ~AliOADBMultSelection();
    
    //Getters
    Long_t GetNEstimators () { return fSelection->GetNEstimators(); }
    
    TH1F *GetCalibHisto(Long_t iEst)
    {
        return (TH1F*) fCalibList->At(iEst);
    }
    TH1F *GetCalibHisto(TString lCalibHistoName)
    {
        return ((TH1F*)fCalibList->FindObject(lCalibHistoName)); 
    }
    
    void AddCalibHisto (TH1F * var) {
        fCalibList->Add(var);
    }
    AliMultSelectionCuts * GetEventCuts()                            { return fEventCuts;     }
    void                   SetEventCuts (AliMultSelectionCuts * var) { fEventCuts = var;      }
    AliMultSelection     * GetMultSelection()                        { return fSelection;     }
    void                   SetMultSelection(AliMultSelection * var ) { fSelection = var;      }
    
    
    // Make it browsable
    void Browse(TBrowser *b);
    virtual Bool_t IsFolder() const { return kTRUE; }
    
private:
    TList *fCalibList; // Calibration Histograms
    AliMultSelectionCuts * fEventCuts; // EventCuts
    AliMultSelection     * fSelection; // Definition of Estimators
    
    ClassDef(AliOADBMultSelection, 1)
    
    
};

#endif
