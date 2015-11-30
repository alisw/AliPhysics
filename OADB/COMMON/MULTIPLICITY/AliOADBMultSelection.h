#ifndef ALIOADBMULTSELECTION_H
#define ALIOADBMULTSELECTION_H

#include <TNamed.h>
#include <AliMultSelection.h>
class TBrowser;
class TH1F;
class TList; 
class AliMultSelectionCuts;
class AliMultEstimator;
class TMap;

class AliOADBMultSelection : public TNamed {
    
public:
    AliOADBMultSelection();
    AliOADBMultSelection(const AliOADBMultSelection& o);
    AliOADBMultSelection(const char * name, const char * title = "OADB for multiplicity and centrality selection");
    AliOADBMultSelection& operator=(const AliOADBMultSelection& o);
    ~AliOADBMultSelection();
    
    //Getters
    Long_t GetNEstimators () { return fSelection ? fSelection->GetNEstimators() : 0; }
    
    TH1F *GetCalibHisto(Long_t iEst) const;
    TH1F *GetCalibHisto(const TString& lCalibHistoName) const;
    
    void AddCalibHisto (TH1F * var);
    AliMultSelectionCuts * GetEventCuts() const                           { return fEventCuts;     }
    void                   SetEventCuts (AliMultSelectionCuts * var);
    AliMultSelection     * GetMultSelection() const                       { return fSelection;     }
    void                   SetMultSelection(AliMultSelection * var );
    
    // Make it browsable
    void Browse(TBrowser *b);
    virtual Bool_t IsFolder() const { return kTRUE; }
    void Dissociate();
    
    //Use internal map
    void Setup();
    TH1F* FindHisto(AliMultEstimator* e);
    void Print(Option_t* option="") const;
    
private:
    TList *fCalibList; // Calibration Histograms
    AliMultSelectionCuts * fEventCuts; // EventCuts
    AliMultSelection     * fSelection; // Definition of Estimators
    TMap*                  fMap; //! Map estimator to histogram
    ClassDef(AliOADBMultSelection, 1)
    
    
};

#endif
