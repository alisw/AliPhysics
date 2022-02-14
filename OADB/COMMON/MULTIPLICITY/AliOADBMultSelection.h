#ifndef ALIOADBMULTSELECTION_H
#define ALIOADBMULTSELECTION_H

#include <TNamed.h>
#include <TMap.h>
#include <AliMultSelection.h>
class TBrowser;
class TH1F;
class TProfile;
class TList; 
class AliMultSelectionCuts;
class AliMultEstimator;

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
  TProfile *GetCalibHistoVtx(Long_t iVar) const;
  TProfile *GetCalibHistoVtx(const TString& lVarHistoName) const;
  
  void AddCalibHisto (TH1F * var);
  void AddCalibHistoVtx (TProfile * varVtx);
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
  TList *fCalibListVtxZ; // Calibration Histograms for vertex-Z
  AliMultSelectionCuts * fEventCuts; // EventCuts
  AliMultSelection     * fSelection; // Definition of Estimators
  TMap*                  fMap; //! Map estimator to histogram
  TMap*                  fMapVtxZ; //! Map raw var to histogram
  ClassDef(AliOADBMultSelection, 2)
  //2 - Add Vertex-Z
  
  
};

#endif
