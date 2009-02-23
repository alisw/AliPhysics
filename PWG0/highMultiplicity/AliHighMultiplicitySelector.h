/* $Id$ */

#ifndef AliHighMultiplicitySelector_H
#define AliHighMultiplicitySelector_H

#include "AliSelectorRL.h"

class TH1F;
class TH2F;
class TH1;
class TH2;
class TNtuple;
class TGraph;

//
// Selector that produces plots needed for the high multiplicity analysis with the
// pixel detector
//

class AliHighMultiplicitySelector : public AliSelectorRL {
  public:
    AliHighMultiplicitySelector();
    virtual ~AliHighMultiplicitySelector();

    virtual void    Init(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual Bool_t  Notify();
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    void WriteHistograms(const char* filename = "highmult.root");
    void ReadHistograms(const char* filename = "highmult.root");
    void DrawHistograms();
    void JPRPlots();
    void Ntrigger(Bool_t relative = kTRUE);
    TGraph* IntFractRate();
    void Contamination();
    void Contamination2();
    void Contamination3();
    void MBComparison();

 protected:
    void MakeGraphs(const char* title, TH1* xSection, TH2* fMvsL, Int_t limit);
    void MakeGraphs2(const char* title, TH1* xSection, TH2* fMvsL);

    TH1* GetXSectionCut(TH1* xSection, TH2* multVsLayer, Int_t cut);
    TH1* GetTriggerEfficiency(TH2* multVsLayer, Int_t cut);

    TH1F* fChipsLayer1;   // fired chips in layer 1
    TH1F* fChipsLayer2;   // fired chips in layer 2

    TH2F* fL1vsL2;        // layer1 hits vs. layer2 hits
    TH2F* fMvsL1;         // mult. vs. hits in layer 1
    TH2F* fMvsL2;         // mult. vs. hits in layer 2

    TH2F* fChipsFired;    // chips fired, module number vs. chip number

    TNtuple* fPrimaryL1;  // multiplicity vs. number of fired chips vs. number of chips fired by primaries
    TNtuple* fPrimaryL2;  // multiplicity vs. number of fired chips vs. number of chips fired by primaries

    TH1F*    fClusterZL1; // number of clusters as funtion of z in layer 1
    TH1F*    fClusterZL2; // number of clusters as funtion of z in layer 2

    TH2F* fClvsL1;        // number of cluster vs. number of fired chips in layer 1
    TH2F* fClvsL2;        // number of cluster vs. number of fired chips in layer 2

    Bool_t fCentralRegion; // only consider the central region (two central modules)

 private:
    AliHighMultiplicitySelector(const AliHighMultiplicitySelector&);
    AliHighMultiplicitySelector& operator=(const AliHighMultiplicitySelector&);

  ClassDef(AliHighMultiplicitySelector, 0);
};

#endif
