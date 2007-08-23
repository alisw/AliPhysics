/* $Id$ */

#ifndef AliHighMultiplicitySelector_H
#define AliHighMultiplicitySelector_H

#include "AliSelectorRL.h"

class TH1F;
class TH2F;
class TH1;
class TH2;
class TNtuple;

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
    void Ntrigger();

 protected:
    void MakeGraphs(const char* title, TH1* xSection, TH2* fMvsL, Int_t limit);
    void MakeGraphs2(const char* title, TH1* xSection, TH2* fMvsL);

    TH1* GetXSectionCut(TH1* xSection, TH2* multVsLayer, Int_t cut);
    TH1* GetTriggerEfficiency(TH2* multVsLayer, Int_t cut);

    TH1F* fChipsLayer1;
    TH1F* fChipsLayer2;

    TH2F* fL1vsL2;
    TH2F* fMvsL1;
    TH2F* fMvsL2;

    TH2F* fChipsFired;

    TNtuple* fPrimaryL1;  // multiplicity vs. number of fired chips vs. number of chips fired by primaries
    TNtuple* fPrimaryL2;  // multiplicity vs. number of fired chips vs. number of chips fired by primaries

    TH1F*    fClusterZL1;
    TH1F*    fClusterZL2;

    TH2F* fClvsL1;
    TH2F* fClvsL2;

    Bool_t centralRegion;

 private:
    AliHighMultiplicitySelector(const AliHighMultiplicitySelector&);
    AliHighMultiplicitySelector& operator=(const AliHighMultiplicitySelector&);

  ClassDef(AliHighMultiplicitySelector, 0);
};

#endif
