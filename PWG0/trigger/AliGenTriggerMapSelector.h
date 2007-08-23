/* $Id$ */

#ifndef AliGenTriggerMapSelector_H
#define AliGenTriggerMapSelector_H

#include "AliSelectorRL.h"

class TH1F;
class TH2F;
class TH1;
class TH2;
class TNtuple;

class AliGenTriggerMapSelector : public AliSelectorRL {
  public:
    AliGenTriggerMapSelector();
    virtual ~AliGenTriggerMapSelector();

    virtual void    Init(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual Bool_t  Notify();
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    void WriteHistograms(const char* filename = "triggerMap.root");
    void ReadHistograms(const char* filename = "triggerMap.root");

 protected:
    TH2F* fChipsFired;

    TNtuple* fTracklets;  // vertex vs. chip_l1 vs. chip_l2 for all tracklets

 private:
    AliGenTriggerMapSelector(const AliGenTriggerMapSelector&);
    AliGenTriggerMapSelector& operator=(const AliGenTriggerMapSelector&);

  ClassDef(AliGenTriggerMapSelector, 0);
};

#endif
