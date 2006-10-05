/* $Id$ */

#ifndef ALIDNDETAANALYSISSELECTORMC_H
#define ALIDNDETAANALYSISSELECTORMC_H

#include "AliSelectorRL.h"

class TH3F;
class TH1F;
class dNdEtaAnalysis;

class AlidNdEtaAnalysisMCSelector : public AliSelectorRL {
  public:
    AlidNdEtaAnalysisMCSelector();
    virtual ~AlidNdEtaAnalysisMCSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual void    SlaveTerminate();
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    Terminate();

 protected:

 private:
    dNdEtaAnalysis* fdNdEtaAnalysis;      // contains the intermediate histograms (on each slave)

    TH3F* fVertex;  //! vertex of counted particles
    TH1F* fPartEta; //! counted particles as function of eta
    TH1F* fPartPt; //! counted particles as function of pt
    Int_t fEvents;  //! number of processed events

    AlidNdEtaAnalysisMCSelector(const AlidNdEtaAnalysisMCSelector&);
    AlidNdEtaAnalysisMCSelector& operator=(const AlidNdEtaAnalysisMCSelector&);

    ClassDef(AlidNdEtaAnalysisMCSelector, 0);
};

#endif
