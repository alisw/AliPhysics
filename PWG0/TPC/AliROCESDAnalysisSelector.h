/* $Id$ */

#ifndef AliROCESDAnalysisSelector_H
#define AliROCESDAnalysisSelector_H

#include "AliSelector.h"

class AliESDfriend;
class TH2F;

// this is an empty selector that can be used to create an analysis

class AliROCESDAnalysisSelector : public AliSelector {
  public:
    AliROCESDAnalysisSelector();
    virtual ~AliROCESDAnalysisSelector();

    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    AliESDfriend* fESDfriend;
    
    TH2F* fhQmaxVsRow;

 private:
    AliROCESDAnalysisSelector(const AliROCESDAnalysisSelector&);
    AliROCESDAnalysisSelector& operator=(const AliROCESDAnalysisSelector&);

  ClassDef(AliROCESDAnalysisSelector, 0);
};

#endif
