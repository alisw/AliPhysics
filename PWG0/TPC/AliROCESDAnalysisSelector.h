/* $Id$ */

#ifndef AliROCESDAnalysisSelector_H
#define AliROCESDAnalysisSelector_H

#include "AliSelector.h"

//#include "TPC/AliTPCClusterHistograms.h"

class AliTPCClusterHistograms;
class AliESDfriend;

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

    AliTPCClusterHistograms* fClusterHistograms;

 private:
    AliROCESDAnalysisSelector(const AliROCESDAnalysisSelector&);
    AliROCESDAnalysisSelector& operator=(const AliROCESDAnalysisSelector&);

  ClassDef(AliROCESDAnalysisSelector, 0);
};

#endif
