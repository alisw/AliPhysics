/* $Id$ */

#ifndef AliROCESDAnalysisSelector_H
#define AliROCESDAnalysisSelector_H

#include "AliSelector.h"

class AliTPCClusterHistograms;
class AliESDfriend;

// 
// TODO explain this
//

class AliROCESDAnalysisSelector : public AliSelector {
  public:
    enum { kTPCSectors = 72 };
  
    AliROCESDAnalysisSelector();
    virtual ~AliROCESDAnalysisSelector();

    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    AliESDfriend* fESDfriend;

    AliTPCClusterHistograms* fClusterHistograms[kTPCSectors];

 private:
    AliROCESDAnalysisSelector(const AliROCESDAnalysisSelector&);
    AliROCESDAnalysisSelector& operator=(const AliROCESDAnalysisSelector&);

  ClassDef(AliROCESDAnalysisSelector, 0);
};

#endif
