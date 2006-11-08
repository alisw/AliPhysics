/* $Id$ */

#ifndef AliROCRawAnalysisSelector_H
#define AliROCRawAnalysisSelector_H

#include "AliSelector.h"


// 
// TODO explain this
//

class AliROCRawAnalysisSelector : public AliSelector {
  public:
  
    AliROCRawAnalysisSelector();
    virtual ~AliROCRawAnalysisSelector();

    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:

    

 private:
    AliROCRawAnalysisSelector(const AliROCRawAnalysisSelector&);
    AliROCRawAnalysisSelector& operator=(const AliROCRawAnalysisSelector&);

  ClassDef(AliROCRawAnalysisSelector, 0);
};

#endif
