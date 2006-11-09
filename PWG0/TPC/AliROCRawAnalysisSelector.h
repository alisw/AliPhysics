/* $Id$ */

#ifndef AliROCRawAnalysisSelector_H
#define AliROCRawAnalysisSelector_H

#include "TSelector.h"

// 
// TODO explain this
//

class AliRawEvent;
class TTree;

class AliROCRawAnalysisSelector : public TSelector {
  public:
  
    AliROCRawAnalysisSelector();
    virtual ~AliROCRawAnalysisSelector();

    virtual Int_t   Version() const {return 1;}
    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    
    AliRawEvent* fRawEvent;
    TTree*       fTree;

 private:
    AliROCRawAnalysisSelector(const AliROCRawAnalysisSelector&);
    AliROCRawAnalysisSelector& operator=(const AliROCRawAnalysisSelector&);

  ClassDef(AliROCRawAnalysisSelector, 0);
};

#endif
