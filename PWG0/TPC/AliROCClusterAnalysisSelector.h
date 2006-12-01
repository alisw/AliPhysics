/* $Id$ */

#ifndef AliROCClusterAnalysisSelector_H
#define AliROCClusterAnalysisSelector_H

#include "AliSelectorRL.h"

class AliTPCClusterHistograms;

class TObjArray;

// 
// TODO explain this
//

class AliROCClusterAnalysisSelector : public AliSelectorRL {
  public:
    enum { kTPCSectors = 72, kTPCHists = kTPCSectors * 2 };
  
    AliROCClusterAnalysisSelector();
    virtual ~AliROCClusterAnalysisSelector();

    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    Int_t           ProcessEvent(Long64_t entry, Bool_t detailedHistogram=kFALSE, const Char_t* label="");


 protected:

    AliTPCClusterHistograms* fClusterHistograms[kTPCHists]; // 0..71 histograms created with all clusters, 72..143 without edges

 private:

    Int_t      fNMaxObjectsToSave;
    TObjArray* fObjectsToSave;


    AliROCClusterAnalysisSelector(const AliROCClusterAnalysisSelector&);
    AliROCClusterAnalysisSelector& operator=(const AliROCClusterAnalysisSelector&);

  ClassDef(AliROCClusterAnalysisSelector, 0);
};

#endif
