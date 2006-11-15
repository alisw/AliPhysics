/* $Id$ */

#ifndef AliROCESDAnalysisSelector_H
#define AliROCESDAnalysisSelector_H

#include "AliSelector.h"

class AliTPCClusterHistograms;
class AliESD;
class AliESDfriend;
class AliTPCseed;

class TObjArray;

// 
// TODO explain this
//

class AliROCESDAnalysisSelector : public AliSelector {
  public:
    enum { kTPCSectors = 72, kTPCHists = kTPCSectors * 2 };
  
    AliROCESDAnalysisSelector();
    virtual ~AliROCESDAnalysisSelector();

    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    Int_t           ProcessEvent(Bool_t detailedHistogram=kFALSE, const Char_t* label="");

    Bool_t          AcceptTrack(const AliTPCseed* track, Int_t minRowsIncluded=0);

 protected:
    AliESDfriend* fESDfriend;  // ESD friend pointer

    AliTPCClusterHistograms* fClusterHistograms[kTPCHists]; // 0..71 histograms created with all clusters, 72..143 without edges

 private:

    TObjArray* fObjectsToSave;

    Int_t fMinNumberOfRowsIsTrack;

    AliROCESDAnalysisSelector(const AliROCESDAnalysisSelector&);
    AliROCESDAnalysisSelector& operator=(const AliROCESDAnalysisSelector&);

  ClassDef(AliROCESDAnalysisSelector, 0);
};

#endif
