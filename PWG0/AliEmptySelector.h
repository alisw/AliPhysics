/* $Id$ */

#ifndef AliEmptySelector_H
#define AliEmptySelector_H

#include "AliSelector.h"

// this is an empty selector that can be used to create an analysis

class AliEmptySelector : public AliSelector {
  public:
    AliEmptySelector();
    virtual ~AliEmptySelector();

    virtual Bool_t  Process(Long64_t entry);

 private:
    AliEmptySelector(const AliEmptySelector&);
    AliEmptySelector& operator=(const AliEmptySelector&);

  ClassDef(AliEmptySelector, 0);
};

#endif
