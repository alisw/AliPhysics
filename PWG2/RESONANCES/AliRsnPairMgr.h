//
// *** Class AliRsnPairMgr ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNMVPAIRMGR_H
#define ALIRSNMVPAIRMGR_H

#include "AliRsnPair.h"

class AliRsnPairMgr : public TNamed
{
  public:
    AliRsnPairMgr(const char*name="default");

    ~AliRsnPairMgr();

    void          AddPair ( AliRsnPair *pair );
    TObjArray     *GetPairs() { return &fPairs; }
    void          PrintPairs();

  private:

    TObjArray       fPairs;                 // array of pairs
    
    ClassDef ( AliRsnPairMgr, 1 )
};

#endif
