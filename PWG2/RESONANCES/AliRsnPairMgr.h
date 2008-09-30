//
// Class AliRsnPairMgr
//
// A collection of pairs for an analysis.
// The function of this collection is just for purposes of well-sorting
// the analyzed pairs into upper-level groups, in the case of a wide
// analysis containing many resonances at once, or different settings for the same one.
//
// Each PairMgr will result in a separate list of histograms, which
// can be seen as a folder in the output file, whose name is given by this object.
//
// author: M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNPAIRMGR_H
#define ALIRSNPAIRMGR_H

#include "AliRsnPair.h"

class AliRsnPairMgr : public TNamed
{
  public:

    AliRsnPairMgr(const char*name = "default");
    ~AliRsnPairMgr();

    void       AddPair(AliRsnPair *pair);
    TObjArray *GetPairs() { return &fPairs; }
    void       PrintPairs();

  private:

    TObjArray  fPairs;  // array of pairs

    ClassDef(AliRsnPairMgr, 1)
};

#endif
