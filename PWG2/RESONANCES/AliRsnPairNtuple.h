//
// *** Class AliRsnPairNtuple ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef AliRsnPairNtuple_H
#define AliRsnPairNtuple_H

#include "AliRsnPair.h"

class TList;
class TNtuple;

class AliRsnPairNtuple : public AliRsnPair
{
  public:

    AliRsnPairNtuple(const char *name = "default", AliRsnPairDef *def = 0);
    AliRsnPairNtuple(const AliRsnPairNtuple &copy);
    AliRsnPairNtuple& operator=(const AliRsnPairNtuple&);
    ~AliRsnPairNtuple();

    Bool_t       AddValue(AliRsnValue*const val);
    void         GenerateNtuple(const char *prefix = "", TList *list = 0);
    virtual void Compute();
    virtual void Init(const char *prefix, TList *list);

  private:

    TClonesArray  fValues;  // single values computed from analyzed objects
    TNtuple      *fNtuple;  // ntuple computed with values

    ClassDef(AliRsnPairNtuple, 2)
};

#endif

