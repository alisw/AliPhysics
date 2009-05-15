//
// Class AliRsnVManager
//
// Base "manager" class.
// It is built in order to manage a list of objects which share
// the same level in the work flow of the analysis.
// This base class contains simply the list of "child" objects
// and the methods to add objects to the list or retrieve the list.
//
// author     : M. Vala       [martin.vala@cern.ch]
// revised by : A. Pulvirenti [alberto.pulvirenti@ct.infn.it]
//

#ifndef AliRsnVManager_H
#define AliRsnVManager_H

#include <TNamed.h>

class AliRsnVManager : public TNamed
{
  public:

    AliRsnVManager(const char*name = "default");
    ~AliRsnVManager();

    virtual void        Add(TObject *pair);
            TObjArray*  GetArray() {return &fArray;}
            Int_t       GetEntries() {return fArray.GetEntries();}
            Int_t       GetEntriesFast() {return fArray.GetEntriesFast();}
    virtual void        PrintArray() const;
    virtual void        Print(Option_t *opt = "") const;

  protected:

    TObjArray  fArray;  // the managed array

    ClassDef(AliRsnVManager, 1)
};

#endif
