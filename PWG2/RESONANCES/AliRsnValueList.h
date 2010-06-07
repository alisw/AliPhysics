//
// Class AliRsn Fcn
//
// This class defines a base classe to implement a typical computation
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from
// this class and override the virtual functions defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNVALUELIST_H
#define ALIRSNVALUELIST_H

#include <TClonesArray.h>
#include <TArrayD.h>
#include <TNamed.h>

class AliRsnValue;
class AliRsnPairParticle;
class AliRsnPairDef;

class AliRsnValueList : public TNamed
{

  public:

    AliRsnValueList();
    AliRsnValueList(const AliRsnValueList &copy);
    virtual ~AliRsnValueList() { Clear(""); }
    const AliRsnValueList& operator=(const AliRsnValueList &copy);

    virtual void Clear(Option_t* /*option = ""*/) {fValueList.Delete();}
    Int_t        GetNumberOfValues() {return fValueList.GetEntries();}

    void         AddValue(const AliRsnValue * const axis);
    Bool_t       Eval(TObject * const obj, const AliRsnPairDef *pairDef = 0x0);
    Double_t     GetValue(Int_t i) {if (i>=0 && i<fArray.GetSize()) return fArray[i]; return 0.0;}

  protected:

    TClonesArray  fValueList;  // list of computators
    TArrayD       fArray;      // list of values

    // ROOT dictionary
    ClassDef(AliRsnValueList, 3)
};

#endif
