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

#ifndef ALIRsnFunctionNew_H
#define ALIRsnFunctionNew_H

#include <TArrayD.h>
#include <TString.h>

//#include "AliRsnCut.h"
#include "AliRsnFunctionDef.h"
#include "AliRsnPairParticle.h"

class TH1D;
class TH2D;
class AliRsnEvent;

class AliRsnFunctionNew : public TObject
{

  public:

    AliRsnFunctionNew();
    AliRsnFunctionNew(AliRsnFunctionDef *fd);
    AliRsnFunctionNew(const AliRsnFunctionNew &copy);
    virtual ~AliRsnFunctionNew() {Clear();}
    virtual void Clear(Option_t *option = "");

    void               SetFunctionDef(AliRsnFunctionDef *def) {fFcnDef = def;}
    AliRsnFunctionDef* GetFunctionDef() {return fFcnDef;}
    TString            GetFcnName();
    TString            GetFcnTitle();

    // working routines
    TList*   Init(const char *histoName, const char *histoTitle);
    void     Init(const char *histoName, const char *histoTitle, TList *tgt);
    Bool_t   Fill(AliRsnPairParticle *pair, AliRsnPairDef *ref);

  private:

    const AliRsnFunctionNew& operator=(const AliRsnFunctionNew &copy);

    AliRsnFunctionDef  *fFcnDef;                               // definitions
    TH1D               *fHistoTot;                             // integrated histogram
    TH2D               *fHistoBin[AliRsnFunctionDef::kBinMax]; // binned histograms

    // ROOT dictionary
    ClassDef(AliRsnFunctionNew, 1)
};

#endif
