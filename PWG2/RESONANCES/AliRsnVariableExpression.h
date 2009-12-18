//
// AliRsnVariableExpresion class
// is used
// to help AliRsnExpresion class
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNVARIABLEEXPRESSION_H
#define ALIRSNVARIABLEEXPRESSION_H

#include "AliRsnExpression.h"

class AliRsnVariableExpression: public AliRsnExpression
{
  public:
    AliRsnVariableExpression(TString a) : AliRsnExpression() { fVname = a;  };
    ~AliRsnVariableExpression() {}
    virtual Bool_t    Value(TObjArray& pgm);
    virtual TString    Unparse() const { return fVname; }

    ClassDef(AliRsnVariableExpression, 1);    // Class to define a variable expression
};

#endif
