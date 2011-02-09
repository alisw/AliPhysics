//
// AliRsnExpresion class is used to
// handle operators &|!
// in AliRsnCut
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNEXPRESSION_H
#define ALIRSNEXPRESSION_H

#include <TObject.h>

class TObjArray;
#include "AliRsnCutSet.h"
class AliRsnVariableExpression;

class AliRsnExpression : public TObject {

public:

   // operators for complex cut expressions
   enum ECutOp {
      kOpAND = 1, // AND '&'
      kOpOR,      // OR '|'
      kOpNOT      // Unary negation '!'
   };

   AliRsnExpression() : fVname(0), fArg1(0), fArg2(0), fOperator(0)  {}
   AliRsnExpression(TString exp);
   virtual    ~AliRsnExpression();
   AliRsnExpression(const AliRsnExpression& exp);
   AliRsnExpression&    operator= (const AliRsnExpression& exp);

   virtual Bool_t     Value(TObjArray & vars);
   virtual TString     Unparse() const;

   void SetCutSet(AliRsnCutSet* const theValue) { fgCutSet = theValue; }
   AliRsnCutSet* GetCutSet() const { return fgCutSet; }


   TString                     fVname;   // Variable name
   static AliRsnCutSet        *fgCutSet; // global cutset

private:
   AliRsnExpression*           fArg1;         // left argument
   AliRsnExpression*           fArg2;         // right argument
   Int_t                       fOperator;     // operator

   AliRsnExpression(int op, AliRsnExpression* a);
   AliRsnExpression(int op, AliRsnExpression* a, AliRsnExpression* b);

   TObjArray*    Tokenize(TString str) const;
   static AliRsnExpression*    Element(TObjArray &st, Int_t &i);
   static AliRsnExpression*    Primary(TObjArray &st, Int_t &i);
   static AliRsnExpression*    Expression(TObjArray &st, Int_t &i);

   ClassDef(AliRsnExpression, 1);    // Class to evaluate an expression
};

#endif
