//
// AliRsnExpresion class is used to
// handle operators &|!
// in AliRsnCut
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliLog.h"

#include "AliRsnVariableExpression.h"
#include "AliRsnExpression.h"

ClassImp(AliRsnExpression)

AliRsnCutSet *AliRsnExpression::fgCutSet = 0;

//______________________________________________________________________________
AliRsnExpression::AliRsnExpression(TString exp) :
   TObject(),
   fVname(""),
   fArg1(0x0),
   fArg2(0x0),
   fOperator(0)
{
   // Default constructor
   TObjArray* tokens = Tokenize(exp);

   Int_t i = -1;
   AliRsnExpression* e = Expression(*tokens, i);
   // Copy !!!
   fArg1 = e->fArg1; e->fArg1 = 0;
   fArg2 = e->fArg2; e->fArg2 = 0;
   fOperator = e->fOperator;
   fVname = e->fVname;
   delete e;
   delete tokens;
}

//______________________________________________________________________________
AliRsnExpression::~AliRsnExpression()
{
   if (fArg1) delete fArg1;
   if (fArg2) delete fArg2;
}

AliRsnExpression::AliRsnExpression(const AliRsnExpression & exp) : TObject(exp),
   fVname(exp.fVname),
   fArg1(exp.fArg1),
   fArg2(exp.fArg2),
   fOperator(exp.fOperator)
{}

//______________________________________________________________________________
AliRsnExpression& AliRsnExpression::operator= (const AliRsnExpression& e)
{
   // AliRsnExpression assignment operator.

   if (this != &e) {
      TObject::operator= (e);
      fArg1 = e.fArg1;
      fArg2 = e.fArg2;
      fOperator = e.fOperator;
      fVname = e.fVname;
   }
   return *this;
}

//______________________________________________________________________________
AliRsnExpression::AliRsnExpression(int op, AliRsnExpression* a, AliRsnExpression* b) :
   TObject(),
   fVname(""),
   fArg1(a),
   fArg2(b),
   fOperator(op)
{
   // Create a new expression
}

//______________________________________________________________________________
AliRsnExpression::AliRsnExpression(int op, AliRsnExpression* a) :
   TObject(),
   fVname(""),
   fArg1(0),
   fArg2(a),
   fOperator(op)
{
   // Create a unary expression.
}

//______________________________________________________________________________
Bool_t AliRsnExpression::Value(TObjArray &vars)
{
   //  Evaluate the expression
   if (fArg2 == 0 && fVname.IsNull()) {
      AliError("Expression undefined.");
      return kFALSE;
   }
   if (!fArg2) {
      AliError("Argument 2 is required.");
      return kFALSE;
   }


//   AliDebug(AliLog::kDebug,Form("fOperator %d",fOperator));

   switch (fOperator) {

      case kOpOR :
         return fArg1->Value(vars) || fArg2->Value(vars);

      case kOpAND :
         return fArg1->Value(vars) && fArg2->Value(vars);

      case kOpNOT :
         return !(fArg2->Value(vars));

      case 0 : {
//       Int_t indexx = fgCutSet->GetIndexByCutName ( fVname.Data() );
         AliDebug(AliLog::kDebug, Form("Vname %s", fVname.Data()));
//       return fgCutSet->GetBoolValue ( indexx );
         return fgCutSet->GetBoolValue(fVname.Atoi());
      }

      default:
         AliError("Illegal operator in expression!");

   }
   return kFALSE;
}


//______________________________________________________________________________
TString AliRsnExpression::Unparse() const
{
   // Unparse the expression

   TString opVals[4] = { "", "&", "|", "!" };
   if (fArg2 == 0 && fVname.IsNull()) {
      AliError("Expression undefined.");
      return "Error";
   }

   if (fArg2 == 0 && !fVname.IsNull()) return fVname;

   if (fArg1 == 0 && fArg2) {
      return opVals[fOperator] + fArg2->Unparse();
   }
   return "(" + fArg1->Unparse() + " " + opVals[fOperator] + " " + fArg2->Unparse() + ")";
}

//______________________________________________________________________________
TObjArray* AliRsnExpression::Tokenize(TString str) const
{
   // tokenize the expression

   // Remove spaces
   TString str1;
   for (Int_t i = 0; i < str.Length(); i++) {
      if (str[i] == ' ') continue;
      str1.Append(str[i]);
   }
   // get variable tokens
   TObjArray* valtok = str1.Tokenize("!&|()");
   // put all variables together
   Int_t nvt = valtok->GetEntriesFast();
   TString sumval;
   for (Int_t i = 0; i < nvt; i++) {
      TObjString* val = (TObjString*) valtok->At(i);
      sumval.Append(val->String());
   }
   // get the operator tokens
   TObjArray* optok = str1.Tokenize(sumval.Data());
   // put all operator in one string
   TString operators;
   Int_t nopt = optok->GetEntriesFast();
   for (Int_t i = 0; i < nopt; i++) {
      TObjString* val1 = (TObjString*) optok->At(i);
      operators.Append(val1->String());
   }
   // add more room to be safe
   TObjString* blank = new TObjString(" ");
   operators.Append(" ");
   valtok->AddLast(blank);
   // Now put var. and oper. together
   TObjArray* tokens = new TObjArray(valtok->GetEntriesFast() + operators.Length());
   int io = 0, iv = 0;
   int index = 0;
   while (1) {
      TString so = operators[io];
      int indexO = str1.Index(so, index);
      TString val2 = ((TObjString*) valtok->At(iv))->String();
      int indexV = str1.Index(val2, index);
      if ((indexO < indexV || indexV < 0) && indexO >= 0) {
         tokens->AddLast(new TObjString(so));
         index += so.Length();
         io++;
      }
      if ((indexV < indexO || indexO < 0) && indexV >= 0) {
         tokens->AddLast(new TObjString(val2));
         index += val2.Length();
         iv++;
      }
      if (index >= str1.Length()) break;
   }

// //  Debug -> Print the tokens
//   Int_t nt = tokens->GetEntriesFast();
//   for ( Int_t i=0; i<nt; i++ )
//   {
//     TObjString* val3 = ( TObjString* ) tokens->At ( i );
//     AliInfo ( Form ( "%d %s",i,val3->String().Data() ) );
//   }
//
//
   delete valtok;
   delete optok;

   return tokens;
}


//______________________________________________________________________________
AliRsnExpression* AliRsnExpression::Element(TObjArray &st, Int_t &i)
{
   // create an element

   AliRsnExpression* result = 0;

   Int_t nt = st.GetEntriesFast();
   TString token = "@";
   TObjString* valt;
   // next token
   if (i < nt - 1) {
      i++;
      valt = (TObjString*) st.At(i);
      token = valt->String();
   }
   // token type
   char ttok = (token[0] != '|' && token[0] != '&' &&
                token[0] != '!' && token[0] != '(' && token[0] != ')') ? 'w' : token[0];
   switch (ttok) {
      case 'w' : {
         result = new AliRsnVariableExpression(token);
         break;
      }
      case '(' :
         result = Expression(st, i);
         // next token
         if (i < nt - 1) {
            i++;
            valt = (TObjString*) st.At(i);
            token = valt->String();
         }
         if (token[0] != ')') {
            //       i--; // push back
            AliErrorGeneral("AliRsnExpression::Element", "Mismatched parenthesis.");
            delete result;
            result = new AliRsnExpression;
         }
         break;
      default:
         i--; // push back
         AliErrorGeneral("AliRsnExpression::Element", Form("Unexpected symbol on input. %s", token.Data()));
         //if (result) delete result;
         result = new AliRsnExpression;
   }
   return result;
}

//______________________________________________________________________________
AliRsnExpression* AliRsnExpression::Primary(TObjArray &st, Int_t &i)
{
   // create a primary

   Int_t nt = st.GetEntriesFast();
   TString token = "@";
   TObjString* valt;
   // next token
   if (i < nt - 1) {
      i++;
      valt = (TObjString*) st.At(i);
      token = valt->String();
   }

   switch (token[0]) {
      case '!' :
         return new AliRsnExpression(kOpNOT, Primary(st, i));
      default:
         i--; // push back
         return Element(st, i);
   }
}

//______________________________________________________________________________
AliRsnExpression* AliRsnExpression::Expression(TObjArray &st, Int_t &i)
{
   // create an expression

   AliRsnExpression* result = 0;
   Bool_t done = kFALSE;
   TString token;
   TObjString* valt;

   static int stack = 0;
   stack++;
   Int_t nt = st.GetEntriesFast();

   result = Primary(st, i);
//   cout <<"i "<<i<< "Primary " << result->Unparse() << endl;
   while (! done) {
      // next token
      if (i < nt - 1) i++;
      else break;
      valt = (TObjString*) st.At(i);
      token = valt->String();
      switch (token[0]) {
         case '&' :
            result = new AliRsnExpression(kOpAND, result, Primary(st, i));
//   cout <<"i "<<i<< " Expression AND " << result->Unparse() << endl;
            break;
         case '|' :
            result = new AliRsnExpression(kOpOR, result, Primary(st, i));
//   cout <<"i "<<i<< " Expression OR " << result->Unparse() << endl;
            break;
         default:
            done = kTRUE;
            i--; // push back
            break;
      }
   }
   stack--;
   if (stack == 0 && !token.IsNull() && token[0] == ')') {
      AliErrorGeneral("AliRsnExpression::Expression", "To many closing parenthesis.");
      delete result;
      result = new AliRsnExpression;
   } else if (stack == 0 && i < nt - 1) {
      AliErrorGeneral("AliRsnExpression::Expression", Form("Unexpected symbol on input. %s", token.Data()));
      delete result;
      result = new AliRsnExpression;
   }
   return result;
}
