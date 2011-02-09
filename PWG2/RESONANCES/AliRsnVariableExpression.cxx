//
// AliRsnVariableExpresion class
// is used
// to help AliRsnExpresion class
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"

#include "AliRsnCut.h"
#include "AliRsnVariableExpression.h"

ClassImp(AliRsnVariableExpression)

//______________________________________________________________________________
Bool_t AliRsnVariableExpression::Value(TObjArray& /*pgm*/)
{

//   Int_t indexx = fgCutSet->GetIndexByCutName ( fVname.Data() );
   AliDebug(AliLog::kDebug, Form("Vname %s", fVname.Data()));
//   return fgCutSet->GetBoolValue ( indexx );

   return fgCutSet->GetBoolValue(fVname.Atoi());
}

