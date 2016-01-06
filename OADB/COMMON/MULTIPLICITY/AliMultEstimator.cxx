/**********************************************
 *
 * This is the basic estimator definition class
 *
 *  Estimator definitions have to contain 
 *  the following: 
 *
 *  --- Estimator Definition (logic)
 *  --- Estimator Value
 *  --- Estimator Mean Value
 *  --- Estimator Percentile
 *
 **********************************************/

#include "AliMultInput.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "TFolder.h"
#include "TObjString.h"
#include "TBrowser.h"
#include "TFormula.h"
#include "RVersion.h"

ClassImp(AliMultEstimator);
//________________________________________________________________
AliMultEstimator::AliMultEstimator() :
  TNamed(), fDefinition(""), fIsInteger(kFALSE), fValue(0), fMean(0), fPercentile(0), fFormula(0),
fkUseAnchor(kFALSE), fAnchorPoint(0), fAnchorPercentile(100.0)
{
  // Constructor
  
}
AliMultEstimator::AliMultEstimator(const char * name, const char * title, TString lInitDef):
TNamed(name,title), fDefinition(""), fIsInteger(kFALSE), fValue(0), fMean(0), fPercentile(0), fFormula(0),
fkUseAnchor(kFALSE), fAnchorPoint(0), fAnchorPercentile(100.0)
{
    //Named, titled, definition constructor
    fDefinition=lInitDef;
}
//________________________________________________________________
AliMultEstimator::AliMultEstimator(const AliMultEstimator& e)
: TNamed(e),
fDefinition(e.fDefinition),
fIsInteger(e.fIsInteger),
fValue(e.fValue),
fMean(e.fMean),
fPercentile(e.fPercentile),
fFormula(0),
fkUseAnchor(e.fkUseAnchor),
fAnchorPoint(e.fAnchorPoint),
fAnchorPercentile(e.fAnchorPercentile)
{
  if (e.fFormula) fFormula = new TFormula(*e.fFormula);
}
//________________________________________________________________
AliMultEstimator& AliMultEstimator::operator=(const AliMultEstimator& e)
{
    if (this == &e) return *this;
    SetName(e.GetName());
    SetTitle(e.GetTitle());
    fDefinition  = e.fDefinition;
    fIsInteger   = e.fIsInteger;
    fValue       = e.fValue;
    fMean        = e.fMean;
    fPercentile  = e.fPercentile;
    
    if (fFormula) delete fFormula;
    fFormula = 0;
    if (e.fFormula) fFormula = new TFormula(*e.fFormula);
    
    //Anchor point configs
    fkUseAnchor         = e.fkUseAnchor;
    fAnchorPoint        = e.fAnchorPoint;
    fAnchorPercentile   = e.fAnchorPercentile;
    
    return *this;
}
//________________________________________________________________
void AliMultEstimator::Set(const AliMultEstimator* e)
{
    if (!e) {
        fValue = 0;
        fMean  = 0;
        fPercentile = 0;
    }
    fValue       = e->fValue;
    fMean        = e->fMean;
    fPercentile  = e->fPercentile;
    
    //Anchor point configs
    fkUseAnchor         = e->fkUseAnchor;
    fAnchorPoint        = e->fAnchorPoint;
    fAnchorPercentile   = e->fAnchorPercentile;
}
//________________________________________________________________
AliMultEstimator::~AliMultEstimator(){
  // destructor
  if (fFormula) delete fFormula;   
}
//________________________________________________________________
Float_t AliMultEstimator::GetZ() const {
    //Don't die for zero
    Double_t lReturnVal = 0;
    if ( TMath::Abs( fValue ) > 1e-6 ){
        lReturnVal = fValue/fMean;
    }
    return lReturnVal; 
}
//________________________________________________________________
void AliMultEstimator::SetupFormula(const AliMultInput* lInput)
{
    TString expr = fDefinition;
    Int_t   nVar = lInput->GetNVariables();
    for (Int_t i = 0; i < nVar; i++) {
        TString repl(Form("[%d]", i));
        TString lVarName = lInput->GetVariable(i)->GetName();
        //IMPORTANT: this is necessary as names may have a common component!
        //Example: fAmplitude_V0A and fAmplitude_V0AEq
        //Required in syntax: parenthesis around all variables
        lVarName.Append (")");
        lVarName.Prepend("(");
        expr.ReplaceAll(lVarName, repl);
    }
    fFormula = new TFormula(Form("e%s", GetName()), expr);
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,4)
    fFormula->Optimize();
#endif
}
//________________________________________________________________
Float_t AliMultEstimator::Evaluate(const AliMultInput* lInput)
{
    if (!fFormula) return fValue = 0;
    for (Int_t i = 0; i < lInput->GetNVariables(); i++) {
        AliMultVariable* v = lInput->GetVariable(i);
        fFormula->SetParameter(i, v->IsInteger() ?
                               v->GetValueInteger() :
                               v->GetValue());
    }
    return fValue = fFormula->Eval(0);
}
