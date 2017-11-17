/**********************************************
 *
 * Class designed to store all multiplicity
 * variables in a TClonesArray
 *
 * Instancing an empty AliMultInput class
 * will allow you to store standard variables
 *
 **********************************************/

#include "THashList.h"
#include "AliLog.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include <TROOT.h>

ClassImp(AliMultInput);

AliMultInput::AliMultInput(const char * name, const char * title)
  : TNamed(name, title)
  , fVariableList(nullptr)
{
    fVariableList = new THashList;
}

AliMultInput::AliMultInput(const AliMultInput& o)
  : TNamed(o)
  , fVariableList(nullptr)
{
    fVariableList = new THashList;
    TIter next(o.fVariableList);
    AliMultVariable* v = nullptr;
    while ((v = static_cast<AliMultVariable*>(next())))
      AddVariable(v);
}

AliMultInput& AliMultInput::operator=(const AliMultInput& o)
{
    if (&o == this)
      return *this;
    SetName(o.GetName());
    SetTitle(o.GetTitle());
    if (!fVariableList) fVariableList = new THashList;
    fVariableList->Clear();
    TIter next(o.fVariableList);
    AliMultVariable* v = nullptr;
    while ((v = static_cast<AliMultVariable*>(next())))
      AddVariable(v);
    return *this;
}

AliMultInput::~AliMultInput()
{
  SafeDelete(fVariableList);
}
void AliMultInput::AddVariable ( AliMultVariable *lVar )
{
    if (!lVar)
      return;

    if (!fVariableList)
      fVariableList = new THashList;

    //Protect against double-declaration
    if ( fVariableList->FindObject( lVar->GetName() ) ){
      Printf("===========================================================================" );
      Printf("                          !!!  WARNING !!!                                 " );
      Printf("Variable named %s already exists, you're doing a double-declaration!",lVar->GetName() );
      Printf("AddVariable call exiting without doing anything... Please check your logic!");
      Printf("===========================================================================" );
      return;
    }
    fVariableList->Add(lVar);
}

AliMultVariable* AliMultInput::GetVariable (const TString& lName) const
{
    if (!fVariableList) return 0;
    return static_cast<AliMultVariable*>(fVariableList->FindObject(lName));
}

AliMultVariable* AliMultInput::GetVariable (Long_t iIdx) const
{
    if (!fVariableList) return nullptr;
    if (iIdx < 0 || iIdx >= GetNVariables()) return nullptr;
    return static_cast<AliMultVariable*>(fVariableList->At(iIdx));
}

void AliMultInput::Clear(Option_t* option)
{
    TIter next(fVariableList);
    AliMultVariable* var = nullptr;
    while ((var = static_cast<AliMultVariable*>(next()))) {
        var->Clear(option);
    }
}

void AliMultInput::Set(const AliMultInput* other)
{
    TIter next(fVariableList);
    AliMultVariable* var = nullptr;
    while ((var = static_cast<AliMultVariable*>(next()))) {
        AliMultVariable* ovar = other->GetVariable(var->GetName());
        var->Set(ovar);
    }
}
void AliMultInput::Print(Option_t* option) const
{
    Printf("%s: %s/%s %ld variables", ClassName(),
           GetName(), GetTitle(), GetNVariables());
    gROOT->IndentLevel();
    Printf("Variables");
    gROOT->IncreaseDirLevel();
    TIter next(fVariableList);
    TObject* o = nullptr;
    while ((o = next())) {
        gROOT->IndentLevel();
        o->Print(option);
    }
    gROOT->DecreaseDirLevel();
}
Bool_t AliMultInput::SetValue(TString name, Float_t val)
{
  AliMultVariable *v = GetVariable(name);
  if (!v) {
    AliWarningF("variable '%s' not found", name.Data());
    return kFALSE;
  }
  v->SetValue(val);
  return kTRUE;
}
Bool_t AliMultInput::SetValue(TString name, Int_t val)
{
  AliMultVariable *v = GetVariable(name);
  if (!v) {
    AliWarningF("variable '%s' not found", name.Data());
    return kFALSE;
  }
  v->SetValueInteger(val);
  return kTRUE;
}
Bool_t AliMultInput::IncrementValue(TString name, Int_t increment)
{
  AliMultVariable *v = GetVariable(name);
  if (!v) {
    AliWarningF("variable '%s' not found", name.Data());
    return kFALSE;
  }
  v->SetValueInteger(v->GetValueInteger()+increment);
  return kTRUE;
}
Int_t AliMultInput::GetValueInteger(TString name) const
{
  AliMultVariable *v = GetVariable(name);
  if (!v) {
    AliWarningF("variable '%s' not found", name.Data());
    return -1;
  }
  return v->GetValueInteger();
}
Float_t AliMultInput::GetValue(TString name) const
{
  AliMultVariable *v = GetVariable(name);
  if (!v) {
    AliWarningF("variable '%s' not found", name.Data());
    return -1.0f;
  }
  return v->GetValue();
}
