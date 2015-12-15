/**********************************************
 *
 * Class designed to store all multiplicity
 * variables in a TClonesArray
 *
 * Instancing an empty AliMultInput class
 * will allow you to store standard variables
 *
 **********************************************/

#include "TList.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include <TROOT.h>

ClassImp(AliMultInput);

AliMultInput::AliMultInput() :
  TNamed(), fNVars(0), fVariableList(0x0)
{
  // Constructor
    fVariableList = new TList();
}

AliMultInput::AliMultInput(const char * name, const char * title):
TNamed(name,title), fNVars(0), fVariableList(0x0)
{
  // Constructor
    fVariableList = new TList();
}

AliMultInput::AliMultInput(const AliMultInput& o)
: TNamed(o), fNVars(0), fVariableList(0x0)
{
    // Constructor
    fVariableList = new TList();
    TIter next(o.fVariableList);
    AliMultVariable* v  = 0;
    while ((v = static_cast<AliMultVariable*>(next())))  AddVariable(v);
}

AliMultInput& AliMultInput::operator=(const AliMultInput& o)
{
    if (&o == this) return *this;
    SetName(o.GetName());
    SetTitle(o.GetTitle());
    if (!fVariableList) fVariableList = new TList();
    fVariableList->Clear();
    fNVars = 0;
    TIter next(o.fVariableList);
    AliMultVariable* v  = 0;
    while ((v = static_cast<AliMultVariable*>(next())))  AddVariable(v);
    return *this;
}

AliMultInput::~AliMultInput(){
  // destructor
  
}
void AliMultInput::AddVariable ( AliMultVariable *lVar )
{
    if (!lVar) return;
    if (!fVariableList) {
        fVariableList = new TList;
        fNVars = 0;
    }
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
    fNVars++;
}

AliMultVariable* AliMultInput::GetVariable (const TString& lName) const
{
    if (!fVariableList) return 0;
    return static_cast<AliMultVariable*>(fVariableList->FindObject(lName));
}

AliMultVariable* AliMultInput::GetVariable (Long_t iIdx) const
{
    if (!fVariableList) return 0;
    if (iIdx < 0 || iIdx >= fNVars) return 0;
    return static_cast<AliMultVariable*>(fVariableList->At(iIdx));
}

void AliMultInput::Clear(Option_t* option)
{
    TIter next(fVariableList);
    AliMultVariable* var = 0;
    while ((var = static_cast<AliMultVariable*>(next()))) {
        var->Clear(option);
    }
}

void AliMultInput::Set(const AliMultInput* other)
{
    TIter next(fVariableList);
    AliMultVariable* var = 0;
    while ((var = static_cast<AliMultVariable*>(next()))) {
        AliMultVariable* ovar = other->GetVariable(var->GetName());
        var->Set(ovar);
    }
}
void AliMultInput::Print(Option_t* option) const
{
    Printf("%s: %s/%s %ld variables", ClassName(),
           GetName(), GetTitle(), fNVars);
    gROOT->IndentLevel();
    Printf("Variables");
    gROOT->IncreaseDirLevel();
    TIter next(fVariableList);
    TObject* o = 0;
    while ((o = next())) {
        gROOT->IndentLevel();
        o->Print(option);
    }
    gROOT->DecreaseDirLevel();
}
