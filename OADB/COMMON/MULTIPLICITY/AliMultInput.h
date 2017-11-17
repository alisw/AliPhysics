#ifndef AliMultInput_H
#define AliMultInput_H
#include <TNamed.h>
#include <TList.h>
#include "AliMultVariable.h"

class AliMultInput : public TNamed {

public:
    AliMultInput();
    AliMultInput(const char * name, const char * title = "MultInput");
    AliMultInput(const AliMultInput& o);
    AliMultInput& operator=(const AliMultInput& o);
    virtual ~AliMultInput();

    void     AddVariable ( AliMultVariable *lVar );
    AliMultVariable* GetVariable (const TString& lName) const;
    AliMultVariable* GetVariable (Long_t iIdx) const;
    Long_t GetNVariables         () const { return fVariableList ? fVariableList->GetEntries() : 0; }
    void Clear(Option_t* option="");
    void Set(const AliMultInput* other);
    void Print(Option_t* option="") const;

    Bool_t SetValue(TString name, Float_t val);
    Bool_t SetValue(TString name, Int_t val);
    Bool_t IncrementValue(TString name, Int_t increment=1);

    Int_t   GetValueInteger(TString name) const;
    Float_t GetValue(TString name) const;

 private:
    TList *fVariableList; //List containing all AliMultVariables

    ClassDef(AliMultInput, 2);
};
#endif
