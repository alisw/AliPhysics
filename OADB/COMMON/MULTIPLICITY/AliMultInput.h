#ifndef AliMultInput_H
#define AliMultInput_H

#include <iostream>
#include "TNamed.h"

using namespace std;

class AliMultInput : public TNamed {
    
public:
    AliMultInput();
    AliMultInput(const char * name, const char * title = "MultInput");
    ~AliMultInput();

    void     AddVariable ( AliMultVariable *lVar ) { fVariableList->Add(lVar); fNVars++; }
    AliMultVariable* GetVariable (TString lName) { return ((AliMultVariable*)fVariableList->FindObject(lName.Data())); }
    AliMultVariable* GetVariable (Long_t iIdx) { return ((AliMultVariable*)fVariableList->At(iIdx)); }
    Long_t GetNVariables () { return fNVars; }
    
private:
    Long_t fNVars;
    TList *fVariableList; //List containing all AliMultVariables
    
    ClassDef(AliMultInput, 1)
};
#endif
