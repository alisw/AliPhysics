#ifndef AliMultInput_H
#define AliMultInput_H
#include <TNamed.h>
#include "AliMultVariable.h"

class AliMultInput : public TNamed {
    
public:
    AliMultInput();
    AliMultInput(const char * name, const char * title = "MultInput");
    AliMultInput(const AliMultInput& o);
    AliMultInput& operator=(const AliMultInput& o);
    ~AliMultInput();

    void     AddVariable ( AliMultVariable *lVar );
    AliMultVariable* GetVariable (const TString& lName) const;
    AliMultVariable* GetVariable (Long_t iIdx) const;
    Long_t GetNVariables         () const { return fNVars; }
    void Clear(Option_t* option="");
    void Set(const AliMultInput* other);
    void Print(Option_t* option="") const;
    
private:
    Long_t fNVars;
    TList *fVariableList; //List containing all AliMultVariables
    
    ClassDef(AliMultInput, 1)
};
#endif
