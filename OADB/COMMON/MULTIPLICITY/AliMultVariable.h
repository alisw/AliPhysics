#ifndef AliMultVariable_H
#define AliMultVariable_H

#include <iostream>
#include "TNamed.h"

using namespace std;

class AliMultVariable : public TNamed {
    
public:
    AliMultVariable();
    AliMultVariable(const char * name, const char * title = "Mult Variable");
    ~AliMultVariable();
    
    void     SetValue ( Float_t lVal ) { fValue = lVal; }
    Float_t GetValue ()  { return fValue; }
    //(specialized) use with care
    Float_t& GetRValue () { return fValue; }
    
    void     SetValueInteger ( Int_t lVal ) { fValueInteger = lVal; }
    Int_t GetValueInteger () { return fValueInteger; }
    //(specialized) use with care
    Int_t& GetRValueInteger () { return fValueInteger; }
    
    void     SetMean ( Float_t lVal ) { fMean = lVal; }
    Float_t GetMean () { return fMean; }
    
    void     SetIsInteger( Bool_t lVal ){ fIsInteger = lVal; }
    Bool_t  IsInteger() { return fIsInteger; } 
    
private:
    Bool_t fIsInteger; //If Integer
    Float_t fValue;    //Variable value
    Int_t   fValueInteger; //Variable value if integer
    Float_t fMean;     //Variable mean
    
    ClassDef(AliMultVariable, 1)
};
#endif
