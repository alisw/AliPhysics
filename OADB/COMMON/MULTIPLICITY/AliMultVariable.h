#ifndef AliMultVariable_H
#define AliMultVariable_H
#include <TNamed.h>

class AliMultVariable : public TNamed {
    
public:
    AliMultVariable();
    AliMultVariable(const char * name, const char * title = "Mult Variable");
    AliMultVariable(const AliMultVariable& o)
      : TNamed(o),
        fIsInteger(o.fIsInteger),
        fValue(o.fValue),
        fValueInteger(o.fValueInteger),
        fMean(o.fMean)
    {}
    AliMultVariable& operator=(const AliMultVariable& o)
    {
        if (&o == this) return *this;
        SetName(o.GetName());
        SetTitle(o.GetTitle());
        fIsInteger = o.fIsInteger;
        fValue = o.fValue;
        fValueInteger = o.fValueInteger;
        fMean = o.fMean;
        return *this;
    }
    void Set(const AliMultVariable* other)
    {
        if (!other) { if (IsInteger()) SetValueInteger(0); else SetValue(0.); }
        if (IsInteger())
            SetValueInteger(other->GetValueInteger());
        else
            SetValue(other->GetValue());
    }
    ~AliMultVariable();
    
    void Clear(Option_t*)
    {
        fValue = 0;
        fValueInteger = 0;
    }
    
    void     SetValue ( Float_t lVal ) { fValue = lVal; }
    Float_t GetValue () const  { return fValue; }
    //(specialized) use with care
    Float_t& GetRValue () { return fValue; }
    
    void     SetValueInteger ( Int_t lVal ) { fValueInteger = lVal; }
    Int_t GetValueInteger () const { return fValueInteger; }
    //(specialized) use with care
    Int_t& GetRValueInteger () { return fValueInteger; }
    
    void     SetMean ( Float_t lVal ) { fMean = lVal; }
    Float_t GetMean () const { return fMean; }
    
    void     SetIsInteger( Bool_t lVal ){ fIsInteger = lVal; }
    Bool_t  IsInteger() const { return fIsInteger; }
    
    void Print(Option_t* option="") const;
    
private:
    Bool_t fIsInteger; //If Integer
    Float_t fValue;    //Variable value
    Int_t   fValueInteger; //Variable value if integer
    Float_t fMean;     //Variable mean
    
    ClassDef(AliMultVariable, 1)
};
#endif
