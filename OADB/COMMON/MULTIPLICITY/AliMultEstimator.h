#ifndef AliMultEstimator_H
#define AliMultEstimator_H
#include <TNamed.h>
class AliMultInput;
class TFormula;

class AliMultEstimator : public TNamed {
    
public:
    AliMultEstimator();
    AliMultEstimator(const char * name, const char * title = "Mult Estimator", TString lInitDef = "");
    AliMultEstimator(const AliMultEstimator& e);
    AliMultEstimator& operator=(const AliMultEstimator& e);
    void Set(const AliMultEstimator* other);
    ~AliMultEstimator();

    void     SetDefinition ( const TString& lVal ) { fDefinition = lVal; }
    const TString&  GetDefinition () const { return fDefinition; }

    void     SetIsInteger ( Bool_t lVal ) { fIsInteger = lVal; }
    Bool_t   IsInteger() const { return fIsInteger; }
    
    void     SetPercentile ( Float_t lVal ) { fPercentile = lVal; }
    Float_t GetPercentile () const { return fPercentile; }
    
    void     SetValue ( Float_t lVal ) { fValue = lVal; }
    Float_t GetValue () const { return fValue; }
    
    void     SetMean ( Float_t lVal ) { fMean = lVal; }
    Float_t GetMean () const { return fMean; }
    
    Float_t GetZ () const; //check for zero

    //Pre-processing for speed
    void SetupFormula(const AliMultInput* lInput);
    Float_t Evaluate(const AliMultInput* lInput);
    
private:
    TString fDefinition; //How to evaluate based on AliMultVariables
    Bool_t fIsInteger; //Requires special treatment when calibrating
    
    Float_t fValue;     // estimator value
    Float_t fMean;   // estimator mean value
    Float_t fPercentile;   //Percentile
    TFormula* fFormula; //!
    
    ClassDef(AliMultEstimator, 1)
};
#endif
