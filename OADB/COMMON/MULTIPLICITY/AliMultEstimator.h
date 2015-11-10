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
    
    void     SetUseAnchor( Bool_t lVal ) { fkUseAnchor = lVal; }
    Bool_t   GetUseAnchor() const { return fkUseAnchor; }

    void     SetAnchorPoint( Float_t lVal ){ fAnchorPoint = lVal; }
    Float_t  GetAnchorPoint() const { return fAnchorPoint; }

    void     SetAnchorPercentile( Float_t lVal ){ fAnchorPercentile = lVal; }
    Float_t  GetAnchorPercentile() const { return fAnchorPercentile; }
    
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
    
    //Anchor point definition
    Bool_t  fkUseAnchor;        //Use Anchor Logic (default: No)
    Float_t fAnchorPoint;       //Raw value below which
    Float_t fAnchorPercentile;  //Percentile of X-section at anchor point
    
    ClassDef(AliMultEstimator, 1)
};
#endif
