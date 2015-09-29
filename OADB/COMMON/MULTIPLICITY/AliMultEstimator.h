#ifndef AliMultEstimator_H
#define AliMultEstimator_H

#include <iostream>
#include "TNamed.h"

using namespace std;

class AliMultEstimator : public TNamed {
    
public:
    AliMultEstimator();
    AliMultEstimator(const char * name, const char * title = "Mult Estimator", TString lInitDef = "");
    ~AliMultEstimator();

    void     SetDefinition ( TString lVal ) { fDefinition = lVal; }
    TString  GetDefinition () { return fDefinition; }

    void     SetIsInteger ( Bool_t lVal ) { fIsInteger = lVal; }
    Bool_t   IsInteger () { return fIsInteger; }
    
    void     SetPercentile ( Float_t lVal ) { fPercentile = lVal; }
    Float_t GetPercentile () { return fPercentile; }
    
    void     SetValue ( Float_t lVal ) { fValue = lVal; }
    Float_t GetValue () { return fValue; }
    
    void     SetMean ( Float_t lVal ) { fMean = lVal; }
    Float_t GetMean () { return fMean; }
    
    Float_t GetZ (); //check for zero

private:
    TString fDefinition; //How to evaluate based on AliMultVariables
    Bool_t fIsInteger; //Requires special treatment when calibrating
    
    Float_t fValue;     // estimator value
    Float_t fMean;   // estimator mean value
    Float_t fPercentile;   //Percentile
    
    ClassDef(AliMultEstimator, 1)
};
#endif
