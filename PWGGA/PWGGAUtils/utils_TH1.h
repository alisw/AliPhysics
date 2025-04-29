#pragma once

class TF1;
class TH1;

#include <string>


// (helper) class for exponential extrapolation of TH1 histograms
class TH1_ExponentialInterpolation
{
    public:    
    TH1_ExponentialInterpolation(std::string const &_id,
                                 TH1 const &theTH1,
                                 bool theIntegrate,
                                 bool theUseXtimesExp);
    ~TH1_ExponentialInterpolation();

    TF1 &GetNewTF1_global(std::string const &_id);
    
    double Evaluate(double *x, double *);

    private:
    std::string id;
    TH1 *th1; // this will point to a clone of the histo to be extrapolated    
    bool integrate;  
    bool useXtimesExp;  // fits a function of the f(x) = x * exp([0] + [1]*x)

    TF1 *fCache;
};

class utils_TH1
{
public:


    // get a TF1 exponential defined by:
    /* 
    case theIntegrate = false: 
        - the next two neares bin centers to theX
        - the two coefficients are calculated analytically from the two bin content values
    case theIntegrate = true: 
        - the lower and upper edges of the two nearest bins to theX 
        - the two coefficients are determined from a fit using the two bins  
    */
    static TF1 *GetLocalExponentialTF1(TH1   &theTH1, 
                                       double theX, 
                                       bool   theIntegrate, 
                                       bool   theUseXtimesExp);
    
    static double LocalExponentialInterpolate(TH1   &theTH1, 
                                              double theX, 
                                              bool   theIntegrate = false,
                                              bool   theUseXtimesExp = false);
    /* 
    get globally defined TF1 which is either fitted to the bin centers 
    OR such that the integrals in each in bin agree with their content. 
    if theUseXtimesExp=true: use f(x) = x * exp([0] + [1]*x). Otherwise f(x) = exp([0] + [1]*x) 
    */
    static TF1 &GlobalPieceWiseExponentialInterpolation(std::string const &theName, 
                                                        TH1 const         &theTH1, 
                                                        bool               theIntegrate = false,
                                                        bool               theUseXtimesExp = false);
};