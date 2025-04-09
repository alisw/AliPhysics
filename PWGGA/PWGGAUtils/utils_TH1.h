#pragma once

class TF1;
class TH1;

#include <string>

// (helper) class for exponential extrapolation of TH1 histograms
class TH1_ExponentialInterpolation
{
    public:    
    TH1_ExponentialInterpolation(std::string const &_id,
                             TH1 const &theTH1);
    ~TH1_ExponentialInterpolation();

    TF1 &GetNewTF1(std::string const &_id);
    double Evaluate(double *x, double *);

    private:
    std::string id;
    TH1 *th1; // this will point to a clone of the histo to be extrapolated        
};

class utils_TH1
{
public:

    // get a TF1 exponential defined by theTH1 adjacent bin centers to theX
    static TF1 *GetLocalExponentialTF1(TH1 &theTH1, double theX);
    
    // returns f(theX) with f defined by the exponential that fits theTH1's adjacent bin centers to theX
    static double LocalExponentialInterpolate(TH1 &theTH1, double theX);
    
    static TF1 &GlobalPieceWiseExponentialInterpolation(std::string const &theName, TH1 const &theTH1);
};