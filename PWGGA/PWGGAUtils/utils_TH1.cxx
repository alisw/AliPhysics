#include "utils_TH1.h"

#include "TF1.h"
#include "TH1.h"
#include <iostream>


// ============== (helper) class for exponential extrapolation of TH1 histograms===================
//_________________________________________________________________________________________________
TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(std::string const &_id,
                                                           TH1 const         &theTH1,
                                                           bool               theIntegrate,
                                                           bool               theUseXtimesExp)
    : id(_id), 
      th1(dynamic_cast<TH1*>(theTH1.Clone(Form("%s_clone", theTH1.GetName())))),
      integrate(theIntegrate),
      useXtimesExp(theUseXtimesExp),
      fCache(nullptr)
{
    printf("TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(): created instance with properties:\n\t"
           "id: %s, th1: %s, integrate: %d, useXtimesExp: %d\n",
           id.data(), th1->GetName(), integrate, useXtimesExp);
}

//_________________________________________________________________________________________________
TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation()
{
    printf("TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation(): deleting instance %s\n", 
           id.data());
    delete th1;
}

//_________________________________________________________________________________________________
TF1 &TH1_ExponentialInterpolation::GetNewTF1_global(std::string const &theName)
{
    printf("TH1_ExponentialInterpolation::GetNewTF1_global(): id: %s returning new TF1 with name '%s'.\n",
           id.data(), theName.data());
    TAxis &lXaxis = *th1->GetXaxis();
    return *new TF1(theName.data(),
                    this,
                    &TH1_ExponentialInterpolation::Evaluate,
                    lXaxis.GetXmin(), lXaxis.GetXmax(),
                    0, // nPar
                    "TH1_ExponentialInterpolation",
                    "Evaluate");
}

//_________________________________________________________________________________________________
double TH1_ExponentialInterpolation::Evaluate(double *x, double *)
{
    // check if there is already a local interpol for this x
    bool canUseExisting = fCache && (*x>=fCache->GetXmin() && *x<=fCache->GetXmax());
    if (!canUseExisting){
        // printf("INFO: TH1_ExponentialInterpolation::Evaluate(): Called with x = %f.\n%s Will create a new one.\n",
        //        *x,
        //        fCache 
        //           ? "Cached function available but with wrong range."
        //           : "No cached function available.");
    }
    TF1 *localFunction = canUseExisting 
      ?  fCache
      :  utils_TH1::GetLocalExponentialTF1(*th1, *x, integrate, useXtimesExp);
    double result = localFunction ? localFunction->Eval(*x) : 0.;
    fCache = localFunction;

    // delete localFunction;
    return result;
}
///////////////////// end TH1_ExponentialInterpolation ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//======================= class utils_TH1==========================================================
//_________________________________________________________________________________________________
TF1 *utils_TH1::GetLocalExponentialTF1(TH1   &theTH1, 
                                       double theX, 
                                       bool   theIntegrate, 
                                       bool   theUseXtimesExp){

    TAxis const &lAxis = *theTH1.GetXaxis();    
    auto findFunctionI = [&](){
        int histoBinI = lAxis.FindBin(theX);
        return (theX >= lAxis.GetBinCenter(histoBinI)) 
            ? histoBinI 
            : std::max(histoBinI-1, 0);
    };
    int const iLeftBin = findFunctionI();
    int const iRightBin = iLeftBin+1;

    std::pair<double, double> lRange_edgeToEdge( 
        { lAxis.GetBinLowEdge(iLeftBin), lAxis.GetBinUpEdge(iRightBin) });    
    
    std::pair<double, double> lRange_centerToCenter( 
        { lAxis.GetBinCenter(iLeftBin), lAxis.GetBinCenter(iRightBin) });
    
    std::pair<double, double> const &lFunctionDefineRange = theIntegrate 
        ? lRange_edgeToEdge
        : lRange_centerToCenter;
    
    std::string lFunctionName(Form("%s_localExponential%s%s", 
                                   theTH1.GetName(), 
                                   theUseXtimesExp 
                                     ? "_*x" 
                                     : "",
                                   theIntegrate 
                                     ? "fitted_w/_int_cond" 
                                     : "calc_analyt_through_bin_centers"));
 
    TF1 *lResult = new TF1(lFunctionName.data(),
                           Form("%sexpo(0)", 
                                theUseXtimesExp ? "x*" : ""),  // = [x*] exp([0] + [1]*x) 
                           lFunctionDefineRange.first,
                           lFunctionDefineRange.second);
    
    std::string lFitOptions("QFMN0"); // Q = minimum printing
    if (theIntegrate){
        theTH1.Fit(lResult, 
                   lFitOptions.append(theIntegrate ? "I" : "").data(), 
                    "" /* global fit options I believe */, 
                    lFunctionDefineRange.first, 
                    lFunctionDefineRange.second);
        lResult->SetRange(lRange_centerToCenter.first, lRange_centerToCenter.second);
    } 
    else {
        double x1 = lRange_centerToCenter.first;
        double x2 = lRange_centerToCenter.second;
        double y1 = theTH1.GetBinContent(iLeftBin);
        double y2 = theTH1.GetBinContent(iRightBin);

        if (!(y1&&y2)) { 
            return nullptr;
        }
        // printf("using bin %d %f %f and bin %d %f %f\n", bin, x1, y1, adjBin, x2, y2);
        double b = TMath::Log(y2/y1) / (x2-x1);
        double a = TMath::Log(y1) - b*x1;
        lResult->SetParameters(a,b);
    }
    return lResult;
}

//_________________________________________________________________________________________________
double utils_TH1::LocalExponentialInterpolate(TH1   &theTH1, 
                                              double theX, 
                                              bool   theIntegrate /*= false*/,
                                              bool   theUseXtimesExp /*= false*/){
    TF1* f = GetLocalExponentialTF1(theTH1, theX, theIntegrate, theUseXtimesExp);
    return f ? f->Eval(theX) : 0.;
}

//_________________________________________________________________________________________________
TF1 &utils_TH1::GlobalPieceWiseExponentialInterpolation(std::string const &theName, 
                                                        TH1 const &theTH1, 
                                                        bool theIntegrate /* = false*/,
                                                        bool theUseXtimesExp /* = false*/) // fits a function of the f(x) = x * exp([0] + [1]*x)
{
    printf("utils_TF1::GlobalPieceWiseExponentialInterpolation(): called with theName: %s, theTH1: %s, theIntegrate = %d, theUseXtimesExp = %d\n",
           theName.data(), theTH1.GetName(), theIntegrate, theUseXtimesExp);

    class TH1_ExponentialInterpolation &lGlobalInterpolation = 
        *new class TH1_ExponentialInterpolation(theName, theTH1, theIntegrate, theUseXtimesExp);
    return lGlobalInterpolation.GetNewTF1_global(theName);
}
// end utils_TH1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

