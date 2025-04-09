#include "utils_TH1.h"

#include "TF1.h"
#include "TH1.h"
#include <iostream>

// ============== class TH1_ExponentialInterpolation ===
TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(std::string const &_id,
                                                           TH1 const         &theTH1)
    : id(_id), 
      th1(dynamic_cast<TH1*>(theTH1.Clone(Form("%s_clone", theTH1.GetName()))))
{
    printf("TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(): created instance with properties:\n\t"
        "id: %s, th1: %s\n",
        id.data(), th1->GetName());
}

TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation()
{
    printf("TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation(): deleting instance %s\n", 
           id.data());
    delete th1;
}
    
TF1 &TH1_ExponentialInterpolation::GetNewTF1(std::string const &theName)
{
    printf("TH1_ExponentialInterpolation::GetNewTF1(): id: %s returning new TF1 with name '%s'.\n",
           id.data(), theName.data());
    TAxis &lXaxis = *th1->GetXaxis();
    return *new TF1(theName.data(),
                    this,
                    &TH1_ExponentialInterpolation::Evaluate,
                    lXaxis.GetXmin(), lXaxis.GetXmax(),
                    0,
                    "TH1_ExponentialInterpolation",
                    "Evaluate");
}


double TH1_ExponentialInterpolation::Evaluate(double *x, double *)
{
    TF1 *localFunction = utils_TH1::GetLocalExponentialTF1(*th1, *x);
    double result = localFunction ? localFunction->Eval(*x) : 0.;
    delete localFunction;
    return result;
}
///////////////////// end TH1_ExponentialInterpolation /////////////
////////////////////////////////////////////////////////////////////

TF1 *utils_TH1::GetLocalExponentialTF1(TH1 &theTH1, double theX){
        
    // find two non-zero bins to use for interpolation
    int bin = theTH1.FindBin(theX);
    float y1 = theTH1.GetBinContent(bin);
    if (!y1) { 
        printf("INFO: GetExponentialFunction(): %s at x=%f, bin=%d, has bin content 0. returning nullptr\n",
            theTH1.GetName(), theX, bin);
        return nullptr;
    }
    
    float x1 = theTH1.GetBinCenter(bin);
    int nBins = theTH1.GetNbinsX();
    int adjBin = (bin==1) ? 2 
                        : bin==nBins ? nBins-1 
                                        : theX<=x1 ? bin-1 
                                                : bin+1;
    
    float x2 = theTH1.GetBinCenter(adjBin);
    float y2 = theTH1.GetBinContent(adjBin);
    if (!y2){
        adjBin = adjBin<bin ? bin+1 : bin-1;
        if (!adjBin || adjBin > nBins) {
            std::cout << "ERROR: GetExponentialFunction(): adjBin out of binrange\n";
            return nullptr;
        }
        x2 = theTH1.GetBinCenter(adjBin);
        y2 = theTH1.GetBinContent(adjBin);
        if (!y2) {
            std::cout << "y2 0.\n";
            return nullptr;
        }
    }
    // printf("using bin %d %f %f and bin %d %f %f\n", bin, x1, y1, adjBin, x2, y2);
    
    float b = TMath::Log(y2/y1) / (x2-x1);
    float a = TMath::Log(y1) - b*x1;
    
    // find out where function needs to live
    double xmin = std::min(x1, x2);
    xmin = std::min(xmin, theX);
    double xmax = std::max(x1, x2);
    xmax = std::max(xmax, theX);
    
    TF1* ret = new TF1("tf1", "expo(0)", xmin, xmax); 
    ret->SetParameters(a,b);
    return ret;
}

double utils_TH1::LocalExponentialInterpolate(TH1 &theTH1, double theX){
    TF1* f = GetLocalExponentialTF1(theTH1, theX);
    return f ? f->Eval(theX) : 0.;
}

// =================================================================================================
TF1 &utils_TH1::GlobalPieceWiseExponentialInterpolation(std::string const &theName, TH1 const &theTH1)
{
    printf("utils_TF1::GlobalPieceWiseExponentialInterpolation(): called with theName: %s, theTH1: %s\n",
           theName.data(), theTH1.GetName());

    class TH1_ExponentialInterpolation &lGlobalInterpolation = *new class TH1_ExponentialInterpolation(theName, theTH1);
    return lGlobalInterpolation.GetNewTF1(theName);
}
// end utils_TH1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

