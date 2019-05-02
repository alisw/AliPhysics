#include "TMath.h"
// ExpoGaussExpo function
// to use:

/*
.L gauspl_function.C++g
Double_t xmin = 0;
Double_t xmax = 1;
TF1 *gauspl = new TF1("gauspl", gauspl_function, xmin, xmax, 5, 1);
gauspl = new TF1("gauspl", gauspl_function, xmin, xmax, 5, 1);
gauspl->SetParNames("Constant", "Mean", "Sigma", "xTailStart", "tailsize");
gauspl->SetParameters(1,0.05,0.05,0.1,0.1);
gauspl->SetParLimits(0,0,1e38);
gauspl->SetParLimits(1,0,1e38);
gauspl->SetParLimits(2,0,1e38);
gauspl->SetParLimits(3,0,1e38);
gauspl->SetParLimits(4,0,1);
gauspl->SetNpx(10000)
gauspl->Draw();
*/


Double_t gauspl_function(Double_t* xv, Double_t *p) {
    //p[0] normalization gaus
    //p[1] gaus mean
    //p[2] gaus sigma
    //p[3] start tail here (x)
    //p[4] relative tail size 0=no tail
    Double_t x = xv[0];    
    Double_t c = p[0];
    Double_t mu = p[1];
    Double_t sigma = p[2];
    Double_t x0 = p[3];
    Double_t xn = (x-mu)/sigma;    
    if (x < x0) return c*TMath::Exp(-0.5*(xn*xn));
    Double_t s = p[4];    
    Double_t x0n = (mu-x0)/sigma;    
    Double_t n = -x0*(mu-x0)/(sigma*sigma);
    return c*s*TMath::Exp(-0.5*(x0n*x0n))*TMath::Power(x0,n)*TMath::Power(x,-n) + c*(1.-s)*TMath::Exp(-0.5*(xn*xn));
}
