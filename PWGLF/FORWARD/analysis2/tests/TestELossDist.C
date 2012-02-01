#include <iomanip>
#ifndef __CINT__
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TArrayD.h>
#include <TStyle.h>
// #include <TFitResult.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <THStack.h>
#include <TList.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#else
class TH1;
class TH2;
class TArrayD;
class TRandom;
class TCanvas;
class TF1;
#endif

//___________________________________________________________________
static Double_t landauGaus1(Double_t* xp, Double_t* pp);
static Double_t landauGausN(Double_t* xp, Double_t* pp);
static Double_t landauGausI(Double_t* xp, Double_t* pp);

//====================================================================
/**
 * 
 * 
 * @ingroup pwg2_forward_scripts_tests
 */
struct Function
{
  /**
   * Enumeration of parameters 
   * 
   */
  enum { 
    /// Constant 
    kC, 
    /// MPV of landau
    kDelta, 
    /// Width of landau
    kXi, 
    /// Sigma of gaussian
    kSigma, 
    /// Number of particles
    kN, 
    /// Weights 
    kA 
  };
  /// MP shift
  static const Double_t fgkMPShift;
  /// 1 / sqrt(2 * pi)
  static const Double_t fgkInvRoot2Pi;
  /// Number of sigmas to integrate over
  static const Double_t fgkConvNSigma;
  /// Number of steps in integration 
  static const Double_t fgkConvNSteps;
  /** 
   * Evaluate shifted landay 
   * @f[                                                                        
   *    f'_{L}(x;\Delta,\xi) = f_L(x;\Delta+0.22278298\xi)                      
   * @f]                                                                        
   *                                                                            
   * where @f$ f_{L}@f$ is the ROOT implementation of the Landau                
   * distribution (known to have @f$ \Delta_{p}=-0.22278298@f$ for              
   * @f$\Delta=0,\xi=1@f$.                                                      
   *                                                                            
   * @param x      Where to evaluate @f$ f'_{L}@f$                              
   * @param delta  Most probable value                                          
   * @param xi     The 'width' of the distribution                              
   *                                                                            
   * @return @f$ f'_{L}(x;\Delta,\xi) @f$                                       
   */
  static Double_t  Landau(Double_t x, Double_t delta, Double_t xi)
  {
    return TMath::Landau(x, delta - xi * fgkMPShift, xi);
  }
  /** 
   * Calculate the value of a Landau convolved with a Gaussian                  
   *                                                                            
   * @f[                                                                        
   * f(x;\Delta,\xi,\sigma') = \frac{1}{\sigma' \sqrt{2 \pi}}                   
   *    \int_{-\infty}^{+\infty} d\Delta' f'_{L}(x;\Delta',\xi)                 
   *    \exp{-\frac{(\Delta-\Delta')^2}{2\sigma^2}}                            
   * @f]                                                                        
   *                                                                            
   * where @f$ f'_{L}@f$ is the Landau distribution, @f$ \Delta@f$ the
   * energy loss, @f$ \xi@f$ the width of the Landau, and @f$\sigma@f$
   * is the variance of the Gaussian.
   *                                                                            
   * Note that this function uses the constants fgkConvNSteps and          
   * fgkConvNSigma                                                        
   *                                                                            
   * References: 
   *  - <a href="dx.doi.org/10.1016/0168-583X(84)90472-5">
   *   Nucl.Instrum.Meth.B1:16</a>
   *  - <a href="dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
   *  - <a href="root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">
   *   ROOT implementation</a> 
   *                                                                            
   * @param x         where to evaluate @f$ f@f$                                
   * @param delta     @f$ \Delta@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$           
   * @param xi        @f$ \xi@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$  
   * @param sigma     @f$ \sigma@f$ of Gaussian
   *                                                                            
   * @return @f$ f@f$ evaluated at @f$ x@f$.                                    
   */
  static Double_t  LandauGaus(Double_t x,  Double_t delta, 
			      Double_t xi, Double_t sigma) 
  {
    Double_t deltap = delta - xi * fgkMPShift;
    Double_t xlow   = x - fgkConvNSigma * sigma;
    Double_t xhigh  = x + fgkConvNSigma * sigma;
    Double_t step   = (xhigh - xlow) / fgkConvNSteps;
    Double_t sum    = 0;

    for (Int_t i = 0; i < fgkConvNSteps / 2; i++) { 
      Double_t x1 = xlow  + (i - .5) * step;
      Double_t x2 = xhigh - (i - .5) * step;

      Double_t s1 = 
	TMath::Landau(x1, deltap, xi, kTRUE) * TMath::Gaus(x, x1, sigma);
      Double_t s2 = 
	TMath::Landau(x2, deltap, xi, kTRUE) * TMath::Gaus(x, x2, sigma);
      sum += s1 + s2;
    }
    Double_t ret = step * sum * fgkInvRoot2Pi / sigma;
    return ret;
  }
  /** 
   * Evaluate                                                                   
   * @f[                                                                        
   *    f_i(x;\Delta,\xi,\sigma) = f(x;\Delta_i,\xi_i,\sigma_i)               
   * @f]                                                                        
   * corresponding to @f$ i@f$ particles i.e., with the substitutions           
   * @f{eqnarray*}{                                                             
   *    \Delta    \rightarrow \Delta_i    &=& i(\Delta + \xi\log(i))\\          
   *    \xi       \rightarrow \xi_i       &=& i \xi\\                           
   *    \sigma    \rightarrow \sigma_i    &=& \sqrt{i}\sigma\\
   * @f}                                                                        
   *                                                                            
   * @param x        Where to evaluate                                          
   * @param delta    @f$ \Delta@f$                                              
   * @param xi       @f$ \xi@f$                                                 
   * @param sigma    @f$ \sigma@f$                                              
   * @param i        @f$ i @f$                                                  
   *                                                                            
   * @return @f$ f_i @f$ evaluated        * 
   */
  static Double_t ILandauGaus(Int_t i, Double_t x, Double_t delta, 
			      Double_t xi, Double_t sigma)
  {
    if (i <= 1) return LandauGaus(x, delta, xi, sigma);
    Double_t di      = i;
    Double_t deltai  = i * (delta + xi * TMath::Log(di));
    Double_t xii     = i * xi;
    Double_t sigmai  = TMath::Sqrt(di) * sigma;
    
    if (sigmai < 1e-10) return Landau(x, deltai, xii);
    Double_t ret =  LandauGaus(x, deltai, xii, sigmai);;
    // Info("ILandauGaus", "Fi(%f;%d,%f,%f,%f)->%f",
    //      x, i, deltai, xii, sigmai, ret);
    return ret;
  }
  /** 
   * Numerically evaluate 
   * @f[ 
   *    \left.\frac{\partial f_i}{\partial p_i}\right|_{x}
   * @f] 
   * where @f$ p_i@f$ is the @f$ i^{\mbox{th}}@f$ parameter.  The mapping 
   * of the parameters is given by 
   *
   * - 0: @f$\Delta@f$ 
   * - 1: @f$\xi@f$ 
   * - 2: @f$\sigma@f$ 
   *
   * This is the partial derivative with respect to the parameter of
   * the response function corresponding to @f$ i@f$ particles i.e.,
   * with the substitutions
   * @f[ 
   *    \Delta    \rightarrow \Delta_i    = i(\Delta + \xi\log(i))\\
   *    \xi       \rightarrow \xi_i       = i \xi\\
   *    \sigma    \rightarrow \sigma_i    = \sqrt{i}\sigma\\
   * @f] 
   * 
   * @param x        Where to evaluate 
   * @param ipar     Parameter number 
   * @param dp       @f$ \epsilon\delta p_i@f$ for some value of @f$\epsilon@f$
   * @param delta    @f$ \Delta@f$ 
   * @param xi       @f$ \xi@f$ 
   * @param sigma    @f$ \sigma@f$ 
   * @param i        @f$ i@f$
   * 
   * @return @f$ f_i@f$ evaluated
   */  
  static Double_t IdLandauGausdPar(Int_t    i,     Double_t x, 
				   UShort_t ipar,  Double_t dPar, 
				   Double_t delta, Double_t xi, 
				   Double_t sigma)
  {
    if (dPar == 0) return 0;
    Double_t dp      = dPar;
    Double_t d2      = dPar / 2;
    Double_t deltaI  =  i * (delta + xi * TMath::Log(i));
    Double_t xiI     =  i * xi;
    Double_t si      =  TMath::Sqrt(Double_t(i));
    Double_t sigmaI  =  si*sigma;
    Double_t y1      = 0;
    Double_t y2      = 0;
    Double_t y3      = 0;
    Double_t y4      = 0;
    switch (ipar) {
    case 0: 
      y1 = ILandauGaus(i, x, deltaI+i*dp, xiI, sigmaI);
      y2 = ILandauGaus(i, x, deltaI+i*d2, xiI, sigmaI);
      y3 = ILandauGaus(i, x, deltaI-i*d2, xiI, sigmaI);
      y4 = ILandauGaus(i, x, deltaI-i*dp, xiI, sigmaI);
      break;
    case 1: 
      y1 = ILandauGaus(i, x, deltaI, xiI+i*dp, sigmaI);
      y2 = ILandauGaus(i, x, deltaI, xiI+i*d2, sigmaI);
      y3 = ILandauGaus(i, x, deltaI, xiI-i*d2, sigmaI);
      y4 = ILandauGaus(i, x, deltaI, xiI-i*dp, sigmaI);
      break;
    case 2: 
      y1 = ILandauGaus(i, x, deltaI, xiI, sigmaI+si*dp);
      y2 = ILandauGaus(i, x, deltaI, xiI, sigmaI+si*d2);
      y3 = ILandauGaus(i, x, deltaI, xiI, sigmaI-si*d2);
      y4 = ILandauGaus(i, x, deltaI, xiI, sigmaI-si*dp);
      break;
    default:
      return 0;
    } 
    
    Double_t d0  = y1 - y4;
    Double_t d1  = 2 * (y2 - y3);
    
    Double_t g   = 1/(2*dp) * (4*d1 - d0) / 3;
    
    return g;
  }
  /** 
   * Evaluate                                                                   
   * @f[                                                                        
   *   f_N(x;\Delta,\xi,\sigma') = \sum_{i=1}^N a_i f_i(x;\Delta,\xi,\sigma'a)  
   * @f]                                                                        
   *                                                                            
   * where @f$ f(x;\Delta,\xi,\sigma')@f$ is the convolution of a               
   * Landau with a Gaussian (see LandauGaus).  Note that                        
   * @f$ a_1 = 1@f$, @f$\Delta_i = i(\Delta_1 + \xi\log(i))@f$,                 
   * @f$\xi_i=i\xi_1@f$, and @f$\sigma_i'^2 = \sigma_n^2 + i\sigma_1^2@f$.      
   *                                                                            
   * References:                                                                
   *  - <a href="dx.doi.org/10.1016/0168-583X(84)90472-5">
   *    Nucl.Instrum.Meth.B1:16</a>
   *  - <a href="dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
   * @param x        Where to evaluate @f$ f_N@f$                               
   * @param delta    @f$ \Delta_1@f$                                            
   * @param xi       @f$ \xi_1@f$                                               
   * @param sigma    @f$ \sigma_1@f$                                            
   * @param n        @f$ N@f$ in the sum above.                                 
   * @param a        Array of size @f$ N-1@f$ of the weights @f$ a_i@f$ for     
   *                 @f$ i > 1@f$                                               
   *                                                                            
   * @return @f$ f_N(x;\Delta,\xi,\sigma)@f$                                   
   */  
  static Double_t NLandauGaus(Double_t x, Double_t delta, Double_t xi, 
			      Double_t sigma, Int_t n, Double_t* a)
  {
    Double_t res = LandauGaus(x, delta, xi, sigma);
    for (Int_t i = 2; i <= n; i++) 
      res += a[i-2] * ILandauGaus(i, x, delta, xi, sigma);
    return res;
  }
  /** 
   * Calculate the the estimate of number of particle contribtutions
   * at energy loss @a x 
   * @f[
   *   E(\Delta;\Delta_p,\xi,\sigma,\mathbf{a}) = 
   *   \frac{\sum_{i=1}^N i a_i f_i(x;\Delta_p,\xi,\sigma)}{
   *     \sum_{i=1}^N a_i f_i(x;\Delta_p,\xi,\sigma)}
   * @f]
   * where @f$a_1=1@f$ 
   * 
   * @param x      Energy loss @f$\Delta@f$
   * @param delta  @f$\Delta_{p}@f$
   * @param xi     @f$\xi@f$ 
   * @param sigma  @f$\sigma@f$
   * @param n      Maximum number of particles @f$N@f$ 
   * @param a      Weights @f$a_i@f$ for @f$i=2,...,N@f$ 
   * 
   * @return @f$E(\Delta;\Delta_p,\xi,\sigma,\mathbf{a})@f$
   */
  Double_t NEstimate(Double_t x, Double_t delta, Double_t xi, 
		     Double_t sigma, Int_t n, Double_t* a)
  {
    Double_t num = LandauGaus(x, delta, xi, sigma);
    Double_t den = num;
    for (Int_t i = 2; i <= n; i++) {
      Double_t f =  ILandauGaus(i, x, delta, xi, sigma);
      num        += a[i-2] * f * i;
      den        += a[i-2] * f;
    }
    if (den < 1e-4) return 0;
    return num / den;
  }
  /** 
   * Calculate the partial derivative of the function 
   * @f$E(\Delta;\Delta_p,\xi,\sigma,\mathbf{a})@f$ with respect to
   * one of the parameters @f$\Delta_{p}@f$, @f$\xi@f$, @f$\sigma@f$,
   * or @f$a_i@f$.  Note for that 
   * @f{eqnarray*}
   *   \frac{\partial E}{\partial C}   & \equiv & 0\\
   *   \frac{\partial E}{\partial a_1} & \equiv & 0\\
   *   \frac{\partial E}{\partial N}   & \equiv & 0\\
   * @f{eqnarray*}
   * 
   * @param x       Where to evaluate the derivative
   * @param ipar    Parameter to differentiate relative to 
   * @param dPar    Variation in parameter to use in calculation
   * @param delta   @f$\Delta_{p}@f$ 
   * @param xi      @f$\xi@f$
   * @param sigma   @f$\sigma@f$ 
   * @param n       @f$N@f$ 
   * @param a       @f$\mathbf{a} = (a_2,...,a_N)@f$ 
   * 
   * @return @f$\partial E/\partial p|_{x}@f$  
   */
  Double_t DEstimatePar(Double_t x, Int_t ipar, Double_t delta, Double_t xi, 
			Double_t sigma, UShort_t n, Double_t* a, 
			Double_t dPar=0.001) 
  {
    Double_t dp    = dPar;
    Double_t d2    = dPar / 2;
    Double_t y1    = 0;
    Double_t y2    = 0;
    Double_t y3    = 0;
    Double_t y4    = 0;
    switch (ipar) { 
    case Function::kC: 
    case Function::kN: 
      return 0;
    case Function::kDelta:
      y1 = NEstimate(x, delta+dp, xi, sigma, n, a);
      y2 = NEstimate(x, delta+d2, xi, sigma, n, a);
      y3 = NEstimate(x, delta-d2, xi, sigma, n, a);
      y4 = NEstimate(x, delta-dp, xi, sigma, n, a);
      break;
    case Function::kXi:
      y1 = NEstimate(x, delta, xi+dp, sigma, n, a);
      y2 = NEstimate(x, delta, xi+d2, sigma, n, a);
      y3 = NEstimate(x, delta, xi-d2, sigma, n, a);
      y4 = NEstimate(x, delta, xi-dp, sigma, n, a);
      break;
    case Function::kSigma:
      y1 = NEstimate(x, delta, xi, sigma+dp, n, a);
      y2 = NEstimate(x, delta, xi, sigma+d2, n, a);
      y3 = NEstimate(x, delta, xi, sigma-d2, n, a);
      y4 = NEstimate(x, delta, xi, sigma-dp, n, a);
      break;
    default: {
      Int_t j = ipar-kA;
      if (j+1 > n) return 0;
      Double_t aa = a[j];
      a[j] = aa+dp; y1 = NEstimate(x, delta, xi, sigma, n, a);
      a[j] = aa+d2; y2 = NEstimate(x, delta, xi, sigma, n, a);
      a[j] = aa-d2; y3 = NEstimate(x, delta, xi, sigma, n, a);
      a[j] = aa-dp; y4 = NEstimate(x, delta, xi, sigma, n, a);
      a[j] = aa;
    }
      break;
    }
    Double_t d0  = y1 - y4;
    Double_t d1  = 2 * (y2 - y3);
    
    Double_t g   = 1/(2*dp) * (4*d1 - d0) / 3;    
    return g;
  }
  /** 
   * Generate a TF1 object of @f$ f_I@f$                                        
   *                                                                            
   * @param n        @f$ n@f$ - the number of particles                         
   * @param c        Constant                                                   
   * @param delta    @f$ \Delta@f$                                              
   * @param xi       @f$ \xi_1@f$                                               
   * @param sigma    @f$ \sigma_1@f$                                            
   * @param a        Array of @f$a_i@f$ for @f$i > 1@f$
   * @param xmin     Least value of range                                       
   * @param xmax     Largest value of range                                     
   *                                                                            
   * @return Newly allocated TF1 object                                         
   */  
  static TF1* MakeFunc(Int_t n, Double_t c, Double_t delta, 
		       Double_t xi, Double_t sigma, Double_t* a, 
		       Double_t xmin, Double_t xmax)
  {
    Int_t nPar = kN + n;
    TF1* f = new TF1(Form("landGaus%d",n), &landauGausN, xmin, xmax, nPar);
    f->SetNpx(500);
    f->SetParName(kC,      "C");
    f->SetParName(kDelta,  "#Delta_{p}");
    f->SetParName(kXi,     "#xi");
    f->SetParName(kSigma,  "#sigma");
    f->SetParName(kN,      "N");
    f->SetParameter(kC,     c);
    f->SetParameter(kDelta, delta);
    f->SetParameter(kXi,    xi);
    f->SetParameter(kSigma, sigma);
    f->FixParameter(kN,     n);
    f->SetParLimits(kDelta, xmin, 4);
    f->SetParLimits(kXi,    0.00, 4);
    f->SetParLimits(kSigma, 0.01, 0.11);

    for (Int_t i = 2; i <= n; i++) { 
      Int_t j = i-2;
      f->SetParName(kA+j, Form("a_{%d}", i));
      f->SetParameter(kA+j, a[j]);
      f->SetParLimits(kA+j, 0, 1);
    }
    return f;
  }
  /** 
   * Make a ROOT TF1 function object corresponding to a single 
   * component for @f$ i@f$ particles 
   * 
   * @param i        Number of particles
   * @param c        @f$C@f$ 
   * @param delta    @f$ \Delta@f$                                              
   * @param xi       @f$ \xi_1@f$                                               
   * @param sigma    @f$ \sigma_1@f$                                            
   * @param xmin     Minimum of range of @f$\Delta@f$
   * @param xmax     Maximum of range of @f$\Delta@f$
   * 
   * @return Pointer to newly allocated ROOT TF1 object 
   */
  static TF1* MakeIFunc(Int_t i, Double_t c, Double_t delta, 
			Double_t xi, Double_t sigma,
			Double_t xmin, Double_t xmax)
  {
    Int_t nPar = 5;
    TF1* f = new TF1(Form("landGausI%d",i), &landauGausI, xmin, xmax, nPar);
    f->SetNpx(100);
    f->SetParName(kC,      "C");
    f->SetParName(kDelta,  "#Delta_{p}");
    f->SetParName(kXi,     "#xi");
    f->SetParName(kSigma,  "#sigma");
    f->SetParName(kN,      "N");
    f->SetParameter(kC,     c);
    f->SetParameter(kDelta, delta);
    f->SetParameter(kXi,    xi);
    f->SetParameter(kSigma, sigma);
    f->FixParameter(kN,     i);
    return f;
  }

  // --- Object code -------------------------------------------------
  /** 
   * Constructor 
   * 
   * @param n 
   * @param c 
   * @param delta 
   * @param xi 
   * @param sigma 
   * @param a 
   * @param xmin 
   * @param xmax 
   * 
   * @return 
   */
  Function(Int_t n, Double_t c, Double_t delta, Double_t xi, Double_t sigma, 
	   Double_t* a, Double_t xmin, Double_t xmax)
    : fF(MakeFunc(n, c, delta, xi, sigma, a, xmin, xmax))
  {}
  Function(TF1* f) : fF(f) {}
  /** 
   * Get pointer to ROOT TF1 object 
   * 
   * @return Pointer to TF1 object
   */
  TF1* GetF1() { return fF; }
  /** 
   * Evaluate the function at @a x
   * @f[ 
   *  f_N(x;\Delta,\xi,\sigma) = 
   *     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i)
   * @f] 
   * 
   * @param x Where to evaluate 
   * 
   * @return 
   */
  Double_t Evaluate(Double_t x) const
  {
    return fF->Eval(x);
  }
  /** 
   * Evaluate the function at @a x
   * @f[ 
   *  f_N(x;\Delta,\xi,\sigma) = 
   *     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i)
   * @f] 
   * 
   * This also calculates the error on the point like 
   * @f[
   *   \delta^2 f_N = 
   *    \left(\frac{\partial f_N}{\partial\Delta_p}\right)^2\delta^2\Delta_p+ 
   *    \left(\frac{\partial f_N}{\partial\xi     }\right)^2\delta^2\xi+ 
   *    \left(\frac{\partial f_N}{\partial\sigma  }\right)^2\delta^2\sigma+ 
   *    \sum_{i=2}^N\left(\frac{\partial f_N}{\partial a_i}\right)^2\delta^2a_i
   * @f]
   *
   * @param x Where to evaluate 
   * @param e On return, contains the error @f$\delta f_N@f$ on the
   * function value
   * 
   * @return @f$ f_N(x;\Delta,\xi,\sigma)@f$ 
   */
  Double_t Evaluate(Double_t x, Double_t& e) const
  {
    Double_t delta     = GetDelta();
    Double_t xi        = GetXi();
    Double_t sigma     = GetSigma();
    Double_t dFdDelta2 = 0;
    Double_t dFdXi2    = 0;
    Double_t dFdSigma2 = 0;
    e              = 0;
    for (Int_t i = 1; i <= GetN(); i++) { 
      Double_t a      = GetA(i);
      Double_t dDelta = a*IdLandauGausdPar(i, x, 1, 0.001, delta, xi, sigma); 
      Double_t dXi    = a*IdLandauGausdPar(i, x, 2, 0.001, delta, xi, sigma);
      Double_t dSigma = a*IdLandauGausdPar(i, x, 3, 0.001, delta, xi, sigma);
      Double_t dFda   = a*ILandauGaus(i, x, delta, xi, sigma);
      dFdDelta2 += dDelta * dDelta;
      dFdXi2    += dXi    * dXi;
      dFdSigma2 += dSigma * dSigma;
      e += TMath::Power(dFda * GetEA(i),2);
    }      
    Double_t edelta = GetEDelta();
    Double_t exi    = GetEXi();
    Double_t esigma = GetESigma();
    e += (dFdDelta2 * edelta * edelta + 
	  dFdXi2    * exi    * exi    + 
	  dFdSigma2 * esigma * esigma);
    e = TMath::Sqrt(e);
    return fF->Eval(x);
  }

  /** 
   * Evaluate 
   * @f[ 
   *   E(x;\Delta,\xi,\sigma,\mathbf{a}) = 
   *   \frac{\sum_{i=1}^{n} i a_i f_i(x;\Delta,\xi,\sigma)}{
   *     f_N(x;\Delta,\xi,\sigma)} = 
   *   \frac{\sum_{i=1}^{n} i a_i f(x;\Delta_i,\xi_i,\sigma_i)}{
   *     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i)}
   * @f] 
   *
   * If the denominator is less than @f$10^{-4}@f$, then 0 is returned.
   * 
   * @param x    Where to evaluate @f$ E@f$ 
   * @param maxN Maximum number of weights 
   * 
   * @return @f$ E(x;\Delta,\xi,\sigma,\mathbf{a})@f$
   */  
  Double_t EstimateNParticles(Double_t x, Short_t maxN=-1)
  {
    UShort_t n     = (maxN < 0 ? GetN() : TMath::Min(UShort_t(maxN), GetN()));
    Double_t delta = GetDelta();
    Double_t xi    = GetXi();
    Double_t sigma = GetSigma();
    return NEstimate(x, delta, xi, sigma, n, GetAs());
  }
  /** 
   * Evaluate 
   * @f[ 
   *   E(x;\Delta,\xi,\sigma,\mathbf{a}) = 
   *   \frac{\sum_{i=1}^{n} i a_i f_i(x;\Delta,\xi,\sigma)}{
   *     f_N(x;\Delta,\xi,\sigma,\mathbf{a})} = 
   *   \frac{\sum_{i=1}^{n} i a_i f(x;\Delta_i,\xi_i,\sigma_i)}{
   *     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i)}
   * @f] 
   *
   * If @f$f_N(x;\Delta,\xi,\sigma,\mathbf{a})<10^{-4}@f$ then 0 is
   * returned.
   *
   * This also calculatues the error @f$\delta E@f$ of the
   * function value at @f$x@f$: 
   *
   * @f[
   *   \delta^2 E = 
   *    \left(\frac{\partial E}{\partial\Delta_p}\right)^2\delta^2\Delta_p+ 
   *    \left(\frac{\partial E}{\partial\xi     }\right)^2\delta^2\xi+ 
   *    \left(\frac{\partial E}{\partial\sigma  }\right)^2\delta^2\sigma+ 
   *    \sum_{i=2}^N\left(\frac{\partial E}{\partial a_i}\right)^2\delta^2a_i
   * @f]
   * The partial derivatives are evaluated numerically. 
   * 
   * @param x    Where to evaluate @f$ E@f$ 
   * @param e    On return, @f$\delta E|_{x}@f$ 
   * @param maxN Maximum number of weights 
   * 
   * @return @f$ E(x;\Delta,\xi,\sigma,\mathbf{a})@f$
   */  
  Double_t EstimateNParticles(Double_t x, Double_t& e, Short_t maxN=-1)
  {
    
    UShort_t  n      = (maxN < 0 ? GetN() : TMath::Min(UShort_t(maxN), GetN()));
    Double_t  delta  = GetDelta();
    Double_t  xi     = GetXi();
    Double_t  sigma  = GetSigma();
    Double_t* a      = GetAs();
    Double_t  dDelta = (DEstimatePar(x,Function::kDelta,delta,xi,sigma,n,a,\
				     0.01 * GetEDelta()) * GetEDelta());
    Double_t  dXi    = (DEstimatePar(x,Function::kXi,   delta,xi,sigma,n,a,
				     0.01 * GetEXi())    * GetEXi());
    Double_t  dSigma = (DEstimatePar(x,Function::kSigma,delta,xi,sigma,n,a,
				     0.01 * GetESigma()) * GetESigma());
    e                = dDelta * dDelta + dXi * dXi + dSigma * dSigma;
    for (Int_t i = 2; i <= n; i++) { 
      Double_t dAi = (DEstimatePar(x, Function::kA+i-2,delta,xi,sigma,n,a,
				   0.01 * GetEA(i)) * GetEA(i));
      e += dAi * dAi;
    }
    e = TMath::Sqrt(e);
    return NEstimate(x, GetDelta(), GetXi(), GetSigma(), n, GetAs());
  }
  /** 
   * Estimate the number of particles by calculating the weights 
   * @f[ 
   *   w_i = a_i f_i(x;\Delta_p,\xi,\sigma) 
   * @f] 
   * 
   * and then draw a random number @f$r@f$ in the range @f$[0,\sum_i^N
   * w_i]@f$ and return the @f$i@f$ for which @f$w_{i-1} < @r \leq
   * w_{i}@f$ (notem @f$w_{-1}=0@f$).
   * 
   * @param x     Where to evaluate the component functions 
   * @param maxN  Maximum weight to use 
   * 
   * @return Estimate of the number of particles using the procedure
   * outlined above. 
   */
  Double_t RandomEstimateNParticles(Double_t x, Short_t maxN=-1)
  {
    UShort_t n     = (maxN < 0 ? GetN() : TMath::Min(UShort_t(maxN), GetN()));
    TArrayD  p(n);
    Double_t delta = GetDelta();
    Double_t xi    = GetXi();
    Double_t sigma = GetSigma();
    
    if (Evaluate(x) < 1e-4) return 0;

    // std::cout << "RandomEstimateNParticles(" << x << "):";
    for (Int_t i = 1; i <= n; i++) { 
      Double_t a = GetA(i);
      Double_t f = ILandauGaus(i, x, delta, xi, sigma);
      p[i-1]     = a * f;
      // std::cout << p[i-1] << ",";
    }
    Double_t r = gRandom->Uniform(p.GetSum());
    Double_t l = 0;
    // std::cout << "s=" << p.GetSum() << ",r=" << r << ",";
    for (Int_t i = 1; i <= n; i++) { 
      if (r > l && r <= l + p[i-1]) {
	// std::cout << "l=" << l << ",l+p=" << l+p[i-1] << "-> " << i 
	//  << std::endl;
	return i;
      }
      l += p[i-1];
    }
    return 0;
  }
  /** 
   * Draw the function 
   * 
   * @param option Drawing option 
   */
  void Draw(Option_t* option="") { fF->Draw(option); }
  /** @return @f$C@f$ */
  Double_t GetC() const { return fF->GetParameter(kC); }
  /** @return @f$\Delta_{mp}@f$ */
  Double_t GetDelta() const { return fF->GetParameter(kDelta); }
  /** @return @f$\xi@f$ */
  Double_t GetXi() const { return fF->GetParameter(kXi); }
  /** @return @f$\sigma@f$ */
  Double_t GetSigma() const { return fF->GetParameter(kSigma); }
  /** @return @f$N@f$ */
  UShort_t GetN() const { return fF->GetParameter(kN); }
  /** @return @f$a_{i}@f$ */
  Double_t GetA(Int_t i) const { return i == 1 ? 1 : fF->GetParameter(kA+i-2);}
  /** @return @f$a_2,...,a_N@f$ */
  Double_t* GetAs() const { return &(fF->GetParameters()[kA]);}
  /** @return @f$\Delta C@f$ */
  Double_t GetEC() const { return fF->GetParError(kC); }
  /** @return @f$\Delta \Delta_{mp}@f$ */
  Double_t GetEDelta() const { return fF->GetParError(kDelta); }
  /** @return @f$\Delta \xi@f$ */
  Double_t GetEXi() const { return fF->GetParError(kXi); }
  /** @return @f$\Delta \sigma@f$ */
  Double_t GetESigma() const { return fF->GetParError(kSigma); }
  /** @return @f$\Delta N@f$ */
  UShort_t GetEN() const { return 0; }
  /** @return @f$\Delta a_{i}@f$ */
  Double_t GetEA(Int_t i) const { return i == 1 ? 0 : fF->GetParError(kA+i-2);}
  /** @return @f$\Delta a_2,...,\Delta a_N@f$ */
  Double_t* GetEAs() const { return &(fF->GetParErrors()[kA]);}
  /** The ROOT object */
  TF1* fF;
};
#ifndef __CINT__
const Double_t Function::fgkMPShift    = -0.22278298;
const Double_t Function::fgkInvRoot2Pi = 1. / TMath::Sqrt(2*TMath::Pi());
const Double_t Function::fgkConvNSigma = 5;
const Double_t Function::fgkConvNSteps = 100;
#endif
  
//====================================================================
/**
 * 
 * 
 * 
 * @ingroup pwg2_forward_scripts_tests
 */
struct Fitter
{
  // --- Object code ------------------------------------------------
  Fitter(Int_t nMinus) 
    : fNMinus(nMinus)
  {
  }
  /** 
   * 
   * 
   * @param dist 
   * 
   * @return 
   */
  TF1* FitFirst(TH1* dist)
  {
    fFunctions.Clear();
    
    // Get upper edge 
    Int_t    midBin = dist->GetMaximumBin();
    Double_t midE   = dist->GetBinLowEdge(midBin);
    
    // Get low edge 
    Int_t    minBin = midBin - fNMinus;
    Double_t minE  = dist->GetBinCenter(minBin);
    Double_t maxE  = dist->GetBinCenter(midBin+2*fNMinus);
    Double_t a[]   = { 0 };
    TF1* f = Function::MakeFunc(1, 0.5, midE, 0.07, 0.1, a, minE, maxE);

    // Do the fit 
    dist->Fit(f, "RNQS", "", minE, maxE);
    fFunctions.AddAtAndExpand(f, 0);
    return f;
  }
  /** 
   * 
   * 
   * @param dist 
   * @param n 
   * 
   * @return 
   */
  Function* FitN(TH1* dist, Int_t n)
  {
    if (n == 1) return new Function(FitFirst(dist));
    TF1*        f1 = static_cast<TF1*>(fFunctions.At(0));
    if (!f1) { 
      f1 = FitFirst(dist);
      if (!f1) { 
	::Warning("FitN", "First fit missing");
	return 0;
      }
    }
    
    Double_t delta1 = f1->GetParameter(Function::kDelta);
    Double_t xi1    = f1->GetParameter(Function::kXi);
    Double_t maxEi  = n * (delta1 + xi1 * TMath::Log(n)) + 2 * n * xi1;
    Double_t minE   = f1->GetXmin();
    
    TArrayD a(n-1);
    for (UShort_t i = 2; i <= n; i++) a.fArray[i-2] = (n==2? 0.05 : 0.000001);
    
    TF1* fn = Function::MakeFunc(n, f1->GetParameter(Function::kC), 
				 delta1, xi1, 
				 f1->GetParameter(Function::kSigma), 
				 a.fArray, minE, maxEi);
    
    dist->Fit(fn, "RSNQ", "", minE, maxEi);
    dist->Fit(fn, "RSNQM", "", minE, maxEi);
    fFunctions.AddAtAndExpand(fn, n-1);
    
    return new Function(fn);
  }
  Int_t     fNMinus;
  TObjArray fFunctions;
};

//====================================================================
/** 
 * Utility function 
 * 
 * @param xp Pointer to independent variables
 * @param pp Pointer to parameters 
 * 
 * @return Function evaluated at xp[0]
 * 
 * @ingroup pwg2_forward_scripts_tests
 */  
static Double_t landauGaus1(Double_t* xp, Double_t* pp)
{
  Double_t x     = xp[0];
  Double_t c     = pp[Function::kC];
  Double_t delta = pp[Function::kDelta];
  Double_t xi    = pp[Function::kXi];
  Double_t sigma = pp[Function::kSigma];
    
  return c * Function::LandauGaus(x, delta, xi, sigma);
}
/** 
 * Utility function 
 * 
 * @param xp Pointer to independent variables
 * @param pp Pointer to parameters 
 * 
 * @return Function evaluated at xp[0]
 * 
 * @ingroup pwg2_forward_scripts_tests
 */
static Double_t landauGausN(Double_t* xp, Double_t* pp)
{
  Double_t  x     = xp[0];
  Double_t  c     = pp[Function::kC];
  Double_t  delta = pp[Function::kDelta];
  Double_t  xi    = pp[Function::kXi];
  Double_t  sigma = pp[Function::kSigma];
  Int_t     n     = Int_t(pp[Function::kN]);
  Double_t* a     = &pp[Function::kA];
    
  return c * Function::NLandauGaus(x, delta, xi, sigma, n, a);
}
/** 
 * Utility function 
 * 
 * @param xp Pointer to independent variables
 * @param pp Pointer to parameters 
 * 
 * @return Function evaluated at xp[0]
 * 
 * @ingroup pwg2_forward_scripts_tests
 */
static Double_t landauGausI(Double_t* xp, Double_t* pp)
{
  Double_t  x     = xp[0];
  Double_t  c     = pp[Function::kC];
  Double_t  delta = pp[Function::kDelta];
  Double_t  xi    = pp[Function::kXi];
  Double_t  sigma = pp[Function::kSigma];
  Int_t     i     = Int_t(pp[Function::kN]);
    
  return c * Function::ILandauGaus(i, x, delta, xi, sigma);
}

//====================================================================
/**
 * 
 * 
 * 
 * @ingroup pwg2_forward_scripts_tests
 */
struct Generator
{
  Double_t fDelta;
  Double_t fXi;
  Double_t fSigma;
  Int_t    fMaxN;
  TH1*     fSingle;
  TH1*     fSummed;
  TH1*     fSimple;
  TH2*     fNSum;
  TH1*     fN;
  TArrayD  fRP; // Raw probabilities
  TArrayD  fP; // Normalized probabilities 
  TF1*     fF;
  Int_t    fNObs;

  /** 
   * 
   * 
   * 
   * @return 
   */
  Generator(Double_t delta=.55, Double_t xi=0.06, Double_t sigma=0.02,
	    Int_t maxN=5, const Double_t* a=0)
    : fDelta(delta), fXi(xi), fSigma(sigma), fMaxN(maxN),
      fSingle(0), fSummed(0), fSimple(0), fNSum(0), fRP(), fP(fMaxN), fF(0),
      fNObs(0)
  {
    fMaxN = TMath::Min(fMaxN, 5);
    const Double_t ps[] = { 1, .05, .02, .005, .002 };
    fRP.Set(5, ps);
    if (a) { 
      for (Int_t i = 0; i < maxN; i++) fRP[i] = a[i];
    }

    Double_t sum = 0; 
    for (Int_t i = 0; i < fMaxN; i++) sum += fRP[i];
    for (Int_t i = 0; i < fP.fN; i++) {
      fP[i] = fRP[i] / sum;
      if (i != 0) fP[i] += fP[i-1];
    }

    fSingle = new TH1D("single", "Single signals", 500, 0, 10);
    fSingle->SetXTitle("Signal");
    fSingle->SetYTitle("Events");
    fSingle->SetFillColor(kRed+1);
    fSingle->SetMarkerColor(kRed+1);
    fSingle->SetLineColor(kRed+1);
    fSingle->SetFillStyle(3001);
    fSingle->SetMarkerStyle(21);
    fSingle->SetDirectory(0);
    fSingle->Sumw2();

    fSummed = static_cast<TH1D*>(fSingle->Clone("summed"));
    fSummed->SetTitle("Summed signals");
    fSummed->SetFillColor(kBlue+1);
    fSummed->SetMarkerColor(kBlue+1);
    fSummed->SetLineColor(kBlue+1);
    fSummed->SetDirectory(0);

    fSimple = static_cast<TH1D*>(fSingle->Clone("summed"));
    fSimple->SetTitle("Summed signals");
    fSimple->SetFillColor(kGreen+1);
    fSimple->SetMarkerColor(kGreen+1);
    fSimple->SetLineColor(kGreen+1);
    fSimple->SetDirectory(0);

    fNSum = new TH2D("single", "Observation signals", 
		     fMaxN, .5, fMaxN+.5, 100, 0, 10);
    fNSum->SetMarkerColor(kMagenta+1);
    fNSum->SetLineColor(kMagenta+1);
    fNSum->SetFillColor(kMagenta+1);
    fNSum->SetYTitle("Signal");
    fNSum->SetXTitle("N");
    fNSum->SetZTitle("Events");
    fNSum->SetDirectory(0);
    fNSum->Sumw2();

    fN = new TH1D("n", "Number of particles", fMaxN, .5, fMaxN+.5);
    fN->SetXTitle("N");
    fN->SetYTitle("Events");
    fN->SetMarkerColor(kMagenta+1);
    fN->SetLineColor(kMagenta+1);
    fN->SetFillColor(kMagenta+1);
    fN->SetMarkerStyle(22);
    fN->SetFillStyle(3001);
    fN->SetDirectory(0);

    std::cout << "Probabilities:" << std::setprecision(4) 
	      << "\n  ";
    for (Int_t i = 0; i < fRP.fN; i++) 
      std::cout << std::setw(7) << fRP[i] << " ";
    std::cout << " = " << fRP.GetSum() << "\n  ";
    for (Int_t i = 0; i < fP.fN; i++) 
      std::cout << std::setw(7) << fP[i] << " ";
    std::cout << std::endl;

    Double_t aa[] = { 0 };
    fF = Function::MakeFunc(1, 1, fDelta, fXi, fSigma, aa,  0, 10);
  }
  /** 
   * 
   * 
   */
  void Clear()
  {
    fSingle->Reset();
    fSummed->Reset();
    fNSum->Reset();
    fN->Reset();
  }
  /** 
   * 
   * 
   * 
   * @return 
   */    
  Int_t GenerateN()
  {
    Double_t rndm = gRandom->Uniform();
    Int_t    ret  = 0;
    if      (                 rndm < fP[0]) ret = 1;
    else if (rndm >= fP[0] && rndm < fP[1]) ret = 2;
    else if (rndm >= fP[1] && rndm < fP[2]) ret = 3;
    else if (rndm >= fP[2] && rndm < fP[3]) ret = 4;
    else if (rndm >= fP[3] && rndm < fP[4]) ret = 5;
    // Printf("%f -> %d", rndm, ret);
    fN->Fill(ret);
    return ret;
  }
  /** 
   * 
   * 
   * @param mpv 
   * @param xi 
   * @param sigma 
   * 
   * @return 
   */
  Double_t Generate1Signal()
  {
    Double_t ret = 0;
    if (fF) ret = fF->GetRandom();
    else { 
      Double_t rmpv = gRandom->Gaus(fDelta, fSigma);
      ret           = gRandom->Landau(rmpv, fXi);
    }
    fSingle->Fill(ret);
    fSimple->Fill(gRandom->Landau(fDelta - fXi * Function::fgkMPShift, fXi));
    return ret;
  }
  /** 
   * 
   * 
   * @param n 
   * @param mpv 
   * @param xi 
   * @param sigma 
   * 
   * @return 
   */  
  Double_t GenerateNSignal(Int_t n)
  {
    Double_t ret = 0;
    for (Int_t i = 0; i < n; i++) ret += Generate1Signal();
    fSummed->Fill(ret);
    fNSum->Fill(n, ret);
    return ret;
  }
  /** 
   * 
   * 
   * @param nObs 
   * @param reset 
   */  
  void Generate(Int_t nObs=1000, Bool_t reset=true)
  {
    if (reset) Clear();
    fNObs = nObs;

    for (Int_t i = 0; i < nObs; i++) {
      // if (((i+1) % (nObs/10)) == 0) 
      //	std::cout << "Event " << std::setw(6) << i << std::endl;
      Int_t n = GenerateN();
      GenerateNSignal(n);
    }
    Double_t m = fSingle->GetMaximum();
    fSingle->Scale(1. / m);
    m = fSummed->GetMaximum();
    fSummed->Scale(1. / m);
    m = fSimple->GetMaximum();
    fSimple->Scale(1. / m);

    std::cout << "Resulting probabilities:\n  ";
    for (Int_t i = 1; i <= fN->GetNbinsX(); i++) 
      std::cout << std::setprecision(4) << std::setw(7) 
		<< fN->GetBinContent(i) << " ";
    std::cout << std::endl;
  }
  /** 
   * Pretty print of parameters. 
   * 
   * @param input 
   * @param f 
   * @param i 
   */
  void PrintParameter(Double_t input, TF1* f, Int_t i)
  {
    Double_t val = 0, err = 0;
    if (i < f->GetNpar()) {
      val = f->GetParameter(i);
      err = f->GetParError(i);
    }
    std::cout << "  " << std::setw(16) << f->GetParName(i) << ": "
	      << std::fixed 
	      << std::setprecision(4) << std::setw(6) << input << " -> " 
	      << std::setprecision(4) << std::setw(6) << val   << " +/- "
	      << std::setprecision(5) << std::setw(7) << err   << " [delta:" 
	      << std::setw(3) << int(100 * (input-val) / input) << "%,err:" 
	      << std::setw(3) << int(100 * err / val) << "%]"
	      << std::endl;
  }
  /** 
   * Print the result 
   * 
   * @param f 
   */
  void PrintResult(TF1* f) 
  {
    Double_t chi2 = f->GetChisquare();
    Int_t    ndf  = f->GetNDF();
    std::cout << "Chi^2/NDF = " << chi2 << "/" << ndf << "=" 
	      << chi2/ndf << std::endl;
    PrintParameter(1,      f, Function::kC);
    PrintParameter(fDelta, f, Function::kDelta);
    PrintParameter(fXi,    f, Function::kXi);
    PrintParameter(fSigma, f, Function::kSigma);
    PrintParameter(fMaxN,  f, Function::kN);
    for (Int_t i = 0; i < fMaxN-1; i++) 
      PrintParameter(fRP[i+1], f, Function::kA+i);
  }
  /** 
   * 
   * 
   * @param x 
   * @param y 
   * @param intput 
   * @param f 
   * @param i 
   */
  void DrawParameter(Double_t x, Double_t y, Double_t input, TF1* f, Int_t i)
  {
    Double_t val = 0, err = 0;
    if (i < f->GetNpar()) {
      val = f->GetParameter(i);
      err = f->GetParError(i);
    }

    TLatex* ltx = new TLatex(x, y, f->GetParName(i));
    ltx->SetTextSize(.04);
    ltx->SetTextFont(132);
    ltx->SetTextAlign(11);
    ltx->SetNDC();
    ltx->Draw();
    
    ltx->DrawLatex(x+.05,  y, Form("%5.3f", input));
    ltx->DrawLatex(x+.14, y, "#rightarrow");
    ltx->SetTextAlign(21);
    ltx->DrawLatex(x+.3,  y, Form("%5.3f #pm %6.4f", val, err));
    ltx->SetTextAlign(11);
    ltx->DrawLatex(x+.48,  y, Form("[#Delta=%3d%%, #delta=%3d%%]",
				  int(100 * (input-val)/input), 
				  int(100*err/val)));
  }
  /** 
   * 
   * 
   * @param x 
   * @param y 
   * @param f 
   */
  void DrawResult(Double_t x, Double_t y, TF1* f)
  {
    Double_t chi2 = f->GetChisquare();
    Int_t    ndf  = f->GetNDF();
    TLatex* ltx = new TLatex(x, y, Form("#chi^{2}/NDF=%5.2f/%3d=%6.3f", 
					chi2, ndf, chi2/ndf));
    ltx->SetTextSize(.04);
    ltx->SetNDC();
    ltx->SetTextFont(132);
    ltx->Draw();

    x += .05;
    y -= .035; DrawParameter(x, y, 1,      f, Function::kC);
    y -= .035; DrawParameter(x, y, fDelta, f, Function::kDelta);
    y -= .035; DrawParameter(x, y, fXi,    f, Function::kXi);
    y -= .035; DrawParameter(x, y, fSigma, f, Function::kSigma);
    y -= .035; DrawParameter(x, y, fMaxN,  f, Function::kN);
    for (Int_t i = 0; i < fMaxN-1; i++) {
      y -= .035; DrawParameter(x, y, fRP[i+1], f, Function::kA+i);
    }
  }
  /** 
   * 
   * 
   */
  void Draw(const char* type="png")
  {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    TCanvas* c = new TCanvas("c", "Generator distributions", 1200, 1000);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.01);

    c->Divide(2, 2);

    // --- First pad:  Data (multi, single, landau) + Fit -------------
    TVirtualPad* p = c->cd(1);
    p->SetLogy();
    p->SetBorderSize(0);
    p->SetBorderMode(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    fSummed->Draw();
    fSingle->Draw("same");
    fSimple->Draw("same");

    if (fF) { 
      Double_t max = fF->GetMaximum();
      Double_t a[] = { 0 };
      TF1* fcopy = Function::MakeFunc(1, 
				      fF->GetParameter(Function::kC) /max, 
				      fF->GetParameter(Function::kDelta), 
				      fF->GetParameter(Function::kXi), 
				      fF->GetParameter(Function::kSigma), 
				      a, 0, 10);
      fcopy->SetLineColor(2);
      fcopy->SetLineStyle(2);
      fcopy->Draw("same");
    }

    Fitter*   fitter = new Fitter(4);
    Function* fit    = fitter->FitN(fSummed, 5);
    TF1*      f      = fit->GetF1();
    f->Draw("same");
    f->SetLineColor(kBlack);
    TF1* fst = static_cast<TF1*>(fitter->fFunctions.At(0));
    fst->SetLineWidth(3);
    fst->SetLineStyle(2);
    fst->SetLineColor(kMagenta);
    fst->Draw("same");
    PrintResult(f);
    DrawResult(.25, .95, f);
    TGraphErrors* gr = new TGraphErrors(fSummed->GetNbinsX());
    gr->SetName("fitError");
    gr->SetTitle("Fit w/errors");
    gr->SetFillColor(kYellow+1);
    gr->SetFillStyle(3001);
    for (Int_t j = 1; j <= fSummed->GetNbinsX(); j++) { 
      Double_t x = fSummed->GetBinCenter(j);
      Double_t e = 0;
      Double_t y = fit->Evaluate(x, e);
      gr->SetPoint(j-1, x, y);
      gr->SetPointError(j-1, 0, e);
    }
    gr->Draw("same l e4");

    TLegend* l = new TLegend(.15, .11, .5, .4);
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetTextFont(132);
    l->AddEntry(fSimple, "L: Single-L", "lp");
    l->AddEntry(fSingle, "LG: Single-L#otimesG", "lp");
    l->AddEntry(fSummed,"NLG: Multi-L#otimesG", "lp");
    l->AddEntry(gr, "f_{N}: Fit to NLG");
    l->Draw();

    // --- Second pad:  Data (1, 2, 3, ...) --------------------------
    p = c->cd(2);
    p->SetLogy();
    p->SetBorderSize(0);
    p->SetBorderMode(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetGridx();
    THStack* stack = new THStack(fNSum, "y");
    TIter next(stack->GetHists());
    TH1* h = 0;
    Int_t i = 2;
    Double_t max = 0;
    TLegend* l2 = new TLegend(.6, .5, .95, .95);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->SetFillColor(0);
    l2->SetTextFont(132);
    while ((h = static_cast<TH1*>(next()))) { 
      if (i == 2) max = h->GetMaximum();
      h->SetLineColor(i); 
      h->SetFillColor(i); 
      h->SetFillStyle(3001);
      h->Scale(1/max);
      l2->AddEntry(h, Form("%d particles", i - 1), "f");
      i++;
    }
    stack->Draw("hist nostack");
    stack->GetHistogram()->SetXTitle("Signal");
    stack->GetHistogram()->SetYTitle("Events");
    i = 2;
    for (Int_t j = 1; j <= 5; j++) {
      TF1* fi = Function::MakeIFunc(j, fRP[j-1], fDelta, fXi, fSigma, 0, 10);
      if (j == 1) max = fi->GetMaximum();
      fi->SetParameter(0, fRP[j-1]/max);
      fi->SetLineColor(i);
      fi->SetLineWidth(3);
      fi->SetLineStyle(2);
      fi->Draw("same");

      TF1* fj = Function::MakeIFunc(j, fit->GetC() * fit->GetA(j), 
				    fit->GetDelta(), fit->GetXi(), 
				    fit->GetSigma(), 0, 10);
      fj->SetLineColor(i);
      fj->SetLineWidth(2);
      fj->Draw("same");

      std::cout << "Component " << j << " scale=" << fi->GetParameter(0) 
		<< " (" << fRP[j-1] << ")" << " " << fj->GetParameter(0) 
		<< " (" << fit->GetC() << "*" << fit->GetA(j) << ")" 
		<< std::endl;
      i++;
    }
    TLegendEntry* ee = l2->AddEntry("", "Input f_{i}", "l");
    ee->SetLineStyle(2);
    ee->SetLineWidth(3);
    l2->AddEntry("", "Fit f_{i}", "l");
    l2->Draw();
    // fNSum->Draw("lego2");

    // --- Third pad:  Mean multiplicity ------------------------------
    TH1* nhist1 = new TH1F("nhist1", "NHist", 5, 0.5, 5.5);
    nhist1->SetFillColor(kCyan+1);
    nhist1->SetMarkerColor(kCyan+1);
    nhist1->SetLineColor(kCyan+1);
    nhist1->SetFillStyle(3002);
    TH1* nhist2 = new TH1F("nhist2", "NHist", 5, 0.5, 5.5);
    nhist2->SetFillColor(kYellow+1);
    nhist2->SetMarkerColor(kYellow+1);
    nhist2->SetLineColor(kYellow+1);
    nhist2->SetFillStyle(3002);

    TH2* resp = new TH2F("resp", "Reponse", 100, 0, 10, 6, -.5, 5.5);
    resp->SetXTitle("Signal");
    resp->SetYTitle("# particles");
    for (Int_t j = 0; j < fNObs; j++) { 
      if (((j+1) % (fNObs / 10)) == 0) 
	std::cout << "Event # " << j+1 << std::endl;
      Double_t x = fSummed->GetRandom();
      Double_t y = fit->RandomEstimateNParticles(x);
      nhist1->Fill(y);
      resp->Fill(x, y);
    }
    TGraphErrors* graph = new TGraphErrors(fSummed->GetNbinsX());
    graph->SetName("evalWeighted");
    graph->SetTitle("Evaluated function");
    graph->SetLineColor(kBlack);
    graph->SetFillColor(kYellow+1);
    graph->SetFillStyle(3001);
    for (Int_t j = 1; j <= fSummed->GetNbinsX(); j++) { 
      Double_t x = fSummed->GetBinCenter(j);
      Double_t e = 0;
      Double_t y = fit->EstimateNParticles(x, e);
      nhist2->Fill(y, fSummed->GetBinContent(j));
      graph->SetPoint(j-1, x, y);
      graph->SetPointError(j-1, 0, e);
    }
    nhist1->Scale(1. / nhist1->GetMaximum());
    nhist2->Scale(1. / nhist2->GetMaximum());
    fN->Scale(1. / fN->GetMaximum());

    p = c->cd(3);
    p->SetLogy();
    p->SetBorderSize(0);
    p->SetBorderMode(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    fN->Draw();
    nhist1->Draw("same");
    nhist2->Draw("same");

    TLegend* l3 = new TLegend(.3, .7, .95, .95);
    l3->SetBorderSize(0);
    l3->SetFillStyle(0);
    l3->SetFillColor(0);
    l3->SetTextFont(132);
    l3->AddEntry(fN,  
		 Form("Input n distribution: #bar{m}=%5.3f, s_{m}=%5.3f",
		      fN->GetMean(), fN->GetRMS()), "f");
    l3->AddEntry(nhist1, 
		 Form("Random draw from fit: #bar{m}=%5.3f, s_{m}=%5.3f", 
		      nhist1->GetMean(), nhist1->GetRMS()), "f");
    l3->AddEntry(nhist2, 
		 Form("Weighted evaluation of fit: #bar{m}=%5.3f, s_{m}=%5.3f", 
		      nhist2->GetMean(), nhist2->GetRMS()), "f");
    l3->Draw();
    
    // --- Fourth pad: Reponse function ------------------------------
    p = c->cd(4);
    // p->SetLogy();
    p->SetBorderSize(0);
    p->SetBorderMode(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetGridx();
    TProfile* prof1 = fNSum->ProfileY();
    prof1->SetLineColor(kMagenta+1);
    prof1->Draw();
    TProfile* prof2 = resp->ProfileX();
    prof2->SetLineColor(kCyan+2);
    prof2->Draw("same");
    graph->SetLineWidth(2);
    graph->Draw("same LP E4");

    TLegend* l4 = new TLegend(.2, .11, .8, .4);
    l4->SetBorderSize(0);
    l4->SetFillStyle(0);
    l4->SetFillColor(0);
    l4->SetTextFont(132);
    l4->AddEntry(prof1, "Input distribution of N particles", "l");
    l4->AddEntry(prof2, "Random estimate of N particles", "l");
    l4->AddEntry(graph, "E: #sum_{i}^{N}i a_{i}f_{i}/"
		 "#sum_{i}^{N}a_{i}f_{i} evaluated over NLG", "lf");
    l4->Draw();

    c->cd();
    
    if (type && type[0] != '\0') 
      c->Print(Form("TestELossDist.%s", type));
  }
};

/** 
 * Test the energy loss fits
 * 
 * @param type Output graphics type 
 * 
 * @ingroup pwg2_forward_scripts_tests
 */
void
TestELossDist(const char* type="png")
{
  gRandom->SetSeed(12345);
  const Double_t a1[] = { 1, .05, .02, .005, .002 }; // 7TeV pp
  const Double_t a2[] = { 1, .1, .05, .01, .005 };
  const Double_t a3[] = { 1, .2, .10, .05, .01 };

  Generator* g = new Generator(0.55, 0.06, 0.06, 5, a1);
  g->Generate(100000);
  g->Draw(type);

  
}


