#ifndef ALILANDAUGAUS_H
#define ALILANDAUGAUS_H
/**
 * @file   AliLandauGaus.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Mar 11 08:53:47 2014
 * 
 * @brief  Declaration and implementation of Landau-Gauss distributions. 
 * 
 * @ingroup pwglf_forward 
 */

#include <TObject.h>
#include <TF1.h>
#include <TMath.h>

/** 
 * This class contains static member functions to calculate the energy
 * loss stragling - most notably the N-particle energy loss as a sum
 * of convolutions of a Landau and Gauss distribution.
 *
 * That is, for a single particle we have the function @f$ f(x)@f$: 
 *
 * @f[ 
 * f(x;\Delta_p,\xi,\sigma') = \frac{1}{\sigma' \sqrt{2 \pi}}
 *    \int_{-\infty}^{+\infty} dx' f'_{L}(x',\Delta_p,\xi)
 *    \exp{-\frac{(x-x')^2}{2\sigma'^2}}
 * @f]
 * 
 * where @f$ f'_{L}@f$ is the Landau distribution, @f$\Delta_p@f$ the
 * most probable energy loss, @f$ \xi@f$ the width of the Landau, and
 * @f$ \sigma'^2=\sigma^2-\sigma_n^2 @f$.  Here, @f$\sigma@f$ is the
 * variance of the Gaussian, and @f$\sigma_n@f$ is a parameter
 * modelling noise in the detector.
 *
 * for @f$ i@f$ particles this is modified to 
 *
 * @f[ 
 *    f_i(x;\Delta_{p},\xi,\sigma')=f(x;\Delta_{p,i},\xi_i,\sigma_i')
 * @f] 
 *
 * corresponding to @f$ i@f$ particles i.e., with the substitutions 
 * @f{eqnarray*}{ 
 *    \Delta_p  \rightarrow \Delta_{p,i}&=& i(\Delta_p + \xi\log(i))\\
 *    \xi       \rightarrow \xi_i       &=& i \xi\\
 *    \sigma    \rightarrow \sigma_i    &=& \sqrt{i}\sigma\\
 *    \sigma'^2 \rightarrow \sigma_i'^2 &=& \sigma_n^2 + \sigma_i^2
 * @f} 
 *
 * Because of the convolution with a Gaussian, the most-probable-value
 * @f$\Delta_p'@f$ of the resulting distribution is not really at the
 * Landau most-probable-value @f$\Delta_p@f$.  In fact we find that
 * @f$\Delta_p' > \Delta_p@f$.
 *
 * Ideally, one would find an analytic expression for this shift by
 * solving
 *
 * @f{eqnarray}{
 *   0 &=& \frac{d f_i(x;\Delta_p,\xi,\sigma)}{d\Delta}\\
 *     &=& \frac{d}{d\Delta}\frac{1}{\sigma' \sqrt{2 \pi}}
 *      \int_{-\infty}^{+\infty} dx' f'_{L}(x',\Delta_p,\xi)
 *      \exp{-\frac{(x-x')^2}{2\sigma'^2}}
 * @f}
 *
 * for @f$ x@f$ as a function of @f$\Delta_p,\xi,\sigma,i@f$. However,
 * do to the complex nature of the Landau distribution this is not
 * really feasible.
 *
 * Instead, the shift was studied numerically. Landau-Gauss
 * distributions for @f$i=1,\ldots@f$ where generated with varying
 * @f$\xi@f$ and @f$\sigma@f$.  The distributions was then numerically
 * differentiated and the root @f$\Delta_p'@f$ of that derivative
 * found numerically.  The difference
 * @f$\delta\Delta_p=\Delta_p'-\Delta_p@f$ was then studied as a
 * function of the @f$\sigma,\xi@f$ parameters and an approximate
 * expression was found
 *
 * @f[ 
 *  \delta\Delta_p \approx \frac{c \sigma u}{(1+1/i)^{p u^{3/2}}}
 * @f]
 *
 * where @f$ u=\sigma/\xi@f$.  The parameters @f$c@f$ and @f$p@f$ is
 * found to depend on @f$ u@f$ only weakly, and for practical
 * applications where @f$u\approx1@f$, we set @f$ c=p=1/2@f$. 
 *
 * For the evaluating the full energy loss distribution from f@f$
 * 1+2+\ldots,n@f$ particles, we evaluate
 *
 * @f[ 
 *   f_N(x;\Delta_p,\xi,\sigma')=\sum_{i=1}^N a_i f_i(x;\Delta_p,\xi,\sigma',a)
 * @f] 
 * 
 * where @f$ f(x;\Delta_p,\xi,\sigma')@f$ is the convolution of a
 * Landau with a Gaussian (see LandauGaus), and @f$ a@f$ is a vector of
 * weights for each @f$ f_i@f$. Note that @f$ a_1 = 1@f$.
 *
 * Everything is defined in this header file to make it easy to move
 * this code around. Nothing here's meant to be persistent, so we
 * can easily do that. 
 * 
 * References: 
 *  - <a href="http://dx.doi.org/10.1016/0168-583X(84)90472-5">
 *        Nucl.Instrum.Meth.B1:16</a>
 *  - <a href="http://dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
 *  - <a href="http://root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">
 *        ROOT implementation</a>
 * 
 * @ingroup pwglf_forward 
 */
class AliLandauGaus
{
public:
  /** Enumeration of parameters */
  enum { 
    kC     = 0,
    kDelta, 
    kXi, 
    kSigma, 
    kSigmaN, 
    kN, 
    kA
  };
  /** Enumeration of colors */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Constants 
   */
  //------------------------------------------------------------------
  /** 
   * The shift of the most probable value for the ROOT function
   * TMath::Landau
   * 
   * @return  Shift of TMath::Landau
   */
  static Double_t MPShift() { return -0.22278298; }
  /** 
   * Constant of @f$\sigma@f$ shift 
   *
   * @return Constant of @f$\sigma@f$ shift 
   */
  static Double_t SigmaShiftC() { return 0.5446; }
  /** 
   * Power factor of @f$\sigma@f$ shift 
   *
   * @return Power factor of @f$\sigma@f$ shift 
   */
  static Double_t SigmaShiftP() { return 0.5746; }
  /** 
   * Normalization constant 
   * 
   * @return The landau-gauss normalization constant
   */
  static Double_t InvSq2Pi() { return 1. / TMath::Sqrt(TMath::TwoPi()); }
  /** 
   * How many sigma's of the Gaussian in the Landau, Gaussian
   * convolution to integrate over
   */
  static Double_t NSigma() { return 5; }
  /**
   * Number of steps to do in the Landau, Gaussiam convolution 
   */
  static Int_t NSteps() { return 100; }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Function calculations 
   */
  //------------------------------------------------------------------
  /** 
   * Calculate the shifted Landau
   * @f[
   *    f'_{L}(x;\Delta_p,\xi) = f_L(x;\Delta_p+0.22278298\xi)
   * @f]
   *
   * where @f$ f_{L}@f$ is the ROOT implementation of the Landau
   * distribution (known to have @f$\Delta_{p}=-0.22278298@f$ for
   * @f$\Delta_p=0,\xi=1@f$. 
   *
   * @param x      Where to evaluate @f$ f'_{L}@f$ 
   * @param delta  Most probable value 
   * @param xi     The 'width' of the distribution 
   *
   * @return @f$ f'_{L}(x;\Delta,\xi) @f$
   */
  static Double_t Fl(Double_t x, Double_t delta, Double_t xi);
  //------------------------------------------------------------------
  /** 
   * Calculate the value of a Landau convolved with a Gaussian 
   * 
   * @f[ 
   * f(x;\Delta_p,\xi,\sigma') = \frac{1}{\sigma' \sqrt{2 \pi}}
     *    \int_{-\infty}^{+\infty} dx' f'_{L}(x';\Delta_p,\xi)
   *    \exp{-\frac{(x-x')^2}{2\sigma'^2}}
   * @f]
   * 
   * Note that this function uses the constants NSteps() and
   * NSigma()
   * 
   * @param x         where to evaluate @f$ f@f$
   * @param delta     @f$ \Delta_p@f$ of @f$ f(x;\Delta_p,\xi,\sigma')@f$
   * @param xi        @f$ \xi@f$ of @f$ f(x;\Delta_p,\xi,\sigma')@f$
   * @param sigma     @f$ \sigma@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
   * @param sigma_n   @f$ \sigma_n@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
   * 
   * @return @f$ f@f$ evaluated at @f$ x@f$.  
   */
  static Double_t F(Double_t x, Double_t delta, Double_t xi, 
		    Double_t sigma, Double_t sigma_n);
  //------------------------------------------------------------------
  /** 
   * Evaluate 
   * @f[ 
   *    f_i(x;\Delta_p,\xi,\sigma') = f(x;\Delta_{p,i},\xi_i,\sigma_i')
   * @f] 
   * corresponding to @f$ i@f$ particles.
   * 
   * @param x        Where to evaluate 
   * @param delta    @f$ \Delta@f$ 
   * @param xi       @f$ \xi@f$ 
   * @param sigma    @f$ \sigma@f$ 
   * @param sigma_n  @f$ \sigma_n@f$
   * @param i        @f$ i @f$
   * 
   * @return @f$ f_i @f$ evaluated
   */  
  static Double_t Fi(Double_t x, Double_t delta, Double_t xi, 
		     Double_t sigma, Double_t sigma_n, Int_t i);

  //------------------------------------------------------------------
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
   * - 3: @f$\sigma_n@f$ 
   *
   * @param x        Where to evaluate 
   * @param ipar     Parameter number 
   * @param dp       @f$ \epsilon\delta p_i@f$ for some value of @f$\epsilon@f$
   * @param delta    @f$ \Delta_p@f$ 
   * @param xi       @f$ \xi@f$ 
   * @param sigma    @f$ \sigma@f$ 
   * @param sigma_n  @f$ \sigma_n@f$
   * @param i        @f$ i@f$
   * 
   * @return @f$ f_i@f$ evaluated
   */  
  static Double_t DFidPar(Double_t x, UShort_t ipar, Double_t dp,
			  Double_t delta, Double_t xi, 
			  Double_t sigma, Double_t sigma_n, Int_t i);

  //------------------------------------------------------------------
  /** 
   * Evaluate 
   * @f[ 
   * f_N(x;\Delta_p,\xi,\sigma')=\sum_{i=1}^N a_i f_i(x;\Delta_p,\xi,\sigma',a)
   * @f] 
   * 
   * where @f$ f(x;\Delta_p,\xi,\sigma')@f$ is the convolution of a
     * Landau with a Gaussian (see LandauGaus).
     * 
     * @param x        Where to evaluate @f$ f_N@f$
     * @param delta    @f$ \Delta_1@f$ 
     * @param xi       @f$ \xi_1@f$
     * @param sigma    @f$ \sigma_1@f$ 
     * @param sigma_n  @f$ \sigma_n@f$ 
     * @param n        @f$ N@f$ in the sum above.
     * @param a        Array of size @f$ N-1@f$ of the weights @f$ a_i@f$ for 
     *                 @f$ i > 1@f$ 
     * 
     * @return @f$ f_N(x;\Delta,\xi,\sigma')@f$ 
   */
  static Double_t Fn(Double_t x, Double_t delta, Double_t xi, 
		     Double_t sigma, Double_t sigma_n, Int_t n, 
		     const Double_t* a);
  /** 
   * Get parameters for the @f$ i@f$ particle response.
   *
   * @f{eqnarray*}{ 
   *    \Delta    \rightarrow \Delta_i    &=& i(\Delta + \xi\log(i))\\
   *    \xi       \rightarrow \xi_i       &=& i \xi\\
   *    \sigma    \rightarrow \sigma_i    &=& \sqrt{i}\sigma\\
   *    \sigma'^2 \rightarrow \sigma_i'^2 &=& \sigma_n^2 + \sigma_i^2
   * @f} 
   * 
   *
   * @param i     Number of particles 
   * @param delta Input: single particle, output @f$ i@f$ particle
   * @param xi    Input: single particle, output @f$ i@f$ particle
   * @param sigma Input: single particle, output @f$ i@f$ particle
   */
  static void IPars(Int_t i, Double_t& delta, Double_t& xi, Double_t& sigma);
  /** 
   * Set and check if sigma shift is enabled 
   * 
   * @param val if <0, then only check.  Otherwise set enabled (>0) or not (=0)
   * 
   * @return whether the sigma shift is enabled or not 
   */
  static Bool_t EnableSigmaShift(Short_t val=-1);
  /** 
   * Get the shift of the MPV due to convolution with a Gaussian. 
   *
   * @f[ 
   *  \delta\Delta_p \approx \frac{c \sigma u}{(1+1/i)^{p u^{3/2}}}
   * @f]
   *
   * where @f$ u=\sigma/\xi@f$.  
   * 
   * @param i       Number of particles 
   * @param xi      Landau width @f$\xi@f$
   * @param sigma   Gaussian variance @f$\sigma@f$ 
   * 
   * @return The shift 
   */
  static Double_t SigmaShift(Int_t i, Double_t xi, Double_t sigma);
  /* @} */

  
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Utilities for defining TF1 objects 
   */
  //------------------------------------------------------------------
  /** 
   * Generate a TF1 object of @f$ f_1@f$ 
   * 
   * @param c        Constant
   * @param delta    @f$ \Delta_1@f$ 
   * @param xi       @f$ \xi_1@f$	       
   * @param sigma    @f$ \sigma_1@f$ 	       
   * @param sigma_n  @f$ \sigma_n@f$ 	       
   * @param xmin     Least value of range
   * @param xmax     Largest value of range
   * 
   * @return Newly allocated TF1 object
   */
  static TF1* MakeF1(Double_t c, 
		     Double_t delta, Double_t xi, 
		     Double_t sigma, Double_t sigma_n,
		     Double_t xmin,  Double_t  xmax);
  //------------------------------------------------------------------
  /** 
   * Generate a TF1 object of @f$ f_I@f$ 
   * 
   * @param c        Constant
   * @param delta    @f$ \Delta_1@f$ 
   * @param xi       @f$ \xi_1@f$	       
   * @param sigma    @f$ \sigma_1@f$ 	       
   * @param sigma_n  @f$ \sigma_n@f$ 	       
   * @param i 	     @f$ i@f$ - the number of particles
   * @param xmin     Least value of range
   * @param xmax     Largest value of range
   * 
   * @return Newly allocated TF1 object
   */
  static TF1* MakeFi(Double_t c, 
		     Double_t delta, Double_t xi, 
		     Double_t sigma, Double_t sigma_n,
		     Int_t    i, 
		     Double_t xmin,  Double_t  xmax);
  //------------------------------------------------------------------
  /** 
   * Generate a TF1 object of @f$ f_N@f$ 
   * 
   * @param c         Constant			       
   * @param delta     @f$ \Delta_1@f$ 		       
   * @param xi 	      @f$ \xi_1@f$	       	       
   * @param sigma     @f$ \sigma_1@f$ 	       	       
   * @param sigma_n   @f$ \sigma_n@f$ 	       	       
   * @param n 	      @f$ N@f$ - how many particles to sum to
   * @param a         Array of size @f$ N-1@f$ of the weights @f$ a_i@f$ for 
   *                  @f$ i > 1@f$ 
   * @param xmin      Least value of range  
   * @param xmax      Largest value of range
   * 
   * @return Newly allocated TF1 object
   */
  static TF1* MakeFn(Double_t c, 
		     Double_t delta, Double_t  xi, 
		     Double_t sigma, Double_t  sigma_n,
		     Int_t    n,     const Double_t* a, 
		     Double_t xmin,  Double_t  xmax);
  /** 
   * Make a TF1 object that describes the convoluted energy loss of
   * primary and secondary particles.
   * 
   * @param c1     Primary weight
   * @param delta  Most-probable value 
   * @param xi1    Primary width
   * @param sigma  Gaussian variance 
   * @param c2     Secondary weight
   * @param xi2    Secondary with
   * @param xmin   Least x
   * @param xmax   Largest x
   * 
   * @return Pointer to newly allocated TF1 object
   */
  static TF1* MakeComposite(Double_t c1, 
			    Double_t delta, 
			    Double_t xi1, 
			    Double_t sigma, 
			    Double_t c2,
			    Double_t xi2, 
			    Double_t xmin, 
			    Double_t xmax);
  
  //------------------------------------------------------------------
  /** 
   * Get the color of the @f$ i@f$ particle response 
   * 
   * @param i Particle number 
   * 
   * @return Color 
   */
  static Color_t GetIColor(Int_t i);
  //------------------------------------------------------------------
  /** 
   * Utility function for TF1 definition 
   *
   * @param pp Array of parameters 
   * @param xp Pointer to independent variable 
   *
   * @see AliLandauGaus::F 
   *
   * @return Landau convolved with a Gauss
   */
  static Double_t F1Func(Double_t* xp, Double_t* pp);
  //------------------------------------------------------------------
  /** 
   * Utility function for TF1 definition 
   *
   * @param pp Array of parameters 
   * @param xp Pointer to independent variable 
   *
   * @see AliLandauGaus::Fn 
   *
   * @return Landau convolved with a Gauss
   */
  static Double_t FnFunc(Double_t* xp, Double_t* pp);
  //------------------------------------------------------------------
  /** 
   * Utility function for TF1 definition 
   *
   * @param pp Array of parameters 
   * @param xp Pointer to independent variable 
   *
   * @see AliLandauGaus::Fn 
   *
   * @return Landau convolved with a Gauss
   */
  static Double_t FiFunc(Double_t* xp, Double_t* pp);
  //------------------------------------------------------------------
  /** 
   * Utility function for TF1 definition 
   *
   * @param pp Array of parameters 
   * @param xp Pointer to independent variable 
   *
   * @see AliLandauGaus::F 
   *
   * @return Landau convolved with a Gauss
   */
  static Double_t CompFunc(Double_t* xp, Double_t* pp);
  /* @} */
};
//____________________________________________________________________
inline Bool_t
AliLandauGaus::EnableSigmaShift(Short_t val)
{
  static Bool_t enabled = true;
  if (val >= 0) enabled = val == 1;
  return enabled;
}
//____________________________________________________________________
inline void
AliLandauGaus::IPars(Int_t i, Double_t& delta, Double_t& xi, Double_t& sigma)
{
#ifndef NO_SIGMA_SHIFT
  Double_t dDelta = EnableSigmaShift() ? SigmaShift(i, xi, sigma) : 0;
#else 
  // This is for testing 
  Double_t dDelta = 0;
#endif
  if (i == 1) {
    delta -= dDelta;
    return;
  }
  
  delta = i * (delta + xi * TMath::Log(i)) - dDelta;
  xi    = i * xi;
  sigma = TMath::Sqrt(Double_t(i)) * sigma;
}
//____________________________________________________________________
inline Double_t
AliLandauGaus::SigmaShift(Int_t i, Double_t xi, Double_t sigma)
{
  if (xi <= 0) return 0;
  if (sigma <= 0) return 0;
  const Double_t c = SigmaShiftC();
  const Double_t p = SigmaShiftP();
  const Double_t u = sigma / xi;
  const Double_t q = p*u*TMath::Sqrt(u);
  if (q > 100) return 0;
  return c * sigma / TMath::Power(1+1./i, q);
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::Fl(Double_t x, Double_t delta, Double_t xi)
{
  Double_t deltaP = delta - xi * MPShift();
  return TMath::Landau(x, deltaP, xi, true);
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::F(Double_t x, Double_t delta, Double_t xi,
		 Double_t sigma, Double_t sigmaN)
{
  if (xi <= 0) return 0;

  const Int_t    nSteps = NSteps();
  const Double_t nSigma = NSigma();
  const Double_t deltaP = delta; // - sigma * sigmaShift; // + sigma * mpshift;
  const Double_t sigma2 = sigmaN*sigmaN + sigma*sigma;
  const Double_t sigma1 = sigmaN == 0 ? sigma : TMath::Sqrt(sigma2);
  const Double_t xlow   = x - nSigma * sigma1;
  const Double_t xhigh  = x + nSigma * sigma1;
  const Double_t step   = (xhigh - xlow) / nSteps;
  Double_t       sum    = 0;
  
  for (Int_t i = 0; i <= nSteps/2; i++) { 
    const Double_t x1 = xlow  + (i - .5) * step;
    const Double_t x2 = xhigh - (i - .5) * step;
    sum += Fl(x1, deltaP, xi) * TMath::Gaus(x, x1, sigma1);
    sum += Fl(x2, deltaP, xi) * TMath::Gaus(x, x2, sigma1);
  }
  return step * sum * InvSq2Pi() / sigma1;
}

//____________________________________________________________________
inline Double_t 
AliLandauGaus::Fi(Double_t x, Double_t delta, Double_t xi, 
		  Double_t sigma, Double_t sigmaN, Int_t i)
{
  Double_t deltaI = delta;
  Double_t xiI    = xi;
  Double_t sigmaI = sigma;
  IPars(i, deltaI, xiI, sigmaI);
  if (sigmaI < 1e-10) 
    // Fall back to landau 
    return Fl(x, deltaI, xiI);
  
  return F(x, deltaI, xiI, sigmaI, sigmaN);
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::Fn(Double_t x, Double_t delta, Double_t xi, 
		  Double_t sigma, Double_t sigmaN, Int_t n, 
		  const Double_t* a)
{
  Double_t result = Fi(x, delta, xi, sigma, sigmaN, 1);
  for (Int_t i = 2; i <= n; i++) 
    result += a[i-2] * Fi(x,delta,xi,sigma,sigmaN,i);
  return result;
}

//____________________________________________________________________
inline Double_t 
AliLandauGaus::DFidPar(Double_t x, 
		       UShort_t par,   Double_t dPar, 
		       Double_t delta, Double_t xi, 
		       Double_t sigma, Double_t sigmaN, 
		       Int_t    i)
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
  switch (par) {
  case 0: 
    y1 = Fi(x, deltaI+i*dp, xiI, sigmaI, sigmaN, i);
    y2 = Fi(x, deltaI+i*d2, xiI, sigmaI, sigmaN, i);
    y3 = Fi(x, deltaI-i*d2, xiI, sigmaI, sigmaN, i);
    y4 = Fi(x, deltaI-i*dp, xiI, sigmaI, sigmaN, i);
    break;
  case 1: 
    y1 = Fi(x, deltaI, xiI+i*dp, sigmaI, sigmaN, i);
    y2 = Fi(x, deltaI, xiI+i*d2, sigmaI, sigmaN, i);
    y3 = Fi(x, deltaI, xiI-i*d2, sigmaI, sigmaN, i);
    y4 = Fi(x, deltaI, xiI-i*dp, sigmaI, sigmaN, i);
    break;
  case 2: 
    y1 = Fi(x, deltaI, xiI, sigmaI+si*dp, sigmaN, i);
    y2 = Fi(x, deltaI, xiI, sigmaI+si*d2, sigmaN, i);
    y3 = Fi(x, deltaI, xiI, sigmaI-si*d2, sigmaN, i);
    y4 = Fi(x, deltaI, xiI, sigmaI-si*dp, sigmaN, i);
    break;
  case 3: 
    y1 = Fi(x, deltaI, xiI, sigmaI, sigmaN+dp, i);
    y2 = Fi(x, deltaI, xiI, sigmaI, sigmaN+d2, i);
    y3 = Fi(x, deltaI, xiI, sigmaI, sigmaN-d2, i);
    y4 = Fi(x, deltaI, xiI, sigmaI, sigmaN-dp, i);
    break;
  default:
    return 0;
  } 
  
  Double_t d0  = y1 - y4;
  Double_t d1  = 2 * (y2 - y3);
  
  Double_t g   = 1/(2*dp) * (4*d1 - d0) / 3;
   
  return g;
}


//____________________________________________________________________
inline Color_t
AliLandauGaus::GetIColor(Int_t i) 
{
  const Int_t kColors[] = { kRed+1, 
			    kPink+3, 
			    kMagenta+2, 
			    kViolet+2, 
			    kBlue+1, 
			    kAzure+3, 
			    kCyan+1, 
			    kTeal+2, 
			    kGreen+2, 
			    kSpring+3, 
			    kYellow+2, 
			    kOrange+2 };
  return kColors[((i-1) % 12)];
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::F1Func(Double_t* xp, Double_t* pp) 
{ 
  Double_t x        = xp[0];
  Double_t constant = pp[kC];
  Double_t delta    = pp[kDelta];
  Double_t xi       = pp[kXi];
  Double_t sigma    = pp[kSigma];
  Double_t sigmaN   = pp[kSigmaN];
  
  return constant * F(x, delta, xi, sigma, sigmaN);
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::FiFunc(Double_t* xp, Double_t* pp) 
{
  Double_t x        = xp[0];
  Double_t constant = pp[kC];
  Double_t delta    = pp[kDelta];
  Double_t xi       = pp[kXi];
  Double_t sigma    = pp[kSigma];
  Double_t sigmaN   = pp[kSigmaN];
  Int_t    i        = Int_t(pp[kN]);

  return constant * Fi(x, delta, xi, sigma, sigmaN, i);
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::FnFunc(Double_t* xp, Double_t* pp) 
{ 
  Double_t  x        = xp[0];
  Double_t constant  = pp[kC];
  Double_t delta     = pp[kDelta];
  Double_t xi        = pp[kXi];
  Double_t sigma     = pp[kSigma];
  Double_t sigmaN    = pp[kSigmaN];
  Int_t     n        = Int_t(pp[kN]);
  Double_t* a        = &(pp[kA]);
  
  return constant * Fn(x, delta, xi, sigma, sigmaN, n, a);
}
//____________________________________________________________________
inline Double_t 
AliLandauGaus::CompFunc(Double_t* xp, Double_t* pp) 
{
  Double_t x           = xp[0];
  Double_t cP          = pp[kC];
  Double_t deltaP      = pp[kDelta];
  Double_t xiP         = pp[kXi];
  Double_t sigmaP      = pp[kSigma];
  Double_t cS          = pp[kSigma+1];
  Double_t deltaS      = deltaP; // pp[kSigma+2];
  Double_t xiS         = pp[kSigma+2/*3*/];
  Double_t sigmaS      = sigmaP; // pp[kSigma+4];
  
  return (cP * F(x,deltaP,xiP,sigmaP,0) + 
	  cS * F(x,deltaS,xiS,sigmaS,0));

}
//____________________________________________________________________
inline TF1*
AliLandauGaus::MakeF1(Double_t  c, 
		      Double_t  delta, Double_t xi, 
		      Double_t  sigma, Double_t sigmaN,
		      Double_t  xmin, Double_t xmax)
{
  // Define the function to fit 
  TF1* f = new TF1("landau1", &F1Func, xmin,xmax,kSigmaN+1);

  // Set initial guesses, parameter names, and limits  
  f->SetParameters(c,delta,xi,sigma,sigmaN);
  f->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}");
  f->SetNpx(500);
  f->SetLineColor(GetIColor(1));
  
  return f;
}
  
//____________________________________________________________________
inline TF1*
AliLandauGaus::MakeFn(Double_t  c, 
		      Double_t  delta, Double_t xi, 
		      Double_t  sigma, Double_t sigmaN, Int_t n, 
		      const Double_t* a, 
		      Double_t  xmin, Double_t xmax)
{
  Int_t npar = kN+n;
  TF1*  f    = new TF1(Form("nlandau%d", n), &FnFunc,xmin,xmax,npar);
  f->SetLineColor(GetIColor(n)); 
  f->SetLineWidth(2);
  f->SetNpx(500);
  f->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "N");
  f->SetParameter(kC,      c);       
  f->SetParameter(kDelta,  delta);   
  f->SetParameter(kXi,     xi);      
  f->SetParameter(kSigma,  sigma);   
  f->SetParameter(kSigmaN, sigmaN); 
  f->FixParameter(kN,      n);       
  for (UShort_t i = 2; i <= n; i++) {
    f->SetParameter(kA+i-2, a[i-2]);
    f->SetParName(kA+i-2, Form("a_{%d}", i));
  }
  return f;
}
//____________________________________________________________________
inline TF1*
AliLandauGaus::MakeFi(Double_t  c, 
		      Double_t  delta, Double_t xi, 
		      Double_t  sigma, Double_t sigmaN, Int_t i, 
		      Double_t  xmin, Double_t xmax)
{
  Int_t npar = kN+1;
  TF1*  f    = new TF1(Form("ilandau%d", i), &FiFunc,xmin,xmax,npar);
  f->SetLineColor(GetIColor(i));
  f->SetLineWidth(1);
  f->SetNpx(500);
  f->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "i");
  f->SetParameter(kC,      c);       
  f->SetParameter(kDelta,  delta);   
  f->SetParameter(kXi,     xi);      
  f->SetParameter(kSigma,  sigma);   
  f->SetParameter(kSigmaN, sigmaN); 
  f->FixParameter(kN,      i);       

  return f;
}
//____________________________________________________________________
inline TF1*
AliLandauGaus::MakeComposite(Double_t c1, 
			     Double_t delta, 
			     Double_t xi1, 
			     Double_t sigma, 
			     Double_t c2, 
			     Double_t xi2, 
			     Double_t xmin, 
			     Double_t xmax)
{
  TF1* comp = new TF1("composite", &CompFunc, xmin, xmax, kSigma+1+2);
  comp->SetParNames("C",       "#Delta_{p}",       "#xi",       "#sigma",
		    "C#prime", "#xi#prime");
  comp->SetParameters(c1,     // 0 Primary weight 
		      delta,  // 1 Primary Delta
		      xi1,    // 2 primary Xi
		      sigma,  // 3 primary sigma
		      c2,     // 5 Secondary weight
		      xi2);   // 6 secondary Xi
  return comp;
}
#endif
// Local Variables:
//  mode: C++ 
// End:



