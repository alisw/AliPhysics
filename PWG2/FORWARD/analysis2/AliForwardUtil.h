#ifndef ALIROOT_PWG2_FORWARD_ALIFORWARDUTIL_H
#define ALIROOT_PWG2_FORWARD_ALIFORWARDUTIL_H
#include <TObject.h>
#include <TString.h>
#include <TObjArray.h>
class TH2D;
class TH1I;
class TH1;
class TF1;
class TAxis;
class AliESDEvent;

/** 
 * Utilities used in the forward multiplcity analysis 
 * 
 * @ingroup pwg2_forward_analysis 
 */
class AliForwardUtil : public TObject
{
public:
  //__________________________________________________________________
  /**
   * Number of steps to do in the Landau, Gaussiam convolution 
   */
  static Int_t fgConvolutionSteps;
  //------------------------------------------------------------------
  /** 
   * How many sigma's of the Gaussian in the Landau, Gaussian
   * convolution to integrate over
   */
  static Double_t fgConvolutionNSigma;
  //------------------------------------------------------------------
  /** 
   * Calculate the shifted Landau
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
  static Double_t Landau(Double_t x, Double_t delta, Double_t xi);
  
  //------------------------------------------------------------------
  /** 
   * Calculate the value of a Landau convolved with a Gaussian 
   * 
   * @f[ 
   * f(x;\Delta,\xi,\sigma') = \frac{1}{\sigma' \sqrt{2 \pi}}
   *    \int_{-\infty}^{+\infty} d\Delta' f'_{L}(x;\Delta',\xi)
   *    \exp{-\frac{(\Delta-\Delta')^2}{2\sigma'^2}}
   * @f]
   * 
   * where @f$ f'_{L}@f$ is the Landau distribution, @f$ \Delta@f$ the
   * energy loss, @f$ \xi@f$ the width of the Landau, and 
   * @f$ \sigma'^2=\sigma^2-\sigma_n^2 @f$.  Here, @f$\sigma@f$ is the
   * variance of the Gaussian, and @f$\sigma_n@f$ is a parameter modelling 
   * noise in the detector.  
   *
   * Note that this function uses the constants fgConvolutionSteps and
   * fgConvolutionNSigma
   * 
   * References: 
   *  - <a href="http://dx.doi.org/10.1016/0168-583X(84)90472-5">Nucl.Instrum.Meth.B1:16</a>
   *  - <a href="http://dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
   *  - <a href="http://root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">ROOT implementation</a>
   * 
   * @param x         where to evaluate @f$ f@f$
   * @param delta     @f$ \Delta@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$
   * @param xi        @f$ \xi@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$
   * @param sigma     @f$ \sigma@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
   * @param sigma_n   @f$ \sigma_n@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
   * 
   * @return @f$ f@f$ evaluated at @f$ x@f$.  
   */
  static Double_t LandauGaus(Double_t x, Double_t delta, Double_t xi, 
			     Double_t sigma, Double_t sigma_n);
  
  //------------------------------------------------------------------
  /** 
   * Evaluate 
   * @f[ 
      f_N(x;\Delta,\xi,\sigma') = \sum_{i=1}^N a_i f(x;\Delta_i,\xi_i,\sigma'_i)
     @f] 
   * 
   * where @f$ f(x;\Delta,\xi,\sigma')@f$ is the convolution of a
   * Landau with a Gaussian (see LandauGaus).  Note that 
   * @f$ a_1 = 1@f$, @f$\Delta_i = i(\Delta_1 + \xi\log(i))@f$, 
   * @f$\xi_i=i\xi_1@f$, and @f$\sigma_i'^2 = \sigma_n^2 + i\sigma_1^2@f$. 
   *  
   * References: 
   *  - <a href="http://dx.doi.org/10.1016/0168-583X(84)90472-5">Nucl.Instrum.Meth.B1:16</a>
   *  - <a href="http://dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
   *  - <a href="http://root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">ROOT implementation</a>
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
  static Double_t NLandauGaus(Double_t x, Double_t delta, Double_t xi, 
			      Double_t sigma, Double_t sigma_n, Int_t n, 
			      Double_t* a);
  //__________________________________________________________________
  /** 
   * Structure to do fits to the energy loss spectrum 
   * 
   */
  struct ELossFitter 
  {
    enum { 
      kC     = 0,
      kDelta, 
      kXi, 
      kSigma, 
      kSigmaN, 
      kN, 
      kA
    };
    /** 
     * Constructor 
     * 
     * @param lowCut     Lower cut of spectrum - data below this cuts is ignored
     * @param maxRange   Maximum range to fit to 
     * @param minusBins  The number of bins below maximum to use 
     */
    ELossFitter(Double_t lowCut, Double_t maxRange, UShort_t minusBins); 
    virtual ~ELossFitter();
    /** 
     * Clear internal arrays 
     * 
     */
    void Clear();
    /** 
     * Fit a 1-particle signal to the passed energy loss distribution 
     * 
     * Note that this function clears the internal arrays first 
     *
     * @param dist    Data to fit the function to 
     * @param sigman If larger than zero, the initial guess of the
     *               detector induced noise. If zero or less, then this 
     *               parameter is ignored in the fit (fixed at 0)
     * 
     * @return The function fitted to the data 
     */
    TF1* Fit1Particle(TH1* dist, Double_t sigman=-1);
    /** 
     * Fit a N-particle signal to the passed energy loss distribution 
     *
     * If there's no 1-particle fit present, it does that first 
     *
     * @param dist   Data to fit the function to 
     * @param n      Number of particle signals to fit 
     * @param sigman If larger than zero, the initial guess of the
     *               detector induced noise. If zero or less, then this 
     *               parameter is ignored in the fit (fixed at 0)
     * 
     * @return The function fitted to the data 
     */
    TF1* FitNParticle(TH1* dist, UShort_t n, Double_t sigman=-1);
     

    const Double_t fLowCut;     // Lower cut on data 
    const Double_t fMaxRange;   // Maximum range to fit 
    const UShort_t fMinusBins;  // Number of bins from maximum to fit 1st peak
    TObjArray fFitResults;      // Array of fit results 
    TObjArray fFunctions;       // Array of functions 
  };
      

  //__________________________________________________________________
  /** 
   * Structure to hold histograms 
   *
   * @ingroup pwg2_forward_analysis 
   */
  struct Histos : public TObject
  {	
    /** 
     * Constructor 
     * 
     * 
     */
    Histos() : fFMD1i(0), fFMD2i(0), fFMD2o(0), fFMD3i(0), fFMD3o(0) {}
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    Histos(const Histos& o) 
      : TObject(o), 
	fFMD1i(o.fFMD1i), 
	fFMD2i(o.fFMD2i), 
	fFMD2o(o.fFMD2o), 
	fFMD3i(o.fFMD3i), 
	fFMD3o(o.fFMD3o)
    {}
    /** 
     * Assignement operator 
     * 
     * @return Reference to this 
     */
    Histos& operator=(const Histos&) { return *this;}
    /** 
     * Destructor
     */
    ~Histos();
    /** 
     * Initialize the object 
     * 
     * @param etaAxis Eta axis to use 
     */
    void Init(const TAxis& etaAxis);
    /** 
     * Make a histogram 
     * 
     * @param d        Detector
     * @param r        Ring 
     * @param etaAxis  Eta axis to use
     * 
     * @return Newly allocated histogram 
     */
    TH2D* Make(UShort_t d, Char_t r, const TAxis& etaAxis) const;
    /** 
     * Clear data 
     * 
     * @param option Not used 
     */
    void  Clear(Option_t* option="");
    // const TH2D* Get(UShort_t d, Char_t r) const;
    /** 
     * Get the histogram for a particular detector,ring
     * 
     * @param d Detector 
     * @param r Ring 
     * 
     * @return Histogram for detector,ring or nul 
     */
    TH2D* Get(UShort_t d, Char_t r) const;
    TH2D* fFMD1i; // Histogram for FMD1i
    TH2D* fFMD2i; // Histogram for FMD2i
    TH2D* fFMD2o; // Histogram for FMD2o
    TH2D* fFMD3i; // Histogram for FMD3i
    TH2D* fFMD3o; // Histogram for FMD3o

    ClassDef(Histos,1) 
  };

  //__________________________________________________________________
  struct RingHistos : public TObject
  {
    RingHistos() : fDet(0), fRing('\0'), fName("") {}
    RingHistos(UShort_t d, Char_t r) 
      : fDet(d), fRing(r), fName(TString::Format("FMD%d%c", d, r)) 
    {}
    RingHistos(const RingHistos& o) 
      : TObject(o), fDet(o.fDet), fRing(o.fRing), fName(o.fName)
    {}
    virtual ~RingHistos() {}
    RingHistos& operator=(const RingHistos& o) 
    {
      TObject::operator=(o);
      fDet  = o.fDet;
      fRing = o.fRing;
      fName = o.fName;
      return *this;
    }
    TList* DefineOutputList(TList* d) const;
    TList* GetOutputList(TList* d) const;
    TH1* GetOutputHist(TList* d, const char* name) const;
    Color_t Color() const 
    { 
      return ((fDet == 1 ? kRed : (fDet == 2 ? kGreen : kBlue))
	      + ((fRing == 'I' || fRing == 'i') ? 2 : -2));
    }

    UShort_t fDet; 
    Char_t   fRing;
    TString  fName;

    ClassDef(RingHistos,1) 
  };
    
};

#endif
// Local Variables:
//  mode: C++
// End:

