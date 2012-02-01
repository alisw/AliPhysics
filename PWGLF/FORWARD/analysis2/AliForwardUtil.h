// 
// Utilities used in the forward multiplcity analysis 
// 
//
#ifndef ALIFORWARDUTIL_H
#define ALIFORWARDUTIL_H
/**
 * @file   AliForwardUtil.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:54 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward 
 */
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
 * @ingroup pwg2_forward 
 */
class AliForwardUtil : public TObject
{
public:
  /** 
   * Get the standard color for a ring  
   *
   * @param d Detector
   * @param r Ring 
   * 
   * @return 
   */
  static Color_t RingColor(UShort_t d, Char_t r)
  { 
    return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	    + ((r == 'I' || r == 'i') ? 2 : -3));
  }
  //==================================================================
  /** 
   * @{ 
   * @name Collision/run parameters 
   */
  /**						
   * Defined collision types 
   */
  enum ECollisionSystem {
    kUnknown, 
    kPP, 
    kPbPb,
    kPPb
  };
  //__________________________________________________________________
  /** 
   * Parse a collision system spec given in a string.   Known values are 
   * 
   *  - "pp", "p-p" which returns kPP 
   *  - "PbPb", "Pb-Pb", "A-A", which returns kPbPb
   *  - "pPb", "p-Pb", "pA", p-A" which returns kPPb
   *  - Everything else gives kUnknown 
   * 
   * @param sys Collision system spec 
   * 
   * @return Collision system id 
   */
  static UShort_t ParseCollisionSystem(const char* sys);
  /** 
   * Get a string representation of the collision system 
   * 
   * @param sys  Collision system 
   * - kPP -> "pp"
   * - kPbPb -> "PbPb" 
   * - kPPb -> "pPb"
   * - anything else gives "unknown"
   * 
   * @return String representation of the collision system 
   */
  static const char* CollisionSystemString(UShort_t sys);
  //__________________________________________________________________
  /** 
   * Parse the center of mass energy given as a float and return known 
   * values as a unsigned integer
   * 
   * @param sys  Collision system (needed for AA)
   * @param cms  Center of mass energy * total charge 
   * 
   * @return Center of mass energy per nucleon
   */
  static UShort_t ParseCenterOfMassEnergy(UShort_t sys, Float_t cms);
  /** 
   * Get a string representation of the center of mass energy per nuclean
   * 
   * @param cms  Center of mass energy per nucleon
   * 
   * @return String representation of the center of mass energy per nuclean
   */
  static const char* CenterOfMassEnergyString(UShort_t cms);
  //__________________________________________________________________
  /** 
   * Parse the magnetic field (in kG) as given by a floating point number
   * 
   * @param field  Magnetic field in kG 
   * 
   * @return Short integer value of magnetic field in kG 
   */
  static Short_t ParseMagneticField(Float_t field);
  /** 
   * Get a string representation of the magnetic field
   * 
   * @param field Magnetic field in kG
   * 
   * @return String representation of the magnetic field
   */
  static const char* MagneticFieldString(Short_t field);
  /* @} */

  /** 
   * @{ 
   * @name Energy stragling functions 
   */
  //__________________________________________________________________
  /**
   * Number of steps to do in the Landau, Gaussiam convolution 
   */
  static Int_t fgConvolutionSteps; // Number of convolution steps
  //------------------------------------------------------------------
  /** 
   * How many sigma's of the Gaussian in the Landau, Gaussian
   * convolution to integrate over
   */
  static Double_t fgConvolutionNSigma; // Number of convolution sigmas 
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
   *    f_i(x;\Delta,\xi,\sigma') = f(x;\Delta_i,\xi_i,\sigma_i')
   * @f] 
   * corresponding to @f$ i@f$ particles i.e., with the substitutions 
   * @f{eqnarray*}{ 
   *    \Delta    \rightarrow \Delta_i    &=& i(\Delta + \xi\log(i))\\
   *    \xi       \rightarrow \xi_i       &=& i \xi\\
   *    \sigma    \rightarrow \sigma_i    &=& \sqrt{i}\sigma\\
   *    \sigma'^2 \rightarrow \sigma_i'^2 &=& \sigma_n^2 + \sigma_i^2
   * @f} 
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
  static Double_t ILandauGaus(Double_t x, Double_t delta, Double_t xi, 
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
   * This is the partial derivative with respect to the parameter of
   * the response function corresponding to @f$ i@f$ particles i.e.,
   * with the substitutions
   * @f[ 
   *    \Delta    \rightarrow \Delta_i    = i(\Delta + \xi\log(i))\\
   *    \xi       \rightarrow \xi_i       = i \xi\\
   *    \sigma    \rightarrow \sigma_i    = \sqrt{i}\sigma\\
   *    \sigma'^2 \rightarrow \sigma_i'^2 = \sigma_n^2 + \sigma_i^2
   * @f] 
   * 
   * @param x        Where to evaluate 
   * @param ipar     Parameter number 
   * @param dp       @f$ \epsilon\delta p_i@f$ for some value of @f$\epsilon@f$
   * @param delta    @f$ \Delta@f$ 
   * @param xi       @f$ \xi@f$ 
   * @param sigma    @f$ \sigma@f$ 
   * @param sigma_n  @f$ \sigma_n@f$
   * @param i        @f$ i@f$
   * 
   * @return @f$ f_i@f$ evaluated
   */  
  static Double_t IdLandauGausdPar(Double_t x, UShort_t ipar, Double_t dp,
				   Double_t delta, Double_t xi, 
				   Double_t sigma, Double_t sigma_n, Int_t i);

  //------------------------------------------------------------------
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
			      const Double_t* a);
  /** 
   * Generate a TF1 object of @f$ f_I@f$ 
   * 
   * @param c        Constant
   * @param delta    @f$ \Delta@f$ 
   * @param xi       @f$ \xi_1@f$	       
   * @param sigma    @f$ \sigma_1@f$ 	       
   * @param sigma_n  @f$ \sigma_n@f$ 	       
   * @param i 	     @f$ i@f$ - the number of particles
   * @param xmin     Least value of range
   * @param xmax     Largest value of range
   * 
   * @return Newly allocated TF1 object
   */
  static TF1* MakeILandauGaus(Double_t c, 
			      Double_t delta, Double_t xi, 
			      Double_t sigma, Double_t sigma_n,
			      Int_t    i, 
			      Double_t xmin,  Double_t  xmax);
  /** 
   * Generate a TF1 object of @f$ f_N@f$ 
   * 
   * @param c         Constant			       
   * @param delta     @f$ \Delta@f$ 		       
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
  static TF1* MakeNLandauGaus(Double_t c, 
			      Double_t delta, Double_t  xi, 
			      Double_t sigma, Double_t  sigma_n,
			      Int_t    n,     const Double_t* a, 
			      Double_t xmin,  Double_t  xmax);
			    			    
  //__________________________________________________________________
  /** 
   * Structure to do fits to the energy loss spectrum 
   * 
   * @ingroup pwg2_forward 
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
    /** 
     * Destructor
     * 
     */
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
    /**
     * Get Lower cut on data 
     *
     * @return Lower cut on data 
     */
    Double_t GetLowCut() const { return fLowCut; }
    /**
     * Get Maximum range to fit 
     *
     * @return Maximum range to fit 
     */
    Double_t GetMaxRange() const { return fMaxRange; }
    /**
     * Get Number of bins from maximum to fit 1st peak
     *
     * @return Number of bins from maximum to fit 1st peak
     */
    UShort_t GetMinusBins() const { return fMinusBins; }
    /**
     * Get Array of fit results 
     *
     * @return Array of fit results 
     */
    const TObjArray& GetFitResults() const { return fFitResults; }
    /** 
     * Get Array of fit results  
     *
     * @return Array of fit results 
     */
    TObjArray& GetFitResults() { return fFitResults; }
    /**
     * Get Array of functions 
     *
     * @return Array of functions 
     */
    const TObjArray& GetFunctions() const { return fFunctions; }
    /** 
     * Get Array of functions  
     *
     * @return Array of functions 
     */
    TObjArray& GetFunctions() { return fFunctions; }
  private:
    const Double_t fLowCut;     // Lower cut on data 
    const Double_t fMaxRange;   // Maximum range to fit 
    const UShort_t fMinusBins;  // Number of bins from maximum to fit 1st peak
    TObjArray fFitResults;      // Array of fit results 
    TObjArray fFunctions;       // Array of functions 
  };
  /* @} */
      

  //==================================================================
  /** 
   * @{
   * @name Convenience containers 
   */
  /** 
   * Structure to hold histograms 
   *
   * @ingroup pwg2_forward 
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
  /**
   * Base class for structure holding ring specific histograms
   * 
   * @ingroup pwg2_forward 
   */
  struct RingHistos : public TObject
  {
    /** 
     * Constructor
     * 
     */
    RingHistos() : fDet(0), fRing('\0'), fName("") {}
    /** 
     * 
     * 
     * @param d Detector
     * @param r Ring 
     */
    RingHistos(UShort_t d, Char_t r) 
      : fDet(d), fRing(r), fName(TString::Format("FMD%d%c", d, r)) 
    {}
    /** 
     * Copy constructor
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o) 
      : TObject(o), fDet(o.fDet), fRing(o.fRing), fName(o.fName)
    {}
    /** 
     * 
     */
    virtual ~RingHistos() {}
    /** 
     * Assignement operator
     * 
     * @param o Object to assign from
     * 
     * @return Reference to this
     */
    RingHistos& operator=(const RingHistos& o) 
    {
      if (&o == this) return *this;
      TObject::operator=(o);
      fDet  = o.fDet;
      fRing = o.fRing;
      fName = o.fName;
      return *this;
    }
    /** 
     * Define the outout list in @a d
     * 
     * @param d Where to put the output list
     * 
     * @return Newly allocated TList object or null
     */
    TList* DefineOutputList(TList* d) const;
    /** 
     * Get our output list from the container @a d
     * 
     * @param d where to get the output list from 
     * 
     * @return The found TList or null
     */
    TList* GetOutputList(const TList* d) const;
    /** 
     * Find a specific histogram in the source list @a d
     * 
     * @param d     (top)-container 
     * @param name  Name of histogram
     * 
     * @return Found histogram or null
     */
    TH1* GetOutputHist(const TList* d, const char* name) const;
    /** 
     * 
     * 
     * 
     * @return 
     */
    Color_t Color() const 
    { 
      return AliForwardUtil::RingColor(fDet, fRing);
    }
    const char* GetName() const { return fName.Data(); } 
    UShort_t fDet;   // Detector
    Char_t   fRing;  // Ring
    TString  fName;  // Name

    ClassDef(RingHistos,1) 
  };
  /* @} */
private:
  AliForwardUtil() {}
  AliForwardUtil(const AliForwardUtil& o) : TObject(o) {}
  AliForwardUtil& operator=(const AliForwardUtil&) { return *this; }
  ~AliForwardUtil() {}
  

  ClassDef(AliForwardUtil,1) // Utilities - do not make object
};

#endif
// Local Variables:
//  mode: C++
// End:

