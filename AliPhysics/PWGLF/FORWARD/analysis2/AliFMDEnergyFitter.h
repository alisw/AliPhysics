//
// Class to fit the energy distribution.  
//
#ifndef ALIFMDENERGYFITTER_H
#define ALIFMDENERGYFITTER_H
/**
 * @file   AliFMDEnergyFitter.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:02:23 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_eloss
 */
#include <TNamed.h>
// #include <TH1D.h>
#include <TAxis.h>
#include <TList.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include "AliFMDCorrELossFit.h"
#include "AliForwardUtil.h"
#include "AliLandauGaus.h"
class TH1;
class TH2;
class AliESDFMD;
class TFitResult;
class TF1;
class TArrayD;

/**
 * Class to fit the energy distribution.  
 *
 * @par Input: 
 *    - AliESDFMD object  - from reconstruction
 *
 * @par Output: 
 *    - Lists of histogram - one per ring.  Each list has a number of 
 *      histograms corresponding to the number of eta bins defined.  
 *
 * @par Corrections used: 
 *    - None
 *
 * @image html alice-int-2012-040-eloss_fits.png "Summary of fits"
 * 
 * @ingroup pwglf_forward_algo
 */
class AliFMDEnergyFitter : public TNamed
{
public: 
  /** 
   * Enumeration of parameters 
   */
  enum { 
    /** Index of pre-constant @f$ C@f$ */
    kC		= AliLandauGaus::kC,
    /** Index of most probable value @f$ \Delta_p@f$ */
    kDelta	= AliLandauGaus::kDelta, 
    /** Index of Landau width @f$ \xi@f$ */
    kXi		= AliLandauGaus::kXi, 
    /** Index of Gaussian width @f$ \sigma@f$ */
    kSigma	= AliLandauGaus::kSigma, 
    /** Index of Gaussian additional width @f$ \sigma_n@f$ */
    kSigmaN	= AliLandauGaus::kSigmaN,
    /** Index of Number of particles @f$ N@f$ */
    kN		= AliLandauGaus::kN, 
    /** Base index of particle strengths @f$ a_i@f$ for 
	@f$i=2,\ldots,N@f$ */
    kA		= AliLandauGaus::kA
  };
  /** 
   * Enumeration of residual methods 
   */
  enum EResidualMethod {
    /** Do not calculate residuals */
    kNoResiduals = 0, 
    /** The residuals stored are the difference, and the errors are
	stored in the error bars of the histogram. */
    kResidualDifference, 
    /** The residuals stored are the differences scaled to the error
	on the data */ 
    kResidualScaledDifference, 
    /** The residuals stored are the square difference scale to the
	square error on the data. */
    kResidualSquareDifference
  };

  /**
   * FMD ring bits for skipping 
   */
   enum FMDRingBits { 
     /** FMD1i */
     kFMD1I=0x01,
     /** All of FMD1 */
     kFMD1 =kFMD1I,
     /** FMD2i */
     kFMD2I=0x02,
     /** FMD2o */
     kFMD2O=0x04,
     /** All of FMD2 */
     kFMD2 =kFMD2I|kFMD2O,
     /** FMD3i */
     kFMD3I=0x08,
     /** FMD3o */
     kFMD3O=0x10,
     /** All of FMD3 */
     kFMD3 =kFMD3I|kFMD3O
  };
  /** 
   * Destructor
   */
  virtual ~AliFMDEnergyFitter();
  /** 
   * Default Constructor - do not use 
   */
  AliFMDEnergyFitter();
  /** 
   * Constructor 
   * 
   * @param title Title of object  - not significant 
   */
  AliFMDEnergyFitter(const char* title);

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Setters of options and parameters 
   */
  /** 
   * Set the eta axis to use.  This will force the code to use this
   * eta axis definition - irrespective of whatever axis is passed to
   * the Init member function.  Therefore, this member function can be
   * used to force another eta axis than one found in the correction
   * objects. 
   * 
   * @param nBins  Number of bins 
   * @param etaMin Minimum of the eta axis 
   * @param etaMax Maximum of the eta axis 
   */
  void SetEtaAxis(Int_t nBins, Double_t etaMin, Double_t etaMax);
  /** 
   * Set the eta axis to use.  This will force the code to use this
   * eta axis definition - irrespective of whatever axis is passed to
   * the Init member function.  Therefore, this member function can be
   * used to force another eta axis than one found in the correction
   * objects. 
   * 
   * @param etaAxis Eta axis to use 
   */
  void SetEtaAxis(const TAxis& etaAxis);
  /** 
   * Set the centrality bins.  E.g., 
   * @code 
   * UShort_t n = 12;
   * Double_t bins[] = {  0.,  5., 10., 15., 20., 30., 
   *                     40., 50., 60., 70., 80., 100. };
   * task->GetFitter().SetCentralityBins(n, bins);
   * @endcode
   * 
   * @param nBins Size of @a bins
   * @param bins  Bin limits. 
   */
  void SetCentralityAxis(UShort_t nBins, Double_t* bins);
  /** 
   * Set the low cut used for energy 
   * 
   * @param lowCut Low cut
   */
  void SetLowCut(Double_t lowCut=0.3) { fLowCut = lowCut; }
  Double_t GetLowCut() const { return fLowCut; }
  /** 
   * Set the number of bins to subtract 
   * 
   * @param n 
   */
  void SetFitRangeBinWidth(UShort_t n=4) { fFitRangeBinWidth = n; }
  /** 
   * Whether or not to enable fitting of the final merged result.  
   * Note, fitting takes quite a while and one should be careful not to do 
   * this needlessly 
   * 
   * @param doFit Whether to do the fits or not 
   */
  void SetDoFits(Bool_t doFit=kTRUE) { fDoFits = doFit; }
  /** 
   * Set whether to make the corrections object on the output.  Note,
   * fits should be enable for this to have any effect.
   * 
   * @param doMake If true (false is default), do make the corrections object. 
   */
  void SetDoMakeObject(Bool_t doMake=kTRUE) { fDoMakeObject = doMake; }
  /** 
   * Set how many particles we will try to fit at most to the data
   * 
   * @param n Max number of particle to try to fit 
   */
  void SetNParticles(UShort_t n) { fNParticles = (n<1 ? 1 : (n>7 ? 7 : n)); }
  /** 
   * Set the minimum number of entries each histogram must have 
   * before we try to fit our response function to it
   * 
   * @param n Minimum number of entries
   */
  void SetMinEntries(UShort_t n) { fMinEntries = (n < 1 ? 1 : n); }
  /**
   * Set maximum energy loss to consider 
   *
   * @param x Maximum energy loss to consider 
   */
  void SetMaxE(Double_t x) { fMaxE = x; }
  /**
   * Set number of energy loss bins 
   *
   * @param x Number of energy loss bins 
   */
  void SetNEbins(Int_t x) { fNEbins = x; }
  /** 
   * Set the maximum relative error 
   * 
   * @param e Maximum relative error 
   */
  void SetMaxRelativeParameterError(Double_t e=0.2) { fMaxRelParError = e; }
  /** 
   * Set the maximum @f$ \chi^2/\nu@f$ 
   * 
   * @param c Maximum @f$ \chi^2/\nu@f$ 
   */
  void SetMaxChi2PerNDF(Double_t c=10) { fMaxChi2PerNDF = c; }
  /** 
   * Set the least weight
   * 
   * @param c Least weight
   */
  void SetMinWeight(Double_t c=1e-7) { fMinWeight = c; }
  /**
   * Set wheter to use increasing bin sizes 
   *
   * @param x Wheter to use increasing bin sizes 
   */
  void SetUseIncreasingBins(Bool_t x) { fUseIncreasingBins = x; }
  /** 
   * Set whether to make residuals, and in that case how. 
   *
   * - Square difference: @f$chi_i^2=(h_i - f(x_i))^2/\delta_i^2@f$ 
   * - Scaled difference: @f$(h_i - f(x_i))/\delta_i@f$ 
   * - Difference: @f$(h_i - f(x_i)) \pm\delta_i@f$ 
   *
   * where @f$h_i, x_i, \delta_i@f$ is the bin content, bin center,
   * and bin error for bin @f$i@f$ respectively, and @f$ f@f$ is the
   * fitted function.
   * 
   * @param x Residual method 
   */
  void SetStoreResiduals(EResidualMethod x=kResidualDifference) 
  { 
    fResidualMethod = x; 
  }
  /** 
   * Set the regularization cut @f$c_{R}@f$.  If a @f$\Delta@f$
   * distribution has more entries @f$ N_{dist}@f$ than @f$c_{R}@f$,
   * then we modify the errors of the the distribution with the factor
   * 
   * @f[
   * \sqrt{N_{dist}/c_{R}}
   * @f]
   *
   * to keep the @f$\chi^2/\nu@f$ within resonable limits. 
   *
   * The large residuals @f$chi_i^2=(h_i - f(x_i))^2/\delta_i^2@f$
   * (see also SetStoreResiduals) comes about on the boundary between
   * the @f$N@f$ and @f$N+1@f$ particle contributions, and seems to
   * fall off for larger @f$N@f$. This may indicate that there's a
   * component in the distributions that the function
   *
   * @f[
   *   f(\Delta;\Delta_p,\xi,\sigma,\mathbf{a}) = \sum_i=1^{n} a_i\int
   *   d\Delta' L(\Delta;\Delta',\xi) G(\Delta';\Delta_p,\sigma)
   * @f]
   * 
   * does not capture.   
   *
   * @param cut
   */
  void SetRegularizationCut(Double_t cut=3e6) 
  {
    fRegularizationCut = cut;
  }
  void SetSkips(UShort_t skip) { fSkips = skip; }
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1);
  /**
   * Whether to enable the extra shift in the MPV from @f$ \sigma/\xi@f$ 
   *
   * @param use If true, enable extra shift @f$\delta\Delta_p(\sigma/\xi)@f$  
   */
  void SetEnableDeltaShift(Bool_t use=true);

  /* @} */
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Processing 
   */
  void Init();
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  virtual void CreateOutputObjects(TList* dir);
  /** 
   * Initialise the task
   * 
   * @param etaAxis The eta axis to use.  Note, that if the eta axis
   * has already been set (using SetEtaAxis), then this parameter will be 
   * ignored
   *
   * @param sys Collision system identifier.  If set, then it will
   * define how many MIP peaks we will fit.  For pp we do 3, for pA,
   * Ap, and AA we do 5.
   */
  virtual void SetupForData(const TAxis& etaAxis, UShort_t sys=0);
  /** 
   * Fitter the input AliESDFMD object
   * 
   * @param input     Input 
   * @param cent      Event centrality (or < 0 if not valid)
   * @param empty     Whether the event is 'empty'
   * 
   * @return True on success, false otherwise 
   */
  virtual Bool_t Accumulate(const AliESDFMD& input, 
			    Double_t         cent,
			    Bool_t           empty);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir Where the histograms are  
   */
  virtual void Fit(const TList* dir);
  /** 
   * Generate the corrections object 
   * 
   * @param dir List to analyse 
   */
  void MakeCorrectionsObject(TList* dir);
  /** @} */
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
  /** 
   * Read the parameters from a list - used when re-running the code 
   * 
   * @param list Input list 
   * 
   * @return true if the parameter where read 
   */
  Bool_t ReadParameters(const TCollection* list);
  //==================================================================
  /** 
   * Internal data structure to keep track of the histograms. Objects
   * of this class are streamed, as we create the objects at task
   * initialization time.  We must therefore create a dictionary for
   * this class.  This does _not_ mean that one should actually use
   * this class outside of this (or a derived) class.
   */
  struct RingHistos : public AliForwardUtil::RingHistos
  { 
    typedef AliFMDCorrELossFit::ELossFit ELossFit_t;
    /** 
     * Default CTOR
     */
    RingHistos();
    /** 
     * Constructor
     * 
     * @param d detector
     * @param r ring 
     */
    RingHistos(UShort_t d, Char_t r);
    /** 
     * Copy constructor - not defined
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o){;}
    /** 
     * Assignment operator  - not defined
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this 
     */
    RingHistos& operator=(const RingHistos& o){return *this;}
    /** 
     * Destructor 
     */
    ~RingHistos();
    /** 
     * Make an axis with increasing bins 
     * 
     * @param n    Number of bins 
     * @param min  Minimum 
     * @param max  Maximum
     * 
     * @return An axis with quadratically increasing bin size 
     */
    virtual TArrayD MakeIncreasingAxis(Int_t    n, 
				       Double_t min, 
				       Double_t max) const;
    /** 
     * Make E/E_mip histogram 
     * 
     * @param name    Name of histogram
     * @param title   Title of histogram
     * @param eAxis   @f$\eta@f$ axis
     * @param deMax   Maximum energy loss 
     * @param nDeBins Number energy loss bins 
     * @param incr    Whether to make bins of increasing size
     */
    TH2* Make(const char*  name, 
	      const char*  title, 
	      const TAxis& eAxis, 
	      Double_t     deMax=12, 
	      Int_t        nDeBins=300, 
	      Bool_t       incr=true);
    /** 
     * Define outputs
     * 
     * @param dir 
     */
    virtual void CreateOutputObjects(TList* dir);
    /** 
     * Initialise object 
     * 
     * @param eAxis      Eta axis
     * @param cAxis      Centrality axis 
     * @param maxDE      Max energy loss to consider 
     * @param nDEbins    Number of bins 
     * @param useIncrBin Whether to use an increasing bin size 
     */
    virtual void SetupForData(const TAxis& eAxis, 
			      const TAxis& cAxis,
			      Double_t     maxDE=10, 
			      Int_t        nDEbins=300, 
			      Bool_t       useIncrBin=true);
    /** 
     * Fill histogram 
     * 
     * @param empty  True if event is empty
     * @param eta    @f$ Eta@f$
     * @param icent  Centrality bin (1 based)
     * @param mult   Signal 
     */
    virtual void Fill(Bool_t empty, Double_t eta, Int_t icent, Double_t mult);
    /** 
     * Get the the 2D histogram eloss name from our sub-list of @a dir
     * and call the Fit function described below (with &fBest) as last
     * argument.
     * 
     * @param dir         Output list 
     * @param lowCut      Lower cut 
     * @param nParticles  Max number of convolved landaus to fit
     * @param minEntries  Minimum number of entries 
     * @param minusBins   Number of bins from peak to subtract to 
     *                    get the fit range 
     * @param relErrorCut Cut applied to relative error of parameter. 
     *                    Note, for multi-particle weights, the cut 
     *                    is loosend by a factor of 2 
     * @param chi2nuCut   Cut on @f$ \chi^2/\nu@f$ - 
     *                    the reduced @f$\chi^2@f$ 
     * @param minWeight   Least weight ot consider
     * @param regCut      Regularization cut-off
     * @param residuals   Mode for residual plots
     *
     * @return List of fit parameters 
     */
    virtual TObjArray* Fit(TList*          dir, 
			   Double_t        lowCut, 
			   UShort_t        nParticles,
			   UShort_t        minEntries,
			   UShort_t        minusBins,
			   Double_t        relErrorCut, 
			   Double_t        chi2nuCut,
			   Double_t        minWeight,
			   Double_t        regCut,
			   EResidualMethod residuals) const;
    /** 
     * Get the the 2D histogram @a name from our sub-list of @a
     * dir. Then for each eta slice, try to fit the energu loss
     * distribution up to @a nParticles particle responses.
     *
     * The fitted distributions (along with the functions fitted) are
     * stored in a newly created sublist (<i>name</i>Dists).
     *
     * The fit parameters are also recorded in the newly created sub-list 
     * <i>name</i>Results.  
     *
     * If @a residuals is not equal to kNoResiduals, then the
     * residuals of the fits will be stored in the newly created
     * sub-list <i>name</i>Residuals.
     *
     * A histogram named <i>name</i>Status is also generated and
     * stored in the output list.
     * 
     * @param dir         Output list 
     * @param name        Name of 2D base histogram in list
     * @param lowCut      Lower cut 
     * @param nParticles  Max number of convolved landaus to fit
     * @param minEntries  Minimum number of entries 
     * @param minusBins   Number of bins from peak to subtract to 
     *                    get the fit range 
     * @param relErrorCut Cut applied to relative error of parameter. 
     *                    Note, for multi-particle weights, the cut 
     *                    is loosend by a factor of 2 
     * @param chi2nuCut   Cut on @f$ \chi^2/\nu@f$ - 
     *                    the reduced @f$\chi^2@f$ 
     * @param minWeight   Least weight ot consider
     * @param regCut      Regularization cut-off
     * @param residuals   Mode for residual plots
     * @param scaleToPeak If true, scale distribution to peak value
     * @param best        Optional array to store fits in
     *
     * @return List of fit parameters 
     */
    virtual TObjArray* FitSlices(TList*          dir, 
				 const char*     name,
				 Double_t        lowCut, 
				 UShort_t        nParticles,
				 UShort_t        minEntries,
				 UShort_t        minusBins,
				 Double_t        relErrorCut, 
				 Double_t        chi2nuCut,
				 Double_t        minWeight,
				 Double_t        regCut,
				 EResidualMethod residuals,
				 Bool_t          scaleToPeak=true,
				 TObjArray*      best=0) const;
    /** 
     * Do scaling of histogram before fitting.  This can be
     * overwritten to do some smoothing or the like. By default, this
     * simply scales to the bin width.
     * 
     * @param dist Histogram to scale. 
     */     
    virtual void Scale(TH1* dist) const;
    /** 
     * Fit a signal histogram.  First, the bin @f$ b_{min}@f$ with
     * maximum bin content in the range @f$ [E_{min},\infty]@f$ is
     * found.  Then the fit range is set to the bin range 
     * @f$ [b_{min}-\Delta b,b_{min}+2\Delta b]@f$, and a 1 
     * particle signal is fitted to that.  The parameters of that fit 
     * is then used as seeds for a fit of the @f$ N@f$ particle response 
     * to the data in the range 
     * @f$ [b_{min}-\Delta b,N(\Delta_1+\xi_1\log(N))+2N\xi@f$
     * 
     * @param dist        Histogram to fit 
     * @param lowCut      Lower cut @f$ E_{min}@f$ on signal 
     * @param minEntries  Least number of entries required
     * @param nParticles  Max number @f$ N@f$ of convolved landaus to fit
     * @param minusBins   Number of bins @f$ \Delta b@f$ from peak to 
     *                    subtract to get the fit range 
     * @param relErrorCut Cut applied to relative error of parameter. 
     *                    Note, for multi-particle weights, the cut 
     *                    is loosend by a factor of 2 
     * @param chi2nuCut   Cut on @f$ \chi^2/\nu@f$ - 
     *                    the reduced @f$\chi^2@f$ 
     * @param minWeight   Least weight ot consider
     * @param regCut      Regularization cut-off
     * @param scaleToPeak If true, scale distribution to peak value
     * @param status      On return, contain the status code (0: OK, 1:
     *                    empty, 2: low statistics, 3: fit failed)
     * 
     * @return The best fit function 
     */
    virtual ELossFit_t* FitHist(TH1*      dist,
				Double_t  lowCut, 
				UShort_t  nParticles,
				UShort_t  minEntries,
				UShort_t  minusBins,
				Double_t  relErrorCut, 
				Double_t  chi2nuCut,
				Double_t  minWeight,
				Double_t  regCut,
				Bool_t    scaleToPeak,
				UShort_t& status) const;
    /** 
     * Find the best fit 
     * 
     * @param dist           Histogram 
     * @param relErrorCut    Cut applied to relative error of parameter. 
     *                       Note, for multi-particle weights, the cut 
     *                       is loosend by a factor of 2 
     * @param chi2nuCut      Cut on @f$ \chi^2/\nu@f$ - 
     *                       the reduced @f$\chi^2@f$ 
     * @param minWeightCut   Least valid @f$ a_i@f$ 
     * 
     * @return Best fit 
     */
    virtual ELossFit_t* FindBestFit(const TH1* dist,
				    Double_t   relErrorCut, 
				    Double_t   chi2nuCut,
				    Double_t   minWeightCut) const;
    /** 
     * Calculate residuals of the fit 
     * 
     * @param mode   How to calculate 
     * @param lowCut Lower cut 
     * @param dist   Distribution 
     * @param fit    Function fitted to distribution
     * @param out    Output list to store residual histogram in
     */
    virtual void CalculateResiduals(EResidualMethod  mode,
				    Double_t         lowCut,
				    TH1*             dist, 
				    ELossFit_t*      fit, 
				    TCollection*     out) const;
    /** 
     * Find the best fits.  This assumes that the array fBest has been
     * filled with the best possible fits for each eta bin, and that
     * the fits are placed according to the bin number of the eta bin.
     *
     * This is called by the parent class when generating the corretion 
     * object. 
     * 
     * @param d    Parent list
     * @param obj  Object to add fits to
     * @param eta  Eta axis 
     */
    virtual void FindBestFits(const TList*        d, 
			      AliFMDCorrELossFit& obj,
			      const TAxis&        eta);
    /** 
     * Make a parameter histogram
     * 
     * @param name   Name of histogram.
     * @param title  Title of histogram. 
     * @param eta    Eta axis 
     * 
     * @return 
     */
    TH1* MakePar(const char* name, const char* title, const TAxis& eta) const;
    /** 
     * Make a histogram that contains the results of the fit over the
     * full ring
     * 
     * @param name  Name 
     * @param title Title
     * @param eta   Eta axis 
     * @param low   Least bin
     * @param high  Largest bin
     * @param val   Value of parameter 
     * @param err   Error on parameter 
     * 
     * @return The newly allocated histogram 
     */
    TH1* MakeTotal(const char*  name, 
		    const char*  title, 
		    const TAxis& eta, 
		    Int_t        low, 
		    Int_t        high, 
		    Double_t     val, 
		    Double_t     err) const;
    TH1*                 fEDist;     // Ring energy distribution 
    TH1*                 fEmpty;     // Ring energy dist for empty events
    TH2*                 fHist;      // Two dimension Delta distribution
    // TList*               fEtaEDists; // Energy distributions per eta bin. 
    TList*               fList;
    mutable TObjArray    fBest;
    mutable TClonesArray fFits;
    Int_t                fDebug;
    ClassDef(RingHistos,4);
  };
protected:
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDEnergyFitter(const AliFMDEnergyFitter& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliFMDEnergyFitter& operator=(const AliFMDEnergyFitter& o);

  virtual RingHistos* CreateRingHistos(UShort_t d, Char_t r) const;
  /** 
   * Get the ring histogram container 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Ring histogram container 
   */
  RingHistos* GetRingHistos(UShort_t d, Char_t r) const;
  /** 
   * Check if the detector @a d, ring @a r is listed <i>in</i> the @a
   * skips bit mask.  If the detector/ring is in the mask, return true.
   * 
   * That is, use case is 
   * @code 
   *  for (UShort_t d=1. d<=3, d++) {
   *    UShort_t nr = (d == 1 ? 1 : 2);
   *    for (UShort_t q = 0; q < nr; q++) { 
   *      Char_t r = (q == 0 ? 'I' : 'O');
   *      if (CheckSkips(d, r, skips)) continue; 
   *      // Process detector/ring 
   *    }
   *  }
   * @endcode
   *
   * @param d      Detector
   * @param r      Ring 
   * @param skips  Mask of detector/rings to skip
   * 
   * @return True if detector @a d, ring @a r is in the mask @a skips 
   */
  static Bool_t CheckSkip(UShort_t d, Char_t r, UShort_t skips);

  TList           fRingHistos;        // List of histogram containers
  Double_t        fLowCut;            // Low cut on energy
  UShort_t        fNParticles;        // Number of landaus to try to fit 
  UShort_t        fMinEntries;        // Minimum number of entries
  UShort_t        fFitRangeBinWidth;  // N-bins to subtract from found max
  Bool_t          fDoFits;            // Whether to actually do the fits 
  Bool_t          fDoMakeObject;      // Whether to make corrections object
  TAxis           fEtaAxis;           // Eta axis 
  TAxis           fCentralityAxis;    // Centrality axis 
  Double_t        fMaxE;              // Maximum energy loss to consider 
  Int_t           fNEbins;            // Number of energy loss bins 
  Bool_t          fUseIncreasingBins; // Wheter to use increasing bin sizes 
  Double_t        fMaxRelParError;    // Relative error cut
  Double_t        fMaxChi2PerNDF;     // chi^2/nu cit
  Double_t        fMinWeight;         // Minimum weight value 
  Int_t           fDebug;             // Debug level 
  EResidualMethod fResidualMethod;    // Whether to store residuals (debugging)
  UShort_t        fSkips;             // Rings to skip when fitting 
  Double_t        fRegularizationCut; // When to regularize the chi^2

  ClassDef(AliFMDEnergyFitter,8); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
