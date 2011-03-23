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
 * @ingroup pwg2_forward_eloss
 */
#include <TNamed.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TList.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include "AliFMDCorrELossFit.h"
#include "AliForwardUtil.h"
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
 *
 * @ingroup pwg2_forward_algo
 * @ingroup pwg2_forward_eloss
 */
class AliFMDEnergyFitter : public TNamed
{
public: 
    enum { 
      kC	= AliForwardUtil::ELossFitter::kC,
      kDelta	= AliForwardUtil::ELossFitter::kDelta, 
      kXi	= AliForwardUtil::ELossFitter::kXi, 
      kSigma	= AliForwardUtil::ELossFitter::kSigma, 
      kSigmaN	= AliForwardUtil::ELossFitter::kSigmaN,
      kN	= AliForwardUtil::ELossFitter::kN, 
      kA	= AliForwardUtil::ELossFitter::kA
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

  /** 
   * Initialise the task
   * 
   * @param etaAxis The eta axis to use.  Note, that if the eta axis
   * has already been set (using SetEtaAxis), then this parameter will be 
   * ignored
   */
  void Init(const TAxis& etaAxis);
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
  void SetNParticles(UShort_t n) { fNParticles = (n < 1 ? 1 : (n > 5 ? 5 : n)); }
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
   * Fitter the input AliESDFMD object
   * 
   * @param input     Input 
   * @param cent      Event centrality (or < 0 if not valid)
   * @param empty     Whether the event is 'empty'
   * 
   * @return True on success, false otherwise 
   */
  Bool_t Accumulate(const AliESDFMD& input, 
		    Double_t         cent,
		    Bool_t           empty);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir Where the histograms are  
   */
  void Fit(const TList* dir);
  /** 
   * Generate the corrections object 
   * 
   * @param dir List to analyse 
   */
  void MakeCorrectionsObject(TList* dir);
  
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  void DefineOutput(TList* dir);
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1);
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
protected:
  /** 
   * Internal data structure to keep track of the histograms
   */
  struct RingHistos : public AliForwardUtil::RingHistos
  { 
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
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o);
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this 
     */
    RingHistos& operator=(const RingHistos& o);
    /** 
     * Destructor 
     */
    ~RingHistos();
    /** 
     * Define outputs
     * 
     * @param dir 
     */
    void Output(TList* dir);
    /** 
     * Initialise object 
     * 
     * @param eAxis      Eta axis
     * @param cAxis      Centrality axis 
     * @param maxDE      Max energy loss to consider 
     * @param nDEbins    Number of bins 
     * @param useIncrBin Whether to use an increasing bin size 
     */
    void Init(const TAxis& eAxis, 
	      const TAxis& cAxis,
	      Double_t     maxDE=10, 
	      Int_t        nDEbins=300, 
	      Bool_t       useIncrBin=true);
    /** 
     * Fill histogram 
     * 
     * @param empty  True if event is empty
     * @param ieta   Eta bin (0 based)
     * @param icent  Centrality bin (1 based)
     * @param mult   Signal 
     */
    void Fill(Bool_t empty, Int_t ieta, Int_t icent, Double_t mult);
    /** 
     * Fit each histogram to up to @a nParticles particle responses.
     * 
     * @param dir         Output list 
     * @param eta         Eta axis 
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
     */
    TObjArray* Fit(TList* dir, 
		   const TAxis& eta,
		   Double_t     lowCut, 
		   UShort_t     nParticles,
		   UShort_t     minEntries,
		   UShort_t     minusBins,
		   Double_t     relErrorCut, 
		   Double_t     chi2nuCut) const;
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
     * @param nParticles  Max number @f$ N@f$ of convolved landaus to fit
     * @param minusBins   Number of bins @f$ \Delta b@f$ from peak to 
     *                    subtract to get the fit range 
     * @param relErrorCut Cut applied to relative error of parameter. 
     *                    Note, for multi-particle weights, the cut 
     *                    is loosend by a factor of 2 
     * @param chi2nuCut   Cut on @f$ \chi^2/\nu@f$ - 
     *                    the reduced @f$\chi^2@f$ 
     * 
     * @return The best fit function 
     */
    TF1* FitHist(TH1*     dist,
		 Double_t lowCut, 
		 UShort_t nParticles,
		 UShort_t minusBins,
		 Double_t relErrorCut, 
		 Double_t chi2nuCut) const;
    /** 
     * Find the best fits 
     * 
     * @param d              Parent list
     * @param obj            Object to add fits to
     * @param eta            Eta axis 
     * @param relErrorCut    Cut applied to relative error of parameter. 
     *                       Note, for multi-particle weights, the cut 
     *                       is loosend by a factor of 2 
     * @param chi2nuCut      Cut on @f$ \chi^2/\nu@f$ - 
     *                       the reduced @f$\chi^2@f$ 
     * @param minWeightCut   Least valid @f$ a_i@f$ 
     */
    void FindBestFits(const TList*        d, 
		      AliFMDCorrELossFit& obj,
		      const TAxis&        eta,     
		      Double_t            relErrorCut, 
		      Double_t            chi2nuCut,
		      Double_t            minWeightCut);
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
    AliFMDCorrELossFit::ELossFit* FindBestFit(const TH1* dist,
					      Double_t relErrorCut, 
					      Double_t chi2nuCut,
					      Double_t minWeightCut);
    /** 
     * Check the result of the fit. Returns true if 
     * - @f$ \chi^2/\nu < \max{\chi^2/\nu}@f$
     * - @f$ \Delta p_i/p_i < \delta_e@f$ for all parameters.  Note, 
     *   for multi-particle fits, this requirement is relaxed by a 
     *   factor of 2
     * - @f$ a_{n} > 10^{-7}@f$ when fitting to an @f$ n@f$ 
     *   particle response 
     * 
     * @param r           Result to check
     * @param relErrorCut Cut @f$ \delta_e@f$ applied to relative error 
     *                    of parameter.  
     * @param chi2nuCut   Cut @f$ \max{\chi^2/\nu}@f$ 
     * 
     * @return true if fit is good. 
     */
    Bool_t CheckResult(TFitResult* r,
		       Double_t    relErrorCut, 
		       Double_t    chi2nuCut) const;
    /** 
     * Make an axis with increasing bins 
     * 
     * @param n    Number of bins 
     * @param min  Minimum 
     * @param max  Maximum
     * 
     * @return An axis with quadratically increasing bin size 
     */
    TArrayD MakeIncreasingAxis(Int_t n, Double_t min, Double_t max) const;
    /** 
     * Make E/E_mip histogram 
     * 
     * @param ieta    Eta bin
     * @param eMin    Least signal
     * @param eMax    Largest signal 
     * @param deMax   Maximum energy loss 
     * @param nDeBins Number energy loss bins 
     * @param incr    Whether to make bins of increasing size
     */
    void Make(Int_t ieta, Double_t eMin, Double_t eMax, 
	      Double_t deMax=12, Int_t nDeBins=300, Bool_t incr=true);
    /** 
     * Make a parameter histogram
     * 
     * @param name   Name of histogram.
     * @param title  Title of histogram. 
     * @param eta    Eta axis 
     * 
     * @return 
     */
    TH1D* MakePar(const char* name, const char* title, const TAxis& eta) const;
    /** 
     * Make a histogram that contains the results of the fit over the full ring 
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
    TH1D* MakeTotal(const char* name, 
		    const char* title, 
		    const TAxis& eta, 
		    Int_t low, 
		    Int_t high, 
		    Double_t val, 
		    Double_t err) const;
    TH1D*        fEDist;        // Ring energy distribution 
    TH1D*        fEmpty;        // Ring energy distribution for empty events
    TList        fEtaEDists;    // Energy distributions per eta bin. 
    TList*       fList;
    TClonesArray fFits;
    Int_t        fDebug;
    ClassDef(RingHistos,1);
  };
  /** 
   * Get the ring histogram container 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Ring histogram container 
   */
  RingHistos* GetRingHistos(UShort_t d, Char_t r) const;

  TList    fRingHistos;    // List of histogram containers
  Double_t fLowCut;        // Low cut on energy
  UShort_t fNParticles;    // Number of landaus to try to fit 
  UShort_t fMinEntries;    // Minimum number of entries
  UShort_t fFitRangeBinWidth;// Number of bins to subtract from found max
  Bool_t   fDoFits;        // Whether to actually do the fits 
  Bool_t   fDoMakeObject;  // Whether to make corrections object
  TAxis    fEtaAxis;       // Eta axis 
  TAxis    fCentralityAxis;// Centrality axis 
  Double_t fMaxE;          // Maximum energy loss to consider 
  Int_t    fNEbins;        // Number of energy loss bins 
  Bool_t   fUseIncreasingBins; // Wheter to use increasing bin sizes 
  Double_t fMaxRelParError;// Relative error cut
  Double_t fMaxChi2PerNDF; // chi^2/nu cit
  Double_t fMinWeight;     // Minimum weight value 
  Int_t    fDebug;         // Debug level 
  

  ClassDef(AliFMDEnergyFitter,2); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
