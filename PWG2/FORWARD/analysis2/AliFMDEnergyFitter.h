#ifndef ALIROOT_PWG2_FORWARD_ALIFMDENERGYFITTER_H
#define ALIROOT_PWG2_FORWARD_ALIFMDENERGYFITTER_H
#include <TNamed.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TList.h>
#include <TObjArray.h>
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
 * @ingroup pwg2_forward_analysis 
 */
class AliFMDEnergyFitter : public TNamed
{
public: 
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

  void Init(const TAxis& etaAxis);
  /** 
   * Set the low cut used for energy 
   * 
   * @param lowCut Low cut
   */
  void SetLowCut(Double_t lowCut=0.3) { fLowCut = lowCut; }
  void SetBinsToSubtract(UShort_t n=4) { fBinsToSubtract = n; }
  /** 
   * Whether or not to enable fitting of the final merged result.  
   * Note, fitting takes quite a while and one should be careful not to do 
   * this needlessly 
   * 
   * @param doFit Whether to do the fits or not 
   */
  void DoFits(Bool_t doFit=kTRUE) { fDoFits = doFit; }
  /** 
   * Whether or not to enable fitting of the final merged result.  
   * Note, fitting takes quite a while and one should be careful not to do 
   * this needlessly 
   * 
   * @param doFit Whether to do the fits or not 
   */
  void SetNLandau(UShort_t n) { fNLandau = (n < 1 ? 1 : (n > 5 ? 5 : n)); }
  /** 
   * Whether or not to enable fitting of the final merged result.  
   * Note, fitting takes quite a while and one should be careful not to do 
   * this needlessly 
   * 
   * @param doFit Whether to do the fits or not 
   */
  void SetMinEntries(UShort_t n) { fMinEntries = (n < 1 ? 1 : n); }
  /** 
   * Fitter the input AliESDFMD object
   * 
   * @param input     Input 
   * @param empty     Whether the event is 'empty'
   * 
   * @return True on success, false otherwise 
   */
  Bool_t Accumulate(const AliESDFMD& input, 
		    Bool_t           empty);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param nEvents Number of events 
   */
  void Fit(TList* dir);
  
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
     * @param eAxis 
     */
    void Init(const TAxis& eAxis, 
	      Double_t maxDE=10, 
	      Int_t nDEbins=300, 
	      Bool_t useIncrBin=true);
    /** 
     * Fill histogram 
     * 
     * @param empty  True if event is empty
     * @param ieta   Eta bin
     * @param mult   Signal 
     */
    void Fill(Bool_t empty, Int_t ieta, Double_t mult);
    /** 
     * Scale the histograms to the total number of events 
     * 
     * @param dir     Output list 
     * @param eta     Eta axis 
     * @param lowCut  Lower cut 
     * @param nLandau Max number of convolved landaus to fit
     */
    TObjArray* Fit(TList* dir, const TAxis& eta,
		   Double_t lowCut, UShort_t nLandau,
		   UShort_t minEntries,
		   UShort_t minusBins) const;
    /** 
     * Fit a signal histogram 
     * 
     * @param dist     Historgam to fit 
     * @param lowCut   Lower cut on signal 
     * @param nLandau  Max number of convolved landaus to fit
     * 
     * @return The best fit function 
     */
    TF1* FitHist(TH1* dist,Double_t lowCut, UShort_t nLandau,
		 UShort_t minEntries, UShort_t minusBins) const;
    /** 
     * Fit more Landau 
     * 
     * @param dist     Histogram to fit 
     * @param nLandau  Number of landaus 
     * @param r        Result from 1st Landau fit
     * @param landau1  Function of 1st Landau fit
     * @param minE     Least signal for range 
     * 
     * @return The result 
     */
    TF1* FitMore(TH1*        dist,
		 UShort_t    nLandau,
		 TFitResult& r, 
		 TF1*        landau1,
		 Double_t    minE,
		 Double_t    maxRange) const;
    /** 
     * Fit more Landau 
     * 
     * @param dist     Histogram to fit 
     * @param nLandau  Number of landaus 
     * @param r        Result from 1st Landau fit
     * @param landau1  Function of 1st Landau fit
     * @param minE     Least signal for range 
     * 
     * @return The result 
     */
    TF1* FitMore2(TH1*        dist,
		  UShort_t    nLandau,
		  TFitResult& r, 
		  TF1*        landau1,
		  Double_t    minE,
		  Double_t    maxRange) const;
    /** 
     * Check the result of the fit. Returns true if the reduced 
     * @f$ \chi^2/\nu@f$ is less than 5, and that the relative error 
     * @f$ \Delta p_i/p_i@f$ on each parameter is less than 20 percent. 
     * 
     * @param r Result to check
     * 
     * @return true if fit is good. 
     */
    Bool_t CheckResult(TFitResult& r) const;
    TArrayD MakeIncreasingAxis(Int_t n, Double_t min, Double_t max) const;
    /** 
     * Make E/E_mip histogram 
     * 
     * @param ieta   Eta bin
     * @param eMin   Least signal
     * @param eMax   Largest signal 
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
    TH1D* MakeTotal(const char* name, 
		    const char* title, 
		    const TAxis& eta, 
		    Int_t low, 
		    Int_t high, 
		    Double_t val, 
		    Double_t err) const;
    TH1D*     fEDist;        // Ring energy distribution 
    TH1D*     fEmpty;        // Ring energy distribution for empty events
    TList     fEtaEDists;    // Energy distributions per eta bin. 
    TList*    fList;
    Int_t     fDebug;
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
  UShort_t fNLandau;       // Number of landaus to try to fit 
  UShort_t fMinEntries;    // Minimum number of entries
  UShort_t fBinsToSubtract;// Number of bins to subtract from found max
  Bool_t   fDoFits;        // Wheter to actually do the fits 
  TAxis    fEtaAxis;       // Eta axis 
  Int_t    fDebug;         // Debug level 

  ClassDef(AliFMDEnergyFitter,1); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
