// Object holding the Energy loss fit 'correction'
// 
// These are generated from Monte-Carlo or real ESDs. 
#ifndef ALIFMDCORRELOSSFIT_H
#define ALIFMDCORRELOSSFIT_H
/**
 * @file   AliFMDCorrELossFit.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:01:15 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_eloss
 * 
 */
#include <TObject.h>
#include <TAxis.h>
#include <TObjArray.h>
class TF1;
class TBrowser;

/** 
 * @defgroup pwg2_forward_corr Corrections 
 * 
 * @ingroup pwg2_forward
 */
/** 
 * Object holding the Energy loss fit 'correction'
 * 
 * These are generated from Monte-Carlo or real ESDs. 
 *
 * @ingroup pwg2_forward_corr
 * @ingroup pwg2_forward_eloss
 */
class AliFMDCorrELossFit : public TObject 
{
public:
  /** 
   * POD structure to hold data from fits 
   * 
   * @ingroup pwg2_forward_corr
   */
  struct ELossFit : public TObject 
  {
    Int_t     fN;      // Number of peaks fitted
    UShort_t  fNu;     // Number degrees of freedom
    Double_t  fChi2;   // Chi square from fit
    Double_t  fC;      // Scaling constant 
    Double_t  fDelta;  // Most probable value 
    Double_t  fXi;     // Width parameter of Landau 
    Double_t  fSigma;  // Sigma on folded gaussian 
    Double_t  fSigmaN; // Sigma of detector noise 
    Double_t* fA;      // [fN] Weights 
    Double_t  fEC;     // Error on C 
    Double_t  fEDelta; // Error on Delta 
    Double_t  fEXi;    // Error on Xi
    Double_t  fESigma; // Error on sigma 
    Double_t  fESigmaN;// Error on sigma (noise)
    Double_t* fEA;     // [fN] Error on weights
    Int_t     fQuality;// Assigned quality 
    UShort_t  fDet;    // Detector 
    Char_t    fRing;   // Ring
    UShort_t  fBin;    // Eta bin

    static Double_t fgMaxRelError;  // Global default max relative error
    static Double_t fgLeastWeight;  // Global default least weight 
    static Double_t fgMaxChi2nu;    // Global default maximum reduced chi^2
    /**
     * Default constructor 
     * 
     */
    ELossFit();
    /** 
     * Construct from a function
     * 
     * @param quality Quality flag
     * @param f       Function
     */
    ELossFit(Int_t quality,const TF1& f);
    /** 
     * Constructor with full parameter set
     * 
     * @param quality   Quality flag
     * @param n         @f$ N@f$ - Number of fitted peaks
     * @param chi2      @f$ \chi^2 @f$
     * @param nu        @f$ \nu @f$ - number degrees of freedom
     * @param c         @f$ C@f$ - scale constant
     * @param ec        @f$ \delta C@f$ - error on @f$ C@f$ 
     * @param delta     @f$ \Delta@f$ - Most probable value		  
     * @param edelta    @f$ \delta\Delta@f$ - error on @f$\Delta@f$ 
     * @param xi        @f$ \xi@f$ - width  
     * @param exi       @f$ \delta\xi@f$ - error on @f$\xi@f$ 
     * @param sigma     @f$ \sigma@f$ - Width of Gaussian		   
     * @param esigma    @f$ \delta\sigma@f$ - error on @f$\sigma@f$ 
     * @param sigman    @f$ \sigma_n@f$ - Noise width		  
     * @param esigman   @f$ \delta\sigma_n@f$ - error on @f$\sigma_n@f$ 
     * @param a         Array of @f$ N-1@f$ weights @f$ a_i@f$ for 
     *                  @f$ i=2,\ldots@f$ 
     * @param ea        Array of @f$ N-1@f$ error on the weights @f$ a_i@f$ for 
     *                  @f$ i=2,\ldots@f$ 
     */
    ELossFit(Int_t     quality,UShort_t  n, 
	     Double_t  chi2,   UShort_t  nu, 
	     Double_t  c,      Double_t  ec, 
	     Double_t  delta,  Double_t  edelta, 
	     Double_t  xi,     Double_t  exi,
	     Double_t  sigma,  Double_t  esigma, 
	     Double_t  sigman, Double_t  esigman, 
	     const Double_t* a,const Double_t* ea);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    ELossFit(const ELossFit& o);
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object 
     */
    ELossFit& operator=(const ELossFit& o);
    /** 
     * Destructor 
     */
    ~ELossFit();
    /**
     * @{
     * @name Access to parameters 
     */
    /**
     * @return Number of peaks fitted
     */
    Int_t GetN() const { return fN; }
    /**
     * @return Number degrees of freedom
     */
    UShort_t GetNu() const { return fNu; }
    /**
     * @return Chi square from fit
     */
    Double_t GetChi2() const { return fChi2; }
    /**
     * @return Scaling constant 
     */
    Double_t GetC() const { return fC; }
    /**
     * @return Most probable value 
     */
    Double_t GetDelta() const { return fDelta; }
    /**
     * @return Width parameter of Landau 
     */
    Double_t GetXi() const { return fXi; }
    /**
     * @return Sigma on folded gaussian 
     */
    Double_t GetSigma() const { return fSigma; }
    /**
     * @return Sigma of detector noise 
     */
    Double_t GetSigmaN() const { return fSigmaN; }
    /**
     * @return Weights 
     */
    Double_t* GetAs() const { return fA; }
    /**
     * @return Weights 
     */
    Double_t GetA(UShort_t i) const;    
    /**
     * @return Error on C 
     */
    Double_t GetEC() const { return fEC; }
    /**
     * @return Error on Delta 
     */
    Double_t GetEDelta() const { return fEDelta; }
    /**
     * @return Error on Xi
     */
    Double_t GetEXi() const { return fEXi; }
    /**
     * @return Error on sigma 
     */
    Double_t GetESigma() const { return fESigma; }
    /**
     * @return Error on sigma (noise)
     */
    Double_t GetESigmaN() const { return fESigmaN; }
    /**
     * @return Error on weights
     */
    Double_t* GetEAs() const { return fEA; }
    /**
     * @return Error on weights
     */
    Double_t GetEA(UShort_t i) const;
    /**
     * @return Assigned quality 
     */
    Int_t GetQuality() const { return fQuality; }
    /**
     * @return Detector 
     */
    UShort_t GetDet() const { return fDet; }
    /**
     * @return Ring
     */
    Char_t GetRing() const { return fRing; }
    /**
     * @return Eta bin
     */
    UShort_t GetBin() const { return fBin; }
    /* @} */

    /** 
     * @{ 
     * @name Evaluation 
     */
    /** 
     * Evaluate 
     * @f[ 
     *  f_N(x;\Delta,\xi,\sigma') = 
     *     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i')
     * @f] 
     *
     * (see AliForwardUtil::NLandauGaus) for the maximum @f$ N @f$
     * that fulfills the requirements 
     *
     * @param x           Where to evaluate 
     * @param maxN 	  @f$ \max{N}@f$    
     * 
     * @return @f$ f_N(x;\Delta,\xi,\sigma')@f$ 
     */
    Double_t Evaluate(Double_t x, 
		      UShort_t maxN=999) const;
    /** 
     * Evaluate 
     * @f[ 
     *   f_W(x;\Delta,\xi,\sigma') = 
     *   \frac{\sum_{i=1}^{n} i a_i f_i(x;\Delta,\xi,\sigma')}{
     *     f_N(x;\Delta,\xi,\sigma')} = 
     *   \frac{\sum_{i=1}^{n} i a_i f(x;\Delta_i,\xi_i,\sigma_i')}{
     *     \sum_{i=1}^{n} a_i f(x;\Delta_i,\xi_i,\sigma_i')}
     * @f] 
     * where @f$ n@f$ fulfills the requirements (see FindMaxWeight). 
     *
     * If the denominator is zero, then 1 is returned. 
     *
     * See also AliForwardUtil::ILandauGaus and AliForwardUtil::NLandauGaus
     * for more information on the evaluated functions. 
     * 
     * @param x           Where to evaluate 
     * @param maxN 	  @f$ \max{N}@f$      
     * 
     * @return @f$ f_W(x;\Delta,\xi,\sigma')@f$.  
     */
    Double_t EvaluateWeighted(Double_t x, 
			      UShort_t maxN=9999) const;
    /** 
     * Find the maximum weight to use.  The maximum weight is the
     * largest i for which 
     * 
     * - @f$ i \leq \max{N}@f$ 
     * - @f$ a_i > \min{a}@f$ 
     * - @f$ \delta a_i/a_i > \delta_{max}@f$ 
     * 
     * @param maxRelError @f$ \min{a}@f$ 
     * @param leastWeight @f$ \delta_{max}@f$ 
     * @param maxN        @f$ \max{N}@f$      
     * 
     * @return The largest index @f$ i@f$ for which the above
     * conditions hold.  Will never return less than 1. 
     */
    Int_t FindMaxWeight(Double_t maxRelError=2*fgMaxRelError, 
			Double_t leastWeight=fgLeastWeight, 
			UShort_t maxN=999) const;
    /* @} */
    /** 
     * @{
     * @name TObject Sortable interface 
     */
    /** 
     * Declare this object as sortable 
     * 
     * @return Always true 
     */
    Bool_t IsSortable() const { return kTRUE; }
    /** 
     * Compare to another ELossFit object. 
     * 
     * - +1, if this quality is better (larger) than other objects quality
     * - -1, if this quality is worse (smaller) than other objects quality
     * - +1, if this @f$|\chi^2/\nu-1|@f$ is smaller than the same for other
     * - -1, if this @f$|\chi^2/\nu-1|@f$ is larger than the same for other
     * - 0 otherwise 
     * 
     * @param o Other object to compare to 
     */
    Int_t Compare(const TObject* o) const;
    /* @} */
    /** 
     * @{ 
     * name Auxiliary member functions  
     */
    /** 
     * Information to standard output 
     * 
     * @param option Not used 
     */
    void Print(Option_t* option) const; // *MENU*
    /** 
     * Draw this fit 
     * 
     * @param option Options 
     *  - COMP  Draw components too 
     */
    void Draw(Option_t* option="comp"); // *MENU*
    /** 
     * Browse this object 
     * 
     * @param b Browser
     */
    void Browse(TBrowser* b);
    /** 
     * Get the name of this object 
     * 
     * 
     * @return 
     */
    const Char_t* GetName() const;
    /** 
     * Calculate the quality 
     */
    void CalculateQuality(Double_t maxChi2nu=fgMaxChi2nu, 
			  Double_t maxRelError=fgMaxRelError, 
			  Double_t leastWeight=fgLeastWeight);
    /* @} */
    ClassDef(ELossFit,1); // Result of fit 
  };

  /** 
   * Default constructor 
   */
  AliFMDCorrELossFit();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDCorrELossFit(const AliFMDCorrELossFit& o);
  /** 
   * Destructor 
   */
  virtual ~AliFMDCorrELossFit(); 
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliFMDCorrELossFit& operator=(const AliFMDCorrELossFit& o);

  /** 
   * @{ 
   * @name Set fits 
   */
  /** 
   * Set the fit parameters from a function 
   * 
   * @param d        Detector
   * @param r        Ring 
   * @param eta      Eta 
   * @param quality  Quality flag
   * @param f        Function from fit 
   */  
  Bool_t SetFit(UShort_t d, Char_t r, Double_t eta, Int_t quality, 
		const TF1& f);
  /** 
   * Set the fit parameters from a function 
   * 
   * @param d    Detector
   * @param r    Ring 
   * @param eta  Eta 
   * @param f    ELoss fit result - note, the object will take ownership
   */  
  Bool_t SetFit(UShort_t d, Char_t r, Double_t eta, ELossFit* f);
  /** 
   * Set the fit parameters from a function 
   * 
   * @param d       Detector
   * @param r       Ring 
   * @param etaBin  Eta (bin number, 1->nBins)
   * @param f       ELoss fit result - note, the object will take ownership
   */  
  Bool_t SetFit(UShort_t d, Char_t r, Int_t etaBin, ELossFit* f);
  /** 
   * Set the fit parameters from a function 
   * 
   * @param d         Detector number
   * @param r         Ring identifier 
   * @param eta       Eta value
   * @param quality   Quality flag
   * @param n         @f$ N@f$ - Number of fitted peaks
   * @param chi2      @f$ \chi^2 @f$
   * @param nu        @f$ \nu @f$ - number degrees of freedom
   * @param c         @f$ C@f$ - scale constant
   * @param ec        @f$ \delta C@f$ - error on @f$ C@f$ 
   * @param delta     @f$ \Delta@f$ - most probable value
   * @param edelta    @f$ \delta\Delta@f$ - error on @f$\Delta@f$ 
   * @param xi        @f$ \xi@f$ - Landau width		  
   * @param exi       @f$ \delta\xi@f$ - error on @f$\xi@f$ 
   * @param sigma     @f$ \sigma@f$ - Gaussian width
   * @param esigma    @f$ \delta\sigma@f$ - error on @f$\sigma@f$ 
   * @param sigman    @f$ \sigma_n@f$ - Noise width		  
   * @param esigman   @f$ \delta\sigma_n@f$ - error on @f$\sigma_n@f$ 
   * @param a         Array of @f$ N-1@f$ weights @f$ a_i@f$ for 
   *                  @f$ i=2,\ldots@f$ 
   * @param ea        Array of @f$ N-1@f$ errors on weights @f$ a_i@f$ for 
   *                  @f$ i=2,\ldots@f$ 
   */
  Bool_t SetFit(UShort_t  d,      Char_t    r, Double_t eta, 
		Int_t     quality,UShort_t  n, 
		Double_t  chi2,   UShort_t  nu, 
		Double_t  c,      Double_t  ec, 
		Double_t  delta,  Double_t  edelta, 
		Double_t  xi,     Double_t  exi,
		Double_t  sigma,  Double_t  esigma, 
		Double_t  sigman, Double_t  esigman, 
		Double_t* a,      Double_t* ea);
  /* @} */
  
  /** 
   * @{
   * @name Set and get eta axis
   */
  /** 
   * Set the eta axis to use 
   * 
   * @param axis Eta axis 
   */
  void SetEtaAxis(const TAxis& axis);
  /** 
   * Set the eta axis to use 
   * 
   * @param nBins Number of bins 
   * @param min   Minimum @f$ \eta@f$
   * @param max   maximum @f$ \eta@f$
   */
  void SetEtaAxis(Int_t nBins, Double_t min, Double_t max);
  /** 
   * Get the eta axis used
   * 
   * @return 
   */
  const TAxis& GetEtaAxis() const { return fEtaAxis; }
  /** 
   * Set the low cut used when fitting 
   * 
   * @param cut Cut value 
   */
  void SetLowCut(Double_t cut) { fLowCut = cut; }
  /** 
   * Get the low cut used when fitting 
   * 
   * @return Low cut used for fitting 
   */
  Double_t GetLowCut() const { return fLowCut; }
  /** 
   * Find the eta bin corresponding to the given eta 
   * 
   * @param eta  Eta value 
   * 
   * @return Bin (in @f$[1,N_{bins}]@f$) corresponding to the given
   * eta, or 0 if out of range.
   */
  Int_t FindEtaBin(Double_t eta) const;
  /* @} */

  /**
   * @{						
   * @name Find fits 
   */
  /** 
   * Find the fit corresponding to the specified parameters 
   * 
   * @param d   Detector 
   * @param r   Ring 
   * @param eta Eta value 
   * 
   * @return Fit parameters or null in case of problems 
   */
  ELossFit* FindFit(UShort_t d, Char_t r, Double_t eta) const;
  /** 
   * Find the fit corresponding to the specified parameters 
   * 
   * @param d      Detector 
   * @param r      Ring 
   * @param etabin Eta bin (1 based)
   * 
   * @return Fit parameters or null in case of problems 
   */
  ELossFit* FindFit(UShort_t d, Char_t r, Int_t etabin) const;
  /** 
   * Find the fit corresponding to the specified parameters 
   * 
   * @param d   Detector 
   * @param r   Ring 
   * @param eta Eta value 
   * 
   * @return Fit parameters or null in case of problems 
   */
  ELossFit* GetFit(UShort_t d, Char_t r, Double_t eta) const;
  /** 
   * Find the fit corresponding to the specified parameters 
   * 
   * @param d      Detector 
   * @param r      Ring 
   * @param etabin Eta bin (1 based)
   * 
   * @return Fit parameters or null in case of problems 
   */
  ELossFit* GetFit(UShort_t d, Char_t r, Int_t etabin) const;
  /* @} */

  /**						
   * @{ 
   * @name Miscellaneous
   */
  /** 
   * Get the ring array corresponding to the specified ring
   * 
   * @param d Detector 
   * @param r Ring 
   * 
   * @return Pointer to ring array, or null in case of problems
   */
  TObjArray* GetRingArray(UShort_t d, Char_t r) const;
  /** 
   * Signal that this is a folder
   * 
   * @return Always true 
   */
  Bool_t IsFolder() const { return true; }
  /** 
   * Browse this object 
   * 
   * @param b 
   */
  void Browse(TBrowser* b);
  /** 
   * Draw this object 
   *
   * @param option Options.  Possible values are 
   *  - error Plot error bars 
   *  - relative Plot relative errors
   */
  void Draw(Option_t* option=""); //*MENU*
  /** 
   * Print this object.  
   * 
   * @param option Options 
   *   - R   Print recursive  
   *
   */
  void Print(Option_t* option="R") const; //*MENU*
  /* @} */
protected:
  /** 
   * Get the ring array corresponding to the specified ring
   * 
   * @param d Detector 
   * @param r Ring 
   * 
   * @return Pointer to ring array, or newly created container 
   */
  TObjArray* GetOrMakeRingArray(UShort_t d, Char_t r);

  TObjArray  fRings;    // Array of rings
  TAxis      fEtaAxis;  // Eta axis used
  Double_t   fLowCut;   // Low cut used when fitting 

  ClassDef(AliFMDCorrELossFit,2); 
};

//____________________________________________________________________
inline void 
AliFMDCorrELossFit::SetEtaAxis(Int_t nBins, Double_t min, Double_t max)
{
  fEtaAxis.Set(nBins, min, max);
}
//____________________________________________________________________
inline void 
AliFMDCorrELossFit::SetEtaAxis(const TAxis& e)
{
  fEtaAxis.Set(e.GetNbins(), e.GetXmin(), e.GetXmax());
}
//____________________________________________________________________
inline Double_t
AliFMDCorrELossFit::ELossFit::GetA(UShort_t i) const
{
  if (i <  1)   return 0;
  if (i >  fN)  return 0;
  if (i == 1)   return 1;
  return fA[i-2];
}
//____________________________________________________________________
inline Double_t
AliFMDCorrELossFit::ELossFit::GetEA(UShort_t i) const
{
  if (i <  1)   return 0;
  if (i >  fN)  return 0;
  if (i == 1)   return 1;
  return fEA[i-2];
}


#endif
// Local Variables:
//   mode: C++ 
// End:

