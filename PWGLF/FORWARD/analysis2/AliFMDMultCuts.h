#ifndef ALIFMDMULTCUTS_H
#define ALIFMDMULTCUTS_H
#include <TObject.h>
class TH2;

/**
 * Cuts used when calculating the multiplicity.
 *
 * We can define our cuts in four ways (in order of priorty)
 *
 * - Using a fixed value @f$ v@f$- AliFMDMultCuts:: SetMultCuts
 * - Using a fraction @f$ f@f$ of the most probably value (@f$ \Delta_p@f$)
 *   from the energy loss fits - AliFMDMultCuts::SetMPVFraction
 * - Using some number @f$ n@$ of widths (@f$ \xi@f$) below the most
 *   probable value (@f$ \Delta_p@f$) from the energy loss fits -
 *   possibly including the Gaussian variance (@f$ \sigma@f$) -
 *   AliFMDMultCuts::SetNXi and AliFMDMultCuts::SetIncludeSigma
 * - Using the @f$ x@f$ value for which @f$ P(x>p)@f$ given some cut
 *   value @f$ p@f$
 * - Using the lower fit range of the energy loss fits
 *
 * Which method used depends on the settings of @f$ v@f$, @f$ f@f$,
 * and @f$ n@f$:
 *
 * - if @f$ v > 0@f$ then give @f$ v@f$ 
 * - if @f$ f > 0@f$ then give @f$ f\Delta_p@f$ 
 * - if @f$ n > 0@f$ and @f$\sigma@f$ included then give 
 *   @f$ \Delta_p - n(\xi+\sigma)@f$ 
 * - if @f$ n > 0@f$ then give @f$ \Delta_p - n\xi@f$ 
 * - if @f$ p > 0@f$ then give @f$ x@f$ for which @f$ P(x>p)@f$ 
 * - otherwise, give lower bound on fit range 
 *
 * The member function AliFMDMultCuts::Reset resets all cut values,
 * meaning the lower bound on the fits will be used by default.  This
 * is useful to ensure a fresh start:
 *
 * @code 
 AliFMDMultCuts c;
 c.Reset();
 c.SetNXi(2);
 @endcode 
 *
 * The member function AliFMDMultCuts::GetMethod will return the
 * method identifier for the current method employed
 * (AliFMDMultCuts::EMethod).  Like wise will the method
 * AliFMDMultCuts::GetMethodString give a human readable string of the
 * current method employed.
 */
class AliFMDMultCuts : public TObject 
{
public:
  enum EMethod { 
    kFixed, 
    kMPVFraction, 
    kFitRange, 
    kLandauWidth,
    kProbability
  };
  /** 
   * CTOR
   */
  AliFMDMultCuts();
  /** 
   * Set the cut for specified method.
   * 
   * @param method Method to use 
   * @param cut1   1st cut value
   * @param cut2   2nd cut value (ignored for method!=kFixed)
   * @param cut3   3rd cut value (ignored for method!=kFixed)
   * @param cut4   4th cut value (ignored for method!=kFixed) 
   * @param cut5   5th cut value (ignored for method!=kFixed)
   */
   AliFMDMultCuts(EMethod method, 
		  Double_t cut1, 
		  Double_t cut2=-1, 
		  Double_t cut3=-1, 
		  Double_t cut4=-1, 
		  Double_t cut5=-1);
  /** 
   * Copy CTOR
   * 
   * @param o Oject to copy from
   */
  AliFMDMultCuts(const AliFMDMultCuts& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to asing from 
   * 
   * @return Reference to this object 
   */
  AliFMDMultCuts& operator=(const AliFMDMultCuts& o);
  /** 
   * Reset all cuts to default value. 
   */
  void Reset();
  /** 
   * Get the multiplicity cuts for a specific ring and pseudo-rapidity 
   * 
   * @param d      Detector 
   * @param r      Ring  
   * @param eta    pseudo-rapidity 
   * @param errors Get error on cut
   * 
   * @return Cut value 
   */
  Double_t GetMultCut(UShort_t d, Char_t r, Double_t eta, Bool_t errors) const;
  /** 
   * Get the multiplicity cuts for a specific ring and pseudo-rapidity 
   * 
   * @param d      Detector 
   * @param r      Ring  
   * @param etabin pseudo-rapidity bin
   * @param errors Get error on cut
   * 
   * @return Cut value 
   */
  Double_t GetMultCut(UShort_t d, Char_t r, Int_t etabin, Bool_t errors) const;
  /** 
   * Clear the cuts 
   * 
   */
  void UnsetMultCuts() { SetMultCuts(-1); }
  /** 
   * Set the cuts
   * 
   * @param fmd1i  Value for FMD1i 
   * @param fmd2i  Value for FMD2i 
   * @param fmd2o  Value for FMD2o 
   * @param fmd3i  Value for FMD3i 
   * @param fmd3o  Value for FMD3o 
   */
  void SetMultCuts(Double_t fmd1i, 
		   Double_t fmd2i=-1, 
		   Double_t fmd2o=-1, 
		   Double_t fmd3i=-1, 
		   Double_t fmd3o=-1);
  /** 
   * Set the faction of most probable value 
   * 
   * @param frac Fraction
   */
  void SetMPVFraction(Double_t frac=0) { fMPVFraction = frac; }
  /** 
   * Set the number times the Landau width
   * 
   * @param nXi Number of widths 
   */
  void SetNXi(Double_t nXi) { fNXi = nXi; }
  /** 
   * Set wether to includle the Gaussian sigma in cut calculation
   * 
   * @param in If true, use Gaussian variance 
   */
  void SetIncludeSigma(Bool_t in) { fIncludeSigma = in; }
  /** 
   * Set probability cut.  See
   * AliFMDCorrELossFit::ELossFit::FindProbabilityCut
   * 
   * @param cut Cut value 
   */
  void SetProbability(Double_t cut=1e-5) { fProbability = cut; }
  /** 
   * Set the cut for specified method.
   *
   * Note, that if @a method is kFixed, and only @a cut1 is specified,
   * then the outer rings cut value is increased by 20% relative to @a
   * cut1.
   *
   * Also note, that if @a method is kLandauWidth, and @a cut2 is
   * larger than zero, then @f$\sigma@f$ of the fits are included in
   * the cut value.
   * 
   * @param method Method to use 
   * @param cut1   1st cut value
   * @param cut2   2nd cut value (ignored for method!=kFixed)
   * @param cut3   3rd cut value (ignored for method!=kFixed)
   * @param cut4   4th cut value (ignored for method!=kFixed) 
   * @param cut5   5th cut value (ignored for method!=kFixed)
   */
  void Set(EMethod method, 
	   Double_t cut1, 
	   Double_t cut2=-1, 
	   Double_t cut3=-1, 
	   Double_t cut4=-1, 
	   Double_t cut5=-1);
  /** 
   * Print information
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /** 
   * Fill a histogram with cut values.  The histogram is assumed to
   * have rings on the y axis, and @f$ \eta@f$ on the x axis.
   * 
   * @param h Histogram to fill 
   */
  void FillHistogram(TH2* h) const;
  /** 
   * Define outputs 
   * 
   * @param l      List to add to 
   * @param name   Name 
   */
  void Output(TList* l, const char* name=0) const;
  /** 
   * Read in cuts stored in file 
   * 
   * @param l      List to read from 
   * @param name   Name of possible sub-list
   * 
   * @return true on success
   */
  Bool_t Input(TList* l, const char* name);
  /** 
   * Get a fixed cut value 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Fixed cut value 
   */
  Double_t GetFixedCut(UShort_t d, Char_t r) const;
  /** 
   * Return methid Id
   * 
   * @return method number
   */
  UShort_t GetMethod() const;
  /** 
   * Get a string descriping the method used
   * 
   * @return String 
   */
  const char* GetMethodString() const;
protected:
  Double_t fMultCuts[5];   // Fixed cuts
  Double_t fMPVFraction;   // Most probably value fraction
  Double_t fNXi;           // Times of Landau width
  Bool_t   fIncludeSigma;  // Include Gaussian variance 
  Double_t fProbability;   // Probability cut

  ClassDef(AliFMDMultCuts,4); // Cuts on ESD Mult 
};

#endif
// Local Variables: 
//  mode: C++
// End:
