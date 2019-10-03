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
 *   from the energy loss fits
 * - Using some number @f$ n@f$ of widths (@f$ \xi@f$) below the most
 *   probable value (@f$ \Delta_p@f$) from the energy loss fits 
 * - Using some number @f$ n@f$ of widths (@f$ \xi+\sigma@f$) below the
 *   most probable value (@f$ \Delta_p@f$) from the energy loss fits 
 * - Using the @f$ x@f$ value for which @f$ P(x>p)@f$ given some cut
 *   value @f$ p@f$
 * - Using the lower fit range of the energy loss fits
 *
 * The member function AliFMDMultCuts::Reset resets all cut values,
 * meaning the lower bound on the fits will be used by default.  This
 * is useful to ensure a fresh start:
 *
 * @code 
 AliFMDMultCuts c;
 c.Reset();
 c.Set(AliFMDMultCuts::kLandauWidth,2);
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
    kLandauSigmaWidth,
    kProbability,
    kLandauSigmaAvg    
  };
  enum {
    kBad = 1024
  };
  
  /** 
   * CTOR
   */
  AliFMDMultCuts();
  /** 
   * Set the cut for specified method.
   * 
   * @param method Cut method
   * @param fmd1i  Value for FMD1i
   * @param fmd2i  Value for FMD2i (if < 0, use fmd1i) 
   * @param fmd2o  Value for FMD2o (if < 0, use fmd1i) 
   * @param fmd3i  Value for FMD3i (if < 0, use fmd1i) 
   * @param fmd3o  Value for FMD3o (if < 0, use fmd1i) 
   */
   AliFMDMultCuts(EMethod method, 
		  Double_t fmd1i, 
		  Double_t fmd2i=-1, 
		  Double_t fmd2o=-1, 
		  Double_t fmd3i=-1, 
		  Double_t fmd3o=-1);
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
   * Get the multiplicity cuts for a specific ring and pseudo-rapidity bin
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
   * Set the cuts
   * 
   * @param fmd1i  Value for FMD1i
   * @param fmd2i  Value for FMD2i (if < 0, use fmd1i) 
   * @param fmd2o  Value for FMD2o (if < 0, use fmd1i) 
   * @param fmd3i  Value for FMD3i (if < 0, use fmd1i) 
   * @param fmd3o  Value for FMD3o (if < 0, use fmd1i) 
   *
   * @deprecated Use AliFMDMultCuts::Set 
   */
  void SetMultCuts(Double_t fmd1i, 
		   Double_t fmd2i=-1, 
		   Double_t fmd2o=-1, 
		   Double_t fmd3i=-1, 
		   Double_t fmd3o=-1) 
  { 
    DepSet("SetMultCuts", kFixed,fmd1i,fmd2i,fmd2o,fmd3i,fmd3o); 
  }
  /** 
   * Set the faction of most probable value 
   * 
   * @param frac Fraction
   *
   * @deprecated Use AliFMDMultCuts::Set 
   */
  void SetMPVFraction(Double_t frac=0) 
  { 
    DepSet("SetMPVFraction",kMPVFraction, frac); 
  }
  /** 
   * Set the number times the Landau width
   * 
   * @param nXi Number of widths 
   *
   * @deprecated Use AliFMDMultCuts::Set 
   */
  void SetNXi(Double_t nXi) { DepSet("SetNXi", kLandauWidth, nXi); }
  /** 
   * Set wether to includle the Gaussian sigma in cut calculation
   * 
   * @param in If true, use Gaussian variance 
   *
   * @deprecated Use AliFMDMultCuts::Set with kLandauSigmaWidth
   */
  void SetIncludeSigma(Bool_t in);
  /** 
   * Set probability cut.  See
   * AliFMDCorrELossFit::ELossFit::FindProbabilityCut
   * 
   * @param cut Cut value 
   *
   * @deprecated Use AliFMDMultCuts::Set 
   */
  void SetProbability(Double_t cut=1e-5) 
  { 
    DepSet("SetProbability", kProbability, cut); 
  }
  /** 
   * Set the cut for specified method.
   *
   * Note, that if @a method is kFixed, and only @a fmd1i is specified,
   * then the outer rings cut value is increased by 20% relative to @a
   * fmd1i.
   *
   * Also note, that if @a method is kLandauWidth, and @a cut2 is
   * larger than zero, then @f$\sigma@f$ of the fits are included in
   * the cut value.
   * 
   * @param method Method to use 
   * @param fmd1i  Value for FMD1i
   * @param fmd2i  Value for FMD2i (if < 0, use fmd1i) 
   * @param fmd2o  Value for FMD2o (if < 0, use fmd1i) 
   * @param fmd3i  Value for FMD3i (if < 0, use fmd1i) 
   * @param fmd3o  Value for FMD3o (if < 0, use fmd1i) 
   */
  void Set(EMethod method, 
	   Double_t fmd1i, 
	   Double_t fmd2i=-1, 
	   Double_t fmd2o=-1, 
	   Double_t fmd3i=-1, 
	   Double_t fmd3o=-1);
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
   * Return methid Id
   * 
   * @return method number
   */
  EMethod GetMethod() const { return fMethod; }
  /** 
   * Get a string descriping the method used
   * 
   * @return String 
   */
  const char* GetMethodString(Bool_t latex=false) const;
  /** 
   * helper function to translate a method into a string 
   * 
   * @param method Method identifier 
   * @param latex  IF true, format for TLatex
   * 
   * @return String 
   */
  static const char* Method2String(EMethod method, Bool_t latex);
  /** 
   * Helper function to translate a string into a method identifier 
   * 
   * @param str String
   * 
   * @return Method identifier 
   */
  static EMethod String2Method(const char* str);
protected:
  /** 
   * Set the cut for specified method.
   *
   * @param what   What to set 
   * @param method Cut method
   * @param fmd1i  Value for FMD1i
   * @param fmd2i  Value for FMD2i (if < 0, use fmd1i) 
   * @param fmd2o  Value for FMD2o (if < 0, use fmd1i) 
   * @param fmd3i  Value for FMD3i (if < 0, use fmd1i) 
   * @param fmd3o  Value for FMD3o (if < 0, use fmd1i) 
   */
  void DepSet(const char* what, 
	      EMethod method, 
	      Double_t fmd1i, 
	      Double_t fmd2i=-1, 
	      Double_t fmd2o=-1, 
	      Double_t fmd3i=-1, 
	      Double_t fmd3o=-1);
  /** 
   * Get the cut parameter for a specific ring 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Cut parameter
   */
  Double_t GetCutParam(UShort_t d, Char_t r) const;

  Double_t fCuts[5];   // Fixed cuts
  EMethod  fMethod;

  ClassDef(AliFMDMultCuts,5); // Cuts on ESD Mult 
};

#endif
// Local Variables: 
//  mode: C++
// End:
