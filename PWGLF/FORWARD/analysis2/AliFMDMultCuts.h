#ifndef ALIFMDMULTCUTS_H
#define ALIFMDMULTCUTS_H
#include <TObject.h>

/**
 * Cuts used when calculating the multiplicity 
 * 
 */
class AliFMDMultCuts : public TObject 
{
public:
  /** 
   * CTOR
   */
  AliFMDMultCuts();
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
   * Print information
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /** 
   * Define outputs 
   * 
   * @param l      List to add to 
   * @param name   Name 
   */
  void Output(TList* l, const char* name=0) const;
  /** 
   * Get a fixed cut value 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Fixed cut value 
   */
  Double_t GetFixedCut(UShort_t d, Char_t r) const;
protected:
  Double_t fMultCuts[5];   // Fixed cuts
  Double_t fMPVFraction;   // Most probably value fraction
  Double_t fNXi;           // Times of Landau width
  Bool_t   fIncludeSigma;  // Include Gaussian variance 
   
  ClassDef(AliFMDMultCuts,1); // Cuts on ESD Mult 
};

#endif
// Local Variables: 
//  mode: C++
// End:
