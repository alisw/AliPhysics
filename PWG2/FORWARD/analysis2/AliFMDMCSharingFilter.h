#ifndef ALIFMDMCSHARINGFILTER_H
#define ALIFMDMCSHARINGFILTER_H
#include "AliFMDSharingFilter.h"
class AliMCEvent;

/**
 * Class to do the sharing correction for MC data.
 *
 * @par Input: 
 *    - AliESDFMD object  - from reconstruction
 *
 * @par Output: 
 *    - AliESDFMD object  - copy of input, but with signals merged 
 *
 * @par Corrections used: 
 *    - AliCorrELossFit 
 *
 * @par Histograms: 
 *    - For each ring (FMD1i, FMD2i, FMD2o, FMD3i, FMD3o) the distribution of 
 *      signals before and after the filter.  
 *    - For each ring (see above), an array of distributions of number of 
 *      hit strips for each vertex bin (if enabled - see Init method)
 * 
 *
 * @ingroup pwg2_forward_algo
 * @ingroup pwg2_forward_mc
 */
class AliFMDMCSharingFilter : public AliFMDSharingFilter
{
public: 
  /** 
   * Destructor
   */
  virtual ~AliFMDMCSharingFilter();
  /** 
   * Default Constructor - do not use 
   */
  AliFMDMCSharingFilter();
  /** 
   * Constructor 
   * 
   * @param title Title of object  - not significant 
   */
  AliFMDMCSharingFilter(const char* title);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDMCSharingFilter(const AliFMDMCSharingFilter& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliFMDMCSharingFilter& operator=(const AliFMDMCSharingFilter& o);
  /** 
   * Filter the input AliESDFMD object
   * 
   * @param input     Input (from ESD) - used for eta
   * @param lowFlux   If this is a low-flux event 
   * @param output    Output AliESDFMD object 
   * 
   * @return True on success, false otherwise 
   */
  Bool_t FilterMC(const AliESDFMD&  input, 
		  const AliMCEvent& event,
		  Double_t          vz,
		  AliESDFMD&        output);
  /** 
   * Compare the result of merging to the monte-carlo truth.  This
   * fills the correlation histograms
   * 
   * @param esd  ESD after sharing correction
   * @param mc   MC ESD 
   */
  void CompareResults(const AliESDFMD&  esd, const AliESDFMD&  mc);
		  
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  void DefineOutput(TList* dir);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir     Where the output is 
   * @param nEvents Number of events 
   */
  void ScaleHistograms(TList* dir, Int_t nEvents);
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
protected:
  void StoreParticle(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
		     AliESDFMD& output) const;
  TH2D* fFMD1i;  // ESD-MC correlation 
  TH2D* fFMD2i;  // ESD-MC correlation 
  TH2D* fFMD2o;  // ESD-MC correlation 
  TH2D* fFMD3i;  // ESD-MC correlation 
  TH2D* fFMD3o;  // ESD-MC correlation 
  TH1D* fSumEta; // MC dN/deta 

  ClassDef(AliFMDMCSharingFilter,1); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
