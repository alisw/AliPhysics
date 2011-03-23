//
// Class to do the sharing correction for MC data.
//
#ifndef ALIFMDMCSHARINGFILTER_H
#define ALIFMDMCSHARINGFILTER_H
/**
 * @file   AliFMDMCSharingFilter.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:03:47 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_aod
 */
#include "AliFMDSharingFilter.h"
class AliMCEvent;

/**
 * Class to do the sharing correction for MC data.
 *
 * @par Input: 
 *    - AliESDFMD object  - from reconstruction
 *    - Kinematics
 *    - Track-References
 *
 * @par Output: 
 *    - AliESDFMD object  - copy of input, but with signals merged 
 *
 * @par Corrections used: 
 *    - None
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
 * @ingroup pwg2_forward_aod
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
  AliFMDMCSharingFilter()
  : AliFMDSharingFilter(), 
    fFMD1i(0),
    fFMD2i(0),
    fFMD2o(0),
    fFMD3i(0),
    fFMD3o(0), 
    fSumEta(0)
  {}
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
   * Filter the input kinematics and track references, using 
   * some of the ESD information
   * 
   * @param input   Input ESD event
   * @param event   Input MC event
   * @param vz      Vertex position 
   * @param output  Output ESD-like object
   * @param primary Per-event histogram of primaries 
   *
   * @return True on succes, false otherwise 
   */
  Bool_t FilterMC(const AliESDFMD&  input, 
		  const AliMCEvent& event,
		  Double_t          vz,
		  AliESDFMD&        output,
		  TH2D*             primary);
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
  void ScaleHistograms(const TList* dir, Int_t nEvents);
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
protected:
  /** 
   * Store a particle hit in FMD<i>dr</i>[<i>s,t</i>] in @a output
   * 
   * @param d       Detector
   * @param r       Ring
   * @param s       Sector
   * @param t       Strip
   * @param output  Output ESD object
   */
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
