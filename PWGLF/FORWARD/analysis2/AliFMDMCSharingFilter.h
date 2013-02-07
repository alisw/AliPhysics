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
 * @ingroup pwglf_forward_aod
 */
#include "AliFMDSharingFilter.h"
#include "AliFMDMCTrackDensity.h"

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
 *      hit strips for each vertex bin (if enabled - see SetupForData method)
 * 
 *
 * @ingroup pwglf_forward_algo
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_aod
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
    fTrackDensity(),
    fFMD1i(0),
    fFMD2i(0),
    fFMD2o(0),
    fFMD3i(0),
    fFMD3o(0), 
    fOperComp(0)
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
   * Return the track density calculator 
   * 
   * @return Track density calculator 
   */
  const AliFMDMCTrackDensity& GetTrackDensity() const { return fTrackDensity; }
  /** 
   * Return the track density calculator 
   * 
   * @return Track density calculator 
   */
  AliFMDMCTrackDensity& GetTrackDensity() { return fTrackDensity; }

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
  void CreateOutputObjects(TList* dir);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir     Where the output is 
   * @param output  Output list
   * @param nEvents Number of events 
   */
  void Terminate(const TList* dir, TList* output, Int_t nEvents);
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;

  virtual void SetDebug(Int_t dbg=1);
protected:
  AliFMDMCTrackDensity fTrackDensity;
  TH2D* fFMD1i;      // ESD-MC correlation 
  TH2D* fFMD2i;      // ESD-MC correlation 
  TH2D* fFMD2o;      // ESD-MC correlation 
  TH2D* fFMD3i;      // ESD-MC correlation 
  TH2D* fFMD3o;      // ESD-MC correlation 
  TH2I* fOperComp;   // Operation vs # trackrefs
  ClassDef(AliFMDMCSharingFilter,2); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
