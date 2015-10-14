//
// Class to do the sharing correction.  That is, a filter that merges 
// adjacent strip signals presumably originating from a single particle 
// that impinges on the detector in such a way that it deposite energy 
// into two or more strips. 
//
#ifndef ALIFMDSHARINGFILTER_H
#define ALIFMDSHARINGFILTER_H
/**
 * @file   AliFMDSharingFilter.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:03:57 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include <TNamed.h>
#include <TH2.h>
#include <TList.h>
#include "AliForwardUtil.h"
#include "AliFMDMultCuts.h"
class AliESDFMD;
class TAxis;
class TList;
class TH2;
class AliFMDFloatMap;

/**
 * Class to do the sharing correction.  That is, a filter that merges 
 * adjacent strip signals presumably originating from a single particle 
 * that impinges on the detector in such a way that it deposite energy 
 * into two or more strips. 
 *
 * @image html alice-int-2012-040-share_fraction.png "Energy loss sharing"
 *
 * @par Input: 
 *    - AliESDFMD object  - from reconstruction
 *
 * @par Output: 
 *    - AliESDFMD object  - copy of input, but with signals merged 
 *
 * @par Corrections used: 
 *    - AliFMDCorrELossFit
 *
 * @par Histograms: 
 *    - For each ring (FMD1i, FMD2i, FMD2o, FMD3i, FMD3o) the distribution of 
 *      signals before and after the filter.  
 *    - For each ring (see above), an array of distributions of number of 
 *      hit strips for each vertex bin (if enabled - see SetupForData method)
 * 
 *
 * @ingroup pwglf_forward_algo 
 */
class AliFMDSharingFilter : public TNamed
{
public: 
  /** 
   * Status of a strip 
   * @deprecated Not used
   */
  enum Status { 
    /** Nothing yet */
    kNone             = 1, 
    /** Candidate for merging */
    kCandidate        = 2, 
    /** This was merged into other strip */
    kMergedWithOther  = 3, 
    /** Other strips was merged into this */
    kMergedInto       = 4
  };
  /** 
   * Destructor
   */
  virtual ~AliFMDSharingFilter();
  /** 
   * Default Constructor - do not use 
   */
  AliFMDSharingFilter();
  /** 
   * Constructor 
   * 
   * @param title Title of object  - not significant 
   */
  AliFMDSharingFilter(const char* title);

  /** 
   * @{ 
   * @name Parameters etc.
   */
  /** 
   * If called with a true argument, then merging is wholy disabled 
   * 
   * @param disable If true, disable merging altogether 
   */
  virtual void SetMergingDisabled(Bool_t disable) {fMergingDisabled = disable; }
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  virtual void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /** 
   * Enable use of angle corrected signals in the algorithm 
   * 
   * @param use If true, use angle corrected signals, 
   * otherwise use de-corrected signals.  In the final output, the 
   * signals are always angle corrected. 
   */
  void SetUseAngleCorrectedSignals(Bool_t use) { fCorrectAngles = use; }
   /** 
   * Enable zeroing of signals if below high cut
   * 
   * @param use zero the signals if below sharing cut
   * 
   */
  void SetZeroSharedHitsBelowThreshold(Bool_t use) { fZeroSharedHitsBelowThreshold = use; }
 /** 
   * Enable a simpler merging algorithm
   * 
   * @param use use the simpler algorithm
   * 
   */
  void SetUseSimpleSharing(Bool_t use) { fUseSimpleMerging = use; }
  /** 
   * In case of a simpler merging algorithm allow 3 strips to be 
   * merged
   * 
   * @param use allow three strips
   * 
   */
  void SetAllow3Strips(Bool_t use) { fThreeStripSharing = use; }  
  /**
   * Set whether to ignore the ESD info when angle correcting, this
   * is to counter a known issue where the info in the ESD is incorrect
   * 
   * @param use ignore the ESD info
   */
  void SetIgnoreESDWhenAngleCorrecting(Bool_t use) { fIgnoreESDForAngleCorrection = use; }
  /* @} */

  /** 
   * @{ 
   * @name Processing 
   */
  /** 
   * Initialize 
   * 
   * @param axis Default eta axis from parent task 
   */
  void SetupForData(const TAxis& axis);
  /** 
   * Filter the input AliESDFMD object
   * 
   * @param input     Input 
   * @param lowFlux   If this is a low-flux event 
   * @param output    Output AliESDFMD object 
   * @param zvtx      Vertex position 
   * 
   * @return True on success, false otherwise 
   */
  Bool_t Filter(const AliESDFMD& input, 
		Bool_t           lowFlux, 
		AliESDFMD&       output, 
		Double_t         zvtx);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir     Where the output is 
   * @param output  Output list
   * @param nEvents Number of events 
   */
  virtual void Terminate(const TList* dir, TList* output, Int_t nEvents);  
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  virtual void CreateOutputObjects(TList* dir);
  /* @} */
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  virtual void Print(Option_t* option="") const;

  /** 
   * @{ 
   * @name Cuts
   */
  /** 
   * Get the low cuts 
   * 
   * @return Reference to low cuts
   */
  AliFMDMultCuts& GetLCuts() { return fLCuts; }
  /** 
   * Get the high cuts 
   * 
   * @return Reference to high cuts
   */
  AliFMDMultCuts& GetHCuts() { return fHCuts; }
  /** 
   * Get the low cuts 
   * 
   * @return Reference to low cuts
   */
  const AliFMDMultCuts& GetLCuts() const { return fLCuts; }
  /** 
   * Get the high cuts 
   * 
   * @return Reference to high cuts
   */
  const AliFMDMultCuts& GetHCuts() const { return fHCuts; }
  /** 
   * Set the low cuts 
   * 
   * @param c Cuts object
   */  
  void SetLCuts(const AliFMDMultCuts& c) { fLCuts = c; }
  /** 
   * Set the high cuts 
   * 
   * @param c Cuts object
   */  
  void SetHCuts(const AliFMDMultCuts& c) { fHCuts = c; }
  /* @} */
protected:
  /** 
   * Copy constructor - not implemented
   */
  AliFMDSharingFilter(const AliFMDSharingFilter& o) : TNamed(o) {}
  /** 
   * Assignment operator  - not implemented
   * 
   * @return Reference to this 
   */
  AliFMDSharingFilter& operator=(const AliFMDSharingFilter&){return *this;}
  /** 
   * Internal data structure to keep track of the histograms.  Objects
   * of this class are never streamed.
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
     * Clear this object
     */
    // void Clear(const Option_t* ="") { fNHits = 0; } 
    /** 
     * Increase number of hits 
     * 
     */
    // void Incr() { fNHits++; } 
    /** 
     * Finish off 
     * 
     */
    // void Finish(); 
    /** 
     * Make output 
     * 
     * @param dir where to store 
     */
    void CreateOutputObjects(TList* dir);
    /** 
     * Scale the histograms to the total number of events 
     * 
     * @param nEvents Number of events 
     * @param dir     Where the output is 
     */
    void Terminate(const TList* dir, Int_t nEvents);
    TH1D*     fBefore;          // Distribution of signals before filter
    TH1D*     fAfter;           // Distribution of signals after filter
    TH1D*     fSingle;          // Distribution of 1 signal after filter
    TH1D*     fDouble;          // Distribution of 2 signals after filter
    TH1D*     fTriple;          // Distribution of 3 signals after filter
    TH2D*     fSinglePerStrip;  // Distribution of 1 signal per strip
    TH2D*     fBeforeAfter;     // Correlation of before and after 
    TH2D*     fNeighborsBefore; // Correlation of neighbors 
    TH2D*     fNeighborsAfter;  // Correlation of neighbors 
    TH2D*     fSumESD;          // Summed ESD signal 
    TH2D*     fSum;             // Summed cluster signal 
    TH1D*     fNConsecutive;    // # consecutive strips with signal > low cut
    // ClassDef(RingHistos,4);
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
  /** 
   * Get the signal in a strip 
   * 
   * @param fmd   ESD object
   * @param d     Detector
   * @param r     Ring 
   * @param s     Sector 
   * @param t     Strip
   * 
   * @return The energy signal 
   */
  Double_t SignalInStrip(const AliESDFMD& fmd, 
			 UShort_t d,
			 Char_t   r,
			 UShort_t s,
			 UShort_t t) const;
  /** 
   * Angle correct the signal 
   * 
   * @param mult Angle Un-corrected Signal 
   * @param eta  Pseudo-rapidity 
   * 
   * @return Angle corrected signal 
   */
  Double_t AngleCorrect(Double_t mult, Double_t eta) const;
  /** 
   * Angle de-correct the signal 
   * 
   * @param mult Angle corrected Signal 
   * @param eta  Pseudo-rapidity 
   * 
   * @return Angle un-corrected signal 
   */
  Double_t DeAngleCorrect(Double_t mult, Double_t eta) const;
  /** 
   * Get the high cut.  The high cut is defined as the 
   * most-probably-value peak found from the energy distributions, minus 
   * 2 times the width of the corresponding Landau.
   * 
   * @param d      Detector 
   * @param r      Ring 
   * @param eta    Eta value 
   * @param errors If false, do not show errors 
   * 
   * @return 0 or less on failure, otherwise @f$\Delta_{mp}-n\xi@f$ 
   */
  virtual Double_t GetHighCut(UShort_t d, Char_t r, Double_t eta,
			      Bool_t errors=true) const;
  /**
   * Get the low cut.  Normally, the low cut is taken to be the lower
   * value of the fit range used when generating the energy loss fits.
   * However, if fLowCut is set (using SetLowCit) to a value greater
   * than 0, then that value is used.
   *
   * @param d    Detector
   * @param r    Ring 
   * @param eta  Eta value 
   * 
   * @return 
   */
  virtual Double_t GetLowCut(UShort_t d, Char_t r, Double_t eta) const;
  TList    fRingHistos;      // List of histogram containers
  Bool_t   fCorrectAngles;   // Whether to work on angle corrected signals
  TH2*     fHighCuts;        // High cuts used
  TH2*     fLowCuts;         // High cuts used
  Int_t    fDebug;           // Debug level 
  Bool_t   fZeroSharedHitsBelowThreshold; // Zero shared strip below cut?
  AliFMDMultCuts fLCuts;     // Cuts object for low cuts
  AliFMDMultCuts fHCuts;     // Cuts object for high cuts
  Bool_t   fUseSimpleMerging;// enable simple sharing by HHD
  Bool_t   fThreeStripSharing; //In case of simple sharing allow 3 strips
  Bool_t   fMergingDisabled; // If true, do not merge
  Bool_t   fIgnoreESDForAngleCorrection; // Ignore ESD information when angle correcting
  ClassDef(AliFMDSharingFilter,11); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
