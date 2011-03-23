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
 * @ingroup pwg2_forward_aod
 */
#include <TNamed.h>
#include <TH2.h>
#include <TList.h>
#include "AliForwardUtil.h"
class AliESDFMD;
class TAxis;
class TList;
class TH2;

/**
 * Class to do the sharing correction.  That is, a filter that merges 
 * adjacent strip signals presumably originating from a single particle 
 * that impinges on the detector in such a way that it deposite energy 
 * into two or more strips. 
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
 *      hit strips for each vertex bin (if enabled - see Init method)
 * 
 *
 * @ingroup pwg2_forward_algo 
 * @ingroup pwg2_forward_aod
 */
class AliFMDSharingFilter : public TNamed
{
public: 
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
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDSharingFilter(const AliFMDSharingFilter& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliFMDSharingFilter& operator=(const AliFMDSharingFilter& o);

  /** 
   * Set the low cut used for sharing 
   * 
   * @param lowCut Low cut
   */
  void SetLowCut(Double_t lowCut=0) { fLowCut = lowCut; }
  /** 
   * Reset the low cut for sharing to use the fit range lower cut 
   * 
   */
  void UnsetLowCut() { fLowCut = 0; }
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1) { fDebug = dbg; }

  /** 
   * Enable use of angle corrected signals in the algorithm 
   * 
   * @param use If true, use angle corrected signals, 
   * otherwise use de-corrected signals.  In the final output, the 
   * signals are always angle corrected. 
   */
  void UseAngleCorrectedSignals(Bool_t use) { fCorrectAngles = use; }
  /** 
   * Set the number of landau width to subtract from the most probably
   * value to get the high cut for the merging algorithm.
   * 
   * @param n Number of @f$ \xi@f$ 
   */
  void SetNXi(Short_t n) { fNXi = n; }
  /** 
   * Filter the input AliESDFMD object
   * 
   * @param input     Input 
   * @param lowFlux   If this is a low-flux event 
   * @param output    Output AliESDFMD object 
   * 
   * @return True on success, false otherwise 
   */
  Bool_t Filter(const AliESDFMD& input, 
		Bool_t           lowFlux, 
		AliESDFMD&       output);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir     Where the output is 
   * @param nEvents Number of events 
   */
  virtual void ScaleHistograms(const TList* dir, Int_t nEvents);
  
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  virtual void DefineOutput(TList* dir);
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  virtual void Print(Option_t* option="") const;
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
     * Clear this object
     */
    void Clear(const Option_t* ="") { fNHits = 0; } 
    /** 
     * Increase number of hits 
     * 
     */
    void Incr() { fNHits++; } 
    /** 
     * Finish off 
     * 
     */
    void Finish(); 
    /** 
     * Make output 
     * 
     * @param dir where to store 
     */
    void Output(TList* dir);
    /** 
     * Scale the histograms to the total number of events 
     * 
     * @param nEvents Number of events 
     * @param dir     Where the output is 
     */
    void ScaleHistograms(const TList* dir, Int_t nEvents);
    TH1D*     fBefore;       // Distribution of signals before filter
    TH1D*     fAfter;        // Distribution of signals after filter
    TH1D*     fHits;         // Distribution of hit strips. 
    Int_t     fNHits;        // Number of hit strips per event
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
   * The actual algorithm 
   * 
   * @param mult      The unfiltered signal in the strip
   * @param eta       Psuedo rapidity 
   * @param prevE     Previous strip signal (or 0)
   * @param nextE     Next strip signal (or 0) 
   * @param lowFlux   Whether this is a low flux event 
   * @param d         Detector
   * @param r         Ring 
   * @param s         Sector 
   * @param t         Strip
   * @param usedPrev  Whether the previous strip was used in sharing or not
   * @param usedThis  Wether this strip was used in sharing or not. 
   * 
   * @return The filtered signal in the strip
   */
  Double_t MultiplicityOfStrip(Double_t mult,
			       Double_t eta,
			       Double_t prevE,
			       Double_t nextE,
			       Bool_t   lowFlux,
			       UShort_t d,
			       Char_t   r,
			       UShort_t s,
			       UShort_t t,
			       Bool_t&  usedPrev, 
			       Bool_t&  usedThis) const;
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
   */
  virtual Double_t GetHighCut(UShort_t d, Char_t r, Double_t eta) const;
  /**
   * Get the low cut.  Normally, the low cut is taken to be the lower
   * value of the fit range used when generating the energy loss fits.
   * However, if fLowCut is set (using SetLowCit) to a value greater
   * than 0, then that value is used.
   */
  virtual Double_t GetLowCut() const;

  TList    fRingHistos;    // List of histogram containers
  Double_t fLowCut;        // Low cut on sharing
  Bool_t   fCorrectAngles; // Whether to work on angle corrected signals
  Short_t  fNXi;           // Number of xi's from Delta to stop merging
  Int_t    fDebug;         // Debug level 

  ClassDef(AliFMDSharingFilter,1); //
};

#endif
// Local Variables:
//  mode: C++ 
// End:
