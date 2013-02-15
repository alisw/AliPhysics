// 
// This class collects the event histograms into single histograms, 
// one for each ring in each vertex bin.  
//
#ifndef ALIFMDHISTCOLLECTOR_H
#define ALIFMDHISTCOLLECTOR_H
/**
 * @file   AliFMDHistCollector.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:03:01 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include <TNamed.h>
#include <TList.h>
#include <TArrayI.h>
#include "AliForwardUtil.h"
class AliESDFMD;
class TH2D;

/** 
 * This class collects the event histograms into single histograms, 
 * one for each ring in each vertex bin.  
 *
 * @par Input:
 *   - AliESDFMD object possibly corrected for sharing
 *
 * @par Output:
 *   - 5 RingHistos objects - each with a number of vertex dependent 
 *     2D histograms of the inclusive charge particle density 
 * 
 * @par HistCollector used: 
 *   - AliFMDCorrSecondaryMap
 *
 * @ingroup pwglf_forward_algo
 * @ingroup pwglf_forward_aod
 */
class AliFMDHistCollector : public TNamed
{
public:
  /** 
   * Methods to use when merging overlapping bins @f$b_1@f$, @f$b_2@f$
   * with content @f$c_1@f$, @f$c_2@f$, and errors @f$e_1@f$,
   * @f$e_2@f$ into bin @f$ b@f$ with content @f$c@f$ and error @f$e@f$ 
   */
  enum MergeMethod {
    /**
     * @f[
     *   c = \frac{1}{2}(c_1+c_2) 
     * @f]
     * @f[
     *   e = \sqrt{e_1^2+e_2^2} 
     * @f]
     */
    kStraightMean,       
    /**
     * As above, exept zero's are ignored 
     */
    kStraightMeanNoZero, 
    /** 
     * @f[ 
     *   c = \frac{\frac{c_1}{e_1^2}+\frac{c_2}{e_2^2}}{
     *             \frac{1}{e_1^2}+\frac{1}{e_2^2}}
     * @f]
     * @f[
     *   e = \sqrt{\frac{1}{\frac{1}{e_1^2}+\frac{1}{e_2^2}}}
     * @f]
     */
    kWeightedMean, 
    /** 
     * @f[
     *     c = \left\{\begin{array}{cl}
     *          c_1 & \mbox{if $e_1 < e_2$} \\
     *          c_2 & \mbox{otherwise}\end{array}\right.
     * @f]
     */
    kLeastError,
    /** 
     * Just sum the signals 
     */
    kSum
  };
  /**
   * How to obtain the fiducial cuts 
   */
  enum FiducialMethod { 
    /**
     * Select bins by fixed cut.  Bins with a secondary correction
     * less than the cut is considered as non-valid
     */
    kByCut, 
    /**
     * A bin is considered non-valid, if it is less then twice as
     * large as it's neighbors (in eta)
     */
    kDistance 
  };
  /**
   * FMD ring bits for skipping 
   */
   enum FMDRingBits { 
     kFMD1I=0x11,
     kFMD1 =kFMD1I,
     kFMD2I=0x21,
     kFMD2O=0x22,
     kFMD2 =kFMD2I|kFMD2O,
     kFMD3I=0x31,
     kFMD3O=0x32,
     kFMD3 =kFMD2I|kFMD2O
  };
  /** 
   * Constructor 
   */
  AliFMDHistCollector();
  /** 
   * Constructor 
   * 
   * @param title Name of object
   */
  AliFMDHistCollector(const char* title);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDHistCollector(const AliFMDHistCollector& o);

  /** 
   * Destructor 
   */
  virtual ~AliFMDHistCollector();
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   *
   * @return Reference to this object
   */
  AliFMDHistCollector& operator=(const AliFMDHistCollector&);
  /** 
   * Intialise 
   * 
   * @param vtxAxis  @f$ v_z@f$ axis 
   * @param etaAxis  @f$ \eta@f$ axis 
   */  
  virtual void SetupForData(const TAxis& vtxAxis,
		    const TAxis& etaAxis);
  /** 
   * Do the calculations 
   * 
   * @param hists    Cache of histograms 
   * @param sums     Cache to sum ring histograms in 
   * @param vtxBin   Vertex bin (1 based)
   * @param out      Output histogram
   * 
   * @return true on successs 
   */
  virtual Bool_t Collect(const AliForwardUtil::Histos& hists, 
			 AliForwardUtil::Histos&       sums, 
			 UShort_t                      vtxBin, 
			 TH2D&                         out,
			 TList* 			lout=0x0,
			 Double_t 		        cent=-1.0,
			 TList*       sumsv=0x0);
  /** 
   * Output diagnostic histograms to directory 
   * 
   * @param dir List to write in
   */  
  virtual void CreateOutputObjects(TList* dir);
  /** 
   * Set the merge method 
   * 
   * @param m Method
   */
  void SetMergeMethod(MergeMethod m) { fMergeMethod = m; }
  /** 
   * Set the method for finding the fidicual area of the secondary maps 
   * 
   * @param m Method
   */
  void SetFiducialMethod(FiducialMethod m) { fFiducialMethod = m; }
  /** 
   * Set the number of extra bins (beyond the secondary map border) 
   * to cut away. 
   * 
   * @param n Number of bins 
   */
  void SetNCutBins(UInt_t n=2) { fNCutBins = n; }
  /** 
   * Set the correction cut, that is, when bins in the secondary
   * correction maps have a value lower than this cut, they are
   * considered uncertain and not used
   * 
   * @param cut Cut-off 
   */
  void SetCorrectionCut(Float_t cut=0.5) { fCorrectionCut = cut; }
  /** 
   * Set FMD rings to skip. Argument should be
   * kFirstRingToSkip|kSecondRingToSkip...
   * 
   * @param mask bit pattern
   */
  void SetFMDRingsToSkip(UShort_t mask) { fSkipFMDRings = mask; }
 /** 
   * Set whether to make bg maps or not
   * 
   * @param use make them
   */
  void SetMakeBGHitMaps(Bool_t use) { fBgAndHitMaps = use; }
 
  /** 
   * Set the debug level. The higher the value the more output 
   * 
   * @param dbg Debug level 
   */

  void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
protected:
  /** 
   * Get the first and last eta bin to use for a given ring and vertex 
   * 
   * @param d        Detector
   * @param r        Ring 
   * @param vtxBin   Vertex bin (1 based)
   * @param first    On return, the first eta bin to use 
   * @param last     On return, the last eta bin to use 
   */
  virtual void GetFirstAndLast(UShort_t d, Char_t r, UShort_t vtxBin, 
			       Int_t& first, Int_t& last) const;
  /** 
   * Get the first and last eta bin to use for a given ring and vertex 
   * 
   * @param idx      Ring index as given by GetIdx
   * @param vtxBin   Vertex bin (1 based) 
   * @param first    On return, the first eta bin to use 
   * @param last     On return, the last eta bin to use 
   */
  virtual void GetFirstAndLast(Int_t idx, UShort_t vtxBin, 
			       Int_t& first, Int_t& last) const;
  /** 
   * Get the first eta bin to use for a given ring and vertex 
   * 
   * @param d Detector 
   * @param r Ring 
   * @param v vertex bin (1 based)
   * 
   * @return First eta bin to use, or -1 in case of problems 
   */  
  Int_t GetFirst(UShort_t d, Char_t r, UShort_t v) const; 
  /** 
   * Get the first eta bin to use for a given ring and vertex 
   * 
   * @param idx Ring index as given by GetIdx
   * @param v vertex bin (1 based)
   * 
   * @return First eta bin to use, or -1 in case of problems 
   */  
  Int_t GetFirst(Int_t idx, UShort_t v) const; 
  /** 
   * Get the last eta bin to use for a given ring and vertex 
   * 
   * @param d Detector 
   * @param r Ring 
   * @param v vertex bin (1 based)
   * 
   * @return Last eta bin to use, or -1 in case of problems 
   */  
  Int_t GetLast(UShort_t d, Char_t r, UShort_t v) const;
  /** 
   * Get the last eta bin to use for a given ring and vertex 
   * 
   * @param idx Ring index as given by GetIdx
   * @param v vertex bin (1 based)
   * 
   * @return Last eta bin to use, or -1 in case of problems 
   */  
  Int_t GetLast(Int_t idx, UShort_t v) const; 
  /** 
   * Get the detector and ring from the ring index 
   * 
   * @param idx Ring index 
   * @param d   On return, the detector or 0 in case of errors 
   * @param r   On return, the ring id or '0' in case of errors 
   */
  void GetDetRing(Int_t idx, UShort_t& d, Char_t& r) const;
  /** 
   * Get the ring index from detector number and ring identifier 
   * 
   * @param d Detector
   * @param r Ring identifier 
   * 
   * @return ring index or -1 in case of problems 
   */
  Int_t GetIdx(UShort_t d, Char_t r) const;
  /** 
   * Get the possibly overlapping histogram of eta bin @a e in 
   * detector and ring 
   * 
   * @param d Detector
   * @param r Ring 
   * @param e Eta bin
   * @param v Vertex bin (1 based)
   *
   * @return Overlapping histogram index or -1
   */
  Int_t GetOverlap(UShort_t d, Char_t r, Int_t e, UShort_t v) const;
  /** 
   * Get the possibly overlapping histogram of eta bin @a e in 
   * detector and ring 
   * 
   * @param i Ring index
   * @param e Eta bin
   * @param v Vertex bin (1 based)
   *
   * @return Overlapping histogram index or -1
   */
  Int_t GetOverlap(Int_t i, Int_t e, UShort_t v) const;
  /** 
   * Check if there's an overlapping histogram with this eta bin of
   * the detector and ring
   * 
   * @param d Detector 
   * @param r Ring 
   * @param e eta bin
   * @param v Vertex bin (1 based)
   * 
   * @return True if there's an overlapping histogram 
   */
  Bool_t HasOverlap(UShort_t d, Char_t r, Int_t e, UShort_t v) const;
  /** 
   * Check if there's an overlapping histogram with this eta bin of
   * ring
   * 
   * @param i Ring index
   * @param e eta bin
   * @param v Vertex bin
   * 
   * @return True if there's an overlapping histogram 
   */
  Bool_t HasOverlap(Int_t i, Int_t e, UShort_t v) const;
  /** 
   * Check if we should include the bin in the data range 
   * 
   * @param bg Secondary map histogram
   * @param ie Eta bin
   * @param ip Phi bin
   * 
   * @return True if to be used
   */
  Bool_t CheckCorrection(const TH2D* bg, Int_t ie, Int_t ip) const;
  /** 
   * Merge bins accoring to set method
   * 
   * @param c   Current content
   * @param e   Current error
   * @param oc  Old content
   * @param oe  Old error
   * @param rc  On return, the new content
   * @param re  On return, tne new error
   */
  void MergeBins(Double_t c,   Double_t e, 
		 Double_t oc,  Double_t oe,
		 Double_t& rc, Double_t& re) const;
  

  Int_t       fNCutBins;        // Number of additional bins to cut away
  Float_t     fCorrectionCut;   // Cut-off on secondary corrections 
  TArrayI     fFirstBins;       // Array of first eta bins 
  TArrayI     fLastBins;        // Array of last eta bins 
  Int_t       fDebug;           // Debug level 
  TList*      fList;		// Output list
  TH2D*       fSumRings;        // Sum per ring (on y-axis)
  TH2D*       fCoverage;        // Sum per ring (on y-axis)
  MergeMethod fMergeMethod;     // Merge methiod for overlapping bins 
  FiducialMethod fFiducialMethod; // Fidicual method
  UShort_t    fSkipFMDRings;    // FMD rings to ignore     
  Bool_t      fBgAndHitMaps;    // Make hit/bg maps or not
  
  ClassDef(AliFMDHistCollector,4); // Calculate Nch density 
};

//____________________________________________________________________
inline void
AliFMDHistCollector::GetFirstAndLast(UShort_t d, Char_t r, UShort_t vtxbin, 
				     Int_t& first, Int_t& last) const
{
  GetFirstAndLast(GetIdx(d,r), vtxbin, first, last);
}
//____________________________________________________________________
inline Int_t
AliFMDHistCollector::GetFirst(UShort_t d, Char_t r, UShort_t v) const 
{
  return GetFirst(GetIdx(d,r), v);
}
//____________________________________________________________________
inline Int_t
AliFMDHistCollector::GetLast(UShort_t d, Char_t r, UShort_t v) const 
{
  return GetLast(GetIdx(d, r), v);
}
//____________________________________________________________________
inline Bool_t
AliFMDHistCollector::HasOverlap(UShort_t d, Char_t r, Int_t e, UShort_t v) const
{
  return GetOverlap(d,r,e,v) >= 0;
}
//____________________________________________________________________
inline Bool_t
AliFMDHistCollector::HasOverlap(Int_t i, Int_t e, UShort_t v) const
{
  return GetOverlap(i,e,v) >= 0;
}

#endif
// Local Variables:
//   mode: C++
// End:

