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
class TH2;
class TH2D;
class TH1D;
class TObjArray;

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
    kSum,
    /** 
     * In overlaps, prefer inners, or if both are inners, 
     * do the straight mean 
     */
    kPreferInner,
    /** 
     * In overlaps, prefer outers, or if both are outers (doesn't happen), 
     * do the straight mean 
     */
    kPreferOuter

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
     kFMD1I=0x01,
     kFMD1 =kFMD1I,
     kFMD2I=0x02,
     kFMD2O=0x04,
     kFMD2 =kFMD2I|kFMD2O,
     kFMD3I=0x08,
     kFMD3O=0x10,
     kFMD3 =kFMD3I|kFMD3O
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
   * @param cent     Centrality
   * @param eta2phi  Copy eta coverage to phi acceptance 
   * @param add      If true, add to internal caches
   * 
   * @return true on successs 
   */
  virtual Bool_t Collect(const AliForwardUtil::Histos& hists, 
			 AliForwardUtil::Histos&       sums, 
			 UShort_t                      vtxBin, 
			 TH2D&                         out,
			 Double_t 		       cent=-1.0,
			 Bool_t                        eta2phi=false,
			 Bool_t                        add=true);
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
  MergeMethod GetMergeMethod() const { return fMergeMethod; }
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
   * Set whether to make by-centrality sums for each ring
   * 
   * @param use If true, make by-centrality sums
   */
  void SetMakeCentralitySums(Bool_t use) { fDoByCent = use; }
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
   * Get the detector and ring from the ring index 
   * 
   * @param idx Ring index 
   * @param d   On return, the detector or 0 in case of errors 
   * @param r   On return, the ring id or '0' in case of errors 
   */
  static void GetDetRing(Int_t idx, UShort_t& d, Char_t& r);
  /** 
   * Get the ring index from detector number and ring identifier 
   * 
   * @param d Detector
   * @param r Ring identifier 
   * 
   * @return ring index or -1 in case of problems 
   */
  static Int_t GetIdx(UShort_t d, Char_t r);
  /** 
   * Check if the detector @a d, ring @a r is listed <i>in</i> the @a
   * skips bit mask.  If the detector/ring is in the mask, return true.
   * 
   * That is, use case is 
   * @code 
   *  for (UShort_t d=1. d<=3, d++) {
   *    UShort_t nr = (d == 1 ? 1 : 2);
   *    for (UShort_t q = 0; q < nr; q++) { 
   *      Char_t r = (q == 0 ? 'I' : 'O');
   *      if (CheckSkips(d, r, skips)) continue; 
   *      // Process detector/ring 
   *    }
   *  }
   * @endcode
   *
   * @param d      Detector
   * @param r      Ring 
   * @param skips  Mask of detector/rings to skip
   * 
   * @return True if detector @a d, ring @a r is in the mask @a skips 
   */
  static Bool_t CheckSkip(UShort_t d, Char_t r, UShort_t skips);
  /** 
   * Check the correction
   * 
   * @param m   Fiducial method used
   * @param cut Cut value 
   * @param bg  Secondary map
   * @param ie  @f$\eta@f$ bin
   * @param ip  @f$\varphi@f$ bin
   * 
   * @return true if OK. 
   */
  static Bool_t CheckCorrection(FiducialMethod m, Double_t cut, 
				const TH2D* bg, Int_t ie, Int_t ip);

  /** 
   * Merge bins accoring to set method
   * 
   * @param m   Merging method
   * @param c   Current content
   * @param e   Current error
   * @param oc  Old content
   * @param oe  Old error
   * @param rc  On return, the new content
   * @param re  On return, tne new error
   */
  static void MergeBins(MergeMethod   m, 
			Double_t c,   Double_t e, 
			Double_t oc,  Double_t oe,
			Double_t& rc, Double_t& re);
  
  //==================================================================
  /**
   * Structure to hold per-vertex bin cache of per-ring histograms 
   */
  struct VtxBin : public TObject
  {
    /** 
     * Constructor 
     * 
     * @param index   Index number
     * @param minIpZ  Least @f$IP_{z}@f$
     * @param maxIpZ  Largest @f$IP_{z}@f$
     * @param nCut    Cut on n
     */
    VtxBin(Int_t index=0, Double_t minIpZ=999, Double_t maxIpZ=-999,
	   Int_t nCut=0);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from
     */
    VtxBin(const VtxBin& o);
    /** 
     * Assignment operator
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object
     */    
    VtxBin& operator=(const VtxBin& o);
    /** 
     * Override to give name based on cuts
     * 
     * @return Name
     */
    const Char_t* GetName() const;
    /** 
     * Set up for data
     * 
     * @param coverage    Diagnostics histogram to be filled 
     * @param skip        Skip flags
     * @param fiducial    Fiducial cut method
     * @param cut         Fiducial cut
     * @param l           Parent output list 
     * @param etaAxis     @f$\eta@f$ axis used
     * @param doHitMap    If true, also do a per-ring sum
     * @param storeSecMap If true, store used secondary map
     */
    void SetupForData(TH2*           coverage,
		      UShort_t       skip,
		      FiducialMethod fiducial, 
		      Double_t       cut,
		      TList*         l, 
		      const TAxis&   etaAxis,
		      Bool_t         doHitMap,
		      Bool_t         storeSecMap);
    /** 
     * Process one event in this vertex bin
     * 
     * @param hists      Histograms
     * @param sums       Sum histograms
     * @param out        Per-event output histogram
     * @param sumRings   Sum per ring 
     * @param skipped    Histogram of skipped rings 
     * @param cent       Event centrality
     * @param m          Merging method
     * @param skips      Which rings to skip
     * @param byCent     List (or null) of per centrality sums
     * @param eta2phi    Copy eta coverage to phi acceptance 
     * @param add      If true, add to internal caches
     *
     * @return true on success
     */
    Bool_t Collect(const AliForwardUtil::Histos& hists, 
		   AliForwardUtil::Histos&       sums, 
		   TH2D&                         out,
		   TH2D*                         sumRings,
		   TH1D*                         skipped,
		   Double_t                      cent,
		   MergeMethod                   m,
		   UShort_t                      skips,
		   TList*                        byCent,
		   Bool_t                        eta2phi,
		   Bool_t                        add);
    /** 
     * Check if there's an overlap between detector @a d, ring @a r
     * and some other ring for the given @f$\eta@f$ @a bin.  If so,
     * return the ring index.  If not, return -1.
     * 
     * @param d    Current detector
     * @param r    Current ring
     * @param bin  Current @f$\eta@f$ bin
     * 
     * @return Index of overlapping ring, or -1
     */    
    Int_t GetOverlap(UShort_t d, Char_t r, Int_t bin) const;
    /** 
     * Get the first and last @f$\eta@f$ bin for a detector 
     * 
     * @param d      Current detector 
     * @param r      Current ring	   
     * @param first  On return, the first @f$\eta@f$ bin
     * @param last   On return, the last @f$\eta@f$ bin
     */
    void GetFirstAndLast(UShort_t d, UShort_t r, 
			 Int_t& first, Int_t& last) const {
      GetFirstAndLast(GetIdx(d,r), first, last);
    }
    /** 
     * Get the first and last @f$\eta@f$ bin for a detector 
     * 
     * @param idx    Current ring index
     * @param first  On return, the first @f$\eta@f$ bin
     * @param last   On return, the last @f$\eta@f$ bin
     */
    void GetFirstAndLast(Int_t idx,Int_t& first, Int_t& last) const;
    /** 
     * Get the first @f$\eta@f$ bin
     * 
     * @param idx Ring index (0-4)
     * 
     * @return bin number
     */
    Int_t GetFirst(Int_t idx) const;
    /** 
     * Get the last @f$\eta@f$ bin
     * 
     * @param idx Ring index (0-4)
     * 
     * @return bin number
     */
    Int_t GetLast(Int_t idx) const;
    /** 
     * Get the first @f$\eta@f$ bin
     * 
     * @param d  Detector
     * @param r  Ring
     * 
     * @return bin number
     */
    Int_t GetFirst(UShort_t d, Char_t r) const { return GetFirst(GetIdx(d,r));}
    /** 
     * Get the last @f$\eta@f$ bin
     * 
     * @param d  Detector
     * @param r  Ring
     * 
     * @return bin number
     */
    Int_t GetLast(UShort_t d, Char_t r) const { return GetLast(GetIdx(d,r));}

    Int_t                   fIndex;     // Vertex bin index
    Double_t                fLow;       // Low @f$ ip_z @f$ 
    Double_t                fHigh;      // High @f$ ip_z @f$
    AliForwardUtil::Histos* fHitMap;    // Hit map (optional)
    TArrayI                 fFirstBin;  // Per-ring first bin
    TArrayI                 fLastBin;   // Per-ring last bin
    Int_t                   fNCutBins;  // Number of bins to cut 

    // ClassDef(VtxBin,1); // Vertex bin in histogram collector
  };
  /** 
   * Get a vertex bin
   * 
   * @param ivtx Bin number (1-nVz)
   * 
   * @return Bin or null
   */
  VtxBin* GetVtxBin(Int_t ivtx);
  /** 
   * Get a vertex bin
   * 
   * @param ivtx Bin number (1-nVz)
   * 
   * @return Bin or null
   */
  const VtxBin* GetVtxBin(Int_t ivtx) const;

  Int_t       fNCutBins;        // Number of additional bins to cut away
  Float_t     fCorrectionCut;   // Cut-off on secondary corrections 
  Int_t       fDebug;           // Debug level 
  TList*      fList;		// Output list
  TH2D*       fSumRings;        // Sum per ring (on y-axis)
  TH2D*       fCoverage;        // Sum per ring (on y-axis)
  TH1D*       fSkipped;         // Skipped rings
  MergeMethod fMergeMethod;     // Merge methiod for overlapping bins 
  FiducialMethod fFiducialMethod; // Fidicual method
  UShort_t    fSkipFMDRings;    // FMD rings to ignore     
  Bool_t      fBgAndHitMaps;    // Make hit/bg maps or not
  TObjArray*  fVtxList;         //! Per-vertex list
  TList*      fByCent;          // By centrality sums
  Bool_t      fDoByCent;        // Whether to do by centrality sum
  ClassDef(AliFMDHistCollector,6); // Calculate Nch density 
};


#endif
// Local Variables:
//   mode: C++
// End:

