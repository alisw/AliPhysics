// 
// This class collects the event histograms into single histograms, 
// one for each ring in each vertex bin.  
//
#ifndef ALIFMDHISTCOLLECTOR_H
#define ALIFMDHISTCOLLECTOR_H
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
 * @ingroup pwg2_forward_algo
 */
class AliFMDHistCollector : public TNamed
{
public:
  /** 
   * Constructor 
   */
  AliFMDHistCollector() 
    : fNCutBins(0), fCorrectionCut(0), fFirstBins(), fLastBins(), fDebug(0)
  {}
  /** 
   * Constructor 
   * 
   * @param title Name of object
   */
  AliFMDHistCollector(const char* title)
    : TNamed("fmdHistCollector", title), 
      fNCutBins(2), 
      fCorrectionCut(0.5), 
      fFirstBins(1), 
      fLastBins(1), 
      fDebug(0)
  {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDHistCollector(const AliFMDHistCollector& o)
    : TNamed(o), 
      fNCutBins(o.fNCutBins), fCorrectionCut(o.fCorrectionCut),
      fFirstBins(o.fFirstBins), fLastBins(o.fLastBins), fDebug(o.fDebug) 
  {}

  /** 
   * Destructor 
   */
  virtual ~AliFMDHistCollector() {}
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
   * @param vtxAxis  Vertex axis 
   */  
  virtual void Init(const TAxis& vtxAxis);
  /** 
   * Do the calculations 
   * 
   * @param hists    Cache of histograms 
   * @param vtxBin   Vertex bin (1 based)
   * @param out      Output histogram
   * 
   * @return true on successs 
   */
  virtual Bool_t Collect(AliForwardUtil::Histos& hists, UShort_t vtxBin, 
			 TH2D& out);
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
   * Set the debug level.  The higher the value the more output 
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


  Int_t       fNCutBins;        // Number of additional bins to cut away
  Float_t     fCorrectionCut;   // Cut-off on secondary corrections 
  TArrayI     fFirstBins;       // Array of first eta bins 
  TArrayI     fLastBins;        // Array of last eta bins 
  Int_t       fDebug;           // Debug level 

  ClassDef(AliFMDHistCollector,1); // Calculate Nch density 
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

