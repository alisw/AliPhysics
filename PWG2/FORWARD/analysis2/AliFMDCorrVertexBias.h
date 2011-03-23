//
// This class contains the secondary correction and the double hit
// correction used in low-flux events.
//
#ifndef ALIFMDCORRVERTEXBIAS_H
#define ALIFMDCORRVERTEXBIAS_H
/**
 * @file   AliFMDCorrVertexBias.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:01:56 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_corr
 */
#include <TObject.h>
#include <TObjArray.h>
#include <TAxis.h>
class TH2D;

/**
 * This class contains the correction for the bias introduced by
 * different vertex bins
 *
 * The correction is given by 
 * @f[
 *   b_{v}(\eta,\varphi) = \frac{1/N_{t}\sum_i^{N_{tv}} N_{ch,i,primary}}{
 *                               1/N_{v}\sum_i^{N_{v}} N_{ch,i,primary}}
 * @f]
 *
 * where @f$N_{ch,i,primary}@f$ is the number of primary particles in
 * the given @f$(\eta,\varphi)@f$, and where the denominator sum runs
 * over all events with a vertex within the given vertex bin, and the
 * sum in the numerator runs over only events that have a valid
 * trigger and reconstructed vertex.  @f$ N_{t}@f$ is the number of
 * events with a valid trigger (but not necessarily a valid vertex).
 * The vertex information used @f$v@f$ is in all cases the MC truth
 * vertex
 *
 * These are generated from Monte-Carlo truth and ESD information. 
 *
 * @ingroup pwg2_forward_corr
 */
class AliFMDCorrVertexBias : public TObject 
{
public:
  /** 
   * Default constructor 
   */
  AliFMDCorrVertexBias();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDCorrVertexBias(const AliFMDCorrVertexBias& o);
  /**
   * Destructor 
   * 
   */
  virtual ~AliFMDCorrVertexBias();
  /** 
   * @{ 
   * @name Get corrections and parameters 
   */
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliFMDCorrVertexBias& operator=(const AliFMDCorrVertexBias& o);
  /** 
   * Get the vertex bias correction @f$ b_{v}@f$ 
   *
   * @param r  Ring type (I or O)
   * @param v  Primary interaction point @f$z@f$ coordinate
   * 
   * @return The correction @f$ b_{v}@f$ 
   */
  TH2D* GetCorrection(Char_t r, Double_t v) const;
  /** 
   * Get the vertex bias correction @f$ b_{v}@f$ 
   *
   * @param r  Ring type (I or O)
   * @param b  Bin corresponding to the primary interaction point 
   *           @f$z@f$ coordinate (1 based)
   * 
   * @return The correction @f$ b_{v}@f$ 
   */
  TH2D* GetCorrection(Char_t r, UShort_t b) const;
  /** 
   * Get the vertex axis used 
   * 
   * @return vertex axis 
   */
  const TAxis& GetVertexAxis() const { return fVertexAxis; }
  /* @} */

  /** 
   * @{ 
   * @name Set corrections and parameters 
   */
  /** 
   * Set the vertex bias correction @f$ b_{v}(\eta,\varphi)@f$.
   * Note, that the object takes ownership of the passed pointer.
   * 
   * @param r    Ring type (I or O)
   * @param v    Primary interaction point @f$z@f$ coordinate  
   * @param h    @f$ b_{v}(\eta,\varphi)@f$ 
   * 
   * @return true if operation succeeded 
   */
  Bool_t SetCorrection(Char_t r, Double_t v, TH2D* h);
  /** 
   * Set the vertex bias correction @f$ b_{v}(\eta,\varphi)@f$ 
   * Note, that the object takes ownership of the passed pointer.
   * 
   * @param r    Ring type (I or O)
   * @param b    Bin corresponding to the primary interaction point 
   *             @f$z@f$ coordinate  (1 based)
   * @param h    @f$ b_{v}(\eta,\varphi)@f$ 
   * 
   * @return true if operation succeeded 
   */
  Bool_t SetCorrection(Char_t r, UShort_t b, TH2D* h);
  /** 
   * Set the vertex axis to use 
   * 
   * @param axis Vertex axis 
   */
  void SetVertexAxis(const TAxis& axis);
  /** 
   * Set the vertex axis to use 
   * 
   * @param nBins Number of bins
   * @param min   Minimum	  
   * @param max   Maximum	  
   */
  void SetVertexAxis(Int_t nBins, Double_t min, Double_t max);
  /* @} */

  /** 
   * @{ 
   * @name Auxiliary member functions 
   */
  /** 
   * Declare this as a folder
   * 
   * @return Always true 
   */
  Bool_t IsFolder() const { return true; }
  /** 
   * Browse this object in the browser
   * 
   * @param b 
   */
  void Browse(TBrowser* b);
  /** 
   * Print this object 
   * 
   * @param option 
   */  
  void Print(Option_t* option="R") const; //*MENU*
  /* @} */
protected:
  /** 
   * Find the vertex bin that corresponds to the passed vertex 
   * 
   * @param vertex The interaction points @f$z@f$-coordinate 
   * 
   * @return Vertex bin in @f$[1,N_{\mbox{vertex}}]@f$ or negative if 
   * out of range 
   */
  Int_t FindVertexBin(Double_t vertex) const;
  /** 
   * Get the vertex array corresponding to the specified ring
   * 
   * @param v vertex bin (1 based)
   * 
   * @return Pointer to vertex array, or null in case of problems
   */
  TObjArray* GetVertexArray(UShort_t v) const;
  /** 
   * Get the vertex array corresponding to the specified ring
   * 
   * @param v vertex bin (1 based)
   * 
   * @return Pointer to vertex array, or newly created container 
   */
  TObjArray* GetOrMakeVertexArray(UShort_t v);

  TObjArray fVertexArray;    // Array of per-ring, per-vertex 2nd map
  TAxis     fVertexAxis;     // The vertex axis 
  ClassDef(AliFMDCorrVertexBias,1); // 
};

//____________________________________________________________________
inline void 
AliFMDCorrVertexBias::SetVertexAxis(Int_t nBins, Double_t min, Double_t max)
{
  fVertexAxis.Set(nBins, min, max);
}
//____________________________________________________________________
inline void 
AliFMDCorrVertexBias::SetVertexAxis(const TAxis& e)
{
  fVertexAxis.Set(e.GetNbins(), e.GetXmin(), e.GetXmax());
}
#endif
// Local Variables:
//  mode: C++
// End:
