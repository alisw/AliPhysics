//
// This class contains the acceptance correction due to dead channels
//
//
#ifndef ALICENTRALCORRACCEPTANCE_H
#define ALICENTRALCORRACCEPTANCE_H
/**
 * @file   AliCentralCorrAcceptance.h
 * @author Hans Hjersing Dalsgaard 
 * @date   Wed Mar 23 13:58:33 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_corr
 */
#include <TObject.h>
#include <TObjArray.h>
#include <TAxis.h>
class TH1D;

/**
 * This class contains the acceptance correction due to dead channels
 *
 * These are generated from the on-line dead channel calculations 
 *
 * @ingroup pwg2_forward_corr
 */
class AliCentralCorrAcceptance : public TObject 
{
public:
  /** 
   * Default constructor 
   */
  AliCentralCorrAcceptance();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliCentralCorrAcceptance(const AliCentralCorrAcceptance& o);
  /**
   * Destructor 
   * 
   */
  virtual ~AliCentralCorrAcceptance();
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
  AliCentralCorrAcceptance& operator=(const AliCentralCorrAcceptance& o);
  /** 
   * Get the acceptance correction @f$ a_{r,v}@f$ 
   *
   * @param v  Primary interaction point @f$z@f$ coordinate
   * 
   * @return The correction @f$ a_{r,v}@f$ 
   */
  TH1D* GetCorrection(Double_t v) const;
  /** 
   * Get the acceptance correction @f$ a_{r,v}@f$ 
   *
   * @param b  Bin corresponding to the primary interaction point 
   *           @f$z@f$ coordinate (1 based)
   * 
   * @return The correction @f$ a_{r,v}@f$ 
   */
  TH1D* GetCorrection(UShort_t b) const;
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
   * Set the acceptance correction @f$ a_{r,v}(\eta)@f$.
   * Note, that the object takes ownership of the passed pointer.
   * 
   * @param v    Primary interaction point @f$z@f$ coordinate  
   * @param h    @f$ a_{r,v}(\eta)@f$ 
   * 
   * @return true if operation succeeded 
   */
  Bool_t SetCorrection(Double_t v, TH1D* h);
  /** 
   * Set the acceptance correction @f$ a_{r,v}(\eta)@f$ 
   * Note, that the object takes ownership of the passed pointer.
   * 
   * @param b    Bin corresponding to the primary interaction point 
   *             @f$z@f$ coordinate  (1 based)
   * @param h    @f$ a_{r,v}(\eta)@f$ 
   * 
   * @return true if operation succeeded 
   */
  Bool_t SetCorrection(UShort_t b, TH1D* h);
  /** 
   * Set the vertex axis to use 
   * 
   * @param axis Vertex axis 
   */
  void SetVertexAxis(const TAxis& axis);
  /** 
   * Set the vertex axis to use 
   * 
   * @param nBins   Number of bins
   * @param min     Minimum
   * @param max     Maximum
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

  TObjArray fArray;          // Array of per-vertex acceptance corr
  TAxis     fVertexAxis;     // The vertex axis 
  ClassDef(AliCentralCorrAcceptance,1); // Acceptance correction due to dead areas
};

//____________________________________________________________________
inline void 
AliCentralCorrAcceptance::SetVertexAxis(Int_t nBins, Double_t min, 
					   Double_t max)
{
  fVertexAxis.Set(nBins, min, max);
}
//____________________________________________________________________
inline void 
AliCentralCorrAcceptance::SetVertexAxis(const TAxis& e)
{
  fVertexAxis.Set(e.GetNbins(), e.GetXmin(), e.GetXmax());
}
#endif
// Local Variables:
//  mode: C++
// End:
