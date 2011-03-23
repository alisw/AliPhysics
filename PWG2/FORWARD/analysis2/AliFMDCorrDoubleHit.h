//
// This class contains the secondary correction and the double hit
// correction used in low-flux events.
//
#ifndef ALIFMDCORRDOUBLEHIT_H
#define ALIFMDCORRDOUBLEHIT_H
/**
 * @file   AliFMDCorrDoubleHit.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:00:50 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_corr
 * 
 */
#include <TObject.h>
#include <TObjArray.h>
class TH1D;

/**
 * This class contains double hit correction used in low-flux events.
 *
 *
 * The double hit correction is given by 
 * @f[
 *   h_{r}(\eta) = \frac{\sum_i N_{i,strips hit}(\eta)}{
 *                       \sum_i N_{i,total hits}(\eta)}
 * @f]
 *
 * where @f$ N_{i,strips hit}(\eta)@f$ is the number of strips in the
 * @f$\eta@f$ bin that had one or more hits in event @f$i@f$, and
 * @f$N_{i,total hits}(\eta)@f$ is the total number hits in the
 * @f$\eta@f$ bin. 
 *
 * These are generated from Monte-Carlo truth information. 
 *
 * @ingroup pwg2_forward_corr
 * 
 */
class AliFMDCorrDoubleHit : public TObject 
{
public:
  /** 
   * Default constructor 
   */
  AliFMDCorrDoubleHit();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDCorrDoubleHit(const AliFMDCorrDoubleHit& o);
  /**
   * Destructor 
   * 
   */
  virtual ~AliFMDCorrDoubleHit();
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliFMDCorrDoubleHit& operator=(const AliFMDCorrDoubleHit& o);
  /** 
   * @{ 
   * @name Get corrections and parameters 
   */
  /** 
   * Get the double hit correction @f$ h_{r}(\eta)@f$ 
   * 
   * @param d Detector number 
   * @param r Ring identifier 
   * 
   * @return @f$ h_{r}(\eta)@f$ 
   */
  TH1D* GetCorrection(UShort_t d, Char_t r) const;
  /* @} */

  /** 
   * @{ 
   * @name Set corrections and parameters 
   */
  /** 
   * Set the double hit correction @f$ h_{r}(\eta)@f$. Note, that the
   * object takes ownership of the passed pointer.
   * 
   * @param d  Detector number (1-3)
   * @param r  Ring identifier (I or O)
   * @param h  @f$ h_{r}(\eta)@f$ 
   * 
   * @return true if operation succeeded 
   */
  Bool_t SetCorrection(UShort_t d, Char_t r, TH1D* h);
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
   * Get the index corresponding to the given ring 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Index (0 based) or negative in case of errors
   */
  Int_t GetRingIndex(UShort_t d, Char_t r) const;

  TObjArray fCorrections; // Array of per-ring double hit corr. 
  ClassDef(AliFMDCorrDoubleHit,1); // 
};
#endif
// Local Variables:
//  mode: C++
// End:
