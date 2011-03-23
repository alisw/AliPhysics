//
// See implementation or Doxygen comments for more information
//
#ifndef ALIAODCENTRALMULT_H
#define ALIAODCENTRALMULT_H
/**
 * @file   AliAODCentralMult.h
 * @author Hans Hjersing Dalsgaard
 * @date   Wed Mar 23 13:57:24 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_aod
 * 
 */
#include <TObject.h>
#include <TH2D.h>
class TBrowser;
/**
 * Class that contains the central multiplicity data per event 
 *
 * This class contains a histogram of 
 * @f[
 *   \frac{d^2N_{ch}}{d\eta d\phi}\quad,
 * @f]
 * as well as a trigger mask for each analysed event.  
 * 
 * The eta acceptance of the event is stored in the underflow bins of
 * the histogram.  So to build the final histogram, one needs to
 * correct for this acceptance (properly weighted by the events), and
 * the vertex efficiency.  This simply boils down to defining a 2D
 * histogram and summing the event histograms in that histogram.  One
 * should of course also do proper book-keeping of the accepted event.
 *
 * @ingroup pwg2_forward
 * @ingroup pwg2_forward_aod
 */
class AliAODCentralMult : public TObject
{
public:
  
  /** 
   * Default constructor 
   * 
   * Used by ROOT I/O sub-system - do not use
   */
  AliAODCentralMult();
  /** 
   * Constructor 
   * 
   * @param isMC Whether this was from MC or not 
   */
  AliAODCentralMult(Bool_t isMC);
  /** 
   * Destructor 
   */
  virtual ~AliAODCentralMult() {} // Destructor 
  /** 
   * Initialize 
   * 
   * @param etaAxis  Pseudo-rapidity axis
   */
  void Init(const TAxis& etaAxis);
  /** 
   * Get the @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   *
   * @return @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   */  
  const TH2D& GetHistogram() const { return fHist; } // Get histogram 
  /** 
   * Get the @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   *
   * @return @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   */  
  TH2D& GetHistogram() { return fHist; } // Get histogram 
  /**
   * Clear Object between events
   *
   * @param opt option (not used)
   *
   */
  void Clear(Option_t* opt);
  /** 
   * browse this object 
   * 
   * @param b Browser 
   */
  void Browse(TBrowser* b);
  /** 
   * This is a folder 
   * 
   * @return Always true
   */
  Bool_t IsFolder() const { return kTRUE; } // Always true 
  /** 
   * Print content 
   * 
   * @param option Passed verbatim to TH2::Print 
   */
  void Print(Option_t* option="") const;

  /** 
   * Get the name of the object 
   * 
   * @return Name of object 
   */
  const Char_t* GetName() const { return (fIsMC ? "CentralClustersMC" : "CentralClusters"); }

protected: 
  Bool_t  fIsMC;     // Whether this is from MC 
  TH2D    fHist;     // Histogram of d^2N_{ch}/(deta dphi) for this event
  
  ClassDef(AliAODCentralMult,1); // AOD forward multiplicity 
};

#endif
// Local Variables:
//  mode: C++
// End:

