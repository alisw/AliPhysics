#ifndef ALIROOT_PWG2_FORWARD_ANALYSIS2_ALIFMDMCDENSITYCALCULATOR_H
#define ALIROOT_PWG2_FORWARD_ANALYSIS2_ALIFMDMCDENSITYCALCULATOR_H
#include "AliFMDDensityCalculator.h"
#include <TList.h>
#include "AliForwardUtil.h"
class AliMCEvent;
class TH2D;
class TH1D;

/** 
 * This class calculates the inclusive charged particle density
 * in each for the 5 FMD rings based on the MC truth.
 *
 * @par Input:
 *   - AliMCEvent  MC truth event infromation
 *
 * @par Output:
 *   - None
 *
 * @par Corrections used: 
 *   - None
 *
 * @ingroup pwg2_forward 
 */
class AliFMDMCDensityCalculator : public AliFMDDensityCalculator
{
public:
  /** 
   * Constructor 
   */
  AliFMDMCDensityCalculator() : AliFMDDensityCalculator() {}
  /** 
   * Constructor 
   * 
   * @param name Name of object
   */
  AliFMDMCDensityCalculator(const char* name) 
   : AliFMDDensityCalculator(name)
  {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDMCDensityCalculator(const AliFMDMCDensityCalculator& o)
   : AliFMDDensityCalculator(o) 
  {}
  /** 
   * Destructor 
   */
  virtual ~AliFMDMCDensityCalculator() {}
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  AliFMDMCDensityCalculator& operator=(const AliFMDMCDensityCalculator& o);
  /** 
   * Do the calculations 
   * 
   * @param fmd      AliESDFMD object (possibly) corrected for sharing
   * @param hists    Histogram cache
   * @param vtxBin   Vertex bin 
   * @param lowFlux  Low flux flag. 
   * 
   * @return true on successs 
   */
  virtual Bool_t Calculate(const AliESDFMD& fmd, 
			   AliForwardUtil::Histos& hists, 
			   UShort_t vtxBin, Bool_t lowFlux);
  /** 
   * Calculate the charged particle density from the MC track references. 
   * 
   * @param event  MC event
   * @param hists  Histograms to fill
   * @param vz     Interaction z coordinate @f$ v_z@f$
   * @param vtxBin bin corresponding to @f$ v_z@f$
   * 
   * @return true on success
   */
  virtual Bool_t CalculateMC(const AliMCEvent& event, 
			     AliForwardUtil::Histos& hists, 
			     Double_t vz,
			     UShort_t vtxBin);			     
protected:

  ClassDef(AliFMDMCDensityCalculator,1); // Calculate Nch density 
};

#endif
// Local Variables:
//   mode: C++
// End:

