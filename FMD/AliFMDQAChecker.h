#ifndef ALIFMDQACHECKER_H
#define ALIFMDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#include "AliQACheckerBase.h"
class TFile; 
class TH1F; 
class TH1I; 

/** @class AliFMDQAChecker 
    @brief Quality assurance checker for the FMD */
class AliFMDQAChecker : public AliQACheckerBase 
{
public:
  /** Constructor */
  AliFMDQAChecker() 
    : AliQACheckerBase("FMD","FMD Quality Assurance Checker") 
  {}          
  /** Destructor */
  virtual ~AliFMDQAChecker() {}
  /** 
   * Member function called to do the actual checking
   * 
   * @param rv   Array of return values. 
   * @param what What to check 
   * @param list Array of arrays of histograms.  There's one arrat for
   *             each 'specie'
   * @param t    Reconstruction parameters - not used. 
   */
  void Check(Double_t* rv, AliQAv1::ALITASK_t what, TObjArray** list, const AliDetectorRecoParam* t);
private:
  ClassDef(AliFMDQAChecker,0)  // Yves? what to do? 
};

#endif // AliFMDQAChecker_H
// Local Variables:
//  mode: c++
// End:
