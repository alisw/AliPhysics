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
    : AliQACheckerBase("FMD","FMD Quality Assurance Checker") ,
      fDoScale(false)
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
  void Check(Double_t* rv, AliQAv1::ALITASK_t what, 
	     TObjArray** list, const AliDetectorRecoParam* t);
  /** 
   * Make output images.  This is overridden relative to the base
   * class so that we can set the log(y) scale and put everything on
   * the same axis. 
   * 
   * @param list  List of specie array of histograms 
   * @param task  What to show 
   * @param mode  Mode 
   */
  void  MakeImage(TObjArray** list, 
		  AliQAv1::TASKINDEX_t task, 
		  AliQAv1::MODE_t mode);
  void SetDoScale(Bool_t on=true) { fDoScale = on; }

protected:
  /** 
   * Check one histogram 
   * 
   * @param specie 
   * @param hist 
   * 
   * @return 
   */
  Double_t CheckOne(AliQAv1::ALITASK_t          what,
		    AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  Double_t CheckRaw(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  Double_t CheckSim(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  Double_t CheckESD(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  Double_t CheckRec(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;

  Bool_t fDoScale;
  ClassDef(AliFMDQAChecker,0)  // Yves? what to do? 
};

#endif // AliFMDQAChecker_H
// Local Variables:
//  mode: c++
// End:
