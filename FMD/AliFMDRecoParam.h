#ifndef ALIFMDRECOPARAM_H
#define ALIFMDRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//
//
// FMD reconstruction parameters
//
//

#include "AliDetectorRecoParam.h"

class AliFMDRecoParam : public AliDetectorRecoParam
{
public: 
  AliFMDRecoParam(Float_t noiseFactor=3, 
		  Bool_t angleCorrect=kTRUE,
		  Bool_t sharingCorrect=kFALSE);
  virtual ~AliFMDRecoParam() {}
  /** 
   * Whether to do angle of passage correction 
   * 
   * @return @c true if we're to do angle of passage correction
   */
  Bool_t  AngleCorrect()   const { return fAngleCorrect; }
  /** 
   * Get the noise suppression factor
   * 
   * @return The number of noise levels away from the pedestal 
   *         that are suppressed. 
   */
  Float_t NoiseFactor()    const { return fNoiseFactor; }
  /** 
   * Whether to do the sharing correction.  A single particle may
   * traverse more than one strip due to it's incident angle.  In that
   * case, part of it's signal is put in one strip, and another in
   * it's adjacent strip.  The sharing correction corrects for this
   * and adds the signal of the two strips into a single strip. 
   * 
   * @return @c true if the reconstruction should also do the sharing
   * correction. 
   */
  Bool_t  SharingCorrect() const { return fSharingCorrect; }

  void SetAngleCorrect(Bool_t doit) { fAngleCorrect = doit; }
  void SetSharingCorrect(Bool_t doit) { fSharingCorrect = doit; }
  void SetNoiseFactor(Float_t f) { fNoiseFactor = f; }

  static AliFMDRecoParam* GetLowFluxParam();
  static AliFMDRecoParam* GetHighFluxParam();
  static AliFMDRecoParam* GetParam(AliRecoParam::EventSpecie_t specie);
private:
  Float_t fNoiseFactor;    // Noise suppression factor 
  Bool_t  fAngleCorrect;   // Whether to do angle correction or not
  Bool_t  fSharingCorrect; // Whether to do sharing correction or not
  
  ClassDef(AliFMDRecoParam, 1)
};


#endif
// Local Variables:
//  mode: C++ 
// End:
