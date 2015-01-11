#ifndef ALIFMDOFFLINETRIGGER_H
#define ALIFMDOFFLINETRIGGER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
// 
// This class implements the FMD offline trigger as requested for 
// ALICE first physics.  
//
//


#include "TObject.h"

class AliESDFMD;
//____________________________________________________________________
/** @brief FMD offline trigger class
    @ingroup FMD_rec
*/
class AliFMDOfflineTrigger : public TObject 
{
public:
  /** CTOR */
  AliFMDOfflineTrigger();
  /** DTOR */
  ~AliFMDOfflineTrigger() {}
  /** Copy ctor 
      @param o Object to copy from  */
  AliFMDOfflineTrigger(const AliFMDOfflineTrigger& o);
  
  AliFMDOfflineTrigger& operator=(const AliFMDOfflineTrigger& o);
    
  Bool_t ASideHasHit(AliESDFMD* fmd);
  Bool_t CSideHasHit(AliESDFMD* fmd);
  void   SetLowCut(Float_t lowcut) {fLowCut = lowcut;}
  void   SetHitCut(Float_t hitcut) {fHitCut = hitcut;}
private:
  
  Float_t fLowCut;
  Float_t fHitCut;
  
  ClassDef(AliFMDOfflineTrigger, 1) 
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


