#ifndef ALIEMCALTRIGGERDCSCONFIG_H
#define ALIEMCALTRIGGERDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 


Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "TObject.h"
#include "TClonesArray.h"

class AliEMCALTriggerSTUDCSConfig;
class AliEMCALTriggerTRUDCSConfig;

class AliEMCALTriggerDCSConfig : public TObject 
{
public:
	
	         AliEMCALTriggerDCSConfig();
	virtual ~AliEMCALTriggerDCSConfig();
	
	void                         SetTRUArr(TClonesArray* const ta)             { fTRUArr    = ta; }
	inline void                  SetSTUObj(AliEMCALTriggerSTUDCSConfig* so, Bool_t isDCAL = false);
  
	TClonesArray*                GetTRUArr()                 const             { return fTRUArr;  }
	inline AliEMCALTriggerSTUDCSConfig* GetSTUDCSConfig(Bool_t isDCAL = false) const;
	AliEMCALTriggerTRUDCSConfig* GetTRUDCSConfig(Int_t iTRU) const             { return (AliEMCALTriggerTRUDCSConfig*)fTRUArr->At(iTRU); }
	
private:

	AliEMCALTriggerDCSConfig(const AliEMCALTriggerDCSConfig &cd);            // Not implemented
	AliEMCALTriggerDCSConfig &operator=(const AliEMCALTriggerDCSConfig &cd); // Not implemented

	TClonesArray*                fTRUArr;   // TRU array
	AliEMCALTriggerSTUDCSConfig* fSTUObj;   // STU
	AliEMCALTriggerSTUDCSConfig* fSTUDCAL;  // STU of DCAL

	ClassDef(AliEMCALTriggerDCSConfig,2)  //
};

void AliEMCALTriggerDCSConfig::SetSTUObj(AliEMCALTriggerSTUDCSConfig* so, Bool_t isDCAL) {
	if(isDCAL)
	  fSTUDCAL = so;
	else
	  fSTUObj  = so;
}

AliEMCALTriggerSTUDCSConfig* AliEMCALTriggerDCSConfig::GetSTUDCSConfig(Bool_t isDCAL) const{
  if(isDCAL)
    return fSTUDCAL;
  return fSTUObj;
}

#endif

