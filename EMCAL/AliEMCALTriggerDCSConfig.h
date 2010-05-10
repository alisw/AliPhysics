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
	void                         SetSTUObj(AliEMCALTriggerSTUDCSConfig* so)    { fSTUObj    = so; }
  
	TClonesArray*                GetTRUArr()                 const             { return fTRUArr;  }
	AliEMCALTriggerSTUDCSConfig* GetSTUDCSConfig(          ) const             { return (AliEMCALTriggerSTUDCSConfig*)fSTUObj;           }
	AliEMCALTriggerTRUDCSConfig* GetTRUDCSConfig(Int_t iTRU) const             { return (AliEMCALTriggerTRUDCSConfig*)fTRUArr->At(iTRU); }
	
private:

	AliEMCALTriggerDCSConfig(const AliEMCALTriggerDCSConfig &cd);            // Not implemented
	AliEMCALTriggerDCSConfig &operator=(const AliEMCALTriggerDCSConfig &cd); // Not implemented

	TClonesArray*                fTRUArr; //
	AliEMCALTriggerSTUDCSConfig* fSTUObj; //

	ClassDef(AliEMCALTriggerDCSConfig,1)  //
};
#endif

