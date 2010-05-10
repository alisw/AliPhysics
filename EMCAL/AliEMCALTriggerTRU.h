#ifndef ALIEMCALTRIGGERTRU_H
#define ALIEMCALTRIGGERTRU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <AliEMCALTriggerBoard.h>

class AliEMCALTriggerSTU;
class AliEMCALDigit;
class AliEMCALTriggerTRUDCSConfig;

class AliEMCALTriggerTRU : public AliEMCALTriggerBoard 
{
public:
	
	               AliEMCALTriggerTRU();
 	               AliEMCALTriggerTRU(AliEMCALTriggerTRUDCSConfig* dcsConf, const TVector2& rSize, Int_t mapType);
	virtual       ~AliEMCALTriggerTRU();
	
	virtual Int_t  L0();
	virtual void   SetADC(Int_t channel, Int_t bin, Int_t sig );
	virtual void   SaveRegionADC(Int_t iTRU, Int_t iEvent);
	virtual void   Reset();
	virtual void   ShowFastOR(Int_t timewindow, Int_t chan = -1);
	
private:
	                    AliEMCALTriggerTRU(const AliEMCALTriggerTRU& rhs);
	         AliEMCALTriggerTRU& operator=(const AliEMCALTriggerTRU& rhs);
	
	AliEMCALTriggerTRUDCSConfig* fDCSConfig;

	Int_t         fADC[96][256]; //! FIXME: Check the maximum number of samples
	
	ClassDef(AliEMCALTriggerTRU,1)
};
 
#endif
