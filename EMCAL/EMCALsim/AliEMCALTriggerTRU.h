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
	virtual Int_t  L0v0(int mask, int pattern);
	virtual Int_t  L0v1(int mask, int pattern);
	virtual void   SetADC(Int_t channel, Int_t bin, Int_t sig );
	virtual void   SaveRegionADC(Int_t iTRU, Int_t iEvent);
	virtual void   Reset();
	virtual void   ShowFastOR(Int_t timewindow, Int_t chan = -1);
	virtual void   GetL0Region(const int time, Int_t arr[][4]);
	virtual Int_t  GetL0Time() const {return fL0Time;}
	
private:
	                    AliEMCALTriggerTRU(const AliEMCALTriggerTRU& rhs);
	         AliEMCALTriggerTRU& operator=(const AliEMCALTriggerTRU& rhs);
	
	AliEMCALTriggerTRUDCSConfig* fDCSConfig; // DCS config

	Int_t         fADC[96][256]; //! FIXME: Check the maximum number of samples
	Int_t         fL0Time;       // Time when the L0 is issued
	
	ClassDef(AliEMCALTriggerTRU,1)
};
 
#endif
