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
class AliEMCALCalibData;

class AliEMCALTriggerTRU : public AliEMCALTriggerBoard 
{
public:
	
	               AliEMCALTriggerTRU();
 	               AliEMCALTriggerTRU(AliEMCALCalibData *calibData, const TVector2& rSize, Int_t mapType);
	virtual       ~AliEMCALTriggerTRU();
	
	virtual Int_t  L0v0(); // space sum has reached a max
	virtual Int_t  L0v1(); // one of the 4 FastOR in the patch has reached a max
	virtual Int_t  L0v2(); // activity trigger
	
	virtual void   PeakFinder(const Int_t idx[], Int_t nfastor, Int_t start, Int_t nup, Int_t ndown, Int_t& npeaks); 
	virtual void   SetADC(Int_t channel, Int_t bin, Int_t sig );
	
	virtual void   SaveRegionADC(Int_t iTRU, Int_t iEvent);
//	virtual void   Scan();
	virtual void   Reset();
	virtual void   Peaks(Int_t arr[96][2]);
	virtual void   ShowFastOR(Int_t timewindow, Int_t chan = -1);
	
private:
	                    AliEMCALTriggerTRU(const AliEMCALTriggerTRU& rhs);
	         AliEMCALTriggerTRU& operator=(const AliEMCALTriggerTRU& rhs);
	
	Int_t         fADC[96][256]; //! FIXME: Check the maximum number of samples
	
	ClassDef(AliEMCALTriggerTRU,1)
};
 
#endif
