#ifndef ALIEMCALTRIGGERELECTRONICS_H
#define ALIEMCALTRIGGERELECTRONICS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
EMCal trigger electronics manager L0/L1
can handle both simulated digits and raw data
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#include "TClonesArray.h"

class AliRawReader;
class AliEMCALTriggerDCSConfig;
class TTree;
class AliEMCALTriggerData;
class AliEMCALTriggerSTU;
class AliESDVZERO;
class AliEMCALTriggerTRU;

class AliEMCALTriggerElectronics : public TObject 
{
public:
			       AliEMCALTriggerElectronics(const AliEMCALTriggerDCSConfig* dcsConfig = 0x0); // ctor
	virtual       ~AliEMCALTriggerElectronics();                                   // dtor
	
	virtual void   Digits2Trigger(const TClonesArray* digits, const Int_t V0M[], AliEMCALTriggerData* data);	
	virtual void   Reset();  
	
	virtual AliEMCALTriggerTRU* GetTRU( Int_t iTRU ) {return (AliEMCALTriggerTRU*)fTRU->At(iTRU);}
	virtual AliEMCALTriggerSTU* GetSTU(            ) {return                      fSTU;          }
	
private:

	AliEMCALTriggerElectronics(const AliEMCALTriggerElectronics& other);            // Not implemented
	AliEMCALTriggerElectronics& operator=(const AliEMCALTriggerElectronics& other); // Not implemented

	TClonesArray*        fTRU; // 32 TRU
	AliEMCALTriggerSTU*  fSTU; //  1 STU
		
  ClassDef(AliEMCALTriggerElectronics,1)
};

#endif
