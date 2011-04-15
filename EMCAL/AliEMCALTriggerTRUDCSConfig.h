#ifndef ALIEMCALTRIGGERTRUDCSCONFIG_H
#define ALIEMCALTRIGGERTRUDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
 
 
 
 Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
 Author: Jiri Kral, JYU
*/

#include "TObject.h"

class AliEMCALTriggerTRUDCSConfig : public TObject 
{
public:

	AliEMCALTriggerTRUDCSConfig();
	virtual ~AliEMCALTriggerTRUDCSConfig() {}

	void    SetSELPF( UInt_t pf)              { fSELPF  = pf; }  
	void    SetL0SEL( UInt_t la)              { fL0SEL  = la; }  
	void    SetL0COSM(UInt_t lc)              { fL0COSM = lc; }
	void    SetGTHRL0(UInt_t lg)              { fGTHRL0 = lg; }
	void    SetMaskReg(UInt_t msk, Int_t pos) { fMaskReg[pos] = msk; }

	UInt_t   GetSELPF()                       const { return fSELPF;  }
	UInt_t   GetL0SEL()                       const { return fL0SEL;  }
	UInt_t   GetL0COSM()                      const { return fL0COSM; }
	UInt_t   GetGTHRL0()                      const { return fGTHRL0; }
	UInt_t   GetMaskReg(Int_t pos) const { return fMaskReg[pos]; }
	
protected:

	AliEMCALTriggerTRUDCSConfig(const AliEMCALTriggerTRUDCSConfig &cd);
	AliEMCALTriggerTRUDCSConfig &operator=(const AliEMCALTriggerTRUDCSConfig &cd);

private:
	
	UInt_t   fSELPF;                         // PeakFinder setup
	UInt_t   fL0SEL;                         // L0 Algo selection
	UInt_t   fL0COSM;                        // 2x2
	UInt_t   fGTHRL0;                        // 4x4
	UInt_t   fMaskReg[6];                    // 6*16 = 96 mask bits per TRU
	
	ClassDef(AliEMCALTriggerTRUDCSConfig,2)  //
};
#endif
