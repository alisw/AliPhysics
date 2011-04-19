#ifndef ALIEMCALTRIGGERDATA_H
#define ALIEMCALTRIGGERDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
EMCal trigger data container: can be used independently of the data stream (simulation or raw data)
for transient storage of trigger data
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerTypes.h"

#include <TObject.h>
#include <TVector2.h>
#include <TClonesArray.h>

class AliEMCALTriggerData : public TObject 
{

public:
	             AliEMCALTriggerData();
	virtual     ~AliEMCALTriggerData();

	virtual void SetMode(Int_t i) {fMode = i;}
	
	virtual void SetL0Trigger(      Int_t i, Int_t j, Int_t k) { fL0Trigger[i][j] = k; }
	virtual void SetL0Region(       Int_t i, const Int_t**& region);
	virtual void SetL1Region(       Int_t i,       Int_t**& region);
	
	virtual void SetPatches(TriggerType_t type, Int_t i, const TClonesArray& patches);
	
	virtual void SetL1GammaThreshold(Int_t  v) {fL1GammaThreshold = v;}
	virtual void SetL1JetThreshold(  Int_t  v) {  fL1JetThreshold = v;}
	virtual void SetL1V0(            Int_t* v) {for (int i = 0; i < 2; i++) fL1V0[i] = v[i];}
	virtual void SetL1FrameMask(     Int_t  v) {     fL1FrameMask = v;}            
	virtual void SetL1TriggerType(   Int_t* v) {for (int i = 0; i < 8; i++) fL1TriggerType[i] = v[i];}

	virtual void SetL1DataDecoded(   Int_t  v) {   fL1DataDecoded = v;} 
	
	virtual void          GetL0Trigger(  Int_t i, Int_t j, Int_t& k  ) const {   k = fL0Trigger[i][j];}
	virtual Int_t         GetL0Trigger(  Int_t i, Int_t j            ) const {return fL0Trigger[i][j];}
	
	virtual void          GetPatches(TriggerType_t type, Int_t i, TClonesArray& patches) const;
	virtual TClonesArray* GetPatches(TriggerType_t type, Int_t i                       ) const;

	virtual void          GetL1Region( Int_t i, Int_t arr[][64]             ) const;
	
	virtual Int_t         GetL1GammaThreshold()           const {return fL1GammaThreshold;}
	virtual Int_t         GetL1JetThreshold()             const {return   fL1JetThreshold;}
	virtual void          GetL1V0(            Int_t  v[]) const {for (int i = 0; i < 2; i++) v[i] = fL1V0[i];}
	virtual Int_t         GetL1FrameMask(               ) const {return fL1FrameMask;}            
	virtual void          GetL1TriggerType(   Int_t  v[]) const {for (int i = 0; i < 8; i++) v[i] = fL1TriggerType[i];}
	
	virtual Int_t         GetL1DataDecoded(             ) const {return fL1DataDecoded;}
	
	virtual Int_t         GetMode() const {return fMode;}
	
	virtual void          Scan() const;
	virtual void          Reset();
	
private:

    AliEMCALTriggerData(const AliEMCALTriggerData& rhs);            // NOT implemented
	AliEMCALTriggerData& operator=(const AliEMCALTriggerData& rhs); // NOT implemented
	
	Int_t                      fMode;          //
	
	Int_t               fL0Trigger[2][32];      //
	
	TClonesArray*       fL0Patches[2];          // array of patches  
	
	Int_t               fL0Region[32][24][4];   // from F-ALTRO data only
	
	TClonesArray*  fL1GammaPatches[2];          // array of patches  
	TClonesArray*    fL1JetPatches[2];          // array of patches
	
	Int_t                fL1Region[2][48][64];  // STU FastOR

	Int_t        fL1GammaThreshold;             //
	Int_t          fL1JetThreshold;             //	
	
	Int_t                    fL1V0[2];          //
	Int_t             fL1FrameMask;             //
	Int_t           fL1TriggerType[8];          //
	
	Int_t           fL1DataDecoded;
	
	ClassDef(AliEMCALTriggerData,2)
};

#endif
