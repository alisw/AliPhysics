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
	
	virtual void SetL1GammaThreshold(int i, int  v) {fL1GammaThreshold[i] = v;}
	virtual void SetL1JetThreshold(  int i, int  v) {  fL1JetThreshold[i] = v;}
	virtual void SetL1V0(            Int_t* v) {for (int i = 0; i < 2; i++) fL1V0[i] = v[i];}
	virtual void SetL1FrameMask(     Int_t  v) {     fL1FrameMask = v;}            
	virtual void SetL1TriggerType(   Int_t* v) {for (int i = 0; i < 15; i++) fL1TriggerType[i] = v[i];}

	virtual void SetL1DataDecoded(   Int_t  v) {   fL1DataDecoded = v;} 
	virtual void SetL1RawData(       Int_t  v) {       fL1RawData = v;}
        virtual void SetMedian(          Int_t  v) {          fMedian = v;}

	virtual void          GetL0Trigger(  Int_t i, Int_t j, Int_t& k  ) const {   k = fL0Trigger[i][j];}
	virtual Int_t         GetL0Trigger(  Int_t i, Int_t j            ) const {return fL0Trigger[i][j];}
	
	virtual Int_t         GetL1GammaThreshold(int i)           const {return fL1GammaThreshold[i];}
	virtual Int_t         GetL1JetThreshold(  int i)           const {return   fL1JetThreshold[i];}
	virtual void          GetL1V0(            Int_t  v[]) const {for (int i = 0; i < 2; i++) v[i] = fL1V0[i];}
	virtual Int_t         GetL1FrameMask(               ) const {return fL1FrameMask;}            
	virtual void          GetL1TriggerType(   Int_t  v[]) const {for (int i = 0; i < 15; i++) v[i] = fL1TriggerType[i];}
	
	virtual Int_t         GetL1DataDecoded(             ) const {return fL1DataDecoded;}
        virtual Int_t             GetL1RawData(             ) const {return     fL1RawData;}
	virtual Int_t                GetMedian(             ) const {return        fMedian;}

	virtual Int_t         GetMode() const {return fMode;}
	
	virtual void          Scan() const;
	virtual void          Reset();
	
private:

    AliEMCALTriggerData(const AliEMCALTriggerData& rhs);            // NOT implemented
	AliEMCALTriggerData& operator=(const AliEMCALTriggerData& rhs); // NOT implemented
	
	Int_t                    fMode;             // Simulation/Raw
	
	Int_t               fL0Trigger[2][52];      // Triggering TRU

	Int_t        fL1GammaThreshold[2];          // L1-g threshold
	Int_t          fL1JetThreshold[2];          // L1-j threshold
	
	Int_t                    fL1V0[2];          // V0 charges
	Int_t             fL1FrameMask;             // Frame mask
	Int_t           fL1TriggerType[19];         // Trigger type
	
	Int_t           fL1DataDecoded;             // Raw data decoded
	Int_t               fL1RawData;             // Raw data
        Int_t                  fMedian;             // Median
        
	ClassDef(AliEMCALTriggerData,3)
};

#endif
