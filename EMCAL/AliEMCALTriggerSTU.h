#ifndef ALIEMCALTRIGGERSTU_H
#define ALIEMCALTRIGGERSTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <AliEMCALTriggerBoard.h>

class TTree;

class AliEMCALTriggerSTU : public AliEMCALTriggerBoard 
{
public:
	
	              AliEMCALTriggerSTU();
	              AliEMCALTriggerSTU(AliEMCALCalibData *calibData, const TVector2& rsize);
	virtual      ~AliEMCALTriggerSTU();
	
	virtual void  FetchFOR(Int_t i, Int_t** Region, const TVector2* rSize);
	virtual void  BuildMap(Int_t i, Int_t**    Map, const TVector2* rSize);
	virtual void  PrintADC(L1TriggerType_t type, TVector2& pos, TVector2& idx);
	virtual void  L1(L1TriggerType_t type);//, TTree& treeV0);
	virtual void  PatchGenerator(const TClonesArray* lpos, Int_t val);
	virtual const Int_t* V0() const {return fV0M;}
	virtual void  SetV0Multiplicity(Int_t M[], Int_t n) { for (Int_t i=0;i<n;i++) fV0M[i] = M[i]; } // for raw data rec
	virtual void  V0Multiplicity(TTree& treeV0);
	virtual void  Reset();
	
protected:
		            AliEMCALTriggerSTU(const AliEMCALTriggerSTU& rhs);
	                AliEMCALTriggerSTU& operator=(const AliEMCALTriggerSTU& rhs);

private:
	
	virtual Float_t ThresholdFromV0(L1TriggerType_t type);
	
	        Int_t   fV0M[2]; //! 0/1: V0C/V0A

	ClassDef(AliEMCALTriggerSTU,1)
};
 
#endif
