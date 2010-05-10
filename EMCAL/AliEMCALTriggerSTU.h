#ifndef ALIEMCALTRIGGERSTU_H
#define ALIEMCALTRIGGERSTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <AliEMCALTriggerBoard.h>

class TTree;
class AliEMCALTriggerSTUDCSConfig;

class AliEMCALTriggerSTU : public AliEMCALTriggerBoard 
{
public:
	
	              AliEMCALTriggerSTU();
	              AliEMCALTriggerSTU(AliEMCALTriggerSTUDCSConfig *dcsConf, const TVector2& rsize);
	virtual      ~AliEMCALTriggerSTU();
	
	virtual void  FetchFOR(Int_t i, Int_t** Region, const TVector2* rSize);
	virtual void  BuildMap(Int_t i, Int_t**    Map, const TVector2* rSize);
	virtual void  PrintADC(L1TriggerType_t type, TVector2& pos, TVector2& idx);
	virtual void  L1(L1TriggerType_t type);//, TTree& treeV0);
	virtual void  PatchGenerator(const TClonesArray* lpos, Int_t val);
	virtual const Int_t* V0() const {return fV0M;}
	virtual void  SetV0Multiplicity(const Int_t M[], Int_t n);
	virtual void  Reset();
	virtual Int_t GetThreshold(L1TriggerType_t type);
	
protected:
		            AliEMCALTriggerSTU(const AliEMCALTriggerSTU& rhs);
	                AliEMCALTriggerSTU& operator=(const AliEMCALTriggerSTU& rhs);

private:
	
	        Int_t   fV0M[2]; //! 0/1: V0C/V0A
		  Int_t   fGammaTh;
		  Int_t   fJetTh;
	
	AliEMCALTriggerSTUDCSConfig* fDCSConfig;

	ClassDef(AliEMCALTriggerSTU,1)
};
 
#endif
