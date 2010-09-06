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
	
	virtual void  Build(TString& str, Int_t i, Int_t** Map, const TVector2* rSize);
	virtual void  PrintADC(TriggerType_t type, TVector2& pos, TVector2& idx);
	virtual void  L1(TriggerType_t type);
	virtual void  PatchGenerator(const TClonesArray* lpos, Int_t val);
	
	virtual void  ComputeThFromV0(const Int_t M[]);
	
	virtual void  SetThreshold(TriggerType_t type, Int_t v);
	
	virtual Int_t GetThreshold(TriggerType_t type);
	virtual Int_t GetRawData() const;

	virtual void  Reset();
	
protected:
	
				  AliEMCALTriggerSTU(const AliEMCALTriggerSTU& rhs);
				  AliEMCALTriggerSTU& operator=(const AliEMCALTriggerSTU& rhs);

private:
	
		  Int_t   fGammaTh;
		  Int_t   fJetTh;
	
	AliEMCALTriggerSTUDCSConfig* fDCSConfig;

	ClassDef(AliEMCALTriggerSTU,1)
};
 
#endif
