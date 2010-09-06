#ifndef ALIEMCALTRIGGERRAWDIGITMAKER_H
#define ALIEMCALTRIGGERRAWDIGITMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <vector>
#include <TObject.h>

class AliEMCALGeometry;
class AliRawReader;
class AliCaloRawStreamV3;
class AliEMCALTriggerSTURawStream;
class AliCaloRawAnalyzerFakeALTRO;
class AliCaloBunchInfo;
class TClonesArray;
class AliEMCALTriggerDCSConfigDB;
class AliEMCALTriggerData;

class AliEMCALTriggerRawDigitMaker : public TObject 
{
	
public:
	         AliEMCALTriggerRawDigitMaker();
	virtual ~AliEMCALTriggerRawDigitMaker();
	
	virtual void SetIO(AliRawReader* reader, AliCaloRawStreamV3& in, AliEMCALTriggerSTURawStream& inSTU, TClonesArray* digits, AliEMCALTriggerData* data);
	virtual void Add(const std::vector<AliCaloBunchInfo> &bunchlist);
	virtual void PostProcess();
	virtual void Reset();
	
private:
	
    AliEMCALTriggerRawDigitMaker(const AliEMCALTriggerRawDigitMaker& rhs);            // NOT implemented
	AliEMCALTriggerRawDigitMaker& operator=(const AliEMCALTriggerRawDigitMaker& rhs); // NOT implemented
	
protected:
	
	AliEMCALGeometry*            fGeometry;
	AliRawReader*                fRawReader;
	AliCaloRawStreamV3*          fCaloRawStream;
	AliEMCALTriggerSTURawStream* fSTURawStream;
	TClonesArray*                fRawDigits;
	AliCaloRawAnalyzerFakeALTRO* fRawAnalyzer;
	AliEMCALTriggerDCSConfigDB*  fDCSConfig;
	AliEMCALTriggerData*         fTriggerData;
	
	Int_t						 fRawDigitIndex[3072];
	
	ClassDef(AliEMCALTriggerRawDigitMaker,1)
};
 
#endif
