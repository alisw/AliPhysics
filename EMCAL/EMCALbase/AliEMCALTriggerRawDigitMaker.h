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
  enum {kMaxDigitIndex=5952};
	         AliEMCALTriggerRawDigitMaker();
	virtual ~AliEMCALTriggerRawDigitMaker();
	
	virtual void SetIO(AliRawReader* reader, AliCaloRawStreamV3& in, AliEMCALTriggerSTURawStream& inSTU, TClonesArray* digits, TClonesArray* data);
	virtual void Add(const std::vector<AliCaloBunchInfo> &bunchlist);
	virtual void PostProcess();
	virtual void Reset();
	
protected:
	
	AliEMCALGeometry*            fGeometry;      // Geometry
	AliRawReader*                fRawReader;     // Raw reader
	AliCaloRawStreamV3*          fCaloRawStream; // Calo raw stream
	AliEMCALTriggerSTURawStream* fSTURawStream;  // STU raw stream
	TClonesArray*                fRawDigits;     // Raw digits
	AliCaloRawAnalyzerFakeALTRO* fRawAnalyzer;   // Raw analyzer
	AliEMCALTriggerDCSConfigDB*  fDCSConfig;     // DCS config
	TClonesArray*                fTriggerData;   // Trigger data
	
	Int_t                        fRawDigitIndex[kMaxDigitIndex]; // Raw digit indexes

private:
	
        AliEMCALTriggerRawDigitMaker(const AliEMCALTriggerRawDigitMaker& rhs);            // NOT implemented
        AliEMCALTriggerRawDigitMaker& operator=(const AliEMCALTriggerRawDigitMaker& rhs); // NOT implemented	
	
	ClassDef(AliEMCALTriggerRawDigitMaker,2)
};
 
#endif
