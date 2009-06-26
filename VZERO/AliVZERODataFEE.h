#ifndef ALIVZERODATAFEE_H
#define ALIVZERODATAFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
class TMap;

#include <TObject.h>
#include <TString.h>
#include <TMap.h>

class AliVZERODataFEE : public TObject {
public:
	//enum {kNAliases=64};
	AliVZERODataFEE();
	AliVZERODataFEE(Int_t nRun, UInt_t startTime, UInt_t endTime);
	AliVZERODataFEE(const AliVZERODataFEE &dataFEE);
	AliVZERODataFEE& operator= (const AliVZERODataFEE &dataFEE);
	virtual ~AliVZERODataFEE();
	
	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() const {return fRun;}
	Int_t GetStartTime() const {return fStartTime;}
	Int_t GetEndTime() const {return fEndTime;}

	void ProcessData(TMap& aliasMap);
	void Init();
	void PrintAliases();
	
	TMap * GetParameters() const {return fParameters;};
	
	enum { kNCIUBoards = 8, kNCIUParam = 13, kNChannelParam = 8, kNCCIUParam = 19, kNAliases  = kNChannelParam*8*kNCIUBoards +kNCIUParam*kNCIUBoards + kNCCIUParam };

private:

	Int_t fRun;       // Run number
	Int_t fStartTime; // Start time
	Int_t fEndTime;   // End time
	TString fAliasNames[kNAliases];	// aliases for DCS data
	Bool_t fIsProcessed; // bool to know processing status

	TString GetFEEParamName(Int_t iParam);
	
	TMap * fParameters;
	
	ClassDef( AliVZERODataFEE, 1 )  
	
};


#endif // ALIVZERODATAFEE_H
