#ifndef ALIVZERODATAFEE_H
#define ALIVZERODATAFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
class TMap;
class TString;

#include <TObject.h>

//
// Class AliVZERODataFEE
// ---------------------
// Used to process the TMap of DCS values comming from the shuttle.
// It stores into a TMap the FEE parameters for the given run number
//


class AliVZERODataFEE : public TObject {
public:
	//enum {kNAliases=64};
	AliVZERODataFEE();
	AliVZERODataFEE(Int_t nRun, UInt_t startTime, UInt_t endTime);
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
	
	AliVZERODataFEE(const AliVZERODataFEE & /*dataFEE*/); // Not implemented
	AliVZERODataFEE& operator= (const AliVZERODataFEE &/*dataFEE*/); // Not implemented

	Int_t fRun;       // Run number
	Int_t fStartTime; // Start time
	Int_t fEndTime;   // End time
	TString fAliasNames[kNAliases];	// aliases for DCS data
	Bool_t fIsProcessed; // bool to know processing status
	TMap * fParameters;  // TMap holding the FEE parameters

	TString GetFEEParamName(Int_t iParam);
	
	
	ClassDef( AliVZERODataFEE, 1 )  
	
};


#endif // ALIVZERODATAFEE_H

