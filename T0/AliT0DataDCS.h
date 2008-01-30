#ifndef AliT0DataDCS_H
#define AliT0DataDCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TMap.h>
#include <TClonesArray.h>
#include <TGraph.h>

// AliT0DataDCS class
// fetching T0 data points from DCS, calculating mean values for the run
// and storing the result to Reference DB

class AliT0DataDCS : public TObject {
public:
	enum {kNAliases=32};

	AliT0DataDCS();
	AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliT0DataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() const {return fRun;}
	Int_t GetStartTime() const {return fStartTime;}
	Int_t GetEndTime() const {return fEndTime;}

	Bool_t ProcessData(TMap& aliasMap);

private:
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr);

	Int_t fRun;	// Run number
	UInt_t fStartTime;	// Start time
	UInt_t fEndTime;	// End time
	Float_t fScalerMean[32];	// Mean value of T0 scaler counts from the entire run 
	TString fAliasNames[kNAliases];		// T0 data points aliases  
	Bool_t fIsProcessed;			// status - was processing data successful
	ClassDef(AliT0DataDCS, 2);

protected:
};

#endif
