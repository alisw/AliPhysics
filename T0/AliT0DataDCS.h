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
	enum {kNAliases=184, kHV=12, kLV=2, kCFD=12, kScalers=32, kTRM=10, kDRM=5}; //184

	AliT0DataDCS();
	AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery);
	AliT0DataDCS(const AliT0DataDCS & data);
  	AliT0DataDCS& operator=(const AliT0DataDCS & data);
	~AliT0DataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	void SetStartTimeDCSQuery(Int_t startTimeDCSQuery) {fStartTimeDCSQuery = startTimeDCSQuery;}
	void SetEndTimeDCSQuery(Int_t endTimeDCSQuery) {fEndTimeDCSQuery = endTimeDCSQuery;}
	Int_t GetRun() const {return fRun;}
	Int_t GetStartTime() const {return fStartTime;}
	Int_t GetEndTime() const {return fEndTime;}
	Int_t GetStartTimeDCSQuery() const {return fStartTimeDCSQuery;}
  	Int_t GetEndTimeDCSQuery() const {return fEndTimeDCSQuery;}	

	Bool_t ProcessData(TMap& aliasMap);

	const char* GetAliasName(Int_t pos) const {return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
private:
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr)const;

	Int_t fRun;	// Run number
	UInt_t fStartTime;	// Start time
	UInt_t fEndTime;	// End time
	Int_t fStartTimeDCSQuery; // Begin time DCSQuery
  	Int_t fEndTimeDCSQuery;   // End time DCSQuery
	Float_t fHViA[kHV];
	Float_t fHVvA[kHV]; 
	Float_t fLViA[kLV]; 
	Float_t fLVvA[kLV];
	Float_t fHViC[kHV];
        Float_t fHVvC[kHV];
        Float_t fLViC[kLV];
        Float_t fLVvC[kLV];
 	Float_t fCFDtA[kCFD];
	Float_t fCFDwA[kCFD];
	Float_t fCFDtC[kCFD];
        Float_t fCFDwC[kCFD];
        Float_t fScalerMean[kScalers];  // Mean value of T0 scaler counts from the entire run
        Float_t fScalerSecMean[kScalers];
	Float_t fTRM[kTRM];
	Float_t fDRM[kDRM];
	Float_t fAtten;	
	TString fAliasNames[kNAliases];		// T0 data points aliases  
	Bool_t fIsProcessed;			// status - was processing data successful
	ClassDef(AliT0DataDCS, 3);

protected:
};

#endif
