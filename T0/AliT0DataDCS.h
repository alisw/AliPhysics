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

	Int_t fRun;			// Run number
	UInt_t fStartTime;		// Start time
	UInt_t fEndTime;		// End time
	Int_t fStartTimeDCSQuery; 	// Begin time DCSQuery
  	Int_t fEndTimeDCSQuery;   	// End time DCSQuery
	Float_t fHViA[kHV];		// Mean value of HV current in uA on A-side
	Float_t fHVvA[kHV]; 		// Mean value of HV voltage in V on A-side
	Float_t fLViA[kLV]; 		// Mean value of LV current in uA on A-side
	Float_t fLVvA[kLV];		// Mean value of LV voltage in V on A-side
	Float_t fHViC[kHV];		// Mean value of HV current in uA on C-side
        Float_t fHVvC[kHV];		// Mean value of HV voltage in V on C-side
        Float_t fLViC[kLV];		// Mean value of LV current in uA on C-side
        Float_t fLVvC[kLV];		// Mean value of LV voltage in V on C-side
 	Float_t fCFDtA[kCFD];		// Mean threshold on CFD in V on A-side	
	Float_t fCFDwA[kCFD];		// Mean walk on CFD in V on A-side
	Float_t fCFDtC[kCFD];		// Mean threshold on CFD in V on C-side
        Float_t fCFDwC[kCFD];		// Mean walk on CFD in V on C-side
        Float_t fScalerMean[kScalers];  // Mean value of T0 scaler counts from the entire run
        Float_t fScalerSecMean[kScalers]; // Mean value of T0 scaler counts per second
	Float_t fTRM[kTRM];		// Mean temperature on TRM in degrees of Celsius
	Float_t fDRM[kDRM];		// Mean temperature on DRM in degrees of Celsius
	Float_t fAtten;			// Laser amplitude in MIPs
	TString fAliasNames[kNAliases];		// T0 data points aliases  
	Bool_t fIsProcessed;			// status - was processing data successful
	ClassDef(AliT0DataDCS, 3);

protected:
};

#endif
