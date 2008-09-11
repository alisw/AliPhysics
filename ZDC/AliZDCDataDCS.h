#ifndef AliZDCDataDCS_H
#define AliZDCDataDCS_H

////////////////////////////////////////////////
//  Class for ZDC DCS data                    //
////////////////////////////////////////////////


#include <TMap.h>
#include <TClonesArray.h>
#include <TGraph.h>

class AliZDCDataDCS : public TObject {
public:
	enum {kNAliases=28, kNGraphs=24};

	AliZDCDataDCS();
	AliZDCDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliZDCDataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	void SetCalibData(Float_t *val) 
	     {for(Int_t i=0; i<kNGraphs; i++) fCalibData[i] = val[i];}
	void SetCalibData(Int_t i, Float_t val) {fCalibData[i] = val;} 
	//
	Int_t GetRun() {return fRun;}
	Int_t GetStartTime() {return fStartTime;}
	Int_t GetEndTime() {return fEndTime;}
	Float_t GetCalibData() {return *fCalibData;}
	Float_t GetCalibData(Int_t i) {return fCalibData[i];}

	void ProcessData(TMap& aliasMap);

	const char* GetAliasName(UInt_t pos)
			{return pos<kNAliases ? fAliasNames[pos].Data() : 0;}

	const TGraph* GetGraph(UInt_t pos)
			{return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos)) : 0;}

	Double_t Eval(int pos, Double_t time)
			{return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos))->Eval(time) : -1;}

	void Draw(const Option_t* option);


private:
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr);
	void CreateGraph(int i, int dim, const Double_t *x, const Double_t *y);

	Int_t fRun;		// Run number
	UInt_t fStartTime;	// Start of run time
	UInt_t fEndTime;	// End of run time

	TString fAliasNames[kNAliases];	// Name of the aliases provided by the DCS
	Double_t fCalibData[kNAliases];	// Array containing calibration data
	Double_t fTimeStamp[kNAliases];	// Array containing time stamps

	TClonesArray fGraphs;		// Array containing PTM HV graphics

	Bool_t fIsProcessed;		// Flag set when data are processed

	ClassDef(AliZDCDataDCS, 3);
};

#endif
