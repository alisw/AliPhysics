#ifndef AliZDCDataDCS_H
#define AliZDCDataDCS_H

#include <TMap.h>
#include <TClonesArray.h>
#include <TGraph.h>

class AliZDCDataDCS : public TObject {
public:
	enum {kNAliases=26, kNGraphs=22};

	AliZDCDataDCS();
	AliZDCDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliZDCDataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() {return fRun;}
	Int_t GetStartTime() {return fStartTime;}
	Int_t GetEndTime() {return fEndTime;}

	void ProcessData(TMap& aliasMap, Float_t *CalibData);

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

	Int_t fRun;
	UInt_t fStartTime;
	UInt_t fEndTime;

	TString fAliasNames[kNAliases];
	TClonesArray fGraphs;
	
	Float_t fCalibData[kNAliases];

	Bool_t fIsProcessed;

	ClassDef(AliZDCDataDCS, 2);
};

#endif
