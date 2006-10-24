#ifndef AliSTARTDataDCS_H
#define AliSTARTDataDCS_H

#include <TMap.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>

class AliSTARTDataDCS : public TObject {
public:
	enum {kNAliases=6, kNHistos=3, kNGraphs=3, kNFunctions=2};
	enum {kHvMin=0, kHvMax=20};

	AliSTARTDataDCS();
	AliSTARTDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliSTARTDataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() {return fRun;}
	Int_t GetStartTime() {return fStartTime;}
	Int_t GetEndTime() {return fEndTime;}

	void ProcessData(TMap& aliasMap);

	const char* GetAliasName(UInt_t pos)
			{return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
	const TH1F* GetHisto(UInt_t pos) 
			{return pos<kNHistos ? fHv[pos] : 0;}

	Float_t GetMean(UInt_t pos) {return pos<kNHistos ? fMean[pos] : 0;}
	Float_t GetWidth(UInt_t pos) {return pos<kNHistos ? fWidth[pos] : 0;}

	const TGraph* GetGraph(UInt_t pos)
			{return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos)) : 0;}

	const TF1* GetFunction() {return fFunc;}

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

	Float_t fMean[kNHistos];
	Float_t fWidth[kNHistos];

	TString fAliasNames[kNAliases];
	TH1F *fHv[kNHistos];
	TClonesArray fGraphs;
	TF1 *fFunc;

	Bool_t fIsProcessed;

	ClassDef(AliSTARTDataDCS, 2);
};

#endif
