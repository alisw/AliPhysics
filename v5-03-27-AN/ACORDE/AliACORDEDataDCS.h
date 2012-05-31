/****************************************************
 AliACORDEDataDCS class
 	Pointer to the DCS objects
 Author: Pedro Gonzalez (CIEMAT, Madrid)
 Last updated: March 21 th 2010
	Fixing coding violatgion
	Mario Rodriguez (FCFM-BUAP, Puebla MX)

*****************************************************/
#ifndef AliACORDEDATADCS_H
#define AliACORDEDATADCS_H

#include <TClonesArray.h>
#include <TGraph.h>
#include <TMap.h>
class AliACORDEDataDCS : public TObject {
public:
	enum {kNAliases=60, kNHistos=60, kNGraphs=60, kNFunctions=2};
	enum {kHvMin=0, kHvMax=20};

	AliACORDEDataDCS();
	AliACORDEDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
        AliACORDEDataDCS(const AliACORDEDataDCS & data);
        AliACORDEDataDCS& operator=(const AliACORDEDataDCS & data);

	~AliACORDEDataDCS();

	void SetRun(Int_t run) {fRun = run;} 
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() const {return fRun;}
	Int_t GetStartTime() const {return fStartTime;}
	Int_t GetEndTime() const {return fEndTime;}

	void ProcessData(TMap& aliasMap);

	const char* GetAliasName(UInt_t pos)
			{return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
	TH1F* GetHisto(UInt_t pos) const {return pos<kNHistos ? fHv[pos] : 0;}

	Float_t GetMean(UInt_t pos) const {return pos<kNHistos ? fMean[pos] : 0;}
	Float_t GetWidth(UInt_t pos) const {return pos<kNHistos ? fWidth[pos] : 0;}

	const TGraph* GetGraph(UInt_t pos) {return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos)) : 0;}

	const TF1* GetFunction() const {return fFunc;}

	Double_t Eval(int pos, Double_t time)
			{return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos))->Eval(time) : -1;}

	void Draw(const Option_t* option);


private:
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr);
	void CreateGraph(int i, int dim, const Double_t *x, const Double_t *y);

	Int_t fRun; // # of Run
	UInt_t fStartTime; // Start time of the Run
	UInt_t fEndTime; // End time of the Run

	Float_t fMean[kNHistos]; // Mean of hits distribution for ACORDE
	Float_t fWidth[kNHistos]; // Width of the hits dist. for ACORDE

	TString fAliasNames[kNAliases]; // Alias names for ACORDE Data Points
	TH1F *fHv[kNHistos]; // High Voltage values
	TClonesArray fGraphs; // Clones of plots for ACORDE
	TF1 *fFunc; // Funtion for ACORDE DP

	Bool_t fIsProcessed; // Boolean flag to know if the vent was processed

	ClassDef(AliACORDEDataDCS, 1);
};

#endif
