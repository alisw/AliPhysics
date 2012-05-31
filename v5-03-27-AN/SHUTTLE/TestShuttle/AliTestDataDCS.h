#ifndef ALITESTDATADCS_H
#define ALITESTDATADCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

//
// Class to simulate DCS Data point
// for a Test Preprocessor
//

#include <TClonesArray.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TF1.h>

class TMap;

class AliTestDataDCS : public TObject {
public:
	enum {kNAliases=6, kNHistos=3, kNGraphs=3, kNFunctions=2};
	enum {kHvMin=0, kHvMax=20};

	AliTestDataDCS();
	AliTestDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliTestDataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() const {return fRun;}
	Int_t GetStartTime() const {return fStartTime;}
	Int_t GetEndTime() const {return fEndTime;}

	void ProcessData(TMap& aliasMap);

	const char* GetAliasName(UInt_t pos)
			{return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
	const TH1F* GetHisto(UInt_t pos) const 
			{return pos<kNHistos ? fHv[pos] : 0;}

	Float_t GetMean(UInt_t pos) const {return pos<kNHistos ? fMean[pos] : 0;}
	Float_t GetWidth(UInt_t pos) const {return pos<kNHistos ? fWidth[pos] : 0;}

	const TGraph* GetGraph(UInt_t pos)
			{return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos)) : 0;}

	const TF1* GetFunction() const {return fFunc;}

	Double_t Eval(int pos, Double_t time)
			{return pos<kNGraphs ? ((TGraph*) fGraphs.UncheckedAt(pos))->Eval(time) : -1;}

	void Draw(const Option_t* option);


private:
	AliTestDataDCS(const AliTestDataDCS& other);
	AliTestDataDCS& operator= (const AliTestDataDCS& other);
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr);
	void CreateGraph(int i, int dim, const Double_t *x, const Double_t *y);

	Int_t fRun;   // run number
	UInt_t fStartTime;  // run start time
	UInt_t fEndTime;    // run end time

	Float_t fMean[kNHistos];   // mean of the histograms
	Float_t fWidth[kNHistos];  // width of the histograms

	TString fAliasNames[kNAliases];  // aliases of the DPs
	TH1F *fHv[kNHistos];   // histograms
	TClonesArray fGraphs;  // graphs
	TF1 *fFunc;   // functions

	Bool_t fIsProcessed;  // flag indicating whether the DCS data have been processed

	ClassDef(AliTestDataDCS, 2);
};

#endif
