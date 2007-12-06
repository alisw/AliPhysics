#ifndef AliT0DataDCS_H
#define AliT0DataDCS_H

#include <TMap.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>

class AliT0DataDCS : public TObject {
public:
	enum {kNAliases=6};


	AliT0DataDCS();
	AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliT0DataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	Int_t GetRun() {return fRun;}
	Int_t GetStartTime() {return fStartTime;}
	Int_t GetEndTime() {return fEndTime;}

	Bool_t ProcessData(TMap& aliasMap);

	//const char* GetAliasName(UInt_t pos)
	//		{return pos<kNAliases ? fAliasNames[pos].Data() : 0;}

	void CalcScalerMean(Float_t* t0_scal);


private:
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr);

	Int_t fRun;
	UInt_t fStartTime;
	UInt_t fEndTime;


	TString fAliasNames[kNAliases];

	Bool_t fIsProcessed;
	void SetScalerMean(Int_t channel, Float_t val) {fScalerMean[channel]=val;}
	Float_t  GetScalerMean(Int_t channel)        const {return fScalerMean[channel];}
	Float_t* GetScalerMean()          const {return (float*)fScalerMean;}


	ClassDef(AliT0DataDCS, 2);

protected:
	Float_t     fScalerMean[24];	
};

#endif
