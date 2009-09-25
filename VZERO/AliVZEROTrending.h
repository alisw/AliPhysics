#ifndef ALIVZEROTRENDING_H
#define ALIVZEROTRENDING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */


// 
// Class AliVZEROTrending
// ---------------------------
// 
//  class used in QA to publish variables evolution versus time in AMORE. 
//  These histo are the one which will be looked at by QA Shifter
// 


#include <TH1.h>

class TGraph;

class AliVZEROTrending  : public TH1 {
public:
	AliVZEROTrending();
	AliVZEROTrending(const char* name, const char* title);
	virtual ~AliVZEROTrending();
	AliVZEROTrending(const AliVZEROTrending &trend);
		
	Double_t * GetTime(){return fTime;};
	Double_t * GetChannel(Int_t i){return fData[i];};
	Double_t  GetLastTime(){return fTime[fNEntries-1];};
	Double_t  GetLastChannel(Int_t i){return fData[i][fNEntries];};
	UInt_t GetNEntries(){return fNEntries;};
	void AddEntry(Double_t * data, UInt_t time);
	void PrintEntry(UInt_t entry);	

private:
	
	AliVZEROTrending& operator= (const AliVZEROTrending & /*trend*/); // Not implemented
	enum{kDataSize = 500};
	Double_t fData[8][kDataSize];
	Double_t fTime[kDataSize];
	UInt_t fNEntries;
	
	
	ClassDef( AliVZEROTrending, 1 )  
	
};

#endif // ALIVZEROTRENDING_H


