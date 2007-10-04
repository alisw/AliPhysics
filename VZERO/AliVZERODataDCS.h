#ifndef AliVZERODataDCS_H
#define AliVZERODataDCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h" 
#include "TString.h"

#include <TClonesArray.h>
#include <TH1F.h>

class TMap;
class TH2F;
class TGraph;
class TF1;

// AliVZERODataDCS class
// main aim is to process DCS data
// in order to obtain the data to be stored in the OCDB

class AliVZERODataDCS : public TObject {
public:
  enum {kNAliases=64,kNGraphs=64};
  enum {kHvMin=0, kHvMax=2000};
  
  AliVZERODataDCS();
  AliVZERODataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
  ~AliVZERODataDCS();
  
  void SetRun(Int_t run) {fRun = run;}
  void SetStartTime(Int_t startTime) {fStartTime = startTime;}
  void SetEndTime(Int_t endTime) {fEndTime = endTime;}
  Int_t GetRun() const {return fRun;}
  Int_t GetStartTime() const {return fStartTime;}
  Int_t GetEndTime() const {return fEndTime;}
  
  void ProcessData(TMap& aliasMap);
  
  const char* GetAliasName(Int_t pos) const 
    {return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
  
  void Draw(const Option_t* option);
  
  Float_t* GetMeanHV()   const {return (float*)fMeanHV;}
  Float_t* GetWidthHV()   const {return (float*)fWidthHV;}

private:
  void Init();
  void Introduce(UInt_t numAlias, const TObjArray* aliasArr) const;
  void CreateGraph(int i, int dim, const Double_t *x, const Double_t *y);
    
  Int_t fRun;       // Run number
  Int_t fStartTime; // start time
  Int_t fEndTime;   // end time
  
  
  TString fAliasNames[kNAliases];        // aliases for DCS data
  TClonesArray fGraphs;		// Array containing  graphics
  TH1F *fHv[kNAliases];
  Float_t fMeanHV[kNAliases]; 
  Float_t fWidthHV[kNAliases];
 
  Bool_t fIsProcessed;                   // bool to know processing status
  
  ClassDef(AliVZERODataDCS, 2);
};

#endif
