#ifndef AliADDataDCS_H
#define AliADDataDCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h> 
#include <TClonesArray.h>

class TMap;
class TH2F;
class TGraph;
class TF1;
class TString;
class TH1F;

// AliADDataDCS class
// main aim is to process DCS data
// in order to obtain the data to be stored in the OCDB

class AliADDataDCS : public TObject {
public:
  enum {kNAliases=213,kNGraphs=32,kNHvChannel=16,kNLvChannel=4,kNCIUBoards = 2};
  enum {kHvMin=0, kHvMax=3000};
  
  AliADDataDCS();
  AliADDataDCS(Int_t nRun, UInt_t timeCreated, UInt_t timeCompleted, UInt_t daqStart, UInt_t daqEnd, UInt_t ctpStart, UInt_t ctpEnd);
  ~AliADDataDCS();
  
  void SetRun(Int_t run) {fRun = run;}
  void SetStartTime(Int_t startTime) {fStartTime = startTime;}
  void SetEndTime(Int_t endTime) {fEndTime = endTime;}
  void SetDaqStartTime(Int_t startTime) {fDaqStartTime = startTime;}
  void SetDaqEndTime(Int_t endTime) {fDaqEndTime = endTime;}
  Int_t GetRun() const {return fRun;}
  Int_t GetStartTime() const {return fStartTime;}
  Int_t GetEndTime() const {return fEndTime;}
  Int_t GetDaqStartTime() const {return fDaqStartTime;}
  Int_t GetDaqEndTime() const {return fDaqEndTime;}
  
  Bool_t ProcessData(TMap& aliasMap);
  
  const char* GetAliasName(Int_t pos) const 
    {return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
  
  void Draw(const Option_t* option);
  
  Float_t* GetMeanHV()    const {return (float*)fMeanHV;}
  Float_t* GetWidthHV()   const {return (float*)fWidthHV;}
  Bool_t * GetDeadMap()   const {return (bool*)fDeadChannel;}
  TMap * GetFEEParameters() const {return fFEEParameters;};
  TClonesArray * GetGraphs() const {return fGraphs;};
  
private:
  AliADDataDCS(const AliADDataDCS&); // Not implemented
  AliADDataDCS& operator=(const AliADDataDCS&); // Not implemented

  void Init();
  void Introduce(UInt_t numAlias, const TObjArray* aliasArr) const;
  void CreateGraph(int i, int dim, const Double_t *x, const Double_t *y);
    
  Int_t fRun;       // Run number
  Int_t fStartTime; // start time (time created)
  Int_t fEndTime;   // end time (time completed)
  UInt_t fDaqStartTime; // DAQ start time
  UInt_t fDaqEndTime;   // DAQ end time
  UInt_t fCtpStartTime; // CTP start time
  UInt_t fCtpEndTime;   // CTP end time
  
  TString fAliasNames[kNAliases];        // aliases for DCS data
  TClonesArray *fGraphs;		         // Array containing  graphics
  TH1F *fHv[kNHvChannel];                  // High Voltage histograms
  Float_t fMeanHV[kNHvChannel];            // High Voltage mean values
  Float_t fWidthHV[kNHvChannel];           // High Voltage widths
  Bool_t  fDeadChannel[kNHvChannel];       // Dead Map 
  TMap * fFEEParameters;  // TMap holding the FEE parameters
    
  Bool_t fIsProcessed;                   // bool to know processing status
  
  ClassDef(AliADDataDCS, 7);
};

#endif
