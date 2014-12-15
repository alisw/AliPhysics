#ifndef AliVZERODataDCS_H
#define AliVZERODataDCS_H

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

// AliVZERODataDCS class
// main aim is to process DCS data
// in order to obtain the data to be stored in the OCDB

class AliVZERODataDCS : public TObject {
public:
  enum {kNAliases=240,kNGraphs=64,kNHvChannel=64,kNLvChannel=16,kNCIUBoards = 8};
  enum {kHvMin=0, kHvMax=3000};
  
  AliVZERODataDCS();
  AliVZERODataDCS(Int_t nRun, UInt_t timeCreated, UInt_t timeCompleted, UInt_t daqStart, UInt_t daqEnd, UInt_t ctpStart, UInt_t ctpEnd);
  ~AliVZERODataDCS();
  
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
  
// Getter of Offline Channel number as used in aliroot (defined by aliroot 
// numbering convention) from DCS Channel number

    Int_t      GetOfflineChannel(Int_t channel)  const
      { Int_t  fOfflineChannel[64] = {32, 33, 34, 35, 36, 37, 38, 39, 
                                      40, 41, 42, 43, 44, 45, 46, 47, 
			              48, 49, 50, 51, 52, 53, 54, 55, 
			              56, 57, 58, 59, 60, 61, 62, 63,
			               0,  1,  2,  3,  4,  5,  6,  7, 
			               8,  9, 10, 11, 12, 13, 14, 15,
			              16, 17, 18, 19, 20, 21, 22, 23, 
			              24, 25, 26, 27, 28, 29, 30, 31};
               return fOfflineChannel[channel]; }

private:
  AliVZERODataDCS(const AliVZERODataDCS&); // Not implemented
  AliVZERODataDCS& operator=(const AliVZERODataDCS&); // Not implemented

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
  TClonesArray fGraphs;		         // Array containing  graphics
  TH1F *fHv[kNHvChannel];                  // High Voltage histograms
  Float_t fMeanHV[kNHvChannel];            // High Voltage mean values
  Float_t fWidthHV[kNHvChannel];           // High Voltage widths
  Bool_t  fDeadChannel[kNHvChannel];       // Dead Map 
  TMap * fFEEParameters;  // TMap holding the FEE parameters of Time Resolution
    
  Bool_t fIsProcessed;                   // bool to know processing status
  
  ClassDef(AliVZERODataDCS, 7);
};

#endif
