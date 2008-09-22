#ifndef AliTOFDataDCS_H
#define AliTOFDataDCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h" 
//#include "TString.h"

class TMap;
class TClonesArray;
//class TH2F;
//class TGraph;
//class TF1;
class TString;
class AliTOFFormatDCS;

// AliTOFDataDCS class
// main aim is to process DCS data
// in order to obtain the data to be stored in the OCDB

class AliTOFDataDCS : public TObject {
public:
  enum {kNAliases=360, kNHV=90};
  
  AliTOFDataDCS();
  AliTOFDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery );
  AliTOFDataDCS(const AliTOFDataDCS & data);
  AliTOFDataDCS& operator=(const AliTOFDataDCS & data);
  ~AliTOFDataDCS();
  
  void SetRun(Int_t run) {fRun = run;}
  void SetStartTime(Int_t startTime) {fStartTime = startTime;}
  void SetEndTime(Int_t endTime) {fEndTime = endTime;}
  void SetStartTimeDCSQuery(Int_t startTimeDCSQuery) {fStartTimeDCSQuery = startTimeDCSQuery;}
  void SetEndTimeDCSQuery(Int_t endTimeDCSQuery) {fEndTimeDCSQuery = endTimeDCSQuery;}
  Int_t GetRun() const {return fRun;}
  Int_t GetStartTime() const {return fStartTime;}
  Int_t GetEndTime() const {return fEndTime;}
  Int_t GetStartTimeDCSQuery() const {return fStartTimeDCSQuery;}
  Int_t GetEndTimeDCSQuery() const {return fEndTimeDCSQuery;}
  
  Bool_t ProcessData(TMap& aliasMap);
  
  const char* GetAliasName(Int_t pos) const 
    {return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
  
  void Draw(const Option_t* /*option*/);
  
  AliTOFFormatDCS* GetHVvp(Int_t pos) const
    {return pos<kNHV ? fHVvp[pos] : 0;}
  AliTOFFormatDCS* GetHVvn(Int_t pos) const 
    {return pos<kNHV ? fHVvn[pos] : 0;}
  AliTOFFormatDCS* GetHVip(Int_t pos) const 
    {return pos<kNHV ? fHVip[pos] : 0;}
  AliTOFFormatDCS* GetHVin(Int_t pos) const 
    {return pos<kNHV ? fHVin[pos] : 0;}

  void SetFDRFlag(Bool_t flag) {fFDR = flag;}
  Bool_t GetFDRFlag() const {return fFDR;}

private:
  void Init();
  void Introduce(UInt_t numAlias, const TObjArray* aliasArr) const;
  void CreateHisto(int nbin);
  
  Int_t fRun;       // Run number
  Int_t fStartTime; // start time
  Int_t fEndTime;   // end time  
  Int_t fStartTimeDCSQuery; // start time DCSQuery
  Int_t fEndTimeDCSQuery;   // end time DCSQuery
  
  TString fAliasNames[kNAliases];        // aliases for DCS data
  AliTOFFormatDCS *fHVvp[kNHV];          // HV voltages, positive ch
  AliTOFFormatDCS *fHVvn[kNHV];          // HV voltages, negative ch
  AliTOFFormatDCS *fHVip[kNHV];          // HV currents, positive ch
  AliTOFFormatDCS *fHVin[kNHV];          // HV currents, negative ch
  
  Bool_t fIsProcessed;                   // bool to know processing status
  Bool_t fFDR;                   // bool to know whether we are in a FDR run
  
  ClassDef(AliTOFDataDCS, 5);
};

#endif
