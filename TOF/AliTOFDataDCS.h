#ifndef AliTOFDataDCS_H
#define AliTOFDataDCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h" 
#include "TString.h"

class TMap;
class TClonesArray;
class TH2F;
class TGraph;
class TF1;
class AliTOFFormatDCS;

// AliTOFDataDCS class
// main aim is to process DCS data
// in order to obtain the data to be stored in the OCDB

class AliTOFDataDCS : public TObject {
public:
  enum {kNAliases=10514, kNHV=90, kNLV=576, 
	kNLV33=72, kNLV50=72, kNLV48=72, 
	kNFEEthr=1152, kNFEEtfeac=576, kNFEEttrm=6840, 
	kNFunctions=0, kNT=1, kNP=1};
  
  AliTOFDataDCS();
  AliTOFDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime);
  AliTOFDataDCS(const AliTOFDataDCS & data);
  AliTOFDataDCS& operator=(const AliTOFDataDCS & data);
  ~AliTOFDataDCS();
  
  void SetRun(Int_t run) {fRun = run;}
  void SetStartTime(Int_t startTime) {fStartTime = startTime;}
  void SetEndTime(Int_t endTime) {fEndTime = endTime;}
  Int_t GetRun() const {return fRun;}
  Int_t GetStartTime() const {return fStartTime;}
  Int_t GetEndTime() const {return fEndTime;}
  
  Bool_t ProcessData(TMap& aliasMap);
  
  const char* GetAliasName(Int_t pos) const 
    {return pos<kNAliases ? fAliasNames[pos].Data() : 0;}
  
  
  Float_t* GetT() const;
  Float_t* GetP() const;
  
  Float_t GetT(Int_t i) const {return fT[i];}
  Float_t GetP(Int_t i) const {return fP[i];}
  Float_t GetSlopeT() const {return fT[0];}
  Float_t GetInterceptT() const {return fT[1];}
  Float_t GetMaxT()const {return fT[2];}
  Float_t GetSlopeP() const {return fP[0];}
  Float_t GetInterceptP() const {return fP[1];}
  Float_t GetMaxP() const {return fP[2];}
 
  void SetSlopeT(Float_t slopeT) {fT[0]=slopeT;}
  void SetInterceptT(Float_t interceptT) {fT[1]=interceptT;}
  void SetMaxT(Float_t maxT) {fT[2]=maxT;}
  void SetSlopeP(Float_t slopeP) {fP[0]=slopeP;}
  void SetInterceptP(Float_t interceptP) {fP[1]=interceptP;}
  void SetMaxP(Float_t maxP) {fP[2]=maxP;}
  
  void Draw(const Option_t* option);
  
  AliTOFFormatDCS* GetHVvp(Int_t pos) const
    {return pos<kNHV ? fHVvp[pos] : 0;}
  AliTOFFormatDCS* GetHVvn(Int_t pos) const 
    {return pos<kNHV ? fHVvn[pos] : 0;}
  AliTOFFormatDCS* GetHVip(Int_t pos) const 
    {return pos<kNHV ? fHVip[pos] : 0;}
  AliTOFFormatDCS* GetHVin(Int_t pos) const 
    {return pos<kNHV ? fHVin[pos] : 0;}
  AliTOFFormatDCS* GetLVv(Int_t pos) const 
    {return pos<kNLV ? fLVv[pos] : 0;}
  AliTOFFormatDCS* GetLVi(Int_t pos) const 
    {return pos<kNLV ? fLVi[pos] : 0;}
  AliTOFFormatDCS* GetLVv33(Int_t pos) const 
    {return pos<kNLV ? fLVv33[pos] : 0;}
  AliTOFFormatDCS* GetLVi33(Int_t pos) const 
    {return pos<kNLV ? fLVi33[pos] : 0;}
  AliTOFFormatDCS* GetLVv50(Int_t pos) const 
    {return pos<kNLV ? fLVv50[pos] : 0;}
  AliTOFFormatDCS* GetLVi50(Int_t pos) const 
    {return pos<kNLV ? fLVi50[pos] : 0;}
  AliTOFFormatDCS* GetLVv48(Int_t pos) const 
    {return pos<kNLV ? fLVv48[pos] : 0;}
  AliTOFFormatDCS* GetLVi48(Int_t pos) const 
    {return pos<kNLV ? fLVi48[pos] : 0;}
  AliTOFFormatDCS* GetFEEthr(Int_t pos) const
    {return pos<kNFEEthr ? fFEEthr[pos] : 0;}
  AliTOFFormatDCS* GetFEEtfeac(Int_t pos) const
    {return pos<kNFEEtfeac ? fFEEtfeac[pos] : 0;}
  AliTOFFormatDCS* GetFEEttrm(Int_t pos) const
    {return pos<kNFEEttrm ? fFEEttrm[pos] : 0;}

private:
  void Init();
  void Introduce(UInt_t numAlias, const TObjArray* aliasArr) const;
  void CreateHisto(int nbin);
  
  Int_t fRun;       // Run number
  Int_t fStartTime; // start time
  Int_t fEndTime;   // end time
  
  
  TString fAliasNames[kNAliases];        // aliases for DCS data
  AliTOFFormatDCS *fHVvp[kNHV];          // HV voltages, positive ch
  AliTOFFormatDCS *fHVvn[kNHV];          // HV voltages, negative ch
  AliTOFFormatDCS *fHVip[kNHV];          // HV currents, positive ch
  AliTOFFormatDCS *fHVin[kNHV];          // HV currents, negative ch
  AliTOFFormatDCS *fLVv[kNLV];           // LV fea voltages
  AliTOFFormatDCS *fLVi[kNLV];           // LV fea currents
  AliTOFFormatDCS *fLVv33[kNLV33];       // LV 3.3 V voltages
  AliTOFFormatDCS *fLVi33[kNLV33];       // LV 3.3 V currents
  AliTOFFormatDCS *fLVv50[kNLV50];       // LV 5.0 V voltages
  AliTOFFormatDCS *fLVi50[kNLV50];       // LV 5.0 V currents
  AliTOFFormatDCS *fLVv48[kNLV48];       // LV 48 V voltages
  AliTOFFormatDCS *fLVi48[kNLV48];       // LV 48 V currents
  AliTOFFormatDCS *fFEEthr[kNFEEthr];    // FEE thresholds
  AliTOFFormatDCS *fFEEtfeac[kNFEEtfeac];// FEE feac temperatures
  AliTOFFormatDCS *fFEEttrm[kNFEEttrm];  // FEE trm temperatures
  Float_t fT[3];                         // environment temperature
  Float_t fP[3];                         // environment pressure
  
  Bool_t fIsProcessed;                   // bool to know processing status
  
  ClassDef(AliTOFDataDCS, 2);
};

#endif
