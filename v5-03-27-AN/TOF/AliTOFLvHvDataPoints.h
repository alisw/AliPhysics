#ifndef ALITOFLVHVDATAPOINTS_H
#define ALITOFLVHVDATAPOINTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:  $ */

/////////////////////////////////////////////////////////
//                                                     //
// AliTOFLvHvDataPoints class                          //
// main aim to introduce                               //
// the aliases for the TOF LV and HV DCS data points   //
// to be then stored in the OCDB, and to process them. //
//                                                     //
/////////////////////////////////////////////////////////

#include "TObject.h"

class TMap;
class TClonesArray;
class TString;
class TH1C;

class AliTOFDCSmaps;

class AliTOFLvHvDataPoints : public TObject {
public:
  enum {kNsectors=18, kNplates=5, kNddl=72, kNpads=18*91*96, kNmaxDataPoints=77777};
  
  AliTOFLvHvDataPoints();
  AliTOFLvHvDataPoints(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery );
  AliTOFLvHvDataPoints(const AliTOFLvHvDataPoints & data);
  AliTOFLvHvDataPoints& operator=(const AliTOFLvHvDataPoints & data);
  ~AliTOFLvHvDataPoints();
  
  void SetRun(Int_t run) {fRun = run;}
  void SetStartTime(Int_t startTime) {fStartTime = startTime;}
  void SetEndTime(Int_t endTime) {fEndTime = endTime;}
  void SetStartTimeDCSQuery(Int_t startTimeDCSQuery) {fStartTimeDCSQuery = startTimeDCSQuery;}
  void SetEndTimeDCSQuery(Int_t endTimeDCSQuery) {fEndTimeDCSQuery = endTimeDCSQuery;}
  void SetNSecondsBeforeEOR(Int_t nSecondsBeforeEOR) {fNSecondsBeforeEOR = nSecondsBeforeEOR;}
  Int_t GetRun() const {return fRun;}
  Int_t GetStartTime() const {return fStartTime;}
  Int_t GetEndTime() const {return fEndTime;}
  Int_t GetStartTimeDCSQuery() const {return fStartTimeDCSQuery;}
  Int_t GetEndTimeDCSQuery() const {return fEndTimeDCSQuery;}
  Int_t GetNSecondsBeforeEOR() const {return fNSecondsBeforeEOR;}
  Bool_t ProcessData(TMap& aliasMap);

  
  const char* GetAliasNameXLV(Int_t pos) const 
    {return pos<kNddl ? fAliasNamesXLVmap[pos].Data() : 0;}
  
  const char* GetAliasNameXHV(Int_t pos1, Int_t pos2) const 
    {return pos1<kNsectors&&pos2<kNplates ? fAliasNamesXHVmap[pos1][pos2].Data() : 0;}
  
  void Draw(const Option_t* /*option*/);
  void DrawHVandLVMap(Int_t index);
  void DrawHVMap(Int_t index);
  void DrawLVMap(Int_t index);
  
  void SetFDRFlag(Bool_t flag) {fFDR = flag;}
  Bool_t GetFDRFlag() const {return fFDR;}

  Int_t GetNumberOfHVandLVmaps() const { return fNumberOfHVandLVmaps; };
  AliTOFDCSmaps * GetHVandLVmapAtSOR() const { return fMap[0]; };
  AliTOFDCSmaps * GetHVandLVmapAtEOR() ;
  AliTOFDCSmaps * GetHVandLVmap(Int_t index) const { if (index>=fNumberOfHVandLVmaps) return 0x0; else return fMap[index]; };
  Int_t GetNumberOfLVmaps() const { return fNumberOfLVdataPoints; };
  AliTOFDCSmaps * GetLVmap(Int_t index) const { if (index>=fNumberOfLVdataPoints) return 0x0; else return fLVDataPoints[index]; };
  Int_t GetNumberOfHVmaps() const { return fNumberOfHVdataPoints; };
  AliTOFDCSmaps * GetHVmap(Int_t index) const { if (index>=fNumberOfHVdataPoints) return 0x0; else return fHVDataPoints[index]; };

private:
  void Init();
  void CreateHisto(int nbin);
  void FillHVarrayPerDataPoint(Int_t sector, Int_t plate, UInt_t baseWord, Short_t *array) const;
  void FillLVarrayPerDataPoint(Int_t nDDL, UInt_t baseWord, Short_t *array) const;
  void GetStripsConnectedToFEAC(Int_t nDDL, Int_t nFEAC, Int_t *iStrip, Int_t &firstPadX, Int_t &lastPadX) const;

  Bool_t ReadHVDataPoints(TMap& aliasMap);
  Bool_t ReadLVDataPoints(TMap& aliasMap);

  Bool_t MergeMaps();
  Bool_t MergeHVmap();
  Bool_t MergeLVmap();

  Int_t InsertHVDataPoint(AliTOFDCSmaps *object);
  Int_t FindHVdpIndex(Int_t z) const;
  Int_t InsertLVDataPoint(AliTOFDCSmaps *object);
  Int_t FindLVdpIndex(Int_t z) const;
  Int_t fRun;       // Run number
  Int_t fStartTime; // start time
  Int_t fEndTime;   // end time  
  Int_t fStartTimeDCSQuery; // start time DCSQuery
  Int_t fEndTimeDCSQuery;   // end time DCSQuery
  
  Bool_t fIsProcessed; // bool to know processing status
  Bool_t fFDR;         // bool to know whether we are in a FDR run

  AliTOFDCSmaps *fLVDataPoints[kNmaxDataPoints]; // LV status map VS time
  AliTOFDCSmaps *fHVDataPoints[kNmaxDataPoints]; // HV status map VS time
  AliTOFDCSmaps *fMap[kNmaxDataPoints]; // LV&&HV status map VS time

  Int_t fNumberOfLVdataPoints; // number of found LV status dps
  Int_t fNumberOfHVdataPoints; // number of found HV status dps
  Int_t fNumberOfHVandLVmaps; // number of found LV&&HV status maps
  
  TString fAliasNamesXLVmap[kNddl]; // aliases for LV map
  TString fAliasNamesXHVmap[kNsectors][kNplates]; // aliases for HV map

  AliTOFDCSmaps *fStartingLVmap; // starting value for LV map
  AliTOFDCSmaps *fStartingHVmap; // starting value for HV map
  
  TH1C * fHisto; // histogram

  Int_t fNSecondsBeforeEOR; // time window to choose if a run ended correctly or not [s]

  ClassDef(AliTOFLvHvDataPoints, 3);
};

#endif
