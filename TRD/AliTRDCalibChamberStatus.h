#ifndef ALITRDCALIBCHAMBERSTATUS_H
#define ALITRDCALIBCHAMBERSTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalibChamberStatus.h 34340 2009-08-20 07:48:28Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for online calibration                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#include <TCanvas.h>
#include <TH2.h>
#endif

class AliRawReader;

class AliTRDCalChamberStatus;
class AliRawReader;
class AliTRDCalDCS;


struct eventHeaderStruct;

class AliTRDCalibChamberStatus : public TObject {

public:

  AliTRDCalibChamberStatus();
  AliTRDCalibChamberStatus(const AliTRDCalibChamberStatus &ped);
  virtual ~AliTRDCalibChamberStatus();

  AliTRDCalibChamberStatus& operator = (const  AliTRDCalibChamberStatus &source);

  void ProcessEvent(AliRawReader    *rawReader, Int_t nevents_physics);
  void ProcessEvent3(AliRawReader    *rawReader, Int_t nevents_physics);
  
  void Init();
  void AnalyseHisto();
  void CheckEORStatus(AliTRDCalDCS *calDCS);

  void Add(AliTRDCalibChamberStatus *calibChamberStatus);

  Int_t GetNumberEventNotEmpty() const { return fCounterEventNotEmpty; };
  
  THnSparseI *GetSparseI()       const {return fHnSparseI;};
  THnSparseI *GetSparseHCM()     const {return fHnSparseHCM;};
  // for fDebugLevel>0
  THnSparseI *GetSparseEvtDet()  const {return fHnSparseEvtDet;};
  THnSparseI *GetSparseDebug()   const {return fHnSparseDebug;};
  THnSparseI *GetSparseMCM()     const {return fHnSparseMCM;};
  

  AliTRDCalChamberStatus *GetCalChamberStatus() const {return fCalChamberStatus;};

  void  DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
  
  Bool_t TestEventHisto(Int_t nevent);

  // Plot
  TH2D *PlotSparseI(Int_t sm, Int_t side);    // Plot fStatus for sm 
  TH2F *MakeHisto2DSmPlEORStatus(AliTRDCalDCS *calDCS, Int_t sm, Int_t pl);
  TCanvas *PlotHistos2DSmEORStatus(AliTRDCalDCS *calDCS,Int_t sm, const Char_t *name);

  // Debug
  void     SetDebugLevel(Short_t level)  { fDebugLevel = level;   }

 private:

  Int_t fDetector;                           //  Current detector
  Int_t fNumberOfTimeBins;                   //  Current number of time bins
  Int_t fCounterEventNotEmpty;               //  Counter Events Not Empty
  
  AliTRDCalChamberStatus *fCalChamberStatus; //  AliTRDCalChamberStatus result
  
  THnSparseI *fHnSparseI;                    //  THnSparse for entries in half chambers
  THnSparseI *fHnSparseHCM;                  //  THnSparse for DCS half chamber status
  
  // for fDebugLevel>0
  THnSparseI *fHnSparseEvtDet;    //  THnSparse for entries in half chambers per events
  THnSparseI *fHnSparseDebug;     //  THnSparse for half chambers satuts
  THnSparseI *fHnSparseMCM;       //  THnSparse for DCS MCM status

  TCanvas *fC1;

  Short_t     fDebugLevel;                   // Flag for debugging

  ClassDef(AliTRDCalibChamberStatus,1)
    
};
#endif