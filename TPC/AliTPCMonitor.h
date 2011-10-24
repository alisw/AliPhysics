
#ifndef ALITPCMONITOR_H
#define ALITPCMONITOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitor class
//// 
//// Main class for the TPC raw data Monitor.
//// The Monitor can handle rootified data, files and online streams in DATE format.
//// The monitor GUI is started by the macro TPCMonitor.C
//// 
//// Author: Stefan Kniege, IKF, Frankfurt=
////         Jens Wiechula, Uni Tuebingen (Jens.Wiechula@cern.ch)
////
/////////////////////////////////////////////////////////////////////////

#include "AliTPCMonitorConfig.h"
#include "TString.h"

class TH1F;
class TH1D;
class TH1;
class TH2F;
class TH2S;
class TCanvas;
class TH3S;
class AliTPCMonitorMappingHandler;
class AliTPCMonitorFFT;
class AliTPCMonitorConfig;
class AliRawReader;
class AliAltroRawStreamV3;

class AliTPCMonitor : public AliTPCMonitorConfig {
  
public:
  
  AliTPCMonitor(const char* name, const char* title);
  AliTPCMonitor(const  AliTPCMonitor &monitor);
  AliTPCMonitor& operator= (const AliTPCMonitor& monitor);
  
  virtual ~AliTPCMonitor();
  
  Int_t         CheckEqId(Int_t secid, Int_t eqid);
  TCanvas*      CreateCanvas(const char* name);
  void          CreateHistos();
  
  void          DeleteHistos();
  void          DisableFit(Int_t val) { fDisableFit =val; }
  void          DrawHists(Int_t histos);
  void          DrawRMSMap();
  void          DumpHeader(AliRawReader*          reader  ) const ;
  
  void          ExecPad() ;
  void          ExecRow() ;
  Int_t         ExecProcess();
  void          ExecTransform();
  
  void          FillGlobal(Int_t sector);
  void          FillHistsPadPlane();
  
  static double Gamma4(const double* x, const double* par);
  Int_t         GetChannelsProc()  const  { return fChannelIter     ;}
  Int_t         GetEventID()       const  { return fEventNumber     ;}
  TH1*          GetHisto(char* histname);
  Int_t         GetRCUPatch(Int_t runid, Int_t eqid) const;
  
  Int_t         GetPadAtX(Float_t xval, Int_t row, Int_t padmax) const ;
  Int_t         GetPadAtX(Float_t xval, Int_t row) const ;
  void          GetXY( Double_t& xval , Double_t& yval , Int_t padmax, Int_t row , Int_t pad) const ;
  
  Int_t         IsLastEvent()      const  { return fLastEv          ;}
  
  Int_t         ProcessEvent();
  
  void          SetEventID(Int_t val)     { fEventNumber =val;}
  void          SetMirror(Int_t val)      { fMirror=val;}
  void          SetVerbose(Int_t val)     { fVerb = val;}
  void          SetMappingHandler(AliTPCMonitorMappingHandler* val ) { fMapHand  = val;}
  void          ShowSel(const Int_t* compval);
  void          SetEqIds();
  
  void          ResizeCanv();
  void          ResetHistos();
  void          ResetArrays();
  Int_t         ReadDataNew(    Int_t secid);
    
  void          WriteHistos() ;
  void          Write10bitChannel();

  void          SetupMonitoringTable(const char* table);
  
 private:
  
    // stats for size of arrays and histograms /////////////////////////////////////////////////
  Int_t**      fPad;                      // array to store channel adc in time
  Float_t*     fPadMapHw;                 // array to store mapping of hardware address and channel number
  Int_t**      fPadMapRCU;                // store multiple information for hardware address for debugging purpose (read out errors occur)
  
  
    // histograms to be used      ////////////////////////////////////////////////////////////////
  TH2F*        fHistIROC;                 // histo for max adc IROC
  TH2F*        fHistOROC;                 // histo for max adc OROC
  TH2S*        fHistIROCIndex;            // histo for channel number in each bin IROC
  TH2S*        fHistOROCIndex;            // histo for channel number in each bin OROC
  TH2F*        fHistIROCTime;             // histo for peaking time in each bin IROC
  TH2F*        fHistOROCTime;             // histo for peaking time in each bin OROC
  TH2F*        fHistIROCClone;            // clone histo for max adc IROC (backup when selecting component)
  TH2F*        fHistOROCClone;            // clone histo for max adc OROC (backup when selecting component)
  TH2F*        fHistIROCRMS;              // histo for RMS IROC
  TH2F*        fHistOROCRMS;              // histo for RMS OROC
  TH2F*        fHistIROCBASE;             // histo for baseline IROC
  TH2F*        fHistOROCBASE;             // histo for baseline OROC
  TH2F*        fHistIROCSUM;              // histo for adc sum IROC
  TH2F*        fHistOROCSUM;              // histo for adc sum OROC
  
  TH2F*        fHistChannelTime;          // histo for adc(channel,time)
  TH1F*        fHistAddrMapIndex;         // histo for channel(hardware address)
  TH1F*        fHistAddrMaxAdc;           // histo for max-adc(hardware address)
  TH1F*        fHistAddrBaseMean;         // histo for baseline(hardware address)
  TH1F*        fHistAddrMaxAdcX;          // histo for max-adc-xposition(hardware address)
  TH1F*        fHistAddrAdcSum;           // histo for adc-sum(hardware address)
  TH1F*        fHistAddrBaseRms;          // histo for baselinr-rms(hardware address)
  
  TH1F*        fHistDistrSumIROC;         // distribution of adc sum for all channels IROC
  TH1F*        fHistDistrMaxIROC;         // distribution of adc max for all channels OROC
  TH1F*        fHistDistrSumOROC;         // distribution of adc sum for all channels OROC
  TH1F*        fHistDistrMaxOROC;         // distribution of adc max for all channels OROC
  
  TH2F*        fHistDistrBase2dIROC;      // distribution of baseline vs rms for all channels IROC
  TH2F*        fHistDistrBase2dOROC;      // distribution of baseline vs rms for all channels OROC
  TH1D*        fHistDistrBaseRmsIROC;     // projection of fHistDistrBase2dIROC on Rms
  TH1D*        fHistDistrBaseMeanIROC;    // projection of fHistDistrBase2dIROC on Mean
  TH1D*        fHistDistrBaseRmsOROC;     // projection of fHistDistrBase2dOROC on Rms
  TH1D*        fHistDistrBaseMeanOROC;    // projection of fHistDistrBase2dOROC on Mean
  
  TH2S*        fHistGlobalMaxA;           // global histogramm for max adc of all sectors Side A
  TH2S*        fHistGlobalMaxC;           // global histogramm for max adc of all sectors Side C
  
  TObjArray*   fHistList;                 // array to store all histogram
  
  
    // row and pad settings
  Int_t        fkNRowsIroc;                // number of rows in IROC
  Int_t        fkNRowsOroc;                // number of rows in OROC
  
  Int_t        fkNPadsIroc;                // number of pads in IROC
  Int_t        fkNPadsOroc;                // number of pads in OROC
  
  Int_t        fkNPadMinIroc;              // min for pad (y-axis) representation in 2D histogram IROC
  Int_t        fkNPadMinOroc;              // min for pad (y-axis) representation in 2D histogram OROC
  Int_t        fkNPadMaxIroc;              // max for pad (y-axis) representation in 2D histogram IROC
  Int_t        fkNPadMaxOroc;              // max for pad (y-axis) representation in 2D histogram IROC
  
  Int_t        fVerb;                     // verbose flag
  
  Int_t        fLastEv;                   // flag for last event
  
  Int_t        fEventNumber;              // current event number (ID)
  Int_t        fEventNumberOld;           // previous event number
  
  Int_t        fDisableFit;               // flag to disable fit of maximum  peak
  
  Int_t        fExecGlob;                 // flag to add Executable to global histo
  Int_t        fExecPlaneMax;             // flag to add executable for plane view ( max iroc/oroc)
  Int_t        fExecPadIrocRms;           // flag to add executable for rms plane view IROC
  Int_t        fExecPadOrocRms;           // flag to add executable for rms plane view OROC
  
  Int_t        fRunId;                    // current run id
  Int_t        fEqId;                     // currend equipment id
  
  Int_t        fPadUsedRoc;               // last ROC canvas the mouse pointed at
  Int_t        fPadUsedHwAddr;            // hwaddress for last pad the mouse pointed at
  
  Int_t        fGdcId;                    // current GDC id
  Int_t        fLdcId;                    // current LDC id
  Int_t        fLdcIdOld;                 // previous LDC id
  
  Int_t**      fMapEqidsSec;              // mapping rcu-patch(equipmentid)
  Int_t*       fMapEqidsRcu;              // mapping equipmentid(sector,rcu-patch)
  
  Int_t        fMirror;                   // mirror x position on C-side
  Int_t        fChannelIter;              // counter for channels read
  
  
  AliTPCMonitorMappingHandler*  fMapHand;            // mapping handler
  
  AliRawReader*                 fRawReader;         // reader for ROOT format
  
  //for monitoring table
  const Char_t **fkMonTable;                              //! table to pass to the raw reader
  TString  fMonTableString;                        //! string that keep the table definition
  TObjArray *fMonTableArray;                       //! tokenized array of the table string
  Bool_t   fMonTableChanged;                       //! indicates that the table changed and a new raw reader instance is needed.
  
  static const Int_t       fgkHwMaskFEC            ;                          // mask for fec in hardware address
  static const Int_t       fgkHwMaskBranch         ;                          // mask for branch in hardware address
  static const Int_t       fgkHwMaskFECChannel     ;                          // mask for fec channel  in hardware address
  static const Int_t       fgkHwMaskAltroChannel   ;                          // mask for altro channel in hardware address
  static const Int_t       fgkHwMaskAltroChip      ;                          // mask for altro chip  in hardware address
  static const Int_t       fgkHwMaskRCU            ;                          // not part of the trailer added afterwards
  
  
  ClassDef(AliTPCMonitor,1);
}; 

#endif
