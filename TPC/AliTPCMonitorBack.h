
#ifndef ALITPCMONITOR_H
#define ALITPCMONITOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitor class
// 
// Main class for TPC Monitor
// Monitor can handle rootified data, files and online streams in DATE format.
// The monitor GUI is started by the macro TPCMonitor.C
// 
// Author: Stefan Kniege, IKF, Frankfurt
//       
//
/////////////////////////////////////////////////////////////////////////

 
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <istream>
#include <ostream>
#include "TStyle.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TFormula.h"
#include <string>
#include "TROOT.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TGButtonGroup.h"
#include "TGButton.h"
#include "TGTextBuffer.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TH3S.h"
#include "AliTPCMonitorMappingHandler.h"
#include "AliTPCMonitorDateFile.h"
#include "AliTPCMonitorDateMonitor.h"
#include "AliTPCMonitorDateFormat.h"
#include "AliTPCMonitorAltro.h"
#include "AliTPCMonitorFFT.h"
#include "AliTPCMonitorConfig.h"
#include "AliSignalProcesor.h"
#include "AliRawReaderRoot.h"
#include "AliRawReader.h"
#include "TGMsgBox.h"
#include "TNamed.h"
#include "TObject.h" 
#include "TDirectory.h"
#include "AliLog.h"
#include "RQ_OBJECT.h"

using namespace std;

class AliTPCMonitor : public AliTPCMonitorConfig {
    
 public:
    
    AliTPCMonitor(char* name, char* title);
    virtual ~AliTPCMonitor();
    
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
    Int_t        kNRowsIroc;                // number of rows in IROC
    Int_t        kNRowsOroc;                // number of rows in OROC
    
    Int_t        kNPadsIroc;                // number of pads in IROC
    Int_t        kNPadsOroc;                // number of pads in OROC
    
    Int_t        kNPadMinIroc;              // min for pad (y-axis) representation in 2D histogram IROC
    Int_t        kNPadMinOroc;              // min for pad (y-axis) representation in 2D histogram OROC
    Int_t        kNPadMaxIroc;              // max for pad (y-axis) representation in 2D histogram IROC
    Int_t        kNPadMaxOroc;              // max for pad (y-axis) representation in 2D histogram IROC

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
    
    AliRawReaderRoot*             fReaderROOT;         // reader for ROOT format
    AliTPCMonitorDateFile*        fReaderDATE;         // reader for DATE files
    AliTPCMonitorDateMonitor*     fReaderDATEMon;      // reader for DATE monitoring
    

    Int_t         CheckEqId(Int_t secid, Int_t eqid);
    TCanvas*      CreateCanvas(char* name);
    void          CreateHistos();
    
    void          DeleteHistos();
    void          DisableFit(Int_t val) { fDisableFit =val; }
    void          DrawHists(Int_t histos);
    void          DrawRMSMap();
    void          DumpHeader(AliRawReaderRoot*          reader  );
    void          DumpHeader(AliTPCMonitorDateFormat*   DateForm);
    
    void          ExecPad() ;
    void          ExecRow() ;
    Int_t         ExecProcess();
    void          ExecTransform();
    
    void          FillGlobal(Int_t sector);
    void          FillHistsDecode( AliTPCMonitorAltro* altro , Int_t rcu_patch, Int_t id=0);
    void          FillHistsPadPlane();
    
    static double Gamma4(double* x, double* par);
    Int_t         GetChannelsProc()   { return fChannelIter     ;}
    Int_t         GetEventID()        { return fEventNumber     ;}
    TH1*          GetHisto(char* histname);
    Int_t         GetRCUPatch(Int_t runid, Int_t eqid);
    
    Int_t         GetPadAtX(Float_t xval, Int_t row, Int_t padmax);
    Int_t         GetPadAtX(Float_t xval, Int_t row);
    void          GetXY( double& xval , double& yval , Int_t rowmax, Int_t row , Int_t pad);
    
    Int_t         IsLastEvent()       { return fLastEv          ;}
    
    Int_t         ProcessEvent(); 
    
    void          SetEventID(Int_t val)     { fEventNumber =val;}
    void          SetMirror(Int_t val)      { fMirror=val;}
    void          SetVerbose(Int_t val)     { fVerb = val;}
    void          SetMappingHandler(AliTPCMonitorMappingHandler* val ) { fMapHand  = val;}
    void          ShowSel(Int_t* comp_val);
    void          SetEqIds();
    
    void          ResizeCanv();
    void          ResetHistos();
    void          ResetArrays();
    Int_t         ReadData(    Int_t secid);
    Int_t         ReadDataDATE(Int_t secid ,Int_t format); 
    Int_t         ReadDataROOT(Int_t secid ); 
    
    void          WriteHistos() ;
    void          Write10bitChannel();
    

 
   
    
    
 private:
     
    ClassDef(AliTPCMonitor,1); 
}; 

#endif
