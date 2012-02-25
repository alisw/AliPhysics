#ifndef ALITPCMONITORCONFIG_H
#define ALITPCMONITORCONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorConfig class
////
//// Configuration handler class for AliTPCMonitor.
//// The class reads and stores basic configurations
//// for the AliTPCMonitor class. The values can be changed
//// online or written to the configuration file  AliTPCMonitorConfig.txt 
//// 
//// Author: Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TString.h>

class AliTPCMonitorConfig: public TNamed
{
 public :
    
    AliTPCMonitorConfig(const Char_t* name,const Char_t* title);
    AliTPCMonitorConfig(const  AliTPCMonitorConfig &config);
    AliTPCMonitorConfig& operator= (const AliTPCMonitorConfig& config);
    virtual ~AliTPCMonitorConfig();
    
    Float_t  GetButtonXSize()                        const { return fButtonXSize;}
    Float_t  GetButtonYSize()                        const { return fButtonYSize;}
    Float_t  GetButtonXFirst1()                      const { return fButtonFirstX1;}
    Float_t  GetButtonXFirst2()                      const { return fButtonFirstX2;}
    Float_t  GetButtonYFirst()                       const { return fButtonFirstY;}
    Int_t    GetMainXSize()                          const { return fMainXSize;}
    Int_t    GetMainYSize()                          const { return fMainYSize;}
    Int_t    GetBorderXSize()                        const { return fBorderXSize;}
    Int_t    GetBorderYSize()                        const { return fBorderYSize;}
    Int_t    GetCanvasXOffset()                      const { return fCanvasXOffset;}
    Int_t    GetCanvasXSize()                        const { return fCanvasXSize;}
    Int_t    GetCanvasYSize()                        const { return fCanvasYSize;}
    Int_t    GetCanvasXSpace()                       const { return fCanvasXSpace;}
    Int_t    GetCanvasYSpace()                       const { return fCanvasYSpace;}
    
    const Float_t* GetComponentSelection()                 const { return fComponents;}
    
    Int_t    GetEventProcessed()                     const { return fEventProcessed  ;}
    
  
    Int_t    GetFormat()                             const { return fFormat      ;}
    const Char_t*  GetFile()                         const { return fFileCurrent.Data();}
    
    Int_t    GetFitPulse()                           const { return fFitPulse     ;}
    const Char_t*  GetLastProcFile();
    Int_t    GetMaxHwAddr()                          const { return fMaxHwAddr            ; } 
    
    Int_t    GetLastSector()                         const { return fSectorLast;}
    Int_t    GetLastSectorDisplayed()                const { return fSectorLastDisplayed;}
    
    Int_t    GetNextEventID()                        const { return fEventNextID      ;}
    Int_t    GetNumOfChannels()                      const { return fNumOfChannels        ; }
    
    Int_t    GetPedestals()                          const { return fPedestals            ; }
    Int_t    GetProcNextEvent()                      const { return fEventNext;}
 
    Int_t    GetProcOneSector()                      const { return fProcOneSector;}
    
    Int_t    GetRangeBaseMin()                       const { return fRangeBaseMin;}
    Int_t    GetRangeBaseMax()                       const { return fRangeBaseMax;}
    
    Int_t    GetRangeMaxAdcMin()                     const { return fRangeMaxAdcMin;}
    Int_t    GetRangeMaxAdcMax()                     const { return fRangeMaxAdcMax;}
    
    Int_t    GetRangeSumMin()                        const { return fRangeSumMin;}
    Int_t    GetRangeSumMax()                        const { return fRangeSumMax;}


    Int_t    GetSectorFilled(Int_t sector,Int_t side)const { return fSectorArr[sector+side*18]   ;}
    Int_t    GetSectorFilled(Int_t sector)           const { return fSectorArr[sector]   ;}
    
    Int_t    GetSamplingFrequency()                  const { return  fSamplingFreq;}
    
    Int_t    GetTimeBins()                           const { return fTimeBins             ; }
    
    Int_t    GetWrite10Bit()                         const { return fWrite10Bit ;}

    

    void     SetBaseConfig(float*  ConfArr);
        
    void     SetLastProcFile(const Char_t* val);
    
    void     SetEventProcessed(Int_t val)            {        fEventProcessed=val;}
 
    void     SetFitPulse( Int_t val)                 {        fFitPulse =val;}
    
    void     SetNumOfChannels(Int_t val)             {        fNumOfChannels = val  ; }
    
    void     SetPedestals(Int_t val)                 {        fPedestals = val      ; }
    
    void     SetTimeBins(Int_t val)                  {        fTimeBins  = val      ; }
    
    void     SetFile(const Char_t* val)              { fFileCurrent=val; }
        
    void     SetLastSector(Int_t val)                { fSectorLast = val;}

    void     SetLastSectorDisplayed(Int_t val)       { fSectorLastDisplayed = val;}
    
    void     SetNextEventID(Int_t val )              {        fEventNextID = val;}
    
    void     SetProcNextEvent(Int_t val)             {        fEventNext = val ;}
  
    void     SetProcOneSector(Int_t val)             {        fProcOneSector = val ;}
    
    void     SetRangeMax( Int_t min, Int_t max)      { fRangeMaxAdcMin=min,  fRangeMaxAdcMax=max;}
    void     SetRangeBase(Int_t min, Int_t max)      { fRangeBaseMin  =min,  fRangeBaseMax  =max;}
    void     SetRangeSum( Int_t min, Int_t max)      { fRangeSumMin   =min,  fRangeSumMax  =max;}
    
    void     SetSectorFilled(Int_t sector,Int_t side){        fSectorArr[sector+side*18] =1;}
    void     SetSectorFilled(Int_t sector)           {        fSectorArr[sector] =1;}
    
    void     SetComponentSelection(float* val)       { for(Int_t i=0;i<10;i++) fComponents[i] = val[i];}
    
    void     SetFormat(Int_t val)                    {         fFormat = val;}
    
    void     SetWrite10Bit(Int_t val)                {        fWrite10Bit =val;}
    void     SetMainSize(Int_t mainx,   Int_t mainy,  
			 Int_t borderx ,Int_t bordery );
    
    void     PrintConfig();

    void     ReadConfig(const Char_t* nameconf);
    void     ResetSectorArray()                      { for(Int_t i=0;i<36; i++) fSectorArr[i]=0;}
    
 private: 
    
    // Data Format  0: DATA 1: ROOT
    Int_t    fFormat;                                                   // Format of the processed file/stream  
    
    Int_t    fSector;                                                   // Currently processed sector 
    Int_t    fSectorLast;                                               // Previously processed sector
    Int_t    fSectorLastDisplayed;                                      // Last displayed sector
    Int_t    fSectorArr[36];                                            // Array of processed sectors
    
    // Current and Last Files and Dirs
    TString  fFileLast;                                                 // Name of last processed file/stream
    Int_t    fFileLastSet ;                                             // Flag showing if last file name was set
    
    TString  fFileCurrent;                                              // Current file/stream  name
    
    Int_t    fEventNext;                                                // Process next event -> do not stay in current event                         
    Int_t    fEventNextID;                                              // Next event ID to be processed (if event id does not exist search for next existing event)                                             
    
    Int_t    fEventProcessed;                                           // Flag to show if event was read in

    // Ranges for determination of ADC max , Baseline and  ADC Sum 
    Int_t    fRangeMaxAdcMin   ;                                        // Min timebin of range to determine max.  adc value
    Int_t    fRangeMaxAdcMax   ;                                        // Max timebin of range to determine max.  adc value

    Int_t    fRangeBaseMin   ;                                          // Min timebin of range to determine baseline
    Int_t    fRangeBaseMax   ;                                          // Max timebin of range to determine basline

    Int_t    fRangeSumMin  ;                                            // Min timebin of range to determine adc sum
    Int_t    fRangeSumMax ;                                             // Max timebin of range to determine adc sum

    // Canvas Size for Monitor Canvases
    Int_t    fCanvasXSize;                                              // Canvas size in x  ( set to fCanvasMainSize )
    Int_t    fCanvasYSize;                                              // Canvas size in y  ( set to fCanvasMainSize )
    Int_t    fCanvasXSpace;                                             // Canvas size in x + border size
    Int_t    fCanvasYSpace;                                             // Canvas size in y + border size
    Int_t    fCanvasXOffset;                                            // Canvas x offset (main window)
    Int_t    fCanvasMainSize;                                           // Canvas size in x and y
    
    // Size of Main frame and Border (depending on Window Manager)     
    Int_t    fMainXSize;                                                // Main window size x
    Int_t    fMainYSize;                                                // Main window size y

    Int_t    fBorderXSize ;                                             // Canvas border size x
    Int_t    fBorderYSize ;                                             // Canvas border size y
    
    // Buttonsize;
    Float_t  fButtonXSize  ;                                            // Button size x
    Float_t  fButtonYSize  ;                                            // Button size y
    Float_t  fButtonFirstX1;                                            // Pos of first button row in x  
    Float_t  fButtonFirstX2;                                            // Pos of second button row in x
    Float_t  fButtonFirstY ;                                            // Position of first button in y
    
    Int_t    fWrite10Bit ;                                              // Flag to write 10 bit data words to file
    
    // Arr to Store Selected components to be displayed
    Float_t  fComponents[10];                                           // Array of components to be selected for display
    
    // Sampling Freq required for FFT 
    Int_t    fSamplingFreq;                                             // Sampling frequency for data taking
    
    
    Int_t    fPedestals     ;                                           // Version for pedestal calculation  
    Int_t    fNumOfChannels ;                                           // Maximum number of channels
    Int_t    fTimeBins      ;                                           // Number of timebins to be displayed
    Int_t    fMaxHwAddr     ;                                           // Max value of hardware addresses

    Int_t    fFitPulse     ;                                            // Flag for fitting pulse around max adc    

    Int_t    fProcOneSector ;                                           // Flag for processing only the specified sector for the next event

 
    ClassDef(AliTPCMonitorConfig,1);  
};
#endif 
