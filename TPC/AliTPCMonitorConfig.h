#ifndef ALITPCMONITORCONFIG_H
#define ALITPCMONITORCONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorConfig class
//
// Configuration handler class for AliTPCMonitor
// 
// Author: Stefan Kniege, IKF, Frankfurt
//       
//
/////////////////////////////////////////////////////////////////////////



#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include "TNamed.h"
#include "TObject.h"
#include "TSystem.h" 
#include "AliLog.h" 
using namespace std;

class AliTPCMonitorConfig: public TNamed
{
 public :
    
    AliTPCMonitorConfig(Char_t* name,Char_t* title);
    virtual ~AliTPCMonitorConfig();
    
    // Data Format  0: DATA 1: ROOT
    Int_t    fFormat;
    
    Int_t    fSector;                                                   // Currently processed sector 
    Int_t    fSectorLast;                                               // Previously processed sector
    Int_t    fSectorLastDisplayed;                                      // Last displayed sector
    Int_t*   fSectorArr;                                                // Array of processed sectors
    
    // Current and Last Files and Dirs
    Char_t*  fFileLast;                                                 // Name of last processed file/stream
    Int_t    fFileLastSet ;                                             // Flag showing if last file name was set
    
    Char_t*  fFileCurrent;                                              // Current file/stream  name
    
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
    
    // Size of Main frame and Border (depending on Window Manager)      // Main window size x
    Int_t    fMainXSize;                                                // Main window size y
    Int_t    fMainYSize;

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
    Float_t* fComponents;                                               // Array of components to be selected for display
    
    // Sampling Freq required for FFT 
    Int_t    fSamplingFreq;                                             // Sampling frequency for data taking
    
    
    Int_t    fPedestals     ;                                           // Version for pedestal calculation  
    Int_t    fNumOfChannels ;                                           // Maximum number of channels
    Int_t    fTimeBins      ;                                           // Number of timebins to be displayed
    Int_t    fMaxHwAddr     ;                                           // Max value of hardware addresses

    Int_t    fFitPulse     ;                                            // Flag for fitting pulse around max adc    

    Int_t    fProcOneSector ;
 
    
    Float_t  GetButtonXSize()                        { return fButtonXSize;}
    Float_t  GetButtonYSize()                        { return fButtonYSize;}
    Float_t  GetButtonXFirst1()                      { return fButtonFirstX1;}
    Float_t  GetButtonXFirst2()                      { return fButtonFirstX2;}
    Float_t  GetButtonYFirst()                       { return fButtonFirstY;}
    Int_t    GetMainXSize()                          { return fMainXSize;}
    Int_t    GetMainYSize()                          { return fMainYSize;}
    Int_t    GetBorderXSize()                        { return fBorderXSize;}
    Int_t    GetBorderYSize()                        { return fBorderYSize;}
    Int_t    GetCanvasXOffset()                      { return fCanvasXOffset;}
    Int_t    GetCanvasXSize()                        { return fCanvasXSize;}
    Int_t    GetCanvasYSize()                        { return fCanvasYSize;}
    Int_t    GetCanvasXSpace()                       { return fCanvasXSpace;}
    Int_t    GetCanvasYSpace()                       { return fCanvasYSpace;}
    
    Float_t* GetComponentSelection()                 { return fComponents;}
    
    Int_t    GetEventProcessed()                     { return fEventProcessed  ;}
    
  
    Int_t    GetFormat()                             { return fFormat      ;}
    Char_t*  GetFile()                               { return fFileCurrent;}
    
    Int_t    GetFitPulse()                           { return fFitPulse     ;}
    Char_t*  GetLastProcFile();
    Int_t    GetMaxHwAddr()                          { return fMaxHwAddr            ; } 
    
    Int_t    GetLastSector()                         { return fSectorLast;}
    Int_t    GetLastSectorDisplayed()                { return fSectorLastDisplayed;}
    
    Int_t    GetNextEventID()                        { return fEventNextID      ;}
    Int_t    GetNumOfChannels()                      { return fNumOfChannels        ; }
    
    Int_t    GetPedestals()                          { return fPedestals            ; }
    Int_t    GetProcNextEvent()                      { return fEventNext;}
 
    Int_t    GetProcOneSector()                      { return fProcOneSector;}
    
    Int_t    GetRangeBaseMin()                       { return fRangeBaseMin;}
    Int_t    GetRangeBaseMax()                       { return fRangeBaseMax;}
    
    Int_t    GetRangeMaxAdcMin()                     { return fRangeMaxAdcMin;}
    Int_t    GetRangeMaxAdcMax()                     { return fRangeMaxAdcMax;}
    
    Int_t    GetRangeSumMin()                        { return fRangeSumMin;}
    Int_t    GetRangeSumMax()                        { return fRangeSumMax;}


    Int_t    GetSectorFilled(Int_t sector,Int_t side){ return fSectorArr[sector+side*18]   ;}
    Int_t    GetSectorFilled(Int_t sector)           { return fSectorArr[sector]   ;}
    
    Int_t    GetSamplingFrequency()                  { return  fSamplingFreq;}
    
    Int_t    GetTimeBins()                           { return fTimeBins             ; }
    
    Int_t    GetWrite10Bit()                         { return fWrite10Bit ;}

    

    void     SetBaseConfig(float*  ConfArr);
        
    void     SetLastProcFile(Char_t* val);
    
    void     SetEventProcessed(Int_t val)            {        fEventProcessed=val;}
 
    void     SetFitPulse( Int_t val)                 {        fFitPulse =val;}
    
    void     SetNumOfChannels(Int_t val)             {        fNumOfChannels = val  ; }
    
    void     SetPedestals(Int_t val)                 {        fPedestals = val      ; }
    
    void     SetTimeBins(Int_t val)                  {        fTimeBins  = val      ; }
    
    void     SetFile(Char_t* val)                      { sprintf(fFileCurrent,val);}
        
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

    void     ReadConfig(Char_t* nameconf);
    void     ResetSectorArray()                      { for(Int_t i=0;i<36; i++) fSectorArr[i]=0;}
    
 
    ClassDef(AliTPCMonitorConfig,1);  
};
#endif 
