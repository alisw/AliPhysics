#ifndef ALIMONITORPLOT_H
#define ALIMONITORPLOT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TString.h>
#include <Gtypes.h>


class TH1;


class AliMonitorPlot : public TNamed {
public:
  AliMonitorPlot();
  AliMonitorPlot(const AliMonitorPlot& plot);
  AliMonitorPlot& operator = (const AliMonitorPlot& plot);
  AliMonitorPlot(const char* name, const char* title);
  virtual ~AliMonitorPlot() {};

  virtual void    SetReference(TH1* ref) = 0;
  virtual void    SetReference(AliMonitorPlot* ref) = 0;
  void            SetDescription(TString description) 
                    {fDescription = description;};
  TString         GetDescription() const {return fDescription;};

  virtual void    Update() = 0;
  virtual void    Add(AliMonitorPlot* plot) = 0;
  virtual void    Reset() = 0;
  virtual void    ResetList() = 0;

  static void     SetDrawRef(Bool_t drawRef = kTRUE); // *MENU*
  static Bool_t   GetDrawRef() {return fgDrawRef;};
  static void     SetThreshold(Float_t threshold); // *MENU*
  static Float_t  GetThreshold() {return fgThreshold;};

  virtual void    DrawEvent(Int_t number = 1); // *MENU*
  virtual void    DrawSum(Int_t number); // *MENU*
  virtual void    DrawRun(); // *MENU*

  virtual Bool_t  CompareEvent(Int_t number = 1);
  virtual Bool_t  CompareSum(Int_t number);
  virtual Bool_t  CompareRun();

protected:
  virtual Bool_t  ComparePlot() = 0;
  virtual Bool_t  GetEvent(Int_t number = 1) = 0;
  virtual Bool_t  GetSum(Int_t number) = 0;
  virtual Bool_t  GetRun() = 0;
  virtual void    DrawPlot() = 0;

  TString         fDescription;  // description of the monitor histogram
  Int_t           fNumberOfEvents;// number of monitored events

  static Bool_t   fgDrawRef;     //! draw the reference histogram or not
  static Float_t  fgThreshold;   //! threshold for the comparison to the reference histogram

  static Color_t  fgColorData;   //! color of the data histogram
  static Color_t  fgColorRef;    //! color of the reference histogram
  static Color_t  fgColorCompare;//! color of the comparison histogram

  ClassDef(AliMonitorPlot, 1)   // base class for plots used to monitor the quality of the recorded data
};
 

#endif









