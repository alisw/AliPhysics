#ifndef ALIMONITORHISTO_H
#define ALIMONITORHISTO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitorPlot.h"
#include <TH1.h>
#include <TList.h>
#include <TString.h>


class AliMonitorHisto : public AliMonitorPlot {
public:
  enum ENorm {kNormNone, kNormEvents, kNormEntries, kNormIntegral};

  AliMonitorHisto();
  AliMonitorHisto(const AliMonitorHisto& histo);
  AliMonitorHisto(TH1* histo, ENorm norm = kNormIntegral);
  virtual ~AliMonitorHisto();

  static void     SetNHistosMax(Int_t nHistosMax) {fgNHistosMax = nHistosMax;};
  static Int_t    GetNHistosMax() {return fgNHistosMax;};

  virtual void    SetReference(TH1* ref);
  virtual void    SetReference(AliMonitorPlot* ref);

  void            Fill(Axis_t x);
  void            Fill(Axis_t x, Axis_t y);
  void            Fill(Axis_t x, Axis_t y, Stat_t w);
  void            ScaleErrorBy(Double_t factor);

  virtual void    Update();
  virtual void    Update(TH1* histo);
  virtual void    Add(AliMonitorPlot* plot);
  virtual void    Reset();
  virtual void    ResetList();

protected:
  void            Scale(Int_t nEvents);
  virtual Bool_t  ComparePlot();
  virtual Bool_t  GetEvent(Int_t number = 1);
  virtual Bool_t  GetSum(Int_t number);
  virtual Bool_t  GetRun();
  virtual void    DrawPlot();

  TH1*            fHisto;        //! the histogram that is currently filled with data
  TList           fHistoList;    // list of histograms of last fNHistos events
  Int_t           fNHistos;      // number of buffered histograms
  static Int_t    fgNHistosMax;  //! maximal number of buffered histograms
  TH1*            fHistoRun;     // sum of histograms for the current run

  TH1*            fHistoDraw;    //! the normalized histogram, used for comparison to a reference histogram
  TH1*            fHistoRef;     //! the reference histogram for comparison
  TH1*            fHistoCompare; //! the result of the comparison to the reference histogram, only bins with large deviation are set

  ENorm           fNorm;         // type of normalization

  ClassDef(AliMonitorHisto, 2)   // histogram for monitoring the quality of the recorded data
};
 

#endif









