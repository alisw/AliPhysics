#ifndef ALIMONITORTREND_H
#define ALIMONITORTREND_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitorPlot.h"
#include "TArrayD.h"


class AliMonitorTrend : public AliMonitorPlot {
public:
  AliMonitorTrend();
  AliMonitorTrend(const AliMonitorTrend& trend);
  AliMonitorTrend& operator =(const AliMonitorTrend& trend);
  AliMonitorTrend(const char* name, const char* title,
		  const char* label, Double_t min = 0, Double_t max = 0);
  virtual ~AliMonitorTrend();

  virtual void    SetReference(TH1* ref);
  virtual void    SetReference(AliMonitorPlot* ref);

  void            Fill(Double_t x);

  virtual void    Update();
  virtual void    Add(AliMonitorPlot* plot);
  virtual void    Reset();
  virtual void    ResetList();

  Double_t        GetMean() const;
  Double_t        GetSigma() const;
protected:
  virtual Bool_t  ComparePlot();
  virtual Bool_t  GetEvent(Int_t number = 1);
  virtual Bool_t  GetSum(Int_t number);
  virtual Bool_t  GetRun();
  virtual void    DrawPlot();

  TH1*            CreateHisto(Int_t nBins);

  TString         fLabel;        // label of the y axis
  Double_t        fMin;          // minimal y axis value
  Double_t        fMax;          // maximal y axis value
  TArrayD         fData;         // monitored values
  static Int_t    fgIncSize;  //! amount by which the TArrayD is increased

  TH1*            fHistoDraw;    //! the histogram for the trend, used for comparison to a reference
  Double_t        fRefMean;      //! mean reference value
  Double_t        fRefSigma;     //! standard deviation of the reference value
  TH1*            fHistoCompare; //! the result of the comparison to the reference, only bins with large deviation are set

  ClassDef(AliMonitorTrend, 1)   // histogram for monitoring the quality of the recorded data (time dependent value)
};
 

#endif









