#ifndef ALIMONITOR_H
#define ALIMONITOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TFolder.h>
#include <TTree.h>
#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliMonitorHisto.h"
#include "AliMonitorTrend.h"


class AliMonitor : public TObject {
public:
  virtual void     CreateHistos(TFolder* folder) = 0;
  virtual void     CreateBranches(TTree* tree);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader) = 0;

protected:
  TFolder*         fFolder;    // sub folder for monitor histograms

  AliMonitorHisto* CreateHisto1(const char* name, const char* title,
				Int_t xBins, Double_t xMin, Double_t xMax,
				const char* xTitle, const char* yTitle,
				AliMonitorHisto::ENorm norm);
  AliMonitorHisto* CreateHisto2(const char* name, const char* title,
				Int_t xBins, Double_t xMin, Double_t xMax,
				Int_t yBins, Double_t yMin, Double_t yMax,
				const char* xTitle, const char* yTitle,
				const char* zTitle,
				AliMonitorHisto::ENorm norm);
  AliMonitorTrend* CreateTrend(const char* name, const char* title,
			       const char* label, 
			       Double_t min = 0, Double_t max = 0);

  ClassDef(AliMonitor, 0)   // base class for the creation and filling of monitor histograms
};
 

#endif









