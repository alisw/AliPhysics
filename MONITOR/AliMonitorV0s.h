#ifndef ALIMONITORV0S_H
#define ALIMONITORV0S_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitor.h"


class AliMonitorV0s : public AliMonitor {
public:
  AliMonitorV0s();
  virtual ~AliMonitorV0s() {};

  virtual void     CreateHistos(TFolder* folder);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader, AliESD* esd);

private:
  AliMonitorV0s(const AliMonitorV0s& monitor);
  AliMonitorV0s& operator = (const AliMonitorV0s& monitor);

  AliMonitorHisto* fRadius;             // radius of V0 vertices
  AliMonitorHisto* fMassK0;             // invariant mass distribution of V0s for pi+ pi- hypothesis
  AliMonitorHisto* fMassLambda;         // invariant mass distribution of V0s for p pi- hypothesis
  AliMonitorHisto* fMassAntiLambda;     // invariant mass distribution of V0s for anti-p pi+ hypothesis

  ClassDef(AliMonitorV0s, 0)   // creation and filling of monitor histograms for correlations between TPC and ITS
};
 

#endif









