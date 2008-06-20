#ifndef ALITPCCALIBCOSMIC_H
#define ALITPCCALIBCOSMIC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
class TH1F;
class TH2F;
class TList;
class AliESDEvent;

#include "TTreeStream.h"


class AliTPCcalibCosmic:public AliTPCcalibBase {
public:
  AliTPCcalibCosmic(); 
  AliTPCcalibCosmic(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibCosmic();
  
  virtual void     Process(AliESDEvent *event);
  virtual Long64_t Merge(TCollection *li);
  
  
public:
  static void BinLogX(TH1 * h);   // method for correct histogram binning

  TH1F  *fHistNTracks;
  TH1F  *fClusters;
  TH2F  *fModules;
  TH1F  *fHistPt;
  TH1F  *fPtResolution;
  TH2F  *fDeDx;

  AliTPCcalibCosmic(const AliTPCcalibCosmic&); 
  AliTPCcalibCosmic& operator=(const AliTPCcalibCosmic&); 

  ClassDef(AliTPCcalibCosmic, 1); 
};

#endif

