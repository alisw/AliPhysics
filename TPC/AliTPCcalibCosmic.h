#ifndef ALITPCCALIBCOSMIC_H
#define ALITPCCALIBCOSMIC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "AliTPCCalPad.h"
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
  void             FindPairs(AliESDEvent *event);
  Bool_t           IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1);
  void             SetGainMap(AliTPCCalPad *GainMap){fGainMap = GainMap;};
public:
  static void BinLogX(TH1 * h);   // method for correct histogram binning
  AliTPCCalPad *fGainMap;         //! gain map from Krypton calibration
  TH1F  *fHistNTracks;
  TH1F  *fClusters;
  TH2F  *fModules;
  TH1F  *fHistPt;
  TH1F  *fPtResolution;
  TH2F  *fDeDx;
  // cuts
  //
  Float_t fCutMaxD;     // maximal distance in rfi ditection
  Float_t fCutTheta;    // maximal distance in theta ditection
  Float_t fCutMinDir;   // direction vector products
private:
  AliTPCcalibCosmic(const AliTPCcalibCosmic&); 
  AliTPCcalibCosmic& operator=(const AliTPCcalibCosmic&); 

  ClassDef(AliTPCcalibCosmic, 1); 
};

#endif

