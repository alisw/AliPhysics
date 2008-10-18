#ifndef ALITPCCALIBTIME_H
#define ALITPCCALIBTIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "TH2F.h"
#include "TF1.h"
#include "TArrayD.h"
#include "TObjArray.h"

class TH1F;
class TH3F;
class TH2F;
class THnSparse;
class TList;
class AliESDEvent;
class AliESDtrack;

#include "TTreeStream.h"

 
class AliTPCcalibTime:public AliTPCcalibBase {
public:
  AliTPCcalibTime(); 
  AliTPCcalibTime(const Text_t *name, const Text_t *title, ULong64_t TriggerMask, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeDeDx, Int_t deltaIntegrationTimeVdrift);
  virtual ~AliTPCcalibTime();
  
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze();
  //
  void                   ProcessCosmic(AliESDEvent *event);
  Bool_t                 IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1);
  //
  THnSparse *                 GetHistVdrift(){return (THnSparse*) fHistVdrift;};
  THnSparse *                 GetHistDeDxVsTgl(){return (THnSparse*) fHistDeDxTgl;};
  THnSparse *                 GetHistDeDx(){return (THnSparse*) fHistDeDx;};

  

private:

  ULong64_t fTriggerMask;               // select certain trigger within one run

  THnSparse * fHistDeDxTgl;             // dEdx vs. dip angle vs time histogram
  THnSparse * fHistDeDx;                // dEdx vs. time histogram (cosmics: all particles on Fermi plateau)
  THnSparse * fHistVdrift;              // drift velocity vs time histogram

  Float_t fIntegrationTimeDeDx;         // required statistics for each dEdx time bin
  Float_t fIntegrationTimeVdrift;       // required statistics for each Vdrift time bin

  // cuts
  //
  Float_t fCutMaxD;     // maximal distance in rfi ditection
  Float_t fCutTheta;    // maximal distance in theta ditection
  Float_t fCutMinDir;   // direction vector products
 
  AliTPCcalibTime(const AliTPCcalibTime&); 
  AliTPCcalibTime& operator=(const AliTPCcalibTime&); 

  ClassDef(AliTPCcalibTime, 1); 
};

#endif


