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
class AliTPCcalibLaser;

#include "TTreeStream.h"
#include "TMap.h"
 
class AliTPCcalibTime:public AliTPCcalibBase {
public:
  AliTPCcalibTime(); 
  AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeDeDx, Int_t deltaIntegrationTimeVdrift);
  virtual ~AliTPCcalibTime();
  
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze();
  //
  static Bool_t           IsLaser(AliESDEvent *event);
  void                   ProcessLaser (AliESDEvent *event);
  void                   ProcessCosmic(AliESDEvent *event);
  Bool_t                 IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1);
  //
  THnSparse *                 GetHistVdrift(){return (THnSparse*) fHistVdrift;};
  THnSparse *                 GetHistDeDxVsTgl(){return (THnSparse*) fHistDeDxTgl;};
  THnSparse *                 GetHistDeDx(){return (THnSparse*) fHistDeDx;};
  THnSparse *                 GetHistVdriftLaserA(Int_t index=1){return (THnSparse*) fHistVdriftLaserA[index];};
  THnSparse *                 GetHistVdriftLaserC(Int_t index=1){return (THnSparse*) fHistVdriftLaserC[index];};
  TMap      *                 GetMapDz(){return (TMap *) fMapDz;};

  
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
private:
  void ResetCurrent();                  // reset current values

  THnSparse * fHistDeDxTgl;             // dEdx vs. dip angle vs time histogram
  THnSparse * fHistDeDx;                // dEdx vs. time histogram (cosmics: all particles on Fermi plateau)
  THnSparse * fHistVdrift;              // drift velocity vs time histogram

  Float_t fIntegrationTimeDeDx;         // required statistics for each dEdx time bin
  Float_t fIntegrationTimeVdrift;       // required statistics for each Vdrift time bin

  AliTPCcalibLaser * fLaser;            //! laser calibration
  //
  // current information
  //
  Float_t fDz;          //! current delta z
  Float_t fdEdx;        //! current dEdx
  Float_t fdEdxRatio;   //! current dEdx ratio
  Float_t fTl;          //! current tan(lambda)
  
  // cuts
  //
  Float_t fCutMaxD;     // maximal distance in rfi ditection
  Float_t fCutMaxDz;     // maximal distance in z ditection
  Float_t fCutTheta;    // maximal distance in theta ditection
  Float_t fCutMinDir;   // direction vector products
 
  AliTPCcalibTime(const AliTPCcalibTime&); 
  AliTPCcalibTime& operator=(const AliTPCcalibTime&); 


  // laser histo
  THnSparse * fHistVdriftLaserA[3];		//NEW! Histograms for V drift from laser
  THnSparse * fHistVdriftLaserC[3];		//NEW! Histograms for V drift from laser
  // DELTA Z histo
  TMap      * fMapDz;			//NEW! Tmap of V drifts for different triggers

  Int_t    fTimeBins;			//NEW! Bins time
  Double_t fTimeStart;			//NEW! Start time
  Double_t fTimeEnd;			//NEW! End time
  Int_t    fPtBins;			//NEW! Bins pt
  Double_t fPtStart;			//NEW! Start pt
  Double_t fPtEnd;			//NEW! End pt
  Int_t    fVdriftBins;			//NEW! Bins vdrift
  Double_t fVdriftStart;		//NEW! Start vdrift
  Double_t fVdriftEnd;			//NEW! End vdrift
  Int_t    fRunBins;			//NEW! Bins run
  Double_t fRunStart;			//NEW! Start run
  Double_t fRunEnd;			//NEW! End run
  Int_t    fBinsVdrift[4];		//NEW! Bins for vdrift
  Double_t fXminVdrift[4];		//NEW! Xmax for vdrift
  Double_t fXmaxVdrift[4];		//NEW! Xmin for vdrift
  ClassDef(AliTPCcalibTime, 1); 
};

#endif


