#ifndef ALITPCCALIBTIMEGAIN_H
#define ALITPCCALIBTIMEGAIN_H

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
class TGraphErrors;
class AliESDEvent;
class AliESDtrack;
class AliTPCcalibLaser;

#include "TTreeStream.h"


class AliTPCcalibTimeGain:public AliTPCcalibBase {
public:
  AliTPCcalibTimeGain(); 
  AliTPCcalibTimeGain(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeGain);
  virtual ~AliTPCcalibTimeGain();
  //
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze();
  //
  void                   ProcessCosmicEvent(AliESDEvent *event);
  void                   ProcessBeamEvent(AliESDEvent *event);
  //
  void                   CalculateBetheAlephParams(TH2F *hist, Double_t * ini);
  static void            BinLogX(THnSparse *h, Int_t axisDim);
  static void            BinLogX(TH1 *h);
  //
  THnSparse *            GetHistGainTime(){return (THnSparse*) fHistGainTime;};
  TH2F      *            GetHistDeDxTotal(){return (TH2F*) fHistDeDxTotal;};
  //
  TGraphErrors *         GetGraphGainVsTime(){return fGainVsTime;};
  //
  void SetMIP(Float_t MIP){fMIP = MIP;};  
  void SetLowerTrunc(Float_t LowerTrunc){fLowerTrunc = LowerTrunc;};
  void SetUpperTrunc(Float_t UpperTrunc){fUpperTrunc = UpperTrunc;};
  void SetUseShapeNorm(Bool_t UseShapeNorm){fUseShapeNorm = UseShapeNorm;};
  void SetUsePosNorm(Bool_t UsePosNorm){fUsePosNorm = UsePosNorm;};
  void SetUsePadNorm(Int_t UsePadNorm){fUsePadNorm = UsePadNorm;};
  void SetIsCosmic(Bool_t IsCosmic){fIsCosmic = IsCosmic;};
  void SetLowMemoryConsumption(Bool_t LowMemoryConsumption){fLowMemoryConsumption = LowMemoryConsumption;};

private:
  //
  THnSparse    * fHistGainTime;            // dEdx vs. time, type, Driftlength, momentum P
  TGraphErrors * fGainVsTime;              // multiplication factor vs. time
  TH2F         * fHistDeDxTotal;           // dEdx vs. momentum for quality assurance
  //
  Float_t fIntegrationTimeDeDx;         // required statistics for each dEdx time bin
  //
  Float_t fMIP;                         // rough MIP position in order to have scaleable histograms
  //
  Float_t fLowerTrunc;                  // lower truncation of dE/dx ; at most 5%
  Float_t fUpperTrunc;                  // upper truncation of dE/dx ; ca. 70%
  Bool_t  fUseShapeNorm;                // use empirical correction of dependencies
  Bool_t  fUsePosNorm;                  // charge correction (analytical?)
  Int_t   fUsePadNorm;                  // normalization of pad geometries
  //
  Bool_t  fIsCosmic;                    // kTRUE if the analyzed runs contain cosmic events
  Bool_t  fLowMemoryConsumption;        // set this option kTRUE if the momenta information should not be stored in order to save memory
  //
  AliTPCcalibTimeGain(const AliTPCcalibTimeGain&); 
  AliTPCcalibTimeGain& operator=(const AliTPCcalibTimeGain&); 
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}

  ClassDef(AliTPCcalibTimeGain, 1); 
};

#endif


