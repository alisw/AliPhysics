#ifndef ALITPCCALIBTIMEGAIN_H
#define ALITPCCALIBTIMEGAIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "TH2F.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TArrayD.h"
#include "TObjArray.h"
#include "AliSplineFit.h"

class TH1F;
class TH3F;
class TH2F;
class TList;
class TGraphErrors;
class AliESDEvent;
class AliESDtrack;
class AliTPCcalibLaser;
class AliTPCseed;

#include "TTreeStream.h"


class AliTPCcalibTimeGain:public AliTPCcalibBase {
public:
  AliTPCcalibTimeGain(); 
  AliTPCcalibTimeGain(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeGain);
  virtual ~AliTPCcalibTimeGain();
  //
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           AnalyzeRun(Int_t minEntries);
  //
  void                   ProcessCosmicEvent(AliESDEvent *event);
  void                   ProcessBeamEvent(AliESDEvent *event);
  //
  void                   CalculateBetheAlephParams(TH2F *hist, Double_t * ini);
  static void            BinLogX(THnSparse *h, Int_t axisDim);
  static void            BinLogX(TH1 *h);
  //
  THnSparse *            GetHistGainTime() const {return (THnSparse*) fHistGainTime;};
  TH2F      *            GetHistDeDxTotal() const {return (TH2F*) fHistDeDxTotal;};
  //
  TGraphErrors *         GetGraphGainVsTime(Int_t runNumber = 0, Int_t minEntries = 2000);
  static AliSplineFit *  MakeSplineFit(TGraphErrors * graph);
  TGraphErrors *         GetGraphAttachment(Int_t minEntries, Int_t nmaxBin, Float_t fracLow=0.1, Float_t fracUp=0.9);
  //
  void SetMIP(Float_t MIP){fMIP = MIP;};
  void SetUseMax(Bool_t UseMax){fUseMax = UseMax;};
  void SetLowerTrunc(Float_t LowerTrunc){fLowerTrunc = LowerTrunc;};
  void SetUpperTrunc(Float_t UpperTrunc){fUpperTrunc = UpperTrunc;};
  void SetUseShapeNorm(Bool_t UseShapeNorm){fUseShapeNorm = UseShapeNorm;};
  void SetUsePosNorm(Bool_t UsePosNorm){fUsePosNorm = UsePosNorm;};
  void SetUsePadNorm(Int_t UsePadNorm){fUsePadNorm = UsePadNorm;};
  void SetIsCosmic(Bool_t IsCosmic){fIsCosmic = IsCosmic;};
  void SetLowMemoryConsumption(Bool_t LowMemoryConsumption){fLowMemoryConsumption = LowMemoryConsumption;};
  void SetUseCookAnalytical(Bool_t UseCookAnalytical){fUseCookAnalytical = UseCookAnalytical;};
  //
  void SetCutMinCrossRows(Int_t crossRows){fCutCrossRows = crossRows;};
  void SetCutMaxEta(Float_t maxEta){fCutEtaWindow = maxEta;};
  void SetCutRequireITSrefit(Bool_t requireItsRefit = kFALSE){fCutRequireITSrefit = requireItsRefit;};
  void SetCutMaxDcaXY(Float_t maxXY){fCutMaxDcaXY = maxXY;};
  void SetCutMaxDcaZ(Float_t maxZ){fCutMaxDcaZ = maxZ;};
  //
  void    SetMinTPCsignalN(Float_t minSignalN) { fMinTPCsignalN=minSignalN; }
  Float_t GetMinTPCsignalN() const             { return fMinTPCsignalN;     }
  //
  void SetMinMomentumMIP(Float_t minMom = 0.4){fMinMomentumMIP = minMom;};
  void SetMaxMomentumMIP(Float_t maxMom = 0.6){fMaxMomentumMIP = maxMom;};
  void SetAlephParameters(Float_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j];};

  static void SetMergeEntriesCut(Double_t entriesCut){fgMergeEntriesCut = entriesCut;}

  Double_t GetEntries() const {return fHistGainTime->GetEntries();}

private:
  static Double_t fgMergeEntriesCut;  //maximal number of entries for merging  -can be modified via setter

  Float_t GetTPCdEdx(AliTPCseed * seed);   // wrapper for CookdEdxNorm or analytical
  //
  THnSparse    * fHistGainTime;            // dEdx vs. time, type, Driftlength, momentum P
  TGraphErrors * fGainVsTime;              // multiplication factor vs. time
  TH2F         * fHistDeDxTotal;           // dEdx vs. momentum for quality assurance
  //
  Float_t fIntegrationTimeDeDx;         // required statistics for each dEdx time bin
  //
  Float_t fMIP;                         // rough MIP position in order to have scaleable histograms
  //
  // track cuts
  //
  Int_t   fCutCrossRows;                // minimum number of crossed rows 
  Float_t fCutEtaWindow;                // maximum eta of tracks
  Bool_t  fCutRequireITSrefit;          // if ITSrefit should be required (dangerous in cpass0)
  Float_t fCutMaxDcaXY;                 // max dca_xy (only TPConly resolution is guaranteed!)
  Float_t fCutMaxDcaZ;                  // max dca_z  (dangerous if vDrift is not calibrated)
  Float_t fMinTPCsignalN;               // minimum number of PID clusters
  //
  // definition of MIP window
  //
  Float_t fMinMomentumMIP;              // minimum momentum of MIP region, e.g. 400 MeV
  Float_t fMaxMomentumMIP;              // maximum momentum of MIP region, e.g. 600 MeV
  Float_t fAlephParameters[5];          // parameters for equalization in MIP window, parameter set should be =1 at MIP
  //
  //
  Bool_t  fUseMax;                      // true: use max charge for dE/dx calculation, false: use total charge for dE/dx calculation
  Float_t fLowerTrunc;                  // lower truncation of dE/dx ; at most 5%
  Float_t fUpperTrunc;                  // upper truncation of dE/dx ; ca. 70%
  Bool_t  fUseShapeNorm;                // use empirical correction of dependencies
  Bool_t  fUsePosNorm;                  // charge correction (analytical?)
  Int_t   fUsePadNorm;                  // normalization of pad geometries
  Bool_t  fUseCookAnalytical;           // true if CookdEdxAnalytical should be used
  //
  Bool_t  fIsCosmic;                    // kTRUE if the analyzed runs contain cosmic events
  Bool_t  fLowMemoryConsumption;        // set this option kTRUE if the momenta information should not be stored in order to save memory
  //
  AliTPCcalibTimeGain(const AliTPCcalibTimeGain&); 
  AliTPCcalibTimeGain& operator=(const AliTPCcalibTimeGain&); 
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}

  ClassDef(AliTPCcalibTimeGain, 3);
};

#endif


