#ifndef ALITPCPREPROCESSOROFFLINE_H
#define ALITPCPREPROCESSOROFFLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//
//    Class to create OCDB entries - processing the results of the OFFLINE calibration
//


#include "TNamed.h"
class TObjArray;
class AliTPCcalibTime;
class AliTPCcalibTimeGain;
class AliTPCcalibGainMult;
class AliTPCROC;
class AliTPCParam;
class TPad;
class AliCDBRunRange;
class AliCDBStorage;

class AliTPCPreprocessorOffline:public TNamed { 
public:
  AliTPCPreprocessorOffline();
  virtual ~AliTPCPreprocessorOffline();
  void UpdateOCDBDrift(Int_t ustartRun, Int_t uendRun, AliCDBStorage* storage);
  void UpdateOCDBGain(Int_t  startRunNumber, Int_t endRunNumber, AliCDBStorage* storage);
  void UpdateDriftParam(AliTPCParam *param, TObjArray *const arr, Int_t lstartRun);

  //
  // v drift part
  //
  void GetRunRange(AliTPCcalibTime* const timeDrift);
  void CalibTimeVdrift(const Char_t* file, Int_t ustartRun, Int_t uendRun,AliCDBStorage* ocdbStorage=0x0);
  void AddHistoGraphs(  TObjArray * vdriftArray, AliTPCcalibTime * const timeDrift, Int_t minEntries);
  void AddAlignmentGraphs(  TObjArray * vdriftArray, AliTPCcalibTime * const timeDrift);
  void AddLaserGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift);
  void SetDefaultGraphDrift(TGraph *graph, Int_t color, Int_t style);
  void MakeDefaultPlots(TObjArray * const arr, TObjArray *picArray);
  void SetMaxVDriftCorr(Double_t maxVDriftCorr=0.03) {fMaxVdriftCorr=maxVDriftCorr;};
  Bool_t ValidateTimeDrift();
  //
  // Gain part
  //
  void CalibTimeGain(const Char_t* fileName, Int_t startRunNumber, Int_t endRunNumber,  AliCDBStorage* ocdbStorage);
  void ReadGainGlobal(const Char_t* fileName="CalibObjectsTrain1.root");
  void MakeQAPlot(Float_t  FPtoMIPratio);
  Bool_t AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit = 500, Float_t FPtoMIPratio = 1.43); 
  Bool_t AnalyzeAttachment(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesFit = 2000);
  Bool_t AnalyzePadRegionGain();
  Bool_t AnalyzeGainMultiplicity();
  Bool_t AnalyzeGainChamberByChamber();
  void SetTimeGainRange(Double_t minGain=2.0, Double_t maxGain = 3.0) 
       {fMinGain = minGain; fMaxGain = maxGain;};
  Bool_t ValidateTimeGain();
  //
  // Alignment time part
  //
  void  MakeChainTime();
  void  MakePrimitivesTime();
  //  static void RegisterPrimitiveTimes();
  void  CreateAlignTime(TString fstring, TVectorD paramC);  
  void  MakeFitTime();
  static Double_t EvalAt(Double_t phi, Double_t refX, Double_t theta, Int_t corr, Int_t ptype);
  static Double_t EvalAtPar(Double_t phi, Double_t snp, Double_t refX, Double_t theta, Int_t corr, Int_t ptype, Int_t nstep);

  // event/track counter related setters and getters
  Int_t GetNeventsVdrift() const {return fNeventsVdrift;}
  Int_t GetNtracksVdrift() const {return fNtracksVdrift;}
  Int_t GetMinEventsVdrift() const {return fMinEventsVdrift;}
  Int_t GetMinTracksVdrift() const {return fMinTracksVdrift;}
  void SetMinEventsVdrift(Int_t min) {fMinEventsVdrift=min;}
  void SetMinTracksVdrift(Int_t min) {fMinTracksVdrift=min;}

  //
  // QA drawing part
  //
  static void SetPadStyle(TPad *pad, Float_t mx0, Float_t mx1, Float_t my0, Float_t my1);
  static void PrintArray(TObjArray *array);
  TChain *GetAlignTree(){return fAlignTree;}
  //
  // graph filtering part
  //
  static TGraphErrors* FilterGraphMedianAbs(TGraphErrors * graph, Float_t cut,Double_t &medianY);
  static TGraphErrors* FilterGraphDrift(TGraphErrors * graph, Float_t errSigmaCut, Float_t medianCutAbs);
  static TGraphErrors * MakeGraphFilter0(THnSparse *hisN, Int_t itime, Int_t ival, Int_t minEntries, Double_t offset=0);

  //
  void SwitchOnValidation(Bool_t val = kTRUE) {fSwitchOnValidation = val;} 
  Bool_t IsSwitchOnValidation() { return fSwitchOnValidation; } 

  //
  Int_t GetStatus();
  enum ECalibStatusBit { kCalibFailedTimeDrift =0x0001,
                         kCalibFailedTimeGain  =0x0002 };

private:
  Int_t fMinEntries;                      // minimal number of entries for fit
  Int_t fStartRun;                         // start Run - used to make fast selection in THnSparse
  Int_t fEndRun;                           // end   Run - used to make fast selection in THnSparse
  Int_t fStartTime;                        // fStartTime - used to make fast selection in THnSparse
  Int_t fEndTime;                          // fEndTime   - used to make fast selection in THnSparse
  AliCDBStorage*  fOCDBstorage;            // OCDB storage
  TObjArray * fVdriftArray;               // array with output calibration graphs
  AliTPCcalibTime * fTimeDrift;           // input data to construct calibration graphs
  TGraphErrors * fGraphMIP;                // graph time dependence of MIP
  TGraphErrors * fGraphCosmic;             // graph time dependence at Plateu
  TGraphErrors * fGraphAttachmentMIP;      // graph time dependence of attachment (signal vs. mean driftlength)
  AliSplineFit * fFitMIP;                  // fit of dependence - MIP
  AliSplineFit * fFitCosmic;               // fit of dependence - Plateu
  TObjArray    * fGainArray;               // array to be stored in the OCDB
  AliTPCcalibTimeGain * fGainMIP;          // calibration component for MIP
  AliTPCcalibTimeGain * fGainCosmic;       // calibration component for cosmic
  AliTPCcalibGainMult * fGainMult;         // calibration component for pad region gain equalization and multiplicity dependence

  TChain   *fAlignTree;        //alignment tree
  //
  Bool_t fSwitchOnValidation;  // flag to switch on validation of OCDB parameters
  Float_t fMinGain;   	       // min gain
  Float_t fMaxGain;            // max gain
  Float_t fMaxVdriftCorr;      // max v-drift correction
  Int_t fNtracksVdrift;           // n tracks used for v drift determination
  Int_t fMinTracksVdrift;         // minimum numner of tracks for v drift determination
  Int_t fNeventsVdrift;           // number of events used for drift calibration
  Int_t fMinEventsVdrift;         // minimum number of events used for drift calibration

  Int_t fCalibrationStatus;       // status of calibration, each set bit signifies a failure in a component (see ECalibStatusBit)

private:
  AliTPCPreprocessorOffline& operator=(const AliTPCPreprocessorOffline&); // not implemented
  AliTPCPreprocessorOffline(const AliTPCPreprocessorOffline&); // not implemented
  ClassDef(AliTPCPreprocessorOffline,1)
};

#endif
