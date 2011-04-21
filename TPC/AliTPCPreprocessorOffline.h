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

class AliTPCPreprocessorOffline:public TNamed { 
public:
  AliTPCPreprocessorOffline();
  virtual ~AliTPCPreprocessorOffline();
  void UpdateOCDBDrift(Int_t ustartRun, Int_t uendRun, const char* storagePath);
  void UpdateOCDBGain(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateDriftParam(AliTPCParam *param, TObjArray *const arr, Int_t startRun);

  //
  // v drift part
  //
  void GetRunRange(AliTPCcalibTime* const timeDrift);
  void CalibTimeVdrift(const Char_t* file, Int_t ustartRun, Int_t uendRun,TString ocdbStorage="");
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
  void CalibTimeGain(const Char_t* fileName, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage);
  void ReadGainGlobal(const Char_t* fileName="CalibObjectsTrain1.root");
  void MakeQAPlot(Float_t  FPtoMIPratio);
  Bool_t AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit = 500, Float_t FPtoMIPratio = 1.43); 
  Bool_t AnalyzeAttachment(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesFit = 2000);
  Bool_t AnalyzePadRegionGain();
  Bool_t AnalyzeGainMultiplicity();
  void SetTimeGainRange(Double_t minGain=2.0, Double_t maxGain = 3.0) 
       {fMinGain = minGain; fMaxGain = maxGain;};
  Bool_t ValidateTimeGain();
  //
  // Alignment time part
  //
  void  MakeChainTime();
  void  MakePrimitivesTime();
  void  CreateAlignTime(TString fstring, TVectorD paramC);  
  void  MakeFitTime();
  static Double_t EvalAt(Double_t phi, Double_t refX, Double_t theta, Int_t corr, Int_t ptype);

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


private:
  Int_t fMinEntries;                      // minimal number of entries for fit
  Int_t startRun;                         // start Run - used to make fast selection in THnSparse
  Int_t endRun;                           // end   Run - used to make fast selection in THnSparse
  Int_t startTime;                        // startTime - used to make fast selection in THnSparse
  Int_t endTime;                          // endTime   - used to make fast selection in THnSparse
  TString  ocdbStorage;                   // path to the OCDB storage
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

private:
  AliTPCPreprocessorOffline& operator=(const AliTPCPreprocessorOffline&); // not implemented
  AliTPCPreprocessorOffline(const AliTPCPreprocessorOffline&); // not implemented
  ClassDef(AliTPCPreprocessorOffline,1)
};

#endif
