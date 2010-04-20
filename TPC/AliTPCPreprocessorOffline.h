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
class AliTPCROC;
class AliTPCParam;
class TPad;
class AliCDBRunRange;

class AliTPCPreprocessorOffline:public TNamed { 
public:
  AliTPCPreprocessorOffline();
  virtual ~AliTPCPreprocessorOffline();
  void MakeAll(){;}
  void UpdateOCDBDrift(Int_t ustartRun, Int_t uendRun, const char* storagePath);
  void UpdateOCDBGain(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);
  void UpdateDriftParam(AliTPCParam *param, TObjArray *arr, Int_t startRun);

  //
  // v drift part
  //
  void GetRunRange(AliTPCcalibTime* timeDrift);
  void CalibTimeVdrift(Char_t* file, Int_t ustartRun, Int_t uendRun,TString ocdbStorage="");
  void AddHistoGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift, Int_t minEntries);
  void AddAlignmentGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift);
  void AddLaserGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift);
  void SetDefaultGraphDrift(TGraph *graph, Int_t color, Int_t style);
  void MakeDefaultPlots(TObjArray * arr, TObjArray *picArray);
  //
  // Gain part
  //
  void CalibTimeGain(Char_t* fileName, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage);
  void ReadGainGlobal(Char_t* fileName="CalibObjectsTrain1.root");
  void MakeQAPlot(Float_t  FPtoMIPratio);
  Bool_t AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit = 500, Float_t FPtoMIPratio = 1.43); 
  //
  // QA drawing part
  //
  static void SetPadStyle(TPad *pad, Float_t mx0, Float_t mx1, Float_t my0, Float_t my1);
  static void PrintArray(TObjArray *array);
  //
  // graph filtering part
  //
  static TGraphErrors* FilterGraphMedianAbs(TGraphErrors * graph, Float_t cut,Double_t &medianY);
  static TGraphErrors* FilterGraphDrift(TGraphErrors * graph, Float_t errSigmaCut, Float_t medianCutAbs);
  static TGraphErrors * MakeGraphFilter0(THnSparse *hisN, Int_t itime, Int_t ival, Int_t minEntries, Double_t offset=0);

public:
  Int_t kMinEntries;                      // minimal number of entries for fit
  Int_t startRun;                         // start Run - used to make fast selection in THnSparse
  Int_t endRun;                           // end   Run - used to make fast selection in THnSparse
  Int_t startTime;                        // startTime - used to make fast selection in THnSparse
  Int_t endTime;                          // endTime   - used to make fast selection in THnSparse
  TString  ocdbStorage;                   // path to the OCDB storage
  TObjArray * fVdriftArray;
  AliTPCcalibTime * fTimeDrift;
  TGraphErrors * fGraphMIP;                // graph time dependence of MIP
  TGraphErrors * fGraphCosmic;             // graph time dependence at Plateu
  AliSplineFit * fFitMIP;                  // fit of dependence - MIP
  AliSplineFit * fFitCosmic;               // fit of dependence - Plateu
  TObjArray    * fGainArray;               // array to be stored in the OCDB
  AliTPCcalibTimeGain * fGainMIP;          // calibration component for MIP
  AliTPCcalibTimeGain * fGainCosmic;       // calibration component for cosmic

private:
  AliTPCPreprocessorOffline& operator=(const AliTPCPreprocessorOffline&); // not implemented
  AliTPCPreprocessorOffline(const AliTPCPreprocessorOffline&); // not implemented
  ClassDef(AliTPCPreprocessorOffline,1)
};

#endif
