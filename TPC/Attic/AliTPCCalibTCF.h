#ifndef ALI_TPC_CALIB_TCF_H
#define ALI_TPC_CALIB_TCF_H


/* Copyright(c) 2007-08, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                             */

///////////////////////////////////////////////////////////////////////////////
//                        Class AliTPCCalibTCF                               //
// Class for Extraction and test of TCF parameters needed by the ALTRO chip  //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include <exception>


class TObject;
class TGraph;
class AliRawReader;
class AliTPCRawStreamV3;
class TNtuple;
class TH1F;
class TH2F;

class AliTPCCalibTCF : public TNamed {

public:

  AliTPCCalibTCF();
  AliTPCCalibTCF(Int_t gateWidth, Int_t Sample, Int_t pulseLength, Int_t lowPulseLim, Int_t upPulseLim, Double_t rmsLim, Double_t ratioIntLim);
  AliTPCCalibTCF(const AliTPCCalibTCF &sig);
  virtual ~AliTPCCalibTCF();
  
  AliTPCCalibTCF& operator = (const  AliTPCCalibTCF &source);

 
  void ProcessRawFileV3(const char *nameRawFile, const char *nameFileOut);
  void ProcessRawEventV3(AliRawReader *rawReader,AliTPCRawStreamV3 *rawStream, const char *nameFileOut);

  void MergeHistoPerSector(const char *nameFileIn);

  void AnalyzeRootFile(const char *nameFileIn, Int_t minNumPulse=1, Int_t histStart=1, Int_t histEnd=1000000);
  Int_t AnalyzePulse(TH1F * const hisIn, Double_t *coefZ, Double_t *coefP); 

  void TestTCFonRootFile(const char *nameFileIn, const char *nameFileTCF, Int_t nPulseMin=0, Int_t plotFlag=0, Int_t lowKey=1, Int_t upKey=1000000);
  void TestTCFonRawFile(const char *nameRawFile, const char *nameFileOut, const char *nameFileTCF, Int_t nPulseMin=0, Int_t plotFlag=0, bool bUseHLTOUT=false);

  Int_t DumpTCFparamToFilePerSector(const char *nameFileTCFPerSec, const char *nameMappingFile="$ALICE_ROOT/TPC/Calib/tpcMapping.root");  
  Int_t DumpTCFparamToFilePerPad(const char *nameFileTCFPerPad,const char *nameFileTCFPerSec, const char *nameMappingFile="$ALICE_ROOT/TPC/Calib/tpcMapping.root");  
 
  TH2F *PlotOccupSummary2Dhist(const char *nameFileIn, Int_t side=0);
  void PlotOccupSummary(const char *nameFile, Int_t side=0, Int_t nPulseMin=0); 

  void PlotQualitySummary(const char *nameFileQuality, const char *plotSpec="widthRed:maxUndershot", const char *cut="maxUndershot<0.1&&maxUndershot>-40&&widthRed>0&&widthRed<100", const char *pOpt="LEGO2Z");

  void PrintPulseThresholds();

  void MergeHistoPerFile(const char *fileNameIn, const char *fileSum, Int_t mode=0);
  void MergeToOneFile(const char *nameFileSum);

private:
  
  // tresholds for proper pulse finder (Analyze functions)
  Int_t fGateWidth;     // expected Gate fluctuation length
  Int_t fSample;        // expected usefull signal length
  Int_t fPulseLength;   // needed pulselength for TC characterisation
  Int_t fLowPulseLim;   // lower pulse height limit
  Int_t fUpPulseLim;    // upper pulse height limit
  Double_t fRMSLim;     // signal RMS limit
  Double_t fRatioIntLim;// ratio of signal-integral/pulse-integral limit

  Int_t FitPulse(TNtuple *dataTuple, Double_t *coefZ, Double_t *coefP);
  static void FitFcn(Int_t &nPar, Double_t *grad, Double_t &f, Double_t * const par, Int_t iflag);

  Double_t* ExtractPZValues(Double_t *param);
  Int_t Equalization(TNtuple *dataTuple, Double_t *coefZ, Double_t *coefP);

  Int_t FindCorTCFparam(TH1F * const hisIn, const char *nameFileTCF, Double_t *coefZ, Double_t *coefP);
  Double_t *GetQualityOfTCF(TH1F *hisIn, Double_t *coefZ, Double_t *coefP,Int_t plotFlag=0); 

  TNtuple *ApplyTCFilter(TH1F * const hisIn, Double_t * const coefZ, Double_t * const coefP, Int_t plotFlag=0);

  ClassDef(AliTPCCalibTCF,1);

};
#endif
