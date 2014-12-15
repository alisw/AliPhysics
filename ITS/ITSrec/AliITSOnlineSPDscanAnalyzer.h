#ifndef ALIITSONLINESPDSCANANALYZER_H
#define ALIITSONLINESPDSCANANALYZER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// This class is used in the detector algorithm framework //
// to process the data stored in special container files  //
// (see AliITSOnlineSPDscan). For instance, minimum       //
// threshold values can be extracted.                     //
////////////////////////////////////////////////////////////

#include <TString.h>
#include <TH1F.h>

class AliITSOnlineSPDscan;
class AliITSOnlineCalibrationSPDhandler;
class TGraph;
class TH2F;
class TArrayI;

class AliITSOnlineSPDscanAnalyzer {

 public:
  AliITSOnlineSPDscanAnalyzer(const Char_t *fileName, AliITSOnlineCalibrationSPDhandler *handler, Bool_t readFromGridFile=kFALSE);
  AliITSOnlineSPDscanAnalyzer(const AliITSOnlineSPDscanAnalyzer& handle);
  ~AliITSOnlineSPDscanAnalyzer();

  AliITSOnlineSPDscanAnalyzer& operator=(const AliITSOnlineSPDscanAnalyzer& handle);

  Bool_t     IsChipPresent(UInt_t hs, UInt_t chipNr);
  Bool_t     IsOverWriteSet() const {return fOverWrite;}
  void       SetCalibHandler(AliITSOnlineCalibrationSPDhandler * const handler) {fHandler=handler;}
  void       SetParam(const Char_t *pname, const Char_t *pval);
  void       ReadParamsFromLocation(const Char_t *dirName);

  UInt_t     GetType() const {return fType;}
  UInt_t     GetDacId() const {return fDacId;}
  
  Int_t      GetDelay(UInt_t hs, UInt_t chipNr);
  Int_t      GetMinTh(UInt_t hs, UInt_t chipNr);
  TArrayI    GetMeanTh(UInt_t hs, UInt_t chipNr);
  
  Int_t      GetNrNoisyUnima(UInt_t hs, UInt_t chipNr);

  Bool_t     ProcessUniformity();
  Bool_t     ProcessDeadPixels();
  Bool_t     ProcessNoisyPixels();

  Bool_t     ProcessNrTriggers();

  AliITSOnlineSPDscan* GetOnlineScan() const {return fScanObj;}
  UInt_t     GetRouterNr();
  Bool_t     GetHalfStavePresent(UInt_t hs);

  TGraph*    GetNrTriggersG();
  TGraph*    GetMeanMultiplicityG(UInt_t hs, UInt_t chipNr);
  TGraph*    GetHitEventEfficiencyG(UInt_t hs, UInt_t chipNr);
  TH2F*      GetHitMapTot(UInt_t step);
  TH2F*      GetPhysicalHitMapTot(UInt_t step);
  TH2F*      GetHitMapChip(UInt_t step, UInt_t hs, UInt_t chip);

  Float_t    GetTPeff() const {return fTPeff;}
  TH1F*      GetTPeffHS() const {return fTPeffHS;}
  TH1F*      GetTPeffChip(UInt_t hs) const {return fTPeffChip[hs];}
  Float_t    GetDeadPixel() const {return fDeadPixel;}
  TH1F*      GetDeadPixelHS() const {return fDeadPixelHS;}
  TH1F*      GetDeadPixelChip(UInt_t hs) const {return fDeadPixelChip[hs];}
  Float_t    GetNoisyPixel() const {return fNoisyPixel;}
  TH1F*      GetNoisyPixelHS() const {return fNoisyPixelHS;}
  TH1F*      GetNoisyPixelChip(UInt_t hs) const {return fNoisyPixelChip[hs];}

 private:
  UInt_t               fType;           // calib type
  UInt_t               fDacId;          // dac id
  TString              fFileName;       // container file name
  enum                 calibvals{kMINTH,kMEANTH,kDAC,kUNIMA,kNOISE,kDELAY};  // calib types

  AliITSOnlineSPDscan               *fScanObj;  // container obj
  AliITSOnlineCalibrationSPDhandler *fHandler;  // calib helper obj
  Bool_t     fbModuleScanned[240];        // is module used in scan?

  TGraph*    fMeanMultiplicity[6][11];   // mean mult graphs
  TGraph*    fHitEventEfficiency[6][11]; // hit event graphs
  TGraph*    fTriggers;                  // trigger graph

  // uniformity scan analysis:
  Float_t    fTPeff;                     // number of good pixels [%] (for full router)
  TH1F*      fTPeffHS;                   // 6 bin histogram, number good pixels [%] (for each hs)
  TH1F*      fTPeffChip[6];              // 10 bin histograms, number good pixels [%] (for each chip)
  Float_t    fDeadPixel;                 // number of dead pixels [%] (for full router)
  TH1F*      fDeadPixelHS;               // 6 bin histogram, number dead pixels [%] (for each hs)
  TH1F*      fDeadPixelChip[6];          // 10 bin histograms, number dead pixels [%] (for each chip)
  Float_t    fNoisyPixel;                // number of 'noisy' pixels [%] (for full router)
  TH1F*      fNoisyPixelHS;              // 6 bin histogram, number 'noisy' pixels [%] (for each hs)
  TH1F*      fNoisyPixelChip[6];         // 10 bin histograms, number 'noisy' pixels [%] (for each chip)
  
  void       Init(Bool_t readFromGridFile=kFALSE);                     // init

  void       CreateUniformityHistograms(); // method to create all histograms to be filled by 'ProcessUniformity'
  void       DeleteUniformityHistograms(); // method to delete all histograms used by uniformity scan analysis

  Bool_t     ProcessMeanMultiplicity();  // process mean mult
  Bool_t     ProcessHitEventEfficiency();// process hit event eff

  Int_t      FindLastMinThDac(UInt_t hs, UInt_t chipNr);  // dac value where fMinIncreaseFromBaseLine reached
  Int_t      FindClosestLowerStep(Float_t dacValueInput); // step closest (lower) to a dacvalue 
  Float_t    GetCompareLine(UInt_t step, UInt_t hs, UInt_t chipNr, Float_t basePar2); // line to compare mean mult with

  // dead noisy parameters:
  Bool_t     fOverWrite;          // overWrite old dead/noisy or just add new ones to it
  // noise scan parameters:
  Float_t    fNoiseThreshold;     // ratio of allowed hits/triggers
  UInt_t     fNoiseMinimumEvents; // minimum events required to find noisy pixels
  // min th scan parameters:
  UInt_t     fMinNrStepsBeforeIncrease; // min nr of steps required before fMinIncreaseFromBaseLine reached
  Float_t    fMinIncreaseFromBaseLine;  // min increase of mean mult from base line
  UInt_t     fStepDownDacSafe;          // nr of steps down to put threshold result (to be on the safe side)
  Float_t    fMaxBaseLineLevel;         // maximum value for the base line to compare with



};

#endif


