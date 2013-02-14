/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// Class for spectrum correction
// Subtraction of hadronic background, Unfolding of the data and
// Renormalization done here
// For more information see the implementation file
//
#ifndef ALIHFEBEAUTYSPECTRUM_H
#define ALIHFEBEAUTYSPECTRUM_H

#include "AliHFECorrectSpectrumBase.h"



class TGraphErrors;
class TObject;
class TH1;
class TF1;
class TList;
class TObjArray;
class AliCFContainer;
class AliHFEcontainer;
class AliCFDataGrid;
class AliCFEffGrid;
class AliHFEInclusiveSpectrumQA;
class AliHFEBeautySpectrumQA;

class AliHFEBeautySpectrum : public AliHFECorrectSpectrumBase{
 public:
  
  enum{
    kElecBgSources = 9,
    kBgLevels = 3,
    kBgPtBins = 44,
    kSignalPtBins = 35,
    kCentrality = 12
  };
  
  AliHFEBeautySpectrum(const char* name);
  ~AliHFEBeautySpectrum();
  
  
  virtual Bool_t Init(const AliHFEcontainer *datahfecontainer, const AliHFEcontainer *mchfecontainer, const AliHFEcontainer *bghfecontainer=0x0, const AliHFEcontainer */*v0hfecontainer*/ = 0x0,AliCFContainer */*photoniccontainerD*/ = 0x0);
  virtual Bool_t Correct(Bool_t subtractcontamination=kTRUE, Bool_t /*subtractphotonic*/=kFALSE);
  
  AliCFDataGrid *SubtractBackground(Bool_t setBackground = kFALSE);
  
  AliCFDataGrid *CorrectParametrizedEfficiency(AliCFDataGrid* const bgsubpectrum = 0x0);
  
  THnSparse *Unfold(AliCFDataGrid* const bgsubpectrum = 0x0);
  void UnfoldBG(AliCFDataGrid* const bgsubpectrum);
  AliCFDataGrid *CorrectForEfficiency(AliCFDataGrid* const bgsubpectrum = 0x0);
  
  void SetPbPbAnalysis(Bool_t isPbPb = kFALSE) { fBeamType=(Char_t) isPbPb; };
  void SetEtaSyst(Bool_t etaSyst = kTRUE) { fEtaSyst = etaSyst; };
  
  void SetParameterizedEff(AliCFContainer *container, AliCFContainer *containermb, AliCFContainer *containeresd, AliCFContainer *containeresdmb, Int_t *dimensions);
  
  void SetNumberOfMCEvents(Int_t nEvents) { fNMCEvents = nEvents; };
  void SetNumberOfMC2Events(Int_t nEvents) { fNMCbgEvents = nEvents; };
  void SetUnSetCorrelatedErrors(Bool_t unsetcorrelatederrors) {fUnSetCorrelatedErrors = unsetcorrelatederrors;};
  void SetFillMoreCorrelationMatrix(Bool_t fillMoreCorrelationMatrix) { fFillMoreCorrelationMatrix = fillMoreCorrelationMatrix;}
  
  void SetNCentralityBinAtTheEnd(Int_t nCentralityBinAtTheEnd) {fNCentralityBinAtTheEnd = nCentralityBinAtTheEnd; };
  void SetLowHighBoundaryCentralityBinAtTheEnd(Int_t low, Int_t high, Int_t i) { fLowBoundaryCentralityBinAtTheEnd[i] = low; fHighBoundaryCentralityBinAtTheEnd[i] = high;};
  
  void CallInputFileForBeauty2ndMethod();
  void SetInputFileForBeauty2ndMethod(const char *filenameb = "BSpectrum2ndmethod.root"){fkBeauty2ndMethodfilename = filenameb; };
  void SetBeautyAnalysis2ndMethod(Bool_t beauty2ndmethod) { fBeauty2ndMethod = beauty2ndmethod; }
  void SetIPEffCombinedSamples(Bool_t ipEffCombinedSamples) { fIPEffCombinedSamples = ipEffCombinedSamples; }
  void SetHadronEffbyIPcut(THnSparseF* hsHadronEffbyIPcut) { fHadronEffbyIPcut = hsHadronEffbyIPcut;};
  void SetNonHFEsyst(Bool_t syst){ fNonHFEsyst = syst; };
  
  void SetDumpToFile(Bool_t dumpToFile) { fDumpToFile=dumpToFile; }; 
  
  void SetDebugLevel(Int_t debugLevel, Bool_t writeToFile = kFALSE) { fDebugLevel = debugLevel; fWriteToFile = writeToFile; };
  void SetUnfoldBG() { fUnfoldBG = kTRUE; };
  
  
  AliCFDataGrid* GetRawBspectra2ndMethod();
  AliCFDataGrid* GetCharmBackground();
  AliCFDataGrid* GetNonHFEBackground(Int_t decay = 0, Int_t source = 0);
  THnSparse* GetCharmWeights(Int_t centBin);
  THnSparse* GetBeautyIPEff(Bool_t isMCpt);
  THnSparse* GetPIDxIPEff(Int_t source);
  void CalculateNonHFEsyst();
  
  void EnableIPanaHadronBgSubtract() { fIPanaHadronBgSubtract = kTRUE; };
  void EnableIPanaCharmBgSubtract() { fIPanaCharmBgSubtract = kTRUE; };
  void EnableIPanaNonHFEBgSubtract() { fIPanaNonHFEBgSubtract = kTRUE; };
  void EnableIPParameterizedEff() { fIPParameterizedEff = kTRUE; };
  
 protected:
  void AddTemporaryObject(TObject *cont);
  void ClearObject(TObject *o);
  
 private:
  AliHFEBeautySpectrum(const AliHFEBeautySpectrum &ref);
  AliHFEBeautySpectrum &operator=(const AliHFEBeautySpectrum &ref);
  virtual void Copy(TObject &o) const;
  
  AliHFEBeautySpectrumQA *fQA; // Beauty QA
  
  //shift some of the following to QA or Base class:
  TList *fTemporaryObjects;     // Emulate garbage collection
  AliCFDataGrid *fBackground;   // Background Grid
  TF1 *fEfficiencyTOFPIDD;       // TOF PID efficiency parameterized
  TF1 *fEfficiencyesdTOFPIDD;    // TOF PID efficiency parameterized
  TF1 *fEfficiencyIPCharmD;      // IP efficiency parameterized for charm
  TF1 *fEfficiencyIPBeautyD;     // IP efficiency parameterized for beauty 
  TF1 *fEfficiencyIPBeautyesdD;  // IP efficiency parameterized for beauty for esd
  TF1 *fEfficiencyIPConversionD; // IP efficiency parameterized for conversion
  TF1 *fEfficiencyIPNonhfeD;     // IP efficiency parameterized for nonhfe 
  TF1 *fNonHFEbg;     // Efficiency Function
  THnSparseF *fWeightCharm;     // Weight for charm bg
  
  AliCFContainer *fConvSourceContainer[kElecBgSources][kBgLevels]; //container for conversion electrons, divided into different photon sources
  AliCFContainer *fNonHFESourceContainer[kElecBgSources][kBgLevels]; //container for non-HF electrons, divided into different sources
  
  Bool_t fInclusiveSpectrum;     // Inclusive Spectrum
  Bool_t fDumpToFile;           // Write Result in a file
  
  Bool_t fUnSetCorrelatedErrors;    // Unset correlated errors
  
  Bool_t fIPanaHadronBgSubtract;     // Hadron background subtraction
  Bool_t fIPanaCharmBgSubtract;      // Charm background subtraction 
  Bool_t fIPanaNonHFEBgSubtract;     // nonHFE background subtraction
  Bool_t fIPParameterizedEff;        // switch to use parameterized efficiency for ip analysis
  Bool_t fNonHFEsyst;            // choose NonHFE background level (upper, lower, central)
  Bool_t fBeauty2ndMethod;      // 2nd method to get beauty spectrum
  Bool_t fIPEffCombinedSamples; // flag to combine two different samples
  
  Int_t fNMCEvents;             // Number of MC Events
  Int_t fNMCbgEvents;       // Number of BG MC Events
  
  Int_t fNCentralityBinAtTheEnd;// Number of centrality class at the end
  Int_t fLowBoundaryCentralityBinAtTheEnd[20];  // Boundary of the bins
  Int_t fHighBoundaryCentralityBinAtTheEnd[20];  // Boundary of the bins
  Bool_t fFillMoreCorrelationMatrix;             // For low stats to have reasonable errors
  
  THnSparseF *fHadronEffbyIPcut;// container for hadron efficiency by IP cut
  TH1D *fEfficiencyCharmSigD; // charm IP cut eff from signal enhanced MC
  TH1D *fEfficiencyBeautySigD; // beauty IP cut eff from signal enhanced MC
  TH1D *fEfficiencyBeautySigesdD; // beauty IP cut eff from signal enhanced MC for esd
  TH1D *fConversionEff;     // conversion IP cut eff
  TH1D *fNonHFEEff;         // nonhfe IP cut eff
  TH1D *fCharmEff;          // charm IP cut eff
  TH1D *fBeautyEff;         // beauty IP cut eff
  TH1D *fConversionEffbgc;      // conversion IP cut eff
  TH1D *fNonHFEEffbgc;          // nonhfe IP cut eff
  TH1D *fBSpectrum2ndMethod;             // beauty spectrum for 2nd method
  const char *fkBeauty2ndMethodfilename;      // name of file, which contains beauty spectrum for 2ndmethod
  Char_t fBeamType;             // beamtype; default -1; pp =0; PbPb=1
  Bool_t fEtaSyst;              // pp 2.76 TeV (= kTRUE) or 7 TeV (= kFALSE)
  
  
  Int_t fDebugLevel;            // Debug Level
  Bool_t fWriteToFile;           // Write plots to eps files
  Bool_t fUnfoldBG;             // flag to unfold backgroud
  
  ClassDef(AliHFEBeautySpectrum, 1) 
    };
#endif

