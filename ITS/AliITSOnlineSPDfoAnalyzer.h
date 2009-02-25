#ifndef ALIITSONLINESPDFOANALYZER_H
#define ALIITSONLINESPDFOANALYZER_H
/* Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                     // 
// This class is used within the detector algorithm framework //
// to analyze FO scan data. It intends to find the best DAC   //
// values to get the best FO trigger efficiency               //
////////////////////////////////////////////////////////////////

class THnSparse;
class TObject;
class TArrayI;
class AliITSOnlineSPDfo;
class AliITSOnlineSPDfoInfo;

class AliITSOnlineSPDfoAnalyzer {
  
 public:
  AliITSOnlineSPDfoAnalyzer(const TString fileName,  Bool_t readFromGridFile=kFALSE);
  AliITSOnlineSPDfoAnalyzer(const AliITSOnlineSPDfoAnalyzer& foan);
  ~AliITSOnlineSPDfoAnalyzer();
  
  AliITSOnlineSPDfoAnalyzer& operator=(const AliITSOnlineSPDfoAnalyzer& handle);
  
  void       Init(Bool_t readFromGridFile=kFALSE);     
  
  enum {kNqualityFlags=3};
  
  void ReadParamsFromLocation(const Char_t* dirName);
  Int_t IsSelected(Float_t eff) const;                  // selection quality (0 = best, 1 tight, 2 loose)
  Int_t Select(const AliITSOnlineSPDfoChip *chip) const;
  void WriteToFile(TString outputfile);
  Bool_t IsExisting(TArrayI dacs,Int_t hs, Int_t chip) const;
  Bool_t CorrectPreVTHChioce(const TH1D *h,Int_t &bin) const;
  
  // SETTERS
  void SetGeneralThresholds(Float_t thre[3]); 
  void SetParam(const Char_t *pname, const Char_t *pval);
  void SetNdimensions();
  void BuildTHnSparse(Int_t ihs, Int_t ichip);
  void Process();
  void CheckResults(TString filename, Int_t hs, Int_t ichip, Int_t iqual) const;
  
  // GETTERS
  void GetCanvases(const THnSparse *hn, Int_t ihs, Int_t ichip,Int_t iqual) const;
  TArrayI ChooseDACValues(Int_t ihs, Int_t ichip) const;
  TArrayI GetCentralDACS(Int_t qualityflag, Int_t hs, Int_t chip, const TH1D *hd) const;
  AliITSOnlineSPDfo * GetFOHandler() const {return fFOHandler;}
  
 private:
  TString fFileName;
  Int_t fNdims; // number of dimensions of the histogram (= #DACs in the scan)
  Int_t *fNbins; //[fNdims]  
  Double_t *fXmin; //[fNdims]
  Double_t *fXmax; //[fNdims]  
  AliITSOnlineSPDfo *fFOHandler;
  Float_t fGeneralThresholds[3];
  THnSparse *fNh[3][6][10];      // N-dim histo per chip per half sector per quality flag [0= exact, 1 = within 0.01, 2 = within 0.05]
  Bool_t fHighOccupancyCheck;
};

#endif
