#ifndef ALITPCCALIBKR_H
#define ALITPCCALIBKR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliTPCclusterKr.h"

class TH3F;
class TH1D;

class AliTPCCalibKr : public TObject {

public:
  AliTPCCalibKr();
  AliTPCCalibKr(const AliTPCCalibKr&); // copy constructor
  virtual ~AliTPCCalibKr();

  AliTPCCalibKr& operator=(const AliTPCCalibKr&); 

  //
  void Init();
  Bool_t Process(AliTPCclusterKr *cluster);
  Bool_t Accept(AliTPCclusterKr *cluster);
  Bool_t Update(AliTPCclusterKr *cluster);
  TH3F*  CreateHisto(Int_t chamber);

  const TObjArray* GetHistoKrArray () {return &fHistoKrArray;}  // get calibration object
  TH3F* GetHistoKr(Int_t sector) const;                         // get refernce histogram

  Bool_t IsCSide(Int_t chamber);
  Bool_t IsIROC(Int_t chamber);

  void Analyse();
  static TH1D* ProjectHisto(TH3F* histo3D, const char* name = "_pz", Int_t firstxbin = 0, Int_t lastxbin = 0, Int_t firstybin = 0, Int_t lastybin = 0);

  void SetASide(Bool_t bA = kTRUE) {fASide = bA;} // fill histo only A TPC side
  void SetBSide(Bool_t bC = kTRUE) {fCSide = bC;} // fill histo only C TPC side

  //Merge output objects (needed by PROOF)
  virtual Long64_t Merge(TCollection* list);

  // Setters
  void SetADCOverClustSizeRange(Float_t min=0.0,Float_t max=1.0e9)   {fADCOverClustSizeMin = min ; fADCOverClustSizeMax = max;  }
  void SetMaxADCOverClustADCRange(Float_t min=0.0,Float_t max=1.0e9) {fMaxADCOverClustADCMin = min ; fMaxADCOverClustADCMax = max;  }
  void SetTimeRange(Float_t min=0.0, Float_t max=1.0e9)              {fTimeMin = min ; fTimeMax = max; }
  void SetClustSizeRange(Float_t min=0.0, Float_t max=1.0e9)         {fClustSizeMin = min ; fClustSizeMax = max; }

  void SetTimebinRmsMin(Float_t iroc=0.0,Float_t oroc=0.0)           {fTimebinRmsIrocMin = iroc ; fTimebinRmsOrocMin = oroc; }
  void SetPadRmsMin(Float_t iroc=0.0,Float_t oroc=0.0)               {fPadRmsIrocMin = iroc ; fPadRmsOrocMin = oroc; }
  void SetRowRmsMin(Float_t iroc=0.0,Float_t oroc=0.0)               {fRowRmsIrocMin = iroc ; fRowRmsOrocMin = oroc; }
  void SetClusterPadSize1DMax(Short_t iroc=200,Short_t oroc=200) {fClusterPadSize1DIrocMax = iroc ; fClusterPadSize1DOrocMax = oroc; }
  void SetCurveCoefficient(Float_t iroc=1.0e9,Float_t oroc=1.0e9)    {fCurveCoefficientIroc = iroc ; fCurveCoefficientOroc = oroc; }

  void SetIrocHistogram(Int_t nbins=200,Float_t min=100,Float_t max=6000) {fIrocHistogramNbins = nbins ; fIrocHistogramMin = min ; fIrocHistogramMax = max; }
  void SetOrocHistogram(Int_t nbins=200,Float_t min=100,Float_t max=5500) {fOrocHistogramNbins = nbins ; fOrocHistogramMin = min ; fOrocHistogramMax = max; }

  void SetRadius(UInt_t row=0, UInt_t pad=0) {fRowRadius = row ; fPadRadius = pad; }
  void SetStep(UInt_t row=1, UInt_t pad=1) {fRowStep = (row>=1?row:1) ; fPadStep = (pad>=1?pad:1) ; }

private:

  Bool_t fASide;              //! Only A side
  Bool_t fCSide;              //! Only C side 
  TObjArray fHistoKrArray;    //  Calibration histograms for Kr distribution

  Float_t fADCOverClustSizeMin; // min ADCcluster over Cluster size ratio
  Float_t fADCOverClustSizeMax; // max ADCcluster over Cluster size ratio
  Float_t fMaxADCOverClustADCMin; // min MaxADC over ADCcluster ratio 
  Float_t fMaxADCOverClustADCMax; // max MaxADC over ADCcluster ratio
  Float_t fTimeMin; // min time bin for MaxADC
  Float_t fTimeMax; // max time bin for MaxADC
  Float_t fClustSizeMin; // min cluster size
  Float_t fClustSizeMax; // max cluster size

  Float_t fTimebinRmsIrocMin; // min Timebin RMS for IROCs
  Float_t fPadRmsIrocMin; // min Pad RMS for IROCs
  Float_t fRowRmsIrocMin; // min Row RMS for IROCs
  Short_t fClusterPadSize1DIrocMax; // max size of cluster in pad dir. for IROCs 
  Float_t fCurveCoefficientIroc; // A coefficient in curve function for IROCs

  Float_t fTimebinRmsOrocMin; // min Timebin RMS for OROCs
  Float_t fPadRmsOrocMin; // min Pad RMS for OROCs
  Float_t fRowRmsOrocMin; // min Row RMS for OROCs
  Short_t fClusterPadSize1DOrocMax; // max size of cluster in pad dir. for OROCs 
  Float_t fCurveCoefficientOroc; // A coefficient in curve function for OROCs

  Float_t fIrocHistogramMin; // minimal range of histogram for IROCs
  Float_t fIrocHistogramMax; // maximal range of histogram for IROCs
  Int_t   fIrocHistogramNbins; // number of bins in IROC histogram
  Float_t fOrocHistogramMin; // minimal range of histogram for OROCs
  Float_t fOrocHistogramMax; // maximal range of histogram for OROCs
  Int_t   fOrocHistogramNbins; // number of bins in OROC histogram

  UInt_t fRowRadius; // window size around pad +/-; set to 0 for pad-by-pad calib
  UInt_t fPadRadius; // window size around pad +/-; set to 0 for pad-by-pad calib
  UInt_t fRowStep; // step size; set to 1 for finest granularity
  UInt_t fPadStep; // step size; set to 1 for finest granularity



public:
  ClassDef(AliTPCCalibKr, 4)  // Implementation of the TPC krypton calibration
};

#endif

