#ifndef ALITPCCALIBPULSER_H
#define ALITPCCALIBPULSER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//                  Implementation of the TPC pulser calibration                       //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

#include <TVectorT.h>
#include "AliTPCCalibRawBase.h"
#include <TObjArray.h>
class TH2S;
class TH2F;
class TTreeSRedirector;
class AliTPCCalPad;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCParam;
class AliRawReader;
class AliTPCRawStream;
class TMap;

struct eventHeaderStruct;

class AliTPCCalibPulser : public AliTPCCalibRawBase {

public:
  AliTPCCalibPulser();
  AliTPCCalibPulser(const AliTPCCalibPulser &sig);
  AliTPCCalibPulser(const TMap *config);
  virtual ~AliTPCCalibPulser();
  
  void Reset();
  
  AliTPCCalibPulser& operator = (const  AliTPCCalibPulser &source);
  
  
  virtual Int_t Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
                       const Int_t iTimeBin, const Float_t signal);
  virtual void Analyse();
     //
  AliTPCCalROC* GetCalRocT0 (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocQ  (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocRMS(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocOutliers(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  
  const TObjArray* GetCalPadT0()  const { return &fCalRocArrayT0; }      // get calibration object
  const TObjArray* GetCalPadQ()   const { return &fCalRocArrayQ;  }      // get calibration object
  const TObjArray* GetCalPadRMS() const{ return &fCalRocArrayRMS;}      // get calibration object
  const TObjArray* GetCalPadOutliers() const { return &fCalRocArrayOutliers;}      // get calibration object
  
  TH2S* GetHistoQ  (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  TH2S* GetHistoT0 (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  TH2S* GetHistoRMS(Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  
  TH2F* GetHistoTSec();                                        // mean abs time distribution histogram
  
  Float_t GetMeanTimeSector(Int_t sector) const {return fVMeanTimeSector[sector];}
  const TVectorF* GetMeanTimeSectorArray() const {return &fVMeanTimeSector;}
  
  Short_t GetDebugLevel()     const { return fDebugLevel;    }
    //
  void  SetRangeTime (Int_t firstTimeBin, Int_t lastTimeBin) { fFirstTimeBin=firstTimeBin;   fLastTimeBin=lastTimeBin;  } //Set range in which the pulser signal is expected
    //
  void  SetRangeRefQ  (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsQ   = nBins; fXminQ   = xMin; fXmaxQ   = xMax; }   //Set range for Q reference histograms
  void  SetRangeRefT0 (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsT0  = nBins; fXminT0  = xMin; fXmaxT0  = xMax; }   //Set range for T0 reference histograms
  void  SetRangeRefRMS(Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsRMS = nBins; fXminRMS = xMin; fXmaxRMS = xMax; }   //Set range for T0 reference histograms
  void  SetRangePeakIntegral(Int_t minus, Int_t plus) { fPeakIntMinus=minus; fPeakIntPlus=plus;}
  
  void  SetDebugLevel(Short_t debug=1){ fDebugLevel = debug;}
  
  void  SetIsZeroSuppressed(Bool_t zs=kTRUE){ fIsZeroSuppressed=zs;}
  
  void  SetPedestalDatabase(AliTPCCalPad * const pedestalTPC, AliTPCCalPad * const padNoiseTPC) {fPedestalTPC = pedestalTPC; fPadNoiseTPC = padNoiseTPC;}
  void  SetOutliers(AliTPCCalPad * const outliers)  {fOutliers = outliers;}
  
  Bool_t GetIsZeroSupperssed() const { return fIsZeroSuppressed; }

  Float_t GetPeakIntegralMinus() const {return fPeakIntMinus;}
  Float_t GetPeakIntegralPlus() const {return fPeakIntPlus;}
  
  void Merge(AliTPCCalibPulser * const sig);
  virtual Long64_t Merge(TCollection * const list);
  
  //
  // Test functions
  TObjArray* TestBinning();
  
protected:
  virtual void ResetEvent();
  virtual void EndEvent();
  
private:
    // reference histogram ranges
  Int_t   fNbinsT0;                 //  Number of bins for T0 reference histogram
  Float_t fXminT0;                  //  xmin   of T0 reference histogram
  Float_t fXmaxT0;                  //  xmax   of T0 reference histogram
  Int_t   fNbinsQ;                  //  Number of bins for T0 reference histogram
  Float_t fXminQ;                   //  xmin   of T0 reference histogram
  Float_t fXmaxQ;                   //  xmax   of T0 reference histogram
  Int_t   fNbinsRMS;                //  Number of bins for T0 reference histogram
  Float_t fXminRMS;                 //  xmin   of T0 reference histogram
  Float_t fXmaxRMS;                 //  xmax   of T0 reference histogram
  Int_t   fPeakIntMinus;            //  Peak integral range for COG determination. Bins used before max bin
  Int_t   fPeakIntPlus;             //  Peak integral range for COG determination. Bins used after max bin
  
  Bool_t  fIsZeroSuppressed;        //  if data is zero suppressed
  
  Int_t     fLastSector;            //! Last sector processed

  AliTPCParam *fParam;              //! TPC information
  
  AliTPCCalPad *fPedestalTPC;       //! Pedestal Information
  AliTPCCalPad *fPadNoiseTPC;       //! Pad noise Information whole TPC
  AliTPCCalPad *fOutliers;          //! Outlier information. Those will not be used for calculating the T0
  AliTPCCalROC *fPedestalROC;       //! Pedestal Information for current ROC
  AliTPCCalROC *fPadNoiseROC;       //! Pad noise Information for current ROC
  
  TObjArray fCalRocArrayT0;         //  Array of AliTPCCalROC class for Time0 calibration
  TObjArray fCalRocArrayQ;          //  Array of AliTPCCalROC class for Charge calibration
  TObjArray fCalRocArrayRMS;        //  Array of AliTPCCalROC class for signal width calibration
  TObjArray fCalRocArrayOutliers;  //  Array of AliTPCCalROC class for signal outliers
  
  TObjArray fHistoQArray;           //  Calibration histograms for Charge distribution
  TObjArray fHistoT0Array;          //  Calibration histograms for Time0  distribution
  TObjArray fHistoRMSArray;         //  Calibration histograms for signal width distribution
  
  TH2F *fHMeanTimeSector;           //  Timing distribution per sector
  TVectorF  fVMeanTimeSector;       //  Mean time per sector from analysis of fHMeanTimeSector
  
  TObjArray fPadTimesArrayEvent;    //! Pad Times for the event, before mean Time0 corrections
  TObjArray fPadQArrayEvent;        //! Charge for the event, only needed for debugging streamer
  TObjArray fPadRMSArrayEvent;      //! Signal width for the event, only needed for debugging streamer
  TObjArray fPadPedestalArrayEvent; //! Signal width for the event, only needed for debugging streamer
  
  Int_t     fCurrentChannel;         //! current channel processed
  Int_t     fCurrentSector;          //! current sector processed
  Int_t     fCurrentRow;             //! current row processed
  Int_t     fCurrentPad;             //! current pad processed
  Float_t   fMaxPadSignal;           //! maximum bin of current pad
  Int_t     fMaxTimeBin;             //! time bin with maximum value
  TVectorF  fPadSignal;              //! signal of current Pad
  Float_t   fPadPedestal;            //! Pedestal Value of current pad
  Float_t   fPadNoise;               //! Noise Value of current pad
  
  TVectorF  fVTime0Offset;          //!  Time0 Offset from preprocessing for each sector;
  TVectorF  fVTime0OffsetCounter;   //!  Time0 Offset from preprocessing for each sector;
  
  
  void   FindPedestal(Float_t part=.6);
  void FindPulserSignal(TVectorD &param, Float_t &qSum);
  
  TH2S* GetHisto(Int_t sector, TObjArray *arr,
                 Int_t nbinsY, Float_t ymin, Float_t ymax,
                 const Char_t *type, Bool_t force);
  
  
  AliTPCCalROC* GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) const;
  
  TVectorF* GetPadTimesEvent(Int_t sector, Bool_t force=kFALSE);
  
  Bool_t IsEdgePad(Int_t sector, Int_t row, Int_t pad);
  
  void ResetPad();
  void ProcessPad();
  
  
  //debug
  TVectorF* GetPadInfoEvent(Int_t sector, TObjArray *arr, Bool_t force=kFALSE);
  TVectorF* GetPadQEvent(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetPadRMSEvent(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetPadPedestalEvent(Int_t sector, Bool_t force=kFALSE);

  
  ClassDef(AliTPCCalibPulser,5)           //Implementation of the TPC pulser calibration
};



#endif

