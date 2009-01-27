#ifndef ALITPCCALIBCE_H
#define ALITPCCALIBCE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//             Implementation of the TPC Central Electrode calibration                //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////

#include <TVectorT.h>
#include "AliTPCCalibRawBase.h"
class TH1S;
class TH2S;
class TH1F;
class TTreeSRedirector;
class AliTPCCalPad;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCParam;
class AliRawReader;
class AliTPCRawStream;
class AliTPCRawStreamFast;
class TGraph;
class TMap;

struct eventHeaderStruct;

class AliTPCCalibCE : public AliTPCCalibRawBase {
  
public:
  AliTPCCalibCE();
  AliTPCCalibCE(const AliTPCCalibCE &sig);
  AliTPCCalibCE(const TMap *config);
  virtual ~AliTPCCalibCE();
  
  AliTPCCalibCE& operator = (const  AliTPCCalibCE &source);
  
  virtual Int_t Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
                       const Int_t iTimeBin, const Float_t signal);
  virtual void Analyse();
  
    //
  AliTPCCalROC* GetCalRocT0  (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocT0Err(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocQ   (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocRMS(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocOutliers(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  
  const TObjArray* GetCalPadT0()    const { return &fCalRocArrayT0; }      // get calibration object
  const TObjArray* GetCalPadT0Err() const { return &fCalRocArrayT0Err; }      // get calibration object
  const TObjArray* GetCalPadQ()     const { return &fCalRocArrayQ;  }      // get calibration object
  const TObjArray* GetCalPadRMS()   const { return &fCalRocArrayRMS;}      // get calibration object
  const TObjArray* GetCalPadOutliers() const { return &fCalRocArrayOutliers;}      // get calibration object
  
  TH2S* GetHistoQ  (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  TH2S* GetHistoT0 (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  TH2S* GetHistoRMS(Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  
  Float_t GetMeanT0rms() const {return fMeanT0rms;}
  Float_t GetMeanQrms() const {return fMeanQrms;}
  Float_t GetMeanRMSrms() const {return fMeanRMSrms;}
  
  Int_t   GetPeakDetectionMinus() const {return fPeakDetMinus;}
  Int_t   GetPeakDetectionPlus()  const {return fPeakDetPlus;}
  Int_t   GetPeakIntRangeMinus() const {return fPeakIntMinus;}
  Int_t   GetPeakIntRangePlus()  const {return fPeakIntPlus;}
  Float_t GetNnoiseThresholdMax() const {return fNoiseThresholdMax;}
  Float_t GetNnoiseThresholdSum() const {return fNoiseThresholdSum;}
  
  TH1S* GetHistoTmean(Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
  
    //needed here to merge ClibCE objects
  TObjArray* GetParamArrayPol1(Int_t sector, Bool_t force=kFALSE);
  TObjArray* GetParamArrayPol2(Int_t sector, Bool_t force=kFALSE);
  
//    TObjArray*  GetTMeanArrayEvent(){ return &fTMeanArrayEvent; }
//    TObjArray*  GetQMeanArrayEvent(){ return &fQMeanArrayEvent; }
  TVectorF* GetTMeanEvents(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetQMeanEvents(Int_t sector, Bool_t force=kFALSE);
  
  TVectorD*   GetEventTimes()     { return &fVEventTime;      }
  TVectorD*   GetEventIds()       { return &fVEventNumber;    }
  
  //
  void  SetRangeRefQ  (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsQ   = nBins; fXminQ   = xMin; fXmaxQ   = xMax; }   //Set range for Q reference histograms
  void  SetRangeRefT0 (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsT0  = nBins; fXminT0  = xMin; fXmaxT0  = xMax; }   //Set range for T0 reference histograms
  void  SetRangeRefRMS(Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsRMS = nBins; fXminRMS = xMin; fXmaxRMS = xMax; }   //Set range for T0 reference histograms
  //
  void  SetRangePeakDetection(Int_t minus, Int_t plus) { fPeakDetMinus=minus; fPeakDetPlus=plus;}
  void  SetRangePeakIntegral(Int_t minus, Int_t plus) { fPeakIntMinus=minus; fPeakIntPlus=plus;}
  void  SetNnoiseThresholdMax(Float_t n) {fNoiseThresholdMax=n;}
  void  SetNnoiseThresholdSum(Float_t n) {fNoiseThresholdSum=n;}
  //
  void  SetTimeStampEvent(Double_t timestamp){ fTimeStamp = timestamp; }
  void  SetRunNumber(Double_t eventnumber){ fRunNumber = eventnumber; }
  void  SetEventInfo(Double_t runNumber, Double_t timestamp, Double_t eventId){ fRunNumber=runNumber; fTimeStamp=timestamp; fEventId=eventId;}
  //
  void  SetPedestalDatabase(AliTPCCalPad *pedestalTPC, AliTPCCalPad *padNoiseTPC) {fPedestalTPC = pedestalTPC; fPadNoiseTPC = padNoiseTPC;}
  void  SetIsZeroSuppressed(Bool_t zs=kTRUE) { fIsZeroSuppressed=zs; }
  void  SetSecRejectRatio(Float_t ratio) { fSecRejectRatio=ratio; }
  //Getters
  Int_t GetNeventsProcessed() const { return fNevents; }
  
  Bool_t GetIsZeroSuppressed() const { return fIsZeroSuppressed; }
  
  Float_t  GetSecRejectRatio() const { return fSecRejectRatio; }

  const TVectorF *GetTime0Side(Int_t side=0) const {return (side==0)?&fVTime0SideA:&fVTime0SideC;}
  Float_t GetPeakIntegralMinus() const {return fPeakIntMinus;}
  Float_t GetPeakIntegralPlus() const {return fPeakIntPlus;}
  
  
  void Merge(AliTPCCalibCE *ce);
  
  TGraph *MakeGraphTimeCE(Int_t sector, Int_t xVariable=0, Int_t fitType=0, Int_t fitParameter=0);
    
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
  Int_t   fPeakDetMinus;               //  Consecutive timebins on rising edge to be regarded as a signal
  Int_t   fPeakDetPlus;                //  Consecutive timebins on falling edge to be regarded as a signal
  Int_t   fPeakIntMinus;            //  Peak integral range for COG determination. Bins used before max bin
  Int_t   fPeakIntPlus;             //  Peak integral range for COG determination. Bins used after max bin
  Float_t fNoiseThresholdMax;       //  Analysis Treshold for signal finding: Max>fNoiseThresholdMax*PadNoise
  Float_t fNoiseThresholdSum;       //  Analysis Treshold for signal finding: Sum>fNoiseThresholdSum*PadNoise
  
  Bool_t  fIsZeroSuppressed;        //  If data is Zero Suppressed -> Don't subtrakt pedestals!
  
  Int_t     fLastSector;            //! Last sector processed
  
  Float_t   fSecRejectRatio;        //! Needed percentage of signals in one chamber. Below it will be rejected
                                      //  This is neede if we do not process a laser event
  
  AliTPCParam *fParam;              //! TPC information
  
  AliTPCCalPad *fPedestalTPC;       //! Pedestal Information whole TPC
  AliTPCCalPad *fPadNoiseTPC;       //! Pad noise Information whole TPC
  AliTPCCalROC *fPedestalROC;       //! Pedestal Information for current ROC
  AliTPCCalROC *fPadNoiseROC;       //! Pad noise Information for current ROC
  
  TObjArray fCalRocArrayT0;         //  Array of AliTPCCalROC class for Time0 calibration
  TObjArray fCalRocArrayT0Err;      //  Array of AliTPCCalROC class for the error (rms) of Time0 calibration
  TObjArray fCalRocArrayQ;          //  Array of AliTPCCalROC class for Charge calibration
  TObjArray fCalRocArrayRMS;        //  Array of AliTPCCalROC class for signal width calibration
  TObjArray fCalRocArrayOutliers;   //  Array of AliTPCCalROC class for signal outliers
  
  TObjArray fHistoQArray;           //  Calibration histograms for Charge distribution
  TObjArray fHistoT0Array;          //  Calibration histograms for Time0  distribution
  TObjArray fHistoRMSArray;         //  Calibration histograms for signal width distribution
  
  Float_t   fMeanT0rms;             // mean of the rms of all pad T0  fits, used as error estimation of T0 results
  Float_t   fMeanQrms;              // mean of the rms of all pad Q   fits, used as error estimation of Q results
  Float_t   fMeanRMSrms;            // mean of the rms of all pad TMS fits, used as error estimation of RMS results
  
  TObjArray fHistoTmean;            //! Calibration histograms of the mean CE position for all sectors
  
  TObjArray fParamArrayEventPol1;   //  Store mean arrival time parameters for each sector event by event from global plane fit
  TObjArray fParamArrayEventPol2;   //  Store mean arrival time parameters for each sector event by event from global parabola fit
  TObjArray fTMeanArrayEvent;       //  Store mean arrival time for each sector event by event
  TObjArray fQMeanArrayEvent;       //  Store mean arrival Charge for each sector event by event
  TVectorD  fVEventTime;            //  Timestamps of the events
  TVectorD  fVEventNumber;          //  Eventnumbers of the events
  TVectorF  fVTime0SideA;           //  Mean Time0 for side A for all events
  TVectorF  fVTime0SideC;           //  Mean Time0 for side C for all events
  Double_t  fTimeStamp;             //! Timestamp of the current event
  Double_t  fEventId;               //! Event Id of the current event
  Double_t  fRunNumber;             //! Run Number of the current event
  Double_t  fOldRunNumber;          //! Old Run Number
  
  TObjArray fPadTimesArrayEvent;    //! Pad Times for the event, before mean Time0 corrections
  TObjArray fPadQArrayEvent;        //! Charge for the event, only needed for debugging streamer
  TObjArray fPadRMSArrayEvent;      //! Signal width for the event, only needed for debugging streamer
  TObjArray fPadPedestalArrayEvent; //! Signal width for the event, only needed for debugging streamer
  
  Int_t     fCurrentChannel;        //! current channel processed
  Int_t     fCurrentSector;         //! current sector processed
  Int_t     fCurrentRow;            //! current row processed
  Float_t   fMaxPadSignal;          //! maximum bin of current pad
  Int_t     fMaxTimeBin;            //! time bin with maximum value
  TVectorF  fPadSignal;             //! signal of current Pad
  Float_t   fPadPedestal;           //! Pedestal Value of current pad
  Float_t   fPadNoise;              //! Noise Value of current pad
  
  TVectorD  fVTime0Offset;          //!  Time0 Offset for each sector;
  TVectorD  fVTime0OffsetCounter;   //!  Time0 Offset counter for each sector;
  TVectorD  fVMeanQ;                //!  Mean Q for each sector;
  TVectorD  fVMeanQCounter;         //!  Mean Q counter for each sector;
  
  Float_t   fCurrentCETimeRef;      //! Time refernce of the current sector
  
  void   FindPedestal(Float_t part=.6);
  void   UpdateCETimeRef(); //Get the time reference of the last valid measurement in sector
  void   FindCESignal(TVectorD &param, Float_t &qSum, const TVectorF maxima);
  void   FindLocalMaxima(TVectorF &maxima);
  Bool_t IsPeak(Int_t pos, Int_t tminus, Int_t tplus) const;
  
  TH2S* GetHisto(Int_t sector, TObjArray *arr,
                 Int_t nbinsY, Float_t ymin, Float_t ymax,
                 const Char_t *type, Bool_t force);
  TH1S* GetHisto(Int_t sector, TObjArray *arr,
                 const Char_t *type, Bool_t force);
  
  AliTPCCalROC* GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) const;
  
  TVectorF* GetVectSector(Int_t sector, TObjArray *arr, UInt_t size, Bool_t force=kFALSE) const;
  TVectorF* GetPadTimesEvent(Int_t sector, Bool_t force=kFALSE);
  
  TObjArray* GetParamArray(Int_t sector, TObjArray *arr, Bool_t force=kFALSE) const;
  
  virtual void EndEvent();
  virtual void ResetEvent();
  void ResetPad();
  void ProcessPad();
  //debug
  TVectorF* GetPadQEvent(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetPadRMSEvent(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetPadPedestalEvent(Int_t sector, Bool_t force=kFALSE);
  
  ClassDef(AliTPCCalibCE,8)  //Implementation of the TPC Central Electrode calibration
};

#endif
