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
#include <THnSparse.h>

#include "AliTPCCalibRawBase.h"
class TH1S;
#include "TObjArray.h"
class TH2S;
class TH1F;
class TTreeSRedirector;
class AliTPCCalPad;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCParam;
class AliRawReader;
class TGraph;
class TMap;
class TCollection;

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
  virtual void ProcessBunch(const Int_t sector, const Int_t row, const Int_t pad,
                            const Int_t length, const UInt_t startTimeBin, const UShort_t* signal);
  
  virtual void Analyse();
  void AnalyseTrack();
  
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
  
  const TVectorD*   GetEventTimes()  const   { return &fVEventTime;      }
  const TVectorD*   GetEventIds()    const   { return &fVEventNumber;    }
  
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
  void  SetEventInfo(UInt_t runNumber,UInt_t timestamp, UInt_t eventId){ fRunNumber=runNumber; fTimeStamp=timestamp; fEventId=eventId;}
  //
  void  SetPedestalDatabase(AliTPCCalPad * const pedestalTPC, AliTPCCalPad * const padNoiseTPC) {fPedestalTPC = pedestalTPC; fPadNoiseTPC = padNoiseTPC;}
  void  SetIsZeroSuppressed(Bool_t zs=kTRUE) { fIsZeroSuppressed=zs; }
  void  SetSecRejectRatio(Float_t ratio) { fSecRejectRatio=ratio; }

  void SetProcessOld(Bool_t process=kTRUE) {fProcessOld=process;}
  void SetProcessNew(Bool_t process=kTRUE) {fProcessNew=process; if (process&&!fHnDrift) CreateDVhist(); }
  //Getters
  Int_t GetNeventsProcessed() const { return fNevents; }
  
  Bool_t GetIsZeroSuppressed() const { return fIsZeroSuppressed; }
  
  Float_t  GetSecRejectRatio() const { return fSecRejectRatio; }

  const TVectorF *GetTime0Side(Int_t side=0) const {return (side==0)?&fVTime0SideA:&fVTime0SideC;}
  Float_t GetPeakIntegralMinus() const {return fPeakIntMinus;}
  Float_t GetPeakIntegralPlus() const {return fPeakIntPlus;}
  
  
  void Merge(AliTPCCalibCE * const ce);
  virtual Long64_t Merge(TCollection * const list);
  
  TGraph *MakeGraphTimeCE(Int_t sector, Int_t xVariable=0, Int_t fitType=0, Int_t fitParameter=0);

  //
  // New functions using also the laser tracks
  //
  Bool_t IsEdgePad(Int_t sector, Int_t row, Int_t pad) const;
  
  void FindLocalMaxima(TObjArray * const arrObj, Double_t timestamp, Int_t burst);
  Int_t FindLaserTrackID(Int_t sector,Int_t row, const Double_t *peakpos,Double_t &mindist, const Double_t *peakposloc, Int_t &itrackMin2);
  
  const THnSparseI *GetHnDrift() const {return fHnDrift;}
  const TObjArray& GetArrHnDrift() const {return fArrHnDrift;}
  const TVectorD&  GetTimeBursts() const {return fTimeBursts;}
  const TObjArray  *GetArrFitGraphs() const {return fArrFitGraphs;}

  virtual void DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
  
  static AliTPCCalibCE *ReadFromFile(const Char_t *filename);
  
protected:
  virtual void EndEvent();
  virtual void ResetEvent();
  
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
  Double_t  fEventId;               //! Event Id of the current event
  UInt_t  fOldRunNumber;          //! Old Run Number
  
  TObjArray fPadTimesArrayEvent;    //! Pad Times for the event, before mean Time0 corrections
  TObjArray fPadQArrayEvent;        //! Charge for the event, only needed for debugging streamer
  TObjArray fPadRMSArrayEvent;      //! Signal width for the event, only needed for debugging streamer
  TObjArray fPadPedestalArrayEvent; //! Signal width for the event, only needed for debugging streamer
  
  Int_t     fCurrentChannel;        //! current channel processed
  Int_t     fCurrentSector;         //! current sector processed
  Int_t     fCurrentRow;            //! current row processed
  Float_t   fMaxPadSignal;          //! maximum bin of current pad
  Int_t     fMaxTimeBin;            //! time bin with maximum value
  Float_t   fPadSignal[1024];       //! signal of current Pad
  Float_t   fPadPedestal;           //! Pedestal Value of current pad
  Float_t   fPadNoise;              //! Noise Value of current pad
  
  TVectorD  fVTime0Offset;          //!  Time0 Offset for each sector;
  TVectorD  fVTime0OffsetCounter;   //!  Time0 Offset counter for each sector;
  TVectorD  fVMeanQ;                //!  Mean Q for each sector;
  TVectorD  fVMeanQCounter;         //!  Mean Q counter for each sector;
  
  Float_t   fCurrentCETimeRef;      //! Time refernce of the current sector
  
  // new part of the algorithm
  Bool_t      fProcessOld;             // Whether to use the old algorithm
  Bool_t      fProcessNew;             // Whether to use the new algorithm
  Bool_t      fAnalyseNew;             //! Whether to analyse the new part of the algorithm.
                                       //In the DA this needs to be switched off, in the Preprocessor on...
  enum {kHnBinsDV=5};
  THnSparseI *fHnDrift;                //! Histogram digits for each pad and timebin for several timestamps
  TObjArray   fArrHnDrift;             // array of sparse histograms for each burst
  TVectorD    fTimeBursts;             //  time stamps of bursts
  UInt_t      fBinsLastAna[100];       // number of bin in the THnSparse during the last analysis
  UShort_t    fPeaks[14];               //! Peak position: 4 laser layers and CE
  UShort_t    fPeakWidths[14];          //! Peak window widths
  TObjArray  *fArrFitGraphs;           // Fit resut graphs for each parameter
  UInt_t      fEventInBunch;           //! event in current bunch
  
  
  //
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
  
  void ResetPad();
  void ProcessPad();

  // new part of the algorithm
  void CreateDVhist();
  
  void   FindLaserLayers();
  Bool_t IsPeakInRange(UShort_t timebin, Int_t roc) const;

  TObjArray *SetupMeasured();
  void ResetMeasured(TObjArray * const arr);
  
  void AddCEtoIdeal(TObjArray *arr);

  void CalculateDV(TObjArray * const arrIdeal, TObjArray * const arrMeasured, Int_t burst);
  Double_t SetBurstHnDrift();
  //debug
  TVectorF* GetPadQEvent(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetPadRMSEvent(Int_t sector, Bool_t force=kFALSE);
  TVectorF* GetPadPedestalEvent(Int_t sector, Bool_t force=kFALSE);
  
  ClassDef(AliTPCCalibCE,10)  //Implementation of the TPC Central Electrode calibration
};

//Inline functions
//_____________________________________________________________________
inline Bool_t AliTPCCalibCE::IsPeakInRange(UShort_t timebin, Int_t roc) const
{
  //
  // Check whether timebin is in the range of a laser layer
  //
  Int_t side=(roc/18)%2;
  Int_t add=7*side;
//   return kTRUE;
  if (fPeaks[13]<2) return kTRUE; //not determined yet
  for (Int_t i=add; i<add+7; ++i){
    if (TMath::Abs((Short_t)timebin-(Short_t)fPeaks[i])<(Short_t)fPeakWidths[i]) return kTRUE;
  }
  return kFALSE;
}

#endif
