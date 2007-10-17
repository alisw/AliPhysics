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
class TObjArray;
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
class TGraph;
struct eventHeaderStruct;

class AliTPCCalibCE : public TObject {

public:
    AliTPCCalibCE();
    AliTPCCalibCE(const AliTPCCalibCE &sig);
    virtual ~AliTPCCalibCE();

    AliTPCCalibCE& operator = (const  AliTPCCalibCE &source);


    Bool_t ProcessEvent(AliTPCRawStream *rawStream);
    Bool_t ProcessEvent(AliRawReader    *rawReader);
    Bool_t ProcessEvent(eventHeaderStruct   *event);

    Int_t Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
	       const Int_t iTimeBin, const Float_t signal);
    void Analyse();
    //
    AliTPCCalROC* GetCalRocT0  (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
    AliTPCCalROC* GetCalRocQ   (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
    AliTPCCalROC* GetCalRocRMS(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
    AliTPCCalROC* GetCalRocOutliers(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector

    const TObjArray* GetCalPadT0()  const { return &fCalRocArrayT0; }      // get calibration object
    const TObjArray* GetCalPadQ()   const { return &fCalRocArrayQ;  }      // get calibration object
    const TObjArray* GetCalPadRMS() const { return &fCalRocArrayRMS;}      // get calibration object
    const TObjArray* GetCalPadOutliers() const { return &fCalRocArrayOutliers;}      // get calibration object

    TH2S* GetHistoQ  (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
    TH2S* GetHistoT0 (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
    TH2S* GetHistoRMS(Int_t sector, Bool_t force=kFALSE);           // get refernce histogram

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

    Short_t GetDebugLevel()     const { return fDebugLevel;    }
    //
    void  SetRangeTime (Int_t firstTimeBin, Int_t lastTimeBin) { fFirstTimeBin=firstTimeBin;   fLastTimeBin=lastTimeBin;  } //Set range in which the pulser signal is expected
    //
    void  SetRangeRefQ  (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsQ   = nBins; fXminQ   = xMin; fXmaxQ   = xMax; }   //Set range for Q reference histograms
    void  SetRangeRefT0 (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsT0  = nBins; fXminT0  = xMin; fXmaxT0  = xMax; }   //Set range for T0 reference histograms
    void  SetRangeRefRMS(Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsRMS = nBins; fXminRMS = xMin; fXmaxRMS = xMax; }   //Set range for T0 reference histograms
    //
    void  SetTimeStampEvent(Double_t timestamp){ fTimeStamp = timestamp; }
    void  SetRunNumber(Double_t eventnumber){ fRunNumber = eventnumber; }

    void  SetEventInfo(Double_t runNumber, Double_t timestamp, Double_t eventId){ fRunNumber=runNumber; fTimeStamp=timestamp; fEventId=eventId;}

    void  SetOldRCUformat(Bool_t format=kTRUE){ fOldRCUformat = format; }

    void  SetDebugLevel(Short_t debug=1){ fDebugLevel = debug;}

    void  SetPedestalDatabase(AliTPCCalPad *pedestalTPC, AliTPCCalPad *padNoiseTPC) {fPedestalTPC = pedestalTPC; fPadNoiseTPC = padNoiseTPC;}

    Int_t GetFirstTimeBin()   const { return fFirstTimeBin;  }
    Int_t GetLastTimeBin()    const { return fLastTimeBin;   }

    Int_t GetNeventsProcessed() const { return fNevents; }

    void Merge(AliTPCCalibCE *ce);

    TGraph *MakeGraphTimeCE(Int_t sector, Int_t xVariable=0, Int_t fitType=0, Int_t fitParameter=0);

    void DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);

private:
    Int_t fFirstTimeBin;              //  First Time bin needed for analysis
    Int_t fLastTimeBin;               //  Last Time bin needed for analysis

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

    Int_t     fLastSector;            //! Last sector processed

    Bool_t  fOldRCUformat;            //! Should we use the old RCU format for data reading

    AliTPCROC   *fROC;                //! ROC information
    AliTPCParam *fParam;              //! TPC information

    AliTPCCalPad *fPedestalTPC;       //! Pedestal Information whole TPC
    AliTPCCalPad *fPadNoiseTPC;       //! Pad noise Information whole TPC
    AliTPCCalROC *fPedestalROC;       //! Pedestal Information for current ROC
    AliTPCCalROC *fPadNoiseROC;       //! Pad noise Information for current ROC

    TObjArray fCalRocArrayT0;         //  Array of AliTPCCalROC class for Time0 calibration
    TObjArray fCalRocArrayQ;          //  Array of AliTPCCalROC class for Charge calibration
    TObjArray fCalRocArrayRMS;        //  Array of AliTPCCalROC class for signal width calibration
    TObjArray fCalRocArrayOutliers;   //  Array of AliTPCCalROC class for signal outliers

    TObjArray fHistoQArray;           //  Calibration histograms for Charge distribution
    TObjArray fHistoT0Array;          //  Calibration histograms for Time0  distribution
    TObjArray fHistoRMSArray;         //  Calibration histograms for signal width distribution

    TObjArray fHistoTmean;            //! Calibration histograms of the mean CE position for all sectors

    TObjArray fParamArrayEventPol1;   //  Store mean arrival time parameters for each sector event by event from global plane fit
    TObjArray fParamArrayEventPol2;   //  Store mean arrival time parameters for each sector event by event from global parabola fit
    TObjArray fTMeanArrayEvent;       //  Store mean arrival time for each sector event by event
    TObjArray fQMeanArrayEvent;       //  Store mean arrival Charge for each sector event by event
    TVectorD  fVEventTime;            //  Timestamps of the events
    TVectorD  fVEventNumber;          //  Eventnumbers of the events
    Int_t     fNevents;               //  Event counter
    Double_t  fTimeStamp;             //! Timestamp of the current event
    Double_t  fEventId;               //! Event Id of the current event
    Double_t  fRunNumber;             //! Run Number of the current event
    Double_t  fOldRunNumber;          //! Old Run Number

    TObjArray fPadTimesArrayEvent;    //! Pad Times for the event, before mean Time0 corrections
    TObjArray fPadQArrayEvent;        //! Charge for the event, only needed for debugging streamer
    TObjArray fPadRMSArrayEvent;      //! Signal width for the event, only needed for debugging streamer
    TObjArray fPadPedestalArrayEvent; //! Signal width for the event, only needed for debugging streamer

    Int_t     fCurrentChannel;         //! current channel processed
    Int_t     fCurrentSector;          //! current sector processed
    Int_t     fCurrentRow;             //! current row processed
    Float_t   fMaxPadSignal;           //! maximum bin of current pad
    Int_t     fMaxTimeBin;             //! time bin with maximum value
    TVectorF  fPadSignal;              //! signal of current Pad
    Float_t   fPadPedestal;            //! Pedestal Value of current pad
    Float_t   fPadNoise;               //! Noise Value of current pad

    TVectorD  fVTime0Offset;          //!  Time0 Offset for each sector;
    TVectorD  fVTime0OffsetCounter;   //!  Time0 Offset counter for each sector;
    TVectorD  fVMeanQ;                //!  Mean Q for each sector;
    TVectorD  fVMeanQCounter;         //!  Mean Q counter for each sector;
    //debugging
//    Int_t fEvent;
    TTreeSRedirector *fDebugStreamer;  //! debug streamer

    Short_t fDebugLevel;              // debug level
    //! debugging

    void   FindPedestal(Float_t part=.6);
    void   FindCESignal(TVectorD &param, Float_t &qSum, const TVectorF maxima);
    void   FindLocalMaxima(TVectorF &maxima);
    Bool_t IsPeak(Int_t pos, Int_t tminus, Int_t tplus) const;

    TH2S* GetHisto(Int_t sector, TObjArray *arr,
		   Int_t nbinsY, Float_t ymin, Float_t ymax,
		   Char_t *type, Bool_t force);
    TH1S* GetHisto(Int_t sector, TObjArray *arr,
		   Char_t *type, Bool_t force);

    AliTPCCalROC* GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) const;

    TVectorF* GetVectSector(Int_t sector, TObjArray *arr, UInt_t size, Bool_t force=kFALSE) const;
    TVectorF* GetPadTimesEvent(Int_t sector, Bool_t force=kFALSE);

    TObjArray* GetParamArray(Int_t sector, TObjArray *arr, Bool_t force=kFALSE) const;

    void ResetEvent();
    void ResetPad();
    void ProcessPad();
    void EndEvent();


    //debug
    TVectorF* GetPadQEvent(Int_t sector, Bool_t force=kFALSE);
    TVectorF* GetPadRMSEvent(Int_t sector, Bool_t force=kFALSE);
    TVectorF* GetPadPedestalEvent(Int_t sector, Bool_t force=kFALSE);

    ClassDef(AliTPCCalibCE,2)  //Implementation of the TPC Central Electrode calibration

};



#endif

