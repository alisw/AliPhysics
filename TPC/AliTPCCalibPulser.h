#ifndef ALITPCCALIBPULSER_H
#define ALITPCCALIBPULSER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TVectorT.h>
#include <TObject.h>
#include <TObjArray.h>

class TObjArray;
class TH2S;
class TTreeSRedirector;
class AliTPCCalPad;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCParam;
class AliRawReader;
class AliTPCRawStream;
struct eventHeaderStruct;

class AliTPCCalibPulser : public TObject {

public:
    AliTPCCalibPulser();
    AliTPCCalibPulser(const AliTPCCalibPulser &sig);
    virtual ~AliTPCCalibPulser();

    AliTPCCalibPulser& operator = (const  AliTPCCalibPulser &source);


    Bool_t ProcessEvent(AliTPCRawStream *rawStream);
    Bool_t ProcessEvent(AliRawReader    *rawReader);
    Bool_t ProcessEvent(eventHeaderStruct   *event);

    Int_t Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
	       const Int_t iTimeBin, const Float_t signal);
    void Analyse();
    //
    AliTPCCalROC* GetCalRocT0 (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
    AliTPCCalROC* GetCalRocQ  (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
    AliTPCCalROC* GetCalRocRMS(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
    AliTPCCalROC* GetCalRocOutliers(Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector

    const TObjArray* GetCalPadT0() { return &fCalRocArrayT0; }      // get calibration object
    const TObjArray* GetCalPadQ()  { return &fCalRocArrayQ;  }      // get calibration object
    const TObjArray* GetCalPadRMS(){ return &fCalRocArrayRMS;}      // get calibration object
    const TObjArray* GetCalPadOutliers(){ return &fCalRocArrayOutliers;}      // get calibration object

    TH2S* GetHistoQ  (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
    TH2S* GetHistoT0 (Int_t sector, Bool_t force=kFALSE);           // get refernce histogram
    TH2S* GetHistoRMS(Int_t sector, Bool_t force=kFALSE);           // get refernce histogram

    Short_t GetDebugLevel()     const { return fDebugLevel;    }
    //
    void  SetRangeTime (Int_t firstTimeBin, Int_t lastTimeBin) { fFirstTimeBin=firstTimeBin;   fLastTimeBin=lastTimeBin;  } //Set range in which the pulser signal is expected
    //
    void  SetRangeRefQ  (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsQ   = nBins; fXminQ   = xMin; fXmaxQ   = xMax; }   //Set range for Q reference histograms
    void  SetRangeRefT0 (Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsT0  = nBins; fXminT0  = xMin; fXmaxT0  = xMax; }   //Set range for T0 reference histograms
    void  SetRangeRefRMS(Int_t nBins, Float_t xMin, Float_t xMax){ fNbinsRMS = nBins; fXminRMS = xMin; fXmaxRMS = xMax; }   //Set range for T0 reference histograms

    void  SetOldRCUformat(Bool_t format=kTRUE){ fOldRCUformat = format; }

    void  SetDebugLevel(Short_t debug=1){ fDebugLevel = debug;}

    void  SetPedestalDatabase(AliTPCCalPad *pedestalTPC, AliTPCCalPad *padNoiseTPC) {fPedestalTPC = pedestalTPC; fPadNoiseTPC = padNoiseTPC;}

    Int_t GetFirstTimeBin()   const { return fFirstTimeBin;  }
    Int_t GetLastTimeBin()    const { return fLastTimeBin;   }

    void Merge(AliTPCCalibPulser *sig);

    void DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
    //
    // Test functions
    TObjArray* TestBinning();

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

    AliTPCCalPad *fPedestalTPC;       //! Pedestal Information
    AliTPCCalPad *fPadNoiseTPC;       //! Pad noise Information whole TPC
    AliTPCCalROC *fPedestalROC;       //! Pedestal Information for current ROC
    AliTPCCalROC *fPadNoiseROC;       //! Pad noise Information for current ROC
//    Bool_t fBpedestal;                //! are we running with pedestal substraction


    TObjArray fCalRocArrayT0;         //  Array of AliTPCCalROC class for Time0 calibration
    TObjArray fCalRocArrayQ;          //  Array of AliTPCCalROC class for Charge calibration
    TObjArray fCalRocArrayRMS;        //  Array of AliTPCCalROC class for signal width calibration
    TObjArray fCalRocArrayOutliers;  //  Array of AliTPCCalROC class for signal outliers

    TObjArray fHistoQArray;           //  Calibration histograms for Charge distribution
    TObjArray fHistoT0Array;          //  Calibration histograms for Time0  distribution
    TObjArray fHistoRMSArray;         //  Calibration histograms for signal width distribution

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

    TVectorF  fVTime0Offset;          //!  Time0 Offset from preprocessing for each sector;
    TVectorF  fVTime0OffsetCounter;   //!  Time0 Offset from preprocessing for each sector;

    //debugging
    Int_t fEvent;
    TTreeSRedirector *fDebugStreamer;  //! debug streamer

    Short_t fDebugLevel;
    //! debugging

    void   FindPedestal(Float_t part=.6);
    void FindPulserSignal(TVectorD &param, Float_t &qSum);

    TH2S* GetHisto(Int_t sector, TObjArray *arr,
		   Int_t nbinsY, Float_t ymin, Float_t ymax,
		   Char_t *type, Bool_t force);


    AliTPCCalROC* GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force);

    TVectorF* GetPadTimesEvent(Int_t sector, Bool_t force=kFALSE);

    void ResetEvent();
    void ResetPad();
    void ProcessPad();
    void EndEvent();


    //debug
    TVectorF* GetPadInfoEvent(Int_t sector, TObjArray *arr, Bool_t force=kFALSE);
    TVectorF* GetPadQEvent(Int_t sector, Bool_t force=kFALSE);
    TVectorF* GetPadRMSEvent(Int_t sector, Bool_t force=kFALSE);
    TVectorF* GetPadPedestalEvent(Int_t sector, Bool_t force=kFALSE);


public:


  ClassDef(AliTPCCalibPulser,1)           //Implementation of the TPC pulser calibration
};



#endif

