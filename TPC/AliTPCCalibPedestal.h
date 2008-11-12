#ifndef ALITPCCALIBPEDESTAL_H
#define ALITPCCALIBPEDESTAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include <TObject.h>
#include <TObjArray.h>

class TArrayF;
class TH2F;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCRawStream;
class AliTPCRawStreamFast;
class AliRawReader;
class AliTPCAltroMapping;
class TMap;

struct eventHeaderStruct;

class AliTPCCalibPedestal : public TObject {

public:
  AliTPCCalibPedestal();
  AliTPCCalibPedestal(const AliTPCCalibPedestal &ped);
  AliTPCCalibPedestal(const TMap *config);
  virtual ~AliTPCCalibPedestal();

  AliTPCCalibPedestal& operator = (const  AliTPCCalibPedestal &source);

  Bool_t ProcessEventFast(AliTPCRawStreamFast *rawStreamFast);
  Bool_t ProcessEventFast(AliRawReader        *rawReader);

  Bool_t ProcessEvent(AliTPCRawStream *rawStream);
  Bool_t ProcessEvent(AliRawReader    *rawReader);
  Bool_t ProcessEvent(eventHeaderStruct   *event);

  Int_t  Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
	        const Int_t iTimeBin, const Float_t signal);
  void   Analyse();
  //
  AliTPCAltroMapping **GetAltroMapping() { return fMapping; };
  void  SetAltroMapping(AliTPCAltroMapping **mapp) { fMapping = mapp; };
  //
  AliTPCCalROC* GetCalRocPedestal (Int_t sector, Bool_t force=kFALSE);  // get calibration object - sector
  AliTPCCalROC* GetCalRocSigma(Int_t sector, Bool_t force=kFALSE);        // get calibration object - sector
  const TObjArray* GetCalPadPedestal (){return &fCalRocArrayPedestal;}  // get calibration object
  const TObjArray* GetCalPadSigma(){return &fCalRocArraySigma;}             // get calibration object

  AliTPCCalROC* GetCalRocMean (Int_t sector, Bool_t force=kFALSE);      // get calibration object - sector
  AliTPCCalROC* GetCalRocRMS(Int_t sector, Bool_t force=kFALSE);        // get calibration object - sector
  const TObjArray* GetCalPadMean (){return &fCalRocArrayMean;}          // get calibration object
  const TObjArray* GetCalPadRMS(){return &fCalRocArrayRMS;}             // get calibration object
  
  TH2F* GetHistoPedestal  (Int_t sector, Bool_t force=kFALSE);          // get refernce histogram
  void  DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
  //
  void  SetTimeAnalysis(Bool_t time = kTRUE);                  // Use ONLY in TPCPEDESTALda on LDC for one sector!
  void  AnalyseTime(Int_t nevents);                            // Makes sense only in TPCPEDESTALda on LDC!
  TArrayF **GetTimePedestals()  const { return fTimeSignal; }  // Get array with time dependent pedestals (for one sector!)
  //
  Int_t   GetFirstTimeBin() const { return fFirstTimeBin; }
  Int_t   GetLastTimeBin()  const { return fLastTimeBin;  }
  Int_t   GetAdcMin()       const { return fAdcMin;       }
  Int_t   GetAdcMax()       const { return fAdcMax;       }
  Float_t GetAnaMeanDown()  const { return fAnaMeanDown;  }
  Float_t GetAnaMeanUp()    const { return fAnaMeanUp;    }
  
  void  SetRangeTime(Int_t tMin, Int_t tMax){ fFirstTimeBin=tMin; fLastTimeBin=tMax; }  // Set time bin range that is used for the pedestal calibration
  void  SetRangeAdc (Int_t aMin, Int_t aMax){ fAdcMin=aMin; fAdcMax=aMax; }  // Set adc range for the pedestal calibration
  void  SetAnalysisTruncationRange(Float_t down, Float_t up) {fAnaMeanDown=down; fAnaMeanUp=up;}    //Set range for truncated mean analysis of the channel information

  void  Merge(AliTPCCalibPedestal *ped);

  Bool_t TestEvent();  // Test the fast approach to fill histogram - used for test purposes

private:

  Int_t fFirstTimeBin;              //  First Time bin needed for analysis
  Int_t fLastTimeBin;               //  Last Time bin needed for analysis

  Int_t fAdcMin;                    //  min adc channel of pedestal value
  Int_t fAdcMax;                    //  max adc channel of pedestal value

  Float_t fAnaMeanDown;             // Truncated mean channel analysis - lower cut
  Float_t fAnaMeanUp;               // Truncated mean channel analysis - upper cut
  
  Bool_t  fTimeAnalysis;            //! Should we use the time dependent analysis? ONLY ON LDC!

  AliTPCROC *fROC;                  //! ROC information
  AliTPCAltroMapping **fMapping;    //! Altro Mapping object

  TObjArray fCalRocArrayPedestal;   //  Array of AliTPCCalROC class for pedestal values from gaus fit
  TObjArray fCalRocArraySigma;      //  Array of AliTPCCalROC class for noise values from gaus fit

  TObjArray fHistoPedestalArray;    //  Calibration histograms for Pedestal distribution

  TArrayF **fTimeSignal;            //! Arrays which hold time dependent signals
  
  TObjArray fCalRocArrayMean;       //  Array of AliTPCCalROC class for pedestal values, simple mean
  TObjArray fCalRocArrayRMS;        //  Array of AliTPCCalROC class for noise values, simple rms

  TH2F* GetHisto(Int_t sector, TObjArray *arr,
		 Int_t nbinsY, Float_t ymin, Float_t ymax,
		 const Char_t *type, Bool_t force);

  AliTPCCalROC* GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force);

public:
  ClassDef(AliTPCCalibPedestal, 6)  // Implementation of the TPC pedestal and noise calibration
};



#endif

