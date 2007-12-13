#ifndef ALITPCDATAQA_H
#define ALITPCDATAQA_H
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
class AliTPCCalPad; 
struct eventHeaderStruct;

class AliTPCdataQA : public TObject {

public:
  AliTPCdataQA();
  AliTPCdataQA(const AliTPCdataQA &ped);
  virtual ~AliTPCdataQA();

  AliTPCdataQA& operator = (const  AliTPCdataQA &source);
 void  DumpToFile(const Char_t *filename, const Char_t *dir="", const Bool_t append=kFALSE);
  //
  Bool_t ProcessEventFast(AliTPCRawStreamFast *rawStreamFast);
  Bool_t ProcessEventFast(AliRawReader        *rawReader);
  Bool_t ProcessEvent(AliTPCRawStream *rawStream);
  Bool_t ProcessEvent(AliRawReader    *rawReader);
  Bool_t ProcessEvent(eventHeaderStruct   *event);

  Int_t  Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
	        const Int_t iTimeBin, const Float_t signal);
  void   Analyse();
  //
  //
  AliTPCCalPad *GetMaxCharge(){ return fMaxCharge;}
  AliTPCCalPad *GetOverThreshold0(){ return fOverThreshold0;}
  AliTPCCalPad *GetOverThreshold5(){ return fOverThreshold5;}
  AliTPCCalPad *GetOverThreshold10(){ return fOverThreshold10;}
  AliTPCCalPad *GetOverThreshold20(){ return fOverThreshold20;}
  AliTPCCalPad *GetOverThreshold30(){ return fOverThreshold30;}

  //
  AliTPCAltroMapping **GetAltroMapping() { return fMapping; };
  void  SetAltroMapping(AliTPCAltroMapping **mapp) { fMapping = mapp; };
  //
  //
  Int_t GetFirstTimeBin() const { return fFirstTimeBin; }
  Int_t GetLastTimeBin()  const { return fLastTimeBin;  }
  Int_t GetAdcMin()       const { return fAdcMin;       }
  Int_t GetAdcMax()       const { return fAdcMax;       }
  void  SetRangeTime(Int_t tMin, Int_t tMax){ fFirstTimeBin=tMin; fLastTimeBin=tMax; }  // Set time bin range that is used for the pedestal calibration
  void  SetRangeAdc (Int_t aMin, Int_t aMax){ fAdcMin=aMin; fAdcMax=aMax; }  // Set adc range for the pedestal calibration

  void  SetOldRCUformat(Bool_t format=kTRUE) { fOldRCUformat = format; }
  

private:

  Int_t fFirstTimeBin;              //  First Time bin needed for analysis
  Int_t fLastTimeBin;               //  Last Time bin needed for analysis
  Int_t fAdcMin;                    //  min adc channel of pedestal value
  Int_t fAdcMax;                    //  max adc channel of pedestal value
  Bool_t  fOldRCUformat;            //! Should we use the old RCU format for data reading

  AliTPCROC *fROC;                  //! ROC information
  AliTPCAltroMapping **fMapping;    //! Altro Mapping object
  //
  //
  AliTPCCalPad * fMaxCharge;       // max charge
  AliTPCCalPad * fOverThreshold0;   // number of digits over threshold
  AliTPCCalPad * fOverThreshold5;   // number of digits over threshold
  AliTPCCalPad * fOverThreshold10;   // number of digits over threshold
  AliTPCCalPad * fOverThreshold20;   // number of digits over threshold
  AliTPCCalPad * fOverThreshold30;   // number of digits over threshold
  Int_t          fEventCounter;      // event Counter

public:
  ClassDef(AliTPCdataQA, 1)  // Implementation of the TPC pedestal and noise calibration
};



#endif

